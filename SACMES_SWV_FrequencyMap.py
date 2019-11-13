
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

                                    ########################
                                    ### Import Libraries ###
                                    ########################


#---Clear mac terminal memory---#
import os
import matplotlib
matplotlib.use('TkAgg')
os.system("clear && printf '\e[3J'")

#---Import Modules---#
import sys
import time
import datetime
try:
    import Tkinter as tk
    from Tkinter.ttk import *
    from Tkinter import *
    from Tknter import filedialog, Menu

except ImportError: # Python 3
    import tkinter as tk
    from tkinter.ttk import *
    from tkinter import *
    from tkinter import filedialog, Menu
    from tkinter.messagebox import showinfo


from matplotlib import style
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import csv
from pylab import *
from numpy import *
from scipy.interpolate import *
from scipy.integrate import simps
from scipy.signal import *
from itertools import *
from math import log10, floor
from decimal import Decimal
from operator import truediv
import threading
from threading import Thread
from queue import Queue
style.use('ggplot')

#---Filter out error warnings---#
import warnings
warnings.simplefilter('ignore', np.RankWarning)         #numpy polyfit_deg warning
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd") #RuntimeWarning


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

#-- file handle variable --#
handle_variable = ''    # default handle variable is nothing
e_var = 'single'        # default input file is 'Multichannel', or a single file containing all electrodes
PHE_method = 'Abs'      # default PHE Extraction is difference between absolute max/min

#------------------------------------------------------------#

InputFrequencies = [10,20,30,40,50,60,70,80,90,100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,550,600,650,700,750,800,850,900,950,1000]  # frequencies initially displayed in Frequency Listbox

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

                                    ########################
                                    ### Global Variables ###
                                    ########################

########################################
### Polynomial Regression Parameters ###
########################################
sg_window = 5           ### Savitzky-Golay window (in mV range), must be odd number (increase signal:noise)
sg_degree = 1           ### Savitzky-Golay polynomial degree
polyfit_deg = 15        ### degree of polynomial fit

cutoff_frequency = 50          ### frequency that separates 'low' and 'high'
                               ### frequencies for regression analysis and
                               ### smoothing manipulation

#############################
### Checkpoint Parameters ###
#############################
key = 0                 ### SkeletonKey
search_lim = 15         ### Search limit (sec)
PoisonPill = False      ### Stop Animation variable
FoundFilePath = False   ### If the user-inputted file is found
ExistVar = False        ### If Checkpoints are not met ExistVar = True
AlreadyInitiated = False    ### indicates if the user has already initiated analysis
HighAlreadyReset = False    ### If data for high frequencies has been reset
LowAlreadyReset = False      ### If data for low frequencies has been reset
analysis_complete = False    ### If analysis has completed, begin PostAnalysis

##################################
### Data Extraction Parameters ###
##################################
delimiter = 1               ### default delimiter is a space; 2 = tab
column_index = -2           ### column index for list_val.
                              # list_val = column_index + 3
                              # defauly column is the second (so index = 1)

###############
### Styling ###
###############
HUGE_FONT = ('Verdana', 18)
LARGE_FONT = ('Verdana', 11)
MEDIUM_FONT = ('Verdnana', 10)
SMALL_FONT = ('Verdana', 8)


                        ########################
                        ### Global Functions ###
                        ########################


##############################
### Retrieve the file name ###
##############################
def _retrieve_file(file, electrode, frequency):

    try:
        if e_var == 'single':
            filename = '%s%dHz_.txt' % (handle_variable, frequency)
            filename2 = '%s%dHz_%d.txt' % (handle_variable, frequency, file)
            filename3 = '%s%dHz__%d.txt' % (handle_variable, frequency, file)


        elif e_var == 'multiple':
            filename = 'E%s_%s%sHz_.txt' % (electrode,handle_variable,frequency)
            filename2 = 'E%s_%s%sHz_%d.txt' % (electrode,handle_variable,frequency,file)
            filename3 = 'E%s_%s%sHz__%d.txt' % (electrode,handle_variable,frequency,file)

        return filename, filename2, filename3

    except:
        print('\nError in _retrieve_file\n')



def ReadData(myfile, electrode):
    global delimiter

    ###############################################################
    ### Get the index value of the data depending on if the     ###
    ### electrodes are in the same .txt file or separate files  ###
    ###############################################################
    if e_var == 'single':
        list_val = (electrode*3) + column_index
    elif e_var == 'multiple':
        list_val = column_index + 3

    #####################
    ### Read the data ###
    #####################

    #---Preallocate Potential and Current lists---#
    with open(myfile,'r',encoding='utf-8') as mydata:
        variables = len(mydata.readlines())
        potentials = ['hold']*variables

        ### key: potential; value: current ##
        data_dict = {}

        currents = [0]*variables

    #---Extract data and dump into lists---#
    with open(myfile,'r',encoding='utf-8') as mydata:
        list_num = 0
        for line in mydata:
            check_split = line.split(delimiter)
            check_split = check_split[0]
            check_split = check_split.replace(',','')
            try:
                check_split = float(check_split)
                check_split = True
            except:
                check_split = False

            if check_split:
                #---Currents---#
                current_value = line.split(delimiter)
                current_value = current_value[list_val]                      # list_val is the index value of the given electrode
                current_value = current_value.replace(',','')
                current_value = float(current_value)
                current_value = current_value*1000000
                currents[list_num] = current_value

                #---Potentials---#
                potential_value = line.split(delimiter)[0]
                potential_value = potential_value.strip(',')
                potential_value = float(potential_value)
                potentials[list_num] = potential_value
                data_dict.setdefault(potential_value, []).append(current_value)
                list_num = list_num + 1

    currents = [abs(value) for value in currents]

    ### if there are 0's in the list (if the preallocation added to many)
    ### then remove them
    cut_value = 0
    for value in potentials:
        if value == 'hold':
            cut_value += 1


    if cut_value > 0:
        potentials = potentials[:-cut_value]
        currents = currents[:-cut_value]

    #######################
    ### Return the data ###
    #######################
    return potentials, currents, data_dict

#######################################
### Retrieve the column index value ###
#######################################
def _get_listval(electrode):

    if e_var == 'single':
        list_val = electrode

    elif e_var == 'multiple':
        list_val = 1

    return list_val


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#



                        ###############################################
                        ########   Graphical User Interface     #######
                        ###############################################


#############################################
### Class that contains all of the frames ###
#############################################
class MainWindow(tk.Tk):

    #--- Initialize the GUI ---#
    def __init__(self,master=None,*args, **kwargs):
        global container, Plot, frame_list, PlotValues, ShowFrames


        #tk.Tk.__init__(self, *args, **kwargs)
        self.master = master
        self.master.wm_title('SACMES')

        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)

        #--- Create a frame for the UI ---#
        container = tk.Frame(self.master,relief='flat',bd=5)
        container.grid(row=0,rowspan=11,padx=10, sticky = 'nsew')         ## container object has UI frame in column 0
        container.rowconfigure(0, weight=1)              ## and PlotContainer (visualization) in column 1
        container.columnconfigure(0, weight=1)


        #--- Raise the frame for initial UI ---#
        ShowFrames = {}                                 # Key: frame handle / Value: tk.Frame object
        frame = InputFrame(container, self.master)
        ShowFrames[InputFrame] = frame
        frame.grid(row=0, column=0, sticky = 'nsew')
        self.show_frame(InputFrame)


        self._create_toolbar()

        #--- High and Low Frequency Dictionary ---#
        HighLowList = {}

    #--- Function to visualize different frames ---#
    def _create_toolbar(self):

        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)

        #################
        ### Edit Menu ###
        #################
        editmenu = Menu(menubar, tearoff=0)
        editmenu.add_separator()
        editmenu.add_command(label="Adjust Current Extraction", command=lambda: self.extraction_adjustment_frame())
        self.delimiter_value = IntVar()
        self.delimiter_value.set(1)

        menubar.add_cascade(label="Settings", menu=editmenu)


    def extraction_adjustment_frame(self):
        global delimiter

        win = tk.Toplevel()
        win.wm_title("Data Extraction Settings")

        l = tk.Label(win, text="Differential Currents are in Column:")
        l.grid(row=0, column=0)

        self.list_val_entry = tk.Entry(win, width=7)
        self.list_val_entry.insert(END,2)
        self.list_val_entry.grid(row=1,column=0,pady=5)

        apply_list_val = ttk.Button(win, text="Apply", command=lambda: self.get_list_val())
        apply_list_val.grid(row=2, column=0,pady=3)


        l = tk.Label(win, text="Delimiter between data columns:")
        l.grid(row=3, column=0)

        self.space_delimiter = tk.Radiobutton(win, text='Space',variable = self.delimiter_value, value = 1)
        self.space_delimiter.grid(row=4,column=0,pady=5)

        self.tab_delimiter = tk.Radiobutton(win, text = 'Tab',variable = self.delimiter_value, value = 2)
        self.tab_delimiter.grid(row=5, column=0,pady=3)

        delimiter = self.delimiter_value.get()

        exit = ttk.Button(win, text="Exit", command= lambda: win.destroy())
        exit.grid(row=6, column=0,pady=3)

    def get_list_val(self):
        global column_index

        column_index = int(self.list_val_entry.get())
        column_index = column_index - 4

    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()

    def onExit(self):
        self.master.destroy()
        self.master.quit()
        quit()


#############################################################################
### Initial Input frame that is displayed when the program is initialized ###
#############################################################################

class InputFrame(tk.Frame):
    def __init__(self, parent, controller):
        global figures, SaveBox, ManipulateFrequenciesFrame

        tk.Frame.__init__(self, parent)             # initialize the frame
        self.parent = parent
        self.controller = controller

        ##############################################
        ### Pack all of the widgets into the frame ###
        ##############################################

        self.SelectFilePath = ttk.Button(self, style = 'Off.TButton', text = 'Select File Path', command = lambda: self.FindFile(parent))
        self.SelectFilePath.grid(row=0,column=0,columnspan=4)

        self.NoSelectedPath = tk.Label(self, text = 'No File Path Selected', font = MEDIUM_FONT, fg = 'red')
        self.PathWarningExists = False               # tracks the existence of a warning label

        ImportFileLabel = tk.Label(self, text = 'Import File Label', font=LARGE_FONT).grid(row=2,column=0,columnspan=2)
        self.ImportFileEntry = tk.Entry(self)
        self.ImportFileEntry.grid(row=3,column=0,columnspan=2,pady=5)
        self.ImportFileEntry.insert(END, handle_variable)

        #--- File Handle Input ---#
        HandleLabel = tk.Label(self, text='Exported File Handle:', font=LARGE_FONT)
        HandleLabel.grid(row=2,column=2,columnspan=2)
        self.filehandle = ttk.Entry(self)
        now = datetime.datetime.now()
        hour = str(now.hour)
        day = str(now.day)
        month = str(now.month)
        year = str(now.year)
        self.filehandle.insert(END, 'DataExport_%s_%s_%s.txt' % (year, month, day))
        self.filehandle.grid(row=3,column=2,columnspan=2,pady=5)
        EmptyLabel = tk.Label(self, text = '',font=LARGE_FONT).grid(row=4,rowspan=2,column=0,columnspan=10)


        #---File Limit Input---#
        numFileLabel = tk.Label(self, text='Number of Files:', font=LARGE_FONT)
        numFileLabel.grid(row=5,column=0,columnspan=2)
        self.numfiles = ttk.Entry(self, width=7)
        self.numfiles.insert(END, '1')
        self.numfiles.grid(row=6,column=0,columnspan=2)

        IntervalLabel = tk.Label(self, text='Analysis Interval:', font=LARGE_FONT)
        IntervalLabel.grid(row=5,column=2,columnspan=2)
        self.Interval = ttk.Entry(self, width=7)
        self.Interval.insert(END, '10')
        self.Interval.grid(row=6,column=2,columnspan=2)


    ##################################
    ### Select and Edit Electrodes ###
    ##################################

        self.ElectrodeListboxFrame = tk.Frame(self)                   # create a frame to pack in the Electrode box and
        self.ElectrodeListboxFrame.grid(row=10,column=0,columnspan=2,padx=10,pady=10, sticky = 'nsew')

        #--- parameters for handling resize ---#
        self.ElectrodeListboxFrame.rowconfigure(0, weight=1)
        self.ElectrodeListboxFrame.rowconfigure(1, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(0, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(1, weight=1)

        electrodes = [1,2,3,4,5,6,7,8]
        self.ElectrodeListExists = False
        self.ElectrodeLabel = tk.Label(self.ElectrodeListboxFrame, text='Select Electrodes:', font=LARGE_FONT)
        self.ElectrodeLabel.grid(row=0,column=0, sticky = 'nswe')
        self.ElectrodeCount = Listbox(self.ElectrodeListboxFrame, relief='groove', exportselection=0, width=10, font=LARGE_FONT, height=6, selectmode = 'multiple', bd=3)
        self.ElectrodeCount.bind('<<ListboxSelect>>',self.ElectrodeCurSelect)
        self.ElectrodeCount.grid(row=1,column=0,columnspan=2,sticky='nswe')
        for electrode in electrodes:
            self.ElectrodeCount.insert(END, electrode)

        self.scrollbar = Scrollbar(self.ElectrodeListboxFrame, orient="vertical")
        self.scrollbar.config(width=10,command=self.ElectrodeCount.yview)
        self.scrollbar.grid(row=1,column=1,sticky='nse')
        self.ElectrodeCount.config(yscrollcommand=self.scrollbar.set)

        #--- Option to have data for all electrodes in a single file ---#
        self.SingleElectrodeFile = ttk.Button(self.ElectrodeListboxFrame, text='Multipotentiostat', style = 'On.TButton', command =  lambda: self.ElectrodeSelect('Multipotentiostat'))
        self.SingleElectrodeFile.grid(row=2,column=0)

        #--- Option to have data for each electrode in a separate file ---#
        self.MultipleElectrodeFiles = ttk.Button(self.ElectrodeListboxFrame,  text='Multiplex', style = 'Off.TButton',command = lambda: self.ElectrodeSelect('Multiplex'))
        self.MultipleElectrodeFiles.grid(row=2,column=1)


    #####################################################
    ### Select and Edit Frequencies for Data Analysis ###
    #####################################################

        self.ListboxFrame = tk.Frame(self)                   # create a frame to pack in the frequency box and scrollbar
        self.ListboxFrame.grid(row=10,column=2,columnspan=2,padx=10,pady=10, sticky='nsew')
        frequencies = InputFrequencies

        #-- resize ---#
        self.ListboxFrame.rowconfigure(0, weight=1)
        self.ListboxFrame.rowconfigure(1, weight=1)
        self.ListboxFrame.columnconfigure(0, weight=1)

        self.FrequencyLabel = tk.Label(self.ListboxFrame, text='Select Frequencies', font= LARGE_FONT)
        self.FrequencyLabel.grid(row=0,padx=10)

        #--- If more than 5 frequencies are in the listbox, add a scrollbar as to not take up too much space ---#
        if len(InputFrequencies) > 5:
            self.ScrollBarVal = True
        else:
            self.ScrollBarVal = False

        #--- Variable to check if the frequency_list contains frequencies ---#
        self.FrequencyListExists = False

        #--- ListBox containing the frequencies given on line 46 (InputFrequencies) ---#
        self.FrequencyList = Listbox(self.ListboxFrame, relief='groove', exportselection=0, width=5, font=LARGE_FONT, height = 5, selectmode='multiple', bd=3)
        self.FrequencyList.bind('<<ListboxSelect>>',self.FrequencyCurSelect)
        self.FrequencyList.grid(row=1,padx=10,sticky='nswe')
        for frequency in frequencies:
            self.FrequencyList.insert(END, frequency)

        #--- Scroll Bar ---#
        if self.ScrollBarVal:
            self.scrollbar = Scrollbar(self.ListboxFrame, orient="vertical")
            self.scrollbar.config(width=10,command=self.FrequencyList.yview)
            self.scrollbar.grid(row=1,sticky='nse')
            self.FrequencyList.config(yscrollcommand=self.scrollbar.set)

        ManipulateFrequencies = tk.Button(self.ListboxFrame, text = 'Edit', font = MEDIUM_FONT, command = lambda: ManipulateFrequenciesFrame.tkraise()).grid(row=2,column=0,columnspan=4)

        ###########################################################
        ### Frame for adding/deleting frequencies from the list ###
        ###########################################################

        ManipulateFrequenciesFrame = tk.Frame(self, width=10, bd = 3, relief = 'groove')
        ManipulateFrequenciesFrame.grid(row=10,column=2,columnspan=2,padx=10,pady=10, sticky='nsew')

        ManipulateFrequencyLabel = tk.Label(ManipulateFrequenciesFrame, text = 'Enter Frequency(s)', font = MEDIUM_FONT)
        ManipulateFrequencyLabel.grid(row=0,column=0,columnspan=4)

        self.FrequencyEntry = tk.Entry(ManipulateFrequenciesFrame, width=8)
        self.FrequencyEntry.grid(row=1,column=0,columnspan=4)

        AddFrequencyButton = tk.Button(ManipulateFrequenciesFrame, text='Add', font = MEDIUM_FONT, command = lambda: self.AddFrequency()).grid(row=2,column=0)
        DeleteFrequencyButton = tk.Button(ManipulateFrequenciesFrame, text='Delete', font = MEDIUM_FONT, command = lambda: self.DeleteFrequency()).grid(row=2,column=1)
        ClearFrequencyButton = tk.Button(ManipulateFrequenciesFrame, text='Clear', font = MEDIUM_FONT, command = lambda: self.Clear()).grid(row=3,column=0,columnspan=2)

        ReturnButton = tk.Button(ManipulateFrequenciesFrame, text = 'Return', font = MEDIUM_FONT, command = lambda: self.Return()).grid(row=4,column=0,columnspan=2)

        ManipulateFrequenciesFrame.rowconfigure(0, weight=1)
        ManipulateFrequenciesFrame.rowconfigure(1, weight=1)
        ManipulateFrequenciesFrame.rowconfigure(2, weight=1)
        ManipulateFrequenciesFrame.rowconfigure(3, weight=1)
        ManipulateFrequenciesFrame.rowconfigure(4, weight=1)
        ManipulateFrequenciesFrame.columnconfigure(0, weight=1)
        ManipulateFrequenciesFrame.columnconfigure(1, weight=1)


        ############################################################
        ### Adjustment of Visualization Parameters: xstart, xend ###
        ############################################################

        #--- Create a frame that will contain all of the widgets ---#
        AdjustmentFrame = tk.Frame(self, relief = 'groove', bd=3)
        AdjustmentFrame.grid(row=11,column=0,columnspan=4,pady=10)
        AdjustmentFrame.rowconfigure(0, weight=1)
        AdjustmentFrame.rowconfigure(1, weight=1)
        AdjustmentFrame.rowconfigure(2, weight=1)
        AdjustmentFrame.rowconfigure(3, weight=1)
        AdjustmentFrame.rowconfigure(4, weight=1)
        AdjustmentFrame.columnconfigure(0, weight=1)
        AdjustmentFrame.columnconfigure(1, weight=1)
        AdjustmentFrame.columnconfigure(2, weight=1)
        AdjustmentFrame.columnconfigure(3, weight=1)

        #--- Y Limit Adjustment Variables ---#
        self.y_limit_parameter_label = tk.Label(AdjustmentFrame, text = 'Select Y Limit Parameters',font=LARGE_FONT)
        self.y_limit_parameter_label.grid(row=0,column=0,columnspan=4,pady=5,padx=5)

        #--- Raw Data Minimum Parameter Adjustment ---#
        self.raw_data_min_parameter_label = tk.Label(AdjustmentFrame, text = 'Raw Min. Factor',font=MEDIUM_FONT)
        self.raw_data_min_parameter_label.grid(row=1,column=0)
        self.raw_data_min = tk.Entry(AdjustmentFrame, width=5)
        self.raw_data_min.insert(END, '0.5')                   # initial minimum is set to 0.5*minimum current (baseline) of file 1
        self.raw_data_min.grid(row=2,column=0,padx=5,pady=2,ipadx=2)

        #--- Raw Data Maximum Parameter Adjustment ---#
        self.raw_data_max_parameter_label = tk.Label(AdjustmentFrame, text = 'Raw Max. Factor',font=MEDIUM_FONT)
        self.raw_data_max_parameter_label.grid(row=3,column=0)
        self.raw_data_max = tk.Entry(AdjustmentFrame, width=5)
        self.raw_data_max.insert(END, '2')                      # initial adjustment is set to 2x the max current (Peak Height) of file 1
        self.raw_data_max.grid(row=4,column=0,padx=5,pady=2,ipadx=2)

        #--- Peak Height/AUC Minimum Parameter Adjustment ---#
        self.data_min_parameter_label = tk.Label(AdjustmentFrame, text = 'Data Min. Factor',font=MEDIUM_FONT)
        self.data_min_parameter_label.grid(row=1,column=1)
        self.data_min = tk.Entry(AdjustmentFrame, width=5)
        self.data_min.insert(END, '0.5')                   # initial minimum is set to 0.5*minimum current (baseline) of file 1
        self.data_min.grid(row=2,column=1,padx=5,pady=2,ipadx=2)

        #--- Peak Height/AUC Maximum Parameter Adjustment ---#
        self.data_max_parameter_label = tk.Label(AdjustmentFrame, text = 'Data Max. Factor',font=MEDIUM_FONT)
        self.data_max_parameter_label.grid(row=3,column=1)
        self.data_max = tk.Entry(AdjustmentFrame, width=5)
        self.data_max.insert(END, '2')                      # initial adjustment is set to 2x the max current (Peak Height) of file 1
        self.data_max.grid(row=4,column=1,padx=5,pady=2,ipadx=2)

        #--- Select Analysis Method---#
        Options = ['Peak Height Extraction','Area Under the Curve']
        OptionsLabel = tk.Label(self, text='Select Data to be Plotted', font=LARGE_FONT)
        self.PlotOptions = Listbox(self, relief='groove', exportselection=0, font=LARGE_FONT, height=len(Options), selectmode='single', bd=3)
        self.PlotOptions.bind('<<ListboxSelect>>', self.SelectPlotOptions)
        OptionsLabel.grid(row=12,column=0,columnspan=4)
        self.PlotOptions.grid(row=13,column=0,columnspan=4)
        for option in Options:
            self.PlotOptions.insert(END, option)

        #--- Warning label for if the user does not select an analysis method ---#
        self.NoOptionsSelected = tk.Label(self, text = 'Select a Data Analysis Method', font = MEDIUM_FONT, fg='red')   # will only be added to the grid (row 16) if they dont select an option
        self.NoSelection = False

        #--- Ask the User if they want to export the data to a .txt file ---#
        self.SaveVar = BooleanVar()
        self.SaveVar.set(False)
        self.SaveBox = Checkbutton(self, variable=self.SaveVar, onvalue=True, offvalue=False, text="Export Data").grid(row=17,column=0,columnspan=4)

        #--- Quit Button ---#
        self.QuitButton = ttk.Button(self, width=9, text='Quit Program',command=lambda: quit())
        self.QuitButton.grid(row=18,column=0,columnspan=2,pady=10,padx=10)

        #--- Button to Initialize Data Analysis --#
        StartButton = ttk.Button(self, width=9, text='Initialize', command = lambda: self.CheckPoint())
        StartButton.grid(row=18,column=2,columnspan=2, pady=10, padx=10)

        for row in range(18):
            row += 1
            self.rowconfigure(row, weight = 1)

        self.columnconfigure(0, weight = 1)
        self.columnconfigure(1, weight = 1)
        self.columnconfigure(2, weight = 1)
        self.columnconfigure(3, weight = 1)

        ### Raise the initial frame for Electrode and Frequency Selection ###
        self.ListboxFrame.tkraise()
        self.ElectrodeListboxFrame.tkraise()


    #################################################
    ### Functions to track Selections and Entries ###
    #################################################

    def AddFrequency(self):
        Frequencies = self.FrequencyEntry.get()
        self.FrequencyEntry.delete(0,END)

        if Frequencies is not None:
            FrequencyList = Frequencies.split(' ')
            for frequency in FrequencyList:
                if int(frequency) not in InputFrequencies:
                    InputFrequencies.append(int(frequency))
            InputFrequencies.sort()

            self.FrequencyList.delete(0,1)
            self.FrequencyList.delete(0,END)

            for frequency in InputFrequencies:
                self.FrequencyList.insert(END,frequency)


    def DeleteFrequency(self):
        Frequencies = self.FrequencyEntry.get()
        self.FrequencyEntry.delete(0,END)

        if Frequencies is not None:
            FrequencyList = Frequencies.split(' ')

            for Frequency in FrequencyList:

                Frequency = int(Frequency)
                if Frequency in InputFrequencies:
                    InputFrequencies.remove(Frequency)

                self.FrequencyList.delete(0,END)

                for frequency in InputFrequencies:
                    self.FrequencyList.insert(END,int(frequency))

    def Clear(self):
        global InputFrequencies

        self.FrequencyList.delete(0, tk.END)
        InputFrequencies = []

    def Return(self):
        self.ListboxFrame.tkraise()
        self.FrequencyEntry.delete(0,tk.END)


    def ElectrodeSelect(self, variable):
        global e_var

        if variable == 'Multiplex':
            e_var = 'multiple'

            self.SingleElectrodeFile['style'] = 'Off.TButton'
            self.MultipleElectrodeFiles['style'] = 'On.TButton'

        elif variable == 'Multipotentiostat':
            e_var = 'single'

            self.SingleElectrodeFile['style'] = 'On.TButton'
            self.MultipleElectrodeFiles['style'] = 'Off.TButton'


    def FindFile(self, parent):
        global FilePath, ExportPath, FoundFilePath, NoSelectedPath

        try:

            ### prompt the user to select a  ###
            ### directory for  data analysis ###
            FilePath = filedialog.askdirectory(parent = parent)
            FilePath = ''.join(FilePath + '/')


            ### Path for directory in which the    ###
            ### exported .txt file will be placed  ###
            ExportPath = FilePath.split('/')

            #-- change the text of the find file button to the folder the user chose --#
            DataFolder = '%s/%s' % (ExportPath[-3],ExportPath[-2])

            self.SelectFilePath['style'] = 'On.TButton'
            self.SelectFilePath['text'] = DataFolder



            del ExportPath[-1]
            del ExportPath[-1]
            ExportPath = '/'.join(ExportPath)
            ExportPath = ''.join(ExportPath + '/')

            ## Indicates that the user has selected a File Path ###
            FoundFilePath = True

            if self.PathWarningExists:
                self.NoSelectedPath['text'] = ''
                self.NoSelectedPath.grid_forget()

        except:
            print('\n\nInputPage.FindFile: Could Not Find File Path\n\n')


    #--- Analysis Method ---#
    def SelectPlotOptions(self, evt):
        global SelectedOptions
        SelectedOptions = str((self.PlotOptions.get(self.PlotOptions.curselection())))

    def SelectXaxisOptions(self, evt):
        global XaxisOptions
        XaxisOptions = str((self.XaxisOptions.get(self.XaxisOptions.curselection())))

    #--- Electrode Selection ---#
    def ElectrodeCurSelect(self, evt):
        global electrode_count, electrode_list, electrode_dict

        electrode_list = [self.ElectrodeCount.get(idx) for idx in self.ElectrodeCount.curselection()]
        electrode_list = [int(electrode) for electrode in electrode_list]
        electrode_count = len(electrode_list)

        index = 0
        electrode_dict = {}
        for electrode in electrode_list:
            electrode_dict[electrode] = index
            index += 1

        if electrode_count is 0:
            self.ElectrodeListExists = False
            self.ElectrodeLabel['fg'] = 'red'

        elif electrode_count is not 0:
            self.ElectrodeListExists = True
            self.ElectrodeLabel['fg'] = 'black'

    #--- Frequency Selection ---#
    def FrequencyCurSelect(self, evt):
        global frequency_list, frequency_dict


        frequency_list = [self.FrequencyList.get(idx) for idx in self.FrequencyList.curselection()]

        if len(frequency_list) is not 0:

            self.FrequencyListExists = True
            self.FrequencyLabel['fg'] = 'black'

            for var in frequency_list:
                var = int(var)

            #--- Frequency Dictionary ---#
            frequency_dict = {}
            count = 0
            for frequency in frequency_list:
                frequency = int(frequency)
                frequency_dict[frequency] = count
                count += 1

        elif len(frequency_list) is 0:
            self.FrequencyListExists = False
            self.FrequencyLabel['fg'] = 'red'


    #--- Functions to switch frames and plots ---#
    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()

    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame):
        frame.tkraise()

    #####################################################################
    ### Check to see if the user has filled out all  required fields: ###
    ### Electrodes, Frequencies, Analysis Method, and File Path. If   ###
    ### they have, initialize the program                             ###
    #####################################################################
    def CheckPoint(self):
        global mypath, Option, SelectedOptions, FileHandle, AlreadyInitiated

        try:
            #--- check to see if the data analysis method has been selected by the user ---#
            Option = SelectedOptions

            #--- If a data analysis method was selected and a warning label was already created, forget it ---#
            if self.NoSelection:
                self.NoSelection = False
                self.NoOptionsSelected.grid_forget()
        except:
            #--- if no selection was made, create a warning label telling the user to select an analysis method ---#
            self.NoSelection = True
            self.NoOptionsSelected.grid(row=14,column=0,columnspan=4)


        #########################################################
        ### Initialize Canvases and begin tracking animation  ###
        #########################################################
        try:
            mypath = FilePath                       # file path
            FileHandle = str(self.filehandle.get()) # handle for exported .txt file
            FileHandle = ''.join(ExportPath + FileHandle)

            if self.PathWarningExists:
                self.NoSelectedPath.grid_forget()
                self.PathWarningExists = False

        except:
            #-- if the user did not select a file path for data analysis, raise a warning label ---#
            if not FoundFilePath:
                self.NoSelectedPath.grid(row=1,column=0,columnspan=4)
                self.PathWarningExists = True

        if not self.FrequencyListExists:
            self.FrequencyLabel['fg'] = 'red'
        elif self.FrequencyListExists:
            self.FrequencyLabel['fg'] = 'black'

        if not self.ElectrodeListExists:
            self.ElectrodeLabel['fg'] = 'red'
        elif self.ElectrodeListExists:
            self.ElectrodeLabel['fg'] = 'black'


        if not self.PathWarningExists:
            if not self.NoSelection:
                if self.FrequencyListExists:
                    initialize = self.InitializeSetup()



                else:
                    print('Could Not Start Program')

    def InitializeSetup(self):
        global FileHandle, text_file_export, starting_file, handle_variable, Interval, SaveVar, numFiles, q, delimiter, max_data,min_data,min_raw,max_raw

        #---Get the User Input and make it globally accessible---#

        q = Queue()

        ### Set the delimiter value for data columns ###
        if delimiter == 1:
            delimiter = ' '
        elif delimiter == 2:
            delimiter = '   '

        starting_file = 1

        SaveVar = self.SaveVar.get()                        # tracks if text file export has been activated
        handle_variable = self.ImportFileEntry.get()        # string handle used for the input file
        numFiles = int(self.numfiles.get())


        min_raw = float(self.raw_data_min.get())            # raw data y limit adjustment variables
        max_raw = float(self.raw_data_max.get())
        min_data = float(self.data_min.get())            # raw data y limit adjustment variables
        max_data = float(self.data_max.get())

        #############################################################
        ### Interval at which the program searches for files (ms) ###
        #############################################################
        Interval = self.Interval.get()

        #################################
        ### Initiate .txt File Export ###
        #################################

        #--- If the user has indicated that text file export should be activated ---#
        if SaveVar:
            print('Initializing Text File Export')
            text_file_export = True

        else:
            text_file_export = None
            print('Text File Export Deactivated')

        ## set the resizeability of the container ##
        ## frame to handle PlotContainer resize   ##
        container.columnconfigure(1, weight=1)

        ################################################################
        ### If all checkpoints have been met, initialize the program ###
        ################################################################
        if FoundFilePath:
            checkpoint = CheckPoint(self.parent, self.controller)


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#

####################################
### Checkpoint TopLevel Instance ###
####################################
class CheckPoint():
    def __init__(self, parent, controller):

        #-- Check to see if the user's settings are accurate
        #-- Search for the presence of the files. If they exist,
        #-- initialize the functions and frames for Real Time Analysis

        self.win = tk.Toplevel()
        self.win.wm_title("CheckPoint")

        title = tk.Label(self.win, text = 'Searching for files...',font=HUGE_FONT).grid(row=0,column=0,columnspan=2,pady=10,padx=10,sticky='news')

        self.parent = parent
        self.win.transient(self.parent)
        self.win.attributes('-topmost', 'true')
        self.controller = controller

        row_value = 0
        column_value = 0

        self.label_dict = {}
        self.already_verified = {}

        self.electrode_frame = tk.Frame(self.win)
        self.electrode_frame.grid(row=1,column=0,columnspan=2)

        for electrode in electrode_list:


            electrode_label = tk.Label(self.electrode_frame, text = 'E%s' % electrode,font=HUGE_FONT)
            electrode_label.grid(row=row_value,column=column_value,pady=5,padx=5)
            self.label_dict[electrode] = electrode_label
            self.already_verified[electrode] = False

            if column_value == 1:
                column_value = 0
                row_value += 1
            else:
                column_value = 1

        row_value += 1

        self.stop = tk.Button(self.win, text = 'Stop', command = self.stop)
        self.stop.grid(row=row_value, column=0,columnspan=2,pady=5)
        self.StopSearch = False

        self.num = 0
        self.count = 0
        self.analysis_count = 0
        self.electrode_limit = electrode_count - 1
        self.frequency_limit = len(frequency_list) - 1

        root.after(50,self.verify)

    def verify(self):

        self.electrode = electrode_list[self.num]

        frequency = frequency_list[0]

        filename, filename2, filename3 = _retrieve_file(1,self.electrode,frequency)

        myfile = mypath + filename               ### path of your file
        myfile2 = mypath + filename2               ### path of your file
        myfile3 = mypath + filename3               ### path of your file

        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

        except:
            try:
                mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                myfile = myfile2
            except:
                try:
                    mydata_bytes = os.path.getsize(myfile3)    ### retrieves the size of the file in bytes
                    myfile = myfile3
                except:
                    mydata_bytes = 1

        if mydata_bytes > 1000:
            if not self.already_verified[self.electrode]:
                self.already_verified[self.electrode] = True
                if not self.StopSearch:
                    self.label_dict[self.electrode]['fg'] = 'green'
                    self.analysis_count += 1

            if self.analysis_count == electrode_count:
                root.after(10,self.proceed)


        if self.num < self.electrode_limit:
            self.num += 1
        else:
            self.num = 0

        if self.analysis_count < electrode_count:
            if not self.StopSearch:
                root.after(10,self.verify)

    def proceed(self):
        global wait_time, track, initialize, post_analysis

        self.win.destroy()

        ##############################
        ### Syncronization Classes ###
        ##############################
        track = Track()

        ######################################################
        ### Matplotlib Canvas, Figure, and Artist Creation ###
        ######################################################
        initialize = InitializeFigureCanvas()

        ############################
        ### Post Analysis Module ###
        ############################
        #post_analysis = PostAnalysis(self.parent, self.controller)
        #ShowFrames[PostAnalysis] = post_analysis
        #post_analysis.grid(row=0, column=0, sticky = 'nsew')

        ################################################
        ### Initialize the RealTimeManipulationFrame ###
        ################################################
        frame = RealTimeManipulationFrame(container, self)
        ShowFrames[RealTimeManipulationFrame] = frame
        frame.grid(row=0, column=0, sticky='nsew')


        #---When initliazed, raise the Start Page and the plot for electrode one---#
        self.show_frame(RealTimeManipulationFrame)              # raises the frame for real-time data manipulation
        self.show_plot(PlotValues[0])           # raises the figure for electrode 1


    def stop(self):
        self.StopSearch = True
        self.win.destroy()

    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame):
        frame.tkraise()

    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()



                                        #############################################################
                                        #############################################################
                                        ### Creation of Matplotlib Canvas, Figures, Axes, Artists ###
                                        ### and all other decorators (e.g. axis labels, titles)   ###
                                        #############################################################
                                        #############################################################


class InitializeFigureCanvas():
    def __init__(self):
        global text_file_export, SaveVar, file_list, electrode_count, anim, Frame, FrameReference, FileHandle, PlotContainer,frequency_list, data_list, plot_list, EmptyPlots, figures, frame_list, Plot, PlotFrames, PlotValues

        ##############################################
        ### Generate global lists for data storage ###
        ##############################################

        self.length = len(frequency_list)
        electrode_count = int(electrode_count)

        #--- Animation list ---#
        anim = []

        #--- file list ---#
        file_list = [0]*numFiles

        #--- Figure lists ---#
        figures = []

        ############################################
        ### Create global lists for data storage ###
        ############################################
        data_list = [0]*electrode_count                             # Peak Height/AUC data (after smoothing and polynomial regression)

        for num in range(electrode_count):
            data_list[num] = [0]*self.length                        # a data list for each eletrode
            for count in range(self.length):                        # a data list for each frequency for that electrode
                data_list[num][count] = [0]*numFiles


        #--- Lists of Frames and Artists ---#
        plot_list = []
        frame_list = []

        ######################################################
        ### Create a figure and artists for each electrode ###
        ######################################################
        for num in range(electrode_count):
            electrode = electrode_list[num]
            figure = self.MakeFigure(electrode)
            figures.append(figure)


        #####################################################
        ### Create a frame for each electrode and embed   ###
        ### within it the figure containing its artists   ###
        #####################################################

        PlotFrames = {}                # Dictionary of frames for each electrode
        PlotValues = []                # create a list of frames

        #--- Create a container that can be created and destroyed when Start() or Reset() is called, respectively ---#
        PlotContainer = tk.Frame(container, relief = 'groove', bd = 3)
        PlotContainer.grid(row=0,column=1, sticky = 'nsew')
        PlotContainer.rowconfigure(0, weight=1)
        PlotContainer.columnconfigure(0, weight=1)

        frame_count = 0
        for electrode_frame in frame_list:                # Iterate through the frame of each electrode

            #--- create an instance of the frame and append it to the global frame dictionary ---#
            FrameReference = VisualizationFrame(electrode_frame, frame_count, PlotContainer, self)            # PlotContainer is the 'parent' frame
            FrameReference.grid(row=0,column=0,sticky='nsew')      # sticky must be 'nsew' so it expands and contracts with resize
            PlotFrames[electrode_frame] = FrameReference

            frame_count += 1

        #--- Create a list containing the Frame objects for each electrode ---#
        for reference, frame in PlotFrames.items():
            PlotValues.append(frame)


        #################################
        ### Initiate .txt File Export ###
        #################################

        #--- If the user has indicated that text file export should be activated ---#
        if SaveVar:
            print('Initializing Text File Export')
            text_file_export = TextFileExport()

        else:
            text_file_export = None
            print('Text File Export Deactivated')



    ############################################
    ### Create the figure and artist objects ###
    ############################################
    def MakeFigure(self, electrode):
        global EmptyPlots, plot_list, frame_list

        try:
            ##########################################
            ### Setup the Figure for voltammograms ###
            ##########################################
            fig, ax = plt.subplots(nrows=2,ncols=1,squeeze=False,figsize=(9,4.5))    ## figsize=(width, height)
            plt.subplots_adjust(bottom=0.2,hspace=0.6,wspace=0.3)         ### adjust the spacing between subplots
            #---Set the electrode index value---#
            if e_var == 'single':
                list_val = (electrode*3)
            elif e_var == 'multiple':
                list_val = 3
            #######################
            ### Set axis labels ###
            #######################
            ax[0,0].set_ylabel('Current (µA)',fontweight='bold')
            ax[0,0].set_xlabel('Voltage (V)',fontweight='bold')

            ax[1,0].set_ylabel('Charge (uC)',fontweight='bold')
            ax[1,0].set_xlabel('Frequency (Hz)',fontweight='bold')
            ##########################################
            ### Set suplot axes for each frequency ###
            ##########################################
            electrode_plot = []

            max_frequency = frequency_list[-1]
            ax[1,0].set_xscale('log')
            #################################################################################
            #################################################################################
            ###       Analyze the first file and create the Y limits of the subplots      ###
            ###               depending on the data range of the first file               ###
            #################################################################################

            self.InitializeSubplots(ax, electrode)

            #################################################################################
            #################################################################################

            #---Initiate the subplots---#
            # this assigns a Line2D artist object to the artist object (Axes)
            smooth, = ax[0,0].plot([],[],'ko',Markersize=2)
            regress, = ax[0,0].plot([],[],'r-')
            charge, = ax[1,0].plot([],[],'ko',MarkerSize=1)

            #--- shading for AUC ---#
            verts = [(0,0),*zip([],[]),(0,0)]
            poly = Polygon(verts, alpha = 0.5)
            ax[0,0].add_patch(poly)

            #####################################################
            ### Create a list of the primitive artists        ###
            ### (Line2D objects) that will be returned        ###
            ### to ElectrochemicalAnimation to be visualized  ###
            #####################################################

            # this is the list that will be returned as _drawn_artists to the Funcanimation class
            plots = [smooth,regress,charge,poly]

            #--- And append that list to keep a global reference ---#
            electrode_plot.append(plots)        # 'plots' is a list of artists that are passed to animate
            electrode_frame = 'Electrode %s' % str(electrode)
            if electrode_frame not in frame_list:
                frame_list.append(electrode_frame)

            #--- Create empty plots to return to animate for initializing---#
            EmptyPlots = [smooth,regress,charge]

            plot_list.append(plots)        # 'plot_list' is a list of lists containing 'plots' for each electrode

            #-- Return both the figure and the axes to be stored as global variables --#
            return fig, ax


        except:
            print('Error in MakeFigure')


    #####################################################################################
    ### Initalize Y Limits of each figure depending on the y values of the first file ###
    #####################################################################################
    def InitializeSubplots(self,ax,electrode):

        self.list_val = _get_listval(electrode)

        try:
            print('Initialize: Electrode %s' % str(electrode))
            frequency = frequency_list[0]
            filename, filename2, filename3 = _retrieve_file(1,electrode,frequency)

            myfile = mypath + filename               ### path of your file
            myfile2 = mypath + filename2               ### path of your file
            myfile3 = mypath + filename3               ### path of your file

            try:
                mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

            except:
                try:
                    mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                    myfile = myfile2
                except:
                    try:
                        mydata_bytes = os.path.getsize(myfile3)    ### retrieves the size of the file in bytes
                        myfile = myfile3
                    except:
                        mydata_bytes = 1


            if mydata_bytes > 1000:
                print('Found File %s' % myfile)
                self.RunInitialization(myfile,ax,electrode)

            else:
                return False


        except:
            print('could not find file for electrode %d' % electrode)
            #--- If search time has not met the search limit keep searching ---#
            root.after(1000, self.InitializeSubplots, ax, electrode)


    def RunInitialization(self, myfile, ax, electrode):
        global high_xstart, high_xend, low_xstart, low_xend

        try:
            avalue = 0
            #########################
            ### Retrieve the data ###
            #########################
            potentials, currents, data_dict = ReadData(myfile, electrode)

            ##########################################
            ### Set the x axes of the voltammogram ###
            ##########################################
            MIN_POTENTIAL = min(potentials)
            MAX_POTENTIAL = max(potentials)

            #-- Reverse voltammogram to match the 'Texas' convention --#
            ax[0,0].set_xlim(MAX_POTENTIAL,MIN_POTENTIAL)

            #######################################
            ### Get the high and low potentials ###
            #######################################

            #-- set the local variables to the global ---#
            xstart = max(potentials)
            xend = min(potentials)

            low_xstart = xstart
            high_xstart = xstart
            low_xend = xend
            high_xend = xend

            cut_value = 0
            for value in potentials:
                if value == 0:
                    cut_value += 1


            if cut_value > 0:
                potentials = potentials[:-cut_value]
                currents = currents[:-cut_value]

            adjusted_potentials = [value for value in potentials if xend <= value <= xstart]

            #########################################
            ### Savitzky-Golay smoothing          ###
            #########################################
            smooth_currents = savgol_filter(currents, 15, sg_degree)
            data_dict = dict(zip(potentials,smooth_currents))

            #######################################
            ### adjust the smooth currents to   ###
            ### match the adjusted potentials   ###
            #######################################
            adjusted_currents = []
            for potential in adjusted_potentials:
                adjusted_currents.append(data_dict[potential])

            ######################
            ### Polynomial fit ###
            ######################
            polynomial_coeffs = np.polyfit(adjusted_potentials,adjusted_currents,polyfit_deg)
            eval_regress = np.polyval(polynomial_coeffs,adjusted_potentials).tolist()
            regression_dict = dict(zip(eval_regress, adjusted_potentials))      # dictionary with current: potential

            fit_half = round(len(eval_regress)/2)
            min1 = min(eval_regress[:-fit_half])
            min2 = min(eval_regress[fit_half:])
            max1 = max(eval_regress[:-fit_half])
            max2 = max(eval_regress[fit_half:])

            linear_fit = np.polyfit([regression_dict[min1],regression_dict[min2]],[min1,min2],1)
            linear_regression = polyval(linear_fit,[regression_dict[min1],regression_dict[min2]]).tolist()

            Peak_Height = max(max1,max2)-min(min1,min2)


            if SelectedOptions == 'Area Under the Curve':
                AUC_index = 1
                AUC = 0

                AUC_potentials = [abs(potential) for potential in adjusted_potentials]
                AUC_min = min(adjusted_currents)
                AUC_currents = [Y - AUC_min for Y in adjusted_currents]

                while AUC_index <= len(AUC_currents) - 1:
                    AUC_height = (AUC_currents[AUC_index] + AUC_currents[AUC_index - 1])/2
                    AUC_width = AUC_potentials[AUC_index] - AUC_potentials[AUC_index - 1]
                    AUC += (AUC_height * AUC_width)
                    AUC_index += 1

            #--- calculate the baseline current ---#
            minimum_current = min(min1,min2)
            maximum_current = max(max1,max2)
            peak_current = maximum_current - minimum_current

            ## Reverse voltammogram to match the 'Texas' convention ##
            ax[0,0].set_xlim(MAX_POTENTIAL,MIN_POTENTIAL)

            ## set the limits of the lovric plot ##
            ax[1,0].set_ylim(0,max_data*.05)
            ax[1,0].set_xlim(int(frequency_list[0]),int(frequency_list[-1]))

            if SelectedOptions == 'Peak Height Extraction':
                ax[0,0].set_ylim(min_raw*minimum_current,10*max_raw*max(max1,max2))         # voltammogram

            elif SelectedOptions == 'Area Under the Curve':
                ax[0,0].set_ylim(0,max_raw*max(max1,max2))


            return True

        except:
            print('\n\nError in RunInitialization\n\n')




                #############################################################
                #############################################################
                ###              END OF INITIATION FUNCTIONS              ###
                #############################################################
                #############################################################





#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#




############################################################
### Frame displayed during experiment with widgets and   ###
### functions for Real-time Data Manipulation            ###
############################################################
class RealTimeManipulationFrame(tk.Frame):

    def __init__(self, parent, controller):
        global PlotValues, container, high_xstart_entry, low_xstart_entry, high_xend_entry, low_xend_entry, ShowFrames

        tk.Frame.__init__(self, parent)         # Initialize the frame


        #######################################
        #######################################
        ### Pack the widgets into the frame ###
        #######################################
        #######################################

        #################################################
        ### Nested Frame for Real-Time adjustment     ###
        ### of voltammogram, polynomial regression,   ###
        ### and savitzky-golay Smoothing              ###
        #################################################

        RegressionFrame = tk.Frame(self,relief='groove',bd=5)
        RegressionFrame.grid(row=7,column=0,columnspan=4,pady=5,padx=5,ipadx=3, sticky='ns')
        RegressionFrame.rowconfigure(0, weight=1)
        RegressionFrame.rowconfigure(1, weight=1)
        RegressionFrame.rowconfigure(2, weight=1)
        RegressionFrame.columnconfigure(0, weight=1)
        RegressionFrame.columnconfigure(1, weight=1)

        #--- Title ---#
        self.RegressionLabel = tk.Label(RegressionFrame, text = 'Real Time Analysis Manipulation', font=LARGE_FONT)
        self.RegressionLabel.grid(row=0,column=0,columnspan=4,pady=5,padx=5)


        ###################################################################
        ### Real Time Manipulation of Savitzky-Golay Smoothing Function ###
        ###################################################################
        self.SmoothingLabel = tk.Label(RegressionFrame, text = 'Savitzky-Golay Window (mV)', font = LARGE_FONT)
        self.SmoothingLabel.grid(row=1,column=0,columnspan=4,pady=1)
        self.SmoothingEntry = tk.Entry(RegressionFrame, width=10)
        self.SmoothingEntry.grid(row=2,column=0,columnspan=4,pady=3)
        self.SmoothingEntry.insert(END, sg_window)

        #--- Check for the presence of high and low frequencies ---#
        if frequency_list[-1] > 50:
            self.High = True
        else:
            self.High = False
        if frequency_list[0] <= 50:
            self.Low = True
        else:
            self.Low = False
        ###################################################
        ### If a frequency <= 50Hz exists, grid a frame ###
        ### for low frequency data manipulation         ###
        ###################################################
        if self.Low is True:
            LowParameterFrame = tk.Frame(RegressionFrame)
            LowParameterFrame.grid(row=3,column=0,columnspan=4, sticky='nsew')
            LowParameterFrame.rowconfigure(0, weight=1)
            LowParameterFrame.rowconfigure(1, weight=1)
            LowParameterFrame.rowconfigure(2, weight=1)
            LowParameterFrame.columnconfigure(0, weight=1)
            LowParameterFrame.columnconfigure(1, weight=1)
            ShowFrames['LowParameterFrame'] = LowParameterFrame

            #--- points discarded at the beginning of the voltammogram, xstart ---#
            self.low_xstart_label = tk.Label(LowParameterFrame, text = 'xstart (V)', font=MEDIUM_FONT).grid(row=0,column=0)
            self.low_xstart_entry = tk.Entry(LowParameterFrame, width=7)
            self.low_xstart_entry.insert(END, str(low_xstart))
            self.low_xstart_entry.grid(row=1,column=0)
            low_xstart_entry = self.low_xstart_entry

            #--- points discarded at the beginning of the voltammogram, xend ---#
            self.low_xend_label = tk.Label(LowParameterFrame, text = 'xend (V)', font=MEDIUM_FONT).grid(row=0,column=1)
            self.low_xend_entry = tk.Entry(LowParameterFrame, width=7)
            self.low_xend_entry.insert(END, str(low_xend))
            self.low_xend_entry.grid(row=1,column=1)
            low_xend_entry = self.low_xend_entry

        ##################################################
        ### If a frequency > 50Hz exists, grid a frame ###
        ### for high frequency data manipulation       ###
        ##################################################
        if self.High is True:
            HighParameterFrame = tk.Frame(RegressionFrame)
            HighParameterFrame.grid(row=3,column=0,columnspan=4, sticky='nsew')
            HighParameterFrame.rowconfigure(0, weight=1)
            HighParameterFrame.rowconfigure(1, weight=1)
            HighParameterFrame.rowconfigure(2, weight=1)
            HighParameterFrame.columnconfigure(0, weight=1)
            HighParameterFrame.columnconfigure(1, weight=1)
            ShowFrames['HighParameterFrame'] = HighParameterFrame

            #--- points discarded at the beginning of the voltammogram, xstart ---#
            self.high_xstart_label = tk.Label(HighParameterFrame, text = 'xstart (V)', font=MEDIUM_FONT).grid(row=0,column=0)
            self.high_xstart_entry = tk.Entry(HighParameterFrame, width=7)
            self.high_xstart_entry.insert(END, str(high_xstart))
            self.high_xstart_entry.grid(row=1,column=0)
            high_xstart_entry = self.high_xstart_entry

            #--- points discarded at the beginning of the voltammogram, xend ---#
            self.high_xend_label = tk.Label(HighParameterFrame, text = 'xend (V)', font=MEDIUM_FONT).grid(row=0,column=1)
            self.high_xend_entry = tk.Entry(HighParameterFrame, width=7)
            self.high_xend_entry.insert(END, str(high_xend))
            self.high_xend_entry.grid(row=1,column=1)
            high_xend_entry = self.high_xend_entry

        ############################################################
        ### If both high and low frequencies are being analyzed, ###
        ### create buttons to switch between the two             ###
        ############################################################
        if self.High is True:
            if self.Low is True:
                self.SelectLowParameters = ttk.Button(RegressionFrame, style = 'Off.TButton', text = 'f <= 50Hz', command = lambda: self.show_frame('LowParameterFrame'))
                self.SelectLowParameters.grid(row=4,column=0,pady=5,padx=5)

                self.SelectHighParameters = ttk.Button(RegressionFrame, style = 'On.TButton', text = 'f > 50Hz', command = lambda: self.show_frame('HighParameterFrame'))
                self.SelectHighParameters.grid(row=4,column=1,pady=5,padx=5)


        #--- Button to apply adjustments ---#
        self.AdjustParameterButton = tk.Button(RegressionFrame, text = 'Apply Adjustments', font=LARGE_FONT, command = lambda: self.AdjustParameters())
        self.AdjustParameterButton.grid(row=5,column=0,columnspan=4,pady=10,padx=10)


        #---Buttons to switch between electrode frames---#
        frame_value = 0
        row_value = 8
        column_value = 0
        for value in PlotValues:
            Button = ttk.Button(self, text=frame_list[frame_value], command = lambda frame_value=frame_value: self.show_plot(PlotValues[frame_value]))
            Button.grid(row=row_value,column=column_value,pady=2,padx=5)

            ## allows .grid() to alternate between
            ## packing into column 1 and column 2
            if column_value == 1:
                column_value = 0
                row_value += 1

            ## if gridding into the 1st column,
            ## grid the next into the 2nd column
            else:
                column_value += 1
            frame_value += 1
        row_value += 1


        #--- Start ---#
        StartButton = ttk.Button(self, text='Start', style='Fun.TButton', command = lambda: self.SkeletonKey())
        StartButton.grid(row=row_value, column=0, pady=5, padx=5)

        #--- Reset ---#
        Reset = ttk.Button(self, text='Reset', style='Fun.TButton', command = lambda: self.Reset())
        Reset.grid(row=row_value, column=1,pady=5, padx=5)
        row_value += 1

        #--- Quit ---#
        QuitButton = ttk.Button(self, text='Quit Program',command=lambda: quit())
        QuitButton.grid(row=row_value,column=0,columnspan=4,pady=5)

        for row in range(row_value):
            row += 1
            self.rowconfigure(row, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)



                                        ###################################################
                                        ###################################################
                                        ###### Real Time Data Manipulation Functions ######
                                        ###################################################
                                        ###################################################



    #################################################
    ### Adjustment of points discarded at the     ###
    ### beginning and end of Regression Analysis  ###
    #################################################
    def AdjustParameters(self):
        #--- Adjusts the parameters used to visualize the raw voltammogram, smoothed currents, and polynomial fit
        global low_xstart, high_xstart, low_xend, high_xend, sg_window

        ###############################################
        ### Polynomical Regression Range Parameters ###
        ###############################################

        if self.Low:
            #--- parameters for frequencies equal or below 50Hz ---#
            low_xstart = float(self.low_xstart_entry.get())          # xstart/xend adjust the points at the start and end of the voltammogram/smoothed currents, respectively
            low_xend = float(self.low_xend_entry.get())
        if self.High:
            #--- parameters for frequencies above 50Hz ---#
            high_xstart = float(self.high_xstart_entry.get())
            high_xend = float(self.high_xend_entry.get())

        #######################################
        ### Savitzky-Golay Smoothing Window ###
        #######################################
        sg_window = float(self.SmoothingEntry.get())
        print('\n\n\nAdjustParamaters: SG_Window (mV) %d\n\n\n' % sg_window)


    ########################################################
    ### Function to Reset and raise the user input frame ###
    ########################################################
    def Reset(self):
        global key, PoisonPill, AlreadyInitiated, AlreadyReset

        key = 0
        PoisonPill = True
        AlreadyInitiated = False # reset the start variable
        AlreadyReset = True

        # Raise the initial user input frame
        self.show_frame(InputFrame)
        self.close_frame(RealTimeManipulationFrame)

    ##########################################################
    ### Function to raise frame to the front of the canvas ###
    ##########################################################
    def show_frame(self, cont):

        frame = ShowFrames[cont]            # Key: frame handle / Value: tk.Frame object
        frame.tkraise()                     # raise the frame objext

        if cont == 'LowParameterFrame':
            self.SelectLowParameters['style'] = 'On.TButton'
            self.SelectHighParameters['style'] = 'Off.TButton'

        elif cont == 'HighParameterFrame':
            self.SelectLowParameters['style'] = 'Off.TButton'
            self.SelectHighParameters['style'] = 'On.TButton'

    ###################################################
    ### Function to start returning visualized data ###
    ###################################################
    def SkeletonKey(self):
        global key, PoisonPill, data_analysis, extrapolate, AlreadyInitiated

        if not AlreadyInitiated:

            ######################################################################
            ### Initialize Animation (Visualization) for each electrode figure ###
            ######################################################################
            fig_count = 0                   # index value for the frame
            for figure in figures:
                fig, self.ax = figure
                electrode = electrode_list[fig_count]
                anim.append(ElectrochemicalAnimation(fig, self.ax, electrode, resize_interval = None, fargs=None))
                fig_count += 1

            AlreadyInitiated = True

            #--- reset poison pill variables --#
            PoisonPill = False

            if key == 0:                                # tells Generate() to start data analysis
                key += 100
        else:
            print('\n\nProgram has already been initiaed\n\n')



    ######################################################
    ### Function to raise frame for specific electrode ###
    ######################################################
    def show_plot(self, frame):
        frame.tkraise()

    #####################################
    ### Destory the frames on Reset() ###
    #####################################
    def close_frame(self, cont):
        frame = ShowFrames[cont]
        frame.grid_forget()

        for value in PlotValues:
            value.destroy()

        PlotContainer.destroy()




#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#




#########################################################
### Electrode Frame Class for data visualization      ###
### displayed next to the RealTimeManipulationFrame   ###
###                                                   ###
###    Embeds a canvas within the tkinter             ###
###    MainWindow containing figures that             ###
###    visualize the data for that electrode          ###
#########################################################

class VisualizationFrame(tk.Frame):
    def __init__(self, electrode, count, parent, controller):
        global FrameFileLabel

        tk.Frame.__init__(self, parent)

        #--- for resize ---#
        self.columnconfigure(0, weight = 2)
        self.columnconfigure(1, weight = 1)
        self.rowconfigure(2, weight=2)

        ElectrodeLabel = tk.Label(self, text='%s' % electrode ,font=HUGE_FONT)
        ElectrodeLabel.grid(row=0,column=0,pady=5,sticky='n')

        FrameFileLabel = tk.Label(self, text = '', font=MEDIUM_FONT)
        FrameFileLabel.grid(row=0,column=1,pady=3,sticky='ne')

        #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
        fig, ax = figures[count]                                                # Call the figure and artists for the electrode
        canvas = FigureCanvasTkAgg(fig, self)                                         # and place the artists within the frame
        canvas.draw()                                                           # initial draw call to create the artists that will be blitted
        canvas.get_tk_widget().grid(row=1,columnspan=2,pady=6,ipady=5,sticky='news')          # does not affect size of figure within plot container

                                        #############################################################
                                        #############################################################
                                        ###                   End of GUI Classes                  ###
                                        #############################################################
                                        #############################################################





#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


    ##########################################################################
    ##########################################################################
    ###   ANIMATION FUNCTION TO HANDLE ALL DATA ANALYSIS AND VISUALIZATION ###
    ##########################################################################
    ##########################################################################


class ElectrochemicalAnimation():
    def __init__(self, fig, ax, electrode, generator = None, func = None, resize_interval = None, fargs = None):

        self.electrode = electrode                               # Electrode for this class instance
        self.num = electrode_dict[self.electrode]                # Electrode index value
        self.spacer = ''.join(['       ']*self.electrode)        # Spacer value for print statements
        self.list_val = (self.electrode*3) + column_index        # Electrode column index value
        self.file = starting_file                                # Starting File
        self.index = 0                                           # File Index Value
        self.ax = ax                                             # Figure Axes object
        self.count = 0                                           # Frequency index value
        self.frequency_limit = len(frequency_list) - 1           # ' -1 ' so it matches the index value

        self.frequency_axis = []
        self.charge_axis = []

        ### Lists for sample rate (time passed)  ###
        ### and file count for each electrode    ###
        self.sample_list = []
        self.file_list = []



        ##############################
        ## Set the generator object ##
        ##############################
        if generator is not None:
            self.generator = generator
        else:
            self.generator = self._raw_generator

        ################################
        ## Set the animation function ##
        ################################
        if func is not None:
            self._func = func
        else:
            self._func = self._animate

        if resize_interval is not None:
            self.resize_interval = resize_interval
        else:
            self.resize_interval = None

        self.resize_limit = self.resize_interval        # set the first limit

        if fargs:
            self._args = fargs
        else:
            self._args = ()

        self._fig = fig

        # Disables blitting for backends that don't support it.  This
        # allows users to request it if available, but still have a
        # fallback that works if it is not.
        self._blit = fig.canvas.supports_blit


        # Instead of starting the event source now, we connect to the figure's
        # draw_event, so that we only start once the figure has been drawn.
        self._first_draw_id = fig.canvas.mpl_connect('draw_event', self._start)

        # Connect to the figure's close_event so that we don't continue to
        # fire events and try to draw to a deleted figure.
        self._close_id = self._fig.canvas.mpl_connect('close_event', self._stop)

        self._setup_blit()


    def _start(self, *args):


        # Starts interactive animation. Adds the draw frame command to the GUI
        # andler, calls show to start the event loop.

        # First disconnect our draw event handler
        self._fig.canvas.mpl_disconnect(self._first_draw_id)
        self._first_draw_id = None  # So we can check on save

        # Now do any initial draw
        self._init_draw()

        ### Create a thread to analyze obtain the file from a Queue
        ### and analyze the data.

        class _threaded_animation(threading.Thread):

            def __init__(self, Queue):
                #global PoisonPill

                threading.Thread.__init__(self)     # initiate the thread

                self.q = Queue

                #-- set the poison pill event for Reset --#
                self.PoisonPill = Event()
                PoisonPill = self.PoisonPill             # global reference

                self.file = 1

                root.after(10,self.start)                       # initiate the run() method

            def run(self):

                while True:
                    try:
                        task = self.q.get(block=False)

                    except:
                        break
                    else:
                        if not PoisonPill:
                            root.after(Interval,task)

                if not analysis_complete:
                    if not PoisonPill:
                        root.after(10, self.run)


        threaded_animation = _threaded_animation(Queue = q)

        self._step()


    def _stop(self, *args):
        # On stop we disconnect all of our events.
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._fig.canvas.mpl_disconnect(self._close_id)

    def _setup_blit(self):
        # Setting up the blit requires: a cache of the background for the
        # axes
        self._blit_cache = dict()
        self._drawn_artists = []
        self._resize_id = self._fig.canvas.mpl_connect('resize_event',
                                                       self._handle_resize)
        self._post_draw(True)

    def _blit_clear(self, artists, bg_cache):
        # Get a list of the axes that need clearing from the artists that
        # have been drawn. Grab the appropriate saved background from the
        # cache and restore.
        axes = {a.axes for a in artists}
        for a in axes:
            if a in bg_cache:
                a.figure.canvas.restore_region(bg_cache[a])


    #######################################################################
    ### Initialize the drawing by returning a sequence of blank artists ###
    #######################################################################
    def _init_draw(self):

        self._drawn_artists = EmptyPlots

        for a in self._drawn_artists:
            a.set_animated(self._blit)


    def _handle_resize(self, *args):
        # On resize, we need to disable the resize event handling so we don't
        # get too many events. Also stop the animation events, so that
        # we're paused. Reset the cache and re-init. Set up an event handler
        # to catch once the draw has actually taken place.

        #################################################
        ### Stop the event source and clear the cache ###
        #################################################
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._blit_cache.clear()
        self._init_draw()
        self._resize_id = self._fig.canvas.mpl_connect('draw_event',
                                                       self._end_redraw)


    def _end_redraw(self, evt):
        # Now that the redraw has happened, do the post draw flushing and
        # blit handling. Then re-enable all of the original events.
        self._post_draw(True)
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._resize_id = self._fig.canvas.mpl_connect('resize_event',
                                                       self._handle_resize)

    def _draw_next_frame(self, framedata, fargs = None):
        # Breaks down the drawing of the next frame into steps of pre- and
        # post- draw, as well as the drawing of the frame itself.
        self._pre_draw(framedata)
        self._draw_frame(framedata, fargs)
        self._post_draw(False)


    def _pre_draw(self, framedata):
        # Perform any cleaning or whatnot before the drawing of the frame.
        # This default implementation allows blit to clear the frame.
        self._blit_clear(self._drawn_artists, self._blit_cache)

    ###########################################################################
    ### Retrieve the data from _animation and blit the data onto the canvas ###
    ###########################################################################
    def _draw_frame(self, framedata, fargs):


        self._drawn_artists = self._func(framedata, *self._args)

        if self._drawn_artists is None:
            raise RuntimeError('The animation function must return a '
                               'sequence of Artist objects.')
        self._drawn_artists = sorted(self._drawn_artists,
                                     key=lambda x: x.get_zorder())

        for a in self._drawn_artists:
            a.set_animated(self._blit)


    def _post_draw(self, redraw):
        # After the frame is rendered, this handles the actual flushing of
        # the draw, which can be a direct draw_idle() or make use of the
        # blitting.

        if redraw:

            # Data plots #
            self._fig.canvas.draw()

        elif self._drawn_artists:

            self._blit_draw(self._drawn_artists, self._blit_cache)


    # The rest of the code in this class is to facilitate easy blitting
    def _blit_draw(self, artists, bg_cache):
        # Handles blitted drawing, which renders only the artists given instead
        # of the entire figure.
        updated_ax = []
        for a in artists:
            # If we haven't cached the background for this axes object, do
            # so now. This might not always be reliable, but it's an attempt
            # to automate the process.
            if a.axes not in bg_cache:
                bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            a.axes.draw_artist(a)
            updated_ax.append(a.axes)

        # After rendering all the needed artists, blit each axes individually.
        for ax in set(updated_ax):
            ax.figure.canvas.blit(ax.bbox)


    ## callback that is called every 'interval' ms ##
    def _step(self):
        global file_list, sample_list, analysis_complete

        ### look for the file here ###
        frequency = int(frequency_list[self.count])

        filename, filename2, filename3 = _retrieve_file(1,self.electrode,frequency)

        myfile = mypath + filename               ### path of your file
        myfile2 = mypath + filename2               ### path of your file
        myfile3 = mypath + filename3               ### path of your file

        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

        except:
            try:
                mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                myfile = myfile2
            except:
                try:
                    mydata_bytes = os.path.getsize(myfile3)    ### retrieves the size of the file in bytes
                    myfile = myfile3
                except:
                    mydata_bytes = 1

        #################################################################
        #### If the file meets the size requirement, analyze the data ###
        #################################################################
        if mydata_bytes > 1000:
            #self.FileLabel['text'] = 'Found %s' % filename
            print('%s%d: Queueing %s' % (self.spacer,self.electrode,filename))
            q.put(lambda: self._run_analysis(myfile,frequency))

        else:
            if not PoisonPill:
                root.after(100,self._step)

    def _check_queue():

        while True:
            try:
                print('%sChecking Queue' % self.spacer)
                task = q.get(block=False)
            except:
                print('%sQueue Empty' % self.spacer)
                break
            else:
                if not PoisonPill:
                    root.after(1,self.task)

        if not analysis_complete:
            if not PoisonPill:
                root.after(5, self._check_queue)

    def _run_analysis(self,myfile,frequency):
        #global post_analysis

        #######################################################
        ### Perform the next iteration of the data analysis ###
        #######################################################
        try:
            framedata = self.generator(myfile, frequency)
            self._draw_next_frame(framedata)
            track.tracking(self.file,frequency)

        except StopIteration:
            return False

        ##########################################################################
        ### if the resize limit has been reached, resize and redraw the figure ###
        ##########################################################################
        if self.resize_limit is not None:
            if self.file == self.resize_limit:

                # Dont redraw if this is the already the last file #

                if self.resize_limit < numFiles:

                    ###############################################################
                    ### If this is the last frequency, move onto the next limit ###
                    ###############################################################
                    if self.count == self.frequency_limit:
                        self.resize_limit = self.resize_limit + self.resize_interval

                        ### If the resize limit is above the number of files (e.g.
                        ### going out of bounds for the last resize event) then
                        ### readjust the final interval to the number of files
                        if self.resize_limit >= numFiles:
                            self.resize_limit = numFiles

                    ############################################################
                    ### 'if' statement used to make sure the plots dont get  ###
                    ### erased when there are no more files to be visualized ###
                    ############################################################
                    try:
                        self._redraw_figures()
                    except:
                        print('\nCould not redraw figure\n')



        ##################################################################
        ### If the function has analyzed each frequency for this file, ###
        ### move onto the next file and reset the frequency index      ###
        ##################################################################
        if self.count == self.frequency_limit:



            #########################################################################
            ### If the function has analyzed the final final, remove the callback ###
            #########################################################################
            if self.file == numFiles:
                print('\n%sFILE %s.\n%sElectrode %d\n%sData Analysis Complete\n' % (self.spacer,str(self.file),self.spacer,self.electrode,self.spacer))

                #post_analysis._analysis_finished()

            else:
                self.file += 1
                self.index += 1
                print('%smoving onto file %s\n' % (self.spacer,str(self.file)))
                self.count = 0
                root.after(1, self._step)



        ##########################################################
        ### Elif the function has not analyzed each frequency  ###
        ### for this file, move onto the next frequency        ###
        ##########################################################
        elif self.count < self.frequency_limit:
            self.count += 1

            root.after(1, self._step)


    def _raw_generator(self, myfile, frequency):

        ########################################
        ### Polynomical Regression Range (V) ###
        ########################################
        #--- if the frequency is equal or below cutoff_frequency, use the low freq parameters ---#
        if frequency <= cutoff_frequency:
            xstart = low_xstart
            xend = low_xend

        #--- if the frequency is above cutoff_frequency, use the high freq parameters ---#
        else:
            xstart = high_xstart
            xend = high_xend

        ###################################
        ### Retrieve data from the File ###
        ###################################
        potentials, currents, data_dict = ReadData(myfile, self.electrode)
        print(myfile)

        cut_value = 0
        for value in potentials:
            if value == 0:
                cut_value += 1

        if cut_value > 0:
            potentials = potentials[:-cut_value]
            currents = currents[:-cut_value]


        ################################################################
        ### Adjust the potentials depending on user-input parameters ###
        ################################################################
        adjusted_potentials = [value for value in potentials if xend <= value <= xstart]

        #########################################
        ### Savitzky-Golay Smoothing          ###
        #########################################
        min_potential = min(potentials)            # find the min potential
        sg_limit = sg_window/1000                  # mV --> V

        # shift all values positive
        sg_potentials = [x - min_potential for x in potentials]

        # find how many points fit within the sg potential window
        # this will be how many points are included in the rolling average
        sg_range = len([x for x in sg_potentials if x <= sg_limit])

        #--- Savitzky-golay Window must be greater than the range ---#
        if sg_range <= sg_degree:
            sg_range = sg_degree + 1

        #-- if the range is even, make it odd --#
        if sg_range % 2 == 0:
            sg_range = sg_range + 1

        # Apply the smoothing function and create a dictionary pairing
        # each potential with its corresponding current
        try:
            smooth_currents = savgol_filter(currents, sg_range, sg_degree)
            data_dict = dict(zip(potentials,smooth_currents))
        except ValueError:
            smooth_currents = savgol_filter(currents, 15, sg_degree)
            data_dict = dict(zip(potentials,smooth_currents))


        #######################################
        ### adjust the smooth currents to   ###
        ### match the adjusted potentials   ###
        #######################################
        adjusted_currents = []
        for potential in adjusted_potentials:
            adjusted_currents.append(data_dict[potential])

        ######################
        ### Polynomial fit ###
        ######################
        polynomial_coeffs = np.polyfit(adjusted_potentials,adjusted_currents,polyfit_deg)

        #############################
        ### Polynomial Regression ###
        #############################
        eval_regress = np.polyval(polynomial_coeffs,adjusted_potentials).tolist()
        regression_dict = dict(zip(eval_regress, adjusted_potentials))      # dictionary with current: potential

        ###############################################
        ### Absolute Max/Min Peak Height Extraction ###
        ###############################################
        #-- If the user selects 'Absolute Max/Min' in the 'Peak Height Extraction Settings'
        #-- within the Settings toolbar this analysis method will be used for PHE
        fit_half = round(len(eval_regress)/2)

        min1 = min(eval_regress[:fit_half])
        min2 = min(eval_regress[fit_half:])
        max1 = max(eval_regress[:fit_half])
        max2 = max(eval_regress[fit_half:])


        ################################################################
        ### If the user selected Peak Height Extraction, analyze PHE ###
        ################################################################

        if SelectedOptions == 'Peak Height Extraction':
            Peak_Height = max(max1,max2)-min(min1,min2)
            data = Peak_Height


        ########################################################
        ### If the user selected AUC extraction, analyze AUC ###
        ########################################################

        elif SelectedOptions == 'Area Under the Curve':
            ##################################
            ### Integrate Area Under the   ###
            ### Curve using a Riemmann Sum  ###
            ##################################
            AUC_index = 1
            AUC = 0

            AUC_potentials = [abs(potential) for potential in adjusted_potentials]

            #--- Find the minimum value and normalize it to 0 ---#
            AUC_min = min(adjusted_currents)
            AUC_currents = [Y - AUC_min for Y in adjusted_currents]

            #--- Midpoint Riemann Sum ---#
            while AUC_index <= len(AUC_currents) - 1:
                AUC_height = (AUC_currents[AUC_index] + AUC_currents[AUC_index - 1])/2
                AUC_width = AUC_potentials[AUC_index] - AUC_potentials[AUC_index - 1]
                AUC += (AUC_height * AUC_width)
                AUC_index += 1

            data = AUC

        #######################################
        ### Save the data into global lists ###
        #######################################

        data_list[self.num][self.count][self.index] = data

        frequency = frequency_list[self.count]
        self.frequency_axis.append(int(frequency))

        charge = (data/frequency) * 100000
        self.charge_axis.append(Peak_Height/frequency)

        #####################################################
        ### Return data to the animate function as 'args' ###
        #####################################################

        return potentials, adjusted_potentials, smooth_currents, adjusted_currents, eval_regress


    def _animate(self, framedata, *args):

        if key > 0:
            while True:

                potentials, adjusted_potentials, smooth_currents, adjusted_currents, regression = framedata

                print('\n%s%d: %dHz\n%s_animate' % (self.spacer,self.electrode,frequency_list[self.count],self.spacer))


                ################################################################
                ### Acquire the artists for this electrode at this frequency ###
                ### and get the data that will be visualized                 ###
                ################################################################
                plots = plot_list[self.num]

                ##########################
                ### Visualize the data ###
                ##########################

                #--- Peak Height ---#

                data = data_list[self.num][self.count][:len(self.file_list)]                     # 'num' is the electrode index value

                ####################################################
                ### Set the data of the artists to be visualized ###
                ####################################################
                plots[0].set_data(potentials,smooth_currents)          # Smooth current voltammogram
                plots[1].set_data(adjusted_potentials,regression)      # Regression voltammogram
                plots[2].set_data(self.frequency_axis,self.charge_axis)


                print('returning plots!')
                return plots


        else:
            file = 1
            EmptyPlots = framedata
            time.sleep(0.1)
            print('\n Yielding Empty Plots in Animation \n')
            return EmptyPlots




                                        ##############################
                                        ##############################
                                        ### END OF ANIMATION CLASS ###
                                        ##############################
                                        ##############################



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



                        ###############################################################################
                        ###############################################################################
                        ###### Classes and Functions for Real-Time Tracking and Text File Export ######
                        ###############################################################################
                        ###############################################################################

class Track():
    def __init__(self):

        self.track_list = [1]*numFiles

    def tracking(self, file, frequency):

        index = file - 1

        if self.track_list[index] == electrode_count:

            if SaveVar:
                text_file_export.ExportData(file,frequency)


            self.track_list[index] = 1

        else:
            self.track_list[index] += 1



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#



                                        ################################################
                                        ################################################
                                        ### Functions for Real Time Text File Export ###
                                        ################################################
                                        ################################################




##################################
### Real-Time Text File Export ###
##################################
class TextFileExport():

    ###############################
    ### Initialize the .txt file ###
    ###############################
    def __init__(self):
        global tracker

        self.TextFileHandle = FileHandle
        tracker = 0

        TxtList = []
        TxtList.append('Frequency(Hz)')

        E_count = 1
        for electrode in frame_list:
            if SelectedOptions == 'Peak Height Extraction':
                TxtList.append('PeakHeight_E%d(uA)' % (E_count))
            elif SelectedOptions == 'Area Under the Curve':
                TxtList.append('AUC_E%d' % (E_count))
            TxtList.append('Charge_E%d(uC)' % (E_count))
            E_count += 1

        TxtList.append('Avg.PeakHeight(uA)')
        TxtList.append('Standard_Deviation')
        TxtList.append('Avg.Charge(uC)')

        with open(FileHandle,'w+',encoding='utf-8', newline = '') as input:
            writer = csv.writer(input, delimiter = ' ')
            writer.writerow(TxtList)

    #################################################################
    ### Write the data from the current file into the Export File ###
    #################################################################
    def ExportData(self, file,frequency):

        list = []
        index = file - 1
        try:

            list.append(str(frequency))
            count = frequency_dict[int(frequency)]

            # Peak Height / AUC
            for num in range(electrode_count):
                list.append(data_list[num][count][index])
                list.append(data_list[num][count][index]/int(frequency))

            # Average Peak Height / AUC
            value = 0
            for num in range(electrode_count):
                value += data_list[num][count][index]

            average = value/electrode_count
            list.append(average)

            # Standard Deviation of a Sample across all electrodes
            std_list = []
            for num in range(electrode_count):
                std_list.append(data_list[num][count][index])

            std_list = [(value - average)**2 for value in std_list]

            standard_deviation = sqrt(sum(std_list)/(electrode_count - 1))
            list.append(standard_deviation)
            list.append(average/int(frequency))


            #--- Write the data into the .txt file ---#
            with open(FileHandle,'a',encoding='utf-8', newline = '') as input:
                writer = csv.writer(input, delimiter = ' ')
                writer.writerow(list)
            with open(FileHandle,'r',encoding='utf-8', newline = '') as filecontents:
                filedata =  filecontents.read()
            filedata = filedata.replace('[','')
            filedata = filedata.replace('"','')
            filedata = filedata.replace(']','')
            filedata = filedata.replace(',','')
            filedata = filedata.replace('\'','')
            with open(FileHandle,'w',encoding='utf-8', newline = '') as output:
                output.write(filedata)


        except:
            print('\n\n','ERROR IN TEXT FILE EXPORT','\n\n')
            time.sleep(3)


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#


                    ############################################
                    ### Initialize GUI to start the program  ###
                    ############################################

if __name__ == '__main__':

    root = tk.Tk()
    app = MainWindow(root)
    style = ttk.Style()
    style.configure('On.TButton', foreground = 'blue', font = LARGE_FONT, relief = 'raised', border = 100)
    style.configure('Off.TButton', foreground = 'black', font = MEDIUM_FONT, relief = 'sunken', border = 5)


    while True:
        #--- initiate the mainloop ---#
        try:
            root.mainloop()
        #--- escape scrolling error ---#
        except UnicodeDecodeError:
            pass


                    #*########################################*#
                    #*############ End of Program ############*#
                    #*########################################*#
