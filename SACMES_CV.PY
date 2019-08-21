
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
import scipy
from scipy import *
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
handle_variable = 'CV_repetitive_100mVs'    # default handle variable is nothing
e_var = 'single'        # default input file is 'Multichannel', or a single file containing all electrodes
PHE_method = 'Abs'      # default PHE Extraction is difference between absolute max/min

#------------------------------------------------------------#

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


#############################
### Checkpoint Parameters ###
#############################
key = 0                 ### SkeletonKey
search_lim = 15         ### Search limit (sec)
PoisonPill = False      ### Stop Animation variable
FoundFilePath = False   ### If the user-inputted file is found
ExistVar = False        ### If Checkpoints are not met ExistVar = True
AlreadyInitiated = False    ### indicates if the user has already initiated analysis
AlreadyReset = False
ExportAltered = False   ### Indicates if user has manually altered the export path

analysis_complete = False    ### If analysis has completed, begin PostAnalysis

##################################
### Data Extraction Parameters ###
##################################
delimiter = 1               ### default delimiter is a space; 2 = tab
column_index = -2           ### column index for list_val.
                              # list_val = column_index + 3
                              # defauly column is the second (so index = 1)
mass_transport = 'surface'

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



############################
#### Update Global Lists ###
############################
def _update_global_lists(file):
    global file_list

    if file not in file_list:
        file_list.append(file)

        if file != numFiles:
            FileLabel['text'] = file + 1

#######################################
### Retrieve the column index value ###
#######################################
def _get_listval(electrode):

    if e_var == 'single':
        list_val = electrode

    elif e_var == 'multiple':
        list_val = 1

    return list_val

##############################
### Retrieve the file name ###
##############################
def _retrieve_file(file, electrode):

    if e_var == 'single':
        filename = '%s_%d.txt' % (handle_variable, file)
        filename2 = '%s__%d.txt' % (handle_variable, file)

    elif e_var == 'multiple':
        filename = 'E%s_%s_%d.txt' % (electrode,handle_variable,file)
        filename2 = 'E%s_%s__%d.txt' % (electrode,handle_variable,file)

    return filename, filename2


def ReadData(myfile, electrode):
    try:
        ###############################################################
        ### Get the index value of the data depending on if the     ###
        ### electrodes are in the same .txt file or separate files  ###
        ###############################################################
        list_val = _get_listval(electrode)

        #####################
        ### Read the data ###
        #####################

        #---Preallocate Potential and Current lists---#
        with open(myfile,'r',encoding='utf-8') as mydata:
            variables = len(mydata.readlines())
            potentials = ['hold']*variables

            ### key: potential; value: current ##
            currents = [0]*variables
            correlation_list = []

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
                    list_num = list_num + 1



        ### if there are 0's in the list (if the preallocation added to many)
        ### then remove them
        cut_value = 0
        for value in potentials:
            if value == 'hold':
                cut_value += 1
        if cut_value > 0:
            potentials = potentials[:-cut_value]
            currents = currents[:-cut_value]

        min_potential = min(potentials)
        max_potential = max(potentials)


        ################################################################
        ### Determine the Number of Segments and Create a dictionary ###
        ### correlating the segments to their index values within    ###
        ### text file                                                ###
        ################################################################
        list_num = 0            # index value
        step_list = []          # step potential list
        segment_number = 1      # segment number
        segment_dictionary = {}
        segment_dictionary[segment_number] = ([],{},{})
        forward_dictionary = {}    # dictionary for forward scans
        reverse_dictionary = {}     # dictionary for reverse scans

        #-- find the initial step potential to determine the  --#
        #-- direction for forward and reverse segments       --#
        initial_step = potentials[1] - potentials[0]
        if initial_step < 0:
            initial_step = 'Negative'
        else:
            initial_step = 'Positive'


        #-- place the potentials and currents into dictionaries --#
        #-- correlated to their respective segment              --#
        for value in potentials:
            if list_num > 0:
                step = potentials[list_num] - potentials[list_num-1]

                if step < 0:
                    step = 'Negative'
                    step_list.append(step)


                else:
                    step = 'Positive'
                    step_list.append(step)

            if list_num > 1:
                #-- if the step potentials do not match, the sign has  --#
                #-- reversed and a new segment has begun               --#

                #-- the step_list is one index behind the segment dictionary --#
                if step_list[list_num-1] != step_list[list_num-2]:

                    segment_number += 1
                    segment_dictionary[segment_number] = ([],{},{})

                    #-- append the inflection point to both segments --#
                    segment_dictionary[segment_number][0].append(list_num-1)
                    segment_dictionary[segment_number][1][potentials[list_num-1]] = list_num-1
                    segment_dictionary[segment_number][2].setdefault(currents[list_num-1],[]).append(value)

                    if step == initial_step:
                        if segment_number not in forward_dictionary:
                            forward_dictionary[segment_number] = {}
                        forward_dictionary[segment_number][potentials[list_num-1]] = list_num-1

                    else:
                        if segment_number not in reverse_dictionary:
                            reverse_dictionary[segment_number] = {}
                        reverse_dictionary[segment_number][potentials[list_num-1]] = list_num-1


                if step == initial_step:
                    if segment_number not in forward_dictionary:
                        forward_dictionary[segment_number] = {}
                    forward_dictionary[segment_number][value] = list_num

                else:
                    if segment_number not in reverse_dictionary:
                        reverse_dictionary[segment_number] = {}
                    reverse_dictionary[segment_number][value] = list_num

            segment_dictionary[segment_number][0].append(list_num)
            segment_dictionary[segment_number][1][value] = list_num
            segment_dictionary[segment_number][2].setdefault(currents[list_num],[]).append(value)

            list_num += 1

        #-- Return the data --#
        return potentials, currents, segment_dictionary, forward_dictionary, reverse_dictionary

    except:
        print('\n\nError in Read Data\n\n')

def ClosestVoltageEstimation(voltage, potentials):

    ###################################################################
    ### Find the closest potential value to the user-inputted value ###
    ###################################################################
    # shift all potentials positive with the minimum potential = 0 #
    min_voltage = min(potentials)
    selection_potentials = [x - min_voltage for x in potentials]
    adjusted_voltage = voltage - min_voltage

    value = 0
    selection_dict = {}
    for x in selection_potentials:
        selection = x - adjusted_voltage
        selection_dict[selection] = value
        value += 1

    #-- Voltage 1 --#
    Voltage_Index = selection_dict[min(selection_dict, key=abs)]
    Voltage = potentials[Voltage_Index]

    return Voltage

#######################################################################
### Function to raise frame of column 0.                            ###
### (RealTimeManipulationFrame,InputFrame,SetupFrame,PostAnalysis)  ###
#######################################################################
def show_frame(cont):

    frame = ShowFrames[cont]
    frame.tkraise()

######################################################
### Function to raise frame for specific electrode ###
######################################################
def show_plot(frame):
    frame.tkraise()


#####################################
### Destory the frames on Reset() ###
#####################################
def close_frame(cont, destroy):
    frame = ShowFrames[cont]
    frame.grid_forget()

    if destroy:
        close_visualization()

def close_visualization():

    # close all matplotlib figures
    plt.close('all')

    # destory the frames holding the figures
    for frame in PlotValues:
        frame.destroy()

    # destory the container holding those frames
    PlotContainer.destroy()


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

        self._create_toolbar(container)

    #--- Function to visualize different frames ---#
    def _create_toolbar(self, parent):
        global forward_boolean,reverse_boolean

        menubar = tk.Menu(self.master)
        self.master.config(menu=menubar)
        self.parent = parent

        #################
        ### Edit Menu ###
        #################

        self.transport_menu = tk.Menu(menubar)
        self.solution = self.transport_menu.add_command(label = "    Solution Phase", command = lambda: self.mass_transport('solution'))
        self.surface = self.transport_menu.add_command(label="✓ Surface Bound", command = lambda: self.mass_transport('surface'))
        menubar.add_cascade(label='Mass Transport', menu=self.transport_menu)

        self.visualization_menu = tk.Menu(menubar)

        forward_boolean = tk.BooleanVar()
        forward_boolean.set(True)
        self.visualize_forward = self.visualization_menu.add_checkbutton(label='Forward Segment',onvalue=True,offvalue=False,variable=forward_boolean)

        reverse_boolean = tk.BooleanVar()
        reverse_boolean.set(True)
        self.visualize_reverse = self.visualization_menu.add_checkbutton(label='Reverse Segment',onvalue=True,offvalue=False,variable=reverse_boolean)
        menubar.add_cascade(label='Segment Visualization',menu=self.visualization_menu)

    def mass_transport(self, transport):
        global mass_transport

        mass_transport = transport

        # Change the labels of the menu options #
        if mass_transport == 'solution':
            self.transport_menu.entryconfigure(0, label='✓ Solution Phase')
            self.transport_menu.entryconfigure(1,label='    Surface Bound')
        elif mass_transport == 'surface':
            self.transport_menu.entryconfigure(0, label='    Solution Phase')
            self.transport_menu.entryconfigure(2,label='✓ Surface Bound')

        print('\nmass transport value = %s' % mass_transport)


    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()

    def onExit(self):
        self.master.destroy()
        self.master.quit()
        quit()

class InputFrame(tk.Frame):                         # first frame that is displayed when the program is initialized
    def __init__(self, parent, controller):
        global figures, SaveBox

        self.parent = parent
        self.controller = controller

        tk.Frame.__init__(self, parent)             # initialize the frame

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
        numFileLabel.grid(row=5,column=0,columnspan=2,pady=4)
        self.numfiles = ttk.Entry(self, width=7)
        self.numfiles.insert(END, '20')
        self.numfiles.grid(row=6,column=0,columnspan=2,pady=6)

        #--- Analysis interval for event callback in ElectrochemicalAnimation ---#
        IntervalLabel = tk.Label(self, text='Analysis Interval (ms):', font=LARGE_FONT)
        IntervalLabel.grid(row=5,column=2,columnspan=2,pady=4)
        self.Interval = ttk.Entry(self, width=7)
        self.Interval.insert(END, '10')
        self.Interval.grid(row=6,column=2,columnspan=2,pady=6)

        #---Scan Rate Variable---#
        ScanLabel = tk.Label(self, text='Scan Rate (V/s):', font=LARGE_FONT)
        ScanLabel.grid(row=7,column=0,columnspan=2)
        self.scan_rate = ttk.Entry(self, width=7)
        self.scan_rate.insert(END, '0.1')
        self.scan_rate.grid(row=8,column=0,columnspan=2)

        self.resize_label = tk.Label(self, text='Resize Interval', font = LARGE_FONT)
        self.resize_label.grid(row=7,column=2,columnspan=2)
        self.resize_entry = tk.Entry(self, width = 7)
        self.resize_entry.insert(END,'200')
        self.resize_entry.grid(row=8,column=2,columnspan=2)

        ##################################
        ### Select and Edit Electrodes ###
        ##################################

        self.ElectrodeListboxFrame = tk.Frame(self)                   # create a frame to pack in the Electrode box and
        self.ElectrodeListboxFrame.grid(row=10,column=0,columnspan=4,padx=10,pady=10, sticky = 'ns')

        #--- parameters for handling resize ---#
        self.ElectrodeListboxFrame.rowconfigure(0, weight=1)
        self.ElectrodeListboxFrame.rowconfigure(1, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(0, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(1, weight=1)

        electrodes = [1,2,3,4,5,6,7,8]
        self.ElectrodeListExists = False
        self.ElectrodeLabel = tk.Label(self.ElectrodeListboxFrame, text='Select Electrodes:', font=LARGE_FONT)
        self.ElectrodeLabel.grid(row=0,column=0,columnspan=2, sticky = 'nswe')
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
        self.SingleElectrodeFile = ttk.Button(self.ElectrodeListboxFrame, text='Multichannel', style = 'On.TButton', command =  lambda: self.ElectrodeSelect('Multichannel'))
        self.SingleElectrodeFile.grid(row=2,column=0)

        #--- Option to have data for each electrode in a separate file ---#
        self.MultipleElectrodeFiles = ttk.Button(self.ElectrodeListboxFrame,  text='Multiplex', style = 'Off.TButton',command = lambda: self.ElectrodeSelect('Multiplex'))
        self.MultipleElectrodeFiles.grid(row=2,column=1)

        #--- Select Analysis Method---#
        Options = ['Peak Height Extraction','Area Under the Curve']
        OptionsLabel = tk.Label(self, text='Select Data to be Plotted', font=LARGE_FONT)
        self.PlotOptions = Listbox(self, relief='groove', exportselection=0, font=LARGE_FONT, height=len(Options), selectmode='single', bd=3)
        self.PlotOptions.bind('<<ListboxSelect>>', self.SelectPlotOptions)
        OptionsLabel.grid(row=11,column=0,columnspan=4)
        self.PlotOptions.grid(row=12,column=0,columnspan=4)
        for option in Options:
            self.PlotOptions.insert(END, option)

        #--- Warning label for if the user does not select an analysis method ---#
        self.NoOptionsSelected = tk.Label(self, text = 'Select a Data Analysis Method', font = MEDIUM_FONT, fg='red')   # will only be added to the grid (row 16) if they dont select an option
        self.NoSelection = False

        ############################################################
        ### Adjustment of Visualization Parameters: xstart, xend ###
        ############################################################

        #--- Ask the User if they want to export the data to a .txt file ---#
        self.SaveVar = BooleanVar()
        self.SaveVar.set(False)
        self.SaveBox = Checkbutton(self, variable=self.SaveVar, onvalue=True, offvalue=False, text="Export Data").grid(row=17,column=0,columnspan=4,pady=10)

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
        self.ElectrodeListboxFrame.tkraise()

    def ElectrodeSelect(self, variable):
        global e_var

        if variable == 'Multiplex':
            e_var = 'multiple'

            self.SingleElectrodeFile['style'] = 'Off.TButton'
            self.MultipleElectrodeFiles['style'] = 'On.TButton'

        elif variable == 'Multichannel':
            e_var = 'single'

            self.SingleElectrodeFile['style'] = 'On.TButton'
            self.MultipleElectrodeFiles['style'] = 'Off.TButton'

    #--- Analysis Method ---#
    def SelectPlotOptions(self, evt):
        global SelectedOptions
        SelectedOptions = str((self.PlotOptions.get(self.PlotOptions.curselection())))


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

    #--- Electrode Selection ---#
    def ElectrodeCurSelect(self, evt):
        ###################################################
        ## electrode_list: list; ints                    ##
        ## electrode_dict: dict; {electrode: index}      ##
        ## electrode_count: int                          ##
        ###################################################
        global electrode_count, electrode_list, electrode_dict, frame_list, PlotValues

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
        global mypath, Option, SelectedOptions, FileHandle, AlreadyInitiated, delimeter

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


        if not self.ElectrodeListExists:
            self.ElectrodeLabel['fg'] = 'red'
        elif self.ElectrodeListExists:
            self.ElectrodeLabel['fg'] = 'black'


        if not self.PathWarningExists:
            if not self.NoSelection:
                self.StartProgram()

        else:
            print('Could Not Start Program')


    ########################################################################
    ### Function To Initialize Data Acquisition, Analysis, and Animation ###
    ########################################################################

    def StartProgram(self):
        global FileHandle, text_file_export, starting_file, scan_rate, handle_variable, track, Interval, e_var, resize_interval, CV_min, CV_max, data_min, data_max, mypath, electrode_count, SaveVar, track, numFiles, frames, generate, figures, Plot, frame_list, PlotValues, anim, q, delimiter

        #---Get the User Input and make it globally accessible---#

        numFiles = int(self.numfiles.get())     # file limit

        q = Queue()

        ### Set the delimiter value for data columns ###
        if delimiter == 1:
            delimiter = ' '
        elif delimiter == 2:
            delimiter = '   '

        starting_file = 1

        SaveVar = self.SaveVar.get()                        # tracks if text file export has been activated
        resize_interval = int(self.resize_entry.get())      # interval at which xaxis of plots resizes
        handle_variable = self.ImportFileEntry.get()        # string handle used for the input file
        scan_rate = float(self.scan_rate.get())

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

        row_value = 1
        self.label_dict = {}
        self.already_verified = {}
        for electrode in electrode_list:
            electrode_label = tk.Label(self.win, text = 'E%s' % electrode,fg='red',font=LARGE_FONT)
            electrode_label.grid(row=row_value,column=0,pady=5,padx=5)
            frame = tk.Frame(self.win, relief='groove',bd=5)
            frame.grid(row = row_value,column=1,pady=5,padx=5)
            self.label_dict[electrode] = electrode_label
            self.already_verified[electrode] = {}
            row_value += 1

        self.stop = tk.Button(self.win, text = 'Stop', command = self.stop)
        self.stop.grid(row=row_value, column=0,columnspan=2,pady=5)
        self.StopSearch = False

        self.num = 0
        self.analysis_count = 0
        self.analysis_limit = electrode_count
        self.electrode_limit = electrode_count - 1

        root.after(50,self.verify)

    def verify(self):

        self.electrode = electrode_list[self.num]

        filename, filename2 = _retrieve_file(1,self.electrode)

        myfile = mypath + filename               ### path of your file
        myfile2 = mypath + filename2               ### path of your file

        try:
            print('Looking for %s' % myfile)
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

        except:
            try:
                print('Looking for %s' % myfile2)
                mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                myfile = myfile2
            except:
                mydata_bytes = 1

        if mydata_bytes > 1000:

            #-- if in multipotentiostat settings, try and see if a column
            #-- exists for the current electrode
            if e_var == 'single':
                list_val = _get_listval(self.electrode)
                try:
                    potentials, currents, segment_dictionary, forward, reverse = ReadData(myfile, self.electrode)

                    if not self.already_verified[self.electrode]:
                        self.already_verified[self.electrode] = True
                        if not self.StopSearch:
                            self.label_dict[self.electrode]['fg'] = 'green'
                            self.analysis_count += 1

                except:
                    if not self.StopSearch:
                        root.after(200,self.verify)

            elif e_var == 'multiple':
                if not self.already_verified[self.electrode]:
                    self.already_verified[self.electrode] = True
                    if not self.StopSearch:
                        self.label_dict[self.electrode]['fg'] = 'green'
                        self.analysis_count += 1

            if self.analysis_count == self.analysis_limit:
                root.after(10,self.proceed)


        if self.num < self.electrode_limit:
            self.num += 1
        else:
            self.num = 0

        if self.analysis_count < self.analysis_limit:
            if not self.StopSearch:
                root.after(10,self.verify)

    def proceed(self):
        global track, initialize

        self.win.destroy()

        ##############################
        ### Syncronization Classes ###
        ##############################
        track = Track()

        ######################################################
        ### Matplotlib Canvas, Figure, and Artist Creation ###
        ######################################################
        initialize = InitializeFigureCanvas()


        ############################################################
        ### Initialize the Initial Parameter Visualization Frame ###
        ############################################################
        frame = SetupVisualizationFrame(container, self)
        ShowFrames[SetupVisualizationFrame] = frame
        frame.grid(row=0,column=1,sticky='nsew')

        ####################################################
        ### Initialize the Initial Parameter Setup Frame ###
        ####################################################
        frame = SetupFrame(container, self)
        ShowFrames[SetupFrame] = frame
        frame.grid(row=0, column=0, sticky='nsew')


        #---When initliazed, raise the Start Page and the plot for electrode one---#
        self.show_frame(SetupFrame)                  # raises the frame for initial parameter setup
        self.show_frame(SetupVisualizationFrame)     # raises the setup visualization frame


    def stop(self):
        self.StopSearch = True
        self.win.destroy()

    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame):
        frame.tkraise()

    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()


#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#

class SetupFrame(tk.Frame):
    def __init__(self, parent, controller):
        global low_voltage,high_voltage,CV_min,CV_max,data_min,data_max,V1_min,V1_max,V2_min,V2_max

        self.parent = parent
        self.controller = controller

        tk.Frame.__init__(self, parent)

        self.potentials = setup_dictionary['potentials']
        self.currents = setup_dictionary['currents']
        self.segment_dictionary = setup_dictionary['segments']
        self.forward = setup_dictionary['forward']
        self.reverse = setup_dictionary['reverse']
        self.potential_dict = setup_dictionary['self.potential_dict']
        self.fig, self.ax = setup_dictionary['figure']


        self.max_voltage = max(self.potentials)
        self.min_voltage = min(self.potentials)

        self.reverse_artist, = self.ax.plot([],[],'o',color='0.75',MarkerSize=1,label='Reverse',zorder=1)
        self.adjusted_reverse_artist, = self.ax.plot([],[],'k-',MarkerSize=0.1,zorder=2)
        self.adjusted_reverse_artist.set(alpha=0.9)          # increase its transparency

        self.forward_artist, = self.ax.plot([],[],'o',color='0.75',MarkerSize=1,label='Forwards',zorder=1)
        self.adjusted_forward_artist, = self.ax.plot([],[],'k-',MarkerSize=0.1,zorder=2)
        self.adjusted_forward_artist.set(alpha=0.9)        # increase its transparency

        self.setup_artists = [self.reverse_artist,self.adjusted_reverse_artist,self.forward_artist,self.adjusted_forward_artist]

        self.bg_cache = {}

        if not AlreadyReset:
            low_voltage = self.min_voltage
            high_voltage = self.max_voltage


        #############
        ### Title ###
        #############
        self.Title = tk.Label(self,text = 'Initial Parameter Setup', font=LARGE_FONT).grid(row=0,column=0,columnspan=2,pady=10)

        if text_file_export:
            self.AlterExportPathButton = ttk.Button(self, text = 'Alter Export Path', style = 'Off.TButton',command = lambda: self.AlterExportPath())
            self.AlterExportPathButton.grid(row=1,column=0,columnspan=2,pady=10)

        ##############################
        ### Forwards Segment Setup ###
        ##############################
        self.ForwardsSegmentLabel = tk.Label(self, text = 'Select Forward\nSegment').grid(row=2,column=0,padx=10,pady=5)
        self.ForwardsSegmentSelection = Listbox(self, relief='groove', exportselection=0, width=10, font=LARGE_FONT, height=6, selectmode = 'single', bd=3)
        self.ForwardsSegmentSelection.bind('<<ListboxSelect>>',self.ForwardsSelection)
        self.ForwardsSegmentSelection.grid(row=3,column=0,padx=10,sticky='nswe')
        self.forward_selection_dict = {}
        count = 0
        for segment in self.forward:
            self.forward_selection_dict[segment] = count
            self.ForwardsSegmentSelection.insert(END, segment)
            count += 1

        # row 4, column 1
        self.NoForwardsSelection = tk.Label(self, text = 'Select Forward Segment',fg='red')
        self.ForwardsSelectionExists = False
        self.ForwardsWarning = False

        ############################
        ### Reverse Segment Setup ###
        ############################
        self.ReverseSegmentLabel = tk.Label(self, text = 'Select Reverse\nSegment').grid(row=2,column=1,padx=10,pady=5)
        self.ReverseSegmentSelection = Listbox(self, relief='groove', exportselection=0, width=10, font=LARGE_FONT, height=6, selectmode = 'single', bd=3)
        self.ReverseSegmentSelection.bind('<<ListboxSelect>>',self.ReverseSelection)
        self.ReverseSegmentSelection.grid(row=3,column=1,padx=10,sticky='nswe')
        self.reverse_selection_dict = {}
        count = 0
        for segment in self.reverse:
            self.reverse_selection_dict[segment] = count
            self.ReverseSegmentSelection.insert(END, segment)
            count += 1

        # row 4, column 0
        self.NoReverseSelection = tk.Label(self, text = 'Select Reverse Segment',fg='red')
        self.ReverseSelectionExists = False
        self.ReverseWarning = False



        ###############################################
        ### Setup Frame for Voltage Range Selection ###
        ###############################################

        #--- Container Frame ---#
        self.VoltageSelectionFrame = tk.Frame(self, relief = 'groove',bd=4)
        self.VoltageSelectionFrame.grid(row=5,column=0,columnspan=2,pady=10)
        #-----------------------#

        #--- Reverse Voltage Adjustment ---#
        self.ReverseFrame = tk.Frame(self.VoltageSelectionFrame)
        self.ReverseFrame.grid(row=1,column=0,columnspan=2,sticky='news')

        self.LowVoltageLabel = tk.Label(self.ReverseFrame, text = 'Low Voltage', font=MEDIUM_FONT).grid(row=0,column=0,pady=2,padx=2)
        self.LowVoltage = tk.Entry(self.ReverseFrame, width=7)
        self.LowVoltage.insert(END,low_voltage)
        self.LowVoltage.grid(row=1,column=0,padx=2,pady=2)

        self.HighVoltageLabel = tk.Label(self.ReverseFrame, text = 'High Voltage', font=MEDIUM_FONT).grid(row=0,column=1,pady=2,padx=2)
        self.HighVoltage = tk.Entry(self.ReverseFrame, width=7)
        self.HighVoltage.insert(END,high_voltage)
        self.HighVoltage.grid(row=1,column=1,padx=2,pady=2)

        self.ApplyVoltages = ttk.Button(self.VoltageSelectionFrame, text='Apply Voltages', command = lambda: self.ApplyVoltage())
        self.ApplyVoltages.grid(row=5,column=0,columnspan=2,pady=5,padx=5)

        ##################################################
        ### Frame to select two unique voltages to be  ###
        ### tracked over the course of the experiment  ###
        ##################################################

        #--- Container Frame ---#
        self.SelectUniqueVoltages = tk.Frame(self, relief='groove',bd=4)
        self.SelectUniqueVoltages.grid(row=6,column=0,columnspan=2,pady=10)

        self.SelectUniqueVoltagesTitle = tk.Label(self.SelectUniqueVoltages, text = 'Voltage Probes',font=LARGE_FONT).grid(row=0,column=0,columnspan=2,padx=5,pady=5)
        #-----------------------#

        self.Voltage1Label = tk.Label(self.SelectUniqueVoltages, text = 'Voltage 1', font=MEDIUM_FONT).grid(row=1,column=0,padx=5,pady=5)
        self.Voltage1 = tk.Entry(self.SelectUniqueVoltages, width=7)
        self.Voltage1.grid(row=2,column=0,padx=5,pady=5)

        self.Voltage1Error = tk.Label(self.SelectUniqueVoltages,text='Select A Voltage',fg='red',font=SMALL_FONT)
        self.Voltage1Warning = False

        if not AlreadyReset:
            self.Voltage1Exists = False
        else:
            self.Voltage1Exists = True

        self.Voltage2Label = tk.Label(self.SelectUniqueVoltages, text = 'Voltage 2', font=MEDIUM_FONT).grid(row=1,column=1,padx=5,pady=5)
        self.Voltage2 = tk.Entry(self.SelectUniqueVoltages, width=7)
        self.Voltage2.grid(row=2,column=1,padx=5,pady=5)

        self.Voltage2Error = tk.Label(self.SelectUniqueVoltages,text='Select A Voltage',fg='red',font=SMALL_FONT)
        self.Voltage2Warning = False
        if not AlreadyReset:
            self.Voltage2Exists = False
        else:
            self.Voltage2Exists = True



        if AlreadyReset:
            self.Voltage1.insert(END,Voltage1)
            self.Voltage2.insert(END,Voltage2)

        self.SelectVoltages = tk.Button(self.SelectUniqueVoltages, text = 'Select Voltages',font=MEDIUM_FONT, command = lambda: self.SelectVoltagesFunction())
        self.SelectVoltages.grid(row=5,column=0,columnspan=2,pady=5,padx=5)


        #--- Create a frame that will contain all of the widgets ---#
        AdjustmentFrame = tk.Frame(self, relief = 'groove', bd=3)
        AdjustmentFrame.grid(row=8,column=0,columnspan=2,pady=15)
        AdjustmentFrame.rowconfigure(0, weight=1)
        AdjustmentFrame.rowconfigure(1, weight=1)
        AdjustmentFrame.rowconfigure(2, weight=1)
        AdjustmentFrame.rowconfigure(3, weight=1)
        AdjustmentFrame.rowconfigure(4, weight=1)
        AdjustmentFrame.columnconfigure(0, weight=1)
        AdjustmentFrame.columnconfigure(1, weight=1)
        AdjustmentFrame.columnconfigure(2, weight=1)
        AdjustmentFrame.columnconfigure(3, weight=1)

        if not AlreadyReset:
            CV_min = 2
            CV_max = 2
            data_min = 1
            data_max = 1
            V1_min = 2
            V1_max = 2
            V2_min = 2
            V2_max = 2


        #--- Y Limit Adjustment Variables ---#
        self.y_limit_parameter_label = tk.Label(AdjustmentFrame, text = 'Select Y Limit Parameters',font=LARGE_FONT)
        self.y_limit_parameter_label.grid(row=0,column=0,columnspan=2,pady=5,padx=5)

        #--- Cyclic Voltammogram Minimum Parameter Adjustment ---#
        self.CV_min_parameter_label = tk.Label(AdjustmentFrame, text = 'Raw Min. Factor',font=MEDIUM_FONT)
        self.CV_min_parameter_label.grid(row=1,column=0)
        self.CV_min = tk.Entry(AdjustmentFrame, width=5)
        self.CV_min.insert(END, CV_min)                   # initial minimum is set to minimum current - (2*minimum_current) (baseline) of file 1
        self.CV_min.grid(row=2,column=0,padx=5,pady=2,ipadx=2)

        #--- Cyclic Voltammogram Parameter Adjustment ---#
        self.CV_max_parameter_label = tk.Label(AdjustmentFrame, text = 'Raw Max. Factor',font=MEDIUM_FONT)
        self.CV_max_parameter_label.grid(row=3,column=0)
        self.CV_max = tk.Entry(AdjustmentFrame, width=5)
        self.CV_max.insert(END, CV_max)                      # initial max is max_current + (2*max_currente)
        self.CV_max.grid(row=4,column=0,padx=5,pady=2,ipadx=2)

        #--- Raw Data Minimum Parameter Adjustment ---#
        self.data_min_parameter_label = tk.Label(AdjustmentFrame, text = 'Data Min. Factor',font=MEDIUM_FONT)
        self.data_min_parameter_label.grid(row=1,column=1)
        self.data_min = tk.Entry(AdjustmentFrame, width=5)
        self.data_min.insert(END, data_min)                   # initial minimum is set to minimum current - (2*minimum_current) (baseline) of file 1
        self.data_min.grid(row=2,column=1,padx=5,pady=2,ipadx=2)

        #--- Raw Data Maximum Parameter Adjustment ---#
        self.data_max_parameter_label = tk.Label(AdjustmentFrame, text = 'Data Max. Factor',font=MEDIUM_FONT)
        self.data_max_parameter_label.grid(row=3,column=1)
        self.data_max = tk.Entry(AdjustmentFrame, width=5)
        self.data_max.insert(END, data_max)                      # initial max is max_current + (2*max_currente)
        self.data_max.grid(row=4,column=1,padx=5,pady=2,ipadx=2)

        #--- Voltage 1 Minimum Parameter Adjustment ---#
        self.V1_min_parameter_label = tk.Label(AdjustmentFrame, text = 'V1 Min. Factor',font=MEDIUM_FONT)
        self.V1_min_parameter_label.grid(row=1,column=2)
        self.V1_min = tk.Entry(AdjustmentFrame, width=5)
        self.V1_min.insert(END, V1_min)                   # initial minimum is set to minimum current - (2*minimum_current) (baseline) of file 1
        self.V1_min.grid(row=2,column=2,padx=5,pady=2,ipadx=2)

        #--- Voltage 1 Maximum Parameter Adjustment ---#
        self.V1_max_parameter_label = tk.Label(AdjustmentFrame, text = 'V1 Max. Factor',font=MEDIUM_FONT)
        self.V1_max_parameter_label.grid(row=3,column=2)
        self.V1_max = tk.Entry(AdjustmentFrame, width=5)
        self.V1_max.insert(END, V1_max)                      # initial max is max_current + (2*max_currente)
        self.V1_max.grid(row=4,column=2,padx=5,pady=2,ipadx=2)

        #--- Voltage 2 Minimum Parameter Adjustment ---#
        self.V2_min_parameter_label = tk.Label(AdjustmentFrame, text = 'V2 Min. Factor',font=MEDIUM_FONT)
        self.V2_min_parameter_label.grid(row=1,column=3)
        self.V2_min = tk.Entry(AdjustmentFrame, width=5)
        self.V2_min.insert(END, V1_min)                   # initial minimum is set to minimum current - (2*minimum_current) (baseline) of file 1
        self.V2_min.grid(row=2,column=3,padx=5,pady=2,ipadx=2)

        #--- Voltage 2 Maximum Parameter Adjustment ---#
        self.V2_max_parameter_label = tk.Label(AdjustmentFrame, text = 'V2 Max. Factor',font=MEDIUM_FONT)
        self.V2_max_parameter_label.grid(row=3,column=3)
        self.V2_max = tk.Entry(AdjustmentFrame, width=5)
        self.V2_max.insert(END, V2_max)                      # initial max is max_current + (2*max_currente)
        self.V2_max.grid(row=4,column=3,padx=5,pady=2,ipadx=2)
        self.ResetButton = tk.Button(self, text = 'Reset', font=LARGE_FONT, width=7, command = lambda: self.Reset())
        self.ResetButton.grid(row=10,column=0,pady=5,padx=5)

        if AlreadyReset:
            self.ReverseSegmentSelection.select_set(self.reverse_selection_dict[reverse_segment]) #This only sets focus on the first item.
            self.ReverseSegmentSelection.event_generate("<<ListboxSelect>>")

            self.ForwardsSegmentSelection.select_set(self.forward_selection_dict[forward_segment])
            self.ForwardsSegmentSelection.event_generate("<<ListboxSelect>>")

            self.ReverseSelection(None)
            self.ForwardsSelection(None)

        #####################################
        ### Initialize Real Time Analysis ###
        #####################################
        self.Start = tk.Button(self, text='Start', font=LARGE_FONT, width=7, command = lambda: self.Checkpoint())
        self.Start.grid(row=10,column=1,pady=5,padx=5)

    def Reset(self):

        # Raise the initial user input frame
        show_frame(InputFrame)
        close_frame(SetupFrame, False)
        close_frame(SetupVisualizationFrame, False)


    def ReverseSelection(self, evt):
        global reverse_segment

        self.reverse_segment = int(self.ReverseSegmentSelection.get(self.ReverseSegmentSelection.curselection()))
        self.index_list, self.potential_dict, self.data_dict = self.segment_dictionary[self.reverse_segment]
        reverse_segment = self.reverse_segment

        self.reverse_potentials = [self.potentials[index] for index in self.index_list]
        self.reverse_currents = [self.currents[index] for index in self.index_list]

        self.adjusted_reverse_potentials = [potential for potential in self.reverse_potentials if low_voltage <= potential <= high_voltage]
        self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_reverse_potentials]
        self.adjusted_reverse_currents = [self.currents[index] for index in self.adjusted_index_list]

        self.adjusted_reverse_artist.set_data(self.adjusted_reverse_potentials,self.adjusted_reverse_currents)
        self.reverse_artist.set_data(self.reverse_potentials,self.reverse_currents)

        if not self.ReverseSelectionExists:
            self.ReverseSelectionExists = True

        if self.reverse_segment is not None:
            if self.Voltage1Exists:
                self.check_V1(Voltage1)

            if self.Voltage2Exists:
                self.check_V2(Voltage2)

        self.blit_data(self.setup_artists)


    def ForwardsSelection(self, evt):
        global forward_segment

        self.forward_segment = int(self.ForwardsSegmentSelection.get(self.ForwardsSegmentSelection.curselection()))
        forward_segment = self.forward_segment

        self.index_list, self.potential_dict,self.data_dict = self.segment_dictionary[self.forward_segment]

        self.forward_potentials = [self.potentials[index] for index in self.index_list]
        self.forward_currents = [self.currents[index] for index in self.index_list]

        self.adjusted_forward_potentials = [potential for potential in self.forward_potentials if low_voltage <= potential <= high_voltage]
        self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_forward_potentials]
        self.adjusted_forward_currents = [self.currents[index] for index in self.adjusted_index_list]

        self.adjusted_forward_artist.set_data(self.adjusted_forward_potentials,self.adjusted_forward_currents)
        self.forward_artist.set_data(self.forward_potentials,self.forward_currents)

        if not self.ForwardsSelectionExists:
            self.ForwardsSelectionExists = True

        if self.Voltage1Exists:
            self.check_V1(Voltage1)

        if self.Voltage2Exists:
            self.check_V2(Voltage2)

        self.blit_data(self.setup_artists)

    def ChangeParameterFrame(self, segment):

        if segment == 'Reverse':
            self.ReverseSelect['style'] = 'On.TButton'
            self.ForwardsSelect['style'] = 'Off.TButton'

            self.ReverseFrame.tkraise()

        elif segment == 'Forwards':
            self.ReverseSelect['style'] = 'Off.TButton'
            self.ForwardsSelect['style'] = 'On.TButton'

            self.ForwardsFrame.tkraise()

    def ApplyVoltage(self):
        global low_voltage, high_voltage

        low_voltage = float(self.LowVoltage.get())
        high_voltage = float(self.HighVoltage.get())

        if self.ForwardsSelectionExists:
            self.ForwardsSelection(None)

        if self.ReverseSelectionExists:

            self.ReverseSelection(None)

        self.blit_data(self.setup_artists)

    def SelectVoltagesFunction(self):
        global Voltage1, Voltage2

        Voltage1 = float(self.Voltage1.get())
        Voltage2 = float(self.Voltage2.get())

        if Voltage1 is not None:
            try:
                self.check_V1(Voltage1)

            except:
                print('\nCould not retrieve voltage 1')

        if Voltage2 is not None:
            try:
                self.check_V2(Voltage2)
            except:
                print('\nCould not retrieve voltage 2')

        if self.ForwardsSelectionExists:
            if self.ReverseSelectionExists:
                self.blit_data(self.setup_artists)

    def check_V1(self, voltage):
        global Voltage1

        Voltage1 = ClosestVoltageEstimation(voltage,self.potentials)

        passkey = True
        if self.ForwardsSelectionExists:
            if Voltage1 not in self.forward_potentials:
                passkey = False

        if self.ReverseSelectionExists:
            if Voltage1 not in self.reverse_potentials:
                passkey = False

        if passkey:
            self.Voltage1Exists = True

            if self.Voltage1Warning:
                self.Voltage1Error.grid_forget()

            return True

        else:
            self.Voltage1Error['text'] = 'Voltage 1 not in\nboth segments'
            self.Voltage1Error.grid(row=4,column=0)
            self.Voltage1Warning = True

            if self.Voltage1Exists:
                self.Voltage1Exists = False

            return False

    def check_V2(self, voltage):
        global Voltage2

        Voltage2 = ClosestVoltageEstimation(voltage,self.potentials)

        passkey = True
        if self.ForwardsSelectionExists:
            if Voltage2 not in self.forward_potentials:
                passkey = False

        if self.ReverseSelectionExists:
            if Voltage2 not in self.reverse_potentials:
                passkey = False

        if passkey:
            self.Voltage2Exists = True

            if self.Voltage2Warning:
                self.Voltage2Error.grid_forget()

            return True

        else:
            self.Voltage2Error['text'] = 'Voltage 2 not in\nboth segments'
            self.Voltage2Error.grid(row=4,column=1)
            self.Voltage2Warning = True

            if self.Voltage2Exists:
                self.Voltage2Exists = False

            return False

    def blit_data(self, artists):

        axes = {a.axes for a in artists}
        for a in axes:
            if a in self.bg_cache:
                a.figure.canvas.restore_region(self.bg_cache[a])

        self._drawn_artists = artists
        self._drawn_artists = sorted(self._drawn_artists,
                                     key=lambda x: x.get_zorder())
        for a in self._drawn_artists:
            a.set_animated(True)

        updated_ax = []
        for a in self._drawn_artists:
            # If we haven't cached the background for this axes object, do
            # so now. This might not always be reliable, but it's an attempt
            # to automate the process.
            if a.axes not in self.bg_cache:
                self.bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            a.axes.draw_artist(a)
            updated_ax.append(a.axes)

        # After rendering all the needed artists, blit each axes individually.
        for ax in set(updated_ax):
            ax.figure.canvas.blit(ax.bbox)

    def AlterExportPath(self):
        global FileHandle

        ### Path for directory in which the    ###
        ### exported .txt file will be placed  ###
        ExportPath = filedialog.askdirectory(parent = self.parent)
        ExportPath = ''.join(ExportPath + '/')

        FileHandle = FileHandle.split('/')
        FileHandle = FileHandle[-1]
        FileHandle = ''.join(ExportPath + FileHandle)

        DataFolder = FileHandle.split('/')
        DataFolder = ''.join(DataFolder[-3]+'/'+DataFolder[-2])

        self.AlterExportPathButton['style'] = 'On.TButton'
        self.AlterExportPathButton['text'] = DataFolder

    def Checkpoint(self):

        passkey = True

        if self.ForwardsSelectionExists:
            if self.ForwardsWarning:
                self.NoForwardsSelection.grid_forget()
        else:
            self.NoForwardsSelection.grid(row=4,column=0)
            self.ForwardsWarning = True

            passkey = False

        if self.ReverseSelectionExists:
            if self.ReverseWarning:
                self.NoReverseSelection.grid_forget()
        else:
            self.NoReverseSelection.grid(row=4,column=1)
            self.ReverseWarning = True

            passkey = False


        if self.Voltage1Exists:
            if self.Voltage1Warning:
                self.Voltage1Error.grid_forget()
        else:
            self.Voltage1Error.grid(row=4,column=0)
            self.Voltage1Warning = True

            passkey = False

        if self.Voltage2Exists:
            if self.Voltage2Warning:
                self.Voltage2Error.grid_forget()
        else:
            self.Voltage2Error.grid(row=4,column=1)
            self.Voltage2Warning = True

            passkey = False

        if passkey:
            self.Initialize()

    def Initialize(self):
        global CV_min, CV_max, data_min, post_analysis, data_max,V1_min,V1_max,V2_min,V2_max,PlotValues, PlotFrames, PlotContainer, text_file_export

        #-- CV Plot Y Limit Adjustment Values --#
        CV_min = float(self.CV_min.get())            # raw data y limit adjustment variables
        CV_max = float(self.CV_max.get())

        #-- Data (PHE/AUC) Plot Y Limit Adjustment Values --#
        data_min = float(self.data_min.get())               # raw data y limit adjustment variables
        data_max = float(self.data_max.get())

        #-- Voltage 1 Plot Y Limit Adjustment Values --#
        V1_min = float(self.V1_min.get())               # raw data y limit adjustment variables
        V1_max = float(self.V1_max.get())

        #-- Voltage 1 Plot Y Limit Adjustment Values --#
        V2_min = float(self.V2_min.get())               # raw data y limit adjustment variables
        V2_max = float(self.V2_max.get())

        ######################################################
        ### Create a figure and artists for each electrode ###
        ######################################################
        for num in range(electrode_count):
            electrode = electrode_list[num]
            figure = initialize.MakeFigure(electrode)
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

        #--- If the user has indicated that text file export should be activated ---#
        if text_file_export:
            print('Initializing Text File Export')
            text_file_export = TextFileExport()




        ################################################
        ### Initialize the RealTimeManipulationFrame ###
        ################################################
        frame = RealTimeManipulationFrame(container, self)
        ShowFrames[RealTimeManipulationFrame] = frame
        frame.grid(row=0,column=0,sticky='news')

        ############################
        ### Post Analysis Module ###
        ############################
        #post_analysis = PostAnalysis(self.parent, self.controller)
        #ShowFrames[PostAnalysis] = post_analysis
        #post_analysis.grid(row=0, column=0, sticky = 'nsew')

        show_frame(RealTimeManipulationFrame)
        show_plot(PlotValues[0])




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
        global PlotValues, container, ShowFrames, FileLabel

        tk.Frame.__init__(self, parent)         # Initialize the frame


        #######################################
        #######################################
        ### Pack the widgets into the frame ###
        #######################################
        #######################################

        #--- Display the file number ---#
        FileTitle = tk.Label(self, text = 'File Number', font=MEDIUM_FONT,)
        FileTitle.grid(row=0,column=0,padx=5,pady=5,sticky='nw')
        FileLabel = ttk.Label(self, text = '1', font=LARGE_FONT, style='Fun.TButton')
        FileLabel.grid(row=1,column=0,padx=5,pady=5,sticky='nw')

        #---Buttons to switch between electrode frames---#
        frame_value = 0
        column_value = 0
        row_value = 2
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



    ########################################################
    ### Function to Reset and raise the user input frame ###
    ########################################################
    def Reset(self):
        global key, PoisonPill, AlreadyReset, analysis_complete, AlreadyInitiated

        key = 0
        PoisonPill = True
        AlreadyInitiated = False # reset the start variable
        AlreadyReset = True

        # Raise the initial user input frame
        show_frame(InputFrame)
        close_frame(RealTimeManipulationFrame, True)
        close_frame(SetupFrame, False)
        close_frame(SetupVisualizationFrame, False)

        #post_analysis._reset()

        ## Take resize weight away from the Visualization Canvas
        container.columnconfigure(1, weight=0)

        analysis_complete = False


    ##########################################################
    ### Function to raise frame to the front of the canvas ###
    ##########################################################
    def show_frame(self, cont):

        frame = ShowFrames[cont]            # Key: frame handle / Value: tk.Frame object
        frame.tkraise()                     # raise the frame objext

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
                anim.append(ElectrochemicalAnimation(fig, self.ax, electrode, resize_interval = resize_interval, fargs=None))
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

        # close all matplotlib figures
        plt.close('all')

        # destory the frames holding the figures
        for frame in PlotValues:
            frame.destroy()

        # destory the container holding those frames
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

class SetupVisualizationFrame(tk.Frame):
    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)

        #--- for resize ---#
        self.columnconfigure(0, weight = 2)

        Label = tk.Label(self, text='Setup' ,font=HUGE_FONT)
        Label.grid(row=0,column=0,pady=5,sticky='n')

        fig, ax = setup_dictionary['figure']

        #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
        canvas = FigureCanvasTkAgg(fig, self)                                         # and place the artists within the frame
        canvas.draw()                                                           # initial draw call to create the artists that will be blitted
        canvas.get_tk_widget().grid(row=1,columnspan=1,pady=6,ipady=5,sticky='news')          # does not affect size of figure within plot container


#########################################################
### Electrode Frame Class for data visualization      ###
### displayed next to the RealTimeManipulationFrame   ###
###                                                   ###
###    Embeds a canvas within the tkinter             ###
###    MainWindow containing figures that             ###
###    visualize the data for that electrode          ###
#########################################################

class VisualizationFrame(tk.Frame):
    def __init__(self, electrode, frame_count, parent, controller):
        global FrameFileLabel

        tk.Frame.__init__(self, parent)

        #--- for resize ---#
        self.columnconfigure(0, weight = 2)
        self.columnconfigure(1, weight = 1)
        self.rowconfigure(2, weight=2)

        ElectrodeLabel = tk.Label(self, text='%s' % electrode ,font=HUGE_FONT)
        ElectrodeLabel.grid(row=0,column=0,columnspan=2,pady=5,sticky='n')

        FrameFileLabel = tk.Label(self, text = '', font=MEDIUM_FONT)
        FrameFileLabel.grid(row=0,column=1,pady=3,sticky='ne')

        #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
        fig, ax = figures[frame_count]                                                # Call the figure and artists for the electrode
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




                                        #############################################################
                                        #############################################################
                                        ### Creation of Matplotlib Canvas, Figures, Axes, Artists ###
                                        ### and all other decorators (e.g. axis labels, titles)   ###
                                        #############################################################
                                        #############################################################


class InitializeFigureCanvas():
    def __init__(self):
        global text_file_export, voltage1_data, voltage2_data,electrode_count, anim, Frame, FrameReference, FileHandle, PlotContainer, data_list, peak_data_list, plot_list, EmptyPlots, file_list, figures, frame_list, Plot, PlotFrames, PlotValues

        ##############################################
        ### Generate global lists for data storage ###
        ##############################################

        electrode_count = int(electrode_count)

        #--- Animation list ---#
        anim = []

        #--- Figure lists ---#
        figures = []

        ############################################
        ### Create global lists for data storage ###
        ############################################
        data_list = [0]*electrode_count                             # Peak Height/AUC data (after smoothing and polynomial regression)
        avg_data_list = []                                          # Average Peak Height/AUC across all electrodes for each frequency
        std_data_list = []                                          # standard deviation between electrodes for each frequency
        for num in range(electrode_count):
            data_list[num] = {}
            for count in ['forward','reverse']:
                data_list[num][count] = [0]*numFiles                      # a data list for each eletrode

        #-- Data for peak potentials --#
        peak_data_list = [0]*electrode_count
        for num in range(electrode_count):
            peak_data_list[num] = {}
            for count in ['forward','reverse']:
                peak_data_list[num][count] = [0]*numFiles                      # a data list for each eletrode

        #-- Data Stored for Voltage 1 --#
        voltage1_data = [0]*electrode_count
        for num in range(electrode_count):
            voltage1_data[num] = {}
            for count in ['forward','reverse']:
                voltage1_data[num][count] = [0]*numFiles

        #-- Data Stored for Voltage 2 --#
        voltage2_data = [0]*electrode_count
        for num in range(electrode_count):
            voltage2_data[num] = {}
            for count in ['forward','reverse']:
                voltage2_data[num][count] = [0]*numFiles

        #--- Lists of Frames and Artists ---#
        plot_list = []
        frame_list = []

        #--- Misc Lists ---#
        file_list = []          # Used for len(file_list)

        ######################################
        ### Create an initial setup figure ###
        ######################################
        self.create_setup(electrode_list[0])


    ############################################
    ### Create the figure and artist objects ###
    ############################################
    def MakeFigure(self, electrode):
        global list_val, EmptyPlots, plot_list, ScanRate, frame_list, numFiles

        ########################
        ### Setup the Figure ###
        ########################
        fig, ax = plt.subplots(nrows=3,ncols=2,squeeze=False,figsize=(9,6))    ## figsize=(width, height)
        plt.subplots_adjust(bottom=0.1,hspace=1.1,wspace=0.3)         ### adjust the spacing between subplots

        #---Set the electrode index value---#
        list_val = _get_listval(electrode)

        #######################
        ### Set axis labels ###
        #######################
        #-- Cyclic Voltammogram Plot --#
        ax[0,0].set_ylabel('Current/µA',fontweight='bold')
        ax[0,0].set_xlabel('Potential/V')

        #-- PHE/AUC Plot --#
        ax[0,1].set_xlabel('File Number')
        if SelectedOptions == 'Peak Height Extraction':
            ax[0,1].set_ylabel('Peak Height/µA',fontweight='bold')
        elif SelectedOptions == 'Area Under the Curve':
            ax[0,1].set_ylabel('AUC/µC',fontweight='bold')

        #-- User-specified voltage Plots --#
        ax[1,0].set_xlabel('File Number')   # Forwards Peak Potential Plot
        ax[1,1].set_xlabel('File Number')   # Reverse Peak Potential Plot
        ax[1,0].set_ylabel('Current/µA')    # Forwards Peak Potential Plot
        ax[1,1].set_ylabel('Current/µA')    # Reverse Peak Potential Plot

        #-- User-specified voltage Plots --#
        ax[2,0].set_xlabel('File Number')    # Voltage 1
        ax[2,1].set_xlabel('File Number')    # Voltage 2
        ax[2,0].set_ylabel('Current/µA')  # Voltage 1
        ax[2,1].set_ylabel('Current/µA')  # Voltage 2


        #############################################
        ### Set X Limits of all plots besides the ###
        ### Cyclic Voltammogram plot (ax[0,0])    ###
        #############################################

        #--- if the resize interval is larger than the number of files, ---#
        #--- make the x lim the number of files (& vice versa)          ---#
        if resize_interval > numFiles:
            xlim_factor = numFiles
        elif resize_interval <= numFiles:
            xlim_factor = resize_interval

        #-- Raw PHE/AUC Data --#
        ax[0,1].set_xlim(-0.05,xlim_factor+0.1)

        #-- Peak Potential Forwards and Reverse --#
        ax[1,0].set_xlim(-0.05,xlim_factor+0.1)
        ax[1,1].set_xlim(-0.05,xlim_factor+0.1)

        #-- Voltage Probes 1 and 2 --#
        ax[2,0].set_xlim(-0.05,xlim_factor+0.1)
        ax[2,1].set_xlim(-0.05,xlim_factor+0.1)

        #################################################################################
        #################################################################################
        ###       Analyze the first file and create the Y limits of the subplots      ###
        ###               depending on the data range of the first file               ###
        #################################################################################
        self.InitializeSubplots(ax, electrode)

        #################################################################################
        #################################################################################


        #---Set Subplot Title---#
        ax[0,0].set_title('Cyclic Voltammogram',fontweight='bold')

        if SelectedOptions == 'Peak Height Extraction':
            ax[0,1].set_title('Peak Height Extraction',fontweight='bold')
        elif SelectedOptions == 'Area Under the Curve':
            ax[0,1].set_title('Area Under the Curve',fontweight='bold')

        ax[1,0].set_title('Forward Peak Potential')
        ax[1,1].set_title('Reverse Peak Potential')

        ax[2,0].set_title(''.join(str(Voltage1)+'V'))
        ax[2,1].set_title(''.join(str(Voltage2)+'V'))

        #---Initiate the subplots---#
        # this assigns a Line2D artist object to the artist object (Axes)

        #-- Cyclic Voltammogram Artists --#
        raw_forward, = ax[0,0].plot([],[],'o',color='0.75',MarkerSize=0.5,zorder=1)   # Forwards
        adjusted_forward, = ax[0,0].plot([],[],'k-',zorder=3)
        forward_baseline, = ax[0,0].plot([],[],'k-',MarkerSize=0.5,zorder=4)
        forward_baseline.set_alpha(1)

        #-- shading for Forwards AUC --#
        forward_verts = [(0,0),*zip([],[]),(0,0)]
        forward_poly = Polygon(forward_verts, alpha = 0.5,color='g',zorder=2)
        ax[0,0].add_patch(forward_poly)

        raw_reverse, = ax[0,0].plot([],[],'o',color='0.75',MarkerSize=0.5,zorder=1)     # Reverse
        adjusted_reverse, = ax[0,0].plot([],[],'k-',zorder=3)
        reverse_baseline, = ax[0,0].plot([],[],'k-',MarkerSize=0.5,zorder=4)
        reverse_baseline.set_alpha(1)

        #-- shading for Reverse AUC --#
        reverse_verts = [(0,0),*zip([],[]),(0,0)]
        reverse_poly = Polygon(reverse_verts, alpha = 0.5, color='r',zorder=2)
        ax[0,0].add_patch(reverse_poly)


        #-- PHE/AUC Artists --#
        forward_data, = ax[0,1].plot([],[],'go',Markersize=2,label='Forwards')
        reverse_data, = ax[0,1].plot([],[],'ro',Markersize=2,label='Reverse')

        #-- Peak Potential Forwards --#
        peak_forward, = ax[1,0].plot([],[],'go',MarkerSize=1)

        #-- Peak Potential Reverse --#
        peak_reverse, = ax[1,1].plot([],[],'ro',MarkerSize=1)

        #-- Voltage 1 Artists --#
        volt1_forward, = ax[2,0].plot([],[],'go',MarkerSize=1)
        volt1_reverse, = ax[2,0].plot([],[],'ro',MarkerSize=1)

        #-- Voltage 2 Artists --#
        volt2_forward, = ax[2,1].plot([],[],'go',MarkerSize=1)
        volt2_reverse, = ax[2,1].plot([],[],'ro',MarkerSize=1)

        #ax[0,1].legend(bbox_to_anchor=(0., 0, 1.05, 105), bbox_transform=gcf().transFigure)

        #####################################################
        ### Create a list of the primitive artists        ###
        ### (Line2D objects) that will be returned        ###
        ### to ElectrochemicalAnimation to be visualized  ###
        #####################################################

        # this is the list that will be returned as _drawn_artists to the Funcanimation class
        plots = [raw_forward,adjusted_forward,forward_baseline,forward_poly,raw_reverse,adjusted_reverse,reverse_baseline,reverse_poly,forward_data,reverse_data,volt1_forward,volt1_reverse,volt2_forward,volt2_reverse,peak_forward,peak_reverse]

        #--- And append that list to keep a global reference ---#
        electrode_frame = 'Electrode %s' % str(electrode)
        if electrode_frame not in frame_list:
            frame_list.append(electrode_frame)

        #--- Create empty plots to return to animate for initializing---#
        EmptyPlots = [raw_forward]

        plot_list.append(plots)        # 'plot_list' is a list of lists containing 'plots' for each electrode

        #-- Return both the figure and the axes to be stored as global variables --#
        return fig, ax

    #####################################################################################
    ### Initalize Y Limits of each figure depending on the y values of the first file ###
    #####################################################################################
    def InitializeSubplots(self,ax,electrode):

        self.list_val = _get_listval(electrode)

        try:
            print('Initialize: Electrode %s' % str(electrode))
            filename, filename2 = _retrieve_file(1,electrode)

            myfile = mypath + filename               ### path of your file
            myfile2 = mypath + filename2               ### path of your file

            try:
                mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

            except:

                try:
                    mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                    if mydata_bytes > 1000:
                        myfile = myfile2
                except:
                    mydata_bytes = 1

            if mydata_bytes > 1000:
                print('Found File %s' % myfile)
                self.RunInitialization(myfile,ax, electrode)

            else:
                return False


        except:
            print('could not find file for electrode %d' % electrode)
            #--- If search time has not met the search limit keep searching ---#
            root.after(1000, self.InitializeSubplots, ax, electrode)


    def RunInitialization(self, myfile, ax, electrode):
        global sign_dictionary

        try:
            ###################################
            ### Retrieve data from the File ###
            ###################################
            potentials, currents, segment_dictionary, forward, reverse = ReadData(myfile, electrode)


            ###################################################
            ### Create local data dictionaries that will be ###
            ### passed on to _animate to be visualized      ###
            ###################################################
            segment_data = {}
            adjusted_segment_data = {}
            baseline_data = {}
            vertex_data = {}

            ##########################################
            ### Set the x axes of the voltammogram ###
            ##########################################
            MIN_POTENTIAL = min(potentials)
            MAX_POTENTIAL = max(potentials)
            MIN_CURRENT = min(currents)
            MAX_CURRENT = max(currents)


            ########################################
            ### Iterate through the two segments ###
            ########################################
            V1_data = []
            V2_data = []
            local_data_list = []
            peak_potential_list = {}

            ########################################
            ### Iterate through the two segments ###
            ########################################

            local_dictionary = {}
            for count in ['forward','reverse']:
                try:
                    #-- 1st, forward --#
                    if count == 'forward':
                        self.index_list, self.potential_dict, self.data_dict = segment_dictionary[forward_segment]
                        self.segment_potentials = [potentials[index] for index in self.index_list]
                        self.segment_currents = [currents[index] for index in self.index_list]
                        self.adjusted_potentials = [potential for potential in self.segment_potentials if low_voltage <= potential <= high_voltage]
                        self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_potentials]
                        self.adjusted_currents = [currents[index] for index in self.adjusted_index_list]
                        f_initial = max(self.adjusted_currents)

                    #-- 2nd, reverse --#
                    elif count == 'reverse':
                        self.index_list, self.potential_dict, self.data_dict = segment_dictionary[reverse_segment]
                        self.segment_potentials = [potentials[index] for index in self.index_list]
                        self.segment_currents = [currents[index] for index in self.index_list]
                        self.adjusted_potentials = [potential for potential in self.segment_potentials if low_voltage <= potential <= high_voltage]
                        self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_potentials]
                        self.adjusted_currents = [currents[index] for index in self.adjusted_index_list]
                        r_initial = max(self.adjusted_currents)
                except:
                    print('RunInitialization: Error in Initial Data Adjustment')

            #- if the reverse scan is negative compared
            #- to the forward scan, it is anodic and the
            #- forward is cathodic
            self.sign_dictionary = {}
            if r_initial < f_initial:
                self.sign_dictionary['reverse'] = 'anodic'
                self.sign_dictionary['forward'] = 'cathodic'
            #- and vice versa
            else:
                self.sign_dictionary['forward'] = 'anodic'
                self.sign_dictionary['reverse'] = 'cathodic'

            #- create a global variable -#
            sign_dictionary = self.sign_dictionary

            #-- initial analysis --#
            for count in ['forward','reverse']:

                #-- 1st, forward --#
                if count == 'forward':
                    segment = forward_segment

                #-- 2nd, reverse --#
                elif count == 'reverse':
                    segment = reverse_segment

                ###############################################
                ### Extact the data for the current segment ###
                ###############################################
                try:
                    self.index_list, self.potential_dict, self.data_dict = segment_dictionary[segment]
                    self.segment_potentials = [potentials[index] for index in self.index_list]
                    self.segment_currents = [currents[index] for index in self.index_list]

                    #-- save it to the artist data list --#
                    segment_data[count] = [self.segment_potentials,self.segment_currents]

                    #-- calculate the step potential --#
                    step_potential = self.segment_potentials[1] - self.segment_potentials[0]
                    if step_potential < 0:
                        step = 'Negative'
                    elif step_potential > 0:
                        step = 'Positive'
                except:
                    print('RunInitialization: Error in %s Segment Extraction' % count)

                #############################################
                ### Adjust the potentials and currents to ###
                ### the voltage range chosen by the user  ###
                #############################################
                try:
                    self.adjusted_potentials = [potential for potential in self.segment_potentials if low_voltage <= potential <= high_voltage]
                    self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_potentials]
                    self.adjusted_currents = [currents[index] for index in self.adjusted_index_list]
                    adjusted_segment_data[count] = [self.adjusted_potentials,self.adjusted_currents]
                except:
                    print('RunInitialization: Error in Value Adjustment')

                ################################################
                ### Extract Data for User-Specified Voltages ###
                ################################################

                #-- if the user chose a value that does not --#
                #-- exist, find the next closest value      --#
                V1 = ClosestVoltageEstimation(Voltage1,self.segment_potentials)
                V2 = ClosestVoltageEstimation(Voltage2,self.segment_potentials)

                #-- find the associated index values --#
                Voltage1_Index = self.potential_dict[V1]
                Voltage2_Index = self.potential_dict[V2]

                #-- and the associated currents --#
                Voltage1_Current = currents[Voltage1_Index]
                Voltage2_Current = currents[Voltage2_Index]

                V1_data.append(Voltage1_Current)
                V2_data.append(Voltage2_Current)

                ###########################################################
                ### Extract the Peak Potential and its associated index ###
                ###########################################################
                try:
                    peak_potential, self.peak_index = self.extract_peak(count)
                    vertex1, vertex2 = self.extract_vertex_currents(count)

                    #-- correlate the peak potential to its segment
                    peak_potential_list[count] = peak_potential

                except:
                    print('RunInitialization: Error in Peak Potential Extraction')

                ###############################################################
                ### Extract the Vertices to be used for the lienar baseline ###
                ###############################################################
                try:
                    #-- Split the potentials before and after the peak potential --#
                    vertex1_potentials = self.adjusted_potentials[:self.peak_index]
                    vertex2_potentials = self.adjusted_potentials[self.peak_index:]

                    #-- Find the potentials that correspond to the vertex currents
                    vertex1_potential = self.data_dict[vertex1]        # extract the associated potential
                    vertex2_potential = self.data_dict[vertex2]

                    #-- if multiple potential values have a corresponding current
                    #-- that equals the vertex current, find the closest to the peak
                    vertex1_potential, index1 = self.extract_vertex_potentials(vertex1_potential,vertex1_potentials)
                    vertex2_potential, index2 = self.extract_vertex_potentials(vertex2_potential,vertex2_potentials)

                except:
                    print('\nRunInitialization: Error in Vertex Potential Extraction\n')

                ####################################################################
                ### Adjust the potentials and currents to the extracted vertices ###
                ### and create a linear baseline using these vertices            ###
                ####################################################################
                try:
                    #-- adjust the segment potentials and
                    #-- currents to the bounds of the vertices
                    baseline_currents = self.adjusted_currents[index1:index2]
                    baseline_potentials = self.adjusted_potentials[index1:index2]

                    #-- create an index dictionary for the new baseline values
                    baseline_current_dict = {}
                    for index in range(len(self.adjusted_currents)):
                        baseline_current_dict[self.adjusted_currents[index]] = index
                except:
                    print('RunInitialization: Error in Baseline Value Adjustment')

                #-- if it is a surface bound species --#
                if mass_transport == 'surface':

                    try:
                        proto_baseline = np.linspace(vertex1,vertex2,len(baseline_currents))
                        linear_baseline, baseline_indeces,baseline_potentials,baseline_currents,peak_index_dict = self.extract_baseline(count,proto_baseline,baseline_potentials,baseline_currents)

                        #-- make a new dictionary with key: value
                        #-- as baseline_current: baseline_potential
                        baseline_dict = {}
                        for index in range(len(baseline_currents)):
                            baseline_dict.setdefault(baseline_currents[index],[]).append(baseline_potentials[index])

                        # recalculate the peak potential of the baseline currents #
                        if self.sign_dictionary[count] == 'cathodic':
                            peak_potential = baseline_dict[max(baseline_currents)][0]  # extract the potential of the peak current
                        if self.sign_dictionary[count] == 'anodic':
                            peak_potential = baseline_dict[min(baseline_currents)][0]  # extract the potential of the peak current

                        baseline_peak_index = peak_index_dict[peak_potential]

                        vertices = [(vertex1_potential,vertex1), *zip(baseline_potentials,baseline_currents), ((vertex2_potential,vertex2))]
                        vertex_data[count] = vertices

                    except:
                        print('RunInitialization: Error in Surface Bound Baseline Creation')

                #-- if it is a solution phase species --#
                elif mass_transport == 'solution':

                    try:
                        try:
                            if self.sign_dictionary[count] == 'cathodic':
                                vertex_current = min(vertex1,vertex2)
                            elif self.sign_dictionary[count] == 'anodic':
                                vertex_current = max(vertex1,vertex2)

                            # retrieve the corresponding potential of the vertex current
                            vertex_potential = self.data_dict[vertex_current][0]

                            # retrieve the corresponding index of the vertex current/potential
                            vertex_index = baseline_current_dict[vertex_current]
                        except:
                            print('RunInitialization: Error in Solution Phase Vertex Calculations')

                        ############################################################
                        ### Calculate the Range to be used for baseline creation ###
                        ############################################################
                        try:
                            # set the `range` to 5% of the total points #
                            value_range = round(0.05*len(baseline_currents))

                            # upper limit is the max index value before the index
                            # is out of the bounds of the adjusted currents
                            upper_limit = len(baseline_currents) - vertex_index

                            if vertex_index < value_range:
                                print('Lower Vertex is under range')
                                lower_range = baseline_currents[:vertex_index]
                                lower_range_potentials = baseline_potentials[:vertex_index]
                            else:
                                lower_index = vertex_index - value_range
                                lower_range = baseline_currents[lower_index:vertex_index]
                                lower_range_potentials = baseline_potentials[lower_index:vertex_index]
                            if vertex_index < upper_limit:
                                print('Upper vertex below upper limit')
                                upper_index = vertex_index + value_range
                                upper_range = baseline_currents[vertex_index:upper_index]
                                upper_range_potentials = baseline_potentials[vertex_index:upper_index]

                            else:
                                upper_range = baseline_currents[vertex_index:]
                                upper_range_potentials = baseline_potentials[vertex_index:]
                        except:
                            print('\nRunInitialization: Error in Solution Phase Range Calculation')

                        try:
                            # create a list with the entire range of current values #
                            slope_currents = []
                            for value in lower_range:
                                slope_currents.append(value)
                            for value in upper_range:
                                slope_currents.append(value)

                            # create a list of slope potentials
                            slope_potentials = []
                            for value in lower_range_potentials:
                                slope_potentials.append(value)
                            for value in upper_range_potentials:
                                slope_potentials.append(value)


                            # calculate the slope of the data points surrounding the vertex #
                            # (a ±5% range)
                            slope = (slope_currents[-1]-slope_currents[0])/(slope_potentials[-1]-slope_potentials[0])

                            # using this slope, create a linear baseline using the entire
                            # range of potentials chosen by the user
                            proto_baseline = []
                            for index in range(len(baseline_currents)):
                                proto_baseline.append(slope*(baseline_potentials[index]-vertex_potential)+vertex_current)
                        except:
                            print('RunInitialization: Error in Solution Phase Slope Calculation')

                        linear_baseline,baseline_indeces,baseline_potentials,baseline_currents,peak_index_dict = self.extract_baseline(count,proto_baseline,baseline_potentials,baseline_currents)

                        ####################################
                        ### Calculate the Peak Potential ###
                        ####################################

                        try:

                            #-- make a new dictionary with key: value
                            #-- as baseline_current: baseline_potential
                            baseline_dict = {}
                            for index in range(len(baseline_currents)):
                                baseline_dict.setdefault(baseline_currents[index],[]).append(baseline_potentials[index])

                            # recalculate the peak potential of the baseline currents #
                            if self.sign_dictionary[count] == 'cathodic':
                                peak_potential = baseline_dict[max(baseline_currents)][0]  # extract the potential of the peak current

                            if self.sign_dictionary[count] == 'anodic':
                                peak_potential = baseline_dict[min(baseline_currents)[0]]  # extract the potential of the peak current

                            baseline_peak_index = peak_index_dict[peak_potential]
                            vertices = [(baseline_potentials[0],baseline_currents[0]), *zip(baseline_potentials,baseline_currents), ((baseline_potentials[-1],baseline_currents[-1]))]
                            vertex_data[count] = vertices
                        except:
                            print('RunInitialization: Error in Solution Phase Peak Potential Calculation')

                    except:
                        print('RunInitialization: Error in Solution Phase Baseline Creation')

                baseline_data[count] = [baseline_potentials, linear_baseline]

                ################################################################
                ### If the user selected Peak Height Extraction, analyze PHE ###
                ################################################################

                if SelectedOptions == 'Peak Height Extraction':

                    try:
                        Peak_Height = abs(baseline_currents[baseline_peak_index] - linear_baseline[baseline_peak_index])
                        data = Peak_Height
                        local_data_list.append(data)
                    except:
                        print('RunInitialization: Error in PHE')

                ########################################################
                ### If the user selected AUC extraction, analyze AUC ###
                ########################################################

                elif SelectedOptions == 'Area Under the Curve':
                    try:
                        AUC_index = 1
                        AUC = 0

                        while AUC_index < len(baseline_currents) - 1:
                            AUC_height = (baseline_currents[AUC_index] - linear_baseline[AUC_index])
                            AUC_width = step_potential
                            AUC += (AUC_height * AUC_width)
                            AUC_index += 1

                        AUC = AUC/scan_rate
                        data = AUC
                        local_data_list.append(data)
                    except:
                        print('RunInitialization: Error in AUC')

            ###########################
            ## Set the X and Y Axes ###
            ###########################
            try:
                ### Cyclic Voltammogram Plot ###

                #-- set the x limits as the min and max potentials of the voltammogram
                ax[0,0].set_xlim(MIN_POTENTIAL,MAX_POTENTIAL)
                #-- set the ylimits as the max/min of the currents plus/minus the
                #-- their absolute value times InputFrame global object 'CV_max/min'
                #-- value
                ax[0,0].set_ylim(MIN_CURRENT-(abs(MIN_CURRENT)*CV_min),MAX_CURRENT+(CV_max*abs(MAX_CURRENT)))   # voltammogram

                ### PHE/AUC Plot Y Limits ###

                #- find the max/min of the Peah Height of the
                #- first file for the forward and reverse scans
                #- local_data_list = [Peak_Height_forward,Peak_Height_reverse]
                local_min = min(local_data_list)
                local_max = max(local_data_list)

                #-- set the ylimits as the max/min plus/minus the aboslute value
                #-- of their value times the InputFrame global objext 'data_max/min'
                ax[0,1].set_ylim(local_min-abs(data_min*local_min),local_max+abs(data_max*local_max))  # raw peak height

                ### Peak Potential Plots ###
                forward_peak = peak_potential_list['forward']
                reverse_peak = peak_potential_list['reverse']
                ax[1,0].set_ylim(forward_peak-abs(data_min*forward_peak),forward_peak+abs(data_max*forward_peak))
                ax[1,1].set_ylim(reverse_peak-abs(data_min*reverse_peak),reverse_peak+abs(data_max*reverse_peak))


                ### Probe VoltagePlots ###

                #-- Set Voltage 1 Plot Y Limits as the max/min of the forward and
                #-- reverse currents for voltage 1 plus/minus the aboslute value
                #-- of their value times the InputFrame global objext 'data_max/min'
                ax[2,0].set_ylim(min(V1_data)-(V1_min*abs(min(V1_data))),max(V1_data)+abs((V1_max*max(V1_data))))

                #-- SetVoltage 2 Plot Y Limits as the max/min of the forward and
                #-- reverse currents for voltage 2 plus/minus the aboslute value
                #-- of their value times the InputFrame global objext 'data_max/min'
                ax[2,1].set_ylim(min(V2_data)-(V2_min*abs(min(V2_data))),max(V2_data)+abs((V2_max*max(V2_data))))
            except:
                print('RunInitialization: Error in Axes Creation')

            return True

        except:
            print('\n\nError in RunInitialization\n\n')

    def create_setup(self,electrode):
        global setup_dictionary

        #########################
        ### Retrieve the data ###
        #########################
        filename, filename2 = _retrieve_file(1,electrode)

        myfile = mypath + filename               ### path of your file
        myfile2 = mypath + filename2               ### path of your file
        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

        except:

            try:
                mydata_bytes = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes
                if mydata_bytes > 1000:
                    myfile = myfile2
            except:
                mydata_bytes = 1

        if mydata_bytes > 1000:
            potentials, currents, segment_dictionary, forward, reverse = ReadData(myfile, 1)

            fig = plt.figure(figsize=(9,6))
            ax = fig.add_subplot(111)

            ax.set_xlim(min(potentials),max(potentials))
            ax.set_ylim(min(currents),max(currents))

            setup_figure = (fig, ax)

            value = 0
            self.potential_dict = {}
            for potential in potentials:
                self.potential_dict[potential] = value
                value += 1

            setup_dictionary = {}
            setup_dictionary['potentials'] = potentials
            setup_dictionary['self.potential_dict'] = self.potential_dict
            setup_dictionary['currents'] = currents
            setup_dictionary['segments'] = segment_dictionary
            setup_dictionary['figure'] = setup_figure
            setup_dictionary['forward'] = forward
            setup_dictionary['reverse'] = reverse



    def extract_peak(self, count):
        ### extract the potential of the peak current and then extract
        ### the associated index adjusted to the user-chosen range
        try:
            #-- if the peak is positive (cathodic) --#
            if self.sign_dictionary[count] == 'cathodic':
                peak_potential = self.data_dict[max(self.adjusted_currents)][0]
                peak_index = self.potential_dict[peak_potential] - self.adjusted_index_list[0]

            #-- if the peak is negative (anodic) --#
            elif self.sign_dictionary[count] == 'anodic':
                peak_potential = self.data_dict[min(self.adjusted_currents)][0]
                peak_index = self.potential_dict[peak_potential]  - self.adjusted_index_list[0]
            return peak_potential, peak_index

        except:
            print('RunInitialization: Error in extract_peak()')


    def extract_vertex_currents(self, count):

        try:

            if sign_dictionary[count] == 'cathodic':
                #-- Split the voltammogram before and after the peak  --#
                #-- potential and extract the maxima for each besides --#
                try:
                    vertex1 = min(self.adjusted_currents[:self.peak_index])
                    vertex2 = min(self.adjusted_currents[self.peak_index:])
                except:
                    fit_half = int(len(self.adjusted_currents)/2)
                    vertex1 = min(self.adjusted_currents[:fit_half])
                    vertex2 = min(self.adjusted_currents[fit_half:])

            if sign_dictionary[count] == 'anodic':
                #-- Split the voltammogram before and after the peak  --#
                #-- potential and extract the maxima for each besides --#
                try:
                    vertex1 = max(self.adjusted_currents[:self.peak_index])
                    vertex2 = max(self.adjusted_currents[self.peak_index:])
                except:
                    fit_half = int(len(self.adjusted_currents)/2)
                    vertex1 = max(self.adjusted_currents[:fit_half])
                    vertex2 = max(self.adjusted_currents[fit_half:])

            return vertex1, vertex2

        except:
            print('RunInitialization: Error in extract_vertex_currents()')


    def extract_vertex_potentials(self, vertex_potential, vertex_potentials):

        ### if multiple potential values have a corresponding current     ###
        ### that equals the vertex current, find the closest to the peak  ###
        ### and extract the associated index value                        ###

        if len(vertex_potential) == 1:
            vertex_potential = vertex_potential[0]
            index = self.potential_dict[vertex_potential] - self.adjusted_index_list[0] # and the index value
        else:
            list = []
            for vertex in vertex_potential:
                if min(vertex_potentials) <= vertex <= max(vertex_potentials):
                    list.append(vertex)

            vertex_potential = list[0]
            index = self.potential_dict[vertex_potential] - self.adjusted_index_list[0] # and the index value

        return vertex_potential,index

    def extract_baseline(self,count,proto_baseline,baseline_potentials,baseline_currents):
        try:
            # remove any currents that are above/below the currents of the segment
            linear_baseline = []
            baseline_indeces = []
            for index in range(len(proto_baseline)):
                if self.sign_dictionary[count] == 'cathodic':
                    if proto_baseline[index] <= baseline_currents[index]:
                        linear_baseline.append(proto_baseline[index])
                        baseline_indeces.append(index)
                if self.sign_dictionary[count] == 'anodic':
                    if proto_baseline[index] >= baseline_currents[index]:
                        linear_baseline.append(proto_baseline[index])
                        baseline_indeces.append(index)

            if len(linear_baseline) == 0:
                linear_baseline = proto_baseline
                baseline_indeces = [value for value in range(len(proto_baseline))]

            baseline_potentials = [baseline_potentials[index] for index in baseline_indeces]
            baseline_currents = [baseline_currents[index] for index in baseline_indeces]

            peak_index_dict = {}
            index = 0
            for potential in baseline_potentials:
                peak_index_dict[potential] = index
                index += 1

            return linear_baseline, baseline_indeces, baseline_potentials, baseline_currents,peak_index_dict

        except:
            print('RunInitialization: Error in extract_baseline()')
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
        self.list_val = _get_listval(electrode)
        self.file = starting_file                                # Starting File
        self.index = 0                                           # File Index Value
        self.ax = ax                                             # Figure Axes object

        ### and file count for each electrode    ###
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
            self.resize_interval = 200

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

    def _redraw_figures(self):

        ############################################
        ### Resize raw and normalized data plots ###
        ############################################
        fig, ax = figures[self.num]

        ax[0,1].set_xlim(0,self.resize_limit+0.1)

        ax[1,1].set_xlim(0,self.resize_limit+0.1)
        ax[1,0].set_xlim(0,self.resize_limit+0.1)

        ax[2,0].set_xlim(0,self.resize_limit+0.1)
        ax[2,1].set_xlim(0,self.resize_limit+0.1)


        #####################################################
        ### Set up the new canvas with an idle draw event ###
        #####################################################
        self._post_draw(True)


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
        global file_list, analysis_complete

        if self.file not in self.file_list:
            self.file_list.append(self.file)

        filename, filename2 = _retrieve_file(self.file, self.electrode)

        myfile = mypath + filename                    ### path of raw data .csv file
        myfile2 = mypath + filename2                   ### path of raw data .csv file

        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes

        except:
            mydata_bytes = 1

        try:
            mydata_bytes2 = os.path.getsize(myfile2)    ### retrieves the size of the file in bytes

        except:
            mydata_bytes2 = 1

        #################################################################
        #### If the file meets the size requirement, analyze the data ###
        #################################################################
        if mydata_bytes > 1000:
            #self.FileLabel['text'] = 'Found %s' % filename
            q.put(lambda: self._run_analysis(myfile))

        elif mydata_bytes2 > 1000:
            q.put(lambda: self._run_analysis(myfile2))

        else:
            if not PoisonPill:
                root.after(100,self._step)

    def _check_queue():

        while True:
            try:
                task = q.get(block=False)
            except:
                break
            else:
                if not PoisonPill:
                    root.after(1,self.task)

        if not analysis_complete:
            if not PoisonPill:
                root.after(5, self._check_queue)

    def _run_analysis(self,myfile):
        global post_analysis

        #######################################################
        ### Perform the next iteration of the data analysis ###
        #######################################################
        try:
            framedata = self.generator(myfile)
            self._draw_next_frame(framedata)

        except StopIteration:
            return False

        ##########################################################################
        ### if the resize limit has been reached, resize and redraw the figure ###
        ##########################################################################
        if self.file == self.resize_limit:

            # Dont redraw if this is the already the last file #
            if self.resize_limit < numFiles:

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

        track.tracking(self.file)

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
            root.after(1, self._step)


    def _raw_generator(self, myfile):

        ###################################
        ### Retrieve data from the File ###
        ###################################
        potentials, currents, segment_dictionary, forward, reverse = ReadData(myfile, self.electrode)

        ################################################
        ### Create local data dictionaries for the   ###
        ### cyclic voltammogram that will be passed  ###
        ### on to _animate to be visualized          ###
        ################################################
        segment_data = {}            # Raw CV Potentials and Currents
        adjusted_segment_data = {}   # Adjusted CV Potentials and Currents
        baseline_data = {}           # Linear Baseline Currents
        vertex_data = {}             # Linear Baseline Vertices

        print('%sElectrode Generator %s ' % (self.spacer,self.electrode))

        #-- initial analysis --#
        for count in ['forward','reverse']:

            #-- 1st, forward --#
            if count == 'forward':
                print('\n%sFORWARD' % self.spacer)
                segment = forward_segment
                color = 'g'

            #-- 2nd, reverse --#
            elif count == 'reverse':
                print('\n%sREVERSE' % self.spacer)
                segment = reverse_segment
                color = 'r'

            ###############################################
            ### Extact the data for the current segment ###
            ###############################################
            self.index_list, self.potential_dict, self.data_dict = segment_dictionary[segment]
            self.segment_potentials = [potentials[index] for index in self.index_list]
            self.segment_currents = [currents[index] for index in self.index_list]

            #-- save it to the artist data list --#
            segment_data[count] = [self.segment_potentials,self.segment_currents]

            #-- calculate the step potential --#
            step_potential = self.segment_potentials[1] - self.segment_potentials[0]
            if step_potential < 0:
                step = 'Negative'
            elif step_potential > 0:
                step = 'Positive'

            #############################################
            ### Adjust the potentials and currents to ###
            ### the voltage range chosen by the user  ###
            #############################################
            self.adjusted_potentials = [potential for potential in self.segment_potentials if low_voltage <= potential <= high_voltage]
            self.adjusted_index_list = [self.potential_dict[potential] for potential in self.adjusted_potentials]
            self.adjusted_currents = [currents[index] for index in self.adjusted_index_list]
            adjusted_segment_data[count] = [self.adjusted_potentials,self.adjusted_currents]

            ################################################
            ### Extract Data for User-Specified Voltages ###
            ################################################

            #-- if the user chose a value that does not --#
            #-- exist, find the next closest value      --#
            V1 = ClosestVoltageEstimation(Voltage1,self.segment_potentials)
            V2 = ClosestVoltageEstimation(Voltage2,self.segment_potentials)

            #-- find the associated index values --#
            Voltage1_Index = self.potential_dict[V1]
            Voltage2_Index = self.potential_dict[V2]

            #-- and the associated currents --#
            Voltage1_Current = currents[Voltage1_Index]
            Voltage2_Current = currents[Voltage2_Index]

            #-- save the data --#
            voltage1_data[self.num][count][self.index] = Voltage1_Current
            voltage2_data[self.num][count][self.index] = Voltage2_Current

            #-- retrieve the orientation of the segment --#
            self.sign_dictionary = sign_dictionary

            ###########################################################
            ### Extract the Peak Potential and its associated index ###
            ###########################################################
            try:
                peak_potential, self.peak_index = self.extract_peak(count)
                vertex1, vertex2 = self.extract_vertex_currents(count)

            except:
                print('\n%s_raw_generator: Error in Peak Potential Extraction\n' % self.spacer)

            ###############################################################
            ### Extract the Vertices to be used for the lienar baseline ###
            ###############################################################
            try:
                #-- Split the potentials before and after the peak potential --#
                vertex1_potentials = self.adjusted_potentials[:self.peak_index]
                vertex2_potentials = self.adjusted_potentials[self.peak_index:]

                #-- Find the potentials that correspond to the vertex currents
                vertex1_potential = self.data_dict[vertex1]        # extract the associated potential
                vertex2_potential = self.data_dict[vertex2]

                #-- if multiple potential values have a corresponding current
                #-- that equals the vertex current, find the closest to the peak
                vertex1_potential, index1 = self.extract_vertex_potentials(vertex1_potential,vertex1_potentials,peak_potential)
                vertex2_potential, index2 = self.extract_vertex_potentials(vertex2_potential,vertex2_potentials,peak_potential)

            except:
                print('\n%s_raw_generator: Error in Vertex Potential Extraction\n' % self.spacer)

            ###########################################################################
            ### Adjust the potentials and currents using the extracted vertices     ###
            ### as the boundaries and create a linear baseline using these vertices ###
            ###########################################################################
            try:
                #-- adjust the segment potentials and
                #-- currents to the bounds of the vertices
                baseline_currents = self.adjusted_currents[index1:index2]
                baseline_potentials = self.adjusted_potentials[index1:index2]
                baseline_indeces = [self.potential_dict[potential] for potential in baseline_potentials]

                #-- create an index dictionary for the new baseline values
                baseline_current_dict = {}
                for index in range(len(self.adjusted_currents)):
                    baseline_current_dict[self.adjusted_currents[index]] = index
            except:
                print('\n%s_raw_generator: Error in Baseline Value Adjustment\n'% self.spacer)


            ################################
            ### Create a linear baseline ###
            ################################

            #-- if it is a surface bound species --#
            if mass_transport == 'surface':
                try:
                    proto_baseline = np.linspace(vertex1,vertex2,len(baseline_currents))
                    linear_baseline = []
                    baseline_indeces = []

                    try:
                        #-- make sure the baseline currents to not go out of the
                        #-- boundaries of the current segment
                        for index in range(len(proto_baseline)):
                            if self.sign_dictionary[count] == 'cathodic':
                                if proto_baseline[index] <= baseline_currents[index]:
                                    linear_baseline.append(proto_baseline[index])
                                    baseline_indeces.append(index)
                            if self.sign_dictionary[count] == 'anodic':
                                if proto_baseline[index] >= baseline_currents[index]:
                                    linear_baseline.append(proto_baseline[index])
                                    baseline_indeces.append(index)
                    except:
                        print('\n%s_raw_generator: Error in Surface Phase Linear Baseline Adjustment' % self.spacer)

                    if len(linear_baseline) == 0:
                        linear_baseline = proto_baseline
                        baseline_indeces = [value for value in range(len(proto_baseline))]

                    try:
                        baseline_potentials = [baseline_potentials[index] for index in baseline_indeces]
                        baseline_currents = [baseline_currents[index] for index in baseline_indeces]

                        #-- make a new dictionary with key: value
                        #-- as baseline_current: baseline_potential
                        baseline_dict = {}
                        for index in range(len(baseline_currents)):
                            baseline_dict.setdefault(baseline_currents[index],[]).append(baseline_potentials[index])
                    except:
                        print('\n%s_raw_generator: Error in Surface Phase Baseline Value Adjustment' % self.spacer)

                    #-- create a dictionary correlating potentials to their
                    #-- index value
                    peak_index_dict = {}
                    index = 0
                    for potential in baseline_potentials:
                        peak_index_dict[potential] = index
                        index += 1

                    try:
                        # recalculate the peak potential of the baseline currents #
                        if self.sign_dictionary[count] == 'cathodic':
                            peak_potential = baseline_dict[max(baseline_currents)][0]  # extract the potential of the peak current
                        if self.sign_dictionary[count] == 'anodic':
                            peak_potential = baseline_dict[min(baseline_currents)][0]  # extract the potential of the peak current

                        baseline_peak_index = peak_index_dict[peak_potential]
                    except:
                        print('\n%s_raw_generator: Error in Surface Bound Peak Potential Extraction' % self.spacer)

                    try:
                        vertices = [(vertex1_potential,vertex1), *zip(baseline_potentials,baseline_currents), ((vertex2_potential,vertex2))]
                        vertex_data[count] = vertices
                    except:
                        print('\n%s_raw_generator: Error in Surface Bound Vertices Creation' % self.spacer)

                except:
                    print('\n%s_raw_generator: Error in Surface Bound Baseline Creation\n' % self.spacer)

            #-- if it is a solution phase species --#
            elif mass_transport == 'solution':
                try:
                    try:
                        if self.sign_dictionary[count] == 'cathodic':
                            vertex_current = min(vertex1,vertex2)
                        elif self.sign_dictionary[count] == 'anodic':
                            vertex_current = max(vertex1,vertex2)

                        # retrieve the corresponding potential of the vertex current
                        vertex_potential = self.data_dict[vertex_current][0]

                        # retrieve the corresponding index of the vertex current/potential
                        vertex_index = baseline_current_dict[vertex_current]
                    except:
                        print('\n%s_raw_generator: Error in Solution Phase Vertex Calculation' % self.spacer)

                    ############################################################
                    ### Calculate the Range to be used for baseline creation ###
                    ############################################################
                    try:
                        # set the `range` to 5% of the total points #
                        #print('Baseline Currents:',len(baseline_currents))
                        value_range = round(0.02*len(self.segment_currents))
                        #print('Value Range:',value_range)

                        # upper limit is the max index value before the index
                        # is out of the bounds of the adjusted currents
                        upper_limit = len(baseline_currents) - vertex_index

                        if vertex_index < value_range:
                            lower_range = baseline_currents[:vertex_index]
                            lower_range_potentials = baseline_potentials[:vertex_index]
                        else:
                            lower_index = vertex_index - value_range
                            lower_range = baseline_currents[lower_index:vertex_index]
                            lower_range_potentials = baseline_potentials[lower_index:vertex_index]
                        if vertex_index < upper_limit:
                            upper_index = vertex_index + value_range
                            upper_range = baseline_currents[vertex_index:upper_index]
                            upper_range_potentials = baseline_potentials[vertex_index:upper_index]

                        else:
                            upper_range = baseline_currents[vertex_index:]
                            upper_range_potentials = baseline_potentials[vertex_index:]
                    except:
                        print('\n%s_raw_generator: Error in Solution Phase Range Calculation' % self.spacer)

                    try:
                        # create a list with the entire range of current values #
                        slope_currents = []
                        for value in lower_range:
                            slope_currents.append(value)
                        for value in upper_range:
                            slope_currents.append(value)

                        # create a list of slope potentials
                        slope_potentials = []
                        for value in lower_range_potentials:
                            slope_potentials.append(value)
                        for value in upper_range_potentials:
                            slope_potentials.append(value)

                        # create a corresponding list of slope potentials #

                        # calculate the slope of the data points surrounding the vertex #
                        # (a ±5% range)
                        slope = (slope_currents[-1]-slope_currents[0])/(slope_potentials[-1]-slope_potentials[0])


                    except:
                        print('\n%s_raw_generator: Error in Solution Phase Slope Calculation' % self.spacer)

                    try:
                        # using this slope, create a linear baseline using the entire
                        # range of potentials chosen by the user
                        proto_baseline = []
                        for index in range(len(baseline_currents)):
                            proto_baseline.append(slope*(baseline_potentials[index]-vertex_potential)+vertex_current)

                        # remove any currents that are above/below the currents of the segment
                        linear_baseline = []
                        baseline_indeces = []

                        for index in range(len(proto_baseline)):

                            if self.sign_dictionary[count] == 'cathodic':
                                if proto_baseline[index] <= baseline_currents[index]:
                                    linear_baseline.append(proto_baseline[index])
                                    baseline_indeces.append(index)
                            if self.sign_dictionary[count] == 'anodic':
                                if proto_baseline[index] >= baseline_currents[index]:
                                    linear_baseline.append(proto_baseline[index])
                                    baseline_indeces.append(index)
                        if len(linear_baseline) == 0:
                            linear_baseline = proto_baseline
                            baseline_indeces = [value for value in range(len(proto_baseline))]
                    except:
                        print('\n%s_raw_generator: Error in Solution Phase Linear Baseline Calculation' % self.spacer)

                    try:
                        # adjust the baseline values to those that correspond with the
                        # linear baseline
                        try:
                            baseline_potentials = [baseline_potentials[index] for index in baseline_indeces]
                            baseline_currents = [baseline_currents[index] for index in baseline_indeces]

                            #-- make a new dictionary with key: value
                            #-- as baseline_current: baseline_potential
                            baseline_dict = {}
                            for index in range(len(baseline_currents)):
                                baseline_dict.setdefault(baseline_currents[index],[]).append(baseline_potentials[index])
                        except:
                            print('\n%s_raw_generator: Error in Solution Phase Baseline Value Adjustment' % self.spacer)

                        peak_index_dict = {}
                        index = 0
                        for potential in baseline_potentials:
                            peak_index_dict[potential] = index
                            index += 1

                        # recalculate the peak potential of the baseline currents #
                        if self.sign_dictionary[count] == 'cathodic':
                            peak_potential = baseline_dict[max(baseline_currents)][0]  # extract the potential of the peak current
                        if self.sign_dictionary[count] == 'anodic':
                            peak_potential = baseline_dict[min(baseline_currents)][0]  # extract the potential of the peak current

                        baseline_peak_index = peak_index_dict[peak_potential]
                        vertices = [(baseline_potentials[0],linear_baseline[0]), *zip(baseline_potentials,baseline_currents), ((baseline_potentials[-1],linear_baseline[-1]))]
                    except:
                        print('\n%s_raw_generator: Error in Solultion Phase Peak Potential Calculation' % self.spacer)
                except:
                    print('\n%s_raw_generator: Error in Solution Phase Baseline Creation\n' % self.spacer)

            ################################################################
            ### If the user selected Peak Height Extraction, analyze PHE ###
            ################################################################

            if SelectedOptions == 'Peak Height Extraction':
                try:
                    Peak_Height = abs(baseline_currents[baseline_peak_index] - linear_baseline[baseline_peak_index])
                    data = Peak_Height
                    print('%s%s Peak Height:' % (self.spacer,count),data)
                except:
                    print('\n%s_raw_generator: Error in PHE' % self.spacer)

            ########################################################
            ### If the user selected AUC extraction, analyze AUC ###
            ########################################################

            elif SelectedOptions == 'Area Under the Curve':

                ##################################
                ### Integrate Area Under the   ###
                ### Curve using a Riemmann Sum  ###
                ##################################
                try:
                    AUC_index = 1
                    AUC = 0
                    AUC_width = step_potential

                    while AUC_index < len(baseline_currents) - 1:
                        AUC_height = (baseline_currents[AUC_index] - linear_baseline[AUC_index])
                        AUC += (AUC_height * AUC_width)
                        AUC_index += 1

                    AUC = AUC/scan_rate
                    print('%s%s Reimann Sum calculated AUC:' % (self.spacer,count),AUC)

                    data = abs(AUC)

                except:
                    print('\n%s_raw_generator: Error in AUC\n' % self.spacer)

            ################################
            ### Save the data into lists ###
            ################################

            #-- Local Lists --#
            baseline_data[count] = [baseline_potentials, linear_baseline]
            vertex_data[count] = vertices

            #-- Global Lists --#
            data_list[self.num][count][self.index] = data
            peak_data_list[self.num][count][self.index] = peak_potential


            # * Optional * Export Baseline Data *
            #LocalFileHandle = ''.join(ExportPath +'%s baseline data E%s.txt' % (count,self.electrode))
            #with open(LocalFileHandle,'w+',encoding='utf-8', newline = '') as input:
            #    writer = csv.writer(input, delimiter = ' ')
            #    zip_list = zip(baseline_data[count][0],baseline_data[count][1])
            #    for item in zip_list:
            #        writer.writerow(item)


        ##########################################################
        ### Return data to the animate function as 'framedata' ###
        ##########################################################
        return segment_data,adjusted_segment_data,baseline_data,vertex_data


    ######################################################
    ### Retrieve the processed data from the generator ###
    ### and animate ('blit') it onto the artists       ###
    ######################################################
    def _animate(self, framedata, *args):

        if key > 0:

            #####################################################
            ### Retrieve the data from the generator function ###
            #####################################################
            segment_data, adjusted_segment_data, baseline_data, vertex_data = framedata

            ################################################################
            ### Acquire the artists for this electrode at this frequency ###
            ### and get the data that will be visualized                 ###
            ################################################################
            plots = plot_list[self.num]
            self.forward_boolean = forward_boolean.get()
            self.reverse_boolean = reverse_boolean.get()
            if self.forward_boolean:
                ########################
                ### Forwards Segment ###
                ########################
                #-- raw data --#
                forward_data = segment_data['forward']
                forward_potentials = forward_data[0]
                forward_currents = forward_data[1]

                #-- data adjusted to user chosen voltage range --#
                adjusted_forward_data = adjusted_segment_data['forward']
                adjusted_forward_potentials = adjusted_forward_data[0]
                adjusted_forward_currents = adjusted_forward_data[1]

                #-- linear baseline data --#
                forward_baseline = baseline_data['forward']
                forward_baseline_potentials = forward_baseline[0]
                forward_baseline_currents = forward_baseline[1]

                #-- set the x and y data of the artists to be animated --#
                plots[0].set_data(forward_potentials,forward_currents)
                plots[1].set_data(adjusted_forward_potentials,adjusted_forward_currents)
                plots[2].set_data(forward_baseline_potentials,forward_baseline_currents)
                plots[3].set_xy(vertex_data['forward'])

                #######################
                ### Peak Height/AUC ###
                #######################
                plots[8].set_data(self.file_list,data_list[self.num]['forward'][:len(self.file_list)])   # Forwards

                #################
                ### Voltage 1 ###
                #################
                plots[10].set_data(self.file_list,voltage1_data[self.num]['forward'][:len(self.file_list)])

                #################
                ### Voltage 2 ###
                #################
                plots[12].set_data(self.file_list,voltage2_data[self.num]['forward'][:len(self.file_list)])

                #######################
                ### Peak Potentials ###
                #######################
                plots[14].set_data(self.file_list,peak_data_list[self.num]['forward'][:len(self.file_list)])   # Forwards

            if self.reverse_boolean:
                #######################
                ### Reverse Segment ###
                #######################
                #-- raw data --#
                reverse_data = segment_data['reverse']
                reverse_potentials = reverse_data[0]
                reverse_currents = reverse_data[1]

                #-- data adjusted to user chosen voltage range --#
                adjusted_reverse_data = adjusted_segment_data['reverse']
                adjusted_reverse_potentials = adjusted_reverse_data[0]
                adjusted_reverse_currents = adjusted_reverse_data[1]

                #-- linear baseline data --#
                reverse_baseline = baseline_data['reverse']
                reverse_baseline_potentials = reverse_baseline[0]
                reverse_baseline_currents = reverse_baseline[1]

                #-- set the x and y data of the artists to be animated --#
                plots[4].set_data(reverse_potentials,reverse_currents)
                plots[5].set_data(adjusted_reverse_potentials,adjusted_reverse_currents)
                plots[6].set_data(reverse_baseline_potentials,reverse_baseline_currents)
                plots[7].set_xy(vertex_data['reverse'])

                #######################
                ### Peak Height/AUC ###
                #######################
                plots[9].set_data(self.file_list,data_list[self.num]['reverse'][:len(self.file_list)])   # Reverse

                #################
                ### Voltage 1 ###
                #################
                plots[11].set_data(self.file_list,voltage1_data[self.num]['reverse'][:len(self.file_list)])

                #################
                ### Voltage 2 ###
                #################
                plots[13].set_data(self.file_list,voltage2_data[self.num]['reverse'][:len(self.file_list)])

                #######################
                ### Peak Potentials ###
                #######################
                plots[15].set_data(self.file_list,peak_data_list[self.num]['reverse'][:len(self.file_list)])   # Forwards


            print('\n%s %d %s_animate\n' % (self.spacer,self.electrode,self.spacer))

            return plots


        else:
            file = 1
            EmptyPlots = framedata
            time.sleep(0.1)
            print('\n Yielding Empty Plots in Animation \n')
            return EmptyPlots


    def extract_peak(self, count):
        ### extract the potential of the peak current and then extract
        ### the associated index adjusted to the user-chosen range
        try:
            #-- if the peak is positive (cathodic) --#
            if self.sign_dictionary[count] == 'cathodic':
                peak_potential = self.data_dict[max(self.adjusted_currents)][0]
                peak_index = self.potential_dict[peak_potential] - self.adjusted_index_list[0]

            #-- if the peak is negative (anodic) --#
            elif self.sign_dictionary[count] == 'anodic':
                peak_potential = self.data_dict[min(self.adjusted_currents)][0]
                peak_index = self.potential_dict[peak_potential]  - self.adjusted_index_list[0]
            return peak_potential, peak_index

        except:
            print('\n%s_raw_generator: Error in extract_peak()\n' % self.spacer)


    def extract_vertex_currents(self, count):

        try:

            if sign_dictionary[count] == 'cathodic':
                #-- Split the voltammogram before and after the peak  --#
                #-- potential and extract the maxima for each besides --#
                try:
                    vertex1 = min(self.adjusted_currents[:self.peak_index])
                    vertex2 = min(self.adjusted_currents[self.peak_index:])
                except:
                    fit_half = int(len(self.adjusted_currents)/2)
                    vertex1 = min(self.adjusted_currents[:fit_half])
                    vertex2 = min(self.adjusted_currents[fit_half:])

            if sign_dictionary[count] == 'anodic':
                #-- Split the voltammogram before and after the peak  --#
                #-- potential and extract the maxima for each besides --#
                try:
                    vertex1 = max(self.adjusted_currents[:self.peak_index])
                    vertex2 = max(self.adjusted_currents[self.peak_index:])
                except:
                    fit_half = int(len(self.adjusted_currents)/2)
                    vertex1 = max(self.adjusted_currents[:fit_half])
                    vertex2 = max(self.adjusted_currents[fit_half:])

            return vertex1, vertex2

        except:
            print('\n%s_raw_generator: Error in extract_vertex_currents()\n' % self.spacer)


    def extract_vertex_potentials(self, vertex_potential, vertex_potentials, peak_potential):

        try:
            ### if multiple potential values have a corresponding current     ###
            ### that equals the vertex current, find the closest to the peak  ###
            ### and extract the associated index value                        ###

            if len(vertex_potential) == 1:
                vertex_potential = vertex_potential[0]
                index = self.potential_dict[vertex_potential] - self.adjusted_index_list[0] # and the index value
            else:
                list = []
                for vertex in vertex_potential:
                    if min(vertex_potentials) <= vertex <= max(vertex_potentials):
                        list.append(vertex)

                #-- find the closest to the peak potential --#
                vertex_check_list = {}
                for vertex in list:
                    check_value = abs(peak_potential - vertex)
                    vertex_check_list[check_value] = vertex
                check_value = min(vertex_check_list)
                vertex_potential = vertex_check_list[check_value]
                index = self.potential_dict[vertex_potential] - self.adjusted_index_list[0] # and the index value

            return vertex_potential, index

        except:
            print('\n%s_raw_generator: Error in extract_vertex_potentials()\n' % self.spacer)




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




##################################################
### Post Analysis Module for data manipulation ###
### after the completion of data analysis      ###
##################################################


class PostAnalysis(tk.Frame):
    def __init__(self, parent, container):
        global analysis_complete

        ############################
        ### Class-wide variables ###
        ############################
        self.parent = parent
        self.container = container

        #-- global boolean to control activation of this class --#
        analysis_complete = False
        self.ExportTopLevelExists = False

        #-- once completion value == electrode_count, analysis_complete --#
        #-- will be changed from False to True                          --#
        self.completion_value = 0


        ##########################################
        ### Initialize the Post Analysis Frame ###
        ##########################################
        self._initialize_frame()

    def _initialize_frame(self):

        ###################################################
        ### Initialize the Frame and create its Widgets ###
        ###################################################
        tk.Frame.__init__(self, self.parent)             # initialize the frame


        self.Title = tk.Label(self, text = 'Post Analysis', font=HUGE_FONT).grid(row=0,column=0,columnspan=2)

        DataAdjustmentFrame = tk.Frame(self, relief='groove',bd=3)
        DataAdjustmentFrame.grid(row=1,column=0,columnspan=2,pady=5, ipadx=50, padx=2, sticky = 'ns')

        self.RedrawButton = ttk.Button(DataAdjustmentFrame, text = 'Redraw Figures', command = lambda: self._draw(), width=12)
        self.RedrawButton.grid(row=3,column=0,pady=7)


        DataAdjustmentFrame.columnconfigure(0,weight=1)
        row_value = 3

        self.DataExportFrame = tk.Frame(self,relief='groove',bd=2)
        self.DataExportFrame.grid(row=row_value,column=0,pady=5,ipady=5)

        self.DataExportSettings = tk.Button(self.DataExportFrame, text = 'Data Export Settings', command = lambda: self.DataExportTopLevel)

        #---Buttons to switch between electrode frames---#
        frame_value = 0
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

        ExportSettings = tk.Frame(self)
        ExportSettings.grid(row=row_value,column=0,columnspan=2,pady=5,ipady=10)

        ExportSettingsButton = tk.Button(ExportSettings, text = 'Post Analysis Data Export', command = self.DataExportTopLevel)
        ExportSettingsButton.grid(row=0, column=0,padx=5)

        ExportSettings.columnconfigure(1,weight=1)
        ExportSettings.rowconfigure(1,weight=1)

        row_value += 1

        #--- Reset ---#
        Reset = ttk.Button(self, text='Reset', style='Fun.TButton', command = lambda: self._reset())
        Reset.grid(row=row_value, column=1,pady=5, padx=5)

        #--- Quit ---#
        QuitButton = ttk.Button(self, text='Quit Program',command=lambda: quit())
        QuitButton.grid(row=row_value,column=0,pady=5)

        for row in range(row_value):
            row += 1
            self.rowconfigure(row, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

    def DataExportTopLevel(self):

        self.win = tk.Toplevel()
        self.win.wm_title("Post Analysis Data Export")

        self.ExportTopLevelExists = True

        self.ElectrodeLabel = tk.Label(self.win, text='Select Electrodes:', font=LARGE_FONT)
        self.ElectrodeLabel.grid(row=1,column=0, sticky = 'nswe')
        self.ElectrodeCount = Listbox(self.win, relief='groove', exportselection=0, width=10, font=LARGE_FONT, height=6, selectmode = 'multiple', bd=3)
        self.ElectrodeCount.bind('<<ListboxSelect>>',self.ElectrodeCurSelect)
        self.ElectrodeCount.grid(row=2,column=0,padx=10,sticky='nswe')
        for electrode in electrode_list:
            self.ElectrodeCount.insert(END, electrode)


        ExportData = tk.Button(self.win, text = 'Export Data', command = lambda: self.PostAnalysisDataExport())
        ExportData.grid(row=3,column=0,columnspan=2)

        CloseButton = tk.Button(self.win, text = 'Close', command = lambda: self.win.destroy())
        CloseButton.grid(row=4,column=0,columnspan=2,pady=10)


    def ElectrodeCurSelect(self, evt):
        ###################################################
        ## electrode_list: list; ints                    ##
        ## electrode_dict: dict; {electrode: index}      ##
        ## electrode_count: int                          ##
        ###################################################

        self.electrode_list = [self.ElectrodeCount.get(idx) for idx in self.ElectrodeCount.curselection()]
        self.electrode_list = [int(electrode) for electrode in self.electrode_list]
        electrode_count = len(electrode_list)

        if electrode_count is 0:
            self.ElectrodeListExists = False
            self.ElectrodeLabel['fg'] = 'red'

        elif electrode_count is not 0:
            self.ElectrodeListExists = True
            self.ElectrodeLabel['fg'] = 'black'


    def _analysis_finished(self):
        global analysis_complete

        self.completion_value += 1

        if self.completion_value == electrode_count:
            analysis_complete = True

            #####################################
            ### Raise the Post Analysis Frame ###
            #####################################
            ShowFrames[PostAnalysis].tkraise()

    def _adjust_data(self):

        if SaveVar:
            text_file_export.TxtFileNormalization()

        ### Draw the readjusted data
        self._draw()


    def _draw(self):
        global peak, norm

        for num in range(electrode_count):

            ## get the figure for the electrode ##
            fig, ax = figures[num]

            Xaxis = file_list

            ################################################################
            ### Acquire the artists for this electrode at this frequency ###
            ### and get the data that will be visualized                 ###
            ################################################################
            plots = plot_list[num]

            ##########################
            ### Visualize the data ###
            ##########################

            #--- Peak Height/AUC ---#

            data = data_list[num]

            ### Draw new data ###
            ax[1,0].clear()
            peak = ax[1,0].plot(Xaxis,data,'bo',markersize=1)

            ax[2,0].clear()

            #####################
            ## Set the Y Label ##
            #####################
            ax[0,0].set_ylabel('Current/µA',fontweight='bold')
            if SelectedOptions == 'Peak Height Extraction':
                ax[1,0].set_ylabel('Peak Height/µA',fontweight='bold')
            elif SelectedOptions == 'Area Under the Curve':
                ax[1,0].set_ylabel('AUC/µC',fontweight='bold')
            ax[2,0].set_ylabel('Normalized', fontweight='bold')

            fig.canvas.draw_idle()

    def PostAnalysisDataExport(self):

        post_analysis_export = TextFileExport(electrodes=self.electrode_list)
        post_analysis_export.TxtFileNormalization(electrodes=self.electrode_list)


    def _reset(self):
        global AlreadyInitiated, PoisonPill, key, container

        self.completion_value = 0
        analysis_complete = False

        if self.ExportTopLevelExists is True:
            self.win.destroy()

        key = 0
        PoisonPill = True
        AlreadyInitiated = False # reset the start variable

        # Raise the initial user input frame
        self.show_frame(InputFrame)
        self.close_frame(PostAnalysis)

        ## Take resize weight away from the Visualization Canvas
        container.columnconfigure(1, weight=0)


    #--- Function to switch between visualization frames ---#
    def show_plot(self, frame):
        frame.tkraise()

    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()

    #####################################
    ### Destory the frames on Reset() ###
    #####################################
    def close_frame(self, cont):
        frame = ShowFrames[cont]
        frame.grid_forget()

        # close all matplotlib figures
        plt.close('all')

        # destory the frames holding the figures
        for frame in PlotValues:
            frame.destroy()

        # destory the container holding those frames
        PlotContainer.destroy()



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

    def tracking(self, file):

        index = file - 1

        if self.track_list[index] == electrode_count:

            ### Global File List
            _update_global_lists(file)

            if SaveVar:
                text_file_export.RealTimeExport(file)


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

    ################################
    ### Initialize the .txt file ###
    ################################
    def __init__(self, electrodes=None):

        if electrodes is None:
            self.electrode_list = electrode_list
        else:
            self.electrode_list = electrodes

        self.TextFileHandle = FileHandle

        TxtList = []
        TxtList.append('File')

        for num in self.electrode_list:
            if SelectedOptions == 'Peak Height Extraction':
                TxtList.append('E%s_Forwards_PeakHeight_(µA)' % (str(num)))
                TxtList.append('E%s_Reverse_PeakHeight_(µA)' % (str(num)))
                TxtList.append('E%s_Forward_Peak_Potential_(V)' % (str(num)))
                TxtList.append('E%s_Reverse_Peak_Potential_(V)' % (str(num)))

            elif SelectedOptions == 'Area Under the Curve':
                TxtList.append('E%s_Forwards_AUC_(µC)' % (str(num)))
                TxtList.append('E%s_Reverse_AUC_(µC)' % (str(num)))
                TxtList.append('E%s_Forward_Peak_Potential_(V)' % (str(num)))
                TxtList.append('E%s_Reverse_Peak_Potential_(V)' % (str(num)))

            TxtList.append('E%s_%sV_Forwards_(µA)' % (str(num),str(Voltage1)))
            TxtList.append('E%s_%sV_Forwards_(µA)' % (str(num),str(Voltage2)))
            TxtList.append('E%s_%sV_Reverse_(µA)' % (str(num),str(Voltage1)))
            TxtList.append('E%s_%sV_Reverse_(µA)' % (str(num),str(Voltage2)))

        for electrode in self.electrode_list:
            LocalFileHandle = FileHandle.split('.txt')
            LocalFileHandle = ''.join(LocalFileHandle[0] + '.txt')
            with open(LocalFileHandle,'w+',encoding='utf-8', newline = '') as input:
                writer = csv.writer(input, delimiter = ' ')
                writer.writerow(TxtList)

    #################################################################
    ### Write the data from the current file into the Export File ###
    #################################################################
    def RealTimeExport(self, _file_):

        index = _file_ - 1

        list = []
        list.append(str(_file_))

        #--- Peak Height ---#
        for num in range(electrode_count):

            #-- Forwards and Reverse PHE/AUC --#
            list.append(data_list[num]['forward'][index])
            list.append(data_list[num]['reverse'][index])

            #-- Forward and Reverse Peak Potential --#
            list.append(peak_data_list[num]['forward'][index])
            list.append(peak_data_list[num]['reverse'][index])

            #-- Voltage 1 (Forwards, Reverse) --#
            list.append(voltage1_data[num]['forward'][index])
            list.append(voltage1_data[num]['reverse'][index])

            #-- Voltage 2 (Forwards, Reverse) --#
            list.append(voltage2_data[num]['forward'][index])
            list.append(voltage2_data[num]['reverse'][index])

            LocalFileHandle = FileHandle.split('.txt')
            LocalFileHandle = ''.join(LocalFileHandle[0] + '.txt')

        print('\n%s TextFileExport' % LocalFileHandle,list)
        #--- Write the data into the .txt file ---#
        with open(LocalFileHandle,'a',encoding='utf-8', newline = '') as input:
            writer = csv.writer(input, delimiter = ' ')
            writer.writerow(list)



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
