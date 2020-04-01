
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
# 01/09/2020

# Need to update the get parameters function due to not all file types
# including the parameters
    # do I add in Settings options for updating these parameters?
#       (DONE) included multiple files types (.DTA, .txt, .csv)

# Need to add the ability to select how many amperograms the user is analyzing
# and adjust the merge dictionary function to accomodate this
#       (DONE)

# Need to update the script to be compatible with the SACMES main script and
# ElectrochemicalAnimation main script

# I could add in the option to calculate the diffusion coefficient using the
# Cotrell equation
    # need electrode area (A), # of electrons (n), and analyte concentration (C)
    # this probaly would be for solution phase reporters
    # Netz says this would be better done with a different equation

# replace ability to average data from multiple files with a settings option

# find where to place the file label update (_update_global_lists)
#       (DONE) I just deleted it. simple.

# 01/16/2020
# could i use a function to remove outliers in poor data sets?
#     03/28/20
#     I could probably utilize the covariance to do this in some way


# 02/07/2020
# I added the ability for multiple Pulses. Now I need to add the functionability
# to all functions that use retrieve_data. Should be able to make the analysis
# more modular - just add back the merger function in the published script.
# - start with create_setup and then move onto CA_Generator
#       (DONE)


# 03/11/2020
# 1. need to figure out where i am going to export data if there are two
# separate pulses in two separate folders being analyzed
#          (DONE) Exporting to the first file path selected
# 2. need to setup the SetupFrameCheckpoint function
# 3. need to add the ability to update the file_label widget in RTMF
#           (DONE)
# 4. WTF am I going to do with updating xstart and xend for different pulses
#           (DONE) set up a master dictionary that will hold variables for
#           both pulses (arbitrary number)
# 5. right now analysis_complete and FoundFilePath have no purpose

# 03/12/2020
# 1. Should add the ability to set x/yscale of axes (e.g. log-log)
# 2. How should i handle creating xstart and xend for multiple pulses?
# 3 .******* my use of GlobalAccessFunction isnt working. I cant seem to use it
# as an object that is part of self.controller. Not only that, I dont even really
# know what self.controller is. Need to figure out a mode of inheritance to allow
# all classes to access the functions without using it as an object. Alternatively
# I can just repeat the functions in the required classes but that would just be
# bad programming.
#           (DONE) by allowing every class to inherit from the masterclass

# 03/18/20
# 1. Leaving off at the exponential analysis classes. Analysis of daniels data
# isnt working.
#       - the analysis of optimal parameters (popt) and decay constant (k)
#         are currently in different functions of the same class but I am
#         currently using only one object (self.controller.exponential_func)
#         for the __call__ method.
#       (DONE) separated the curve_fitting (__call__) and application of the
#              exponential decay parameters given by the curve_fitting - func()

# 03/19/20
# extract_fit on 2630 needs to include the value for PathSelection
# (DONE)

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

#---Clear mac terminal memory---#
import os
os.system("clear && printf '\e[3J'")
import matplotlib
matplotlib.use('TkAgg')

#---Import Modules---#
import sys
import time as t
import datetime
try:
    import Tkinter as tk
    from Tkinter.ttk import *
    from Tkinter import *
    from Tknter import filedialog

except ImportError: # Python 3
    import tkinter as tk
    from tkinter.ttk import *
    from tkinter import *
    from tkinter import filedialog

from matplotlib import style
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import numpy as np
import csv
from pylab import *
from numpy import *
from scipy.interpolate import *
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy.signal import *
from itertools import *
import math
from math import log10, floor
from decimal import Decimal
from operator import truediv
import threading
from threading import Thread, Event
from collections import Counter
from queue import Queue
style.use('ggplot')

#---Filter out error warnings---#
import warnings
warnings.simplefilter('ignore', np.RankWarning)         #numpy polyfit_deg warning
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd") #RuntimeWarning


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
                        ####################
                        ### Master class ###
                        ####################
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

class Controller():

    def __init__(self, master):

        self.master = master

        self.Interval = 10
        self.numFiles = 20
        self.handle_variable = ''
        self.AlreadyInitiated = False
        self.FoundFilePath = False   ### If the user-inputted file is found

        now = datetime.datetime.now()
        day = str(now.day)
        month = str(now.month)
        year = str(now.year)
        self.ExportFileHandle = 'DataExport_%s_%s_%s.txt' % (year, month, day)

        self.delimiter_value = IntVar()
        self.delimiter_value.set(1)

        self.extension_value = IntVar()
        self.extension_value.set(1)

        self.current_column_index = 1
        self.spacing_index = 1
        self.time_column_index = 0

        #####################################
        ### Initialize the User Interface ###
        #####################################
        app = MainWindow(master=self.master, controller=self)

    ##############################
    ### Retrieve the file name ###
    ##############################
    def _retrieve_file(self, file, electrode, handle_variable, extension, mypath):

        try:

            if self.e_var == 'single':
                filename = '%s%s' % (handle_variable,extension)
                filename2 = '%s_%s' % (handle_variable,extension)
                filename3 = '%s_%d%s' % (handle_variable,file, extension)
                filename4 = '%s__%d%s' % (handle_variable,file, extension)
                filename5 = '%s_#%d%s' % (handle_variable,file, extension)
                filename6 = '%s__#%d%s' % (handle_variable,file, extension)


            elif self.e_var == 'multiple':
                filename = 'E%s_%s%s' % (electrode,handle_variable, extension)
                filename2 = 'E%s_%s_%s' % (electrode,handle_variable, extension)
                filename3 = 'E%s_%s_%d%s' % (electrode,handle_variable, file, extension)
                filename4 = 'E%s_%s__%d%s' % (electrode,handle_variable, file, extension)
                filename5 = 'E%s_%s_#%d%s' % (electrode,handle_variable, file, extension)
                filename6 = 'E%s_%s__#%d%s' % (electrode,handle_variable, file, extension)


            myfile = mypath + filename               ### path of your file
            myfile2 = mypath + filename2
            myfile3 = mypath + filename3
            myfile4 = mypath + filename4
            myfile5 = mypath + filename5
            myfile6 = mypath + filename6

            try:
                ### retrieves the size of the file in bytes
                mydata_bytes = os.path.getsize(myfile)
            except:
                try:
                    mydata_bytes = os.path.getsize(myfile2)
                    myfile = myfile2
                except:
                    try:
                        mydata_bytes = os.path.getsize(myfile3)
                        myfile = myfile3
                    except:
                        try:
                            mydata_bytes = os.path.getsize(myfile4)
                            myfile = myfile4
                        except:
                            try:
                                mydata_bytes = os.path.getsize(myfile5)
                                myfile = myfile5
                            except:
                                try:
                                    mydata_bytes = os.path.getsize(myfile6)
                                    myfile = myfile6
                                except:
                                    mydata_bytes = 1

            return myfile, mydata_bytes

        except:
            print('\n%sError in _retrieve_file\n')

    #####################################
    ### Read data from retrieved file ###
    #####################################
    def _read_data(self, myfile, electrode, spacer=None):
        try:
            ###############################################################
            ### Get the index value of the data depending on if the     ###
            ### electrodes are in the same .txt file or separate files  ###
            ###############################################################
            list_val = self._get_list_val(electrode)

            #####################
            ### Read the data ###
            #####################
            try:
                #---Preallocate Potential and Current lists---#
                with open(myfile,'r',encoding='utf-8') as mydata:
                    self.encoding = 'utf-8'

                    variables = len(mydata.readlines())
                    time = ['hold']*variables
                    currents = [0]*variables
            except:
                try:
                    #---Preallocate Potential and Current lists---#
                    with open(myfile,'r',encoding='utf-16') as mydata:
                        self.encoding = 'utf-16'

                        variables = len(mydata.readlines())
                        time = ['hold']*variables
                        currents = [0]*variables
                except:
                    print('\n%s_read_data: Could not determine encoding\n' % spacer)

            try:
                with open(myfile,'r',encoding=self.encoding) as mydata:
                    list_num = 0
                    for line in mydata:
                        #if list_num < self.pulse_length:

                        #---Extract data and dump into lists---#
                        line = line.strip()
                        check_split_list = line.split(self.delimiter)
                        check_split = check_split_list[0]
                        check_split = check_split.replace(',','')

                        try:
                            check_split = float(check_split)
                            check_split = True
                        except:
                            check_split = False

                        if check_split:

                            #--- extract the time ---#
                            seconds = check_split_list[self.time_column_index]
                            seconds = seconds.replace(',','')
                            seconds = float(seconds)

                            if seconds > 0:
                                time[list_num] = seconds

                                #---Currents---#
                                current_value = check_split_list[list_val]
                                current_value = current_value.replace(',','')
                                current_value = float(current_value)
                                current_value = abs(current_value)
                                current_value = current_value*1000000
                                currents[list_num] = current_value

                                list_num += 1


                        #--- if this is the last file for this pulse, break ---#
                        else:
                            pass

            except:
                print('\n%sError in read_data' % spacer)

            ### if there are 0's in the list (if the preallocation added to many)
            ### then remove them
            cut_value = 0
            for value in time:
                if value == 'hold':
                    cut_value += 1

            if cut_value > 0:
                time = time[:-cut_value]
                currents = currents[:-cut_value]

            return time, currents


        except:
            print('\n%sError in read_data' % spacer)

    #######################################################
    ### Adjust the data to match the user's constraints ###
    #######################################################
    def _apply_adjustment(self, time_list, current_list, value, spacer=None):
        try:

            zip_list = set(zip(time_list, current_list))

            # create a key,value dictionary with format {time: current}
            dict = {}
            for (time, current) in zip_list:
                dict[time] = current

            #-- values for adjustment of range used to
            xstart = self.pulse_adjustments[value]['xstart']
            xend = self.pulse_adjustments[value]['xend']

            adjusted_time = [time for time in time_list if xstart <= time <= xend ]
            adjusted_currents = [dict[time] for time in adjusted_time]

            return adjusted_time, adjusted_currents

        except:
            print('\n%sError in _apply_adjustment\n' % spacer)

    def _decay_analysis(self, time, currents, spacer=None):
            #####################################################################
            # Use a curve fitting function to extrapolate the                   #
            # paramers of the exponential decay                                 #
            #                                                                   #
            # popt: array = optimal parameters returned from least squares fit  #
            #      Monoexponential = [a k c]                                    #
            #      Biexponential = [a k1 b k2]                                  #
            # pcov: 2D array = estimated covariance of popt                     #
            #####################################################################
            popt, pcov, C = self.set_fit(time, currents, spacer=spacer)
            #print('\n\n%sOptimizatized Parameters' % spacer,popt)
            #print('%sCovariance:' % spacer,pcov,'\n\n')

            if self.exponential_str == 'monoexponential':
                # popt = [a,k,c)

                # Get the time constant, lambda, with (1/k) from y = A * e^(-kt) + C
                # *1000 for seconds --> milliseconds
                fitted_lifetime = [(1/popt[1])*1000]

            elif self.exponential_str == 'biexponential':
                # popt = [a,k1,b,k2)
                fitted_lifetime = [((1/popt[1])*1000),((1/popt[3])*1000)]

            print('%sFitted Lifetime (ms):' % spacer,fitted_lifetime)

            # use the curve fitting function to approximate values
            # fot A,k, and C in y = A * e^(-kt) + C
            fitted_data = []
            for x in time:
                fit = self.exponential_analysis.func(x, *popt, e = C)
                fitted_data.append(fit)

            return fitted_data, fitted_lifetime


    ####################################################################
    ### Using the data provided, calculate the lifetime of the decay ###
    ####################################################################
    def set_fit(self, time, currents, spacer=None):
        try:
            #-- create the range of points that will be used for regression --#
            adjusted_time, adjusted_currents, adjusted_dict = self.exponential_adjustment(time, currents,spacer=spacer)
            try:
                #-- find the optimized parameters for exponential decay --#
                popt, pcov, C = self.exponential_analysis(adjusted_time, adjusted_currents)

            except:
                try:
                    #-- if fail, call the other exponential analysis method --#
                    popt, pcov, C = self.alternate_analysis(adjusted_time, adjusted_currents)


                    #-- upon failure of currently chosen method switch the --#
                    #-- default method to the alternate and vice versa     --#
                    if self.exponential_str == 'monoexponential':
                        self.monoexponential.deselect()
                        self.biexponential.select()
                        self.exponential_str = 'biexponential'

                    elif self.exponential_str == 'biexponential':
                        self.monoexponential.select()
                        self.biexponential.deselect()
                        self.exponential_str = 'monoexponential'

                    print('\nChosen analysis method failed.\nSwitching to alternate method: %s\n' % self.exponential_str)

                    #-- switch the analysis method --#
                    alternate = self.exponential_analysis
                    self.exponential_analysis = self.alternate_analysis
                    self.alternate_analysis = alternate

                except:
                    print('\nSet Fit. Both Exponential Analysis Methods Failed.\n')

            return popt, pcov, C

        except:
            print('\n%sError in set_fit\n' % spacer)

    ####################################################
    ### Adjust the time data to fit the range set by ###
    ### the user for exponential decay analysis      ###
    ####################################################
    def exponential_adjustment(self, time_data, current_data, spacer=None):

        try:
            try:
                adjusted_time = [time for time in time_data if self.exp_low <= time <= self.exp_high]
            except:
                print('\nCould not apply exponential adjustments to time data\n')

            ###################################################
            ### Create a new dictionary that will correlate ###
            ### the currents to the adjusted time data      ###
            ###################################################
            dict = {x:y for (x,y) in zip(time_data, current_data)}
            adjusted_currents = [dict[time] for time in adjusted_time]
            adjusted_dict = {current: time for (time,current) in zip(adjusted_time,adjusted_currents)}

            return adjusted_time, adjusted_currents, adjusted_dict

        except:
            print('\n%sError in exponential_adjustment\n' % spacer)

    ################################
    ### Monoexponential Analysis ###
    ################################
    class monoexponential_analysis():
        def __init__(self, controller):

            self.controller = controller

        def __call__(self, time, currents):
            try:
                decay = sorted(set(zip(time,currents)))

                #-- extract the initial parameters for exponential decay --#
                A1, C1, K1 = self.controller.extract_fit(decay,spacer='spacer')

                #-- use the initial parameters to
                popt, pcov = curve_fit(self.func, time, currents, p0 = (A1, K1, C1))

                return popt, pcov, C1

            except:
                print('\nError in Monoexponential Analysis\n')

        def func(self, t, a, k, c, e = None):
            try:
                return (a * np.exp(-t * k) + c)
            except:
                print('\n%sError in Monoexponential func\n' % spacer)

    ##############################
    ### Biexponential Analysis ###
    ##############################
    class biexponential_analysis():
        def __init__(self, controller):

            self.controller = controller

        def __call__(self, time, currents):
            try:
                ###############################################
                ### Separate the data into two zip files to ###
                ### fit to a biexponential function         ###
                ###############################################
                half = math.ceil(0.5 * len(time))
                decay_1 = sorted(set(zip(time[:half],currents[:half])))
                decay_2 = sorted(set(zip(time[half:],currents[half:])))

                #-- extract the parameters for biexponential decay --#
                A, C1, K1 = self.controller.extract_fit(decay_1,spacer='spacer')

                B, C2, K2 = self.controller.extract_fit(decay_2,spacer='spacer')

                popt, pcov = curve_fit(self.func, time, currents,
                            p0 = (A, K1, B, K2))

                return popt, pcov, C2
            except:
                print('\nError in biexponential_analysis\n')

        def func(self, t, a, k1, b, k2, e = None):

            ##################################################################
            ### return parameters fitted to the curve with a minimization  ###
            ### function and based upon the suggestion of the variables p0 ###
            ##################################################################
            try:
                if e is None:
                    e = self.controller.constant

                return a * np.exp(-t * k1) + b * np.exp(-t * k2) + e
            except:
                print('\n%sError in biexponential_func\n' % spacer)

    ##################################################################
    ### Extract the parameters used for exponential decay analysis ###
    ##################################################################
    def extract_fit(self, curve, spacer=None):

        try:
            time_data = []
            currents = []
            dict = {}
            for item in curve:
                time, current = item
                time_data.append(time)
                currents.append(current)
                dict[current] = time

            ### get a initial assumption for the A value in y = Ae^(-kt) + C
            A = currents[5]

            ### get an initial assumption for the C value in y = Ae^(-kt) + C
            C = currents[-20]
            self.constant = C

            ### get the y value for the half-life
            half_life = 0.33*(float(A))

            ### find the closest real value to the theoretical half life
            closest_current = min(currents, key = lambda x:abs(x-half_life))
            closest_time = dict[closest_current]

            ## linearize both sides  with the natural log and extrapolate K
            ## in exponential equation y = Ae^(-kt) + C
            #log[(y-C)/A]/t = [log((y-c) - log(A)]/t
            try:
                K = (math.log(A/(closest_current - C))/closest_time)
            except:
                K = (math.log(A) - math.log(closest_current))/closest_time
            return A, C, K
        except:
            print('\n%sError in extract_fit\n' % spacer)

    def apply_offset(self, currents, alt):

        if alt:
            try:
                self.offset = float(self.select_offset.get())
            except:
                pass

        if min(currents) <= 0:
            if self.offset is not None:
                if abs(min(currents)) > self.offset:
                    self.offset = (min(currents) * -1) + 0.001
                    if alt:
                        self.select_offset.insert(END,self.offset)
                    else:
                        self.offset_label['text'] = self.offset
            else:
                self.offset = (min(currents) * -1) + 0.001
                if alt:
                    self.select_offset.insert(END,self.offset)
                else:
                    self.offset_label['text'] = self.offset

        if self.offset is not None:
            currents = [current+self.offset for current in currents]

        return currents

    ###########################################################################
    ### if multiple pulses are being analyzed, merge the data into a single ###
    ### dictionary. Overlapping data will be overwritten by one of the      ###
    ### dictionaries.                                                       ###
    ###########################################################################
    def merge_two_dicts(self, x, y):

        #-- any keys (time) that overlap will have items  --#
        #-- (currents) be overridden with the that from y --#
        z = x.copy()
        z.update(y)    # modifies z with y's keys and values & returns None
        dict = {x: y for (x,y) in z.items()}

        return dict

    ##############################################################
    ### Retrieve the column index value for current extraction ###
    ##############################################################
    def _get_list_val(self, electrode):
        spacer = ''.join(['       ']*electrode)
        try:
            if self.e_var == 'single':
                list_val = self.current_column_index + (electrode-1)*self.spacing_index

            elif self.e_var == 'multiple':
                list_val = self.current_column_index

            return list_val
        except:
            print('\n%sError in _get_list_val' % spacer)

    def update(self, file):

        if file not in self.file_list:
            self.file_list.append(file)

            if file <= self.numFiles:
                self.FileNumber['text'] = file

    ################################
    ### Raise a particular frame ###
    ################################
    def show_frame(self, cont):

        frame = self.ShowFrames[cont]
        frame.tkraise()

    #####################################
    ### Destory the frames on Reset() ###
    #####################################
    def close_frame(self, cont, destroy):
        frame = self.ShowFrames[cont]
        frame.grid_forget()

        if destroy:
            self.close_visualization()

    ######################################
    ### Close all visualization frames ###
    ######################################
    def close_visualization(self):

        # close all matplotlib figures
        plt.close('all')

        # destory the container holding those frames
        self.PlotContainer.destroy()

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#


#############################################
### Class that contains all of the frames ###
#############################################
class MainWindow():

    #--- Initialize the GUI ---#
    def __init__(self, master=None, controller=None, *args, **kwargs):

        # tk.Tk instance
        self.master = master
        self.controller = controller

        self.master.wm_title('Real-Time E-AB Sensing Platform')

        ## make the window size static
        self.master.resizable(width=False, height=False)

        #--- Create a frame for the UI ---#
        self.controller.container = tk.Frame(self.master,relief='flat',bd=5)
        self.controller.container.grid(row=0,column=0,padx=5,sticky='nsw')

        #--- Raise the frame for initial UI ---#
        # Key: frame handle / Value: tk.Frame object
        self.controller.ShowFrames = {}
        frame = InputFrame(self.controller.container, self.master, self.controller)
        frame.grid(row=0, column=0, sticky='nsw')

        self.controller.ShowFrames[InputFrame] = frame

        self.controller.show_frame(InputFrame)

        self._create_toolbar()

    def _create_toolbar(self):

        self.controller.menubar = tk.Menu(self.master)
        self.master.config(menu=self.controller.menubar)

        #################
        ### Edit Menu ###
        #################
        self.controller.editmenu = Menu(self.controller.menubar, tearoff=0)
        self.controller.menubar.add_cascade(label="Settings", menu=self.controller.editmenu)

        #######################################
        ### General Settings Drop Down Menu ###
        #######################################
        self.controller.editmenu.add_separator()
        self.controller.editmenu.add_command(label="Customize File Format",
                    command=lambda: self.extraction_adjustment_frame())

        #-- settings menu to adjust the minumum number of bytes a file must
        #-- have in order to pass the checkpoint for data analysis
                # important for potentiostats that preallocate hard drive
                # space for data files prior to interrogation
        self.byte_menu = tk.Menu(self.controller.menubar)

        #-- set the initial limit in bytes to filter out preinitialized
        #-- files < 2000b and the index value it is correlated to
        self.controller.byte_limit = 2000
        self.controller.byte_index = 1

        #- set the initial bite index to match the checkbutton
        #- index in the toolbar menu MainWindow.byte_menu


        self.onethousand = self.byte_menu.add_command(label = "   1000",
                    command = lambda: self.set_bytes('1000',0))
        self.twothousand = self.byte_menu.add_command(label = "✓ 2000",
                    command = lambda: self.set_bytes('2000',1))
        self.threethousand = self.byte_menu.add_command(label="   3000",
                    command = lambda: self.set_bytes('3000',2))
        self.fourthousand = self.byte_menu.add_command(label = "   4000",
                    command = lambda: self.set_bytes('4000',3))
        self.fivethousand = self.byte_menu.add_command(label = "   5000",
                    command = lambda: self.set_bytes('5000',4))
        self.fivethousand = self.byte_menu.add_command(label = "   7500",
                    command = lambda: self.set_bytes('7500',5))

        self.controller.editmenu.add_cascade(label='Byte Limit', menu=self.byte_menu)

        #####################################
        ### Pulse Settings Drop Down Menu ###
        #####################################
        self.pulse_menu = tk.Menu(self.controller.menubar)
        self.controller.menubar.add_cascade(label = 'Pulse Settings',menu=self.pulse_menu)

        #-- Average Multiple Pulses --#
        self.controller.pulse_averaging = tk.BooleanVar()
        self.controller.pulse_averaging.set(False)
        self.average_pulses = self.pulse_menu.add_checkbutton(label='Average Pulses',
                        onvalue=True,offvalue=False,variable=self.controller.pulse_averaging)

        #-- Average Multiple Steps --#
        self.controller.step_averaging = tk.BooleanVar()
        self.controller.step_averaging.set(False)
        self.average_steps = self.pulse_menu.add_checkbutton(label='Average Steps',
                        onvalue=True,offvalue=False,variable=self.controller.step_averaging)

        #-- Step analysis settings --#
        self.step_analysis = tk.Menu(self.pulse_menu)

        #- step analysis variables -#
        self.controller.step_index = 0            # step index
        self.controller.step_str = 'Forward'    # step analysis

        self.pulse_menu.add_cascade(label='Step Analysis',
                    menu=self.step_analysis)
        self.forward_step = self.step_analysis.add_command(label = "✓ Forward",
                    command = lambda: self.set_step_analysis('Forward',0))
        self.reverse_step = self.step_analysis.add_command(label = "   Reverse",
                    command = lambda: self.set_step_analysis('Reverse',1))
        self.both_steps = self.step_analysis.add_command(label = "   Both",
                    command = lambda: self.set_step_analysis('Both',2))


    def set_bytes(self, bytes, index):

        #-- reset the self.byte_menu widgets --#
        self.byte_menu.entryconfigure(index, label='✓%s' % bytes)
        self.byte_menu.entryconfigure(self.controller.byte_index,
                    label='   %s' % str(self.controller.byte_limit))

        #-- now change the current data being used --#
        self.controller.byte_limit = int(bytes)
        self.controller.byte_index = index

    def extraction_adjustment_frame(self):

        win = tk.Toplevel()
        win.wm_title("Customize File Format")

        #-- new frame --#
        row_value = 0
        local_container = tk.Frame(win, relief='groove',bd=2)
        local_container.grid(row=row_value,column=0,columnspan=2,padx=5,pady=5,ipadx=3)

        local_container_value = 0
        l = tk.Label(local_container, text="Time is in Column:")
        l.grid(row=local_container_value, column=0)

        local_container_value += 1
        self.time_column = tk.Entry(local_container, width=5)
        self.time_column.insert(END,1)
        self.time_column.grid(row=local_container_value,column=0,pady=5)

        local_container_value = 0
        l = tk.Label(local_container, text="Current is in Column:")
        l.grid(row=local_container_value, column=1)

        #-- entry for current column index selection --#
        local_container_value += 1
        self.list_val_entry = tk.Entry(local_container, width=5)
        self.list_val_entry.insert(END,2)
        self.list_val_entry.grid(row=local_container_value,column=1,pady=5)



        local_container_value += 1
        l = tk.Label(local_container, text="Multipotentiostat Settings\nSpace Between Current Columns:")
        l.grid(row=local_container_value, column=0,columnspan=2)

        #-- frameception --#
        local_container_value += 1
        inner_frame = tk.Frame(local_container)
        inner_frame.grid(row=local_container_value,column=0,columnspan=2)
        self.spacing_label = tk.Label(inner_frame,
                        text = '\t         Columns').grid(row=0,column=0)

        self.spacing_val_entry = tk.Entry(inner_frame, width=4)
        self.spacing_val_entry.insert(END,3)
        self.spacing_val_entry.grid(row=0,column=0,pady=1)

        #-- new frame --#
        row_value += 1
        box = tk.Frame(win, relief='groove',bd=2)
        box.grid(row=row_value,column=0,columnspan=2,pady=7)

        box_value = 0
        l = tk.Label(box, text="Delimiter between\ndata columns:")
        l.grid(row=box_value, column=0)

        box_value += 1
        self.space_delimiter = tk.Radiobutton(box, text='Space',
                        variable = self.controller.delimiter_value, value = 1)
        self.space_delimiter.grid(row=box_value,column=0,pady=5)

        box_value += 1
        self.tab_delimiter = tk.Radiobutton(box, text = 'Tab',
                        variable = self.controller.delimiter_value, value = 2)
        self.tab_delimiter.grid(row=box_value, column=0,pady=3)

        box_value += 1
        self.tab_delimiter = tk.Radiobutton(box, text = 'Comma',
                        variable = self.controller.delimiter_value, value = 3)
        self.tab_delimiter.grid(row=box_value, column=0,pady=3)

        box_value = 0
        l = tk.Label(box, text="File Extension")
        l.grid(row=box_value, column=1)

        box_value += 1
        self.txt_value = tk.Radiobutton(box, text='txt',
                        variable = self.controller.extension_value, value = 1)
        self.txt_value.grid(row=box_value,column=1,pady=5)

        box_value += 1
        self.csv_value = tk.Radiobutton(box, text = 'csv',
                        variable = self.controller.extension_value, value = 2)
        self.csv_value.grid(row=box_value, column=1,pady=3)

        box_value += 1
        self.dta_value = tk.Radiobutton(box, text = 'dta',
                        variable = self.controller.extension_value, value = 3)
        self.dta_value.grid(row=box_value, column=1,pady=3)


        row_value += 1
        apply_list_val = ttk.Button(win, text="Apply",
                        command=lambda: self.get_extraction_adjustments())
        apply_list_val.grid(row=row_value, column=0,pady=6)

        exit = ttk.Button(win, text="Exit", command= lambda: win.destroy())
        exit.grid(row=row_value, column=1,pady=3)

    def get_extraction_adjustments(self):

        #-- indeces for the columns containing currents and time --#
        self.controller.current_column_index = int(self.list_val_entry.get()) - 1
        self.controller.time_column_index = int(self.time_column.get()) - 1

        #-- number of columns in between currents for multipotentiostat files -#
        self.controller.spacing_index = int(self.spacing_val_entry.get())


    def set_step_analysis(self, step, index):

        #-- reset the self.step_analysis widgets --#
        self.step_analysis.entryconfigure(index, label='✓%s' % step)
        self.step_analysis.entryconfigure(self.controller.step_index,
                    label='   %s' % str(self.controller.step_str))

        #-- now change the current data being used --#
        self.controller.step_str = step
        self.controller.step_index = index

#--- Frame for User Input ---#
class InputFrame(tk.Frame):

    def __init__(self, parent, master, controller):

        tk.Frame.__init__(self, parent)             # initialize the frame

        self.parent = parent
        self.controller = controller
        self.master = master

        row_value = 0

        ##############################################
        ### Pack all of the widgets into the frame ###
        ##############################################

        #-- File Selection --#

        self.controller.PathSelection = {}

        #-- first pulse --#
        self.controller.PathSelection[1] = {}
        self.controller.PathSelection[1]['found'] = False

        self.SelectFirstFilePath = ttk.Button(self, style = 'Off.TButton',
                    text = 'Select File Path',
                    command = lambda: self.FindFile(parent,1))
        self.SelectFirstFilePath.grid(row=row_value,column=0,columnspan=3,pady=10)
        self.controller.PathSelection[1]['path'] = self.SelectFirstFilePath

        row_value += 3

        self.ImportFirstFileLabel = tk.Label(self, text = 'File Handle',
                    font=LARGE_FONT).grid(row=2,column=0,columnspan=3)
        self.ImportFileEntry = tk.Entry(self)
        self.ImportFileEntry.grid(row=row_value,column=0,columnspan=3,pady=5)
        self.ImportFileEntry.insert(END, self.controller.handle_variable)
        self.controller.PathSelection[1]['input'] = self.ImportFileEntry

        row_value += 1

        #-- second pulse (optional) --#
        self.AddSecondPulse = ttk.Button(self, style = 'Off.TButton',
                    text = 'Add Second Pulse',
                    command = lambda row_value = row_value: self.AddPulse(parent,row_value))
        self.AddSecondPulse.grid(row=row_value,column=0,columnspan=3,pady=10)

        row_value += 1

        self.SelectSecondFilePath = ttk.Button(self, style = 'Off.TButton',
                    text = 'Select Second File Path',
                    command = lambda: self.FindFile(parent,2))

        row_value += 3

        self.ImportSecondFileLabel = tk.Label(self, text = 'File Handle 2',
                    font=LARGE_FONT)
        self.SecondFileEntry = tk.Entry(self)

        row_value += 1

        #--- File Handle Input ---#
        HandleLabel = tk.Label(self, text='Exported File Handle:', font=LARGE_FONT)
        HandleLabel.grid(row=row_value,column=0,columnspan=3)
        row_value += 1

        self.ExportFileHandleEntry = ttk.Entry(self,width=20)
        self.ExportFileHandleEntry.insert(END, self.controller.ExportFileHandle)
        self.ExportFileHandleEntry.grid(row=row_value,column=0,columnspan=3,pady=5)

        row_value += 1

        self.numFiles_label = tk.Label(self, text='Enter File Limit',font=LARGE_FONT)
        self.numFiles_label.grid(row=row_value,column=0,columnspan=3)

        row_value += 1

        self.numFiles_entry = tk.Entry(self, width=10)
        self.numFiles_entry.insert(END, str(self.controller.numFiles))
        self.numFiles_entry.grid(row=row_value,column=0,columnspan=3)

        row_value += 1

        self.sample_interval_label = tk.Label(self, text="Sample Interval (ms)",
                    font=LARGE_FONT).grid(row=row_value,column=0,padx=5,pady=5)
        self.sample_interval_entry = tk.Entry(self, width=6)
        self.sample_interval_entry.insert(END,0.1)
        self.sample_interval_entry.grid(row=row_value+1,column=0,padx=5,pady=5)

        self.pulse_width_label = tk.Label(self, text="Pulse Width (s)",
                    font=LARGE_FONT).grid(row=row_value,column=1,padx=5,pady=5)
        self.pulse_width_entry = tk.Entry(self, width=6)
        self.pulse_width_entry.insert(END,0.1)
        self.pulse_width_entry.grid(row=row_value+1,column=1,padx=5,pady=5)

        row_value += 2

        ##################################
        ### Select and Edit Electrodes ###
        ##################################

        self.ElectrodeListboxFrame = tk.Frame(self)
        self.ElectrodeListboxFrame.grid(row=row_value,column=0,columnspan=2,
                        padx=10,pady=10,ipady=5, sticky = 'nsew')


        #--- parameters for handling resize ---#
        self.ElectrodeListboxFrame.rowconfigure(0, weight=1)
        self.ElectrodeListboxFrame.rowconfigure(1, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(0, weight=1)
        self.ElectrodeListboxFrame.columnconfigure(1, weight=1)

        self.ElectrodeListExists = False
        self.ElectrodeLabel = tk.Label(self.ElectrodeListboxFrame,
                        text='Select Electrodes:', font=LARGE_FONT)
        self.ElectrodeLabel.grid(row=0,column=0,columnspan=2, sticky = 'nswe')
        self.ElectrodeCount = Listbox(self.ElectrodeListboxFrame, relief='groove',
                        exportselection=0, width=10, font=LARGE_FONT,
                        height=6, selectmode = 'multiple', bd=3)
        self.ElectrodeCount.bind('<<ListboxSelect>>',self.ElectrodeCurSelect)
        self.ElectrodeCount.grid(row=1,column=0,columnspan=2,sticky='nswe')
        for electrode in range(16):
            electrode += 1
            self.ElectrodeCount.insert(END, electrode)

        self.scrollbar = Scrollbar(self.ElectrodeListboxFrame, orient="vertical")
        self.scrollbar.config(width=10,command=self.ElectrodeCount.yview)
        self.scrollbar.grid(row=1,column=1,sticky='nse')
        self.ElectrodeCount.config(yscrollcommand=self.scrollbar.set)

        #--- Option to have data for all electrodes in a single file ---#
        self.SingleElectrodeFile = ttk.Button(self.ElectrodeListboxFrame,
                        text='Multichannel', style = 'On.TButton',
                        command =  lambda: self.ElectrodeSelect('Multichannel'))
        self.SingleElectrodeFile.grid(row=2,column=0)

        #--- Option to have data for each electrode in a separate file ---#
        self.MultipleElectrodeFiles = ttk.Button(self.ElectrodeListboxFrame,
                        text='Multiplex', style = 'Off.TButton',
                        command = lambda: self.ElectrodeSelect('Multiplex'))
        self.MultipleElectrodeFiles.grid(row=2,column=1)

        row_value += 1

        start_button = ttk.Button(self, text = 'Initialize',
                    command = lambda: self.CheckPoint())
        start_button.grid(row=row_value,column=0, columnspan=3,pady=5,padx=5)
        row_value += 1

        quit_button = tk.Button(self, text = 'Quit', command = lambda: quit())
        quit_button.grid(row=row_value,column=0,columnspan=3,pady=5,padx=5)
        row_value += 1

        #--- Ask the User if they want to export the data to a .txt file ---#
        self.SaveVar = BooleanVar()
        self.SaveVar.set(False)
        self.SaveBox = Checkbutton(self, variable=self.SaveVar,
                    onvalue=True, offvalue=False,
                    text="Export Data").grid(row=row_value,column=0,columnspan=3)
        row_value += 1

    def AddPulse(self, parent, row_value):

        self.AddSecondPulse.grid_forget()

        self.SelectSecondFilePath.grid(row=row_value,column=0,columnspan=3,pady=10)
        row_value += 3

        self.ImportSecondFileLabel.grid(row=row_value,column=0,columnspan=3)
        row_value += 1

        self.SecondFileEntry.grid(row=row_value,column=0,columnspan=3,pady=5)
        self.SecondFileEntry.insert(END, self.controller.handle_variable)

        self.controller.PathSelection[2] = {}
        self.controller.PathSelection[2]['path'] = self.SelectSecondFilePath
        self.controller.PathSelection[2]['input'] = self.SecondFileEntry
        self.controller.PathSelection[2]['found'] = False


    def FindFile(self, parent, value):

        try:

            ### prompt the user to select a  ###
            ### directory for  data analysis ###
            mypath = filedialog.askdirectory(parent = parent)

            ### Path for directory in which the    ###
            ### exported .txt file will be placed  ###
            ExportPath = mypath.split('/')

            #-- change the text of the find file button to the folder the user chose --#
            DataFolder = '%s/%s' % (ExportPath[-2],ExportPath[-1])

            del ExportPath[-1]  # folder containing data

            ############################################
            ### Create the file path for data export ###
            ############################################
            if value == 1:
                ExportPath = '/'.join(ExportPath)
                ExportFilePath = ''.join(ExportPath + '/')
                self.controller.ExportFilePath = ''.join(ExportFilePath + self.controller.ExportFileHandle)

            self.controller.PathSelection[value]['mypath'] = ''.join(mypath + '/')
            self.controller.PathSelection[value]['found'] = True
            self.controller.PathSelection[value]['path']['style'] = 'On.TButton'
            self.controller.PathSelection[value]['path']['text'] = DataFolder

        except:

            self.controller.PathSelection[value] = False
            self.controller.PathSelection[value]['path']['style'] = 'Off.TButton'
            self.controller.PathSelection[value]['path']['text'] = 'Select File Path'

            print('\n\nInputPage.FindFile: Could Not Find File Path\n\n')

    #--- Electrode Selection ---#
    def ElectrodeCurSelect(self, evt):
        ###################################################
        ## electrode_list: list; ints                    ##
        ## electrode_dict: dict; {electrode: index}      ##
        ## electrode_count: int                          ##
        ###################################################

        self.controller.electrode_list = [self.ElectrodeCount.get(idx) for idx in self.ElectrodeCount.curselection()]
        self.controller.electrode_list = [int(electrode) for electrode in self.controller.electrode_list]
        self.controller.electrode_count = len(self.controller.electrode_list)

        index = 0
        self.controller.electrode_dict = {}
        for electrode in self.controller.electrode_list:
            self.controller.electrode_dict[electrode] = index
            index += 1

        if self.controller.electrode_count is 0:
            self.ElectrodeListExists = False
            self.ElectrodeLabel['fg'] = 'red'

        elif self.controller.electrode_count is not 0:
            self.ElectrodeListExists = True
            self.ElectrodeLabel['fg'] = 'black'

    def ElectrodeSelect(self, variable):

        if variable == 'Multiplex':
            self.controller.e_var = 'multiple'

            self.SingleElectrodeFile['style'] = 'Off.TButton'
            self.MultipleElectrodeFiles['style'] = 'On.TButton'

        elif variable == 'Multichannel':
            self.controller.e_var = 'single'

            self.SingleElectrodeFile['style'] = 'On.TButton'
            self.MultipleElectrodeFiles['style'] = 'Off.TButton'

    #####################################################################
    ### Check to see if the user has filled out all  required fields: ###
    ### Electrodes, Frequencies, Analysis Method, and File Path. If   ###
    ### they have, initialize the program                             ###
    #####################################################################
    def CheckPoint(self):

        #-- retrieve the values from the Customize File Format Topframe --#
        self.controller.delimiter_var = int(self.controller.delimiter_value.get())
        self.controller.extension_var = int(self.controller.extension_value.get())

        if self.controller.delimiter_var == 1:
            self.controller.delimiter = ' '
        elif self.controller.delimiter_var == 2:
            self.controller.delimiter = '\t'
        elif self.controller.delimiter_var == 3:
            self.controller.delimiter = ','

        if self.controller.extension_var == 1:
            self.controller.extension = '.txt'
        elif self.controller.extension_var == 2:
            self.controller.extension = '.csv'
        elif self.controller.extension_var == 3:
            self.controller.extension = '.DTA'

        for value in self.controller.PathSelection:
            self.controller.PathSelection[value]['handle'] = str(self.controller.PathSelection[value]['input'].get())

        self.controller.numFiles = int(self.numFiles_entry.get())
        self.controller.pulse_width = float(self.pulse_width_entry.get())          # seconds
        self.controller.sample_interval = float(self.sample_interval_entry.get())
        self.controller.pulse_length = self.controller.pulse_width/(self.controller.sample_interval/1000)   # milliseconds

        #--- If the user has indicated that text file export should be activated ---#
        self.controller.ExportFileHandle = str(self.ExportFileHandleEntry.get())
        print('Export File Handle: %s' % self.controller.ExportFileHandle)

        ### Path for directory in which the    ###
        ### exported .txt file will be placed  ###
        ExportPath = self.controller.PathSelection[1]['mypath'].split('/')

        del ExportPath[-1]  # blank
        del ExportPath[-1]  # folder containing data

        ############################################
        ### Create the file path for data export ###
        ############################################
        if value == 1:
            ExportPath = '/'.join(ExportPath)
            self.controller.ExportFilePath = ''.join(ExportPath + '/')
            self.controller.ExportFilePath = ''.join(self.controller.ExportFilePath + self.controller.ExportFileHandle)

        #-- boolean; True = activate file export
        self.controller.SaveVar = self.SaveVar.get()

        if self.controller.SaveVar:
            print('Initializing Text File Export')
            self.controller.text_file_export = True

        else:
            self.controller.text_file_export = None
            print('Text File Export Deactivated')

        ## set the resizeability of the container ##
        self.controller.container.columnconfigure(1, weight=1)

        ################################################################
        ### If all checkpoints have been met, initialize the program ###
        ################################################################
        check_ = True
        for value in self.controller.PathSelection:
            if self.controller.PathSelection[value]['found'] is False:
                check_ = False

        if check_:
            checkpoint = Verification(self.parent, self.master, self.controller)


######################################
### Verification TopLevel Instance ###
######################################
class Verification():
    def __init__(self, parent, master, controller):

        #-- Check to see if the user's settings are accurate
        #-- Search for the presence of the files. If they exist,
        #-- initialize the functions and frames for Real Time Analysis
        self.parent = parent
        self.master = master
        self.controller = controller

        #-- create a temporary checkpoint pop up frame --#
        self.win = tk.Toplevel()
        self.win.wm_title("CheckPoint")
        self.win.transient(self.parent)
        self.win.attributes('-topmost', 'true')

        title = tk.Label(self.win, text = 'Searching for files...',
                    font=HUGE_FONT).grid(row=0,column=1,
                    columnspan=2,pady=10,padx=10,sticky='news')

        row_value = 1

        # dictionaries to keep track of which checkpoint objects correlate
        # to which pulse for which electrode
        # dict[pulse][electrode] = object
        self.label_dict = {}
        self.already_verified = {}

        # for value in the number of path selections #
        for value in self.controller.PathSelection:
            pulse_label = tk.Label(self.win, text = self.controller.PathSelection[value]['handle'])
            pulse_label.grid(row=1,column=value,pady=5)

            self.label_dict[value] = {}
            self.already_verified[value] = {}

            row_value = 2
            for electrode in self.controller.electrode_list:

                electrode_label = tk.Label(self.win, text = 'E%s' % electrode,fg='red',font=LARGE_FONT)
                electrode_label.grid(row=row_value,column=value,pady=5,padx=5)

                self.label_dict[value][electrode] = electrode_label
                self.already_verified[value][electrode] = False

                row_value += 1

        self.stop = tk.Button(self.win, text = 'Stop', command = self.stop)
        self.stop.grid(row=row_value, column=0,columnspan=2,pady=5)
        self.StopSearch = False

        self.num = 0
        self.analysis_count = 0
        self.analysis_limit = self.controller.electrode_count * len(self.controller.PathSelection)
        self.electrode_limit = self.controller.electrode_count - 1

        self.pulse_verification = {}
        for value in self.controller.PathSelection:
            self.pulse_verification[value] = {}
            self.pulse_verification[value]['count'] = 0
            self.pulse_verification[value]['bool'] = False
            root.after(50,self.verify,value)

    def verify(self, value):

        self.electrode = self.controller.electrode_list[self.num]

        #######################################################
        ### Perform the next iteration of the data analysis ###
        #######################################################
        myfile, mydata_bytes = self.controller._retrieve_file(1, self.electrode,
                    self.controller.PathSelection[value]['handle'],
                    self.controller.extension, self.controller.PathSelection[value]['mypath'])

        if mydata_bytes > self.controller.byte_limit:

            #-- if in multipotentiostat settings, try and see if a column
            #-- exists for the current electrode

            if self.controller.e_var == 'single':
                check_ = self.verify_multi(myfile)
            else:
                check_ = True

            if check_:
                if not self.already_verified[value][self.electrode]:
                    self.already_verified[value][self.electrode] = True
                    if not self.StopSearch:
                        self.label_dict[value][self.electrode]['fg'] = 'green'
                        self.pulse_verification[value]['count'] += 1
                        self.analysis_count += 1

            if self.analysis_count == self.analysis_limit:
                if not self.StopSearch:
                    self.StopSearch = True
                    self.win.destroy()

                    self.master.after(10,self.proceed)

            if self.pulse_verification[value]['count'] == self.controller.electrode_count:
                self.pulse_verification[value]['bool'] = True

            else:
                if not self.StopSearch:
                    self.master.after(200,self.verify, value)

        if self.num < self.electrode_limit:
            self.num += 1
        else:
            self.num = 0

        if self.analysis_count < self.analysis_limit:
            if not self.StopSearch:
                if not self.pulse_verification[value]['bool']:
                    self.master.after(10,self.verify, value)

    def verify_multi(self, myfile):

        # changing the column index
        #---Set the electrode index value---#
        check_ = False
        with open(myfile,'r',encoding='utf-8') as mydata:
            for line in mydata:
                check_split_list = line.split(self.controller.delimiter)
                # delete any spaces that may come before the first value #
                while True:
                    if check_split_list[0] == '':
                        del check_split_list[0]
                    else:
                        break

                # delete any tabs that may come before the first value #
                while True:
                    if check_split_list[0] == ' ':
                        del check_split_list[0]
                    else:
                        break

                check_split = check_split_list[0]
                check_split = check_split.replace(',','')
                try:
                    check_split = float(check_split)
                    check_split = True
                except:
                    check_split = False

                if check_split:
                    self.controller.total_columns = len(check_split_list)
                    check_ = True
                    break


        if check_:
            list_val = self.controller.current_column + (self.electrode-1)*self.controller.spacing_index

            if list_val > self.controller.total_columns:
                return False

            else:
                return True

        else:
            print('\nverify_multi: could not find a line\nthat began with an integer\n')
            return False

    def stop(self):
        self.StopSearch = True
        self.win.destroy()

    def proceed(self):

        self.win.destroy()

        self.controller.track = Track(self.controller)

        #--- Create a container that can be created and destroyed when Start() or Reset() is called, respectively ---#
        self.controller.PlotContainer = tk.Frame(self.controller.container, relief = 'groove', bd = 3)
        self.controller.PlotContainer.grid(row=0,column=1, sticky = 'nsew')
        self.controller.PlotContainer.rowconfigure(0, weight=1)
        self.controller.PlotContainer.columnconfigure(0, weight=1)

        ####################################################
        ### Initialize the Initial Parameter Setup Frame ###
        ####################################################
        frame = SetupFrame(self.controller.container, self.controller,self.master)
        self.controller.ShowFrames[SetupFrame] = frame
        frame.grid(row=0, column=0, sticky='nsew')

        ############################################################
        ### Initialize the Initial Parameter Visualization Frame ###
        ############################################################
        frame = SetupVisualizationFrame(self.controller.PlotContainer, self.controller, self.master)
        self.controller.ShowFrames[SetupVisualizationFrame] = frame
        frame.grid(row=0,column=0,sticky='nsew')

        #---When initliazed, raise the Start Page and the plot for electrode one---#
        self.controller.show_frame(SetupFrame)                  # raises the frame for initial parameter setup
        self.controller.show_frame(SetupVisualizationFrame)     # raises the setup visualization frame


class SetupVisualizationFrame(tk.Frame):
    def __init__(self, parent, controller, master):

        tk.Frame.__init__(self, parent)

        self.master = master
        self.controller = controller

        #--- for resize ---#
        self.columnconfigure(0, weight = 2)

        Label = tk.Label(self, text='Setup' ,font=HUGE_FONT)
        Label.grid(row=0,column=0,pady=5,sticky='n')

        fig, ax = self.controller.setup_dictionary['figure']

        # and place the artists within the frame
        canvas = FigureCanvasTkAgg(fig, self)
        # initial draw call to create the artists that will be blitted
        canvas.draw()
        # does not affect size of figure within plot container
        canvas.get_tk_widget().grid(row=1,pady=3,ipady=3,sticky='news')


###############################################################
### Class for creating instances of the visualization frame ###
### which contains the matplotlib plots and artists         ###
###############################################################
class VisualizationFrame(tk.Frame):
    def __init__(self, parent, electrode, figure):
        tk.Frame.__init__(self, parent)         # Initialize the frame

        #--- for resize ---#
        self.columnconfigure(0, weight = 2)
        self.columnconfigure(1, weight = 1)
        self.rowconfigure(2, weight=2)

        ElectrodeLabel = tk.Label(self, text='Electrode %d' % electrode,
                    font=HUGE_FONT)
        ElectrodeLabel.grid(row=0,column=0,pady=5,sticky='new')

        # and place the artists within the frame
        canvas = FigureCanvasTkAgg(figure, self)
        # initial draw call to create the artists that will be blitted
        canvas.draw()
        # does not affect size of figure within plot container
        canvas.get_tk_widget().grid(row=2,padx=5,pady=6,ipady=5)

#---------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------#

class SetupFrame(tk.Frame):
    def __init__(self, parent, controller, master):

        tk.Frame.__init__(self, parent)             # initialize the frame

        self.parent = parent
        self.controller = controller
        self.master = master


        #####################################################
        ### Create the initial SetupVisualization Artists ###
        #####################################################
        self.controller.pulse_adjustments = {}
        self.controller.pulse_widgets = {}
        for value in self.controller.PathSelection:
            self.controller.pulse_adjustments[value] = {}

        self.bg_cache = {}

        self.RegressionFrame = tk.Frame(self,relief='groove',bd=5)
        self.create_setup()

        row_value = 0

        #-- Title --#
        self.Title = tk.Label(self,text = 'Initial Parameter Setup',
                    font=LARGE_FONT).grid(row=row_value,column=0,
                    columnspan=2,pady=10,sticky='news')

        row_value += 1

        #-- ability to alter the export path of the analyzed data --#
        if self.controller.text_file_export:
            self.AlterExportPathButton = ttk.Button(self,
                        text = 'Alter Export Path',
                        style = 'Off.TButton',
                        command = lambda: self.AlterExportPath())
            self.AlterExportPathButton.grid(row=row_value,
                        column=0,columnspan=2,pady=10,sticky='news')
            row_value += 1

        #############################################################
        ### Adjustment of points discarded at beginning of pulses ###
        #############################################################
        self.RegressionFrame.grid(row=row_value,column=0,columnspan=2,
                    pady=5,padx=5,ipadx=3,sticky='news')
        row_value += 1

        self.RegressionLabel = tk.Label(self.RegressionFrame,
                    text = 'Regression Analysis Parameters', font=LARGE_FONT)
        self.RegressionLabel.grid(row=0,column=0,columnspan=2,
                    pady=5,padx=5,sticky='news')

        ##############################
        ### Pulse Parameter Frames ###
        ##############################
        self.controller.pulse_widgets['frames'] = {}

        self.FirstParameterFrame = tk.Frame(self.RegressionFrame,
                    relief='groove',bd=2)
        self.FirstParameterFrame.grid(row=2,column=0,columnspan=4)
        self.controller.pulse_widgets['frames'][1] = self.FirstParameterFrame

        if len(self.controller.PathSelection) == 2:

            #-- the first and second parameter frame overlap in row 2 column 0 --#
            self.SecondParameterFrame = tk.Frame(self.RegressionFrame,
                        relief='groove',bd=2)
            self.SecondParameterFrame.grid(row=2,column=0,columnspan=4)
            self.controller.pulse_widgets['frames'][2] = self.SecondParameterFrame

            #########################################################
            ### Create buttons to switch between overlayed frames ###
            #########################################################

            #--- Buttons to switch in between long and short pulse parameter frames ---#
            self.controller.pulse_widgets['buttons'] = {}

            self.SelectFirstParameters = ttk.Button(self.RegressionFrame,
                        text = '%s' % self.controller.PathSelection[1]['handle'],
                        style = '.TButton',
                        command = lambda: self.SwitchParameterFrame(1))
            self.SelectFirstParameters.grid(row=3,column=0,
                        pady=5,padx=5,sticky='nsew')
            self.controller.pulse_widgets['buttons'][1] = self.SelectFirstParameters

            self.SelectSecondParameters = ttk.Button(self.RegressionFrame,
                        text = '%s' % self.controller.PathSelection[2]['handle'],
                        style = 'On.TButton',
                        command = lambda: self.SwitchParameterFrame(2))
            self.SelectSecondParameters.grid(row=3,column=1,
                        pady=5,padx=5,sticky='nsew')
            self.controller.pulse_widgets['buttons'][2] = self.SelectSecondParameters

        #############################################################################
        ### Parameters for points discarded at the beginning and end of the pulse ###
        #############################################################################
        for value in self.controller.PathSelection:

            frame = self.controller.pulse_widgets['frames'][value]

            #--- Title ---#
            self.PulseLabel = tk.Label(frame,text = 'Pulse Parameters',
                        font=MEDIUM_FONT).grid(row=0,column=0,columnspan=2)

            xstart_label = tk.Label(frame, text = 'xstart (s)',
                        font=MEDIUM_FONT).grid(row=1,column=0)
            xstart_entry = tk.Entry(frame, width=7)
            xstart_entry.insert(END, self.controller.pulse_adjustments[value]['xstart'])
            xstart_entry.grid(row=2,column=0)
            self.controller.pulse_adjustments[value]['xstart entry'] = xstart_entry

            #--- points discarded at the end of the pulse ---#
            xend_label = tk.Label(frame, text = 'xend (s)',
                        font=MEDIUM_FONT).grid(row=1,column=1)
            xend_entry = tk.Entry(frame, width=7)
            xend_entry.insert(END, self.controller.pulse_adjustments[value]['xend'])
            xend_entry.grid(row=2,column=1)
            self.controller.pulse_adjustments[value]['xend entry'] = xend_entry

        self.controller.offset_label.grid(row=4,column=0,columnspan=2,pady=5)
        if self.controller.offset is not None:
            self.controller.select_offset.insert(END,self.controller.offset)
        self.controller.select_offset.grid(row=5,column=0,columnspan=2)

        self.AdjustParameterButton = tk.Button(self.RegressionFrame,
                    text = 'Apply Adjustments', font=LARGE_FONT,
                    command = lambda: self.AdjustParameters())
        self.AdjustParameterButton.grid(row=6,column=0,
                    columnspan=2,pady=10,padx=10)

        #################################################################
        ### Frame for determining real-time bounds of exponential fit ###
        #################################################################
        self.ExponentialFitParameterFrame = tk.Frame(self, bd=5, relief='groove')                      # container
        self.ExponentialFitParameterFrame.grid(row=row_value,column=0,
                    columnspan=2,pady=5)
        row_value += 1

        self.ExponentialFrameLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Curve Fit Range').grid(row=0,column=0,columnspan=2,
                    pady=5,sticky='news')

        self.Exponential_LowLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Lower Limit (ms)', font = MEDIUM_FONT).grid(row=1,
                    column=0,sticky='news')
        self.Exponential_LowLimit = tk.Entry(self.ExponentialFitParameterFrame,
                    width=7)
        self.Exponential_LowLimit.grid(row=2,column=0,sticky='news')
        self.Exponential_LowLimit.insert(END, self.controller.exp_low)

        self.Exponential_HighLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Higher Limit (ms)', font = MEDIUM_FONT).grid(row=1,
                    column=1,sticky='news')
        self.Exponential_HighLimit = tk.Entry(self.ExponentialFitParameterFrame,
                    width=7)
        self.Exponential_HighLimit.grid(row=2,column=1,sticky='news')
        self.Exponential_HighLimit.insert(END, self.controller.exp_high)

        #-- set initial exponential analysis variables --#
        self.controller.exponential_analysis = self.controller.monoexponential_analysis(self.controller)
        self.controller.alternate_analysis = self.controller.biexponential_analysis(self.controller)
        self.controller.exponential_str = 'monoexponential'

        self.exp_analysis_label = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Exponential Decay Method', font = MEDIUM_FONT).grid(row=3,
                    column=0,columnspan=2,sticky='news')

        #-- Exponential Analysis Method Selection --#
        self.controller.exponential_var = tk.IntVar(value=1)

        self.controller.monoexponential = tk.Radiobutton(self.ExponentialFitParameterFrame, text='Monoexpoential',
                    variable=self.controller.exponential_var,
                    value=1, command = lambda: self.set_exponential_analysis())
        self.controller.monoexponential.grid(row=4,column=0,columnspan=2,pady=5)

        self.controller.biexponential = tk.Radiobutton(self.ExponentialFitParameterFrame, text='Biexponential',
                    variable=self.controller.exponential_var,
                    value=2,command = lambda: self.set_exponential_analysis())
        self.controller.biexponential.grid(row=5,column=0,columnspan=2,pady=5)


        self.ApplyExponentialParameters = tk.Button(self.ExponentialFitParameterFrame,
                    text = 'Apply Adjustments',
                    command = lambda: self.ApplyExponentialAdjustments())
        self.ApplyExponentialParameters.grid(row=6,column=0,columnspan=2,
                    pady=5,sticky='news')

        ######################################################
        ### Create the initial artists for parameter setup ###
        ######################################################
        self.create_setup_artists()

        ###############################
        ### Start and Reset Buttons ###
        ###############################
        self.ResetButton = tk.Button(self, text = 'Reset', font=LARGE_FONT,
                    width=7, command = lambda: self.Reset())
        self.ResetButton.grid(row=row_value,column=0,columnspan=1,
                    pady=5,padx=5,sticky='news')

        self.Start = tk.Button(self, text='Start', font=LARGE_FONT, width=7,
                    command = lambda: self.SetupFrameCheckpoint())
        self.Start.grid(row=row_value,column=1,columnspan=1,
                    pady=5,padx=5,sticky='news')

    #####################################################
    ### Create the initial SetupVisualization Artists ###
    #####################################################
    def create_setup(self):

        #-- will contain all of the data for the  --#
        #-- first file of the first electrode     --#
        self.controller.setup_dictionary = {}

        #-- set the initial verical offset for currents as None --#
        self.controller.offset_label = tk.Label(self.RegressionFrame, text = 'Current Adjustment (uA)')
        self.controller.select_offset = tk.Entry(self.RegressionFrame, width=4)
        self.controller.offset = None

        #########################
        ### Retrieve the data ###
        #########################
        electrode = self.controller.electrode_list[0]
        self.spacer = ''.join(['       ']*electrode)

        #-- iterate through all of the pulses selected by the user --#
        for value in self.controller.PathSelection:

            self.controller.setup_dictionary[value] =  {}

            myfile, mydata_bytes = self.controller._retrieve_file(1, electrode,
                        self.controller.PathSelection[value]['handle'],
                        self.controller.extension,
                        self.controller.PathSelection[value]['mypath'])

            #-- if the byte size of the file exceeds the checkpoint --#
            #-- limit then proceed to analyze its contents          --#
            if mydata_bytes > self.controller.byte_limit:

                ##############################
                ### Get data from the file ###
                ##############################
                time_list, current_list  = self.controller._read_data(myfile,
                            electrode,spacer=self.spacer)

                self.controller.setup_dictionary[value]['time'] = time_list
                self.controller.setup_dictionary[value]['currents'] = current_list
                self.controller.setup_dictionary[value]['dict']  = {x: y for (x,y) in zip(time_list, current_list)}

                #-- create initial xstart and xend adjustment variables --#
                self.controller.pulse_adjustments[value]['xstart'] = min(time_list)
                self.controller.pulse_adjustments[value]['xend'] = max(time_list)

                #-- and adjust the data (without manipulation)  --#
                #-- to the users specifications                 --#
                time, currents = self.controller._apply_adjustment(time_list,
                            current_list,value,spacer=self.spacer)

        #-- if there are multiple pulses (2), merge the data into a single dict --#
        if len(self.controller.PathSelection) > 1:
            data_dict = self.controller.merge_two_dicts(self.controller.setup_dictionary[1]['dict'],
                        self.controller.setup_dictionary[2]['dict'])

            time, currents = zip(*sorted(data_dict.items()))

        #-- if there are any negative values in the    --#
        #-- currents then raise them up to be positive --#
        currents = self.controller.apply_offset(currents, True)

        #-- save the final time and current values --#
        self.controller.setup_dictionary['time'] = time
        self.controller.setup_dictionary['currents'] = currents

        ####################################################################
        ### Set the initial range of points used for exponenential decay ###
        ####################################################################
        self.controller.exp_low = min(time)
        self.controller.exp_high = max(time)

        ########################################################
        ### Create a figure that will animate the setup data ###
        ########################################################
        self.fig = plt.figure(figsize=(9,6),constrained_layout=True) # (width, height)
        self.ax = self.fig.add_subplot(111)

        self.ax.set_ylabel('Current (µA)',fontweight='bold')
        self.ax.set_xlabel('Time (s)',fontweight='bold')

        self.ax.set_xlim(min(time),max(time))
        self.ax.set_ylim(min(currents),max(currents))

        #-- set the chronoamperogram to a loglog scale
        self.ax.set_yscale('log')
        self.ax.set_xscale('log')

        setup_figure = (self.fig, self.ax)
        self.controller.setup_dictionary['figure'] = setup_figure

    ######################################################
    ### Create the initial artists for parameter setup ###
    ######################################################
    def create_setup_artists(self):

        # loglog amperogram
        self.setup_artists = {}

        self.raw_decay, = self.ax.plot([],[], 'ko', MarkerSize=1)
        self.setup_artists['raw decay'] = self.raw_decay

        # loglog curve fit
        self.curve_fit, = self.ax.plot([],[], 'r-', MarkerSize=1)
        self.setup_artists['curve fit'] = self.curve_fit

        self.fig.canvas.draw()

    ####################################################
    ### Adjustment of setup visualization parameters ###
    ####################################################
    def AdjustParameters(self):
        try:
            #--- parameters for pulse ---#
            local_dict = {}
            for value in self.controller.PathSelection:
                local_dict[value] = {}

                self.controller.pulse_adjustments[value]['xstart'] = float(self.controller.pulse_adjustments[value]['xstart entry'].get())
                self.controller.pulse_adjustments[value]['xend'] = float(self.controller.pulse_adjustments[value]['xend entry'].get())

                ##############################
                ### Get data from the file ###
                ##############################
                time = self.controller.setup_dictionary[value]['time']
                currents = self.controller.setup_dictionary[value]['currents']

                #-- and adjust the data (without manipulation) to the users specifications --#
                time, currents = self.controller._apply_adjustment(time,
                            currents, value, spacer=self.spacer)

                local_dict[value]['time'] = time
                local_dict[value]['currents'] = currents
                local_dict[value]['dict'] = {x: y for (x,y) in zip(time, currents)}

            #-- if there are multiple pulses (2), merge the data into a single dict --#
            if len(self.controller.PathSelection) > 1:

                data_dict = self.controller.merge_two_dicts(local_dict[1]['dict'],
                            local_dict[2]['dict'])

                time, currents = zip(*sorted(data_dict.items()))

            ###########################################################
            ### Apply the user chosen offset OR create a new offset ###
            ### if any of the user's data has negative Y values     ###
            ###########################################################
            currents = self.controller.apply_offset(currents, True)

            ##############################################
            ### Apply the data to a decay function and ###
            ### extract the lifetime of the decay      ###
            ##############################################
            exp_time, exp_currents, dict = self.controller.exponential_adjustment(time,
                        currents, spacer='  ')

            try:
                fitted_data, fitted_lifetime = self.controller._decay_analysis(time,
                            currents, spacer='  ')


            except:
                print('\nAdjust Parameters: could not complete Fit\n')

            ##################################
            ### Set the data of the artist ###
            ##################################
            self.raw_decay.set_data(time, currents)
            self.curve_fit.set_data(time, fitted_data)

            if exp_print:
                #--- Write the data into the .txt file ---#
                with open('/Users/samdcurtis/Desktop/Merged Data.txt','w+',encoding='utf-8', newline = '') as input:
                    writer = csv.writer(input, delimiter = ' ')
                    writer.writerow(['Time','Currents'])
                with open('/Users/samdcurtis/Desktop/Merged Data.txt','a',encoding='utf-8', newline = '') as input:
                    for value in range(len(time)):
                        list = []
                        list.append(time[value])
                        list.append(currents[value])
                        writer = csv.writer(input, delimiter = ' ')
                        writer.writerow(list)

                with open('/Users/samdcurtis/Desktop/Exponential Data.txt','w+',encoding='utf-8', newline = '') as input:
                    writer = csv.writer(input, delimiter = ' ')
                    writer.writerow(['Time','Currents'])
                with open('/Users/samdcurtis/Desktop/Exponential Data.txt','a',encoding='utf-8', newline = '') as input:
                    for value in range(len(exp_time)):
                        list = []
                        list.append(exp_time[value])
                        list.append(exp_currents[value])
                        writer = csv.writer(input, delimiter = ' ')
                        writer.writerow(list)


            ##########################
            ### Animate the artist ###
            ##########################
            self.blit_data([self.raw_decay,self.curve_fit])

        except:
            print('\nError in SetupFrame: Apply Adjustments')

    def ApplyExponentialAdjustments(self):

        self.controller.exp_low = (float(self.Exponential_LowLimit.get()))
        self.controller.exp_high = (float(self.Exponential_HighLimit.get()))

        #-- apply the new adjustments and visualize the data --#
        self.AdjustParameters()

    def SwitchParameterFrame(self, frame):

        self.controller.pulse_widgets['frames'][frame].tkraise()
        self.controller.pulse_widgets['buttons'][frame]['style'] = 'On.TButton'

        if frame == 1:
            alternate = 2
        else:
            alternate = 1

        self.controller.pulse_widgets['buttons'][alternate]['style'] = 'Off.TButton'

    def set_exponential_analysis(self):

        #-- now change the current data being used --#
        self.exponential_value = int(self.controller.exponential_var.get())

        # Monoexponential
        if self.exponential_value == 1:
            self.controller.exponential_str = 'monoexponential'
            self.controller.exponential_analysis = self.controller.monoexponential_analysis(self.controller)
            self.controller.alternate_analysis = self.controller.biexponential_analysis(self.controller)

            self.controller.monoexponential.select()
            self.controller.biexponential.deselect()

        # Biexponential
        elif self.exponential_value == 2:
            self.controller.exponential_str = 'biexponential'
            self.controller.exponential_analysis = self.controller.biexponential_analysis(self.controller)
            self.controller.alternate_analysis = self.controller.monoexponential_analysis(self.controller)

            self.controller.monoexponential.deselect()
            self.controller.biexponential.select()


    def Reset(self):

        # Raise the initial user input frame
        self.controller.show_frame(InputFrame)
        self.controller.close_frame(SetupFrame, False)
        self.controller.close_frame(SetupVisualizationFrame, True)

    #    emptyMenu = Menu(root)
    #    root.config(menu=emptyMenu)

    ### Only for data setup; blit setup frame
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

        try:
            ### Path for directory in which the    ###
            ### exported .txt file will be placed  ###
            ExportPath = filedialog.askdirectory(parent = self.parent)
            ExportPath = ''.join(ExportPath+'/')

            self.controller.ExportFilePath = ''.join(ExportPath + self.controller.ExportFileHandle)

            DataFolder = self.controller.ExportFilePath.split('/')
            DataFolder = ''.join(DataFolder[-3]+'/'+DataFolder[-2])

            self.AlterExportPathButton['style'] = 'On.TButton'
            self.AlterExportPathButton['text'] = DataFolder

        except:
            print('Could not alter export file path')


    def SetupFrameCheckpoint(self):

        passkey = True

        if passkey:
            self.Initialize()

    def Initialize(self):

        #******************************************************************
        #******************************************************************
        # Here I can add in a checkpoint class similar to current version
        # of SACMES
        #*******************************************************************
        #******************************************************************


        #################################################
        ### Initiate the Frame for Data Visualization ###
        #################################################
        InitializeFigureCanvas(self.controller)

        #############################################
        ### Initialize Data Export Class Instance ###
        #############################################
        if self.controller.SaveVar:
            self.controller.export_data = TextFileExport(self.controller)
        else:
            self.controller.export_data = None

        #################################################
        ### Initialize the Real-Time User Input Frame ###
        #################################################
        frame = RealTimeManipulationFrame(self.controller.container, self.controller, self.master)
        frame.grid(row=0,column=0,sticky='nsew')
        self.controller.ShowFrames[RealTimeManipulationFrame] = frame
        frame.tkraise()

        #-- raise the frame for electrode 1 --#
        self.controller.show_frame(self.controller.electrode_list[0])


########################################################
### Create the figures and axes for the loglog plots ###
########################################################
class InitializeFigureCanvas():
    def __init__(self, master):

        self.controller = master

        #####################################################
        ### Create a figure for each electrode and embed  ###
        ### within it the figure containing its artists   ###
        #####################################################
        self.controller.figures = {}    # global figure and axes list
        self.controller.artist_list = {}

        #--- create a figure and axes that will render the data ---#
        for electrode in self.controller.electrode_list:
            fig, ax = self.MakeFigure(electrode)
            self.controller.figures[electrode] = {}
            self.controller.figures[electrode]['figure'] = fig
            self.controller.figures[electrode]['ax'] = ax

            #-- place the figure within a tkinter frame --#
            FrameReference = VisualizationFrame(self.controller.PlotContainer, electrode, fig)
            FrameReference.grid(row=0,column=0,sticky='news')
            FrameReference.tkraise()

            #-- save a reference for this frame to assign it to a button
            #-- in the RealTimeManipulationFrame
            self.controller.ShowFrames[electrode] = FrameReference

            electrode_artists = self.CreateArtists(ax)
            self.controller.artist_list[electrode] = electrode_artists

    def MakeFigure(self, electrode):

        ############################################
        ### Create the figure and axes that will ###
        ### contain the artists to be animated   ###
        ############################################
        fig = plt.figure(constrained_layout = True,figsize=(10.5,6)) # (width,height)

        if self.controller.exponential_str == 'monoexponential':
            ax = fig.add_gridspec(1, 2) # row, column values
        elif self.controller.exponential_str == 'biexponential':
            ax = fig.add_gridspec(1, 3) # row, column values

            #-- create a plot for Fitted Half life v File Number --#
            self.lifetime_plot2 = fig.add_subplot(ax[2])

        self.chronoamperogram = fig.add_subplot(ax[0])
        self.lifetime_plot = fig.add_subplot(ax[1])

        #-- get the data from the first file for this electrode --#
        time, currents, fitted_lifetime = self.RunInitialization(electrode)

        #-- create the artists for data visualization --#
        self.CreateArtists(ax)

        self.chronoamperogram.set_xlim(min(time),max(time))
        self.chronoamperogram.set_ylim(min(currents),max(currents))

        self.chronoamperogram.set_ylabel('Current (µA)',fontweight='bold')
        self.chronoamperogram.set_xlabel('Time (s)',fontweight='bold')

        #-- set the chronoamperogram to a loglog scale --#
        self.chronoamperogram.set_yscale('log')
        self.chronoamperogram.set_xscale('log')

        #-- create a plot for Fitted Half life v File Number --#
        self.lifetime_plot.set_ylim(0.5*fitted_lifetime[0],1.5*fitted_lifetime[0])  # Half Life
        self.lifetime_plot.set_xlim(0,int(self.controller.numFiles+0.1))      # File Number

        self.lifetime_plot.set_ylabel('Half Life (ms)',fontweight='bold')

        if self.controller.exponential_str == 'biexponential':

            #-- change the y title for the first exponential decay --#
            self.lifetime_plot.set_ylabel('Half Life 1 (ms)',fontweight='bold')

            self.lifetime_plot2.set_ylim(0.5*fitted_lifetime[1],1.5*fitted_lifetime[1])  # Half Life
            self.lifetime_plot2.set_xlim(0,int(self.controller.numFiles+0.1))      # File Number

            self.lifetime_plot2.set_ylabel('Half Life 2 (ms)',fontweight='bold')
            self.lifetime_plot2.set_xlabel('File Number',fontweight='bold')

        return fig, ax

    #####################################################################
    ### Get the data from the first file to populate widgets, set the ###
    ### axes of the SetupVisualizationFrame plots, and create initial ###
    ### parameters for data adjustment and regression analysis        ###
    #####################################################################
    def RunInitialization(self, electrode):

        spacer = ''.join([' ']*electrode)

        #-- iterate through all of the pulses selected by the user --#
        local_dictionary = {}

        for value in self.controller.PathSelection:

            local_dictionary[value] = {}

            myfile, mydata_bytes = self.controller._retrieve_file(1, electrode,
                        self.controller.PathSelection[value]['handle'],
                        self.controller.extension,
                        self.controller.PathSelection[value]['mypath'])

            #-- if the byte size of the file exceeds the checkpoint --#
            #-- limit then proceed to analyze its contents          --#
            if mydata_bytes > self.controller.byte_limit:

                ##############################
                ### Get data from the file ###
                ##############################
                time_list, current_list  = self.controller._read_data(myfile,
                            electrode,spacer=spacer)


                local_dictionary[value]  = {x: y for (x,y) in zip(time_list, current_list)}

                #-- adjust the data (without manipulation)  --#
                #-- to the users specifications             --#
                time, currents = self.controller._apply_adjustment(time_list,
                            current_list,value,spacer=spacer)

        #-- if there are multiple pulses (2), merge the data into a single dict --#
        if len(self.controller.PathSelection) > 1:
            data_dict = self.controller.merge_two_dicts(local_dictionary[1],
                        local_dictionary[2])

            time, currents = zip(*sorted(data_dict.items()))

        ###########################################################
        ### Apply the user chosen offset OR create a new offset ###
        ### if any of the user's data has negative Y values     ###
        ###########################################################
        currents = self.controller.apply_offset(currents, False)

        ##############################################
        ### Apply the data to a decay function and ###
        ### extract the lifetime of the decay      ###
        ##############################################

        fitted_data, fitted_lifetime = self.controller._decay_analysis(time,
                    currents, spacer=spacer)

        return time, currents, fitted_lifetime


    def CreateArtists(self, ax):

        ##########################################
        ### Create artists to animate the data ###
        ##########################################
        electrode_artists = {}

        # loglog amperogram
        self.raw_decay, = self.chronoamperogram.plot(1,1, 'ko', MarkerSize=1)
        electrode_artists['raw decay'] = self.raw_decay

        # loglog curve fit
        self.curve_fit, = self.chronoamperogram.plot(1,1, 'r-', MarkerSize=1)
        electrode_artists['curve fit'] = self.curve_fit

        # fitted lifetime (ms)
        self.fitted_lifetime, = self.lifetime_plot.plot(0,0, ' bo', MarkerSize=1)
        electrode_artists['fitted lifetime'] = self.fitted_lifetime

        if self.controller.exponential_str == 'biexponential':
            self.fitted_lifetime2, = self.lifetime_plot2.plot(0,0, ' bo', MarkerSize=1)
            electrode_artists['fitted lifetime 2'] = self.fitted_lifetime2

        self.controller.EmptyPlots = [self.raw_decay]

        return electrode_artists

##############################################################################
### Frame for Real-Time User Input for Data Manipulation and Visualization ###
##############################################################################
class RealTimeManipulationFrame(tk.Frame):

    def __init__(self, parent, controller, master):

        self.master = master
        self.parent = parent
        self.controller = controller

        tk.Frame.__init__(self, parent)          # initiate the Frame

        row_value = 0
        self.controller.FileLabel = tk.Label(self,text = 'File Number',font=LARGE_FONT)
        self.controller.FileLabel.grid(row=row_value,column=0,columnspan=2,
                    ipady=5,pady=5,ipadx=5,padx=5)
        row_value += 1

        self.controller.FileNumber = ttk.Label(self, text = '1',
                    font=LARGE_FONT, style='Fun.TButton')
        self.controller.FileNumber.grid(row=row_value,column=0,columnspan=2,
                    padx=5,ipadx=5,sticky='news')
        row_value += 1

        #######################################################################
        ### Real-Time adjustment of points discarded at beginning of pulses ###
        #######################################################################
        self.RegressionFrame = tk.Frame(self,relief='groove',bd=5)
        self.RegressionFrame.grid(row=row_value,column=0,columnspan=2,
                    pady=5,padx=5,ipadx=3,sticky='news')
        row_value += 1

        #--- Title ---#
        self.RegressionLabel = tk.Label(self.RegressionFrame,
                    text = 'Regression Analysis Parameters', font=LARGE_FONT)
        self.RegressionLabel.grid(row=0,column=0,columnspan=2,
                    pady=5,padx=5,sticky='news')

        #--- Create frames to hold parameters for pulse, respectively ---#
        self.PulseParameterFrame = tk.Frame(self.RegressionFrame,relief='groove',bd=2)
        self.PulseParameterFrame.grid(row=1,column=0,columnspan=2)

        pulse_row_value = 1
        if len(self.controller.PathSelection) == 2:

            self.FirstParameterFrame = tk.Frame(self.RegressionFrame,relief='groove',bd=2)
            self.FirstParameterFrame.grid(row=2,column=0,columnspan=4)
            self.controller.pulse_widgets['frames'][1] = self.FirstParameterFrame

            self.SecondParameterFrame = tk.Frame(self.RegressionFrame,relief='groove',bd=2)
            self.SecondParameterFrame.grid(row=2,column=0,columnspan=4)
            self.controller.pulse_widgets['frames'][2] = self.SecondParameterFrame

            self.SelectFirstParameters = ttk.Button(self.RegressionFrame,
                        text = '%s' % self.controller.PathSelection[1]['handle'],
                        style = '.TButton', command = lambda: self.SwitchParameterFrame(1))
            self.SelectFirstParameters.grid(row=3,column=0,pady=5,padx=5,sticky='nsew')
            self.controller.pulse_widgets['buttons'][1] = self.SelectFirstParameters

            self.SelectSecondParameters = ttk.Button(self.RegressionFrame,
                        text = '%s' % self.controller.PathSelection[2]['handle'],
                        style = 'On.TButton', command = lambda: self.SwitchParameterFrame(2))
            self.SelectSecondParameters.grid(row=3,column=1,pady=5,padx=5,sticky='nsew')
            self.controller.pulse_widgets['buttons'][2] = self.SelectSecondParameters
        else:
            self.FirstParameterFrame = tk.Frame(self.RegressionFrame,relief='groove',bd=2)
            self.FirstParameterFrame.grid(row=2,column=0,columnspan=4)
            self.controller.pulse_widgets['frames'][1] = self.FirstParameterFrame


        for value in self.controller.PathSelection:

            #############################################################################
            ### Parameters for points discarded at the beginning and end of the pulse ###
            #############################################################################

            frame = self.controller.pulse_widgets['frames'][value]

            #--- Title ---#
            self.PulseLabel = tk.Label(frame, text = 'Pulse Parameters',
                        font=MEDIUM_FONT).grid(row=0,column=0,columnspan=2)

            #--- points discarded at the beginning of pulse ---#
            xstart_label = tk.Label(frame, text = 'xstart (s)',
                        font=MEDIUM_FONT).grid(row=1,column=0)
            xstart_entry = tk.Entry(frame, width=7)
            xstart_entry.insert(END, self.controller.pulse_adjustments[value]['xstart'])
            xstart_entry.grid(row=2,column=0)
            self.controller.pulse_adjustments[value]['xstart entry'] = xstart_entry

            #--- points discarded at the end of the pulse ---#
            xend_label = tk.Label(frame, text = 'xend (s)',
                        font=MEDIUM_FONT).grid(row=1,column=1)
            xend_entry = tk.Entry(frame, width=7)
            xend_entry.insert(END, self.controller.pulse_adjustments[value]['xend'])
            xend_entry.grid(row=2,column=1)
            self.controller.pulse_adjustments[value]['xend entry'] = xend_entry

        #-- set the initial verical offset for currents as None --#
        self.controller.offset_header = tk.Label(self.RegressionFrame, text = 'Current Adjustment (uA)')
        self.controller.offset_label = tk.Label(self.RegressionFrame)
        if self.controller.offset is not None:
            self.controller.offset_label['text'] = '%.2f' % self.controller.offset
        self.controller.offset_header.grid(row=4,column=0,columnspan=2,pady=5)
        self.controller.offset_label.grid(row=5,column=0,columnspan=2)

        self.AdjustParameterButton = tk.Button(self.RegressionFrame,
                    text = 'Apply Adjustments', font=LARGE_FONT,
                    command = lambda: self.AdjustParameters())
        self.AdjustParameterButton.grid(row=6,column=0,columnspan=2,pady=10,padx=10)

        #################################################################
        ### Frame for determining real-time bounds of exponential fit ###
        #################################################################
        self.ExponentialFitParameterFrame = tk.Frame(self, bd=5, relief='groove')                      # container
        self.ExponentialFitParameterFrame.grid(row=row_value,column=0,
                    columnspan=2,pady=5)
        row_value += 1

        self.ExponentialFrameLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Curve Fit Range').grid(row=0,column=0,columnspan=2,
                    pady=5,sticky='news')

        self.Exponential_LowLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Lower Limit (ms)', font = MEDIUM_FONT).grid(row=1,
                    column=0,sticky='news')
        self.Exponential_LowLimit = tk.Entry(self.ExponentialFitParameterFrame,
                    width=7)
        self.Exponential_LowLimit.grid(row=2,column=0,sticky='news')
        self.Exponential_LowLimit.insert(END, self.controller.exp_low)

        self.Exponential_HighLabel = tk.Label(self.ExponentialFitParameterFrame,
                    text = 'Higher Limit (ms)', font = MEDIUM_FONT).grid(row=1,
                    column=1,sticky='news')
        self.Exponential_HighLimit = tk.Entry(self.ExponentialFitParameterFrame,
                    width=7)
        self.Exponential_HighLimit.grid(row=2,column=1,sticky='news')
        self.Exponential_HighLimit.insert(END, self.controller.exp_high)

        self.ApplyExponentialParameters = tk.Button(self.ExponentialFitParameterFrame,
                    text = 'Apply Adjustments',
                    command = lambda: self.ApplyExponentialAdjustments())
        self.ApplyExponentialParameters.grid(row=6,column=0,columnspan=2,
                    pady=5,sticky='news')

        #---Buttons to switch between electrode frames---#
        column_value = 0
        if self.controller.electrode_count > 1:
            for electrode in self.controller.electrode_list:
                Button = ttk.Button(self, text='E%d' % electrode,
                            command = lambda electrode=electrode: self.controller.show_frame(electrode))
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

        ###############################
        ### Start and Reset Buttons ###
        ###############################
        self.StartButton = tk.Button(self, text = 'Begin Visualization', command = lambda: self.SkeletonKey())
        self.StartButton.grid(row=row_value,column=0,pady=5)
        self.ResetButton = tk.Button(self, text = 'Reset', command = lambda: self.Reset())
        self.ResetButton.grid(row=row_value,column=1,pady=5)

    #--- Real-Time Adjustment of visualization parameters ---#
    def AdjustParameters(self):
        for value in self.controller.PathSelection:
            self.controller.pulse_adjustments[value]['xstart'] = float(self.controller.pulse_adjustments[value]['xstart entry'].get())
            self.controller.pulse_adjustments[value]['xend'] = float(self.controller.pulse_adjustments[value]['xend entry'].get())

    def ApplyExponentialAdjustments(self):

        self.controller.exp_low = (float(self.Exponential_LowLimit.get()))      # msec --> sec
        self.controller.exp_high = (float(self.Exponential_HighLimit.get()))

    def SwitchParameterFrame(self, frame):

        self.controller.pulse_widgets['frames'][frame].tkraise()
        self.controller.pulse_widgets['buttons'][frame]['style'] = 'On.TButton'

        if frame == 1:
            alternate = 2
        else:
            alternate = 1

        self.controller.pulse_widgets['buttons'][alternate]['style'] = 'Off.TButton'

    def set_exponential_analysis(self):

        #-- now change the current data being used --#
        self.exponential_value = int(self.controller.exponential_var.get())

        # Monoexponential
        if self.exponential_value == 1:
            self.controller.exponential_str = "mono"
            self.controller.exponential_analysis = self.controller.monoexponential_analysis(self.controller)
            self.controller.alternate_analysis = self.controller.biexponential_analysis(self.controller)

            self.controller.monoexponential.select()
            self.controller.biexponential.deselect()

        # Biexponential
        elif self.exponential_value == 2:
            self.controller.exponential_str = "bi"
            self.controller.exponential_analysis = self.controller.biexponential_analysis(self.controller)
            self.controller.alternate_analysis = self.controller.monoexponential_analysis(self.controller)

            self.controller.monoexponential.deselect()
            self.controller.biexponential.select()

    #--- Function to Reset and raise the user input frame ---#
    def Reset(self):

        if self.controller.AlreadyInitiated:
            self.controller.PoisonPill = True

        self.controller.AlreadyInitiated = False # reset the start variable
        self.controller.analysis_complete = False

        ## Take resize weight away from the Visualization Canvas
        self.controller.container.columnconfigure(1, weight=0)

        self.controller.close_frame(RealTimeManipulationFrame, True)
        self.controller.show_frame(InputFrame)


    def SkeletonKey(self):

        #--- Misc Lists ---#
        self.controller.file_list = []

        ### create a list to hold lifetime values
        self.controller.fitted_lifetime_list = {}

        #-- set the initial parameters for exponential function --#
        self.ApplyExponentialAdjustments()

        anim = {}
        if not self.controller.AlreadyInitiated:

            ######################################################################
            ### Initialize Animation (Visualization) for each electrode figure ###
            ######################################################################
            for electrode in self.controller.electrode_list:
                generator = CA_Generator(electrode, self.controller)
                func = generator.animate
                anim[electrode] = ElectrochemicalAnimation(self.master, self.controller,
                            generator, func, electrode, self.controller.figures[electrode]['figure'],
                            frames=self.controller.numFiles,init_func=generator.init_func)

            self.controller.AlreadyInitiated = True

            #--- reset poison pill variables --#
            self.controller.PoisonPill = False

        else:
            print('\n\nProgram has already been initiaed\n\n')

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


    ##########################################################################
    ##########################################################################
    ###   ANIMATION FUNCTION TO HANDLE ALL DATA ANALYSIS AND VISUALIZATION ###
    ##########################################################################
    ##########################################################################


class ElectrochemicalAnimation():
    def __init__(self, master, controller, generator, func, electrode,fig, frames = None,
                        init_func = None,resize_interval = None,
                        Interval = None, redraw = None, fargs = None):

        # object controller
        self.controller = controller

        # tk.Tk instance
        self.master = master

        # Electrode for this class instance
        self.electrode = electrode

        # Spacer value for print statements
        self.spacer = ''.join(['       ']*self.electrode)

        if frames is not None:
            self.frames = frames
        else:
            self.frames = None

        # interval at which generator callbacks are called
        if Interval is not None:
            self.Interval = Interval
        else:
            self.Interval = 10

        # (optional) generator initialization funtion
        if init_func is not None:
            self.init_func = init_func
        else:
            self.init_func = []

        ##############################
        ## Set the generator object ##
        ##############################
        self.generator = generator

        ################################
        ## Set the animation function ##
        ################################
        self._func = func

        # resize interval for matplotlib (animation) plots
        if redraw is not None:
            if resize_interval is not None:
                self.resize_interval = resize_interval
                self.redraw_figures = redraw
            else:
                self.resize_interval = None
        else:
            self.resize_interval = None

        self.resize_limit = self.resize_interval        # set the first limit

        # (optional) arguments to be passed to the animation function
        if fargs:
            self._args = fargs
        else:
            self._args = ()

        # matplotlib figure that contains the artists to be animated
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

        class _threaded_animation(threading.Thread):
            def __init__(self, parent, master, controller, Interval, frames = None):

                threading.Thread.__init__(self)

                self.parent = parent
                self.master = master
                self.controller = controller

                self.Interval = Interval

                if frames is not None:
                    self.frames = frames
                else:
                    self.frames = 1
                self.frame = 1

                # initiate the run() method
                self.master.after(1,self.start)

            def run(self):

                if not self.controller.PoisonPill:
                    if self.frame <= self.frames:
                        self.master.after(self.Interval,self.parent._step)

                        self.frame += 1

                        self.master.after(10, self.run)

        threaded_animation = _threaded_animation(self, self.master, self.controller,
                self.Interval,frames = self.frames)


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

        self._drawn_artists = self.init_func()

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


    def _step(self):

        try:
            framedata = self.generator()

            if framedata is not False:
                self._draw_next_frame(framedata)
            else:
                root.after(100, self._step)
                return False

        except StopIteration:
            return False

                    ##############################
                    ##############################
                    ### END OF ANIMATION CLASS ###
                    ##############################
                    ##############################



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

##############################################################################
### Chronoamperometry Generator to analyze data and send it to be animated ###
##############################################################################

class CA_Generator():
    def __init__(self, electrode, master, fargs = None):

        self.controller = master

        self.electrode = electrode              # electrode number
        self.num = self.controller.electrode_dict[electrode]    # electrode index
        self.file = 1       # file number
        self.index = 0      # file index

        self.local_dict = {}

        #-- Set the artists of the electrode for animation --#
        self.raw_decay = self.controller.artist_list[electrode]['raw decay']
        self.curve_fit = self.controller.artist_list[electrode]['curve fit']
        self.fitted_lifetime = self.controller.artist_list[electrode]['fitted lifetime']
        if self.controller.exponential_str == 'biexponential':
            self.fitted_lifetime2 = self.controller.artist_list[electrode]['fitted lifetime 2']

        #-- path analysis variables --#
        # total number of pulses being analyzed
        self.path_number = len(self.controller.PathSelection)
        self.path_index = 1 # track the index value of the current pulse being analyzed
        self.path_count = 0 # tracks number of pulses that have completed analysis
        # dictionary for switching between pulses for analysis
        self.path_progress = {}
        self.switch_path = {}
        if self.path_number == 2:
            self.path_progress[1] = False
            self.path_progress[2] = False
            self.switch_path[1] = 2
            self.switch_path[2] = 1
        # if you dont have multiple pulses there is no need to switch
        else:
            self.path_progress[1] = False
            self.switch_path[1] = 1

        # Spacer value for print statements
        self.spacer = ''.join(['       ']*self.electrode)

        self.controller.fitted_lifetime_list[self.num] = []

    def init_func(self):

        return self.controller.EmptyPlots

    ###############################################
    ### Generator Function for data acquisition ###
    ###############################################
    def __call__(self):


        if self.file <= self.controller.numFiles:

            self.controller.update(self.file)

            #######################################################
            ### Perform the next iteration of the data analysis ###
            #######################################################

            myfile, mydata_bytes = self.controller._retrieve_file(self.file,
                        self.electrode,self.controller.PathSelection[self.path_index]['handle'],
                        self.controller.extension,self.controller.PathSelection[self.path_index]['mypath'])

            if mydata_bytes > self.controller.byte_limit:

                time, currents = self._generator(myfile, self.path_index)

                #-- Compile data into local dictionary --#
                self.local_dict[self.path_index] = {}

                self.local_dict[self.path_index]['time'] = time
                self.local_dict[self.path_index]['currents'] = currents
                self.local_dict[self.path_index]['dict'] = {x: y for (x,y) in zip(time,currents)}

                #-- after finding the file and compiling the --#
                #-- data increase the checkpoint value       --#
                self.path_progress[self.path_index] = True
                self.path_count += 1

            ################################################################
            ### Once all of the pulses have been analyzed, move forwards ###
            ### with exponential decay analysis and data animation       ###
            ################################################################
            if self.path_count == self.path_number:

                #-- Reset the path_count checkpoint variable --#
                self.path_count = 0

                #-- if there are multiple pulses (2), merge the data into a single dict --#
                if self.path_number > 1:
                    data_dict = self.controller.merge_two_dicts(self.local_dict[1]['dict'],
                                            self.local_dict[2]['dict'])

                    #-- create a correlation dictionary --#
                    time, currents = zip(*sorted(data_dict.items()))

                else:
                    time = self.local_dict[1]['time']
                    currents = self.local_dict[1]['currents']

                #-- apply the vertical offset --#
                currents = self.controller.apply_offset(currents, False)

                ##############################################
                ### Apply the data to a decay function and ###
                ### extract the lifetime of the decay      ###
                ##############################################
                fitted_data, fitted_lifetime = self.controller._decay_analysis(time,
                            currents, spacer=self.spacer)

                self.controller.fitted_lifetime_list[self.num].append(fitted_lifetime)
                self.controller.track(self.file)

                return time, currents, fitted_data

            else:
                #-- switch paths (if applicable) --#
                path_index = self.switch_path[self.path_index]

                #-- if this path has not completed     --#
                #-- analysis then commit to the switch --#
                if not self.path_progress[path_index]:
                    self.path_index = path_index

                return False

        else:
            return False


    #--- Generator function to return data to analyze ---#
    def _generator(self, myfile, value):

        try:

            print('\n%sElectrode %d: Analyzing File %d\n' % (self.spacer,
                        self.electrode,self.file))

            ##############################
            ### Get data from the file ###
            ##############################
            time, currents  = self.controller._read_data(myfile,
                        self.electrode,spacer=self.spacer)

            #-- and adjust the data (without manipulation) to the users specifications --#
            time, currents = self.controller._apply_adjustment(time,
                        currents, value, spacer=self.spacer)

            return time, currents

        except:
            print('\n%sError in _generator\n' % self.spacer)


    def animate(self, framedata):

        try:
            print('\n%sanimate\n' % self.spacer)
            time_data, current_data, fitted_data = framedata

            # create a list of artists to be animated
            animation_list = []

            # set the chronoamperogram
            self.raw_decay.set_data(time_data,current_data)     # raw curve
            animation_list.append(self.raw_decay)

            self.curve_fit.set_data(time_data,fitted_data)  # set the curve fit as an artist
            animation_list.append(self.curve_fit)

            # set the half life plots
            file_list = range(1,self.file+1)
            if self.controller.exponential_str == 'monoexponential':
                self.fitted_lifetime.set_data(file_list, self.controller.fitted_lifetime_list[self.num])
                animation_list.append(self.fitted_lifetime)

            elif self.controller.exponential_str == 'biexponential':

                fitted_lifetime_one = [x[0] for x in self.controller.fitted_lifetime_list[self.num]]
                fitted_lifetime_two = [x[1] for x in self.controller.fitted_lifetime_list[self.num]]

                self.fitted_lifetime.set_data(file_list, fitted_lifetime_one)
                self.fitted_lifetime2.set_data(file_list, fitted_lifetime_two)

                animation_list.append(self.fitted_lifetime)
                animation_list.append(self.fitted_lifetime2)

            ################################
            ### Move on to the next file ###
            ################################
            if self.file < self.controller.numFiles:
                self.file += 1

            ### return artists to be animated
            return animation_list

        except:
            return self.controller.EmptyPlots
            print('\nError in animate\n')

class Track():
    def __init__(self, master):

        self.controller = master

        self.track_list = [1]*self.controller.numFiles

    def __call__(self, file):

        index = file - 1

        if self.track_list[index] == self.controller.electrode_count:

            if self.controller.SaveVar:
                self.controller.export_data.RealTimeExport(file)

        else:
            self.track_list[index] += 1

##################################
### Real-Time Text File Export ###
##################################
class TextFileExport():

    ###############################
    ### Initialize the .txt file ###
    ###############################
    def __init__(self, master):

        self.controller = master

        try:
            list = []
            list.append('File')
            for electrode in self.controller.electrode_list:
                list.append('E%s_Lifetime(ms)' % str(electrode))
                if self.controller.exponential_str == 'biexponential':
                    list.append('E%s_Lifetime_2(ms)' % str(electrode))


            #--- Write the data into the .txt file ---#
            with open(self.controller.ExportFilePath,'w+',encoding='utf-8', newline = '') as input:
                writer = csv.writer(input, delimiter = ' ')
                writer.writerow(list)

        except:
            print('\n\n','ERROR IN TEXT FILE EXPORT','\n\n')
            time.sleep(0.5)

    def RealTimeExport(self, file):

        index = file - 1

        list = [file]

        for electrode in self.controller.electrode_list:
            num = self.controller.electrode_dict[electrode]
            for item in self.controller.fitted_lifetime_list[num][index]:
                list.append(item)

        #--- Write the data into the .txt file ---#
        with open(self.controller.ExportFilePath,'a',encoding='utf-8', newline = '') as input:
            writer = csv.writer(input, delimiter = ' ')
            writer.writerow(list)
        with open(self.controller.ExportFilePath,'r',encoding='utf-8', newline = '') as filecontents:
            filedata =  filecontents.read()
        filedata = filedata.replace('[','')
        filedata = filedata.replace('"','')
        filedata = filedata.replace(']','')
        filedata = filedata.replace(',','')
        filedata = filedata.replace('\'','')
        with open(self.controller.ExportFilePath,'w',encoding='utf-8', newline = '') as output:
            output.write(filedata)




#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

                    ############################################
                    ### Initialize GUI to start the program  ###
                    ############################################

if __name__ == '__main__':

    #-- Global Editing Variables for managers --#
    exp_print = False

    ###############
    ### Styling ###
    ###############
    HUGE_FONT = ('Verdana', 18)
    LARGE_FONT = ('Verdana', 11)
    MEDIUM_FONT = ('Verdnana', 10)
    SMALL_FONT = ('Verdana', 8)

    root = tk.Tk()
    app = Controller(master=root)

    #--- initiate the mainloop ---#
    try:
        root.mainloop()

        style = ttk.Style()
        style.configure('On.TButton', foreground = 'blue', font = LARGE_FONT, relief = 'raised', border = 100)
        style.configure('Off.TButton', foreground = 'black', font = MEDIUM_FONT, relief = 'sunken', border = 5)

    #--- escape scrolling error ---#
    except UnicodeDecodeError:
        pass

                    #*########################################*#
                    #*############ End of Program ############*#
                    #*########################################*#
