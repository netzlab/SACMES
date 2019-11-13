# chronoamp script #


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

#---Clear mac terminal memory---#
import os
import matplotlib
matplotlib.use('TkAgg')
os.system("clear && printf '\e[3J'")

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


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

handle_variable = ''
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

                                ########################
                                ### Global Variables ###
                                ########################

Interval = 100
file_limit = 52
point_limit = 20
short_aa = 10**(-3)         # points removed from the the beginning of the short pulse
short_zz = 10**(-2)         # points removed from the the end of the short pulse
long_aa = 10**(-2)            # points removed from the the beginning of the long pulse
long_zz = 1                # points removed from the the end of the long pulse

AlreadyInitiated = False

#--- Styling ---#
HUGE_FONT = ('Verdana', 18)
LARGE_FONT = ('Verdana', 12)
MEDIUM_FONT = ('Verdnana', 10)
SMALL_FONT = ('Verdana', 8)


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#



#############################################
### Class that contains all of the frames ###
#############################################
class RealTimeAnalysis(tk.Tk):

    #--- Initialize the GUI ---#
    def __init__(self, *args, **kwargs):
        global container, PlotContainer, fig, ax, ShowFrames

        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self,'Real-Time E-AB Sensing Platform')
        tk.Tk.resizable(self,width=False, height=False)             ## make the window size static

        #--- Create a frame for the UI ---#
        container = tk.Frame(self,relief='flat',bd=5).grid(row=0,column=0,padx=5,sticky='nsw')

        #--- Raise the frame for initial UI ---#
        ShowFrames = {}                                 # Key: frame handle / Value: tk.Frame object
        frame = StartPage(container, self)
        ShowFrames[StartPage] = frame
        frame.grid(row=0, column=0, sticky='nsw')
        self.show_frame(StartPage)


    #--- Function to visualize different frames ---#
    def show_frame(self, cont):

        frame = ShowFrames[cont]
        frame.tkraise()


    def show_plot(self, frame):
        frame.tkraise()


#--- Frame for User Input ---#
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        global point_limit_variable, file_limit_entry, point_limit_entry

        tk.Frame.__init__(self, parent)             # initialize the frame

        self.SelectFilePath = ttk.Button(self, style = 'Off.TButton', text = 'Select File Path', command = lambda: self.FindFile(parent))
        self.SelectFilePath.grid(row=0,column=0,columnspan=3)

        self.NoSelectedPath = tk.Label(self, text = 'No File Path Selected', font = MEDIUM_FONT, fg = 'red')
        self.PathWarningExists = False               # tracks the existence of a warning label

        ImportFileLabel = tk.Label(self, text = 'Data File Handle', font=LARGE_FONT).grid(row=2,column=0,columnspan=3)
        self.ImportFileEntry = tk.Entry(self)
        self.ImportFileEntry.grid(row=3,column=0,columnspan=3,pady=5)
        self.ImportFileEntry.insert(END, handle_variable)

        #--- File Handle Input ---#
        HandleLabel = tk.Label(self, text='Exported File Handle:', font=LARGE_FONT)
        HandleLabel.grid(row=4,column=0,columnspan=3)
        self.filehandle = ttk.Entry(self)
        now = datetime.datetime.now()
        hour = str(now.hour)
        day = str(now.day)
        month = str(now.month)
        year = str(now.year)
        self.filehandle.insert(END, 'DataExport_%s_%s_%s.txt' % (year, month, day))
        self.filehandle.grid(row=5,column=0,columnspan=3,pady=5)

        file_limit_label = tk.Label(self, text='Enter File Limit',font=LARGE_FONT)
        file_limit_label.grid(row=7,column=0,columnspan=3)

        file_limit_entry = tk.Entry(self, width=10)
        file_limit_entry.insert(END, str(file_limit))
        file_limit_entry.grid(row=10,column=0,columnspan=3)

        point_limit_label = tk.Label(self, text='Enter Number of Points to Average', font=LARGE_FONT)
        point_limit_label.grid(row=15,column=0,columnspan=3)

        point_limit_entry = tk.Entry(self, width=10)
        point_limit_entry.insert(END, str(point_limit))
        point_limit_entry.grid(row=20,column=0,columnspan=3)

        start_button = ttk.Button(self, text = 'Initiate Visualization', command = self.StartProgram).grid(row=25,column=0,columnspan=3,pady=5,padx=5)
        quit_button = tk.Button(self, text = 'Quit', command = lambda: quit()).grid(row=30,column=0,columnspan=3,pady=5,padx=5)

        #--- Ask the User if they want to export the data to a .txt file ---#
        self.SaveVar = BooleanVar()
        self.SaveVar.set(False)
        self.SaveBox = Checkbutton(self, variable=self.SaveVar, onvalue=True, offvalue=False, text="Export Data").grid(row=35,column=0,columnspan=3)


    def FindFile(self, parent):
        global mypath, ExportPath, FoundFilePath, NoSelectedPath

        try:

            ### prompt the user to select a  ###
            ### directory for  data analysis ###
            mypath = filedialog.askdirectory(parent = parent)
            mypath = ''.join(mypath + '/')


            ### Path for directory in which the    ###
            ### exported .txt file will be placed  ###
            ExportPath = mypath.split('/')

            #-- change the text of the find file button to the folder the user chose --#
            DataFolder = '%s/%s' % (ExportPath[-3],ExportPath[-2])

            self.SelectFilePath['style'] = 'On.TButton'
            self.SelectFilePath['text'] = DataFolder



            del ExportPath[-1]
            del ExportPath[-1]
            ExportPath = '/'.join(ExportPath)
            ExportPath = ''.join(ExportPath + '/')
            print(ExportPath)

            ## Indicates that the user has selected a File Path ###
            FoundFilePath = True

            if self.PathWarningExists:
                self.NoSelectedPath['text'] = ''
                self.NoSelectedPath.grid_forget()

        except:
            print('\n\nInputPage.FindFile: Could Not Find File Path\n\n')



    def StartProgram(self):
        global initiate_figures, FileHandle, PlotContainer, handle_variable, export_data, SaveVar

        #--- Initiate the Frame for Data Visualization ---#
        initiate_figures = InitiateFigures()
        handle_variable = self.ImportFileEntry.get()
        SaveVar = self.SaveVar.get()
        FileHandle = self.filehandle.get()


        export_data = TextFileExport()


        FileHandle = str(self.filehandle.get()) # handle for exported .txt file
        FileHandle = ''.join(ExportPath + FileHandle)
        print('FileHandle',FileHandle)

        ##############################################################################
        ### Frame for Real-Time User Input for Data Manipulation and Visualization ###
        ##############################################################################
        class RealTimeInputFrame(tk.Frame):

            def __init__(self, parent, controller):
                global pulse_label, file_label

                tk.Frame.__init__(self, parent)          # initiate the Frame

                file_header = tk.Label(self,text = 'File Number',font=LARGE_FONT)
                file_header.grid(row=0,column=0,ipady=5,pady=5,ipadx=5,padx=5,)
                file_label = ttk.Label(self, text = '1', font=LARGE_FONT, style='Fun.TButton')
                file_label.grid(row=1,column=0,padx=10,ipadx=10)

                pulse_header = tk.Label(self, text= 'Point Number',font=LARGE_FONT)
                pulse_header.grid(row=0,column=1,ipady=5,pady=5,ipadx=5,padx=5)
                pulse_label = ttk.Label(self, text = '1', font=LARGE_FONT, style='Fun.TButton')
                pulse_label.grid(row=1,column=1,padx=10,ipadx=10)


                #######################################################################
                ### Real-Time adjustment of points discarded at beginning of pulses ###
                #######################################################################
                RegressionFrame = tk.Frame(self,relief='groove',bd=5)
                RegressionFrame.grid(row=2,column=0,columnspan=4,pady=5,padx=5,ipadx=3)

                #--- Title ---#
                self.RegressionLabel = tk.Label(RegressionFrame, text = 'Regression Analysis Parameters', font=LARGE_FONT)
                self.RegressionLabel.grid(row=0,column=0,columnspan=4,pady=5,padx=5)

                #--- Create frames to hold parameters for short and long pulses, respectively ---#
                self.ShortParameterFrame = tk.Frame(RegressionFrame,relief='groove',bd=2)             # Frame to hold the parameters for frequencies <= 50Hz
                self.ShortParameterFrame.grid(row=1,column=0,columnspan=4)
                self.LongParameterFrame = tk.Frame(RegressionFrame,relief='groove',bd=2)            # Frame to hold the parameters for frequencies > 50Hz
                self.LongParameterFrame.grid(row=1,column=0,columnspan=4)

                #--- Buttons to switch in between long and short pulse parameter frames ---#
                self.SelectShortParameters = ttk.Button(RegressionFrame, text = 'Short Parameters', style = 'Off.TButton', command = lambda: self.SwitchParameterFrame('Short'))
                self.SelectShortParameters.grid(row=2,column=0,pady=5,padx=5,sticky='nsew')
                self.SelectLongParameters = ttk.Button(RegressionFrame, text = 'Long Parameters', style = 'On.TButton', command = lambda: self.SwitchParameterFrame('Long'))
                self.SelectLongParameters.grid(row=2,column=1,pady=5,padx=5,sticky='nsew')

                ##################################################################################
                ### Parameters for points discarded at the beginning and end of the short pulse ###
                ##################################################################################
                #--- Title ---#
                self.ShortLabel = tk.Label(self.ShortParameterFrame, text = 'Short Pulse Parameters', font=MEDIUM_FONT).grid(row=0,column=0,columnspan=2)

                #--- points discarded at the beginning of pulse ---#
                self.short_aa_label = tk.Label(self.ShortParameterFrame, text = 'short_aa (s)', font=MEDIUM_FONT).grid(row=1,column=0)
                self.short_aa_entry = tk.Entry(self.ShortParameterFrame, width=7)
                self.short_aa_entry.insert(END, str(short_aa))
                self.short_aa_entry.grid(row=2,column=0)
                short_aa_entry = self.short_aa_entry

                #--- points discarded at the end of the pulse ---#
                self.short_zz_label = tk.Label(self.ShortParameterFrame, text = 'short_zz (s)', font=MEDIUM_FONT).grid(row=1,column=1)
                self.short_zz_entry = tk.Entry(self.ShortParameterFrame, width=7)
                self.short_zz_entry.insert(END, str(short_zz))
                self.short_zz_entry.grid(row=2,column=1)
                short_zz_entry = self.short_zz_entry

                ##################################################################################
                ### Parameters for points discarded at the beginning and end of the long pulse ###
                ##################################################################################
                #--- Frame Title ---#
                self.LongLabel = tk.Label(self.LongParameterFrame, text = 'Long Pulse Parameters', font=MEDIUM_FONT).grid(row=0,column=0,columnspan=2)

                #--- points discarded at the beginning of the pulse ---#
                self.long_aa_label = tk.Label(self.LongParameterFrame, text = 'long_aa (s)', font=MEDIUM_FONT).grid(row=1,column=0)
                self.long_aa_entry = tk.Entry(self.LongParameterFrame, width=7)
                self.long_aa_entry.insert(END, str(long_aa))
                self.long_aa_entry.grid(row=2,column=0)
                long_aa_entry = self.long_aa_entry

                #--- points discarded at the end of the pulse ---#
                self.long_zz_label = tk.Label(self.LongParameterFrame, text = 'long_zz (s)', font=MEDIUM_FONT).grid(row=1,column=1)
                self.long_zz_entry = tk.Entry(self.LongParameterFrame, width=7)
                self.long_zz_entry.insert(END, str(long_zz))
                self.long_zz_entry.grid(row=2,column=1)
                long_zz_entry = self.long_zz_entry

                #--- Button to apply adjustments ---#
                self.AdjustParameterButton = tk.Button(RegressionFrame, text = 'Apply Adjustments', font=LARGE_FONT, command = lambda: self.AdjustParameters())
                self.AdjustParameterButton.grid(row=3,column=0,columnspan=4,pady=10,padx=10)

                #################################################################
                ### Frame for determining real-time bounds of exponential fit ###
                #################################################################
                self.ExponentialFitParameterFrame = tk.Frame(self, bd=5, relief='groove')                      # container
                self.ExponentialFitParameterFrame.grid(row=3,column=0,columnspan=3,pady=5)

                self.ExponentialFrameLabel = tk.Label(self.ExponentialFitParameterFrame, text = 'Curve Fit Range').grid(row=0,column=0,columnspan=2,pady=5)

                self.Exponential_LowLabel = tk.Label(self.ExponentialFitParameterFrame, text = 'Lower Limit (ms)', font = MEDIUM_FONT).grid(row=1,column=0)
                self.Exponential_LowLimit = tk.Entry(self.ExponentialFitParameterFrame, width=7)
                self.Exponential_LowLimit.grid(row=2,column=0)
                self.Exponential_LowLimit.insert(END, 4)

                self.Exponential_HighLabel = tk.Label(self.ExponentialFitParameterFrame, text = 'Higher Limit (ms)', font = MEDIUM_FONT).grid(row=1,column=1)
                self.Exponential_HighLimit = tk.Entry(self.ExponentialFitParameterFrame, width=7)
                self.Exponential_HighLimit.grid(row=2,column=1)
                self.Exponential_HighLimit.insert(END, 110)

                self.ApplyExponentialParameters = tk.Button(self.ExponentialFitParameterFrame, text = 'Apply Adjustments',command = lambda: self.ApplyExponentialAdjustments())
                self.ApplyExponentialParameters.grid(row=3,column=0,columnspan=2,pady=5)

                ###############################
                ### Start and Reset Buttons ###
                ###############################
                self.StartButton = tk.Button(self, text = 'Begin Visualization', command = lambda: self.StartProgram()).grid(row=4,column=0,pady=5)
                self.ResetButton = tk.Button(self, text = 'Reset', command = lambda: self.Reset()).grid(row=4,column=1,pady=5)


            #--- Real-Time Adjustment of visualization parameters ---#
            def AdjustParameters(self):
                #--- Adjusts the parameters used to visualize the raw voltammogram, smoothed currents, and polynomial fit
                global short_aa, long_aa, short_zz, long_zz

                #--- parameters for short pulse ---#
                short_aa = float(self.short_aa_entry.get())               # aa/zz adjust the points at the start and end of the polynomial fit, respectively
                short_zz = float(self.short_zz_entry.get())

                #--- parameters for long pulse ---#
                long_aa = float(self.long_aa_entry.get())
                long_zz = float(self.long_zz_entry.get())

            def SwitchParameterFrame(self, frame):

                if frame == 'Short':
                    self.ShortParameterFrame.tkraise()

                    self.SelectShortParameters['style'] = 'On.TButton'
                    self.SelectLongParameters['style'] = 'Off.TButton'

                if frame == 'Long':
                    self.LongParameterFrame.tkraise()

                    self.SelectShortParameters['style'] = 'Off.TButton'
                    self.SelectLongParameters['style'] = 'On.TButton'

            def ApplyExponentialAdjustments(self):
                global exp_high, exp_low

                exp_low = (int(self.Exponential_LowLimit.get()))/1000       # msec --> sec
                exp_high = (int(self.Exponential_HighLimit.get()))/1000

            def StartProgram(self):
                global chronoamperometry_analysis, AlreadyInitiated, exp_low, exp_high

                AlreadyInitiated = True

                exp_low = int(self.Exponential_LowLimit.get())/1000
                exp_high = int(self.Exponential_HighLimit.get())/1000

                #--- begin the animation of the chronoamperometry data ---#
                chronoamperometry_analysis = ChronoamperometryAnalysis()


            #--- Function to Reset and raise the user input frame ---#
            def Reset(self):
                global PoisonPill, ThreadQueue

                if AlreadyInitiated:
                    PoisonPill.set()

                    chronoamperometry_analysis.Reset()

                self.close_frame(RealTimeInputFrame)
                self.show_frame(StartPage)

            def show_frame(self, cont):

                frame = ShowFrames[cont]
                frame.tkraise()

            def close_frame(self, cont):

                frame = ShowFrames[cont]
                frame.grid_forget()

                PlotContainer.destroy()


        #--- Initiate the Real-Time User Input Frame ---#
        RT_InputFrame = RealTimeInputFrame(container, self)
        RT_InputFrame.grid(row=0,column=0,sticky='nsw')
        ShowFrames[RealTimeInputFrame] = RT_InputFrame
        RT_InputFrame.tkraise()



########################################################
### Create the figures and axes for the loglog plots ###
########################################################
class InitiateFigures():
    def __init__(self):
        global figures, container, PlotContainer

        figures = {}    # global figure and axes list

        #--- create a figure and axes that will render the data ---#
        self.MakeFigure()

        class PlotContainer(tk.Frame):
            def __init__(self, parent, controller):

                tk.Frame.__init__(self, parent)             # initialize the frame

        #--- Create a container that can be created and destroyed when Start() or Reset() is called, respectively ---#
        PlotContainer = PlotContainer(container, self)
        PlotContainer.relief = 'flat'
        PlotContainer.grid(row=0,column=1)

        ###############################################################
        ### Class for creating instances of the visualization frame ###
        ###############################################################
        class VisualizationFrame(tk.Frame):

            def __init__(self, parent, controller):
                tk.Frame.__init__(self, parent)         # Initialize the frame

                #--- Voltammogram, Raw Peak Height, and Normalized Figure and Artists ---#
                fig = figures['figure']
                canvas = FigureCanvasTkAgg(fig, self)                                         # and place the artists within the frame
                canvas.draw()                                                                 # initial draw call to create the artists that will be blitted
                canvas.get_tk_widget().grid(row=2,padx=5,pady=6,ipady=5)         # does not affect size of figure within plot container


        FrameReference = VisualizationFrame(PlotContainer, self)
        FrameReference.grid(row=0,column=0)
        FrameReference.tkraise()

    def MakeFigure(self):
        global figures

        fig, ax = plt.subplots(nrows=1,ncols=2)
        plt.subplots_adjust(bottom=0.1,hspace=0.6,wspace=0.3)         ### adjust the spacing between subplots

        #-- set the limits of the loglog axes ---#
        ax[0].set_ylim(5*10**-9,10**-4)           # Amperes
        ax[0].set_xlim(4*10**-4,5*10**-1)         # seconds

        ax[0].set_ylabel('Current (A)',fontweight='bold')
        ax[0].set_xlabel('Time (s)',fontweight='bold')

        #-- set the chronoamperogram to a loglog scale
        ax[0].set_yscale('log')
        ax[0].set_xscale('log')

        #--- create a plot for Fitted Half life v File Number
        ax[1].set_ylim(0,15)                                    # Half Life
        ax[1].set_xlim(0,int(file_limit_entry.get())+0.1)      # File Number

        ax[1].set_ylabel('Half Life (ms)',fontweight='bold')
        ax[1].set_xlabel('File Number',fontweight='bold')

        figures['figure'] = fig
        figures['axes'] = ax

        print('Made Figures')





#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#


    ##########################################################################
    ##########################################################################
    ###   ANIMATION FUNCTION TO HANDLE ALL DATA ANALYSIS AND VISUALIZATION ###
    ##########################################################################
    ##########################################################################


class ElectrochemicalAnimation():
    def __init__(self, fig, ax, generator = None, func = None, resize_interval = None, fargs = None):
        global PoisonPill

        print('ASDFASF')

        self.file = 1                                            # Starting File
        self.index = 1                                           # File Index Value
        self.ax = ax                                             # Figure Axes object

        ### Lists for sample rate (time passed)  ###
        self.sample_list = []
        self.file_list = []

        #-- set the poison pill event for Reset --#
        self.PoisonPill = Event()
        PoisonPill = self.PoisonPill             # global reference


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

        print('START')
        # First disconnect our draw event handler
        self._fig.canvas.mpl_disconnect(self._first_draw_id)
        self._first_draw_id = None  # So we can check on save

        # Now do any initial draw
        self._init_draw()

        ### Create a thread to analyze obtain the file from a Queue
        ### and analyze the data.

        class _threaded_animation(threading.Thread):

            def __init__(self, Queue):
                print('THREADED INIT')
                #global PoisonPill

                threading.Thread.__init__(self)     # initiate the thread

                self.q = Queue

                #-- set the poison pill event for Reset --#
                self.PoisonPill = Event()
                PoisonPill = self.PoisonPill             # global reference

                self.file = 1

                root.after(10,self.start)                       # initiate the run() method

            def run(self):

                print('RUNNING')
                while True:
                    try:
                        print('Getting Task')
                        task = self.q.get(block=False)
                        print('Got Task')

                    except:
                        break
                    else:
                        if not PoisonPill:
                            print('Perorming Task')
                            root.after(Interval,task)

                if not PoisonPill:
                    print('Running Again')
                    root.after(10, self.run)

        self.q = Queue()
        threaded_animation = _threaded_animation(Queue = self.q)

        self._step()


    def _stop(self, *args):
        # On stop we disconnect all of our events.
        self._fig.canvas.mpl_disconnect(self._resize_id)
        self._fig.canvas.mpl_disconnect(self._close_id)

    def _setup_blit(self):
        # Setting up the blit requires: a cache of the background for the
        # axes
        print('BLIT')
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

        self._drawn_artists = plots['EmptyPlots']
        print(plots['EmptyPlots'])

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
        self._post_draw(False)
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

        self._blit_draw(self._drawn_artists, self._blit_cache)


    # The rest of the code in this class is to facilitate easy blitting
    def _blit_draw(self, artists, bg_cache):
        # Handles blitted drawing, which renders only the artists given instead
        # of the entire figure.
        updated_ax = []
        print('BLIT DRAW')
        try:
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
        except:
            print('BLIT FAIL')


    ## callback that is called every 'interval' ms ##
    def _step(self):

        print('STEP')

        if self.file not in self.file_list:
            self.file_list.append(self.file)

        #--- get the file handle ---#
        short = 'Short_%s_#%d_#%d.DTA' % (handle_variable,self.file,self.index)
        long = 'Long_%s_#%d_#%d.DTA' % (handle_variable,self.file,self.index)

        short_file = mypath + short
        long_file = mypath + long

        print('SHORT %s' % short_file)
        print('LONG %s' % long_file)

        try:
            mydata_bytes_short = os.path.getsize(short_file)    ### retrieves the size of the file in bytes
            mydata_bytes_long = os.path.getsize(long_file)    ### retrieves the size of the file in bytes
            print('Found FILES')
            mydata_bytes = True

        except:
            mydata_bytes = False


        #################################################################
        #### If the file meets the size requirement, analyze the data ###
        #################################################################
        if mydata_bytes:
            print('Putting task')
            self.q.put(self.run_analysis())


        else:
            if not PoisonPill:
                root.after(100,self._step)

    def run_analysis(self):

        print('Running Analysis')
        try:
            framedata = next(self.generator(file=self.file))
            self._draw_next_frame(framedata)
            self.file += 1
            root.after(10,self._step)

        except:
            print('ERROR IN RUN ANALYSIS')

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

        if not PoisonPill:
            root.after(5, self._check_queue)


#############################################
### Data Analysis and Visualization class ###
#############################################
class ChronoamperometryAnalysis():
    def __init__(self):
        global anim

        self.initialize()

        #########################################
        ### Create a thread for Data Analysis ###
        #########################################

        anim = []
        self.anim = ElectrochemicalAnimation(figures['figure'],figures['axes'], func=self.animate, generator=self.frames)

        anim.append(self.anim)

    ####################################################################
    ### Analze the first file to retrieve the experiments parameters ###
    ### and create a dictionary containing their values              ###
    ####################################################################
    def initialize(self):
        global anim, ParameterPipe

        self.file = 1
        self.index = 1

        self.exp_low = exp_low
        self.exp_high = exp_high

        self.ParameterPipe = Queue()
        ParameterPipe = self.ParameterPipe

        self.volt_dict = {}     # dictionary for voltage parameters

        self.short_dict = self.get_parameters('short')
        self.long_dict = self.get_parameters('long')

        self.file_limit = int(file_limit_entry.get())
        self.point_limit = int(point_limit_entry.get())

        ### create a list to hold lifetime values
        self.lifetime_list = []
        self.fitted_lifetime_list = []

        ### prepare figures and artists to be visualized
        self.prepare_figure()

    def Reset(self):

        self.file = 1

    ######################################
    ### Create artists for each figure ###
    ### 1. raw exponential decay       ###
    ### 2. fitted exponential decay    ###
    ### 3. extrapolated lifetime (ms)  ###
    ######################################
    def prepare_figure(self):
        global plots

        #--- create a primtive artist that will contain the visualized data ---#
        ax = figures['axes']

        # Data Artists
        self.raw_decay, = ax[0].plot(1,1, 'ko', MarkerSize=1)        # loglog amperogram
        self.curve_fit, = ax[0].plot(1,1, 'r-', MarkerSize=1)   # loglog curve fit
        self.fitted_lifetime, = ax[1].plot(0,0, ' bo', MarkerSize=1) # fitted lifetime (ms)

        # empty arists for empty init functions
        self.EmptyPlots, = ax[0].plot(1,1, 'ro', MarkerSize=1)

        plots = {}
        plots['line'] = self.raw_decay        # global reference
        plots['curve fit'] = self.curve_fit
        plots['fitted lifetime'] = self.fitted_lifetime
        plots['EmptyPlots'] = [self.EmptyPlots]

    ### Analyze the first file and extrapolate parameters ###
    ### 1. voltage_1: 1st potential step (int; V)
    ### 2. voltage_2: 2nd potential step (int; V)
    ### 3. pulse_width: width of voltage_1 (int; V)
    ### 4. sample_rate: sample rate (int; s)
    ### 5. variables: # of points per pulse (int)
    def get_parameters(self, pulse):

        # dicionary for voltage parameters #

        if pulse == 'short':
            filehandle = 'Short_%s_#1_#1.DTA' % handle_variable
        elif pulse == 'long':
            filehandle = 'Long_%s_#1_#1.DTA' % handle_variable

        myfile = mypath + filehandle

        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes
        except:
            mydata_bytes = 1

        #---if the file meets the size requirement, analyze data---#
        if mydata_bytes > 1000:

            #---Preallocate Potential and Current lists---#
            with open(myfile,'r',encoding='utf-8') as mydata:
                voltages = []
                dictionary = {}
                val = 0
                for line in mydata:
                    line = ' '.join(line.split())
                    line = line.strip()

                    #--- Step One Potential ---#
                    if line.startswith('VSTEP1'):
                        voltage_1 = ' '.join(line.split())          # remove repeated white space
                        voltage_1 = voltage_1.split(' ')            # then separate each string
                        voltage_1 = voltage_1[2]                    # and grab the potential from the 3rd column
                        volt_1 = voltage_1
                        voltage_1 = voltage_1.split('E')
                        exponent = voltage_1[1]
                        voltage_1 = voltage_1[0]
                        voltage_1 = float(voltage_1) * (10**float(exponent))
                        self.volt_dict[voltage_1] = volt_1

                    #--- Step Two Potential ---#
                    if line.startswith('VSTEP2'):
                        voltage_2 = ' '.join(line.split())
                        voltage_2 = voltage_2.split(' ')
                        voltage_2 = voltage_2[2]
                        volt_2 = voltage_2
                        voltage_2 = voltage_2.split('E')
                        exponent = voltage_2[1]
                        voltage_2 = voltage_2[0]
                        voltage_2 = float(voltage_2) * (10**float(exponent))
                        self.volt_dict[voltage_2] = volt_2

                    #--- Step One Pulse Width ---#
                    if line.startswith('TSTEP1'):
                        pulse_width = ' '.join(line.split())
                        pulse_width = pulse_width.split(' ')
                        pulse_width = pulse_width[2]
                        pulse_width = pulse_width.split('E')
                        exponent = pulse_width[1]
                        pulse_width = pulse_width[0]
                        pulse_width = float(pulse_width) * (10**float(exponent))
                        dictionary['pulse width'] = pulse_width

                    #--- Sample Rate ---#
                    if line.startswith('SAMPLETIME'):
                        sample_rate = ' '.join(line.split())
                        sample_rate = sample_rate.split(' ')
                        sample_rate = sample_rate[2]

                        sample_rate = sample_rate.split('E')
                        exponent = sample_rate[1]
                        sample_rate = sample_rate[0]
                        sample_rate = float(sample_rate) * (10**float(exponent))
                        dictionary['sample period'] = sample_rate

                    val += 1

            #--- extrapolate the number of files in one pulse ---#
            variables = math.ceil(int((pulse_width/sample_rate)))
            dictionary['variables'] = variables

            #--- extrapolate the reduction potential --#
            reduction_voltage = max(abs(voltage_1),abs(voltage_2))            # reduction voltage
            if voltage_1 < 0:
                reduction_voltage = reduction_voltage * -1
            dictionary['reduction voltage'] = reduction_voltage

            return dictionary

        ### If the file has not been found, keep searching
        else:
            print('could not initialize')
            root.after(200, self.initialize)

    ####################################################################
    ### Using the data provided, calculate the lifetime of the decay ###
    ####################################################################
    def set_fit(self, time_data, current_data):

        A1, C1, K1, A2, C2, K2, adjusted_time, adjusted_currents = self.extract_parameters(time_data, current_data)

        print(A1, C1, K1, A2, C2, K2)
        print('Initial Constant value:',C2)

        popt, pcov = curve_fit(self.biexponential_func, adjusted_time, adjusted_currents, p0 = (A1, K1, A2, K2))

        print(popt)
        print(pcov)

        return popt, pcov, C2


    def extract_parameters(self, time_data, current_data):

        zip_list = sorted(set(zip(time_data, current_data)))

        # create a key,value dictionary with format {current: time}
        initial_dict = {}
        for item in zip_list:
            time, current = item
            initial_dict[time] = current

        ## Adjust the time data to fit the parameters set by the user
        try:
            adjusted_time = [time for time in time_data if exp_low <= time <= exp_high]
        except:
            print('\nCould not apply exponential adjustments to time data\n')

        ## create a new dictionary that will correlate
        ## the currents to the adjusted time_data
        adjusted_currents = []
        adjusted_dict = {}
        for item in zip_list:
            time, current = item
            if time in adjusted_time:
                adjusted_currents.append(current)
                adjusted_dict[current] = time

        ###########################
        ### Monoexponential Fit ###
        ###########################

        ### get a initial assumption for the A value in y = Ae^(-kt) + C
        initial_A = adjusted_currents[5]

        ### get an initial assumption for the C value in y = Ae^(-kt) + C
        initial_C = adjusted_currents[-5]

        ### get the y value for the half-life
        half_life = 0.33*(float(initial_A))

        ### find the closest real value to the theoretical half life
        closest_current = min(adjusted_currents, key = lambda x:abs(x-half_life))
        closest_time = adjusted_dict[closest_current]

        ## linearize both sides  with the natural log and extrapolate K
        ## in exponential equation y = Ae^(-kt) + C
        initial_K = (math.log(initial_A) - math.log(closest_current - initial_C))/closest_time
        print('\nMonoexponential Fit: Time Constant', initial_K)

        ### Separate the data into two zip files to ###
        ### fit to a biexponential function         ###

        half = math.ceil(0.5 * len(adjusted_time))

        decay_1 = sorted(set(zip(adjusted_time[:half],adjusted_currents[:half])))
        decay_2 = sorted(set(zip(adjusted_time[half:],adjusted_currents[half:])))
        print('\nDecay 1')
        A1, C1, K1 = self.extract_fit(decay_1)
        print('\nDecay 2')
        A2, C2, K2 = self.extract_fit(decay_2)


        print('extract parameters: returning values')
        return A1, C1, K1, A2, C2, K2, adjusted_time, adjusted_currents
        ### Get initial parameters for each decay ###


        ## return a parameters fitted to the curve with a minimization function
        ## and based upon the suggestion of the variables p0

    ### Curve is a zip list of (time, current)
    def extract_fit(self, curve):

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
            print('Initial Amplitude:',A)

            ### get an initial assumption for the C value in y = Ae^(-kt) + C
            C = currents[-20]
            print('Initial Constant:',C)
            self.constant = C

            ### get the y value for the half-life
            half_life = 0.33*(float(A))
            print('Initial Half Life:',half_life)

            ### find the closest real value to the theoretical half life
            closest_current = min(currents, key = lambda x:abs(x-half_life))
            print('closest current:',closest_current)
            closest_time = dict[closest_current]
            print('Closest Time:',closest_time)

            ## linearize both sides  with the natural log and extrapolate K
            ## in exponential equation y = Ae^(-kt) + C
            try:
                K = (math.log(A) - math.log(closest_current - C))/closest_time
            except:
                print('Constant > Closest Current')
                K = (math.log(A) - math.log(closest_current))/closest_time

            print('\nExtract Fit: Time Constant', K)

            return A, C, K
        except:
            print('couldt perform extraction')

    #######################
    ### Monoexponential ###
    #######################
    def exponential_func(self, t, a, k, c):
        return (a * np.exp(-t * k) + c)

    #####################
    ### Biexponential ###
    #####################
    def biexponential_func(self, t, a, b, c, d, e = None):

        if e is None:
            e = self.constant

        return a * np.exp(-t * b) + c * np.exp(-t * d) + e

    ###############################################
    ### Generator Function for data acquisition ###
    ###############################################
    def frames(self, file = None):
        global file_label, plots

        if file is None:
            self.file = 1
        else:
            self.file = file

        print('File Recieved from Queue: %s' % str(self.file))

        #try:
        file_label['text'] = str(self.file)

        ## Yield a list containing data
        ## avergaed over 'point_limit' number of
        ## chronoamperomegrams
        averaged_data = self.data_analysis()

        yield averaged_data


    def data_analysis(self):
        try:
            self.total_list = []
            self.index = 1

            #--- start the thread to analyze the data ---#
            #thread = threading.Thread(target = self.background_analyze()).start().join()          - NOPE
            self.analyze()

            ##########################################################
            ### Average the data from the merged data dictionaries ###
            ##########################################################
            sums = Counter()
            counters = Counter()
            for itemset in self.total_list:
                sums.update(itemset)
                counters.update(itemset.keys())

            try:
                averaged_data = {x: float(sums[x])/counters[x] for x in sums.keys()}
            except:
                averaged_data = None

            # return the average list

            return averaged_data

        except:
            print('\n\nChronoamperometryAnalysis.data_analysis: Error\n\n')
            t.sleep(5)


    def analyze(self):

        ###########################################################
        ### Iterate through each file that will used to create  ###
        ### a more accurate dataset of averaged data            ###
        ###########################################################
        while self.index <= self.point_limit:
            merged_dictionary = next(self.analysis_generator())
            self.total_list.append(merged_dictionary)
            t.sleep(0.01)


    #--- Generator function to return data to analyze ---#
    def analysis_generator(self):
        #--- Generator function that will yield data from all the potential steps
        #--- that will be averaged per file (titration point)
        while True:
            fig = figures['figure']

            print('\n\n\n POINT %s \n\n\n ' %str(self.index))

            ##############################
            ### Get data from the file ###
            ##############################
            short_time, short_currents  = self.retrieve_data('short')
            long_time, long_currents = self.retrieve_data('long')

            short_dictionary = {}
            count = 0
            for value in short_time:
                short_dictionary[value] = short_currents[count]
                count += 1

            long_dictionary = {}
            count = 0
            for value in long_time:
                long_dictionary[value] = long_currents[count]
                count += 1

            ###########################################
            ### Merge the two dictionaries into one ###
            ###########################################
            merged_dictionary = self.merge_two_dicts(short_dictionary, long_dictionary)

            t.sleep(0.001)

            self.index += 1

            yield merged_dictionary


    def retrieve_data(self, pulse):

        # variable used to tell retrive_data that the first
        # file has been found and that it should continue to
        # iterate through pulse_width/sample_rate number of files,
        # also known as 'variables'
        self.variable = False

        #--- get the file handle ---#
        if pulse == 'short':
            dictionary = self.short_dict
            variables = dictionary['variables']
            filehandle = 'Short_%s_#%d_#%d.DTA' % (handle_variable,self.file,self.index)
            aa = short_aa
            zz = short_zz

        elif pulse == 'long':
            dictionary = self.long_dict
            variables = dictionary['variables']
            filehandle = 'Long_%s_#%d_#%d.DTA' % (handle_variable,self.file,self.index)
            aa = long_aa
            zz = long_zz

        time = [0]*variables
        currents = [0]*variables

        myfile = mypath + filehandle


        ################################################################
        ### Create a loop in which the thread will continue to       ###
        ### search for the file for  'search_lim' amount of seconds  ###
        ################################################################
        try:
            mydata_bytes = os.path.getsize(myfile)    ### retrieves the size of the file in bytes
        except:
            print('\n\nChronoamperometryAnalysis.analyze: Could not find %s\n\n' % myfile)
            mydata_bytes = 1


        #---if the file meets the size requirement, analyze data---#
        if mydata_bytes > 1000:

            try:
                with open(myfile,'r',encoding='utf-8') as mydata:
                    val = 0
                    for line in mydata:
                        line = ' '.join(line.split())
                        line = line.strip()

                        if self.variable:
                            #print(line)
                            line = line.split(' ')

                            #--- extract the time ---#
                            seconds = line[1]

                            try:
                                #-- if its an exceptionally small number, it will
                                #-- be split with an 'E' exponent
                                seconds = seconds.split('E')
                                exponent = float(seconds[1])
                                seconds = float(seconds[0])
                                seconds = seconds * (10**exponent)

                            except:
                                seconds = seconds[0]
                            time[val] = seconds

                            #--- extract the current ---#
                            current = line[3]
                            current = current.split('E')
                            exponent = float(current[1])
                            current = float(current[0])
                            current = (current) * (10**(exponent))
                            #print('current: %s' % str(current))
                            currents[val] = current

                            #print('val: %s' % str(val))

                            #--- if this is the last file for this pulse, break ---#
                            if val == variables - 1:
                                self.variable = False
                                break

                            val += 1

                        #--- Find the first line and qery
                        elif line.startswith('1 '):
                            reduction_voltage = dictionary['reduction voltage']
                            reduction_voltage = self.volt_dict[reduction_voltage]
                            if reduction_voltage in line:
                                print('\n\n')
                                line = line.split(' ')
                                #print(line)
                                #t.sleep(1)

                                #--- extract the time ---#
                                seconds = line[1]
                                try:
                                    seconds = seconds.split('E')
                                    exponent = float(seconds[1])
                                    seconds = float(seconds[0])
                                    seconds = seconds * (10**exponent)

                                except:
                                    seconds = seconds[0]
                                time[val] = seconds

                                #--- extract the current ---#
                                current = line[3]
                                current = current.split('E')
                                exponent = float(current[1])
                                current = float(current[0])
                                current = (current) * (10**(exponent))
                                #print('current: %s' % str(current))
                                currents[val] = current

                                #print('start val: %s' % str(val))
                                #t.sleep(1)

                                val += 1

                                self.variable = True

                    time, currents = self._apply_adjustment(time, currents, aa, zz)

                    return time, currents

            except:
                print('Error in retrieve data')

        else:
            root.after(200, self.retrieve_data)

    def _apply_adjustment(self, time_list, current_list, aa, zz):

        time_list = [float(time) for time in time_list]

        zip_list = set(zip(time_list, current_list))

        # create a key,value dictionary with format {time: current}
        dict = {}
        for item in zip_list:
            time, current = item
            dict[time] = current

        adjusted_time = [time for time in time_list if aa <= time <= zz ]
        adjusted_currents = [dict[time] for time in adjusted_time]

        return adjusted_time, adjusted_currents

    def animate(self, framedata):

        #try:
        chronoamp_data = framedata
        time_data, current_data = zip(*sorted(chronoamp_data.items()))
        time_data = [abs(float(item)) for item in time_data]
        current_data = [abs(item) for item in current_data]

        # use a curve fitting function to extrpolate the
        # paramers of the exponential decay

        #####################################################################
        # popt: array = optimal parameters returned from least squares fit  #
        # pcov: 2D array = estimated covariance of popt                     #
        #####################################################################
        popt, pcov, C2 = self.set_fit(time_data, current_data)
        print(popt)

        # Get the time constant, lambda, with (1/k) from y = A * e^(-kt) + C
        fitted_lifetime = (1/popt[1])*1000                  # popt[1] = k
        print('animate..\n k=',popt[1])
        print('lifetime = ',fitted_lifetime)
        self.fitted_lifetime_list.append(fitted_lifetime)

        # set the chronoamperogram
        self.raw_decay.set_data(time_data,current_data)     # raw curve


        # use the curve fitting function to approximate values
        # fot A,k, and C in y = A * e^(-kt) + C
        fitted_data = []
        for x in time_data:
            fit = self.biexponential_func(x, *popt, e = C2)     # use the curve fit parameters
            fitted_data.append(fit)                     # for each point in the chronoamperogram

        self.curve_fit.set_data(time_data,fitted_data)  # set the curve fit as an artist

        # set the half life plots
        file_list = range(1,self.file+1)
        self.fitted_lifetime.set_data(file_list,self.fitted_lifetime_list) # checking curve fit


        #export any relevant dat
        if SaveVar:
            export_data._export_data([self.file,fitted_lifetime])

        ### return artists to be animated
        return self.raw_decay, self.curve_fit, self.fitted_lifetime


    def merge_two_dicts(self, x, y):            # any keys that overlap will be overridden with the items from y
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z



##################################
### Real-Time Text File Export ###
##################################
class TextFileExport():

    ###############################
    ### Initialize the .txt file ###
    ###############################
    def __init__(self):
        self.TextFileHandle = FileHandle

        try:

            list = ['File','Lifetime (ms)']

            #--- Write the data into the .txt file ---#
            with open(FileHandle,'w+',encoding='utf-8', newline = '') as input:
                writer = csv.writer(input, delimiter = ' ')
                writer.writerow(list)

        except:
            print('\n\n','ERROR IN TEXT FILE EXPORT','\n\n')
            time.sleep(0.5)

    def _export_data(self, data):

        #--- Write the data into the .txt file ---#
        with open(FileHandle,'a',encoding='utf-8', newline = '') as input:
            writer = csv.writer(input, delimiter = ' ')
            writer.writerow(data)
        with open(FileHandle,'r',encoding='utf-8', newline = '') as filecontents:
            filedata =  filecontents.read()
        filedata = filedata.replace('[','')
        filedata = filedata.replace('"','')
        filedata = filedata.replace(']','')
        filedata = filedata.replace(',','')
        filedata = filedata.replace('\'','')
        with open(FileHandle,'w',encoding='utf-8', newline = '') as output:
            output.write(filedata)




#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

                    ############################################
                    ### Initialize GUI to start the program  ###
                    ############################################

if __name__ == '__main__':

    root = RealTimeAnalysis()

    #-- Button Styling --#
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
