'''
Similar to IDLs lp_ximovie.pro but in python, 
for quick inspection of t-series.

Just provide a cube with dimensions (nt, ny, nx)
J. de la Cruz Rodriguez (ISP-SU 2019)

'''
import matplotlib
#matplotlib.use("TkAgg")

import sys
if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


class Gui:
    def __init__(self, root, input_data, figsize=8, dpi=100, cmap='gist_gray'):

        self.cmap = cmap
        self.nt, self.ny, self.nx = input_data.shape[:]
        self.d = input_data
        self.aspect = float(self.ny)/self.nx
        self.fsize = (figsize+1.5, self.aspect*figsize)
        print(self.aspect)
        self.autoplay = False
        
        # get min, max
        self.mi = self.d.min()
        self.ma = self.d.max()
        self.ra = self.ma-self.mi
        
        self.firsttime = True
        self.tstep = Tk.IntVar()
        self.tstep.set(0)

        self.Mival = Tk.DoubleVar()
        self.Mival.set(self.mi)
        self.Maval = Tk.DoubleVar()
        self.Maval.set(self.ma)
        
        # Copy instance of Tk
        self.master = root
        self.master.wm_title("Tseries animation")
        
        
        # Create frame
        self.frame = Tk.Frame(self.master)
        self.fig = Figure(figsize=self.fsize, dpi=dpi)
        
        # Figure
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0,column=1, rowspan=35)

        # time slider
        rowcounter = 0
        slider_length = figsize*dpi*0.2

        T_label = Tk.Label(self.frame, text='Time step:').grid(row=rowcounter, column=0, rowspan=1)
        rowcounter+=1
        self.T_entry = Tk.Scale(self.frame,
                                orient='horizontal',
                                length=slider_length,
                                from_=0, 
                                to=self.nt-1,
                                resolution=1,
                                variable=self.tstep,
                                command=self.redraw_from_event)
        
        self.T_entry.grid(row=rowcounter, column=0, columnspan=1)
        rowcounter +=1

        # autoplay
        play_button = Tk.Button(self.frame,text='Play',command=self.handle_play_loop)
        play_button.grid(row=rowcounter, column=0, columnspan=1,sticky='ew')
        rowcounter += 1
        

        # scale min
        Mi_label = Tk.Label(self.frame, text='Image min:').grid(row=rowcounter, column=0, rowspan=1)
        rowcounter += 1
        self.S_min = Tk.Scale(self.frame,
                              orient='horizontal',
                              length=slider_length,
                              from_=self.mi, 
                              to=self.mi+self.ra*0.75,
                              resolution=self.ra*0.02,
                              variable=self.Mival,
                              command=self.updateScaling)
        

        self.S_min.grid(row=rowcounter, column=0, columnspan=1)
        rowcounter+=1
        # scale max
        Ma_label = Tk.Label(self.frame, text='Image max:').grid(row=rowcounter, column=0, rowspan=1)
        rowcounter += 1
        self.S_max = Tk.Scale(self.frame,
                              orient='horizontal',
                              length=slider_length,
                              from_=self.ma-self.ra*0.75, 
                              to=self.ma,
                              resolution=self.ra*0.02,
                              variable=self.Maval,
                              command=self.updateScaling)
        

        self.S_max.grid(row=rowcounter, column=0, columnspan=1)
        rowcounter+=1
        
        # quit button
        quit_button = Tk.Button(self.frame,text='Quit',command=self.quit)
        quit_button.grid(row=rowcounter, column=0, columnspan=1,sticky='ew')
        
        
        # pack frame in main window
        self.frame.pack()
        
        #show default graph
        self.redraw()

    def updateScaling(self, event):
        self.image.set_clim(self.Mival.get(), self.Maval.get())
        self.canvas.draw()
        self.canvas.flush_events()
        
    def redraw_from_event(self, event):
        self.redraw()
            
    def redraw(self):
        tt = self.tstep.get()
        
        if(self.firsttime):
            self.fig.clear()
            ax = self.fig.add_subplot(111)
            self.image = ax.imshow(self.d[tt], cmap=self.cmap, vmax=self.ma, vmin=self.mi,
                                   interpolation='nearest')
            self.firsttime = False
            self.fig.set_tight_layout(True)
        else:
            self.image.set_data(self.d[tt])
            
        self.canvas.draw()
        self.canvas.flush_events()
        #self.image.figure.canvas.draw()
        
    def handle_play_loop(self):
        if(self.autoplay):
            self.autoplay = False
        else:
            self.autoplay = True
            self.play_loop()
    def play_loop(self):
        while(self.autoplay):
            tt =  self.tstep.get() + 1
            if(tt == self.nt): tt = 0
            self.tstep.set(tt)
            self.T_entry.set(tt)
            self.redraw()
            
    def quit(self):
        self.autoplay = False
        print("exiting")
        self.master.destroy()
        
        
def tseries(input_data, figsize=8, dpi=80, cmap='gist_gray'):
    """
    Similar to IDLs lp_ximovie.pro but in python, 
    for quick inspection of t-series.
    
    Just provide a cube with dimensions (nt, ny, nx)
    J. de la Cruz Rodriguez (ISP-SU 2019)
    """
    
    root = Tk.Tk()
    app = Gui(root, input_data, figsize=figsize, dpi=dpi, cmap=cmap)
    root.mainloop()

    del app
    del root
    
