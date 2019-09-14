"""Z-score_tkinter.py: Written by Phil Wilmarth, OHSU.
Applies a Z-transformation to pair-wise quantitative proteomics data.
Requires Python 3.x, numpy, and pandas.

written by Phil Wilmarth 2011-2016, Oregon Health & Science University.

MIT License

Copyright (c) 2019 OHSU and Phillip Wilmarth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.THE SOFTWARE WILL NOT INFRINGE
ANY PATENT, TRADEMARK OR OTHER RIGHTS.
"""
# Python Tkinter program for calculating Z-score transformations
# written by Phil Wilmarth, OHSU, 2012, 2016
# updated for Python3 and uses pandas
#
# standard libraries
import tkinter as tk
from tkinter import ttk
import os
import sys

# scientific stack libraries
import numpy as np
import pandas as pd
import scipy.stats
from scipy.optimize import curve_fit


# create some globals
WINDOW = 301    # sliding window width (use 41, 61, 81, or 101)
TRIM_PC = 5.0 # trim value, in percent (upper AND lower X% data points trimmed)
ZERO_CORR = 50.0 # zero correction to avoid math errors (0.15 for spectral counts)

# status bar class

class StatusBar(tk.Frame):
    """window status bar class with set and clear methods
    borrowed code from "http://effbot.org/tkinterbook/text.htm"
    """
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.label = tk.Label(self, bd=1, padx=5, pady=1, relief=tk.SUNKEN, anchor=tk.W)
        self.label.pack(fill=tk.X)
                
    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()
                
    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()
    
    # end class

class ZScoreGUI:
    """GUI for calculating Z-score transformations of quantitative data.
    """
    def __init__(self):
        # create root window
        self.root = tk.Tk()
        self.root.title("Z-score Calculator")

        # create a toolbar with buttons
        self.toolbar = tk.Frame(self.root)
        tk.Button(self.toolbar, text='Get Data', width=8, command=self.get_data, borderwidth=2,
                  relief=tk.RAISED).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(self.toolbar, text='Compute', width=8, command=self.compute, borderwidth=2,
                  relief=tk.RAISED).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(self.toolbar, text='Clear', width=8, command=self.clear_data, borderwidth=2,
                  relief=tk.RAISED).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(self.toolbar, text='Help', width=8, command=self.print_help, borderwidth=2,
                  relief=tk.RAISED).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(self.toolbar, text='Quit', width=8, command=self.quit_me, borderwidth=2,
                  relief=tk.RAISED).pack(side=tk.LEFT, padx=5, pady=5)
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        
        self.create_defaults_frame()

        # add a text box
        self.text_frame = tk.Frame(self.root)
        self.text_frame.pack()
        self.text = tk.Text(self.text_frame, wrap=tk.NONE, height=40, width=132, padx=5, pady=5)
        self.vscroll = tk.Scrollbar(self.text_frame, command=self.text.yview)
        self.text.configure(yscrollcommand=self.vscroll.set)
        self.hscroll = tk.Scrollbar(self.text_frame, orien=tk.HORIZONTAL, command=self.text.xview)
        self.text.configure(xscrollcommand=self.hscroll.set)
        self.hscroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.vscroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.text.pack(side=tk.LEFT)

        # add a status line
        self.status = StatusBar(self.root)
        self.status.pack(side=tk.BOTTOM, fill=tk.X)
        self.status.set("%s", "Status line")

        # data attributes
        self.data_frame = None
        self.print_help()

        # maybe this gets the window to the top?
        self.root.lift() 
        self.root.attributes('-topmost', True)
        self.root.attributes('-topmost', False)
        self.root.focus_force()
        
        # enter event loop
        self.root.mainloop()
        
    def create_defaults_frame(self):
        """Lets the user change minimum on count, intensity, etc.
        """
        defaults_frame = ttk.Labelframe(self.root, text='Default parameters:')
        defaults_frame.pack(fill=tk.X, expand=tk.YES, padx=5, pady=5)
        top_defaults = ttk.Frame(defaults_frame)
        top_defaults.pack(fill=tk.X, expand=tk.YES)
        bottom_defaults = ttk.Frame(defaults_frame)
        bottom_defaults.pack(fill=tk.X, expand=tk.YES)
        
        #Variables
        self.window = tk.IntVar()
        self.trim_pc = tk.DoubleVar()
        self.zero_corr = tk.DoubleVar()
        self.low = tk.DoubleVar()
        self.med = tk.DoubleVar()
        self.high = tk.DoubleVar()
        
        #Creation
        self.create_entry(top_defaults, 'Sliding window width (odd #): ', self.window).pack(side=tk.LEFT, padx=5, pady=5)
        self.create_entry(top_defaults, 'Trim %: ', self.trim_pc).pack(side=tk.LEFT, padx=5, pady=5)
        self.create_entry(top_defaults, 'Missing data input: ', self.zero_corr).pack(side=tk.LEFT, padx=5, pady=5)
        self.create_entry(bottom_defaults, 'Low p-value: ', self.low).pack(side=tk.LEFT, padx=5, pady=5)
        self.create_entry(bottom_defaults, 'Medium p-value: ', self.med).pack(side=tk.LEFT, padx=5, pady=5)
        self.create_entry(bottom_defaults, 'High p-value: ', self.high).pack(side=tk.LEFT, padx=5, pady=5)
      
        #Set Defaults
        self.window.set(WINDOW)
        self.trim_pc.set(TRIM_PC)
        self.zero_corr.set(ZERO_CORR)
        self.low.set(0.10)
        self.med.set(0.05)
        self.high.set(0.01)
        return
        
    #Functions to help create widgets             
    def create_entry(self, root, label, variable):
        """Creates a text entry widget.
        """
        frame = tk.Frame(root)
        tk.Label(frame, text=label).pack(side=tk.LEFT)
        tk.Entry(frame, textvariable=variable, width=10).pack(side=tk.LEFT)
        return frame 
        
    # toolbar button functions
    def get_data(self):
        """Gets numerical data from the clipboard.
        """
        # get data from the clipboard and show it in the window
        self.data_frame = pd.read_clipboard(thousands=',')
        if len(self.data_frame.columns) != 2:
            self.clear_screen()
            self.text.insert("1.0", "WARNING: data should be 2 columns!")
            return
        self.data_frame[self.data_frame == 0.0] = self.zero_corr.get()
        self.clear_screen()
        self.print_frame()
        self.status.set("%s", "%s data points read from clipboard" % len(self.data_frame))
           
    def compute(self):
        """Computes Ave SpC, log2 ratios, and Z-scores; results to window and clipboard
        """
        self.clear_data()
        if len(self.data_frame) == 0:
           self.text.insert("1.0", 'WARNING: data needs to be loaded first!')
           return
        
        # add computed columns
        cols = self.data_frame.columns.values
        self.data_frame['Original'] = np.arange(1, len(self.data_frame)+1) # original order
        temp = self.data_frame[[cols[0], cols[1]]]
        self.data_frame['AveAB'] = temp.mean(axis=1) # average of A and B
        self.data_frame['Log2(B/A)'] = np.log2(self.data_frame[cols[1]]/self.data_frame[cols[0]]) # Log2 of B/A ratio
        self.data_frame['FC'] = self.add_FC(self.data_frame[cols[0]], self.data_frame[cols[1]]) # fold-change (B vs A)
        self.data_frame.sort_values(by='AveAB', ascending=False, inplace=True)     # sort descending by average
        self.data_frame.reset_index(drop=True, inplace=True) # pandas does row indexing on index and (mostly) by label
        
        zs = []
        for i, row in enumerate(self.data_frame['Log2(B/A)']):
            zs.append(self.sliding_zscore(self.data_frame['Log2(B/A)'], i, self.window.get())) # compute sliding-window Z-score
        self.data_frame['Z-Score'] = zs
        
        # put back in original order
        self.data_frame.sort_values(by='Original', inplace=True)
        self.data_frame.reset_index(drop=True, inplace=True) # restore index back to original order
        self.data_frame.drop('Original', axis=1, inplace=True)
        
        # histogram, fit Gaussian, and compute p-values
        self.fit_Gaussian()
        self.p_values()
        self.BH_correction()
        self.data_frame['candidate'] = self.data_frame['FDR'].map(self.set_candidates)
        
        # print table to console and to clipboard
        self.root.clipboard_clear()
        self.print_frame()
        self.status.set("%s", "computed %s z-scores" % len(self.data_frame))

    def clear_screen(self):
        """Clears the window
        """
        self.text.delete("1.0", tk.END)
            
    def clear_data(self):
        """Clears the window contents and clipboard
        """
        self.text.delete("1.0", tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append("")
        self.status.set("%s", "Data cleared")

    def print_help(self):
        """Prints some help info to the window
        """
        self.clear_screen()
        help_text = \
"""Computes "Ave Value", "Log2 Ratio(A/B)", "Fold-change", "Z-score",
"p-values", and "FDR" (Benjamini-Hochberg correction) of quant values
for A and B, where A and B are two conditions being compared. The average
value computed so that proteins can be sorted from highest to lowest
abundance. Log2 of the ratio of counts between condition A and condition B
may have abundance-dependent bias that the sliding-window Z-score is
designed to remove. The Z-transformation uses trimmed means and standard
deviations (the largest and smallest 5% of ratios are trimmed in the protein
sliding window). The input data is two quantitative columns for A and B from 
Excel that have been pasted onto the clipboard. The data can have optional 
header text or not.
 
The two columns of quantitative data should be selected in Excel and 
copied to the clipboard. Once the data has been copied, click
the "Get Data" button. Data will be read and echoed in the text window.
The status line at the bottom will show the number of data points read.
Click the "Compute" button to have the computed quantities calculated 
and displayed. The computed values are also written to the clipboard for 
pasting back into Excel. "Clear" clears the clipboard contents, internal 
data structures and the screen. Note: "Clear" does not need to be pressed 
to process more data. Just overwrite the clipboard contents by pasting in 
new pairs of quantitative columns and click "Get Data" again to update internal 
data structures. "Help" prints this text. "Quit" ends the program and closes
the window.

Written by Phil Wilmarth, OHSU, 2012-6."""
        self.text.insert("1.0", help_text)

    def quit_me(self):
        """Quits the GUI.
        """
        self.status.set("%s", "Bye")
        self.root.withdraw()
        self.root.update_idletasks()
        self.root.quit()
        return
            
    def print_frame(self):
        """prints data to window and clipboard
        """
        self.text.insert("1.0", self.data_frame.to_string(index=False))
        self.data_frame.to_clipboard(index=False)

    def add_FC(self, A, B):
        """Adds a traditional fold-change column.
        """
        fc = []
        for a, b in zip(A, B):
            if b > a:
                FC = b/float(a)
            else:
                FC = -1 * a/float(b)
            fc.append(FC)
        return fc

    def sliding_zscore(self, vector, index, window=101):
        """Computes the trimmed Z-score in a symmetric window.
        written by Phil Wilmarth, OHSU, 2012.
        """
        # vector needs to be longer than window width
        if len(vector) < window:
            print('...WARNING vector is shorter than sliding window.')
            return None
        
        # test that index is in valid range
        if index < 0 or index > len(vector):
            print('...WARNING index out of range.')
            return None
        
        # test for beginning, middle, or end of vector and set window
        trim = int((window-1) * (self.trim_pc.get()/100.0))
        if index < window//2:
            vec_win = vector[:window]
        elif index > (len(vector) - (window//2)):
            vec_win = vector[-window:]
        else:
            vec_win = vector[index-window//2:index+window//2+1]
    
        # compute and return Z-score
        vec_win.sort_values(ascending=False, inplace=True)
        vec_trim = vec_win[trim:window-trim]
        mean = vec_trim.mean()
        stdev = vec_trim.std()
        try:
            z_score = (vector.iloc[index] - mean)/stdev
        except ZeroDivisionError:
            z_score = np.nan    # this was set to equal a space
        return z_score
        
    def Gaussian(self, x, amp, mean, sigma):
        """Basic Gaussian function."""
        return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

    def fit_Gaussian(self):
        """Hisotgrams Z-scores and fits a Gaussian to histogram."""        
        # histogram the Z-scores
        bins = np.linspace(-3.1, 3.1, 63)
        centers = bins + 0.05
        centers = centers[:-1]
        hist, bin_edges = np.histogram(self.data_frame['Z-Score'], bins)
        
        # fit the histogram with a Gaussian
        amp_guess = 0.01 * hist.sum() / np.sqrt(2 * np.pi)
        mean_guess = 0.0
        sigma_guess = 1.0
        guesses = [amp_guess, mean_guess, sigma_guess]
        params, covar = curve_fit(self.Gaussian, centers, hist, p0=guesses)
        self.mean = params[1]
        self.sigma = params[2]
        return

    def p_values(self):
        """Computes p-values of Z-scores."""
        # uses cumulative distribution function for 2-tailed probabilities
        Z = self.data_frame['Z-Score'].abs()
        cdf = scipy.stats.norm(self.mean, self.sigma).cdf(Z) 
        self.data_frame['p-value'] = 2 * (1.0 - cdf)
        return
        
    def BH_correction(self):
        """Computes a Benjamini-Hochberg multiple-testing correction."""
        p_frame = self.data_frame['p-value'].to_frame()
        p_frame['original'] = np.arange(len(p_frame))
        p_frame.sort_values(by='p-value', ascending=True, inplace=True)
        p_frame.reset_index(drop=True, inplace=True)
        bh_values = []
        total_tests = len(p_frame)
        prev_value = 0.0
        for i, p_value in enumerate(p_frame['p-value']):
            bh_value = (p_value * total_tests) / float(i + 1)
            bh_value = min(bh_value, 1.0)
            bh_value = max(bh_value, prev_value)
            prev_value = bh_value
            bh_values.append(bh_value)
            
        # save values and put back in original order
        p_frame['FDR'] = bh_values
        p_frame.sort_values(by='original', ascending=True, inplace=True)
        p_frame.reset_index(drop=True, inplace=True)
        self.data_frame['FDR'] = p_frame['FDR']
        return

    def set_candidates(self, fdr):
        """Label candidates according to ranges of p-values."""
        if fdr >= self.low.get():
            label = 'no'
        elif self.low.get() > fdr >= self.med.get():
            label = 'low'
        elif self.med.get() > fdr >= self.high.get():
            label = 'med'
        elif self.high.get() > fdr:
            label = 'high'
        else:
            label = None
        return label
            

# MAIN program starts here

if __name__ == '__main__':
    print(os.getcwd())
    zscore = ZScoreGUI()
    df = zscore.data_frame
    print('at end of program')

