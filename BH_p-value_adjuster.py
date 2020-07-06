"""BH_p-value_adjuster.py: Written by Phil Wilmarth, OHSU.
Applies a Benjamini-Hochberg correcion to a list of p-values.
adapted from:
http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
Copyright 2013, Oregon Health & Science University.
All Rights Reserved.

Permission to use, copy, modify, and distribute any part of this program
for non-profit scientific research or educational use, without fee, and
without a written agreement, is hereby granted, provided that the above
copyright notice, and this license agreement appear in all copies.
Inquiries regarding use of this software in commercial products or for
commercial purposes should be directed to:

Technology & Research Collaborations, Oregon Health & Science University,
2525 SW 1st Ave, Suite 120, Portland, OR 97210
Ph: 503-494-8200, FAX: 503-494-4729, Email: techmgmt@ohsu.edu.

IN NO EVENT SHALL OREGON HEALTH & SCIENCE UNIVERSITY BE LIABLE TO ANY
PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE.  THE
SOFTWARE IS PROVIDED "AS IS", AND OREGON HEALTH &SCIENCE UNIVERSITY HAS
NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, OR ENHANCEMENTS.
OREGON HEALTH & SCIENCE UNIVERSITY MAKES NO REPRESENTATIONS NOR EXTENDS
WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A
PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE
ANY PATENT, TRADEMARK OR OTHER RIGHTS.
"""
#
# Python Tkinter program for doing BH adjusted p-values
# written by Phil Wilmarth, OHSU, 2013
#
from tkinter import *
import os
import sys
import math
#
# create some globals
#
original = []
adjusted = []
#
# status bar class
#
class StatusBar(Frame):
    """window status bar class with set and clear methods
    borrowed code from "http://effbot.org/tkinterbook/text.htm"
    """
    def __init__(self, master):
        Frame.__init__(self, master)
        self.label = Label(self, bd=1, padx=5, pady=1, relief=SUNKEN, anchor=W)
        self.label.pack(fill=X)
    #            
    def set(self, format, *args):
        self.label.config(text=format % args)
        self.label.update_idletasks()
    #            
    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()
    #
    # end class
    #
#
# toolbar button functions
#
def get_data():
    """Gets numerical data from the clipboard
    """
    global original
    header = []
    #
    # get data from the clipboard and show it in the window
    #
    contents = root.clipboard_get()
    lines = contents.splitlines()
    original = []
    for i, line in enumerate(lines):
        line = line.replace(',', '').strip() # remove any commas
        temp = line.split('\t')
        if len(temp) != 1:
            clear_screen()
            text.insert("1.0", "WARNING: data should be a single column!")
            return
        try:
            original.append([i, float(temp[0])])
        except ValueError:
            header = temp
    clear_screen()
    print_data(header, original)
#    print 'length of xy data:', len(xy_data)
    status.set("%s", "%s lines read from clipboard" % len(lines))
#        
def compute():
    """Computes BH adjusted p-values
    """
    global original
    clear_data()
    if len(original) == 0:
       text.insert("1.0", 'WARNING: data needs to be loaded first!')
       return
    #
    # sort by increasing p-value and compute adjusted p-values
    #
    prev_bh_value = 0.0
    num_total_tests = float(len(original))
    original = [[x[-1]]+x for x in original]    # prefix p-values to lists
    original.sort()     # sort by increasing p-values
    for i, orig_list in enumerate(original):
        p_value = orig_list[-1]
        bh_value = (p_value * num_total_tests) / float(i+1)
        bh_value = min(bh_value, 1.0)
        bh_value = max(bh_value, prev_bh_value)     # makes BH values monotonic
        prev_bh_value = bh_value
        orig_list.append(bh_value)
##        if i < 10:
##            print 'p-val', p_value
##            print 'BH', (p_value * num_total_tests) / float(i+1)
##            print i, p_value, bh_value
    #
    # put back in original order
    #
    original = [x[1:] for x in original]
    original.sort()
    #
    # print Z-scores vector to console and to clipboard
    #
    root.clipboard_clear()
    print_results(original)
    status.set("%s", "%s p-values were adjusted" % len(original))
#
def clear_screen():
    """Clears the window
    """
    text.delete("1.0", END)
#        
def clear_data():
    """Clears the window contents and clipboard
    """
    text.delete("1.0", END)
##    original = []
    root.clipboard_clear()
    root.clipboard_append("")
    status.set("%s", "Data cleared")
#
def print_help():
    """Prints some help info to the window
    """
    clear_screen()
    help_text = \
"""Computes Benjanini-Hochberg adjusted p-values. P-values are input
via the clipboard from an Excel sheet. P-values should be displayed
many decimal places or using scientific notation. P=values do not have
to be sorted. The number of tests will be the number of p-values.
Output is original p-values and adjusted p-values displayed in window
and written to the clipboard.

Written by Phil Wilmarth, OHSU, 2013."""
    text.insert("1.0", help_text)

def quit_me():
    status.set("%s", "Bye")
    root.destroy()
#
# other functions
#
def print_data(header, original):
    """echos the data from the clipboard to the window
    """
    text.configure(tabs=("2.5c", NUMERIC))   # set tabs for numbers
    if len(original) == 0:
        text.insert("1.0", 'WARNING: Clipboard was empty!')
    else:
        if len(header) == 0:
            string = 'Row  \t %s\n' % ('p-values')
        else:
            string = 'Row  \t %s\n' % (header[0])
        text.insert("1.0", string)
        for row in original:
            text.insert(CURRENT, '%d\t%0.8f\n' % (row[0], row[1]))
        text.insert(CURRENT, '\n')                
#
def print_results(original):
    """prints data to window and clipboard
    """
    text.configure(tabs=("3.1c", NUMERIC, "6.5c", NUMERIC))   # set tabs for numbers
    text.insert("1.0", 'p-value\t   BH_adjusted\n')
    root.clipboard_append('p-value\tBH_adjusted\r')
    for row in original:
        text.insert(CURRENT, '%0.8f\t%0.8f\n' % (row[1], row[2]))
        root.clipboard_append('%0.10f\t%0.10f\r' % (row[1], row[2]))
#
#
# MAIN program starts here
#
# create root window
root = Tk()
root.title("BH p-value adjuster")

# create a toolbar
toolbar = Frame(root)

b = Button(toolbar, text='Get Data', width=8, command=get_data,
           borderwidth=2, relief=RAISED)
b.pack(side=LEFT, padx=5, pady=5)

b = Button(toolbar, text='Compute', width=8, command=compute,
           borderwidth=2, relief=RAISED)
b.pack(side=LEFT, padx=5, pady=5)

b = Button(toolbar, text='Clear', width=8, command=clear_data,
           borderwidth=2, relief=RAISED)
b.pack(side=LEFT, padx=5, pady=5)

b = Button(toolbar, text='Help', width=8, command=print_help,
           borderwidth=2, relief=RAISED)
b.pack(side=LEFT, padx=5, pady=5)

b = Button(toolbar, text='Quit', width=8, command=quit_me,
           borderwidth=2, relief=RAISED)
b.pack(side=LEFT, padx=5, pady=5)

toolbar.pack(side=TOP, fill=X)

# add a text box
text_frame = Frame(root)
text_frame.pack()
text = Text(text_frame, height=40, padx=5, pady=5)
scroll = Scrollbar(text_frame, command=text.yview)
text.configure(yscrollcommand=scroll.set)
text.pack(side=LEFT)
scroll.pack(side=RIGHT, fill=Y)

# add a status line
status = StatusBar(root)
status.pack(side=BOTTOM, fill=X)
status.set("%s", "Status line")

root.mainloop()
#
# end when user hits quit button or closes window
#
