#!/usr/bin/env python2.7
'''
Interactive dialog to select file(s).
Created on 29 Jun 2015

:author: yaswant.pradhan
'''

from os import getcwd
from Tkinter import Tk
# import os.path as path


def pickfiles(initial_dir=None, extension=None):
    """
    Dialog picker for one or more files.

    Parameters
    ----------
    initial_dir : str
        initial directory to select files (default ``os.getcwd()``)
    extension : str or list of str, optional
        list of file extensions (default '*.*')

    Returns
    -------
    str or list of str
        selected filename[s]

    """
    from tkFileDialog import askopenfilenames
    # remove the empty application window that pops up behind the file dialogs
    Tk().withdraw()

    # parse input parameters
    initial_dir = [initial_dir, getcwd()][initial_dir is None]
    extension = [extension, '*.*'][extension is None]
    if isinstance(extension, list) is False:
        extension = extension.split()

    file_type = zip(['file type'] * len(extension), extension)
    files = askopenfilenames(
        initialdir=initial_dir, filetypes=file_type, title='Pick file[s]')

    if files:
        out = [f for f in files] if len(files) > 1 else files[0]
        return out
    else:
        print 'Warning! Nothing picked, None returned.'


# def pickfile(initialdir=None, filetype=[('All types', '*')]):
#    """
#    Single file selector dialog.
#
#    Parameters
#    ----------
#    initialdir : str
#        initial directory to select files
#    filetype : list, optional
#        list of tuples in ('file type', 'file extension') format
#
#    Returns
#    -------
#    str
#        selected filename
#
#    Raises
#    ------
#    ValueError
#        When nothing selected from the dialog box
#    """
#
#    from tkFileDialog import askopenfilename
#
#    # remove the empty application window that pops up behind the file dialogs
#    Tk().withdraw()
#
#    # set initial directory (OS independent)
#    home = path.expanduser("~")
#    if not initialdir:
#        initialdir = home
#    if initialdir and not path.exists(initialdir):
#        print "** {0} doesnot exist. Please select a valid path **".\
#            format(initialdir)
#        initialdir = home
#
#    fi = askopenfilename(initialdir=initialdir,
#                         filetypes=filetype,
#                         title='Pick a file')
#    if not fi:
#        raise ValueError("Nothingpicked.")
#    return fi
#

# def pickfiles(initialdir=None, filetype=[('All types', '*')]):
#    """
#    Multiple file selector dialog.
#
#    Parameters
#    ----------
#    initialdir : str
#        initial directory to select files
#    filetype : list, optional
#        list of tuples in ('file type', 'file extension') format
#
#    Returns
#    -------
#    tuple of str
#        Selected filenames
#
#    Raises
#    ------
#    ValuError
#        when nothing selected from the dialog box.
#    """
#
#    from tkFileDialog import askopenfilenames
#
#    # remove the empty application window that pops up behind the file dialogs
#    Tk().withdraw()
#
#    # set initial directory (OS independent)
#    home = path.expanduser("~")
#    if not initialdir:
#        initialdir = home
#    if initialdir and not path.exists(initialdir):
#        print "** {0} doesnot exist. Please select a valid path **".\
#            format(initialdir)
#        initialdir = home
#    files = askopenfilenames(initialdir=initialdir,
#                             filetypes=filetype,
#                             title='Pick file(s)')
#    if not files:
#        raise UserWarning("Nothing picked.")
#    return files
#

if __name__ == '__main__':
    pass
