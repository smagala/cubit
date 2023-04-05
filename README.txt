CUBIT README
============

CUBIT (CUBIT Utility BioInformatics Toolkit) is a colleciton of
bioinformatics utility scripts developed to work with influenza virus
sequences.  Many of the scripts have broader applicability, though a few
are specific to influenza.

CUBIT on the Web
----------------

CUBIT does not yet have a web presence, but it will soon.

Installing CUBIT
----------------

Please see the INSTALL.txt file.

Using CUBIT
-----------

CUBIT is composed of many small programs.  Specific help documents are
included for a few.  To find out some basic information about any of the
command line programs, open a DOS prompt in you CUBIT and run
`<program> -h`.

Users of the original james/james-tool package can still access most of
the tools from the command line, or from BioEdit, just like you have
always done.

For a small number of tools, a graphical interface has been added.
These have unique icons on the executable file once CUBIT is installed
(e.g. AA Diff).  To run these tools, double click the icon.  If you
commonly use the tool, you can make a shortcut to it by right-clicking
and selecting `Create Shortcut`.  This shortcut can be moved to your
desktop for easy access.

BioEdit Integration
-------------------

CUBIT is most useful when installed with the BioEdit package on Windows.
It has been tested with BioEdit 7.0.5.2.  An `accApp.ini` file should
have been distributed script with the cubit.zip file, or you can find a
copy in the `source/cubit_source.zip/misc/` folder.  Putting this file
in your `BioEdit/apps/` directory will make several of the command line
tools available as Accessory Applications in the BioEdit menu.

This file cannot currently be auto-installed, as it may overwrite a
previously customized `accApp.ini` file.  If you have a custom
accApp.ini file, you can append the entries for CUBIT to it to make them
available in BioEdit without overwriting your previous customizations.

Development
-----------

To develop CUBIT, you will need to install the following:

* Python 2.5
* wxPython 2.8 w/unicode support
* pyExcelerator 0.6.3a
* markdown

If you plan to compile for Windows, you should also have:

* py2exe 0.6.9
* Inno Setup 5.2.3

The CUBIT source is split into two files - `cubit_source.zip` and
`cubit_extras_source.zip`.  Only CDC influenza users need extras.  Unzip
these in the same directory.

Re-compiling CUBIT for Windows
------------------------------

* If you added new GUI tools, make sure you modify `get_windows()` in
  `setup.py`.
* Double-click `setup.py`.  It is optimized to run `setup.py py2exe`
  when double-clicked.

If there are no errors, this will generate two .zip files and two .exe
auto-installer files.

Platforms
---------

CUBIT is designed to be cross-platform.  It is known to work on:

* Windows 2000
* Windows XP
* Ubuntu Linux 7.10

It is expected to work with no modification on:

* Windows Vista
* Most *nix platforms, including:
    * Mac OS X
    * Any modern Linux distributions

Distributing CUBIT
------------------

All of the tools in `cubit.zip` can be freely distributed.  They are
licensed under the GNU GPL, and source code is included.  Feel free to
distribute within the Flu Division, to other groups in the CDC, or to
outside users.  Please also distribute `accApp.ini` and this file to
help new users get set up.

The files in cubit_extras.zip are not licensed for distribution.  These
files are specific to the Influenza Sequencing Activity, are unlikely to
be useful to anyone else, and should not be distributed.

Copyright
---------

Various parts of CUBIT are copyrighted by Dr. James Smagala, the regents
of the Univerisity of Colorado, or by the Centers for Disease Control
and Prevention.
