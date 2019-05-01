There are two classes with hold the information
TInfo -> Temperature
LInfo -> Led info

They are filled with readOCDB macro which for a set of run numbers
given per period will make an TObjArray of the info objects and store
to either the tempinfo.root or ledinfo.root files.

Then there are two macros plotOCDB_Temperature.C and plotOCDB_LED.C
which define two classes TDraw and LDraw to plot information.
The macro that does this efficiently is readOCDB.C which defines
a function for every period. If you want to add a period follow the
structure there. Note that the function read_all should then be 
extended by the call for the new period.

Then the createTree.C macro will be used to make trees with the relevant
info stored as TCalInfo and TCalCell
The TCalCell class relates Led and temperature information, and should
be used to plot the led/ledmon vs T for various runs.

A macro called anaTree exists, that shows how to do that.

In case of questions please ask C.Loizides.


