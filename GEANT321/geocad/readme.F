*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1995/10/24 10:20:45  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.32  by  S.Giani
*-- Author :
*
*        INFORMATION ABOUT THE GEANT-SET-INTERFACE (CADINT 1.20)
*
C    The interface is available in the CMZ patch CADINT. (CAR file CADINT.car)
C    An updated GXINT as well as command file CDF are also required.
C
C    You must recompile gxint and cdf, then you copy the interface code to your
C    directory in CMZ (YTOC cadint.car.)
C
C    You have to specify the name of the volume where you want to
C    start your conversion and the number of division instances to be written
C    into the SET file. This is useful in a case of many divisions.
C
C    Selection of volumes to export into SET is done with the visibility
C    attribute command (SATT volume SEEN x). Volume names will
C    be written into SET with an index to distinguish different instances
C    of a volume. In a case of many divisions, the indexing of slices will be
C    reset in each division.
C
C    The colors of the volumes will be transformed inte SET file as well.
C
C    When you export GEANT geometry to a CAD system,
C    be a bit careful not to have
C    too many volumes, because this will fill up the CAD system.  (Divided
C    volumes which  contains 100's of equal daughters are
C    therefore better to be SATT volume SEEN -2, or you can define the number
C    of division instances).
C
C    A good strategy is to make several SET files,
C    one for each sub part of the detector.
C    Material information is currently written into .mat file as well as
C    the tree information.
C
C    Remember that you must specify the drawing parameters again if you
C    want to draw your detector after you have used the interface command.
C
C    Contact Nils Hoimyr/CN/CE or Jouko Vuoskoski/CN/AS if you have any
C    questions regarding the interface.
 
 
