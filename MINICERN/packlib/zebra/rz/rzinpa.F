*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:47:24  mclareni
* Zebra
*
*
#include "zebra/pilot.h"
      SUBROUTINE RZINPA(CHPATH,IXDIV,LSUP,JBIAS,KEY,ICYCLE,CHOPT)
*
************************************************************************
*
*      Same as RZIN. In addition a pathname may be specified
* Input:
*   CHPATH  Pathname of the directory where the information has
*           to be read
*   IXDIV   Index of the division to receive the data structure
*           IXDIV = 0 means division 2 of the primary store
*   *LSUP*
*   JBIAS   JBIAS < 1 : LSUP is the supporting bank and JBIAS is the link
*           bias specifying where the data structure has to be introduced
*           into this bank, i.e.  the data structure will be connected to
*           LQ(LSUP+JBIAS).
*           JBIAS =  1 :  LSUP  is the  supporting link,  i.e.   the data
*           structure is connected to LSUP (top level data structure)
*           JBIAS = 2 : Stand alone data structure, no connection.
*   KEY     Keyword vector of the information to be read
*   ICYCLE  Cycle number of the key to be read
*           ICYCLE > highest cycle number means read the highest cycle
*           ICYCLE = 0 means read the lowest cycle
*           ICYCLE =  -1,-2,...  means  read the  highest cycle  minus 1,
*           2,...
*   CHOPT   Character variable specifying the options selected.
*           data structure
*             default
*                   Same as 'D' below
*             'A'   Read continuation of the previously read data structure
*                   with identifier KEYU,ICYCLE
*                   Given that option implies that the record was written with
*                   the same option by a call to RZOUT.
*             'C'   Provide   information   about   the   cycle   numbers
*                   associated with KEY.
*                   The  total number  of  cycles  and the  cycle  number
*                   identifiers of the 19 highest  cycles are returned in
*                   IQUEST(50) and IQUEST(51..89) respectively
*             'D'   Read the  Data structure  with the  (key,cycle)  pair
*                   specified.
*             'N'   Read the neighbouring. keys (i.e. those preceding and
*                   following KEY).
*                   The  key-vectors of  the previous  and  next key  are
*                   available   respectively   as   IQUEST(31..35)    and
*                   IQUEST(41..45), see below.
* Output:
*   *LSUP*  For JBIAS = 1  or 2,  LSUP contains the entry  address to the
*           data structure
*           In any case IQUEST(11) returns the entry address
*
* Called by <USER>
*
*  Author  : R.Brun DD/US/PD
*  Written : 12.04.86
*  Last mod: 12.04.86
*
************************************************************************
#include "zebra/rzcl.inc"
#include "zebra/rzclun.inc"
#include "zebra/rzch.inc"
      CHARACTER*(*) CHOPT,CHPATH
      DIMENSION KEY(*)
*
*-----------------------------------------------------------------------
*
*
      CALL RZCDIR(CHWOLD,'R')
      ISAVE=0
      CALL RZCDIR(CHPATH,' ')
      IF(IQUEST(1).NE.0)GO TO 99
      CALL RZIN(IXDIV,LSUP,JBIAS,KEY,ICYCLE,CHOPT)
      IQSAVE=IQUEST(1)
      CALL RZCDIR(CHWOLD,' ')
      IQUEST(1)=IQSAVE
      ISAVE=1
*
  99  RETURN
      END
