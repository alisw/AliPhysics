#ifndef ALITOFCONSTANTS_H
#define ALITOFCONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

////////////////////////////////////////////////////////////////////////
//
// AliTOFConstants class
//
// This class serves to group constants needed by TOF detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
//
//
// Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTOFConstants {
 public:
    // return number of chambers
    static const Int_t fgkNStripA = 15; // number of strips in A type module 
    static const Int_t fgkNStripB = 19; // number of strips in B type module 
    static const Int_t fgkNStripC = 20; // number of strips in C type module 
    static const Int_t fgkNpadX   = 48; // Number of pads in a strip along the X direction
    static const Int_t fgkNpadZ   =  2; // Number of pads in a strip along the Z direction
    static const Int_t fgkPadXSector =
      (fgkNStripA + 2*fgkNStripB + 2*fgkNStripC)*fgkNpadX*fgkNpadZ;
    static const Int_t fgkNSectors   =  18;
    static const Int_t fgkNPlates    =  5;

// if two signals ar eseparated less than fgkTimeDiff, they are merged
// and considered as one     
    static const Int_t fgkTimeDiff   =  25000; // time in ps, 
    
    
 private:
    AliTOFConstants(){}
    virtual ~AliTOFConstants(){}

    ClassDef(AliTOFConstants, 0)             // TOF global constants 
};
	
#endif
