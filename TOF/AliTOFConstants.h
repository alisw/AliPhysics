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
// Author: Jiri Chudoba (CERN), F. Pierella
//
////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTOFConstants {
 public:
    // return number of chambers
    static const Int_t fgkNStripA; // number of strips in A type module 
    static const Int_t fgkNStripB; // number of strips in B type module 
    static const Int_t fgkNStripC; // number of strips in C type module 
    static const Int_t fgkNpadX;   // Number of pads in a strip along the X direction
    static const Int_t fgkNpadZ; // Number of pads in a strip along the Z direction
    static const Int_t fgkPadXSector;
    static const Int_t fgkNSectors;
    static const Int_t fgkNPlates;

    static const Float_t fgkrmin; // inner radius of the TOF detector (cm)
    static const Float_t fgkrmax; // outer radius of the TOF detector (cm)
    static const Int_t fgkmaxtoftree;    // numer of geom. levels: 
                              // 1 - sector, 2 - module(plate), 3 - strip, 4 - padZ, 5 - padX
    static const Int_t fgkmaxNstrip;   //20 - max. number of strips, A - 15, B - 19, C - 20
    static const Int_t fgkPadXStrip; // number of pads per strip
    static const Float_t fgkzlenA;// length (cm) of the A module, need for generation of add. noise
    static const Float_t fgkzlenB;// length (cm) of the B module, need for generation of add. noise
    static const Float_t fgkzlenC;// length (cm) of the C module, need for generation of add. noise
    static const Float_t fgkXPad;  //size of a pad in the x direction (cm)
    static const Float_t fgkZPad;  //size of a pad in the z direction (cm)
    static const Float_t fgkMaxhZtof;//max.half z-size of TOF (cm)
    static const Float_t fgkSigmaForTail1;// sigma for simulation of tails in TDC (1)
    static const Float_t fgkSigmaForTail2;// sigma for simulation of tails in TDC (2)

// if two signals ar eseparated less than fgkTimeDiff, they are merged
// and considered as one     
    static const Int_t fgkTimeDiff; // time in ps, 
    // speed of light (used in reconstruction) given in 10^9 m/s
    static const Float_t fgkSpeedOfLight; //
    // mass values for pi/K/p/e/mu in [GeV/c^2]
    // used in reconstruction
    static const Float_t fgkPionMass; // charged pion mass
    static const Float_t fgkKaonMass; // charged kaon mass
    static const Float_t fgkProtonMass; // proton mass   
    static const Float_t fgkElectronMass; // electron mass
    static const Float_t fgkMuonMass; // muon mass

 private:
    AliTOFConstants(){}
    virtual ~AliTOFConstants(){}

    ClassDef(AliTOFConstants, 0)             // TOF global constants 
};
	
#endif
