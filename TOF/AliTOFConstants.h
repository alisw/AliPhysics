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
    static const Int_t fgkNStripA = 15; // number of strips in A type module 
    static const Int_t fgkNStripB = 19; // number of strips in B type module 
    static const Int_t fgkNStripC = 20; // number of strips in C type module 
    static const Int_t fgkNpadX   = 48; // Number of pads in a strip along the X direction
    static const Int_t fgkNpadZ   =  2; // Number of pads in a strip along the Z direction
    static const Int_t fgkPadXSector =
      (fgkNStripA + 2*fgkNStripB + 2*fgkNStripC)*fgkNpadX*fgkNpadZ;
    static const Int_t fgkNSectors   =  18;
    static const Int_t fgkNPlates    =  5;

    static const Float_t fgkrmin     = 370.; // inner radius of the TOF detector (cm)
    static const Float_t fgkrmax     = 399.; // outer radius of the TOF detector (cm)
    static const Int_t fgkmaxtoftree = 5;    // numer of geom. levels: 
                              // 1 - sector, 2 - module(plate), 3 - strip, 4 - padZ, 5 - padX
    static const Int_t fgkmaxNstrip  = 20;   //20 - max. number of strips, A - 15, B - 19, C - 20
    static const Int_t fgkPadXStrip  = fgkNpadX*fgkNpadZ; // number of pads per strip
    static const Float_t fgkzlenA    = 106.0;// length (cm) of the A module, need for generation of add. noise
    static const Float_t fgkzlenB    = 141.0;// length (cm) of the B module, need for generation of add. noise
    static const Float_t fgkzlenC    = 177.5;// length (cm) of the C module, need for generation of add. noise
    static const Float_t fgkXPad     = 2.5;  //size of a pad in the x direction (cm)
    static const Float_t fgkZPad     = 3.5;  //size of a pad in the z direction (cm)
    static const Float_t fgkMaxhZtof = 371.5;//max.half z-size of TOF (cm)
    static const Float_t fgkSigmaForTail1= 2.;// sigma for simulation of tails in TDC (1)
    static const Float_t fgkSigmaForTail2=0.5;// sigma for simulation of tails in TDC (2)

// if two signals ar eseparated less than fgkTimeDiff, they are merged
// and considered as one     
    static const Int_t fgkTimeDiff   =  25000; // time in ps, 
    // speed of light (used in reconstruction) given in 10^9 m/s
    static const Float_t fgkSpeedOfLight   =  0.299792458; //
    // mass values for pi/K/p/e/mu in [GeV/c^2]
    // used in reconstruction
    static const Float_t fgkPionMass     = 0.13957; // charged pion mass
    static const Float_t fgkKaonMass     = 0.49368; // charged kaon mass
    static const Float_t fgkProtonMass   = 0.93827; // proton mass   
    static const Float_t fgkElectronMass = 0.00051; // electron mass
    static const Float_t fgkMuonMass     = 0.10566; // muon mass

 private:
    AliTOFConstants(){}
    virtual ~AliTOFConstants(){}

    ClassDef(AliTOFConstants, 0)             // TOF global constants 
};
	
#endif
