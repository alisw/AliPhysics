/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// Geometry parameter for the TRD
//

const Int_t   kNsect   = 18;      // Number of sectors in the full detector
const Int_t   kNplan   = 6;       // Number of planes of the TRD
const Int_t   kNcham   = 5;       // Number of chambers in z-direction

const Float_t kRmin    = 294.0;   // r-dimensions of the TRD
const Float_t kRmax    = 368.0;

const Float_t kZmax1   = 378.35;  // z-dimensions of the TRD
const Float_t kZmax2   = 302.0;

const Float_t kSheight =  74.0;   // Height of the TRD-volume in spaceframe (BTR1-3)
const Float_t kSwidth1 =  99.613; // Lower width of the TRD-volume in spaceframe (BTR1-3)
const Float_t kSwidth2 = 125.707; // Upper width of the TRD-volume in spaceframe (BTR1-3)
const Float_t kSlenTR1 = 751.0;   // Length of the TRD-volume in spaceframe (BTR1)
const Float_t kSlenTR2 = 313.5;   // Length of the TRD-volume in spaceframe (BTR2)
const Float_t kSlenTR3 = 159.5;   // Length of the TRD-volume in spaceframe (BTR3)

const Float_t kCheight =  11.0;   // Height of the chambers
const Float_t kCspace  =   1.6;   // Vertical spacing of the chambers
const Float_t kCaframe =   2.675; // Height of the aluminum frame
const Float_t kCathick =   1.0;   // Thickness of the aluminum frame
const Float_t kCcthick =   1.0;   // Thickness of the carbon frame

const Float_t kCwidcha = (kSwidth2 - kSwidth1) / kSheight * (kCheight + kCspace);
const Float_t kCcframe = kCheight - kCaframe;

// Thicknesses of the the material layers
const Float_t kSeThick = 0.02;                // Radiator seal
const Float_t kRaThick = 4.8;                 // Radiator
const Float_t kPeThick = 0.20;                // PE-layer in the radiator
const Float_t kMyThick = 0.005;               // Mylar-layer
const Float_t kXeThick = 3.5;                 // Gas mixture
const Float_t kDrThick = 3.0;                 // Drift region
const Float_t kAmThick = kXeThick - kDrThick; // Amplification region

const Float_t kCuThick = 0.001;               // Pad plane
const Float_t kSuThick = 0.06;                // HEXCEL+G10 support structure (= 0.31% X0)
const Float_t kFeThick = 0.0044;              // FEE + signal lines (= 0.31% X0)
const Float_t kCoThick = 0.02;                // PE of the cooling device
const Float_t kWaThick = 0.01;                // Cooling water

//  Position of the material layers
const Float_t kSeZpos  = -4.1525;             // Radiator seal
const Float_t kRaZpos  = -1.7425;             // Radiator
const Float_t kPeZpos  =  0.0000;             // PE-layer in the radiator
const Float_t kMyZpos  =  0.6600;             // Mylar-layer
const Float_t kDrZpos  =  2.1625;             // Drift region
const Float_t kAmZpos  =  4.1125;             // Amplification region

const Float_t kCuZpos  = -1.3370;             // Pad plane
const Float_t kSuZpos  =  0.0000;             // Support structure
const Float_t kFeZpos  =  1.3053;             // FEE + signal lines
const Float_t kCoZpos  =  1.3175;             // PE of the cooling device
const Float_t kWaZpos  =  1.3325;             // Cooling Water

//_____________________________________________________________________________
