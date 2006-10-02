/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////
// This class is for the Silicon Drift Detector, SDD, specific geometry.
// It is being replaced by AliITSsegmentationSDD class. This file also
// constains classes derived from AliITSgeomSDD which do nothing but
// initilize this one with predefined values.
////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>
#include <TShape.h>

#include "AliITSgeomSDD.h"

ClassImp(AliITSgeomSDD)
AliITSgeomSDD::AliITSgeomSDD():
TObject(),
fPeriod(0.0),
fDvelocity(0.0),
fNAnodesL(0),
fNAnodesR(0),
fAnodeXL(0.0),
fAnodeXR(0.0),
fAnodeLowEdgeL(0),
fAnodeLowEdgeR(0),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
//    const Float_t kDx = 3.500;//cm. (Geant 3.12 units) Orthonormal to y and z
//    const Float_t kDy = 0.014;//cm. (Geant 3.12 units) Radialy from the Beam
//    const Float_t kDz = 3.763;//cm. (Geant 3.12 units) Allong the Beam Pipe

    fPeriod        = 0.0;
    fDvelocity     = 0.0;
    fNAnodesL      = 0;
    fNAnodesR      = 0;
    fAnodeXL       = 0.0;
    fAnodeXR       = 0.0;
    fAnodeLowEdgeL = 0;
    fAnodeLowEdgeR = 0;
}
//________________________________________________________________________
AliITSgeomSDD::AliITSgeomSDD(const Float_t *box,Float_t per,Float_t vel,
			     Float_t axL,Float_t axR,
			     Int_t nAL,Float_t *leL,
			     Int_t nAR,Float_t *leR):
TObject(),
fPeriod(0.0),
fDvelocity(0.0),
fNAnodesL(0),
fNAnodesR(0),
fAnodeXL(0.0),
fAnodeXR(0.0),
fAnodeLowEdgeL(0),
fAnodeLowEdgeR(0),
fName(),
fTitle(),
fMat(),
fDx(0.0),
fDy(0.0),
fDz(0.0){
////////////////////////////////////////////////////////////////////////
//    Standard constructor
////////////////////////////////////////////////////////////////////////
    fPeriod        = 0.0;
    fDvelocity     = 0.0;
    fNAnodesL      = 0;
    fNAnodesR      = 0;
    fAnodeXL       = 0.0;
    fAnodeXR       = 0.0;
    fAnodeLowEdgeL = 0;
    fAnodeLowEdgeR = 0;
    ResetSDD(box,per,vel,axL,axR,nAL,leL,nAR,leR);
}
//________________________________________________________________________
void AliITSgeomSDD::ResetSDD(const Float_t *box,Float_t per,Float_t vel,
			     Float_t axL,Float_t axR,
			     Int_t nAL,Float_t *leL,
			     Int_t nAR,Float_t *leR){
////////////////////////////////////////////////////////////////////////
//    Standard Filler
////////////////////////////////////////////////////////////////////////
    Int_t i;

    fPeriod        = per;
    fDvelocity     = vel;
    fNAnodesL      = nAL;
    fNAnodesR      = nAR;
    fAnodeXL       = axL;
    fAnodeXR       = axR;
    if(fAnodeLowEdgeL!=0) delete fAnodeLowEdgeL;
    fAnodeLowEdgeL = new Float_t[fNAnodesL];
    if(fAnodeLowEdgeR!=0) delete fAnodeLowEdgeR;
    fAnodeLowEdgeR = new Float_t[fNAnodesR];
    for(i=0;i<fNAnodesL;i++) fAnodeLowEdgeL[i] = leL[i];
    for(i=0;i<fNAnodesR;i++) fAnodeLowEdgeR[i] = leR[i];
    fName="ActiveSDD";
    fTitle="Active volume of SDD";
    fMat="SDD Si Det";
    fDx=box[0];
    fDy=box[1];
    fDz=box[2];
}
//________________________________________________________________________
AliITSgeomSDD::~AliITSgeomSDD(){
// Destructor

    if(fAnodeLowEdgeL!=0) delete [] fAnodeLowEdgeL;
    if(fAnodeLowEdgeR!=0) delete [] fAnodeLowEdgeR;
    fPeriod    = 0.0;
    fDvelocity = 0.0;
    fAnodeXL   = 0.0;
    fAnodeXR   = 0.0;
    fNAnodesL  = 0;
    fNAnodesR  = 0;
    fAnodeLowEdgeL = 0;
    fAnodeLowEdgeR = 0;
}
//________________________________________________________________________
AliITSgeomSDD::AliITSgeomSDD(AliITSgeomSDD &source) : TObject(source),
fPeriod(source.fPeriod),
fDvelocity(source.fDvelocity),
fNAnodesL(source.fNAnodesL),
fNAnodesR(source.fNAnodesR),
fAnodeXL(source.fAnodeXL),
fAnodeXR(source.fAnodeXR),
fAnodeLowEdgeL(0),
fAnodeLowEdgeR(0),
fName(source.fName),
fTitle(source.fTitle),
fMat(source.fMat),
fDx(source.fDx),
fDy(source.fDy),
fDz(source.fDz){
    // Copy constructor  
  fAnodeLowEdgeL = new Float_t[fNAnodesL];
  fAnodeLowEdgeR = new Float_t[fNAnodesR];
  for(Int_t i=0;i<fNAnodesL;i++) fAnodeLowEdgeL[i] = source.fAnodeLowEdgeL[i];
  for(Int_t i=0;i<fNAnodesR;i++) fAnodeLowEdgeR[i] = source.fAnodeLowEdgeR[i];
 
}
//________________________________________________________________________
AliITSgeomSDD& AliITSgeomSDD::operator=(AliITSgeomSDD &source){
    // = operator

  this->~AliITSgeomSDD();
  new(this) AliITSgeomSDD(source);
  return *this;

}
//______________________________________________________________________
void AliITSgeomSDD::Local2Det(Float_t xl,Float_t zl,Int_t &a,Int_t &t,Int_t &s){
// Give the local detector coordinate it returns the anode number, time
// bucket, and detector side.
    Int_t i;

    if(xl>0) {
	if(zl<fAnodeLowEdgeR[0]) i=-1;
	else for(i=0;i<fNAnodesR;i++) if(zl<fAnodeLowEdgeR[i]) break;
	a = i;
	s = 1;
    } else if(xl<0){
	if(zl<fAnodeLowEdgeL[0]) i=-1;
	else for(i=0;i<fNAnodesL;i++) if(zl<fAnodeLowEdgeL[i]) break;
	a = i;
	s = 0;
    } else { // x==0.
	if(zl<fAnodeLowEdgeR[0]) i=-1;
	else for(i=0;i<fNAnodesR;i++) if(zl<fAnodeLowEdgeR[i]) break;
	a = i;
	if(zl<fAnodeLowEdgeL[0]) i=-1;
	else for(i=0;i<fNAnodesL;i++) if(zl<fAnodeLowEdgeL[i]) break;
	s = -i;
    } // end if
    t = (Int_t)TMath::Abs((GetAnodeX(a,s)-xl)/fDvelocity/fPeriod);
    return;
}
//______________________________________________________________________
void AliITSgeomSDD::Det2Local(Int_t a,Int_t t,Int_t s,Float_t &xl,Float_t &zl){
// Give the anode number, time bucket, and detector side, it returns the
// local detector coordinate.

    zl = 0.5*GetAnodeZ(a,s);
    if(s==0){
	xl = GetAnodeX(a,s)+(t+0.5)*fPeriod*fDvelocity;
    } else { // s==1
	xl = GetAnodeX(a,s)-(t+0.5)*fPeriod*fDvelocity;
    } // end if s==0;
    return;
}
//______________________________________________________________________
void AliITSgeomSDD::Print(ostream *os) const {
////////////////////////////////////////////////////////////////////////
// Standard output format for this class.
////////////////////////////////////////////////////////////////////////
    Int_t i;
#if defined __GNUC__ 
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << "TBRIK" << " ";
    *os << setprecision(16) << GetDx() << " ";
    *os << setprecision(16) << GetDy() << " ";
    *os << setprecision(16) << GetDz() << " ";
    *os << setprecision(16) << fPeriod << " ";
    *os << setprecision(16) << fDvelocity << " ";
    *os << fNAnodesL << " ";
    *os << fNAnodesR << " ";
    *os << fAnodeXL << " ";
    *os << fAnodeXR << " ";
    for(i=0;i<fNAnodesL;i++) *os <<setprecision(16)<<fAnodeLowEdgeL[i]<< " ";
    for(i=0;i<fNAnodesR;i++) *os <<setprecision(16)<<fAnodeLowEdgeR[i]<< " ";
    *os << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//______________________________________________________________________
void AliITSgeomSDD::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
    Int_t i;
    char shp[20];

    *is >> shp;
    *is >> fDx >> fDy >> fDz;
    fName="AcrtiveSDD";
    fTitle="Active volulme of SDD";
    fMat="SDD Si Det";
    *is >> fPeriod >> fDvelocity >> fNAnodesL >> fNAnodesR;
    *is >> fAnodeXL >> fAnodeXR;
    if(fAnodeLowEdgeL!=0) delete fAnodeLowEdgeL;
    this->fAnodeLowEdgeL = new Float_t[fNAnodesL];
    if(fAnodeLowEdgeR!=0) delete fAnodeLowEdgeR;
    this->fAnodeLowEdgeR = new Float_t[fNAnodesR];
    for(i=0;i<fNAnodesL;i++) *is >> fAnodeLowEdgeL[i];
    for(i=0;i<fNAnodesR;i++) *is >> fAnodeLowEdgeR[i];
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITSgeomSDD &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSDD &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}

//======================================================================

ClassImp(AliITSgeomSDD256)

AliITSgeomSDD256::AliITSgeomSDD256() : AliITSgeomSDD(){
    // Default Constructor
}
//----------------------------------------------------------------------
AliITSgeomSDD256::AliITSgeomSDD256(Int_t npar,const Float_t *par) : 
    AliITSgeomSDD(){
////////////////////////////////////////////////////////////////////////
//    constructor
/*
Pads for probe cards in ALICE-D2       /05.03.2000/
(X,Y) coordinates are quoted in microns and referred to the centers of 
bonding pads. (0, 0) corrispond to the center of the detector.
Convention: 
left is for negative X, right is for positive X;
DOWN half is for negative Y, UP half is for positive Y.

Detector size: X= 87588 micrometers; Y= 72500 micrometers. 
Detector corners:
  LEFT UP: (-43794, 36250)
RIGHT UP:  (43794, 36250)
 LEFT DOWN: (-43794, -36250)
RIGHT DOWN:  (43794, -36250)

	Drift cathodes (n-side)

cathode #0 (Ubias) 
(-1477, 0), pad size (150, 60)
   (875, 0), pad size (150, 60)
(-36570, 0), pad size (200, 70)
(-37570, 0), pad size (200, 70)
(36570, 0), pad size (200, 70)
(37570, 0), pad size (200, 70)

cathode #1 DOWN half 
(-1477, -120), pad size (150, 60)
   (875, -120), pad size (150, 60)
(-36570, -120), pad size (200, 70)
(-37570, -120), pad size (200, 70)
(36570, -120), pad size (200, 70)
(37570, -120), pad size (200, 70)

cathode #2 DOWN half 
(-1477, -240), pad size (150, 60)
   (875, -240), pad size (150, 60)
(-36570, -240), pad size (200, 70)
(-37570, -240), pad size (200, 70)
(36570, -240), pad size (200, 70)
(37570, -240), pad size (200, 70)

cathode #3 DOWN half 
(-1477, -360), pad size (150, 60)
   (875, -360), pad size (150, 60)
(-36570, -360), pad size (200, 70)
(-37570, -360), pad size (200, 70)
(36570, -360), pad size (200, 70)
(37570, -360), pad size (200, 70)
.....................................
......................................
......................................
cathode #30 DOWN half
(-1477, -3600), pad size (150, 60)
   (875, -3600), pad size (150, 60)
(-36570, -3600), pad size (200, 70)
(-37570, -3600), pad size (200, 70)
(36570, -3600), pad size (200, 70)
(37570, -3600), pad size (200, 70)
...................................
cathode #60 DOWN half
(-1477, -7200), pad size (150, 60)
   (875, -7200), pad size (150, 60)
(-36570, -7200), pad size (200, 70)
(-37570, -7200), pad size (200, 70)
(36570, -7200), pad size (200, 70)
(37570, -7200), pad size (200, 70)
....................................
cathode #90 DOWN half
(-1477, -10800), pad size (150, 60)
   (875, -10800), pad size (150, 60)
(-36570, -10800), pad size (200, 70)
(-37570, -10800), pad size (200, 70)
(36570, -10800), pad size (200, 70)
(37570, -10800), pad size (200, 70)
....................................
cathode #120 DOWN half
(-1477, -14400), pad size (150, 60)
   (875, -14400), pad size (150, 60)
(-36570, -14400), pad size (200, 70)
(-37570, -14400), pad size (200, 70)
(36570, -14400), pad size (200, 70)
(37570, -14400), pad size (200, 70)
....................................
cathode #150 DOWN half
(-1477, -18000), pad size (150, 60)
   (875, -18000), pad size (150, 60)
(-36570, -18000), pad size (200, 70)
(-37570, -18000), pad size (200, 70)
(36570, -18000), pad size (200, 70)
(37570, -18000), pad size (200, 70)
....................................
cathode #180 DOWN half
(-1477, -21600), pad size (150, 60)
   (875, -21600), pad size (150, 60)
(-36570, -21600), pad size (200, 70)
(-37570, -21600), pad size (200, 70)
(36570, -21600), pad size (200, 70)
(37570, -21600), pad size (200, 70)
....................................
cathode #210 DOWN half
(-1477, -25200), pad size (150, 60)
   (875, -25200), pad size (150, 60)
(-36570, -25200), pad size (200, 70)
(-37570, -25200), pad size (200, 70)
(36570, -25200), pad size (200, 70)
(37570, -25200), pad size (200, 70)
....................................
cathode #240 DOWN half
(-1477, -28800), pad size (150, 60)
   (875, -28800), pad size (150, 60)
(-36570, -28800), pad size (200, 70)
(-37570, -28800), pad size (200, 70)
(36570, -28800), pad size (200, 70)
(37570, -28800), pad size (200, 70)
....................................
cathode #270 DOWN half
(-1477, -32400), pad size (150, 60)
   (875, -32400), pad size (150, 60)
(-36570, -32400), pad size (200, 70)
(-37570, -32400), pad size (200, 70)
(36570, -32400), pad size (200, 70)
(37570, -32400), pad size (200, 70)
....................................
cathode #290 DOWN half
(-1477, -34800), pad size (150, 60)
   (875, -34800), pad size (150, 60)
(-36570, -34800), pad size (200, 70)
(-37570, -34800), pad size (200, 70)
(36570, -34800), pad size (200, 70)
(37570, -34800), pad size (200, 70)
___________________________________________________
cathode #1 UP half 
(-1477, 120), pad size (150, 60)
   (875, 120), pad size (150, 60)
(-36570, 120), pad size (200, 70)
(-37570, 120), pad size (200, 70)
(36570, 120), pad size (200, 70)
(37570, 120), pad size (200, 70)

cathode #2 UP half 
(-1477, 240), pad size (150, 60)
   (875, 240), pad size (150, 60)
(-36570, 240), pad size (200, 70)
(-37570, 240), pad size (200, 70)
(36570, 240), pad size (200, 70)
(37570, 240), pad size (200, 70)

cathode #3 UP half 
(-1477, 360), pad size (150, 60)
   (875, 360), pad size (150, 60)
(-36570, 360), pad size (200, 70)
(-37570, 360), pad size (200, 70)
(36570, 360), pad size (200, 70)
(37570, 360), pad size (200, 70)
.....................................
......................................
......................................
cathode #30 UP half
(-1477, 3600), pad size (150, 60)
   (875, 3600), pad size (150, 60)
(-36570, 3600), pad size (200, 70)
(-37570, 3600), pad size (200, 70)
(36570, 3600), pad size (200, 70)
(37570, 3600), pad size (200, 70)
......................................
cathode #60 UP half
(-1477, 7200), pad size (150, 60)
   (875, 7200), pad size (150, 60)
(-36570, 7200), pad size (200, 70)
(-37570, 7200), pad size (200, 70)
(36570, 7200), pad size (200, 70)
(37570, 7200), pad size (200, 70)
......................................
cathode #90 UP half
(-1477, 10800), pad size (150, 60)
   (875, 10800), pad size (150, 60)
(-36570, 10800), pad size (200, 70)
(-37570, 10800), pad size (200, 70)
(36570, 10800), pad size (200, 70)
(37570, 10800), pad size (200, 70)
......................................
cathode #120 UP half
(-1477, 14400), pad size (150, 60)
   (875, 14400), pad size (150, 60)
(-36570, 14400), pad size (200, 70)
(-37570, 14400), pad size (200, 70)
(36570, 14400), pad size (200, 70)
(37570, 14400), pad size (200, 70)
......................................
cathode #150 UP half
(-1477, 18000), pad size (150, 60)
   (875, 18000), pad size (150, 60)
(-36570, 18000), pad size (200, 70)
(-37570, 18000), pad size (200, 70)
(36570, 18000), pad size (200, 70)
(37570, 18000), pad size (200, 70)
......................................
cathode #180 UP half
(-1477, 21600), pad size (150, 60)
   (875, 21600), pad size (150, 60)
(-36570, 21600), pad size (200, 70)
(-37570, 21600), pad size (200, 70)
(36570, 21600), pad size (200, 70)
(37570, 21600), pad size (200, 70)
......................................
cathode #210 UP half
(-1477, 25200), pad size (150, 60)
   (875, 25200), pad size (150, 60)
(-36570, 25200), pad size (200, 70)
(-37570, 25200), pad size (200, 70)
(36570, 25200), pad size (200, 70)
(37570, 25200), pad size (200, 70)
......................................
cathode #240 UP half
(-1477, 28800), pad size (150, 60)
   (875, 28800), pad size (150, 60)
(-36570, 28800), pad size (200, 70)
(-37570, 28800), pad size (200, 70)
(36570, 28800), pad size (200, 70)
(37570, 28800), pad size (200, 70)
......................................
cathode #270 UP half
(-1477, 32400), pad size (150, 60)
   (875, 32400), pad size (150, 60)
(-36570, 32400), pad size (200, 70)
(-37570, 32400), pad size (200, 70)
(36570, 32400), pad size (200, 70)
(37570, 32400), pad size (200, 70)
......................................
cathode #290 UP half
(-1477, 34800), pad size (150, 60)
   (875, 34800), pad size (150, 60)
(-36570, 34800), pad size (200, 70)
(-37570, 34800), pad size (200, 70)
(36570, 34800), pad size (200, 70)
(37570, 34800), pad size (200, 70)
	Injectors (n-side)

Central line injectors (DOWN half)
(-1237, -660), pad size (150, 65)
(1115, -660), pad size (150, 65)
(37890, -660), pad size (100, 74)

Middle line injectors (DOWN half)
(-1237, -17460), pad size (150, 80)
(1115, -17460), pad size (150, 80)
(37890, -17460), pad size (100, 74)

Bottom line injectors (DOWN half)
(-1237, -34020), pad size (150, 80)
(1115, -34020), pad size (150, 80)
(37890, -34020), pad size (100, 74)
___________________________________________
Central line injectors (UP half)
(-1237, 660), pad size (150, 65)
(1115, 660), pad size (150, 65)
(37890, 660), pad size (100, 74)

Middle line injectors (UP half)
(-1237, 17460), pad size (150, 80)
(1115, 17460), pad size (150, 80)
(37890, 17460), pad size (100, 74)

Bottom line injectors (UP half)
(-1237, 34020), pad size (150, 80)
(1115, 34020), pad size (150, 80)
(37890, 34020), pad size (100, 74)

Drift cathodes and injectors of p-side have the bonding pads with the same 
coordinates as for the n-side (when looking through the masks)

	Cathodes of the collection zone (n-side)

cathode #291 (-40 V) DOWN half
(-38220, -35055), pad size (120, 160)
(38190, -34992), pad size (120, 145)

GRID cathode (-15 V) DOWN half
(-37988, -35085), pad size (144, 210)
  (37988, -35085), pad size (144, 210)

cathode #292 (-30 V) DOWN half
(-38245, -35290), pad size (100, 170)
(38210, -35242), pad size (150, 215)

cathode #293 (-15 V) DOWN half
(-38055, -35460), pad size (690, 70)
(36488, -35460), pad size (3805, 70)

n+ bulk contact (GND) DOWN half
(-38300, -36050), pad size (1000, 395)
(38300, -36050), pad size (1000, 395)

bonding pad of the last integrated resistor DOWN half
/it has to be connected to the GND/
(-38190, -35620) pad size (160, 110)
________________________________________________ 
cathode #291 (-40 V) UP half
(-38220, 35055), pad size (120, 160)
(38190, 34992), pad size (120, 145)

GRID cathode (-15 V) UP half
(-37988, 35085), pad size (144, 210)
  (37988, 35085), pad size (144, 210)

cathode #292 (-30 V) UP half
(-38245, 35290), pad size (100, 170)
(38210, 35242), pad size (150, 215)

cathode #293 (-15 V) UP half
(-38055, 35460), pad size (690, 70)
(36488, 35460), pad size (3805, 70)

n+ bulk contact (GND) UP half
(-38300, 36050), pad size (1000, 395)
(38300, 36050), pad size (1000, 395)

bonding pad of the last integrated resistor UP half
/it has to be connected to the GND/
(-38190, 35620) pad size (160, 110)

Cathodes of the collection zone (p-side)

cathode #291 (-40 V) DOWN half
(-38215, -35055), pad size (120, 160)
(38190, -34992), pad size (120, 145)

cathode W1 (-60 V) DOWN half
(-38000, -35110), pad size (140, 240)
 (38000, -35110), pad size (140, 240)

cathode W2 (-80 V) DOWN half
( 0, -35090), pad size (75600, 110)

cathode #292 (-40 V) DOWN half
(-38220, -35290), pad size (150, 170)
(38210, -35242), pad size (150, 215)

p+ bulk contact (GND) DOWN half
(-38300, -36050), pad size (1000, 395)
(38300, -36050), pad size (1000, 395)

It is necessary to connect cathode #291 to cathode #292 in order to
close the integrated divider to p+ bulk contact (GND).

_______________________________________________
cathode #291 (-40 V) UP half
(-38215, 35055), pad size (120, 160)
(38190, 34992), pad size (120, 145)

cathode W1 (-60 V) UP half
(-38000, 35110), pad size (140, 240)
 (38000, 35110), pad size (140, 240)

cathode W2 (-80 V) UP half
( 0, 35090), pad size (75600, 110)

cathode #292 (-40 V) UP half
(-38220, 35290), pad size (150, 170)
(38210, 35242), pad size (150, 215)

p+ bulk contact (GND) UP half
(-38300, 36050), pad size (1000, 395)
(38300, 36050), pad size (1000, 395)

It is necessary to connect cathode #291 to cathode #292 in order to
close the integrated divider to p+ bulk contact (GND).

	Anodes (n-side)
There are 256 anodes to be bonded to the inputs of front-end electronics. In 
addition there are 2 anodes (one at the left edge and one at the right edge 
of the anode array) that have to be bonded to the ground. I call these 2 
anodes #L and #R. The pitch of all anodes is 294 micrometers.

		DOWN half anodes
#L             (-37779, -35085), pad size (184, 140)
#1             (-37485, -35085), pad size (184, 140)
.........................................
.........................................
#256.............(37485, -35085), pad size (184, 140)
#R              (37779, -35085), pad size (184, 140)
_____________________________________________
		UP half anodes
#L             (-37779, 35085), pad size (184, 140)
#1             (-37485, 35085), pad size (184, 140)
.........................................
.........................................
#256.............(37485, 35085), pad size (184, 140)
#R              (37779, 35085), pad size (184, 140)
*/
////////////////////////////////////////////////////////////////////////
//   const Float_t kDxyz[]   = {3.6250,0.01499,4.3794};//cm. (Geant 3.12 units)
                                      // Size of sensitive region of detector
    const Float_t kPeriod   = 25.0E-09; // 40 MHz sampling frequence
    const Float_t kVelocity = 5.46875E+03; // cm/s drift velocity
    const Float_t kAnodeXL  = -3.5085; // cm location in x of anodes left side
    const Float_t kAnodeXR  =  3.5085; // cm location in x of anodes right side
    const Int_t   kNAnodes  = 256;  // nuber of anodes connected
    const Float_t kAnodePitch = 0.0294; // cm
    const Float_t kAnodesZ  = -3.7485; // cm Starting location of anodes in z
    Float_t anodeLowEdges[kNAnodes+1];
    Int_t i;

    if(npar<3){
	Error("AliITSgeomSDD256","npar=%d<3. array par must be [3] or greater",
	      npar);
	return;
    } // end if
   anodeLowEdges[0] = kAnodesZ;
    for(i=0;i<kNAnodes;i++)anodeLowEdges[i+1] = kAnodePitch+anodeLowEdges[i];
    AliITSgeomSDD::ResetSDD(par,kPeriod,kVelocity,kAnodeXL,kAnodeXR,
			    kNAnodes+1,anodeLowEdges,
			    kNAnodes+1,anodeLowEdges);
}
//________________________________________________________________________
ostream &operator<<(ostream &os,AliITSgeomSDD256 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSDD256 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}

//======================================================================

ClassImp(AliITSgeomSDD300)

AliITSgeomSDD300::AliITSgeomSDD300() : AliITSgeomSDD(){
////////////////////////////////////////////////////////////////////////
//    default constructor
////////////////////////////////////////////////////////////////////////
    const Float_t kDxyz[] = {3.500,0.014,3.763};//cm.
    const Float_t kPeriod = 25.0E-09; // 40 MHz
    const Float_t kVelocity = 5.46875E+3; // cm/s
    const Int_t kNAnodes = 300; // number of anodes
    const Float_t kAnodeXL = -3.500; // cm
    const Float_t kAnodeXR = +3.500; // cm
    const Float_t kAnodesZ = -3.75; // cm
    Float_t anodeLowEdges[kNAnodes+1];
    const Float_t kanode = 0.0250;// cm anode separation.
    Int_t i;

    anodeLowEdges[0] = kAnodesZ;
    for(i=0;i<kNAnodes;i++)anodeLowEdges[i+1] = kanode+anodeLowEdges[i];
    AliITSgeomSDD::ResetSDD(kDxyz,kPeriod,kVelocity,kAnodeXL,kAnodeXR,
			    kNAnodes+1,anodeLowEdges,
			    kNAnodes+1,anodeLowEdges);
}
//________________________________________________________________________
ostream &operator<<(ostream &os,AliITSgeomSDD300 &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////

    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITSgeomSDD300 &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////

    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
