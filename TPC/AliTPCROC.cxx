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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Geometry        class for a single ROC                                   //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCROC.h"
#include "TMath.h"

ClassImp(AliTPCROC)


AliTPCROC* AliTPCROC::fgInstance = 0;




//_ singleton implementation __________________________________________________
AliTPCROC* AliTPCROC::Instance()
{
  //
  // Singleton implementation
  // Returns an instance of this class, it is created if neccessary
  //
  if (fgInstance == 0){
    fgInstance = new AliTPCROC();
    fgInstance->Init();    
  }
  return fgInstance;
}




void AliTPCROC::Init(){
  //
  // initialize static variables
  //
  if (AliTPCROC::fNSectorsAll>0) return;
  fNSectorsAll =72;
  fNSectors[0] =36;
  fNSectors[1] =36;
  //
  fNRows[0]= 63;
  fNRows[1]= 96;
  //
  // number of pads in padrow
  fNPads[0] = new UInt_t[fNRows[0]];
  fNPads[1] = new UInt_t[fNRows[1]];  
  //
  // padrow index in array
  //
  fRowPosIndex[0] = new UInt_t[fNRows[0]];
  fRowPosIndex[1] = new UInt_t[fNRows[1]];
  //
  // inner sectors
  //
  UInt_t index =0;
  for (UInt_t irow=0; irow<fNRows[0];irow++){
    UInt_t npads = (irow==0) ? 68 : 2 *Int_t(Double_t(irow)/3. +33.67);
    fNPads[0][irow] = npads;
    fRowPosIndex[0][irow] = index;
    index+=npads;
  }
  fNChannels[0] = index;
  //
  index =0;
  Double_t k1 = 10.*TMath::Tan(10*TMath::DegToRad())/6.;
  Double_t k2 = 15.*TMath::Tan(10*TMath::DegToRad())/6.;
  for (UInt_t irow=0; irow<fNRows[1];irow++){    
    UInt_t npads = (irow<64) ? 
      2*Int_t(k1*Double_t(irow)+37.75):
      2*Int_t(k2*Double_t(irow-64)+56.66);
    fNPads[1][irow] = npads;
    fRowPosIndex[1][irow] = index;
    index+=npads;
  }
  fNChannels[1] = index;
  SetGeometry();
}




void AliTPCROC::SetGeometry()
{
  //
  //set ROC geometry parameters
  //
  const  Float_t kInnerRadiusLow = 83.65;
  const  Float_t kInnerRadiusUp  = 133.3;
  const  Float_t kOuterRadiusLow = 133.5;
  const  Float_t kOuterRadiusUp  = 247.7;
  const  Float_t kInnerFrameSpace = 1.5;
  const  Float_t kOuterFrameSpace = 1.5;
  const  Float_t kInnerWireMount = 1.2;
  const  Float_t kOuterWireMount = 1.4;
  const  Float_t kZLength =250.;
  const  UInt_t   kNRowLow = 63;
  const  UInt_t   kNRowUp1 = 64;
  const  UInt_t   kNRowUp2 = 32;
  const  UInt_t   kNRowUp  = 96;
  const  Float_t kInnerAngle = 20; // 20 degrees
  const  Float_t kOuterAngle = 20; // 20 degrees
  //
  //  pad     parameters
  // 
  const Float_t  kInnerPadPitchLength = 0.75;
  const Float_t  kInnerPadPitchWidth = 0.40;
  const Float_t  kInnerPadLength = 0.75;
  const Float_t  kInnerPadWidth = 0.40;
  const Float_t  kOuter1PadPitchLength = 1.0;
  const Float_t  kOuterPadPitchWidth = 0.6;
  const Float_t  kOuter1PadLength = 1.0;
  const Float_t  kOuterPadWidth = 0.6;
  const Float_t  kOuter2PadPitchLength = 1.5;
  const Float_t  kOuter2PadLength = 1.5;  

  //
  //wires default parameters
  //
//   const UInt_t    kNInnerWiresPerPad = 3;
//   const UInt_t    kInnerDummyWire = 2;
//   const Float_t  kInnerWWPitch = 0.25;
//   const Float_t  kRInnerFirstWire = 84.475;
//   const Float_t  kRInnerLastWire = 132.475;
//   const Float_t  kInnerOffWire = 0.5;
//   const UInt_t    kNOuter1WiresPerPad = 4;
//   const UInt_t    kNOuter2WiresPerPad = 6;
//   const Float_t  kOuterWWPitch = 0.25;  
//   const Float_t  kROuterFirstWire = 134.225;
//   const Float_t  kROuterLastWire = 246.975;
//   const UInt_t    kOuterDummyWire = 2;
//   const Float_t  kOuterOffWire = 0.5;
  //
  //set sector parameters
  //
  fInnerRadiusLow = kInnerRadiusLow;
  fOuterRadiusLow = kOuterRadiusLow;
  fInnerRadiusUp  = kInnerRadiusUp;
  fOuterRadiusUp  = kOuterRadiusUp;  
  fInnerFrameSpace = kInnerFrameSpace;
  fOuterFrameSpace = kOuterFrameSpace;
  fInnerWireMount  = kInnerWireMount;
  fOuterWireMount  = kOuterWireMount;
  fZLength         = kZLength;
  fInnerAngle      =  TMath::DegToRad()*kInnerAngle;
  fOuterAngle      =  TMath::DegToRad()*kOuterAngle;

  fNRowLow       = kNRowLow;
  fNRowUp1      = kNRowUp1;
  fNRowUp2       = kNRowUp2;
  fNRowUp        = kNRowUp;
  //
  //set pad parameter
  //
  fInnerPadPitchLength = kInnerPadPitchLength;
  fInnerPadPitchWidth  = kInnerPadPitchWidth;
  fInnerPadLength      = kInnerPadLength;
  fInnerPadWidth       = kInnerPadWidth;
  fOuter1PadPitchLength = kOuter1PadPitchLength; 
  fOuter2PadPitchLength = kOuter2PadPitchLength;
  fOuterPadPitchWidth   = kOuterPadPitchWidth;
  fOuter1PadLength      = kOuter1PadLength;
  fOuter2PadLength      = kOuter2PadLength;
  fOuterPadWidth        = kOuterPadWidth; 

  //
  //set wire parameters
  //
  // SetInnerNWires(kNInnerWiresPerPad);
  //   SetInnerDummyWire(kInnerDummyWire);
  //   SetInnerOffWire(kInnerOffWire);
  //   SetOuter1NWires(kNOuter1WiresPerPad);
  //   SetOuter2NWire(kNOuter2WiresPerPad);
  //   SetOuterDummyWire(kOuterDummyWire);
  //   SetOuterOffWire(kOuterOffWire);
  //   SetInnerWWPitch(kInnerWWPitch);
  //   SetRInnerFirstWire(kRInnerFirstWire);
  //   SetRInnerLastWire(kRInnerLastWire);
  //   SetOuterWWPitch(kOuterWWPitch);
  //   SetROuterFirstWire(kROuterFirstWire);
  //   SetROuterLastWire(kROuterLastWire);  
} 




//_____________________________________________________________________________
AliTPCROC::AliTPCROC()
          :TObject(), 
           fNSectorsAll(0),
	   fInnerRadiusLow(0.),
	   fInnerRadiusUp(0.),
	   fOuterRadiusUp(0.),
	   fOuterRadiusLow(0.),
	   fInnerFrameSpace(0.),
	   fOuterFrameSpace(0.),
	   fInnerWireMount(0.),
	   fOuterWireMount(0.),
	   fZLength(0.),
	   fInnerAngle(0.),
	   fOuterAngle(0.),
	   fNInnerWiresPerPad(0),
	   fInnerWWPitch(0.),
	   fInnerDummyWire(0),
	   fInnerOffWire(0.),
	   fRInnerFirstWire(0.),
	   fRInnerLastWire(0.),
	   fLastWireUp1(0.),
	   fNOuter1WiresPerPad(0),
	   fNOuter2WiresPerPad(0),
	   fOuterWWPitch(0.),
	   fOuterDummyWire(0),
	   fOuterOffWire(0),
	   fROuterFirstWire(0.),
	   fROuterLastWire(0),
	   fInnerPadPitchLength(0.),
	   fInnerPadPitchWidth(0.),
	   fInnerPadLength(0.),
	   fInnerPadWidth(0.),
	   fOuter1PadPitchLength(0.),
	   fOuter2PadPitchLength(0),
	   fOuterPadPitchWidth(0),
	   fOuter1PadLength(0.),
	   fOuter2PadLength(0),
	   fOuterPadWidth(0),
	   fNRowLow(0),
	   fNRowUp1(0),
	   fNRowUp2(0),
	   fNRowUp(0),
	   fNtRows(0)
{
  //
  // Default constructor
  for (UInt_t i=0;i<2;i++){
    fNSectors[i]  = 0;
    fNRows[i]     = 0;
    fNChannels[i] = 0;
    fNPads[i]     = 0;
    fRowPosIndex[i]= 0;
  }
}


//_____________________________________________________________________________
AliTPCROC::AliTPCROC(const AliTPCROC &roc)
          :TObject(roc),
           fNSectorsAll(0),
	   fInnerRadiusLow(0.),
	   fInnerRadiusUp(0.),
	   fOuterRadiusUp(0.),
	   fOuterRadiusLow(0.),
	   fInnerFrameSpace(0.),
	   fOuterFrameSpace(0.),
	   fInnerWireMount(0.),
	   fOuterWireMount(0.),
	   fZLength(0.),
	   fInnerAngle(0.),
	   fOuterAngle(0.),
	   fNInnerWiresPerPad(0),
	   fInnerWWPitch(0.),
	   fInnerDummyWire(0),
	   fInnerOffWire(0.),
	   fRInnerFirstWire(0.),
	   fRInnerLastWire(0.),
	   fLastWireUp1(0.),
	   fNOuter1WiresPerPad(0),
	   fNOuter2WiresPerPad(0),
	   fOuterWWPitch(0.),
	   fOuterDummyWire(0),
	   fOuterOffWire(0),
	   fROuterFirstWire(0.),
	   fROuterLastWire(0),
	   fInnerPadPitchLength(0.),
	   fInnerPadPitchWidth(0.),
	   fInnerPadLength(0.),
	   fInnerPadWidth(0.),
	   fOuter1PadPitchLength(0.),
	   fOuter2PadPitchLength(0),
	   fOuterPadPitchWidth(0),
	   fOuter1PadLength(0.),
	   fOuter2PadLength(0),
	   fOuterPadWidth(0),
	   fNRowLow(0),
	   fNRowUp1(0),
	   fNRowUp2(0),
	   fNRowUp(0),
	   fNtRows(0)

{
  //
  // AliTPCROC copy constructor
  //
  fNSectorsAll = roc.fNSectorsAll;
  fNSectors[0] = roc.fNSectors[0];
  fNSectors[1] = roc.fNSectors[1];
  fNRows[0]    = roc.fNRows[0];
  fNRows[1]    = roc.fNRows[1];
  fNChannels[0]= roc.fNChannels[0];
  fNChannels[1]= roc.fNChannels[1];
  //
  // number of pads in padrow
  fNPads[0] = new UInt_t[fNRows[0]];
  fNPads[1] = new UInt_t[fNRows[1]];  
  //
  // padrow index in array
  //
  fRowPosIndex[0] = new UInt_t[fNRows[0]];
  fRowPosIndex[1] = new UInt_t[fNRows[1]];
  //
  for (UInt_t irow =0; irow<fNRows[0];irow++){
    fNPads[0][irow]       = roc.fNPads[0][irow];
    fRowPosIndex[0][irow] = roc.fRowPosIndex[0][irow];
  }
  for (UInt_t irow =0; irow<fNRows[1];irow++){
    fNPads[1][irow]       = roc.fNPads[1][irow];
    fRowPosIndex[1][irow] = roc.fRowPosIndex[1][irow];
  }
}
//____________________________________________________________________________
AliTPCROC & AliTPCROC::operator =(const AliTPCROC & roc)
{
  //
  // assignment operator - dummy
  //
  fZLength = roc.fZLength;
  return (*this);
}
//_____________________________________________________________________________
AliTPCROC::~AliTPCROC()
{
  //
  // AliTPCROC destructor
  //
  delete [] fNPads[0];
  delete [] fNPads[1];
  delete [] fRowPosIndex[0];
  delete [] fRowPosIndex[1];
  fgInstance = 0x0;
  
}

