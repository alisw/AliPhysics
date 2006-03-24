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
//  Calibration base class for a single ROC                                  //
//  Contains one float value per pad                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalROC.h"
#include "TMath.h"

ClassImp(AliTPCCalROC)
  Int_t  AliTPCCalROC::fgNSectorsAll =0;
Int_t  AliTPCCalROC::fgNSectors[2]={0,0};
Int_t  AliTPCCalROC::fgNRows[2]={0,0};
Int_t *AliTPCCalROC::fgNPads[2]={0,0};
Int_t *AliTPCCalROC::fgRowPosIndex[2] ={0,0}; 
Int_t  AliTPCCalROC::fgNChannels[2]={0,0};

void AliTPCCalROC::Init(){
  //
  // initialize static variables
  //
  if (AliTPCCalROC::fgNSectorsAll>0) return;
  fgNSectorsAll =72;
  fgNSectors[0] =36;
  fgNSectors[1] =36;
  //
  fgNRows[0]= 63;
  fgNRows[1]= 96;
  //
  // number of pads in padrow
  fgNPads[0] = new Int_t[fgNRows[0]];
  fgNPads[1] = new Int_t[fgNRows[1]];  
  //
  // padrow index in array
  //
  fgRowPosIndex[0] = new Int_t[fgNRows[0]];
  fgRowPosIndex[1] = new Int_t[fgNRows[1]];
  //
  // inner sectors
  //
  Int_t index =0;
  for (Int_t irow=0; irow<fgNRows[0];irow++){
    Int_t npads = (irow==0) ? 68 : 2 *Int_t(Double_t(irow)/3. +33.67);
    fgNPads[0][irow] = npads;
    fgRowPosIndex[0][irow] = index;
    index+=npads;
  }
  fgNChannels[0] = index;
  //
  index =0;
  Double_t k1 = 10.*TMath::Tan(10*TMath::DegToRad())/6.;
  Double_t k2 = 15.*TMath::Tan(10*TMath::DegToRad())/6.;
  for (Int_t irow=0; irow<fgNRows[1];irow++){    
    Int_t npads = (irow<64) ? 
      2*Int_t(k1*Double_t(irow)+37.75):
      2*Int_t(k2*Double_t(irow-64)+56.66);
    fgNPads[1][irow] = npads;
    fgRowPosIndex[1][irow] = index;
    index+=npads;
  }
  fgNChannels[1] = index;
}


//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC():TObject()
{
  //
  // Default constructor
  //
  fSector       = -1;
  fIndex        = 0;
  fData         = 0;
}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(Int_t sector):TObject()
{
  //
  // Constructor that initializes a given sector
  //
  Init();
  fSector = sector;
  fIndex  = (sector<fgNSectors[0]) ? 0:1;      
  fData = new Float_t[fgNChannels[fIndex]];
}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(const AliTPCCalROC &c):TObject(c)
{
  //
  // AliTPCCalROC copy constructor
  //
  fSector = c.fSector;
  fIndex  = c.fIndex;
  Int_t nchannels =  fgNChannels[fIndex];
  fData   = new Float_t[nchannels];
  for (Int_t idata = 0; idata< nchannels; idata++) fData[idata] = c.fData[idata];
}

//_____________________________________________________________________________
AliTPCCalROC::~AliTPCCalROC()
{
  //
  // AliTPCCalROC destructor
  //

  if (fData) {
    delete [] fData;
    fData = 0;
  }
}

