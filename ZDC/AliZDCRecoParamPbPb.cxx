/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
// Class with ZDC reconstruction parameters                                  //
// Origin: Chiara.Oppedisano@to.infn.it                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH1D.h>

#include "AliZDCRecoParam.h"
#include "AliZDCRecoParamPbPb.h"

ClassImp(AliZDCRecoParamPbPb)

//_____________________________________________________________________________
AliZDCRecoParamPbPb::AliZDCRecoParamPbPb() :
  AliZDCRecoParam(),
  fhNpartDist(0x0),
  fhbDist(0x0),
  fClkCenter(0)
{
  //
  //Default constructor
}
//_____________________________________________________________________________
AliZDCRecoParamPbPb::AliZDCRecoParamPbPb(TH1D *hNpart, TH1D *hb, Float_t clkCent) :
  AliZDCRecoParam(),
  fhNpartDist(hNpart),
  fhbDist(hb),
  fClkCenter(clkCent)
{
  //
  //Standard constructor
}

//______________________________________________________________________________
AliZDCRecoParamPbPb::AliZDCRecoParamPbPb(const AliZDCRecoParamPbPb &oldrecopar) :
  AliZDCRecoParam(),
  fhNpartDist(0x0),
  fhbDist(0x0),
  fClkCenter(oldrecopar.fClkCenter)
{
  //Copy constructor
  if(oldrecopar.fhNpartDist){
    fhNpartDist = new TH1D(*oldrecopar.fhNpartDist);
    fhNpartDist->SetDirectory(0);
  }
  if(oldrecopar.fhbDist){
      fhbDist = new TH1D(*oldrecopar.fhbDist);
      fhbDist->SetDirectory(0);
  }
}

//_____________________________________________________________________________
AliZDCRecoParamPbPb &AliZDCRecoParamPbPb::operator =(const AliZDCRecoParamPbPb &recpar)
{
  // Equal operator.
  this->~AliZDCRecoParamPbPb();
  new(this) AliZDCRecoParamPbPb(recpar);
  return *this;  
 
}
 
//_____________________________________________________________________________
AliZDCRecoParamPbPb::~AliZDCRecoParamPbPb()
{
  // destructor

  if(fhNpartDist) delete fhNpartDist;
  if(fhbDist)     delete fhbDist;
}

//_____________________________________________________________________________
void AliZDCRecoParamPbPb::SetGlauberMCDist()
{
  // Setting Glauber MC distributions
  // from histos file stored in $ALICE_ROOT/ZDC
  TFile * fileHistos = TFile::Open("$ALICE_ROOT/ZDC/GlauberMCHistos.root");
  //
  fhNpartDist = (TH1D*) fileHistos->Get("hDist");
  fhNpartDist->SetDirectory(0);
  fhbDist = (TH1D*) fileHistos->Get("hbDist");
  fhbDist->SetDirectory(0);
  
  fileHistos->Close();
}

//_____________________________________________________________________________
AliZDCRecoParamPbPb *AliZDCRecoParamPbPb::GetHighFluxParam() 
{
  // Create high flux reco parameter
  TH1::AddDirectory(0);
  TH2::AddDirectory(0);
  //
  TFile * fileHistos = TFile::Open("$ALICE_ROOT/ZDC/GlauberMCHistos.root");
  fileHistos->cd();
  //
  TH1D* hDist = (TH1D*) fileHistos->Get("hDist");
  hDist->SetDirectory(0);
  //
  TH1D* hbDist = (TH1D*) fileHistos->Get("hbDist");
  hbDist->SetDirectory(0);
  
  AliZDCRecoParamPbPb* zdcRecoParam = new AliZDCRecoParamPbPb(hDist, hbDist, 0.1);
  //
  fileHistos->Close();
	      
  return zdcRecoParam;
  
}
