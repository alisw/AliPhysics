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
void AliZDCRecoParamPbPb::SetGlauberMCDist(Float_t beamEnergy)
{
  // Setting Glauber MC distributions
  // from histos file stored in $ALICE_ROOT/ZDC
  TH1::AddDirectory(0);
  TH2::AddDirectory(0);
  
  TFile *fileGlauberMC =  TFile::Open("$ALICE_ROOT/ZDC/GlauberMCDist.root");
  if(!fileGlauberMC) {
    AliError((" Opening file $ALICE_ROOT/ZDC/SpectatorSignal.root failed\n"));
    return;
  }
  
  Float_t sqrtS = 2*beamEnergy;
  //
  if(TMath::Abs(sqrtS-5500) < 100.){
    AliDebug(2, " ZDC -> Looking for energy5500 in file $ALICE_ROOT/ZDC/GlauberMCDist.root");
    fileGlauberMC->cd("energy5500");
    fileGlauberMC->GetObject("energy5500/hNpartGlauber;1", fhNpartDist);
    if(!fhNpartDist) AliError("  PROBLEM!!! Can't get Glauber MC Npart distribution from file GlauberMCDist.root\n");
    fileGlauberMC->GetObject("energy5500/hbGlauber;1", fhbDist);
    if(!fhbDist) AliError("  PROBLEM!!! Can't get Glauber MC b distribution from file GlauberMCDist.root\n");
  }
  else if(TMath::Abs(sqrtS-2760) < 100.){
    AliDebug(2, " ZDC -> Looking for energy2760 in file $ALICE_ROOT/ZDC/GlauberMCDist.root");
    fileGlauberMC->cd("energy2760");
    fileGlauberMC->GetObject("energy2760/hNpartGlauber;1", fhNpartDist);
    if(!fhNpartDist) AliError("  PROBLEM!!! Can't get Glauber MC Npart distribution from file GlauberMCDist.root\n");
    fileGlauberMC->GetObject("energy2760/hbGlauber;1", fhbDist);
    if(!fhbDist) AliError("  PROBLEM!!! Can't get Glauber MC b distribution from file GlauberMCDist.root\n");
  }
  else AliError(Form(" No AliZDCRecoParam provided for Pb-Pb @ sqrt(s) = %1.0f GeV\n", sqrtS));
  //
  fhNpartDist->SetDirectory(0);
  fhbDist->SetDirectory(0);
  
  fileGlauberMC->Close();
}

//_____________________________________________________________________________
AliZDCRecoParamPbPb *AliZDCRecoParamPbPb::GetHighFluxParam(Float_t beamEnergy) 
{
  // Create high flux reco parameter
  TH1::AddDirectory(0);
  TH2::AddDirectory(0);
  //
  TFile *fileGlauberMC =  TFile::Open("$ALICE_ROOT/ZDC/GlauberMCDist.root");
  if(!fileGlauberMC) printf(" AliZDCRecoParamPbPb::GetHighFluxParam() ERROR opening file $ALICE_ROOT/ZDC/SpectatorSignal.root\n");
  
  TH1D *hNpartDist, *hbDist;
  if(TMath::Abs(beamEnergy-5500)<100.){
    fileGlauberMC->cd("energy5500");
    fileGlauberMC->GetObject("energy5500/hNpartGlauber;1", hNpartDist);
    if(!hNpartDist) printf("  AliZDCRecoParamPbPb::GetHighFluxParam() PROBLEM!!! Can't get Glauber MC Npart distribution from file GlauberMCDist.root\n");
    fileGlauberMC->GetObject("energy5500/hbGlauber;1", hbDist);
    if(!hbDist) printf("  AliZDCRecoParamPbPb::GetHighFluxParam() PROBLEM!!! Can't get Glauber MC b distribution from file GlauberMCDist.root\n");
  }
  else if(TMath::Abs(beamEnergy-2760)<100.){
    fileGlauberMC->cd("energy2760");
    fileGlauberMC->GetObject("energy2760/hNpartGlauber;1", hNpartDist);
    if(!hNpartDist) printf("  PROBLEM!!! Can't get Glauber MC Npart distribution from file GlauberMCDist.root\n");
    fileGlauberMC->GetObject("energy2760/hbGlauber;1", hbDist);
    if(!hbDist) printf("  AliZDCRecoParamPbPb::GetHighFluxParam() PROBLEM!!! Can't get Glauber MC b distribution from file GlauberMCDist.root\n");
  }
  else printf(" No AliZDCRecoParam provided for Pb-Pb @ sqrt(s) = %1.0f GeV\n", beamEnergy);
  //
  hNpartDist->SetDirectory(0);
  hbDist->SetDirectory(0);

  AliZDCRecoParamPbPb* zdcRecoParam = new AliZDCRecoParamPbPb(hNpartDist, hbDist, 0.1);
  //
  fileGlauberMC->Close();
	      
  return zdcRecoParam;
  
}
