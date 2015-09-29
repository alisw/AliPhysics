// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

///  @file   AliHLTLumiRegComponent.cxx
///  @author Davide Caffarri <davide.caffarri@cern.ch>
///  @date   January 2015
///  @brief  Luminous Region calculation component

////////////////////////////////////////////////////////////
//                                                        //
// a Luminous Region determination component for the HLT  //
//                                                        //
////////////////////////////////////////////////////////////

#include "AliHLTLumiRegComponent.h"
#include "TRegexp.h" 
#include "TString.h" 
#include "AliSysInfo.h" 
#include "AliESDVertex.h"
#include "TF1.h"
#include "AliHLTMessage.h"

using namespace std;

ClassImp( AliHLTLumiRegComponent )
AliHLTLumiRegComponent::AliHLTLumiRegComponent()
:
fPushBackPeriodLHC(60),
fPushBackPeriodDQM(300),
fLastPushBackTime(0),
fEventSpecie(0)
{
  for(Int_t i=0; i<3; i++){
    fPrimaryLHC[i] = 0x0;
    fPrimaryDQM[i] = 0x0;
  }
  
  for(Int_t i=0; i<2; i++){
    fPrimaryDefMultLHC[i] = 0x0;
    fPrimaryDefMultDQM[i] = 0x0;
  }
  
}

AliHLTLumiRegComponent::AliHLTLumiRegComponent( const AliHLTLumiRegComponent& )
:AliHLTProcessor()
{
    // see header file for class documentation
    HLTFatal( "copy constructor untested" );
}

AliHLTLumiRegComponent& AliHLTLumiRegComponent::operator=( const AliHLTLumiRegComponent& )
{
    // see header file for class documentation
    HLTFatal( "assignment operator untested" );
    return *this;
}

AliHLTLumiRegComponent::~AliHLTLumiRegComponent()
{
  // see header file for class documentation
  
  for(Int_t i=0; i<3; i++){
    delete fPrimaryLHC[i];
    delete fPrimaryDQM[i];
  }

}

AliHLTComponentDataType AliHLTLumiRegComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;
}

void AliHLTLumiRegComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list )
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS); //spd vertex
  list.push_back(kAliHLTDataTypeFlatESDVertex|kAliHLTDataOriginITS); // ITSSA primary vertex
  }

void AliHLTLumiRegComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //8000 is a temporary number. Find the way to
  //estimate 3 serialized TH2F.
  constBase = 80000;
  inputMultiplier = 1.;

}

AliHLTComponent *AliHLTLumiRegComponent::Spawn(){

return new AliHLTLumiRegComponent;

}


int AliHLTLumiRegComponent::DoInit(int argc, const char **argv){

  int iResults=0;
  TString strVtx = "";
  
  for(Int_t iHisto=0; iHisto<3; iHisto++){
    
    if (iHisto == 0) {
      strVtx = "X";
      fPrimaryLHC[iHisto] = new TH1F(Form("h%sTRKVtxLHC", strVtx.Data()),Form("h%sTRKVtxLHC", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDQM[iHisto] = new TH1F(Form("h%sTRKVtxDQM", strVtx.Data()),Form("h%sTRKVtxDQM", strVtx.Data()), 200,-0.7,0.7);
      
      fPrimaryDefMultLHC[iHisto] = new TH1F (Form("h%sTRKDefMultLHC", strVtx.Data()),Form("h%sTRKDefMultLHC", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDefMultDQM[iHisto] = new TH1F (Form("h%sTRKDefMultDQM", strVtx.Data()),Form("h%sTRKDefMultDQM", strVtx.Data()), 200,-0.7,0.7);
    }
    
    if (iHisto == 1) {
      strVtx = "Y";
      fPrimaryLHC[iHisto] = new TH1F(Form("h%sTRKVtxLHC", strVtx.Data()),Form("h%sTRKVtxLHC", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDQM[iHisto] = new TH1F(Form("h%sTRKVtxDQM", strVtx.Data()),Form("h%sTRKVtxDQM", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDefMultLHC[iHisto] = new TH1F (Form("h%sTRKDefMultLHC", strVtx.Data()),Form("h%sTRKDefMultLHC", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDefMultDQM[iHisto] = new TH1F (Form("h%sTRKDefMultDQM", strVtx.Data()),Form("h%sTRKDefMultDQM", strVtx.Data()), 200,-0.7,0.7);
    }
      
    if (iHisto == 2) {
      strVtx = "Z";
      fPrimaryLHC[iHisto] = new TH1F(Form("h%sTRKVtxLHC", strVtx.Data()),Form("h%sTRKVtxLHC", strVtx.Data()), 200,-20,20);
      fPrimaryDQM[iHisto] = new TH1F(Form("h%sTRKVtxDQM", strVtx.Data()),Form("h%sTRKVtxDQM", strVtx.Data()), 200,-20,20);
    }
  }
  
  fPushBackPeriodLHC = 60;
  fPushBackPeriodDQM = 300;
  fLastPushBackTime = 0;
  
  AliCDBEntry *grpEntry = AliCDBManager::Instance()->Get("/GRP/GRP/Data/");
  if(!grpEntry) {
    ::Error("AliHLTLumiRegComponent","Cannot get AliCDBEntry");
    return 0;
  }
  const AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(grpEntry->GetObject());
  if(!grpData) {
    ::Error("AliHLTLumiRegComponent","Cannot get AliGRPObject");
    return 0;
  }
  TString beamType(grpData->GetBeamType());
  
  if(beamType == "A-A") fEventSpecie = AliHLTLumiRegComponent::kPbPb;
  else
    if((beamType == "p-A")||(beamType == "A-p")) fEventSpecie = AliHLTLumiRegComponent::kpPb;
    else fEventSpecie = AliHLTLumiRegComponent::kpp;

  return 1;
}

Int_t AliHLTLumiRegComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){

  AliSysInfo::AddStamp("AliHLTLumiRegComponent::DoEvent.Start");
  Int_t vertexITSSATRKok=0;
  AliESDVertex vertexTracks;
  AliFlatESDVertex *vtxFlat;
  
  const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDVertex|kAliHLTDataOriginITS);
  if (pBlock) {
   // fpBenchmark->AddInput(pBlock->fSize);
    vtxFlat =  reinterpret_cast<AliFlatESDVertex*>( pBlock->fPtr );
    if (!vtxFlat) {
      return 0;
    }
    if (vtxFlat->GetNContributors()>0)
      vertexITSSATRKok=1;
    
    if (!vertexITSSATRKok) {
      HLTInfo("Vertex Tracks not found, trying the SPD one...");
      return 0;
      //SPD vertex!!
    }
  }
  
  Double_t vtxPos[3] = {0.,0.,0.};
  Int_t nContrib = vtxFlat->GetNContributors();
  vtxFlat->GetXYZ(vtxPos);
  
  for (Int_t iCoord =0; iCoord<3; iCoord++){
    fPrimaryLHC[iCoord]->Fill(vtxPos[iCoord]);
    fPrimaryDQM[iCoord]->Fill(vtxPos[iCoord]);
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kpp){
    if ((nContrib>30) && (nContrib<50)){
      for (Int_t iCoord =0; iCoord<2; iCoord++){
        fPrimaryDefMultLHC[iCoord]->Fill(vtxPos[iCoord]);
        fPrimaryDefMultDQM[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kpPb) {
    if (nContrib>100){
      for (Int_t iCoord =0; iCoord<2; iCoord++){
        fPrimaryDefMultLHC[iCoord]->Fill(vtxPos[iCoord]);
        fPrimaryDefMultDQM[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kPbPb){
    if (nContrib>500){
      for (Int_t iCoord =0; iCoord<2; iCoord++){
        fPrimaryDefMultLHC[iCoord]->Fill(vtxPos[iCoord]);
        fPrimaryDefMultDQM[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }
  
  if(fPushBackPeriodLHC>0) {
    TDatime time;
    if ((fLastPushBackTime<0) || ((int)time.Get()-fLastPushBackTime < fPushBackPeriodLHC)) {
      return 0; 
    }
    else {
      fLastPushBackTime = (int)time.Get();
      Float_t meanVtx[3] ={0.,0.,0.};
      Float_t sigmaVtx[3] ={0.,0.,0.};
      Float_t meanLR[3] ={0.,0.,0.};
      Float_t sigmaLR[3] ={0.,0.,0.};
      Int_t fitVtxResults = FitPositions(fPrimaryLHC, meanVtx, sigmaVtx);
      Int_t lumiRegResults = LuminousRegionExtraction(fPrimaryDefMultLHC, meanLR, sigmaLR);
      
      if (lumiRegResults == 1) {
        HLTWarning("Luminous Region determination ok");
        for (Int_t iCoord =0; iCoord<3; iCoord++){
          fPrimaryLHC[iCoord]->Reset();
          if ((iCoord ==0)||(iCoord==1)) fPrimaryDefMultLHC[iCoord]->Reset();
        }
      }
      
      if (lumiRegResults == 2) HLTWarning("Problems in the luminous region fit, using unconvoluted sigma");
       if (lumiRegResults == 0) HLTWarning("Problems in the luminous region fit, returning 0");
    }
  }
  
  if(fPushBackPeriodDQM>0) {
    TDatime time;
    if ((fLastPushBackTime<0) || ((int)time.Get()-fLastPushBackTime < fPushBackPeriodDQM)) {
      return 0; 
    }
    else {
      fLastPushBackTime = (int)time.Get();
      Float_t meanVtx[3] ={0.,0.,0.};
      Float_t sigmaVtx[3] ={0.,0.,0.};
      Float_t meanLR[3] ={0.,0.,0.};
      Float_t sigmaLR[3] ={0.,0.,0.};
      Int_t fitVtxResults = FitPositions(fPrimaryDQM, meanVtx, sigmaVtx);
      Int_t lumiRegResults = LuminousRegionExtraction(fPrimaryDefMultDQM, meanLR,sigmaLR);
      
      if (lumiRegResults == 1) {
        HLTWarning("Luminous Region determination ok");
        for (Int_t iCoord =0; iCoord<3; iCoord++){
          fPrimaryDQM[iCoord]->Reset();
          if ((iCoord ==0)||(iCoord==1)) fPrimaryDefMultDQM[iCoord]->Reset();
        }
      }
      
      if (lumiRegResults == 2) HLTWarning("Problems in the luminous region fit, using unconvoluted sigma");
      if (lumiRegResults == 0) HLTWarning("Problems in the luminous region fit, returning 0");
      
    }
  }

  PushBack( (TH1F*) fPrimaryLHC[0], kAliHLTDataTypeHistogram);
  // push results to LHC interface or DQM
  return 0;
  
}


Int_t AliHLTLumiRegComponent::FitPositions(TH1F *histos[], Float_t* mean, Float_t* sigma){
  
  Int_t fitVtxResults1 = 0, fitVtxResults2 = 0,fitVtxResults3 = 0;
  fitVtxResults1 = FitHistos(histos[0], mean[0], sigma[0], -1, 1);
  fitVtxResults2 = FitHistos(histos[1], mean[1], sigma[1], -1, 1);
  fitVtxResults3 = FitHistos(histos[2], mean[2], sigma[2], -12, 12);
  
  return fitVtxResults1||fitVtxResults2||fitVtxResults3;
  
}


Int_t AliHLTLumiRegComponent::LuminousRegionExtraction(TH1F* histos[], Float_t* meanLR, Float_t* sigmaLR){
  
  Float_t meanMult = 38.;
  Float_t p2 = 1.3;
  Float_t resolVtx = 0.04;
  Float_t lumiregsquared = 0;
  
  FitHistos(histos[0], meanLR[0], sigmaLR[0], -0.4, 0.4);
  FitHistos(histos[1], meanLR[1], sigmaLR[1], -0.2, 0.7);
  
  if (AliHLTLumiRegComponent::kpp){
    lumiregsquared = (sigmaLR[0]*sigmaLR[0] - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
    
    if (lumiregsquared < 0 || lumiregsquared < 1E-5){
      HLTWarning(Form("Difficult luminous region determination X, keep convoluted sigma"));
      return 2;
    }
    
    if (lumiregsquared > 0 && lumiregsquared < 0.0005) {
      sigmaLR[0] = TMath::Sqrt(lumiregsquared);
      return 1;
    }
    
    lumiregsquared = (sigmaLR[1]*sigmaLR[1] - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
    
    if (lumiregsquared < 0 || lumiregsquared < 1E-5){
      HLTWarning(Form("Difficult luminous region determination Y, keep convoluted sigma"));
      return 2;
    }
    
    if (lumiregsquared > 0 && lumiregsquared < 0.0005) {
      sigmaLR[1] = TMath::Sqrt(lumiregsquared);
      return 1;
    }
  }
  
  return 2;
  
}

Int_t AliHLTLumiRegComponent::FitHistos(TH1F *hVtx, Float_t &mean, Float_t &sigma, Float_t rangelow, Float_t rangeup){
  
  Int_t isFitok = 1;
  Int_t nEntries =  hVtx->GetEntries();
  if (nEntries < 50) {
    HLTWarning("too few entries for fitting");
    //do something!
  }
  
  hVtx->Fit("gaus", "M", "", rangelow, rangeup);
  TF1 *vtxFunct = hVtx->GetFunction("gaus");
  if (!vtxFunct) {
    HLTError("No fit function!"); 
    return 0;
  }
  mean = vtxFunct->GetParameter(1);
  sigma = vtxFunct->GetParameter(2);
  
  delete vtxFunct;
  return isFitok;
}

