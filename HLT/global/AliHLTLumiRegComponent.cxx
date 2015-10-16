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
AliHLTLumiRegComponent::AliHLTLumiRegComponent():
AliHLTProcessor(),
fEventSpecie(0)
{
  for(Int_t i=0; i<3; i++){
    fPrimary[i] = 0x0;
  }
  
  for(Int_t i=0; i<2; i++){
    fPrimaryDefMult[i] = 0x0;
  }
  
}

AliHLTLumiRegComponent::AliHLTLumiRegComponent( const AliHLTLumiRegComponent& ):
AliHLTProcessor(),
fEventSpecie(0)
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
    delete fPrimary[i];
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
      fPrimary[iHisto] = new TH1F(Form("h%sTRKVtx", strVtx.Data()),Form("h%sTRKVtx", strVtx.Data()), 200,-0.7,0.7);
      
      fPrimaryDefMult[iHisto] = new TH1F (Form("h%sTRKDefMult", strVtx.Data()),Form("h%sTRKDefMult", strVtx.Data()), 200,-0.7,0.7);
    }
    
    if (iHisto == 1) {
      strVtx = "Y";
      fPrimary[iHisto] = new TH1F(Form("h%sTRKVtx", strVtx.Data()),Form("h%sTRKVtx", strVtx.Data()), 200,-0.7,0.7);
      fPrimaryDefMult[iHisto] = new TH1F (Form("h%sTRKDefMult", strVtx.Data()),Form("h%sTRKDefMult", strVtx.Data()), 200,-0.7,0.7);
    }
      
    if (iHisto == 2) {
      strVtx = "Z";
      fPrimary[iHisto] = new TH1F(Form("h%sTRKVtx", strVtx.Data()),Form("h%sTRKVtx", strVtx.Data()), 200,-20,20);
      fPrimaryDefMult[iHisto] = new TH1F (Form("h%sTRKDefMult", strVtx.Data()),Form("h%sTRKDefMult", strVtx.Data()), 200,-20.,20.);
    }
  }
  
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

  if (!IsDataEvent())
  {
    //on EOR push unconditionally
    const AliHLTComponentBlockData* pBlock =
      GetFirstInputBlock(kAliHLTDataTypeEOR|kAliHLTDataOriginAny);
    if (pBlock)
    {
      PushAndReset();
    }
    return 0;
  }

  const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDVertex|kAliHLTDataOriginITS);
  if (!pBlock) return 0;
  
  // fpBenchmark->AddInput(pBlock->fSize);
  vtxFlat =  reinterpret_cast<AliFlatESDVertex*>( pBlock->fPtr );
  if (!vtxFlat) return 0;
  int nContrib = vtxFlat->GetNContributors();
  if (nContrib>0)
    vertexITSSATRKok=1;

  if (!vertexITSSATRKok) {
    HLTInfo("Vertex Tracks not found, trying the SPD one...");
    return 0;
    //SPD vertex!!
  }
  
  Double_t vtxPos[3] = {0.,0.,0.};
  vtxFlat->GetXYZ(vtxPos);
  
  for (Int_t iCoord =0; iCoord<3; iCoord++){
    fPrimary[iCoord]->Fill(vtxPos[iCoord]);
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kpp){
    if ((nContrib>30) && (nContrib<50)){
      for (Int_t iCoord =0; iCoord<3; iCoord++){
        fPrimaryDefMult[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kpPb) {
    if (nContrib>100){
      for (Int_t iCoord =0; iCoord<3; iCoord++){
        fPrimaryDefMult[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }
  
  if (fEventSpecie == AliHLTLumiRegComponent::kPbPb){
    if (nContrib>500){
      for (Int_t iCoord =0; iCoord<3; iCoord++){
        fPrimaryDefMult[iCoord]->Fill(vtxPos[iCoord]);
      }
    }
  }

  PushAndReset();

  return 0;
}

int AliHLTLumiRegComponent::PushAndReset()
{
  HLTInfo("pushing and resetting");
  for (Int_t iCoord =0; iCoord<3; iCoord++){
    if (PushBack(fPrimary[iCoord], kAliHLTDataTypeHistogram | kAliHLTDataOriginAny) > 0)
      fPrimary[iCoord]->Reset();
    
    if (PushBack(fPrimaryDefMult[iCoord], kAliHLTDataTypeHistogram | kAliHLTDataOriginAny) > 0)
      fPrimaryDefMult[iCoord]->Reset();
  }
  return 0;
}

Int_t AliHLTLumiRegComponent::FitPositions(TH1F *histos[], Float_t* mean, Float_t* sigma){
  
  Int_t fitVtxResults1 = 0, fitVtxResults2 = 0,fitVtxResults3 = 0;
  fitVtxResults1 = AliHLTLumiRegComponent::FitHistos(histos[0], mean[0], sigma[0], -1, 1);
  fitVtxResults2 = AliHLTLumiRegComponent::FitHistos(histos[1], mean[1], sigma[1], -1, 1);
  fitVtxResults3 = AliHLTLumiRegComponent::FitHistos(histos[2], mean[2], sigma[2], -12, 12);
  
  return fitVtxResults1||fitVtxResults2||fitVtxResults3;
  
}


Int_t AliHLTLumiRegComponent::LuminousRegionExtraction(TH1F* histos[], Float_t* meanLR, Float_t* sigmaLR){
  
  Float_t meanMult = 38.;
  Float_t p2 = 1.3;
  Float_t resolVtx = 0.04;
  Float_t lumiregsquaredX = 0;
  Float_t lumiregsquaredY = 0;
  int rc = 0;
  
  AliHLTLumiRegComponent::FitHistos(histos[0], meanLR[0], sigmaLR[0], -0.4, 0.4);
  AliHLTLumiRegComponent::FitHistos(histos[1], meanLR[1], sigmaLR[1], -0.2, 0.7);
  
  if (AliHLTLumiRegComponent::kpp){
    lumiregsquaredX = (sigmaLR[0]*sigmaLR[0] - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
    lumiregsquaredY = (sigmaLR[1]*sigmaLR[1] - ((resolVtx*resolVtx)/TMath::Power(meanMult, p2)));
    
    if (lumiregsquaredX > 0 && lumiregsquaredX < 0.0005) {
      sigmaLR[0] = TMath::Sqrt(lumiregsquaredX);
    }
    
    if (lumiregsquaredY > 0 && lumiregsquaredY < 0.0005) {
      sigmaLR[1] = TMath::Sqrt(lumiregsquaredY);
    }
    
    if (lumiregsquaredX < 0 || lumiregsquaredX < 1E-5){
      Printf("Difficult luminous region determination X, keep convoluted sigma");
      rc = 1;
    }
    
    if (lumiregsquaredY < 0 || lumiregsquaredY < 1E-5){
      Printf("Difficult luminous region determination Y, keep convoluted sigma");
      rc += 2;
    }
    
  }
  
  //rc: 0 on success, 1 on bad X fit, 2 on bad Y fit, 3 on both bad
  return rc;
  
}

Int_t AliHLTLumiRegComponent::FitHistos(TH1F *hVtx, Float_t &mean, Float_t &sigma, Float_t rangelow, Float_t rangeup){
  
  Int_t isFitok = 1;
  Int_t nEntries =  hVtx->GetEntries();
  if (nEntries < 50) {
    Printf("too few entries for fitting");
    //do something!
  }
  
  hVtx->Fit("gaus", "M", "", rangelow, rangeup);
  TF1 *vtxFunct = hVtx->GetFunction("gaus");
  if (!vtxFunct) {
    Printf("No fit function!"); 
    return 1;
  }
  mean = vtxFunct->GetParameter(1);
  sigma = vtxFunct->GetParameter(2);
  
  delete vtxFunct;
  return isFitok;
}

