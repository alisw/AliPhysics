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

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Class container for identified charged hadron spectra: TOF            //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

#define LOG_NO_INFO
#define LOG_NO_DEBUG
#include "AliLog.h"
#include "TObject.h"
#include "AliAnTOFtrack.h"
#include <iostream>



ClassImp(AliAnTOFtrack)

//________________________________________________________________________
AliAnTOFtrack::AliAnTOFtrack() :
fTrkMask(0), 
fTPCPIDMask(0),
fTrkCutMask(0),
fDCAXYIndex(0), 
fDCAZIndex(0), 
fLength(-999), 
fLengthRatio(-999), 
fTOFTime(-999), 
fTOFMismatchTime(-999), 
fT0TrkTime(-999), 
fTOFchan(-999), 
fEta(-999), 
fPhi(-999),
fPt(-999), 
fNTOFClusters(-999),
fTPCSignal(-999) 
// fT0TrkSigma(-999)
{
  
  //
  // standard constructur which should be used
  //
  AliInfo("**** CONSTRUCTOR CALLED ****");
  for(Int_t i = 0; i < kExpSpecies; i++){
    fTOFExpTime[i] = -999;
    fTOFExpSigma[i] = -999;
  }
  
  AliInfo("**** END OF CONSTRUCTOR ****");
}

//________________________________________________________________________
AliAnTOFtrack::~AliAnTOFtrack(){//Destructor
  AliInfo("**** DESTRUCTOR CALLED ****");
  
  
  AliInfo("**** END OF DESTRUCTOR ****");
  
}

///////////////////
///Cut functions///
///////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::PassStdCut(){
  for(Int_t i = 0; i < nCuts; i++) if(!GetMaskBit(fTrkCutMask, CutStdIndexInMask[i])) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::PassCut(const Int_t cut){
  if(cut < -1 || cut >= nCutVars) AliFatal("requested cut is out of bound");
  //Always apply standard cuts except if requiring one different cut
  if(cut == -1) return PassStdCut();;
  
  Int_t cutindex = -1;  
  Int_t sum = 0;
  for(Int_t i = 0; i < nCuts; i++){//Identify if the cut Type is of the Loose kind
    if(cut == sum){
      cutindex = i;
      break;
    }
    sum += CutIndex[i]-1;//One cut is not necessary
  }
  
  //Finally check the cuts
  for(Int_t i = 0; i < nCuts; i++){
    if(i == cutindex) continue;
    else if(!GetMaskBit(fTrkCutMask, CutStdIndexInMask[i])) return kFALSE;
  }
  if(cutindex == -1) if(!GetMaskBit(fTrkCutMask, cut)) return kFALSE; //Cut type must be of the Tight kind
  
  return kTRUE;
  
}

/////////////////////////
//////PID functions//////
/////////////////////////

//________________________________________________________________________
Float_t AliAnTOFtrack::GetDeltaT(const UInt_t id){
  if(id > kExpSpecies) AliFatal("Index required is out of bound");
  return fTOFTime - fTOFExpTime[id] - fT0TrkTime;
}

//________________________________________________________________________
Float_t AliAnTOFtrack::GetDeltaSigma(const UInt_t id, const UInt_t hypo) {
  if(id > kExpSpecies || hypo > kExpSpecies) AliFatal("Index required is out of bound");
  return GetDeltaT(hypo)/fTOFExpSigma[id];
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCElectron(){
  if(GetMaskBit(fTPCPIDMask, kIsTPCElectron)) return kTRUE; //1.5 sigma cut for Electrons in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCPiKP(const UInt_t i){
  if(i >= 3) AliFatal("Wrong index required");
  if(GetMaskBit(fTPCPIDMask, kIsTPCPion + i)) return kTRUE; //5 sigma cut for Pi/K/P in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCPiKP(){
  for(Int_t i = 0; i < 3; i++) if(IsTPCPiKP(i)) return kTRUE; //5 sigma cut for Pi/K/P in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::ConsistentTPCTOF(){
  if(!IsTPCPiKP()) return kFALSE;
  for (Int_t i = 0; i < 3; i++) {
    if(GetMaskBit(fTPCPIDMask, kIsTPCPion + i) && TMath::Abs(GetDeltaSigma(kpi + i, kpi + i)) < 5) return kTRUE;
  }
  return kFALSE;
}

///////////////////////
////Track functions////
///////////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsNegative(){
  if(GetMaskBit(fTrkMask, kNegTrk)) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Double_t AliAnTOFtrack::GetDCA(const Bool_t xy){
  Double_t u, d;
  GetBinnedDCA(d, u, xy);
  return d + (u-d)/2;
}

////////////////////////
/////TEST functions/////
////////////////////////

//________________________________________________________________________
void AliAnTOFtrack::TestDCAXYBinning(){
  Int_t dim = static_cast<Int_t>(TMath::Power(2., 8.*sizeof(fDCAXYIndex)));
  Double_t range = (2.*fDCAXYRange);
  std::cout<<"fDCAXYIndex has "<<dim<<" possibilities"<<std::endl;
  std::cout<<"DCAXY range is ["<<-fDCAXYRange<<","<<fDCAXYRange<<"] -> "<<range<<" bin width = "<<range/dim<<std::endl;
  std::cout<<"Variable dimension is "<<dim<<" while actual bin required are "<<kDCAXYBins<<std::endl;
  
  dim = static_cast<Int_t>(TMath::Power(2., 8.*sizeof(fDCAZIndex)));
  range = (2.*fDCAZRange);
  std::cout<<"fDCAZIndex has "<<dim<<" possibilities"<<std::endl;
  std::cout<<"DCAZ range is ["<<-fDCAZRange<<","<<fDCAZRange<<"] -> "<<range<<" bin width = "<<range/dim<<std::endl;
  std::cout<<"Variable dimension is "<<dim<<" while actual bin required are "<<kDCAZBins<<std::endl;
  
  for (Int_t var = 0; var <=20; ++var) {
    //       Double_t rnd = 2.*fDCAXYRange*gRandom->Rndm()-1;
    Double_t rnd = 2.*fDCAXYRange*(Double_t)var/20.-fDCAXYRange;
    ComputeDCABin(rnd, kTRUE);
    
    Double_t binlow, binup;
    GetBinnedDCA(binlow, binup, 1);
    
    std::cout<<rnd<<" Index: "<<fDCAXYIndex<<"   ["<<binlow<<" ; "<<binup<<"] --> Diff --> ["<< binlow - rnd <<" ; "<<binup - rnd<<"]"<<std::endl;
    //       cout<<rnd<<"  "<<TMath::Abs((rnd+fDCAXYRange)*kDCAXYBins/(2*fDCAXYRange))<<" Index: "<<fDCAXYIndex<<"   ["<<binlow<<" ; "<<binup<<"] -->  ["<< binlow - rnd <<" ; "<<binup - rnd<<"]"<<endl;
  }
  
  //     for (Int_t c = 0 ; c <= kDCAXYBins; ++c) cout<<c<<"  "<<-fDCAXYRange+2.*fDCAXYRange*c/kDCAXYBins<<endl;
  
}

///////////////////////////
///DCA Utility functions///
///////////////////////////

//________________________________________________________________________
void AliAnTOFtrack::ComputeDCABin(const Double_t dca, const Bool_t xy){
  if(xy){
    if(dca > fDCAXYRange){//Overflow
      fDCAXYIndex = kDCAXYBins+2;
      return;
    }
    else if(dca <= -fDCAXYRange){//Underflow
      fDCAXYIndex = kDCAXYBins+1;
      return;
    }
    fDCAXYIndex = static_cast<UShort_t>(TMath::Ceil((dca+fDCAXYRange)/((2*fDCAXYRange)/kDCAXYBins)) -1);
  }
  else{
    if(dca > fDCAZRange){//Overflow
      fDCAZIndex = kDCAZBins+2;
      return;
    }
    else if(dca <= -fDCAZRange){//Underflow
      fDCAZIndex = kDCAZBins+1;
      return;
    }
    fDCAZIndex = static_cast<UShort_t>(TMath::Ceil((dca+fDCAZRange)/((2*fDCAZRange)/kDCAZBins)) -1);
  }
  
};

//________________________________________________________________________
void AliAnTOFtrack::ComputeDCABin(const Double_t dcaxy, const Double_t dcaz){
  
  ComputeDCABin(dcaxy, (Bool_t) kTRUE);
  ComputeDCABin(dcaz, (Bool_t) kFALSE);
};


//////////////////////////
///T0 Utility functions///
//////////////////////////

//________________________________________________________________________
Double_t AliAnTOFtrack::GetT0Resolution(const Double_t TOFsigma) const {
  return TMath::Sqrt(TMath::Power(fTOFExpSigma[0], 2) - TMath::Power(TOFsigma, 2) - TMath::Power(15./GetMomentum(), 2));
}

////////////////////////
///T0 Usage functions///
////////////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF(const Bool_t exclusive){
  if(GetMaskBit(fTrkMask, kT0_0)){
    if(exclusive && !GetMaskBit(fTrkMask, kT0_1) && !GetMaskBit(fTrkMask, kT0_2)) return kTRUE;
    else if(exclusive) return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0A(const Bool_t exclusive){
  if(GetMaskBit(fTrkMask, kT0_1)){
    if(exclusive && !GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_2)) return kTRUE;
    else if(exclusive) return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0C(const Bool_t exclusive){
  if(GetMaskBit(fTrkMask, kT0_2)){
    if(exclusive && !GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_1)) return kTRUE;
    else if(exclusive) return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0Fill(){
  if(!GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_1) && !GetMaskBit(fTrkMask, kT0_2)) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0A(){
  if(IsT0TOF() && IsT0A()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0A(){
  if(IsT0TOF() || IsT0A()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0C(){
  if(IsT0TOF() && IsT0C()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0C(){
  if(IsT0TOF() || IsT0C()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0A_T0C(){
  if(IsT0A() && IsT0C()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0A_T0C(){
  if(IsT0A() || IsT0C()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0A_T0C(){
  if(IsT0TOF() && IsT0A() && IsT0C()) return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0A_T0C(){
  if(IsT0TOF() || IsT0A() || IsT0C()) return kTRUE;
  return kFALSE;
}
