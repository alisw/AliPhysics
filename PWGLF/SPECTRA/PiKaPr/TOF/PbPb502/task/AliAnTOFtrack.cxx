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
#include "Riostream.h"



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
