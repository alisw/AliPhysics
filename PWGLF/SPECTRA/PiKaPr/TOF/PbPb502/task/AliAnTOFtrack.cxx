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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                                                                       //
// Class container for identified charged hadron spectra: TOF            //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


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

Float_t AliAnTOFtrack::GetDeltaT(const UInt_t id){
  if(id > kExpSpecies) AliFatal("Index required is out of bound");
  return fTOFTime - fTOFExpTime[id] - fT0TrkTime;
}

Float_t AliAnTOFtrack::GetDeltaSigma(const UInt_t id, const UInt_t hypo) {
  if(id > kExpSpecies || hypo > kExpSpecies) AliFatal("Index required is out of bound");
  return GetDeltaT(hypo)/fTOFExpSigma[id];
}
