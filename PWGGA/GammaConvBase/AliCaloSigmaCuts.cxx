/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Authors: Steven Merkel, Adrian mechler                                  *
* Version 1.0                                                             *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling all kinds of selection cuts for
// Sigma+ analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliCaloSigmaCuts.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TPDGCode.h"
#include "TDatabasePDG.h"
#include "AliAODMCParticle.h"


class iostream;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliCaloSigmaCuts)
/// \endcond


const char* AliCaloSigmaCuts::fgkCutNames[AliCaloSigmaCuts::kNCuts] = {
  "FilterBitCut", //0
  "NClusterTPCCut", //1
  "Chi2TPCCut", //2
  "NClusterITSCut", //3
  "Chi2ITSCut", //4
  "DCAXYCut", //5
  "DCAZCut", //6
  "NSigmaTPCCut", //7
  "NSigmaTOFCut", //8
  "PionMassLowerCut", //9
  "PionMassUpperCut", //10
  "PodolanskiCut", //11
  "OpeningAngleCut", //12
  "BackgroundEstimation", //13
};

//________________________________________________________________________
AliCaloSigmaCuts::AliCaloSigmaCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fCutString(NULL),
  fCutStringRead(""),
  fAmenterosCut(NULL),
  fFilterBit(0),
  fNClusterTPC(0),
  fChi2TPC(0),
  fNClusterITS(0),
  fChi2ITS(0),
  fDCAXY(0),
  fDCAZ(0),
  fNSigmaTPC(0),
  fNSigmaTOF(0),
  fMaxPionMass(0),
  fMinPionMass(0),
  fMaxAlpha(0),
  fMinAlpha(0),
  fQt(0),
  fMaxOpeningAngle(0),
  fMinOpeningAngle(0),
  fBackgroundestimation(0)
  
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//________________________________________________________________________
AliCaloSigmaCuts::AliCaloSigmaCuts(const AliCaloSigmaCuts &ref) :
  AliAnalysisCuts(ref),
  fCutString(NULL),
  fCutStringRead(""),
  fAmenterosCut(NULL),
  fFilterBit(ref. fFilterBit),
  fNClusterTPC(ref. fNClusterTPC),
  fChi2TPC(ref. fChi2TPC),
  fNClusterITS(ref. fNClusterITS),
  fChi2ITS(ref. fChi2ITS),
  fDCAXY(ref. fDCAXY),
  fDCAZ(ref. fDCAZ),
  fNSigmaTPC(ref. fNSigmaTPC),
  fNSigmaTOF(ref. fNSigmaTOF),
  fMaxPionMass(ref.fMaxPionMass),
  fMinPionMass(ref.fMinPionMass),
  fMaxAlpha(ref.fMaxAlpha),
  fMinAlpha(ref.fMinAlpha),
  fQt(ref.fQt),
  fMaxOpeningAngle(ref.fMaxOpeningAngle),
  fMinOpeningAngle(ref.fMinOpeningAngle),
  fBackgroundestimation(ref.fBackgroundestimation)

{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//________________________________________________________________________
AliCaloSigmaCuts::~AliCaloSigmaCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }
  if(fAmenterosCut != NULL){
    delete fAmenterosCut;
    fAmenterosCut = NULL;
  }
}

//________________________________________________________________________
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloSigmaCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  // AliInfo(Form("Set Meson Cutnumber: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  //cout << "Updating cut  " << fgkCutNames[cutID] << " (" << cutID << ") to " << value << endl;
  switch (cutID) {
  case kFilterBit:
    if( SetFilterBitCut(value)) {
      fCuts[kFilterBit] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kNTPCCluster:
    if( SetNClusterTPCCut(value)) {
      fCuts[kNTPCCluster] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kChi2TPC:
    if( SetChi2TPCCut(value)) {
      fCuts[kChi2TPC] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kNITSCluster:
    if( SetNClusterITSCut(value)) {
      fCuts[kNITSCluster] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kChi2ITS:
    if( SetChi2ITSCut(value)) {
      fCuts[kChi2ITS] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kMinDCAXY:
    if( SetDCAXYCut(value)) {
      fCuts[kMinDCAXY] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kMinDCAZ:
    if( SetDCAZCut(value)) {
      fCuts[kMinDCAZ] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kNSigmaTPC:
    if( SetNSigmaTPCCut(value)) {
      fCuts[kNSigmaTPC] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kNSigmaTOF:
    if( SetNSigmaTOFCut(value)) {
      fCuts[kNSigmaTOF] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;                  
  case kPionMassLower:
    if( SetMinPionMassCut(value)) {
      fCuts[kPionMassLower] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kPionMassUpper:
    if( SetMaxPionMassCut(value)) {
      fCuts[kPionMassUpper] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kAmenterosCut:
    if( SetAmenterosCut(value)) {
      fCuts[kAmenterosCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kOpeningAngleCut:
    if( SetOpeningAngleCut(value)) {
      fCuts[kOpeningAngleCut] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kBackgroundEstimation:
    if( SetBackgroundEstimation(value)) {
      fCuts[kBackgroundEstimation] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
  case kNCuts:
    cout << "Error:: Cut id out of range"<< endl;
    return kFALSE;
  }

  cout << "Error:: Cut id " << cutID << " not recognized "<< endl;
  return kFALSE;

}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetFilterBitCut(Int_t FilterBitCut){
  // Set Cut
  switch(FilterBitCut){
  case 0:
    fFilterBit = 1<<0;  //TPC Only. 1<<0
    break;
  case 1:
    fFilterBit = 1<<4;  // TPC+ITS 1<<4
    break;  
  case 2:
    fFilterBit = 1<<8;   //global hybrids
    break;    
  default:
    cout<<"Warning: FilterBitCut not defined"<<FilterBitCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetNClusterTPCCut(Int_t NClusterTPCCut){
  // Set Cut
  switch(NClusterTPCCut){
  case 0:
    fNClusterTPC = 0;
    break;
  case 1:
    fNClusterTPC = 80;
    break;  
  default:
    cout<<"Warning: NClusterTPCCut not defined"<<NClusterTPCCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetChi2TPCCut(Int_t Chi2TPCCut){
  // Set Cut
  switch(Chi2TPCCut){
  case 0:
    fChi2TPC = 10.;
    break;
  case 1:
    fChi2TPC = 4.;
    break;  
  default:
    cout<<"Warning: Chi2TPCCut not defined"<<Chi2TPCCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetNClusterITSCut(Int_t NClusterITSCut){
  // Set Cut
  switch(NClusterITSCut){
  case 0:
    fNClusterITS = 0;
    break;
  case 1:
    fNClusterITS = 2;
    break;  
  default:
    cout<<"Warning: NClusterITSCut not defined"<<NClusterITSCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetChi2ITSCut(Int_t Chi2ITSCut){
  // Set Cut
  switch(Chi2ITSCut){
  case 0:
    fChi2ITS = 10.;
    break;
  case 1:
    fChi2ITS = 3.;
    break;  
  default:
    cout<<"Warning: Chi2ITSCut not defined"<<Chi2ITSCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetDCAXYCut(Int_t DCAXYCut){
  // Set Cut
  switch(DCAXYCut){
  case 0:
    fDCAXY = 0.;
    break;
  case 1:
    fDCAXY = 0.05;
    break;  
  default:
    cout<<"Warning: DCAXYCut not defined"<<DCAXYCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetDCAZCut(Int_t DCAZCut){
  // Set Cut
  switch(DCAZCut){
  case 0:
    fDCAZ = 0.;
    break;
  case 1:
    fDCAZ = 0.05;
    break;  
  default:
    cout<<"Warning: DCAZCut not defined"<<DCAZCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetNSigmaTPCCut(Int_t NSigmaTPCCut){
  // Set Cut
  switch(NSigmaTPCCut){
  case 0:
    fNSigmaTPC = 5.;
    break;
  case 1:
    fNSigmaTPC = 3.;
    break;  
  default:
    cout<<"Warning: NSigmaTPCCut not defined"<<NSigmaTPCCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetNSigmaTOFCut(Int_t NSigmaTOFCut){
  // Set Cut
  switch(NSigmaTOFCut){
  case 0:
    fNSigmaTOF = 5.;
    break;
  case 1:
    fNSigmaTOF = 3.;
    break;  
  default:
    cout<<"Warning: NSigmaTOFCut not defined"<<NSigmaTOFCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetMinPionMassCut(Int_t PionMinMassCut){
  // Set Cut
  switch(PionMinMassCut){
  case 0:
    fMinPionMass = 0.;
    break;
  case 1:
    fMinPionMass = 0.118;
    break;  
  default:
    cout<<"Warning: PionMinMassCut not defined"<<PionMinMassCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetMaxPionMassCut(Int_t PionMaxMassCut){
  // Set Cut
  switch(PionMaxMassCut){
  case 0:
    fMaxPionMass = 0.5;
    break;
  case 1:
    fMaxPionMass = 0.148;
    break;  
  default:
    cout<<"Warning: PionMaxMassCut not defined"<<PionMaxMassCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetAmenterosCut(Int_t AmenterosCut){
  // Set Cut
  switch(AmenterosCut){
  case 0:
    fMaxAlpha = 10.;
    fMinAlpha = 0.;
    fQt = 10.;
    break;
  case 1:
    fMaxAlpha = 1.02;
    fMinAlpha = 0.2;
    fQt = 0.3;
    break;
  default:
    cout<<"Warning: AmenterosCut not defined"<<AmenterosCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetOpeningAngleCut(Int_t OpeningAngleCut){
  // Set Cut
  switch(OpeningAngleCut){
  case 0:
    fMaxOpeningAngle = 4.;
    fMinOpeningAngle = 0.;
    break;
  case 1:
    fMaxOpeningAngle = 0.35;
    fMinOpeningAngle = 0.05;
    break;
  default:
    cout<<"Warning: OpeningAngleCut not defined"<<OpeningAngleCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetBackgroundEstimation(Int_t BackgroundEstimation){
  // Set Cut
  switch(BackgroundEstimation){
  case 0:
    fBackgroundestimation = 0;
    break;
  case 1:
    fBackgroundestimation = 1;
    break;  

  default:
    cout<<"Warning: BackgroundEstimation not defined"<<BackgroundEstimation<<endl;
    return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
TString AliCaloSigmaCuts::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}


//_____________________________________________________________________________
// Pion Mass Cut on sigma daughters
Bool_t AliCaloSigmaCuts::PionIsSelectedByMassCut(Double_t pionMass){
  if(fMaxPionMass){
    if(fMinPionMass){
      if(fMaxPionMass > pionMass && fMinPionMass < pionMass)
        return kTRUE;
      else
        return kFALSE;
    }    
  }
  else{
    return kTRUE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
// Pion Mass Cut on sigma daughters
Bool_t AliCaloSigmaCuts::ArmenterosLikeQtCut(Double_t alpha, Double_t qT){
  if(fMinAlpha){
    if(fMaxAlpha){
      if(fQt){
        if(fMaxAlpha > alpha && fMinAlpha < alpha && qT < (sqrt(fQt*fQt*(1-(((alpha-(fMaxAlpha/2.))*(alpha-(fMaxAlpha/2.)))/(((fMaxAlpha-fMinAlpha)/2.)*((fMaxAlpha-fMinAlpha)/2.)))))))
          return kTRUE;
        else
          return kFALSE;
      }  
    }    
  }
  else{
    return kTRUE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
// Openingangle Cut on sigma daughters
Bool_t AliCaloSigmaCuts::SigmaDaughtersOpeningangleCut(Double_t openingangle){
  if(fMinOpeningAngle){
    if(fMaxOpeningAngle){
      if(fMaxOpeningAngle > openingangle && fMinOpeningAngle < openingangle)
        return kTRUE;
      else
        return kFALSE;
    }    
  }
  else{
    return kTRUE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
// Openingangle Cut on sigma daughters
Bool_t AliCaloSigmaCuts::TrackIsSelected(AliAODTrack* track, AliPIDResponse* fPIDResponse){
   
    if(!(track->TestFilterBit(fFilterBit))) return kFALSE;
    if(!(track->GetTPCNcls())) return kFALSE;
    if(!(track->GetTPCchi2())) return kFALSE;
    if((track->GetTPCNcls()) < fNClusterTPC || (track->GetTPCchi2perCluster()) > fChi2TPC) return kFALSE;

    if(!(track->GetITSNcls())) return kFALSE;
    if(!(track->GetITSchi2())) return kFALSE;
    if((track->GetITSNcls()) < fNClusterITS || (track->GetITSchi2()) > fChi2ITS) return kFALSE;
    
    if(!(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) return kFALSE;
    if(!(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton))) return kFALSE;
    if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) > fNSigmaTPC || (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton))) > fNSigmaTOF ) return kFALSE;

    Float_t dcaXY = 0.0, dcaZ = 0.0;
    track->GetImpactParameters(dcaXY,dcaZ);

    if(dcaXY < fDCAXY || dcaZ < fDCAZ) return kFALSE;

    return kTRUE;
  
}
