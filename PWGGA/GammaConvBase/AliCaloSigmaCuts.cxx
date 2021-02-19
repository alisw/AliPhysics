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
  "PIDVariationCut", //1
  "FilterBitCut", //2
  "NClusterTPCCut", //3
  "Chi2TPCCut", //4
  "NClusterITSCut", //5
  "Chi2ITSCut", //6
  "InvertChi2ITSCut", //7
  "DCAXYCut", //8
  "DCAZCut", //9
  "LowPNSigmaTPCCut", //10
  "NSigmaTPCCut", //11
  "NSigmaTOFCut", //12
  "PionMassLowerCut", //13
  "PionMassUpperCut", //14
  "PodolanskiCut", //15
  "OpeningAngleCut", //16
  "BackgroundEstimation", //17
};

//________________________________________________________________________
AliCaloSigmaCuts::AliCaloSigmaCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fCutString(NULL),
  fCutStringRead(""),
  fAmenterosCut(NULL),
  fPIDVariation(0),
  fFilterBit(0),
  fNClusterTPC(0),
  fChi2TPC(0),
  fNClusterITS(0),
  fChi2ITS(0),
  fInvertChi2ITS(0),
  fDCAXY(0),
  fDCAZ(0),
  fLowPNSigmaTPC(0),
  fNSigmaTPC(0),
  fNSigmaTOF(0),
  fMaxPionMass(0),
  fMinPionMass(0),
  fMaxAlpha(0),
  fMinAlpha(0),
  fQt(0),
  fMaxOpeningAngle(0),
  fMinOpeningAngle(0),
  fBackgroundestimation(0),
  fHistDEDx(0),
  fHistTOFBeta(0),
  fHistTPCSignal(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}
//________________________________________________________________________
AliCaloSigmaCuts::AliCaloSigmaCuts(const AliCaloSigmaCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fCutString(NULL),
  fCutStringRead(""),
  fAmenterosCut(NULL),
  fPIDVariation(ref. fPIDVariation),
  fFilterBit(ref. fFilterBit),
  fNClusterTPC(ref. fNClusterTPC),
  fChi2TPC(ref. fChi2TPC),
  fNClusterITS(ref. fNClusterITS),
  fChi2ITS(ref. fChi2ITS),
  fInvertChi2ITS(ref. fInvertChi2ITS),
  fDCAXY(ref. fDCAXY),
  fDCAZ(ref. fDCAZ),
  fLowPNSigmaTPC(ref. fLowPNSigmaTPC),
  fNSigmaTPC(ref. fNSigmaTPC),
  fNSigmaTOF(ref. fNSigmaTOF),
  fMaxPionMass(ref.fMaxPionMass),
  fMinPionMass(ref.fMinPionMass),
  fMaxAlpha(ref.fMaxAlpha),
  fMinAlpha(ref.fMinAlpha),
  fQt(ref.fQt),
  fMaxOpeningAngle(ref.fMaxOpeningAngle),
  fMinOpeningAngle(ref.fMinOpeningAngle),
  fBackgroundestimation(ref.fBackgroundestimation),
  fHistDEDx(ref.fHistDEDx),
  fHistTOFBeta(ref.fHistTOFBeta),
  fHistTPCSignal(ref.fHistTPCSignal)
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
void AliCaloSigmaCuts::InitCutHistograms(TString name){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);

  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }

  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("SigmaPlusCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  fHistDEDx = new TH2F("fHistDEDx", "fHistDEDx;#it{p};d#it{E}/d#it{x}", 200,0.,10.,200,1.,201.);
  fHistograms->Add(fHistDEDx);
  fHistTOFBeta = new TH2F("fHistTOFBeta", "fHistTOFBeta;#it{p};#beta", 200,0.,10.,130,0.1,1.3);
  fHistograms->Add(fHistTOFBeta);
  fHistTPCSignal = new TH2F("fHistTPCSignal", "fHistTPCSignal;#it{p};#sigma_{TPC}", 200, 0., 10., 60, -3., 3.);
  fHistograms->Add(fHistTPCSignal);
  TH1::AddDirectory(kTRUE);
  return;
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
  case kPIDVariation:
    if( SetPIDVariationCut(value)) {
      fCuts[kPIDVariation] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;
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
  case kInvertChi2ITS:
    if( SetInvertChi2ITSCut(value)) {
      fCuts[kInvertChi2ITS] = value;
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
  case kLowPNSigmaTPC:
    if( SetLowPNSigmaTPCCut(value)) {
      fCuts[kLowPNSigmaTPC] = value;
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
Bool_t AliCaloSigmaCuts::SetPIDVariationCut(Int_t PIDVariationCut){
  // Set Cut
  switch(PIDVariationCut){
  case 0:
    fPIDVariation = 0;  //TPC Only. 
    break;
  case 1:
    fPIDVariation = 1;  // TPC+TOF
    break;  
  case 2:
    fPIDVariation = 2;   //TPC with TOF veto
    break;   
  default:
    cout<<"Warning: PIDVariationCut not defined"<<PIDVariationCut<<endl;
    return kFALSE;
  }
  return kTRUE;
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
  case 2:
    fNClusterTPC = 50;
    break;
  case 3:
    fNClusterTPC = 65;
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
   case 2:
    fChi2TPC = 3.;
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
  case 2:
    fChi2ITS = 1.5;
    break;    
  default:
    cout<<"Warning: Chi2ITSCut not defined"<<Chi2ITSCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetInvertChi2ITSCut(Int_t InvertChi2ITSCut){
  // Set Cut
  switch(InvertChi2ITSCut){
  case 0:
    fInvertChi2ITS = 0;
    break;
  case 1:
    fInvertChi2ITS = 1;
    break;  
  default:
    cout<<"Warning: InvertChi2ITSCut not defined"<<InvertChi2ITSCut<<endl;
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
  case 2:
    fDCAXY = 0.02;
    break;  
  case 3:
    fDCAXY = 0.01;
    break;  
  case 4:
    fDCAXY = 0.015;
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
  case 2:
    fDCAZ = 0.02;
    break;
  case 3:
    fDCAZ = 0.01;
    break;  
  case 4:
    fDCAZ = 0.015;
    break;            
  default:
    cout<<"Warning: DCAZCut not defined"<<DCAZCut<<endl;
    return kFALSE;
  }
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliCaloSigmaCuts::SetLowPNSigmaTPCCut(Int_t LowPNSigmaTPCCut){
  // Set Cut
  switch(LowPNSigmaTPCCut){
  case 0:
    fLowPNSigmaTPC = 5.;
    break;
  case 1:
    fLowPNSigmaTPC = 2.;
    break;  
  default:
    cout<<"Warning: LowPNSigmaTPCCut not defined"<<LowPNSigmaTPCCut<<endl;
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
  case 2:
    fNSigmaTPC = 2.;
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
  case 2:
    fNSigmaTOF = 2.;
    break;    
  case 3:
    fNSigmaTOF = 5.;
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
  case 2:
    fMinPionMass = 0.125;
    break; 
  case 3:
    fMinPionMass = 0.11;
    break;
  case 4:
    fMinPionMass = 0.12;
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
  case 2:
    fMaxPionMass = 0.145;
    break;
  case 3:
    fMaxPionMass = 0.15;
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
  case 2:
    fMaxOpeningAngle = 0.5;
    fMinOpeningAngle = 0.05;
    break; 
  case 3:
    fMaxOpeningAngle = 0.3;
    fMinOpeningAngle = 0.05;
    break;  
  case 4:
    fMaxOpeningAngle = 0.25;
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
  if(!(track->GetTPCchi2perCluster())) return kFALSE;
  if((track->GetTPCNcls()) < fNClusterTPC || (track->GetTPCchi2perCluster()) > fChi2TPC) return kFALSE;

  if(!(track->GetITSNcls())) return kFALSE;
  if(!(track->GetITSchi2())) return kFALSE;
  if(fInvertChi2ITS == 0){
    if((track->GetITSNcls()) < fNClusterITS || (track->GetITSchi2()) > fChi2ITS) return kFALSE;
  }
  if(fInvertChi2ITS == 1){
    if((track->GetITSNcls()) < fNClusterITS || (track->GetITSchi2()) < fChi2ITS || (track->GetITSchi2()) > 10.) return kFALSE;
  }

  //dE/dx and TOF Beta Plot
  Double_t tpcSignal = track->GetTPCsignal();
  if(fHistDEDx) fHistDEDx->Fill(track->P(), tpcSignal);
  if(track->GetTOFsignal()){
    const float len = track->GetIntegratedLength();
    const float tim = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->GetTPCmomentum());
    const float beta = len / (tim * (2.99792457999999984e-02));
    if(fHistTOFBeta) fHistTOFBeta->Fill(track->P(), beta);
  }
  
  if(fHistTPCSignal) fHistTPCSignal->Fill(track->P(), fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
  
  if(!(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) return kFALSE;
  if(track->P() <= 0.1) return kFALSE;
  if(track->P() <= 0.8){
    if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) > fLowPNSigmaTPC ) return kFALSE;
  }
  if(track->P() > 0.8){
    if(fPIDVariation == 0){
      if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) > fNSigmaTPC) return kFALSE; 
    }
    if(fPIDVariation == 1){
      if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) > fNSigmaTPC || (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton))) > fNSigmaTOF ) return kFALSE; 
    }
    if(fPIDVariation == 2){
      if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) > fNSigmaTPC) return kFALSE; 
      if((TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton))) <= fNSigmaTPC){
        if((fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTOF,track)) == AliPIDResponse::kDetPidOk){
          if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) > fNSigmaTOF) return kFALSE; 
        }
      }
    }
  }
  return kTRUE;  
}
//_____________________________________________________________________________
// Openingangle Cut on sigma daughters
Bool_t AliCaloSigmaCuts::TrackIsSelectedByDCACut(AliAODTrack* track){
   
  Float_t dcaXY = 0.0, dcaZ = 0.0;
  track->GetImpactParameters(dcaXY,dcaZ);
  if(dcaXY < fDCAXY || dcaZ < fDCAZ) return kFALSE;

  return kTRUE;
  
}
