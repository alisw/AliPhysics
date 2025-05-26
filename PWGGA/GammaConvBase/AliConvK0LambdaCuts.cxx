/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.   *
 *                                                                          *
 * Authors: Friederike Bock                                                 *
 * Version 1.0                                                              *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ***************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling photon selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConvK0LambdaCuts.h"

#include "AliGAKFVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliDataFile.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliTRDTriggerAnalysis.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class iostream;
using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliConvK0LambdaCuts)
/// \endcond

const char* AliConvK0LambdaCuts::fgkCutNames[AliConvK0LambdaCuts::kNCuts] = {
  "V0FinderType",           // 0
  "EtaCut",                 // 1
  "MinRCut",                // 2
  "EtaForPhiCut",           // 3
  "MinPhiCut",              // 4
  "MaxPhiCut",              // 5
  "SinglePtCut",            // 6
  "ClsTPCCut",              // 7
  "ededxSigmaCut",          // 8
  "pidedxSigmaCut",         // 9
  "piMomdedxSigmaCut",      // 10
  "piMaxMomdedxSigmaCut",   // 11
  "LowPRejectionSigmaCut",  // 12
  "TOFelectronPID",         // 13
  "ITSelectronPID",         // 14 -- new ITS PID
  "TRDelectronPID",         // 15 -- new TRD PID
  "QtMaxCut",               // 16
  "Chi2GammaCut",           // 17
  "PsiPair",                // 18
  "DoPhotonAsymmetryCut",   // 19
  "CosinePointingAngle",    // 20
  "SharedElectronCuts",     // 21
  "RejectToCloseV0s",       // 22
  "DcaRPrimVtx",            // 23
  "DcaZPrimVtx",            // 24
  "EventPlane"              // 25
};

//________________________________________________________________________
AliConvK0LambdaCuts::AliConvK0LambdaCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(0),
  fMaxR(0.),
  fMinR(0.),
  fEtaCut(0.),
  fEtaCutMin(0.),
  fSinglePtCutPos(0.),
  fSinglePtCutNeg(0.),
  fMinClsTPC(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0.),
  fLineCutZRSlopeMin(0.),
  fLineCutZValueMin(0.),
  fChi2CutV0(0.),
  fPIDnSigmaAbove(0.),
  fPIDnSigmaBelow(0.),
  fTofPIDnSigmaAbove(0.),
  fTofPIDnSigmaBelow(0.),
  fDoQtGammaSelection(0),
  fDo2DQt(false),
  fQtMax(0.),
  fQtPtMax(0.),
  fUseOnFlyV0Finder(false),
  fCosPAngleCut(0.),
  fDCAZPrimVtxCut(0.),
  fDCAPrimVtxCut(0.),
  fCutString(nullptr),
  fCutStringRead(""),
  fIsHeavyIon(0),
  fExcludeMinR(0.),
  fExcludeMaxR(0.),
  fHistoCutIndex(nullptr),
  fHistoArmenterosbefore(nullptr),
  fHistoArmenterosafter(nullptr)
{
  InitPIDResponse();
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());
}

//________________________________________________________________________
AliConvK0LambdaCuts::AliConvK0LambdaCuts(const AliConvK0LambdaCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fMaxR(0.),
  fMinR(0.),
  fEtaCut(0.),
  fEtaCutMin(0.),
  fSinglePtCutPos(0.),
  fSinglePtCutNeg(0.),
  fMinClsTPC(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0.),
  fLineCutZRSlopeMin(0.),
  fLineCutZValueMin(0.),
  fChi2CutV0(0.),
  fPIDnSigmaAbove(0.),
  fPIDnSigmaBelow(0.),
  fTofPIDnSigmaAbove(0.),
  fTofPIDnSigmaBelow(0.),
  fDoQtGammaSelection(0),
  fDo2DQt(false),
  fQtMax(0.),
  fQtPtMax(0.),
  fUseOnFlyV0Finder(false),
  fCosPAngleCut(0.),
  fDCAZPrimVtxCut(0.),
  fDCAPrimVtxCut(0.),
  fCutString(nullptr),
  fCutStringRead(""),
  fIsHeavyIon(0),
  fExcludeMinR(0.),
  fExcludeMaxR(0.),
  fHistoCutIndex(nullptr),
  fHistoArmenterosbefore(nullptr),
  fHistoArmenterosafter(nullptr)
{
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());
}


//________________________________________________________________________
AliConvK0LambdaCuts::~AliConvK0LambdaCuts() {
  // Destructor
  if(fHistograms != NULL){
    delete fHistograms;
  }
}

//________________________________________________________________________
void AliConvK0LambdaCuts::InitCutHistograms(TString name, bool preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);


  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }

  if(fDoLightOutput==2) {
      AliInfo("Minimal output chosen");
      return;
  }

  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("K0Lambda_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  // // IsPhotonSelected
  fHistoCutIndex=new TH2F(Form("K0LambdaCuts %s",GetCutNumber().Data()),"IsPhotonSelected",11,-0.5,10.5, 300, 0, 30);
  fHistoCutIndex->GetXaxis()->SetBinLabel(1,"in");
  fHistoCutIndex->GetXaxis()->SetBinLabel(2,"onfly status");
  fHistoCutIndex->GetXaxis()->SetBinLabel(3,"eta");
  fHistoCutIndex->GetXaxis()->SetBinLabel(4,"radius");
  fHistoCutIndex->GetXaxis()->SetBinLabel(5,"min pt");
  fHistoCutIndex->GetXaxis()->SetBinLabel(6,"TPC Clus");
  fHistoCutIndex->GetXaxis()->SetBinLabel(7,"dEdX");
  fHistoCutIndex->GetXaxis()->SetBinLabel(8,"chi2");
  fHistoCutIndex->GetXaxis()->SetBinLabel(9,"Cos(PA)");
  fHistoCutIndex->GetXaxis()->SetBinLabel(10,"DCA");
  fHistoCutIndex->GetXaxis()->SetBinLabel(11,"out");
  fHistograms->Add(fHistoCutIndex);

  if(!fDoLightOutput){
    fHistoArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,100,0,1.);
    fHistograms->Add(fHistoArmenterosbefore);

    fHistoArmenterosafter=new TH2F(Form("Armenteros_after %s",GetCutNumber().Data()),"Armenteros_after",200,-1,1,100,0,1.);
    fHistograms->Add(fHistoArmenterosafter);
  }

  TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
bool AliConvK0LambdaCuts::InitPIDResponse(){
  // Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;

  }

  return kFALSE;
}

bool AliConvK0LambdaCuts::DoV0ReaderTypeCut(bool status) const{
  if(status == fUseOnFlyV0Finder){
    return true;
  }
  return false;
}

bool AliConvK0LambdaCuts::DoEtaCut(double etaV0) const{
  if(std::abs(etaV0) > fEtaCut) return false;
  return true;
}

bool AliConvK0LambdaCuts::DoRCut(AliAODv0 *v0) const{
  double tDecayVertexV0[3]; 
  v0->GetXYZ(tDecayVertexV0); 
  double RV0 = sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

  if(RV0 < fMinR) return false;
  else if(RV0 >fMaxR) return false;
  return true;
}

bool AliConvK0LambdaCuts::DoSinglePtCut(const AliAODTrack *trNeg, const AliAODTrack *trPos) const{
  if(trNeg->Pt() < fSinglePtCutNeg){
    return false;
  } else if (trPos->Pt() < fSinglePtCutPos){
    return false;
  }
  return true;
}

bool AliConvK0LambdaCuts::DoTPCClusCut(const AliAODTrack *trNeg, const AliAODTrack *trPos) const{
  if(trNeg->GetTPCClusterInfo(2,1) < fMinClsTPC) {
    return false;
  } else if (trNeg->GetTPCClusterInfo(2,1) < fMinClsTPC) {
    return false;
  }
  return true;
}

bool AliConvK0LambdaCuts::DodEdXCut(int mode, const AliAODTrack *trNeg, const AliAODTrack *trPos) const{
  if(mode == 2){ // K0s, both are pions
    if(fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kPion) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kPion) > fPIDnSigmaAbove){
      return false;
    } else if(fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kPion) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kPion) > fPIDnSigmaAbove){
      return false;
    }
    return true;
  } else if (mode == 3){ // Lambda: Pos track is proton, negative is pion
    if(fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kPion) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kPion) > fPIDnSigmaAbove){
      return false;
    } else if(fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kProton) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kProton) > fPIDnSigmaAbove){
      return false;
    }
    return true;
  } else if (mode == 4){ // Anti-Lambda: Neg track is proton, positive is pion
    if(fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kProton) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trNeg,AliPID::kProton) > fPIDnSigmaAbove){
      return false;
    } else if(fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kPion) < fPIDnSigmaBelow || fPIDResponse->NumberOfSigmasTPC(trPos,AliPID::kPion) > fPIDnSigmaAbove){
      return false;
    }
    return true;
  }
  return true;
}


bool AliConvK0LambdaCuts::DoChi2Cut(double chi2) const{
  if(chi2 > fChi2CutV0){
    return false;
  }
  return true;
}

bool AliConvK0LambdaCuts::DoCosPACut(double cosPA) const{
  if(cosPA < fCosPAngleCut){
    return false;
  }
  return true;
}

bool AliConvK0LambdaCuts::DoDCAToPrimVtxCut(double dca) const{
  if(dca > fDCAPrimVtxCut){
    return false;
  }
  return true;
}


bool AliConvK0LambdaCuts::IsK0sLambdaAccepted(AliAODv0 *v0, int fDoProcessK0){
  const AliAODTrack *ntrack=(AliAODTrack*) v0->GetDaughter(1);
  const AliAODTrack *ptrack=(AliAODTrack*) v0->GetDaughter(0);

  if(!ntrack || !ptrack){
    return false;
  }

  // Fill some histograms before
  if(fHistoArmenterosbefore)fHistoArmenterosbefore->Fill(v0->AlphaV0(),v0->PtArmV0());

  double cutindex = 0.;
  fHistoCutIndex->Fill(cutindex, v0->Pt());

  cutindex++;
  // Now do the cuts
  if(!DoV0ReaderTypeCut(v0->GetOnFlyStatus())){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  if(!DoEtaCut(v0->Eta())){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;
  
  if(!DoRCut(v0)){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  if(!DoSinglePtCut(ntrack, ptrack)){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  if(!DoTPCClusCut(ntrack, ptrack)){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  if(!DodEdXCut(fDoProcessK0, ntrack, ptrack)){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  if(!DoChi2Cut(v0->Chi2V0())){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  // if(!DoCosPACut(v0->GetV0CosineOfPointingAngle())){ // not a member of 
  //   fHistoCutIndex->Fill(cutindex, v0->Pt());
  //   return false;
  // }
  cutindex++;

  if(!DoDCAToPrimVtxCut(v0->DcaV0ToPrimVertex())){
    fHistoCutIndex->Fill(cutindex, v0->Pt());
    return false;
  }
  cutindex++;

  fHistoCutIndex->Fill(cutindex, v0->Pt());
  if(fHistoArmenterosafter)fHistoArmenterosafter->Fill(v0->AlphaV0(),v0->PtArmV0());
  return true;
}


///________________________________________________________________________
bool AliConvK0LambdaCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)
  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
bool AliConvK0LambdaCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  AliInfo(Form("Set Photoncut Number: %s",analysisCutSelection.Data()));
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
///________________________________________________________________________
bool AliConvK0LambdaCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID
  switch (cutID) {

    case kv0FinderType:
      if( SetV0Finder(value)) {
        fCuts[kv0FinderType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ketaCut:
      if( SetEtaCut(value)) {
        fCuts[ketaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kRCut:
      if( SetRCut(value)) {
        fCuts[kRCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kEtaForPhiSector:
      cout << "nothing to be done for now" << endl;
      return kTRUE;
    case kMinPhiSector:
      if( SetMinPhiSectorCut(value)) {
        fCuts[kMinPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kMaxPhiSector:
      if( SetMaxPhiSectorCut(value)) {
        fCuts[kMaxPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ksinglePtCut:
      if( SetSinglePtCut(value)) {
        fCuts[ksinglePtCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kclsTPCCut:
      if( SetTPCClusterCut(value)) {
        fCuts[kclsTPCCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kededxSigmaCut:
      if( SetTPCdEdxCut(value)) {
        fCuts[kededxSigmaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kpidedxSigmaCut:
      cout << "nothing to be done for now" << endl;
      return kTRUE;

    case kpiMomdedxSigmaCut:
      if( SetMinMomPiondEdxCut(value)) {
        fCuts[kpiMomdedxSigmaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kpiMaxMomdedxSigmaCut:
      cout << "nothing to be done for now" << endl;
      return kTRUE;

    case kLowPRejectionSigmaCut:
      cout << "nothing to be done for now" << endl;
      return kTRUE;

    case kTOFelectronPID:
      if( SetTOFElectronPIDCut(value)) {
        fCuts[kTOFelectronPID] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kQtMaxCut:
      if( SetQtMaxCut(value)) {
        fCuts[kQtMaxCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kchi2GammaCut:
      if( SetChi2V0Cut(value)) {
        fCuts[kchi2GammaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPsiPair:
      if( SetPsiPairCut(value)) {
        fCuts[kPsiPair] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kdoPhotonAsymmetryCut:
      if( SetAsymmetryCut(value)) {
        fCuts[kdoPhotonAsymmetryCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kCosPAngle:
      if( SetCosPAngleCut(value)) {
        fCuts[kCosPAngle] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kElecShare:
      cout << "nothing to be done for now" << endl;
      return kTRUE;
    case kToCloseV0s:
      cout << "nothing to be done for now" << endl;
      return kTRUE;

    case kDcaRPrimVtx:
      if( SetDCAPrimVtxCut(value)) {
        fCuts[kDcaRPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDcaZPrimVtx:
      if( SetDCAZPrimVtxCut(value)) {
        fCuts[kDcaZPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kInPlaneOutOfPlane:
      cout << "nothing to be done for now" << endl;
      return kTRUE;

    case kITSelectronPID:
      if( SetITSElectronPIDCut(value)) {
        fCuts[kITSelectronPID] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kTRDelectronPID:
      if( SetTRDElectronPIDCut(value)) {
        fCuts[kTRDelectronPID] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }
  AliError("Cut id %d not recognized");
  return kFALSE;
}
///________________________________________________________________________
void AliConvK0LambdaCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetV0Finder(Int_t v0FinderType){   // Set Cut
  switch (v0FinderType){
  case 0:  // on fly V0 finder
    cout << "have chosen onfly V0" << endl;
    fUseOnFlyV0Finder=kTRUE;
    break;
  case 1:  // offline V0 finder
    cout << "have chosen offline V0" << endl;
    fUseOnFlyV0Finder=kFALSE;
    break;
  default:
    AliError(Form(" v0FinderType not defined %d",v0FinderType));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetEtaCut(Int_t etaCut){   // Set Cut

  //Set Standard LineCutZValues
  fLineCutZValueMin = -2;
  fLineCutZValue = 7.;

  switch(etaCut){
  case 0: // 0.9
    fEtaCut     = 0.9;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 1:  // 0.6
    fEtaCut     = 0.6;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 2: // 0.5
    fEtaCut     = 0.5;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 3: // d - 0.8
    fEtaCut     = 0.8;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  default:
    AliError(Form(" EtaCut not defined %d",etaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetRCut(Int_t RCut){
  // Set Cut
  switch(RCut){
  case 0:
    fMinR=0;
    fMaxR = 180.;
    break;
  case 1:
    fMinR=2.8;
    fMaxR = 180.;
    break;
  case 2:
    fMinR=5.;
    fMaxR = 180.;
    break;
  case 3:
    fMaxR = 5.;
    fMinR = 100.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 4:
    fMaxR = 5.;
    fMinR = 70.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 5:
    fMaxR = 5.;
    fMinR = 30.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 6:
    fMaxR = 0.;
    fMinR = 30.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  default:
    AliError("RCut not defined");
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetEtaForPhiCut(Int_t etaPhiCut) {

  cout << "Nothing to be selected";

  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetMinPhiSectorCut(Int_t minPhiCut) {

  cout << "SetMinPhiSectorCut nothing to be selected" << endl;
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetMaxPhiSectorCut(Int_t maxPhiCut) {

  cout << "SetMaxPhiSectorCut nothing to be selected" << endl;

  return kTRUE;
}


///________________________________________________________________________
bool AliConvK0LambdaCuts::SetSinglePtCut(Int_t singlePtCut){   // Set Cut
  switch(singlePtCut){
  case 0: // 0.050 GeV
    fSinglePtCutPos = 0.050;
    fSinglePtCutNeg = 0.050;
    break;
  case 1:  // 0.100 GeV
    fSinglePtCutPos = 0.1;
    fSinglePtCutNeg = 0.1;
    break;
  case 2:  // 0.150 GeV
    fSinglePtCutPos = 0.15;
    fSinglePtCutNeg = 0.15;
    break;
  case 3:  // 0.2 GeV
    fSinglePtCutPos = 0.2;
    fSinglePtCutNeg = 0.2;
    break;
  case 4:  // 0.3 GeV
    fSinglePtCutPos = 0.3;
    fSinglePtCutNeg = 0.3;
    break;
  
  default:
    AliError(Form("singlePtCut not defined %d",singlePtCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetTPCClusterCut(Int_t clsTPCCut){   // Set Cut
  switch(clsTPCCut){
  case 0: // 0
    fMinClsTPC= 0.;
    break;
  case 1:  // 60
    fMinClsTPC= 60.;
    break;
  case 2:  // 80
    fMinClsTPC= 80.;
    break;
  case 3:  // 100
    fMinClsTPC= 100.;
    break;
  default:
    AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetTPCdEdxCut(Int_t ededxSigmaCut){   // Set Cut
  switch(ededxSigmaCut){
  case 0:
    fPIDnSigmaBelow=-10;
    fPIDnSigmaAbove=10;
    break;
  case 1: //
    fPIDnSigmaBelow=-1;
    fPIDnSigmaAbove=1;
    break;
  case 2:
    fPIDnSigmaBelow=-2;
    fPIDnSigmaAbove=2;
    break;
  case 3:
    fPIDnSigmaBelow=-3;
    fPIDnSigmaAbove=3;
    break;
  case 4:
    fPIDnSigmaBelow=-4;
    fPIDnSigmaAbove=4;
    break;
  case 5:
    fPIDnSigmaBelow=-5;
    fPIDnSigmaAbove=5;
    break;
  default:
    AliError("TPCdEdxCutElectronLine not defined");
    return kFALSE;

  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetMinMomPiondEdxCut(Int_t piMomdedxSigmaCut){   // Set Cut
  cout << "nothing to be done here SetMinMomPiondEdxCut" << endl;
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
  // Set Cut
  cout << "nothing to be done in SetTOFElectronPIDCut" << endl;
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetITSElectronPIDCut(Int_t ITSelectronPID){
  // Set Cut
  cout << "nothing to be done in SetITSElectronPIDCut" << endl;
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetTRDElectronPIDCut(Int_t TRDelectronPID){
  // Set Cut
  cout << "Not cut to be selected in SetTRDElectronPIDCut" << endl;
  return kTRUE;
}


///________________________________________________________________________
bool AliConvK0LambdaCuts::SetQtMaxCut(Int_t QtMaxCut){   // Set Cut
  switch(QtMaxCut){
  case 0: //
    fQtMax=1.;
    fDoQtGammaSelection=0;
    fDo2DQt=kFALSE;
    break;
  case 1:
    fQtPtMax=0.11;
    fQtMax=0.040;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  default:
    AliError(Form("Warning: QtMaxCut not defined %d",QtMaxCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetChi2V0Cut(Int_t chi2Cut){   // Set Cut

  switch(chi2Cut){
  case 0: // 100
    fChi2CutV0 = 100.;
    break;
  case 1:  // 50
    fChi2CutV0 = 50.;
    break;
  case 2:  // 30
    fChi2CutV0 = 30.;
    break;
  case 3:
    fChi2CutV0 = 200.;
    break;
  case 4:
    fChi2CutV0 = 100000.;
    break;
  default:
    AliError(Form("Warning: Chi2Cut not defined %d",chi2Cut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetPsiPairCut(Int_t psiCut) {
  cout << "SetPsiPairCut nothing to be selected" << endl;

  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetAsymmetryCut(Int_t doPhotonAsymmetryCut){
  cout << "SetAsymmetryCut nothing to be selected" << endl;
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetCosPAngleCut(Int_t cosCut) {

  switch(cosCut){
  case 0:
    fCosPAngleCut = -1;
    break;
  case 1:
    fCosPAngleCut = 0;
    break;
  case 2:
    fCosPAngleCut = 0.5;
    break;
  case 3:
    fCosPAngleCut = 0.75;
    break;
  case 4:
    fCosPAngleCut = 0.85;
    break;
  case 5:
    fCosPAngleCut = 0.88;
    break;
  case 6:
    fCosPAngleCut = 0.9;
    break;
  case 7:
    fCosPAngleCut = 0.95;
    break;
  case 8:
    fCosPAngleCut = 0.98;
    break;
  case 9:
    fCosPAngleCut = 0.99;
    break;
  case 10://a
    fCosPAngleCut = 0.995;
    break;
  case 11://b
    fCosPAngleCut = 0.985;
    break;
  case 12://c
    fCosPAngleCut = 0.996;
    break;
  case 13://d
    fCosPAngleCut = 0.997;
    break;
  case 14://e
    fCosPAngleCut = 0.998;
    break;
  case 15://f
    fCosPAngleCut = 0.999;
    break;
  default:
    AliError(Form("Cosine Pointing Angle cut not defined %d",cosCut));
    return kFALSE;
  }

  return kTRUE;
}


///________________________________________________________________________
bool AliConvK0LambdaCuts::SetDCAZPrimVtxCut(Int_t DCAZPrimVtx){
  // Set Cut
  switch(DCAZPrimVtx){
  case 0:  //
    fDCAZPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCAZPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCAZPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCAZPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCAZPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCAZPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCAZPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCAZPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCAZPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCAZPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCAZPrimVtx not defined "<<DCAZPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
bool AliConvK0LambdaCuts::SetDCAPrimVtxCut(Int_t DCARPrimVtx){
  // Set Cut
  switch(DCARPrimVtx){
  case 0:  //
    fDCAPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCAPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCAPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCAPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCAPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCAPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCAPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCAPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCAPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCAPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCARPrimVtx not defined "<<DCARPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
TString AliConvK0LambdaCuts::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}
