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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Check basic detector results at ESD level                            //
//   - Geometrical efficiency                                           //
//   - Tracking efficiency                                              //
//   - PID efficiency                                                   //
//   - Refit efficiency                                                 //
//                                                                      //
// Author                                                               //
//   Alex Bercuci <A.Bercuci@gsi.de>                                    //
//   Ionut Arsene <i.c.arsene@gsi.de>                                   //
//                                                                      //
//     The analysis task fills AliCFContainer objects defined in the    //
//   configuration macro using the AddCFContainer() method.             //
//   The CF containers can be filled with any of the variables defined  //
//   in ETrdCfVariables and at any of the steps defined in ETrdCfSteps. //
//     To define a new variable one needs to:                           //
//   1. Add an entry in the ETrdCfVariables enumeration                 //
//   2. Add the corresponding variable name (and in the correct order)  //
//      in fgkVarNames                                                  //
//   3. Define how the variable is filled in one of the Fill functions: //
//      FillEventInfo(), FillTrackInfo(), FillTrackletInfo(),           //
//      FillTrackletSliceInfo().                                        //
//     To define a new step one needs to:                               //
//   1. Add an entry in the ETrdCfSteps                                 //
//   2. Add the corresponding name in fgkStepNames                      //
//   3. Define the track level cuts for this step in IsTrackSelected()  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TH3S.h>
#include <TH3F.h>
#include <TProfile2D.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TParticle.h>
#include <TTimeStamp.h>
#include <TRandom.h>
#include <TString.h>

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDEvent.h"
#include "AliESDkink.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliESDpid.h"
#include "AliExternalTrackParam.h"

#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliTrackReference.h"
#include "AliMultiplicity.h"
#include "AliCFContainer.h"

#include "AliTRDcheckESD.h"
#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliTRDcheckESD)

const Float_t AliTRDcheckESD::fgkxTPC = 290.;
const Float_t AliTRDcheckESD::fgkxTOF = 365.;
const Char_t* AliTRDcheckESD::fgkVarNames[AliTRDcheckESD::kNTrdCfVariables] = {
    "vtxZ", "multiplicity", "trigger", "BC", "TOFBC", "DCAxy", "DCAz", "charge", "OuterParam rad.", "phiVtx", "phi", 
    "etaVtx", "eta", "pt", "ptTRD", "P", "PTRD", "TRDchi2", "tracklets", "clusters", "TrdQuality", 
    "TrdBudget", "TOFchi2", "Qtot0", "ClustersPerRows", "Clusters/tracklet", "TrdP", "TrdPloss", "layer", "slice", "PH0"
}; 
const Char_t* AliTRDcheckESD::fgkStepNames[AliTRDcheckESD::kNSteps] = {"TPC", "TRD", "TOF", "TOFin", "TOFout"};  

FILE* AliTRDcheckESD::fgFile = NULL;

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD():
  AliAnalysisTaskSE()
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fESDpid(new AliESDpid)
  ,fHistos(NULL)
  ,fReferenceTrackFilter(NULL)
  ,fPhysSelTriggersEnabled(kFALSE)
  ,fUserEnabledTriggers("")
  ,fNAssignedTriggers(0)
{
  //
  // Default constructor
  //
  SetNameTitle("TRDcheckESD", "Check TRD @ ESD level");
  SetMC(kTRUE);
}

//____________________________________________________________________
AliTRDcheckESD::AliTRDcheckESD(char* name):
  AliAnalysisTaskSE(name)
  ,fStatus(0)
  ,fNRefFigures(0)
  ,fESD(NULL)
  ,fMC(NULL)
  ,fESDpid(new AliESDpid)
  ,fHistos(NULL)
  ,fReferenceTrackFilter(NULL)
  ,fPhysSelTriggersEnabled(kFALSE)
  ,fUserEnabledTriggers("")
  ,fNAssignedTriggers(0)
{
  //
  // Default constructor
  //
  SetMC(kTRUE);
  SetTitle("Check TRD @ ESD level");
  DefineOutput(1, TObjArray::Class());
}

//____________________________________________________________________
AliTRDcheckESD::~AliTRDcheckESD()
{
  // Destructor
  if(fHistos && !(AliAnalysisManager::GetAnalysisManager() && AliAnalysisManager::GetAnalysisManager()->IsProofMode())){
    if(fHistos->IsOwner()) fHistos->Delete();
    delete fHistos;
    fHistos = NULL;
  }
}

//____________________________________________________________________
void AliTRDcheckESD::FillEventInfo(Double_t* values) {
  //
  // Fill event information
  //
  values[kEventVtxZ] = fESD->GetPrimaryVertex()->GetZ();
  values[kEventBC] = fESD->GetBunchCrossNumber();
  
  const AliMultiplicity* mult=fESD->GetMultiplicity();
  Double_t itsNTracklets = mult->GetNumberOfTracklets();
  values[kEventMult] = itsNTracklets;
}

//____________________________________________________________________
void AliTRDcheckESD::FillTrackInfo(Double_t* values, AliESDtrack* esdTrack) {
  //
  // Fill track information
  //
  Float_t dcaxy,dcaz;
  esdTrack->GetImpactParameters(dcaxy,dcaz);
  values[kTrackDCAxy]  = dcaxy;
  values[kTrackDCAz]   = dcaz;  
  values[kTrackCharge] = esdTrack->Charge();
  values[kTrackPt]     = esdTrack->Pt();
  values[kTrackPhi]    = esdTrack->Phi();
  values[kTrackEta]    = esdTrack->Eta();
  values[kTrackP]      = esdTrack->P();
  values[kTrackTrdTracklets] = esdTrack->GetTRDntracklets();
  values[kTrackTrdClusters]  = esdTrack->GetTRDncls();
  values[kTrackTrdChi2]      = esdTrack->GetTRDchi2()/(esdTrack->GetTRDntracklets()>0 ? esdTrack->GetTRDntracklets() : 1.0);
  values[kTrackTrdQuality]   = esdTrack->GetTRDQuality();
  values[kTrackTRDBudget]    = -1.0*esdTrack->GetTRDBudget();
  values[kTrackTOFBC]        = esdTrack->GetTOFBunchCrossing(fESD->GetMagneticField());
  values[kTrackTOFchi2]      = esdTrack->GetTOFchi2();
  const AliExternalTrackParam *out=esdTrack->GetOuterParam();
  Double_t p[3];
  if(out->GetXYZ(p)) 
    values[kTrackOuterParamRadius] = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  else 
    values[kTrackOuterParamRadius] = 0.0;
}

//____________________________________________________________________
void AliTRDcheckESD::FillTrackletInfo(Double_t* values, AliESDtrack* esdTrack, Int_t iPlane,
                                      Double_t* localSagitaPhi, Double_t localMom[][3], Bool_t* localMomGood) {
  //
  // Fill TRD tracklet info
  //
  values[kTrackletClustersVsRows] = esdTrack->GetTRDtrkltClCross(iPlane);
  values[kTrackletClusters]       = esdTrack->GetTRDtrkltOccupancy(iPlane);
  values[kTrackletQtot]           = esdTrack->GetTRDslice(iPlane, 0);
  values[kTrackletP]              = esdTrack->GetTRDmomentum(iPlane);
  values[kTrackPlossTRDlayer]     = 1000.0*(esdTrack->P() - values[kTrackletP]);    // p loss in MeV    
  values[kTrackletLayer]          = iPlane;
  values[kTrackPhiTRD]            = localSagitaPhi[iPlane];
  values[kTrackPtTRD]         = (localMomGood[iPlane] ? TMath::Sqrt(localMom[iPlane][0]*localMom[iPlane][0]+
	                                                            localMom[iPlane][1]*localMom[iPlane][1]) : values[kTrackPt]);
  values[kTrackPTRD]          = (localMomGood[iPlane] ? TMath::Sqrt(values[kTrackPtTRD]*values[kTrackPtTRD]+
                                                                    localMom[iPlane][2]*localMom[iPlane][2]) : values[kTrackP]);	
  values[kTrackEtaTRD] = values[kTrackPTRD]-(localMomGood[iPlane] ? localMom[iPlane][2] : esdTrack->Pz());
  values[kTrackEtaTRD] = (TMath::Abs(values[kTrackPTRD])>1.0e-8 ? (values[kTrackPTRD]+(localMomGood[iPlane] ? localMom[iPlane][2] : esdTrack->Pz()))/values[kTrackEtaTRD] : 0.0);
  values[kTrackEtaTRD] = (values[kTrackEtaTRD]>1.0e-8 ? 0.5*TMath::Log(values[kTrackEtaTRD]) : -999.);
}

//____________________________________________________________________
void AliTRDcheckESD::FillTrackletSliceInfo(Double_t* values, AliESDtrack* esdTrack, Int_t iSlice) {
  //
  // Fill TRD tracklet info
  //
  values[kTrackletPHslice] = esdTrack->GetTRDslice(Int_t(values[kTrackletLayer]), iSlice);
  values[kTrackletSlice] = iSlice;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::IsTrackSelected(AliESDtrack* track, Double_t* /*values*/, Int_t step) {
  //
  // Select tracks at each step
  //  
  Bool_t referenceFilter = fReferenceTrackFilter->IsSelected(track);
  if(step==kTPCreference) {    // reference track filter
    return referenceFilter;
  }
  if(step==kTRD) {    // TRD reference track filter
    return (referenceFilter && Int_t(track->GetTRDntracklets()>0));
  }
  if(step==kTOF) {    // TRD+(TOFout || TOFpid) request
    return (referenceFilter && Int_t(track->GetTRDntracklets())>0 && 
           ((track->GetStatus() & AliESDtrack::kTOFout) || (track->GetStatus() & AliESDtrack::kTOFpid))); 
  }
  if(step==kTOFin) {    // TOFin request
    return (referenceFilter && (track->GetStatus() & AliESDtrack::kTOFin)); 
  }
  if(step==kTOFout) {    // TOFout request
    return (referenceFilter && (track->GetStatus() & AliESDtrack::kTOFout)); 
  }
  return kFALSE;
}

//____________________________________________________________________
void AliTRDcheckESD::UserCreateOutputObjects()
{	
  //
  // Create Output Containers (TObjectArray containing 1D histograms)
  //
  Histos();
  PostData(1, fHistos);
}

//____________________________________________________________________
void AliTRDcheckESD::MakeSummaryFromCF(Double_t* trendValues, const Char_t* /*triggerName*/, Bool_t /*useIsolatedBC*/, Bool_t /*cutTOFbc*/){
  //
  // Draw summary plots for the ESDcheck task using the CF container
  //

  cout << "Make summary from CF" << endl;
  TCanvas *cOut=0x0;
  if(gROOT->FindObject("trackingSummary")) delete gROOT->FindObject("trackingSummary");
  cOut = new TCanvas("trackingSummary", "Tracking summary for the ESD task", 1600, 1200);
  cOut->cd();
  //PlotTrackingSummaryFromCF(trendValues, triggerName, useIsolatedBC, cutTOFbc);
  PlotTrackingSummaryFromCF(trendValues);
  cOut->SaveAs("trackingSummary.gif");
  
  if(gROOT->FindObject("pidSummary")) delete gROOT->FindObject("pidSummary");
  cOut = new TCanvas("pidSummary", "PID summary for the ESD task", 1600, 1200);
  cOut->cd();
  //PlotPidSummaryFromCF(trendValues, triggerName, useIsolatedBC, cutTOFbc);
  PlotPidSummaryFromCF(trendValues);
  cOut->SaveAs("pidSummary.gif");

  if(gROOT->FindObject("centSummary")) delete gROOT->FindObject("centSummary");
  cOut = new TCanvas("centSummary", "Centrality summary for the ESD task", 1600, 1200);
  cOut->cd();
  //PlotCentSummaryFromCF(trendValues, triggerName, useIsolatedBC, cutTOFbc);
  PlotCentSummaryFromCF(trendValues);
  cOut->SaveAs("centSummary.gif");
  
  PlotOtherSummaryFromCF(trendValues);
  
  if(trendValues)
    for(Int_t i=0;i<50;++i) cout << "trend #" << i << " :: " << trendValues[i] << endl;
}


//____________________________________________________________________
void AliTRDcheckESD::UserExec(Option_t *){
  //
  // Run the Analysis
  //
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fMC = MCEvent();

  if(!fESD){
    AliError("ESD event missing.");
    return;
  }
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) return;
  
  if(!fPhysSelTriggersEnabled) {
    InitializeCFContainers();
    fPhysSelTriggersEnabled = kTRUE;
  }
    
  UInt_t isSelected = AliVEvent::kAny;
  if(inputHandler){
    if(inputHandler->GetEventSelection()) {
      isSelected = inputHandler->IsEventSelected();
    }
  }
  if(!isSelected) return;

  TString triggerClasses = fESD->GetFiredTriggerClasses();
  //cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++triggers fired:  " << triggerClasses.Data() << endl;
  TObjArray* triggers = triggerClasses.Tokenize(" ");
  TObjArray* userTriggers = fUserEnabledTriggers.Tokenize(";");
  if(triggers->GetEntries()<1) {delete triggers; delete userTriggers; return;}
  Bool_t hasGoodTriggers = kFALSE;
  Int_t triggerIndices[kNMaxAssignedTriggers] = {0};
  Int_t nTrigFired=0;
  Bool_t trigAlreadyChecked = kFALSE;
  Bool_t trigSelected = kFALSE;
  Int_t trigIdx = 0;
  for(Int_t i=0; i<triggers->GetEntries(); ++i) {
    TString trigStr=triggers->At(i)->GetName();
    if(!trigStr.Contains("NOTRD") && !trigStr.Contains("MUON")) hasGoodTriggers = kTRUE;    // check wheter TRD was read out in this event
    if(trigStr.Contains("NOTRD")) continue;
    if(trigStr.Contains("MUON")) continue;
    if(i>=kNMaxAssignedTriggers) continue;
    
    hasGoodTriggers = kTRUE;
    // enable the "All triggers" bit
    trigIdx = 1;
    trigAlreadyChecked = kFALSE;
    for(Int_t k=0;k<nTrigFired;++k) 
      if(triggerIndices[k]==trigIdx) {
	trigAlreadyChecked = kTRUE;
	break;
    }
    if(!trigAlreadyChecked) triggerIndices[nTrigFired++] = trigIdx;  
    
    trigSelected = kFALSE;
    // check whether this trigger matches any of the user defined trigger families
    for(Int_t j=0;j<userTriggers->GetEntries();++j) {
      TString userTrigStr=userTriggers->At(j)->GetName();
      if(trigStr.Contains(userTrigStr.Data())) {
	trigSelected = kTRUE;
	trigIdx = GetTriggerIndex(userTrigStr.Data(), kFALSE);
	trigAlreadyChecked = kFALSE;
	for(Int_t k=0;k<nTrigFired;++k) 
	  if(triggerIndices[k]==trigIdx) {
	    trigAlreadyChecked = kTRUE;
	    break;
	  }
	if(!trigAlreadyChecked) {    // add trigger to the list of enabled triggers only if it was not added already
	  triggerIndices[nTrigFired++] = trigIdx;  
	}
      }
    }
    
    trigIdx = GetTriggerIndex(trigStr.Data(), kFALSE);
    if(trigIdx>0) trigSelected = kTRUE;
    if(trigIdx==-1) trigIdx=1;  
    trigAlreadyChecked = kFALSE;
    for(Int_t k=0;k<nTrigFired;++k) 
      if(triggerIndices[k]==trigIdx) {
        trigAlreadyChecked = kTRUE;
	break;
      }
    if(!trigAlreadyChecked) {
      triggerIndices[nTrigFired++]=1;  // 0-assigned to all other triggers
    }
  }  // end loop over triggers  
  
  if(!trigSelected && hasGoodTriggers) {
    triggerIndices[nTrigFired++]=2;
  }
  
  TH1F* hTrig = (TH1F*)fHistos->FindObject("hTriggerDefs");
  for(Int_t i=0; i<nTrigFired; ++i)
    hTrig->Fill(triggerIndices[i]);

  if(!hasGoodTriggers) {
    PostData(1, fHistos);
    delete triggers;
    delete userTriggers;
    return;
  }
  
  Int_t* trigFiredIdx=new Int_t[nTrigFired];
  for(Int_t i=0;i<nTrigFired;++i) trigFiredIdx[i] = triggerIndices[i];
  
  // Get MC information if available
  //AliStack * fStack = NULL;
  if(HasMC()){
    if(!fMC){ 
      AliWarning("MC event missing");
      SetMC(kFALSE);
    } else {
      if(!fMC->Stack()){
        AliWarning("MC stack missing");
        SetMC(kFALSE);
      }
    }
  }
  
  Double_t values[kNTrdCfVariables];      // array where the CF container variables are stored
  for(Int_t i=0;i<kNTrdCfVariables; ++i) values[i] = -999.;
  FillEventInfo(values);  
  
  Int_t multLimits[6] = {0, 700, 1400, 2100, 2800, 3500};
  Int_t centralityClass = 0;
  for(Int_t iCent=0; iCent<5; ++iCent) {
    if(values[kEventMult]>=multLimits[iCent] && values[kEventMult]<multLimits[iCent+1])
      centralityClass=iCent+1;
  }
  if(centralityClass == 0) return;
  
  // radius of TRD entrance plane in each layer
  Double_t rTRD[6] = {298.0, 311.0, 324.0, 337.0, 350.0, 363.0};
  
  AliESDtrack *esdTrack(NULL);
  for(Int_t itrk = 0; itrk < fESD->GetNumberOfTracks(); itrk++){
    esdTrack = fESD->GetTrack(itrk);
    
    Bool_t stepSelections[kNSteps]; 
    for(Int_t is=0;is<kNSteps;++is) {
      stepSelections[is] = IsTrackSelected(esdTrack, values, is);
    }
    if(!stepSelections[0]) continue;
    
    FillTrackInfo(values, esdTrack);
    
    // find position and momentum of the track at entrance in TRD
    const AliExternalTrackParam *outerParam = esdTrack->GetOuterParam();
    Double_t localCoord[6][3] = {{0.0}};
    Bool_t localCoordGood[6];
    for(Int_t il=0;il<6;++il) 
      localCoordGood[il] = (outerParam ? outerParam : esdTrack)->GetXYZAt(rTRD[il], fESD->GetMagneticField(), localCoord[il]);  
    Double_t localMom[6][3] = {{0.0}};
    Bool_t localMomGood[6];
    for(Int_t il=0; il<6; ++il) 
      localMomGood[il] = (outerParam ? outerParam : esdTrack)->GetPxPyPzAt(rTRD[il], fESD->GetMagneticField(), localMom[il]);
    Double_t localSagitaPhi[6] = {-999.};
    for(Int_t il=0; il<6; ++il) 
      localSagitaPhi[il] = (localCoordGood[il] ? TMath::ATan2(localCoord[il][1], localCoord[il][0]) : -999.);
    if(!localMomGood[0]) continue;
    
    // fill tracklet values such that the TRD local coordinates are filled
    FillTrackletInfo(values, esdTrack, 0, localSagitaPhi, localMom, localMomGood);
        
    for(Int_t itrig=0; itrig<nTrigFired; ++itrig) {
      values[kEventTrigger] = Double_t(trigFiredIdx[itrig]);
      
      // check if cf needs tracklet or slice info
      FillGlobalTrackContainers(values, stepSelections, itrig);
      
      for(Int_t iPlane=0; iPlane<6; iPlane++) {
        FillTrackletInfo(values, esdTrack, iPlane, localSagitaPhi, localMom, localMomGood);
	if(values[kTrackletQtot]>20.0) FillTrdTrackletContainers(values, stepSelections, itrig);
	
	for(Int_t iSlice=0; iSlice<8; iSlice++) {
	  FillTrackletSliceInfo(values, esdTrack, iSlice);
	  if(values[kTrackletPHslice]>20.0) FillTrdSliceContainers(values, stepSelections, itrig);
	}  // end loop over slices
      }  // end loop over TRD layers
    }  // end loop over triggers
  }  // end loop over tracks
  
  delete triggers;
  delete userTriggers;
  delete [] trigFiredIdx;  
  
  PostData(1, fHistos);
}

//____________________________________________________________________
TObjArray* AliTRDcheckESD::Histos()
{
  // Retrieve histograms array if already build or build it
  if(!fHistos) {
    fHistos = new TObjArray();
    fHistos->SetOwner(kTRUE);
  }
  
  TH1* h = 0x0;
  // Trigger definitions
  if(!(h=(TH1F*)gROOT->FindObject("hTriggerDefs"))) {
    h = new TH1F("hTriggerDefs", "Trigger definitions", kNMaxAssignedTriggers, 0.5, 0.5+Float_t(kNMaxAssignedTriggers));
  }
  else h->Reset();
  fHistos->Add(h);
    
  return fHistos;
}


//__________________________________________________________________________________________________________
void AliTRDcheckESD::InitializeCFContainers() {
  //
  //  Initialize the CF container
  //
  AliAnalysisManager* man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) return;

  GetTriggerIndex("All triggers", kTRUE);
  GetTriggerIndex("Not specified triggers", kTRUE);

  AliPhysicsSelection* physSel = (AliPhysicsSelection*)inputHandler->GetEventSelection();
  const TList* trigList = (physSel ? physSel->GetCollisionTriggerClasses() : 0x0);
  const TList* bgTrigList = (physSel ? physSel->GetBGTriggerClasses() : 0x0);
  
  // Add collision triggers from PhysicsSelection
  if(trigList) {
    for(Int_t it=0; it<trigList->GetEntries(); ++it) {
      TString trigName = trigList->At(it)->GetName();
      TObjArray* arr = trigName.Tokenize(" ");
      trigName = arr->At(0)->GetName();
      trigName.Remove(0,1);
      TObjArray* arr2 = trigName.Tokenize(",");
      for(Int_t jt=0; jt<arr2->GetEntries(); ++jt) {
        // Assign an index into the trigger histogram and the CF container for this trigger
        GetTriggerIndex(arr2->At(jt)->GetName(), kTRUE);
      }
      delete arr;
    }
  }
  // Add background triggers from PhysicsSelection
  if(bgTrigList) {
    for(Int_t it=0; it<bgTrigList->GetEntries(); ++it) {
      TString trigName = bgTrigList->At(it)->GetName();
      TObjArray* arr = trigName.Tokenize(" ");
      trigName = arr->At(0)->GetName();
      trigName.Remove(0,1);
      TObjArray* arr2 = trigName.Tokenize(",");
      for(Int_t jt=0; jt<arr2->GetEntries(); ++jt) {
        // Assign an index into the trigger histogram and the CF container for this trigger
        GetTriggerIndex(arr2->At(jt)->GetName(), kTRUE);
      }
      delete arr;
    }
  }
    
  // Add user enabled triggers
  TObjArray* arr = fUserEnabledTriggers.Tokenize(";");
  for(Int_t it=0; it<arr->GetEntries(); ++it) {
    GetTriggerIndex(arr->At(it)->GetName(), kTRUE);
  }
  delete arr;
}


//__________________________________________________________________________________________________________
void AliTRDcheckESD::AddCFContainer(const Char_t* name, const Char_t* title,
				    Int_t nSteps, Int_t* steps, 
				    Int_t nVars, UInt_t* vars, TArrayD* binLimits) {
  //
  // Add a CF container
  //
  if(!fHistos) {
    fHistos = new TObjArray();
    fHistos->SetOwner(kTRUE);
  }
  // get number of bins for each variable
  Int_t* nBins = new Int_t[nVars];
  for(Int_t iv=0;iv<nVars;++iv)
    nBins[iv] = binLimits[iv].GetSize()-1;  
  // create the CF container
  AliCFContainer* cf = new AliCFContainer(name, title, nSteps, nVars, nBins);
  // set CF container variable binning and name
  for(Int_t iv=0;iv<nVars;++iv) {
    cf->SetBinLimits(iv, binLimits[iv].GetArray());
    cf->SetVarTitle(iv, fgkVarNames[vars[iv]]);
  }
  // set the step names 
  for(Int_t is=0; is<nSteps; ++is) {
    cf->SetStepTitle(is, fgkStepNames[steps[is]]);
    for(Int_t iv=0;iv<nVars;++iv) cf->GetAxis(iv, is)->SetUniqueID(vars[iv]);
  }
  fHistos->Add(cf);
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::FillTrdSliceContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig) {
  //
  // fill TRD slice info
  //
  if(!fHistos) return;
  for(Int_t i=0;i<fHistos->GetEntries();++i) {
    TString objType = fHistos->At(i)->IsA()->GetName();
    if(!objType.Contains("AliCFContainer")) continue;
    AliCFContainer* cf = (AliCFContainer*)fHistos->At(i);
    TString varNames="";
    for(Int_t ivar=0;ivar<cf->GetNVar();++ivar) {
      varNames += cf->GetVarTitle(ivar); varNames += ";";
    }
    //if(cf->GetVar(fgkVarNames[kTrackletSlice])<0 && cf->GetVar(fgkVarNames[kTrackletPHslice])<0) continue;
    if(!varNames.Contains(fgkVarNames[kTrackletSlice]) && !varNames.Contains(fgkVarNames[kTrackletPHslice])) continue;
    //if((cf->GetVar(fgkVarNames[kEventTrigger])<0 && itrig==0) || (cf->GetVar(fgkVarNames[kEventTrigger])>=0))
    if((!varNames.Contains(fgkVarNames[kEventTrigger]) && itrig==0) || varNames.Contains(fgkVarNames[kEventTrigger]))
      FillCFContainer(cf, values, stepSelections);
  }
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::FillTrdTrackletContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig) {
  //
  // fill global track info
  //
  if(!fHistos) return;
  for(Int_t i=0;i<fHistos->GetEntries();++i) {
    TString objType = fHistos->At(i)->IsA()->GetName();
    if(!objType.Contains("AliCFContainer")) continue;
    AliCFContainer* cf = (AliCFContainer*)fHistos->At(i);
    TString varNames="";
    for(Int_t ivar=0;ivar<cf->GetNVar();++ivar) {
      varNames += cf->GetVarTitle(ivar); varNames += ";";
    }
    //if(cf->GetVar(fgkVarNames[kTrackletSlice])>=0 || cf->GetVar(fgkVarNames[kTrackletPHslice])>=0) continue;
    if(varNames.Contains(fgkVarNames[kTrackletSlice]) || varNames.Contains(fgkVarNames[kTrackletPHslice])) continue;
    //if(cf->GetVar(fgkVarNames[kTrackletLayer])<0 && cf->GetVar(fgkVarNames[kTrackletQtot])<0) continue;
    if(!varNames.Contains(fgkVarNames[kTrackletLayer]) && !varNames.Contains(fgkVarNames[kTrackletQtot])) continue;
    //if((cf->GetVar(fgkVarNames[kEventTrigger])<0 && itrig==0) || (cf->GetVar(fgkVarNames[kEventTrigger])>=0))
    if((!varNames.Contains(fgkVarNames[kEventTrigger]) && itrig==0) || varNames.Contains(fgkVarNames[kEventTrigger]))
      FillCFContainer(cf, values, stepSelections);
  }
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::FillGlobalTrackContainers(Double_t* values, Bool_t* stepSelections, Int_t itrig) {
  //
  // fill global track info
  //
  if(!fHistos) return;
  for(Int_t i=0;i<fHistos->GetEntries();++i) {
    TString objType = fHistos->At(i)->IsA()->GetName();
    if(!objType.Contains("AliCFContainer")) continue;
    AliCFContainer* cf = (AliCFContainer*)fHistos->At(i);
    TString varNames="";
    for(Int_t ivar=0;ivar<cf->GetNVar();++ivar) {
      varNames += cf->GetVarTitle(ivar); varNames += ";";
    }
    /*if(cf->GetVar(fgkVarNames[kTrackletLayer])>=0 || 
       cf->GetVar(fgkVarNames[kTrackletSlice])>=0 || 
       cf->GetVar(fgkVarNames[kTrackletQtot])>=0 || 
       cf->GetVar(fgkVarNames[kTrackletPHslice])>=0) continue;*/
    if(varNames.Contains(fgkVarNames[kTrackletLayer]) || 
       varNames.Contains(fgkVarNames[kTrackletSlice]) || 
       varNames.Contains(fgkVarNames[kTrackletQtot]) || 
       varNames.Contains(fgkVarNames[kTrackletPHslice])) continue;
    /*if((cf->GetVar(fgkVarNames[kEventTrigger])<0 && itrig==0) || 
       (cf->GetVar(fgkVarNames[kEventTrigger])>=0))*/
    if((!varNames.Contains(fgkVarNames[kEventTrigger]) && itrig==0) || 
       varNames.Contains(fgkVarNames[kEventTrigger]))
      FillCFContainer(cf, values, stepSelections);
  }
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::FillCFContainer(AliCFContainer* cf, Double_t* values, Bool_t* stepSelections) {
  //
  // Fill CF container
  //
  Double_t* cfValues=new Double_t[cf->GetNVar()];
  for(Int_t iv=0;iv<cf->GetNVar();++iv)
    cfValues[iv] = values[cf->GetAxis(iv,0)->GetUniqueID()];
    
  for(Int_t istep=0; istep<cf->GetNStep(); ++istep) {
    TString stepTitle = cf->GetStepTitle(istep);
    Int_t stepNo = -1;
    for(Int_t is=0;is<kNSteps;++is) 
      if(!stepTitle.CompareTo(fgkStepNames[is])) {
	stepNo = is;
	break;
      }
    if(stepSelections[stepNo]) cf->Fill(cfValues, istep);
  }  // end loop over steps

  delete [] cfValues;
}

//____________________________________________________________________
Bool_t AliTRDcheckESD::Load(const Char_t *file, const Char_t *dir, const Char_t *name)
{
  //
  // Load data from performance file
  //
  if(!TFile::Open(file)){
    AliWarning(Form("Couldn't open file %s.", file));
    return kFALSE;
  }
  if(dir){
    if(!gFile->cd(dir)){
      AliWarning(Form("Couldn't cd to %s in %s.", dir, file));
      return kFALSE;
    }
  }
  TObjArray *o(NULL);
  const Char_t *tn=(name ? name : GetName());
  if(!(o = (TObjArray*)gDirectory->Get(tn))){
    AliWarning(Form("Missing histogram container %s.", tn));
    return kFALSE;
  }
  fHistos = (TObjArray*)o->Clone(GetName());
    
  TH1F* trigHist = (TH1F*)fHistos->FindObject("hTriggerDefs");
  for(Int_t i=1;i<=trigHist->GetXaxis()->GetNbins();++i) {
    if(trigHist->GetXaxis()->GetBinLabel(i)[0]!='\0') ++fNAssignedTriggers;
  }
  gFile->Close();
  return kTRUE;
}

//_______________________________________________________
Bool_t AliTRDcheckESD::PutTrendValue(const Char_t *name, Double_t val)
{
// Dump trending value to default file

  if(!fgFile){
    fgFile = fopen("TRD.Performance.txt", "at");
  }
  fprintf(fgFile, "%s_%s %f\n", GetName(), name, val);
  return kTRUE;
}

//____________________________________________________________________
void AliTRDcheckESD::Terminate(Option_t *)
{
  // Steer post-processing 
  if(!fHistos){
    fHistos = dynamic_cast<TObjArray *>(GetOutputData(1));
    if(!fHistos){
      AliError("Histogram container not found in output");
      return;
    }
  }
}

//____________________________________________________________________
Int_t AliTRDcheckESD::Pdg2Idx(Int_t pdg) const
{
  //
  // Helper function converting PDG code into AliPID index
  //
  switch(pdg){
  case kElectron: 
  case kPositron: return AliPID::kElectron;  
  case kMuonPlus:
  case kMuonMinus: return AliPID::kMuon;  
  case kPiPlus: 
  case kPiMinus: return AliPID::kPion;  
  case kKPlus: 
  case kKMinus: return AliPID::kKaon;
  case kProton: 
  case kProtonBar: return AliPID::kProton;
  } 
  return -1;
}

  
//________________________________________________________
void AliTRDcheckESD::Process2D(TH2 * const h2, TGraphErrors **g)
{
  //
  // Do the processing
  //

  Int_t n = 0;
  if((n=g[0]->GetN())) for(;n--;) g[0]->RemovePoint(n);
  if((n=g[1]->GetN())) for(;n--;) g[1]->RemovePoint(n);
  TF1 f("fg", "gaus", -3.,3.);
  for(Int_t ibin = 1; ibin <= h2->GetNbinsX(); ibin++){
    Double_t x = h2->GetXaxis()->GetBinCenter(ibin);
    TH1D *h = h2->ProjectionY("py", ibin, ibin);
    if(h->GetEntries()<100) continue;
    //AdjustF1(h, f);

    h->Fit(&f, "QN");
    Int_t ip = g[0]->GetN();
    g[0]->SetPoint(ip, x, f.GetParameter(1));
    g[0]->SetPointError(ip, 0., f.GetParError(1));
    g[1]->SetPoint(ip, x, f.GetParameter(2));
    g[1]->SetPointError(ip, 0., f.GetParError(2));
  }
  return;
}
//____________________________________________________________________
void AliTRDcheckESD::PrintStatus(ULong_t status)
{
// Dump track status to stdout

  printf("ITS[i(%d) o(%d) r(%d)] TPC[i(%d) o(%d) r(%d) p(%d)] TRD[i(%d) o(%d) r(%d) p(%d) s(%d)] HMPID[o(%d) p(%d)]\n"
    ,Bool_t(status & AliESDtrack::kITSin)
    ,Bool_t(status & AliESDtrack::kITSout)
    ,Bool_t(status & AliESDtrack::kITSrefit)
    ,Bool_t(status & AliESDtrack::kTPCin)
    ,Bool_t(status & AliESDtrack::kTPCout)
    ,Bool_t(status & AliESDtrack::kTPCrefit)
    ,Bool_t(status & AliESDtrack::kTPCpid)
    ,Bool_t(status & AliESDtrack::kTRDin)
    ,Bool_t(status & AliESDtrack::kTRDout)
    ,Bool_t(status & AliESDtrack::kTRDrefit)
    ,Bool_t(status & AliESDtrack::kTRDpid)
    ,Bool_t(status & AliESDtrack::kTRDStop)
    ,Bool_t(status & AliESDtrack::kHMPIDout)
    ,Bool_t(status & AliESDtrack::kHMPIDpid)
  );
}

//____________________________________________________________________
TH1D* AliTRDcheckESD::Proj2D(TH2* hist, TH1* mpvErr, TH1* widthErr, TH1* chi2) {
  //
  // project the PH vs Slice 2D-histo into a 1D histo with Landau MPV and widths
  //
  
  TH1D* hProjection = (TH1D*)hist->ProjectionX(Form("hProjection_%f", gRandom->Rndm()));
  hProjection->Reset();
  
  TF1* fitLandau = new TF1("landauFunc","landau",20.,3000.);
  TH1D *hD;
  for(Int_t iBin=1;iBin<=hist->GetXaxis()->GetNbins();iBin++) {
    if(gROOT->FindObject("projection"))
      delete gROOT->FindObject("projection");
    hD = (TH1D*)hist->ProjectionY("projection",iBin,iBin);
    //hD->Rebin(4);
    if(hD->Integral()>10) {
      fitLandau->SetParameter(1, hD->GetBinCenter(hD->GetMaximumBin()));
      //fitLandau->SetParLimits(1, 0.2*hD->GetBinCenter(hD->GetMaximumBin()), 3.0*hD->GetBinCenter(hD->GetMaximumBin()));
      fitLandau->SetParameter(0, 1000.);
      //fitLandau->SetParLimits(0, 1., 10000000.);
      fitLandau->SetParameter(2, 0.5*hD->GetBinCenter(hD->GetMaximumBin()));
      //fitLandau->SetParLimits(2, 0.01*hD->GetBinCenter(hD->GetMaximumBin()), 1.0*hD->GetRMS());
      hD->Fit(fitLandau, "Q0", "", hD->GetXaxis()->GetXmin(), hD->GetXaxis()->GetXmax());
      hD->Fit(fitLandau, "Q0", "", hD->GetXaxis()->GetXmin(), hD->GetXaxis()->GetXmax());
      hProjection->SetBinContent(iBin, fitLandau->GetParameter(1));
      hProjection->SetBinError(iBin, fitLandau->GetParameter(2));
      if(mpvErr) {
	mpvErr->SetBinContent(iBin, fitLandau->GetParameter(1));
	mpvErr->SetBinError(iBin, fitLandau->GetParError(1));
      }
      if(widthErr) {
	widthErr->SetBinContent(iBin, fitLandau->GetParameter(2));
	widthErr->SetBinError(iBin, fitLandau->GetParError(2));
      }
      if(chi2) {
	chi2->SetBinContent(iBin, (fitLandau->GetNDF()>0 ? fitLandau->GetChisquare()/Double_t(fitLandau->GetNDF()) : 0.0));
      }
    }
    else{
      hProjection->SetBinContent(iBin, 0);
      hProjection->SetBinError(iBin, 0);
    }
  }
  return hProjection;
}

//____________________________________________________________________
TH2F* AliTRDcheckESD::Proj3D(TH3* hist, TH2* accMap, Int_t zbinLow, Int_t zbinHigh, Float_t &entries) {
  //
  //  Project a 3D histogram to a 2D histogram in the Z axis interval [zbinLow,zbinHigh] 
  //  Return the 2D histogram and also the number of entries into this projection (entries)

  Int_t nBinsX = hist->GetXaxis()->GetNbins();   // X and Y axis bins are assumed to be all equal
  Float_t minX = hist->GetXaxis()->GetXmin();
  Float_t maxX = hist->GetXaxis()->GetXmax();
  Int_t nBinsY = hist->GetYaxis()->GetNbins();
  Float_t minY = hist->GetYaxis()->GetXmin();
  Float_t maxY = hist->GetYaxis()->GetXmax();
  Int_t nBinsZ = hist->GetZaxis()->GetNbins();  // Z axis bins (pt) might have different widths

  TH2F* projHisto = (TH2F*)gROOT->FindObject("projHisto");
  if(projHisto) 
    projHisto->Reset();
  else
    projHisto = new TH2F("projHisto", "projection", nBinsX, minX, maxX, nBinsY, minY, maxY);

  for(Int_t iZ=1; iZ<=nBinsZ; iZ++) {
    if(iZ<zbinLow) continue;
    if(iZ>zbinHigh) continue;
    for(Int_t iX=1; iX<=nBinsX; iX++) {
      for(Int_t iY=1; iY<=nBinsY; iY++) {
        if(accMap) {
          if(accMap->GetBinContent(iX,iY)>0.1)
            projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
        }
        else    // no acc. cut 
          projHisto->SetBinContent(iX, iY, projHisto->GetBinContent(iX, iY)+hist->GetBinContent(iX,iY,iZ));
        // count only the entries which are inside the acceptance map
        if(accMap) {
          if(accMap->GetBinContent(iX,iY)>0.1)
            entries+=hist->GetBinContent(iX,iY,iZ);
        }
        else    // no acc. cut
          entries+=hist->GetBinContent(iX,iY,iZ);
      }
    }
  }
  return projHisto;
}

//____________________________________________________________________
void AliTRDcheckESD::CheckActiveSM(TH1D* phiProj, Bool_t activeSM[18]) {
  //
  // Check the active super-modules
  //
  Double_t entries[18] = {0.0};
  Double_t smPhiLimits[19];
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  for(Int_t phiBin=1; phiBin<=phiProj->GetXaxis()->GetNbins(); ++phiBin) {
    Double_t phi = phiProj->GetBinCenter(phiBin);
    Int_t sm = -1;
    for(Int_t ism=0; ism<18; ++ism) 
      if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1]) sm = ism;
    if(sm==-1) continue;
    entries[sm] += phiProj->GetBinContent(phiBin);
  }
  Double_t avEntries = Double_t(phiProj->Integral())/18.0;
  for(Int_t ism=0; ism<18; ++ism) 
    if(entries[ism]>0.5*avEntries) activeSM[ism] = kTRUE;
}


//__________________________________________________________________________________________________
TH1F* AliTRDcheckESD::EfficiencyFromPhiPt(AliCFContainer* cf, Int_t minNtrkl, Int_t maxNtrkl, 
					  Int_t stepNom, Int_t stepDenom, Int_t var) {
  //
  // Use the CF container to extract the efficiency vs pt (other variables beside pt also posible)
  //
  Int_t varTrackPhi = cf->GetVar(fgkVarNames[kTrackPhiTRD]);
  Int_t otherVar = cf->GetVar(fgkVarNames[var]);  
  Int_t trdStepNumber = cf->GetStep(fgkStepNames[kTRD]);
  Int_t tpcStepNumber = cf->GetStep(fgkStepNames[kTPCreference]);
    
  TH1D* phiProj = (TH1D*)cf->Project(trdStepNumber, varTrackPhi);
  Bool_t activeSM[18] = {kFALSE};
  CheckActiveSM(phiProj, activeSM); delete phiProj;
  Double_t smPhiLimits[19];
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  
  TH2D* hNomPhiVar=0x0;
  TH2D* hDenomPhiVar=0x0;
  
  if((stepNom!=tpcStepNumber) &&
     (minNtrkl>-1 && minNtrkl<7 && maxNtrkl>-1 && maxNtrkl<7)) {
    cf->SetRangeUser(cf->GetVar(fgkVarNames[kTrackTrdTracklets]), Double_t(minNtrkl), Double_t(maxNtrkl));
    hNomPhiVar = (TH2D*)cf->Project(stepNom, otherVar, varTrackPhi);
    cf->SetRangeUser(cf->GetVar(fgkVarNames[kTrackTrdTracklets]), 0.0,6.0);
  }
  else
    hNomPhiVar = (TH2D*)cf->Project(stepNom, otherVar, varTrackPhi);
  if((stepDenom!=tpcStepNumber) &&
     (minNtrkl>-1 && minNtrkl<7 && maxNtrkl>-1 && maxNtrkl<7)) {
    cf->SetRangeUser(cf->GetVar(fgkVarNames[kTrackTrdTracklets]), Double_t(minNtrkl), Double_t(maxNtrkl));
    hDenomPhiVar = (TH2D*)cf->Project(stepDenom, otherVar, varTrackPhi);
    cf->SetRangeUser(cf->GetVar(fgkVarNames[kTrackTrdTracklets]), 0.0,6.0);
  } 
  else
    hDenomPhiVar = (TH2D*)cf->Project(stepDenom, otherVar, varTrackPhi);
    
  TH1F* hEff = new TH1F(Form("hEff%s_%d_%d_%f", fgkVarNames[var], stepNom, stepDenom, gRandom->Rndm()), "", 
			hNomPhiVar->GetXaxis()->GetNbins(), hNomPhiVar->GetXaxis()->GetXbins()->GetArray());
  for(Int_t ib=1;ib<=hNomPhiVar->GetXaxis()->GetNbins();++ib)
    hEff->GetXaxis()->SetBinLabel(ib, hNomPhiVar->GetXaxis()->GetBinLabel(ib));
  
  for(Int_t ivar=1; ivar<=hEff->GetXaxis()->GetNbins(); ++ivar) {
    Double_t nom = 0.0; Double_t denom = 0.0;
    Double_t eff = 0.0; Double_t err = 0.0;
    for(Int_t iphi=1; iphi<=hNomPhiVar->GetYaxis()->GetNbins(); ++iphi) {
      Double_t phi = hNomPhiVar->GetYaxis()->GetBinCenter(iphi);
      Bool_t isActive = kFALSE;
      for(Int_t ism=0; ism<18; ++ism) 
        if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1] && activeSM[ism]) 
	  isActive = kTRUE;
      if(!isActive) continue;
      nom += hNomPhiVar->GetBinContent(ivar, iphi);
      denom += hDenomPhiVar->GetBinContent(ivar, iphi);
    }
    eff = (denom>0.001 ? nom/denom : 0.0);
    err = (denom>0.001 && (denom-nom)>0.001 && nom>0.001 ? (TMath::Sqrt(nom*(denom-nom)/denom/denom/denom)) : 0.0);
    hEff->SetBinContent(ivar, eff);
    hEff->SetBinError(ivar, err);
  }   // end loop over pt bins
  delete hNomPhiVar; delete hDenomPhiVar;
  return hEff;
}


//____________________________________________________________________
TH1F* AliTRDcheckESD::EfficiencyTRD(TH3* tpc3D, TH3* trd3D, Bool_t useAcceptance) {
  //
  // Calculate the TRD-TPC matching efficiency as function of pt
  //
  
  if(!tpc3D || !trd3D) return NULL;
  Int_t nBinsZ = trd3D->GetZaxis()->GetNbins();
  // project everything on the eta-phi map to obtain an acceptance map
  Float_t nada = 0.;
  TH2F *trdAcc = (useAcceptance ? (TH2F*)Proj3D(trd3D, 0x0, 1, nBinsZ, nada)->Clone(Form("trdAcc%f", gRandom->Rndm())) : 0x0);
  TH1D *phiProj = (trdAcc ? trdAcc->ProjectionY(Form("phiProj%f", gRandom->Rndm())) : 0x0);
  
  // prepare the acceptance map
  Bool_t activeSM[18] = {kFALSE};
  Double_t smPhiLimits[19];
  for(Int_t ism=0; ism<=18; ++ism) smPhiLimits[ism] = -TMath::Pi() + (2.0*TMath::Pi()/18.0)*ism;
  if(phiProj) {
    CheckActiveSM(phiProj, activeSM);   // get the active SMs
    trdAcc->Reset();
    // Put 1 entry in every bin which belongs to an active SM
    for(Int_t iY=1; iY<=trdAcc->GetYaxis()->GetNbins(); ++iY) {
      Double_t phi = trdAcc->GetYaxis()->GetBinCenter(iY);
      Bool_t isActive = kFALSE;
      for(Int_t ism=0; ism<18; ++ism) {
        if(phi>=smPhiLimits[ism] && phi<smPhiLimits[ism+1] && activeSM[ism]) {
	  isActive = kTRUE;
        }
      }
      if(!isActive) continue;
      for(Int_t iX=1; iX<=trdAcc->GetXaxis()->GetNbins(); ++iX) 
        if(trdAcc->GetXaxis()->GetBinCenter(iX)>=-0.85 && trdAcc->GetXaxis()->GetBinCenter(iX)<=0.85) trdAcc->SetBinContent(iX, iY, 1.0);
    }  // end for over Y(phi) bins
  }  // end if phiProj
    
  // get the bin limits from the Z axis of 3D histos
  Float_t *ptBinLimits = new Float_t[nBinsZ+1];
  for(Int_t i=1; i<=nBinsZ; i++) {
    ptBinLimits[i-1] = trd3D->GetZaxis()->GetBinLowEdge(i);
  }
  ptBinLimits[nBinsZ] = trd3D->GetZaxis()->GetBinUpEdge(nBinsZ);
  
  TH1F *efficiency = new TH1F(Form("eff%d", Int_t(1000000.0*gRandom->Rndm())), "TRD-TPC matching efficiency", nBinsZ, ptBinLimits);
  
  // loop over Z bins
  Bool_t effGood = kFALSE;
  for(Int_t i=1; i<=nBinsZ; i++) {
    Float_t tpcEntries = 0.0; Float_t trdEntries = 0.0;
    Proj3D(tpc3D, trdAcc, i, i, tpcEntries);
    Proj3D(trd3D, trdAcc, i, i, trdEntries);
    Float_t ratio = 0;
    if(tpcEntries>0) ratio = trdEntries/tpcEntries;
    Float_t error = 0;
    if(tpcEntries>0 && trdEntries>0 && (tpcEntries-trdEntries)>=0.0) 
      error = TMath::Sqrt(trdEntries*(tpcEntries-trdEntries)/tpcEntries/tpcEntries/tpcEntries);
    if(ratio>0.001) {
      efficiency->SetBinContent(i,ratio);
      efficiency->SetBinError(i,error);
      effGood = kTRUE;
    }
  }     // end loop over Z bins
  if(!effGood) return 0x0;
  
  return efficiency;
}


//__________________________________________________________________________________________________
void AliTRDcheckESD::PlotCentSummaryFromCF(Double_t* /*trendValues*/, const Char_t* /*triggerName*/, Bool_t /*useIsolatedBC*/, Bool_t /*cutTOFbc*/) {
  //
  // Make the centrality summary figure from the CF container 
  // 
  if(!fHistos) return;
  AliCFContainer* matchPt=(AliCFContainer*)fHistos->FindObject("MatchingPt");
  if(!matchPt) return;
  AliCFContainer* centCF=(AliCFContainer*)fHistos->FindObject("CentralityCF");
  if(!centCF) return;
  AliCFContainer* clustersCF=(AliCFContainer*)fHistos->FindObject("ClustersCF");
  if(!clustersCF) return;
  
  TLatex* lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);

  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001); gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  TVirtualPad* pad=0x0;
  
  if(gROOT->FindObject("rangeEffPt")) delete gROOT->FindObject("rangeEffPt");
  TH2F* rangeEffPt=new TH2F("rangeEffPt", "",10,0.,10.,10,0.,1.3);
  rangeEffPt->SetStats(kFALSE);
  SetStyle(rangeEffPt->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  
  Int_t padsForEffs[5] = {0,3,6,1,4};
  for(Int_t iCent=1; iCent<6; ++iCent) {
    pad = ((TVirtualPad*)l->At(padsForEffs[iCent-1])); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02); pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    rangeEffPt->Draw();
    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.7, rangeEffPt->GetXaxis()->GetXmax(), 0.7);
    line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.9, rangeEffPt->GetXaxis()->GetXmax(), 0.9);
    
    matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kEventMult]), Double_t(iCent), Double_t(iCent), kTRUE);
       
    matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), +1.0, +1.0);  
    TH1F* hEffPosAll = EfficiencyFromPhiPt(matchPt,  0,  6, 1, 0);
    TH1F* hEffPosTrk4 = EfficiencyFromPhiPt(matchPt, 4,  4, 1, 0);  
    TH1F* hEffPosTrk5 = EfficiencyFromPhiPt(matchPt, 5,  5, 1, 0);
    TH1F* hEffPosTrk6 = EfficiencyFromPhiPt(matchPt, 6,  6, 1, 0);
     
    matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);  
    TH1F* hEffNegAll = EfficiencyFromPhiPt(matchPt, 0, 6, 1, 0);  
    TH1F* hEffNegTrk4 = EfficiencyFromPhiPt(matchPt, 4, 4, 1, 0);  
    TH1F* hEffNegTrk5 = EfficiencyFromPhiPt(matchPt, 5, 5, 1, 0);  
    TH1F* hEffNegTrk6 = EfficiencyFromPhiPt(matchPt, 6, 6, 1, 0);
    matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackTrdTracklets]), 0.0, 6.0);  
    matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), -1.0, +1.0);  
    
    SetStyle(hEffPosAll,  1, kRed, 1, 24, kRed, 1);
    SetStyle(hEffPosTrk4, 1, kRed, 1, 25, kRed, 1);
    SetStyle(hEffPosTrk5, 1, kRed, 1, 26, kRed, 1);
    SetStyle(hEffPosTrk6, 1, kRed, 1, 27, kRed, 1);
    SetStyle(hEffNegAll,  1, kBlue, 1, 24, kBlue, 1);
    SetStyle(hEffNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
    SetStyle(hEffNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
    SetStyle(hEffNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
    hEffPosAll->Draw("same");
    hEffNegAll->Draw("same");
    hEffPosTrk4->Draw("same");
    hEffNegTrk4->Draw("same");
    hEffPosTrk5->Draw("same");
    hEffNegTrk5->Draw("same");
    hEffPosTrk6->Draw("same");
    hEffNegTrk6->Draw("same");    
        
    TLegend* leg=new TLegend(0.18, 0.7, 0.77, 0.89);
    if(iCent==1) {
      leg->SetFillColor(0);
      leg->SetNColumns(2);
      leg->SetMargin(0.1);
      leg->SetBorderSize(0);
      leg->AddEntry(hEffPosAll,  "pos. (#geq 1 tracklet)", "p");
      leg->AddEntry(hEffNegAll,  "neg. (#geq 1 tracklet)", "p");
      leg->AddEntry(hEffPosTrk4, "pos. (4 tracklets)", "p");
      leg->AddEntry(hEffNegTrk4, "neg. (4 tracklets)", "p");
      leg->AddEntry(hEffPosTrk5, "pos. (5 tracklets)", "p");
      leg->AddEntry(hEffNegTrk5, "neg. (5 tracklets)", "p");
      leg->AddEntry(hEffPosTrk6, "pos. (6 tracklets)", "p");     
      leg->AddEntry(hEffNegTrk6, "neg. (6 tracklets)", "p");
      leg->Draw();
    }
    lat->DrawLatex(0.2, 1.32, Form("%.0f < SPD tracklets < %.0f", matchPt->GetAxis(matchPt->GetVar(fgkVarNames[kEventMult]),0)->GetBinLowEdge(iCent), matchPt->GetAxis(matchPt->GetVar(fgkVarNames[kEventMult]),0)->GetBinUpEdge(iCent)));
  }   // end for loop over multiplicity classes
  
  // Reset the modified user ranges of the CF container
  matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kEventMult]), 0., 3500.);
     
  // Cluster distributions in all multiplicity classes
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeNcls")) delete gROOT->FindObject("rangeNcls");
  TH2F* rangeNcls = new TH2F("rangeNcls", "", 10, 0.0, 199.9, 10, 0.0, 1.199);
  SetStyle(rangeNcls->GetXaxis(), "# TRD clusters", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNcls->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeNcls->SetStats(kFALSE);
  rangeNcls->Draw();
    
  TH1D* hNcls[6]={0x0};
  TLegend* legCls=new TLegend(0.7, 0.75, 0.97, 0.97);
  legCls->SetBorderSize(0);
  legCls->SetFillColor(0);
  legCls->SetMargin(0.15);
  
  for(Int_t iCent=0; iCent<6; ++iCent) {
    if(iCent>0)
      clustersCF->SetRangeUser(clustersCF->GetVar(fgkVarNames[kEventMult]), Double_t(iCent), Double_t(iCent), kTRUE);
    hNcls[iCent] = (TH1D*)clustersCF->Project(0, clustersCF->GetVar(fgkVarNames[kTrackTrdClusters]));
    if(!hNcls[iCent]) continue;
    
    hNcls[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
    Double_t maximum = hNcls[iCent]->GetMaximum();
    if(maximum>1.0)
      hNcls[iCent]->Scale(1.0/maximum);
    hNcls[iCent]->SetStats(kFALSE);
    hNcls[iCent]->SetTitle("");
    hNcls[iCent]->SetLineWidth(2);
    
    if(hNcls[iCent]->Integral()>0.01) {
      hNcls[iCent]->Draw("same");
      legCls->AddEntry(hNcls[iCent], (iCent==0 ? "all centralities" : Form("%.0f < SPD tracklets < %.0f", clustersCF->GetAxis(clustersCF->GetVar(fgkVarNames[kEventMult]),0)->GetBinLowEdge(iCent), 
									   clustersCF->GetAxis(clustersCF->GetVar(fgkVarNames[kEventMult]),0)->GetBinUpEdge(iCent))), "l");
    }
  }
  legCls->Draw();
  clustersCF->SetRangeUser(clustersCF->GetVar(fgkVarNames[kEventMult]), 0.0, 6.0, kTRUE);
  
  // Qtot distributions
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.02); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  if(gROOT->FindObject("rangeQtot")) delete gROOT->FindObject("rangeQtot");
  TH2F* rangeQtot = new TH2F("rangeQtot", "", 10, 0.0, 9.999, 10, 0.0, 1.199);
  SetStyle(rangeQtot->GetXaxis(), "Q_{tot} (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtot->GetYaxis(), "entries (a.u.)", 0.07, 0.8, kTRUE, 0.05);
  rangeQtot->SetStats(kFALSE);
  rangeQtot->Draw();
  
  TH1D* hQtot[6+1]={0x0};
  TLegend* leg2=new TLegend(0.6, 0.7, 0.9, 0.97);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  
  for(Int_t iCent=0; iCent<6; ++iCent) {
    if(iCent>0)
      centCF->SetRangeUser(centCF->GetVar(fgkVarNames[kEventMult]), Double_t(iCent), Double_t(iCent), kTRUE);
    
    hQtot[iCent] = (TH1D*)centCF->Project(0, centCF->GetVar(fgkVarNames[kTrackletQtot]));
    if(!hQtot[iCent]) continue;
    hQtot[iCent]->SetBinContent(1, 0);
    
    Double_t maximum = hQtot[iCent]->GetMaximum();
    if(maximum>1.0)
      hQtot[iCent]->Scale(1.0/maximum);
    hQtot[iCent]->SetLineColor(iCent<4 ? iCent+1 : iCent+2);
    hQtot[iCent]->SetStats(kFALSE);
    hQtot[iCent]->SetTitle("");
    hQtot[iCent]->SetLineWidth(2);
    if(hQtot[iCent]->Integral()>0.01) {
      hQtot[iCent]->Draw(iCent==0 ? "" : "same");
      leg2->AddEntry(hQtot[iCent], (iCent==0 ? "all centralities" : Form("%.0f < SPD tracklets < %.0f", centCF->GetAxis(centCF->GetVar(fgkVarNames[kEventMult]),0)->GetBinLowEdge(iCent), 
									   centCF->GetAxis(centCF->GetVar(fgkVarNames[kEventMult]),0)->GetBinUpEdge(iCent))), "l");
    }
  }
  leg2->Draw();
  centCF->SetRangeUser(centCF->GetVar(fgkVarNames[kEventMult]), 0.0, 5.0, kTRUE);
}



//_________________________________________________________________
void AliTRDcheckESD::PlotTrackingSummaryFromCF(Double_t* trendValues, const Char_t* /*triggerName*/, Bool_t /*useIsolatedBC*/, Bool_t /*cutTOFbc*/) {
  //
  //  Plot tracking summary
  //
  //  trendValues will be filled with trending variables
  //  trendValues[0] : TPC-TRD matching efficiency for positive tracks in the range 1.0<pt<3.0 GeV/c
  //  trendValues[1] : statistical error of trendValues[0]
  //  trendValues[2] : TPC-TRD matching efficiency for negative tracks in the range 1.0<pt<3.0 GeV/c
  //  trendValues[3] : statistical error of trendValues[2]
  //  trendValues[4] : TRD-TOF matching efficiency for positive tracks in the range 1.0<pt<3.0 GeV/c
  //  trendValues[5] : statistical error of trendValues[4]
  //  trendValues[6] : TRD-TOF matching efficiency for negative tracks in the range 1.0<pt<3.0 GeV/c
  //  trendValues[7] : statistical error of trendValues[6]
  //  trendValues[8] : Average number of TRD tracklets per track in the range 1.0<pt<3.0 GeV/c
  //  trendValues[9] : statistical error of trendValues[8]
  //  trendValues[10]: Average number of TRD clusters per track in the range 1.0<p<3.0 GeV/c
  //  trendValues[11]: statistical error of trendValues[10]
  // 
  if(!fHistos) return;
  AliCFContainer* matchPhiEta=(AliCFContainer*)fHistos->FindObject("MatchingPhiEta");
  if(!matchPhiEta) return;
  AliCFContainer* matchPt=(AliCFContainer*)fHistos->FindObject("MatchingPt");
  if(!matchPt) return;
  AliCFContainer* clustersCF=(AliCFContainer*)fHistos->FindObject("ClustersCF");
  if(!clustersCF) return;
  AliCFContainer* bcCF=(AliCFContainer*)fHistos->FindObject("BunchCrossingsCF");
  if(!bcCF) return;
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);
  lat->SetTextFont(42);
  
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  
  // eta-phi distr. for positive TPC tracks
  TVirtualPad* pad = ((TVirtualPad*)l->At(0)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  TH2D* hTPCrefPos = 0x0; TH2D* hTRDrefPos = 0x0; TH2D* hTOFrefPos = 0x0;
  TH2D* hTPCrefNeg = 0x0; TH2D* hTRDrefNeg = 0x0; TH2D* hTOFrefNeg = 0x0;
  matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackTrdTracklets]), 0.0, 6.0);
  matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackCharge]), +1.0, +1.0);      // positive charges
  hTPCrefPos = (TH2D*)matchPhiEta->Project(0, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  hTRDrefPos = (TH2D*)matchPhiEta->Project(1, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  hTOFrefPos = (TH2D*)matchPhiEta->Project(2, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);      // negative charges
  hTPCrefNeg = (TH2D*)matchPhiEta->Project(0, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  hTRDrefNeg = (TH2D*)matchPhiEta->Project(1, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  hTOFrefNeg = (TH2D*)matchPhiEta->Project(2, matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]));
  matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackCharge]), -1.0, +1.0);      // reset charge cut
    
  if(gROOT->FindObject("rangeEtaPhi")) delete gROOT->FindObject("rangeEtaPhi");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi", "", 10, -0.99, +0.99, 10, -3.15, 3.15);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);  
  
  //----------------------------------------------
  // eta-phi efficiency for positive TRD tracks
  pad = ((TVirtualPad*)l->At(0)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  TH2D* hTRDeffPos = (hTRDrefPos ? (TH2D*)hTRDrefPos->Clone("hTRDeffPos") : 0x0);
  if(hTRDeffPos) {
    hTRDeffPos->Reset();
    hTRDeffPos->SetStats(kFALSE);
    hTRDeffPos->Divide(hTRDrefPos, hTPCrefPos);
    hTRDeffPos->SetMaximum(1.0);
    hTRDeffPos->Draw("samecolz");
    lat->DrawLatex(-0.9, 3.3, "TPC-TRD matching for positive tracks");
    DrawTRDGrid();
  }
  
  //----------------------------------------------
  // eta-phi efficiency for negative TRD tracks
  pad = ((TVirtualPad*)l->At(3)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  TH2D* hTRDeffNeg = (hTRDrefNeg ? (TH2D*)hTRDrefNeg->Clone("hTRDeffNeg") : 0x0);
  if(hTRDeffNeg) {
    hTRDeffNeg->Reset();
    hTRDeffNeg->SetStats(kFALSE);
    hTRDeffNeg->Divide(hTRDrefNeg, hTPCrefNeg);
    hTRDeffNeg->SetMaximum(1.0);
    hTRDeffNeg->Draw("samecolz");
    lat->DrawLatex(-0.9, 3.3, "TPC-TRD matching for negative tracks");
    DrawTRDGrid();  
  }
  
  //----------------------------------------------
  // eta-phi TRD-TOF matching efficiency for positive tracks
  pad = ((TVirtualPad*)l->At(1)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  TH2D* hTOFeffPos = (hTOFrefPos ? (TH2D*)hTOFrefPos->Clone("hTOFeffPos") : 0x0);
  if(hTOFeffPos) {
    hTOFeffPos->Reset();
    hTOFeffPos->SetStats(kFALSE);
    hTOFeffPos->Divide(hTOFrefPos, hTRDrefPos);
    hTOFeffPos->SetMaximum(1.0);
    hTOFeffPos->Draw("samecolz");
    lat->DrawLatex(-0.9, 3.3, "TRD-TOF matching for positive tracks");
    DrawTRDGrid();
  }
  
  //----------------------------------------------
  // eta-phi TRD-TOF matching efficiency for negative tracks
  pad = ((TVirtualPad*)l->At(4)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  rangeEtaPhi->Draw();
  
  TH2D* hTOFeffNeg = (hTOFrefNeg ? (TH2D*)hTOFrefNeg->Clone("hTOFeffNeg") : 0x0);
  if(hTOFeffNeg) {
    hTOFeffNeg->Reset();
    hTOFeffNeg->SetStats(kFALSE);
    hTOFeffNeg->Divide(hTOFrefNeg, hTRDrefNeg);
    hTOFeffNeg->SetMaximum(1.0);
    hTOFeffNeg->Draw("samecolz");
    lat->DrawLatex(-0.9, 3.3, "TRD-TOF matching for negative tracks");
    DrawTRDGrid();
  }
  
  if(hTRDrefPos) delete hTRDrefPos; if(hTPCrefPos) delete hTPCrefPos; if(hTOFrefPos) delete hTOFrefPos;
  if(hTRDrefNeg) delete hTRDrefNeg; if(hTPCrefNeg) delete hTPCrefNeg; if(hTOFrefNeg) delete hTOFrefNeg;
  
  // switch to the Pt cf container
  matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), +1.0, +1.0);  
  TH1F* hTRDEffPtPosAll = EfficiencyFromPhiPt(matchPt, 0, 6, 1, 0);
  TH1F* hTOFEffPtPosAll = EfficiencyFromPhiPt(matchPt, 0, 6, 2, 1);
  TH1F* hTRDEffPtPosTrk4 = EfficiencyFromPhiPt(matchPt, 4, 4, 1, 0);
  TH1F* hTOFEffPtPosTrk4 = EfficiencyFromPhiPt(matchPt, 4, 4, 2, 1);
  TH1F* hTRDEffPtPosTrk5 = EfficiencyFromPhiPt(matchPt, 5, 5, 1, 0);
  TH1F* hTOFEffPtPosTrk5 = EfficiencyFromPhiPt(matchPt, 5, 5, 2, 1);
  TH1F* hTRDEffPtPosTrk6 = EfficiencyFromPhiPt(matchPt, 6, 6, 1, 0);
  TH1F* hTOFEffPtPosTrk6 = EfficiencyFromPhiPt(matchPt, 6, 6, 2, 1);
  
  matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);  
  TH1F* hTRDEffPtNegAll = EfficiencyFromPhiPt(matchPt, 0, 6, 1, 0);
  TH1F* hTOFEffPtNegAll = EfficiencyFromPhiPt(matchPt, 0, 6, 2, 1);
  TH1F* hTRDEffPtNegTrk4 = EfficiencyFromPhiPt(matchPt, 4, 4, 1, 0);
  TH1F* hTOFEffPtNegTrk4 = EfficiencyFromPhiPt(matchPt, 4, 4, 2, 1);
  TH1F* hTRDEffPtNegTrk5 = EfficiencyFromPhiPt(matchPt, 5, 5, 1, 0);
  TH1F* hTOFEffPtNegTrk5 = EfficiencyFromPhiPt(matchPt, 5, 5, 2, 1);
  TH1F* hTRDEffPtNegTrk6 = EfficiencyFromPhiPt(matchPt, 6, 6, 1, 0);
  TH1F* hTOFEffPtNegTrk6 = EfficiencyFromPhiPt(matchPt, 6, 6, 2, 1);
  matchPt->SetRangeUser(matchPt->GetVar(fgkVarNames[kTrackCharge]), -1.0, +1.0);  
  
  
  TF1* funcConst = new TF1("constFunc", "[0]", 1.0, 3.0);
  if(trendValues) {
    if(hTRDEffPtPosAll && hTRDEffPtPosAll->Integral()>0.1) {
      hTRDEffPtPosAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[0] = funcConst->GetParameter(0);
      trendValues[1] = funcConst->GetParError(0);
    }
  }
  if(trendValues) { 
    if(hTRDEffPtNegAll && hTRDEffPtNegAll->Integral()>0.1) {
      hTRDEffPtNegAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[2] = funcConst->GetParameter(0);
      trendValues[3] = funcConst->GetParError(0);
    }
  }
  if(trendValues) { 
    if(hTOFEffPtPosAll && hTOFEffPtPosAll->Integral()>0.1) {
      hTOFEffPtPosAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[4] = funcConst->GetParameter(0);
      trendValues[5] = funcConst->GetParError(0);
    }
  }
  if(trendValues) { 
    if(hTOFEffPtNegAll && hTOFEffPtNegAll->Integral()>0.1) {
      hTOFEffPtNegAll->Fit(funcConst, "Q0ME", "goff", 1.0, 3.0);
      trendValues[6] = funcConst->GetParameter(0);
      trendValues[7] = funcConst->GetParError(0);
    }
  }
  
  //---------------------------------------------------------
  // TPC-TRD matching efficiency vs pt
  pad = ((TVirtualPad*)l->At(6)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  if(gROOT->FindObject("rangeEffPt2")) delete gROOT->FindObject("rangeEffPt2");
  TH2F* rangeEffPt=new TH2F("rangeEffPt2", "",10,0.,10.,10,0.,1.4);
  rangeEffPt->SetStats(kFALSE);
  SetStyle(rangeEffPt->GetXaxis(), "p_{T} [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEffPt->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  rangeEffPt->Draw();
  lat->DrawLatex(0.2, 1.44, "TPC-TRD matching efficiency");
  //++++++++++++++++++
  TLine line;
  line.SetLineStyle(2);
  line.SetLineWidth(2);
  line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.7, rangeEffPt->GetXaxis()->GetXmax(), 0.7);
  line.DrawLine(rangeEffPt->GetXaxis()->GetXmin(), 0.9, rangeEffPt->GetXaxis()->GetXmax(), 0.9);
  TLegend* leg=new TLegend(0.2, 0.7, 0.7, 0.89);
  leg->SetNColumns(2);
  leg->SetMargin(0.15);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  SetStyle(hTRDEffPtPosAll, 1, kRed, 1, 24, kRed, 1);
  SetStyle(hTRDEffPtNegAll, 1, kBlue, 1, 24, kBlue, 1);
  SetStyle(hTRDEffPtPosTrk4, 1, kRed, 1, 25, kRed, 1);
  SetStyle(hTRDEffPtNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
  SetStyle(hTRDEffPtPosTrk5, 1, kRed, 1, 26, kRed, 1);
  SetStyle(hTRDEffPtNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
  SetStyle(hTRDEffPtPosTrk6, 1, kRed, 1, 27, kRed, 1);
  SetStyle(hTRDEffPtNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
  if(hTRDEffPtPosAll) {hTRDEffPtPosAll->Draw("same"); leg->AddEntry(hTRDEffPtPosAll, "pos. (#geq 1 tracklet)", "p");}
  if(hTRDEffPtNegAll) {hTRDEffPtNegAll->Draw("same"); leg->AddEntry(hTRDEffPtNegAll, "neg. (#geq 1 tracklet)", "p");}
  hTRDEffPtPosTrk4->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk4, "pos. (4 tracklets)", "p");
  hTRDEffPtNegTrk4->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk4, "neg. (4 tracklets)", "p");
  hTRDEffPtPosTrk5->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk5, "pos. (5 tracklets)", "p");
  hTRDEffPtNegTrk5->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk5, "neg. (5 tracklets)", "p");
  hTRDEffPtPosTrk6->Draw("same"); leg->AddEntry(hTRDEffPtPosTrk6, "pos. (6 tracklets)", "p");
  hTRDEffPtNegTrk6->Draw("same"); leg->AddEntry(hTRDEffPtNegTrk6, "neg. (6 tracklets)", "p");
  leg->Draw();
  
  //---------------------------------------------------------
  // TRD-TOF matching efficiency vs pt
  pad = ((TVirtualPad*)l->At(7)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  rangeEffPt->Draw();
  lat->DrawLatex(0.2, 1.44, "TRD-TOF matching efficiency");
  SetStyle(hTOFEffPtPosAll, 1, kRed, 1, 24, kRed, 1);
  SetStyle(hTOFEffPtPosTrk4, 1, kRed, 1, 25, kRed, 1);
  SetStyle(hTOFEffPtPosTrk5, 1, kRed, 1, 26, kRed, 1);
  SetStyle(hTOFEffPtPosTrk6, 1, kRed, 1, 27, kRed, 1);
  SetStyle(hTOFEffPtNegAll, 1, kBlue, 1, 24, kBlue, 1);
  SetStyle(hTOFEffPtNegTrk4, 1, kBlue, 1, 25, kBlue, 1);
  SetStyle(hTOFEffPtNegTrk5, 1, kBlue, 1, 26, kBlue, 1);
  SetStyle(hTOFEffPtNegTrk6, 1, kBlue, 1, 27, kBlue, 1);
  if(hTOFEffPtPosAll) hTOFEffPtPosAll->Draw("same"); 
  hTOFEffPtPosTrk4->Draw("same"); 
  hTOFEffPtPosTrk5->Draw("same"); 
  hTOFEffPtPosTrk6->Draw("same"); 
  if(hTOFEffPtNegAll) hTOFEffPtNegAll->Draw("same"); 
  hTOFEffPtNegTrk4->Draw("same"); 
  hTOFEffPtNegTrk5->Draw("same"); 
  hTOFEffPtNegTrk6->Draw("same");  
    
  //-----------------------------------------------------
  // <ntracklets> vs (phi,eta)
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  rangeEtaPhi->Draw();
  lat->DrawLatex(-0.9, 3.3, "TRD <N_{tracklets}>");
  
  TH3D* hNtracklets = (TH3D*)matchPhiEta->Project(kTRD, matchPhiEta->GetVar(fgkVarNames[kTrackPhiTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackEtaTRD]), matchPhiEta->GetVar(fgkVarNames[kTrackTrdTracklets]));
  
  TProfile2D* hNtrackletsProf = hNtracklets->Project3DProfile();
  delete hNtracklets;
  if(hNtrackletsProf) {
    hNtrackletsProf->SetStats(kFALSE);
    hNtrackletsProf->SetMinimum(0.);
    hNtrackletsProf->SetMaximum(6.);
    hNtrackletsProf->Draw("samecolz");
    DrawTRDGrid();
  }
  
  // calculate the trend value for tracklets/track
  TH2D* hNtrackletsVsP = (TH2D*)matchPt->Project(kTRD, matchPt->GetVar(fgkVarNames[kTrackPt]), matchPt->GetVar(fgkVarNames[kTrackTrdTracklets]));
  if(trendValues &&  hNtrackletsVsP && hNtrackletsVsP->GetEntries()>0.1) {
    TProfile* hNtrackletsVsPprof = hNtrackletsVsP->ProfileX("hNtrackletsVsPprof");
    hNtrackletsVsPprof->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[8] = funcConst->GetParameter(0);
    trendValues[9] = funcConst->GetParError(0);
    delete hNtrackletsVsP;
  }
    
  //--------------------------------------------------------------
  // Nclusters per TRD track vs momentum
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.12);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  
  if(gROOT->FindObject("rangeNclsP")) delete gROOT->FindObject("rangeNclsP");
  TH2F* rangeNclsP = new TH2F("rangeNclsP", "", 10, 0.0, 11.99, 10, 0.0, 199.0);
  SetStyle(rangeNclsP->GetXaxis(), "p [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeNclsP->GetYaxis(), "#clusters", 0.07, 0.8, kTRUE, 0.05);
  rangeNclsP->SetStats(kFALSE);
  rangeNclsP->Draw();
  lat->DrawLatex(1.0, 205., "TRD Clusters / track");
  
  TH2D* hNclsVsP = (TH2D*)clustersCF->Project(0, clustersCF->GetVar(fgkVarNames[kTrackP]), clustersCF->GetVar(fgkVarNames[kTrackTrdClusters]));
  if(hNclsVsP) {
    hNclsVsP->SetStats(kFALSE);
    hNclsVsP->Draw("samecolz");
  }
    
  if(trendValues && hNclsVsP && hNclsVsP->GetEntries()>10) {
    TProfile* hNclsVsPprof = hNclsVsP->ProfileX("hNclsVsPprof");
    hNclsVsPprof->Fit(funcConst, "QME0", "goff", 1.0, 3.0);
    trendValues[10] = funcConst->GetParameter(0);
    trendValues[11] = funcConst->GetParError(0);
  }
  
  //--------------------------------------------------------------
  // TRD-TPC and TOF-TRD matching efficiency vs bunch crossing
  pad = ((TVirtualPad*)l->At(8)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.02);
  pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);

  TH1F* hTRDEffBC = EfficiencyFromPhiPt(bcCF, -1, -1, 1, 0, kEventBC);
  TH1F* hTOFEffBC = EfficiencyFromPhiPt(bcCF, -1, -1, 2, 1, kEventBC);
   
  if(gROOT->FindObject("rangeBC")) delete gROOT->FindObject("rangeBC");
  TH2F* rangeBC = new TH2F("rangeBC", "", 10, -0.5, 3499.5, 10, 0.0, 1.4);
  rangeBC->SetStats(kFALSE);
  SetStyle(rangeBC->GetXaxis(), "Bunch crossing", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeBC->GetYaxis(), "efficiency", 0.07, 0.8, kTRUE, 0.05);
  rangeBC->Draw();
  
  TLegend* legBC=new TLegend(0.8, 0.7, 0.95, 0.89);
  legBC->SetBorderSize(0);
  legBC->SetMargin(0.15);
  legBC->SetFillColor(0);
  if(hTRDEffBC) {
    hTRDEffBC->SetStats(kFALSE);
    SetStyle(hTRDEffBC, 1, kRed, 2, 24, kRed, 1); legBC->AddEntry(hTRDEffBC, "TPC-TRD", "p");
    SetStyle(hTOFEffBC, 1, kBlue, 2, 24, kBlue, 1); legBC->AddEntry(hTOFEffBC, "TRD-TOF", "p");
    hTRDEffBC->Draw("same");
    hTOFEffBC->Draw("same");
    legBC->Draw();
    lat->DrawLatex(200., 1.44, "Matching efficiency at 1<p_{T}<3 GeV/c");
  }
  
  delete funcConst;
}



//_________________________________________________________________
void AliTRDcheckESD::PlotPidSummaryFromCF(Double_t* trendValues, const Char_t* /*triggerName*/, Bool_t /*useIsolatedBC*/, Bool_t /*cutTOFbc*/) {
  //
  // PID summary
  //
  //  trendValues will be filled with trending variables
  //  trendValues[12] : PH plateau height from slices times 0.002
  //  trendValues[13] : statistical error of trendValues[12]
  //  trendValues[14] : PH slope from slices times 0.002
  //  trendValues[15] : statistical error of trendValues[14]
  //  trendValues[16] : Landau MPV of tracklet Qtot distribution at p=1GeV/c times 0.002
  //  trendValues[17] : Landau width of tracklet Qtot distribution at p=1GeV/c times 0.002
  //  trendValues[18] : PH plateau height from slices
  //  trendValues[19] : statistical error of trendValues[19]
  //  trendValues[20] : PH slope from slices
  //  trendValues[21] : statistical error of trendValues[20]
  //  trendValues[22] : Landau MPV of tracklet Qtot distribution at p=1GeV/c
  //  trendValues[23] : Landau width of tracklet Qtot distribution at p=1GeV/c
  //
  if(!fHistos) return;
  AliCFContainer* qtotCF = (AliCFContainer*)fHistos->FindObject("QtotCF");
  if(!qtotCF) return;
  AliCFContainer* phCF = (AliCFContainer*)fHistos->FindObject("PulseHeightCF");
  if(!phCF) return;
  AliCFContainer* centCF = (AliCFContainer*)fHistos->FindObject("CentralityCF");
  if(!centCF) return;   
  
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.07);
  lat->SetTextColor(2);
  lat->SetTextFont(42);
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
    
  if(gROOT->FindObject("rangeEtaPhi2")) delete gROOT->FindObject("rangeEtaPhi2");
  TH2F* rangeEtaPhi = new TH2F("rangeEtaPhi2", "", 10, -0.99, +0.99, 10, -3.15, 3.15);
  SetStyle(rangeEtaPhi->GetXaxis(), "#eta", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeEtaPhi->GetYaxis(), "detector #varphi", 0.07, 0.8, kTRUE, 0.05);
  rangeEtaPhi->SetStats(kFALSE);  
  
  // eta-phi distr. for <Qtot> in layer 0
  TVirtualPad* pad;
  TProfile2D* hProf2D;
  TH2D* hQtotMPV;
  TH1D* tempProj1D;
  TH1D* hqtot = (TH1D*)qtotCF->Project(0, qtotCF->GetVar(fgkVarNames[kTrackletQtot]));
  qtotCF->SetRangeUser(qtotCF->GetVar(fgkVarNames[kTrackletQtot]), 0.0, 3000.);
  for(Int_t iLayer=0; iLayer<6; ++iLayer) {
    pad = ((TVirtualPad*)l->At((iLayer<3 ? iLayer*3 : (iLayer-3)*3+1))); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.1); pad->SetBottomMargin(0.15);
    pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
    rangeEtaPhi->Draw();
    qtotCF->SetRangeUser(qtotCF->GetVar(fgkVarNames[kTrackletLayer]), Double_t(iLayer), Double_t(iLayer));
    TH3D* hQtotEtaPhi = (TH3D*)qtotCF->Project(0, qtotCF->GetVar(fgkVarNames[kTrackPhiTRD]), qtotCF->GetVar(fgkVarNames[kTrackEtaTRD]), qtotCF->GetVar(fgkVarNames[kTrackletQtot]));  
    hQtotEtaPhi->Rebin3D(2,2,1);
    //hProf2D = (hQtotEtaPhi ? hQtotEtaPhi->Project3DProfile() : 0x0);
    
    hQtotMPV = new TH2D(Form("QtotMPV_layer%d",iLayer),"",hQtotEtaPhi->GetYaxis()->GetNbins(),
                        hQtotEtaPhi->GetYaxis()->GetXmin(), hQtotEtaPhi->GetYaxis()->GetXmax(), 
                        hQtotEtaPhi->GetXaxis()->GetNbins(), hQtotEtaPhi->GetXaxis()->GetXmin(), hQtotEtaPhi->GetXaxis()->GetXmax());
    for(Int_t iphi=1; iphi<hQtotEtaPhi->GetYaxis()->GetNbins(); ++iphi) {
       for(Int_t ieta=1; ieta<hQtotEtaPhi->GetXaxis()->GetNbins(); ++ieta) {
          tempProj1D = hQtotEtaPhi->ProjectionZ("tempProj",ieta,ieta,iphi,iphi);
          Int_t maxBin = tempProj1D->GetMaximumBin();
          if(maxBin<3) continue;
          if(maxBin>tempProj1D->GetXaxis()->GetNbins()-3) continue;
          Double_t mpv = tempProj1D->GetBinContent(maxBin)*tempProj1D->GetBinCenter(maxBin)+
                                      tempProj1D->GetBinContent(maxBin-1)*tempProj1D->GetBinCenter(maxBin-1)+
                                      tempProj1D->GetBinContent(maxBin+1)*tempProj1D->GetBinCenter(maxBin+1)+
                                      tempProj1D->GetBinContent(maxBin-2)*tempProj1D->GetBinCenter(maxBin-2)+
                                      tempProj1D->GetBinContent(maxBin+2)*tempProj1D->GetBinCenter(maxBin+2);
          Double_t norm = tempProj1D->GetBinContent(maxBin-2)+tempProj1D->GetBinContent(maxBin-1)+tempProj1D->GetBinContent(maxBin)+tempProj1D->GetBinContent(maxBin+1)+tempProj1D->GetBinContent(maxBin+2);
          if(norm>0.5) mpv = mpv / (norm);
          delete tempProj1D;
          hQtotMPV->SetBinContent(iphi,ieta,mpv);
       }
   }
   if(hQtotMPV) {
      hQtotMPV->SetStats(kFALSE);
      hQtotMPV->SetMinimum(500.);
      hQtotMPV->SetMaximum((hqtot->GetMean()<10 ? 4.0 : 1500.));
      hQtotMPV->Draw("samecolz");
   }
    
    if(hQtotEtaPhi) delete hQtotEtaPhi;    
    if(hProf2D) {
      hProf2D->SetName(Form("Qtot_layer%d",iLayer));
      hProf2D->SetStats(kFALSE);
      hProf2D->SetMinimum(500.);
      hProf2D->SetMaximum((hqtot->GetMean()<10 ? 4.0 : 1500.));
      hProf2D->Draw("samecolz");
    }
    lat->DrawLatex(-0.9, 3.3, Form("TRD Q_{tot} MPV Layer %d", iLayer));
    DrawTRDGrid();
  }
  qtotCF->SetRangeUser(qtotCF->GetVar(fgkVarNames[kTrackletLayer]), 0.0, 5.0);
  qtotCF->SetRangeUser(qtotCF->GetVar(fgkVarNames[kTrackletQtot]), 0.0, 10000.);
  // PH versus slice number
  pad = ((TVirtualPad*)l->At(2)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.03); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  
  if(gROOT->FindObject("rangePHslice")) delete gROOT->FindObject("rangePHslice");
  TH2F* rangePHslice=new TH2F("rangePHslice", "", 8, -0.5, 7.5, 10, 0.0, (hqtot->GetMean()<10.0 ? 6.0 : 3000.));
  rangePHslice->SetStats(kFALSE);
  SetStyle(rangePHslice->GetXaxis(), "slice", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangePHslice->GetYaxis(), "PH", 0.07, 0.8, kTRUE, 0.05);
  rangePHslice->Draw();
  
  TF1* funcPol1 = new TF1("funcPol1", "[0]+[1]*x", 2.9, 6.4);
  
  TH2D* hPH = (TH2D*)phCF->Project(0, phCF->GetVar(fgkVarNames[kTrackletSlice]), phCF->GetVar(fgkVarNames[kTrackletPHslice]));
  TH1D* hSliceErr = new TH1D(Form("hSliceErr%f", gRandom->Rndm()), "", hPH->GetXaxis()->GetNbins(), hPH->GetXaxis()->GetXbins()->GetArray());
  TH1D* hLandauFit = Proj2D(hPH, hSliceErr);
  hPH->SetStats(kFALSE);
  hPH->Draw("samecolz");
  const Double_t kQx = 0.002;
  if(trendValues) {
    hSliceErr->Fit(funcPol1, "QME0", "goff", 2.9, 6.4);
    trendValues[12] = kQx*funcPol1->GetParameter(0);  // PH plateau
    trendValues[13] = kQx*funcPol1->GetParError(0);   // PH plateau
    trendValues[14] = kQx*funcPol1->GetParameter(1);  // PH slope
    trendValues[15] = kQx*funcPol1->GetParError(1);   // PH slope
    trendValues[18] = funcPol1->GetParameter(0);  // PH plateau
    trendValues[19] = funcPol1->GetParError(0);   // PH plateau
    trendValues[20] = funcPol1->GetParameter(1);  // PH slope
    trendValues[21] = funcPol1->GetParError(1);   // PH slope
  }
  hLandauFit->SetLineWidth(2);
  hLandauFit->SetLineStyle(2);
  hLandauFit->Draw("same");
  
  delete funcPol1; delete hSliceErr;
  // Qtot vs P
  pad = ((TVirtualPad*)l->At(5)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.03); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  
  if(gROOT->FindObject("rangeQtotP")) delete gROOT->FindObject("rangeQtotP");
  TH2F* rangeQtotP = new TH2F("rangeQtotP", "", 10, 0.0, 11.99, 10, 0.0, (hqtot->GetMean()<10.0 ? 11.99 : 5999.));
  SetStyle(rangeQtotP->GetXaxis(), "P [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtotP->GetYaxis(), "Q_{tot}", 0.07, 0.8, kTRUE, 0.05);
  rangeQtotP->SetStats(kFALSE);
  rangeQtotP->Draw();
  
  Int_t pVar = centCF->GetVar(fgkVarNames[kTrackP]);
  if(pVar<0) pVar = centCF->GetVar(fgkVarNames[kTrackletP]);
  TH2D* hQtotP = (TH2D*)centCF->Project(0, pVar, centCF->GetVar(fgkVarNames[kTrackletQtot]));
  TH1D* mpvErr=new TH1D("mpvErr", "Landau MPV error vs. P", hQtotP->GetXaxis()->GetNbins(), hQtotP->GetXaxis()->GetXbins()->GetArray());  
  TH1D* widthErr=new TH1D("widthErr", "Landau width error vs. P", hQtotP->GetXaxis()->GetNbins(), hQtotP->GetXaxis()->GetXbins()->GetArray());
  TH1D* landauChi2=new TH1D("landauChi2", "Landau fit #chi^{2} vs. P", hQtotP->GetXaxis()->GetNbins(), hQtotP->GetXaxis()->GetXbins()->GetArray());  
  if(hQtotP)
    for(Int_t i=1; i<=hQtotP->GetXaxis()->GetNbins(); ++i) 
      hQtotP->SetBinContent(i, 1, 0.0);  
  TH1D* hQtotProj = (hQtotP ? Proj2D(hQtotP, mpvErr, widthErr, landauChi2) : 0x0);
  //landauChi2->Scale(0.001);
  if(hQtotProj) SetStyle(hQtotProj, 2, kBlue, 2, 1, kBlue, 1);
  if(trendValues && hQtotProj && hQtotProj->GetEntries()>2) {
    trendValues[16] = kQx*hQtotProj->GetBinContent(hQtotProj->FindBin(1.0));   // Landau MPV at 1GeV/c
    trendValues[17] = kQx*hQtotProj->GetBinError(hQtotProj->FindBin(1.0));     // Landau width at 1 GeV/c
    trendValues[22] = hQtotProj->GetBinContent(hQtotProj->FindBin(1.0));   // Landau MPV at 1GeV/c
    trendValues[23] = hQtotProj->GetBinError(hQtotProj->FindBin(1.0));     // Landau width at 1 GeV/c
  }
  if(hQtotP) {
    hQtotP->SetStats(kFALSE);
    hQtotP->Draw("samecolz");
    hQtotProj->Draw("same");
  }
  
  // Qtot vs P (fit results)
  pad = ((TVirtualPad*)l->At(8)); pad->cd();
  pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
  pad->SetTopMargin(0.03); pad->SetBottomMargin(0.15);
  pad->SetGridx(kFALSE); pad->SetGridy(kFALSE);
  pad->SetLogz();
  
  if(gROOT->FindObject("rangeQtotPfit")) delete gROOT->FindObject("rangeQtotPfit");
  TH2F* rangeQtotPfit = new TH2F("rangeQtotPfit", "", 100, 0.0, 11.99, 100, 0.0, (hqtot->GetMean()<10.0 ? 6.0 : 3999.));
  SetStyle(rangeQtotPfit->GetXaxis(), "P [GeV/c]", 0.07, 0.8, kTRUE, 0.05);
  SetStyle(rangeQtotPfit->GetYaxis(), "Q_{tot}", 0.07, 0.8, kTRUE, 0.05);
  rangeQtotPfit->SetStats(kFALSE);
  rangeQtotPfit->Draw();
  
  if(mpvErr) SetStyle(mpvErr, 1, kBlue, 2, 1, kBlue, 1);
  if(widthErr) SetStyle(widthErr, 2, kRed, 2, 1, kRed, 1);
  if(mpvErr) {
    mpvErr->SetStats(kFALSE);
    mpvErr->Draw("same");
  }
  if(widthErr) {
    widthErr->SetStats(kFALSE);
    widthErr->Draw("same");
  }
  TLegend* leg=new TLegend(0.2,0.6,0.5,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry("mpvErr","Landau MPV","l");
  leg->AddEntry("widthErr","Landau width","l");
  leg->Draw();
}


//__________________________________________________________________________________________________
void AliTRDcheckESD::PlotOtherSummaryFromCF(Double_t* /*trendValues*/) {
  //
  // Plot additional QA 
  //
  if(!fHistos) return;
  AliCFContainer* matchPhiEta = (AliCFContainer*)fHistos->FindObject("MatchingPhiEta");
  AliCFContainer* trdChi2 = (AliCFContainer*)fHistos->FindObject("trdChi2");
  AliCFContainer* trdBudget = (AliCFContainer*)fHistos->FindObject("trdBudget");
  AliCFContainer* ploss = (AliCFContainer*)fHistos->FindObject("Ploss");
  AliCFContainer* clusters = (AliCFContainer*)fHistos->FindObject("clustersPerTracklet");  
  AliCFContainer* clsRows = (AliCFContainer*)fHistos->FindObject("clustersVsRows");
    
  TLatex *lat=new TLatex();
  lat->SetTextSize(0.06);
  lat->SetTextColor(2);
  lat->SetNDC();
  lat->SetTextFont(42);
  TCanvas* c1 = new TCanvas("ESDsummary", "ESD summary 1", 1600, 1200);
  gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
  gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
  gPad->Divide(3,3,0.,0.);
  TList* l=gPad->GetListOfPrimitives();
  TVirtualPad* pad=0x0;
  
  // matching as a function of trigger class
  if(matchPhiEta) {
    matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);
    TH1F* hTRDEffTriggerNeg = EfficiencyFromPhiPt(matchPhiEta, -1, -1, 1, 0, kEventTrigger);
    matchPhiEta->SetRangeUser(matchPhiEta->GetVar(fgkVarNames[kTrackCharge]), +1.0, +1.0);
    TH1F* hTRDEffTriggerPos = EfficiencyFromPhiPt(matchPhiEta, -1, -1, 1, 0, kEventTrigger);
      
    pad = ((TVirtualPad*)l->At(0)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.01);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    hTRDEffTriggerNeg->SetStats(kFALSE);
    SetStyle(hTRDEffTriggerNeg->GetYaxis(), "efficiency", 0.06, 1.0, kTRUE, 0.06);
    hTRDEffTriggerNeg->GetXaxis()->SetRange(1,fNAssignedTriggers);
    hTRDEffTriggerPos->GetXaxis()->SetRange(1,fNAssignedTriggers);
    SetStyle(hTRDEffTriggerNeg, 1, 2, 2, 20, 2, 1);
    SetStyle(hTRDEffTriggerPos, 1, 4, 2, 20, 4, 1);
    hTRDEffTriggerNeg->Draw();
    hTRDEffTriggerPos->Draw("same");
    TLegend* legEff=new TLegend(0.5,0.5,0.7,0.7);
    legEff->SetFillColor(0);
    legEff->SetBorderSize(0);
    legEff->AddEntry(hTRDEffTriggerNeg, "negatives", "l");
    legEff->AddEntry(hTRDEffTriggerPos, "positives", "l");
    legEff->Draw();
    lat->DrawLatex(0.2, 0.95, "TPC-TRD matching efficiency");
  }
  
  if(trdChi2) {
    // Track TRD chi2 vs (eta,phi)
    TH3D* trdChi23D = (TH3D*)trdChi2->Project(0, trdChi2->GetVar(fgkVarNames[kTrackEtaTRD]), 
		                               trdChi2->GetVar(fgkVarNames[kTrackPhiTRD]),
				  	       trdChi2->GetVar(fgkVarNames[kTrackTrdChi2]));
    trdChi23D->SetName("trdChi23D");
    TProfile2D* prof2DChi2 = trdChi23D->Project3DProfile("yx");
    prof2DChi2->SetName("prof2DChi2");
    pad = ((TVirtualPad*)l->At(3)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    prof2DChi2->SetStats(kFALSE);
    prof2DChi2->SetTitle("");
    SetStyle(prof2DChi2->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(prof2DChi2->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
    prof2DChi2->SetMaximum(2.9);
    prof2DChi2->Draw("colz");
    lat->DrawLatex(0.2, 0.95, "TRD #chi^{2}");
    DrawTRDGrid();
  
    // Track TRD chi2 vs pt and charge
    trdChi2->SetRangeUser(trdChi2->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);
    TH2D* trdChi2VsPtNeg = (TH2D*)trdChi2->Project(0, trdChi2->GetVar(fgkVarNames[kTrackPt]),
                                                    trdChi2->GetVar(fgkVarNames[kTrackTrdChi2]));
    trdChi2VsPtNeg->SetName("trdChi2VsPtNeg");
    TProfile* trdChi2VsPtNegProf = trdChi2VsPtNeg->ProfileX();
    trdChi2VsPtNegProf->SetName("trdChi2VsPtNegProf");
    trdChi2->SetRangeUser(trdChi2->GetVar(fgkVarNames[kTrackCharge]), 1.0, 1.0);
    TH2D* trdChi2VsPtPos = (TH2D*)trdChi2->Project(0, trdChi2->GetVar(fgkVarNames[kTrackPt]),
                                                    trdChi2->GetVar(fgkVarNames[kTrackTrdChi2]));
    trdChi2VsPtPos->SetName("trdChi2VsPtPos");
    TProfile* trdChi2VsPtPosProf = trdChi2VsPtPos->ProfileX();
    trdChi2VsPtPosProf->SetName("trdChi2VsPtPosProf");
    pad = ((TVirtualPad*)l->At(6)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.01);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    trdChi2VsPtNegProf->SetStats(kFALSE);
    trdChi2VsPtNegProf->SetTitle("");
    SetStyle(trdChi2VsPtNegProf->GetXaxis(), "p_{T} (GeV/c)", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(trdChi2VsPtNegProf->GetYaxis(), "<TRD #chi^{2}>", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(trdChi2VsPtNegProf, 1, 2, 2, 20, 2, 1);
    SetStyle(trdChi2VsPtPosProf, 1, 4, 2, 20, 4, 1);
    trdChi2VsPtNegProf->Draw();
    trdChi2VsPtPosProf->Draw("same");
    lat->DrawLatex(0.2, 0.95, "TRD #chi^{2}");
  }
  
  if(trdBudget) {
    // Track TRD budget vs (eta,phi)
    TH3D* trdBudget3D = (TH3D*)trdBudget->Project(0, trdBudget->GetVar(fgkVarNames[kTrackEtaTRD]), 
		                                   trdBudget->GetVar(fgkVarNames[kTrackPhiTRD]),
						   trdBudget->GetVar(fgkVarNames[kTrackTRDBudget]));
    trdBudget3D->SetName("trdBudget3D");
    TProfile2D* prof2DBudget = trdBudget3D->Project3DProfile("yx");
    prof2DBudget->SetName("prof2DBudget");
    pad = ((TVirtualPad*)l->At(1)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    prof2DBudget->SetStats(kFALSE);
    prof2DBudget->SetTitle("");
    SetStyle(prof2DBudget->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(prof2DBudget->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
    prof2DBudget->Draw("colz");
    lat->DrawLatex(0.2, 0.95, "TRD budget");
    DrawTRDGrid();
  }
  
  if(ploss) {
    // momentum loss
    ploss->SetRangeUser(ploss->GetVar(fgkVarNames[kTrackCharge]), -1.0, -1.0);
    TH2D* plossLayerNeg = (TH2D*)ploss->Project(0, ploss->GetVar(fgkVarNames[kTrackletLayer]), 
		                                   ploss->GetVar(fgkVarNames[kTrackPlossTRDlayer]));
    plossLayerNeg->SetName("plossLayerNeg");
    TProfile* plossLayerNegProf = plossLayerNeg->ProfileX();
    plossLayerNegProf->SetName("plossLayerNegProf");
    ploss->SetRangeUser(ploss->GetVar(fgkVarNames[kTrackCharge]), +1.0, +1.0);
    TH2D* plossLayerPos = (TH2D*)ploss->Project(0, ploss->GetVar(fgkVarNames[kTrackletLayer]), 
		                                   ploss->GetVar(fgkVarNames[kTrackPlossTRDlayer]));
    plossLayerPos->SetName("plossLayerPos");
    TProfile* plossLayerPosProf = plossLayerPos->ProfileX();
    plossLayerPosProf->SetName("plossLayerPosProf");
    ploss->SetRangeUser(ploss->GetVar(fgkVarNames[kTrackCharge]), -1.5, +1.5);
    
    ploss->SetRangeUser(ploss->GetVar(fgkVarNames[kTrackletLayer]), 0.0, 0.0);
    TH3D* ploss3Dl0 = (TH3D*)ploss->Project(0, ploss->GetVar(fgkVarNames[kTrackEtaTRD]), 
		                             ploss->GetVar(fgkVarNames[kTrackPhiTRD]), 
		                             ploss->GetVar(fgkVarNames[kTrackPlossTRDlayer]));
    ploss3Dl0->SetName("ploss3Dl0");
    TProfile2D* plossEtaPhiL0Prof = ploss3Dl0->Project3DProfile("yx");
    plossEtaPhiL0Prof->SetName("plossEtaPhiL0Prof");
    ploss->SetRangeUser(ploss->GetVar(fgkVarNames[kTrackletLayer]), 5.0, 5.0);
    TH3D* ploss3Dl5 = (TH3D*)ploss->Project(0, ploss->GetVar(fgkVarNames[kTrackEtaTRD]), 
		                             ploss->GetVar(fgkVarNames[kTrackPhiTRD]), 
		                             ploss->GetVar(fgkVarNames[kTrackPlossTRDlayer]));
    ploss3Dl5->SetName("ploss3Dl5");
    TProfile2D* plossEtaPhiL5Prof = ploss3Dl5->Project3DProfile("yx");
    plossEtaPhiL5Prof->SetName("plossEtaPhiL5Prof");
    pad = ((TVirtualPad*)l->At(4)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    plossEtaPhiL0Prof->SetStats(kFALSE);
    plossEtaPhiL0Prof->SetTitle("");
    SetStyle(plossEtaPhiL0Prof->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(plossEtaPhiL0Prof->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
    plossEtaPhiL0Prof->SetMaximum(80.0);
    plossEtaPhiL0Prof->SetMinimum(-20.0);
    plossEtaPhiL0Prof->Draw("colz");
    lat->DrawLatex(0.2, 0.95, "P_{loss} at layer 0");
    DrawTRDGrid();
    pad = ((TVirtualPad*)l->At(7)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    plossEtaPhiL5Prof->SetStats(kFALSE);
    plossEtaPhiL5Prof->SetTitle("");
    SetStyle(plossEtaPhiL5Prof->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(plossEtaPhiL5Prof->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
    plossEtaPhiL5Prof->SetMaximum(80.0);
    plossEtaPhiL5Prof->SetMinimum(-20.0);
    plossEtaPhiL5Prof->Draw("colz");
    lat->DrawLatex(0.2, 0.95, "P_{loss} at layer 5");
    DrawTRDGrid();
    pad = ((TVirtualPad*)l->At(2)); pad->cd();
    pad->SetLeftMargin(0.15); pad->SetRightMargin(0.01);
    pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
    plossLayerNegProf->SetStats(kFALSE);
    plossLayerNegProf->SetTitle("");
    SetStyle(plossLayerNegProf->GetYaxis(), "#Delta P (MeV/c)", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(plossLayerNegProf->GetXaxis(), "TRD layer", 0.06, 1.0, kTRUE, 0.06);
    SetStyle(plossLayerNegProf, 1, 2, 2, 20, 2, 1);
    SetStyle(plossLayerPosProf, 1, 4, 2, 20, 4, 1);
    plossLayerNegProf->GetYaxis()->SetRangeUser(TMath::Min(plossLayerNegProf->GetMinimum(),plossLayerPosProf->GetMinimum()-5.0),
                                                TMath::Max(plossLayerNegProf->GetMaximum(),plossLayerPosProf->GetMaximum())+5.0);
    plossLayerNegProf->Draw();  
    plossLayerPosProf->Draw("same");
    lat->DrawLatex(0.2, 0.95, "P_{loss} vs layer");
  }
  
  // clusters/tracklet and clusters/crossed rows
  TH3D* clustersEtaPhi[6]={0x0};
  TH3D* clsRowsEtaPhi[6]={0x0};
  for(Int_t il=0;il<6;++il) {
    if(clusters) {
      clusters->SetRangeUser(clusters->GetVar(fgkVarNames[kTrackletLayer]), Double_t(il), Double_t(il));
      clustersEtaPhi[il]=(TH3D*)clusters->Project(0, clusters->GetVar(fgkVarNames[kTrackEtaTRD]),
                                                   clusters->GetVar(fgkVarNames[kTrackPhiTRD]),
                                                   clusters->GetVar(fgkVarNames[kTrackletClusters]));
      clustersEtaPhi[il]->SetName(Form("clustersEtaPhi%d",il));
    }
    if(clsRows) {
      clsRows->SetRangeUser(clsRows->GetVar(fgkVarNames[kTrackletLayer]), Double_t(il), Double_t(il));
      clsRowsEtaPhi[il]=(TH3D*)clsRows->Project(0, clsRows->GetVar(fgkVarNames[kTrackEtaTRD]),
                                                 clsRows->GetVar(fgkVarNames[kTrackPhiTRD]),
                                                 clsRows->GetVar(fgkVarNames[kTrackletClustersVsRows]));
      clsRowsEtaPhi[il]->SetName(Form("clsRowsEtaPhi%d",il));
    }
  }
    
  lat->SetTextSize(0.05);
  Int_t layerPads[6] = {0, 2, 4, 1, 3, 5};
  
  TCanvas* c2=0x0;
  if(clusters) {
    c2 = new TCanvas("ESDsummary2", "ESD summary 2", 1600, 1200);
    gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
    gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
    gPad->Divide(3,2,0.,0.);
    l=gPad->GetListOfPrimitives();
    for(Int_t il=0;il<6;++il) {
      TProfile2D* clustersEtaPhiProf = clustersEtaPhi[il]->Project3DProfile("yx");
      pad = ((TVirtualPad*)l->At(layerPads[il])); pad->cd();
      pad->SetLeftMargin(0.15); pad->SetRightMargin(0.1);
      pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
      clustersEtaPhiProf->SetStats(kFALSE);
      clustersEtaPhiProf->SetTitle("");
      SetStyle(clustersEtaPhiProf->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
      SetStyle(clustersEtaPhiProf->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
      clustersEtaPhiProf->SetMaximum(30.);
      clustersEtaPhiProf->Draw("colz");
      lat->DrawLatex(0.2, 0.95, Form("Clusters / tracklet, layer %d", il));
      DrawTRDGrid();
    }
  }
  
  TCanvas* c3=0x0;
  if(clsRows) {
    c3 = new TCanvas("ESDsummary3", "ESD summary 3", 1600, 1200);
    gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.001);
    gPad->SetLeftMargin(0.001); gPad->SetRightMargin(0.001);
    gPad->Divide(3,2,0.,0.);
    l=gPad->GetListOfPrimitives();
    for(Int_t il=0;il<6;++il) {
      TProfile2D* clsRowsEtaPhiProf = clsRowsEtaPhi[il]->Project3DProfile("yx");
      pad = ((TVirtualPad*)l->At(layerPads[il])); pad->cd();
      pad->SetLeftMargin(0.15); pad->SetRightMargin(0.13);
      pad->SetTopMargin(0.06); pad->SetBottomMargin(0.15);
      clsRowsEtaPhiProf->SetStats(kFALSE);
      clsRowsEtaPhiProf->SetTitle("");
      SetStyle(clsRowsEtaPhiProf->GetXaxis(), "#eta", 0.06, 1.0, kTRUE, 0.06);
      SetStyle(clsRowsEtaPhiProf->GetYaxis(), "#varphi (rad.)", 0.06, 1.0, kTRUE, 0.06);
      clsRowsEtaPhiProf->Draw("colz");
      lat->DrawLatex(0.2, 0.95, Form("Clusters / crossed rows, layer %d", il));
      DrawTRDGrid();
    }
  }
  
  if(matchPhiEta || trdChi2 || trdBudget || ploss) c1->SaveAs("esdSummary1.gif");
  if(clusters) c2->SaveAs("esdSummary2.gif");
  if(clsRows) c3->SaveAs("esdSummary3.gif");
}


//__________________________________________________________________________________________________
void AliTRDcheckESD::DrawTRDGrid() {
  //
  //   Draw a grid of lines showing the TRD supermodule and stack structure in (eta,phi) coordinates.
  //   The canvas on which to draw must already exist.
  //
  TLine line;
  line.SetLineColor(2);
  line.SetLineWidth(1);
  line.SetLineStyle(2);
  for(Int_t i=-9; i<=9; ++i) {
    line.DrawLine(-0.92, 2.0*TMath::Pi()/18.0*i, +0.92, 2.0*TMath::Pi()/18.0*i);
  }
  line.DrawLine(-0.85, -3.15, -0.85, 3.15);
  line.DrawLine(-0.54, -3.15, -0.54, 3.15);
  line.DrawLine(-0.16, -3.15, -0.16, 3.15);
  line.DrawLine(+0.16, -3.15, +0.16, 3.15);
  line.DrawLine(+0.54, -3.15, +0.54, 3.15);
  line.DrawLine(+0.85, -3.15, +0.85, 3.15);
}

//_________________________________________________________________
void AliTRDcheckESD::SetStyle(TH1* hist, 
			      Int_t lineStyle, Int_t lineColor, Int_t lineWidth, 
			      Int_t markerStyle, Int_t markerColor, Int_t markerSize) {
  //
  // Set style settings for histograms
  //
  if(!hist) return;
  hist->SetLineStyle(lineStyle);
  hist->SetLineColor(lineColor);
  hist->SetLineWidth(lineWidth);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerColor(markerColor);
  hist->SetMarkerSize(markerSize);
}

//____________________________________________________________________
void AliTRDcheckESD::SetStyle(TAxis* axis, const Char_t* title, Float_t titleSize, Float_t titleOffset, Bool_t centerTitle, 
                              Float_t labelSize) {
  //
  // Set style settings for axes
  //
  if(!axis) return;
  axis->SetTitle(title);
  axis->SetTitleSize(titleSize);
  axis->SetTitleOffset(titleOffset); 
  axis->CenterTitle(centerTitle);
  axis->SetLabelSize(labelSize);
  axis->SetTitleFont(42);
  axis->SetLabelFont(42);
  axis->SetNdivisions(507);
}

//____________________________________________________________________
void AliTRDcheckESD::FindIsolatedBCs(TH1D* bcHist, Bool_t isIsolated[3500]) {
  //
  // Find the isolated bunch crossings
  //
  Int_t isolationSize = 5;      // number of free bunches in both directions
  for(Int_t bcBin=1; bcBin<=bcHist->GetXaxis()->GetNbins(); ++bcBin) {
    Int_t bc = TMath::Nint(bcHist->GetBinCenter(bcBin));
    if(bc<-0.001 || bc>3499.01) {
      isIsolated[bc] = kFALSE;
      continue;
    }
    Double_t entries = bcHist->GetBinContent(bcBin);
    if(entries<0.001) {
      isIsolated[bc] = kFALSE;
      continue;     // no entries
    }
        
    // check isolation
    isIsolated[bc] = kTRUE;
    for(Int_t ibc = TMath::Max(1,bcBin-isolationSize); ibc<=TMath::Min(3499, bcBin+isolationSize); ++ibc) {
      if(ibc==bcBin) continue;
      if(bcHist->GetBinContent(ibc)>0.01) {
        isIsolated[bc] = kFALSE;
        break;
      }
    }
  }   // end loop over BC bins
  
  cout << "Isolated bunches: " << endl;
  for(Int_t ibc=0; ibc<3500; ++ibc) 
    if(isIsolated[ibc]) cout << "BC #" << ibc << endl; 
}


//__________________________________________________________________________________________________
Int_t AliTRDcheckESD::GetTriggerIndex(const Char_t* name, Bool_t createNew/*=kTRUE*/) {
  //
  //  Return the index of trigger "name" in the trigger histogram.
  //  If the index for this trigger does not exist yet, then assign one if createNew is set to TRUE 
  //
  //cout << "GetTriggerIndex for " << name << endl;
  TH1F* triggerHist = (TH1F*)fHistos->FindObject("hTriggerDefs");
  TString nameStr=name;
  for(Int_t i=1; i<=triggerHist->GetXaxis()->GetNbins(); ++i) {
    if(!nameStr.CompareTo(triggerHist->GetXaxis()->GetBinLabel(i))) {
      //cout << "       index found: " << i << endl;
      return i;
    }
  }
  if(createNew) {
    triggerHist->GetXaxis()->SetBinLabel(fNAssignedTriggers+1, name);
    for(Int_t i=1;i<fHistos->GetEntries();++i) {
      TString objType = fHistos->At(i)->IsA()->GetName();
      if(!objType.Contains("AliCFContainer")) continue;
      AliCFContainer* cf=(AliCFContainer*)fHistos->At(i);
      Int_t trigVar = cf->GetVar(fgkVarNames[kEventTrigger]);
      if(trigVar>=0)
	for(Int_t istep=0;istep<cf->GetNStep();++istep) 
	  cf->GetAxis(trigVar, istep)->SetBinLabel(fNAssignedTriggers+1, name);
    }  // end loop over histograms and CFs
        
    ++fNAssignedTriggers;
    return fNAssignedTriggers+1;
  }
  else {
    return -1;
  }
}

//__________________________________________________________________________________________________
void AliTRDcheckESD::PrintTriggers() const {
  //
  //  Print the available triggers for this run
  //
  if(!fHistos) {
    cout << "Warning in AliTRDcheckESD::PrintTriggers(): No file loaded!" << endl;
    return;
  }
  TH1F* hTriggers = (TH1F*)fHistos->FindObject("hTriggerDefs");
  cout << "Triggers found in this run" << endl;
  cout << "==========================" << endl;
  cout << "Name   Index   Entries    " << endl;
  for(Int_t it=1; it<hTriggers->GetXaxis()->GetNbins(); ++it) {
    if(hTriggers->GetXaxis()->GetBinLabel(it)[0]!='\0') {
      cout << hTriggers->GetXaxis()->GetBinLabel(it) << "  " << hTriggers->GetXaxis()->GetBinCenter(it) << "  " << hTriggers->GetBinContent(it) << endl;
    }
  }
}


//__________________________________________________________________________________________________
Int_t AliTRDcheckESD::GetTriggerCounter(const Char_t* triggerName) const {
  //
  // Get the number of events for a given trigger name
  //
  if(!fHistos) {
    cout << "Warning in AliTRDcheckESD::PrintTriggers(): No file loaded!" << endl;
    return -1;
  }
  TH1F* hTriggers = (TH1F*)fHistos->FindObject("hTriggerDefs");
  Int_t counter = -1;
  for(Int_t it=1; it<hTriggers->GetXaxis()->GetNbins(); ++it) {
    TString trgString = hTriggers->GetXaxis()->GetBinLabel(it);
    if(!trgString.CompareTo(triggerName)) 
      counter = (Int_t)hTriggers->GetBinContent(it);
  }
  if(counter<0) {cout << "AliTRDcheckESD::GetTriggerCounter()  Trigger not found !!";}
  return counter;
}


//__________________________________________________________________________________________________________
Int_t AliTRDcheckESD::GetNAssignedTriggers() {
  //
  // Return the number of assigned triggers
  //
  return fNAssignedTriggers;
}
