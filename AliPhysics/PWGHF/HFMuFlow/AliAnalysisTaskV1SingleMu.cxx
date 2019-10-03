/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskV1SingleMu
/// Analysis task for v1 of single muons in the spectrometer with the scalar prduct method
/// The output is a list of histograms and CF containers.
/// The macro class can run on AODs or ESDs.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Audrey Francisco and Jacopo Margutti from AliAnalysisTaskDimu
//-----------------------------------------------------------------------------

#define AliAnalysisTaskV1SingleMu_cxx

#include "AliAnalysisTaskV1SingleMu.h"

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TArrayD.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TLatex.h"

// STEER includes
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliVHeader.h"
#include "AliAODMCHeader.h"
#include "AliStack.h"
#include "AliEventplane.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliFlowEvent.h"
#include "AliAnalysisTaskZDCEP.h"

// CORRFW includes
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

// PWG includes
#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMuonAnalysisOutput.h"
#include "AliUtilityMuonAncestor.h"

#include <fstream>

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskV1SingleMu) // Class implementation in ROOT context
/// \endcond

//________________________________________________________________________
AliAnalysisTaskV1SingleMu::AliAnalysisTaskV1SingleMu() :
AliAnalysisTaskSE(),
fUtilityMuonAncestor(0x0),
fNPtBins(160),
fNEtaBins(3),
fNSPBins(500),
fNCentBins(50),
fNPhiBins(50),
fHarmonic(1),
fMergeableCollection(0x0),
fZdcCalibOnly(kTRUE),
fSparse(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskV1SingleMu::AliAnalysisTaskV1SingleMu(const char *name) :
AliAnalysisTaskSE(name),
// fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
// fMuonTrackCuts(new AliMuonTrackCuts(cuts)),
fUtilityMuonAncestor(0x0),
fNPtBins(160),
fNEtaBins(3),
fNSPBins(500),
fNCentBins(50),
fNPhiBins(50),
fHarmonic(1),
fMergeableCollection(0x0),
fZdcCalibOnly(kTRUE),
fSparse(0x0)
{
  //
  /// Constructor.
  //
  DefineInput(1,AliFlowEventSimple::Class());
  DefineOutput(1, AliMergeableCollection::Class());
}


//________________________________________________________________________
AliAnalysisTaskV1SingleMu::~AliAnalysisTaskV1SingleMu()
{
  //
  /// Destructor
  //
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fMergeableCollection;
  }
  delete fSparse;
  // delete fAODEvent;
  // delete fESDEvent;
}
//________________________________________________________________________
TObject* AliAnalysisTaskV1SingleMu::GetMergeableObject ( TString identifier, TString objectName )
{
  /// Get mergeable object

  TObject* obj = fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
  if ( obj ) return obj;


  if ( objectName == "MuSparse" ) obj = fSparse->Clone(objectName.Data());
  else if ( objectName == "nevents" ) {
    TH1* histo = new TH1D(objectName.Data(),objectName.Data(),1,0.5,1.5);
    obj = histo;
  }
  else if( (objectName == "hNormQA")|| (objectName == "hNormQB") ){
    TH1* histo = new TH1D(objectName.Data(),objectName.Data(),100,0.,1);
    obj = histo;
  }
  else if(objectName == "hScalProdQAQB"){
    TProfile* prof =new TProfile(objectName.Data(),objectName.Data(),300,0.,100.);
    obj = prof;
  }
  else {
    AliError(Form("Unknown object %s\n",objectName.Data()));
  }

  fMergeableCollection->Adopt(identifier, obj);
  AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  return obj;
}

//___________________________________________________________________________
void AliAnalysisTaskV1SingleMu::UserCreateOutputObjects()
{
  //
  /// Create outputs
  //

  //
  Int_t nCentBins = fNCentBins;
  Double_t centMin = 0., centMax = 100.;
  TString centName("Centrality"), centTitle("Centrality"), centUnits("\%");

  Int_t nPtBins = fNPtBins;
  Double_t ptMin = 0., ptMax = 80.;
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");

  Int_t nEtaBins = fNEtaBins;
  Double_t etaMin = -4.5, etaMax = -2.;
  TString etaName("Eta"), etaTitle("#eta"), etaUnits("");

  Int_t nPhiBins = fNPhiBins;
  Double_t phiMin = 0., phiMax = 2.*TMath::Pi();
  TString phiTitle("#phi"), phiUnits("rad");

  Int_t nChargeBins = 2;
  Double_t chargeMin = -2, chargeMax = 2.;
  TString chargeName("Charge"), chargeTitle("charge"), chargeUnits("e");

  Int_t nSPBins = fNSPBins;
  Double_t SPMin = -2, SPMax = 2.;
  TString SPNameA("Scalar product with A"), SPTitleA("SP A"), SPUnits("");
  TString SPNameB("Scalar product with C"), SPTitleB("SP C");

  Int_t nbins[kNvars] = {nCentBins, nPtBins, nEtaBins, nChargeBins, nPhiBins, nSPBins, nSPBins};
  Double_t xmin[kNvars] = {centMin, ptMin, etaMin, chargeMin, phiMin, SPMin, SPMin};
  Double_t xmax[kNvars] = {centMax, ptMax, etaMax, chargeMax, phiMax, SPMax, SPMax};
  TString axisTitle[kNvars] = {centTitle, ptTitle, etaTitle, chargeTitle, phiTitle, SPTitleA, SPTitleB};
  TString axisUnits[kNvars] = {centUnits, ptUnits, etaUnits, chargeUnits, phiUnits, SPUnits, SPUnits};

  fSparse = new THnSparseF("MuSparse","Sparse for muons",kNvars,nbins);

  TString histoTitle = "";
  for ( Int_t idim = 0; idim<kNvars; idim++ ) {
    histoTitle = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
    histoTitle.ReplaceAll("()","");
    fSparse->GetAxis(idim)->SetTitle(histoTitle.Data());

    Double_t array[nbins[idim]+1];
    for ( Int_t ibin=0; ibin<=nbins[idim]; ibin++ ) array[ibin] = xmin[idim] + ibin * (xmax[idim]-xmin[idim])/nbins[idim] ;
    fSparse->SetBinEdges(idim, array);
  }

  // histo = new TH1D("hNormQA", "QA flow vector norm;QAx;QAy",100,0.,1);
  // fMergeableCollection->Adopt(identifier, obj);

  // histo = new TH1D("hNormQB", "QB flow vector norm;QBx;QBy",100,0.,1);

  // TProfile* hScalProdDenom = new TProfile("hScalProdQAQB","hScalProdQAQB",100,0.,100.);

  // for ( Int_t jep=ch+1; jep<fChargeKeys->GetEntries(); jep++ ) {
  //   TString auxCharge = fChargeKeys->At(jep)->GetName();
  //   currTitle = Form("cos(#psi_{%s} - #psi_{%s})",currCharge.Data(),auxCharge.Data());
  //   histo = new TH1D(Form("hResoSub3%s%s",currCharge.Data(),auxCharge.Data()), currTitle.Data(), 100, -1.,1.);
  //   histo->GetXaxis()->SetTitle(currTitle.Data());
  // }

  // currTitle = Form("cos(#psi_{1} - #psi_{2})");
  // histo = new TH1D(Form("hResoSub2%s",currCharge.Data()), Form("%s: %s", currCharge.Data(), currTitle.Data()), 100, -1.,1.);
  // histo->GetXaxis()->SetTitle(currTitle.Data());
  fMergeableCollection = new AliMergeableCollection(GetOutputSlot(1)->GetContainer()->GetName());

  fMuonEventCuts.Print("mask");
  fMuonTrackCuts.Print("mask");

  PostData(1,fMergeableCollection);
  fUtilityMuonAncestor = new AliUtilityMuonAncestor();
}

//________________________________________________________________________
void AliAnalysisTaskV1SingleMu::UserExec ( Option_t * /*option*/ )
{
  //
  /// Fill output objects
  //
  // fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  // if ( ! fAODEvent )
  //   fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  // if ( ! fAODEvent && ! fESDEvent ) {
  //   AliError ("AOD or ESD event not found. Nothing done!");
  //   return;
  // }

  if ( ! fMuonEventCuts.IsSelected(fInputHandler) ) return;

  if ( fZdcCalibOnly && !IsZDCCalibrated(aod->GetRunNumber()) ) return;

  // //
  // // Global event info
  // //
  const TObjArray* selectTrigClasses = fMuonEventCuts.GetSelectedTrigClassesInEvent(fInputHandler);
  // Int_t nSelTrigClasses = selectTrigClasses->GetEntries();

  Double_t containerInput[kNvars];
  Double_t centrality = fMuonEventCuts.GetCentrality(InputEvent());
  if(centrality <5. || centrality > 40.) return;
  containerInput[kHCentrality] = centrality;

  // //
  // // Flow info
  // //
  AliFlowVector vQarray[2];
  Double_t QA[2]={0.,0.};
  Double_t QB[2]={0.,0.};
  // get ZDC Q-vectors
  AliFlowEvent* anEvent = dynamic_cast<AliFlowEvent*>(GetInputData(1));
  if(anEvent) {
      // Get Q vectors for the subevents
      anEvent->GetZDC2Qsub(vQarray);
      QA[0] = vQarray[0].X(); // ZNA
      QA[1] = vQarray[0].Y(); // ZNA
      QB[0] = vQarray[1].X(); // ZNC
      QB[1] = vQarray[1].Y(); // ZNC
  } else {
      AliWarning("WARNING: FlowEvent not found !!! \n");
  }

  if( (QA[0]*QA[0]+QA[1]*QA[1])<10.E-6) {
    AliWarning("WARNING : FlowEvent did not compute the flow vectors for ZDC");
    return;
  }

 //Resolution
  // for ( auto& trigClass : selTrigClasses ) {
    // TString identifier = "Qnorm";
    // static_cast<TH1*>(GetMergeableObject(identifier,"hNormQA"))->Fill(TMath::Sqrt(QA[0]*QA[0]+QA[1]*QA[1]));
    // static_cast<TH1*>(GetMergeableObject(identifier,"hNormQB"))->Fill(TMath::Sqrt(QB[0]*QB[0]+QB[1]*QB[1]));
    // static_cast<TProfile*>(GetMergeableObject(identifier,"hScalProdQAQB"))->Fill(centrality,QA[0]*QB[0] + QA[1]*QB[1]);
  // }

  std::vector<TString> selTrigClasses;
  TIter nextTrig(selectTrigClasses);
  TObject* obj = NULL;
  while ( (obj = nextTrig()) ) selTrigClasses.push_back(obj->GetName());
  for ( auto& trigClass : selTrigClasses ) {
    TString identifier = Form("/%s",trigClass.Data());
    static_cast<TH1*>(GetMergeableObject(identifier,"hNormQA"))->Fill(TMath::Sqrt(QA[0]*QA[0]+QA[1]*QA[1]));
    static_cast<TH1*>(GetMergeableObject(identifier,"hNormQB"))->Fill(TMath::Sqrt(QB[0]*QB[0]+QB[1]*QB[1]));
    static_cast<TProfile*>(GetMergeableObject(identifier,"hScalProdQAQB"))->Fill(centrality,QA[0]*QB[0] + QA[1]*QB[1]);
  } // loop on selected trigger classes

  AliVParticle* track = 0x0;

  Int_t nSteps = MCEvent() ? 2 : 1;
  for ( Int_t istep = 0; istep<nSteps; ++istep ) {
    std::vector<TString> selTrigClasses;
    if ( istep == kStepReconstructed ) {
      TIter nextTrig(selectTrigClasses);
      TObject* obj = NULL;
      while ( (obj = nextTrig()) ) selTrigClasses.push_back(obj->GetName());
    }
    else selTrigClasses.push_back("generated");

    Int_t nTracks = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetNTracks(InputEvent()) : MCEvent()->GetNumberOfTracks();
    //loop on tracks
    // Int_t nSelected = 0;
    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      track = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetTrack(itrack,InputEvent()) : MCEvent()->GetTrack(itrack);

      // In case of MC we usually ask that the particle is a muon
      // However, in W or Z simulations, Pythia stores both the initial muon
      // (before ISR, FSR and kt kick) and the final state one.
      // The first muon is of course there only for information and should be rejected.
      // The Pythia code for initial state particles is 21
      // When running with POWHEG, Pythia puts the hard process input of POWHEG in the stack
      // with state 21, and then re-add it to stack before applying ISR, FSR and kt kick.
      // This muon produces the final state muon, and its status code is 11
      // To avoid all problems, keep only final state muon (status code <10)
      // FIXME: is the convention valid for other generators as well?
      Bool_t isSelected = ( istep == kStepReconstructed ) ? fMuonTrackCuts.IsSelected(track) : ( TMath::Abs(track->PdgCode()) == 13 && AliAnalysisMuonUtility::GetStatusCode(track) < 10 );
      if ( ! isSelected ) {
        continue;
      }

      containerInput[kHvarPt]         = track->Pt();
      containerInput[kHvarEta]        = track->Eta();
      containerInput[kHvarCharge]     = track->Charge()/3.;
      containerInput[kHvarPhi]        = track->Phi();
      // containerInput[kHvarVz]         = ( istep == kStepReconstructed && !fUseMCKineForRecoTracks ) ? ipVz : ipVzMC;

      // Check with Jacopo
      // containerInput[kHvarCosPhi]     = TMath::Cos(fHarmonic*track->Phi());
      // containerInput[kHvarSinPhi]     = TMath::Sin(fHarmonic*track->Phi());
      containerInput[kHvarV1QA]       = TMath::Cos(fHarmonic*track->Phi())*QA[0]+TMath::Sin(fHarmonic*track->Phi())*QA[1];
      containerInput[kHvarV1QB]       = TMath::Cos(fHarmonic*track->Phi())*QB[0]+TMath::Sin(fHarmonic*track->Phi())*QB[1];
      // containerInput[kHvarOdd]        = TMath::Cos(fHarmonic*track->Phi())*(QB[0]-QA[0])+TMath::Sin(fHarmonic*track->Phi())*(QB[1]-QA[1]);
      // containerInput[kHvarEven]       = TMath::Cos(fHarmonic*track->Phi())*(QB[0]+QA[0])+TMath::Sin(fHarmonic*track->Phi())*(QB[1]+QA[1]);



      for ( auto& trigClass : selTrigClasses ) {
        TString identifier = Form("/%s",trigClass.Data());
        static_cast<TH1*>(GetMergeableObject(identifier, "nevents"))->Fill(1.);
        if ( istep == kStepReconstructed && ! fMuonTrackCuts.TrackPtCutMatchTrigClass(track, fMuonEventCuts.GetTrigClassPtCutLevel(trigClass.Data())) ) continue;
        // TString identifier = Form("/%s/%f",trigClass.Data(),centrality);
        static_cast<THnSparse*>(GetMergeableObject(identifier, "MuSparse"))->Fill(containerInput,1.);
      } // loop on selected trigger classes
    } // loop on tracks

  }// loop on container steps

  PostData(1, fMergeableCollection);
}
//________________________________________________________________________
void AliAnalysisTaskV1SingleMu::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts.SetRun(fInputHandler);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskV1SingleMu::IsZDCCalibrated(Int_t run){
  Bool_t calib = kFALSE;
  Int_t dRun15hPos[] = {246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424};
    for (Int_t i=0; i<40; i++) {
      if(run==dRun15hPos[i]) calib = kTRUE;
  }
  return calib;
}
//________________________________________________________________________
Int_t AliAnalysisTaskV1SingleMu::GetParticleType ( AliVParticle* track )
{
  //
  /// Get particle type from matched MC track
  //

  // CAVEAT: here the order matters
  // The muon ancestor class tracks the particle up to the first ancestor
  // So, for example the case
  // W+ -> b -> mu+ -> pi- (from scattering in absorber) -> mu-
  // will be flagged as:
  // - kWbosonMu
  // - kBeautyMu
  // - kSecondaryMu
  // since it is all of them.
  // The user can have one or the other by changing the order of the if conditions
  // in the following
  if ( fUtilityMuonAncestor->IsUnidentified(track,MCEvent()) ) return kUnidentified;
  if ( fUtilityMuonAncestor->IsHadron(track,MCEvent()) ) return kRecoHadron;
  if ( fUtilityMuonAncestor->IsSecondaryMu(track,MCEvent()) ) return kSecondaryMu;
  if ( fUtilityMuonAncestor->IsDecayMu(track,MCEvent()) ) return kDecayMu;
  if ( fUtilityMuonAncestor->IsBeautyMu(track,MCEvent()) ) return kBeautyMu;
  if ( fUtilityMuonAncestor->IsCharmMu(track,MCEvent()) ) return kCharmMu;
  if ( fUtilityMuonAncestor->IsWBosonMu(track,MCEvent()) ) return kWbosonMu;
  if ( fUtilityMuonAncestor->IsZBosonMu(track,MCEvent()) ) return kZbosonMu;
  if ( fUtilityMuonAncestor->IsQuarkoniumMu(track,MCEvent()) ) return kQuarkoniumMu;

  return kDecayMu;
}

//________________________________________________________________________
void AliAnalysisTaskV1SingleMu::Terminate(Option_t *) {
  // //
  // /// Draw some histograms at the end.
  // //

  fMergeableCollection = static_cast<AliMergeableCollection*>(GetOutputData(1));

  if ( ! fMergeableCollection ) return;

  // AliMuonAnalysisOutput muonOut(fOutputList);
  // TList* trigClasses = fMergeableCollection->CreateListOfKeys(0);
  // TIter nextClass(trigClasses);


  // //////////////////////
  // // Event statistics //
  // //////////////////////
  // printf("\nTotal analyzed events:\n");
  // TString evtSel = Form("trigger:%s", trigClassName.Data());
  // muonOut.GetCounterCollection()->PrintSum(evtSel.Data());
  // printf("Physics selected analyzed events:\n");
  // evtSel = Form("trigger:%s/selected:yes", trigClassName.Data());
  // muonOut.GetCounterCollection()->PrintSum(evtSel.Data());

  // TString countPhysSel = "any";
  // if ( physSel.Contains(fPhysSelKeys->At(kPhysSelPass)->GetName()) ) countPhysSel = "yes";
  // else if ( physSel.Contains(fPhysSelKeys->At(kPhysSelReject)->GetName()) ) countPhysSel="no";
  // countPhysSel.Prepend("selected:");
  // printf("Analyzed events vs. centrality:\n");
  // evtSel = Form("trigger:%s/%s", trigClassName.Data(), countPhysSel.Data());
  // muonOut.GetCounterCollection()->Print("centrality",evtSel.Data(),kTRUE);

  // if ( ! outFilename.IsNull() ) {
  //   printf("\nWriting output file %s\n", outFilename.Data());
  //   TFile* file = TFile::Open(outFilename.Data(), "RECREATE");
  //   outList.Write();
  //   file->Close();
  // }
}
