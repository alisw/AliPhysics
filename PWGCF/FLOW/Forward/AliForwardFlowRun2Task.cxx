
//
// Calculate flow in the forward and central regions using the Q cumulants method.
//
// Inputs:
//  - AliAODEvent
//
// Outputs:
//  - AnalysisResults.root or forward_flow.root
//
#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TList.h>
#include <THn.h>

#include "AliLog.h"
#include "AliForwardFlowRun2Task.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"

#include "AliAODForwardMult.h"
#include "AliAODCentralMult.h"
#include "AliAODEvent.h"

#include "AliForwardUtil.h"

#include "AliVVZERO.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"

#include "AliAnalysisFilter.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliForwardSecondariesTask.h"
using namespace std;
ClassImp(AliForwardFlowRun2Task)
#if 0
; // For emacs
#endif

//_____________________________________________________________________
AliForwardFlowRun2Task::AliForwardFlowRun2Task() : AliAnalysisTaskSE(),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fAnalysisList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  useEvent(true)
  {
  //
  //  Default constructor
  //
  }

//_____________________________________________________________________
  AliForwardFlowRun2Task::AliForwardFlowRun2Task(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fAnalysisList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  useEvent(true)
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //
    DefineOutput(1, TList::Class());
  }

//_____________________________________________________________________
  void AliForwardFlowRun2Task::UserCreateOutputObjects()
  {
    //
    //  Create output objects
    //

    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested
    fEventCuts.AddQAplotsToList(fOutputList);

    TRandom r = TRandom();              // random integer to use for creation of samples (used for error bars).
                                          // Needs to be created here, otherwise it will draw the same random number.

    fAnalysisList    = new TList();
    fEventList = new TList();
    fAnalysisList   ->SetName("Analysis");
    fEventList->SetName("EventInfo");

    fEventList->Add(new TH1D("Centrality","Centrality",fSettings.fCentBins,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin",
       20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram
       fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));
       fEventList->Add(new TH1D("EventCuts_FMD","EventCuts_FMD",3,0,3));
       fEventList->Add(new TH1D("dNdeta","dNdeta",fSettings.fNDiffEtaBins,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge));

    fAnalysisList->Add(new TList());
    fAnalysisList->Add(new TList());
    fAnalysisList->Add(new TList());
    static_cast<TList*>(fAnalysisList->At(0))->SetName("Reference");
    static_cast<TList*>(fAnalysisList->At(1))->SetName("Differential");
    static_cast<TList*>(fAnalysisList->At(2))->SetName("AutoCorrection");

    fOutputList->Add(fAnalysisList);
    fOutputList->Add(fEventList);

    // do analysis from v_2 to a maximum of v_5
    Int_t fMaxMoment = 5;
    Int_t dimensions = 5;

    Int_t dbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
    Int_t rbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
    Double_t xmin[5] = {0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0, 0};
    Double_t xmax[5] = {10,fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 100, static_cast<Double_t>(fSettings.kSinphi1phi2phi3p+1)};

    static_cast<TList*>(fAnalysisList->At(2))->Add(new THnD("fQcorrfactor", "fQcorrfactor", dimensions, rbins, xmin, xmax)); //(eta, n)
    static_cast<TList*>(fAnalysisList->At(2))->Add(new THnD("fpcorrfactor","fpcorrfactor", dimensions, dbins, xmin, xmax)); //(eta, n)

    // create a THn for each harmonic
    for (Int_t n = 2; n <= fMaxMoment; n++) {

      static_cast<TList*>(fAnalysisList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      static_cast<TList*>(fAnalysisList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));
      // The THn has dimensions [random samples, vertex position, eta, centrality, kind of variable to store]
      // set names
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(4)->SetName("identifier");
    }
    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardFlowRun2Task::UserExec(Option_t *)
  {
  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //
  AliAODForwardMult* aodfmult = 0;
  TH2D centraldNdedp;
  TH2D forwarddNdedp;
  AliMCEvent* fAODMC;
  if (fSettings.mc){
    fAODMC = this->MCEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      return;
    }
    if (fSettings.use_primaries) forwarddNdedp = TH2D("forwarddNdedp","forwarddNdedp",400,-4,6,400,0,2*TMath::Pi());
    else forwarddNdedp = TH2D("forwarddNdedp","forwarddNdedp",200,-4,6,20,0,2*TMath::Pi());

    if (fSettings.esd){
      AliStack* stack = fAODMC->Stack();
    }
  }
  else{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

    //..check if I have AOD
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      return;
    }
    aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
    forwarddNdedp = aodfmult->GetHistogram();
  }

  //..AliEventCuts selection
  if (!fSettings.esd && !fEventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    return;
  }

  // Get centrality
  AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  double cent = MultSelection->GetMultiplicityPercentile(fSettings.centrality_estimator);

  // Get vertex
  double zvertex = fAOD->GetPrimaryVertex()->GetZ();

if (fSettings.useSPD) centraldNdedp = TH2D("centraldNdedp","centraldNdedp",400,-2.5,2.5,400,0,2*TMath::Pi());
else centraldNdedp = TH2D("centraldNdedp","centraldNdedp",400,-1.5,1.5,400,0,2*TMath::Pi());



  FillHists(centraldNdedp, forwarddNdedp);

  if (!fSettings.ExtraEventCutFMD(forwarddNdedp, cent, fSettings.mc)) {
    useEvent = false;
  }
  static_cast<TH1D*>(fEventList->FindObject("FMDHits"))->Fill(0.5);

  if (useEvent){
    static_cast<TH1D*>(fEventList->FindObject("FMDHits"))->Fill(1.5);

    UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

    static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(cent);
    static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(zvertex);

    //AliForwardQCumulantRun2 calculator = AliForwardQCumulantRun2();
    AliForwardGenericFramework calculator = AliForwardGenericFramework();
    calculator.fSettings = fSettings;

    calculator.CumulantsAccumulate(centraldNdedp,fOutputList, cent, zvertex,"central",true,true);
    calculator.CumulantsAccumulate(forwarddNdedp, fOutputList, cent, zvertex,"forward",false,true);
    calculator.saveEvent(fOutputList, cent, zvertex,  randomInt);
    calculator.reset();

    PostData(1, fOutputList);
  }
  return;
}


void AliForwardFlowRun2Task::FillHists(TH2D centraldNdedp, TH2D forwarddNdedp){
  TH1D* dNdeta = static_cast<TH1D*>(fEventList->FindObject("dNdeta"));

  if (fSettings.mc && fSettings.esd && !fSettings.use_primaries){
    Int_t nTracks = this->MCEvent()->Stack()->GetNtrack();

    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
        AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
      if (p->Charge() == 0) continue;

      if (AliTrackReference *ref = this->IsHitFMD(p)) {
        forwarddNdedp.Fill(p->Eta(),p->Phi(),1);
        dNdeta->Fill(p->Eta(),1);
      }
      if (AliTrackReference *ref = this->IsHitTPC(p)) {
        centraldNdedp.Fill(p->Eta(),p->Phi(),1);
        dNdeta->Fill(p->Eta(),1);
      }
    }
  }

  else if (fSettings.use_primaries && fSettings.mc){
    Int_t nTracksMC   = this->MCEvent()->GetNumberOfTracks();

    for (Int_t iTr = 0; iTr < nTracksMC; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
      if (!p->IsPhysicalPrimary()) continue;
      if (p->Charge() == 0) continue;

      if ( (p->Eta() < -1.7 && p->Eta() > -3.4) || (p->Eta() > 1.7 && p->Eta() < 5.0) ){
        forwarddNdedp.Fill(p->Eta(),p->Phi());
        dNdeta->Fill(p->Eta(),1);
      }
      if (p->Pt()>=0.2 && p->Pt()<=5) {
        if (fabs(p->Eta()) < 1.1) {
          centraldNdedp.Fill(p->Eta(),p->Phi());
          dNdeta->Fill(p->Eta(),1);
        }
      }
    }
  }
  else if (fSettings.fFlowFlags & fSettings.kSPD){
    AliAODTracklets* aodTracklets = fAOD->GetTracklets();

    for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
      centraldNdedp.Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i), 1);
      dNdeta->Fill(aodTracklets->GetEta(i),1);
    }
  }
  else{
    Int_t  iTracks(fAOD->GetNumberOfTracks());
    for(Int_t i(0); i < iTracks; i++) {

    // loop  over  all  the  tracks
      AliAODTrack* track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
      if (track->TestFilterBit(fSettings.tracktype)){
        if (track->Pt() >= 0.2 && track->Pt() <= 5){
          centraldNdedp.Fill(track->Eta(),track->Phi(), 1);
          dNdeta->Fill(track->Eta(),1);
        }
      }
    }
    for (Int_t etaBin = 1; etaBin <= forwarddNdedp.GetNbinsX(); etaBin++) {
      for (Int_t phiBin = 1; phiBin <= forwarddNdedp.GetNbinsY(); phiBin++) {
         dNdeta->Fill(forwarddNdedp.GetXaxis()->GetBinCenter(etaBin),forwarddNdedp.GetBinContent(etaBin, phiBin));
      }
    }
  }
}


AliTrackReference* AliForwardFlowRun2Task::IsHitFMD(AliMCParticle* p) {
  //std::cout << "p->GetNumberOfTrackReferences() = " << p->GetNumberOfTrackReferences() << std::endl;
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    //std::cout << "ref->DetectorId() = " << ref->DetectorId() << std::endl;
    //std::cout << "AliTrackReference::kFMD = " << AliTrackReference::kFMD << std::endl;
    if (!ref || AliTrackReference::kFMD != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}

AliTrackReference* AliForwardFlowRun2Task::IsHitTPC(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) {
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kTPC != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}

//_____________________________________________________________________
void AliForwardFlowRun2Task::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF
