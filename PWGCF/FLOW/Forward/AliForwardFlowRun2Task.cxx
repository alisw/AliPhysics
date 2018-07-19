
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

using namespace std;
ClassImp(AliForwardFlowRun2Task)
#if 0
; // For emacs 
#endif

//_____________________________________________________________________
AliForwardFlowRun2Task::AliForwardFlowRun2Task() : AliAnalysisTaskSE(),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fStdQCList(0), 
  fGFList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut()
  {
  //
  //  Default constructor
  //
  }

//_____________________________________________________________________
  AliForwardFlowRun2Task::AliForwardFlowRun2Task(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0),           // input event
  fOutputList(0),    // output list
  fStdQCList(0), 
  fGFList(0),
  fEventList(0),
  fRandom(0),
  fSettings(),
  fEventCuts(),
  fMultTOFLowCut(),
  fMultTOFHighCut(),
  fMultCentLowCut()
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

  fStdQCList = new TList(); 
  fGFList = new TList();
  fStdQCList->SetName("StdQC"); 
  fGFList->SetName("GF");

  fEventList = new TList();

  fEventList->Add(new TH1D("Centrality","Centrality",10,0,100));
  fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
  fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin", 
     20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram 
  fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));

  fEventList->SetName("EventInfo");

  fStdQCList->Add(new TList());
  fStdQCList->Add(new TList());
  static_cast<TList*>(fStdQCList->At(0))->SetName("Reference");
  static_cast<TList*>(fStdQCList->At(1))->SetName("Differential");  

  fGFList->Add(new TList());
  fGFList->Add(new TList());   
  fGFList->Add(new TList());   
  static_cast<TList*>(fGFList->At(0))->SetName("Reference");
  static_cast<TList*>(fGFList->At(1))->SetName("Differential"); 
  static_cast<TList*>(fGFList->At(2))->SetName("AutoCorrection"); 

    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fQcorrfactor", "fQcorrfactor", 1, -6.0, 6.0)); //(eta, n)
    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fpcorrfactor", "fpcorrfactor", fSettings.fNDiffEtaBins, -6.0, 6.0)); //(eta, n)

    fOutputList->Add(fStdQCList);
    fOutputList->Add(fGFList);

    fOutputList->Add(fEventList);

    // do analysis to a maximum of v_5
    Int_t fMaxMoment = 5;
    // create a THn for each harmonic
    for (Int_t n = 2; n <= fMaxMoment; n++) {

      Int_t dimensions = 5;

      Int_t dbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
      Int_t rbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins, fSettings.kSinphi1phi2phi3p+1} ;
      Double_t xmin[5] = {0,fSettings.fZVtxAcceptanceLowEdge, -6.0, 0, 0};
      Double_t xmax[5] = {10,fSettings.fZVtxAcceptanceUpEdge, 6, 100, (double)(fSettings.kSinphi1phi2phi3p+1)};


      static_cast<TList*>(fGFList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      static_cast<TList*>(fGFList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));
      static_cast<TList*>(fStdQCList->At(0))->Add(new THnD(Form("cumuRef_v%d", n), Form("cumuRef_v%d", n), dimensions, rbins, xmin, xmax));
      static_cast<TList*>(fStdQCList->At(1))->Add(new THnD(Form("cumuDiff_v%d", n),Form("cumuDiff_%d", n), dimensions, dbins, xmin, xmax));

      // The THn has dimensions [random samples, vertex position, eta, centrality, kind of variable to store]
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(0))   ->FindObject(Form("cumuRef_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(0))->FindObject(Form("cumuRef_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fGFList->At(1))   ->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(4)->SetName("identifier");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(0)->SetName("samples");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(1)->SetName("vertex");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(2)->SetName("eta");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(3)->SetName("cent");
      static_cast<THnD*>(static_cast<TList*>(fStdQCList->At(1))->FindObject(Form("cumuDiff_v%d", n)))->GetAxis(4)->SetName("identifier");
    }
    PostData(1, fOutputList);
  }


//_____________________________________________________________________
void AliForwardFlowRun2Task::UserExec(Option_t */*option*/)
  {
  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //

  //..check if I have AOD
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      return;
    }

  //..AliEventCuts selection
    if(!fEventCuts.AcceptEvent(fInputEvent)) {
      PostData(1, fOutputList);
      return;
    }  
    
    AliAODVertex* aodVtx = fAOD->GetPrimaryVertex();
    
    float v0Centr = 0;

  // Get detector objects
    AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
    AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));

  //AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters")); // only exists if created by user from ESDs
  TH2D spddNdedp = TH2D("spddNdedp","spddNdedp",400,-1.5,1.5,400,0,2*TMath::Pi()); // Histogram to contain the central tracks

  /*
    AliAODTracklets* aodTracklets = fAOD->GetTracklets();

    for (Int_t i = 0; i < aodTracklets->GetNumberOfTracklets(); i++) {
      spddNdedp.Fill(aodTracklets->GetEta(i),aodTracklets->GetPhi(i), 1);
    }
  */




  Int_t  iTracks(fAOD->GetNumberOfTracks());
  for(Int_t i(0); i < iTracks; i++) {

  // loop  over  all  the  tracks
    AliAODTrack* track = static_cast<AliAODTrack *>(fAOD->GetTrack(i));
    if (track->TestFilterBit(kHybrid)){
      if (track->Pt() >= 0.2 && track->Pt() <= 5){
        spddNdedp.Fill(track->Eta(),track->Phi(), 1);
      }
    }
  }


  //const AliAODTracklets* spdmult = fAOD->GetMultiplicity();

  TH2D& forwarddNdedp = aodfmult->GetHistogram(); // also known as dNdetadphi

  if(!fAOD) return;              
  //AliMultSelection* MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  Float_t lPerc = v0Centr; 

 /* if ( MultSelection ) {
    lPerc = MultSelection->GetMultiplicityPercentile("V0M");
    //Quality check
    Int_t lEvSelCode = MultSelection->GetEvSelCode();
    if( lEvSelCode > 0 ) lPerc = lEvSelCode; // also if lEvSelCode > 200, the event is probably useless
    //disregard!
  } 
  else
  {
    //If this happens, re-check if AliMultSelectionTask ran before your task!
    AliInfo("Didn't find MultSelection!"); 
  }*/

  UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

  static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(lPerc);
  static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(aodVtx->GetZ());

  //AliForwardQCumulantRun2 calculator = AliForwardQCumulantRun2();
  AliForwardGenericFramework calculator = AliForwardGenericFramework();
  calculator.fSettings = fSettings;

  calculator.CumulantsAccumulate(spddNdedp,fOutputList, lPerc, aodVtx->GetZ(),"central");

  if (calculator.useEvent) calculator.CumulantsAccumulate(forwarddNdedp, fOutputList, lPerc, aodVtx->GetZ(),"forward");
  if (calculator.useEvent) calculator.saveEvent(fOutputList, lPerc, aodVtx->GetZ(),  randomInt);
  calculator.reset();

  PostData(1, fOutputList); 

  return;
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
