
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
#include "TGraph.h"
#include "AliLog.h"
#include "AliForwardMCClosure.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"
#include "AliAODMCHeader.h"
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
#include "AliMCParticle.h"
#include "AliForwardSecondariesTask.h"
#include "AliGenEventHeaderTunedPbPb.h"
using namespace std;
ClassImp(AliForwardMCClosure)
#if 0
; // For emacs 
#endif

//_____________________________________________________________________
AliForwardMCClosure::AliForwardMCClosure() : AliAnalysisTaskSE(),
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
  AliForwardMCClosure::AliForwardMCClosure(const char* name) : AliAnalysisTaskSE(name),
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
  void AliForwardMCClosure::UserCreateOutputObjects()
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

  fEventList->Add(new TH1D("Centrality","Centrality",fSettings.fCentBins,0,100));
  fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
  fEventList->Add(new TH2F("hOutliers","Maximum #sigma from mean N_{ch} pr. bin", 
     20, 0., 100., 500, 0., 5.)); //((fFlags & kMC) ? 15. : 5. // Sigma <M> histogram 
  fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));
  fEventList->Add(new TH1D("dNdeta","dNdeta",fSettings.fNDiffEtaBins,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge));
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

    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fQcorrfactor", "fQcorrfactor", 1, -4.0, 6.0)); //(eta, n)
    static_cast<TList*>(fGFList->At(2))->Add(new TH1F("fpcorrfactor", "fpcorrfactor", fSettings.fNDiffEtaBins, -4.0, 6.0)); //(eta, n)

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
      Double_t xmin[5] = {0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0, 0};
      Double_t xmax[5] = {10,fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 100, static_cast<Double_t>(fSettings.kSinphi1phi2phi3p+1)};


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
void AliForwardMCClosure::UserExec(Option_t */*option*/)
  {
  //
  //  Analyses the event with use of the helper class AliForwardQCumulantRun2
  //
  //  Parameters:
  //   option: Not used
  //
  AliMCEvent* fAOD = this->MCEvent();
  AliStack* stack = fAOD->Stack();

  if(!fAOD) {
        std::cout << "no aod" << std::endl;
    return;
  }
  /*
  if (!stack) {
    std::cout << "no stack" << std::endl;
    return;
  }
  */
    /*
    if(!fEventCuts.AcceptEvent(fInputEvent)) {
      PostData(1, fOutputList);
      return;
    } 
    */ 

  // Disregard events without reconstructed vertex
  Float_t zvertex = fAOD->GetPrimaryVertex()->GetZ();
  if (!(TMath::Abs(zvertex) > 0)) {
    return;
  }

    AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    double cent = MultSelection->GetMultiplicityPercentile("SPDTracklets");
//double cent = 10;



/*


  AliAODMCHeader* fAODMCHeader = static_cast<AliAODMCHeader*>(fAOD->FindListObject(AliAODMCHeader::StdBranchName()));
  Double_t impactParam[] = { 0.00,  3.72,  5.23,  7.31,  8.88, 10.20, 
          11.38, 12.47, 13.50, 14.51, 16.679};
  Double_t centrality[]  = { 0.,    5.,   10.,   20.,   30.,   40., 
          50.,   60.,   70.,   80.,  100.};

  Int_t nPoints = sizeof(impactParam)/sizeof(Double_t);
  TGraph* fImpactParToCent = new TGraph(nPoints, impactParam, centrality);

    double fCent = 0;
    if (fAODMCHeader){
  AliGenEventHeaderTunedPbPb* header = 
    dynamic_cast<AliGenEventHeaderTunedPbPb*>(fAODMCHeader->GetCocktailHeader(0));
  if (header) fCent = header->GetCentrality();
}

    double cent1 = 0;


if (fAODMCHeader){
  Double_t b = fAODMCHeader->GetImpactParameter();
  cent1 = fImpactParToCent->Eval(b);
}

std::cout << fCent << std::endl;
std::cout << cent1 << std::endl;

*/

  //AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(fAOD->FindListObject("CentralClusters")); // only exists if created by user from ESDs
  TH2D spddNdedp = TH2D("spddNdedp","spddNdedp",400,-1.5,1.5,400,0,2*TMath::Pi()); // Histogram to contain the central tracks
  TH2D centralDiff = TH2D("central","central",400,-1.5,1.5,400,0,2*TMath::Pi()); // Histogram to contain the central tracks
  //TH2D forwarddNdedp = TH2D("forwarddNdedp","forwarddNdedp",400,-4,6,400,0,2*TMath::Pi()); // also known as dNdetadphi
  TH2D forwarddNdedp = TH2D("forwarddNdedp","forwarddNdedp",200,-4,6,20,0,2*TMath::Pi()); // also known as dNdetadphi

  Int_t nTracks   = fAOD->GetNumberOfTracks();
  //Int_t nTracks   = stack->GetNtrack();


std::cout << "noTracks = " << nTracks << std::endl;

      TH1D* dNdeta = static_cast<TH1D*>(fEventList->FindObject("dNdeta"));

  Bool_t truth_run = true;
  // AODs
  /*
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
      //if (!p->IsPhysicalPrimary()) continue;
      if (fabs(p->Eta())>1.6 && p->Charge() != 0) {
        forwarddNdedp.Fill(p->Eta(),p->Phi(),1);
        dNdeta->Fill(p->Eta(),1);
      }
      if (fabs(p->Eta()) <= 1.1 && p->Charge() != 0) {        
        //if (p->Pt()>=0.2 && p->Pt()<5) 
        spddNdedp.Fill(p->Eta(),p->Phi(),1);
        dNdeta->Fill(p->Eta(),1);
      }
    }
    */
  
  //ESDs w. primary
    /*
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
      if (!p->IsPhysicalPrimary()) continue;
      if (p->Charge() == 0) continue;

        if ( (p->Eta() < -1.7 && p->Eta() > -3.4) || (p->Eta() > 1.7 && p->Eta() < 5.0) ){
          forwarddNdedp.Fill(p->Eta(),p->Phi(),1);
          dNdeta->Fill(p->Eta(),1);
        }
        if (p->Pt()>=0.2 && p->Pt()<=5) {
          if (fabs(p->Eta()) < 1.1) {
            spddNdedp.Fill(p->Eta(),p->Phi(),1);
            dNdeta->Fill(p->Eta(),1);
          }
        }
    }
    */
    //ESDs w. secondary
    
    for (Int_t iTr = 0; iTr < nTracks; iTr++) {
      AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
      if (p->Charge() == 0) continue;

      if (AliTrackReference *ref = this->IsHitFMD(p)){
          forwarddNdedp.Fill(p->Eta(),p->Phi(),1);
          dNdeta->Fill(p->Eta(),1);
      }
      if (AliTrackReference *ref = this->IsHitTPC(p)) {
        if (p->Pt()>=0.2 && p->Pt()<=5) {
            spddNdedp.Fill(p->Eta(),p->Phi(),1);
            dNdeta->Fill(p->Eta(),1);
        }
      }
    }

        //if (p->Pt()>0.2 && p->Pt() < 2) centralDiff.Fill(p->Eta(),p->Phi(),1);




    bool useEvent = kTRUE;
  if (nTracks < 10) useEvent = kFALSE;
  //if (!fSettings.ExtraEventCutFMD(forwarddNdedp, cent, true)) useEvent = kFALSE;

  if (useEvent){

    UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

    static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(cent);
    static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(zvertex);

    //AliForwardQCumulantRun2 calculator = AliForwardQCumulantRun2();
    AliForwardGenericFramework calculator = AliForwardGenericFramework();
    calculator.fSettings = fSettings;

    calculator.CumulantsAccumulate(spddNdedp,fOutputList, cent, zvertex,"central",true,true);
    //calculator.CumulantsAccumulate(centralDiff,fOutputList, cent, zvertex,"central",false,true);

    if (calculator.useEvent) calculator.CumulantsAccumulate(forwarddNdedp, fOutputList, cent, zvertex,"forward",false,true);
    if (calculator.useEvent) calculator.saveEvent(fOutputList, cent, zvertex,  randomInt);
    calculator.reset();
  }
  PostData(1, fOutputList); 

  return;
}



AliTrackReference* AliForwardMCClosure::IsHitFMD(AliMCParticle* p) {
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

AliTrackReference* AliForwardMCClosure::IsHitTPC(AliMCParticle* p) {
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



AliMCParticle* AliForwardMCClosure::GetMother(AliMCParticle* p) {
  // Recurses until the mother IsPhysicalPrimary
  // Return NULL if no mother was found
  AliMCEvent* event = this->MCEvent();
  // GetLabel() is the index on the Stack!
  // event->Stack()->IsPhysicalPrimary(p->GetLabel());
  Bool_t isPP = this->IsRedefinedPhysicalPrimary(p);
  // Return this particle if it is stable
  if (isPP) {
    return p;
  }
  else {
    // No stable particle found and no mother left !?
    if (p->GetMother() < 0) {
      return 0x0;
    }
    AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
    return GetMother(ancestor);
  }
}

/*
Bool_t AliForwardMCClosure::AddMotherIfFirstTimeSeen(AliMCParticle* p, std::vector<Int_t> v){

  //Checking if v contains elements (is empty):
  if(v.empty()){
     return false;
  } 
  Int_t x = p->GetLabel();
  if(std::find(v.begin(), v.end(), x) != v.end()) {
      /* v contains x */

//    return true;
  //} else {
      /* v does not contain x */
    //return false;
  //}

//}

Bool_t AliForwardMCClosure::IsRedefinedPhysicalPrimary(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // Is this a pi0 which was produced as a primary particle?
  if (TMath::Abs(p->PdgCode()) == 111 &&
      p->GetLabel() < event->GetNumberOfPrimaries()) {
    std::cout << "found a pi0" << std::endl;
    return true;
  }
  // Is it a Physical Primary by the standard definition?
  Bool_t isPPStandardDef = p->IsPhysicalPrimary();
  AliMCParticle *pi0Candidate = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Check if this is a primary originating from a pi0
  if (isPPStandardDef && pi0Candidate) {
    if (TMath::Abs(pi0Candidate->PdgCode()) == 111){/*pi0*/
    std::cout << "found wrong pi0Candidate" << std::endl;
      return false; // Don't allow stable particles stemming from pi0!
    }
  }
  return isPPStandardDef;
}




//_____________________________________________________________________
void AliForwardMCClosure::Terminate(Option_t */*option*/)
{
  return;
}


//_____________________________________________________________________
//
//
// EOF