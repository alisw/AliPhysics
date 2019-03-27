
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
#include "AliMCEvent.h"

#include "AliForwardFlowUtil.h"

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
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

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
  centralDist(),
  refDist(),
  forwardDist(),
  fSettings(),
  fUtil(),
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
  centralDist(),
  refDist(),
  forwardDist(),
  fSettings(),
  fUtil(),
  useEvent(true)
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //

  // Rely on validation task for event and track selection
  DefineInput(1, AliForwardTaskValidation::Class());
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________
void AliForwardFlowRun2Task::UserCreateOutputObjects()
  {
    //
    //  Create output objects
    //
    bool saveAutoAdd = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);

    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested

    TRandom r = TRandom();              // random integer to use for creation of samples (used for error bars).
                                        // Needs to be created here, otherwise it will draw the same random number.

    fAnalysisList    = new TList();
    fEventList       = new TList();
    fAnalysisList   ->SetName("Analysis");
    fEventList      ->SetName("EventInfo");

    fEventList->Add(new TH1D("Centrality","Centrality",fSettings.fCentBins,0,100));
    fEventList->Add(new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge));
    fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));
    fEventList->Add(new TH1F("dNdeta","dNdeta",100 /*fSettings.fNDiffEtaBins*/,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge));

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
    Int_t ptnmax =  (fSettings.doPt ? 10 : 0);

    // create a THn for each harmonic
    for (Int_t n = 2; n <= fMaxMoment; n++) {
      for (Int_t ptn = 0; ptn <= ptnmax; ptn++){

        static_cast<TList*>(fAnalysisList->At(0))->Add(new THnD(Form("cumuRef_v%d_pt%d", n,ptn), Form("cumuRef_v%d_pt%d", n,ptn), dimensions, rbins, xmin, xmax));
        static_cast<TList*>(fAnalysisList->At(1))->Add(new THnD(Form("cumuDiff_v%d_pt%d", n,ptn),Form("cumuDiff_%d_pt%d", n,ptn), dimensions, dbins, xmin, xmax));
        // The THn has dimensions [random samples, vertex position, eta, centrality, kind of variable to store]
        // set names
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(0)->SetName("samples");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(1)->SetName("vertex");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(2)->SetName("eta");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(3)->SetName("cent");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))   ->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(4)->SetName("identifier");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(0)->SetName("samples");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(1)->SetName("vertex");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(2)->SetName("eta");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(3)->SetName("cent");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))   ->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(4)->SetName("identifier");
      }
    }

    PostData(1, fOutputList);
    TH1::AddDirectory(saveAutoAdd);
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
  // Get the event validation object
   AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
   if (!ev_val->IsValidEvent()){
      PostData(1, this->fOutputList);
     return;
   }

  if (!fSettings.esd){
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(InputEvent());
    fUtil.fAODevent = aodevent;
    if(!aodevent) throw std::runtime_error("Not AOD as expected");
  }
  if (fSettings.mc) fUtil.fMCevent = this->MCEvent();

  fUtil.fevent = fInputEvent;
  fUtil.fSettings = fSettings;

  // Make centralDist
  Int_t   centralEtaBins = (fSettings.useITS ? 200 : 400);
  Int_t   centralPhiBins = (fSettings.useITS ? 20 : 400);
  Double_t centralEtaMin = (fSettings.useSPD ? -2.5 : fSettings.useITS ? -4 : -1.5);
  Double_t centralEtaMax = (fSettings.useSPD ? 2.5 : fSettings.useITS ? 6 : 1.5);

  // Make refDist
  Int_t   refEtaBins = (((fSettings.ref_mode & fSettings.kITSref) || (fSettings.ref_mode & fSettings.kFMDref)) ? 200 : 400);
  Int_t   refPhiBins = (((fSettings.ref_mode & fSettings.kITSref) || (fSettings.ref_mode & fSettings.kFMDref)) ? 20  : 400);
  Double_t refEtaMin = ((fSettings.ref_mode & fSettings.kSPDref) ? -2.5 
                           : ((fSettings.ref_mode & fSettings.kITSref) || (fSettings.ref_mode & fSettings.kFMDref)) ? -4 
                           : -1.5);
  Double_t refEtaMax = ((fSettings.ref_mode & fSettings.kSPDref) ?  2.5 
                           : ((fSettings.ref_mode & fSettings.kITSref) || (fSettings.ref_mode & fSettings.kFMDref)) ? 6 
                           : 1.5);


  TH2D centralDist_tmp = TH2D("c","",centralEtaBins,centralEtaMin,centralEtaMax,centralPhiBins,0,2*TMath::Pi());
  centralDist_tmp.SetDirectory(0);

  TH2D refDist_tmp = TH2D("r","",refEtaBins,refEtaMin,refEtaMax,refPhiBins,0,2*TMath::Pi());
  refDist_tmp.SetDirectory(0);

  TH2D forwardTrRef  ("ft","",200,-4,6,20,0,TMath::TwoPi());
  TH2D forwardPrim  ("fp","",400,-4,6,400,0,TMath::TwoPi());
  forwardTrRef.SetDirectory(0);
  forwardPrim.SetDirectory(0);
  forwardDist = (fSettings.use_primaries_fwd ? &forwardPrim : &forwardTrRef);

  centralDist = &centralDist_tmp;
  centralDist->SetDirectory(0);
  refDist = &refDist_tmp;
  refDist->SetDirectory(0);

  TH1F* dNdeta = static_cast<TH1F*>(fEventList->FindObject("dNdeta"));

  dNdeta->SetDirectory(0);
  fUtil.dodNdeta = kTRUE;
  fUtil.dNdeta = dNdeta;
  fUtil.FillData(refDist,centralDist,forwardDist);


  Double_t zvertex = fUtil.GetZ();
  Double_t cent = fUtil.GetCentrality(fSettings.centrality_estimator);

  if (fSettings.makeFakeHoles) fUtil.MakeFakeHoles(*forwardDist);

  static_cast<TH1D*>(fEventList->FindObject("Centrality"))->Fill(cent);
  static_cast<TH1D*>(fEventList->FindObject("Vertex"))->Fill(zvertex);
  AliForwardGenericFramework calculator = AliForwardGenericFramework();
  calculator.fSettings = fSettings;

  calculator.CumulantsAccumulate(*refDist, fOutputList, cent, zvertex,"central",true,false);
  calculator.CumulantsAccumulate(*centralDist, fOutputList, cent, zvertex,"central",false,true);  
  calculator.CumulantsAccumulate(*forwardDist, fOutputList, cent, zvertex,"forward",false,true);

  Int_t ptnmax =  (fSettings.doPt ? 10 : 0);

  TH1F pthist = TH1F("pthist", "", ptnmax, fSettings.minpt, fSettings.maxpt);

  for (Int_t ptn = 0; ptn <=ptnmax; ptn ++ ){
    
    fUtil.fSettings.minpt = pthist.GetXaxis()->GetBinLowEdge(ptn+1);
    fUtil.fSettings.maxpt = pthist.GetXaxis()->GetBinUpEdge(ptn+1);

    centralDist->Reset();

    // Fill centralDist
    if (fSettings.useSPD) fUtil.FillFromTracklets(centralDist);
    else  fUtil.FillFromTracks(centralDist, fSettings.tracktype); //(fSettings.useTPC)

    UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

    calculator.CumulantsAccumulate(*centralDist, fOutputList, cent, zvertex,"central",false,true);  
    calculator.saveEvent(fOutputList, cent, zvertex,  randomInt, ptn);    
    calculator.fpvector->Reset();
  }

  calculator.reset();

  PostData(1, fOutputList);
  return;
}


//_____________________________________________________________________
void AliForwardFlowRun2Task::Terminate(Option_t */*option*/)
{
  std::cout << "Terminating" << '\n';
  return;
}


//_______________________________________________________________
