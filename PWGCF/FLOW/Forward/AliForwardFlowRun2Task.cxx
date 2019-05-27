
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
#include <TList.h>
#include <THn.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliMCEvent.h"

#include "AliForwardFlowRun2Task.h"
#include "AliForwardQCumulantRun2.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardFlowUtil.h"
#include "AliAODForwardMult.h"

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
  fStorage(nullptr),
  fSettings(),
  fUtil(),
  fCalculator()
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
  fStorage(nullptr),
  fSettings(),
  fUtil(),
  fCalculator()
  {
  //
  //  Constructor
  //
  //  Parameters:
  //   name: Name of task
  //

  // Rely on validation task for event and track selection
  DefineInput(1, AliForwardTaskValidation::Class());

  DefineOutput(1, AliForwardFlowResultStorage::Class());
}

//_____________________________________________________________________
void AliForwardFlowRun2Task::UserCreateOutputObjects()
  {
    //
    //  Create output objects
    //
    //bool saveAutoAdd = TH1::AddDirectoryStatus();

    fOutputList = new TList();          // the final output list
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested

    TRandom r = TRandom();              // random integer to use for creation of samples (used for error bars).
                                        // Needs to be created here, otherwise it will draw the same random number.

    fAnalysisList    = new TList();
    //fEventList       = new TList();
    fAnalysisList   ->SetName("Analysis");
    //fEventList      ->SetName("EventInfo");

    //fCent = new TH1D("Centrality","Centrality",fSettings.fCentBins,0,60);
    //fEventList->Add(fCent);
    //fVertex = new TH1D("Vertex","Vertex",fSettings.fNZvtxBins,fSettings.fZVtxAcceptanceLowEdge,fSettings.fZVtxAcceptanceUpEdge);
    //fEventList->Add(fVertex);
    //fCent->SetDirectory(0);
    //fVertex->SetDirectory(0);

    //fEventList->Add(new TH1D("FMDHits","FMDHits",100,0,10));

    // fdNdeta = new TH2D("dNdeta","dNdeta",200 /*fSettings.fNDiffEtaBins*/,fSettings.fEtaLowEdge,fSettings.fEtaUpEdge,fSettings.fCentBins,0,60);
    // fdNdeta->SetDirectory(0);
    // fOutputList->Add(fdNdeta);

    fAnalysisList->Add(new TList());
    fAnalysisList->Add(new TList());
    //fAnalysisList->Add(new TList());
    static_cast<TList*>(fAnalysisList->At(0))->SetName("Reference");
    static_cast<TList*>(fAnalysisList->At(1))->SetName("Differential");
    //static_cast<TList*>(fAnalysisList->At(2))->SetName("AutoCorrection");

    fOutputList->Add(fAnalysisList);
    //fOutputList->Add(fEventList);

    // do analysis from v_2 to a maximum of v_5
    Int_t fMaxMoment = 4;
    Int_t dimensions = 5;

    Int_t dbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, static_cast<Int_t>(fSettings.kW4Four)} ;
    Int_t rbins[5] = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins, static_cast<Int_t>(fSettings.kW4Four)} ;
    Double_t xmin[5] = {0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0, 1};
    Double_t xmax[5] = {10,fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 60, static_cast<Double_t>(fSettings.kW4Four)+1};

    //static_cast<TList*>(fAnalysisList->At(2))->Add(new THnF("fQcorrfactor", "fQcorrfactor", dimensions, rbins, xmin, xmax)); //(eta, n)
    //static_cast<TList*>(fAnalysisList->At(2))->Add(new THnF("fpcorrfactor","fpcorrfactor", dimensions, dbins, xmin, xmax)); //(eta, n)
    Int_t ptnmax =  (fSettings.doPt ? 10 : 0);

    // create a THn for each harmonic
    for (Int_t n = 2; n <= fMaxMoment; n++) {
      for (Int_t ptn = 0; ptn <= ptnmax; ptn++){

        static_cast<TList*>(fAnalysisList->At(0))->Add(new THnD(Form("cumuRef_v%d_pt%d", n,ptn), Form("cumuRef_v%d_pt%d", n,ptn), dimensions, rbins, xmin, xmax));
        static_cast<TList*>(fAnalysisList->At(1))->Add(new THnD(Form("cumuDiff_v%d_pt%d", n,ptn),Form("cumuDiff_v%d_pt%d", n,ptn), dimensions, dbins, xmin, xmax));
        // The THn has dimensions [random samples, vertex position, eta, centrality, kind of variable to store]
        // set names
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(0)->SetName("samples");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(1)->SetName("vertex");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(2)->SetName("eta");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(3)->SetName("cent");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(0))->FindObject(Form("cumuRef_v%d_pt%d", n,ptn)))->GetAxis(4)->SetName("identifier");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(0)->SetName("samples");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(1)->SetName("vertex");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(2)->SetName("eta");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(3)->SetName("cent");
        static_cast<THnD*>(static_cast<TList*>(fAnalysisList->At(1))->FindObject(Form("cumuDiff_v%d_pt%d", n,ptn)))->GetAxis(4)->SetName("identifier");
      }
    }

  // Make centralDist
  Int_t   centralEtaBins = (fSettings.useITS ? 200 : 400);
  Int_t   centralPhiBins = (fSettings.useITS ? 20 : 400);
  Double_t centralEtaMin = (fSettings.useSPD ? -2.5 : fSettings.useITS ? -4 : -1.5);
  Double_t centralEtaMax = (fSettings.useSPD ? 2.5 : fSettings.useITS ? 6 : 1.5);

  // Make refDist
  Int_t   refEtaBins = (((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 200 : 400);
  Int_t   refPhiBins = (((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 20  : 400);
  Double_t refEtaMin = ((fSettings.ref_mode & fSettings.kSPDref) ? -2.5 
                             : ((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? -4 
                             : -1.5);
  Double_t refEtaMax = ((fSettings.ref_mode & fSettings.kSPDref) ?  2.5 
                             : ((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 6 
                             : 1.5);

  centralDist = new TH2D("c","",centralEtaBins,centralEtaMin,centralEtaMax,centralPhiBins,0,2*TMath::Pi());
  centralDist ->SetDirectory(0);
  refDist     = new TH2D("r","",refEtaBins,refEtaMin,refEtaMax,refPhiBins,0,2*TMath::Pi());
  refDist     ->SetDirectory(0);
  forwardDist = new TH2D("ft","",200,-4,6,20,0,TMath::TwoPi());
  forwardDist ->SetDirectory(0);

  fStorage = new AliForwardFlowResultStorage(fSettings.fileName, fOutputList);

  PostData(1, fStorage);
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
  //forwardDist = 0;


  fCalculator.fSettings = fSettings;
  fUtil.fSettings = fSettings;

  // Get the event validation object
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent()){
  //  PostData(1, this->fOutputList);
    PostData(1, fStorage);

    return;
  }

  if (!fSettings.esd){
    AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(InputEvent());
    fUtil.fAODevent = aodevent;
    if(!aodevent) throw std::runtime_error("Not AOD as expected");
  }
  if (fSettings.mc) fUtil.fMCevent = this->MCEvent();

  fUtil.fevent = fInputEvent;

  Double_t cent = fUtil.GetCentrality(fSettings.centrality_estimator);
  if (cent > 60.0){
    //PostData(1, fOutputList);
    return;
  }


  fUtil.FillData(refDist,centralDist,forwardDist);
  
  // dNdeta
  // for (Int_t etaBin = 1; etaBin <= centralDist->GetNbinsX(); etaBin++) {
  //   Double_t eta = centralDist->GetXaxis()->GetBinCenter(etaBin);
  //   for (Int_t phiBin = 1; phiBin <= centralDist->GetNbinsX(); phiBin++) {
  //     fdNdeta->Fill(eta,cent,centralDist->GetBinContent(etaBin,phiBin));
  //   }
  // }
  // for (Int_t etaBin = 1; etaBin <= forwardDist->GetNbinsX(); etaBin++) {
  //   Double_t eta = forwardDist->GetXaxis()->GetBinCenter(etaBin);
  //   for (Int_t phiBin = 1; phiBin <= forwardDist->GetNbinsX(); phiBin++) {
  //     fdNdeta->Fill(eta,cent,forwardDist->GetBinContent(etaBin,phiBin));
  //   }
  // }

  Double_t zvertex = fUtil.GetZ();

  //if (fSettings.makeFakeHoles) fUtil.MakeFakeHoles(*forwardDist);

  //fCent->Fill(cent);
  //fVertex->Fill(zvertex);
  
  if (fSettings.a5){
    fCalculator.CumulantsAccumulate(forwardDist, fOutputList, cent, zvertex,kTRUE, true,false);
    fCalculator.CumulantsAccumulate(centralDist, fOutputList, cent, zvertex,kFALSE,true,false);
  }
  else{
    if (fSettings.ref_mode & fSettings.kFMDref) fCalculator.CumulantsAccumulate(refDist, fOutputList, cent, zvertex,kTRUE,true,false);
    else fCalculator.CumulantsAccumulate(refDist, fOutputList, cent, zvertex,kFALSE,true,false);
  }
  fCalculator.CumulantsAccumulate(forwardDist, fOutputList, cent, zvertex,kTRUE,false,true);
  fCalculator.CumulantsAccumulate(centralDist, fOutputList, cent, zvertex,kFALSE,false,true);  

  UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);
  fCalculator.saveEvent(fOutputList, cent, zvertex,  randomInt, 0);   

  fCalculator.reset();

  centralDist->Reset();
  
  if (!(fSettings.ref_mode & fSettings.kFMDref)) refDist->Reset();
  if (fSettings.mc && fSettings.use_primaries_fwd) forwardDist->Reset();

  //PostData(1, fOutputList);    
  PostData(1, fStorage);


  return;
}


//_____________________________________________________________________
void AliForwardFlowRun2Task::Terminate(Option_t */*option*/)
{
  std::cout << "Terminating" << '\n';
  return;
}


//_____________________________________________________________________
