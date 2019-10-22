
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
  fAnalysisList   ->SetName("cumulants");

  fOutputList->Add(fAnalysisList);

  // do analysis from v_2 to a maximum of v_5
  constexpr Int_t dimensions = 7;
  Int_t ptnmax =  (fSettings.doPt ? 5 : 0);



  // create a THn for each harmonic
  Int_t rbins[dimensions] = {4,ptnmax + 1,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins, static_cast<Int_t>(fSettings.rW4Four)} ; // n, pt, s, zvtx,eta,cent,kind
  Double_t rmin[dimensions] = {0,0, 0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0, 1};
  Double_t rmax[dimensions] = {4,double(ptnmax+1),double(fSettings.fnoSamples),fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 100, static_cast<Double_t>(fSettings.kW4Four)+1};

  fAnalysisList->Add(new THnD("reference", "reference", dimensions, rbins, rmin, rmax));

  Int_t std_dbins[dimensions] = {4,ptnmax + 1,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, static_cast<Int_t>(fSettings.dW4Four)} ;
  Int_t mixed_dbins[dimensions] = {4,ptnmax + 1,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins, static_cast<Int_t>(fSettings.dWTwoTwoD)} ;
  Double_t dmin[dimensions] = {0,0, 0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0, 1};
  Double_t dmax[dimensions] = {4,double(ptnmax+1),double(fSettings.fnoSamples),fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 100, static_cast<Double_t>(fSettings.dW4Four)+1};
  Double_t dmax_mixed[dimensions] = {4,double(ptnmax+1),double(fSettings.fnoSamples),fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, 100, static_cast<Double_t>(fSettings.dWTwoTwoD)+1};

  fAnalysisList->Add(new THnD("standard_differential","standard_differential", dimensions, std_dbins, dmin, dmax));
  fAnalysisList->Add(new THnD("mixed_differential","mixed_differential", dimensions, mixed_dbins, dmin, dmax_mixed));

  // The THn has dimensions [n, pt, random samples, vertex position, eta, centrality, kind of variable to store]
  // set names
  THnD* reference = static_cast<THnD*>(fAnalysisList->At(0));
  reference->GetAxis(0)->SetName("n");
  reference->GetAxis(1)->SetName("pt");
  reference->GetAxis(2)->SetName("samples");
  reference->GetAxis(3)->SetName("vertex");
  reference->GetAxis(4)->SetName("eta");
  reference->GetAxis(5)->SetName("cent");
  reference->GetAxis(6)->SetName("identifier");
  
  THnD* standard_differential = static_cast<THnD*>(fAnalysisList->At(1));
  standard_differential->GetAxis(0)->SetName("n");
  standard_differential->GetAxis(1)->SetName("pt");
  standard_differential->GetAxis(2)->SetName("samples");
  standard_differential->GetAxis(3)->SetName("vertex");
  standard_differential->GetAxis(4)->SetName("eta");
  standard_differential->GetAxis(5)->SetName("cent");
  standard_differential->GetAxis(6)->SetName("identifier");

  THnD* mixed_differential = static_cast<THnD*>(fAnalysisList->At(2));
  mixed_differential->GetAxis(0)->SetName("n");
  mixed_differential->GetAxis(1)->SetName("pt");
  mixed_differential->GetAxis(2)->SetName("samples");
  mixed_differential->GetAxis(3)->SetName("vertex");
  mixed_differential->GetAxis(4)->SetName("eta");
  mixed_differential->GetAxis(5)->SetName("cent");
  mixed_differential->GetAxis(6)->SetName("identifier");

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
  //fCalculator.fUtil = fUtil;
  Bool_t isgoodrun = kFALSE;
  if (!fSettings.mc){
    isgoodrun = fUtil.IsGoodRun(fInputEvent->GetRunNumber());
  }
  if (!fSettings.mc) fSettings.nua_runnumber = fUtil.GetNUARunNumber(fInputEvent->GetRunNumber());



  // Get the event validation object
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent() || !isgoodrun){
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
  // if (cent > 60.0){
  //   //PostData(1, fOutputList);
  //   return;
  // }


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
    fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE, true,false);
    fCalculator.CumulantsAccumulate(centralDist, cent, zvertex,kFALSE,true,false);
  }
  else{
    if (fSettings.ref_mode & fSettings.kFMDref) fCalculator.CumulantsAccumulate(refDist, cent, zvertex,kTRUE,true,false);
    else fCalculator.CumulantsAccumulate(refDist, cent, zvertex,kFALSE,true,false);
  }
  fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE,false,true);
  fCalculator.CumulantsAccumulate(centralDist, cent, zvertex,kFALSE,false,true);  

  UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);
  fCalculator.saveEvent(fOutputList, cent, zvertex,  randomInt, 0);

  fCalculator.reset();

  centralDist->Reset();
  
  if (!(fSettings.ref_mode & fSettings.kFMDref) || (fSettings.mc && fSettings.esd)) refDist->Reset();
  if ((fSettings.mc && fSettings.use_primaries_fwd) || (fSettings.mc && fSettings.esd)) {
    forwardDist->Reset();
    refDist->Reset();
  }

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
