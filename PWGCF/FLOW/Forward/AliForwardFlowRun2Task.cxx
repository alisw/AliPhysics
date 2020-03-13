
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

#include "AliAODEvent.h"

#include "AliForwardFlowRun2Task.h"
#include "AliForwardGenericFramework.h"
#include "AliForwardFlowUtil.h"

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
  fCalculator(2)
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
  fCalculator(2)
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
  std::cout << "void AliForwardFlowRun2Task::UserCreateOutputObjects()"<< std::endl;

  fWeights = AliForwardWeights();
  fWeights.fSettings = this->fSettings;

  if (fSettings.nua_file != "") fWeights.connectNUA(); else fWeights.fSettings.doNUA = kFALSE;
  if (fSettings.nue_file != "") fWeights.connectNUE(); else fWeights.fSettings.doNUE = kFALSE;
  if (fSettings.sec_file != "") fWeights.connectSec(); else fWeights.fSettings.sec_corr = kFALSE;
  if (fSettings.sec_cent_file != "") fWeights.connectSecCent();
  this->fSettings = fWeights.fSettings;


  if (fSettings.etagap){
    fCalculator = AliForwardGenericFramework(2);
  }
  else{
    fCalculator = AliForwardGenericFramework(3);
  }

  fOutputList = new TList();          // the final output list
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them if requested

  TRandom r = TRandom();              // random integer to use for creation of samples (used for error bars).
                                      // Needs to be created here, otherwise it will draw the same random number.

  fAnalysisList    = new TList();
  fAnalysisList   ->SetName("cumulants");

  fOutputList->Add(fAnalysisList);


  fReferenceList    = new TList();
  fReferenceList   ->SetName("reference");
  fAnalysisList->Add(fReferenceList);

  fStandardList    = new TList();
  fStandardList   ->SetName("standard");
  fAnalysisList->Add(fStandardList);

  fMixedList    = new TList();
  fMixedList   ->SetName("mixed");
  fAnalysisList->Add(fMixedList);


  // do analysis from v_2 to a maximum of v_4
  constexpr Int_t dimensions = 5;
  Int_t ptnmax =  (fSettings.doPt ? 5 : 0);

  // create a THn for each harmonic
  Int_t rbins[dimensions]        = {3,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins} ; // n, pt, s, zvtx,eta,cent
  Int_t dbins[dimensions]        = {3,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins} ;
  Int_t negonly_bins[dimensions] = {3,fSettings.fnoSamples, fSettings.fNZvtxBins, 14, fSettings.fCentBins};
  Int_t non_rbins[dimensions-1]  = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, fSettings.fCentBins} ; // n, pt, s, zvtx,eta,cent
  Int_t non_dbins[dimensions-1]  = {fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, fSettings.fCentBins} ; // n, pt, s, zvtx,eta,cent
  
  Double_t non_min[dimensions-1] = {0, fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0};
  Double_t dmin[dimensions]      = {0, 0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0};

  Double_t dmax[dimensions]      = {3,double(fSettings.fnoSamples), fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, double(fSettings.fCentUpEdge)};
  Double_t negonly_max[dimensions]={3,double(fSettings.fnoSamples), fSettings.fZVtxAcceptanceUpEdge, 0, double(fSettings.fCentUpEdge)};
  Double_t non_max[dimensions-1]  = {double(fSettings.fnoSamples),fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, double(fSettings.fCentUpEdge)};


  if ((fSettings.normal_analysis || fSettings.second_analysis) || fSettings.SC_analysis){
    fCalculator.cumu_rW2     = new THnD("cumu_rW2",     "cumu_rW2",     dimensions-1,non_rbins,non_min,non_max);
    fCalculator.cumu_rW2Two  = new THnD("cumu_rW2Two" , "cumu_rW2Two" , dimensions,rbins,dmin,dmax); 
    TList* list_rW2     = new TList(); list_rW2    ->SetName("rW2"    ); list_rW2    ->Add(fCalculator.cumu_rW2);     fReferenceList->Add(list_rW2    );
    TList* list_rW2Two  = new TList(); list_rW2Two ->SetName("rW2Two" ); list_rW2Two ->Add(fCalculator.cumu_rW2Two);  fReferenceList->Add(list_rW2Two );  
  }

  if ((fSettings.normal_analysis || fSettings.second_analysis) || fSettings.decorr_analysis){
    fCalculator.cumu_dW2B    = new THnD("cumu_dW2B"   , "cumu_dW2B"   , dimensions-1, non_dbins, non_min, non_max);
    fCalculator.cumu_dW2TwoB = new THnD("cumu_dW2TwoB", "cumu_dW2TwoB", dimensions, dbins, dmin, dmax);
    TList* list_dW2B    = new TList(); list_dW2B   ->SetName("dW2B"   ); list_dW2B   ->Add(fCalculator.cumu_dW2B   ); fStandardList->Add(list_dW2B   );
    TList* list_dW2TwoB = new TList(); list_dW2TwoB->SetName("dW2TwoB"); list_dW2TwoB->Add(fCalculator.cumu_dW2TwoB); fStandardList->Add(list_dW2TwoB);
  }

  if (fSettings.SC_analysis || fSettings.normal_analysis){
    fCalculator.cumu_dW4     = new THnD("cumu_dW4"    , "cumu_dW4"    , dimensions-1, non_dbins, non_min, non_max);
    TList* list_dW4     = new TList(); list_dW4    ->SetName("dW4"    ); list_dW4    ->Add(fCalculator.cumu_dW4    ); fStandardList->Add(list_dW4    );
  }

  if (fSettings.normal_analysis) {
    fCalculator.cumu_rW4     = new THnD("cumu_rW4"    , "cumu_rW4"    , dimensions-1,non_rbins,non_min,non_max);
    fCalculator.cumu_rW4Four = new THnD("cumu_rW4Four", "cumu_rW4Four", dimensions,rbins,dmin,dmax);
    TList* list_rW4     = new TList(); list_rW4    ->SetName("rW4"    ); list_rW4    ->Add(fCalculator.cumu_rW4);     fReferenceList->Add(list_rW4    );
    TList* list_rW4Four = new TList(); list_rW4Four->SetName("rW4Four"); list_rW4Four->Add(fCalculator.cumu_rW4Four); fReferenceList->Add(list_rW4Four);


    fCalculator.cumu_dW4Four = new THnD("cumu_dW4Four", "cumu_dW4Four", dimensions, dbins, dmin, dmax);    
    TList* list_dW4Four = new TList(); list_dW4Four->SetName("dW4Four"); list_dW4Four->Add(fCalculator.cumu_dW4Four); fStandardList->Add(list_dW4Four);
  }

  if (fSettings.decorr_analysis){
    fCalculator.cumu_dW2A    = new THnD("cumu_dW2A"   , "cumu_dW2A"   , dimensions-1, non_dbins, non_min, non_max);
    fCalculator.cumu_dW2TwoA = new THnD("cumu_dW2TwoA", "cumu_dW2TwoA", dimensions, dbins, dmin, dmax);
    TList* list_dW2A    = new TList(); list_dW2A   ->SetName("dW2A"   ); list_dW2A   ->Add(fCalculator.cumu_dW2A   ); fStandardList->Add(list_dW2A   );
    TList* list_dW2TwoA = new TList(); list_dW2TwoA->SetName("dW2TwoA"); list_dW2TwoA->Add(fCalculator.cumu_dW2TwoA); fStandardList->Add(list_dW2TwoA);


    fCalculator.cumu_dW2TwoTwoD = new THnD("cumu_dW2TwoTwoD", "cumu_dW2TwoTwoD", dimensions, negonly_bins, dmin, negonly_max) ;
    TList* list_dW2TwoTwoD = new TList(); list_dW2TwoTwoD->SetName("dW2TwoTwoD"); list_dW2TwoTwoD->Add(fCalculator.cumu_dW2TwoTwoD); fMixedList->Add(list_dW2TwoTwoD);
    fCalculator.cumu_dW2TwoTwoN = new THnD("cumu_dW2TwoTwoN", "cumu_dW2TwoTwoN", dimensions, negonly_bins, dmin, negonly_max) ;
    TList* list_dW2TwoTwoN = new TList(); list_dW2TwoTwoN->SetName("dW2TwoTwoN"); list_dW2TwoTwoN->Add(fCalculator.cumu_dW2TwoTwoN); fMixedList->Add(list_dW2TwoTwoN);
  }

  if (fSettings.SC_analysis){
    fCalculator.cumu_dW4FourTwo  = new THnD("cumu_dW4FourTwo" , "cumu_dW4FourTwo" , dimensions-1, non_dbins, non_min, non_max) ;
    fCalculator.cumu_dW4ThreeTwo = new THnD("cumu_dW4ThreeTwo", "cumu_dW4ThreeTwo", dimensions-1, non_dbins, non_min, non_max) ;
    TList* list_dW4FourTwo  = new TList(); list_dW4FourTwo ->SetName("dW4FourTwo" ); list_dW4FourTwo ->Add(fCalculator.cumu_dW4FourTwo ); fMixedList->Add(list_dW4FourTwo );
    TList* list_dW4ThreeTwo = new TList(); list_dW4ThreeTwo->SetName("dW4ThreeTwo"); list_dW4ThreeTwo->Add(fCalculator.cumu_dW4ThreeTwo); fMixedList->Add(list_dW4ThreeTwo);

    fCalculator.cumu_wSC = new THnD("cumu_wSC", "cumu_wSC", dimensions-1, non_dbins, non_min, non_max) ;
    TList* list_wSC = new TList(); list_wSC->SetName("wSC"); list_wSC->Add(fCalculator.cumu_wSC); fMixedList->Add(list_wSC);


    fCalculator.cumu_dW22TwoTwoN   = new THnD("cumu_dW22TwoTwoN"  , "cumu_d22WTwoTwoN"  , dimensions-1, non_dbins, non_min, non_max) ;
    fCalculator.cumu_dW22TwoTwoD   = new THnD("cumu_dW22TwoTwoD"  , "cumu_d22WTwoTwoD"  , dimensions-1, non_dbins, non_min, non_max) ;
    TList* list_dW22TwoTwoN   = new TList(); list_dW22TwoTwoN  ->SetName("dW22TwoTwoN"  ); list_dW22TwoTwoN  ->Add(fCalculator.cumu_dW22TwoTwoN  ); fMixedList->Add(list_dW22TwoTwoN  );
    TList* list_dW22TwoTwoD   = new TList(); list_dW22TwoTwoD  ->SetName("dW22TwoTwoD"  ); list_dW22TwoTwoD  ->Add(fCalculator.cumu_dW22TwoTwoD  ); fMixedList->Add(list_dW22TwoTwoD  );
  }

  // Make centralDist
  Int_t   centralEtaBins = (fSettings.useITS ? 200 : 300);
  Int_t   centralPhiBins = (fSettings.useITS ? 20 : 300);
  Double_t centralEtaMin = (fSettings.useSPD ? -2.5 : fSettings.useITS ? -4 : -1.2);
  Double_t centralEtaMax = (fSettings.useSPD ? 2.5 : fSettings.useITS ? 6 : 1.2);

  // Make refDist
  Int_t   refEtaBins = (((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 200 : 300);
  Int_t   refPhiBins = (((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 20  : 300);
  Double_t refEtaMin = ((fSettings.ref_mode & fSettings.kSPDref) ? -2.5 
                             : ((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? -4 
                             : -1.2);
  Double_t refEtaMax = ((fSettings.ref_mode & fSettings.kSPDref) ?  2.5 
                             : ((fSettings.ref_mode & fSettings.kITSref) | (fSettings.ref_mode & fSettings.kFMDref)) ? 6 
                             : 1.2);

  centralDist = new TH2D("c","",centralEtaBins,centralEtaMin,centralEtaMax,centralPhiBins,0,2*TMath::Pi());
  centralDist ->SetDirectory(0);
  refDist     = new TH2D("r","",refEtaBins,refEtaMin,refEtaMax,refPhiBins,0,2*TMath::Pi());
  refDist     ->SetDirectory(0);
  forwardDist = new TH2D("ft","",200,-4,6,20,0,TMath::TwoPi());
  forwardDist ->SetDirectory(0);

  fStorage = new AliForwardFlowResultStorage(fSettings.fileName, fOutputList);
  
  fCalculator.fSettings = fSettings;
  fUtil.fSettings = fSettings;

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
  
  if (fSettings.doNUA) fSettings.nua_runnumber = fUtil.GetNUARunNumber(fInputEvent->GetRunNumber());

  // Get the event validation object
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent()){
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
  if (cent > Double_t(fSettings.fCentUpEdge)) return;

  fUtil.FillData(refDist,centralDist,forwardDist);
  if (fSettings.makeFakeHoles) fUtil.MakeFakeHoles(*forwardDist);

  Double_t zvertex = fUtil.GetZ();


  if (!fSettings.etagap){
    fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE,true,false);
    fCalculator.CumulantsAccumulate(refDist, cent, zvertex,kFALSE,true,false);
  }
  else{
    if (fSettings.ref_mode & fSettings.kFMDref) {
      fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE,true,false);
    }
    else {
      fCalculator.CumulantsAccumulate(refDist, cent, zvertex,kFALSE,true,false);
    }
  }
  
  fCalculator.CumulantsAccumulate(centralDist, cent, zvertex,kFALSE,false,true);  
  fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE, false,true);  

  UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

  fCalculator.saveEvent(cent, zvertex,  randomInt, 0);

  fCalculator.reset();
  centralDist->Reset();
  
  if (!(fSettings.ref_mode & fSettings.kFMDref)) refDist->Reset();
  if ((fSettings.mc && fSettings.use_primaries_fwd) || (fSettings.mc && fSettings.esd)) {
    forwardDist->Reset();
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
