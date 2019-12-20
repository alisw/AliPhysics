
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
#include "AliForwardQCumulantRun2.h"
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
  fCalculator(),
  fCentCounter()
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
  fCalculator(),
  fCentCounter()
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
  constexpr Int_t dimensions = 6;
  Int_t ptnmax =  (fSettings.doPt ? 5 : 0);

  // create a THn for each harmonic
  Int_t rbins[dimensions] = {3,ptnmax + 1,fSettings.fnoSamples, fSettings.fNZvtxBins, fSettings.fNRefEtaBins, 
                             fSettings.fCentBins} ; // n, pt, s, zvtx,eta,cent
  Double_t dmin[dimensions]     = {0,0, 0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0};

  Double_t dmax[dimensions]     = {3,double(ptnmax+1),double(fSettings.fnoSamples),
                                   fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, double(fSettings.fCentUpEdge)};
  Int_t dbins[dimensions]   = {3,ptnmax + 1,fSettings.fnoSamples, 
                                   fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, 
                                   fSettings.fCentBins} ;


  Double_t sc_dmin[dimensions]     = {0,0, 0,fSettings.fZVtxAcceptanceLowEdge, fSettings.fEtaLowEdge, 0};

  Double_t sc_dmax[dimensions]     = {1,double(ptnmax+1),double(fSettings.fnoSamples),
                                   fSettings.fZVtxAcceptanceUpEdge, fSettings.fEtaUpEdge, double(fSettings.fCentUpEdge)};
  Int_t sc_dbins[dimensions]   = {1,ptnmax + 1,fSettings.fnoSamples, 
                                   fSettings.fNZvtxBins, fSettings.fNDiffEtaBins, 
                                   fSettings.fCentBins} ;


                                   

  fCalculator.cumu_rW2     = new THnD("cumu_rW2",     "cumu_rW2",     dimensions,rbins,dmin,dmax);
  fCalculator.cumu_rW2Two  = new THnD("cumu_rW2Two" , "cumu_rW2Two" , dimensions,rbins,dmin,dmax); 
  fCalculator.cumu_rW4     = new THnD("cumu_rW4"    , "cumu_rW4"    , dimensions,rbins,dmin,dmax);;
  fCalculator.cumu_rW4Four = new THnD("cumu_rW4Four", "cumu_rW4Four", dimensions,rbins,dmin,dmax);;  



  TList* list_rW2     = new TList(); list_rW2    ->SetName("rW2"    ); list_rW2    ->Add(fCalculator.cumu_rW2);     fReferenceList->Add(list_rW2    );
  TList* list_rW2Two  = new TList(); list_rW2Two ->SetName("rW2Two" ); list_rW2Two ->Add(fCalculator.cumu_rW2Two);  fReferenceList->Add(list_rW2Two );
  TList* list_rW4     = new TList(); list_rW4    ->SetName("rW4"    ); list_rW4    ->Add(fCalculator.cumu_rW4);     fReferenceList->Add(list_rW4    );
  TList* list_rW4Four = new TList(); list_rW4Four->SetName("rW4Four"); list_rW4Four->Add(fCalculator.cumu_rW4Four); fReferenceList->Add(list_rW4Four);



  fCalculator.cumu_dW2A    = new THnD("cumu_dW2A"   , "cumu_dW2A"   , dimensions, dbins, dmin, dmax); // multiplicity for all particles in subevent A (note subevent A can also be the entire event)
  fCalculator.cumu_dW2TwoA = new THnD("cumu_dW2TwoA", "cumu_dW2TwoA", dimensions, dbins, dmin, dmax); // <w2*two>
  fCalculator.cumu_dW2B    = new THnD("cumu_dW2B"   , "cumu_dW2B"   , dimensions, dbins, dmin, dmax); // multiplicity for all particles in subevent B (note subevent B can NOT be the entire event)
  fCalculator.cumu_dW2TwoB = new THnD("cumu_dW2TwoB", "cumu_dW2TwoB", dimensions, dbins, dmin, dmax); // <w2*two>  Int_t kW4          = 3; // <w4>
  fCalculator.cumu_dW4     = new THnD("cumu_dW4"    , "cumu_dW4"    , dimensions, dbins, dmin, dmax);
  fCalculator.cumu_dW4Four = new THnD("cumu_dW4Four", "cumu_dW4Four", dimensions, dbins, dmin, dmax);

  TList* list_dW2A    = new TList(); list_dW2A   ->SetName("dW2A"   ); list_dW2A   ->Add(fCalculator.cumu_dW2A   ); fStandardList->Add(list_dW2A   );
  TList* list_dW2TwoA = new TList(); list_dW2TwoA->SetName("dW2TwoA"); list_dW2TwoA->Add(fCalculator.cumu_dW2TwoA); fStandardList->Add(list_dW2TwoA);
  TList* list_dW2B    = new TList(); list_dW2B   ->SetName("dW2B"   ); list_dW2B   ->Add(fCalculator.cumu_dW2B   ); fStandardList->Add(list_dW2B   );
  TList* list_dW2TwoB = new TList(); list_dW2TwoB->SetName("dW2TwoB"); list_dW2TwoB->Add(fCalculator.cumu_dW2TwoB); fStandardList->Add(list_dW2TwoB);
  TList* list_dW4     = new TList(); list_dW4    ->SetName("dW4"    ); list_dW4    ->Add(fCalculator.cumu_dW4    ); fStandardList->Add(list_dW4    );
  TList* list_dW4Four = new TList(); list_dW4Four->SetName("dW4Four"); list_dW4Four->Add(fCalculator.cumu_dW4Four); fStandardList->Add(list_dW4Four);

  fCalculator.cumu_dW4FourTwo  = new THnD("cumu_dW4FourTwo" , "cumu_dW4FourTwo" , dimensions, sc_dbins, sc_dmin, sc_dmax) ;
  fCalculator.cumu_dW4ThreeTwo = new THnD("cumu_dW4ThreeTwo", "cumu_dW4ThreeTwo", dimensions, sc_dbins, sc_dmin, sc_dmax) ;

  fCalculator.cumu_dWTwoTwoN   = new THnD("cumu_dWTwoTwoN"  , "cumu_dWTwoTwoN"  , dimensions, dbins, dmin, dmax) ; // Numerator of R_{n,n; 2}
  fCalculator.cumu_dWTwoTwoD   = new THnD("cumu_dWTwoTwoD"  , "cumu_dWTwoTwoD"  , dimensions, dbins, dmin, dmax) ; // Denominator of R_{n,n; 2}

  TList* list_dW4FourTwo  = new TList(); list_dW4FourTwo ->SetName("dW4FourTwo" ); list_dW4FourTwo ->Add(fCalculator.cumu_dW4FourTwo ); fMixedList->Add(list_dW4FourTwo );
  TList* list_dW4ThreeTwo = new TList(); list_dW4ThreeTwo->SetName("dW4ThreeTwo"); list_dW4ThreeTwo->Add(fCalculator.cumu_dW4ThreeTwo); fMixedList->Add(list_dW4ThreeTwo);
  TList* list_dWTwoTwoN   = new TList(); list_dWTwoTwoN  ->SetName("dWTwoTwoN"  ); list_dWTwoTwoN  ->Add(fCalculator.cumu_dWTwoTwoN  ); fMixedList->Add(list_dWTwoTwoN  );
  TList* list_dWTwoTwoD   = new TList(); list_dWTwoTwoD  ->SetName("dWTwoTwoD"  ); list_dWTwoTwoD  ->Add(fCalculator.cumu_dWTwoTwoD  ); fMixedList->Add(list_dWTwoTwoD  );


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
  fCalculator.fSettings = fSettings;
  fUtil.fSettings = fSettings;
  fCalculator.fUtil = fUtil;

  Bool_t isgoodrun = kTRUE;
  if (!fSettings.mc){
    isgoodrun = fUtil.IsGoodRun(fInputEvent->GetRunNumber());
  }
  if (fSettings.doNUA) fSettings.nua_runnumber = fUtil.GetNUARunNumber(fInputEvent->GetRunNumber());

  // Get the event validation object
  AliForwardTaskValidation* ev_val = dynamic_cast<AliForwardTaskValidation*>(this->GetInputData(1));
  if (!ev_val->IsValidEvent() || !isgoodrun){
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
  //Int_t centBin = fCalculator.cumu_rW2->GetAxis(5)->FindBin(cent);


  fUtil.FillData(refDist,centralDist,forwardDist);

  Double_t zvertex = fUtil.GetZ();


  /*
  if (fSettings.a5){
    fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE, true,false);
    fCalculator.CumulantsAccumulate(centralDist, cent, zvertex,kFALSE,true,false);
  }
  else{
  */
  if (fSettings.ref_mode & fSettings.kFMDref) {
    fCalculator.CumulantsAccumulate(forwardDist, cent, zvertex,kTRUE,true,true);
  }
  else {
    fCalculator.CumulantsAccumulate(refDist, cent, zvertex,kFALSE,true,false);
  }
  
  fCalculator.CumulantsAccumulate(centralDist, cent, zvertex,kFALSE,false,true);  

  UInt_t randomInt = fRandom.Integer(fSettings.fnoSamples);

  fCalculator.saveEvent(cent, zvertex,  randomInt, 0);

  fCalculator.reset();
  centralDist->Reset();
  
  if (!(fSettings.ref_mode & fSettings.kFMDref)) refDist->Reset();
  if ((fSettings.mc && fSettings.use_primaries_fwd) || (fSettings.mc && fSettings.esd)) {
    forwardDist->Reset();
  }
  /*
  fSettings.track_sample++;
  if (fSettings.track_sample == fSettings.fnoSamples) fSettings.track_sample = 0;
  fCentCounter[centBin-1]++;
  if (fCentCounter[centBin-1] == 10) fCentCounter[centBin-1] = 0;
  */
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
