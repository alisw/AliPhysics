/**********************************************
 *
 *   Class meant to perform calibration and
 *   write an OADB file containing histos and
 *   Event selection criteria used
 *
 *********************************************
 *
 * --- Revised version ---
 *
 *  - David Dobrigkeit Chinellato
 *  - Alberica Toia
 *  - Tatiana Drozhzhova
 *
 **********************************************/

#include "AliMultSelectionCuts.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "TList.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TArrayL64.h"
#include "TArrayF.h"

ClassImp(AliMultSelectionCalibrator);

AliMultSelectionCalibrator::AliMultSelectionCalibrator() : TNamed(),
fInput(0), fSelection(0), lDesiredBoundaries(0), lNDesiredBoundaries(0),
fRunToUseAsDefault(-1), fMaxEventsPerRun(1e+9), fCheckTriggerType(kFALSE),
fTrigType(AliVEvent::kAny), fPrefilterOnly(kFALSE), fFiredTrigString(""),
fNRunRanges(0), fRunRangesMap(), fMultSelectionList(0),
fInputFileName(""), fBufferFileName("buffer.root"),
fOutputFileName(""), fMultSelectionCuts(0), fCalibHists(0)
{
  // Constructor
  Bool_t lVarDoAutocalib[] = {
    kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
    kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE
  };
  
  for(Int_t iVar = 0; iVar<AliMultInput::kNVariables; iVar++)
  VarDoAutocalib[iVar] = lVarDoAutocalib[iVar];
  
  // Create Event Selector
  fMultSelectionCuts = new AliMultSelectionCuts();
  fMultSelectionCuts -> Print();
  
  //Basic I/O for MultSelection framework
  fInput     = new AliMultInput();
  
  //List of objects to use as AliMultSelection
  fMultSelectionList = new TList();
  
  //Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists -> SetOwner(kTRUE);
  
  //Initialize
  for( Int_t iv=0;iv<1000;iv++){ fFirstRun[iv] = 0; }
  for( Int_t iv=0;iv<1000;iv++){ fLastRun[iv] = 0; }
}

AliMultSelectionCalibrator::AliMultSelectionCalibrator(const char * name, const char * title):
TNamed(name,title),
fInput(0), fSelection(0), lDesiredBoundaries(0), lNDesiredBoundaries(0),
fRunToUseAsDefault(-1), fMaxEventsPerRun(1e+9), fCheckTriggerType(kFALSE),
fTrigType(AliVEvent::kAny), fPrefilterOnly(kFALSE), fFiredTrigString(""),
fNRunRanges(0), fRunRangesMap(), fMultSelectionList(0),
fInputFileName(""), fBufferFileName("buffer.root"),
fOutputFileName(""), fMultSelectionCuts(0), fCalibHists(0)
{
  // Named Constructor
  Bool_t lVarDoAutocalib[] = {
    kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
    kTRUE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kFALSE,
    kFALSE, kFALSE,
    kTRUE, kTRUE, kTRUE, kTRUE, kTRUE,
    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE, kFALSE,
    kFALSE, kFALSE, kFALSE
  };
  
  for(Int_t iVar = 0; iVar<AliMultInput::kNVariables; iVar++)
  VarDoAutocalib[iVar] = lVarDoAutocalib[iVar];
  
  // Create Event Selector
  fMultSelectionCuts = new AliMultSelectionCuts();
  
  //Set Standard Cuts (warning: this is pp-like!)
  fMultSelectionCuts -> SetVzCut (10.0);
  fMultSelectionCuts -> SetTriggerCut(kTRUE);
  fMultSelectionCuts -> SetINELgtZEROCut(kTRUE);
  fMultSelectionCuts -> SetTrackletsVsClustersCut(kTRUE);
  fMultSelectionCuts -> SetRejectPileupInMultBinsCut(kTRUE);
  fMultSelectionCuts -> SetVertexConsistencyCut(kTRUE);
  fMultSelectionCuts -> SetNonZeroNContribs(kFALSE);
  fMultSelectionCuts -> SetIsNotAsymmetricInVZERO(kFALSE);
  fMultSelectionCuts -> SetIsNotIncompleteDAQ(kFALSE);
  fMultSelectionCuts -> SetHasGoodVertex2016(kFALSE);
  
  //Basic I/O for MultSelection framework
  fInput     = new AliMultInput();
  
  //List of objects to use as AliMultSelection
  fMultSelectionList = new TList();
  
  //Make sure the TList owns its objects
  fCalibHists = new TList();
  fCalibHists -> SetOwner(kTRUE);
  
  //Initialize
  for( Int_t iv=0;iv<1000;iv++){ fFirstRun[iv] = 0; }
  for( Int_t iv=0;iv<1000;iv++){ fLastRun[iv] = 0; }
}
//________________________________________________________________
AliMultSelectionCalibrator::~AliMultSelectionCalibrator() {
  // Destructor
  
  if ( fInput ) {
    delete fInput;
    fInput = 0x0;
  }
  if ( fSelection ) {
    delete fSelection;
    fSelection = 0x0;
  }
  if ( fMultSelectionCuts ) {
    delete fMultSelectionCuts;
    fMultSelectionCuts = 0x0;
  }
  if ( fCalibHists ) {
    delete fCalibHists;
    fCalibHists = 0x0;
  }
}

//________________________________________________________________
void AliMultSelectionCalibrator::AddRunRange ( Int_t lFirst, Int_t lLast, AliMultSelection *lMultSelProvided ){
  //Add mapping : all runs in range go to current value of fNRunRanges
  //Ease of access
  fFirstRun[fNRunRanges] = lFirst;
  fLastRun [fNRunRanges] = lLast;
  for( Int_t iRun = lFirst; iRun<=lLast; iRun++){
    fRunRangesMap.insert( std::pair<int,int>(iRun,fNRunRanges));
  }
  //Add the provided AliMultSelection object among those to be used
  fMultSelectionList -> Add ( lMultSelProvided ) ;
  
  AliInfoF("Added Run Range #%i (%i - %i)", (Int_t)fNRunRanges, lFirst, lLast) ;
  AliInfoF("Added AliMultSelection object. Existing objects: %i", fMultSelectionList->GetEntries()) ;
  fNRunRanges++;
}
//________________________________________________________________
Bool_t AliMultSelectionCalibrator::Calibrate() {
  // Function meant to generate calibration OADB
  //
  // --- input : fInputFileName, containing a TTree object
  // --- output: fOutputFileName, containing OABD object
  //
  // Steps involved:
  //  (1) Set up basic I/O
  //  (2) Detect Runs From Input File
  //  (3) Determine Averages
  //     (3a) Create run-by-run buffer files with averages
  //  (4) Determine Quantile Boundaries
  //     (4a) Create run-by-run buffer files (requires averages)
  //     (4b) Compute Quantile Boundaries for all estimators
  //  (4) Save Quantiles + AliMultSelectionCuts to OADB File
  
  
  cout<<"=== STARTING CALIBRATION PROCEDURE ==="<<endl;
  cout<<" * Input File.....: "<<fInputFileName.Data()<<endl;
  cout<<" * Output File....: "<<fOutputFileName.Data()<<endl;
  cout<<endl;
  cout<<" * Event Selection Peformed: "<<endl;
  fMultSelectionCuts -> Print();
  cout<<endl;
  
  // STEP 1: Basic I/O
  cout<<"(1) Opening File"<<endl;
  
  //Open File
  TFile *fInputFile = TFile::Open( fInputFileName.Data(), "READ");
  if(!fInputFile) {
    AliWarningF("File %s not found!", fInputFileName.Data() );
    return kFALSE;
  }
  //Locate TTree object
  TTree* fTree = (TTree*)fInputFile->FindObjectAny("fTreeEvent");
  if(!fTree) {
    AliWarning("fTreeEvent object not found!" );
    return kFALSE;
  }
  
  //Event Selection Variables
  Bool_t fEvSel_IsNotPileupInMultBins      = kFALSE ;
  Bool_t fEvSel_Triggered                  = kFALSE ;
  Bool_t fEvSel_INELgtZERO                 = kFALSE ;
  Bool_t fEvSel_PassesTrackletVsCluster    = kFALSE ;
  Bool_t fEvSel_HasNoInconsistentVertices  = kFALSE ;
  Bool_t fEvSel_IsNotAsymmetricInVZERO     = kFALSE ;
  Bool_t fEvSel_IsNotIncompleteDAQ         = kFALSE ;
  Bool_t fEvSel_HasGoodVertex2016          = kFALSE ;
  Int_t fRunNumber;
  
  //FIXME/CAUTION: non-zero if using tree without that branch
  Int_t fnContributors = 1000;
  
  UInt_t fEvSel_TriggerMask; //! save full info for checking later
  
  //SetBranchAddresses for event Selection Variables
  //(multiplicity related will be done automatically!)
  fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
  fTree->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
  fTree->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
  fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
  fTree->SetBranchAddress("fEvSel_TriggerMask",&fEvSel_TriggerMask);
  fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
  fTree->SetBranchAddress("fRunNumber",&fRunNumber);
  fTree->SetBranchAddress("fnContributors", &fnContributors);
  fTree->SetBranchAddress("fEvSel_IsNotAsymmetricInVZERO", &fEvSel_IsNotAsymmetricInVZERO);
  fTree->SetBranchAddress("fEvSel_IsNotIncompleteDAQ", &fEvSel_IsNotIncompleteDAQ);
  fTree->SetBranchAddress("fEvSel_HasGoodVertex2016", &fEvSel_HasGoodVertex2016);
  
  TString *fFiredTriggerClasses = new TString();
  fTree->SetBranchAddress("fFiredTriggerClasses",&fFiredTriggerClasses);
  
  //============================================================
  // Auto-configure Input
  //============================================================
  
  Bool_t lAutoDiscover = kFALSE;
  
  if ( fInput->GetNVariables() < 1 ){
    cout<<"Error: No Input Variables configured!"<<endl;
    cout<<"The simplest way to get rid of this problem is to remember to call SetupStandardInput()!"<<endl;
    return kFALSE; //failure to calibrate
  }
  
  if ( fMultSelectionList->GetEntries() == 0 ){
    AliInfo("===============================================");
    AliInfo(" Calibrator invoked without run mappings");
    AliInfo(" Auto run discovery mode will be used!");
    AliInfo("===============================================");
    lAutoDiscover = kTRUE;
  }
  
  if ( lAutoDiscover && !fSelection ) {
    cout<<"Error: no default AliMultSelection defined!"<<endl;
    cout<<"The simplest way to get rid of this problem is to remember to call SetMultSelection(...)!"<<endl;
    return kFALSE; //failure to calibrate
  }
  
  //Binding to input variables
  for(Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++) {
    if( !fInput->GetVariable(iVar)->IsInteger() ) {
      fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValue());
    } else {
      fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValueInteger());
    }
  }
  
  Long64_t lNEv = fTree->GetEntries();
  cout<<"(1) File opened, event count is "<<lNEv<<endl;
  
  cout<<"(2) Creating buffer, computing averages"<<endl;
  const int lMax = 1000;
  const int lMaxQuantiles = 10000;
  Int_t lRunNumbers[lMaxQuantiles];
  Long_t lRunStats[lMaxQuantiles];
  
  for( Int_t ix=0; ix<1000;ix++){
    lRunNumbers[ix] = 0;
    lRunStats[ix] = 0;
  }
  
  Int_t lNRuns = 0;
  Bool_t lNewRun = kTRUE;
  Int_t lThisRunIndex = -1;
  //Buffer file with run-by-run TTree objects needed for later processing
  
  TFile *fOutput = new TFile (fBufferFileName.Data(), "RECREATE");
  TTree *sTree[lMaxQuantiles];
  cout<<"Creating Trees..."<<endl;
  //N.B. No need to Exceed Run Ranges in Calibration Code here!
  Int_t lNTrees = 0;
  if( !lAutoDiscover ){
    lNTrees = fNRunRanges;
    lNRuns = fNRunRanges;
  }else{
    lNTrees = lMax;
  }
  for(Int_t iRun=0; iRun<lNTrees; iRun++) {
    sTree[iRun] = new TTree(Form("sTree%i",iRun),Form("sTree%i",iRun));
    
    //useful for debugging / cross-checking
    sTree[iRun]->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    
    for( Int_t iQvar = 0; iQvar<fInput->GetNVariables(); iQvar++) {
      if( !fInput->GetVariable(iQvar)->IsInteger() ) {
        sTree[iRun]->Branch(Form("%s", fInput->GetVariable(iQvar)->GetName()  ),
                            &fInput->GetVariable(iQvar)->GetRValue(),Form("%s/F",fInput->GetVariable(iQvar)->GetName()));
      } else {
        sTree[iRun]->Branch(Form("%s", fInput->GetVariable(iQvar)->GetName()  ),
                            &fInput->GetVariable(iQvar)->GetRValueInteger(),Form("%s/I",fInput->GetVariable(iQvar)->GetName()));
      }
    }
  }
  
  //const int lNEstimators = fSelection->GetNEstimators();
  
  const int lNEstimators = 50; //this is the MAX VALUE!
  
  //For computing average values of estimators
  Double_t lAvEst[lNEstimators][lMax];
  //For computing extreme values (useful for integer calibration mode)
  Double_t lMaxEst[lNEstimators][lMax];
  Double_t lMinEst[lNEstimators][lMax];
  
  //Sanity check. If insane, add kNoCalib histogram
  Bool_t lInsane[lNEstimators][lMax];
  
  //Index of first value above anchor point threshold
  Long64_t lAnchorEst[lNEstimators][lMax];
  
  //Vertex Z profiles
  TProfile *hCalibVtx[1000][AliMultInput::kNVariables];
  for(Int_t ii=0; ii<1000; ii++){
    for(Int_t jj=0; jj<AliMultInput::kNVariables; jj++){
      hCalibVtx[ii][jj]=0x0; //init to null to be sure
    }
  }
  
  cout<<"Creating histograms"<<endl;
  if( !lAutoDiscover ){
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
      for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
        hCalibVtx[iRun][iVar] = new TProfile(Form("hCalibVtx_%i_%s",fFirstRun[iRun],AliMultInput::VarName[iVar].Data()),"",100, -20, 20);
        hCalibVtx[iRun][iVar] ->SetDirectory(0);
      }
    }
  }
  
  for(Long_t iEst=0; iEst<lNEstimators; iEst++) {
    for(Long_t iRun=0; iRun<lMax; iRun++) lAvEst[iEst][iRun] = 0;
    for(Long_t iRun=0; iRun<lMax; iRun++) lMaxEst[iEst][iRun] = -1e+3;
    for(Long_t iRun=0; iRun<lMax; iRun++) lMinEst[iEst][iRun] = 1e+6; //not more than a million, I hope?
    for(Long_t iRun=0; iRun<lMax; iRun++) lAnchorEst[iEst][iRun] = -1; //invalid index
    for(Long_t iRun=0; iRun<lMax; iRun++) lInsane[iEst][iRun] = kFALSE; //we're nice people. We assume no insanity unless there's proof otherwise
  }
  
  //Add Timer
  TStopwatch* timer = new TStopwatch();
  timer->Start ( kTRUE );
  
  //Compute events-per-hour performance metric
  Double_t lEventsPerSecond = 0;
  
  AliMultVariable *lVtxZLocalPointer = fInput -> GetVariable("fEvSel_VtxZ");
  if( ! lVtxZLocalPointer ) cout<<"PROBLEM singling out vtx-Z variable"<<endl; 
  
  for(Long64_t iEv = 0; iEv<fTree->GetEntries(); iEv++) {
    if ( iEv % 100000 == 0 ) {
      Double_t complete = 100. * ( double ) ( iEv ) / ( double ) ( fTree->GetEntries() );
      cout << "\r" << "Event # " << iEv << "/" << fTree->GetEntries() << " (" << complete << "%, Time Left: ";
      timer->Stop();
      Double_t time = timer->RealTime();
      
      //events per hour:
      lEventsPerSecond = ( ( Double_t ) ( iEv ) ) /time;
      
      timer->Start ( kFALSE );
      Double_t secondsperstep = time / ( Double_t ) ( iEv+1 );
      Double_t secondsleft = ( Double_t ) ( fTree->GetEntries()-iEv-1 ) * secondsperstep;
      Long_t minutesleft = ( Long_t ) ( secondsleft / 60. );
      secondsleft = ( Double_t ) ( ( Long_t ) ( secondsleft ) % 60 );
      cout << minutesleft << "min " << secondsleft << "s, working at "<<lEventsPerSecond<<" Events/s..." << flush;
    }
    fTree->GetEntry(iEv);
    //Perform Event selection
    Bool_t lSaveThisEvent = kTRUE; //let's be optimistic
    
    //Apply trigger mask (will only work if not kAny)
    Bool_t isSelected = 0;
    isSelected = fEvSel_TriggerMask & fTrigType;
    if(!isSelected && fCheckTriggerType) lSaveThisEvent = kFALSE;
    
    if(fFiredTrigString.EqualTo("")==kFALSE) {
      if (fFiredTriggerClasses->Contains( fFiredTrigString.Data() ) == kFALSE ) lSaveThisEvent = kFALSE;
    }
    
    //Check Selections as they are in the fMultSelectionCuts Object
    if( fMultSelectionCuts->GetTriggerCut()    && ! fEvSel_Triggered  ) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO ) lSaveThisEvent = kFALSE;
    //if( TMath::Abs( lVtxZLocalPointer->GetValue() ) > fMultSelectionCuts->GetVzCut()      ) lSaveThisEvent = kFALSE;
    //ADD ME HERE: Tracklets Vs Clusters Cut?
    if( fMultSelectionCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins    ) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster  ) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetNonZeroNContribs()          &&  fnContributors < 1 ) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetIsNotAsymmetricInVZERO()    && ! fEvSel_IsNotAsymmetricInVZERO) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetIsNotIncompleteDAQ()        && ! fEvSel_IsNotIncompleteDAQ) lSaveThisEvent = kFALSE;
    if( fMultSelectionCuts->GetHasGoodVertex2016()         && ! fEvSel_HasGoodVertex2016) lSaveThisEvent = kFALSE;
    Int_t lIndex = -1;
    if ( !lAutoDiscover ){
      //Consult map for run range equivalency
      if ( fRunRangesMap.find( fRunNumber ) != fRunRangesMap.end() ) {
        lIndex = fRunRangesMap[ fRunNumber ];
        //Create vertex-Z calibration object if it does not exist
//        for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
//          cout<<"Creating hCalibVtx["<<lIndex<<"]["<<iVar<<"]"<<endl;
//          if( !hCalibVtx[lIndex][iVar] ){
//            hCalibVtx[lIndex][iVar] = new TProfile(Form("hCalibVtx_%i_%s",fRunNumber,AliMultInput::VarName[iVar].Data()),"",100, -20, 20);
//            hCalibVtx[lIndex][iVar] ->SetDirectory(0);
//          }
//        }
      }else{
        lSaveThisEvent = kFALSE;
      }
    }else{
      lNewRun = kTRUE;
      for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        if( lRunNumbers[iRun] == fRunNumber ) {
          lNewRun = kFALSE;
          lIndex = iRun;
        }
      }
      if( lNewRun == kTRUE ) {
        cout<<endl<<"(Autodiscover) New Run Found: "<<fRunNumber<<", added as #"<<lNRuns<<" (so far: "<<lNRuns<<" runs)"<<endl;
        //Adding vertex-Z profiles for this run
        for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
          cout<<"Creating hCalibVtx["<<lNRuns<<"]["<<iVar<<"]"<<endl;
          hCalibVtx[lNRuns][iVar] = new TProfile(Form("hCalibVtx_%i_%s",fRunNumber,AliMultInput::VarName[iVar].Data()),"",100, -20, 20);
          hCalibVtx[lNRuns][iVar] ->SetDirectory(0);
        }
        
        //Add to Map
        fRunRangesMap.insert( std::pair<int,int>(fRunNumber,lNRuns));
        lRunNumbers[lNRuns] = fRunNumber;
        lIndex = lNRuns;
        lNRuns++;
        fNRunRanges++;
      }
    }
    if ( lSaveThisEvent ) {
      if( sTree[lIndex]->GetEntries()<fMaxEventsPerRun ){
        
        //Fill vertex-Z profiles for all Z
        for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
          if(VarDoAutocalib[iVar]){
            AliMultVariable *v1 = fInput->GetVariable(iVar);
            if(!v1) cout<<"PROBLEM finding this variable"<<endl;
            if(!hCalibVtx[lIndex][iVar]) cout<<"Undefined vertex-Z calib histo at indices "<<lIndex<<", "<<iVar<<", please debug"<<endl;
            hCalibVtx[lIndex][iVar]->Fill( lVtxZLocalPointer->GetValue(), v1->IsInteger() ? v1->GetValueInteger() : v1->GetValue() );
          }
        }
        
        //Only fill if this passes the vertex cut (default behaviour)
        if ( TMath::Abs( lVtxZLocalPointer->GetValue() ) < fMultSelectionCuts->GetVzCut() ) sTree [ lIndex ] -> Fill();
        
      }
    }
  }
  cout<<endl;
  
  //Write buffer to file
  for(Int_t iRun=0; iRun<lNRuns; iRun++) sTree[iRun]->Write();
  for(Int_t iRun=0; iRun<lNRuns; iRun++) lRunStats[iRun] = sTree[iRun]->GetEntries();
  
  //Write copies of the vertex-Z histos for cross-checks at this stage already
  for(Int_t iRun=0; iRun<lNRuns; iRun++) {
    for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
      if(lRunStats[iRun] < 5.0e+5 ) hCalibVtx[iRun][iVar]->Rebin(2);
      if(lRunStats[iRun] < 1.5e+5 ) hCalibVtx[iRun][iVar]->Rebin(2);
      if(VarDoAutocalib[iVar])
        hCalibVtx[iRun][iVar] -> Write(Form("Debug_VertexZ_RunIdx%i_VarIdx%i",iRun,iVar));
    }
  }
  
  if(!lAutoDiscover){
    cout<<"(3) Inspect Run Ranges and corresponding statistics: "<<endl;
    for(Int_t iRun = 0; iRun<fNRunRanges; iRun++) {
      cout<<" --- Range #"<<iRun<<", ("<<fFirstRun[iRun]<<" - "<<fLastRun[iRun]<<"), N(events) = "<<sTree[iRun]->GetEntries()<<endl;
    }
    cout<<endl;
  }else{
    cout<<"(3) Inspect Runs and corresponding statistics: "<<endl;
    for(Int_t iRun = 0; iRun<fNRunRanges; iRun++) {
      cout<<" --- Run #"<<iRun<<", (#"<<lRunNumbers[iRun]<<"), N(events) = "<<sTree[iRun]->GetEntries()<<endl;
    }
    cout<<endl;
  }
  
  if( fPrefilterOnly ){
    cout<<"Won't do anything else! But you should have filtered trees now for debugging..."<<endl;
    fOutput->Write();
    return 0;
  }
  
  //FIXME Receive as parameter from the test macro
  Double_t lNrawBoundaries[1000];
  Double_t lMiddleOfBins[1000];
  
  for( Long_t lB=1; lB<lNDesiredBoundaries; lB++) {
    //place squarely at the middle to ensure it's all fine
    lMiddleOfBins[lB-1] = 0.5*(lDesiredBoundaries[lB]+lDesiredBoundaries[lB-1]);
  }
  
  //Histograms to store calibration information
  TH1F *hCalib[1000][lNEstimators];
  
  cout<<"(4) Look at average values"<<endl;
  for(Int_t iRun=0; iRun<fNRunRanges; iRun++) {
    //Contextualize AliMultSelection for this run
    if ( !lAutoDiscover ) fSelection = (AliMultSelection*) fMultSelectionList->At(iRun);
    
    // Calibration pre-optimization and setup
    fSelection->Setup ( fInput );
    
    const Int_t lNEstimatorsThis = fSelection->GetNEstimators();
    
    const Long64_t ntot = (Long64_t) sTree[iRun]->GetEntries();
    if ( !lAutoDiscover ){
      cout<<"--- Processing run range "<<fFirstRun[iRun]<<"-"<<fLastRun[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<"), with "<<ntot<<" events..."<<endl;
    }else{
      cout<<"--- Processing run "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<"), with "<<ntot<<" events..."<<endl;
    }
    if( ntot < 1 ){ cout<<"Sample empty! Skipping..."<<endl; continue; }
    sTree[iRun]->SetEstimate(ntot+1);
    //Cast Run Number into drawing conditions
    
    fSelection->Setup ( fInput );
    
    cout<<"--- Initializing. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    fInput->ClearVtxZ();
    cout<<"--- Cleared. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++){
      //cout<<"--- Adding var: "<<AliMultInput::VarName[iVar].Data()<<", pointer "<<hCalibVtx[iRun][iVar]<<endl;
      fInput->AddVtxZ( hCalibVtx[iRun][iVar] );
    }
    cout<<"--- Setting up input maps for vertex-Z correction..."<<endl;
    fInput->SetupAutoVtxZCorrection();
    cout<<"--- Finished initalization. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    
    timer->Reset();
    timer->Start ( kTRUE );
    lEventsPerSecond = 0;
    
    for( Long64_t iEntry=0; iEntry<ntot; iEntry++) {
      if ( iEntry % 10000 == 0 ) {
        cout << "\r" << "--- Looping over tree [ "<<Form("%.3f",100.*iEntry/ntot)<<"% ] [ ETA: "<< flush;
        timer->Stop();
        Double_t time = timer->RealTime();
        
        //events per hour:
        lEventsPerSecond = ( ( Double_t ) ( iEntry ) ) /time;
        
        timer->Start ( kFALSE );
        Double_t secondsperstep = time / ( Double_t ) ( iEntry+1 );
        Double_t secondsleft = ( Double_t ) ( ntot-iEntry-1 ) * secondsperstep;
        Long_t minutesleft = ( Long_t ) ( secondsleft / 60. );
        secondsleft = ( Double_t ) ( ( Long_t ) ( secondsleft ) % 60 );
        cout << minutesleft << "min " << secondsleft << "s ] [ Speed: "<<lEventsPerSecond<<" events/s ]" << flush;
      }
      sTree[iRun]->GetEntry(iEntry); //Get all variables
      fSelection->Evaluate ( fInput ); //Evaluate based on current variables
      for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        Double_t lThisVal = fSelection->GetEstimator(iEst)->GetValue();
        lAvEst[iEst][iRun] += lThisVal;
        if( lThisVal < lMinEst[iEst][iRun] ) {
          lMinEst[iEst][iRun] = lThisVal;
        }
        if( lThisVal > lMaxEst[iEst][iRun] ) {
          lMaxEst[iEst][iRun] = lThisVal;
        }
      }
    }
    cout<<endl;
    
    for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      if( sTree[iRun]->GetEntries() < 1 ) {
        lAvEst[iEst][iRun] = -1;
      } else {
        lAvEst[iEst][iRun] /= ( (Double_t) (sTree[iRun]->GetEntries()) );
      }
      cout<<" Min = "<<lMinEst[iEst][iRun]<<", Max = "<<lMaxEst[iEst][iRun]<<", Av = "<<lAvEst[iEst][iRun]<<endl;
      
      if ( TMath::Abs( lMinEst[iEst][iRun] - lMaxEst[iEst][iRun] ) < 1e-6 ){
        lInsane[iEst][iRun] = kTRUE; //No valid information to do calibration, please be careful !
      }
    }
  }
  //might be needed
  Long64_t lAcceptedEvents;
  
  //=========================================
  // Determine Calibration Information
  //=========================================
  
  //Open Run 3 - compliant calibration object, create structure
  TString lRun3File = fOutputFileName.Data();
  lRun3File.Prepend("Run3-");
  TFile *fRun3 = new TFile(lRun3File.Data(), "RECREATE");
  //Do not use custom objects, ROOT compliance only
  
  fRun3->cd();
  //Master list to store all calibration objects
  TList *lR3List = new TList();
  TH1F *hR3EventSelection[500];
  TProfile *hR3CalibDataVtx[500][AliMultInput::kNVariables];
  TH1F *hR3DefData[500][lNEstimators];
  TH1F *hR3CalibData[500][lNEstimators];
  TH1F *hR3EventSelectionDefault;
  TProfile *hR3CalibDataVtxDefault[AliMultInput::kNVariables];
  TH1F *hR3DefDataDefault[lNEstimators];
  TH1F *hR3CalibDataDefault[lNEstimators];
  
  //Open output OADB file, generate everything within loop
  TFile * f = new TFile (fOutputFileName.Data(), "recreate");
  AliOADBContainer * oadbContMS = new AliOADBContainer("MultSel");
  
  AliOADBMultSelection * oadbMultSelection = 0x0;
  AliMultSelectionCuts * cuts = 0x0;
  AliMultSelection     * fsels = 0x0;
  
  //Actual Calibration Histograms
  TH1F * hCalibData[lNEstimators];
  TProfile * hCalibDataVtx[AliMultInput::kNVariables];
  
  cout<<"(5) Generate Boundaries through a loop in all desired estimators"<<endl;
  for(Int_t iRun=0; iRun<fNRunRanges; iRun++) {
    cout<<"--- Initializing. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    fInput->ClearVtxZ();
    cout<<"--- Cleared. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    for(Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++){
      cout<<"--- Adding var: "<<AliMultInput::VarName[iVar].Data()<<", pointer "<<hCalibVtx[iRun][iVar]<<" for hCalibVtx["<<iRun<<"]["<<iVar<<"]"<<endl;
      if(hCalibVtx[iRun][iVar]) {
        fInput->AddVtxZ( hCalibVtx[iRun][iVar] );
      }else{
        cout<<"--- Access to hCalibVtx["<<iRun<<"]["<<iVar<<"] failed, presuming that it does not exist!"<<endl;
      }
    }
    cout<<"--- Setting up input maps for vertex-Z correction..."<<endl;
    fInput->SetupAutoVtxZCorrection();
    cout<<"--- Finished initalization. Number of elements: "<<fInput->GetNVtxZ()<<endl;
    
    //Contextualize AliMultSelection for this run
    if ( !lAutoDiscover ) fSelection = (AliMultSelection*) fMultSelectionList->At(iRun);
    
    const Int_t lNEstimatorsThis = fSelection->GetNEstimators();
    
    const Int_t ntot = (Int_t) sTree[iRun]->GetEntries();
    if ( !lAutoDiscover ){
      cout<<"--- Processing run range "<<fFirstRun[iRun]<<"-"<<fLastRun[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<"), with "<<ntot<<" events..."<<endl;
    }else{
      cout<<"--- Processing run "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<"), with "<<ntot<<" events..."<<endl;
    }
    if( ntot < 1 ){ cout<<"Sample empty! Skipping..."<<endl; continue; }
    sTree[iRun]->SetEstimate(ntot+1);
    
    
    //Memory allocation: this may be dangerous in computers with little memory, but it will make the process faster
    //TArrayL64 index(ntot);
    //TArrayF values(ntot);
    
    Int_t **index = new Int_t*[lNEstimatorsThis];
    Float_t **values = new Float_t*[lNEstimatorsThis];
    Long_t lAcceptedEventArray[lNEstimatorsThis];
    
    for(Int_t ii=0; ii<lNEstimatorsThis; ii++){
      index[ii] = new Int_t[ntot];
      values[ii] = new Float_t[ntot];
      lAcceptedEventArray[ii] = 0;
    }
    
    timer->Reset();
    timer->Start ( kTRUE );
    lEventsPerSecond = 0;
    
    cout << "\r" << "--- Looping over tree [ 0% ] "<< flush;
    for( Long64_t iEntry=0; iEntry<ntot; iEntry++) {
      if ( iEntry % 10000 == 0 ) {
        cout << "\r" << "--- Looping over tree [ "<<Form("%.3f",100.*iEntry/ntot)<<"% ] [ ETA: "<< flush;
        timer->Stop();
        Double_t time = timer->RealTime();
        
        //events per hour:
        lEventsPerSecond = ( ( Double_t ) ( iEntry ) ) /time;
        
        timer->Start ( kFALSE );
        Double_t secondsperstep = time / ( Double_t ) ( iEntry+1 );
        Double_t secondsleft = ( Double_t ) ( ntot-iEntry-1 ) * secondsperstep;
        Long_t minutesleft = ( Long_t ) ( secondsleft / 60. );
        secondsleft = ( Double_t ) ( ( Long_t ) ( secondsleft ) % 60 );
        cout << minutesleft << "min " << secondsleft << "s ] [ Speed: "<<lEventsPerSecond<<" events/s ]" << flush;
      }
      sTree[iRun]->GetEntry(iEntry); //Get all variables
      fSelection->Evaluate ( fInput ); //Evaluate based on current variables
      for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        values[iEst][iEntry] = fSelection->GetEstimator(iEst)->GetValue();
        if( !fSelection->GetEstimator(iEst)->GetUseAnchor() ){
          lAcceptedEventArray[iEst]++;
        }else{
          if(values[iEst][iEntry] > fSelection->GetEstimator(iEst)->GetAnchorPoint())
            lAcceptedEventArray[iEst]++;
        }
      }
    }
    cout<<endl;
    
    //Cast Run Number into drawing conditions
    for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      lAcceptedEvents = lAcceptedEventArray[iEst];
      cout<<"--- Sorting estimator "<<fSelection->GetEstimator(iEst)->GetName()<<"..."<<flush;
      
      Int_t *index1 = index[iEst];
      Float_t *values1 = values[iEst];
      TMath::Sort(ntot, values1, index1 );
      cout<<" Done! Getting Boundaries... "<<flush;
      //Special override in case anchored estimator
      
      lNrawBoundaries[0] = 0.0; //Defined OK even if anchored
      //Overwrite lower boundary in case this has a negative minimum...
      if ( lMinEst[iEst][iRun] < 0 ) {
        lNrawBoundaries[0] = lMinEst[iEst][iRun];
        cout<<"Min Value Override, Negative..."<<flush;
      }
      
      for( Long_t lB=1; lB<lNDesiredBoundaries; lB++) {
        Long64_t position = (Long64_t) ( 0.01 * ((Double_t)(ntot)* lDesiredBoundaries[lB] ) );
        if( fSelection->GetEstimator(iEst)->GetUseAnchor() && ntot != 0 ){
          //Make sure index position lAnchorEst corresponds to lAnchorPercentile
          Double_t lAnchorPercentile = (Double_t) fSelection->GetEstimator(iEst)->GetAnchorPercentile();
          Double_t lFractionAccepted = (((Double_t) lAcceptedEvents )/((Double_t) ntot));
          Double_t lScalingFactor    = lFractionAccepted/((0.01)*lAnchorPercentile);
          //Make sure: if AnchorPercentile requested, cut at AnchorPoint
          position = (Long64_t) ( ( 0.01 * ((Double_t)(ntot)* lDesiredBoundaries[lB] ) ) * lScalingFactor );
          if(position > ntot-1 ) position = ntot-1; //protection !
        }
        //sTree[iRun]->GetEntry( index[position] );
        //Calculate the estimator with this input, please
        //fSelection->Evaluate ( fInput );
        //fSelection->PrintInfo();
        lNrawBoundaries[lB] = values[iEst][index[iEst][position]];
        //cout<<sTree[iRun]->GetV1()<<" "<<lNrawBoundaries[lB]<<endl;
      }
      //Cross-check correct rejection of anything beyond anchor point
      if( fSelection->GetEstimator(iEst)->GetUseAnchor() && ntot != 0 ){
        for( Long_t lB=0; lB<lNDesiredBoundaries-1; lB++) {
          if (lNrawBoundaries[lB+1]>fSelection->GetEstimator(iEst)->GetAnchorPoint()){
            if(lNrawBoundaries[lB]<fSelection->GetEstimator(iEst)->GetAnchorPoint()){
              //This is the threshold, should actually be identical to anchor point please
              lNrawBoundaries[lB] = fSelection->GetEstimator(iEst)->GetAnchorPoint();
            }
          }
        }
      }
      
      cout<<" Done! Saving... "<<endl;
      
      if( lInsane[iEst][iRun] == kFALSE) {
        //Create a sane calibration histogram
        //Should not be the source of excessive memory consumption...
        //...but can be rearranged if needed!
        hCalib[iRun][iEst] = new TH1F(Form("hCalib_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()),"",lNDesiredBoundaries-1,lNrawBoundaries);
        hCalib[iRun][iEst]->SetDirectory(0);
        hCalib[iRun][iEst]->SetBinContent(0,100.5); //Just in case correction functions screw up the values ...
        for(Long_t ibin=1; ibin<hCalib[iRun][iEst]->GetNbinsX()+1; ibin++){
          hCalib[iRun][iEst] -> SetBinContent(ibin, lMiddleOfBins[ibin-1]);
          
          //override in case anchored!
          if( fSelection->GetEstimator(iEst)->GetUseAnchor() ){
            if ( hCalib[iRun][iEst]->GetBinCenter(ibin) < fSelection->GetEstimator(iEst)->GetAnchorPoint() ){
              //Override, this is useless!
              //Alberica's recommendation: outside of user range to be sure!
              hCalib[iRun][iEst] -> SetBinContent(ibin, 100.5);
            }
          }
        }
      }else{
        hCalib[iRun][iEst] = new TH1F(Form("hCalib_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()),"",1,0,1);
        hCalib[iRun][iEst]->SetDirectory(0);
        //There was insufficient information to generate a meaningful calibration for this estimator!
        hCalib[iRun][iEst]->SetBinContent(0,AliMultSelectionCuts::kNoCalib);
        hCalib[iRun][iEst]->SetBinContent(1,AliMultSelectionCuts::kNoCalib);
        hCalib[iRun][iEst]->SetBinContent(2,AliMultSelectionCuts::kNoCalib);
      }
      //==== End Floating Point Calibration Engine ====
    }

    //Write OADB object
    if ( !lAutoDiscover ){
      cout<<"--- Processing run range "<<fFirstRun[iRun]<<"-"<<fLastRun[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<")..."<<endl;
    }else{
      cout<<"--- Processing run "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<fNRunRanges<<")..."<<endl;
    }
    oadbMultSelection = new AliOADBMultSelection();
    cuts              = new AliMultSelectionCuts();
    cuts = fMultSelectionCuts;
    fsels             = new AliMultSelection( fSelection );
    
    oadbMultSelection->SetEventCuts    (cuts );
    oadbMultSelection->SetMultSelection(fsels);
    
    for ( Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
      hCalibDataVtx[iVar] = (TProfile*) hCalibVtx[iRun][iVar]->Clone( Form("hCalibVtx_%s",AliMultInput::VarName[iVar].Data() ) );
      hCalibDataVtx[iVar]->SetDirectory(0);
      if( VarDoAutocalib[iVar] ) oadbMultSelection->AddCalibHistoVtx( hCalibDataVtx[iVar] );
    }
    
    for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      //Average values
      fsels->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );
      // hCalibVtx[iRun][iEst] = (TProfile*)gDirectory->Get(Form("hCalibVtx_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()));
      //Beware similar names! Will be saved ...
      hCalibData[iEst] = (TH1F*) hCalib[iRun][iEst]->Clone( Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
      hCalibData[iEst]->SetDirectory(0);
      oadbMultSelection->AddCalibHisto( hCalibData[iEst] );
    }
    cout<<"=================================================================================="<<endl;
    if ( !lAutoDiscover ){
      cout<<"AliMultSelection Object to be saved for run range ("<<fFirstRun[iRun]<<"-"<<fLastRun[iRun]<<")"<<endl;
    }else{
      cout<<"AliMultSelection Object to be saved for run "<<lRunNumbers[iRun]<<""<<endl;
    }
    fsels->PrintInfo();
    cuts->Print();
    cout<<"=================================================================================="<<endl;
    //Protection against saving a calibration object that has been acquired
    //with insufficient statistics
    if ( lRunStats[iRun] > 1000){
      if ( !lAutoDiscover ) {
        oadbContMS->AppendObject(oadbMultSelection, fFirstRun[iRun], fLastRun[iRun] );
      }else{
        oadbContMS->AppendObject(oadbMultSelection, lRunNumbers[iRun], lRunNumbers[iRun] );
      }
    }
    
    TString lR3ListName = lAutoDiscover?Form("single_%i_", lRunNumbers[iRun]):Form("range_%i_%i_", fFirstRun[iRun], fLastRun[iRun]);
    
    //Run 3 calibration object save
    fRun3->cd();
    
    hR3EventSelection[iRun] = new TH1F(Form("%s_hEventSelection",lR3ListName.Data()), "Event selection requirements", 10,0,10);
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(1,"Phys. Sel.");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(2,"INEL > 0");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(3,"Vtx-Z cut");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(4,"SPD Pileup in Mult bins");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(5,"SPD/Track vtx consistency");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(6,"Tracklets vs clusters");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(7,"Non zero contrib to vtx");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(8,"DEPRECATED");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(9,"Incomplete DAQ");
    hR3EventSelection[iRun]->GetXaxis()->SetBinLabel(10,"DEPRECATED");
    
    hR3EventSelection[iRun]->SetBinContent(1, cuts->GetTriggerCut() );
    hR3EventSelection[iRun]->SetBinContent(2, cuts->GetINELgtZEROCut() );
    hR3EventSelection[iRun]->SetBinContent(3, cuts->GetVzCut() );
    hR3EventSelection[iRun]->SetBinContent(4, cuts->GetRejectPileupInMultBinsCut() );
    hR3EventSelection[iRun]->SetBinContent(5, cuts->GetVertexConsistencyCut() );
    hR3EventSelection[iRun]->SetBinContent(6, cuts->GetTrackletsVsClustersCut() );
    hR3EventSelection[iRun]->SetBinContent(7, cuts->GetNonZeroNContribs() );
    hR3EventSelection[iRun]->SetBinContent(8, cuts->GetIsNotAsymmetricInVZERO() );
    hR3EventSelection[iRun]->SetBinContent(9, cuts->GetIsNotIncompleteDAQ() );
    hR3EventSelection[iRun]->SetBinContent(10, cuts->GetHasGoodVertex2016() );
    
    lR3List -> Add(hR3EventSelection[iRun]);
    
    //Variable calibration
    for ( Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
      hR3CalibDataVtx[iRun][iVar] = (TProfile*) hCalibVtx[iRun][iVar]->Clone( Form("%s_hVtx_%s",lR3ListName.Data(),AliMultInput::VarName[iVar].Data() ) );
      hR3CalibDataVtx[iRun][iVar]->SetDirectory(0);
      if( VarDoAutocalib[iVar] ) lR3List -> Add(hR3CalibDataVtx[iRun][iVar]);
    }
    //Estimator definition and properties
    for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      hR3DefData[iRun][iEst] = new TH1F( Form("%s_hEstimatorDef_%s",lR3ListName.Data(),fSelection->GetEstimator(iEst)->GetName()) , "", 5,0,5 );
      hR3DefData[iRun][iEst]->SetTitle(fSelection->GetEstimator(iEst)->GetDefinition());
      
      hR3DefData[iRun][iEst]->GetXaxis()->SetBinLabel(1,"Is Integer");
      hR3DefData[iRun][iEst]->GetXaxis()->SetBinLabel(2,"Mean value");
      hR3DefData[iRun][iEst]->GetXaxis()->SetBinLabel(3,"Is anchored");
      hR3DefData[iRun][iEst]->GetXaxis()->SetBinLabel(4,"Anchor point");
      hR3DefData[iRun][iEst]->GetXaxis()->SetBinLabel(5,"Anchor percentile");
      
      hR3DefData[iRun][iEst]->SetBinContent(1, fSelection->GetEstimator(iEst)->IsInteger());
      hR3DefData[iRun][iEst]->SetBinContent(2, fSelection->GetEstimator(iEst)->GetMean());
      hR3DefData[iRun][iEst]->SetBinContent(3, fSelection->GetEstimator(iEst)->GetUseAnchor());
      hR3DefData[iRun][iEst]->SetBinContent(4, fSelection->GetEstimator(iEst)->GetAnchorPoint());
      hR3DefData[iRun][iEst]->SetBinContent(5, fSelection->GetEstimator(iEst)->GetAnchorPercentile());
      
      lR3List -> Add(hR3DefData[iRun][iEst]);
    }
    //Estimator calibration
    for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      hR3CalibData[iRun][iEst] = (TH1F*) hCalib[iRun][iEst]->Clone( Form("%s_hMultSelCalib_%s",lR3ListName.Data(),fSelection->GetEstimator(iEst)->GetName()) );
      lR3List -> Add(hR3CalibData[iRun][iEst]);
    }
        
    Bool_t lThisIsReference = kFALSE;
    if(!lAutoDiscover){
      if ( fFirstRun[iRun] <= fRunToUseAsDefault && fRunToUseAsDefault <= fLastRun[iRun]) lThisIsReference = kTRUE;
    }else{
      if ( lRunNumbers[iRun] == fRunToUseAsDefault ) lThisIsReference = kTRUE;
    }
    if( lThisIsReference ){
      //========================================================================
      //DEFAULT OADB Object saving procedure STARTS here
      oadbMultSelection = new AliOADBMultSelection("Default");
      cuts              = new AliMultSelectionCuts();
      cuts = fMultSelectionCuts;
      fsels             = new AliMultSelection    ( fSelection         );
      
      const Int_t lNEstimatorsThis = fSelection->GetNEstimators();
      
      //Default Stuff
      TH1F * hDummy[lNEstimatorsThis];
      TProfile * hDummyVtxZ[AliMultInput::kNVariables];
      for ( Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
        hDummyVtxZ[iVar] = (TProfile*) hCalibVtx[iRun][iVar]->Clone( Form("hCalibVtx_%s",AliMultInput::VarName[iVar].Data() ) );
        hDummyVtxZ[iVar] -> SetDirectory(0);
        if( VarDoAutocalib[iVar] ) oadbMultSelection->AddCalibHistoVtx( hDummyVtxZ[iVar] );
      }
      
      for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        //Get Meaningful means
        fsels->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );
        //Clone last histogram ...
        hDummy[iEst]= (TH1F*) hCalib[iRun][iEst]->Clone( Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
        hDummy[iEst]->SetDirectory(0);
      }
      
      cout<<"=================================================================================="<<endl;
      cout<<" Detected that this particular run / run range is special, will save it as default"<<endl;
      cout<<" AliMultSelection Object to be saved (DEFAULT)"<<endl;
      fsels->PrintInfo();
      cout<<"=================================================================================="<<endl;
      
      oadbMultSelection->SetEventCuts        ( cuts  );
      oadbMultSelection->SetMultSelection    ( fsels );
      for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        oadbMultSelection->AddCalibHisto( hDummy[iEst] );
      }
      oadbContMS->AddDefaultObject(oadbMultSelection);
      //DEFAULT OADB Object saving procedure ENDS here
      //========================================================================
      
      TString lR3ListName = "Default";
      
      //Run 3 calibration object save
      fRun3->cd();
      
      hR3EventSelectionDefault = new TH1F("default_hEventSelection", "Event selection requirements", 10,0,10);
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(1,"Phys. Sel.");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(2,"INEL > 0");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(3,"Vtx-Z cut");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(4,"SPD Pileup in Mult bins");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(5,"SPD/Track vtx consistency");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(6,"Tracklets vs clusters");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(7,"Non zero contrib to vtx");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(8,"DEPRECATED");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(9,"Incomplete DAQ");
      hR3EventSelectionDefault->GetXaxis()->SetBinLabel(10,"DEPRECATED");
      
      hR3EventSelectionDefault->SetBinContent(1, cuts->GetTriggerCut() );
      hR3EventSelectionDefault->SetBinContent(2, cuts->GetINELgtZEROCut() );
      hR3EventSelectionDefault->SetBinContent(3, cuts->GetVzCut() );
      hR3EventSelectionDefault->SetBinContent(4, cuts->GetRejectPileupInMultBinsCut() );
      hR3EventSelectionDefault->SetBinContent(5, cuts->GetVertexConsistencyCut() );
      hR3EventSelectionDefault->SetBinContent(6, cuts->GetTrackletsVsClustersCut() );
      hR3EventSelectionDefault->SetBinContent(7, cuts->GetNonZeroNContribs() );
      hR3EventSelectionDefault->SetBinContent(8, cuts->GetIsNotAsymmetricInVZERO() );
      hR3EventSelectionDefault->SetBinContent(9, cuts->GetIsNotIncompleteDAQ() );
      hR3EventSelectionDefault->SetBinContent(10, cuts->GetHasGoodVertex2016() );
      
      lR3List -> Add(hR3EventSelectionDefault);
      
      //Variable calibration
      for ( Int_t iVar=0; iVar<AliMultInput::kNVariables; iVar++) {
        hR3CalibDataVtxDefault[iVar] = (TProfile*) hCalibVtx[iRun][iVar]->Clone( Form("default_hVtx_%s",AliMultInput::VarName[iVar].Data() ) );
        hR3CalibDataVtxDefault[iVar]->SetDirectory(0);
        if( VarDoAutocalib[iVar] ) lR3List -> Add(hR3CalibDataVtxDefault[iVar]);
      }
      //Estimator definition and properties
      for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        hR3DefDataDefault[iEst] = new TH1F( Form("default_hEstimatorDef_%s",fSelection->GetEstimator(iEst)->GetName()) , "", 5,0,5 );
        hR3DefDataDefault[iEst]->SetTitle(fSelection->GetEstimator(iEst)->GetDefinition());
        
        hR3DefDataDefault[iEst]->GetXaxis()->SetBinLabel(1,"Is Integer");
        hR3DefDataDefault[iEst]->GetXaxis()->SetBinLabel(2,"Mean value");
        hR3DefDataDefault[iEst]->GetXaxis()->SetBinLabel(3,"Is anchored");
        hR3DefDataDefault[iEst]->GetXaxis()->SetBinLabel(4,"Anchor point");
        hR3DefDataDefault[iEst]->GetXaxis()->SetBinLabel(5,"Anchor percentile");
        
        hR3DefDataDefault[iEst]->SetBinContent(1, fSelection->GetEstimator(iEst)->IsInteger());
        hR3DefDataDefault[iEst]->SetBinContent(2, fSelection->GetEstimator(iEst)->GetMean());
        hR3DefDataDefault[iEst]->SetBinContent(3, fSelection->GetEstimator(iEst)->GetUseAnchor());
        hR3DefDataDefault[iEst]->SetBinContent(4, fSelection->GetEstimator(iEst)->GetAnchorPoint());
        hR3DefDataDefault[iEst]->SetBinContent(5, fSelection->GetEstimator(iEst)->GetAnchorPercentile());
        
        lR3List -> Add(hR3DefDataDefault[iEst]);
      }
      //Estimator calibration
      for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
        hR3CalibDataDefault[iEst] = (TH1F*) hCalib[iRun][iEst]->Clone( Form("default_hMultSelCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
        lR3List -> Add(hR3CalibDataDefault[iEst]);
      }
           
      
    }
    //Cleanup
    for(int ii = 0; ii < lNEstimatorsThis; ii++){
      delete[] index[ii];
      delete[] values[ii];
    }
    delete[] index;
    delete[] values;
  }
  
  if( fRunToUseAsDefault < 0 ){
    //========================================================================
    //DEFAULT OADB Object saving procedure STARTS here
    oadbMultSelection = new AliOADBMultSelection("Default");
    cuts              = new AliMultSelectionCuts();
    cuts = fMultSelectionCuts;
    fsels             = new AliMultSelection    ( fSelection         );
    
    const Int_t lNEstimatorsThis = fSelection->GetNEstimators();
    
    cout<<"=================================================================================="<<endl;
    cout<<" AliMultSelection Object to be saved (DEFAULT)"<<endl;
    cout<<" Warning: this corresponds to the last calibrated run!"<<endl;
    fsels->PrintInfo();
    cout<<"=================================================================================="<<endl;
    
    //Default Stuff
    TH1F * hDummy[lNEstimatorsThis];
    TProfile * hDummyVtxZ[lNEstimatorsThis];
    for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
      //Clone last histogram ...
      hDummy[iEst]= (TH1F*) hCalib[fNRunRanges-1][iEst]->Clone( Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
      hDummy[iEst]->SetDirectory(0);
    }
    
    oadbMultSelection->SetEventCuts        ( cuts  );
    oadbMultSelection->SetMultSelection    ( fsels );
    for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++){
      oadbMultSelection->AddCalibHisto( hDummy[iEst] );
      // oadbMultSelection->AddCalibHistoVtx( hDummyVtxZ[iEst] );
    }
    oadbContMS->AddDefaultObject(oadbMultSelection);
    //DEFAULT OADB Object saving procedure ENDS here
    //========================================================================
  }
  
  cout<<"All done, will write OADB..."<<endl;
  fRun3->Write();
  lR3List->Write("multSel",1);
  f->cd();
  oadbContMS->Write();
  cout<<" Done!"<<endl;
  return kTRUE;
}
//________________________________________________________________
Float_t AliMultSelectionCalibrator::MinVal( Float_t A, Float_t B ) {
  if( A < B ) {
    return A;
  }
  else {
    return B;
  }
}
//________________________________________________________________
void AliMultSelectionCalibrator::SetupStandardInput() {
  //============================================================
  // --- Definition of Variables for estimators ---
  //============================================================
  cout<<"===================================="<<endl;
  cout<<"==== Setting up standard input..."<<endl;
  cout<<"===================================="<<endl;
  
  //Create input variables in AliMultInput Class
  AliMultVariable *fVariable[AliMultInput::kNVariables];
  
  //More compact than explicitly writing everything
  for(Int_t ii=0; ii<AliMultInput::kNVariables; ii++){
    fVariable[ii] = new AliMultVariable(AliMultInput::VarName[ii].Data());
    if( AliMultInput::VarIsInteger[ii] ) fVariable[ii] -> SetIsInteger(kTRUE);
    fInput->AddVariable( fVariable[ii] );
  }
}



