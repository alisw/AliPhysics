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
#include "TFile.h"
#include "TStopwatch.h"
#include "TArrayL64.h"

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

    for(Long64_t iEv = 0; iEv<fTree->GetEntries(); iEv++) {

        if ( iEv % 100000 == 0 ) {
            Double_t complete = 100. * ( double ) ( iEv ) / ( double ) ( fTree->GetEntries() );
            cout << "Event # " << iEv << "/" << fTree->GetEntries() << " (" << complete << "%, Time Left: ";
            timer->Stop();
            Double_t time = timer->RealTime();

            //events per hour:
            lEventsPerSecond = ( ( Double_t ) ( iEv ) ) /time;

            timer->Start ( kFALSE );
            Double_t secondsperstep = time / ( Double_t ) ( iEv+1 );
            Double_t secondsleft = ( Double_t ) ( fTree->GetEntries()-iEv-1 ) * secondsperstep;
            Long_t minutesleft = ( Long_t ) ( secondsleft / 60. );
            secondsleft = ( Double_t ) ( ( Long_t ) ( secondsleft ) % 60 );
            cout << minutesleft << "min " << secondsleft << "s, working at "<<lEventsPerSecond<<" Events/s..." << endl;
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
        if( TMath::Abs( lVtxZLocalPointer->GetValue() ) > fMultSelectionCuts->GetVzCut()      ) lSaveThisEvent = kFALSE;
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
                cout<<"(Autodiscover) New Run Found: "<<fRunNumber<<", added as #"<<lNRuns<<" (so far: "<<lNRuns<<" runs)"<<endl;

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
                sTree [ lIndex ] -> Fill();
            }
        }

    }

    //Write buffer to file
    for(Int_t iRun=0; iRun<lNRuns; iRun++) sTree[iRun]->Write();

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

    // STEP 4: Actual determination of boundaries...
    Long64_t *index;
    Double_t *lValues;

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
        sTree[iRun]->SetEstimate(ntot+1);
        //Cast Run Number into drawing conditions
        for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
            lRunStats[iRun] = sTree[iRun]->Draw(fSelection->GetEstimator(iEst)->GetDefinition(),"","goff");
            lValues = sTree[iRun]->GetV1();
            cout<<"--- Calculating averages: "<<flush;
            for( Long64_t iEntry=0; iEntry<ntot; iEntry++) {
                Float_t lThisVal = lValues[iEntry]; //Test
                lAvEst[iEst][iRun] += lThisVal;
                if( lThisVal < lMinEst[iEst][iRun] ) {
                    lMinEst[iEst][iRun] = lThisVal;
                }
                if( lThisVal > lMaxEst[iEst][iRun] ) {
                    lMaxEst[iEst][iRun] = lThisVal;
                }
            }
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

    //Open output OADB file, generate everything within loop
    TFile * f = new TFile (fOutputFileName.Data(), "recreate");
    AliOADBContainer * oadbContMS = new AliOADBContainer("MultSel");

    AliOADBMultSelection * oadbMultSelection = 0x0;
    AliMultSelectionCuts * cuts = 0x0;
    AliMultSelection     * fsels = 0x0;

    //Actual Calibration Histograms
    TH1F * hCalibData[lNEstimators];

    cout<<"(5) Generate Boundaries through a loop in all desired estimators"<<endl;
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
        sTree[iRun]->SetEstimate(ntot+1);
        // Memory allocation: don't repeat it per estimator! only per run
        TArrayL64 index(ntot);
        //Cast Run Number into drawing conditions
        for(Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
            if( ! ( fSelection->GetEstimator(iEst)->IsInteger() ) ) {
                //==== Floating Point Calibration Engine ====
                lRunStats[iRun] = sTree[iRun]->Draw(fSelection->GetEstimator(iEst)->GetDefinition(),"","goff");
                cout<<"--- Sorting estimator "<<fSelection->GetEstimator(iEst)->GetName()<<"..."<<flush;

                TMath::Sort(ntot,sTree[iRun]->GetV1(), index.GetArray() );
                cout<<" Done! Getting Boundaries... "<<flush;

                //Special override in case anchored estimator
                if( fSelection->GetEstimator(iEst)->GetUseAnchor() ){
                    cout<<"Anchoring... "<<flush;
                    //Require determination of index after which values are to be discarded
                    //Count fraction of accepted
                    TString lCondition = fSelection->GetEstimator(iEst)->GetDefinition();
                    lCondition.Append(Form("> %.10f",fSelection->GetEstimator(iEst)->GetAnchorPoint() ) );
                    lAcceptedEvents = sTree[iRun]->Draw(fSelection->GetEstimator(iEst)->GetDefinition(),lCondition.Data(),"goff");
                    lRunStats[iRun] = lAcceptedEvents;
                }
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
                    //cout<<"Position requested: "<<position<<flush;
                    sTree[iRun]->GetEntry( index[position] );
                    //Calculate the estimator with this input, please
                    fSelection->Evaluate ( fInput );
                    //fSelection->PrintInfo();
                    lNrawBoundaries[lB] = fSelection->GetEstimator(iEst)->GetValue();
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
            } else {
                //==== Integer Value Calibration Engine ====
                //Procedure: Create histogram to be filled
                cout<<"Integer calibration engine started!"<<endl;
                const Long_t lNBins    = lMaxEst[iEst][iRun]-lMinEst[iEst][iRun]+1;
                Float_t lLowEdge = lMinEst[iEst][iRun]-0.5;
                Float_t lHighEdge= lMaxEst[iEst][iRun]+0.5;
                cout<<"Inspect: "<<lNBins<<", low "<<lLowEdge<<", high "<<lHighEdge<<endl;
                if( sTree[iRun]->GetEntries() < 1 ) {
                    //Case of an empty run!
                    hCalib[iRun][iEst] = new TH1F(Form("hCalib_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()),"",1,0,1);
                    hCalib[iRun][iEst]->SetDirectory(0);
                } else {
                    TH1F *hTemporary = new TH1F("hTemporary", "", lNBins, lMinEst[iEst][iRun]-0.5, lMaxEst[iEst][iRun]+0.5 );
                    //hTemporary->SetDirectory(0);
                    lRunStats[iRun] = sTree[iRun]->Draw(Form("%s>>hTemporary",fSelection->GetEstimator(iEst)->GetDefinition().Data()),"","goff");
                    cout<<"entries = "<<lRunStats[iRun]<<endl;
                    //In memory now: histogram with content, please normalize to unity
                    hTemporary->Scale(1./((double)(lRunStats[iRun])));

                    Float_t lBoundaries[lNBins+1]; //to store cumulative function
                    lBoundaries[0] = 0;
                    for(Long_t iB=1; iB<hTemporary->GetNbinsX()+1; iB++) {
                        lBoundaries[iB] = lBoundaries[iB-1]+hTemporary->GetBinContent(iB);
                    }
                    //This won't follow what was requested (it cannot, mathematically)
                    hCalib[iRun][iEst] = new TH1F(Form("hCalib_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()),"",lNBins,lLowEdge,lHighEdge);
                    hCalib[iRun][iEst]->SetDirectory(0);
                    for(Long_t ibin=1; ibin<hCalib[iRun][iEst]->GetNbinsX()+1; ibin++) hCalib[iRun][iEst] -> SetBinContent(ibin, 100.0-50.0*(lBoundaries[ibin-1]+lBoundaries[ibin]));
                    //Enough info for calibration determined...
                    delete hTemporary;
                    hTemporary = 0x0;
                }
            }
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
        for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
            //Average values
            fsels->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );

            //Beware similar names! Will be saved ...
            hCalibData[iEst] = (TH1F*) hCalib[iRun][iEst]->Clone( Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
            oadbMultSelection->AddCalibHisto( hCalibData[iEst] );
            hCalibData[iEst]->SetDirectory(0);
        }
        cout<<"=================================================================================="<<endl;
        if ( !lAutoDiscover ){
            cout<<"AliMultSelection Object to be saved for run range "<<fFirstRun[iRun]<<"-"<<fLastRun[iRun]<<")"<<endl;
        }else{
            cout<<"AliMultSelection Object to be saved for run "<<lRunNumbers[iRun]<<")"<<endl;
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
            for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) oadbMultSelection->AddCalibHisto( hDummy[iEst] );
            oadbContMS->AddDefaultObject(oadbMultSelection);
            //DEFAULT OADB Object saving procedure ENDS here
            //========================================================================
        }
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
        for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) {
            //Clone last histogram ...
            hDummy[iEst]= (TH1F*) hCalib[fNRunRanges-1][iEst]->Clone( Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
            hDummy[iEst]->SetDirectory(0);
        }

        oadbMultSelection->SetEventCuts        ( cuts  );
        oadbMultSelection->SetMultSelection    ( fsels );
        for ( Int_t iEst=0; iEst<lNEstimatorsThis; iEst++) oadbMultSelection->AddCalibHisto( hDummy[iEst] );
        oadbContMS->AddDefaultObject(oadbMultSelection);
        //DEFAULT OADB Object saving procedure ENDS here
        //========================================================================
    }

    cout<<"All done, will write OADB..."<<endl;

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

    //Create input variables in AliMultInput Class
    //V0 related
    AliMultVariable *fAmplitude_V0A         = new AliMultVariable("fAmplitude_V0A");
    AliMultVariable *fAmplitude_V0A1        = new AliMultVariable("fAmplitude_V0A1");
    AliMultVariable *fAmplitude_V0A2        = new AliMultVariable("fAmplitude_V0A2");
    AliMultVariable *fAmplitude_V0A3        = new AliMultVariable("fAmplitude_V0A3");
    AliMultVariable *fAmplitude_V0A4        = new AliMultVariable("fAmplitude_V0A4");
    AliMultVariable *fAmplitude_V0C         = new AliMultVariable("fAmplitude_V0C");
    AliMultVariable *fAmplitude_V0C1        = new AliMultVariable("fAmplitude_V0C1");
    AliMultVariable *fAmplitude_V0C2        = new AliMultVariable("fAmplitude_V0C2");
    AliMultVariable *fAmplitude_V0C3        = new AliMultVariable("fAmplitude_V0C3");
    AliMultVariable *fAmplitude_V0C4        = new AliMultVariable("fAmplitude_V0C4");
    AliMultVariable *fAmplitude_V0Apartial = new AliMultVariable("fAmplitude_V0Apartial");
    AliMultVariable *fAmplitude_V0Cpartial = new AliMultVariable("fAmplitude_V0Cpartial");
    AliMultVariable *fAmplitude_V0AEq      = new AliMultVariable("fAmplitude_V0AEq");
    AliMultVariable *fAmplitude_V0CEq      = new AliMultVariable("fAmplitude_V0CEq");
    AliMultVariable *fAmplitude_OnlineV0A  = new AliMultVariable("fAmplitude_OnlineV0A");
    AliMultVariable *fAmplitude_OnlineV0C  = new AliMultVariable("fAmplitude_OnlineV0C");
    //SPD Related
    AliMultVariable *fnSPDClusters         = new AliMultVariable("fnSPDClusters");
    AliMultVariable *fnSPDClusters0        = new AliMultVariable("fnSPDClusters0");
    AliMultVariable *fnSPDClusters1        = new AliMultVariable("fnSPDClusters1");
    fnSPDClusters->SetIsInteger( kTRUE );
    fnSPDClusters0->SetIsInteger( kTRUE );
    fnSPDClusters1->SetIsInteger( kTRUE );
    //AD Related
    AliMultVariable *fMultiplicity_ADA     = new AliMultVariable("fMultiplicity_ADA");
    AliMultVariable *fMultiplicity_ADC     = new AliMultVariable("fMultiplicity_ADC");

    AliMultVariable *fRefMultEta5     = new AliMultVariable("fRefMultEta5");
    fRefMultEta5->SetIsInteger( kTRUE );
    AliMultVariable *fRefMultEta8     = new AliMultVariable("fRefMultEta8");
    fRefMultEta8->SetIsInteger( kTRUE );
    AliMultVariable *fnTracklets     = new AliMultVariable("fnTracklets");
    fnTracklets->SetIsInteger( kTRUE );
    AliMultVariable *fnTracklets08     = new AliMultVariable("fnTracklets08");
    fnTracklets08->SetIsInteger( kTRUE );
    AliMultVariable *fnTracklets15     = new AliMultVariable("fnTracklets15");
    fnTracklets15->SetIsInteger( kTRUE );

    //ZDC Related
    AliMultVariable *fZncEnergy = new AliMultVariable("fZncEnergy");
    AliMultVariable *fZpcEnergy = new AliMultVariable("fZpcEnergy");
    AliMultVariable *fZnaEnergy = new AliMultVariable("fZnaEnergy");
    AliMultVariable *fZpaEnergy = new AliMultVariable("fZpaEnergy");
    AliMultVariable *fZem1Energy = new AliMultVariable("fZem1Energy");
    AliMultVariable *fZem2Energy = new AliMultVariable("fZem2Energy");

    AliMultVariable *fZnaTower = new AliMultVariable("fZnaTower");
    AliMultVariable *fZncTower = new AliMultVariable("fZncTower");
    AliMultVariable *fZpaTower = new AliMultVariable("fZpaTower");
    AliMultVariable *fZpcTower = new AliMultVariable("fZpcTower");

    //Fired or not booleans (stored as integer for compatibility)
    AliMultVariable *fZnaFired = new AliMultVariable("fZnaFired");
    fZnaFired->SetIsInteger(kTRUE);
    AliMultVariable *fZncFired = new AliMultVariable("fZncFired");
    fZncFired->SetIsInteger(kTRUE);
    AliMultVariable *fZpaFired = new AliMultVariable("fZpaFired");
    fZpaFired->SetIsInteger(kTRUE);
    AliMultVariable *fZpcFired = new AliMultVariable("fZpcFired");
    fZpcFired->SetIsInteger(kTRUE);

    //Track counters (now useable as AliMultVariables as well)
    AliMultVariable *fNTracks =                  new AliMultVariable("fNTracks");
    fNTracks->SetIsInteger(kTRUE);
    AliMultVariable *fNTracksGlobal2015 =        new AliMultVariable("fNTracksGlobal2015");
    fNTracksGlobal2015->SetIsInteger(kTRUE);
    AliMultVariable *fNTracksGlobal2015Trigger = new AliMultVariable("fNTracksGlobal2015Trigger");
    fNTracksGlobal2015Trigger->SetIsInteger(kTRUE);
    AliMultVariable *fNTracksITSsa2010 =         new AliMultVariable("fNTracksITSsa2010");
    fNTracksITSsa2010->SetIsInteger(kTRUE);

    AliMultVariable *fNTracksINELgtONE =       new AliMultVariable("fNTracksINELgtONE");
    AliMultVariable *fNPartINELgtONE   =       new AliMultVariable("fNPartINELgtONE");

    //vertex-Z
    AliMultVariable *fEvSel_VtxZ = new AliMultVariable("fEvSel_VtxZ");

    AliMultVariable *fMC_NPart =         new AliMultVariable("fMC_NPart");
    fMC_NPart->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NColl =         new AliMultVariable("fMC_NColl");
    fMC_NColl->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchV0A =         new AliMultVariable("fMC_NchV0A");
    fMC_NchV0A->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchV0C =         new AliMultVariable("fMC_NchV0C");
    fMC_NchV0C->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchEta05 =         new AliMultVariable("fMC_NchEta05");
    fMC_NchEta05->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchEta08 =         new AliMultVariable("fMC_NchEta08");
    fMC_NchEta08->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchEta10 =         new AliMultVariable("fMC_NchEta10");
    fMC_NchEta10->SetIsInteger(kTRUE);
    AliMultVariable *fMC_NchEta14 =         new AliMultVariable("fMC_NchEta14");
    fMC_NchEta14->SetIsInteger(kTRUE);

    //Add to AliMultInput Object
    fInput->AddVariable( fAmplitude_V0A );
    fInput->AddVariable( fAmplitude_V0A1 );
    fInput->AddVariable( fAmplitude_V0A2 );
    fInput->AddVariable( fAmplitude_V0A3 );
    fInput->AddVariable( fAmplitude_V0A4 );
    fInput->AddVariable( fAmplitude_V0C );
    fInput->AddVariable( fAmplitude_V0C1 );
    fInput->AddVariable( fAmplitude_V0C2 );
    fInput->AddVariable( fAmplitude_V0C3 );
    fInput->AddVariable( fAmplitude_V0C4 );
    fInput->AddVariable( fAmplitude_V0Apartial );
    fInput->AddVariable( fAmplitude_V0Cpartial );
    fInput->AddVariable( fAmplitude_V0AEq );
    fInput->AddVariable( fAmplitude_V0CEq );
    fInput->AddVariable( fAmplitude_OnlineV0A );
    fInput->AddVariable( fAmplitude_OnlineV0C );
    fInput->AddVariable( fMultiplicity_ADA );
    fInput->AddVariable( fMultiplicity_ADC );
    fInput->AddVariable( fnSPDClusters );
    fInput->AddVariable( fnSPDClusters0 );
    fInput->AddVariable( fnSPDClusters1 );
    fInput->AddVariable( fnTracklets   );
    fInput->AddVariable( fnTracklets08   );
    fInput->AddVariable( fnTracklets15   );
    fInput->AddVariable( fRefMultEta5  );
    fInput->AddVariable( fRefMultEta8  );
    fInput->AddVariable( fZncEnergy );
    fInput->AddVariable( fZpcEnergy );
    fInput->AddVariable( fZnaEnergy );
    fInput->AddVariable( fZpaEnergy );
    fInput->AddVariable( fZem1Energy );
    fInput->AddVariable( fZem2Energy );
    fInput->AddVariable( fZnaTower );
    fInput->AddVariable( fZncTower );
    fInput->AddVariable( fZpaTower );
    fInput->AddVariable( fZpcTower );
    fInput->AddVariable( fZnaFired );
    fInput->AddVariable( fZncFired );
    fInput->AddVariable( fZpaFired );
    fInput->AddVariable( fZpcFired );
    fInput->AddVariable( fNTracks                  );
    fInput->AddVariable( fNTracksGlobal2015        );
    fInput->AddVariable( fNTracksGlobal2015Trigger );
    fInput->AddVariable( fNTracksITSsa2010         );
    fInput->AddVariable( fNTracksINELgtONE );
    fInput->AddVariable( fNPartINELgtONE         );
    fInput->AddVariable( fEvSel_VtxZ  );

    fInput->AddVariable( fMC_NPart );
    fInput->AddVariable( fMC_NColl );
    fInput->AddVariable( fMC_NchV0A );
    fInput->AddVariable( fMC_NchV0C );
    fInput->AddVariable( fMC_NchEta05 );
    fInput->AddVariable( fMC_NchEta08 );
    fInput->AddVariable( fMC_NchEta10 );
    fInput->AddVariable( fMC_NchEta14 );
    //============================================================

}
