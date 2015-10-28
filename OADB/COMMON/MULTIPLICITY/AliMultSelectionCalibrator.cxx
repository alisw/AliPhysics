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
#include "AliESDEvent.h"
#include "TList.h"
#include "TFile.h"
#include "TStopwatch.h"

ClassImp(AliMultSelectionCalibrator);

AliMultSelectionCalibrator::AliMultSelectionCalibrator() :
    TNamed(), fInputFileName(""), fBufferFileName("buffer.root"),
    fOutputFileName(""), fInput(0), fSelection(0), fMultSelectionCuts(0), fCalibHists(0),
    lNDesiredBoundaries(0), lDesiredBoundaries(0)
{
    // Constructor

    // Create Event Selector
    fMultSelectionCuts = new AliMultSelectionCuts();
    fMultSelectionCuts -> Print();

    //Basic I/O for MultSelection framework
    fInput     = new AliMultInput();
    fSelection = new AliMultSelection();

    //Make sure the TList owns its objects
    fCalibHists = new TList();
    fCalibHists -> SetOwner(kTRUE);
}

AliMultSelectionCalibrator::AliMultSelectionCalibrator(const char * name, const char * title):
    TNamed(name,title), fInputFileName(""), fBufferFileName("buffer.root"),
    fOutputFileName(""), fInput(0), fSelection(0), fMultSelectionCuts(0), fCalibHists(0),
    lNDesiredBoundaries(0), lDesiredBoundaries(0)
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

    //Basic I/O for MultSelection framework
    fInput     = new AliMultInput();
    fSelection = new AliMultSelection();

    //Make sure the TList owns its objects
    fCalibHists = new TList();
    fCalibHists -> SetOwner(kTRUE);

}
AliMultSelectionCalibrator::~AliMultSelectionCalibrator() {
    // Destructor

    if ( fMultSelectionCuts ) {
        delete fMultSelectionCuts;
        fMultSelectionCuts = 0x0;
    }

    //Make sure the TList owns its objects
    fCalibHists -> SetOwner(kTRUE);
}
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
    //     (3b) Compute Averages: <V0A>, <V0C>, <V0Apartial>, <V0Cpartial>
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
    Float_t fEvSel_VtxZ                      = 10.0 ;
    Int_t fRunNumber;

    //SetBranchAddresses for event Selection Variables
    //(multiplicity related will be done automatically!)
    fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTree->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTree->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTree->SetBranchAddress("fEvSel_VtxZ",&fEvSel_VtxZ);
    fTree->SetBranchAddress("fRunNumber",&fRunNumber);

    //============================================================
    // --- Definition of Variables for estimators ---
    //============================================================
    // -> only this part needs changing for any additional
    //    variables that may be required for estimators
    //============================================================

    //Create input variables in AliMultInput Class
    //V0 related
    AliMultVariable *fAmplitude_V0A        = new AliMultVariable("fAmplitude_V0A");
    AliMultVariable *fAmplitude_V0C        = new AliMultVariable("fAmplitude_V0C");
    AliMultVariable *fAmplitude_V0Apartial = new AliMultVariable("fAmplitude_V0Apartial");
    AliMultVariable *fAmplitude_V0Cpartial = new AliMultVariable("fAmplitude_V0Cpartial");
    AliMultVariable *fAmplitude_V0AEq      = new AliMultVariable("fAmplitude_V0AEq");
    AliMultVariable *fAmplitude_V0CEq      = new AliMultVariable("fAmplitude_V0CEq");
    AliMultVariable *fAmplitude_OnlineV0A  = new AliMultVariable("fAmplitude_OnlineV0A");
    AliMultVariable *fAmplitude_OnlineV0C  = new AliMultVariable("fAmplitude_OnlineV0C");
    //SPD Related
    AliMultVariable *fnSPDClusters         = new AliMultVariable("fnSPDClusters");
    fnSPDClusters->SetIsInteger( kTRUE );
    //AD Related
    AliMultVariable *fMultiplicity_ADA     = new AliMultVariable("fMultiplicity_ADA");
    AliMultVariable *fMultiplicity_ADC     = new AliMultVariable("fMultiplicity_ADC");

    AliMultVariable *fRefMultEta5     = new AliMultVariable("fRefMultEta5");
    fRefMultEta5->SetIsInteger( kTRUE );
    AliMultVariable *fRefMultEta8     = new AliMultVariable("fRefMultEta8");
    fRefMultEta8->SetIsInteger( kTRUE );
    AliMultVariable *fnTracklets     = new AliMultVariable("fnTracklets");
    fnTracklets->SetIsInteger( kTRUE );
    
    //Add to AliMultInput Object
    fInput->AddVariable( fAmplitude_V0A );
    fInput->AddVariable( fAmplitude_V0C );
    fInput->AddVariable( fAmplitude_V0Apartial );
    fInput->AddVariable( fAmplitude_V0Cpartial );
    fInput->AddVariable( fAmplitude_V0AEq );
    fInput->AddVariable( fAmplitude_V0CEq );
    fInput->AddVariable( fAmplitude_OnlineV0A );
    fInput->AddVariable( fAmplitude_OnlineV0C );
    fInput->AddVariable( fMultiplicity_ADA );
    fInput->AddVariable( fMultiplicity_ADC );
    fInput->AddVariable( fnSPDClusters );
    fInput->AddVariable( fnTracklets   );
    fInput->AddVariable( fRefMultEta5  );
    fInput->AddVariable( fRefMultEta8  );

    //============================================================

    //Binding to input variables
    for(Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++) {
        if( !fInput->GetVariable(iVar)->IsInteger() ) {
            fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValue());
        } else {
            fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValueInteger());
        }
    }

    //============================================================
    // --- Definition of Estimators ---
    //============================================================
    // -> only this part needs changing for any additional
    //    estimators that use known variables
    //============================================================

    AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    AliMultEstimator *fEstV0A = new AliMultEstimator("V0A", "", "(fAmplitude_V0A)");
    AliMultEstimator *fEstV0C = new AliMultEstimator("V0C", "", "(fAmplitude_V0C)");

    AliMultEstimator *fEstOnlineV0M = new AliMultEstimator("OnlineV0M", "", "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)");
    AliMultEstimator *fEstOnlineV0A = new AliMultEstimator("OnlineV0A", "", "(fAmplitude_OnlineV0A)");
    AliMultEstimator *fEstOnlineV0C = new AliMultEstimator("OnlineV0C", "", "(fAmplitude_OnlineV0C)");

    AliMultEstimator *fEstADM = new AliMultEstimator("ADM", "", "(fMultiplicity_ADA)+(fMultiplicity_ADC)");
    AliMultEstimator *fEstADA = new AliMultEstimator("ADA", "", "(fMultiplicity_ADA)");
    AliMultEstimator *fEstADC = new AliMultEstimator("ADC", "", "(fMultiplicity_ADC)");

    //Integer estimators
    AliMultEstimator *fEstnSPDClusters = new AliMultEstimator("SPDClusters", "", "(fnSPDClusters)");
    fEstnSPDClusters->SetIsInteger(kTRUE);
    AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
    fEstnSPDTracklets->SetIsInteger(kTRUE);
    AliMultEstimator *fEstRefMultEta5 = new AliMultEstimator("RefMult05", "", "(fRefMultEta5)");
    fEstRefMultEta5->SetIsInteger(kTRUE);
    AliMultEstimator *fEstRefMultEta8 = new AliMultEstimator("RefMult08", "", "(fRefMultEta8)");
    fEstRefMultEta8->SetIsInteger(kTRUE);
    
    fSelection -> AddEstimator( fEstV0M );
    fSelection -> AddEstimator( fEstV0A );
    fSelection -> AddEstimator( fEstV0C );
    fSelection -> AddEstimator( fEstOnlineV0M );
    fSelection -> AddEstimator( fEstOnlineV0A );
    fSelection -> AddEstimator( fEstOnlineV0C );
    fSelection -> AddEstimator( fEstADM );
    fSelection -> AddEstimator( fEstADA );
    fSelection -> AddEstimator( fEstADC );
    fSelection -> AddEstimator( fEstnSPDClusters  );
    fSelection -> AddEstimator( fEstnSPDTracklets );
    fSelection -> AddEstimator( fEstRefMultEta5 );
    fSelection -> AddEstimator( fEstRefMultEta8 );
    
    //============================================================

    Long64_t lNEv = fTree->GetEntries();
    cout<<"(1) File opened, event count is "<<lNEv<<endl;
    
    cout<<"(2) Creating buffer, computing averages"<<endl;
    const int lMax = 1000;
    const int lMaxQuantiles = 10000;
    Int_t lRunNumbers[lMaxQuantiles];
    Long_t lRunStats[lMaxQuantiles];
    Int_t lNRuns = 0;
    Bool_t lNewRun = kTRUE;
    Int_t lThisRunIndex = -1;
    //Buffer file with run-by-run TTree objects needed for later processing

    TFile *fOutput = new TFile (fBufferFileName.Data(), "RECREATE");
    TTree *sTree[lMaxQuantiles];
    cout<<"Creating Trees..."<<endl;
    for(Int_t iRun=0; iRun<lMax; iRun++) {
        sTree[iRun] = new TTree(Form("sTree%i",iRun),Form("sTree%i",iRun));
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

    const int lNEstimators = fSelection->GetNEstimators();
    //For computing average values of estimators
    Double_t lAvEst[lNEstimators][lMax];
    //For computing extreme values (useful for integer calibration mode)
    Double_t lMaxEst[lNEstimators][lMax];
    Double_t lMinEst[lNEstimators][lMax];

    for(Long_t iEst=0; iEst<lNEstimators; iEst++) {
        for(Long_t iRun=0; iRun<lMax; iRun++) lAvEst[iEst][iRun] = 0;
        for(Long_t iRun=0; iRun<lMax; iRun++) lMaxEst[iEst][iRun] = -1e+3;
        for(Long_t iRun=0; iRun<lMax; iRun++) lMinEst[iEst][iRun] = 1e+6; //not more than a million, I hope?
    }

    //Add Timer
    TStopwatch* timer = new TStopwatch();
    timer->Start ( kTRUE );

    //Compute events-per-hour performance metric
    Double_t lEventsPerSecond = 0;

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

        lNewRun = kTRUE;
        fTree->GetEntry(iEv); //Look at next event
        //May not be the smartest procedure but will work
        for(Int_t iRun=0; iRun<lNRuns; iRun++) {
            if( lRunNumbers[iRun] == fRunNumber ) {
                lNewRun = kFALSE;
                lThisRunIndex = iRun;
            }
        }
        if( lNewRun == kTRUE ) {
            cout<<"(2) New Run Found: "<<fRunNumber<<", added as #"<<lNRuns<<" (so far: "<<lNRuns<<" runs)"<<endl;
            lRunNumbers[lNRuns] = fRunNumber;
            lThisRunIndex = lNRuns;
            lNRuns++;
        }

        //Perform Event selection
        Bool_t lSaveThisEvent = kTRUE; //let's be optimistic

        //Check Selections as they are in the fMultSelectionCuts Object
        if( fMultSelectionCuts->GetTriggerCut()    && ! fEvSel_Triggered  ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO ) lSaveThisEvent = kFALSE;
        if( TMath::Abs(fEvSel_VtxZ) > fMultSelectionCuts->GetVzCut()      ) lSaveThisEvent = kFALSE;
        //ADD ME HERE: Tracklets Vs Clusters Cut?
        if( fMultSelectionCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins    ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster  ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices) lSaveThisEvent = kFALSE;
        
        if ( lSaveThisEvent ) {
            sTree [lThisRunIndex] -> Fill();
        }
        if(lNRuns>lMax) {
            lNRuns = lMax;
            AliWarningF("Exceeded maximum allowed number of runs to quantile! (Nruns now = %i)",lNRuns );
            AliWarningF("Will continue using only %i runs.", lMax);
            break;
        }
    }
/*
    cout<<"Inspect average estimator values: "<<endl;
    for(Long_t iRun=0; iRun<lNRuns; iRun++) {
        for(Long_t iEst=0; iEst<lNEstimators; iEst++) {
            cout<<"iRun: "<<iRun<<", iEst: "<<iEst<<", average = "<<lAvEst[iEst][iRun]/((Double_t)sTree[iRun]->GetEntries())<<endl;
            if( sTree[iRun]->GetEntries() > 0 ) {
                lAvEst[iEst][iRun] =  lAvEst[iEst][iRun] / ((Double_t)sTree[iRun]->GetEntries());
            } else {
                lAvEst[iEst][iRun] = -1;
            }
        }
    }
*/
    //Write buffer to file
    for(Int_t iRun=0; iRun<lNRuns; iRun++) sTree[iRun]->Write();

    cout<<"(3) Inspect List of Runs and their corresponding statistics passing cuts: "<<endl;
    for(Int_t iRun = 0; iRun<lNRuns; iRun++) {
        cout<<" --- "<<lRunNumbers[iRun]<<", N(events) = "<<sTree[iRun]->GetEntries()<<endl;
    }
    cout<<endl;

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
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        const Long64_t ntot = (Long64_t) sTree[iRun]->GetEntries();
        cout<<"--- Processing run number "<<lRunNumbers[iRun]<<"/"<<lNRuns<<", with "<<ntot<<" events..."<<endl;
        sTree[iRun]->SetEstimate(ntot+1);
        //Cast Run Number into drawing conditions
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
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
        }
    }
    
    cout<<"(5) Generate Boundaries through a loop in all desired estimators"<<endl;
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        const Long64_t ntot = (Long64_t) sTree[iRun]->GetEntries();
        cout<<"--- Processing run number "<<lRunNumbers[iRun]<<"/"<<lNRuns<<", with "<<ntot<<" events..."<<endl;
        sTree[iRun]->SetEstimate(ntot+1);
        //Cast Run Number into drawing conditions
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            if( ! ( fSelection->GetEstimator(iEst)->IsInteger() ) ) {
                //==== Floating Point Calibration Engine ====
                lRunStats[iRun] = sTree[iRun]->Draw(fSelection->GetEstimator(iEst)->GetDefinition(),"","goff");
                //This will make sure we use only a projection of the TTree!
                index = new Long64_t[ntot];

                cout<<"--- Sorting estimator "<<fSelection->GetEstimator(iEst)->GetName()<<"..."<<flush;
                TMath::Sort(ntot,sTree[iRun]->GetV1(),index);
                cout<<" Done! Getting Boundaries..."<<endl;
                lNrawBoundaries[0] = 0.0;
                for( Long_t lB=1; lB<lNDesiredBoundaries; lB++) {
                    Long64_t position = (Long64_t) ( 0.01 * ((Double_t)(ntot)* lDesiredBoundaries[lB] ) );
                    //cout<<"Position requested: "<<position<<flush;
                    sTree[iRun]->GetEntry( index[position] );
                    //Calculate the estimator with this input, please
                    fSelection->Evaluate ( fInput );
                    lNrawBoundaries[lB] = fSelection->GetEstimator(iEst)->GetValue();
                }
                //Should not be the source of excessive memory consumption...
                //...but can be rearranged if needed!
                hCalib[iRun][iEst] = new TH1F(Form("hCalib_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName()),"",lNDesiredBoundaries-1,lNrawBoundaries);
                hCalib[iRun][iEst]->SetDirectory(0);
                for(Long_t ibin=1; ibin<hCalib[iRun][iEst]->GetNbinsX()+1; ibin++) hCalib[iRun][iEst] -> SetBinContent(ibin, lMiddleOfBins[ibin-1]);
                //Cleanup: Delete index variable
                delete[] index;
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
                    lRunStats[iRun] = sTree[iRun]->Draw(Form("%s>>hTemporary",fSelection->GetEstimator(iEst)->GetDefinition().Data()),"","goff");

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
    }

    cout<<"(6) Write OADB"<<endl;

    cout<<"Inspect current fSelection"<<endl;
    fSelection->PrintInfo();

    TFile * f = new TFile (fOutputFileName.Data(), "recreate");
    AliOADBContainer * oadbContMS = new AliOADBContainer("MultSel");
    AliOADBMultSelection * oadbMultSelection = new AliOADBMultSelection("Default");
    AliMultSelectionCuts * cuts              = new AliMultSelectionCuts;
    AliMultSelection     * fsels             = new AliMultSelection ( fSelection );

    cout<<"=================================================================================="<<endl; 
    cout<<"AliMultSelection Object to be saved (DEFAULT)"<<endl;
    fsels->PrintInfo();
    cout<<"=================================================================================="<<endl; 

    //Default Stuff
    TH1F * hDummy[lNEstimators];
    for ( Int_t iEst=0; iEst<lNEstimators; iEst++) {
        hDummy[iEst]= new TH1F (Form("hCalib_000000_%s",fSelection->GetEstimator(iEst)->GetName()), "hdummy", 100, 0, 10);
        hDummy[iEst]->SetDirectory(0);
    }

    oadbMultSelection->SetEventCuts        ( cuts  );
    oadbMultSelection->SetMultSelection    ( fsels );
    for ( Int_t iEst=0; iEst<lNEstimators; iEst++) oadbMultSelection->AddCalibHisto( hDummy[iEst] );
    oadbContMS->AddDefaultObject(oadbMultSelection);

    //Actual Calibration Histograms
    TH1F * hCalibData[lNEstimators];

    //Loop over existing runs and write objects as needed
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        cout<<"Processing run number "<<lRunNumbers[iRun]<<endl;
        oadbMultSelection = new AliOADBMultSelection();
        cuts              = new AliMultSelectionCuts;
        fsels             = new AliMultSelection( fSelection );
        //cout<<"Dump"<<endl;
        //fsels->PrintInfo();
        //TODO FIXME: Copy configurations
        //cuts = fMultSelectionCuts;

        oadbMultSelection->SetEventCuts    (cuts );
        oadbMultSelection->SetMultSelection(fsels);
        for ( Int_t iEst=0; iEst<lNEstimators; iEst++) {
            //Average values
            fsels->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );

            hCalibData[iEst] = (TH1F*) hCalib[iRun][iEst]->Clone(Form("hCalib_%i_%s",lRunNumbers[iRun], fSelection->GetEstimator(iEst)->GetName()) );
            oadbMultSelection->AddCalibHisto( hCalibData[iEst]);
            hCalibData[iEst]->SetDirectory(0);
        }
        cout<<"=================================================================================="<<endl; 
        cout<<"AliMultSelection Object to be saved for run "<<lRunNumbers[iRun]<<": "<<endl;
        fsels->PrintInfo();
	cout<<"=================================================================================="<<endl; 
        oadbContMS->AppendObject(oadbMultSelection, lRunNumbers[iRun] ,lRunNumbers[iRun] );
    }
    cout<<"Write OADB..."<<endl;
    //pre-write dump
    //fsels->PrintInfo();

    oadbContMS->Write();
    cout<<" Done!"<<endl;
    return kTRUE;
}

Float_t AliMultSelectionCalibrator::MinVal( Float_t A, Float_t B ) {
    if( A < B ) {
        return A;
    }
    else {
        return B;
    }
}
