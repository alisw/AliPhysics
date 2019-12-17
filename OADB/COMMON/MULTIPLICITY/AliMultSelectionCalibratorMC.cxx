/**********************************************
 *
 *   Class meant to perform calibration and
 *   write an OADB file containing histos and
 *   Event selection criteria used
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
#include "AliMultSelectionCalibratorMC.h"
#include "AliESDEvent.h"
#include "TList.h"
#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"
#include "TStopwatch.h"

ClassImp(AliMultSelectionCalibratorMC);

AliMultSelectionCalibratorMC::AliMultSelectionCalibratorMC() :
    TNamed(), fInputFileNameData(""), fInputFileNameOADB(""), fInputFileNameMC(""),
    fBufferFileNameData("buffer.root"  ),
    fBufferFileNameMC  ("bufferMC.root"),
    fDebugFileName  ("debug.root"),
    fOutputFileName(""), fInput(0), fSelection(0), fMultSelectionCuts(0), fCalibHists(0),
    lNDesiredBoundaries(0), lDesiredBoundaries(0), fRunToUseAsDefault(-1),
    fkUseQuadraticMapping(kFALSE), fNV0MCutoffs(0), fV0MCutoffs()
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

AliMultSelectionCalibratorMC::AliMultSelectionCalibratorMC(const char * name, const char * title):
    TNamed(), fInputFileNameData(""), fInputFileNameOADB(""), fInputFileNameMC(""),
    fBufferFileNameData("buffer.root"  ),
    fBufferFileNameMC  ("bufferMC.root"),
    fDebugFileName  ("debug.root"),
    fOutputFileName(""), fInput(0), fSelection(0), fMultSelectionCuts(0), fCalibHists(0),
    lNDesiredBoundaries(0), lDesiredBoundaries(0), fRunToUseAsDefault(-1),
    fkUseQuadraticMapping(kFALSE), fNV0MCutoffs(0), fV0MCutoffs()
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

    //Basic I/O for MultSelection framework
    fInput     = new AliMultInput();
    fSelection = new AliMultSelection();

    //Make sure the TList owns its objects
    fCalibHists = new TList();
    fCalibHists -> SetOwner(kTRUE);

}
AliMultSelectionCalibratorMC::~AliMultSelectionCalibratorMC() {
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
Bool_t AliMultSelectionCalibratorMC::Calibrate() {
    // Function meant to generate calibration OADB for MC using scaling

    cout<<"=== STARTING MC CALIBRATION PROCEDURE ==="<<endl;
    cout<<" * Input File, Data.......: "<<fInputFileNameData.Data()<<endl;
    cout<<" * Input File, Data OADB..: "<<fInputFileNameOADB.Data()<<endl;
    cout<<" * Input File, MC.........: "<<fInputFileNameMC.Data()<<endl;
    cout<<" * Output File, MC OADB...: "<<fOutputFileName.Data()<<endl;
    cout<<endl;

    // STEP 1: Basic I/O
    cout<<"(2) Opening OADB file for data..."<<endl;
    TFile * fOADB = new TFile (fInputFileNameOADB.Data());

    //Acquire OADB Container
    AliOADBContainer * oadbContMS = (AliOADBContainer*) fOADB->Get("MultSel");

    //Pointer to selection parameters
    AliOADBMultSelection * oadbMultSelection = 0x0;

    //= (AliOADBMultSelection* )oadbContMS->GetObject(run, "Default");

    cout<<"(1) Opening Data and MC TTree files..."<<endl;

    //Open File
    TFile *fInputFileData = TFile::Open( fInputFileNameData.Data(), "READ");
    if(!fInputFileData) {
        AliWarningF("File %s not found!", fInputFileNameData.Data() );
        return kFALSE;
    }
    //Locate TTree object
    TTree* fTree = (TTree*)fInputFileData->FindObjectAny("fTreeEvent");
    if(!fTree) {
        AliWarning("fTreeEvent object not found!" );
        return kFALSE;
    }
    fTree->SetName("fTreeEventData");

    //Open File
    TFile *fInputFileMC = TFile::Open( fInputFileNameMC.Data(), "READ");
    if(!fInputFileMC) {
        AliWarningF("File %s not found!", fInputFileNameMC.Data() );
        return kFALSE;
    }
    //Locate TTree object
    TTree* fTreeMC = (TTree*)fInputFileMC->FindObjectAny("fTreeEvent");
    if(!fTreeMC) {
        AliWarning("fTreeEvent object not found!" );
        return kFALSE;
    }
    fTreeMC->SetName("fTreeEventMC");

    //Event Selection Variables
    Bool_t fEvSel_IsNotPileupInMultBins      = kFALSE ;
    Bool_t fEvSel_Triggered                  = kFALSE ;
    Bool_t fEvSel_INELgtZERO                 = kFALSE ;
    Bool_t fEvSel_PassesTrackletVsCluster    = kFALSE ;
    Bool_t fEvSel_HasNoInconsistentVertices  = kFALSE ;
    //FIXME/CAUTION: non-zero if using tree without that branch
    Int_t fnContributors = 1000;
    Int_t fRunNumber;

    //SetBranchAddresses for event Selection Variables
    fTree->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTree->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTree->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTree->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTree->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTree->SetBranchAddress("fRunNumber",&fRunNumber);
    fTree->SetBranchAddress("fnContributors", &fnContributors);

    fTreeMC->SetBranchAddress("fEvSel_IsNotPileupInMultBins",&fEvSel_IsNotPileupInMultBins);
    fTreeMC->SetBranchAddress("fEvSel_PassesTrackletVsCluster",&fEvSel_PassesTrackletVsCluster);
    fTreeMC->SetBranchAddress("fEvSel_HasNoInconsistentVertices",&fEvSel_HasNoInconsistentVertices);
    fTreeMC->SetBranchAddress("fEvSel_Triggered",&fEvSel_Triggered);
    fTreeMC->SetBranchAddress("fEvSel_INELgtZERO",&fEvSel_INELgtZERO);
    fTreeMC->SetBranchAddress("fRunNumber",&fRunNumber);
    fTreeMC->SetBranchAddress("fnContributors", &fnContributors);

    //============================================================
    // Auto-configure Input
    //============================================================

    if ( fInput->GetNVariables() < 1 ) {
        cout<<"Error: No Input Variables configured!"<<endl;
        cout<<"The simplest way to get rid of this problem is to remember to call SetupStandardInput()!"<<endl;
        return kFALSE; //failure to calibrate
    }
    //Data bindings
    for(Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++) {
        if( !fInput->GetVariable(iVar)->IsInteger() ) {
            fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValue());
        } else {
            fTree->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValueInteger());
        }
    }

    //MC bindings
    for(Long_t iVar=0; iVar<fInput->GetNVariables(); iVar++) {
        if( !fInput->GetVariable(iVar)->IsInteger() ) {
            fTreeMC->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValue());
        } else {
            fTreeMC->SetBranchAddress(fInput->GetVariable(iVar)->GetName(),&fInput->GetVariable(iVar)->GetRValueInteger());
        }
    }

    //WARNING/FIXME- This will not allow for event selections to change in the middle of the processed dataset
    //This will shortly be replaced with a proper selection
    //in which the OADB is queried and a run range mapping object is created
    fTree->GetEntry ( 0 ) ;
    oadbMultSelection = (AliOADBMultSelection* )oadbContMS->GetObject(fRunNumber, "Default");
    if(!oadbMultSelection) {
        AliWarningF("OADB object for run %i not found! Aborting...", fRunNumber );
        return kFALSE;
    }
    AliInfoF("OADB object acquired for run %i, setting up",fRunNumber);
    fSelection         = oadbMultSelection->GetMultSelection();
    fMultSelectionCuts = oadbMultSelection->GetEventCuts();
    AliInfo("Setup complete!");
    cout<<" * Event Selection Read from OADB: "<<endl;
    fMultSelectionCuts -> Print();


    //============================================================
    // Calibration pre-optimization and setup
    //============================================================
    //Pre-optimize and create TFormulas
    fSelection->Setup ( fInput );
    //============================================================

    Long64_t lNEv = fTree->GetEntries();
    cout<<"(1) Data File opened, event count is "<<lNEv<<endl;

    cout<<"(2) Creating buffer, computing averages"<<endl;
    const int lMax = 1000;
    const int lMaxQuantiles = 10000;
    Int_t lRunNumbers[lMaxQuantiles];
    Long_t lRunStats[lMaxQuantiles];
    Int_t lNRuns = 0;
    Bool_t lNewRun = kTRUE;
    Int_t lThisRunIndex = -1;
    //Buffer file with run-by-run TTree objects needed for later processing

    Int_t iVtxZ_index = -1;
    TString lTempStr;

    TFile *fOutput = new TFile (fBufferFileNameData.Data(), "RECREATE");
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
            //For future reference
            lTempStr = fInput->GetVariable(iQvar)->GetName();
            if ( lTempStr.EqualTo("fEvSel_VtxZ") ) iVtxZ_index = iQvar;
        }
    }


    const int lNEstimators = fSelection->GetNEstimators();
    //For computing average values of estimators
    Double_t lAvEst[lNEstimators][lMax];
    //For computing extreme values (useful for integer calibration mode)
    Double_t lMaxEst[lNEstimators][lMax];
    Double_t lMinEst[lNEstimators][lMax];

    //Index of first value above anchor point threshold
    Long64_t lAnchorEst[lNEstimators][lMax];

    for(Long_t iEst=0; iEst<lNEstimators; iEst++) {
        for(Long_t iRun=0; iRun<lMax; iRun++) lAvEst[iEst][iRun] = 0;
        for(Long_t iRun=0; iRun<lMax; iRun++) lMaxEst[iEst][iRun] = -1e+3;
        for(Long_t iRun=0; iRun<lMax; iRun++) lMinEst[iEst][iRun] =  1e+6; //not more than a million, I hope?
        for(Long_t iRun=0; iRun<lMax; iRun++) lAnchorEst[iEst][iRun] = -1; //invalid index
    }

    //Add Timer
    TStopwatch* timer = new TStopwatch();
    timer->Start ( kTRUE );

    //Compute events-per-hour performance metric
    Double_t lEventsPerSecond = 0;

    //Setup Map: We want sTree and sTreeMC to share the same index <-> Run mapping...
    std::map<int, int> fRunMap;

    //==============================================================================
    // Data Loop for Run Number determination + Filtering
    //==============================================================================

    fOutput->cd();
    for(Long64_t iEv = 0; iEv<fTree->GetEntries(); iEv++) {

        if ( iEv % 100000 == 0 ) {
            Double_t complete = 100. * ( double ) ( iEv ) / ( double ) ( fTree->GetEntries() );
            cout << "(Data Loop 1) Event # " << iEv << "/" << fTree->GetEntries() << " (" << complete << "%, Time Left: ";
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
            //Add Map so that sTreeMC[i] matches run number with sTree[i] later
            fRunMap.insert( std::pair<int,int>(fRunNumber,lNRuns)); //[fRunNumber -> lNRuns] for later, please
            lRunNumbers[lNRuns] = fRunNumber;
            lThisRunIndex = lNRuns;
            lNRuns++;
        }

        //Perform Event selection
        Bool_t lSaveThisEvent = kTRUE; //let's be optimistic

        Float_t lLocalVtxZ = fInput->GetVariable(iVtxZ_index)->GetValue();

        //Check Selections as they are in the fMultSelectionCuts Object
        if( fMultSelectionCuts->GetTriggerCut()    && ! fEvSel_Triggered  ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO ) lSaveThisEvent = kFALSE;
        if( TMath::Abs(lLocalVtxZ) > fMultSelectionCuts->GetVzCut()      ) lSaveThisEvent = kFALSE;
        //ADD ME HERE: Tracklets Vs Clusters Cut?
        if( fMultSelectionCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins    ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster  ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices) lSaveThisEvent = kFALSE;
	if( fMultSelectionCuts->GetNonZeroNContribs()          &&  fnContributors < 1 ) lSaveThisEvent = kFALSE;

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

    //Write buffer to file
    for(Int_t iRun=0; iRun<lNRuns; iRun++) sTree[iRun]->Write();

    //==============================================================================
    // MC Loop for Run Number determination + Filtering
    //==============================================================================

    //Stop, reset Timer
    timer->Stop();
    timer->Start(kTRUE);

    TFile *fOutputMC = new TFile (fBufferFileNameMC.Data(), "RECREATE");
    TTree *sTreeMC[lMaxQuantiles];
    cout<<"Creating Trees..."<<endl;
    for(Int_t iRun=0; iRun<lMax; iRun++) {
        sTreeMC[iRun] = new TTree(Form("sTreeMC%i",iRun),Form("sTreeMC%i",iRun));
        for( Int_t iQvar = 0; iQvar<fInput->GetNVariables(); iQvar++) {
            if( !fInput->GetVariable(iQvar)->IsInteger() ) {
                sTreeMC[iRun]->Branch(Form("%s", fInput->GetVariable(iQvar)->GetName()  ),
                                      &fInput->GetVariable(iQvar)->GetRValue(),Form("%s/F",fInput->GetVariable(iQvar)->GetName()));
            } else {
                sTreeMC[iRun]->Branch(Form("%s", fInput->GetVariable(iQvar)->GetName()  ),
                                      &fInput->GetVariable(iQvar)->GetRValueInteger(),Form("%s/I",fInput->GetVariable(iQvar)->GetName()));
            }
        }
    }

    for(Long64_t iEv = 0; iEv<fTreeMC->GetEntries(); iEv++) {

        if ( iEv % 100000 == 0 ) {
            Double_t complete = 100. * ( double ) ( iEv ) / ( double ) ( fTreeMC->GetEntries() );
            cout << "(MC Loop 1) Event # " << iEv << "/" << fTreeMC->GetEntries() << " (" << complete << "%, Time Left: ";
            timer->Stop();
            Double_t time = timer->RealTime();

            //events per hour:
            lEventsPerSecond = ( ( Double_t ) ( iEv ) ) /time;

            timer->Start ( kFALSE );
            Double_t secondsperstep = time / ( Double_t ) ( iEv+1 );
            Double_t secondsleft = ( Double_t ) ( fTreeMC->GetEntries()-iEv-1 ) * secondsperstep;
            Long_t minutesleft = ( Long_t ) ( secondsleft / 60. );
            secondsleft = ( Double_t ) ( ( Long_t ) ( secondsleft ) % 60 );
            cout << minutesleft << "min " << secondsleft << "s, working at "<<lEventsPerSecond<<" Events/s..." << endl;
        }

        fTreeMC->GetEntry(iEv); //Look at next event

        //Perform Event selection
        Bool_t lSaveThisEvent = kTRUE; //let's be optimistic

        Float_t lLocalVtxZ = fInput->GetVariable(iVtxZ_index)->GetValue();

        //Check Selections as they are in the fMultSelectionCuts Object
        if( fMultSelectionCuts->GetTriggerCut()    && ! fEvSel_Triggered  ) lSaveThisEvent = kFALSE; //FIXME
        if( fMultSelectionCuts->GetINELgtZEROCut() && ! fEvSel_INELgtZERO ) lSaveThisEvent = kFALSE;
        if( TMath::Abs(lLocalVtxZ) > fMultSelectionCuts->GetVzCut()      ) lSaveThisEvent = kFALSE;
        //ADD ME HERE: Tracklets Vs Clusters Cut?
        if( fMultSelectionCuts->GetRejectPileupInMultBinsCut() && ! fEvSel_IsNotPileupInMultBins    ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetTrackletsVsClustersCut()    && ! fEvSel_PassesTrackletVsCluster  ) lSaveThisEvent = kFALSE;
        if( fMultSelectionCuts->GetVertexConsistencyCut()      && ! fEvSel_HasNoInconsistentVertices) lSaveThisEvent = kFALSE;
	if( fMultSelectionCuts->GetNonZeroNContribs()          &&  fnContributors < 1 ) lSaveThisEvent = kFALSE;

        //Consult map for run range equivalency
        Int_t lIndex = -1;
        if ( fRunMap.find( fRunNumber ) != fRunMap.end() ) {
            lIndex = fRunMap[ fRunNumber ];
        } else {
            lSaveThisEvent = kFALSE;
        }
        if ( lSaveThisEvent ) {
            sTreeMC [ lIndex ] -> Fill();
        }

        if(lNRuns>lMax) {
            lNRuns = lMax;
            AliWarningF("Exceeded maximum allowed number of runs to quantile! (Nruns now = %i)",lNRuns );
            AliWarningF("Will continue using only %i runs.", lMax);
            break;
        }
    }

    //Write buffer to file
    for(Int_t iRun=0; iRun<lNRuns; iRun++) sTreeMC[iRun]->Write();

    cout<<"(3) Inspect List of Runs and their corresponding statistics passing cuts: "<<endl;
    for(Int_t iRun = 0; iRun<lNRuns; iRun++) {
        cout<<" --- "<<lRunNumbers[iRun]<<", N_{events}^{data} = "<<sTree[iRun]->GetEntries()<<
            ", N_{events}^{MC} = "<<sTreeMC[iRun]->GetEntries()<<endl;
    }
    cout<<endl;

    //==============================================================================
    // Determine min and max values for all estimators for 2D correlation plot
    //==============================================================================

    Double_t *lValues;

    //Determine max SPD tracklets: requires finding which estimator is SPD tracklets
    Int_t iEstSPDtracklets     = -1;
    Int_t iEstSPDclusters      = -1;
    Int_t iEstSPDtrackletscorr = -1;
    Int_t iEstCL0  = -1;
    Int_t iEstCL1  = -1;

    //FIXME: This can be done in a loop, will be adjusted later (but this will work as is)

    TString lEstName;
    for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
        lEstName = fSelection->GetEstimator(iEst)->GetName();
        if ( lEstName.EqualTo("SPDTracklets") ) {
            cout<<"-> Found SPD tracklets at estimator index position "<<iEst<<"..."<<endl;
            iEstSPDtracklets = iEst;
        }
        if ( lEstName.EqualTo("SPDClusters") ) {
            cout<<"-> Found SPD clusters at estimator index position "<<iEst<<"..."<<endl;
            iEstSPDclusters = iEst;
        }
        if ( lEstName.EqualTo("SPDClustersCorr") ) {
            cout<<"-> Found SPD clusters (vz corrected) at estimator index position "<<iEst<<"..."<<endl;
            iEstSPDtrackletscorr = iEst;
        }
        if ( lEstName.EqualTo("CL0") ) {
            cout<<"-> Found CL0 at estimator index position "<<iEst<<"..."<<endl;
            iEstCL0 = iEst;
        }
        if ( lEstName.EqualTo("CL1") ) {
            cout<<"-> Found CL1 at estimator index position "<<iEst<<"..."<<endl;
            iEstCL1 = iEst;
        }
    }

    cout<<"(4) Look at average values"<<endl;
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        const Long64_t ntot = (Long64_t) sTree[iRun]->GetEntries();
        cout<<"(4) --- Processing run number "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<lNRuns<<"), with "<<ntot<<" events..."<<endl;
        sTree[iRun]->SetEstimate(ntot+1);
        //Cast Run Number into drawing conditions
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            lRunStats[iRun] = sTree[iRun]->Draw(fSelection->GetEstimator(iEst)->GetDefinition(),"","goff");
            lValues = sTree[iRun]->GetV1();
            cout<<"(4) --- Calculating values for "<<fSelection->GetEstimator(iEst)->GetName()<<":\t"<<flush;
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

    //==============================================================================
    // Create all 2D Correlation plots for scaling factors
    //==============================================================================

    Float_t lGlobalAxisScaling = 1.35;
    Float_t lRebinningSPD = 10;
    Float_t lRebinningEstimator = 10;
    TFile *fOutputDebug = new TFile(fDebugFileName.Data(), "RECREATE");

    //Create 2D correlation plots as needed
    TH2F *l2dTrackletVsEstimatorData [1000] [lNEstimators];
    TH2F *l2dTrackletVsEstimatorMC   [1000] [lNEstimators];
    for(Int_t ii=0; ii<1000; ii++)
        for(Int_t jj=0; jj<lNEstimators; jj++){
            l2dTrackletVsEstimatorData[ii][jj] = 0x0;
            l2dTrackletVsEstimatorMC[ii][jj] = 0x0;
        }

    cout<<"(5) Creating histograms..."<<endl;

    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            //May need later adjustment: axis scaling factor
            Float_t lSPDRange = (Int_t) ( lGlobalAxisScaling*lMaxEst[iEstSPDtracklets][iRun]+1.0 );
            Int_t lNbinsSPD = (Int_t) ( lSPDRange / lRebinningSPD );
            if ( !fSelection->GetEstimator(iEst)->IsInteger() ) {
                l2dTrackletVsEstimatorData[iRun][iEst] = new TH2F(Form("l2dTrackletVsEstimatorData_%i_%s",lRunNumbers[iRun],
                        fSelection->GetEstimator(iEst)->GetName()), "",
                        200, 0.0, lGlobalAxisScaling*lMaxEst[iEst][iRun],
                        lNbinsSPD, -0.5, -0.5 + lSPDRange );
                l2dTrackletVsEstimatorMC[iRun][iEst] = new TH2F(Form("l2dTrackletVsEstimatorMC_%i_%s",lRunNumbers[iRun],
                        fSelection->GetEstimator(iEst)->GetName()), "",
                        200, 0.0, lGlobalAxisScaling*lMaxEst[iEst][iRun],
                        lNbinsSPD, -0.5, -0.5 + lSPDRange );
            } else {
                Int_t lNbinsEst = (Int_t) (lGlobalAxisScaling*lMaxEst[iEst][iRun]+1.0);
                Float_t lLowEdgeEst = lMinEst[iEst][iRun]-0.5;
                Float_t lHighEdgeEst = lMinEst[iEst][iRun]-0.5+lNbinsEst;
                //Rebin, please!
                lNbinsEst = lNbinsEst/lRebinningEstimator;
                l2dTrackletVsEstimatorData[iRun][iEst] = new TH2F(Form("l2dTrackletVsEstimatorData_%i_%s",lRunNumbers[iRun],
                        fSelection->GetEstimator(iEst)->GetName()), "",
                        lNbinsEst, lLowEdgeEst, lHighEdgeEst,
                        lNbinsSPD, -0.5, -0.5 + lSPDRange );
                l2dTrackletVsEstimatorMC[iRun][iEst] = new TH2F(Form("l2dTrackletVsEstimatorMC_%i_%s",lRunNumbers[iRun],
                        fSelection->GetEstimator(iEst)->GetName()), "",
                        lNbinsEst, lLowEdgeEst, lHighEdgeEst ,
                        lNbinsSPD, -0.5, -0.5 + lSPDRange );
            }
        }
    }

    //Prepare 2D correlations
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            //Data Loop
            cout<<"Filling run "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<lNRuns<<"): Data..."<<flush;
            TString lEvalMe = "";
            lEvalMe.Append(Form("fnTracklets:%s>>l2dTrackletVsEstimatorData_%i_%s",((TString)fSelection->GetEstimator(iEst)->GetDefinition()).Data(),lRunNumbers[iRun], fSelection->GetEstimator(iEst)->GetName()));
            sTree[iRun] -> Draw( lEvalMe.Data() , "", "goff" );
            //MC loop
            cout<<"Done! MC..."<<flush;
            TString current = fSelection->GetEstimator(iEst)->GetDefinition();
            lEvalMe = "";
            //Ignore fired condition ONLY in MC
            current.ReplaceAll("fZnaFired", "1");
            current.ReplaceAll("fZncFired", "1");
            current.ReplaceAll("fZpaFired", "1");
            current.ReplaceAll("fZpcFired", "1");
            lEvalMe.Append(Form("fnTracklets:%s>>l2dTrackletVsEstimatorMC_%i_%s", current.Data(), lRunNumbers[iRun], fSelection->GetEstimator(iEst)->GetName()));
            sTreeMC[iRun] -> Draw( lEvalMe.Data() , "", "goff" );
            cout<<"Done!"<<endl;
        }
    }

    //Prepare 1D fits
    TProfile *profdata[1000][lNEstimators];
    TProfile *profmc[1000][lNEstimators];

    TF1 *fitdata[1000][lNEstimators];
    TF1 *fitmc[1000][lNEstimators];

    for(Int_t ii=0; ii<1000; ii++)
        for(Int_t jj=0; jj<lNEstimators; jj++){
            profdata[ii][jj]=0x0;
            profmc[ii][jj]=0x0;
            fitdata[ii][jj]=0x0;
            fitmc[ii][jj]=0x0;
        }

    TString fFormula = "[0]*x";

    //Experimental: make use of improved fitting procedure
    if(fkUseQuadraticMapping) fFormula = "[0]+[1]*TMath::Power(x,[2])";

    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {

            //Check if anchored, disregard anchored range if the case
            Double_t lLowestX=0.0;
            if(fSelection->GetEstimator(iEst)->GetUseAnchor()){
                lLowestX=fSelection->GetEstimator(iEst)->GetAnchorPoint(); //Remove lowest
            }

            cout<<"---> At Run "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<lNRuns<<"), estimator "<<fSelection->GetEstimator(iEst)->GetName()<<", fit range "<<lMaxEst[iEst][iRun]<<endl;
            if(!l2dTrackletVsEstimatorData[iRun][iEst]) cout<<"Null pointer for l2dTrackletVsEstimatorData"<<endl;
            if(!l2dTrackletVsEstimatorMC[iRun][iEst]) cout<<"Null pointer for l2dTrackletVsEstimatorMC"<<endl;
            cout<<"Profdata: setting up..."<<endl;
            profdata[ iRun ][ iEst ] = l2dTrackletVsEstimatorData[iRun][iEst]->ProfileX(Form("profdata_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ) ) ;
            cout<<"Profmc: setting up..."<<endl;
            profmc[ iRun ][ iEst ] = l2dTrackletVsEstimatorMC[iRun][iEst]->ProfileX(Form("profmc_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ) ) ;

            fitdata[iRun][iEst] = new TF1(Form("fitdata_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ), fFormula.Data(), lLowestX, lMaxEst[iEst][iRun]);
            fitmc[iRun][iEst] = new TF1(Form("fitmc_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ), fFormula.Data(), lLowestX, lMaxEst[iEst][iRun]);

            if(fkUseQuadraticMapping){
                fitdata[iRun][iEst]->SetParameter(0,0);
                fitdata[iRun][iEst]->SetParameter(1,1);
                fitdata[iRun][iEst]->SetParameter(2,1);
                fitmc[iRun][iEst]->SetParameter(0,0);
                fitmc[iRun][iEst]->SetParameter(1,1);
                fitmc[iRun][iEst]->SetParameter(2,1);
            }

            //Adjust range if needed
            //fitdata[iRun][iEst] -> SetRange(0,15000);
            //fitmc  [iRun][iEst] -> SetRange(0,15000);

            /*
             Power law logic

             ydata = Adata + Bdata * xdata ^ Cdata (1)
             ymc = Amc + Bmc * xmc ^ Cmc (2)

             From (1), isolate xdata:

             xdata = TMath::Power((ydata - Adata)/Bdata, 1./Cdata)

             ydata -> ymc = Amc + Bmc * xmc ^ Cmc

             */
            //Die hard fitting
            //TVirtualFitter::SetMaxIterations(1000000);

            //Guess initial parameters more wisely: linear approximation, please

            Double_t lIncline   = 1e-3;
            Double_t lInclineMC = 1e-3;

            if( !profdata[iRun][iEst]) cout<<"Null pointer!"<<endl;
            if( !profmc[iRun][iEst]) cout<<"Null pointer!"<<endl;
            if( TMath::Abs(lAvEst[iEst][iRun])>1e-3 ){
                lIncline   = profdata[iRun][iEst]->GetBinContent( profdata[iRun][iEst]->FindBin(lAvEst[iEst][iRun]) ) / lAvEst[iEst][iRun];
                lInclineMC = profmc  [iRun][iEst]->GetBinContent( profmc  [iRun][iEst]->FindBin(lAvEst[iEst][iRun]) ) / lAvEst[iEst][iRun];
            }

            if( !fitdata[iRun][iEst]) cout<<"Null pointer!"<<endl;
            if( !fitmc[iRun][iEst]) cout<<"Null pointer!"<<endl;

            fitdata[iRun][iEst]->SetParameter(0,lIncline);
            fitmc[iRun][iEst]->SetParameter(0,lIncline);
            //fitdata[iRun][iEst]->SetParameter(1,lMaxEst[iEst][iRun]*5);
            //fitmc[iRun][iEst]->SetParameter(1,lMaxEst[iEst][iRun]*5);
            //fitdata[iRun][iEst]->SetParameter(2,0.0);
            //fitmc[iRun][iEst]->SetParameter(2,0.0);

            //remember to not be silly...
            TString lEstName = fSelection->GetEstimator(iEst)->GetName();
            if( !lEstName.Contains("SPD") &&
               !lEstName.Contains("CL0") &&
               !lEstName.Contains("CL1") ){
                cout<<"Fit DATA: "<<endl;
                profdata[iRun][iEst] -> Fit( Form("fitdata_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ), "IREM0" );
                cout<<"Fit MONTE CARLO: "<<endl;
                profmc[iRun][iEst] -> Fit( Form("fitmc_%i_%s",lRunNumbers[iRun],fSelection->GetEstimator(iEst)->GetName() ), "IREM0" );
            }
        }
    }

    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        cout<<"Run Report: "<<lRunNumbers[iRun]<<" ("<<iRun<<"/"<<lNRuns<<"): "<<endl;
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            cout<<"Estimator: "<<fSelection->GetEstimator(iEst)->GetName()<<" Data: "<<fitdata[iRun][iEst]->GetParameter(0)<<" MC: "<<fitmc[iRun][iEst]->GetParameter(0)<<endl;
        }
    }

    //Create Histograms for storing scaling factors ...
    TH1F *hScaleFactor[lNEstimators];
    Float_t lScaleFactors     [lNEstimators][lNRuns];
    Float_t lScaleFactorsError[lNEstimators][lNRuns];
    for(Int_t iEst=0; iEst<lNEstimators; iEst++){
        hScaleFactor[iEst] = new TH1F(Form("hScaleFactor_%s",fSelection->GetEstimator(iEst)->GetName()),"",lNRuns,0,lNRuns);
        for(Int_t iRun=0; iRun<lNRuns; iRun++){
            hScaleFactor[iEst]->GetXaxis()->SetBinLabel(iRun+1, Form("%i",lRunNumbers[iRun]) );
            lScaleFactors[iEst][iRun] = -1;
            lScaleFactorsError[iEst][iRun] = -1e-6;
            if ( sTreeMC[iRun]->GetEntries() > 10 && sTree[iRun]->GetEntries() > 10 ){
                if( TMath::Abs(fitdata[iRun][iEst]->GetParameter(0))>1e-6 ){
                    Float_t lkAerr = fitmc[iRun][iEst]->GetParError(0);
                    Float_t lkBerr = fitdata[iRun][iEst]->GetParError(0);
                    Float_t lkA = fitmc[iRun][iEst]->GetParameter(0);
                    Float_t lkB = fitdata[iRun][iEst]->GetParameter(0);

                    //Central value
                    lScaleFactors[iEst][iRun] = lkA / lkB;

                    //Standard Error Propagation
                    Float_t errorfromtop = lkAerr*lkAerr / (lkB*lkB) ;
                    Float_t errorfrombottom = ((lkA*lkA)/(lkB*lkB*lkB*lkB)) * lkBerr * lkBerr;
                    lScaleFactorsError[iEst][iRun] = TMath::Sqrt( errorfromtop + errorfrombottom );

                }
            }
            hScaleFactor[iEst]->SetBinContent(iRun+1, lScaleFactors[iEst][iRun]);
            hScaleFactor[iEst]->SetBinError(iRun+1, lScaleFactorsError[iEst][iRun]);
        }
    }

    cout<<"Write output..."<<endl;

    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            l2dTrackletVsEstimatorData[iRun][iEst]->Write(l2dTrackletVsEstimatorData[iRun][iEst]->GetName(), TObject::kOverwrite);
            profdata[iRun][iEst]->Write(profdata[iRun][iEst]->GetName(), TObject::kOverwrite);
            fitdata[iRun][iEst]->Write(fitdata[iRun][iEst]->GetName(), TObject::kOverwrite);
        }
    }
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        for(Int_t iEst=0; iEst<lNEstimators; iEst++) {
            l2dTrackletVsEstimatorMC[iRun][iEst]->Write(l2dTrackletVsEstimatorMC[iRun][iEst]->GetName(), TObject::kOverwrite);
            profmc[iRun][iEst]->Write(profmc[iRun][iEst]->GetName(), TObject::kOverwrite);
            fitmc[iRun][iEst]->Write(fitmc[iRun][iEst]->GetName(), TObject::kOverwrite);
        }
    }

    cout<<"Write Debug file..."<<endl;
    fOutputDebug->Write();
    fOutputDebug->Close();

    //Save scaling factors as pre-pend to estimator evaluation for the OADB!

    //Create Stuff
    TFile * f = new TFile (fOutputFileName.Data(), "recreate");
    AliOADBContainer * oadbContMSout = new AliOADBContainer("MultSel");

    AliOADBMultSelection * oadbMultSelectionout = 0x0;
    AliMultSelectionCuts * cuts              = 0x0;
    AliMultSelection     * fsels             = 0x0;

    AliOADBMultSelection * oadbMultSelectionoutdef = 0x0;
    AliMultSelectionCuts * cutsdef              = 0x0;
    AliMultSelection     * fselsdef             = 0x0;

    //Default Stuff
    TH1F * hDummy[lNEstimators];
    for ( Int_t iEst=0; iEst<lNEstimators; iEst++) {
        hDummy[iEst]= new TH1F (Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()), "hdummy", 1, 0, 1);
        hDummy[iEst]->SetDirectory(0);
    }

    //Actual Calibration Histograms
    TH1F * hCalibData[lNEstimators];

    //Loop over existing runs and write objects as needed
    for(Int_t iRun=0; iRun<lNRuns; iRun++) {
        cout<<"Processing run number "<<lRunNumbers[iRun]<<endl;
        if ( sTreeMC[iRun]->GetEntries() > 10 && sTree[iRun]->GetEntries() > 10 ){
            oadbMultSelectionout = new AliOADBMultSelection();
            cuts              = new AliMultSelectionCuts();
            cuts = fMultSelectionCuts;
            fsels             = new AliMultSelection( fSelection );

            //Write in scaling factors!
            TString lTempDef;
            for(Int_t iEst=0; iEst<lNEstimators; iEst++){
                lTempDef = fsels->GetEstimator( iEst )->GetDefinition();
                lTempDef.ReplaceAll("fZnaFired", "1");
                lTempDef.ReplaceAll("fZncFired", "1");
                lTempDef.ReplaceAll("fZpaFired", "1");
                lTempDef.ReplaceAll("fZpcFired", "1");

                //remember to not be silly...
                TString lEstName = fsels->GetEstimator(iEst)->GetName();
                //Construction of estimator re-definition
                if( !lEstName.Contains("SPD") &&
                   !lEstName.Contains("CL0") &&
                   !lEstName.Contains("CL1") ){
                    if(!fkUseQuadraticMapping){
                        lTempDef.Prepend(Form("%.10f*(",lScaleFactors[iEst][iRun] ));
                        lTempDef.Append(")"); //don't forget parentheses...
                    }else{
                        //Experimental quadratic fit
                        TString lTemporary = lTempDef.Data();

                        /*
                         //Substitution logic:
                         xdata = TMath::Power((ydata - Adata)/Bdata, 1./Cdata)
                         ydata -> ymc = Amc + Bmc * xmc ^ Cmc

                         full formula:

                         xdata = TMath::Power(( Amc + Bmc * TMath::Power(xmc,Cmc) - Adata )/Bdata, 1./Cdata);
                         */

                        //lTempDef = Form(TMath::Power(( (Amc + Bmc * TMath::Power(xmc,Cmc)) - Adata )/Bdata, 1./Cdata),
                        lTempDef = Form("(TMath::Power(( (%.10f + %.10f * TMath::Power(%s,%.10f)) - %.10f )/%.10f, 1./%.10f))",
                                        fitmc[iRun][iEst]->GetParameter(0),
                                        fitmc[iRun][iEst]->GetParameter(1),
                                        lTempDef.Data(),
                                        fitmc[iRun][iEst]->GetParameter(2),
                                        fitdata[iRun][iEst]->GetParameter(0),
                                        fitdata[iRun][iEst]->GetParameter(1),
                                        fitdata[iRun][iEst]->GetParameter(2));
                        lTempDef.ReplaceAll("ESTIMATOR",lTemporary.Data());

                        //Trick: add anchor-point-like logic to avoid low-multiplicity issues
                        if(fsels->GetEstimator(iEst)->GetUseAnchor()){
                            lTempDef.Prepend(Form("(%s>%.10f)?(",fsels->GetEstimator(iEst)->GetDefinition().Data(),0.5*fsels->GetEstimator(iEst)->GetAnchorPoint()));
                            lTempDef.Append("):0.0");
                        }

                        //Use V0M high-cutoff if provided
                        if( lEstName.Contains("V0M"))
                            if ( fV0MCutoffs.find( lRunNumbers[iRun] ) != fV0MCutoffs.end() ) {
                                Float_t lV0MCutoff = fV0MCutoffs[ lRunNumbers[iRun] ];
                                cout<<"-> Run "<<lRunNumbers[iRun]<<": will use V0M cutoff of "<<lV0MCutoff<<endl;
                                lTempDef.Prepend(Form("(((fAmplitude_V0A)+(fAmplitude_V0C))<%.5f)?(", lV0MCutoff));
                                lTempDef.Append("):0.0"); //reject otherwise
                            }
                    }

                    cout<<"================================================================================"<<endl;
                    cout<<" Quadratic fit print obtained for estimator "<<fsels->GetEstimator( iEst )->GetName()<<endl;
                    cout<<lTempDef.Data()<<endl;
                    cout<<"================================================================================"<<endl;
                }

                fsels->GetEstimator( iEst )->SetDefinition ( lTempDef.Data() );
            }

            oadbMultSelectionout->SetEventCuts    (cuts );
            oadbMultSelectionout->SetMultSelection(fsels);
            for ( Int_t iEst=0; iEst<lNEstimators; iEst++) {
                //Average values
                fsels->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );

                //Protect from crashes in the general case
                hCalibData[iEst] = (TH1F*) hDummy[iEst]->Clone(Form("hCalib_%s",fSelection->GetEstimator(iEst)->GetName()) );
                oadbMultSelectionout->AddCalibHisto( hCalibData[iEst]);
                hCalibData[iEst]->SetDirectory(0);
            }
            cout<<"=================================================================================="<<endl;
            cout<<"AliMultSelection Object to be saved for run "<<lRunNumbers[iRun]<<": "<<endl;
            fsels->PrintInfo();
            cuts->Print();
            cout<<"=================================================================================="<<endl;
            oadbContMSout->AppendObject(oadbMultSelectionout, lRunNumbers[iRun] ,lRunNumbers[iRun] );

            //Check if this corresponds to the reference run !
            Bool_t lThisIsReference = kFALSE;
            if( lRunNumbers[iRun] == fRunToUseAsDefault ) lThisIsReference = kTRUE;
            if( lThisIsReference ){
                cout<< "Reference Run found! Will save..."<<endl;
                //========================================================================
                //DEFAULT OADB Object saving procedure STARTS here
                oadbMultSelectionoutdef = new AliOADBMultSelection("Default");
                cutsdef              = new AliMultSelectionCuts();
                cutsdef = fMultSelectionCuts;
                fselsdef             = new AliMultSelection( fSelection );

                //Write in scaling factors!
                TString lTempDef;
                for(Int_t iEst=0; iEst<lNEstimators; iEst++){
                    lTempDef = fselsdef->GetEstimator( iEst )->GetDefinition();
                    //Construction of estimator re-definition
                    TString lEstName = fsels->GetEstimator(iEst)->GetName();
                    if( !lEstName.Contains("SPD") &&
                       !lEstName.Contains("CL0") &&
                       !lEstName.Contains("CL1") ){
                        if(!fkUseQuadraticMapping){
                            lTempDef.Prepend(Form("%.10f*(",lScaleFactors[iEst][iRun] ));
                            lTempDef.Append(")"); //don't forget parentheses...
                        }else{
                            //Experimental quadratic fit
                            TString lTemporary = lTempDef.Data();

                            /*
                             //Substitution logic:
                             xdata = TMath::Power((ydata - Adata)/Bdata, 1./Cdata)
                             ydata -> ymc = Amc + Bmc * xmc ^ Cmc

                             full formula:

                             xdata = TMath::Power(( Amc + Bmc * TMath::Power(xmc,Cmc) - Adata )/Bdata, 1./Cdata);
                             */

                            lTempDef = Form("TMath::Power(( (%.10f + %.10f * TMath::Power(%s,%.10f)) - %.10f )/%.10f, 1./%.10f)",
                                            fitmc[iRun][iEst]->GetParameter(0),
                                            fitmc[iRun][iEst]->GetParameter(1),
                                            lTempDef.Data(),
                                            fitmc[iRun][iEst]->GetParameter(2),
                                            fitdata[iRun][iEst]->GetParameter(0),
                                            fitdata[iRun][iEst]->GetParameter(1),
                                            fitdata[iRun][iEst]->GetParameter(2));
                            lTempDef.ReplaceAll("ESTIMATOR",lTemporary.Data());

                            //Trick: add anchor-point-like logic to avoid low-multiplicity issues
                            if(fsels->GetEstimator(iEst)->GetUseAnchor()){
                                lTempDef.Prepend(Form("(%s>%.10f)?(",fselsdef->GetEstimator(iEst)->GetDefinition().Data(),0.5*fselsdef->GetEstimator(iEst)->GetAnchorPoint()));
                                lTempDef.Append("):0.0");
                            }

                            //Use V0M high-cutoff if provided
                            if( lEstName.Contains("V0M"))
                                if ( fV0MCutoffs.find( lRunNumbers[iRun] ) != fV0MCutoffs.end() ) {
                                    Float_t lV0MCutoff = fV0MCutoffs[ lRunNumbers[iRun] ];
                                    cout<<"-> Run "<<lRunNumbers[iRun]<<": will use V0M cutoff of "<<lV0MCutoff<<endl;
                                    lTempDef.Prepend(Form("(((fAmplitude_V0A)+(fAmplitude_V0C))<%.5f)?(", lV0MCutoff));
                                    lTempDef.Append("):0.0"); //reject otherwise
                                }
                        }

                        cout<<"================================================================================"<<endl;
                        cout<<" Quadratic fit print obtained for estimator "<<fsels->GetEstimator( iEst )->GetName()<<endl;
                        cout<<lTempDef.Data()<<endl;
                        cout<<"================================================================================"<<endl;
                    }

                    //if ZxxFired included in the estimator, ignore it
                    lTempDef.ReplaceAll("fZnaFired", "1");
                    lTempDef.ReplaceAll("fZncFired", "1");
                    lTempDef.ReplaceAll("fZpaFired", "1");
                    lTempDef.ReplaceAll("fZpcFired", "1");

                    fselsdef->GetEstimator( iEst )->SetDefinition ( lTempDef.Data() );
                }

                oadbMultSelectionoutdef->SetEventCuts    (cutsdef );
                oadbMultSelectionoutdef->SetMultSelection(fselsdef);
                for ( Int_t iEst=0; iEst<lNEstimators; iEst++) {
                    //Average values
                    fselsdef->GetEstimator(iEst)->SetMean( lAvEst[iEst][iRun] );

                    //Protect from crashes in the general case
                    //hCalibData[iEst] = (TH1F*) hDummy[iEst]->Clone(Form("hCalib_%s",fselsdef->GetEstimator(iEst)->GetName()) );
                    oadbMultSelectionoutdef->AddCalibHisto( hDummy[iEst]);
                    hDummy[iEst]->SetDirectory(0);
                }
                cout<<"=================================================================================="<<endl;
                cout<<" Detected that this particular run / run range is special, will save it as default"<<endl;
                cout<<" AliMultSelection Object to be saved (DEFAULT)"<<endl;
                fselsdef->PrintInfo();
                cout<<"=================================================================================="<<endl;

                //for ( Int_t iEst=0; iEst<lNEstimators; iEst++) oadbMultSelectionoutdef->AddCalibHisto( hDummy[iEst] );
                oadbContMSout->AddDefaultObject(oadbMultSelectionoutdef);
                //DEFAULT OADB Object saving procedure ENDS here
                //========================================================================
            }
        }else{
            cout<<"Insufficient statistics in either MC or data in run #"<<lRunNumbers[iRun]<<", skipping!"<<endl;
        }
    }
    cout<<"Write OADB..."<<endl;
    //pre-write dump
    //fsels->PrintInfo();

    oadbContMSout->Write();
    f->Write();
    f->Close();

    cout<<" Done!"<<endl;
    return kTRUE;
}
//________________________________________________________________
Float_t AliMultSelectionCalibratorMC::MinVal( Float_t A, Float_t B ) {
    if( A < B ) {
        return A;
    }
    else {
        return B;
    }
}
//________________________________________________________________
void AliMultSelectionCalibratorMC::SetupStandardInput() {
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
    AliMultVariable *fRefMultEta14     = new AliMultVariable("fRefMultEta14");
    fRefMultEta14->SetIsInteger( kTRUE );
    AliMultVariable *fnTracklets     = new AliMultVariable("fnTracklets");
    fnTracklets->SetIsInteger( kTRUE );

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

    //vertex-Z
    AliMultVariable *fEvSel_VtxZ = new AliMultVariable("fEvSel_VtxZ");

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
    fInput->AddVariable( fRefMultEta5  );
    fInput->AddVariable( fRefMultEta8  );
    fInput->AddVariable( fRefMultEta14  );
    fInput->AddVariable( fEvSel_VtxZ  );
    //============================================================

}

//________________________________________________________________
void AliMultSelectionCalibratorMC::AddV0MCutoff ( Int_t lRunNumber, Float_t lV0MCutoff ){
    //Add mapping : all runs in range go to current value of fNRunRanges
    //Ease of access
    fV0MCutoffs.insert( std::pair<int,float>(lRunNumber,lV0MCutoff));
    fNV0MCutoffs++;
    AliInfoF("Added V0M cutoff for run #%i at %.5f", (Int_t)lRunNumber, (Float_t) lV0MCutoff) ;
}

