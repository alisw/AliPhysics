//------------------------------------
// Configuration macro:
//
// Configure JETAN FastJet analysis.
//
// Modified by: K. Read
//
//------------------------------------

AliJetFinder*  ConfigJetAnalysis()
{
    //
    // Configuration goes here
    // 
    printf("========================== \n");
    printf("ConfigJetAnalysisFastJet() \n");
    printf("========================== \n");

    Bool_t kInputIsESD = kTRUE;     //uncomment for input ESD
  //Bool_t kInputIsESD = kFALSE;    //uncomment for input AODs
    Bool_t kFollowsFilter = kTRUE;  //uncomment if follows ESD filter task
  //Bool_t kFollowsFilter = kFALSE; //uncomment if no ESD filter task

    //Alternatively, select input via anaInputData environment variable.
    if (gSystem->Getenv("anaInputData")){
      TString kInputData = gSystem->Getenv("anaInputData");
      if( kInputData == "AOD" ){
        kInputIsESD = kFALSE;
        kFollowsFilter = kFALSE;
      }
    }

    // Define the grids
    AliJetGrid *grid = new AliJetGrid(419,119,0.,2*TMath::Pi(),-0.9,0.9); 
    grid->SetGridType(1);
    grid->InitParams(80.*TMath::Pi()/180,190.*TMath::Pi()/180,-0.7,0.7); 
    grid->SetMatrixIndexes();
    grid->SetIndexIJ();
    AliJetGrid *grid2 = new AliJetGrid(131,95,80.*TMath::Pi()/180.,190.*TMath::Pi()/180.,-0.7,0.7); 
    grid2->SetGridType(0);
    grid2->SetMatrixIndexes();
    grid2->SetIndexIJ();

    // Define reader header
    if(kInputIsESD && !kFollowsFilter) AliJetESDReaderHeader *jrh = new AliJetESDReaderHeader();
    else            AliJetAODReaderHeader *jrh = new AliJetAODReaderHeader();
    jrh->SetComment("Testing");
    if(kInputIsESD && !kFollowsFilter) jrh->SetReadSignalOnly(kFALSE);

    // Detector options: 0 = Charged particles only (MomentumArray)
    //                   1 = Charged particles only (UnitArray)
    //                   2 = Neutral cells only (UnitArray)
    //                   3 = Charged particles + neutral cells (UnitArray)
    jrh->SetDetector(3);
    //jrh->SetDebug(-1);
    //jrh->SetFiducialEta(-0.7,0.7);
    //jrh->SetFiducialPhi(80.*TMath::Pi()/180,190.*TMath::Pi()/180);
    jrh->SetPtCut(0.1);
    jrh->SetFiducialEta(-0.9,0.9);    //fiducial range used by AliJetFillUnitArrayTracks
    jrh->SetFiducialPhi(0,2*TMath::Pi()); //fiducial range used by AliJetFillUnitArrayTracks
    
    // Define reader and set its header
    if(kInputIsESD && !kFollowsFilter) AliJetESDReader *er = new AliJetESDReader();
    else            AliJetAODReader *er = new AliJetAODReader();
    er->SetReaderHeader(jrh);
    er->SetTPCGrid(grid);
    er->SetEMCalGrid(grid2);
    er->SetApplyMIPCorrection(kFALSE);

    // Define jet header
    AliFastJetHeaderV1 *jh=new AliFastJetHeaderV1();
    Double_t R=0.4;
    Double_t Rbkg=0.2;

    // AliFastJetHeaderV1 *jh=new AliFastJetHeaderV1();
    jh->SetComment("Fast jet code with default parameters");
    //jh->SetDebug(-1);
    //jh->SetBGMode(1); //Do BG Subtraction
    jh->SetBGMode(0); //No BG Subtraction.  Store AOD track refs.
    jh->SetRparam(R); // setup parameters
    jh->SetRparamBkg(Rbkg); // setup parameters
    jh->SetPtMin(0.2);
    //jh->SetGhostEtaMax(0.9);
    jh->SetGhostArea(0.01);
    jh->SetGhostEtaMax(0.7);//used to set the rap_min and rap_max, that are then used by FJ
    jh->SetPhiRange(80.*TMath::Pi()/180+R,190.*TMath::Pi()/180-R);//used in AliFastJetFinder for the range

    // Define jet finder. Set its header and reader
    jetFinder = new AliFastJetFinder();
    jetFinder->SetJetHeader(jh);
    jetFinder->SetJetReader(er);
    //jetFinder->SetPlotMode(kTRUE);

    printf("============================== \n");
    printf("END ConfigJetAnalysisFastJet() \n");
    printf("============================== \n");
 
    return jetFinder;
}
