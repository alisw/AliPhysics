//#include "AliDielectron.h"
//#include "PWGDQ/dielectron/macrosLMEE/LMEECutLib.C"
#include<string>

void      InitHistograms(AliDielectron *die, Bool_t doPairing);
void      InitCF(AliDielectron* die, Int_t cutDefinition);
TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max);
TVectorD *GetVector(Int_t var);
enum {kMee=0, kMee500, kPtee, kP2D, kRuns, kPhiV, kOpAng, kOpAng2, kEta2D, kEta3D, kSigmaEle, kSigmaOther, kTPCdEdx, kCent, kPhi2D};

TString names=("all;electrons;highMult;midMult;lowMult");
//TString names=("all;electrons;lowPt;midLowPt;midPt;midHighPt;highPt");
//TString names = ("all;electrons;lowPt;midLowPt;midPt;midHighPt;highPt;highMult;midMult;lowMult");
//TString names=("all;electrons");
TObjArray *arrNames = names.Tokenize(";");
const Int_t nDie = arrNames->GetEntries();
Bool_t MCenabled = kFALSE; //Needed for LMEEcutlib
Bool_t isQAtask = kTRUE;
Int_t selectedPID = -1;
Bool_t pairCuts = kTRUE;

AliDielectron* Config_acapon(TString cutDefinition, Bool_t hasMC=kFALSE, Bool_t isESD=kFALSE, Bool_t SDDstatus =kFALSE, Bool_t doPairing = kTRUE, Bool_t doMixing = kTRUE)
{

    //Setup the instance of AliDielectron
    LMEECutLib*  LMcutlib = new LMEECutLib(SDDstatus);

    //Task name
    TString name = cutDefinition;

    //Init AliDielectron
    AliDielectron *die = new AliDielectron(Form("%s",name.Data()), Form("AliDielectron with cuts: %s",name.Data()));
    //die->SetHasMC(hasMC);
    MCenabled=hasMC;

    // deactivate pairing to check track cuts or run with loose pid cuts:
    if(!doPairing){
        die->SetNoPairing();
    }
 
    die->SetPreFilterUnlikeOnly(kTRUE);

    cout << "cutDefinition = " << cutDefinition << endl;
    // Setup Analysis Selection
    if(cutDefinition == "all"){
        selectedPID = LMEECutLib::kAllSpecies;
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "electrons"){
        selectedPID = LMEECutLib::kElectrons;
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "highMult"){
        selectedPID = LMEECutLib::kHighMult;
        die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "midMult"){
        selectedPID = LMEECutLib::kMidMult;
        die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "lowMult"){
        selectedPID = LMEECutLib::kLowMult;
        die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "lowPt"){
        selectedPID = LMEECutLib::kLowPt;
        //die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "midLowPt"){
        selectedPID = LMEECutLib::kMidLowPt;
        //die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "midPt"){
        selectedPID = LMEECutLib::kMidPt;
        //die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "midHighPt"){
        selectedPID = LMEECutLib::kMidHighPt;
        //die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else if(cutDefinition == "highPt"){
        selectedPID = LMEECutLib::kHighPt;
        //die->GetEventFilter().AddCuts( LMcutlib->GetCentralityCuts(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetKineCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetTrackCutsAna(selectedPID) );
        die->GetTrackFilter().AddCuts( LMcutlib->GetPIDCutsAna(selectedPID) );
        if(pairCuts){
            //die->GetPairPreFilter().AddCuts( LMcutlib->GetPairCutsPre(selectedPID) );
            die->GetPairFilter().AddCuts( LMcutlib->GetPairCutsAna(selectedPID) );
        }
    }
    else{
        cout << " =============================== " << endl;
        cout << " ==== INVALID CONFIGURATION ==== " << endl;
        cout << " cutDefinition = " << cutDefinition << endl;
        cout << " =============================== " << endl;
    }

    die->SetUseKF(kFALSE);

    AliDielectronMixingHandler* mix = 0x0;
    if(doMixing){
        mix = LMcutlib->GetMixingHandler(selectedPID);
        die->SetMixingHandler(mix);
    }

    InitHistograms(die, doPairing);

    return die;
}

//______________________________________________________________________________________

void InitHistograms(AliDielectron *die, Bool_t doPairing)
{
    //Define histogram names based on cut value, in order to avoid mem. warning error

    //Setup histogram Manager
    AliDielectronHistos *histos = new AliDielectronHistos(die->GetName(),die->GetTitle());

    //Initialise histogram classes
    histos->SetReservedWords("Track;Pair;Pre;RejTrack;RejPair");//;Track_Legs");//RejPair;RejTrack");

    //Event class
    histos->AddClass("Event");

    //Track classes
    //0,1: +- ev1, 2,3: +- ev2
    for(Int_t i = 0; i < 2; ++i){
        histos->AddClass(Form("Track_%s",AliDielectron::TrackClassName(i)));
    }

    //Pair classes
    // to fill also mixed event histograms loop until 10
    // fgkPairClassNames[11] = {
    //  "ev1+_ev1+",  "ev1+_ev1-",  "ev1-_ev1-",  // 0-2 (same event)
    //  "ev1+_ev2+",  "ev1-_ev2+",  "ev2+_ev2+",  // 3-4 (+5)
    //  "ev1+_ev2-",  "ev1-_ev2-",                // 6-7
    //  "ev2+_ev2-",  "ev2-_ev2-",  "ev1+_ev1-_TR"
    // };
    for(Int_t i = 0; i < 3; ++i){
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(i)));
        // Legs of final Pairs. Both charges together. No duplicate entries.
        histos->AddClass(Form("Track_Legs_%s",AliDielectron::PairClassName(i))); // not TrackClassName, see 'AliDielectron::FillHistograms(...)'
    }

    //ME and track rot
    if(die->GetMixingHandler()){
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(3)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(4)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(6)));
        histos->AddClass(Form("Pair_%s",AliDielectron::PairClassName(7)));
    }


    TH1::AddDirectory(kFALSE);
    //add histograms to event class
    histos->UserHistogram("Event","nEvents","",1,0.,1.,AliDielectronVarManager::kNevents);
    histos->UserHistogram("Event","Centrality","",100,0,100,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","nESDTracks","",500,0,500,AliDielectronVarManager::kNTrk);
    histos->UserHistogram("Event","zVertexPrimary","",122,-11,11,AliDielectronVarManager::kZvPrim);
    histos->UserHistogram("Event","NVtxContrib","Number of Vertex Contributor;N of Vertex Contributors;N of events",200,-0.5,199.5,AliDielectronVarManager::kNVtxContrib);
    //------ Num. tracks -----/
    histos->UserHistogram("Event","Accepted tracks","",50,0,50,AliDielectronVarManager::kNacc);
    histos->UserHistogram("Event","Ntracks","Number of tracks;N of tracks;N of events", 100, -0.5, 99.5,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","NtracksVsVtxZ","N tracks vs VtxZ;Vertex Z [cm];N of tracks",150,-15,15,50,-0.5,49.5,AliDielectronVarManager::kZvPrim,AliDielectronVarManager::kTracks);
    histos->UserHistogram("Event","RefMultTPConly","",300,0,300,AliDielectronVarManager::kRefMultTPConly);
    histos->UserHistogram("Event","SPD clusters vs tracklets","",75,0,150,30,0,60,AliDielectronVarManager::kNaccTrcklts10,AliDielectronVarManager::kITSLayerFirstCls);

    //--------- V0 plots ------------------------//
    histos->UserHistogram("Event","MultV0","Multiplicity V0;V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kMultV0);
    histos->UserHistogram("Event","EqMultV0","Equalized Multiplicity V0;Equalized V0M amplitude",4000,-0.5,3999.5,AliDielectronVarManager::kEqMultV0);
    histos->UserHistogram("Event","ChMultV0","Charged Multiplicity V0;Charged V0M amplitude",1000,-0.5,999.5,AliDielectronVarManager::kVZEROchMult);
    histos->UserHistogram("Event","CentralityV0M","Centrality V0;V0M percentile",300,-50,250,AliDielectronVarManager::kCentralityNew);
    histos->UserHistogram("Event","CentralityV0Mzoomed","Centrality V0 zoomed;V0M percentile",200,0,2,AliDielectronVarManager::kCentralityNew);


    //add histograms to Track classes
    histos->UserHistogram("Track","Pt",";Pt [GeV];#tracks",100,0,5.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Track","Px",";Px [GeV];#tracks",100,0,5.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Track","Py",";Py [GeV];#tracks",100,0,5.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Track","Pz",";Pz [GeV];#tracks",100,0,5.,AliDielectronVarManager::kPz);
    histos->UserHistogram("Track","P_PIn",";p (GeV/c);p_{in} (GeV/c)",100,0,10,AliDielectronVarManager::kPIn);
    // Eta and Phi
    histos->UserHistogram("Track","Eta","",200,-2,2,AliDielectronVarManager::kEta);
    histos->UserHistogram("Track","Phi","",120,0.,TMath::TwoPi(),AliDielectronVarManager::kPhi);
    histos->UserHistogram("Track","Eta_Phi","",100,-1,1,120,0,TMath::TwoPi(),AliDielectronVarManager::kEta,AliDielectronVarManager::kPhi);

    //Check ITS hits per run
    /*histos->UserHistogram("Track", "nITS vs run","", GetVector(kRuns), BinsToVector(6, -0.5, 6.5), AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kNclsITS);
    histos->UserHistogram("Track", "ITSclusterMap vs run","", GetVector(kRuns), BinsToVector(70, 0, 70), AliDielectronVarManager::kRunNumber, AliDielectronVarManager::kITSclusterMap);*/

    if(isQAtask){
        // DCA
        histos->UserHistogram("Track","dXY","",400,-2.,2.,AliDielectronVarManager::kImpactParXY);
        histos->UserHistogram("Track","dZ" ,"",600,-4.,4.,AliDielectronVarManager::kImpactParZ);
        histos->UserHistogram("Track","dXY_dZ","",100,-1.,1.,150,-3.,3.,AliDielectronVarManager::kImpactParXY,AliDielectronVarManager::kImpactParZ);
        histos->UserHistogram("Track","SPD clusters vs. tracklets",";tracklets;SPD clusters",
                            150,0,150,6,0,6,AliDielectronVarManager::kNTrk ,AliDielectronVarManager::kITSLayerFirstCls);
        histos->UserHistogram("Track","DCA_{xy} vs p_T","",300,0,5,100,0,0.5,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY);
        histos->UserHistogram("Track","DCA_{Z} vs p_T","",300,0,5,100,-1,1,AliDielectronVarManager::kPt,AliDielectronVarManager::kImpactParXY);

        //Track cut variables for trackQA
        //ITS
        histos->UserHistogram("Track","ITSnCls",";ITS number clusters;#tracks",6,-0.5,6.5,AliDielectronVarManager::kNclsITS);
        histos->UserHistogram("Track","ITSnClsClusterMap",";ITS cluster map;#tracks",100, 0.0, 100.0,AliDielectronVarManager::kITSclusterMap);
        histos->UserHistogram("Track","ITSchi2",";ITS chi2/Cl;#tracks",110,0.,11.,AliDielectronVarManager::kITSchi2Cl);
        histos->UserHistogram("Track","nITSshared","#shared ITS clusters", 7, 0, 7, AliDielectronVarManager::kNclsSITS);
        histos->UserHistogram("Track","fracITSshared","frac. shared ITS clusters", 120, 0,  1.2, AliDielectronVarManager::kNclsSITS);

        //TPC
        histos->UserHistogram("Track","TPCnCls",";TPC number clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNclsTPC);
        histos->UserHistogram("Track","TPCchi2",";TPC chi2/Cl;#tracks",100,0.,10.,AliDielectronVarManager::kTPCchi2Cl);
        histos->UserHistogram("Track","nTPCshared","#shared TPC clusters", 170, 0, 170, AliDielectronVarManager::kNclsSTPC);
        histos->UserHistogram("Track","NclsSFracTPC",";TPC fraction of shared clusters;#tracks",200,0,10.,AliDielectronVarManager::kNclsSFracTPC);
        histos->UserHistogram("Track","TPCnCrossed",";TPC findable clusters;#tracks",170,-0.5,169.5,AliDielectronVarManager::kNFclsTPCr);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindable",";TPC crossed rows over findable clusters;#tracks",240,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
        //Repeats to check for binning issues
        histos->UserHistogram("Track","TPCcrossedRowsOverFindableThin",";TPC crossed rows over findable clusters;#tracks",360,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindableWide",";TPC crossed rows over findable clusters;#tracks",120,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindableMed",";TPC crossed rows over findable clusters;#tracks",200,0.,1.2,AliDielectronVarManager::kNFclsTPCfCross);

        //Quality
        histos->UserHistogram("Track","TPCclsDiff",";TPC cluster difference;#tracks",200,0,20.,AliDielectronVarManager::kTPCclsDiff);
        histos->UserHistogram("Track","TPCsignalN",";TPC number PID clusters;#tracks",160,-0.5,159.5,AliDielectronVarManager::kTPCsignalN);



        // ITS
        histos->UserHistogram("Track","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                            GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
        histos->UserHistogram("Track","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                            GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
        histos->UserHistogram("Track","ITSnSigmaEle_Eta",";Eta;n#sigma_{ele}^{ITS}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaEle);
        histos->UserHistogram("Track","ITSnSigmaEle_Phi",";Phi;n#sigma_{ele}^{ITS}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaEle);
        histos->UserHistogram("Track","ITSnSigmaEle_Cent", ";Centrality;n#sigma_{ele}^{ITS}",
                            GetVector(kCent), GetVector(kSigmaEle), AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kITSnSigmaEle);
        histos->UserHistogram("Track","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
        histos->UserHistogram("Track","ITSnSigmaPio_Eta",";p (GeV/c);n#sigma_{pion}^{ITS}",
                            GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kITSnSigmaPio);
        histos->UserHistogram("Track","ITSnSigmaPio_Phi",";p (GeV/c);n#sigma_{pion}^{ITS}",
                            GetVector(kPhi2D), GetVector(kSigmaOther), AliDielectronVarManager::kPhi,AliDielectronVarManager::kITSnSigmaPio);
        histos->UserHistogram("Track","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
        histos->UserHistogram("Track","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);
        // TPC
        histos->UserHistogram("Track","TPC_dEdx_P",";p (GeV/c);TPC signal (arb units)",
                            GetVector(kP2D), GetVector(kTPCdEdx), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCsignal);
        histos->UserHistogram("Track","TPCnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{TPC}",
                            GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaEle_Eta",";Eta;n#sigma_{ele}^{TPC}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaEle_Phi",";Phi;n#sigma_{ele}^{TPC}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaEle_Nacc",";N_{acc}; n#sigma_{ele}^{TPC}",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kNacc,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly",";N_{TPC ref}; n#sigma_{ele}^{TPC}",
                            BinsToVector(100,0.,5000.), GetVector(kSigmaEle), AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaEle_RunNumber",";run;n#sigma_{ele}^{TPC}",
                            GetVector(kRuns), GetVector(kSigmaEle), AliDielectronVarManager::kRunNumber,AliDielectronVarManager::kTPCnSigmaEle);
        histos->UserHistogram("Track","TPCnSigmaPio_P",";p (GeV/c);n#sigma_{ele}^{TPC}",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaPio);
        histos->UserHistogram("Track","TPCnSigmaPio_Eta",";Eta;n#sigma_{ele}^{TPC}",
                            GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaPio);
        histos->UserHistogram("Track","TPCnSigmaPio_Phi",";Phi;n#sigma_{ele}^{TPC}",
                            GetVector(kPhi2D), GetVector(kSigmaOther), AliDielectronVarManager::kPhi,AliDielectronVarManager::kTPCnSigmaPio);

            
        histos->UserHistogram("Track","TPC_dEdx_Eta",";Eta;TPC signal (arb units)",
                              GetVector(kEta2D), GetVector(kTPCdEdx), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal);
        histos->UserHistogram("Track","TPC_dEdx_Eta_P",";Eta;TPC signal (arb units);p (GeV/c)",
                              GetVector(kEta3D), GetVector(kTPCdEdx), GetVector(kP2D),
                              AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kP);
        histos->UserHistogram("Track","TPC_dEdx_P_RunNumber",";p (GeV/c);TPC signal (arb units);run",
                              GetVector(kP2D), GetVector(kTPCdEdx), GetVector(kRuns), 
                              AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal,AliDielectronVarManager::kRunNumber);
        histos->UserHistogram("Track","TPCnSigmaEle_P_dEdx",";p (GeV/c);n#sigma_{ele}^{TPC};TPC signal (arb units)",GetVector(kP2D), GetVector(kSigmaEle), GetVector(kTPCdEdx), AliDielectronVarManager::kP, AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kTPCsignal);
        histos->UserHistogram("Track","TPC nsig electron vs eta",";eta;p_{in} (GeV/c);TPC nSigam elec vs. eta",
                              BinsToVector(100.,0.,1.), GetVector(kSigmaEle),
                              AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle);
        
        histos->UserHistogram("Track","TPCnSigmaEle_Eta_P",";Eta;n#sigma_{ele}^{TPC};p_{in} (GeV/c)",
                                  GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kP2D),
                                  AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kPIn);
        histos->UserHistogram("Track","TPCnSigmaEle_Eta_Nacc",";Eta;n#sigma_{ele}^{TPC};N_{acc}",
                              GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                              AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
        histos->UserHistogram("Track","TPCnSigmaEle_Eta_RefMultTPConly",";Eta;n#sigma_{ele}^{TPC};N_{TPC ref}",
                                  GetVector(kEta3D), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                                  AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRefMultTPConly);
        histos->UserHistogram("Track","TPCnSigmaEle_Eta_RunNumber",";Eta;n#sigma_{ele}^{TPC};run",
                              GetVector(kEta3D), GetVector(kSigmaEle), GetVector(kRuns),
                              AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
        histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_RunNumber",";N_{TPC ref};n#sigma_{ele}^{TPC};run",
                              BinsToVector(100,0.,5000.), GetVector(kSigmaEle), GetVector(kRuns),
                              AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kRunNumber);
        histos->UserHistogram("Track","TPCnSigmaEle_RefMultTPConly_Nacc",";N_{TPC ref};n#sigma_{ele}^{TPC};N_{acc}",
                              BinsToVector(100,0.,5000.), GetVector(kSigmaEle), BinsToVector(100,0.,5000.),
                              AliDielectronVarManager::kRefMultTPConly,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kNacc);
        
        histos->UserHistogram("Track","TPCnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{TPC}",
                              GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaKao);
        histos->UserHistogram("Track","TPCnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{TPC}",
                              GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaPro);
        histos->UserHistogram("Track","TPCnSigmaKao_Eta",";Eta;n#sigma_{kaon}^{TPC}",
                              GetVector(kEta2D), GetVector(kSigmaOther), AliDielectronVarManager::kEta,AliDielectronVarManager::kTPCnSigmaKao);
      
          
            
        // TOF
        histos->UserHistogram("Track","TOFbeta_P",";p (GeV/c);TOF beta",
                            GetVector(kP2D), BinsToVector(120,0.,1.2) ,AliDielectronVarManager::kP,AliDielectronVarManager::kTOFbeta);
        histos->UserHistogram("Track","TOFnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{elec}^{TOF}",
                            GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaEle);
        histos->UserHistogram("Track","TOFnSigmaEle_Eta",";Eta;n#sigma_{elec}^{TOF}",
                            GetVector(kEta2D), GetVector(kSigmaEle), AliDielectronVarManager::kEta,AliDielectronVarManager::kTOFnSigmaEle);
        histos->UserHistogram("Track","TOFnSigmaEle_Phi",";Phi;n#sigma_{elec}^{TOF}",
                            GetVector(kPhi2D), GetVector(kSigmaEle), AliDielectronVarManager::kPhi,AliDielectronVarManager::kTOFnSigmaEle);
        histos->UserHistogram("Track","TOFnSigmaPio_P",";p (GeV/c);TOF number of sigmas Pions",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaPio);
        histos->UserHistogram("Track","TOFnSigmaKao_P",";p (GeV/c);TOF number of sigmas Kaons",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaKao);
        histos->UserHistogram("Track","TOFnSigmaPro_P",";p (GeV/c);TOF number of sigmas Protons",
                            GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kTOFnSigmaPro);
        histos->UserHistogram("Track","TOFnSigmaEle_Cent", ";Centrality;n#sigma_{ele}^{ITS}", 
                            GetVector(kCent), GetVector(kSigmaEle), AliDielectronVarManager::kCentralityNew, AliDielectronVarManager::kTOFnSigmaEle);
          
        // 2D-PID
        histos->UserHistogram("Track","PIn_TPCnSigmaEle_ITSnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};n#sigma_{ele}^{ITS}",
                              50,0.,2.5, 160,-12.,20., 150,-10.,20.,
                            AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kITSnSigmaEle);
        histos->UserHistogram("Track","PIn_TPCnSigmaEle_TOFnSigmaEle",";p_{in} (GeV/c);n#sigma_{ele}^{TPC};TOF number of sigmas Electrons",
                              50,0.,2.5, 160,-12.,20., 50,-5.,5.,
                            AliDielectronVarManager::kP,AliDielectronVarManager::kTPCnSigmaEle,AliDielectronVarManager::kTOFnSigmaEle);


        histos->UserHistogram("Track","TPCcrossedRows_TPCnCls",";TPC number clusters;TPC crossed rows",
                            160,-0.5,159.5, 160,-0.5,159.5, AliDielectronVarManager::kNclsTPC,AliDielectronVarManager::kNFclsTPCr);
        histos->UserHistogram("Track","TPCcrossedRows_P",";Pt [GeV];TPC crossed rows",
                            GetVector(kP2D), BinsToVector(160,-0.5,159.5), AliDielectronVarManager::kP,AliDielectronVarManager::kNFclsTPCr);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindable_P",";P [GeV];TPC crossed rows over findable",
                            GetVector(kP2D), BinsToVector(120,0.,1.2), AliDielectronVarManager::kP,AliDielectronVarManager::kNFclsTPCfCross);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Eta",";Eta;TPC crossed rows over findable",
                            100,-1,1, 120,0.,1.2, AliDielectronVarManager::kEta,AliDielectronVarManager::kNFclsTPCfCross);
        histos->UserHistogram("Track","TPCcrossedRowsOverFindable_Phi",";Phi;TPC crossed rows over findable",
                            120,0.,TMath::TwoPi(), 120,0.,1.2, AliDielectronVarManager::kPhi,AliDielectronVarManager::kNFclsTPCfCross);
    } 
  
    if(doPairing){
        //add histograms to Pair classes
        histos->UserHistogram("Pair","InvMass","",500,0.,5.,AliDielectronVarManager::kM);
        histos->UserHistogram("Pair","PairPt","",160,0.,8., AliDielectronVarManager::kPt);
        histos->UserHistogram("Pair","Rapidity","",200,-2.,2.,AliDielectronVarManager::kY);
        histos->UserHistogram("Pair","OpeningAngle","",240,0.,TMath::Pi(),AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","PhiV","", GetVector(kPhiV), AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","dXY abs (sqrt)",""    ,200 ,0,2.0 , AliDielectronVarManager::kPairDCAabsXY); 
        histos->UserHistogram("Pair","dZ abs (sqrt)",""     ,500 ,0,5.0 , AliDielectronVarManager::kPairDCAabsZ); 
        histos->UserHistogram("Pair","dXY sigma (sqrt)",""  ,2000,0,20.0, AliDielectronVarManager::kPairDCAsigXY); 
        histos->UserHistogram("Pair","dZ sigma (sqrt)",""   ,2000,0,20.0, AliDielectronVarManager::kPairDCAsigZ); 
        histos->UserHistogram("Pair","dXY abs (linear)",""  ,100 ,0,1.0 , AliDielectronVarManager::kPairLinDCAabsXY); 
        histos->UserHistogram("Pair","dZ abs (linear)",""   ,500 ,0,5.0 , AliDielectronVarManager::kPairLinDCAabsZ); 
        histos->UserHistogram("Pair","dXY sigma (linear)","",2000,0,20.0, AliDielectronVarManager::kPairLinDCAsigXY); 
        histos->UserHistogram("Pair","dZ sigma (linear)","" ,2000,0,20.0, AliDielectronVarManager::kPairLinDCAsigZ); 

                //2D and 3D histograms
        histos->UserHistogram("Pair","InvMass_PairPt",";Inv. Mass [GeV];Pair Pt [GeV];#pairs",
                              GetVector(kMee), GetVector(kPtee),
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPt);
        histos->UserHistogram("Pair","Eta_Phi_Pair",";Eta;Phi;#pairs",
                              200,-2.,2, 120,0.,TMath::TwoPi(),
                              AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
        histos->UserHistogram("Pair","InvMass_PairPt_PhivPair",";Inv. Mass [GeV];Pair Pt [GeV];PhiV",
                              GetVector(kMee), GetVector(kPtee), GetVector(kPhiV), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","InvMass_PairPt_OpeningAngle",";Inv. Mass [GeV];Pair Pt [GeV];Opening Angle",
                              GetVector(kMee), GetVector(kPtee), GetVector(kOpAng), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","InvMass_PhivPair_OpeningAngle",";Inv. Mass [GeV];PhiV;Opening Angle",
                              GetVector(kMee500), GetVector(kPhiV), GetVector(kOpAng2), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair, AliDielectronVarManager::kOpeningAngle);

        //opening angle and PhiV
        histos->UserHistogram("Pair","InvMass_OpeningAngle",";Inv. Mass [GeV];Opening Angle;#pairs",
                              GetVector(kMee), GetVector(kOpAng), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","InvMass_PhivPair",";Inv. Mass [GeV];PhiV;#pairs",
                              GetVector(kMee), GetVector(kPhiV), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","PairPt_OpeningAngle",";Pair Pt [GeV];Opening Angle;#pairs",
                              GetVector(kPtee), GetVector(kOpAng), 
                              AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
        histos->UserHistogram("Pair","PairPt_PhivPair",";Pair Pt [GeV];PhiV;#pairs",
                              GetVector(kPtee), GetVector(kPhiV), 
                              AliDielectronVarManager::kPt, AliDielectronVarManager::kPhivPair);
        histos->UserHistogram("Pair","OpeningAngle_PhivPair",";Opening Angle;PhiV;#pairs",
                              GetVector(kOpAng), GetVector(kPhiV), 
                              AliDielectronVarManager::kOpeningAngle, AliDielectronVarManager::kPhivPair);

        //centrality
        histos->UserHistogram("Pair","InvMass_Centrality",";Inv. Mass [GeV];Centrality;#pairs",
                              GetVector(kMee), BinsToVector(102,-1,101), 
                              AliDielectronVarManager::kM, AliDielectronVarManager::kCentrality);
        histos->UserHistogram("Pair","PairPt_Centrality",";Pair Pt [GeV];Centrality;#pairs",
                              GetVector(kPtee), BinsToVector(102,-1,101), 
                              AliDielectronVarManager::kPt, AliDielectronVarManager::kCentrality);
    }//End doMixing histograms

    //add histograms to Track classes
    histos->UserHistogram("Pre","Pt",";Pt [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPt);
    histos->UserHistogram("Pre","Px_pre",";Px [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPx);
    histos->UserHistogram("Pre","Py",";Py [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPy);
    histos->UserHistogram("Pre","Pz",";Pz [GeV];#tracks",200,0,10.,AliDielectronVarManager::kPz);

    histos->UserHistogram("Pre","ITS_dEdx_P",";p (GeV/c);ITS signal (arb units)",
                          GetVector(kP2D), BinsToVector(700,0.,700.), AliDielectronVarManager::kP,AliDielectronVarManager::kITSsignal);
    histos->UserHistogram("Pre","TPC_dEdx_P",";p_{in} (GeV/c);TPC signal (arb units)",
                          GetVector(kP2D), BinsToVector(120,0.,120.), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCsignal);

    histos->UserHistogram("Pre","ITSnSigmaEle_P",";p (GeV/c);n#sigma_{ele}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaEle);
    histos->UserHistogram("Pre","ITSnSigmaPio_P",";p (GeV/c);n#sigma_{pion}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPio);
    histos->UserHistogram("Pre","ITSnSigmaKao_P",";p (GeV/c);n#sigma_{kaon}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaKao);
    histos->UserHistogram("Pre","ITSnSigmaPro_P",";p (GeV/c);n#sigma_{proton}^{ITS}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kP,AliDielectronVarManager::kITSnSigmaPro);

    histos->UserHistogram("Pre","TPCnSigmaEle_P",";p_{in} (GeV/c);n#sigma_{ele}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaEle), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaEle);
    histos->UserHistogram("Pre","TPCnSigmaPio_P",";p_{in} (GeV/c);n#sigma_{pion}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPio);
    histos->UserHistogram("Pre","TPCnSigmaKao_P",";p_{in} (GeV/c);n#sigma_{kaon}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaKao);
    histos->UserHistogram("Pre","TPCnSigmaPro_P",";p_{in} (GeV/c);n#sigma_{proton}^{TPC}",
                          GetVector(kP2D), GetVector(kSigmaOther), AliDielectronVarManager::kPIn,AliDielectronVarManager::kTPCnSigmaPro);

    die->SetHistogramManager(histos);
}



TVectorD *GetVector(Int_t var) 
{
  switch (var) 
  {
    case kPhiV:   return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng:  return AliDielectronHelper::MakeLinBinning(100, 0., TMath::Pi());
    case kOpAng2: return AliDielectronHelper::MakeLinBinning( 50, 0., TMath::Pi()/2.);
    case kPhi2D:  return AliDielectronHelper::MakeLinBinning(100, 0, 2*TMath::Pi());
    case kEta2D:  return AliDielectronHelper::MakeLinBinning(100,-1,1);
    case kEta3D:  return AliDielectronHelper::MakeLinBinning( 50,-1,1);
      
    case kSigmaEle:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-10.,10.);
      else          return AliDielectronHelper::MakeLinBinning( 50, -5., 5.);
    case kSigmaOther:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(100,-20.,20.);
      else          return AliDielectronHelper::MakeLinBinning( 50,-10.,10.);
    case kTPCdEdx:
      if (isQAtask) return AliDielectronHelper::MakeLinBinning(120,  0.,120.);
      else          return AliDielectronHelper::MakeLinBinning( 50, 50.,100.);
      
    case kMee:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.74, 0.78, 0.82, 0.86, 
                                                                   0.90, 0.94, 0.98, 1.02, 1.06, 
                                                                   1.10, 1.30, 1.50, 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 
                                                                   3.10, 3.30, 3.50, 4.00, 4.50, 5.00 
                                                                   ");
    case kMee500: return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 
                                                                   0.10, 0.14, 0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42, 0.46, 
                                                                   0.50 
                                                                   ");
    case kPtee:   return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 
                                                                   3.00, 3.20, 3.40, 3.60, 3.80, 4.00, 4.20, 4.40, 4.60, 4.80, 
                                                                   5.00, 6.00, 7.00, 8.00 
                                                                   ");
    case kP2D:    return AliDielectronHelper::MakeArbitraryBinning("0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 
                                                                   0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 
                                                                   1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 
                                                                   1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 
                                                                   2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80, 
                                                                   4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00 
                                                                   ");
      //2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 
      //2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 
      //3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 
    case kCent: return AliDielectronHelper::MakeLinBinning(10,0.,100.);
    
    //First and lasts bins added for clearer plotting
    case kRuns:   return AliDielectronHelper::MakeArbitraryBinning("265300, 265309, 265332, 265334, 265335, 265336, 265338, 265339, 265342, 265343, 265344, 265377, 265378, 265381, 265383, 265384, 265385, 265387, 265388, 265419, 265420, 265421, 265422, 265424, 265425, 265426, 265427, 265435, 265499, 265500, 265501, 265521, 265525, 265530");

    default: cout << "ERROR: in 'GetVector(...var)' variable for axis range not defined!" << endl;
      break;
  } 
  //if ( var.EqualTo("p_2D"      , kIgnoreCase) ) return AliDielectronHelper::MakeLinBinning(160,0.,8.);
}



TVectorD *BinsToVector(Int_t nbins, Double_t min, Double_t max) {
  return AliDielectronHelper::MakeLinBinning(nbins,min,max);
  //  TVectorD *vec = new TVectorD(nbins+1);
  //
  //  Double_t binwdth = (max-min)/nbins;
  //  for (int i = 0; i < nbins+1; i++) (*vec)[i] = min + i*binwdth;
  //  
  //  return vec;
}



void InitCF(AliDielectron* die, Int_t cutDefinition)
{
  //
  // Setup the CF Manager if needed
  //
  AliDielectronCF *cf=new AliDielectronCF(die->GetName(),die->GetTitle());
  
  //pair variables
  cf->AddVariable(AliDielectronVarManager::kP,100,0.,5.);
  cf->AddVariable(AliDielectronVarManager::kM,200,-0.01,3.99); //20Mev Steps
  cf->AddVariable(AliDielectronVarManager::kPairType,10,0,10);
  
  cf->AddVariable(AliDielectronVarManager::kCentrality,"0.,5.,10.,20.,30.,50.,80.,100.");
  
  //leg variables
  cf->AddVariable(AliDielectronVarManager::kP,160,0.,8.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kITSsignal,350,0.,700.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kTPCsignal,60,0.,120.,kTRUE);
  cf->AddVariable(AliDielectronVarManager::kHaveSameMother,21,-10,10,kTRUE);
  
  //only in this case write MC truth info
  if (MCenabled) { // more elegant: die->GetHasMC() 
    cf->SetStepForMCtruth();
    cf->SetStepsForMCtruthOnly();
    cf->AddVariable(AliDielectronVarManager::kPdgCode,10000,-5000.5,4999.5,kTRUE);
    cf->AddVariable(AliDielectronVarManager::kPdgCodeMother,10000,-5000.5,4999.5,kTRUE);
  }
  
  cf->SetStepsForSignal();
  die->SetCFManagerPair(cf);
}
