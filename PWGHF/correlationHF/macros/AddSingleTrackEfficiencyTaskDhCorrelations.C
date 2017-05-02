

//Specific Single Track efficiency Addtask for the Dh corr analysis
//Last ModifiedFor:  Track cutfile ON/OFF
//Modified on: Feb16, 2016
//Jitendra


// Particle cuts
const Double_t etamin = -0.8;
const Double_t etamax =  0.8;
const Double_t ptmin = 0.3;
const Double_t ptmax = 24.0;
const Double_t phimin = -2*TMath::Pi();
const Double_t phimax = 2*TMath::Pi();
const Double_t thetamin = 0;
const Double_t thetamax = TMath::Pi();
const Double_t zvtxmin = -10.0;
const Double_t zvtxmax =  10.0;
//
// The following cuts are used on a specific container step, not used in standard efficiency evaluation for D-h analyses (and are active only for ESD analysis)
const Int_t mintrackrefsTPC = 5;
const Int_t mintrackrefsITS = 4;
const Int_t mintrackrefsTOF = 0;
const Int_t mintrackrefsMUON = 0;
const Int_t minclustersTPC = 70;
const Int_t minclustersITS = 2;
const Bool_t TPCRefit = kTRUE;
const Bool_t ITSRefit = kFALSE;
const Bool_t ischarged = kTRUE;
const Int_t  fBit = 0;
//const TString centralityEstimator = "ZNA";

//
// Container settings
//
// Container mutliplicity bins
const Float_t multmin_0_20 = 0;
const Float_t multmax_0_20 = 20;
const Float_t multmin_20_50 = 20;
const Float_t multmax_20_50 = 50;
const Float_t multmin_50_102 = 50;
const Float_t multmax_50_102 = 150;
//  Container Pt bins
Double_t ptmin_0_2   = 0.0;
Double_t ptmax_0_2   = 2.0;
Double_t ptmin_2_6   = 2.0;
Double_t ptmax_2_6   = 6.0;
Double_t ptmin_6_8   = 6.0;
Double_t ptmax_6_8   = 8.0;
Double_t ptmin_8_16  = 8.0;
Double_t ptmax_8_16  = 16.0;
Double_t ptmin_16_24 = 16.0;
Double_t ptmax_16_24 = 24.0;
// Container centrality bins
const Float_t centmin_0_10 = 0.;
const Float_t centmax_0_10 = 10.;
const Float_t centmin_10_60 = 10.;
const Float_t centmax_10_60 = 60.;
const Float_t centmin_60_100 = 60.;
const Float_t centmax_60_100 = 100.;

//PID Threshold
const Float_t thresholdPID = 0.8;
AliCFSingleTrackEfficiencyTask *AddSingleTrackEfficiencyTaskDhCorrelations(const Bool_t readAOD = 0, // Flag to read AOD:1 or ESD:0
                                                                           TString suffix="default", // suffix for the output directory
                                                                           AliPID::EParticleType specie=AliPID::kPion, Int_t pdgcode=0, //particle specie
                                                                           ULong64_t triggerMask=AliVEvent::kAnyINT,
                                                                           Bool_t useCentrality = kFALSE,
                                                                           Int_t configuration=AliCFSingleTrackEfficiencyTask::kFast,
                                                                           Int_t usageOfBayesianPID=AliSingleTrackEffCuts::kNoBayesianPID,
                                                                           TString effName="",
                                                                           TString cutObjName="",
									   TString centralityEstimator="ZNA")
{
    
    Info("AliCFSingleTrackEfficiencyTask","SETUP CONTAINER");
    
    //
    // Setting up the container
    //
    // Variables
    const Int_t nvar = 7; // number of variables on the grid: pt, y, phi, theta, zvtx, multiplicity, centrality
    UInt_t nstep = 8;     // number of container steps
    const UInt_t ipt = 0;
    const UInt_t iy  = 1;
    const UInt_t iphi = 2;
    const UInt_t itheta = 3;
    const UInt_t izvtx = 4;
    const UInt_t imult = 5;
    const UInt_t icent = 6;
    //
    // Containter bining
    //   A1. Bins variation by hand for pt = 48 ?
    const Int_t nbinpt_0_2 = 20;  //bins in pt from 0 to 2 GeV
    const Int_t nbinpt_2_6 = 20;   //bins in pt from 2 to 6 GeV
    const Int_t nbinpt_6_8 = 4;   //bins in pt from 6 to 8 GeV
    const Int_t nbinpt_8_16 = 3;  //bins in pt from 8 to 16 GeV
    const Int_t nbinpt_16_24 = 1; //bins in pt from 16 to 24 GeV
    //   A2. Bins variation by hand for other variables
    const Int_t nbin2 = 18; //bins in eta
    const Int_t nbin3 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 9 : 1; //bins in phi
    const Int_t nbin4 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 9 : 1; //bins in theta
    const Int_t nbin5 = 20; //bins in zvtx
    //   A3. Bins for multiplicity
    const Int_t nbinmult = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 24 : 3;  //bins in multiplicity (total number)
    const Int_t nbinmult_0_20 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 10 : 1; //bins in multiplicity between 0 and 20
    const Int_t nbinmult_20_50 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 10 : 1; //bins in multiplicity between 20 and 50
    const Int_t nbinmult_50_102 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 4 : 1; //bins in multiplicity between 50 and 102
    //  A4. Bins for centrality
    const Int_t nbincent = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 16 : 3;  //bins in centrality
    const Int_t nbincent_0_10 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 2 : 1;  //bins in centrality between 0 and 10
    const Int_t nbincent_10_60 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 10 : 1;  //bins in centrality between 10 and 60
    const Int_t nbincent_60_100 = configuration==AliCFSingleTrackEfficiencyTask::kSlow ? 4 : 1;  //bins in centrality between 60 and 100
    
    //arrays for the number of bins in each dimension
    Int_t iBin[nvar];
    iBin[0]=nbinpt_0_2+nbinpt_2_6+nbinpt_6_8+nbinpt_8_16+nbinpt_16_24;
    iBin[1]=nbin2;
    iBin[2]=nbin3;
    iBin[3]=nbin4;
    iBin[4]=nbin5;
    iBin[5]=nbinmult;
    iBin[6]=nbincent;
    
    //arrays for lower bounds :
    Double_t *binLimpT = new Double_t[iBin[0]+1];
    Double_t *binLim2 = new Double_t[iBin[1]+1];
    Double_t *binLim3 = new Double_t[iBin[2]+1];
    Double_t *binLim4 = new Double_t[iBin[3]+1];
    Double_t *binLim5 = new Double_t[iBin[4]+1];
    Double_t *binLimmult = new Double_t[iBin[5]+1];
    Double_t *binLimcent = new Double_t[iBin[6]+1];
    
    // set the pt bins
    for(Int_t i=0; i<=nbinpt_0_2; i++) binLimpT[i]=(Double_t)ptmin_0_2 + (ptmax_0_2-ptmin_0_2)/nbinpt_0_2*(Double_t)i ;
    for(Int_t i=0; i<=nbinpt_2_6; i++) binLimpT[i+nbinpt_0_2]=(Double_t)ptmin_2_6 + (ptmax_2_6-ptmin_2_6)/nbinpt_2_6*(Double_t)i ;
    for(Int_t i=0; i<=nbinpt_6_8; i++) binLimpT[i+nbinpt_0_2+nbinpt_2_6]=(Double_t)ptmin_6_8 + (ptmax_6_8-ptmin_6_8)/nbinpt_6_8*(Double_t)i ;
    for(Int_t i=0; i<=nbinpt_8_16; i++) binLimpT[i+nbinpt_0_2+nbinpt_2_6+nbinpt_6_8]=(Double_t)ptmin_8_16 + (ptmax_8_16-ptmin_8_16)/nbinpt_8_16*(Double_t)i ;
    for(Int_t i=0; i<=nbinpt_16_24; i++) binLimpT[i+nbinpt_0_2+nbinpt_2_6+nbinpt_6_8+nbinpt_8_16]=(Double_t)ptmin_16_24 + (ptmax_16_24-ptmin_16_24)/nbinpt_16_24*(Double_t)i;
    
    // Other Variables
    for(Int_t i=0; i<=nbin2; i++) binLim2[i]=(Double_t)etamin + (etamax-etamin)/nbin2*(Double_t)i ;
    for(Int_t i=0; i<=nbin3; i++) binLim3[i]=(Double_t)phimin + (phimax-phimin)/nbin3*(Double_t)i ;
    for(Int_t i=0; i<=nbin4; i++) binLim4[i]=(Double_t)thetamin + (thetamax-thetamin)/nbin4*(Double_t)i ;
    for(Int_t i=0; i<=nbin5; i++) binLim5[i]=(Double_t)zvtxmin + (zvtxmax-zvtxmin)/nbin5*(Double_t)i ;
    
    // multiplicity bining..
    for(Int_t i=0; i<=nbinmult_0_20; i++) binLimmult[i]=(Double_t)multmin_0_20 + (multmax_0_20-multmin_0_20)/nbinmult_0_20*(Double_t)i ;
    for(Int_t i=0; i<=nbinmult_20_50; i++) binLimmult[i+nbinmult_0_20]=(Double_t)multmin_20_50 + (multmax_20_50-multmin_20_50)/nbinmult_20_50*(Double_t)i ;
    for(Int_t i=0; i<=nbinmult_50_102; i++) binLimmult[i+nbinmult_0_20+nbinmult_20_50]=(Double_t)multmin_50_102 + (multmax_50_102-multmin_50_102)/nbinmult_50_102*(Double_t)i ;
    
    // centrality bining
    for(Int_t i=0; i<=nbincent_0_10; i++) binLimcent[i]=(Double_t)centmin_0_10 + (centmax_0_10-centmin_0_10)/nbincent_0_10*(Double_t)i;
    for(Int_t i=0; i<=nbincent_10_60; i++) binLimcent[i+nbincent_0_10]=(Double_t)centmin_10_60 + (centmax_10_60-centmin_10_60)/nbincent_10_60*(Double_t)i;
    for(Int_t i=0; i<=nbincent_60_100; i++) binLimcent[i+nbincent_0_10+nbincent_10_60]=(Double_t)centmin_60_100 + (centmax_60_100-centmin_60_100)/nbincent_60_100*(Double_t)i;
    
    // Container
    AliCFContainer* container = new AliCFContainer(Form("container%s",suffix.Data()),"container for tracks",nstep,nvar,iBin);
    container -> SetBinLimits(ipt,binLimpT);    // pt
    container -> SetBinLimits(iy,binLim2);      // eta
    container -> SetBinLimits(iphi,binLim3);    // phi
    container -> SetBinLimits(itheta,binLim4);  // theta
    container -> SetBinLimits(izvtx,binLim5);   // Zvtx
    container -> SetBinLimits(imult,binLimmult);// multiplicity
    container -> SetBinLimits(icent,binLimcent);// centrality
    
    // Variable Titles
    container -> SetVarTitle(ipt,"pt");
    container -> SetVarTitle(iy, "y");
    container -> SetVarTitle(iphi,"phi");
    container -> SetVarTitle(itheta, "theta");
    container -> SetVarTitle(izvtx, "Zvtx");
    container -> SetVarTitle(imult, "Multiplicity");
    container -> SetVarTitle(icent, "Centrality");
    
    // Step Titles
    container -> SetStepTitle(0, " MC Particle with Generated Cuts");
    container -> SetStepTitle(1, " MC Particle with Kine Acceptance Cuts");
    container -> SetStepTitle(2, " MC Particle with Track Ref Acceptance Cuts");
    container -> SetStepTitle(3, " Total Reconstructed  Particle ");
    container -> SetStepTitle(4, " Reco Particle With Kine Acceptance Cuts");
    container -> SetStepTitle(5, " Reco Particle to MC True pt particles ");
    container -> SetStepTitle(6, " Reco Particle With Quality Cuts");
    container -> SetStepTitle(7, " Reco PID With Quality Cuts");
    
    
    // SET TLIST FOR QA HISTOS
    TList* qaList = new TList();
    TObjArray* emptyList = new TObjArray(0);
    
    //CREATE THE INTERFACE TO CORRECTION FRAMEWORK USED IN THE TASK
    printf("CREATE INTERFACE AND CUTS\n");
    AliCFManager* man = new AliCFManager();
    
    man->SetNStepEvent(2);
    man->SetEventContainer(container);
    man->SetEventCutsList(0,emptyList);//evtmcList);
    man->SetEventCutsList(1,emptyList);//evtrecoList);
    
    man->SetParticleContainer(container);
    man->SetParticleCutsList(0,emptyList);//mcGenList);
    man->SetParticleCutsList(1,emptyList);//mcKineList);
    man->SetParticleCutsList(2,emptyList);//mcaccList);
    man->SetParticleCutsList(3,emptyList);//evtrecoPureList);
    man->SetParticleCutsList(4,emptyList);//recKineList);
    man->SetParticleCutsList(5,emptyList);
    man->SetParticleCutsList(6,emptyList);
    man->SetParticleCutsList(7,emptyList);
    
    // Simulated particle & event cuts
    AliSingleTrackEffCuts* cuts = new AliSingleTrackEffCuts();
    cuts->SetPtRange(ptmin,ptmax);
    cuts->SetEtaRange(etamin,etamax);
    cuts->SetIsCharged(ischarged);
    cuts->SetMinVtxContr(1);
    cuts->SetMaxVtxZ(zvtxmax);
    cuts->SetNumberOfClusters(mintrackrefsITS,mintrackrefsTPC,mintrackrefsTOF,mintrackrefsMUON);
    cuts->SetTriggerMask(triggerMask);
    cuts->SetIsAOD(readAOD);
    //
    // Pid selection here
    //
    if(pdgcode>0){
        cuts->SetUsePid(true);
        cuts->SetParticleSpecie(specie);
        cuts->SetPdgCode(pdgcode);
        //
        const Int_t nlims=1;
        Float_t plims[nlims+1]={0.,999.}; //TPC limits in momentum [GeV/c]
        Float_t sigmas[nlims]={3.};
        cuts->SetUseTPCPid();
        cuts->SetTPCSigmaPtBins(nlims,plims,sigmas);
        cuts->SetMaximumPTPC(4.);
        //
        const Int_t nlims2=1;
        Float_t plims2[nlims2+1]={0.,999.}; //TPC limits in momentum [GeV/c]
        Float_t sigmas2[nlims2]={3.};
        cuts->SetUseTOFPid();
        cuts->SetTOFSigmaPtBins(nlims2,plims2,sigmas2);
        cuts->SetMaximumPTOF(4.);
        
        if(usageOfBayesianPID>0) {
            cuts->SetUseCombinPID(usageOfBayesianPID);
            if(usageOfBayesianPID==AliSingleTrackEffCuts::kThresholdBayesianProb)
                cuts->SetPIDThreshold(thresholdPID);
        }
        
    }
    
    //  Track Quality cuts via ESD track cuts
    AliESDtrackCuts* QualityCuts = new AliESDtrackCuts();
    if(!effName.Contains(".root")){
        cout<<"*** Track efficiency input file not found/set! Using hardcoded selection... ***"<<endl;
        QualityCuts->SetRequireSigmaToVertex(kFALSE);
        QualityCuts->SetRequireITSRefit(kFALSE);
        QualityCuts->SetRequireTPCRefit(kTRUE);
        QualityCuts->SetMinNClustersITS(3);
        QualityCuts->SetMinNClustersTPC(70);
        QualityCuts->SetMaxChi2PerClusterTPC(4);
        QualityCuts->SetEtaRange(etamin,etamax);
        QualityCuts->SetPtRange(ptmin,ptmax);
        QualityCuts->SetMaxDCAToVertexZ(1.0);
        QualityCuts->SetMaxDCAToVertexXY(0.25);
    }else {
        cout<<"Track efficiency input file is set !"<<endl;
        TFile* fileeff=TFile::Open(effName.Data());
        if(!fileeff->IsOpen())cout<<"*** Not able to open file... ***"<<endl;
        else cout<<".... File Quality Seems Ok !"<<endl;
        AliHFAssociatedTrackCuts* HFTrackcuts = (AliHFAssociatedTrackCuts*)fileeff->Get(cutObjName.Data());
        if(!HFTrackcuts) {::Error("AddTask: cutFile","Wrong objectname or file! Exiting..."); return NULL;}
        QualityCuts = (AliESDtrackCuts*)HFTrackcuts->GetESDTrackCuts();
    }
    
    
    //CREATE THE TASK
    printf("CREATE CF Single track task\n");
    
    AliCFSingleTrackEfficiencyTask *task = new AliCFSingleTrackEfficiencyTask("AliCFSingleTrackEfficiencyTask",QualityCuts,cuts);
    if(readAOD) task->SetFilterBit(kTRUE);
    else task->SetFilterBit(kFALSE);
    task->SetFilterType(fBit);
    //  task->SelectCollisionCandidates(triggerMask);//AliVEvent::kMB);
    if(useCentrality) task->SetUseCentrality(useCentrality,centralityEstimator);
    task->SetConfiguration(configuration);
    task->SetCFManager(man); //here is set the CF manager
    
    //
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTask", "No analysis manager to connect to.");
        return NULL;
    }
    
    // This task requires an ESD or AOD input handler and an AOD output handler.
    // Check this using the analysis manager.
    //===============================================================================
    TString type = mgr->GetInputEventHandler()->GetDataType();
    if (!type.Contains("ESD") && !type.Contains("AOD")) {
        ::Error("AddSingleTrackEfficiencyTask", "AliCFSingleTrackEfficiency task needs the manager to have an ESD or AOD input handler.");
        return NULL;
    }
    
    mgr->AddTask(task);
    printf(" Create the output container\n");
    
    //
    // Create and connect containers for input/output
    //
    // ----- output data -----
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    TString input1name="cchain0";
    TString output2name="HistEventsProcessed", output3name="container",output4name="list",output5name="ESDtrackCuts",output6name="MCtrackCuts";
    outputfile += ":PWGPP_CFSingleTrack";
    //  outputfile += suffix;
    output2name += suffix;
    output3name += suffix;
    output4name += suffix;
    output5name += suffix;
    output6name += suffix;
    
    
    // ------ input data ------
    AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
    // ----- output data -----
    // output TH1I for event counting
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(output2name, TH1I::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    // output Correction Framework Container (for acceptance & efficiency calculations)
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(output3name, AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    // output QA histograms
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(output4name, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    // output ESD track cuts for book keeping
    AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(output5name, AliESDtrackCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    // output event and particle selection cuts for book keeping
    AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(output6name, AliSingleTrackEffCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput1);
    mgr->ConnectOutput(task,2,coutput2);
    mgr->ConnectOutput(task,3,coutput3);
    mgr->ConnectOutput(task,4,coutput4);
    mgr->ConnectOutput(task,5,coutput5);
    
    return task;
}
