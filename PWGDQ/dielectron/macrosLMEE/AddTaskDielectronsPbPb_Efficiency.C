//____________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask *AddTaskDielectronsPbPb_Efficiency (
                                                    Int_t    iMCprod = 0,
                                                    Int_t    iset = 0,
                                                    const char *centrality = "Central",
                                                    Double_t CentralityMin =  0.0,
                                                    Double_t CentralityMax = 10.0,
                                                    const char *magfield = "up",
                                                    const char *file_DetectorResponseMatrices = "DetectorResponseMatrices.root",
                                                    const char *file_WeightsPtDistributions =   "WeightsPtDistributions.root",
                                                    const char *file_CentralityBins =   "CentralityBins.root",
                                                    Double_t DCAxy_param0 = 0.00515869,
                                                    Double_t DCAxy_param1 = 0.01016680,
                                                    Double_t DCAxy_param2 = 1.34489,
                                                    Double_t DCAz_max     = 0.1,
                                                    Int_t    ITS_minNcls     = 4,
                                                    Int_t    TPC_minNcls     = 80,
                                                    Int_t    TPC_nClsdEdx    = 50,
                                                    Int_t    TPC_minCr       = 100,
                                                    Double_t MinCrOverFindableCls = 0.8,
                                                    Double_t MaxGoldenChi2 = 36,
                                                    Double_t MaxTPCchi2 = 4,
                                                    Double_t MaxITSchi2 = 3,
                                                    Double_t MaxFracSharedCls = 0.4,
                                                    Double_t nsigmaTOF_max = 3.0,
                                                    Double_t nsigmaITS_max = 1.0,
                                                    Double_t nsigmaTPC_min = -3.0,
                                                    Double_t nsigmaTPC_max = 3.0,
                                                    Double_t massLim = 0.02,
                                                    Double_t phivLim = 60.0
                                               )  {
    

    //Load Libraries
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    
    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        printf("ERROR: Analysis Manager Not Found!!!\n");
        return NULL;
    }
    
    
    //Retrieve Input Event Handler
    if (!mgr->GetInputEventHandler()) {
        printf("ERROR: No Input Event Handler!!!\n");
        return NULL;
    }
    
    
    //Detector Response Matrices
    TFile *inputDetRespMatrix = TFile::Open (file_DetectorResponseMatrices);
    TH2F *fDetMatrixP       = (TH2F*) inputDetRespMatrix->Get ("fHistoDetResponseMatrix_Momentum");
    TH2F *fDetMatrixTheta   = (TH2F*) inputDetRespMatrix->Get ("fHistoDetResponseMatrix_Theta");
    TH2F *fDetMatrixPhiElec = (TH2F*) inputDetRespMatrix->Get ("fHistoDetResponseMatrix_Phi_Electrons");
    TH2F *fDetMatrixPhiPos  = (TH2F*) inputDetRespMatrix->Get ("fHistoDetResponseMatrix_Phi_Positrons");
    
    //Weights Pt Distributions
    TFile *inputWeights = TFile::Open (file_WeightsPtDistributions);
    TH2F *fWeight0 = (TH2F*) inputWeights->Get ("fHistoWeights(0)");
    TH2F *fWeight1 = (TH2F*) inputWeights->Get ("fHistoWeights(1)");
    TH2F *fWeight2 = (TH2F*) inputWeights->Get ("fHistoWeights(2)");
    TH2F *fWeight3 = (TH2F*) inputWeights->Get ("fHistoWeights(3)");
    TH2F *fWeight4 = (TH2F*) inputWeights->Get ("fHistoWeights(4)");
    TH2F *fWeight5 = (TH2F*) inputWeights->Get ("fHistoWeights(5)");
    
    //Centrality Bins
    TFile *inputCentrBins = TFile::Open (file_CentralityBins);
    TH1F *fHistoCentralityBins = (TH1F*) inputCentrBins->Get ("fHistoCentralityBins");

    
    
    //Task Name,InputBox & Output File
    const char * TaskName =   Form("DielectronTask_Efficiency_MCprod%d_Set%d_%s_MagField_%s",iMCprod,iset,centrality,magfield);
    const char * InputBox =   Form("Input_%d_%s_%s",iset,centrality,magfield);
    const char * OutputFile = Form("Input_DielectronsPbPb_Efficiency_MCProd%d_Set%d_%s_MagField_%s.root",iMCprod,iset,centrality,magfield);
    
    
    //Input Container
    AliAnalysisDataContainer *input = mgr->GetCommonInputContainer();
    
    //Analysis Task
    AliAnalysisTask *task = new AliAnalysisTaskDielectronsPbPb_Efficiency (TaskName);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::SetCentralityRange (CentralityMin,CentralityMax);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::SetDCAparameters   (DCAxy_param0,DCAxy_param1,DCAxy_param2,DCAz_max);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::SetTrackCuts       (ITS_minNcls,TPC_minNcls,TPC_nClsdEdx,TPC_minCr,MinCrOverFindableCls,MaxGoldenChi2,MaxTPCchi2,MaxITSchi2,MaxFracSharedCls);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::SetPIDCuts         (nsigmaTOF_max,nsigmaITS_max,nsigmaTPC_min,nsigmaTPC_max);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::SetPrefilterCuts   (massLim,phivLim);

    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::GetDetectorResponseMatrices (fDetMatrixP,fDetMatrixTheta,fDetMatrixPhiElec,fDetMatrixPhiPos);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::GetWeightsPtDistributions   (fWeight0,fWeight1,fWeight2,fWeight3,fWeight4,fWeight5);
    task -> AliAnalysisTaskDielectronsPbPb_Efficiency::GetCentralityBins           (fHistoCentralityBins);
    mgr -> AddTask(task);
    AliAnalysisDataContainer *output = mgr->CreateContainer(InputBox,TList::Class(),AliAnalysisManager::kOutputContainer,OutputFile);
    mgr -> ConnectInput  (task,0,input);
    mgr -> ConnectOutput (task,1,output);
    
    return task;
}
//____________________________________________________________________________________________________________________________________________________________________________________________________________

