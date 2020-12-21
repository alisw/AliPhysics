///////////////////////////////////////////////////////////////////
//                                                               //
// AddMyTask                                                     //
// Author: Daichi Kawana, Univ. of Tsukuba                       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

 AliAnalysisTaskCaloHFEpPbRun2* AddTaskCaloHFEpPbRun2(
                                  TString name = "name",
                                  TString infoname = "infoname",
                                  TString period,
                                  Bool_t flagMC,
                                  Bool_t flagEMCalCorrection,
                                  Bool_t flagEMCal,
                                  Bool_t flagDCal,
                                  Bool_t flagEG1,
                                  Bool_t flagEG2,
                                  Bool_t flagDG1,
                                  Bool_t flagDG2,
                                  Double_t centmin,
                                  Double_t centmax,
                                  Double_t TrackEtaLow,
                                  Double_t TrackEtaHigh,
                                  Int_t NTPCClust,
                                  Int_t NITSClust,
                                  Int_t NCrossedRow,
                                  Double_t DCAxy,
                                  Double_t DCAz,
                                  Double_t TrackMatchPhi,
                                  Double_t TrackMatchEta,
                                  Double_t NsigmaLow,
                                  Double_t NsigmaHigh,
                                  Double_t M20Low,
                                  Double_t M20High,
                                  Double_t EopLow,
                                  Double_t EopHigh,
                                  Double_t AssoMinpT,
                                  Int_t AssoNTPCClust,
                                  Double_t MassCut
                                  )
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file

    // now we create an instance of your task
    AliAnalysisTaskCaloHFEpPbRun2* task = new AliAnalysisTaskCaloHFEpPbRun2(name.Data());
    if(!task) return 0x0;
    task -> SetMC(flagMC);
    task -> SetEMCalCorrection(flagEMCalCorrection);
    task -> SetRunPeriod(period);
    if(!flagEG2 && !flagEG1 && !flagDG2 && !flagDG1) task -> SelectCollisionCandidates(AliVEvent::kINT7);
    else task -> SelectCollisionCandidates(AliVEvent::kEMCEGA);
    task -> SetClusterTypeEMC(flagEMCal);
    task -> SetClusterTypeDCAL(flagDCal);
    task -> SetEG1(flagEG1);
    task -> SetEG2(flagEG2);
    task -> SetDG1(flagDG1);
    task -> SetDG2(flagDG2);
    task -> SetCentrality(centmin,centmax);

    // TString PathEffHad = "~/cernbox/Analysis/HFE/p-Pb_8TeV/rootfile/EleHadCorrlpPbRun2/EleHadCorrl/MC/LHC17i5b2/180702/EHCorrlOutputStage1.root";
    TString PathEffHad = "alien:///alice/cern.ch/user/d/dkawana/Efficiency/EffHadron/EHCorrlOutputStage1.root";
    TFile *file = TFile::Open(PathEffHad.Data());
    TDirectory *dir = (TDirectory*)file->Get("Dir_kINT7_integ");
    TH1D *hist = (TH1D*)dir->Get("fHadRecoEff");
    if(hist) task -> SetEffHadron(hist);
    else cout << "No hadron efficiency file!" << endl;
    //#########################//
    //Systematic uncertainties //
    //#########################//
    task -> SetTrackEta(TrackEtaLow,TrackEtaHigh);
    task -> SetTrackClust(NTPCClust,NITSClust,NCrossedRow);
    task -> SetDCA(DCAxy,DCAz);
    task -> SetTrackMatch(TrackMatchPhi,TrackMatchEta);
    task -> SetNsigma(NsigmaLow,NsigmaHigh);
    task -> SetM20(M20Low,M20High);
    task -> SetEop(EopLow,EopHigh);
    task -> SetMassPara(AssoMinpT,AssoNTPCClust,MassCut);
    // add your task to the manager
    mgr->AddTask(task);

    TString containerName = mgr->GetCommonFileName();
    containerName += ":PWGHF_hfeCalpA";
    TString SubcontainerName = Form("hfeCalpA");
    SubcontainerName += name;
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(SubcontainerName, TList::Class(),AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);

    /*
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

    // same for the output
    TString Outputname = "MyOutputContainer_";
    Outputname += infoname.Data();
    mgr->ConnectOutput(task,1,mgr->CreateContainer(Outputname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    //mgr->ConnectOutput(task,1,mgr->CreateContainer("MyOutputContainer", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    */

    return task;
}
