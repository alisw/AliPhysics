#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisManager.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "AliAnalysisTaskDeuteronCoalescence.h"
#include <TString.h>
#include <TList.h>
#include <TTree.h>
#endif

//_________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronCoalescence *addTaskDeuteronsCoalescence (Double_t Nch_Transv=7.295, Int_t iShift=0)  {
    
    //Get Input Files
    TFile *inputfileWeights     = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/reshaping_protons_pythia/reshaping_protons.root");
    TFile *inputfileWeightsDiff = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/reshaping_protons_pythia/data_to_pythia_differential.root");
    TFile *inputfileDeutWF      = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/deuteron_wave_function/deuteron_wave_function.root");
    TFile *inputfileSource      = TFile::Open ("alien:///alice/cern.ch/user/a/alcaliva/source_radius_pp/source_radius_pp_collisions.root");
       
    //Get Histograms
    TH1D  *hWeight_toward = (TH1D*) inputfileWeightsDiff -> Get("hWeight_toward");
    TH1D  *hWeight_transv = (TH1D*) inputfileWeightsDiff -> Get("hWeight_transv");
    TH1D  *hWeight_away   = (TH1D*) inputfileWeightsDiff -> Get("hWeight_away");
    TH1D  *hProtWeights   = (TH1D*) inputfileWeights     -> Get("hDataToPythia");
    TF1   *fProtWeights   = (TF1*)  inputfileWeights     -> Get("protonWeight");
    TH1F  *hSourceR0      = (TH1F*) inputfileSource      -> Get("hRzeroMt");
    TF1   *fDeuteronWF    = (TF1*)  inputfileDeutWF      -> Get(Form("deuteronWF[%d]",iShift));

    //Get Analysis Manager
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;
    
    //Get Input Event Handler
    if (!mgr->GetInputEventHandler()) return 0x0;
    
    //File Name
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":deuteron_coalescence";
    
    //Analysis Task
    AliAnalysisTaskDeuteronCoalescence *task = new AliAnalysisTaskDeuteronCoalescence ("task_deuteron_coalescence");
    task -> AliAnalysisTaskDeuteronCoalescence::SetAverageTransverseMultiplicity(Nch_Transv);
    task -> AliAnalysisTaskDeuteronCoalescence::SetReshapingProtons (hProtWeights,fProtWeights);
    task -> AliAnalysisTaskDeuteronCoalescence::SetReshapingProtonsDiff (hWeight_toward,hWeight_transv,hWeight_away);
    task -> AliAnalysisTaskDeuteronCoalescence::SetDeuteronWaveFunc (fDeuteronWF);
    task -> AliAnalysisTaskDeuteronCoalescence::SetSourceSizeRadius (hSourceR0);
    mgr -> AddTask(task);
    mgr -> ConnectInput (task,0,mgr->GetCommonInputContainer());
    mgr -> ConnectOutput(task,1,mgr->CreateContainer("Results", TList::Class(),AliAnalysisManager::kOutputContainer,fileName.Data()));
    mgr -> ConnectOutput(task,2,mgr->CreateContainer("QAPlots",TList::Class(),AliAnalysisManager::kOutputContainer,fileName.Data()));
    
    return task;
}
//_________________________________________________________________________________________________________________________________________________

