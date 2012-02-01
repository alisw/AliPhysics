#include "AliFMDAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "iostream"
#include "AliESDFMD.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliFMDAnaParameters.h"
#include "AliLog.h"
#include "AliFMDDndeta.h"
#include "TDirectory.h"
#include "TProfile2D.h"
ClassImp(AliFMDAnalysisTaskSE)
//
// This task controls the running of the FMD analysis. The current version is made for 
// dN/deta analysis but multiplicity and correlation tasks will be added here as well.
//
// To get the dN/detadphi per selected event as a TH2F* object connect to the 
// output list of this task in the analysis framework and do
//
// TH2F* hFMDdNdetadphi = (TH2F*)list->FindObject("dNdetadphiHistogramTrVtx");
//_____________________________________________________________________
AliFMDAnalysisTaskSE::AliFMDAnalysisTaskSE():
AliAnalysisTaskSE(),
  fListOfHistos(0),
  fSharing("Sharing",kFALSE),
  fDensity("Density",kFALSE),
  fBackground("BackgroundCorrected",kFALSE),
  fDndeta("dNdeta",kFALSE), 
  fBFCorrelation("BFCorrelation",kFALSE), 
  fParams(0),
  fFirstEvent(kTRUE),
  fCentralityLow(0),
  fCentralityHigh(100)
{
  // Default constructor
}
//_____________________________________________________________________
AliFMDAnalysisTaskSE::AliFMDAnalysisTaskSE(const char* name):
  AliAnalysisTaskSE(name),
  fListOfHistos(0),
  fSharing("Sharing",kFALSE),
  fDensity("Density",kFALSE),
  fBackground("BackgroundCorrected",kFALSE),
  fDndeta("dNdeta",kFALSE), 
  fBFCorrelation("BFCorrelation",kFALSE), 
  fParams(0),
  fFirstEvent(kTRUE),
  fCentralityLow(0),
  fCentralityHigh(100)
{
  SetParams(AliFMDAnaParameters::Instance());
  DefineOutput(1, TList::Class());
  // DefineOutput(2, TH2F::Class());
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::UserCreateOutputObjects()
{
// Create the output containers
//
  
  
  fListOfHistos = new TList();
  
  AliESDFMD* fmd = new AliESDFMD();
  AliESDVertex* vertex = new AliESDVertex();
  
  TList* densitylist = new TList();
  
  TList* bgcorlist = new TList();
    
  fSharing.SetFMDData(fmd);
  fSharing.SetVertex(vertex);
  fSharing.SetOutputList(fListOfHistos);
  
  fDensity.Init();
  fDensity.SetOutputList(densitylist);
  fDensity.SetInputESDFMD(fmd) ;
  fDensity.SetInputVertex(vertex);
  
  fBackground.SetInputList(densitylist);
  fBackground.SetOutputList(bgcorlist);
  fBackground.SetHitList(fListOfHistos);

  fDndeta.SetInputList(bgcorlist); 
  fDndeta.SetOutputList(fListOfHistos); 
  fBFCorrelation.SetInputList(bgcorlist); 
  fBFCorrelation.SetOutputList(fListOfHistos); 
  
  fSharing.CreateOutputObjects();
  fDensity.CreateOutputObjects();
  fBackground.CreateOutputObjects();
  fDndeta.CreateOutputObjects();
  fBFCorrelation.CreateOutputObjects();
 
  
  PostData(1, fListOfHistos);
  
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::Init()
{
  std::cout<<"Init"<<std::endl;
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  AliESDEvent* fESD = (AliESDEvent*)InputEvent();
  
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  // Centrality selection - work in progress
  Float_t centrality = 1;
  if( centrality < fCentralityLow || centrality > fCentralityHigh )  return;
  
  //End of centrality selection
  
  if(fFirstEvent) {
    pars->SetParametersFromESD(fESD);
    pars->PrintStatus();
    fFirstEvent = kFALSE;
  }
  
  pars->SetTriggerStatus(fESD);
  fSharing.SetInputESD(fESD);
  
  fSharing.Exec("");
  if(fSharing.GetEventStatus()) {
    fDensity.Exec("");
    if(fDensity.GetEventStatus()) {
      fBackground.Exec("");  
      if(pars->GetRunDndeta())        fDndeta.Exec("");
      if(pars->GetRunBFCorrelation()) fBFCorrelation.Exec("");
    }
    else return;
  }
  else return;
  
 
  PostData(1, fListOfHistos);
     
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::Terminate(Option_t */*option*/)
{
  
  TList* outputList = (TList*)GetOutputData(1);
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  if(outputList) {
    fSharing.SetOutputList(outputList);
    fBackground.SetHitList(outputList);
    fDndeta.SetOutputList(outputList); 
    //fBFCorrelation.SetOutputList(outputList); 
    fSharing.Terminate("");
    fBackground.Terminate("");
    if(fSharing.GetVtxEfficiencyFromData() > 0)
      fDndeta.SetVtxEfficiency(fSharing.GetVtxEfficiencyFromData());
    else
      fDndeta.SetVtxEfficiency(pars->GetVtxSelectionEffFromMC());
    
    AliInfo(Form("Vertex efficiencies:  NSD_data=%f, INEL_data=%f, INEL_mc=%f", 
		 fSharing.GetNSDVtxEfficiencyFromData(),
		 fSharing.GetVtxEfficiencyFromData(),
		 pars->GetVtxSelectionEffFromMC()));
    
    if(fSharing.GetNSDVtxEfficiencyFromData() > 0)
      fDndeta.SetVtxEfficiencyNSD(fSharing.GetNSDVtxEfficiencyFromData());
    else
      fDndeta.SetVtxEfficiencyNSD(pars->GetVtxSelectionEffFromMC());
    
    fDndeta.Terminate("");
    //fBFCorrelation.Terminate("");
    
    AliFMDDndeta t;
    t.SetNbinsToCut(2);
    t.Init(outputList);
    t.GenerateMult(AliFMDDndeta::kMult);
    
    TList* dNdetalist = t.GetMultList(AliFMDDndeta::kMult);
    TList* cloneList = (TList*)dNdetalist->Clone("dNdeta");
    cloneList->SetName("dNdeta");
    outputList->Add(cloneList);
    
    t.GenerateMult(AliFMDDndeta::kMultTrVtx);
    TList* dNdetalist2 = t.GetMultList(AliFMDDndeta::kMultTrVtx);
    TList* cloneList2 = (TList*)dNdetalist2->Clone("dNdetaTrVtx");
    cloneList2->SetName("dNdetaTrVtx");
    outputList->Add(cloneList2);
  
    t.GenerateMult(AliFMDDndeta::kHits);
    TList* dNdetalist3 = t.GetMultList(AliFMDDndeta::kHits);
    TList* cloneList3 = (TList*)dNdetalist3->Clone("Hits");
    cloneList3->SetName("Hits");
    outputList->Add(cloneList3);
    
    t.GenerateMult(AliFMDDndeta::kHitsTrVtx);
    TList* dNdetalist4 = t.GetMultList(AliFMDDndeta::kHits);
    TList* cloneList4 = (TList*)dNdetalist4->Clone("HitsTrVtx");
    cloneList4->SetName("HitsTrVtx");
    outputList->Add(cloneList4);
    
    t.GenerateMult(AliFMDDndeta::kMultNSD);
    TList* dNdetalist5 = t.GetMultList(AliFMDDndeta::kMultNSD);
    TList* cloneList5 = (TList*)dNdetalist5->Clone("MultNSD");
    cloneList5->SetName("MultNSD");
    outputList->Add(cloneList5);
    
    // TFile file("fmd_ana_histos_tmp.root","RECREATE");
    //  fListOfHistos->Write();
    // file.Close();
  }
  else
    AliWarning("no merged output from manager");
  
  
}

//_____________________________________________________________________
void AliFMDAnalysisTaskSE::Print(Option_t* option) const
{
  AliInfo(Form("FMD Single Event Analysis Task\n"
	       "Parameters set to %p", fParams));
  TString opt(option);
  opt.ToLower();
  if (opt.Contains("s")) { 
    fSharing.Print(option);     
    fDensity.Print(option);     
    fBackground.Print(option);  
    fDndeta.Print(option); 
    fBFCorrelation.Print(option); 
  }
  if (opt.Contains("p") && fParams) 
    fParams->Print(option);      
}

//_____________________________________________________________________
//
// EOF
//
// EOF
