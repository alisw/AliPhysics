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
ClassImp(AliFMDAnalysisTaskSE)

//_____________________________________________________________________
AliFMDAnalysisTaskSE::AliFMDAnalysisTaskSE():
AliAnalysisTaskSE(),
  fListOfHistos(0),
  fSharing("Sharing",kFALSE),
  fDensity("Density",kFALSE),
  fBackground("BackgroundCorrected",kFALSE),
  fDndeta("dNdeta",kFALSE), 
  fBFCorrelation("BFCorrelation",kFALSE), 
  fParams(0)
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
  fParams(0)
{
  SetParams(AliFMDAnaParameters::Instance());
  DefineOutput(1, TList::Class());
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
  
  //AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  
  AliESDEvent* fESD = (AliESDEvent*)InputEvent();
  //std::cout<<fESD->GetBeamEnergy()<<"   "<<fESD->GetBeamType()<<"    "<<fESD->GetCurrentL3()<<std::endl;
  fSharing.SetInputESD(fESD);
  
  fSharing.Exec("");
  if(fSharing.GetEventStatus()) {
    fDensity.Exec("");
    if(fDensity.GetEventStatus()) {
      fBackground.Exec("");  
      fDndeta.Exec("");
      fBFCorrelation.Exec("");
    }
    else return;
  }
  else return;
  
  PostData(1, fListOfHistos);
  
  //fListOfHistos = fBackground.GetOutputList();
  
 
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::Terminate(Option_t */*option*/)
{
  
  TList* outputList = (TList*)GetOutputData(1);
  
  if(outputList) {
    fSharing.SetOutputList(outputList);
    fBackground.SetHitList(outputList);
    fDndeta.SetOutputList(outputList); 
    fBFCorrelation.SetOutputList(outputList); 
    fSharing.Terminate("");
    fBackground.Terminate("");
    fDndeta.Terminate("");
    fBFCorrelation.Terminate("");
    
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
