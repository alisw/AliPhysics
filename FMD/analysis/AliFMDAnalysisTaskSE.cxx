#include "AliFMDAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "iostream"
#include "AliESDFMD.h"

ClassImp(AliFMDAnalysisTaskSE)

//_____________________________________________________________________
AliFMDAnalysisTaskSE::AliFMDAnalysisTaskSE():
AliAnalysisTaskSE(),
  fListOfHistos(0),
  fSharing("Sharing",kFALSE),
  fDensity("Density",kFALSE),
  fBackground("BackgroundCorrected",kFALSE),
  fDndeta("dNdeta",kFALSE)
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
  fDndeta("dNdeta",kFALSE)
{
 
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
  
  TObjString* vtxString1 = new TObjString();
  
  TList* densitylist = new TList();
  
  TList* bgcorlist = new TList();
  
  fSharing.SetFMDData(fmd);
  fSharing.SetVertex(vertex);
  
  fDensity.SetOutputList(densitylist);
  fDensity.SetInputESDFMD(fmd) ;
  fDensity.SetInputVertex(vertex);
  
  fBackground.SetInputList(densitylist);
  fBackground.SetOutputList(bgcorlist);
  fBackground.SetHitList(fListOfHistos);
  fBackground.SetOutputVertex(vtxString1);
  
  fDndeta.SetInputVertex(vtxString1);
  fDndeta.SetInputList(bgcorlist); 
  fDndeta.SetOutputList(fListOfHistos); 
  
  fSharing.CreateOutputObjects();
  fDensity.CreateOutputObjects();
  fBackground.CreateOutputObjects();
  fDndeta.CreateOutputObjects();
  
  
  
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
  fSharing.SetInputESD(fESD);
  
  fSharing.Exec("");
  if(fSharing.GetEventStatus()) {
    fDensity.Exec("");
    if(fDensity.GetEventStatus()) {
      fBackground.Exec("");  
      AliMCEvent* mcevent = MCEvent();
      fDndeta.SetMCEvent(mcevent);
      fDndeta.Exec("");
    }
  }
  else
    return;
  
  //fListOfHistos = fBackground.GetOutputList();
  
  PostData(1, fListOfHistos);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSE::Terminate(Option_t */*option*/)
{
  fDndeta.Terminate("");

}
//_____________________________________________________________________
//
// EOF
//
