 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>

#include "AliFMDAnalysisTaskSharing.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliFMDAnaParameters.h"

ClassImp(AliFMDAnalysisTaskSharing)

//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing()
: fDebug(0),
  fESD(0x0),
  fOutputESD(0x0),
  foutputESDFMD(0x0)
{
  // Default constructor
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDEvent::Class());
}
//_____________________________________________________________________
AliFMDAnalysisTaskSharing::AliFMDAnalysisTaskSharing(const char* name):
    AliAnalysisTask(name, "AnalysisTaskFMD"),
    fDebug(0),
    fESD(0x0),
    fOutputESD(0x0),
    foutputESDFMD(0x0)
{
  DefineInput (0, AliESDEvent::Class());
  DefineOutput(0, AliESDEvent::Class());
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::CreateOutputObjects()
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliESDFMD* outputESDFMD      = 0; 
  Int_t nvtxbins = pars->GetNvtxBins();
  
  fOutputESD    = new AliESDEvent();
  fOutputESD->CreateStdContent();
  
  foutputESDFMD = new AliESDFMD();
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::ConnectInputData(Option_t */*option*/)
{
  fESD = (AliESDEvent*)GetInputData(0);
}
//_____________________________________________________________________
void AliFMDAnalysisTaskSharing::Exec(Option_t */*option*/)
{
  AliFMDAnaParameters* pars = AliFMDAnaParameters::Instance();
  
  AliESD* old = fESD->GetAliESDOld();
  if (old) {
    fESD->CopyFromOldESD();
  }
  
  foutputESDFMD->Clear();
  
  fOutputESD->SetPrimaryVertexSPD(fESD->GetPrimaryVertexSPD());
  
  AliESDFMD* fmd = fESD->GetFMDData();
  
  if (!fmd) return;
  
  for(UShort_t det=1;det<=3;det++) {
    Int_t nRings = (det==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for(UShort_t sec =0; sec < nsec;  sec++)  {
	for(UShort_t strip = 0; strip < nstr; strip++) {
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,0.);
	  Float_t mult = fmd->Multiplicity(det,ring,sec,strip);
	  if(mult == AliESDFMD::kInvalidMult) continue;
	  //Sharing algorithm goes here
	  if(mult < 0.75) continue;
	  
	  foutputESDFMD->SetMultiplicity(det,ring,sec,strip,1.);
	  foutputESDFMD->SetEta(det,ring,sec,strip,fmd->Eta(det,ring,sec,strip));
	  
	}
      }
    }
  }
  fOutputESD->SetFMDData(foutputESDFMD);
    
  PostData(0, fOutputESD); 
  
}
//_____________________________________________________________________
//
// EOF
//
