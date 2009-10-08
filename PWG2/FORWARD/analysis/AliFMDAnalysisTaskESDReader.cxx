 
#include <TROOT.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>

#include "AliFMDAnalysisTaskESDReader.h"
#include "AliAnalysisManager.h"
#include "AliESDFMD.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliFMDAnaParameters.h"

ClassImp(AliFMDAnalysisTaskESDReader)

//_____________________________________________________________________
AliFMDAnalysisTaskESDReader::AliFMDAnalysisTaskESDReader()
: fDebug(0),
  fChain(0x0),
  fESD(0x0),
  fOutputESD(0x0)
{
  // Default constructor
  DefineInput (0, TTree::Class());
  DefineOutput(0, AliESDEvent::Class());
 
}
//_____________________________________________________________________
AliFMDAnalysisTaskESDReader::AliFMDAnalysisTaskESDReader(const char* name):
    AliAnalysisTask(name, "AnalysisTaskFMD"),
    fDebug(0),
    fChain(0x0),
    fESD(0x0),
    fOutputESD(0x0)
{
  DefineInput (0, TTree::Class());
  DefineOutput(0, AliESDEvent::Class());
 
}

//_____________________________________________________________________
void AliFMDAnalysisTaskESDReader::ConnectInputData(Option_t */*option*/)
{
  fChain = (TChain*)GetInputData(0);
  fESD = new AliESDEvent();
  fESD->ReadFromTree(fChain);
  
}
//_____________________________________________________________________

void AliFMDAnalysisTaskESDReader::Exec(Option_t */*option*/)
{
  //  std::cout<<fOutputESD<<std::endl;
  fOutputESD = fESD;
  PostData(0, fOutputESD); 
  
}
//_____________________________________________________________________

