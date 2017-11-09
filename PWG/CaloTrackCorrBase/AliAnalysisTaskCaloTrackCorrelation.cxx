/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <cstdlib>

// --- Root ---
#include <TROOT.h>
#include <TInterpreter.h>
#include <TClonesArray.h>
//#include <Riostream.h>
//#include <TObjectTable.h>

// --- Analysis ---
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "AliCaloTrackReader.h"
#include "AliPDG.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskCaloTrackCorrelation) ;
/// \endcond

//________________________________________________________________________
/// Default constructor.
//________________________________________________________________________
AliAnalysisTaskCaloTrackCorrelation::AliAnalysisTaskCaloTrackCorrelation() :
  AliAnalysisTaskSE(),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), 
  fCuts(0x0),
  fFirstEvent(0),
  fLastEvent(0),
  fStoreEventSummary(0)
{
}

//________________________________________________________________________________________
/// Default constructor.
//________________________________________________________________________________________
AliAnalysisTaskCaloTrackCorrelation::AliAnalysisTaskCaloTrackCorrelation(const char* name) :
  AliAnalysisTaskSE(name),
  fAna(0x0),
  fOutputContainer(0x0),
  fConfigName(""), 
  fCuts(0x0),
  fFirstEvent(0),
  fLastEvent(0),
  fStoreEventSummary(0)
{ 
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // will contain cuts or local params
}

//_________________________________________________________________________
/// Destructor.
//_________________________________________________________________________
AliAnalysisTaskCaloTrackCorrelation::~AliAnalysisTaskCaloTrackCorrelation() 
{  
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  
  if (fOutputContainer)
  {
    fOutputContainer->Clear() ;
    delete fOutputContainer ;
  }
  
  if (fAna) delete fAna;
}

//_________________________________________________________________
/// Create the output container, recover it from the maker 
/// (*AliAnaCaloTrackMaker fAna*) pointer.
//_________________________________________________________________
void AliAnalysisTaskCaloTrackCorrelation::UserCreateOutputObjects()
{
  AliDebug(1,"Begin");
  
  // Get list of aod arrays, add each aod array to analysis frame
  TList * list = fAna->FillAndGetAODBranchList(); //Loop the analysis and create the list of branches
  
  AliDebug(1,Form("n AOD branches %d",list->GetEntries()));
  
  // Put the delta AODs in output file, std or delta
  if((fAna->GetReader())->WriteDeltaAODToFile())
  {
    TString deltaAODName = (fAna->GetReader())->GetDeltaAODFileName();
    for(Int_t iaod = 0; iaod < list->GetEntries(); iaod++)
    {
      TClonesArray * array = (TClonesArray*) list->At(iaod);
      if(deltaAODName!="") AddAODBranch("TClonesArray", &array, deltaAODName);//Put it in DeltaAOD file
      else AddAODBranch("TClonesArray", &array);//Put it in standard AOD file
    }
  }
  
  // Histograms container
  OpenFile(1);
  fOutputContainer = fAna->GetOutputContainer();
  
  AliDebug(1,Form("n histograms %d",fOutputContainer->GetEntries()));
  
  fOutputContainer->SetOwner(kTRUE);
  
  AliDebug(1,"End");
  
  PostData(1,fOutputContainer);
}

//___________________________________________________
/// Local Initialization.
/// Call the Init to initialize the configuration of the analysis.
//___________________________________________________
void AliAnalysisTaskCaloTrackCorrelation::LocalInit()
{  
  Init();
}

//______________________________________________
/// Analysis configuration, if provided, and initialization.
//______________________________________________
void AliAnalysisTaskCaloTrackCorrelation::Init()
{
  AliDebug(1,"Begin");
  
  if( fDebug >= 0 )
    (AliAnalysisManager::GetAnalysisManager())->AddClassDebug(this->ClassName(),fDebug);
  
  // Call configuration file if specified
  
  if (fConfigName.Length())
  {
    AliInfo(Form("### Configuration file is %s.C ###", fConfigName.Data()));
    gROOT->LoadMacro(fConfigName+".C");
    fAna = (AliAnaCaloTrackCorrMaker*) gInterpreter->ProcessLine("ConfigAnalysis()");
  }
  
  if(!fAna)
  {
    AliFatal("Analysis maker pointer not initialized, no analysis specified, STOP!");
    return; // coverity
  }
  
  // Add different generator particles to PDG Data Base
  // to avoid problems when reading MC generator particles
  AliPDG::AddParticlesToPdgDataBase();
  
  // Set in the reader the name of the task in case is needed
  (fAna->GetReader())->SetTaskName(GetName());
	
  // Initialise analysis
  fAna->Init();
  
  // Delta AOD
  if((fAna->GetReader())->GetDeltaAODFileName()!="")
    AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile((fAna->GetReader())->GetDeltaAODFileName());
  
  // Selected Trigger
  if(fAna->GetReader()->IsEventTriggerAtSEOn()) fAna->GetReader()->SetEventTriggerMask(GetCollisionCandidates());
  
  AliDebug(1,"End");
}

//______________________________________________________________________
/// Execute analysis for current event.
//______________________________________________________________________
void AliAnalysisTaskCaloTrackCorrelation::UserExec(Option_t */*option*/)
{  
  if ( !fAna->IsEventProcessed() ) return;
  
  if ( (fLastEvent  > 0 && Entry() > fLastEvent )  || 
       (fFirstEvent > 0 && Entry() < fFirstEvent)     ) return ;
  
  AliDebug(1,Form("Begin event %d", (Int_t) Entry()));
    
  // Get the type of data, check if type is correct
  Int_t  datatype = fAna->GetReader()->GetDataType();
  if(datatype != AliCaloTrackReader::kESD && datatype != AliCaloTrackReader::kAOD &&
     datatype != AliCaloTrackReader::kMC)
  {
    AliError("Wrong type of data");
    return ;
  }
  
  fAna->GetReader()->SetInputOutputMCEvent(InputEvent(), AODEvent(), MCEvent());
  
  // Process event
  fAna->ProcessEvent((Int_t) Entry(), CurrentFileName());
  
  PostData(1, fOutputContainer);
  
  AliDebug(1,"End");
  
  //gObjectTable->Print();
}

//_______________________________________________________________________
/// Terminate analysis. Do some plots (plotting not used so far).
//_______________________________________________________________________
void AliAnalysisTaskCaloTrackCorrelation::Terminate(Option_t */*option*/)
{
  // Get merged histograms from the output container
  // Propagate histagrams to maker
  fAna->Terminate((TList*)GetOutputData(1));
  
  // Create cuts/param objects and publish to slot
  fCuts = fAna->GetListOfAnalysisCuts();
  fCuts ->SetOwner(kTRUE);
  
  // Post Data
  PostData(2, fCuts);
}

//__________________________________________________________
/// Put in the output some standard event summary histograms.
//__________________________________________________________
void AliAnalysisTaskCaloTrackCorrelation::FinishTaskOutput()
{
  if ( !fStoreEventSummary ) return ;
    
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  
  AliInputEventHandler *inputH = dynamic_cast<AliInputEventHandler*>(am->GetInputEventHandler());
  
  if (!inputH) return;
  
  TH2F *histStat = dynamic_cast<TH2F*>(inputH->GetStatistics());
  TH2F *histBin0 = dynamic_cast<TH2F*>(inputH->GetStatistics("BIN0"));
  
  if ( histStat ) 
  {
    if ( histStat == histBin0 ) histBin0 = 0 ;
    
    histStat  = (TH2F*) histStat->Clone(Form("%s_%s",histStat->GetName(),"CaloTrackCorr"));
    
    fOutputContainer->Add(histStat);
  }
   
  if ( histBin0 )
  {
    histBin0  = (TH2F*) histBin0->Clone(Form("%s_%s",histBin0->GetName(),"CaloTrackCorr"));

    fOutputContainer->Add(histBin0);
  }
}

