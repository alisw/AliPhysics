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

//
// Create a new AOD starting from the general AOD. This Task can be used also strating 
//from ESD changing the input handler. (Method to be testeted on the grid)
// filtering of the ESD. 
//
// Authors: S. Bufalino (stefania.bufalino@cern.ch)
//          R. Lea      (ramona.lea@cern.ch)
// Based on AliAnalysisTaskESDMuonFilter.cxx  
//
// (see AddFilteredAOD method)
//

#include "AliAnalysisTaskESDNuclExFilterMC.h"

#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAODMCNuclExReplicator.h"
#include "AliAODVertex.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliCodeTimer.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "TF1.h"
//#include "AliPIDResponse.h"

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskESDNuclExFilterMC)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDNuclExFilterMC::AliAnalysisTaskESDNuclExFilterMC(Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Int_t nsigmaTrk1,Int_t nsigmaTrk2, Int_t partType1,Int_t partType2):
  AliAnalysisTaskSE(),
  fTrackFilter(0x0),
  fEnableMuonAOD(kTRUE),
  fEnableDimuonAOD(kTRUE),
  fOnlyMuon(onlyMuon),
  fKeepAllEvents(keepAllEvents),
  fMCMode(mcMode),
  fnSigmaTrk1(nsigmaTrk1),
  fnSigmaTrk2(nsigmaTrk2),
  fpartType1(partType1),
  fpartType2(partType2),
  murep(0x0)
{
  // Default constructor
}

AliAnalysisTaskESDNuclExFilterMC::AliAnalysisTaskESDNuclExFilterMC(const char* name, Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Int_t nsigmaTrk1,Int_t nsigmaTrk2, Int_t partType1,Int_t partType2):
  AliAnalysisTaskSE(name),
  fTrackFilter(0x0),
  fEnableMuonAOD(kTRUE),
  fEnableDimuonAOD(kTRUE),
  fOnlyMuon(onlyMuon),
  fKeepAllEvents(keepAllEvents),
  fMCMode(mcMode),
  fnSigmaTrk1(nsigmaTrk1),
  fnSigmaTrk2(nsigmaTrk2),
  fpartType1(partType1),
  fpartType2(partType2),
  murep(0x0)
  
{
  // Constructor
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilterMC::UserCreateOutputObjects()
{
  
  // Create the output container
  if (fTrackFilter) OutputTree()->GetUserInfo()->Add(fTrackFilter);

}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilterMC::PrintTask(Option_t *option, Int_t indent) const
{
  // Specify how we are configured
  
  AliAnalysisTaskSE::PrintTask(option,indent);
  
  TString spaces(' ',indent+3);
  
  if ( fOnlyMuon ) 
    {
      cout << spaces.Data() << "Keep only muon information " << endl;        
    }
  else 
    {
      cout << spaces.Data() << "Keep all information from standard AOD" << endl;
    }
  
  if ( fKeepAllEvents ) 
    {
      cout << spaces.Data() << "Keep all events, regardless of number of muons" << endl;    
    }
  else 
    {
      cout << spaces.Data() << "Keep only events with at least one muon" << endl;
    }
  
  if ( fMCMode > 0 ) 
    {
      cout << spaces.Data() << "Assuming work on MC data (i.e. will transmit MC branches)" << endl;
    }
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilterMC::AddFilteredAOD(const char* aodfilename, const char* title, Bool_t toMerge)
{
  
  //cout<<"Entro ne ADDFILTETEDAOD"<<endl;

  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) Fatal("UserCreateOutputObjects", "No AOD handler");

  if(aodH){

    AliAODExtension* ext = aodH->AddFilteredAOD(aodfilename,title,toMerge);
    
    if (!ext) return;
    
    if ( fOnlyMuon ) 
      {    
	
	if(!murep)delete murep;
	
	murep = new AliAODMCNuclExReplicator("NuclExReplicator",
					     "remove non interesting tracks",
					     fMCMode,fnSigmaTrk1,fnSigmaTrk2,fpartType1,fpartType2);
	
	
	ext->DropUnspecifiedBranches(); // all branches not part of a FilterBranch call (below) will be dropped
	
	ext->FilterBranch("header",murep);    
	ext->FilterBranch("vertices",murep);    
	ext->FilterBranch("nuclei",murep);  
	ext->FilterBranch("secvertices",murep); //per test
	ext->FilterBranch("daughtertracks",murep);
	
	//cout<<"add filterd aod"<<endl;
	
	if ( fMCMode > 0 ) 
	  {
	    // MC branches will be copied (if present), as they are, but only
	    // for events with at least one muon. 
	    // For events w/o muon, mcparticles array will be empty and mcheader will be dummy
	    // (e.g. strlen(GetGeneratorName())==0)
	    
	    ext->FilterBranch("mcparticles",murep);
	    ext->FilterBranch("mcHeader",murep);
	  }
      }  
    
  }
  //cout<<"fine add filterd"<<endl;
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilterMC::Init()
{

  // Initialization
  if(fEnableMuonAOD) 
    AddFilteredAOD("AliAOD.NuclEx.root", "NuclexFilteredEvents",kTRUE);
}


//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilterMC::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event					    
  
  Long64_t ientry = Entry();
  if(fDebug)printf("Muon Filter: Analysing event # %5d\n", (Int_t) ientry);
  
  //***************************************************

  ConvertESDtoAOD();

  //*************************************************
  
}


void AliAnalysisTaskESDNuclExFilterMC::ConvertESDtoAOD() 

{

  //cout<<"========================> CONVERT ESD TO AOD <============================="<<endl;

  
  AliAODEvent *lAODevent=(AliAODEvent*)InputEvent();

  // Read primary vertex from AOD event 
  // AliAODVertex *primary = *(AODEvent()->GetPrimaryVertex());

  AliAODVertex *primary = lAODevent->GetPrimaryVertex();
  if (fDebug && primary) primary->Print();
  //cout<<"Primary vtx x: "<<primary->GetX()<<" "<<primary->GetY()<<" "<<primary->GetZ()<<endl;
      
  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
 
  if ( handler ){
    
    AliAODExtension *extNuclEx = handler->GetFilteredAOD("AliAOD.NuclEx.root");
    
    if ( extNuclEx ) {				
    
      extNuclEx->SetEvent(lAODevent);
      extNuclEx->SelectEvent();
      extNuclEx->FinishEvent();
      
    }
  }
  
}
//------------------------------------------------
void AliAnalysisTaskESDNuclExFilterMC::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  //  delete murep;
  
  if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}
