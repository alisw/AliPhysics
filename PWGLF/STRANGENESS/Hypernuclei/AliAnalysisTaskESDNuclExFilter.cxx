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

#include "AliAnalysisTaskESDNuclExFilter.h"

#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODMCParticle.h"
#include "AliAODNuclExReplicator.h"
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
#include "AliPIDResponse.h"

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskESDNuclExFilter)
ClassImp(AliAnalysisNonMuonTrackCuts)

////////////////////////////////////////////////////////////////////////

AliAnalysisNonMuonTrackCuts::AliAnalysisNonMuonTrackCuts()
{
  // default ctor 
}

Bool_t AliAnalysisNonMuonTrackCuts::IsSelected(TObject* obj)
{
  // Returns true if the object is a muon track
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);

  ULong_t  status;

  if(track){
    
    status  = (ULong_t)track->GetStatus();

    if(track->GetTPCNcls() > 80 &&
       track->Chi2perNDF() < 5  &&
       track->IsOn(AliAODTrack::kTPCrefit) &&
       track->IsOn(AliAODTrack::kTPCin)    &&
       !track->IsOn(AliAODTrack::kITSpureSA))
      {
	return kTRUE;
      }
  } 
  
  else 
    return kFALSE;
  

}

AliAnalysisNonPrimaryVertices::AliAnalysisNonPrimaryVertices()
{
  // default ctor   
}

Bool_t AliAnalysisNonPrimaryVertices::IsSelected(TObject* obj)
{
  // Returns true if the object is a primary vertex
  
  AliAODVertex* vertex = dynamic_cast<AliAODVertex*>(obj);
  if (vertex)
    {
      if ( vertex->GetType() == AliAODVertex::kPrimary     ||
	   vertex->GetType() == AliAODVertex::kMainSPD     ||
	   vertex->GetType() == AliAODVertex::kPileupSPD   ||
	   vertex->GetType() == AliAODVertex::kPileupTracks||
	   vertex->GetType() == AliAODVertex::kMainTPC )
	{
	  return kTRUE;
	}
    }
  
  //  enum AODVtx_t {kUndef=-1, kPrimary, kKink, kV0, kCascade, kMulti, kMainSPD, kPileupSPD, kPileupTracks,kMainTPC};

  return kFALSE;
  
}

AliAnalysisTaskESDNuclExFilter::AliAnalysisTaskESDNuclExFilter(Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Int_t nsigmaTrk1,Int_t nsigmaTrk2, Int_t partType1,Int_t partType2):
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
  fPIDResponse(0)
{
  // Default constructor
}

AliAnalysisTaskESDNuclExFilter::AliAnalysisTaskESDNuclExFilter(const char* name, Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Int_t nsigmaTrk1,Int_t nsigmaTrk2, Int_t partType1,Int_t partType2):
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
  murep(0),
  fPIDResponse(0)
{
  // Constructor
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilter::UserCreateOutputObjects()
{
  //-----------------------------------------------
  // Particle Identification Setup (new PID object)
  //-----------------------------------------------
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //cout<<"===========================================Manager: "<<man<<endl;
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  //cout<<"ih: "<<inputHandler<<endl;
  
  fPIDResponse = inputHandler->GetPIDResponse();

  // Create the output container
  if (fTrackFilter) OutputTree()->GetUserInfo()->Add(fTrackFilter);
  //cout<<"Sotto"<<endl;
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilter::PrintTask(Option_t *option, Int_t indent) const
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
void AliAnalysisTaskESDNuclExFilter::AddFilteredAOD(const char* aodfilename, const char* title)
{
  
  //cout<<"Entro ne ADDFILTETEDAOD"<<endl;

  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) Fatal("UserCreateOutputObjects", "No AOD handler");
  //cout<<"Add Filterd AOD "<<aodH->AddFilteredAOD(aodfilename,title)<<endl;
  AliAODExtension* ext = aodH->AddFilteredAOD(aodfilename,title);
  //cout<<"Handle inside add filterAOD: "<<aodH<<endl;
  //cout<<"########### ext: "<<ext<<endl;
  
  if (!ext) return;
  
  //cout<<"ONLY MUON?? "<<fOnlyMuon<<endl;

  if ( fOnlyMuon ) 
    {    
    
      //cout<<"Inside fonly muon: "<<endl;

      
      if(!murep)delete murep;

      murep = new AliAODNuclExReplicator("NuclExReplicator",
					 "remove non interesting tracks",
					 new AliAnalysisNonMuonTrackCuts,
					 new AliAnalysisNonPrimaryVertices,
					 fMCMode,fnSigmaTrk1,fnSigmaTrk2,fpartType1,fpartType2);
      
      //cout<<"murep: "<<murep<<endl;
      
      ext->DropUnspecifiedBranches(); // all branches not part of a FilterBranch call (below) will be dropped
      
      // ext->FilterBranch("header",murep);    
      // ext->FilterBranch("tracks",murep);    
      // ext->FilterBranch("vertices",murep);  
      // ext->FilterBranch("dimuons",murep); //per test
      // ext->FilterBranch("AliAODVZERO",murep);
      // ext->FilterBranch("AliAODTZERO",murep);
      
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
  
  //cout<<"fine add filterd"<<endl;
}

//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilter::Init()
{

  //cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sono in INIT"<<endl;
  // Initialization
  if(fEnableMuonAOD) 
    AddFilteredAOD("AliAOD.NuclEx.root", "MuonEvents");
  //cout<<"Fine INIT"<<endl;
  //  if(fEnableDimuonAOD) AddFilteredAOD("AliAOD.Dimuons.root", "DimuonEvents");    
}


//______________________________________________________________________________
void AliAnalysisTaskESDNuclExFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event					    
  
  //cout<<">>>>>>>>>>>>Inside User exec<<<<<<<<<<<<<<<<<"<<endl;

  Long64_t ientry = Entry();
  if(fDebug)printf("Muon Filter: Analysing event # %5d\n", (Int_t) ientry);
  //cout<<"--------------->Enter the user exec<------------------------------"<<endl;
  
  //Check the PIDresponse
  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }

  //***************************************************

  ConvertESDtoAOD();

  //*************************************************
  
  //cout<<"-------------------------------------------------------------------FINE ESD TO AOD CONVERTER!!!"<<endl;
  
  // if(!murep)
  //   delete murep;
  
}





//---- Funziona


void AliAnalysisTaskESDNuclExFilter::ConvertESDtoAOD() 

{
  //cout<<"========================> CONVERT ESD TO AOD <============================="<<endl;

  
  AliVEvent *event = InputEvent();
  //cout<<"VEvent: "<<event<<endl;
  AliAODEvent *lAODevent=(AliAODEvent*)InputEvent();
  //cout<<"AOD Event: "<<event<<endl;
  AliAODHeader* header =lAODevent->GetHeader();
  //cout<<"header :"<<header<<endl;
  Int_t jTracks =  lAODevent->GetNumberOfTracks();
  //cout<<"n jtracks :"<<jTracks<<endl;

  // Read primary vertex from AOD event 
  // AliAODVertex *primary = *(AODEvent()->GetPrimaryVertex());

  AliAODVertex *primary = lAODevent->GetPrimaryVertex();
  if (fDebug && primary) primary->Print();
  //cout<<"Primary vtx x: "<<primary->GetX()<<" "<<primary->GetY()<<" "<<primary->GetZ()<<endl;
  
    
  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  //AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  //cout<<"Mag field Filer: "<<lAODevent->GetMagneticField()<<endl;

  // lAODevent->Print();

  //cout<<"handler inside convert "<<handler<<endl;
  if ( handler ){
    
    //cout<<"Hadler in the loop "<<handler<<endl;

    //cout<<"Inside if handler loop"<<endl;

    AliAODExtension *extNuclEx = handler->GetFilteredAOD("AliAOD.NuclEx.root");
    //  AliAODExtension *extNuclEx = handler->GetFilteredAOD("../1/pass2/AliAOD.root");
    
    //cout<<"extmuon? "<<extNuclEx<<endl;

    if ( extNuclEx ) {				
     //   extNuclEx->Init("");
     extNuclEx->SetEvent(lAODevent);
     extNuclEx->SelectEvent();
     extNuclEx->Print();
     extNuclEx->FinishEvent();
     
     //cout<<"extMuons? "<<extMuons<<endl;

   }
  }


}
//------------------------------------------------
void AliAnalysisTaskESDNuclExFilter::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  //  delete murep;
  
  if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}
