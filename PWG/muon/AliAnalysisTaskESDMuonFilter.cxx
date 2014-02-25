/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///
/// Add the muon tracks to the generic AOD track branch during the
/// filtering of the ESD.
///
///
/// Note that we :
///
///   - completely disable all the branches that are not required by (most) the muon analyses,
///     e.g. cascades, v0s, kinks, jets, etc...
///   - filter severely the tracks (keep only muon tracks) and vertices (keep only primary -including
///     pile-up - vertices) branches
///
/// \see AliAODMuonReplicator
///
/// (see AddFilteredAOD method)
///

#include "AliAnalysisTaskESDMuonFilter.h"

#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisNonMuonTrackCuts.h"
#include "AliAnalysisNonPrimaryVertices.h"
#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODExtension.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMuonReplicator.h"
#include "AliAODVertex.h"
#include "AliCodeTimer.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDMuonTrack.h"
#include "AliESDtrack.h"
#include "AliESDMuonGlobalTrack.h"   // AU
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMultiplicity.h"
#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>

using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskESDMuonFilter)

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter(Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Bool_t withSPDtracklets):
  AliAnalysisTaskSE(),
  fTrackFilter(0x0),
  fEnableMuonAOD(kFALSE),
  fEnableDimuonAOD(kFALSE),
  fOnlyMuon(onlyMuon),
  fKeepAllEvents(keepAllEvents),
  fMCMode(mcMode),
  fWithSPDTracklets(withSPDtracklets)
{
  /// Default constructor
}

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter(const char* name, Bool_t onlyMuon, Bool_t keepAllEvents, Int_t mcMode, Bool_t withSPDtracklets):
  AliAnalysisTaskSE(name),
  fTrackFilter(0x0),
  fEnableMuonAOD(kFALSE),
  fEnableDimuonAOD(kFALSE),
  fOnlyMuon(onlyMuon),
  fKeepAllEvents(keepAllEvents),
  fMCMode(mcMode),
  fWithSPDTracklets(withSPDtracklets)
{
  /// Constructor
}

//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::UserCreateOutputObjects()
{
  /// Create the output container
  if (fTrackFilter) OutputTree()->GetUserInfo()->Add(fTrackFilter);
}

//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::PrintTask(Option_t *option, Int_t indent) const
{
  /// Specify how we are configured
  
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
    if ( fMCMode == 1 )
    {
      cout << spaces.Data() << "  (will write MC information irrespective of whether or not we have reconstructed muons in the event)" << endl;
    }
    else
    {
      cout << spaces.Data() << "  (will write MC information only if we have reconstructed muons in the event)" << endl;
    }
  }
  
  if ( fWithSPDTracklets )
  {
    cout << spaces.Data() << "Will also keep SPD tracklets" << endl;
  }
}

//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::AddFilteredAOD(const char* aodfilename, const char* title)
{
  /// Add an output filtered and replicated aod
  
  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) Fatal("UserCreateOutputObjects", "No AOD handler");

  AliAODExtension* ext = aodH->AddFilteredAOD(aodfilename,title);

  if (!ext) return;
  
  if ( fOnlyMuon ) 
  {    
    
    AliAODMuonReplicator* murep = new AliAODMuonReplicator("MuonReplicator",
                                                           "remove non muon tracks and non primary or pileup vertices",
                                                           new AliAnalysisNonMuonTrackCuts,
                                                           new AliAnalysisNonPrimaryVertices,
                                                           fMCMode,
                                                           kFALSE,
                                                           fWithSPDTracklets);
    
    ext->DropUnspecifiedBranches(); // all branches not part of a FilterBranch call (below) will be dropped
    
    ext->FilterBranch("tracks",murep);    
    ext->FilterBranch("vertices",murep);  
    ext->FilterBranch("dimuons",murep);
    ext->FilterBranch("AliAODVZERO",murep);
    ext->FilterBranch("AliAODTZERO",murep);
    ext->FilterBranch("AliAODZDC",murep);
    
    if ( fWithSPDTracklets )
    {
      ext->FilterBranch("tracklets",murep);
    }
    
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

//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::Init()
{
  /// Initialization
  if(fEnableMuonAOD) AddFilteredAOD("AliAOD.Muons.root", "MuonEvents");
  if(fEnableDimuonAOD) AddFilteredAOD("AliAOD.Dimuons.root", "DimuonEvents");    
}


//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::UserExec(Option_t */*option*/)
{
  /// Execute analysis for current event
  
  Long64_t ientry = Entry();
  if(fDebug)printf("Muon Filter: Analysing event # %5d\n", (Int_t) ientry);
  
  ConvertESDtoAOD();
}

//______________________________________________________________________________
void AliAnalysisTaskESDMuonFilter::ConvertESDtoAOD() 
{
  /// ESD Muon Filter analysis task executed for each event
  
  AliCodeTimerAuto("",0);
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) 
  {
    AliError("Could not get input ESD event");
    return;    
  }
  
  AliMCEventHandler *mcH = static_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    
  // Define arrays for muons
  Double_t pos[3];
  Double_t p[3];
  //  Double_t pid[10];
  
  // has to be changed once the muon pid is provided by the ESD
  //  for (Int_t i = 0; i < 10; pid[i++] = 0.) {}
  //  pid[AliAODTrack::kMuon]=1.;
  
  AliAODHeader* header = AODEvent()->GetHeader();
  AliAODTrack *aodTrack = 0x0;
  AliESDMuonTrack *esdMuTrack = 0x0;
  
  // Access to the AOD container of tracks
  TClonesArray &tracks = *(AODEvent()->GetTracks());
  Int_t jTracks = tracks.GetEntriesFast();
  
  // Read primary vertex from AOD event 
  AliAODVertex *primary = AODEvent()->GetPrimaryVertex();
  if (fDebug && primary) primary->Print();
  
  // ----------------- MUON TRACKS -----------------------------------------------------

  // Loop on muon tracks to fill the AOD track branch
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();

  for (Int_t iTrack=0; iTrack<nMuTracks; ++iTrack) esd->GetMuonTrack(iTrack)->SetESDEvent(esd);
  
  // Update number of positive and negative tracks from AOD event (M.G.)
  Int_t nTracks    = header->GetRefMultiplicity();
  Int_t nPosTracks = header->GetRefMultiplicityPos();
  Int_t nNegTracks = header->GetRefMultiplicityNeg();
  
  // Access to the AOD container of dimuons
  TClonesArray &dimuons = *(AODEvent()->GetDimuons());
  
  Int_t nMuons=0;
  Int_t nDimuons=0;
  Int_t jDimuons=0;
  Int_t nMuonTrack[100];
  UChar_t itsClusMap(0);
  
  for(int imuon=0;imuon<100;imuon++) nMuonTrack[imuon]=0;
  
  for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack)
  {
    esdMuTrack = esd->GetMuonTrack(nMuTrack);
    
    if (!esdMuTrack->ContainTrackerData()) continue;
    
    UInt_t selectInfo(0);
    
    // Track selection
    if (fTrackFilter) {
     	selectInfo = fTrackFilter->IsSelected(esdMuTrack);
     	if (!selectInfo) {
     	  continue;
     	}  
    }
    
    p[0] = esdMuTrack->Px(); 
    p[1] = esdMuTrack->Py(); 
    p[2] = esdMuTrack->Pz();
    
    pos[0] = esdMuTrack->GetNonBendingCoor(); 
    pos[1] = esdMuTrack->GetBendingCoor(); 
    pos[2] = esdMuTrack->GetZ();
    
    if (mcH) mcH->SelectParticle(esdMuTrack->GetLabel()); // to insure that particle's ancestors will be in output MC branches
    
    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdMuTrack->GetUniqueID(), // ID
                                                  esdMuTrack->GetLabel(), // label
                                                  p, // momentum
                                                  kTRUE, // cartesian coordinate system
                                                  pos, // position
                                                  kFALSE, // isDCA
                                                  0x0, // covariance matrix
                                                  esdMuTrack->Charge(), // charge
                                                  itsClusMap, // ITSClusterMap
                                                  //pid, // pid
                                                  primary, // primary vertex
                                                  kFALSE, // used for vertex fit?
                                                  kFALSE, // used for primary vertex fit?
                                                  AliAODTrack::kPrimary,// track type
                                                  selectInfo); 
    
    aodTrack->SetPIDForTracking(AliPID::kMuon);
    aodTrack->SetXYAtDCA(esdMuTrack->GetNonBendingCoorAtDCA(), esdMuTrack->GetBendingCoorAtDCA());
    aodTrack->SetPxPyPzAtDCA(esdMuTrack->PxAtDCA(), esdMuTrack->PyAtDCA(), esdMuTrack->PzAtDCA());
    aodTrack->SetRAtAbsorberEnd(esdMuTrack->GetRAtAbsorberEnd());
    aodTrack->ConvertAliPIDtoAODPID();
    aodTrack->SetChi2perNDF(esdMuTrack->GetChi2() / (2.*esdMuTrack->GetNHit() - 5.));
    aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
    aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
    UInt_t pattern = esdMuTrack->GetHitsPatternInTrigCh();
    AliESDMuonTrack::AddEffInfo(pattern, 0, esdMuTrack->LoCircuit(), (AliESDMuonTrack::EAliTriggerChPatternFlag)0);
    aodTrack->SetMUONtrigHitsMapTrg(pattern);
    aodTrack->SetMUONtrigHitsMapTrk(esdMuTrack->GetHitsPatternInTrigChTrk());
    aodTrack->SetMuonClusterMap(esdMuTrack->GetMuonClusterMap());
    aodTrack->SetMatchTrigger(esdMuTrack->GetMatchTrigger());
    aodTrack->Connected(esdMuTrack->IsConnected());
    primary->AddDaughter(aodTrack);
    
    ++nTracks;
    if (esdMuTrack->Charge() > 0) nPosTracks++;
    else nNegTracks++;
    
    nMuonTrack[nMuons]= jTracks-1;
    ++nMuons;
  }
  
  if(nMuons>=2) 
  { 
    for(int i=0;i<nMuons;i++){
      Int_t index0 = nMuonTrack[i];
      for(int j=i+1;j<nMuons;j++){
        Int_t index1 = nMuonTrack[j];
        new(dimuons[jDimuons++]) AliAODDimuon(tracks.At(index0),tracks.At(index1));
        ++nDimuons;
        if (fDebug > 1){
          AliAODDimuon *dimuon0 = (AliAODDimuon*)dimuons.At(jDimuons-1);
          printf("Dimuon: mass = %f, px=%f, py=%f, pz=%f\n",dimuon0->M(),dimuon0->Px(),dimuon0->Py(),dimuon0->Pz());  
          AliAODTrack  *mu0 = (AliAODTrack*) dimuon0->GetMu(0);
          AliAODTrack  *mu1 = (AliAODTrack*) dimuon0->GetMu(1);
          printf("Muon0 px=%f py=%f pz=%f\n",mu0->Px(),mu0->Py(),mu0->Pz());
          printf("Muon1 px=%f py=%f pz=%f\n",mu1->Px(),mu1->Py(),mu1->Pz());
        }  
      }
    }
  }
  
  
  header->SetRefMultiplicity(nTracks); 
  header->SetRefMultiplicityPos(nPosTracks);
  header->SetRefMultiplicityNeg(nNegTracks);
  header->SetNumberOfMuons(nMuons);
  header->SetNumberOfDimuons(nDimuons);
  
  // ----------------- MFT + MUON TRACKS -----------------------------------------------------

  AliESDMuonGlobalTrack *esdMuGlobalTrack = 0x0;

  // Loop on muon global tracks to fill the AOD track branch. 
  // We won't update the number of total, pos and neg tracks in the event

  Int_t nMuGlobalTracks = esd->GetNumberOfMuonGlobalTracks();

  for (Int_t iTrack=0; iTrack<nMuGlobalTracks; ++iTrack) esd->GetMuonGlobalTrack(iTrack)->SetESDEvent(esd);
  
  Int_t nGlobalMuons=0;
  Int_t nGlobalDimuons=0;
  Int_t jGlobalDimuons=0;
  Int_t nMuonGlobalTrack[100];
  itsClusMap = 0;
  
  for (Int_t iMuon=0; iMuon<100; iMuon++) nMuonGlobalTrack[iMuon]=0;    // position of the i-th muon track in the tracks array of the AOD event  

  for (Int_t nMuTrack=0; nMuTrack<nMuGlobalTracks; ++nMuTrack) {

    esdMuGlobalTrack = esd->GetMuonGlobalTrack(nMuTrack);

    if (!esdMuGlobalTrack->ContainTrackerData()) continue;
    
    UInt_t selectInfo(0);
    
    // Track selection
    if (fTrackFilter) {
      selectInfo = fTrackFilter->IsSelected(esdMuGlobalTrack);
      if (!selectInfo) {
	continue;
      }  
    }
    
    p[0] = esdMuGlobalTrack->Px(); 
    p[1] = esdMuGlobalTrack->Py(); 
    p[2] = esdMuGlobalTrack->Pz();
    
    esdMuGlobalTrack -> GetFirstTrackingPoint(pos);

    if (mcH) mcH->SelectParticle(esdMuGlobalTrack->GetLabel()); // to insure that particle's ancestors will be in output MC branches
    
    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdMuGlobalTrack->GetUniqueID(), // ID
                                                  esdMuGlobalTrack->GetLabel(),    // label
                                                  p,                               // momentum
                                                  kTRUE,                           // cartesian coordinate system
                                                  pos,                             // position
                                                  kFALSE,                          // isDCA
                                                  0x0,                             // covariance matrix
                                                  esdMuGlobalTrack->Charge(),      // charge
                                                  itsClusMap,                      // ITSClusterMap
                                                  primary,                         // origin vertex
                                                  kFALSE,                          // used for vertex fit?
                                                  kFALSE,                          // used for primary vertex fit?
                                                  AliAODTrack::kPrimary,           // track type
                                                  selectInfo);

    Double_t xyAtVertex[2] = {0};
    esdMuGlobalTrack -> GetXYAtVertex(xyAtVertex);
    
    aodTrack->SetIsMuonGlobalTrack(kTRUE);

    aodTrack->SetMFTClusterPattern(esdMuGlobalTrack->GetMFTClusterPattern());
    aodTrack->SetXYAtDCA(xyAtVertex[0], xyAtVertex[1]);
    aodTrack->SetPxPyPzAtDCA(p[0], p[1], p[2]);
    aodTrack->SetRAtAbsorberEnd(esdMuGlobalTrack->GetRAtAbsorberEnd());
    aodTrack->ConvertAliPIDtoAODPID();
    aodTrack->SetChi2perNDF(esdMuGlobalTrack->GetChi2OverNdf());
    aodTrack->SetChi2MatchTrigger(esdMuGlobalTrack->GetChi2MatchTrigger());
    aodTrack->SetHitsPatternInTrigCh(esdMuGlobalTrack->GetHitsPatternInTrigCh());
    UInt_t pattern = esdMuGlobalTrack->GetHitsPatternInTrigCh();
    AliESDMuonTrack::AddEffInfo(pattern, 0, esdMuGlobalTrack->GetLoCircuit(), (AliESDMuonTrack::EAliTriggerChPatternFlag)0);
    aodTrack->SetMUONtrigHitsMapTrg(pattern);
    aodTrack->SetMUONtrigHitsMapTrk(esdMuGlobalTrack->GetHitsPatternInTrigChTrk());
    aodTrack->SetMuonClusterMap(esdMuGlobalTrack->GetMuonClusterMap());
    aodTrack->SetMatchTrigger(esdMuGlobalTrack->GetMatchTrigger());
    aodTrack->Connected(esdMuGlobalTrack->IsConnected());

    printf("Added Muon Global Track %d of %d\n",nMuTrack,nMuGlobalTracks);
    aodTrack->Print();     // to be removed!!!

    primary->AddDaughter(aodTrack);
    
    nMuonGlobalTrack[nGlobalMuons] = jTracks-1;
    ++nGlobalMuons;
  }

  if (nGlobalMuons >= 2) { 
    for (Int_t i=0; i<nGlobalMuons; i++) {
      Int_t index0 = nMuonGlobalTrack[i];
      for (Int_t j=i+1; j<nGlobalMuons; j++) {
        Int_t index1 = nMuonGlobalTrack[j];
        new (dimuons[jGlobalDimuons++]) AliAODDimuon(tracks.At(index0), tracks.At(index1));
        ++nDimuons;
        if (fDebug > 1) {
          AliAODDimuon *dimuon0 = (AliAODDimuon*)dimuons.At(jGlobalDimuons-1);
          printf("Dimuon: mass = %f, px=%f, py=%f, pz=%f\n",dimuon0->M(),dimuon0->Px(),dimuon0->Py(),dimuon0->Pz());  
          AliAODTrack  *mu0 = (AliAODTrack*) dimuon0->GetMu(0);
          AliAODTrack  *mu1 = (AliAODTrack*) dimuon0->GetMu(1);
          printf("Muon0 px=%f py=%f pz=%f\n",mu0->Px(),mu0->Py(),mu0->Pz());
          printf("Muon1 px=%f py=%f pz=%f\n",mu1->Px(),mu1->Py(),mu1->Pz());
        }  
      }
    }
  }
  
  header->SetNumberOfGlobalMuons(nGlobalMuons);
  header->SetNumberOfGlobalDimuons(nGlobalDimuons);
  
  // -----------------------------------------------------------------------

  AliAODHandler* handler = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  
  if ( handler && fEnableMuonAOD && ( (nMuons>0) || fKeepAllEvents ) )
  {
    AliAODExtension *extMuons = handler->GetFilteredAOD("AliAOD.Muons.root");
    if ( extMuons ) extMuons->SelectEvent();
  }
  
  if ( handler && fEnableDimuonAOD && ( (nMuons>1) || fKeepAllEvents )  )
  {
    AliAODExtension *extDimuons = handler->GetFilteredAOD("AliAOD.Dimuons.root");
    if ( extDimuons ) extDimuons->SelectEvent();
  }
  
}
