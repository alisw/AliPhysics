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

/* $Id: AliAnalysisTaskESDMuonFilter.cxx 24535 2008-03-16 22:43:30Z fca $ */
 
#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisTaskESDMuonFilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskESDMuonFilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter():
    AliAnalysisTaskSE(),
    fTrackFilter(0x0)
{
  // Default constructor
}

AliAnalysisTaskESDMuonFilter::AliAnalysisTaskESDMuonFilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilter(0x0)
{
  // Constructor
}

void AliAnalysisTaskESDMuonFilter::UserCreateOutputObjects()
{
// Create the output container
    OutputTree()->GetUserInfo()->Add(fTrackFilter);
}

void AliAnalysisTaskESDMuonFilter::Init()
{
    // Initialization
    if (fDebug > 1) AliInfo("Init() \n");
    // Call configuration file
}


void AliAnalysisTaskESDMuonFilter::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//					    
  Long64_t ientry = Entry();
  printf("Filter: Analysing event # %5d\n", (Int_t) ientry);

  ConvertESDtoAOD();
}

void AliAnalysisTaskESDMuonFilter::ConvertESDtoAOD() {
    // ESD Filter analysis task executed for each event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    AliESD* old = esd->GetAliESDOld();
    
    // set arrays and pointers
    Double_t pos[3];
    Double_t p[3];
    Double_t covVtx[6];
    Double_t covTr[21];
    Double_t pid[10];

    for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
    for (Int_t i = 0; i < 21; i++) covTr [i] = 0.;

    
  // loop over events and fill them
  
  // Multiplicity information needed by the header (to be revised!)
    Int_t nTracks    = esd->GetNumberOfTracks();
    Int_t nPosTracks = 0;
    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) 
	if (esd->GetTrack(iTrack)->GetSign()> 0) nPosTracks++;
    
    // Update the header

    AliAODHeader* header = AODEvent()->GetHeader();
    header->SetRunNumber(esd->GetRunNumber());
    if (old) {
	header->SetBunchCrossNumber(0);
	header->SetOrbitNumber(0);
	header->SetPeriodNumber(0);
	header->SetEventType(0);
	header->SetMuonMagFieldScale(-999.); // FIXME
	header->SetCentrality(-999.);        // FIXME
    } else {
	header->SetBunchCrossNumber(esd->GetBunchCrossNumber());
	header->SetOrbitNumber(esd->GetOrbitNumber());
	header->SetPeriodNumber(esd->GetPeriodNumber());
	header->SetEventType(esd->GetEventType());
	header->SetMuonMagFieldScale(-999.); // FIXME
	header->SetCentrality(-999.);        // FIXME
    }
    
    header->SetTriggerMask(esd->GetTriggerMask()); 
    header->SetTriggerCluster(esd->GetTriggerCluster());
    header->SetMagneticField(esd->GetMagneticField());
    header->SetZDCN1Energy(esd->GetZDCN1Energy());
    header->SetZDCP1Energy(esd->GetZDCP1Energy());
    header->SetZDCN2Energy(esd->GetZDCN2Energy());
    header->SetZDCP2Energy(esd->GetZDCP2Energy());
    header->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));
    header->SetRefMultiplicity(nTracks);
    header->SetRefMultiplicityPos(nPosTracks);
    header->SetRefMultiplicityNeg(nTracks - nPosTracks);
//
//    
    Int_t nV0s      = esd->GetNumberOfV0s();
    Int_t nCascades = esd->GetNumberOfCascades();
    Int_t nKinks    = esd->GetNumberOfKinks();
    Int_t nVertices = nV0s + 2*nCascades /*could lead to two vertices, one V0 and the Xi */+ nKinks + 1 /* = prim. vtx*/;    
    Int_t nJets     = 0;
    Int_t nCaloClus = esd->GetNumberOfCaloClusters();
    Int_t nFmdClus  = 0;
    Int_t nPmdClus  = esd->GetNumberOfPmdTracks();
    
    printf("   NV0=%d  NCASCADES=%d  NKINKS=%d\n", nV0s, nCascades, nKinks);

    AODEvent()->ResetStd(nTracks, nVertices, nV0s+nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);

    AliAODTrack *aodTrack = 0x0;
    
    // Access to the AOD container of vertices
    TClonesArray &vertices = *(AODEvent()->GetVertices());
    Int_t jVertices=0;
    
    // Access to the AOD container of tracks
    TClonesArray &tracks = *(AODEvent()->GetTracks());
    Int_t jTracks=0; 
    
    // Add primary vertex. 
    const AliESDVertex *vtx = esd->GetPrimaryVertex();
    
    vtx->GetXYZ(pos); // position
    vtx->GetCovMatrix(covVtx); //covariance matrix
    
    AliAODVertex * primary = new(vertices[jVertices++])
	AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, AliAODVertex::kPrimary);
    primary->Print();
    
    // muon tracks
    Int_t nMuTracks = esd->GetNumberOfMuonTracks();
    
    printf("Number of Muon Tracks=%d\n",nMuTracks);
    
    for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
	
	AliESDMuonTrack *esdMuTrack = esd->GetMuonTrack(nMuTrack); 
	    
	p[0] = esdMuTrack->Px(); 
	p[1] = esdMuTrack->Py(); 
	p[2] = esdMuTrack->Pz();
	pos[0] = primary->GetX(); 
	pos[1] = primary->GetY(); 
	pos[2] = primary->GetZ();
	
	// has to be changed once the muon pid is provided by the ESD
	for (Int_t i = 0; i < 10; pid[i++] = 0.); pid[AliAODTrack::kMuon]=1.;
	
//	primary->AddDaughter(
	primary->AddDaughter(aodTrack =
	    new(tracks[jTracks++]) AliAODTrack(0, // no ID provided
					       0, // no label provided
					       p,
					       kTRUE,
					       pos,
					       kFALSE,
					       NULL, // no covariance matrix provided
					       (Short_t)-99, // no charge provided
					       0, // no ITSClusterMap
					       pid,
					       primary,
					       kTRUE,  // check if this is right
					       kTRUE,  // not used for vertex fit
					       AliAODTrack::kPrimary)
	    );
  
	aodTrack->ConvertAliPIDtoAODPID();  
	aodTrack->SetHitsPatternInTrigCh(esdMuTrack->GetHitsPatternInTrigCh());
	Int_t track2Trigger = esdMuTrack->GetMatchTrigger();
	aodTrack->SetMatchTrigger(track2Trigger);
	if (track2Trigger) 
	  aodTrack->SetChi2MatchTrigger(esdMuTrack->GetChi2MatchTrigger());
	else 
	  aodTrack->SetChi2MatchTrigger(0.);
    }
        
    tracks.Expand(jTracks); // remove 'empty slots' due to unwritten tracks
  

    return;
}

void AliAnalysisTaskESDMuonFilter::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}
