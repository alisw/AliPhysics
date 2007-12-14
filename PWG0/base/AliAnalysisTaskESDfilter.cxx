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

/* $Id$ */
 
#include <TChain.h>
#include <TFile.h>

#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDkink.h"
#include "AliESDcascade.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskESDfilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter():
    fDebug(0),
    fTree(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeA(0x0),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0)
{
  // Default constructor
}

AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter(const char* name):
    AliAnalysisTask(name, "AnalysisTaskESDfilter"),
    fDebug(0),
    fTree(0x0),
    fESD(0x0),
    fAOD(0x0),
    fTreeA(0x0),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0)
{
  // Default constructor
    DefineInput (0, TChain::Class());
    DefineOutput(0, TTree::Class());
}

void AliAnalysisTaskESDfilter::CreateOutputObjects()
{
// Create the output container
    if (fDebug > 1) AliInfo("CreateOutPutData() \n");
    AliAODHandler* handler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    
    fAOD   = handler->GetAOD();
    fTreeA = handler->GetTree();
    fTreeA->GetUserInfo()->Add(fTrackFilter);
}

void AliAnalysisTaskESDfilter::Init()
{
    // Initialization
    if (fDebug > 1) AliInfo("Init() \n");
    // Call configuration file
}

void AliAnalysisTaskESDfilter::ConnectInputData(Option_t */*option*/)
{
// Connect the input data
//
    if (fDebug > 1) AliInfo("ConnectInputData() \n");
    // Input 
    AliESDInputHandler* esdH = (AliESDInputHandler*) 
	((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
    fESD = (AliESDEvent*) esdH->GetEvent();
    fTree = esdH->GetTree();
}

void AliAnalysisTaskESDfilter::Exec(Option_t */*option*/)
{
// Execute analysis for current event
//
    // ESD Filter analysis task executed for each event
    AliESD* old = fESD->GetAliESDOld();
    
    Long64_t ientry = fTree->GetReadEntry();
    printf("Filter: Analysing event # %5d\n", (Int_t) ientry);

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
    Int_t nTracks    = fESD->GetNumberOfTracks();
    Int_t nPosTracks = 0;
    for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) 
	if (fESD->GetTrack(iTrack)->GetSign()> 0) nPosTracks++;
    
    // Update the header

    AliAODHeader* header = fAOD->GetHeader();
    header->SetRunNumber(fESD->GetRunNumber());
    if (old) {
	header->SetBunchCrossNumber(0);
	header->SetOrbitNumber(0);
	header->SetPeriodNumber(0);
	header->SetEventType(0);
	header->SetMuonMagFieldScale(-999.);
	header->SetCentrality(-999.);
    } else {
	header->SetBunchCrossNumber(fESD->GetBunchCrossNumber());
	header->SetOrbitNumber(fESD->GetOrbitNumber());
	header->SetPeriodNumber(fESD->GetPeriodNumber());
	header->SetEventType(fESD->GetEventType());
	header->SetMuonMagFieldScale(-999.);
	header->SetCentrality(-999.);
    }
    
    header->SetTriggerMask(fESD->GetTriggerMask()); 
    header->SetTriggerCluster(fESD->GetTriggerCluster());
    header->SetMagneticField(fESD->GetMagneticField());
    header->SetZDCN1Energy(fESD->GetZDCN1Energy());
    header->SetZDCP1Energy(fESD->GetZDCP1Energy());
    header->SetZDCN2Energy(fESD->GetZDCN2Energy());
    header->SetZDCP2Energy(fESD->GetZDCP2Energy());
    header->SetZDCEMEnergy(fESD->GetZDCEMEnergy(0),fESD->GetZDCEMEnergy(1));
    header->SetRefMultiplicity(nTracks);
    header->SetRefMultiplicityPos(nPosTracks);
    header->SetRefMultiplicityNeg(nTracks - nPosTracks);
//
//    
    Int_t nV0s      = fESD->GetNumberOfV0s();
    Int_t nCascades = fESD->GetNumberOfCascades();
    Int_t nKinks    = fESD->GetNumberOfKinks();
    Int_t nVertices = nV0s + nCascades + nKinks + 1 /* = prim. vtx*/;
    Int_t nJets     = 0;
    Int_t nCaloClus = fESD->GetNumberOfCaloClusters();
    Int_t nFmdClus  = 0;
    Int_t nPmdClus  = fESD->GetNumberOfPmdTracks();
    
    printf("   NV0=%d  NCASCADES=%d  NKINKS=%d\n", nV0s, nCascades, nKinks);

//    fAOD->ResetStd(nTracks, nVertices, nV0s+nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);

    AliAODTrack *aodTrack;
    
    
    // Array to take into account the tracks already added to the AOD
    Bool_t * usedTrack = NULL;
    if (nTracks>0) {
	usedTrack = new Bool_t[nTracks];
	for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) usedTrack[iTrack]=kFALSE;
    }
    // Array to take into account the V0s already added to the AOD
    Bool_t * usedV0 = NULL;
    if (nV0s>0) {
	usedV0 = new Bool_t[nV0s];
	for (Int_t iV0=0; iV0<nV0s; ++iV0) usedV0[iV0]=kFALSE;
    }
    // Array to take into account the kinks already added to the AOD
    Bool_t * usedKink = NULL;
    if (nKinks>0) {
	usedKink = new Bool_t[nKinks];
	for (Int_t iKink=0; iKink<nKinks; ++iKink) usedKink[iKink]=kFALSE;
    }
    
    // Access to the AOD container of vertices
    TClonesArray &vertices = *(fAOD->GetVertices());
    Int_t jVertices=0;
    
    // Access to the AOD container of tracks
    TClonesArray &tracks = *(fAOD->GetTracks());
    Int_t jTracks=0; 
    
    // Add primary vertex. The primary tracks will be defined
    // after the loops on the composite objects (V0, cascades, kinks)
    const AliESDVertex *vtx = fESD->GetPrimaryVertex();
    
    vtx->GetXYZ(pos); // position
    vtx->GetCovMatrix(covVtx); //covariance matrix
    
    AliAODVertex * primary = new(vertices[jVertices++])
	AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, AliAODVertex::kPrimary);
    primary->Print();


    // Create vertices starting from the most complex objects
    Double_t chi2 = 0.;
    
    // Cascades
    for (Int_t nCascade = 0; nCascade < nCascades; ++nCascade) {
	AliESDcascade *cascade = fESD->GetCascade(nCascade);
	
	cascade->GetXYZ(pos[0], pos[1], pos[2]);

	if (!old) {
	    chi2 = vtx->GetChi2toNDF();
	    cascade->GetPosCovXi(covVtx);
	} else {
	    chi2 = -999.;
	}
	// Add the cascade vertex
	AliAODVertex * vcascade = new(vertices[jVertices++]) AliAODVertex(pos,
									  covVtx,
									  chi2, // = chi2/NDF since NDF = 2*2-3 (AM)
									  primary,
									  AliAODVertex::kCascade);
	
	primary->AddDaughter(vcascade);
	
	// Add the V0 from the cascade. The ESD class have to be optimized...
	// Now we have to search for the corresponding Vo in the list of V0s
	// using the indeces of the positive and negative tracks
	
	Int_t posFromV0 = cascade->GetPindex();
	Int_t negFromV0 = cascade->GetNindex();
	
      
	AliESDv0 * v0 = 0x0;
	Int_t indV0 = -1;
	
	for (Int_t iV0=0; iV0<nV0s; ++iV0) {
	    
	    v0 = fESD->GetV0(iV0);
	    Int_t posV0 = v0->GetPindex();
	    Int_t negV0 = v0->GetNindex();
	    
	    if (posV0==posFromV0 && negV0==negFromV0) {
		indV0 = iV0;
		break;
	    }
	}
	
	AliAODVertex * vV0FromCascade = 0x0;
	
	if (indV0>-1 && !usedV0[indV0] ) {
	    
	    // the V0 exists in the array of V0s and is not used
	    
	    usedV0[indV0] = kTRUE;
	    
	    v0->GetXYZ(pos[0], pos[1], pos[2]);
	    if (!old) {
		chi2 = v0->GetChi2V0();
		v0->GetPosCov(covVtx);
	    } else {
		chi2 = -999.;
	    }

	    vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
								     covVtx,
								     chi2, // = chi2/NDF since NDF = 2*2-3 (AM)
								     vcascade,
								     AliAODVertex::kV0);
	} else {
	    
	    // the V0 doesn't exist in the array of V0s or was used
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " cascade " << nCascade
//		 << " The V0 " << indV0 
//		 << " doesn't exist in the array of V0s or was used!" << endl;
	    
	    cascade->GetXYZ(pos[0], pos[1], pos[2]);
	    
	    if (!old) {
		chi2 = v0->GetChi2V0();
		cascade->GetPosCov(covVtx);
	    } else {
		chi2 = -999.;
	    }

	    vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
								     covVtx,
								     chi2, // = chi2/NDF since NDF = 2*2-3 (AM)
								     vcascade,
								     AliAODVertex::kV0);
	    vcascade->AddDaughter(vV0FromCascade);
	}
	
	// Add the positive tracks from the V0
	
	if (posFromV0>-1 && !usedTrack[posFromV0]) {
	    
	    usedTrack[posFromV0] = kTRUE;
	    
	    AliESDtrack *esdTrack = fESD->GetTrack(posFromV0);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    vV0FromCascade->AddDaughter(aodTrack =
					new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
									   esdTrack->GetLabel(), 
									   p, 
									   kTRUE,
									   pos,
									   kFALSE,
									   covTr, 
									   (Short_t)esdTrack->GetSign(),
									   esdTrack->GetITSClusterMap(), 
									   pid,
									   vV0FromCascade,
									   kTRUE,  // check if this is right
									   kFALSE, // check if this is right
									   AliAODTrack::kSecondary)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " cascade " << nCascade
//		 << " track " << posFromV0 << " has already been used!" << endl;
	}
	
	// Add the negative tracks from the V0
	
	if (negFromV0>-1 && !usedTrack[negFromV0]) {
	    
	    usedTrack[negFromV0] = kTRUE;
	    
	    AliESDtrack *esdTrack = fESD->GetTrack(negFromV0);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    vV0FromCascade->AddDaughter(aodTrack =
					new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
									   esdTrack->GetLabel(),
									   p,
									   kTRUE,
									   pos,
									   kFALSE,
									   covTr, 
									   (Short_t)esdTrack->GetSign(),
									   esdTrack->GetITSClusterMap(), 
									   pid,
									   vV0FromCascade,
									   kTRUE,  // check if this is right
									   kFALSE, // check if this is right
									   AliAODTrack::kSecondary)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " cascade " << nCascade
//		 << " track " << negFromV0 << " has already been used!" << endl;
	}
	
	// Add the bachelor track from the cascade
	
	Int_t bachelor = cascade->GetBindex();
	
	if(bachelor>-1 && !usedTrack[bachelor]) {
	    
	    usedTrack[bachelor] = kTRUE;
	    
	    AliESDtrack *esdTrack = fESD->GetTrack(bachelor);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    vcascade->AddDaughter(aodTrack =
				  new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
								     esdTrack->GetLabel(),
								     p,
								     kTRUE,
								     pos,
								     kFALSE,
								     covTr, 
								     (Short_t)esdTrack->GetSign(),
								     esdTrack->GetITSClusterMap(), 
								     pid,
								     vcascade,
								     kTRUE,  // check if this is right
								     kFALSE, // check if this is right
								     AliAODTrack::kSecondary)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " cascade " << nCascade
//		 << " track " << bachelor << " has already been used!" << endl;
	}
	
	// Add the primary track of the cascade (if any)
	
    } // end of the loop on cascades
   
    // V0s
    
    for (Int_t nV0 = 0; nV0 < nV0s; ++nV0) {
	
	if (usedV0[nV0]) continue; // skip if aready added to the AOD
	
	AliESDv0 *v0 = fESD->GetV0(nV0);
	
	v0->GetXYZ(pos[0], pos[1], pos[2]);

	if (!old) {
	    chi2 = v0->GetChi2V0();
	    v0->GetPosCov(covVtx);
	} else {
	    chi2 = -999.;
	}


	AliAODVertex * vV0 = 
	    new(vertices[jVertices++]) AliAODVertex(pos,
						    covVtx,
						    chi2, // = chi2/NDF since NDF = 2*2-3
						    primary,
						    AliAODVertex::kV0);
	primary->AddDaughter(vV0);
	
	Int_t posFromV0 = v0->GetPindex();
	Int_t negFromV0 = v0->GetNindex();
	
	// Add the positive tracks from the V0
	
	if (posFromV0>-1 && !usedTrack[posFromV0]) {
	    
	    usedTrack[posFromV0] = kTRUE;
	    
	    AliESDtrack *esdTrack = fESD->GetTrack(posFromV0);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    vV0->AddDaughter(aodTrack =
			     new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
								esdTrack->GetLabel(), 
								p, 
								kTRUE,
								pos,
								kFALSE,
								covTr, 
								(Short_t)esdTrack->GetSign(),
								esdTrack->GetITSClusterMap(), 
								pid,
								vV0,
								kTRUE,  // check if this is right
								kFALSE, // check if this is right
								AliAODTrack::kSecondary)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " V0 " << nV0
//		 << " track " << posFromV0 << " has already been used!" << endl;
	}
	
	// Add the negative tracks from the V0
	
	if (negFromV0>-1 && !usedTrack[negFromV0]) {
	    
	    usedTrack[negFromV0] = kTRUE;
	    
	    AliESDtrack *esdTrack = fESD->GetTrack(negFromV0);
	    esdTrack->GetPxPyPz(p);
	    esdTrack->GetXYZ(pos);
	    esdTrack->GetCovarianceXYZPxPyPz(covTr);
	    esdTrack->GetESDpid(pid);
	    
	    vV0->AddDaughter(aodTrack =
			     new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
								esdTrack->GetLabel(),
								p,
								kTRUE,
								pos,
								kFALSE,
								covTr, 
								(Short_t)esdTrack->GetSign(),
								esdTrack->GetITSClusterMap(), 
								pid,
								vV0,
								kTRUE,  // check if this is right
								kFALSE, // check if this is right
								AliAODTrack::kSecondary)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
//	    cerr << "Error: event " << fESD->GetEventNumberInFile() << " V0 " << nV0
//		 << " track " << negFromV0 << " has already been used!" << endl;
	}
	
    } // end of the loop on V0s
    
    // Kinks: it is a big mess the access to the information in the kinks
    // The loop is on the tracks in order to find the mother and daugther of each kink
    
    
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {
	
	
	AliESDtrack * esdTrack = fESD->GetTrack(iTrack);
	
	Int_t ikink = esdTrack->GetKinkIndex(0);
	
	if (ikink && nKinks) {
	    // Negative kink index: mother, positive: daughter
	    
	    // Search for the second track of the kink
	    
	    for (Int_t jTrack = iTrack+1; jTrack<nTracks; ++jTrack) {
		
		AliESDtrack * esdTrack1 = fESD->GetTrack(jTrack);
		
		Int_t jkink = esdTrack1->GetKinkIndex(0);
		
		if ( TMath::Abs(ikink)==TMath::Abs(jkink) ) {
		    
		    // The two tracks are from the same kink
		    
		    if (usedKink[TMath::Abs(ikink)-1]) continue; // skip used kinks
		    
		    Int_t imother = -1;
		    Int_t idaughter = -1;
		    
		    if (ikink<0 && jkink>0) {
			
			imother = iTrack;
			idaughter = jTrack;
		    }
		    else if (ikink>0 && jkink<0) {
			
			imother = jTrack;
			idaughter = iTrack;
		    }
		    else {
//			cerr << "Error: Wrong combination of kink indexes: "
//			     << ikink << " " << jkink << endl;
			continue;
		    }
		    
		    // Add the mother track
		    
		    AliAODTrack * mother = NULL;
		    
		    if (!usedTrack[imother]) {
			
			usedTrack[imother] = kTRUE;
			
			AliESDtrack *esdTrack = fESD->GetTrack(imother);
			esdTrack->GetPxPyPz(p);
			esdTrack->GetXYZ(pos);
			esdTrack->GetCovarianceXYZPxPyPz(covTr);
			esdTrack->GetESDpid(pid);
			
			mother = 
			    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
							       esdTrack->GetLabel(),
							       p,
							       kTRUE,
							       pos,
							       kFALSE,
							       covTr, 
							       (Short_t)esdTrack->GetSign(),
							       esdTrack->GetITSClusterMap(), 
							       pid,
							       primary,
							       kTRUE, // check if this is right
							       kTRUE, // check if this is right
							       AliAODTrack::kPrimary);
			primary->AddDaughter(mother);
			mother->ConvertAliPIDtoAODPID();
		    }
		    else {
//			cerr << "Error: event " << fESD->GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
//			     << " track " << imother << " has already been used!" << endl;
		  }
		    
		    // Add the kink vertex
		    AliESDkink * kink = fESD->GetKink(TMath::Abs(ikink)-1);
		    
		    AliAODVertex * vkink = 
			new(vertices[jVertices++]) AliAODVertex(kink->GetPosition(),
								NULL,
								0.,
								mother,
								AliAODVertex::kKink);
		    // Add the daughter track
		  
		    AliAODTrack * daughter = NULL;
		    
		    if (!usedTrack[idaughter]) {
			
			usedTrack[idaughter] = kTRUE;
			
			AliESDtrack *esdTrack = fESD->GetTrack(idaughter);
			esdTrack->GetPxPyPz(p);
			esdTrack->GetXYZ(pos);
			esdTrack->GetCovarianceXYZPxPyPz(covTr);
			esdTrack->GetESDpid(pid);
			
			daughter = 
			    new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
							       esdTrack->GetLabel(),
							       p,
							       kTRUE,
							       pos,
							       kFALSE,
							       covTr, 
							       (Short_t)esdTrack->GetSign(),
							       esdTrack->GetITSClusterMap(), 
							       pid,
							       vkink,
							       kTRUE, // check if this is right
							       kTRUE, // check if this is right
							       AliAODTrack::kPrimary);
			vkink->AddDaughter(daughter);
			daughter->ConvertAliPIDtoAODPID();
		    }
		    else {
//			cerr << "Error: event " << fESD->GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
//			     << " track " << idaughter << " has already been used!" << endl;
		    }
		}
	    }
	}      
    }
    
  
    // Tracks (primary and orphan)

    printf("NUMBER OF TRACKS %5d\n", nTracks);
    
    for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
	
	
	if (usedTrack[nTrack]) continue;
	AliESDtrack *esdTrack = fESD->GetTrack(nTrack);
	UInt_t selectInfo = 0;
	//
	// Track selection
	if (fTrackFilter) {
	    selectInfo = fTrackFilter->IsSelected(esdTrack);
	    if (!selectInfo) continue;
	}
	
	//
	esdTrack->GetPxPyPz(p);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);
	
	Float_t impactXY, impactZ;
	
	esdTrack->GetImpactParameters(impactXY,impactZ);
	
	if (impactXY<3) {
	    // track inside the beam pipe
	    
	    primary->AddDaughter(aodTrack =
				 new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
								    esdTrack->GetLabel(),
								    p,
								    kTRUE,
								    pos,
								    kFALSE,
								    covTr, 
								    (Short_t)esdTrack->GetSign(),
								    esdTrack->GetITSClusterMap(), 
								    pid,
								    primary,
								    kTRUE, // check if this is right
								    kTRUE, // check if this is right
								    AliAODTrack::kPrimary, 
								    selectInfo)
		);
	    aodTrack->ConvertAliPIDtoAODPID();
	}
	else {
	    // outside the beam pipe: orphan track
	    aodTrack =
		new(tracks[jTracks++]) AliAODTrack(esdTrack->GetID(),
						   esdTrack->GetLabel(),
						   p,
						   kTRUE,
						   pos,
						   kFALSE,
						   covTr, 
						   (Short_t)esdTrack->GetSign(),
						   esdTrack->GetITSClusterMap(), 
						   pid,
						   NULL,
						   kFALSE, // check if this is right
						   kFALSE, // check if this is right
						   AliAODTrack::kOrphan,
						   selectInfo);
	    aodTrack->ConvertAliPIDtoAODPID();
	}	
    } // end of loop on tracks
    
    // muon tracks
    Int_t nMuTracks = fESD->GetNumberOfMuonTracks();
    for (Int_t nMuTrack = 0; nMuTrack < nMuTracks; ++nMuTrack) {
	
	AliESDMuonTrack *esdMuTrack = fESD->GetMuonTrack(nMuTrack);     
	p[0] = esdMuTrack->Px(); 
	p[1] = esdMuTrack->Py(); 
	p[2] = esdMuTrack->Pz();
	pos[0] = primary->GetX(); 
	pos[1] = primary->GetY(); 
	pos[2] = primary->GetZ();
	
	// has to be changed once the muon pid is provided by the ESD
	for (Int_t i = 0; i < 10; pid[i++] = 0.); pid[AliAODTrack::kMuon]=1.;
	
	primary->AddDaughter(
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
    }
  
    delete [] usedTrack;
    delete [] usedV0;
    delete [] usedKink;


    //
    PostData(0, fTreeA);
    return;
}

void AliAnalysisTaskESDfilter::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}

