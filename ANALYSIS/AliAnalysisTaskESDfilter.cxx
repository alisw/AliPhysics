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

/* $Id: AliAnalysisTaskESDfilter.cxx 24535 2008-03-16 22:43:30Z fca $ */
 
#include <TChain.h>
#include <TTree.h>
#include <TList.h>
#include <TArrayI.h>
#include <TRandom.h>
#include <TParticle.h>

#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliStack.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisFilter.h"
#include "AliESDMuonTrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDkink.h"
#include "AliESDcascade.h"
#include "AliESDPmdTrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliMultiplicity.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskESDfilter)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter():
    AliAnalysisTaskSE(),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0),
    fCascadeFilter(0x0),
    fHighPthreshold(0),
    fPtshape(0x0)
{
  // Default constructor
}

AliAnalysisTaskESDfilter::AliAnalysisTaskESDfilter(const char* name):
    AliAnalysisTaskSE(name),
    fTrackFilter(0x0),
    fKinkFilter(0x0),
    fV0Filter(0x0),
    fCascadeFilter(0x0),
    fHighPthreshold(0),
    fPtshape(0x0)
{
  // Constructor
}

void AliAnalysisTaskESDfilter::UserCreateOutputObjects()
{
  //
  // Create Output Objects conenct filter to outputtree
  // 
    OutputTree()->GetUserInfo()->Add(fTrackFilter);
}

void AliAnalysisTaskESDfilter::Init()
{
    // Initialization
    if (fDebug > 1) AliInfo("Init() \n");
    // Call configuration file
}


void AliAnalysisTaskESDfilter::UserExec(Option_t */*option*/)
{
// Execute analysis for current event
//
					    
  Long64_t ientry = Entry();
  
  if (fDebug > 0) {
      printf("Filter: Analysing event # %5d\n", (Int_t) ientry);
      if (fHighPthreshold == 0) AliInfo("detector PID signals are stored in each track");
      if (!fPtshape) AliInfo("detector PID signals are not stored below the pt threshold");
  }
  
  ConvertESDtoAOD();
}

void AliAnalysisTaskESDfilter::ConvertESDtoAOD() {
    // ESD Filter analysis task executed for each event

    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
    AliESD*      old = esd->GetAliESDOld();

    // Fetch Stack for debuggging if available 
    AliStack *pStack = 0;
    AliMCEventHandler *mcH = 0;
    if(MCEvent()){
      pStack = MCEvent()->Stack();
      mcH = (AliMCEventHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler()); 
    }
    // set arrays and pointers
    Float_t posF[3];
    Double_t pos[3];
    Double_t p[3];
    Double_t momPos[3]; 
    Double_t momNeg[3];
    Double_t momBach[3];
    Double_t momPosAtV0vtx[3];
    Double_t momNegAtV0vtx[3];
    Double_t momBachAtCascadeVtx[3];
    Double_t covVtx[6];
    Double_t covTr[21];
    Double_t pid[10];

    for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
    for (Int_t i = 0; i < 21; i++) covTr [i] = 0.;

    
    // loop over events and fill them
    // Multiplicity information needed by the header (to be revised!)
    Int_t nTracks    = esd->GetNumberOfTracks();
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) esd->GetTrack(iTrack)->SetESDEvent(esd);

    //    if (fDebug > 0) printf("-------------------Bo: Number of ESD tracks %d \n",nTracks);

    Int_t nPosTracks = 0;
    
    // Update the header

    AliAODHeader* header = AODEvent()->GetHeader();
    
    header->SetRunNumber(esd->GetRunNumber());
    if (old) {
	header->SetBunchCrossNumber(0);
	header->SetOrbitNumber(0);
	header->SetPeriodNumber(0);
	header->SetEventType(0);
	header->SetMuonMagFieldScale(-999.);
	header->SetCentrality(-999.);       
    } else {
	header->SetBunchCrossNumber(esd->GetBunchCrossNumber());
	header->SetOrbitNumber(esd->GetOrbitNumber());
	header->SetPeriodNumber(esd->GetPeriodNumber());
	header->SetEventType(esd->GetEventType());
	header->SetCentrality(-999.);        // FIXME
    }
    // Trigger
    header->SetFiredTriggerClasses(esd->GetFiredTriggerClasses());
    header->SetTriggerMask(esd->GetTriggerMask()); 
    header->SetTriggerCluster(esd->GetTriggerCluster());
    

    header->SetMagneticField(esd->GetMagneticField());
    header->SetMuonMagFieldScale(esd->GetCurrentDip()/6000.);
    header->SetZDCN1Energy(esd->GetZDCN1Energy());
    header->SetZDCP1Energy(esd->GetZDCP1Energy());
    header->SetZDCN2Energy(esd->GetZDCN2Energy());
    header->SetZDCP2Energy(esd->GetZDCP2Energy());
    header->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));

    
    Float_t diamxy[2]={esd->GetDiamondX(),esd->GetDiamondY()};
    Float_t diamcov[3]; esd->GetDiamondCovXY(diamcov);
    header->SetDiamond(diamxy,diamcov);
//
//
    Int_t nV0s      = esd->GetNumberOfV0s();
    Int_t nCascades = esd->GetNumberOfCascades();
    Int_t nKinks    = esd->GetNumberOfKinks();
    Int_t nVertices = nV0s + nCascades /*V0 wihtin cascade already counted*/+ nKinks + 1 /* = prim. vtx*/;
    Int_t nJets     = 0;
    Int_t nCaloClus = esd->GetNumberOfCaloClusters();
    Int_t nFmdClus  = 0;
    Int_t nPmdClus  = esd->GetNumberOfPmdTracks();
    
    if (fDebug > 0) 
	printf("   NV0=%d  NCASCADES=%d  NKINKS=%d\n", nV0s, nCascades, nKinks);
       
    AODEvent()->ResetStd(nTracks, nVertices, nV0s, nCascades, nJets, nCaloClus, nFmdClus, nPmdClus);


    AliAODTrack   *aodTrack       = 0x0;
    AliAODPid     *detpid         = 0x0;
    Double_t      timezero        = 0;   //TO BE FIXED
    AliAODVertex  *vV0FromCascade = 0x0;
    AliAODv0      *aodV0          = 0x0;
    AliAODcascade *aodCascade     = 0x0;
    
    // RefArray to store the mapping between esd track number and newly created AOD-Track
    TRefArray   *aodTrackRefs = NULL;
    if (nTracks > 0) aodTrackRefs = new TRefArray(nTracks);

    // RefArray to store a mapping between esd V0 number and newly created AOD-Vertex V0
    TRefArray   *aodV0VtxRefs = NULL;
    if (nV0s > 0) aodV0VtxRefs = new TRefArray(nV0s);

    // RefArray to store the mapping between esd V0 number and newly created AOD-V0
    TRefArray   *aodV0Refs = NULL;
    if (nV0s > 0) aodV0Refs = new TRefArray(nV0s);

    



    // Array to take into account the tracks already added to the AOD
    Bool_t * usedTrack = NULL;
    if (nTracks>0) {
	usedTrack = new Bool_t[nTracks];
	for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) usedTrack[iTrack]=kFALSE;
    }
    // Array to take into account the V0s already added to the AOD (V0 within cascades)
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
    TClonesArray &vertices = *(AODEvent()->GetVertices());
    Int_t jVertices=0;
    
    // Access to the AOD container of tracks
    TClonesArray &tracks = *(AODEvent()->GetTracks());
    Int_t jTracks=0; 
    
    // Access to the AOD container of Cascades
    TClonesArray &cascades = *(AODEvent()->GetCascades());
    Int_t jCascades=0;

    // Access to the AOD container of V0s
    TClonesArray &v0s = *(AODEvent()->GetV0s());
    Int_t jV0s=0;
    
    // Add primary vertex. The primary tracks will be defined
    // after the loops on the composite objects (V0, cascades, kinks)
    const AliESDVertex *vtx = esd->GetPrimaryVertex();
    
    vtx->GetXYZ(pos); // position
    vtx->GetCovMatrix(covVtx); //covariance matrix
    
    AliAODVertex * primary = new(vertices[jVertices++])
	AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);
    primary->SetName(vtx->GetName());
    primary->SetTitle(vtx->GetTitle());

    TString vtitle = vtx->GetTitle();
    if (!vtitle.Contains("VertexerTracks")) 
      primary->SetNContributors(vtx->GetNContributors());

    if (fDebug > 0) primary->Print();

    // Create vertices starting from the most complex objects
    Double_t chi2 = 0.;
    
    // Cascades (Modified by A.Maire - February 2009)
    for (Int_t nCascade = 0; nCascade < nCascades; ++nCascade) {

	if (fDebug > 1) 
	printf("\n ******** Cascade number : %d/%d *********\n", nCascade, nCascades);

	// 0- Preparation
	//
	AliESDcascade *esdCascade = esd->GetCascade(nCascade);
		Int_t  idxPosFromV0Dghter  = esdCascade->GetPindex();
		Int_t  idxNegFromV0Dghter  = esdCascade->GetNindex();
		Int_t  idxBachFromCascade  = esdCascade->GetBindex();
	
	AliESDtrack  *esdCascadePos  = esd->GetTrack( idxPosFromV0Dghter);
	AliESDtrack  *esdCascadeNeg  = esd->GetTrack( idxNegFromV0Dghter);
	AliESDtrack  *esdCascadeBach = esd->GetTrack( idxBachFromCascade);

	// Identification of the V0 within the esdCascade (via both daughter track indices)
	AliESDv0 * currentV0   = 0x0;
	Int_t      idxV0FromCascade = -1;
	
	for (Int_t iV0=0; iV0<nV0s; ++iV0) {
	
			 currentV0 = esd->GetV0(iV0);
		Int_t posCurrentV0 = currentV0->GetPindex();
		Int_t negCurrentV0 = currentV0->GetNindex();
	
		if (posCurrentV0==idxPosFromV0Dghter && negCurrentV0==idxNegFromV0Dghter) {
		idxV0FromCascade = iV0;
		break;
		}
	}

	if(idxV0FromCascade < 0){
		printf("Cascade - no matching for the V0 (index V0 = -1) ! Skip ... \n");
		continue;
	}// a priori, useless check, but safer ... in case of pb with tracks "out of bounds"

	if (fDebug > 1) 
		printf("Cascade %d - V0fromCascade ind : %d/%d\n", nCascade, idxV0FromCascade, nV0s);

	AliESDv0 *esdV0FromCascade   = esd->GetV0(idxV0FromCascade);
	

	// 1 - Cascade selection 
	
	//	AliESDVertex *esdPrimVtx = new AliESDVertex(*(esd->GetPrimaryVertex()));
	// 	TList cascadeObjects;
	// 	cascadeObjects.AddAt(esdV0FromCascade, 0);
	// 	cascadeObjects.AddAt(esdCascadePos,    1);
	// 	cascadeObjects.AddAt(esdCascadeNeg,    2);
	// 	cascadeObjects.AddAt(esdCascade,       3);
	// 	cascadeObjects.AddAt(esdCascadeBach,   4);
	// 	cascadeObjects.AddAt(esdPrimVtx,       5);
	// 
	// 	UInt_t selectCascade = 0;
	// 	if (fCascadeFilter) {
	// 	  // selectCascade = fCascadeFilter->IsSelected(&cascadeObjects); 
	// 	  	// FIXME AliESDCascadeCuts to be implemented ...
	// 
	// 		// Here we may encounter a moot point at the V0 level 
	// 		// between the cascade selections and the V0 ones :
	// 		// the V0 selected along with the cascade (secondary V0) may 
	// 		// usually be removed from the dedicated V0 selections (prim V0) ...
	// 		// -> To be discussed !
	// 
	// 	  // this is a little awkward but otherwise the 
	// 	  // list wants to access the pointer (delete it) 
	// 	  // again when going out of scope
	// 	  delete cascadeObjects.RemoveAt(5); // esdPrimVtx created via copy construct
	// 	  esdPrimVtx = 0;
	// 	  if (!selectCascade) 
	// 	    continue;
	// 	}
	// 	else{
	// 	  delete cascadeObjects.RemoveAt(5); // esdPrimVtx created via copy construct
	// 	  esdPrimVtx = 0;
	// 	}

	// 2 - Add the cascade vertex
	
	esdCascade->GetXYZcascade(pos[0], pos[1], pos[2]);
	esdCascade->GetPosCovXi(covVtx);
	chi2 = esdCascade->GetChi2Xi(); 
	
	AliAODVertex *vCascade = new(vertices[jVertices++]) AliAODVertex( pos,
									  covVtx,
									  chi2, // FIXME = Chi2/NDF will be needed
									  primary,
									  nCascade, // id
									  AliAODVertex::kCascade);
	primary->AddDaughter(vCascade);

	if (fDebug > 2) {
		printf("---- Cascade / Cascade Vertex (AOD) : \n");
		vCascade->Print();
	}


	// 3 - Add the bachelor track from the cascade
	
	if (!usedTrack[idxBachFromCascade]) {

	esdCascadeBach->GetPxPyPz(momBach);
	esdCascadeBach->GetXYZ(pos);
	esdCascadeBach->GetCovarianceXYZPxPyPz(covTr);
	esdCascadeBach->GetESDpid(pid);

	    usedTrack[idxBachFromCascade] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadeBach);
	    if(mcH)mcH->SelectParticle(esdCascadeBach->GetLabel());
	    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdCascadeBach->GetID(),
							  esdCascadeBach->GetLabel(), 
							  momBach, 
							  kTRUE,
							  pos,
							  kFALSE, // Why kFALSE for "isDCA" ? FIXME
							  covTr, 
							  (Short_t)esdCascadeBach->GetSign(),
							  esdCascadeBach->GetITSClusterMap(), 
							  pid,
							  vCascade,
							  kTRUE,  // usedForVtxFit = kFALSE ? FIXME
							  vtx->UsesTrack(esdCascadeBach->GetID()),
							  AliAODTrack::kSecondary,
							  selectInfo);
	    aodTrack->SetTPCClusterMap(esdCascadeBach->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdCascadeBach->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadeBach));
	    aodTrackRefs->AddAt(aodTrack,idxBachFromCascade);
	    
	    if (esdCascadeBach->GetSign() > 0) nPosTracks++;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdCascadeBach->GetStatus());
            SetAODPID(esdCascadeBach,aodTrack,detpid,timezero,esd->GetMagneticField());
	}
	else {
	    aodTrack = dynamic_cast<AliAODTrack*>( aodTrackRefs->At(idxBachFromCascade) );
	}

	vCascade->AddDaughter(aodTrack);

	if (fDebug > 4) {
		printf("---- Cascade / bach dghter : \n");
		aodTrack->Print();
	}
	

	// 4 - Add the V0 from the cascade. 
	// = V0vtx + both pos and neg daughter tracks + the aodV0 itself
	//

	if ( !usedV0[idxV0FromCascade] ) {
	// 4.A - if VO structure hasn't been created yet
		
		// 4.A.1 - Create the V0 vertex of the cascade
	
		esdV0FromCascade->GetXYZ(pos[0], pos[1], pos[2]);
		esdV0FromCascade->GetPosCov(covVtx);
		chi2 = esdV0FromCascade->GetChi2V0();  // = chi2/NDF since NDF = 2*2-3 ?
			
		vV0FromCascade = new(vertices[jVertices++]) AliAODVertex(pos,
									covVtx,
									chi2,
									vCascade,
									idxV0FromCascade, //id of ESDv0
									AliAODVertex::kV0);
		// Note:
		//    one V0 can be used by several cascades.
		// So, one AOD V0 vtx can have several parent vtx.
		// This is not directly allowed by AliAODvertex.
		// Setting the parent vtx (here = param "vCascade") doesn't lead to a crash
		// but to a problem of consistency within AODEvent.
		// -> See below paragraph 4.B, for the proposed treatment of such a case.

		// Add the vV0FromCascade to the aodVOVtxRefs
		aodV0VtxRefs->AddAt(vV0FromCascade,idxV0FromCascade);
		
		
		// 4.A.2 - Add the positive tracks from the V0
		
		esdCascadePos->GetPxPyPz(momPos);
		esdCascadePos->GetXYZ(pos);
		esdCascadePos->GetCovarianceXYZPxPyPz(covTr);
		esdCascadePos->GetESDpid(pid);
		
	
		if (!usedTrack[idxPosFromV0Dghter]) {
		usedTrack[idxPosFromV0Dghter] = kTRUE;
	
		UInt_t selectInfo = 0;
		if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadePos);
		if(mcH) mcH->SelectParticle(esdCascadePos->GetLabel());
		aodTrack = new(tracks[jTracks++]) AliAODTrack(  esdCascadePos->GetID(),
								esdCascadePos->GetLabel(), 
								momPos, 
								kTRUE,
								pos,
								kFALSE, // Why kFALSE for "isDCA" ? FIXME
								covTr, 
								(Short_t)esdCascadePos->GetSign(),
								esdCascadePos->GetITSClusterMap(), 
								pid,
								vV0FromCascade,
								kTRUE,  // usedForVtxFit = kFALSE ? FIXME
								vtx->UsesTrack(esdCascadePos->GetID()),
								AliAODTrack::kSecondary,
								selectInfo);
		aodTrack->SetTPCClusterMap(esdCascadePos->GetTPCClusterMap());
		aodTrack->SetTPCSharedMap (esdCascadePos->GetTPCSharedMap());
		aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadePos));
		aodTrackRefs->AddAt(aodTrack,idxPosFromV0Dghter);

		if (esdCascadePos->GetSign() > 0) nPosTracks++;
		aodTrack->ConvertAliPIDtoAODPID();
		aodTrack->SetFlags(esdCascadePos->GetStatus());
		SetAODPID(esdCascadePos,aodTrack,detpid,timezero,esd->GetMagneticField());
		}
		else {
			aodTrack = dynamic_cast<AliAODTrack*>(aodTrackRefs->At(idxPosFromV0Dghter));
		}
		vV0FromCascade->AddDaughter(aodTrack);
	
	
		// 4.A.3 - Add the negative tracks from the V0
		
		esdCascadeNeg->GetPxPyPz(momNeg);
		esdCascadeNeg->GetXYZ(pos);
		esdCascadeNeg->GetCovarianceXYZPxPyPz(covTr);
		esdCascadeNeg->GetESDpid(pid);
		
		
		if (!usedTrack[idxNegFromV0Dghter]) {
		usedTrack[idxNegFromV0Dghter] = kTRUE;
	
		UInt_t selectInfo = 0;
		if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdCascadeNeg);
		if(mcH)mcH->SelectParticle(esdCascadeNeg->GetLabel());
		aodTrack = new(tracks[jTracks++]) AliAODTrack(  esdCascadeNeg->GetID(),
								esdCascadeNeg->GetLabel(),
								momNeg,
								kTRUE,
								pos,
								kFALSE, // Why kFALSE for "isDCA" ? FIXME
								covTr, 
								(Short_t)esdCascadeNeg->GetSign(),
								esdCascadeNeg->GetITSClusterMap(), 
								pid,
								vV0FromCascade,
								kTRUE,  // usedForVtxFit = kFALSE ? FIXME
								vtx->UsesTrack(esdCascadeNeg->GetID()),
								AliAODTrack::kSecondary,
								selectInfo);
		aodTrack->SetTPCClusterMap(esdCascadeNeg->GetTPCClusterMap());
		aodTrack->SetTPCSharedMap (esdCascadeNeg->GetTPCSharedMap());
		aodTrack->SetChi2perNDF(Chi2perNDF(esdCascadeNeg));
		aodTrackRefs->AddAt(aodTrack,idxNegFromV0Dghter);
		
		if (esdCascadeNeg->GetSign() > 0) nPosTracks++;
		aodTrack->ConvertAliPIDtoAODPID();
		aodTrack->SetFlags(esdCascadeNeg->GetStatus());
		SetAODPID(esdCascadeNeg,aodTrack,detpid,timezero,esd->GetMagneticField());
		}
		else {
			aodTrack = dynamic_cast<AliAODTrack*>(aodTrackRefs->At(idxNegFromV0Dghter));
		}
	
		vV0FromCascade->AddDaughter(aodTrack);
	
			
		// 4.A.4 - Add the V0 from cascade to the V0 array

		Double_t  dcaV0Daughters      = esdV0FromCascade->GetDcaV0Daughters();
		Double_t  dcaV0ToPrimVertex   = esdV0FromCascade->GetD( esd->GetPrimaryVertex()->GetX(),
								        esd->GetPrimaryVertex()->GetY(),
								        esd->GetPrimaryVertex()->GetZ() );
		esdV0FromCascade->GetPPxPyPz( momPosAtV0vtx[0],momPosAtV0vtx[1],momPosAtV0vtx[2] ); 
		esdV0FromCascade->GetNPxPyPz( momNegAtV0vtx[0],momNegAtV0vtx[1],momNegAtV0vtx[2] ); 
	
		Double_t dcaDaughterToPrimVertex[2] = { 999., 999.}; // ..[0] = DCA in (x,y) for Pos and ..[1] = Neg
		dcaDaughterToPrimVertex[0] = TMath::Abs(esdCascadePos->GetD(	esd->GetPrimaryVertex()->GetX(),
										esd->GetPrimaryVertex()->GetY(),
										esd->GetMagneticField())        );
		dcaDaughterToPrimVertex[1] = TMath::Abs(esdCascadeNeg->GetD(	esd->GetPrimaryVertex()->GetX(),
										esd->GetPrimaryVertex()->GetY(),
										esd->GetMagneticField())        );
		
		aodV0 = new(v0s[jV0s++]) AliAODv0( vV0FromCascade, 
						   dcaV0Daughters,
						   dcaV0ToPrimVertex, 
						   momPosAtV0vtx, 
						   momNegAtV0vtx, 
						   dcaDaughterToPrimVertex); 
		// set the aod v0 on-the-fly status
		aodV0->SetOnFlyStatus(esdV0FromCascade->GetOnFlyStatus());

		// Add the aodV0 to the aodVORefs
		aodV0Refs->AddAt(aodV0,idxV0FromCascade);
		
	usedV0[idxV0FromCascade] = kTRUE;

	} else { 
	// 4.B - if V0 structure already used

		// Note :
		//    one V0 can be used by several cascades (frequent in PbPb evts) : 
		// same V0 which used but attached to different bachelor tracks
		// -> aodVORefs and aodV0VtxRefs are needed.
		// Goal : avoid a redundancy of the info in "Vertices" and "v0s" clones array.
	
		vV0FromCascade = dynamic_cast<AliAODVertex*>( aodV0VtxRefs->At(idxV0FromCascade) );
		aodV0          = dynamic_cast<AliAODv0*>    ( aodV0Refs   ->At(idxV0FromCascade) );
		 
		// - Treatment of the parent for such a "re-used" V0 :
		// Insert the cascade that reuses the V0 vertex in the lineage chain
		// Before : vV0 -> vCascade1 -> vPrimary
		//  - Hyp : cascade2 uses the same V0 as cascade1
		//  After :  vV0 -> vCascade2 -> vCascade1 -> vPrimary
		
		AliAODVertex *vCascadePreviousParent = dynamic_cast<AliAODVertex*> (vV0FromCascade->GetParent());
			vV0FromCascade->SetParent(vCascade);
			vCascade      ->SetParent(vCascadePreviousParent);

		if(fDebug > 2)	
			printf("---- Cascade / Lineage insertion\n"
				"Parent of V0 vtx                 = Cascade vtx %p\n"
				"Parent of the cascade vtx        = Cascade vtx %p\n"
				"Parent of the parent cascade vtx = Cascade vtx %p\n", 
						static_cast<void*> (vV0FromCascade->GetParent()),
						static_cast<void*> (vCascade->GetParent()),
						static_cast<void*> (vCascadePreviousParent->GetParent()) );
		
	}// end if V0 structure already used

	if (fDebug > 2) {
		printf("---- Cascade / V0 vertex: \n");
		vV0FromCascade->Print();
	}

	if (fDebug > 4) {
		printf("---- Cascade / pos dghter : \n");
			aodTrack->Print();
		printf("---- Cascade / neg dghter : \n");
			aodTrack->Print();
		printf("---- Cascade / aodV0 : \n");
			aodV0->Print();
	}

	// In any case (used V0 or not), add the V0 vertex to the cascade one.
	vCascade->AddDaughter(vV0FromCascade);	
	
		
	// 5 - Add the primary track of the cascade (if any)


	// 6 - Add the cascade to the AOD array of cascades

	Double_t dcaBachToPrimVertexXY = TMath::Abs(esdCascadeBach->GetD(esd->GetPrimaryVertex()->GetX(),
						 		 	 esd->GetPrimaryVertex()->GetY(),
								 	 esd->GetMagneticField())        );

	esdCascade->GetBPxPyPz(momBachAtCascadeVtx[0], momBachAtCascadeVtx[1], momBachAtCascadeVtx[2]);

	aodCascade = new(cascades[jCascades++]) AliAODcascade(  vCascade,
								esdCascade->Charge(),
								esdCascade->GetDcaXiDaughters(),
								-999.,
		// DCAXiToPrimVtx -> needs to be calculated   ----|
		// doesn't exist at ESD level;
		// See AODcascade::DcaXiToPrimVertex(Double, Double, Double)
								dcaBachToPrimVertexXY,
								momBachAtCascadeVtx,
								*aodV0);
	
	if (fDebug > 3) {
		printf("---- Cascade / AOD cascade : \n\n");
		aodCascade->PrintXi(primary->GetX(), primary->GetY(), primary->GetZ());
	}

    } // end of the loop on cascades

    cascades.Expand(jCascades);


    //
    // V0s
    //
    
    for (Int_t nV0 = 0; nV0 < nV0s; ++nV0) {
	
	if (usedV0[nV0]) continue; // skip if already added to the AOD
	
	AliESDv0 *v0 = esd->GetV0(nV0);
	Int_t posFromV0 = v0->GetPindex();
	Int_t negFromV0 = v0->GetNindex();
	
	// V0 selection 
	//
	AliESDVertex *esdVtx   = new AliESDVertex(*(esd->GetPrimaryVertex()));
	AliESDtrack  *esdV0Pos = esd->GetTrack(posFromV0);
	AliESDtrack  *esdV0Neg = esd->GetTrack(negFromV0);
	TList v0objects;
	v0objects.AddAt(v0,                      0);
	v0objects.AddAt(esdV0Pos,                1);
	v0objects.AddAt(esdV0Neg,                2);
	v0objects.AddAt(esdVtx,                  3);
	UInt_t selectV0 = 0;
	if (fV0Filter) {
	  selectV0 = fV0Filter->IsSelected(&v0objects);
	  // this is a little awkward but otherwise the 
	  // list wants to access the pointer (delete it) 
	  // again when going out of scope
	  delete v0objects.RemoveAt(3); // esdVtx created via copy construct
	  esdVtx = 0;
	  if (!selectV0) 
	    continue;
	}
	else{
	  delete v0objects.RemoveAt(3); // esdVtx created via copy construct
	  esdVtx = 0;
	}
    
	v0->GetXYZ(pos[0], pos[1], pos[2]);

	if (!old) {
	    chi2 = v0->GetChi2V0(); // = chi2/NDF since NDF = 2*2-3
	    v0->GetPosCov(covVtx);
	} else {
	    chi2 = -999.;
	    for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
	}


	AliAODVertex * vV0 = 
	  new(vertices[jVertices++]) AliAODVertex(pos,
						  covVtx,
						  chi2,
						  primary,
						  nV0,
						  AliAODVertex::kV0);
	primary->AddDaughter(vV0);
	

	// Add the positive tracks from the V0
	
	esdV0Pos->GetPxPyPz(momPos);
	esdV0Pos->GetXYZ(pos);
	esdV0Pos->GetCovarianceXYZPxPyPz(covTr);
	esdV0Pos->GetESDpid(pid);
	
	if (!usedTrack[posFromV0]) {
	    usedTrack[posFromV0] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdV0Pos);
	    if(mcH)mcH->SelectParticle(esdV0Pos->GetLabel());
	    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdV0Pos->GetID(),
							  esdV0Pos->GetLabel(), 
							  momPos, 
							  kTRUE,
							  pos,
							  kFALSE,
							  covTr, 
							  (Short_t)esdV0Pos->GetSign(),
							  esdV0Pos->GetITSClusterMap(), 
							  pid,
							  vV0,
							  kTRUE,  // check if this is right
							  vtx->UsesTrack(esdV0Pos->GetID()),
							  AliAODTrack::kSecondary,
							  selectInfo);
	    aodTrack->SetTPCClusterMap(esdV0Pos->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdV0Pos->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdV0Pos));
	    aodTrackRefs->AddAt(aodTrack,posFromV0);
	    //	    if (fDebug > 0) printf("-------------------Bo: pos track from original pt %.3f \n",aodTrack->Pt());
	    if (esdV0Pos->GetSign() > 0) nPosTracks++;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdV0Pos->GetStatus());
            SetAODPID(esdV0Pos,aodTrack,detpid,timezero,esd->GetMagneticField());
	}
	else {
	    aodTrack = dynamic_cast<AliAODTrack*>(aodTrackRefs->At(posFromV0));
	    //	    if (fDebug > 0) printf("-------------------Bo pos track from refArray pt %.3f \n",aodTrack->Pt());
	}
	vV0->AddDaughter(aodTrack);
    
	// Add the negative tracks from the V0
	
	esdV0Neg->GetPxPyPz(momNeg);
	esdV0Neg->GetXYZ(pos);
	esdV0Neg->GetCovarianceXYZPxPyPz(covTr);
	esdV0Neg->GetESDpid(pid);
	
	if (!usedTrack[negFromV0]) {
	    usedTrack[negFromV0] = kTRUE;
	    UInt_t selectInfo = 0;
	    if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdV0Neg);
	    if(mcH)mcH->SelectParticle(esdV0Neg->GetLabel());
	    aodTrack = new(tracks[jTracks++]) AliAODTrack(esdV0Neg->GetID(),
							  esdV0Neg->GetLabel(),
							  momNeg,
							  kTRUE,
							  pos,
							  kFALSE,
							  covTr, 
							  (Short_t)esdV0Neg->GetSign(),
							  esdV0Neg->GetITSClusterMap(), 
							  pid,
							  vV0,
							  kTRUE,  // check if this is right
							  vtx->UsesTrack(esdV0Neg->GetID()),
							  AliAODTrack::kSecondary,
							  selectInfo);
	    aodTrack->SetTPCClusterMap(esdV0Neg->GetTPCClusterMap());
	    aodTrack->SetTPCSharedMap (esdV0Neg->GetTPCSharedMap());
	    aodTrack->SetChi2perNDF(Chi2perNDF(esdV0Neg));
	    
	    aodTrackRefs->AddAt(aodTrack,negFromV0);
	    //	    if (fDebug > 0) printf("-------------------Bo: neg track from original pt %.3f \n",aodTrack->Pt());
	    if (esdV0Neg->GetSign() > 0) nPosTracks++;
	    aodTrack->ConvertAliPIDtoAODPID();
	    aodTrack->SetFlags(esdV0Neg->GetStatus());
            SetAODPID(esdV0Neg,aodTrack,detpid,timezero,esd->GetMagneticField());
	}
	else {
	    aodTrack = dynamic_cast<AliAODTrack*>(aodTrackRefs->At(negFromV0));
	    //	    if (fDebug > 0) printf("-------------------Bo neg track from refArray pt %.3f \n",aodTrack->Pt());
	}
	vV0->AddDaughter(aodTrack);


	// Add the V0 the V0 array as well
	
	Double_t  dcaV0Daughters      = v0->GetDcaV0Daughters();
	Double_t  dcaV0ToPrimVertex   = v0->GetD(esd->GetPrimaryVertex()->GetX(),
						 esd->GetPrimaryVertex()->GetY(),
						 esd->GetPrimaryVertex()->GetZ());
	v0->GetPPxPyPz(momPosAtV0vtx[0],momPosAtV0vtx[1],momPosAtV0vtx[2]); 
	v0->GetNPxPyPz(momNegAtV0vtx[0],momNegAtV0vtx[1],momNegAtV0vtx[2]); 
	
	Double_t dcaDaughterToPrimVertex[2] = { 999., 999.}; // ..[0] = DCA in (x,y) for Pos and ..[1] = Neg
	dcaDaughterToPrimVertex[0] = TMath::Abs(esdV0Pos->GetD(  esd->GetPrimaryVertex()->GetX(),
						 		 esd->GetPrimaryVertex()->GetY(),
								 	esd->GetMagneticField()) );
	dcaDaughterToPrimVertex[1] = TMath::Abs(esdV0Neg->GetD(  esd->GetPrimaryVertex()->GetX(),
						 		 esd->GetPrimaryVertex()->GetY(),
								 	esd->GetMagneticField()) );
	
	aodV0 = new(v0s[jV0s++]) AliAODv0(vV0, 
					dcaV0Daughters,
					dcaV0ToPrimVertex,
					momPosAtV0vtx,
					momNegAtV0vtx,
					dcaDaughterToPrimVertex);

	// set the aod v0 on-the-fly status
	aodV0->SetOnFlyStatus(v0->GetOnFlyStatus());
    }//End of loop on V0s 

    v0s.Expand(jV0s);	 

    if (fDebug > 0)   printf("   NAODCascades=%d / NAODV0s=%d\n", jCascades, jV0s);
    // end of V0 parts


    // Kinks: it is a big mess the access to the information in the kinks
    // The loop is on the tracks in order to find the mother and daugther of each kink
    
    
    for (Int_t iTrack=0; iTrack<nTracks; ++iTrack) {
	
	AliESDtrack * esdTrack = esd->GetTrack(iTrack);
	
	Int_t ikink = esdTrack->GetKinkIndex(0);
	
	if (ikink && nKinks) {
	    // Negative kink index: mother, positive: daughter
	    
	    // Search for the second track of the kink
	    
	    for (Int_t jTrack = iTrack+1; jTrack<nTracks; ++jTrack) {
		
		AliESDtrack * esdTrack1 = esd->GetTrack(jTrack);
		
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
		    
		    // Add the mother track if it passed primary track selection cuts
		    
		    AliAODTrack * mother = NULL;
		    
		    UInt_t selectInfo = 0;
		    if (fTrackFilter) {
			selectInfo = fTrackFilter->IsSelected(esd->GetTrack(imother));
			if (!selectInfo) continue;
		    }
		    
		    if (!usedTrack[imother]) {
			
			usedTrack[imother] = kTRUE;
			
			AliESDtrack *esdTrackM = esd->GetTrack(imother);
			esdTrackM->GetPxPyPz(p);
			esdTrackM->GetXYZ(pos);
			esdTrackM->GetCovarianceXYZPxPyPz(covTr);
			esdTrackM->GetESDpid(pid);
			if(mcH)mcH->SelectParticle(esdTrackM->GetLabel());
			mother = 
			    new(tracks[jTracks++]) AliAODTrack(esdTrackM->GetID(),
							       esdTrackM->GetLabel(),
							       p,
							       kTRUE,
							       pos,
							       kFALSE,
							       covTr, 
							       (Short_t)esdTrackM->GetSign(),
							       esdTrackM->GetITSClusterMap(), 
							       pid,
							       primary,
							       kTRUE, // check if this is right
							       vtx->UsesTrack(esdTrack->GetID()),
							       AliAODTrack::kPrimary,
							       selectInfo);
			mother->SetTPCClusterMap(esdTrackM->GetTPCClusterMap());
			mother->SetTPCSharedMap (esdTrackM->GetTPCSharedMap());
			mother->SetChi2perNDF(Chi2perNDF(esdTrackM));

			aodTrackRefs->AddAt(mother, imother);
			
			if (esdTrackM->GetSign() > 0) nPosTracks++;
			mother->SetFlags(esdTrackM->GetStatus());
			mother->ConvertAliPIDtoAODPID();
			primary->AddDaughter(mother);
			mother->ConvertAliPIDtoAODPID();
                        SetAODPID(esdTrackM,mother,detpid,timezero,esd->GetMagneticField());
		    }
		    else {
//			cerr << "Error: event " << esd->GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
//			     << " track " << imother << " has already been used!" << endl;
		  }
		    
		    // Add the kink vertex
		    AliESDkink * kink = esd->GetKink(TMath::Abs(ikink)-1);
		    
		    AliAODVertex * vkink = 
			new(vertices[jVertices++]) AliAODVertex(kink->GetPosition(),
								NULL,
								0.,
								mother,
								esdTrack->GetID(),  // This is the track ID of the mother's track!
								AliAODVertex::kKink);
		    // Add the daughter track
		  
		    AliAODTrack * daughter = NULL;
		    
		    if (!usedTrack[idaughter]) {
			
			usedTrack[idaughter] = kTRUE;
			
			AliESDtrack *esdTrackD = esd->GetTrack(idaughter);
			esdTrackD->GetPxPyPz(p);
			esdTrackD->GetXYZ(pos);
			esdTrackD->GetCovarianceXYZPxPyPz(covTr);
			esdTrackD->GetESDpid(pid);
			selectInfo = 0;
			if (fTrackFilter) selectInfo = fTrackFilter->IsSelected(esdTrackD);
			if(mcH)mcH->SelectParticle(esdTrackD->GetLabel());
			daughter = 
			    new(tracks[jTracks++]) AliAODTrack(esdTrackD->GetID(),
							       esdTrackD->GetLabel(),
							       p,
							       kTRUE,
							       pos,
							       kFALSE,
							       covTr, 
							       (Short_t)esdTrackD->GetSign(),
							       esdTrackD->GetITSClusterMap(), 
							       pid,
							       vkink,
							       kTRUE, // check if this is right
							       vtx->UsesTrack(esdTrack->GetID()),
							       AliAODTrack::kSecondary,
							       selectInfo);
			daughter->SetTPCClusterMap(esdTrackD->GetTPCClusterMap());
			daughter->SetTPCSharedMap (esdTrackD->GetTPCSharedMap());
			aodTrackRefs->AddAt(daughter, idaughter);
			
			if (esdTrackD->GetSign() > 0) nPosTracks++;
			daughter->SetFlags(esdTrackD->GetStatus());
			daughter->ConvertAliPIDtoAODPID();
			vkink->AddDaughter(daughter);
			daughter->ConvertAliPIDtoAODPID();
                       SetAODPID(esdTrackD,daughter,detpid,timezero,esd->GetMagneticField());
		    }
		    else {
//			cerr << "Error: event " << esd->GetEventNumberInFile() << " kink " << TMath::Abs(ikink)-1
//			     << " track " << idaughter << " has already been used!" << endl;
		    }
		}
	    }
	}      
    }
    
  
    // Tracks (primary and orphan)

    if (fDebug > 0) printf("NUMBER OF ESD TRACKS %5d\n", nTracks);
    
    for (Int_t nTrack = 0; nTrack < nTracks; ++nTrack) {
	
	
	if (usedTrack[nTrack]) continue;

	AliESDtrack *esdTrack = esd->GetTrack(nTrack);
	UInt_t selectInfo = 0;
	//
	// Track selection
	if (fTrackFilter) {
	    selectInfo = fTrackFilter->IsSelected(esdTrack);
	    if (!selectInfo && !vtx->UsesTrack(esdTrack->GetID())) continue;
	}
	
	//
	esdTrack->GetPxPyPz(p);
	esdTrack->GetXYZ(pos);
	esdTrack->GetCovarianceXYZPxPyPz(covTr);
	esdTrack->GetESDpid(pid);
	if(mcH)mcH->SelectParticle(esdTrack->GetLabel());
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
								vtx->UsesTrack(esdTrack->GetID()),
								AliAODTrack::kPrimary, 
								selectInfo)
			     );
	aodTrack->SetTPCClusterMap(esdTrack->GetTPCClusterMap());
	aodTrack->SetTPCSharedMap (esdTrack->GetTPCSharedMap());
	aodTrack->SetChi2perNDF(Chi2perNDF(esdTrack));

	aodTrackRefs->AddAt(aodTrack, nTrack);

	
	if (esdTrack->GetSign() > 0) nPosTracks++;
	aodTrack->SetFlags(esdTrack->GetStatus());
	aodTrack->ConvertAliPIDtoAODPID();
	SetAODPID(esdTrack,aodTrack,detpid,timezero,esd->GetMagneticField());
    } // end of loop on tracks
    
    // Update number of AOD tracks in header at the end of track loop (M.G.)
    header->SetRefMultiplicity(jTracks);
    header->SetRefMultiplicityPos(nPosTracks);
    header->SetRefMultiplicityNeg(jTracks - nPosTracks);
    if (fDebug > 0) 
      printf("   NAODTRACKS=%d  NPOS=%d  NNEG=%d\n", jTracks, nPosTracks, jTracks - nPosTracks);
    // Do not shrink the array of tracks - other filters may add to it (M.G)
    //   tracks.Expand(jTracks); // remove 'empty slots' due to unwritten tracks
  
    // Access to the AOD container of PMD clusters
    TClonesArray &pmdClusters = *(AODEvent()->GetPmdClusters());
    Int_t jPmdClusters=0;
  
    for (Int_t iPmd = 0; iPmd < nPmdClus; ++iPmd) {
      // file pmd clusters, to be revised!
      AliESDPmdTrack *pmdTrack = esd->GetPmdTrack(iPmd);
      Int_t nLabel = 0;
      Int_t *label = 0x0;
      Double_t posPmd[3] = { pmdTrack->GetClusterX(), pmdTrack->GetClusterY(), pmdTrack->GetClusterZ()};
      Double_t pidPmd[13] = { 0.}; // to be revised!
      // type not set!
      // assoc cluster not set
      new(pmdClusters[jPmdClusters++]) AliAODPmdCluster(iPmd, nLabel, label, pmdTrack->GetClusterADC(), posPmd, pidPmd);
    }

    // Access to the AOD container of clusters
    TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
    Int_t jClusters=0;
 
    for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {

      AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);

      Int_t id        = cluster->GetID();
      Int_t nLabel    = cluster->GetNLabels();
      TArrayI* labels = cluster->GetLabels();
      Int_t *label = 0;
      if (labels){
	label = (cluster->GetLabels())->GetArray();
	for(int i = 0;i < labels->GetSize();++i){
	  if(mcH)mcH->SelectParticle(label[i]);
	}
      }     

      Float_t energy = cluster->E();
      cluster->GetPosition(posF);
      Char_t ttype = AliAODCluster::kUndef; 

      if (cluster->GetClusterType() == AliESDCaloCluster::kPHOSCluster) {
	ttype=AliAODCluster::kPHOSNeutral;
      } 
      else if (cluster->GetClusterType() == AliESDCaloCluster::kEMCALClusterv1) {
	ttype = AliAODCluster::kEMCALClusterv1;
      }

      
      AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) AliAODCaloCluster(id,
											nLabel,
											label,
											energy,
											posF,
											NULL,
											ttype);
      
      caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
				  cluster->GetClusterDisp(),
				  cluster->GetM20(), cluster->GetM02(),
				  cluster->GetEmcCpvDistance(),  
				  cluster->GetNExMax(),cluster->GetTOF()) ;

      caloCluster->SetPIDFromESD(cluster->GetPid());
      caloCluster->SetNCells(cluster->GetNCells());
      caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
      caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());

      TArrayI* matchedT = 	cluster->GetTracksMatched();
      if (nTracks>0 && matchedT && cluster->GetTrackMatched() >= 0) {	
	for (Int_t im = 0; im < matchedT->GetSize(); im++) {
	    Int_t iESDtrack = matchedT->At(im);;
	    if (aodTrackRefs->At(iESDtrack) != 0) {
		caloCluster->AddTrackMatched((AliAODTrack*)aodTrackRefs->At(iESDtrack));
	    }
	}
      }
      
    } 
    caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
    // end of loop on calo clusters

    // fill EMCAL cell info
    if (esd->GetEMCALCells()) { // protection against missing ESD information
      AliESDCaloCells &esdEMcells = *(esd->GetEMCALCells());
      Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
      
      AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
      aodEMcells.CreateContainer(nEMcell);
      aodEMcells.SetType(AliAODCaloCells::kEMCAL);
      for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
	aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell));
      }
      aodEMcells.Sort();
    }

    // fill PHOS cell info
    if (esd->GetPHOSCells()) { // protection against missing ESD information
      AliESDCaloCells &esdPHcells = *(esd->GetPHOSCells());
      Int_t nPHcell = esdPHcells.GetNumberOfCells() ;
      
      AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
      aodPHcells.CreateContainer(nPHcell);
      aodPHcells.SetType(AliAODCaloCells::kPHOS);
      for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
	aodPHcells.SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell));
      }
      aodPHcells.Sort();
    }

    // tracklets    
    AliAODTracklets &SPDTracklets = *(AODEvent()->GetTracklets());
    const AliMultiplicity *mult = esd->GetMultiplicity();
    if (mult) {
      if (mult->GetNumberOfTracklets()>0) {
	SPDTracklets.CreateContainer(mult->GetNumberOfTracklets());

	for (Int_t n=0; n<mult->GetNumberOfTracklets(); n++) {
	  if(mcH){
	    mcH->SelectParticle(mult->GetLabel(n, 0));
	    mcH->SelectParticle(mult->GetLabel(n, 1));
	  }
	  SPDTracklets.SetTracklet(n, mult->GetTheta(n), mult->GetPhi(n), mult->GetDeltaPhi(n), mult->GetLabel(n, 0),mult->GetLabel(n, 1));
	}
      }
    } else {
      //Printf("ERROR: AliMultiplicity could not be retrieved from ESD");
    }

    delete [] usedKink;
    delete [] usedV0;
    delete [] usedTrack;    
    delete    aodV0Refs;
    delete    aodV0VtxRefs;
    delete    aodTrackRefs;

    return;
}


void AliAnalysisTaskESDfilter::SetAODPID(AliESDtrack *esdtrack, AliAODTrack *aodtrack, AliAODPid *detpid, Double_t timezero, Double_t bfield)
{
  //
  // Setter for the raw PID detector signals
  //

  // Save PID object for candidate electrons
    Bool_t pidSave = kFALSE;
    if (fTrackFilter) {
	Bool_t selectInfo = fTrackFilter->IsSelected("Electrons");
	if (selectInfo)  pidSave = kTRUE;
    }


    // Tracks passing pt cut 
    if(esdtrack->Pt()>fHighPthreshold) {
	pidSave = kTRUE;
    } else {
	if(fPtshape){
	    if(esdtrack->Pt()> fPtshape->GetXmin()){
		Double_t y = fPtshape->Eval(esdtrack->Pt())/fPtshape->Eval(fHighPthreshold);
		if(gRandom->Rndm(0)<1./y){
		    pidSave = kTRUE;
		}//end rndm
	    }//end if p < pmin
	}//end if p function
    }// end else

    if (pidSave) {
      if(!aodtrack->GetDetPid()){// prevent memory leak when calling SetAODPID twice for the same track
	detpid = new AliAODPid();
	SetDetectorRawSignals(detpid,esdtrack,timezero, bfield);
	aodtrack->SetDetPID(detpid);
      }
    }
}

void AliAnalysisTaskESDfilter::SetDetectorRawSignals(AliAODPid *aodpid, AliESDtrack *track, Double_t timezero, Double_t bfield)
{
//
//assignment of the detector signals (AliXXXesdPID inspired)
//
 if(!track){
 AliInfo("no ESD track found. .....exiting");
 return;
 }
 // TPC momentum
 const AliExternalTrackParam *in=track->GetInnerParam();
 if (in) {
   aodpid->SetTPCmomentum(in->GetP());
 }else{
   aodpid->SetTPCmomentum(-1.);
 }


 aodpid->SetITSsignal(track->GetITSsignal());
 aodpid->SetTPCsignal(track->GetTPCsignal());

 //n TRD planes = 6
 Int_t nslices = track->GetNumberOfTRDslices()*6;
 Double_t *trdslices = new Double_t[nslices];
 for(Int_t iSl =0; iSl < track->GetNumberOfTRDslices(); iSl++) {
   for(Int_t iPl =0; iPl<6; iPl++) trdslices[iPl*track->GetNumberOfTRDslices()+iSl] = track->GetTRDslice(iPl,iSl);
 }
 
//TRD momentum
 for(Int_t iPl=0;iPl<6;iPl++){
   Double_t trdmom=track->GetTRDmomentum(iPl);
   aodpid->SetTRDmomentum(iPl,trdmom);
 }

 aodpid->SetTRDsignal(track->GetNumberOfTRDslices()*6,trdslices);
 Double_t times[AliAODPid::kSPECIES]; track->GetIntegratedTimes(times);
 aodpid->SetIntegratedTimes(times);

 aodpid->SetTOFsignal(track->GetTOFsignal()-timezero); // to be fixed
 aodpid->SetHMPIDsignal(track->GetHMPIDsignal());

 //Extrapolate track to EMCAL surface for AOD-level track-cluster matching
 Double_t emcpos[3] = {0.,0.,0.};
 Double_t emcmom[3] = {0.,0.,0.};
 aodpid->SetEMCALPosition(emcpos);
 aodpid->SetEMCALMomentum(emcmom);

 AliExternalTrackParam *outerparam = (AliExternalTrackParam*)track->GetOuterParam();
 if(!outerparam) return;

 //To be replaced by call to AliEMCALGeoUtils when the class becomes available
 Double_t radius = 441.0; //[cm] EMCAL radius +13cm

 Bool_t okpos = outerparam->GetXYZAt(radius,bfield,emcpos);
 Bool_t okmom = outerparam->GetPxPyPzAt(radius,bfield,emcmom);
 if(!(okpos && okmom)) return;

 aodpid->SetEMCALPosition(emcpos);
 aodpid->SetEMCALMomentum(emcmom);

}

Double_t  AliAnalysisTaskESDfilter::Chi2perNDF(AliESDtrack* track)
{
    // Calculate chi2 per ndf for track
    Int_t  nClustersTPC = track->GetTPCNcls();

    if ( nClustersTPC > 5) {
       return (track->GetTPCchi2()/Float_t(nClustersTPC - 5));
    } else {
       return (-1.);
    }
 }



void AliAnalysisTaskESDfilter::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisESDfilter: Terminate() \n");
}

void  AliAnalysisTaskESDfilter::PrintMCInfo(AliStack *pStack,Int_t label){
  if(!pStack)return;
  label = TMath::Abs(label);
  TParticle *part = pStack->Particle(label);
  Printf("########################");
  Printf("%s:%d %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,label,part->GetUniqueID(),part->GetPdgCode(),part->P());
  part->Print();
  TParticle* mother = part;
  Int_t imo = part->GetFirstMother();
  Int_t nprim = pStack->GetNprimary();
  //  while((imo >= nprim) && (mother->GetUniqueID() == 4)) {
  while((imo >= nprim)) {
    mother =  pStack->Particle(imo);
    Printf("Mother %s:%d Label %d UniqueID %d PDG %d P %3.3f",(char*)__FILE__,__LINE__,imo,mother->GetUniqueID(),mother->GetPdgCode(),mother->P());
    mother->Print();
    imo =  mother->GetFirstMother();
  }
  Printf("########################");
}
