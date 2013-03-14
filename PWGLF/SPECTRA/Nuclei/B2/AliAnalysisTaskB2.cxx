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

// Analysis task for B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

#include <AliLog.h>

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>

#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrackCuts.h>
#include <AliExternalTrackParam.h>

#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliMCVertex.h>
#include <AliStack.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TMCProcess.h>

#include <AliTriggerAnalysis.h> // for offline signals
#include <AliCentrality.h>
#include <AliESDUtils.h>

#include <AliESDpid.h>

#include <TH1D.h>
#include <TH2D.h>

#include "AliLnID.h"
#include "AliLnHistoMap.h"
#include "AliAnalysisTaskB2.h"

ClassImp(AliAnalysisTaskB2)

AliAnalysisTaskB2::AliAnalysisTaskB2()
: AliAnalysisTask()
, fSpecies("Proton")
, fPartCode(AliPID::kProton)
, fHeavyIons(0)
, fSimulation(0)
, fMultTrigger(0)
, fCentTrigger(0)
, fTriggerFired(0)
, fGoodVertex(0)
, fPileUpEvent(0)
, fV0AND(0)
, fNoFastOnly(0)
, fMinKNOmult(-10)
, fMaxKNOmult(100000)
, fMinCentrality(0)
, fMaxCentrality(100)
, fNch(0)
, fNtrk(0)
, fMeanNtrk(1)
, fKNOmult(1)
, fMinVx(-1)
, fMaxVx(1)
, fMinVy(-1)
, fMaxVy(1)
, fMinVz(-10)
, fMaxVz(10)
, fMinDCAxy(-1)
, fMaxDCAxy(1)
, fMinDCAz(-2)
, fMaxDCAz(2)
, fMaxNSigma(3)
, fMinEta(-0.8)
, fMaxEta(0.8)
, fMinY(-0.5)
, fMaxY(0.5)
, fMCevent(0)
, fESDevent(0)
, fHistoMap(0)
, fESDtrackCuts(0)
, fLnID(0)
, fMaxNSigmaITS(3)
, fMaxNSigmaTPC(3)
, fMaxNSigmaTOF(3)
, fTrigAna(0)
, fESDpid(0)
, fIsPidOwner(0)
, fTimeZeroType(AliESDpid::kTOF_T0)
, fTOFmatch(0)
, fMinM2(2.)
, fMaxM2(6.)

{
//
// Default constructor
//
	AliLog::SetGlobalLogLevel(AliLog::kFatal);
	
	fTrigAna = new AliTriggerAnalysis();
}

AliAnalysisTaskB2::AliAnalysisTaskB2(const char* name)
: AliAnalysisTask(name,"")
, fSpecies("Proton")
, fPartCode(AliPID::kProton)
, fHeavyIons(0)
, fSimulation(0)
, fMultTrigger(0)
, fCentTrigger(0)
, fTriggerFired(0)
, fGoodVertex(0)
, fPileUpEvent(0)
, fV0AND(0)
, fNoFastOnly(0)
, fMinKNOmult(-10)
, fMaxKNOmult(100000)
, fMinCentrality(-1)
, fMaxCentrality(100)
, fNtrk(0)
, fMeanNtrk(1)
, fKNOmult(1)
, fMinVx(-1)
, fMaxVx(1)
, fMinVy(-1)
, fMaxVy(1)
, fMinVz(-10)
, fMaxVz(10)
, fMinDCAxy(-1)
, fMaxDCAxy(1)
, fMinDCAz(-2)
, fMaxDCAz(2)
, fMaxNSigma(3)
, fMinEta(-0.8)
, fMaxEta(0.8)
, fMinY(-0.5)
, fMaxY(0.5)
, fMCevent(0)
, fESDevent(0)
, fHistoMap(0)
, fESDtrackCuts(0)
, fLnID(0)
, fMaxNSigmaITS(3)
, fMaxNSigmaTPC(3)
, fMaxNSigmaTOF(3)
, fTrigAna(0)
, fESDpid(0)
, fIsPidOwner(0)
, fTimeZeroType(AliESDpid::kTOF_T0)
, fTOFmatch(0)
, fMinM2(2.)
, fMaxM2(6.)

{
//
// Constructor
//
	DefineInput(0, TChain::Class());
	DefineOutput(0, AliLnHistoMap::Class());
	
	//kFatal, kError, kWarning, kInfo, kDebug, kMaxType
	AliLog::SetGlobalLogLevel(AliLog::kFatal);
	
	fTrigAna = new AliTriggerAnalysis();
}

void AliAnalysisTaskB2::ConnectInputData(Option_t *)
{
//
// Connect input data
//
	TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
	if(!tree)
	{
		this->Error("ConnectInputData", "could not read chain from input slot 0");
		return;
	}
	
	// Get the pointer to the existing analysis manager via the static access method.
	AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr)
	{
		this->Error("ConnectInputData", "could not get analysis manager");
		return;
	}
	
	// cast type AliVEventHandler
	AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
	if (!esdH)
	{
		this->Error("ConnectInputData", "could not get ESDInputHandler");
		return;
	}
	
	// Get pointer to esd event from input handler
	fESDevent = esdH->GetEvent();
	
	// PID object for TOF
	fESDpid = esdH->GetESDpid();
	
	if(!fESDpid)
	{ //in case of no Tender attached
		fESDpid = new AliESDpid();
		fIsPidOwner = kTRUE;
	}
	
	if(!fSimulation) return;
	
	AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (mgr->GetMCtruthEventHandler());
	if (!mcH)
	{
		this->Error("ConnectInputData", "could not get AliMCEventHandler");
		return;
	}
	
	fMCevent = mcH->MCEvent();
	if (!fMCevent)
	{
		this->Error("ConnectInputData", "could not get MC fLnEvent");
		return;
	}
}

void AliAnalysisTaskB2::CreateOutputObjects()
{
//
// Create output objects
//
	if(fHistoMap == 0) exit(1); // should be created somewhere else
	
	PostData(0, fHistoMap);
}

AliAnalysisTaskB2::~AliAnalysisTaskB2()
{
//
// Default destructor
//
	delete fHistoMap;
	delete fESDtrackCuts;
	delete fLnID;
	
	delete fTrigAna;
	
	if(fIsPidOwner) delete fESDpid;
}

void AliAnalysisTaskB2::SetParticleSpecies(const TString& species)
{
//
// set the particle species and the AliPID code
//
	fSpecies = species;
	fPartCode = this->GetPidCode(species);
}

void AliAnalysisTaskB2::Exec(Option_t* )
{
//
// Execute analysis for the current event
//
	if(fESDevent == 0)
	{
		this->Error("Exec", "AliESDEvent not available");
		return;
	}
	
	if(fESDtrackCuts == 0)
	{
		this->Error("Exec", "ESD track cuts not set");
		return;
	}
	
	if(fLnID == 0)
	{
		this->Error("Exec", "PID not set");
		return;
	}
	
	fMultTrigger  = kFALSE;
	fCentTrigger  = kFALSE;
	fTriggerFired = kFALSE;
	fGoodVertex   = kFALSE;
	fPileUpEvent  = kFALSE;
	
	// --------- multiplicity and centrality ------------------
	
	fNtrk = AliESDtrackCuts::GetReferenceMultiplicity(fESDevent);
	if(fSimulation) fNch = this->GetChargedMultiplicity(0.5);
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Event_Ntrk"))->Fill(fNtrk);
	
	fKNOmult = fNtrk/fMeanNtrk;
	
	if( (fKNOmult >= fMinKNOmult) && (fKNOmult < fMaxKNOmult) ) fMultTrigger = kTRUE;
	
	if(fHeavyIons)
	{
		AliCentrality* esdCent = fESDevent->GetCentrality();
		
		Float_t centrality = esdCent->GetCentralityPercentile("V0M");
		
		if((centrality >= fMinCentrality) && (centrality < fMaxCentrality)) fCentTrigger = kTRUE;
		
		Float_t v0ScaMult;
		Float_t v0Mult = AliESDUtils::GetCorrV0(fESDevent,v0ScaMult);
		
		((TH1D*)fHistoMap->Get(fSpecies + "_V0_Mult"))->Fill(v0Mult);
		((TH1D*)fHistoMap->Get(fSpecies + "_V0_Scaled_Mult"))->Fill(v0ScaMult);
	}
	
	// ----------------- trigger --------------------
	
	AliInputEventHandler* eventH = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(eventH == 0)
	{
		this->Error("Exec", "could not get AliInputEventHandler");
		return;
	}
	
	UInt_t triggerBits = eventH->IsEventSelected();
	
	if(fHeavyIons)
	{
		fTriggerFired = ( this->IsMB(triggerBits) && fCentTrigger );
	}
	else
	{
		fTriggerFired = this->IsMB(triggerBits);
		if(fNoFastOnly) fTriggerFired = ( fTriggerFired && !this->IsFastOnly(triggerBits) );
		if(fV0AND) fTriggerFired = ( fTriggerFired && this->IsV0AND() );
		
		fTriggerFired = ( fTriggerFired && fMultTrigger );
	}
	
	// --------------------- vertex --------------
	
	const AliESDVertex* vtx = fESDevent->GetPrimaryVertex(); // best primary vertex
	
	if(vtx->GetStatus())
	{
		fGoodVertex=kTRUE;
		
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_YX"))->Fill(vtx->GetX(), vtx->GetY());
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_YZ"))->Fill(vtx->GetZ(), vtx->GetY());
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_XZ"))->Fill(vtx->GetZ(), vtx->GetX());
	}
	
	if(fESDevent->GetPrimaryVertex()->IsFromVertexerZ())
	{
		if(fESDevent->GetPrimaryVertex()->GetDispersion() > 0.02) fGoodVertex=kFALSE;
		if(fESDevent->GetPrimaryVertex()->GetZRes() > 0.25 ) fGoodVertex=kFALSE;
	}
	
	if( (vtx->GetX() <= fMinVx) || (vtx->GetX() >= fMaxVx) ) fGoodVertex=kFALSE;
	if( (vtx->GetY() <= fMinVy) || (vtx->GetY() >= fMaxVy) ) fGoodVertex=kFALSE;
	if( (vtx->GetZ() <= fMinVz) || (vtx->GetZ() >= fMaxVz) ) fGoodVertex=kFALSE;
	
	// -------------------- pile up ------------------
	
	if(fESDevent->IsPileupFromSPDInMultBins()) fPileUpEvent = kTRUE;
	
	// ---------------- event stats ----------------
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(0); // number of events
	
	if(fTriggerFired)
	{
		((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(1); // triggering events
		
		if(vtx->GetStatus())
		{
			((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(3); // valid vertex within a triggering event
			
			if((vtx->GetZ() > fMinVz) && (vtx->GetZ() < fMaxVz))
			{
				((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(4); // Vz
				
				if(   (vtx->GetX() > fMinVx) && (vtx->GetX() < fMaxVx)
				   && (vtx->GetY() > fMinVy) && (vtx->GetY() < fMaxVy))
				{
					((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(5); // Vx, Vy
					
					if(fGoodVertex)
					{
						((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(6); // VertexerZ
						
						if(!fPileUpEvent)
						{
							((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(7); // no pile-up
						}
					}
				}
			}
		}
	}
	
	// ------------- Particles and tracks for this event ---------------
	
	Int_t nParticles = 0;
	if(fSimulation)
	{
		nParticles = this->GetParticles();
	}
	
	this->GetTracks();
	
	// Post the data (input/output slots #0 already used by the base class)
	PostData(0, fHistoMap);
}

Int_t AliAnalysisTaskB2::GetParticles()
{
//
// Get particles from current event
//
	Int_t nParticles = 0;
	
	AliStack* stack = fMCevent->Stack();
	if (!stack)
	{
		AliDebug(AliLog::kWarning, "stack not available");
		return 0;
	}
	
	for (Int_t i = 0; i < fMCevent->GetNumberOfTracks(); ++i)
	{
		TParticle* iParticle = stack->Particle(i);
		
		if(!iParticle) continue;
		
		Int_t pid = fLnID->GetPID(iParticle);
		
		if(pid != fPartCode) continue;
		
		// ------- physical primaries ------------
		
		if(!stack->IsPhysicalPrimary(i)) continue;
		
		TString particle = fSpecies;
		if(iParticle->GetPDG()->Charge() < 0) particle.Prepend("Anti");
		
		Double_t genP   = iParticle->P();
		Double_t genPt  = iParticle->Pt();
		Double_t genY   = iParticle->Y();
		Double_t genPhi = iParticle->Phi();
		Double_t genEta = iParticle->Eta();
		
		// ------- multiplicity and centrality -------------
		
		if( fHeavyIons && !fCentTrigger) continue;
		if(!fHeavyIons && !fMultTrigger) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Pt"))->Fill(genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Y"))->Fill(genY);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Phi"))->Fill(genPhi);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Eta"))->Fill(genEta);
		((TH2D*)fHistoMap->Get(particle + "_Gen_Prim_PtY"))->Fill(genY, genPt);
		((TH2D*)fHistoMap->Get(particle + "_Gen_Prim_EtaY"))->Fill(genY, genEta);
		
		// ------ is within phase space? ----------
		
		if(TMath::Abs(genY) >= fMaxY) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_PhS_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_PhS_Prim_Pt"))->Fill(genPt);
		
		// ------- is from a triggering event? (same as rec.) --------
		
		if(!fTriggerFired)
		{
			((TH1D*)fHistoMap->Get(particle + "_Gen_NoTrig_Prim_Pt"))->Fill(genPt);
			continue;
		}
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Trig_Prim_Pt"))->Fill(genPt);
		
		// ------- is from a triggering event with good vertex? -------------
		
		if(!fESDevent->GetPrimaryVertex()->GetStatus())
		{
			((TH1D*)fHistoMap->Get(particle + "_Gen_NoVtx_Prim_Pt"))->Fill(genPt);
		}
		
		if(!fGoodVertex) continue;
		if(fPileUpEvent) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Nch"))->Fill(fNch);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Vtx_Prim_Pt"))->Fill(genPt);
		
		// ------ is within the geometrical acceptance? ------------
		
		if(TMath::Abs(genEta) >= fMaxEta) continue;
		
		((TH2D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_PtY"))->Fill(genY,genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Pt"))->Fill(genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Phi"))->Fill(genPhi);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Y"))->Fill(genY);
	}
	
	return nParticles;
}

Int_t AliAnalysisTaskB2::GetTracks()
{
//
// Get tracks from current event
//
	using namespace std;
	
	Int_t nTracks = 0;
	
	// --------- trigger, vertex and pile-up -----------
	
	if(!fTriggerFired) return 0;
	if(!fGoodVertex)   return 0;
	if(fPileUpEvent)   return 0;
	
	// ------------ this is a 'good' event -------------
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(2); // analyzed events
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Ana_Event_Ntrk"))->Fill(fNtrk);
	
	if(fSimulation)
	{
		((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Event_Nch_Ntrk"))->Fill(fNtrk, fNch);
	}
	
	const AliESDVertex* vtx = fESDevent->GetPrimaryVertex();
	
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_YX"))->Fill(vtx->GetX(), vtx->GetY());
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_YZ"))->Fill(vtx->GetZ(), vtx->GetY());
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_XZ"))->Fill(vtx->GetZ(), vtx->GetX());
	
	if(fHeavyIons)
	{
		Float_t v0ScaMult;
		Float_t v0Mult = AliESDUtils::GetCorrV0(fESDevent,v0ScaMult);
		
		((TH1D*)fHistoMap->Get(fSpecies + "_Ana_V0_Mult"))->Fill(v0Mult);
		((TH1D*)fHistoMap->Get(fSpecies + "_Ana_V0_Scaled_Mult"))->Fill(v0ScaMult);
	}
	
	// ------- track loop --------
	
	for(Int_t i = 0; i < fESDevent->GetNumberOfTracks(); ++i)
	{
		AliESDtrack* iTrack = fESDevent->GetTrack(i);
		
		if(!iTrack) continue;
		
		Float_t sign = 1;
		TString particle = fSpecies;
		if(iTrack->GetSign()<0 )
		{
			particle.Prepend("Anti");
			sign = -1;
		}
		
		// -------- track cuts ------------------------
		
		Double_t theta = this->GetTheta(iTrack);
		Double_t phi   = this->GetPhi(iTrack);
		
		((TH2D*)fHistoMap->Get(particle + "_Before_Phi_Theta"))->Fill(theta,phi);
		
		if(!fESDtrackCuts->AcceptTrack(iTrack)) continue;  // with next track
		if(fTOFmatch && !this->AcceptTOFtrack(iTrack)) continue; // with next track
		
		// --------------------- end track cuts -----------------------
		
		((TH2D*)fHistoMap->Get(particle + "_After_Phi_Theta"))->Fill(theta, phi);
		
		++nTracks;
		
		// charge
		
		Double_t z = 1;
		if(fPartCode>AliPID::kTriton)  z = 2;
		// impact parameters
		
		Float_t dcaxy, dcaz;
		iTrack->GetImpactParameters(dcaxy, dcaz);
		
		dcaxy *= sign;
		dcaz  *= sign;
		
		Double_t nSigmaVtx = fESDtrackCuts->GetSigmaToVertex(iTrack);
		
		// momentum
		
		Double_t pt      = iTrack->Pt()*z; // pt at DCA
		Double_t y       = this->GetRapidity(iTrack, fPartCode);
		Double_t pITS    = this->GetITSmomentum(iTrack);
		Double_t pTPC    = iTrack->GetTPCmomentum();
		Double_t pTOF    = this->GetTOFmomentum(iTrack);
		Double_t dEdxITS = iTrack->GetITSsignal();
		Double_t dEdxTPC = iTrack->GetTPCsignal();
		Int_t nPointsITS = this->GetITSnPointsPID(iTrack);
		Int_t nPointsTPC = iTrack->GetTPCsignalN();
		
		Double_t beta = 0;
		Double_t mass = 0;
		Double_t m2   = 0;
		
		Double_t simPt  = 0;
		Double_t simPhi = 0;
		Double_t simY   = 0;
		
		// --------- track cuts ------------
		
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_DCAxy"))->Fill(dcaxy);
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_DCAz"))->Fill(dcaz);
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_NSigma"))->Fill(nSigmaVtx);
		
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_ITSchi2PerCls"))->Fill(this->GetITSchi2PerCluster(iTrack));
		
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCncls"))->Fill(iTrack->GetTPCNcls());
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCclsOverF"))->Fill((Double_t)iTrack->GetTPCNcls()/(Double_t)iTrack->GetTPCNclsF());
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCxRows"))->Fill(iTrack->GetTPCCrossedRows());
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCchi2PerCls"))->Fill(iTrack->GetTPCchi2()/iTrack->GetTPCNcls());
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCchi2Global"))->Fill(iTrack->GetChi2TPCConstrainedVsGlobal(fESDevent->GetPrimaryVertex()));
		
		// -------------
		
		((TH2D*)fHistoMap->Get(particle + "_ITS_dEdx_P"))->Fill(pITS, dEdxITS);
		((TH2D*)fHistoMap->Get(particle + "_TPC_dEdx_P"))->Fill(pTPC, dEdxTPC);
		
		if(this->TOFmatch(iTrack))
		{
			beta = this->GetBeta(iTrack);
			m2   = this->GetMassSquare(iTrack);
			mass = TMath::Sqrt(TMath::Abs(m2));
			
			((TH2D*)fHistoMap->Get(particle + "_TOF_Beta_P"))->Fill(pTOF, beta);
			((TH2D*)fHistoMap->Get(particle + "_TOF_Mass_P"))->Fill(pTOF, mass);
		}
		
		// -----------------------------------------------
		//  for pid and efficiency
		// -----------------------------------------------
		
		Int_t simpid = -1;
		TParticle* iParticle = 0;
		
		if(fSimulation)
		{
			iParticle = this->GetParticle(iTrack);
			
			if(iParticle == 0) continue;
			
			simPt  = iParticle->Pt();
			simPhi = iParticle->Phi();
			simY   = iParticle->Y();
			
			simpid = fLnID->GetPID(iParticle);
			
			if(simpid == fPartCode)
			{
				TString simparticle = fSpecies;
				if(this->GetSign(iParticle)<0) simparticle.Prepend("Anti");
				
				((TH2D*)fHistoMap->Get(simparticle + "_Sim_PtY"))->Fill(simY, simPt);
				
				if(TMath::Abs(y) < fMaxY)
				{
					((TH2D*)fHistoMap->Get(simparticle + "_Response_Matrix"))->Fill(pt, simPt);
					
					if(this->IsPhysicalPrimary(iParticle))
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Ntrk"))->Fill(fNtrk);
					}
					
					// for pid
					if(!this->IsFakeTrack(iTrack))
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Pt"))->Fill(simPt);
						
						if(this->IsPhysicalPrimary(iParticle)) // the efficiency is calculated on the primaries
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Pt"))->Fill(simPt);
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Y"))->Fill(simY);
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Phi"))->Fill(simPhi);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_DCAxy_Pt"))->Fill(simPt,dcaxy);
						
							((TH2D*)fHistoMap->Get(simparticle + "_Prim_Response_Matrix"))->Fill(pt, simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Prim_DiffPt_RecPt"))->Fill(pt,simPt-pt);
							
							if( this->TOFmatch(iTrack) )
							{
								((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_M2_P"))->Fill(pTOF, m2);
								((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_M2_Pt"))->Fill(pt, m2);
							}
						}
						else if(this->IsFromWeakDecay(iParticle))
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fdwn_Pt"))->Fill(simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Fdwn_DCAxy_Pt"))->Fill(simPt,dcaxy);
						}
						else // if(this->IsFromMaterial(iParticle)
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Mat_Pt"))->Fill(simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Mat_DCAxy_Pt"))->Fill(simPt,dcaxy);
						}
					}
					else // fake tracks
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Pt"))->Fill(simPt);
						if(this->IsPhysicalPrimary(iParticle))
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Prim_Pt"))->Fill(simPt);
						}
						else if(this->IsFromWeakDecay(iParticle))
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Fdwn_Pt"))->Fill(simPt);
						}
						else
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Mat_Pt"))->Fill(simPt);
						}
					}
				}
			}
		}
		
		// get the pid
		Int_t pid = fLnID->GetPID( fPartCode, pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta, fMaxNSigmaITS, fMaxNSigmaTPC, fMaxNSigmaTOF);
		
		// pid table for prior probabilities (only Bayes)
		
		Int_t offset = AliPID::kDeuteron - 5;
		
		if( fSimulation )
		{
			Int_t sim = simpid > AliPID::kProton ? simpid - offset : simpid;
			Int_t rec = pid > AliPID::kProton ? pid - offset : pid;
			
			if((sim > -1) && (rec > -1)) // pid performance
			{
				((TH2D*)fHistoMap->Get(fSpecies + "_Stats_PID_Table"))->Fill(sim, rec);
				((TH2D*)fHistoMap->Get(fSpecies + "_Stats_PID_Table"))->Fill(9, rec);
			}
			if(sim > -1)
			{
				((TH2D*)fHistoMap->Get(fSpecies + "_Stats_PID_Table"))->Fill(sim, 9);
			}
		}
		
		// for bayes iteration
		if(pid != -1)
		{
			Int_t iBin = pid > AliPID::kProton ? pid - offset : pid;
			((TH1D*)fHistoMap->Get(fSpecies + "_Stats_PID"))->Fill(iBin);
		}
		
		if(pid != fPartCode) continue;
		
		Bool_t goodPid = 0;
		
		if(fSimulation)
		{
			goodPid = ( simpid == pid );
		}
		
		((TH1D*)fHistoMap->Get(particle + "_PID_Ntrk_pTPC"))->Fill(pTPC,fNtrk);
		((TH1D*)fHistoMap->Get(particle + "_PID_Zmult_pTPC"))->Fill(pTPC,fKNOmult);
		
		// pid performance
		((TH2D*)fHistoMap->Get(particle + "_PID_ITSdEdx_P"))->Fill(pITS, dEdxITS);
		((TH2D*)fHistoMap->Get(particle + "_PID_TPCdEdx_P"))->Fill(pTPC, dEdxTPC);
		
		if(this->TOFmatch(iTrack))
		{
			((TH2D*)fHistoMap->Get(particle + "_PID_Beta_P"))->Fill(pTOF, beta);
			((TH2D*)fHistoMap->Get(particle + "_PID_Mass_P"))->Fill(pTOF, mass);
			if(fSimulation && goodPid)
			{
				((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Mass"))->Fill(mass);
			}
		}
		
		// ---------------------------------
		//  results in |y| < fMaxY
		// ---------------------------------
		
		((TH1D*)fHistoMap->Get(particle + "_PID_Y"))->Fill(y);
		((TH2D*)fHistoMap->Get(particle + "_PID_Pt_Y"))->Fill(y, pt);
		
		if(TMath::Abs(y) >= fMaxY) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_PID_Pt"))->Fill(pt);
		((TH1D*)fHistoMap->Get(particle + "_PID_Phi"))->Fill(phi);
		
		if( this->TOFmatch(iTrack) )
		{
			((TH2D*)fHistoMap->Get(particle + "_PID_M2_Pt"))->Fill(pt, m2);
		}
		
		// secondaries
		
		Bool_t m2match = kTRUE;
		
		if( fTOFmatch && (fLnID->GetPidProcedure() > AliLnID::kMaxLikelihood))
		{
			if((m2 < fMinM2) || (m2 >= fMaxM2)) m2match = kFALSE;
		}
		
		if(m2match)
		{
			((TH2D*)fHistoMap->Get(particle + "_PID_DCAxy_Pt"))->Fill(pt, dcaxy);
			((TH2D*)fHistoMap->Get(particle + "_PID_DCAz_Pt"))->Fill(pt, dcaz);
			((TH2D*)fHistoMap->Get(particle + "_PID_NSigma_Pt"))->Fill(pt, nSigmaVtx);
		}
		
		if(fSimulation && goodPid)
		{
			// for unfolding and pid contamination
			if(!this->IsFakeTrack(iTrack))
			{
				((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Pt"))->Fill(simPt);
				
				if( this->TOFmatch(iTrack) )
				{
					((TH2D*)fHistoMap->Get(particle + "_Sim_PID_M2_Pt"))->Fill(pt, m2);
				}
				
				if(this->IsPhysicalPrimary(iParticle))
				{
					((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Prim_Pt"))->Fill(simPt);
				}
				
				if(m2match)
				{
					if(this->IsPhysicalPrimary(iParticle))
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_DCAxy_Pt"))->Fill(pt, dcaxy);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_DCAz_Pt"))->Fill(pt, dcaz);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_NSigma_Pt"))->Fill(pt, nSigmaVtx);
					}
					else if(this->IsFromWeakDecay(iParticle))
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fdwn_DCAxy_Pt"))->Fill(pt, dcaxy);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fdwn_DCAz_Pt"))->Fill(pt, dcaz);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fdwn_NSigma_Pt"))->Fill(pt, nSigmaVtx);
					}
					else // from materials
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Mat_DCAxy_Pt"))->Fill(pt, dcaxy);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Mat_DCAz_Pt"))->Fill(pt, dcaz);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Mat_NSigma_Pt"))->Fill(pt, nSigmaVtx);
					}
				}
			}
			else // fake tracks
			{
				((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Fake_Pt"))->Fill(simPt);
				
				if(m2match)
				{
					if(this->IsPhysicalPrimary(iParticle))
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fake_Prim_DCAxy_Pt"))->Fill(pt, dcaxy);
					}
					else if(this->IsFromWeakDecay(iParticle))
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fake_Fdwn_DCAxy_Pt"))->Fill(pt, dcaxy);
					}
					else // from materials
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fake_Mat_DCAxy_Pt"))->Fill(pt, dcaxy);
					}
				}
			}
		}
	}
	
	return nTracks;
}

void AliAnalysisTaskB2::Terminate(Option_t* )
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

Bool_t AliAnalysisTaskB2::IsV0AND() const
{
//
// signals in both V0A and V0C
//
	return ( fTrigAna->IsOfflineTriggerFired(fESDevent, AliTriggerAnalysis::kV0A) && 
	                 fTrigAna->IsOfflineTriggerFired(fESDevent, AliTriggerAnalysis::kV0C) );
}

Bool_t AliAnalysisTaskB2::IsFastOnly(UInt_t triggerBits) const
{
//
// kFastOnly trigger
//
	return ( (triggerBits&AliVEvent::kFastOnly) == AliVEvent::kFastOnly );
}

Bool_t AliAnalysisTaskB2::IsMB(UInt_t triggerBits) const
{
//
// MB event
//
	return ( (triggerBits&AliVEvent::kMB) == AliVEvent::kMB );
}

Bool_t AliAnalysisTaskB2::TOFmatch(const AliESDtrack* trk) const
{
//
// Check for TOF match signal
//
	return ( trk->IsOn(AliESDtrack::kTOFout) && trk->IsOn(AliESDtrack::kTIME) );
}

Bool_t AliAnalysisTaskB2::AcceptTOFtrack(const AliESDtrack* trk) const
{
//
// Additional checks for TOF match signal
//
	if( !this->TOFmatch(trk) ) return kFALSE;
	if( trk->GetIntegratedLength() < 350) return kFALSE;
	if( trk->GetTOFsignal() < 1e-6) return kFALSE;
	
	return kTRUE;
}

TParticle* AliAnalysisTaskB2::GetParticle(const AliESDtrack* trk) const
{
//
// Particle that left the track
//
	AliStack* stack = fMCevent->Stack();
	
	Int_t label = TMath::Abs(trk->GetLabel()); // if negative then it shares points from other tracks
	if( label >= fMCevent->GetNumberOfTracks() ) return 0;
	
	return stack->Particle(label);
}

Bool_t AliAnalysisTaskB2::IsFakeTrack(const AliESDtrack* trk) const
{
//
// Check if the track shares some clusters with different particles
// (definition changed to label=0? )
//
	return ( trk->GetLabel() < 0 );
}

Bool_t AliAnalysisTaskB2::IsPhysicalPrimary(const TParticle* prt) const
{
//
// Check if the particle is physical primary
//
	AliStack* stack = fMCevent->Stack();
	Int_t index = stack->Particles()->IndexOf(prt);
	
	return stack->IsPhysicalPrimary(index);
}

Bool_t AliAnalysisTaskB2::IsFromMaterial(const TParticle* prt) const
{
//
// Check if the particle is originated at the materials
//
	AliStack* stack = fMCevent->Stack();
	Int_t index = stack->Particles()->IndexOf(prt);
	
	return stack->IsSecondaryFromMaterial(index);
}

Bool_t AliAnalysisTaskB2::IsFromWeakDecay(const TParticle* prt) const
{
//
// Check if the particle comes from a weak decay
//
	AliStack* stack = fMCevent->Stack();
	Int_t index = stack->Particles()->IndexOf(prt);
	
	return stack->IsSecondaryFromWeakDecay(index);
}

Double_t AliAnalysisTaskB2::GetSign(TParticle* prt) const
{
//
// Sign of the particle
//
	TParticlePDG* pdg = prt->GetPDG();
	
	if(pdg != 0) return pdg->Charge();
	
	return 0;
}

Double_t AliAnalysisTaskB2::GetPhi(const AliESDtrack* trk) const
{
//
// Azimuthal angle [0,2pi) using the pt at the DCA point
//
	Double_t px = trk->Px();
	Double_t py = trk->Py();
	
	return TMath::Pi()+TMath::ATan2(-py, -px);
}

Double_t AliAnalysisTaskB2::GetTheta(const AliESDtrack* trk) const
{
//
// Polar angle using the pt at the DCA point
//
	Double_t p  = trk->GetP();
	Double_t pz = trk->Pz();
	
	return (pz == 0) ? TMath::PiOver2() : TMath::ACos(pz/p);
}

Double_t AliAnalysisTaskB2::GetRapidity(const AliESDtrack* trk, Int_t pid) const
{
//
// Rapidity assuming the given particle hypothesis
// and using the momentum at the DCA
//
	Double_t m  = AliPID::ParticleMass(pid);
	Double_t p  = (pid>AliPID::kTriton) ? 2.*trk->GetP() : trk->GetP();
	Double_t e  = TMath::Sqrt(p*p + m*m);
	Double_t pz = (pid>AliPID::kTriton) ? 2.*trk->Pz() : trk->Pz();
	
	if(e <= pz) return 1.e+16;
	
	return 0.5*TMath::Log( (e+pz)/(e-pz) );
}

Double_t AliAnalysisTaskB2::GetITSmomentum(const AliESDtrack* trk) const
{
//
// Momentum for ITS pid
//
	Double_t pDCA = trk->GetP();
	Double_t pTPC = trk->GetTPCmomentum();
	
	return (pDCA+pTPC)/2.;
}

Double_t AliAnalysisTaskB2::GetTOFmomentum(const AliESDtrack* trk) const
{
//
// Momentum for TOF pid
//
	Double_t pDCA = trk->GetP();
	
	const AliExternalTrackParam* param = trk->GetOuterParam();
	Double_t pOut = param ? param->GetP() : trk->GetP();
	
	return (pDCA+pOut)/2.;
}

Double_t AliAnalysisTaskB2::GetBeta(const AliESDtrack* trk) const
{
//
// Velocity
//
	Double_t t = this->GetTimeOfFlight(trk); // ps
	Double_t l = trk->GetIntegratedLength(); // cm
	
	if(t <= 0) return 1.e6; // 1M times the speed of light ;)
	
	return (l/t)/2.99792458e-2;
}

Double_t AliAnalysisTaskB2::GetMassSquare(const AliESDtrack* trk) const
{
//
// Square mass
//
	Double_t p = (fPartCode>AliPID::kTriton) ? 2.*this->GetTOFmomentum(trk) : this->GetTOFmomentum(trk);
	Double_t beta = this->GetBeta(trk);
	
	return p*p*(1./(beta*beta) - 1.);
}

Int_t AliAnalysisTaskB2::GetChargedMultiplicity(Double_t etaMax) const
{
//
// Charged particle multiplicity using ALICE physical primary definition
//
	AliStack* stack = fMCevent->Stack();
	
	Int_t nch = 0;
	//for (Int_t i=0; i < stack->GetNprimary(); ++i)
	for (Int_t i=0; i < stack->GetNtrack(); ++i)
	{
		if(!stack->IsPhysicalPrimary(i)) continue;
		
		TParticle* iParticle = stack->Particle(i);
		
		if(TMath::Abs(iParticle->Eta()) >= etaMax) continue;
		//if (iParticle->Pt() < -1) continue;
		
		TParticlePDG* iPDG = iParticle->GetPDG(); // There are some particles with no PDG
		if (iPDG && iPDG->Charge() == 0) continue;
		
		++nch;
	}
	
	return nch;
}

Double_t AliAnalysisTaskB2::GetTimeOfFlight(const AliESDtrack* trk) const
{
//
// Time of flight associated to the track.
// Adapted from ANALYSIS/AliAnalysisTaskESDfilter.cxx
//
	if(!fESDevent->GetTOFHeader())
	{ //protection in case the pass2 LHC10b,c,d have been processed without tender. 
		Float_t t0spread[10];
		Float_t intrinsicTOFres=100; //ps ok for LHC10b,c,d pass2!! 
		for (Int_t j=0; j<10; j++) t0spread[j] = (TMath::Sqrt(fESDevent->GetSigma2DiamondZ()))/0.03; //0.03 to convert from cm to ps
		fESDpid->GetTOFResponse().SetT0resolution(t0spread);
		fESDpid->GetTOFResponse().SetTimeResolution(intrinsicTOFres);
		
		fESDpid->SetTOFResponse(fESDevent, (AliESDpid::EStartTimeType_t)fTimeZeroType);
	}
	
	if(fESDevent->GetTOFHeader() && fIsPidOwner) fESDpid->SetTOFResponse(fESDevent, (AliESDpid::EStartTimeType_t)fTimeZeroType); //in case of AOD production strating form LHC10e without Tender. 
	
	Double_t timeZero = fESDpid->GetTOFResponse().GetStartTime(trk->P());
	
	return trk->GetTOFsignal()-timeZero;
}

Int_t AliAnalysisTaskB2::GetITSnClusters(const AliESDtrack* trk) const
{
//
// ITS number of clusters
//
	UChar_t map = trk->GetITSClusterMap();
	
	Int_t npoints=0;
	for(Int_t j=0; j<6; j++)
	{
		if(map&(1<<j)) ++npoints;
	}
	
	return npoints;
}

Double_t AliAnalysisTaskB2::GetITSchi2PerCluster(const AliESDtrack* trk) const
{
//
// ITS chi2 per number of clusters
//
	Double_t chi2 = trk->GetITSchi2();
	Int_t ncls = this->GetITSnClusters(trk);
	
	Double_t chi2ncls = (ncls==0) ? 1.e10 : chi2/ncls;
	
	return chi2ncls;
}

Int_t AliAnalysisTaskB2::GetITSnPointsPID(const AliESDtrack* trk) const
{
//
// ITS number of points for PID
//
	UChar_t map = trk->GetITSClusterMap();
	
	Int_t npoints = 0;
	for(Int_t j=2; j<6; j++)
	{
		if(map&(1<<j)) ++npoints;
	}
	
	return npoints;
}

Int_t AliAnalysisTaskB2::GetPidCode(const TString& species) const
{
//
// Return AliPID code of the given species
//
	TString name = species;
	name.ToLower();
	
	if(name == "electron") return AliPID::kElectron;
	if(name == "muon")     return AliPID::kMuon;
	if(name == "pion")     return AliPID::kPion;
	if(name == "kaon")     return AliPID::kKaon;
	if(name == "proton")   return AliPID::kProton;
	if(name == "deuteron") return AliPID::kDeuteron;
	if(name == "triton")   return AliPID::kTriton;
	if(name == "he3")      return AliPID::kHe3;
	if(name == "alpha")    return AliPID::kAlpha;
	
	return -1;
}
