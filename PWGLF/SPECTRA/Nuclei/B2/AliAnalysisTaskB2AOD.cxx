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

// Analysis task for B2 (AOD)
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <AliAnalysisTaskSE.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliAODVertex.h>
#include <AliAODTrack.h>
#include <AliAODMCParticle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TProfile.h>
#include <AliLog.h>

#include "AliLnID.h"
#include "AliLnHistoMap.h"
#include "AliLnAODtrackCuts.h"
#include "AliAnalysisTaskB2AOD.h"

ClassImp(AliAnalysisTaskB2AOD)

AliAnalysisTaskB2AOD::AliAnalysisTaskB2AOD()
: AliAnalysisTaskSE()
, fSpecies("Proton")
, fPartCode(AliPID::kProton)
, fHeavyIons(0)
, fSimulation(0)
, fMultTriggerFired(0)
, fCentTriggerFired(0)
, fTriggerFired(0)
, fGoodVertex(0)
, fPileUpEvent(0)
, fV0AND(0)
, fNoFastOnly(0)
, fNtrkMultTrigger(0)
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
, fAODevent(0)
, fOutputContainer(0)
, fHistoMap(0)
, fTrackCuts(0)
, fLnID(0)
, fMaxNSigmaITS(3)
, fMaxNSigmaTPC(3)
, fMaxNSigmaTOF(3)
, fMinM2(2.)
, fMaxM2(6.)
, fMomentumCorrection(0)
, fMoCpfx(0)

{
//
// Default constructor
//
	AliLog::SetGlobalLogLevel(AliLog::kFatal);
}

AliAnalysisTaskB2AOD::AliAnalysisTaskB2AOD(const char* name)
: AliAnalysisTaskSE(name)
, fSpecies("Proton")
, fPartCode(AliPID::kProton)
, fHeavyIons(0)
, fSimulation(0)
, fMultTriggerFired(0)
, fCentTriggerFired(0)
, fTriggerFired(0)
, fGoodVertex(0)
, fPileUpEvent(0)
, fV0AND(0)
, fNoFastOnly(0)
, fNtrkMultTrigger(0)
, fMinKNOmult(-10)
, fMaxKNOmult(100000)
, fMinCentrality(-1)
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
, fAODevent(0)
, fOutputContainer(0)
, fHistoMap(0)
, fTrackCuts(0)
, fLnID(0)
, fMaxNSigmaITS(3)
, fMaxNSigmaTPC(3)
, fMaxNSigmaTOF(3)
, fMinM2(2.)
, fMaxM2(6.)
, fMomentumCorrection(0)
, fMoCpfx(0)

{
//
// Constructor
//
	DefineOutput(1, TList::Class());
	
	AliLog::SetGlobalLogLevel(AliLog::kFatal);
}

void AliAnalysisTaskB2AOD::UserCreateOutputObjects()
{
//
// Create output objects
//
	if(fHistoMap == 0) AliFatal("no histogram map"); // should be created somewhere else
	
	fOutputContainer = new TList();
	fOutputContainer->SetOwner(kTRUE);
	
	TObjString* key;
	TIter iter(fHistoMap->GetMap());
	while((key = dynamic_cast<TObjString*>(iter.Next()))) fOutputContainer->Add((TH1*)fHistoMap->Get(key));
	
	PostData(1, fOutputContainer);
}

AliAnalysisTaskB2AOD::~AliAnalysisTaskB2AOD()
{
//
// Default destructor
//
	delete fHistoMap;
	delete fOutputContainer;
	delete fTrackCuts;
	delete fLnID;
	
	delete fMoCpfx;
}

void AliAnalysisTaskB2AOD::SetParticleSpecies(const TString& species)
{
//
// set the particle species and the AliPID code
//
	fSpecies = species;
	fPartCode = this->GetPidCode(species);
}

void AliAnalysisTaskB2AOD::UserExec(Option_t* )
{
//
// Execute analysis for the current event
//
	fAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
	
	if (!fAODevent)
	{
		this->Error("Exec", "AOD event not available");
		return;
	}
	
	if(fTrackCuts == 0) AliFatal("track cuts not set");
	
	if(fLnID == 0) AliFatal("PID not set");
	
	// multiplicity and centrality
	
	fNtrk = (fMaxEta > 0.5) ? fAODevent->GetHeader()->GetRefMultiplicityComb08() : fAODevent->GetHeader()->GetRefMultiplicityComb05();
	
	if(fSimulation) fNch = this->GetChargedMultiplicity(fMaxEta);
	
	fKNOmult = fNtrk/fMeanNtrk;
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Event_Ntrk"))->Fill(fNtrk);
	((TH1D*)fHistoMap->Get(fSpecies + "_Event_Zmult"))->Fill(fKNOmult);
	
	fMultTriggerFired = (fNtrk > 0) && (fKNOmult >= fMinKNOmult) && (fKNOmult < fMaxKNOmult);
	
	if(fHeavyIons)
	{
		Double_t centrality = fAODevent->GetHeader()->GetCentrality();

		fCentTriggerFired = (centrality >= fMinCentrality) && (centrality < fMaxCentrality);
	}
	
	// trigger
	
	AliInputEventHandler* eventH = dynamic_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(eventH == 0)
	{
		this->Error("Exec", "could not get AliInputEventHandler");
		return;
	}
	
	UInt_t triggerBits = eventH->IsEventSelected();
	
	if(fHeavyIons)
	{
		fTriggerFired = ( this->IsMB(triggerBits) && fCentTriggerFired );
	}
	else
	{
		fTriggerFired = this->IsMB(triggerBits);
		if(fNoFastOnly)      fTriggerFired = ( fTriggerFired && !this->IsFastOnly(triggerBits) );
		if(fV0AND)           fTriggerFired = ( fTriggerFired && this->IsV0AND(triggerBits) );
		if(fNtrkMultTrigger) fTriggerFired = ( fTriggerFired && fMultTriggerFired );
	}
	
	// vertex
	
	fGoodVertex  = kFALSE;
	
	const AliAODVertex* vtx = fAODevent->GetPrimaryVertex(); // best primary vertex
	
	if(vtx->GetNContributors()>0)
	{
		fGoodVertex=kTRUE;
		
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_YX"))->Fill(vtx->GetX(), vtx->GetY());
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_YZ"))->Fill(vtx->GetZ(), vtx->GetY());
		((TH2D*)fHistoMap->Get(fSpecies + "_Vertex_XZ"))->Fill(vtx->GetZ(), vtx->GetX());
	}
	
	if( (vtx->GetX() <= fMinVx) || (vtx->GetX() >= fMaxVx) ) fGoodVertex=kFALSE;
	if( (vtx->GetY() <= fMinVy) || (vtx->GetY() >= fMaxVy) ) fGoodVertex=kFALSE;
	if( (vtx->GetZ() <= fMinVz) || (vtx->GetZ() >= fMaxVz) ) fGoodVertex=kFALSE;
	
	// pile up
	
	fPileUpEvent = fAODevent->IsPileupFromSPDInMultBins();
	
	// event stats
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(0); // number of events
	
	if(fTriggerFired)
	{
		((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(1); // triggering events
		
		if(vtx->GetNContributors()>0)
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
	
	// Particles and tracks for this event
	
	if(fSimulation) this->GetParticles();
	
	this->GetTracks();
	
	// Post the data (input/output slots #0 already used by the base class)
	PostData(1, fOutputContainer);
}

Int_t AliAnalysisTaskB2AOD::GetParticles()
{
//
// Get particles from current event
//
	Int_t nParticles = 0;
	
	TList* lst = fAODevent->GetList();
	TClonesArray* lop = dynamic_cast<TClonesArray*>(lst->FindObject(AliAODMCParticle::StdBranchName()));
	
	if (!lop)
	{
		AliDebug(AliLog::kWarning, "list of particles not available");
		return 0;
	}
	
	for (Int_t i = 0; i < lop->GetEntries(); ++i)
	{
		AliAODMCParticle* iParticle = dynamic_cast<AliAODMCParticle*>(lop->At(i));
		
		if(!iParticle) continue;
		
		Int_t pid = fLnID->GetPID(iParticle);
		
		if(pid != fPartCode) continue;
		
		// physical primaries
		
		if(!iParticle->IsPhysicalPrimary()) continue;
		
		TString particle = fSpecies;
		if(iParticle->Charge() < 0) particle.Prepend("Anti");
		
		Double_t genP   = iParticle->P();
		Double_t genPt  = iParticle->Pt();
		Double_t genY   = iParticle->Y();
		Double_t genPhi = iParticle->Phi();
		Double_t genEta = iParticle->Eta();
		
		// multiplicity and centrality
		
		if( fHeavyIons && !fCentTriggerFired) continue;
		if( fNtrkMultTrigger && !fMultTriggerFired) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Pt"))->Fill(genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Y"))->Fill(genY);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Phi"))->Fill(genPhi);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Prim_Eta"))->Fill(genEta);
		((TH2D*)fHistoMap->Get(particle + "_Gen_Prim_PtY"))->Fill(genY, genPt);
		((TH2D*)fHistoMap->Get(particle + "_Gen_Prim_EtaY"))->Fill(genY, genEta);
		
		// is within phase space?
		
		if(TMath::Abs(genY) >= fMaxY) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_PhS_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_PhS_Prim_Pt"))->Fill(genPt);
		
		// is from a triggering event? (same as rec.)
		
		if(!fTriggerFired)
		{
			((TH1D*)fHistoMap->Get(particle + "_Gen_NoTrig_Prim_Pt"))->Fill(genPt);
			continue;
		}
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Trig_Prim_Pt"))->Fill(genPt);
		
		// is from a triggering event with good vertex?
		
		if(fAODevent->GetPrimaryVertex()->GetNContributors() <= 0)
		{
			((TH1D*)fHistoMap->Get(particle + "_Gen_NoVtx_Prim_Pt"))->Fill(genPt);
		}
		
		if(!fGoodVertex) continue;
		if(fPileUpEvent) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_Gen_Nch"))->Fill(fNch);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Vtx_Prim_Pt"))->Fill(genPt);
		
		// is within the geometrical acceptance?
		
		if(!fTrackCuts->IsWithinGeoAcceptance(iParticle)) continue;
		
		((TH2D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_PtY"))->Fill(genY,genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_P"))->Fill(genP);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Pt"))->Fill(genPt);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Phi"))->Fill(genPhi);
		((TH1D*)fHistoMap->Get(particle + "_Gen_Acc_Prim_Y"))->Fill(genY);
	}
	
	return nParticles;
}

Int_t AliAnalysisTaskB2AOD::GetTracks()
{
//
// Get tracks from current event
//
	using namespace std;
	
	Int_t nTracks = 0;
	
	// trigger, vertex and pile-up
	
	if(!fTriggerFired) return 0;
	if(!fGoodVertex)   return 0;
	if(fPileUpEvent)   return 0;
	
	// this is a 'good' event
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Stats"))->Fill(2); // analyzed events
	
	((TH1D*)fHistoMap->Get(fSpecies + "_Ana_Event_Ntrk"))->Fill(fNtrk);
	((TH1D*)fHistoMap->Get(fSpecies + "_Ana_Event_Zmult"))->Fill(fKNOmult);
	
	if(fSimulation)
	{
		((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Event_Nch_Ntrk"))->Fill(fNtrk, fNch);
	}
	
	const AliAODVertex* vtx = fAODevent->GetPrimaryVertex();
	
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_YX"))->Fill(vtx->GetX(), vtx->GetY());
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_YZ"))->Fill(vtx->GetZ(), vtx->GetY());
	((TH2D*)fHistoMap->Get(fSpecies + "_Ana_Vertex_XZ"))->Fill(vtx->GetZ(), vtx->GetX());
	
	// track loop
	for(Int_t i = 0; i < fAODevent->GetNumberOfTracks(); ++i)
	{
		AliAODTrack* iTrack = fAODevent->GetTrack(i);
		
		if(!iTrack) continue;
		
		TString particle = fSpecies;
		if(iTrack->Charge() < 0) particle.Prepend("Anti");
		
		// track parameters at the DCA to the primary vertex
		
		AliExternalTrackParam trkparam;
		trkparam.CopyFromVTrack(iTrack);
		
		if(trkparam.GetX() > 3.) continue; // only valid for propagation inside the beam pipe
		
		// impact parameters
		Double_t b[2];
		Double_t bCov[3];
		if(!trkparam.PropagateToDCA(fAODevent->GetPrimaryVertex(), fAODevent->GetMagneticField(), 10000., b, bCov)) continue;
		
		Double_t dcaxy = b[0];
		Double_t dcaz  = b[1];
		
		if(iTrack->Charge() < 0) // in case of asymmetry
		{
			dcaxy = -dcaxy;
			dcaz  = -dcaz;
		}
		
		Double_t nSigmaVtx = fTrackCuts->GetNSigmaToVertex(b, bCov);
		
		// momentum
		Double_t pDCA[3] = { iTrack->Px(), iTrack->Py(), iTrack->Pz() };
		//trkparam.GetPxPyPz(pDCA);
		
		Double_t p  = TMath::Sqrt(pDCA[0]*pDCA[0] + pDCA[1]*pDCA[1] + pDCA[2]*pDCA[2]);
		Double_t pt = TMath::Sqrt(pDCA[0]*pDCA[0] + pDCA[1]*pDCA[1]);
		Double_t pz = pDCA[2];
		
		// track cuts
		
		TVector3 pVect(pDCA[0],pDCA[1],pDCA[2]);
		Double_t eta = pVect.Eta();
		Double_t phi = this->GetPhi(pDCA);
		
		((TH2D*)fHistoMap->Get(particle + "_Before_Phi_Eta"))->Fill(eta,phi);
		
		if(!fTrackCuts->AcceptTrack(iTrack,b, bCov)) continue;  // with next track
		if(!fTrackCuts->IsWithinGeoAcceptance(pDCA)) continue; // with next track
		
		// end track cuts
		
		++nTracks;
		
		((TH2D*)fHistoMap->Get(particle + "_After_Phi_Eta"))->Fill(eta, phi);
		
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_DCAxy"))->Fill(dcaxy);
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_DCAz"))->Fill(dcaz);
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_NSigma"))->Fill(nSigmaVtx);
		//((TH1D*)fHistoMap->Get(particle + "_TrackCuts_ITSchi2PerCls"))->Fill(this->GetITSchi2PerCluster(iTrack));
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCncls"))->Fill(iTrack->GetTPCNcls());
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCxRowsOverF"))->Fill(fTrackCuts->GetNTPCXRowsOverFindable(iTrack));
		((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCxRows"))->Fill(iTrack->GetTPCNCrossedRows());
		//((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCchi2PerCls"))->Fill(iTrack->GetTPCchi2()/iTrack->GetTPCNcls());
		//((TH1D*)fHistoMap->Get(particle + "_TrackCuts_TPCchi2Global"))->Fill(iTrack->GetChi2TPCConstrainedVsGlobal(fAODevent->GetPrimaryVertex()));
		
		// detector signals
		
		Double_t pITS    = p;
		Double_t pTPC    = iTrack->GetTPCmomentum();
		Double_t pTOF    = p;
		Double_t dEdxITS = iTrack->GetITSsignal();
		Double_t dEdxTPC = iTrack->GetTPCsignal();
		Int_t nPointsITS = this->GetITSnPointsPID(iTrack);
		Int_t nPointsTPC = iTrack->GetTPCsignalN();
		Double_t beta    = 0;
		Double_t mass    = 0;
		Double_t m2      = 0;
		
		((TH2D*)fHistoMap->Get(particle + "_ITS_dEdx_P"))->Fill(pITS, dEdxITS);
		((TH2D*)fHistoMap->Get(particle + "_TPC_dEdx_P"))->Fill(pTPC, dEdxTPC);
		
		if(fTrackCuts->TOFmatch())
		{
			beta = this->GetBeta(iTrack);
			m2   = this->GetMassSquared(p, beta);
			mass = TMath::Sqrt(TMath::Abs(m2));
			
			((TH2D*)fHistoMap->Get(particle + "_TOF_Beta_P"))->Fill(pTOF, beta);
			((TH2D*)fHistoMap->Get(particle + "_TOF_Mass_P"))->Fill(pTOF, mass);
		}
		
		// get the pid
		Int_t pid = fLnID->GetPID( fPartCode, pITS, dEdxITS, nPointsITS, pTPC, dEdxTPC, nPointsTPC, pTOF, beta, fMaxNSigmaITS, fMaxNSigmaTPC, fMaxNSigmaTOF);
		
		Int_t offset = AliPID::kDeuteron - 5;
		
		// for bayes iteration
		if(pid != -1)
		{
			Int_t iBin = (pid > AliPID::kProton) ? pid - offset : pid;
			((TH1D*)fHistoMap->Get(fSpecies + "_Stats_PID"))->Fill(iBin);
		}
		
		// fix momentum
		
		if(fPartCode > AliPID::kTriton)
		{
			p  *= 2.;
			pt *= 2.;
			pz *= 2.;
		}
		
		if(fMomentumCorrection)
		{
			pt += this->GetMomentumCorrection(pt);
			p   = TMath::Sqrt(pt*pt + pz*pz);
			
			if(fTrackCuts->TOFmatch())
			{
				m2   = this->GetMassSquared(p, beta);
				mass = TMath::Sqrt(TMath::Abs(m2));
			}
		}
		
		// for pid and efficiency
		
		Double_t simPt  = 0;
		Double_t simPhi = 0;
		Double_t simY   = 0;
		
		Int_t simpid = -1;
		AliAODMCParticle* iParticle = 0;
		
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
				if(iParticle->Charge()<0) simparticle.Prepend("Anti");
				
				((TH2D*)fHistoMap->Get(simparticle + "_Sim_PtY"))->Fill(simY, simPt);
				
				if(TMath::Abs(simY) < fMaxY)
				{
					((TH2D*)fHistoMap->Get(simparticle + "_Response_Matrix"))->Fill(pt, simPt);
					
					if(iParticle->IsPhysicalPrimary())
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Ntrk"))->Fill(fNtrk);
					}
					
					// for pid
					if(!this->IsFakeTrack(iTrack))
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Pt"))->Fill(simPt);
						
						if(iParticle->IsPhysicalPrimary()) // the efficiency is calculated on the primaries
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Pt"))->Fill(simPt);
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Y"))->Fill(simY);
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Phi"))->Fill(simPhi);
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Prim_Rec_Pt"))->Fill(pt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_DCAxy_Pt"))->Fill(simPt,dcaxy);
						
							((TH2D*)fHistoMap->Get(simparticle + "_Prim_Response_Matrix"))->Fill(pt, simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Prim_DiffPt_RecPt"))->Fill(pt,simPt-pt);
							
							if( fTrackCuts->TOFmatch() )
							{
								((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_M2_P"))->Fill(pTOF, m2);
								((TH2D*)fHistoMap->Get(simparticle + "_Sim_Prim_M2_Pt"))->Fill(pt, m2);
							}
						}
						else if(iParticle->IsSecondaryFromWeakDecay())
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fdwn_Pt"))->Fill(simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Fdwn_DCAxy_Pt"))->Fill(simPt,dcaxy);
						}
						else
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Mat_Pt"))->Fill(simPt);
							
							((TH2D*)fHistoMap->Get(simparticle + "_Sim_Mat_DCAxy_Pt"))->Fill(simPt,dcaxy);
						}
					}
					else // fake tracks
					{
						((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Pt"))->Fill(simPt);
						if(iParticle->IsPhysicalPrimary())
						{
							((TH1D*)fHistoMap->Get(simparticle + "_Sim_Fake_Prim_Pt"))->Fill(simPt);
						}
						else if(iParticle->IsSecondaryFromWeakDecay())
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
			
			// pid table for prior probabilities (only Bayes)
			
			Int_t sim = (simpid > AliPID::kProton) ? simpid - offset : simpid;
			Int_t rec = (pid > AliPID::kProton) ? pid - offset : pid;
			
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
		
		// candidate tracks
		
		if(pid != fPartCode) continue;
		
		Bool_t goodPid = 0;
		
		if(fSimulation)
		{
			goodPid = ( simpid == pid );
		}
		
		((TH2D*)fHistoMap->Get(particle + "_PID_Ntrk_pTPC"))->Fill(pTPC,fNtrk);
		((TH2D*)fHistoMap->Get(particle + "_PID_Zmult_pTPC"))->Fill(pTPC,fKNOmult);
		
		// pid performance
		((TH2D*)fHistoMap->Get(particle + "_PID_ITSdEdx_P"))->Fill(pITS, dEdxITS);
		((TH2D*)fHistoMap->Get(particle + "_PID_TPCdEdx_P"))->Fill(pTPC, dEdxTPC);
		
		if(fTrackCuts->TOFmatch())
		{
			((TH2D*)fHistoMap->Get(particle + "_PID_Beta_P"))->Fill(pTOF, beta);
			((TH2D*)fHistoMap->Get(particle + "_PID_Mass_P"))->Fill(pTOF, mass);
			if(fSimulation && goodPid)
			{
				((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Mass"))->Fill(mass);
			}
		}
		
		Double_t y = this->GetRapidity(p, pz, AliPID::ParticleMass(fPartCode));
		
		((TH1D*)fHistoMap->Get(particle + "_PID_Y"))->Fill(y);
		((TH2D*)fHistoMap->Get(particle + "_PID_Pt_Y"))->Fill(y, pt);
		
		//  results in |y| < fMaxY
		
		if(TMath::Abs(y) >= fMaxY) continue;
		
		((TH1D*)fHistoMap->Get(particle + "_PID_Pt"))->Fill(pt);
		((TH1D*)fHistoMap->Get(particle + "_PID_Phi"))->Fill(phi);
		
		if(iTrack->IsOn(AliAODTrack::kTRDin))
		{
			((TH1D*)fHistoMap->Get(particle + "_PID_TRDin_Pt"))->Fill(pt);
			
			if(iTrack->IsOn(AliAODTrack::kTRDout)) ((TH1D*)fHistoMap->Get(particle + "_PID_TRDin_TRDout_Pt"))->Fill(pt);
			if(iTrack->IsOn(AliAODTrack::kTOFout)) ((TH1D*)fHistoMap->Get(particle + "_PID_TRDin_TOFout_Pt"))->Fill(pt);
		}
		
		if(iTrack->IsOn(AliAODTrack::kTOFin))
		{
			((TH1D*)fHistoMap->Get(particle + "_PID_TOFin_Pt"))->Fill(pt);
			
			if(iTrack->IsOn(AliAODTrack::kTOFout)) ((TH1D*)fHistoMap->Get(particle + "_PID_TOFin_TOFout_Pt"))->Fill(pt);
		}
		
		if(fTrackCuts->TOFmatch())
		{
			Double_t dm2 = this->GetM2Difference(beta, p, AliPID::ParticleMass(fPartCode));
			Double_t t   = iTrack->GetTOFsignal()*1.e-3; // ns
			Double_t dt  = t - this->GetExpectedTime(iTrack, AliPID::ParticleMass(fPartCode))*1.e-3;
			
			((TH2D*)fHistoMap->Get(particle + "_PID_M2_Pt"))->Fill(pt, m2);
			((TH2D*)fHistoMap->Get(particle + "_PID_DM2_Pt"))->Fill(pt, dm2);
			((TH2D*)fHistoMap->Get(particle + "_PID_Time_Pt"))->Fill(pt, t);
			((TH2D*)fHistoMap->Get(particle + "_PID_DTime_Pt"))->Fill(pt, dt);
			((TH1D*)fHistoMap->Get(particle + "_PID_TOFmatch_Pt"))->Fill(pt);
		}
		
		// secondaries
		
		Bool_t m2match = kTRUE;
		
		if( fTrackCuts->TOFmatch() && (fLnID->GetPidProcedure() > AliLnID::kMaxLikelihood))
		{
			if ((m2 < fMinM2) || (m2 >= fMaxM2)) m2match = kFALSE;
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
				
				if( fTrackCuts->TOFmatch() )
				{
					((TH2D*)fHistoMap->Get(particle + "_Sim_PID_M2_Pt"))->Fill(pt, m2);
				}
				
				if(iParticle->IsPhysicalPrimary())
				{
					((TH1D*)fHistoMap->Get(particle + "_Sim_PID_Prim_Pt"))->Fill(simPt);
				}
				
				if(m2match)
				{
					if(iParticle->IsPhysicalPrimary())
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_DCAxy_Pt"))->Fill(pt, dcaxy);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_DCAz_Pt"))->Fill(pt, dcaz);
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Prim_NSigma_Pt"))->Fill(pt, nSigmaVtx);
					}
					else if(iParticle->IsSecondaryFromWeakDecay())
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
					if(iParticle->IsPhysicalPrimary())
					{
						((TH2D*)fHistoMap->Get(particle + "_Sim_PID_Fake_Prim_DCAxy_Pt"))->Fill(pt, dcaxy);
					}
					else if(iParticle->IsSecondaryFromWeakDecay())
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

void AliAnalysisTaskB2AOD::Terminate(Option_t* )
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

}

Bool_t AliAnalysisTaskB2AOD::IsV0AND(UInt_t /*triggerBits*/) const
{
//
// signals in both V0A and V0C
//
	//return ( (triggerBits&AliVEvent::kINT7) == AliVEvent::kINT7 );
	
	AliAODVZERO* vzero = fAODevent->GetVZEROData();
	if(vzero == 0) return kFALSE;
	
	if(!vzero->TestBit(AliAODVZERO::kDecisionFilled)) return kFALSE;
	
	return ( (vzero->GetV0ADecision() == AliVVZERO::kV0BB) && (vzero->GetV0CDecision() == AliVVZERO::kV0BB) );
}

Bool_t AliAnalysisTaskB2AOD::IsFastOnly(UInt_t triggerBits) const
{
//
// kFastOnly trigger
//
	return ( (triggerBits&AliVEvent::kFastOnly) == AliVEvent::kFastOnly );
}

Bool_t AliAnalysisTaskB2AOD::IsMB(UInt_t triggerBits) const
{
//
// MB event
//
	return ( (triggerBits&AliVEvent::kMB) == AliVEvent::kMB );
}

AliAODMCParticle* AliAnalysisTaskB2AOD::GetParticle(const AliAODTrack* trk) const
{
//
// Particle that left the track
//
	TList* lst = fAODevent->GetList();
	TClonesArray* lop = dynamic_cast<TClonesArray*>(lst->FindObject(AliAODMCParticle::StdBranchName()));
	
	if (!lop) return 0;
	
	Int_t label = TMath::Abs(trk->GetLabel()); // if negative then it shares points from other tracks
	if( label >= lop->GetEntries() ) return 0;
	
	return dynamic_cast<AliAODMCParticle*>(lop->At(label));
}

Bool_t AliAnalysisTaskB2AOD::IsFakeTrack(const AliAODTrack* trk) const
{
//
// Check if the track shares some clusters with different particles
// (definition changed to label=0? )
//
	return ( trk->GetLabel() < 0 );
}

Double_t AliAnalysisTaskB2AOD::GetPhi(Double_t p[3]) const
{
//
// Azimuthal angle [0,2pi)
//
	Double_t px = p[0];
	Double_t py = p[1];
	
	return TMath::Pi()+TMath::ATan2(-py, -px);
}

Double_t AliAnalysisTaskB2AOD::GetBeta(const AliAODTrack* trk) const
{
//
// Velocity
//
	Double_t t = trk->GetTOFsignal(); // ps
	Double_t l = fTrackCuts->GetIntegratedLength(trk, fPartCode); // cm
	
	if(t <= 0) return 1.e+6; // 1M times the speed of light ;)
	
	return (l/t)/2.99792458e-2;
}

Double_t AliAnalysisTaskB2AOD::GetExpectedTime(const AliAODTrack* trk, Double_t m) const
{
//
// Expected time workaround for the given mass hypothesis (ps)
//
	Double_t p = (fPartCode > AliPID::kTriton) ? 2.*trk->P() : trk->P();
	Double_t beta = p/TMath::Sqrt(p*p + m*m);
	Double_t l = fTrackCuts->GetIntegratedLength(trk, fPartCode);
	
	return l/beta/2.99792458e-2;
}

Double_t AliAnalysisTaskB2AOD::GetMassSquared(Double_t p, Double_t beta) const
{
//
// Mass squared
//
	return p*p*(1./(beta*beta) - 1.);
}

Double_t AliAnalysisTaskB2AOD::GetM2Difference(Double_t beta, Double_t p, Double_t m) const
{
//
// Mass squared difference
//
	Double_t expBeta2 = p*p/(p*p+m*m);
	
	return p*p*(1./(beta*beta)-1./expBeta2);
}

Int_t AliAnalysisTaskB2AOD::GetChargedMultiplicity(Double_t etaMax) const
{
//
// Charged particle multiplicity using ALICE physical primary definition
//
	TList* lst = fAODevent->GetList();
	TClonesArray* lop = dynamic_cast<TClonesArray*>(lst->FindObject(AliAODMCParticle::StdBranchName()));
	
	if (!lop) return 0;
	
	Int_t nch = 0;
	
	for (Int_t i = 0; i < lop->GetEntries(); ++i)
	{
		AliAODMCParticle* iParticle = dynamic_cast<AliAODMCParticle*>(lop->At(i));

		if(!iParticle) continue;
		
		if(!iParticle->IsPhysicalPrimary()) continue;
		
		if(TMath::Abs(iParticle->Eta()) >= etaMax) continue;
		
		if(iParticle->Charge()==0) continue;
		
		++nch;
	}
	
	return nch;
}

Int_t AliAnalysisTaskB2AOD::GetITSnClusters(const AliAODTrack* trk) const
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

Int_t AliAnalysisTaskB2AOD::GetITSnPointsPID(const AliAODTrack* trk) const
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

Int_t AliAnalysisTaskB2AOD::GetPidCode(const TString& species) const
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

Double_t AliAnalysisTaskB2AOD::GetMomentumCorrection(Double_t ptrec) const
{
//
// momentum correction for low pt
//
	if(fMoCpfx == 0) return 0;
	
	return fMoCpfx->Interpolate(ptrec);
}

Double_t AliAnalysisTaskB2AOD::GetRapidity(Double_t p, Double_t pz, Double_t m) const
{
//
// Rapidity
//
	Double_t e  = TMath::Sqrt(p*p + m*m);
	
	if(e <= pz) return 1.e+16;
	
	return 0.5*TMath::Log( (e+pz)/(e-pz) );
}
