/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskCorrHadrons.h 45956 2010-12-10 12:55:37Z agheata $ */
/* AliAnalysisTaskCorrHadronsPP.h

 * This task is written to analyse the properties of jet-particles by two-particle correlation method
 * Author : Minwoo Kim ( minwoo.kim@cern.ch )

 *
 */
#ifndef ALIANALYSISTASKEPCORRPP_H
#define ALIANALYSISTASKEPCORRPP_H

class TH1F;
class TList;
class AliAnalysisManager;
class AliESDtrackCuts;

#include "TObjArray.h"

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "AliESDEvent.h"
#include "AliInputEventHandler.h" // event mixing
#include "AliEventPoolManager.h"  // event mixing
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliVVertex.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"

#include <vector>
using std::vector;

class AliAnalysisTaskEPCorrPP : public AliAnalysisTaskSE {
	public:
		enum{
            kCentBin = 1,
            kZvertBin = 7,
            kpTBin = 7,
            kPID = 1
        };

		enum{
			kTPCOnlyTrackCut = 128,
			kHybridTrackCut = 768
		};

		AliAnalysisTaskEPCorrPP();
		AliAnalysisTaskEPCorrPP(const char *name);
		virtual ~AliAnalysisTaskEPCorrPP();

		virtual void     UserCreateOutputObjects();
		virtual void     UserExec(Option_t *option);
		virtual void     Terminate(Option_t *);

		
	private:
		TList           *fOutput;        // Output list
//		AliESDtrackCuts *fTrackCuts;     // Track cuts
		AliEventPoolManager*  fPoolMgr; //!
		AliPIDResponse *fPIDResponse; //!

		TObjArray*        fMyprimRecoTracks; //
		TObjArray*        fTracksMixing; //

		TObjArray*    AcceptTracksReduced(AliAODEvent *event, Bool_t useCuts);
        float CalculatedPhiStar(float dPhi, float dEta, float Zv, float Zv2, float pT, float pT2, float bSign);
		void SetupForMixing();
		void FillMixedHistos(TObjArray* particle, TObjArray* particleMixed, int cBin, int zBin, float bSign, double weight);
		void DoMixing(double cent, double zvtx, TObjArray* fMyprimRecoTracks, float bSign);
//	    double DeltaPhi(double phi1, double phi2);
		int GetParticleID(AliAODTrack *track, double nsigmaCut);

		Int_t         fMinNumTrack; // AliEventPoolManager(), Size of track buffer for event mixing (number of tracks to fill the pool)
		Int_t         fPoolSize; // AliEventPoolManager(), max number of event to mix
		Int_t         fMinNEventsToMix; //

        double      fNsigmaCut; // nsigmaCut with TOF and TPC

		TH1F            *fHistPt[kPID];        //! Pt spectrum
		TH1F            *fHistEta[kPID];       //! pseudorapidity spectrum
		TH1F            *fHistPhi[kPID];       //! Azimuthal angle
		TH1F            *fHistZvertex;       //! primary z-vertex of each events
		TH1F            *fHistZvertexBin[kZvertBin];       //! primary z-vertex of each events
		TH1F 			*fHistCent; //! centrality by V0M
		TH1F			*fHistCentBin[kCentBin]; //! centrality by V0M

		// Nevt histo for normalisation
		TH1F 			*fHistNevtSame[kCentBin][kZvertBin]; //! Nevt for each categories(cBin, zBin) in same events (! removed 140505)
		TH1F 			*fHistNevtMixed[kCentBin][kZvertBin]; //! Nevt for each categories(cBin, zBin) in mixed events

		// array for naming
//		char hname[10000]; char htit[10000];

		// Histograms for pT,trig 
		TH1F			*fHistPtSame[kCentBin][kZvertBin]; //! pT spectrum of triggered particles in same evt
		TH1F			*fHistPtMixed[kCentBin][kZvertBin]; //! pT spectrum of triggered particles in mixed evt

		// Histograms for 2PC
		// cent : 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 (10 bins)
		// zvert : -8, -6, -4, -2, 0, 2, 4, 6, 8 (8 bins)
		// pT : 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 3,0, 4.0, 5.0, 10.0, 15.0 (10 bins)
		// PID : 0(hadron), 1(pion), 2(kaon), 3(proton)
		TH2F            *fHistdEtadPhiSame[kCentBin][kZvertBin][kpTBin][kpTBin][kPID]; //!
		TH2F            *fHistdEtadPhiMixed[kCentBin][kZvertBin][kpTBin][kpTBin][kPID]; //!

		// Histograms for PID QA
		TH2F            *fHistNsigmaTPCpT[kCentBin];       //! nsigma_{TPC} vs pT 
		TH2F            *fHistNsigmaTOFpT[kCentBin];       //! nsigma_{TOF} vs pT 
		TH2F            *fHistNsigmaITSpT[kCentBin];       //! nsigma_{ITS} vs pT 
		TH2F			*fHistNsigmaITSTPC[kCentBin][kpTBin];       //! nsigma_{ITS} vs nsigma_{TPC}
		TH2F			*fHistNsigmaTPCTOF[kCentBin][kpTBin];       //! nsigma_{TPC} vs nsigma_{TOF}
		TH2F            *fHistNsigmaITSTOF[kCentBin][kpTBin];       //! nsigma_{ITS} vs nsigma_{TOF}

//AliAnalysisTaskCorrHadrons(const AliAnalysisTaskCorrHadrons&); // not implemented
//AliAnalysisTaskCorrHadrons& operator=(const AliAnalysisTaskCorrHadrons&); // not implemented

ClassDef(AliAnalysisTaskEPCorrPP, 1);
};



//_____ Reduced Tracks -- contains only quantities requires for this analysis to reduce memory consumption for event mixing
class AliCorrReducedTrackPP : public AliVParticle // TObject
{
	public:
		AliCorrReducedTrackPP() {}
		AliCorrReducedTrackPP(Int_t partID, Double_t eta, Double_t phi, Double_t pt, Double_t zv, Short_t charge)
			: fParticleIDReduced(partID), fEtaReduced(eta), fPhiReduced(phi), fPtReduced(pt), fZvReduced(zv), fChargeReduced(charge) {}
		~AliCorrReducedTrackPP() {}

		// AliVParticle functions
		virtual Double_t  Px() const { return 0; }
		virtual Double_t  Py() const { return 0; }
		virtual Double_t  Pz() const { return 0; }
		virtual Double_t  Pt() const { return fPtReduced; }
		virtual Double_t  P() const   { return 0; }
		virtual Bool_t        PxPyPz(Double_t[3]) const { return 0; }
		virtual Double_t  Xv() const  { return 0; }
		virtual Double_t  Yv() const  { return 0; }
		virtual Double_t  Zv() const  { return fZvReduced; }
		virtual Bool_t    XvYvZv(Double_t[3]) const { return 0; }
		virtual Double_t  OneOverPt() const   { return 0; }
		virtual Double_t  Phi() const { return fPhiReduced; }
		virtual Double_t  Theta() const   { return 0; }
		virtual Double_t  E() const { return 0; }
		virtual Double_t  M() const   { return 0; }
		virtual Double_t  Eta() const { return fEtaReduced; }
		virtual Double_t  Y() const   { return 0; }
		virtual Short_t   Charge() const  { return fChargeReduced; }

		// void Print() { Printf(Form("Reduced track, eta: %lf phi: %lf pt: %lf p: %lf", fEtaReduced, fPhiReduced, fPtReduced, fPReduced)); }

		//________ PID
		Int_t             GetMyPartID() const { return fParticleIDReduced; }
		virtual Bool_t            IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }
		virtual Int_t         GetLabel() const { return 0; }
		virtual Int_t         PdgCode() const { return 0; }
		virtual const Double_t*   PID() const { return 0; }

	private:

		Int_t     fParticleIDReduced; // particle ID 
		Double_t  fEtaReduced;            // eta
		Double_t  fPhiReduced;        // phi
		Double_t  fPtReduced;         // pT
		Double_t  fZvReduced;			// z vertex
		Short_t   fChargeReduced;     // charge

		ClassDef(AliCorrReducedTrackPP, 1); // reduced track class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};

#endif
