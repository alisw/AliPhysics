/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskDiHadCorrelHighPt_H
#define AliAnalysisTaskDiHadCorrelHighPt_H

#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "AliEventCuts.h"

class AliPIDResponse;
class AliEventPoolManager;
class THnSparse;
class AliAODv0;
class AliAODVertex;

class AliAnalysisTaskDiHadCorrelHighPt : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskDiHadCorrelHighPt();
                                AliAnalysisTaskDiHadCorrelHighPt(const char *name, Bool_t analysisMC);
        virtual                 ~AliAnalysisTaskDiHadCorrelHighPt();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        
        Bool_t 					IsMyGoodPrimaryTrack(const AliAODTrack* aodtrack);
		Bool_t 					IsMyGoodV0AngleLambda(const AliAODv0 *t, AliAODVertex *pv); 
	    Bool_t 					IsMyGoodV0AngleK0(const AliAODv0 *t, AliAODVertex *pv); 
		Bool_t					IsMyGoodV0RapidityK0(const AliAODv0 *t); 
		Bool_t					IsMyGoodV0RapidityLambda(const AliAODv0 *t); 
		Bool_t 					IsMyGoodLifeTimeK0(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0); 
		Bool_t 					IsMyGoodLifeTimeLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0);
        Bool_t                  IsMyGoodLifeTimeAntiLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0);
        Bool_t                  IsMyGoodPID(const AliAODTrack *TrackPos, const AliAODTrack *TrackNeg);
        Bool_t                  IsMyGoodDaughterTrack(const AliAODTrack *t)  ;
	    Bool_t                  IsMyGoodV0(const AliAODv0 *v0,const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg, Int_t oSta);
        Bool_t                  IsMyGoodV0Topology(const AliAODv0 *v0);

        Int_t                   GetOStatus() { return fOStatus; }
        void                    SetOStatus(Int_t stat) {  fOStatus=stat; }
        void                    SetCutsCrosscheck(Bool_t stat) {  fCutsCrosscheck=stat; }

        void                    SetPtTrigMin(Double_t var) {fPtTrigMin=var;}
        void                    SetPtAsocMin(Double_t var) {fPtAsocMin=var;}
        void                    Corelations(TObjArray *triggers, TObjArray *associated, THnSparse * fHistKor, Double_t lPVz,THnSparse* fHistNumOfTrig,Bool_t hh,Bool_t V0h,Float_t perc,TH3F *fHistPtHard, Double_t ptHard,Bool_t hV0);
        void                    CorelationsMixing(TObjArray *triggers, TObjArray *bgTracks, THnSparse * fHistKor, Double_t lPVz,Float_t perc);
        void                    TopologCuts(THnSparse* fHist,Double_t pttrig,Double_t mass,Double_t dcaNeg, Double_t dcaPos,Double_t dcaDau, Double_t V0rad, Double_t cosPA,Double_t lifetime,Double_t massSell,Double_t triggType, Double_t status);
    void                    FillMC(const AliAODv0 *V0,TClonesArray *mcArray,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, TObjArray * selectedMCV0Triggersrec,THnSparse * fHistRecV0, TH3F * fHistMassPtCut,Double_t lPVz, const AliAODTrack * myTrackPos,const AliAODTrack * myTrackNeg,Bool_t status,THnSparse * histPur, TObjArray * selectedMCV0assoc,TH3F * fHistresol);
        void                    CorelationsMixingV0h(TObjArray *bgTracks, TObjArray *assocArray, THnSparse * fHistKor, Double_t lPVz, Float_t perc);
    
        void                    SetCosPAK0(Float_t cosPAK0) { fCosPointAngleK0 = cosPAK0; }
        void                    SetCosPALam(Float_t cosPALam) { fCosPointAngleLam = cosPALam; }
        void                    SetDCAV0Daughters(Float_t dca) { fDCAV0Daughters = dca; }
        void                    SetDCAposDaughter(Float_t dcapos) { fDCAposDaughter = dcapos; }
        void                    SetDCAnegDaughter(Float_t dcaneg) { fDCAnegDaughter = dcaneg; }
        void                    SetNTPCcrossedRows(Int_t Ncr) { fnumOfTPCcrossedRows = Ncr; }
        void                    SetEffAnalysis(Bool_t eff) { fEfficiency = eff; }
        void                    SetPurityCheckAnalysis(Bool_t pur) { fPurityCheck = pur; }
        void                    SetCorrelationsAnalysis(Bool_t corr) { fCorrelations = corr; }
        void                    SetMixedEvents(Int_t nmix) { fMixedEvents = nmix; }
        void                    SetNofPtBinsTrigger(Int_t nbins) { fNumberOfPtBinsTrigger = nbins; }
        void                    SetNofPtBinsAssoc(Int_t nbins) { fNumberOfPtBinsAssoc = nbins; }
        void                    SetV0RadiusCut(Double_t rad) { fV0Radius = rad; }
        void                    SetSigmaCut(Double_t sigma) { fSigmaCut = sigma; }
        void                    SetEtaCut(Double_t eta) { fEtaCut = eta; }
        void                    SetRapidityCut(Double_t rap) { fRapidityCut = rap; }
        void                    SetLifeTimeCutLambda(Double_t lifetime) { fLifeTimeLam = lifetime; }
        void                    SetLifeTimeCutK0(Double_t lifetime) { fLifeTimeK0 = lifetime; }
        void                    SetMassRejectCutK0(Double_t rej) { fMassRejectCutK0 = rej; }
        void                    SetLifeTimeCutLamda(Double_t rej) { fMassRejectCutLam = rej; }
        void                    SetEventPileUpCut(Bool_t cut) { fRejectEventPileUp = cut; }
        void                    SetTrackPileUpCut(Bool_t cut) { fRejectTrackPileUp = cut; }
    
        AliEventCuts            fAliEventCuts;
    
    private:
        AliAODEvent*            fAOD;           		//! input event
        AliPIDResponse*         fPIDResponse;           //!
        TList*                  fOutputList;    		//! output list
        TH3F*					fHistLambdaMassPtCut;	//! 
        TH3F*					fHistK0MassPtCut;		//!
        TH3F*                   fHistAntiLambdaMassPtCut; //!
        TH3F*                   fHistPtHard;                //!
        THnSparse*              fHistKorelacie;             //!
        THnSparse*              fHistdPhidEtaMix;           //!
        TH1D*                   fHistV0Multiplicity;         //!
        TH2D*                   fHistMultVtxz;              //!
        TH1D*                   fHistEtaAssoc;              //!
        TH1D*                   fHistEtaTrigg;              //!
        TH1D*                   fHistPhiAssoc;              //!
        TH1D*                   fHistPhiTrigg;              //!
        TH3D*                   fHistMCPtAs;            //!
        TH3D*                   fHistRCPtAs;            //!
        THnSparse*              fHistNumberOfTriggers;  //!
        THnSparse*              fHistMCKorelacie;       //!
        THnSparse*              fHistMCMixingRec;       //!

        Bool_t          			  fFillMixed;  // enable event mixing (default: ON)
        Int_t           			  fMixingTracks;      // size of track buffer for event mixing
        AliEventPoolManager*          fPoolMgr;         //! event pool manager
        AliEventPool*                 fPool; //!
        Bool_t                        fAnalysisMC; // enable MC study
        Int_t                         fOStatus; //
        Double_t                      fPtTrigMin; //
        Double_t                      fPtAsocMin; //
        Bool_t                        fCutsCrosscheck; //
        Int_t                         fMixedEvents; // number of minimum mixed events
        Double_t                      fV0Radius; // V0 radius
        Double_t                      fSigmaCut; // TPC PID
        Double_t                      fEtaCut; // pseudorapidity cut for primary hadrons
        Double_t                      fRapidityCut; // rapidity cut for V0
        Double_t                      fLifeTimeLam; // LifeTime cut for Lambdas and AnitLambdas
        Double_t                      fLifeTimeK0; // LifeTime cut for K0
        Double_t                      fMassRejectCutK0; // V0 competing rejection for K0
        Double_t                      fMassRejectCutLam;// V0 competing rejection for Lambdas and AntiLambdas
        Bool_t                        fRejectEventPileUp; // enable to use Pile-up cuts
        Bool_t                        fRejectTrackPileUp; // enable to use Bunch-Of Pile-up cuts

        THnSparse*              fHistKorelacieMCrec; //!  
        THnSparse*              fHistNumberOfTriggersGen;  //!
        THnSparse*              fHistNumberOfTriggersRec;  //!
    
        THnSparse*              fHistRecV0; //!
        THnSparse*              fHistGenV0; //!
    
        TH3D*                   fHistMCPtTrigg;            //!
        TH3D*                   fHistRCPtTrigg;            //!
        TH1D*                   fHistSelection;            //!
        TH3D*                   fHistPhiEtaGen;            //!
        TH3D*                   fHistPhiEtaRec;            //!
        TH1F*                   fHistMultipPercentile;     //!
        THnSparse*              fHistTopolCut;             //!
        THnSparse*              fHistTopolCutMC;           //!
        THnSparse*              fHistPurityCheck;          //!
    
        Float_t                 fCosPointAngleK0; //
        Float_t                 fCosPointAngleLam; //
        Float_t                 fDCAV0Daughters; //
        Float_t                 fDCAposDaughter; //
        Float_t                 fDCAnegDaughter; //
        Int_t                   fnumOfTPCcrossedRows; //
        Bool_t                  fEfficiency; //
        Bool_t                  fPurityCheck;//
        Bool_t                  fCorrelations;//
        TH3F*                   fHistPtResolution;//!
        Int_t                   fNumberOfPtBinsTrigger; //
        Int_t                   fNumberOfPtBinsAssoc;//
        TH1F *                  fHistV0MultiplicityK0; //!
        TH1F *                  fHistV0Lam;//!
        TH1F *                  fHistMultiplicityALam;//!
        TH2F *                  fHitsNTracks;//!
        THnSparse *             fHistPhiEta;//!

        AliAnalysisTaskDiHadCorrelHighPt(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented
        AliAnalysisTaskDiHadCorrelHighPt& operator=(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented

        ClassDef(AliAnalysisTaskDiHadCorrelHighPt, 11);
};

class AliV0ChParticle : public AliVParticle
{
    public:
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Bool_t status,Double_t mass)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fRecStatus(status), fMass(mass)
      {
      }
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label, Bool_t status)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fRecStatus(status)
      {
      }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh, Bool_t status)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh), fRecStatus(status)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg, Bool_t status,Double_t mass)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDpos(idpos), fIDneg(idneg), fRecStatus(status), fMass(mass)
    {
    }
     virtual ~AliV0ChParticle() {}
  
      // kinematics
      virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Pt() const { return fpT; }
      virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
      virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  
      virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
      virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
  
      virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Phi()        const { return fPhi; }
      virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
  
  
      virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
      virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
  
      virtual Double_t Eta()        const { return fEta; }
      virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
  
      virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
      virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
      // PID
      virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
      virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
  
      virtual Short_t WhichCandidate()      const { return fCandidate; }
    Int_t   MyLabel()             const { return fLabel; }
    Int_t   GetIDCh()             const { return fIDh; }
    Int_t   GetIDPos()            const { return fIDpos; }
    Int_t   GetIDNeg()            const { return fIDneg; }
    Bool_t  GetRecStatus()        const { return fRecStatus; }
    Double_t GetMass()            const { return fMass; }
  
      private:
      Double_t fEta;      // eta
      Double_t fPhi;      // phi
      Double_t fpT;       // pT
      Short_t fCandidate;   // V0 candidate: 1 - K0, 2 - Lambda, 3 - Antilambda, 4 - ChTrack
      Int_t fLabel;   // Label of MC particles
      Int_t fIDh;   // Label
      Int_t fIDpos;   // Label of possitive charged daughter
      Int_t fIDneg;   // Label of negative charged daughter
      Bool_t fRecStatus;   // reconstruction status
      Double_t fMass; // mass
  
      ClassDef( AliV0ChParticle, 2) // class required for correlatios calculation and event mixing
 };

#endif
