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
class AliESDEvent;
class AliESDtrack;

class AliAnalysisTaskDiHadCorrelHighPt : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskDiHadCorrelHighPt();
                                AliAnalysisTaskDiHadCorrelHighPt(const char *name, Bool_t analysisMC);
        virtual                 ~AliAnalysisTaskDiHadCorrelHighPt();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        
        Bool_t                  IsMyGoodPrimaryTrack(const AliAODTrack* aodtrack);
        Bool_t                  IsMyGoodPrimaryTrackESD(const AliESDtrack* esdtrack);
        Bool_t                  IsMyGoodV0AngleLambda(const AliAODv0 *t, AliAODVertex *pv); 
        Bool_t                  IsMyGoodV0AngleLambdaESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0AngleK0(const AliAODv0 *t, AliAODVertex *pv); 
        Bool_t                  IsMyGoodV0AngleK0ESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0RapidityK0(const AliAODv0 *t); 
        Bool_t                  IsMyGoodV0RapidityK0ESD(const AliESDv0 *t);
        Bool_t                  IsMyGoodV0RapidityLambda(const AliAODv0 *t); 
        Bool_t                  IsMyGoodV0RapidityLambdaESD(const AliESDv0 *t); 
        Bool_t                  IsMyGoodLifeTimeK0(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0); 
        Bool_t                  IsMyGoodLifeTimeLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0);
        Bool_t                  IsMyGoodLifeTimeESD(Double_t x,Double_t y, Double_t z, const AliESDv0 *V0, Double_t masshypotesis);
        Bool_t                  IsMyGoodLifeTimeAntiLambda(Double_t x,Double_t y, Double_t z, const AliAODv0 *V0);
        Bool_t                  IsMyGoodPID(const AliAODTrack *TrackPos, const AliAODTrack *TrackNeg);
        Bool_t                  IsMyGoodDaughterTrack(const AliAODTrack *t)  ;
        Bool_t                  IsMyGoodDaughterTrackESD(const AliESDtrack *t)  ;
        Bool_t                  IsMyGoodV0(const AliAODv0 *v0,const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg, Int_t oSta);
        Bool_t                  IsMyGoodV0ESD(const AliESDv0 *v0,const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg, Int_t oSta);
        Bool_t                  IsMyGoodV0Topology(const AliAODv0 *v0);
        Bool_t                  IsMyGoodV0TopologyESD(const AliESDv0 *v0, const AliESDtrack* myTrackPos, const AliESDtrack* myTrackNeg, Double_t lPVx, Double_t lPVy, Double_t lMagneticField);

        Int_t                   GetOStatus() { return fOStatus; }
        void                    SetOStatus(Int_t stat) {  fOStatus=stat; }
        void                    SetCutsCrosscheck(Bool_t stat) {  fCutsCrosscheck=stat; }

        void                    SetPtTrigMin(Double_t var) {fPtTrigMin=var;}
        void                    SetPtAsocMin(Double_t var) {fPtAsocMin=var;}
        void                    Corelations(TObjArray *triggers, TObjArray *associated, THnSparse * fHistKor, Double_t lPVz,THnSparse* fHistNumOfTrig,Bool_t hh,Bool_t V0h,Float_t perc,TH3F *fHistPtHard, Double_t ptHard,Bool_t hV0,Bool_t isMCGen);
        void                    CorelationsMixing(TObjArray *triggers, TObjArray *bgTracks, THnSparse * fHistKor, Double_t lPVz,Float_t perc);
        void                    TopologCuts(THnSparse* fHist,Double_t pttrig,Double_t mass,Double_t dcaNeg, Double_t dcaPos,Double_t dcaDau, Double_t V0rad, Double_t cosPA,Double_t lifetime,Double_t massSell,Double_t triggType, Double_t status);
        void                    FillMC(const AliVParticle *V0,TClonesArray *mcArray,Int_t pdgV0,Int_t pdgDau1, Int_t pdgDau2,Int_t triggerType, Double_t mass, TObjArray * selectedMCV0Triggersrec,THnSparse * fHistRecV0, TH3F * fHistMassPtCut,Double_t lPVz,THnSparse * histPur, TObjArray * selectedMCV0assoc,TH3F * fHistresol, Double_t posTrackProp[6],Double_t negTrackProp[6]);
        void                    CorelationsMixinghV0(TObjArray *bgTracks, TObjArray *assocArray, THnSparse * fHistKor, Double_t lPVz, Float_t perc);
    
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
        void                    SetTrackPileUpTOFCut(Bool_t cut) { fRejectTOF = cut; }
        void                    SetV0PileUpCut(Bool_t cut) { fRejectV0PileUp = cut; }
        void                    SetMultiplicityEstimator(TString est) { fMultEstimator = est; }
        void                    SetCorrelationsGen(Bool_t correl) { fCorrelationsGen = correl; }
        void                    SetCorrelationsV0h(Bool_t correl) { fV0hCorr = correl; }
        void                    SetCorrelationshh(Bool_t correl) { fhhCorr = correl; }
        void                    SetCorrelationshV0(Bool_t correl) { fhV0Corr = correl; }
        void                    SetFilterBit (Int_t filter) { fFilterBit = filter; }
        void                    SetRemoveHadronsFromV0 (Bool_t rem) { fRemoveHadrFromV0 = rem; }
        void                    SetRemoveLamhFromCascade (Bool_t rem) { fRemoveLamhFromCascade = rem; }
        void                    SetAcceptLambdasFromCascades (Bool_t accept) { fAacceptLambdasFromCasscade = accept; }
        void                    SetAcceptPurePrimHadrons (Bool_t accept) { fPurePrimHadrons = accept; }
        void                    SetAcceptPureV0 (Bool_t accept) { fPureV0 = accept; }
        Float_t                 GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius);
        void                    SetMergingCut(Double_t cut) { fMergingCut = cut; }
        void                    SetNumberOfDeltaPhiBins(Int_t bins) { fNumberOfDeltaPhiBins = bins; }
        void                    SetNumberOfDeltaEtaBins(Int_t bins) { fNumberOfDeltaEtaBins = bins; }
        void                    SetNumberOfEtaBins(Int_t bins) { fNumberOfEtaBins = bins; }
        void                    SetMixing(Bool_t mix) { fMixing = mix; } 
        void                    SetMixingGen (Bool_t mix) { fMixingGen = mix; }
        void                    SetNumberOfVzBins(Int_t nbins) { fNumOfVzBins = nbins; }
        AliAODTrack *           SetAliAODTrack(Double_t theta,Double_t phi, Double_t pt , Short_t charge);
        void                    SetPVCut(Int_t cut) { fPrimaryVertexCut = cut; }
        void                    SetESDtrackCuts(AliESDtrackCuts *const trcuts) { fESDTrackCuts = trcuts; }
        void                    SetAnalysisAOD(Bool_t aodAnalysis) { fAnalysisAOD = aodAnalysis; }
        void                    SetTestPions(Bool_t pions) { fTestPions = pions; }
        void                    SetMergingCutV0(Bool_t cut) { fMergingCutV0 = cut; }
        void                    SetMergingCutDaughters(Bool_t cut) { fMergingCutDaughters = cut; }

        AliEventCuts            fAliEventCuts;
    
    private:
        AliAODEvent*            fAOD;           		//! input event
        AliESDEvent*            fESD;                   //! input event
        AliMCEvent*             fmcEvent;               //! input MC event 
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
        TH3D*                   fHistMCPtAs;            //!
        TH3D*                   fHistRCPtAs;            //!
        THnSparse*              fHistNumberOfTriggers;  //!
        THnSparse*              fHistMCKorelacie;       //!
        THnSparse*              fHistMCMixingRec;       //!
        THnSparse*              fHistMCMixingGen;       //!

        Bool_t          			  fFillMixed;  // enable event mixing (default: ON)
        Int_t           			  fMixingTracks;      // size of track buffer for event mixing
        AliEventPoolManager*          fPoolMgr;         //! event pool manager
        AliEventPool*                 fPool; //!
        AliEventPool*                 fPoolMCGen; //!
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
        Bool_t                        fRejectTrackPileUp; // enable to use Bunch-Of Pile-up cuts for tracks

        THnSparse*              fHistKorelacieMCrec; //!  
        THnSparse*              fHistNumberOfTriggersGen;  //!
        THnSparse*              fHistNumberOfTriggersRec;  //!
    
        THnSparse*              fHistRecV0; //!
        THnSparse*              fHistGenV0; //!
    
        TH3D*                   fHistMCPtTrigg;            //!
        TH3D*                   fHistRCPtTrigg;            //!
        TH1D*                   fHistSelection;            //!
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
        Bool_t                  fRejectTOF; // enable to use Bunch-Of Pile-up cuts for tracks using TOF
        Bool_t                  fRejectV0PileUp;// enable to use Bunch-Of Pile-up cuts for V0
        TString                 fMultEstimator; // enable to change the Multiplicity estimator
        Bool_t                  fCorrelationsGen; // enable to compute the correlation function from generated MC particles
        Bool_t                  fV0hCorr; // enable to run V0-h correlations separately
        Bool_t                  fhhCorr; // enable to run h-h correlations separately
        Bool_t                  fhV0Corr; // enable to run h-V0 correlations separately
        Int_t                   fFilterBit; // enable to vary filter bit for systematic studies
        Bool_t                  fRemoveLamhFromCascade; // enable to remove hh corelations, which are from the same V0
        Bool_t                  fRemoveHadrFromV0; // enable to remove Lamh corelations, which are from the same cascade
        Bool_t                  fAacceptLambdasFromCasscade; // accept Lambdas from cascades to efficiency calculation 
        Bool_t                  fPurePrimHadrons; // anable to accept only pure primary hadrons for MC closure test
        Bool_t                  fPureV0; // anable to accept only pure good ID V0 for MC closure test
        Bool_t                  fAnalysisAOD;// anable to change between ESDs and AODs
        Double_t                fMergingCut; // cut for track spliting/merging
        Int_t                   fNumberOfDeltaPhiBins; // Number of DeltaPhi Bins in correlation function and in mixing
        Int_t                   fNumberOfDeltaEtaBins; // Number of DeltaEta Bins in correlation functionand in mixing
        Int_t                   fNumberOfEtaBins; // Number of Eta bins for efficiency correction
        Bool_t                  fMixing; // enable mixing
        Bool_t                  fMixingGen; //
        Double_t                fNumOfVzBins; // number of PV bins for mixing
        Int_t                   fPrimaryVertexCut; // PV position acceptance
        

        TH1D *                  fHistGenMultiplicity; //! 
        TH2D *                  fHistTPCTracksVsClusters; //!
        TH2D *                  fHistVZeroPercentileTPCMult;//!

        AliESDtrackCuts *       fESDTrackCuts; //
        Bool_t                  fTestPions; // for testing MC trains
        Bool_t                  fMergingCutV0; // enable merrging cut based on V0 
        Bool_t                  fMergingCutDaughters; // enable merrging cut based on V0 daughters

        AliAnalysisTaskDiHadCorrelHighPt(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented
        AliAnalysisTaskDiHadCorrelHighPt& operator=(const AliAnalysisTaskDiHadCorrelHighPt&); // not implemented

        ClassDef(AliAnalysisTaskDiHadCorrelHighPt, 20);
};

class AliV0ChParticle : public AliVParticle
{
    public:
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate,Double_t mass)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(-1),fIDh(-1),fIDpos(-1),fIDneg(-1),fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
      {
      }
      AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label)
        : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label),fIDh(-1),fIDpos(-1),fIDneg(-1),fMass(0),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
      {
      }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh),fIDpos(-1),fIDneg(-1),fMass(0),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Short_t charge, Double_t pz, Double_t energ)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(-1),fIDneg(-1),fMass(0),fCharge(charge), fPz(pz), fEnergie(energ),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh,Short_t charge, Double_t pz, Double_t energ)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(iDh),fIDpos(-1),fIDneg(-1),fMass(0),fCharge(0), fPz(pz), fEnergie(energ),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg,Double_t mass,Double_t phi1, Double_t pt1,Double_t eta1, Double_t char1, Double_t phi2,Double_t pt2,Double_t eta2, Double_t char2)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(phi1),fPt1(pt1),fEta1(eta1),fChar1(char1),fPhi2(phi2),fPt2(pt2),fEta2(eta2),fChar2(char2)
    {
    }
    AliV0ChParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg,Double_t mass)
    : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate), fLabel(label), fIDh(-1),fIDpos(idpos), fIDneg(idneg), fMass(mass),fCharge(0),fPz(0),fEnergie(0),fPhi1(0),fPt1(0),fEta1(0),fChar1(0),fPhi2(0),fPt2(0),fEta2(0),fChar2(0)
    {
    }
     virtual ~AliV0ChParticle() {}
  
      // kinematics
      virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
      virtual Double_t Pz() const { return fPz; }
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
  
  
      virtual Double_t E()          const { return fEnergie; }
      virtual Double_t M()          const { return fMass; }
  
      virtual Double_t Eta()        const { return fEta; }
      virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
  
      virtual Short_t Charge()      const { return fCharge; }
      virtual Int_t   GetLabel()    const { return fLabel; }
      // PID
      virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
      virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
  
      virtual Short_t WhichCandidate()      const { return fCandidate; }
    Int_t   GetIDCh()             const { return fIDh; }
    Int_t   GetIDPos()            const { return fIDpos; }
    Int_t   GetIDNeg()            const { return fIDneg; }
    Short_t GetCharge1()          const { return fChar1; }
    Short_t GetCharge2()          const { return fChar2; }
    Double_t GetPhi1()            const { return fPhi1; }
    Double_t GetPhi2()            const { return fPhi2; }
    Double_t GetPt1()             const { return fPt1; }
    Double_t GetPt2()             const { return fPt2; }
    Double_t GetEta1()            const { return fEta1; }
    Double_t GetEta2()            const { return fEta2; }
  
      private:
      Double_t fEta;      // eta
      Double_t fPhi;      // phi
      Double_t fpT;       // pT
      Short_t fCandidate;   // V0 candidate: 1 - K0, 2 - Lambda, 3 - Antilambda, 4 - ChTrack
      Int_t fLabel;   // Label of MC particles
      Int_t fIDh;   // Label
      Int_t fIDpos;   // Label of possitive charged daughter
      Int_t fIDneg;   // Label of negative charged daughter
      Double_t fMass; // mass
      Short_t fCharge; // charge of the track
      Double_t fPz; //pZ
      Double_t fEnergie; //E
      Double_t fPhi1; // phi of first daughter 
      Double_t fPt1; //pt of first daughter 
      Double_t fEta1; // eta of first daughter  
      Double_t fChar1; // charge of first daughter 
      Double_t fPhi2; // phi of second daughter 
      Double_t fPt2; //pt of second daughter 
      Double_t fEta2; // eta of second daughter 
      Double_t fChar2; // charge of second daughter 
   
      ClassDef( AliV0ChParticle, 5) // class required for correlatios calculation and event mixing
 };

#endif
