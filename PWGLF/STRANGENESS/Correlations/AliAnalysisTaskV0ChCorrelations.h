/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved.  
 * See cxx source for full Copyright notice  
 *
 * AliAnalysisTaskV0ChCorrelations class
 *
 * The task selects candidates for K0s, Lambdas and AntiLambdas (trigger particles)
 * and calculates correlations with charged unidentified particles (associated particles) in phi and eta. 
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 * Edited by Marek Bombara, last update January 2013, Marek.Bombara@cern.ch
 */

#ifndef ALIANALYSISTASKV0CHCORRELATIONS_H
#define ALIANALYSISTASKV0CHCORRELATIONS_H

class TH1F;
class TH1D;
class TH2F;
class THnSparse;
class TList;
class AliPIDResponse;
class AliEventPoolManager;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskV0ChCorrelations : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskV0ChCorrelations(const char *name = "AliAnalysisTaskV0ChCorrelations");
   //AliAnalysisTaskV0ChCorrelations(const AliAnalysisTaskV0ChCorrelations&);            //not implemented
   //AliAnalysisTaskV0ChCorrelations& operator=(const AliAnalysisTaskV0ChCorrelations&); //not implemented 
   virtual ~AliAnalysisTaskV0ChCorrelations();

   // Setting the global variables
   void SetAnalysisMC(Bool_t AnalysisMC = kTRUE) {fAnalysisMC = AnalysisMC;}
   void SetDcaDToPV(Float_t DcaDToPV = 0.5) {fDcaDToPV = DcaDToPV;}
   void SetDcaV0D(Float_t DcaV0D = 0.1) {fDcaV0D = DcaV0D;}

   // Getting the global variables
   Float_t GetDcaDToPV() { return fDcaDToPV; }
   Float_t GetDcaV0D() { return fDcaV0D; }

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     Terminate(Option_t *);

   Bool_t IsMyGoodPrimaryTrack(const AliAODTrack* aodtrack);
   Bool_t IsMyGoodDaughterTrack(const AliAODTrack* aodtrack);
   Bool_t IsMyGoodV0(const AliAODEvent* aod, const AliAODv0* aodv0, const AliAODTrack* tr1, const AliAODTrack* tr2);

private:

   AliAnalysisTaskV0ChCorrelations(const AliAnalysisTaskV0ChCorrelations&);            //not implemented
   AliAnalysisTaskV0ChCorrelations& operator=(const AliAnalysisTaskV0ChCorrelations&); //not implemented 

   Bool_t 		   fAnalysisMC; // enable MC study
   Bool_t          fFillMixed;  // enable event mixing (default: ON)
   Int_t           fMixingTracks;      // size of track buffer for event mixing
   AliEventPoolManager*          fPoolMgr;         //! event pool manager

   Float_t         fDcaDToPV;   // DCA of the daughter to primary vertex
   Float_t         fDcaV0D;     // DCA between daughters

   TList           *fOutput;        // Output list
   AliPIDResponse  *fPIDResponse;   // PID response

   TH2F            *fHistCentVtx; // centrality vs. z vertex - to see statistics for mixing
   TH1F            *fHistMultiMain; // multiplicity of main events
   
   THnSparseF      *fHistMassK0;       // K0 mass
   THnSparseF      *fHistMassLambda;     // Lambda mass
   THnSparseF      *fHistMassAntiLambda;   // AntiLambda mass
   
   THnSparseF	   *fHistdPhidEtaSib;   // dPhi vs. dEta, same event
   THnSparseF	   *fHistdPhidEtaMix;   // dPhi vs. dEta, mixed events
   THnSparseF	   *fHistTrigSib;   // pt of trigger particles, same event
   THnSparseF	   *fHistTrigMix;   // pt of trigger particles involved in mixing

   TH2D	           *fHistMCPtCent;   // pt vs. centrality of MC particles
   TH2D	           *fHistRCPtCent;   // pt vs. centrality of reconstructed tracks
   
   TH1D			   *fHistTemp;   // temporary histogram for debugging
   TH1D			   *fHistTemp2;   // temporary histogram for debugging

   ClassDef(AliAnalysisTaskV0ChCorrelations, 1); // class for V0Ch correlation analysis
};

/*  AliV0ChBasicParticle class contains only quantities 
 *	required for the analysis in order to reduce memory consumption for event mixing.
 */
class AliV0ChBasicParticle : public AliVParticle
{
  public:
    AliV0ChBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t candidate)
      : fEta(eta), fPhi(phi), fpT(pt), fCandidate(candidate)
    {
    }
    virtual ~AliV0ChBasicParticle() {}

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

  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCandidate;   // V0 candidate: 1 - K0sig, 2 - Lamsig, 3 - Alamsig, 4 - K0bg, 5 - Lambg, 6 - Alambg

    ClassDef( AliV0ChBasicParticle, 1); // class required for event mixing
};

#endif

