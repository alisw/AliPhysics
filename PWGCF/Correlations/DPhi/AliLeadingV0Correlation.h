/* Leading Charged Track+V0 Correlation.(Works for Real,Monte Carlo Data)
 *                            Sandun Jayarathna
 *                          University of Houston
 *                      sandun.pahula.hewage@cern.ch
 *****************************************************************************************/

#ifndef ALILEADINGV0CORRELATIONH
#define ALILEADINGV0CORRELATIONH

#include "AliAnalysisTask.h"
#include "TString.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliPID.h"
#include "THnSparse.h"

class TList;
class AliAODEvent;
class AliEventPoolManager;
class AliVParticle;
class AliPIDResponse;
class AliPID;
class AliAODv0;
class AliAnalyseLeadingTrackUE;


#ifndef ALIANALYSISTASKSEH
#include "AliAnalysisTaskSE.h"
#endif

//---------------------------------------------------------------------------------------
class AliLeadingV0Correlation : public AliAnalysisTaskSE {
public:
   AliLeadingV0Correlation();
   AliLeadingV0Correlation(const char *name);
   virtual ~AliLeadingV0Correlation();

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     Terminate(Option_t *);
	
	void SetMaxNEventsInPool(Int_t events){fPoolMaxNEvents=events;}
	void SetMinNTracksInPool(Int_t tracks){fPoolMinNTracks=tracks;}
	void SetMinEventsToMix(Int_t events){fMinEventsToMix=events;}

	void SetPoolPVzBinLimits(Int_t Nzvtxbins,const Double_t *ZvtxBins){
		fNzVtxBins = Nzvtxbins;
		for(int ix = 0;ix<fNzVtxBins+1;ix++){fZvtxBins[ix] = ZvtxBins[ix];}
	}
	
	void SetPoolCentBinLimits(Int_t Ncentbins,const Double_t *CentBins){
		fNCentBins = Ncentbins;
		for(int ix = 0;ix<fNCentBins+1;ix++){fCentBins[ix] = CentBins[ix];}
	}
	
	void SetCollidingSystem(TString system){fcollidingSys = system;}
	void SetTriggerType(TString triggertype){ftriggertype = triggertype;}
	void SetPrimeryVertexCut(Double_t pvzcut){fpvzcut = pvzcut;}
	void SetEatCut(Double_t  TrackEtaCut){fTrackEtaCut = TrackEtaCut;}
	void SetFilterBit(UInt_t  filterBit){fFilterBit = filterBit;}
	void SetMCAnalysis(Bool_t aAnalysisMC){fAnalysisMC=aAnalysisMC;}
	void SetCase(Int_t aCase){fCase=aCase;}
	void SetRemovePileUp(Bool_t aRemovePileUP){fRemovePileUP=aRemovePileUP;}
	void SetRemoveAutoCorr(Bool_t aRemoveAutoCorr){fRemoveAutoCorr=aRemoveAutoCorr;}
	void SetCutRap(Double_t aRapidityCut){fRapidityCut=aRapidityCut;}
	void SetV0Radius(Double_t aV0radius){fV0radius=aV0radius;}
	void SetV0PostoPVz(Double_t aV0PostoPVz){fV0PostoPVz=aV0PostoPVz;}
	void SetV0NegtoPVz(Double_t aV0NegtoPVz){fV0NegtoPVz=aV0NegtoPVz;}
	void SetDCAV0Daughters(Double_t aDCAV0Daughters){fDCAV0Daughters=aDCAV0Daughters;}
	void SetCPAK0(Double_t aCPAK0){fCPAK0=aCPAK0;}
	void SetCPALam(Double_t aCPALam){fCPALam=aCPALam;}
	void SetRejectLamK0(Double_t aRejectLamK0){fRejectLamK0=aRejectLamK0;}
	void SetRejectK0Lam(Double_t aRejectK0Lam){fRejectK0Lam=aRejectK0Lam;}
	void SetSigmaPID(Double_t aSigmaPID){fSigmaPID=aSigmaPID;}
	void SetCTK0(Double_t aCutCTK0){fCutCTK0=aCutCTK0;}
	void SetCTLa(Double_t aCutCTLa){fCutCTLa=aCutCTLa;}
	void SetMassCutK0(Double_t aMassCutK0){fMassCutK0=aMassCutK0;}
	void SetMassCutLa(Double_t aMassCutLambda){fMassCutLa=aMassCutLambda;}
	
private:
	AliLeadingV0Correlation(const  AliLeadingV0Correlation &det);
    AliLeadingV0Correlation&   operator=(const  AliLeadingV0Correlation &det);
	
	Bool_t IsAcseptedDaughterTrack(const AliAODTrack *itrack);
	Bool_t IsAcseptedV0(const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg);
	Bool_t IsAcseptedK0(Double_t v0rad,
						Double_t dcaptp,
						Double_t dcantp,
						Double_t dcav0d,
						Double_t cpa,
						Double_t massLa,
						Double_t massALa);
	Bool_t IsAcseptedLA(Double_t v0rad,
						Double_t dcaptp,
						Double_t dcantp,
						Double_t dcav0d,
						Double_t cpa,
						Double_t massK0);
	Bool_t IsK0InvMass(const Double_t mass)const;
	Bool_t IsLambdaInvMass(const Double_t mass) const;
	Double_t RangePhi(Double_t DPhi);
	Bool_t IsTrackFromV0(AliAODTrack* track);
	void FillCorrelationSibling(Double_t MultipOrCent,
								TObjArray*triggerArray,TObjArray*selectedV0Array,
								THnSparse*triggerHist,THnSparse*associateHist);
	void FillCorrelationMixing(Double_t MultipOrCentMix,Double_t pvxMix,
							   Double_t poolmax,Double_t poolmin,
							   TObjArray*triggerArray,TObjArray*selectedV0Array,
							   THnSparse*triggerHist,THnSparse*associateHist);
	Bool_t IsFeedDownV0(UInt_t lIdxPosV0,
				  UInt_t lIdxNegV0,
				  Double_t lInvMassXiMinus,
				  Double_t lInvMassXiPlus,
				  TH1F*HistXiMinus);
	
	AliAODEvent              * fAODEvent;        //  AOD Event
	AliEventPoolManager      * fPoolMgr;         //  event pool manager for Event Mixing
	AliPIDResponse           * fPIDResponse;     //  PID response
	AliAnalyseLeadingTrackUE * fAnalyseUE;       //  Leading Track Underling Event
	
	Int_t fPoolMaxNEvents;                       // set maximum number of events in the pool
	Int_t fPoolMinNTracks;                       // set minimum number of tracks in the pool
	Int_t fMinEventsToMix;                       // set the minimum number of events want to mix
	Int_t fNzVtxBins;                            // number of z vrtx bins
	Double_t fZvtxBins[100];                     // [fNzVtxBinsDim]
	Int_t fNCentBins;                            // number of centrality bins
	Double_t fCentBins[100];                     // [fNCentBinsDim]
	
	TString         fcollidingSys;               // "PP" or "PbPb"
	TString         ftriggertype;                // "Leading" or "Alltrigs"
	Double_t        fpvzcut;                     // PVz cut of event
    Double_t      	fTrackEtaCut;                // Eta cut on particles
    UInt_t         	fFilterBit;                  // Select tracks from an specific track cut (default 0xFF all track selected)
	Bool_t          fAnalysisMC;                 // MC or Not
	Int_t           fCase;                       // Case number
	Bool_t          fRemovePileUP;               // 1Remove or 0 Not Remove
	Bool_t          fRemoveAutoCorr;             // 1Remove or 0 Not Remove
	
	Double_t        fRapidityCut;                // Rapidity cut V0
	Double_t        fV0radius;                   // Topological selection for systamatics
	Double_t        fV0PostoPVz;                 // Topological selection for systamatics
	Double_t        fV0NegtoPVz;                 // Topological selection for systamatics
	Double_t        fDCAV0Daughters;             // Topological selection for systamatics
	Double_t        fCPAK0;                      // Topological selection for systamatics
	Double_t        fCPALam;                     // Topological selection for systamatics
	Double_t        fRejectLamK0;                // selection for systamatics
	Double_t        fRejectK0Lam;                // selection for systamatics
	Double_t        fSigmaPID;                   // selection for systamatics
	Double_t        fCutCTK0;                    // selection for systamatics
	Double_t        fCutCTLa;                    // selection for systamatics
	
	Double_t        fMassCutK0;                  // selection for systamatics
	Double_t        fMassCutLa;                  // selection for systamatics
	
	Bool_t          fUseChargeHadrons;           // Only pi,k,and proton
	Double_t        fPtMin;                      // 0.15 
	
	
	TList       * fOutputList;                   // Output list
	
	THnSparse   *fHistEventViceGen;
	THnSparse   *fHistEventViceReconst;
	THnSparse	*fHistMCGenK0;
	THnSparse	*fHistMCGenLAM;
	THnSparse	*fHistMCGenALAM;
	THnSparse	*fHistReconstK0;
	THnSparse	*fHistReconstLA;
	THnSparse	*fHistReconstALA;
	THnSparse	*fHistMCAssoK0;
	THnSparse	*fHistMCAssoLA;
	THnSparse	*fHistMCAssoALA;
	
	THnSparse   *fHistReconstSib;
	THnSparse   *fHistReconstMix;
	THnSparse   *fHistReconstSibMC;
	THnSparse   *fHistReconstMixMC;
	THnSparse   *fHistReconstSibMCAssoc;
	THnSparse   *fHistReconstMixMCAssoc;
	THnSparse   *fHistTriggerSib;
	THnSparse   *fHistTriggerMix;
	THnSparse   *fHistTriggerSibMC;
	THnSparse   *fHistTriggerMixMC;

	ClassDef(AliLeadingV0Correlation, 1); 
};
//---------------------------------------------------------------------------------------
class V0Correlationparticle : public AliVParticle
{
public:
    V0Correlationparticle(Float_t eta, 
							Float_t phi, 
							Float_t pt, 
							Short_t candidate):
	  fEta(eta), 
	  fPhi(phi), 
	  fpT(pt), 
	  fCandidate(candidate)
    {
    }
    virtual ~V0Correlationparticle(){}
	
    virtual Double_t Px()                 const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Py()                 const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Pz()                 const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Pt()                 const { return fpT;}
    virtual Double_t P()                  const { AliFatal("Not implemented"); return 0;}
    virtual Bool_t   PxPyPz(Double_t[3])  const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Xv()                 const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Yv()                 const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Zv()                 const { AliFatal("Not implemented"); return 0;}
    virtual Bool_t   XvYvZv(Double_t[3])  const { AliFatal("Not implemented"); return 0;}
    virtual Double_t OneOverPt()          const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Phi()                const { return fPhi;}
    virtual Double_t Theta()              const { AliFatal("Not implemented"); return 0;}
    virtual Double_t E()                  const { AliFatal("Not implemented"); return 0;}
    virtual Double_t M()                  const { AliFatal("Not implemented"); return 0;}
    virtual Double_t Eta()                const { return fEta;}
    virtual Double_t Y()                  const { AliFatal("Not implemented"); return 0;}
    virtual Short_t  Charge()             const { AliFatal("Not implemented"); return 0;}
    virtual Int_t    GetLabel()           const { AliFatal("Not implemented"); return 0;}
    virtual Int_t    PdgCode()            const { AliFatal("Not implemented"); return 0;}
    virtual const    Double_t *PID()      const { AliFatal("Not implemented"); return 0;}
    virtual Short_t  WhichCandidate()     const { return fCandidate;}

	
private:
    Float_t  fEta;            // Eta
    Float_t  fPhi;            // Phi
    Float_t  fpT;             // pT
    Short_t  fCandidate;      // 1-K0,2-Lam,3-Alam
	
	
    ClassDef( V0Correlationparticle, 1);
};

#endif