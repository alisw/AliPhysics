/* Leading Charged Track+V0 Correlation.(Works for Real,Monte Carlo Data)
 *                            Sandun Jayarathna
 *                          University of Houston
 *                      sandun.pahula.hewage@cern.ch
 *****************************************************************************************/

#ifndef ALILEADINGV0CORRELATIONH
#define ALILEADINGV0CORRELATIONH

#include "AliAnalysisTask.h"
#include "AliUEHist.h"
#include "TString.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliPID.h"

class TList;
class TH2;
class AliAODEvent;
class AliEventPoolManager;
class AliEventPool;
class AliVParticle;
class AliPIDResponse;
class AliPID;
class AliAODv0;


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
	void SetTrigger(TString aTriggerMask){fTriggerMask = aTriggerMask;}
	void SetPrimeryVertexCut(Double_t pvzcut){fpvzcut = pvzcut;}
	void SetEatCut(Double_t  TrackEtaCut){fTrackEtaCut = TrackEtaCut;}
	void SetFilterBit(UInt_t  filterBit){fFilterBit = filterBit;}
	void SetTrigPtBinLimits(Double_t trigPtLow,Double_t trigPtHigh){
		ftrigPtLow = trigPtLow;
		ftrigPtHigh = trigPtHigh;
	}
	void SetAssocPtBinLimits(Double_t assocPtLow,Double_t assocPtHigh){
		fassocPtLow = assocPtLow;
		fassocPtHigh = assocPtHigh;
	}
	void SetMCAnalysis(Bool_t aAnalysisMC){fAnalysisMC=aAnalysisMC;}
	void SetCutV0Radius(Double_t aCutV0Radius){fCutV0Radius=aCutV0Radius;}
	void SetCutDCANToP(Double_t aCutDCANegToPV){fCutDCANegToPV=aCutDCANegToPV;}
	void SetCutDCAPToP(Double_t aCutDCAPosToPV){fCutDCAPosToPV=aCutDCAPosToPV;}
	void SetCutDCADaughters(Double_t aCutDCAV0Daughters){fCutDCAV0Daughters=aCutDCAV0Daughters;}
	void SetCutCPA(Double_t aCutV0CosPA){fCutV0CosPA=aCutV0CosPA;}
	void SetCutArmenterosK0s(Double_t aSpecialArmenterosCutK0s){fSpecialArmenterosCutK0s=aSpecialArmenterosCutK0s;}
	void SetUSEPID(TString aUsePID){fUsePID=aUsePID;}
	void SetCutRap(Double_t aRapidityCut){fRapidityCut=aRapidityCut;}
	void SetCutCTauK0(Double_t aCTauK0){fCTauK0=aCTauK0;}
	void SetCutCTauLambda(Double_t aCTauLambda){fCTauLambda=aCTauLambda;}
	
	
private:
	AliLeadingV0Correlation(const  AliLeadingV0Correlation &det);
    AliLeadingV0Correlation&   operator=(const  AliLeadingV0Correlation &det);
	
	AliVParticle* ParticleWithCuts(TObject* obj, Int_t ipart);
	TObjArray* CloneAndReduceTrackList(TObjArray* tracks);
	TObjArray* FindLeadingObjects(TObject *obj);
	TObjArray* FindLeadingObjectsMC(TObject *obj);
	Int_t  NParticles(TObject* obj);
	Double_t RangePhi(Double_t DPhi);
	void FillCorrelations(TObjArray* particles, TObjArray* mixed,TH3F*histo,TH3F*leadinfo);
	void QSortTracks(TObjArray &a, Int_t first, Int_t last);
	Bool_t IsAcseptedPrimaryTrack(const AliAODTrack *itrack);
	Bool_t IsAcseptedDaughterTrack(const AliAODTrack *itrack);
	Bool_t IsK0InvMass(const Double_t mass)const;
	Bool_t IsLambdaInvMass(const Double_t mass) const;
	Bool_t IsTrackNotFromV0(AliAODTrack* track);
	Bool_t IsAcseptedV0(const AliAODEvent*aod,const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg);

    
	AliAODEvent              * fAODEvent;        //! AOD Event
	AliEventPoolManager      * fPoolMgr;         //! event pool manager
	AliEventPoolManager      * fPoolMgrMC;         //! event pool manager
	AliEventPool             * fPoolK0;          //! Pool for event mixing
	AliEventPool             * fPoolLambda;      //! Pool for event mixing
	AliEventPool             * fPoolK0MC;          //! Pool for event mixing
	AliEventPool             * fPoolLambdaMC;      //! Pool for event mixing
	AliPIDResponse           * fPIDResponse;     //  PID response
	
	Int_t fPoolMaxNEvents;                       // set maximum number of events in the pool
	Int_t fPoolMinNTracks;                       // set minimum number of tracks in the pool
	Int_t fMinEventsToMix;                       // set the minimum number of events want to mix
	Int_t fNzVtxBins;                            // number of z vrtx bins
	Double_t fZvtxBins[100];                     // [fNzVtxBinsDim]
	Int_t fNCentBins;                            // number of centrality bins
	Double_t fCentBins[100];                     // [fNCentBinsDim]
	
	TString         fcollidingSys;               // "PP" or "PbPb"
	TString         fTriggerMask;                // "Default"
	Double_t        fpvzcut;                     // PVz cut of event
    Double_t      	fTrackEtaCut;                // Eta cut on particles
    UInt_t         	fFilterBit;                  // Select tracks from an specific track cut (default 0xFF all track selected)
	Double_t        ftrigPtLow;
	Double_t        ftrigPtHigh;
	Double_t        fassocPtLow;
	Double_t        fassocPtHigh;
	Bool_t          fAnalysisMC;
	TString         fUsePID;
	Double_t        fRapidityCut;                // Rapidity cut V0
	
	
	//--- 5 Topological Selections And other quality cuts for V0
	Double_t        fCutV0Radius;
	Double_t        fCutDCANegToPV;
	Double_t        fCutDCAPosToPV;
	Double_t        fCutDCAV0Daughters;
	Double_t        fCutV0CosPA;
	Double_t        fSpecialArmenterosCutK0s;
	Double_t        fCTauK0;
	Double_t        fCTauLambda;
	
	
	TList       * fOutputList;           // Output list
	
	//---------------------MC--------------
	TH1F        *fHistMCPrimaryVertexX;       
	TH1F        *fHistMCPrimaryVertexY;       
	TH1F        *fHistMCPrimaryVertexZ;                       
	TH1F        *fHistMCPtAllK0s;
	TH1F        *fHistMCPtAllLambda;
	TH1F        *fHistMCPtAllAntiLambda;
	TH1F        *fHistMCRapK0s;                  
	TH1F        *fHistMCRapLambda;
	TH1F        *fHistMCRapAntiLambda;           
	TH1F        *fHistMCPtK0s;        
	TH1F        *fHistMCPtLambda;        
	TH1F        *fHistMCPtAntiLambda;
	TH1F        *fHistMCPtLambdaFromSigma;        
	TH1F        *fHistMCPtAntiLambdaFromSigma;        
	TH2F        *fHistPrimRawPtVsYK0s;
	TH2F        *fHistPrimRawPtVsYLambda;
	TH2F        *fHistPrimRawPtVsYAntiLambda;
	//-------------------REAL--------------
	TH1F        *fHistPrimaryVertexX;
	TH1F        *fHistPrimaryVertexY;
	TH1F        *fHistPrimaryVertexZ;
	TH2F        *fHistDcaPosToPrimVertexK0vsMassK0;   
	TH2F        *fHistDcaNegToPrimVertexK0vsMassK0;   
	TH2F        *fHistRadiusV0K0vsMassK0;             
	TH2F        *fHistDecayLengthV0K0vsMassK0;        
	TH2F        *fHistDcaV0DaughtersK0vsMassK0;       
	TH2F        *fHistCosPointAngleK0vsMassK0;        
	TH2F        *fHistDcaPosToPrimVertexLvsMassL;       
	TH2F        *fHistDcaNegToPrimVertexLvsMassL;       
	TH2F        *fHistRadiusV0LvsMassL;                  
	TH2F        *fHistDecayLengthV0LvsMassL;             
	TH2F        *fHistDcaV0DaughtersLvsMassL;          
	TH2F        *fHistCosPointAngleLvsMassL;                
	TH2F        *fHistDcaPosToPrimVertexAntiLvsMass;       
	TH2F        *fHistDcaNegToPrimVertexAntiLvsMass;       
	TH2F        *fHistRadiusV0AntiLvsMass;                  
	TH2F        *fHistDecayLengthV0AntiLvsMass;             
	TH2F        *fHistDcaV0DaughtersAntiLvsMass;          
	TH2F        *fHistCosPointAngleAntiLvsMass;             
	TH1F        *fHistMassK0;        
	TH1F        *fHistMassLambda;        
	TH1F        *fHistMassAntiLambda;                
	TH2F        *fHistPtVsMassK0;        
	TH2F        *fHistPtVsMassLambda;        
	TH2F        *fHistPtVsMassAntiLambda;
	TH2F        *fHistArmenterosPodolanskiK0; 
	TH2F        *fHistArmenterosPodolanskiLambda;
	TH2F        *fHistArmenterosPodolanskiAntiLambda;
	//-----------Associated----------------       
	TH1F        *fHistAsMcPtK0;        
	TH1F        *fHistAsMcPtLambda;        
	TH1F        *fHistAsMcPtAntiLambda;        
	TH1F        *fHistAsMcProdRadiusK0;        
	TH1F        *fHistAsMcProdRadiusLambda;        
	TH1F        *fHistAsMcProdRadiusAntiLambda;        
	TH2F        *fHistAsMcProdRadiusXvsYK0s;        
	TH2F        *fHistAsMcProdRadiusXvsYLambda;        
	TH2F        *fHistAsMcProdRadiusXvsYAntiLambda;        
	TH1F        *fHistPidMcMassK0;        
	TH1F        *fHistPidMcMassLambda;        
	TH1F        *fHistPidMcMassAntiLambda;                             
	TH1F        *fHistAsMcPtLambdaFromSigma;        
	TH1F        *fHistAsMcPtAntiLambdaFromSigma;        
	TH2F        *fHistAsMcSecondaryPtVsRapK0s;        
	TH2F        *fHistAsMcSecondaryPtVsRapLambda;        
	TH2F        *fHistAsMcSecondaryPtVsRapAntiLambda;        
	TH1F        *fHistAsMcSecondaryProdRadiusK0s;        
	TH1F        *fHistAsMcSecondaryProdRadiusLambda;        
	TH1F        *fHistAsMcSecondaryProdRadiusAntiLambda;        
	TH2F        *fHistAsMcSecondaryProdRadiusXvsYK0s;        
	TH2F        *fHistAsMcSecondaryProdRadiusXvsYLambda;        
	TH2F        *fHistAsMcSecondaryProdRadiusXvsYAntiLambda;        
	TH1F        *fHistAsMcSecondaryPtLambdaFromSigma;        
	TH1F        *fHistAsMcSecondaryPtAntiLambdaFromSigma;        
	//-----------Correlation----------------
	TH3F        * fHistSibK0;
	TH3F        * fHistMixK0;
	TH3F        * fHistSibLambda;
	TH3F        * fHistMixLambda;
	TH3F        * fHistSibK0MC;
	TH3F        * fHistMixK0MC;
	TH3F        * fHistSibLambdaMC;
	TH3F        * fHistMixLambdaMC;
	TH3F        * fHistLeadInfo;
	TH3F        * fHistLeadInfoMC;
	TH3F        * fHistLeadInfoMix;
	TH3F        * fHistLeadInfoMixMC;
	
	ClassDef(AliLeadingV0Correlation, 1); 
};
//---------------------------------------------------------------------------------------
class AliLeadingBasicParticle : public AliVParticle
{
  public:
    AliLeadingBasicParticle(Float_t eta, Float_t phi, Float_t pt)
      : fEta(eta), fPhi(phi), fpT(pt)
    {
    }
    virtual ~AliLeadingBasicParticle() {}

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
	
  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT

    ClassDef( AliLeadingBasicParticle, 1); // class required for event mixing
};

#endif

