/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved.  
 * See cxx source for full Copyright notice  
 *
 * AliAnalysisTaskHadronCascadeCorrelations class
 *
 * The task selects candidates for V0 and cascade (associated particles)
 * and calculates correlations with unidentified charged trigger particles in phi and eta. 
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 * Authored by Zhongabo Yin,  Zhongbao.Yin@cern.ch
 */

#ifndef ALIANALYSISTASKLAMBDAK0S_H
#define ALIANALYSISTASKLAMBDAK0S_H
#include "AliVParticle.h"
class TH1F;
class TH1F;
class TH2F;
class THnSparse;
class TList;
class AliPIDResponse;
class AliEventPoolManager;
class AliVParticle;
class TObjArray;
class AliAnalysisUtils;
class AliEventCuts;
#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif


class AliAnalysisTaskLambdaK0s : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLambdaK0s();
  AliAnalysisTaskLambdaK0s(const char *name,
				     Double_t centMin, Double_t centMax,
				     Bool_t effCorr);
  virtual ~AliAnalysisTaskLambdaK0s();
  AliEventCuts *fEventCuts;
  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);

  void SetDCAToPrimVtxCut(Double_t value){fDCAToPrimVtx = value;}//cut on DCA of daughter track to Primary Vertex
  void SetDcaV0DaughtersCut(Double_t value){fDCABetweenDaughters = value;}
   
    void SetDCANegtoPrimVertexMinK0s(Double_t value){fDCANegtoPrimVertexMinK0s = value;}
    void SetDCAPostoPrimVertexMinK0s(Double_t value){fDCAPostoPrimVertexMinK0s = value;;}
    void SetDCANegtoPrimVertexMinLambda(Double_t value){fDCANegtoPrimVertexMinLambda = value;}
    void SetDCAPostoPrimVertexMinLambda(Double_t value){fDCAPostoPrimVertexMinLambda = value; }
    void SetDCANegtoPrimVertexMinAntiLambda(Double_t value){fDCANegtoPrimVertexMinAntiLambda = value; }
    void SetDCAPostoPrimVertexMinAntiLambda(Double_t value){fDCAPostoPrimVertexMinAntiLambda = value; }
    void SetPLtimeK0s( Double_t value){fPLtimeK0s=value;   }
    void SetPLtimeLambda( Double_t value){fPLtimeLambda=value;   }
   
  void SetCPA(Double_t value){fCPA = value;}
  void SetLambdaCPA(Double_t value){fLambdaCPA = value;}
  void SetLambdaDCA(Double_t value){fLambdaDCA = value;}
  void Set2DFiducial(Double_t value){f2DFiducial = value;}
  void SetTPCPidNsigmaCut(Float_t value){fNsigma = value;}
  void SetTPCNcls(Double_t value){fTPCNcls = value;}
  void SetMinCtau(Double_t value){fMinCtau = value;}
  void SetMaxCtau(Double_t value){fMaxCtau = value;}
  void SetMaxEta(Double_t value){fEtaCut = value;}
  void SetV0Eta(Double_t value){fV0Eta = value;}
  void SetV0DaughterEtaCut(Double_t value){fV0DaughterEtaCut = value;}
  void SetV0DaughterPtCut(Double_t value){fV0DaughterPtMinCut = value;}
  void SetAssoPtMin(Double_t value){fPtAssoMin = value;}
  void SetAssoPtMax(Double_t value){fPtAssoMax = value;}
  void SetPtTrigMin(Double_t value){fPtTrigMin = value;}
  void SetPtTrigMax(Double_t value){fPtTrigMax = value;}

void SetMixingTracks(Int_t value){fMixingTracks = value;}
void SetMixingPoolSize(Int_t value){fPoolSize =value ;}
void SetMCpileupgenandrc(Bool_t  gen=kFALSE ,Bool_t rc=kFALSE ){fIsPileUpCuts=gen; fRejectTrackPileUp=rc;}
void SetVtxXMin(Double_t value){fVtxXMin = value;}
void SetVtxYMin(Double_t value){fVtxYMin = value;}
void SetVtxZMin(Double_t value){fVtxZMin = value;}
void SetPrimVertexCut(Double_t value){fPrimVertexCut = value;}

  void SetAnalysisChXi(Bool_t AnalysisChXi = kFALSE)  {fChXi = AnalysisChXi;}
  void SetAnalysisChV0(Bool_t AnalysisChV0 = kFALSE)  {fChV0 = AnalysisChV0;}
  void SetAnalysisMC(Bool_t AnalysisMC = kFALSE) {fAnalysisMC = AnalysisMC;}
  void SetEventMixing(Bool_t eventMixing = kFALSE){fEventMixing = eventMixing;}
  void SetUseHybridGlobalTrack(){fUseHybridGlobalTrack = kFALSE;}

  void SetRapidity(Bool_t Rapidity = kFALSE){fRapidity = Rapidity;}


  
 private:
  AliAnalysisTaskLambdaK0s(const AliAnalysisTaskLambdaK0s&);            //not implemented
  AliAnalysisTaskLambdaK0s& operator=(const AliAnalysisTaskLambdaK0s&); //not implemented 

  Bool_t IsGoodPrimaryTrack(const AliAODTrack* aodtrack);
  Bool_t IsGoodDaughterTrack(const AliAODTrack* aodtrack);
  Bool_t IsGoodV0(AliAODv0 *v1);
  Bool_t IsMCV0Primary(AliAODv0 *v0, Int_t specie = 0); //0 for K0S, 1 for Lambda, 2 for AntiLambda
  Bool_t IsMCV0FromXi(AliAODv0 *v0, Int_t specie = 1); // 1 for Lambda, 2 for AntiLambda
  Bool_t IsMCXiPrimary(AliAODcascade *xi);//////////////////////////////// Cascade 
  Int_t IsMcV0(AliAODv0 *v0, Int_t specie = 0) const; // check whether the V0 is true one
  Int_t GetV0Label(Int_t lab1, Int_t lab2, Int_t specie) const;

  

  void AddQAEvent();
  void AddAnalysisTrk();
  void AddAnalysisK0s();
  void AddAnalysisLambda();
  void AddAnalysisAntiLambda();
  void AddAnalysisXiMinus();
  void AddQAV0Candidates();
  static const int    kNVtxZ = 8;   


  Int_t fMixingTracks; // size of track buffer for event mixing
  Int_t fPoolSize; 
  AliEventPoolManager* fPoolMgr; //! event pool manager

  TList *fEffList; //!
  TH2F * fHistEffEtaPtK0s; //! [kNVtxZ];
  TH2F * fHistEffEtaPtLambda; //! [kNVtxZ];
  TH2F * fHistEffEtaPtAntiLambda; //! [kNVtxZ];
  TH2F * fHistEffEtaPtXiMinus; //! [kNVtxZ];

  TH1F * fHistEffPtXiMinus; //! [kNVtxZ];

  TList *fOutput; //! Output list
  AliPIDResponse  *fPIDResponse;   // PID response
  
  TH2F * fHistArmenterosPodolanski;//!


  TH2F * fHistRCTrkPtEta; //! reconstructed track pt vs eta for eff. study
  TH2F * fHistRCPriTrkPtEta; //! rec. prim. track pt vs eta for eff. study
  TH2F * fHistRCSecTrkPtEta; //! rec. sec. track pt vs eta for contamination 
  TH2F * fHistMCtruthTrkPtEta; //! MC truth track pt vs eta for eff. study
  
  TH2F * fHistRCK0sPt; //! reconstructed K0s pt for eff.
  TH2F * fHistMCtruthK0sPt; //! MC truth K0s pt for eff.

TH2F * fHistRCXiMinusPt; //!
  TH2F * fHistMCtruthXiMinusPt;//!
//---PT XI
TH1F * fHistRCXiPt; //!
 TH1F * HistRCXiPt; //!
 
 TH1F * fHistMCXiMinusPt;//!

//Matrix 
 TH2F * fHistXiPtVSLambdaPt; //!
  TH2F * fHistXiPtVSAntiLambdaPt;//!

//------------------------------------------- 
  TH2F * fHistRCLambdaPt; //! 
  TH2F * fHistMCtruthLambdaPt; //! 
  
  TH2F * fHistallRCLambdaPt; //! 
  TH2F * fHistallDCALambdaPt; //! 
  TH2F * fHistallDLLambdaPt; //! 
  TH2F * fHistallRLambdaPt; //! 
  
  TH2F * fHistRCAntiLambdaPt; //! 
  TH2F * fHistMCtruthAntiLambdaPt; //! 

  TH2F * fHistRCPrimLambdaDCAPt; //!
  TH2F * fHistRCPrimLambdaDLPt; //! 
  TH2F * fHistRCPrimLambdaCTPt; //! 
  TH2F * fHistRCPrimLambdaRPt;//!  
 
  TH2F * fHistRCSecLambdaDCAPt; //! 
  TH2F * fHistRCSecLambdaDLPt; //! 
  TH2F * fHistRCSecLambdaCTPt; //! 
  TH2F * fHistRCSecLambdaRPt; //! 
 
 
  TH2F * fHistRCPrimAntiLambdaDCAPt; //!
  TH2F * fHistRCSecAntiLambdaDCAPt; //!



  TH3F * fHistRCTrk[kNVtxZ]; //! reconstructed track pt, eta, phi for eff. study
  TH3F * fHistRCPriTrk[kNVtxZ]; //! rec. prim. track pt, eta, phi for eff. study 
  TH3F * fHistRCSecTrk[kNVtxZ]; //! rec. sec. track pt, eta, phi for contamination
  TH3F * fHistMCtruthTrk[kNVtxZ]; //! MC truth track pt, eta, phi for eff. study

  TH3F * fHistRCK0s[kNVtxZ]; //! reconstructed K0s pt, eta, phi for eff.
  TH3F * fHistMCtruthK0s[kNVtxZ]; //! MC truth K0s pt, eta, phi for eff. 

  TH3F * fHistRCLambda[kNVtxZ]; //! reconstructed Lambda pt, eta, phi
  TH3F * fHistMCtruthLambda[kNVtxZ]; //! MC truth Lambda pt, eta, phi

  TH3F * fHistRCAntiLambda[kNVtxZ]; //! reconstructed Anti-L pt, eta, phi
  TH3F * fHistMCtruthAntiLambda[kNVtxZ]; //! MC truth Anti-L pt, eta, phi
//-----------------------------------------------------------------------------------------
 TH3F * fHistRCXiMinus[kNVtxZ]; //! reconstructed Xi- pt, eta, phi for eff.
 TH3F * fHistMCtruthXiMinus[kNVtxZ]; //! MC truth Xi- pt, eta, phi for eff.

//---------------------------------------------------------------------------------
  
  //TH2F * fHistRCXiPt; //! reconstructed Xi pt for eff.
  //TH2F * fHistMCtruthXiPt; //! MC truth Xi pt for eff.
    Bool_t fRejectTrackPileUp;
    Bool_t fIsPileUpCuts;
    
  Bool_t fEffCorr; //
  Bool_t fChXi;
  Bool_t fChV0;
  Bool_t fAnalysisMC; // enable MC study 
  Bool_t fEventMixing; // enable event mixing
  Bool_t fUseHybridGlobalTrack;
  Bool_t fRapidity;
  Double_t fV0Eta;



  AliAnalysisUtils * fUtils;

  Double_t fMinCent, fMaxCent; 
  Double_t fPrimVertexCut;
   Double_t        fVtxXMin;
    Double_t        fVtxYMin;
	 Double_t        fVtxZMin; 

  Double_t fDCAToPrimVtx;
  Double_t fDCABetweenDaughters;
  Double_t fPLtimeK0s;
  Double_t fPLtimeLambda;    
    Double_t        fDCANegtoPrimVertexMinK0s;
    Double_t        fDCAPostoPrimVertexMinK0s;
    Double_t        fDCANegtoPrimVertexMinLambda;
    Double_t            fDCAPostoPrimVertexMinLambda;
    Double_t            fDCANegtoPrimVertexMinAntiLambda;
    Double_t            fDCAPostoPrimVertexMinAntiLambda;


    
  Double_t fCPA;
  Double_t fLambdaCPA;
  Double_t fLambdaDCA;
  Double_t f2DFiducial;
  
  Float_t fNsigma; //for TPC dEdx pid
  Double_t fTPCNcls; // number of TPC clusters
  Double_t fMinCtau;
  Double_t fMaxCtau;
  Double_t fEtaCut;
  Double_t fV0DaughterEtaCut;
  Double_t fV0DaughterPtMinCut;
  
  Double_t fPtAssoMin; //V0 pt min
  Double_t fPtAssoMax; // V0 pt max
 
  Double_t fPtTrigMin; //V0 pt min
  Double_t fPtTrigMax; // V0 pt max



  
  Double_t fMassMean[2]; //mass mean for K0s, Lambda 
  Double_t fMassRes[2];  //mass resolution for K0s, Lambda

  Double_t fBestPrimaryVtxPos[3];

  TClonesArray* fMCArray; //! MC array for AOD 
  TObjArray       *selectedK0s;//!
  TObjArray       *selectedLambda;//!
  TObjArray       *selectedAntiLambda;//!
  TObjArray       *selectedTracks;//!
  TObjArray       *assocParticles;//!
  TObjArray       *selectedXiMinus;//!  ////////////////////////////


  ClassDef(AliAnalysisTaskLambdaK0s, 1); // class for HadronCascade correlation analysis
};

class AliV0XiParticleall : public AliVParticle
{
 public:
  AliV0XiParticleall(){}
  AliV0XiParticleall(Double_t pt, Double_t phi, Double_t eta, Double_t mass, 
		  	Short_t candidate)
    :fPt(pt), fPhi(phi), fEta(eta), fMass(mass), 
     fCandidate(candidate)
    {
    }

    virtual ~AliV0XiParticleall() {}

    // kinematics                                             
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Pt() const { return fPt; }
    virtual Double_t P()  const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { return fMass; }

    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

    virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID                       
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

    //    virtual Double_t Dca2PV() const {return fDca2PV;}
    virtual Short_t WhichCandidate() const { return fCandidate; }

 private:
    Double_t fPt;         // pt 
    Double_t fPhi;        // phi
    Double_t fEta;        // eta 
    Double_t fMass;       // inv. mass
    //Double_t fDca2PV;     // DCA of V0 to prim. vertex
    Short_t fCandidate;   // Candidate: 0-K0s candidates, 1-Lambda candidates, 2-AntiLambda candidates,3-Xi

    ClassDef( AliV0XiParticleall, 1); // class required for event mixing
};


#endif

