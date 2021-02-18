/*
 *****************************************************************************************/
#ifndef ALIANALYSISTASKSEPBPBCORRELATIONSJETV2
#define ALIANALYSISTASKSEPBPBCORRELATIONSJETV2

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "TSpline.h"
#include "TGrid.h"
#include "AliEventPoolManager.h"
#include "AliTHn.h"


class AliEventPoolManager;
class AliOADBContainer;
class THnSparse;

//====================================================================================================================================================

class  AliAnalysisTaskSEPbPbCorrelationsJetV2 : public AliAnalysisTaskSE {

 public:
 
  AliAnalysisTaskSEPbPbCorrelationsJetV2();
  AliAnalysisTaskSEPbPbCorrelationsJetV2(const Char_t *name);
  virtual ~AliAnalysisTaskSEPbPbCorrelationsJetV2();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();

  // ------------- Cuts -----------------

  void SetRemovePileup(Bool_t flag) { fRemovePileup = flag; }
  void SetRemovePileup2(Bool_t flag) { fRemovePileup2 = flag; }
  void SetRemovePileup3(Bool_t flag) { fRemovePileup3 = flag; }
  void SetUseRes(Bool_t flag) { fUseRes = flag; }
  void SetAverageRes(Bool_t flag) {fAverageRes = flag; }
  void SetSubDivide(Int_t n1, Int_t n2) { fN1 = n1; fN2 = n2; }
  void SetESEcuts(Double_t low, Double_t high) { fLowESEcut = low; fHighESEcut = high; }
  void SetESEdet(Int_t det) { fESEdet = det; }
  void SetEP(Bool_t flag) { fEP = flag; }
  void SetMinHardPt(Double_t pt) { fMinHardPt = pt; };
  
  // ------------- Analysis -------------

  Float_t GetV0Multiplicity();
  Double_t GetITSMultiplicity();
  Bool_t IsTriggerFired();
  void SetPtBinning(Int_t nBins, Double_t *limits);
  void SetTrigPtBinning(Int_t nBins, Double_t *limits);
  void SetAssocPtBinning(Int_t nBins, Double_t *limits);
  void SetCentBinning(Int_t nBins, Double_t *limits);
  void SetEtaBinning(Int_t nBins, Double_t *limits);
  void SetZvtxBinning(Int_t nBins, Double_t *limits);
  void SetCentMethod(TString method) { fCentMethod = method; }
  
  void FillHistogramsV2(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
			Double_t resA2, Double_t resC2, Double_t resT2,
			Double_t resA2LowESE, Double_t resC2LowESE, Double_t resT2LowESE,
			Double_t resA2HighESE, Double_t resC2HighESE, Double_t resT2HighESE,
			Int_t index);
  void FillHistogramsV3(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
			Double_t resA3, Double_t resC3, Double_t resT3,
			Double_t resA3LowESE, Double_t resC3LowESE, Double_t resT3LowESE,
			Double_t resA3HighESE, Double_t resC3HighESE, Double_t resT3HighESE,
			Int_t index);
  void FillHistogramsdPhidEta(Int_t trIndex, Short_t charge, Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
			      Double_t resA2, Double_t resC2, Double_t resT2,
			      Double_t resA3, Double_t resC3, Double_t resT3);
  void FillHistogramsdPhidEtaMixed(Int_t trIndex, Short_t charge, Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
				   TObjArray *mixedTracks);
  Int_t GetCentBin();
  Int_t GetCentBin(Double_t percentile);
  Int_t GetZvtxBin(Double_t zvtx);

  Double_t CalcCorrectedPhi(Double_t phi, Double_t dPhi) const;

  void ProcessEvent(Double_t percentile, Int_t centBin, Int_t zvtxBin);

  void OpenInfoCalbration(Int_t run);
  Short_t GetVertexZ(Double_t vtxZ) const;
  Bool_t ComputeQ(AliAODEvent* aod, Double_t Zvtx);

  void CalcResolutions(Double_t percentile,
		       Double_t &resA2, Double_t &resC2, Double_t &resT2,
		       Double_t &resA2LowESE, Double_t &resC2LowESE, Double_t &resT2LowESE,
		       Double_t &resA2HighESE, Double_t &resC2HighESE, Double_t &resT2HighESE,
		       Double_t &resA3, Double_t &resC3, Double_t &resT3,
		       Double_t &resA3LowESE, Double_t &resC3LowESE, Double_t &resT3LowESE,
		       Double_t &resA3HighESE, Double_t &resC3HighESE, Double_t &resT3HighESE);
  
 private:

  static const Int_t fNMaxBinsCentrality = 25; //!
  static const Int_t fNMaxBinsPt = 50; //!
  static const Int_t fNMaxBinsAssocPt = 10; //!
  static const Int_t fNMaxBinsEta = 6; //!
  static const Int_t fNMaxBinsZvtx = 20; //!

  static AliOADBContainer *cont; //!
  static AliOADBContainer *contQx2am[14]; //!
  static AliOADBContainer *contQy2am[14]; //!
  static AliOADBContainer *contQx2as[14]; //!
  static AliOADBContainer *contQy2as[14]; //!
  static AliOADBContainer *contQx2cm[14]; //!
  static AliOADBContainer *contQy2cm[14]; //!
  static AliOADBContainer *contQx2cs[14]; //!
  static AliOADBContainer *contQy2cs[14]; //!
  static AliOADBContainer *contQx2trm[14]; //!
  static AliOADBContainer *contQy2trm[14]; //!
  static AliOADBContainer *contQx2trs[14]; //!
  static AliOADBContainer *contQy2trs[14]; //!
  static AliOADBContainer *contQx3am[14]; //!
  static AliOADBContainer *contQy3am[14]; //!
  static AliOADBContainer *contQx3as[14]; //!
  static AliOADBContainer *contQy3as[14]; //!
  static AliOADBContainer *contQx3cm[14]; //!
  static AliOADBContainer *contQy3cm[14]; //!
  static AliOADBContainer *contQx3cs[14]; //!
  static AliOADBContainer *contQy3cs[14]; //!
  static AliOADBContainer *contQx3trm[14]; //!
  static AliOADBContainer *contQy3trm[14]; //!
  static AliOADBContainer *contQx3trs[14]; //!
  static AliOADBContainer *contQy3trs[14]; //!
  
  static TSpline3*    fSplQ2V0A[90];           //! splines for q2 V0A
  static TSpline3*    fSplQ3V0A[90];           //! splines for q3 V0A
  static TSpline3*    fSplQ2V0C[90];           //! splines for q2 V0A
  static TSpline3*    fSplQ3V0C[90];           //! splines for q3 V0A
  static TSpline3*    fSplQ2trk[90];           //! splines for q2 tracklets
  static TSpline3*    fSplQ3trk[90];           //! splines for q3 tracklets
  static TSpline3*    fSplQ2V0AC[90];           //! splines for q2 V0A+V0C
  static TSpline3*    fSplQ3V0AC[90];           //! splines for q3 V0A+V0C

  AliAODEvent *fAOD; //!
  AliEventPoolManager *fPoolMgr; //! event pool manager
  TClonesArray *ftrk; //!
  
  Int_t fNbinsCent; 
  Int_t fNbinsPt;     
  Int_t fNbinsPtTrig; 
  Int_t fNbinsAssocPt; 
  Int_t fNbinsZvtx; 
  
  TAxis *fCentAxis; 
  TAxis *fPtAxis; 
  TAxis *fPtTrigAxis; 
  TAxis *fPtAssocAxis; 
  TAxis *fEtaAxis; 
  TAxis *fZvtxAxis; 

  Bool_t fRemovePileup; 
  Bool_t fRemovePileup2; 
  Bool_t fRemovePileup3; 
  Bool_t fUseRes; 
  Bool_t fAverageRes; 
  Double_t fMinHardPt; 
  
  Int_t fN1; 
  Int_t fN2; 

  Double_t fLowESEcut; 
  Double_t fHighESEcut; 
  Int_t    fESEdet; 

  Bool_t   fEP; 
  
  TProfile *fHistSP2A[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2C[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2T[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2ALowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2CLowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2TLowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2AHighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2CHighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP2THighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3A[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3C[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3T[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3ALowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3CLowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3TLowESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3AHighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3CHighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile *fHistSP3THighESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!

  TProfile2D *fHistSP2AvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile2D *fHistSP2CvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile2D *fHistSP2TvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile2D *fHistSP3AvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile2D *fHistSP3CvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!
  TProfile2D *fHistSP3TvsESE[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!

  TProfile2D *fHistSP2AdPhidEta[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  TProfile2D *fHistSP2CdPhidEta[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  TProfile2D *fHistSP2TdPhidEta[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  //TH3D       *fHistdPhidEta[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  AliTHn     *fHistdPhidEtaPt;
  TProfile2D *fHistSP2AdPhidEtaSS[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  TProfile2D *fHistSP2CdPhidEtaSS[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  TProfile2D *fHistSP2TdPhidEtaSS[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  //TH3D       *fHistdPhidEtaSS[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  AliTHn     *fHistdPhidEtaPt_SS; //!

  TProfile2D *fHistSP3AdPhidEtaSS[fNMaxBinsCentrality][fNMaxBinsZvtx][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  
  //TH3D       *fHistdPhidEtaMixed[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  //TH3D       *fHistdPhidEtaSSMixed[fNMaxBinsCentrality][fNMaxBinsPt][fNMaxBinsAssocPt]; //!
  AliTHn     *fHistdPhidEtaPt_Mixed; //!
  AliTHn     *fHistdPhidEtaPt_Mixed_SS; //!
  AliTHn     *fHistTrig; //!

  TProfile *fHistACv2; //!
  TProfile *fHistATv2; //!
  TProfile *fHistCTv2; //!
  TProfile *fHistACv3; //!
  TProfile *fHistATv3; //!
  TProfile *fHistCTv3; //!
  TProfile *fHistACv2LowESE; //!
  TProfile *fHistATv2LowESE; //!
  TProfile *fHistCTv2LowESE; //!
  TProfile *fHistACv2HighESE; //!
  TProfile *fHistATv2HighESE; //!
  TProfile *fHistCTv2HighESE; //!
  TProfile *fHistACv3LowESE; //!
  TProfile *fHistATv3LowESE; //!
  TProfile *fHistCTv3LowESE; //!
  TProfile *fHistACv3HighESE; //!
  TProfile *fHistATv3HighESE; //!
  TProfile *fHistCTv3HighESE; //!
  TProfile2D *fHistACv2vsESE; //!
  TProfile2D *fHistATv2vsESE; //!
  TProfile2D *fHistCTv2vsESE; //!
  TProfile2D *fHistACv3vsESE; //!
  TProfile2D *fHistATv3vsESE; //!
  TProfile2D *fHistCTv3vsESE; //!

  TProfile *fHistQxT2; //!
  TProfile *fHistQyT2; //!
  TProfile *fHistQxT3; //!
  TProfile *fHistQyT3; //!

  TH2D *fHistv2ESEvsCent; //!
  TH2D *fHistv3ESEvsCent; //!
  
  TH2D *fHistV0Multiplicity; //!
  TH1D *fHistITSMultiplicity; //!
  TH1D *fHistCentrality; //!
  TH1D *fHistEvStat; //!
  TH2D *fHistNtrVsCent; //!
  TH2D *fHistV0CVsCent; //!
  TH2D *fHistESDvsTPC; //!
  TH2D *fPileup1_Before; //!
  TH2D *fPileup1_After; //!
  TH2D *fPileup2_Before; //!
  TH2D *fPileup2_After; //!
  TH2D *fPileup3_Before_Low; //!
  TH2D *fPileup3_After_Low; //!
  TH2D *fPileup3_Before_High; //!
  TH2D *fPileup3_After_High; //!


  TH2D *fHistCentVsZ; //!
  TH2D *fHistCentVsZMixed; //!
  
  TString fCentMethod; //!

  TList *fOutputList; //!
  TList *fOutputList1; //!

  Int_t   fRunN; //!
  Float_t fv0mult; //!
  Float_t fv0multonline; //!
  Double_t fitsmult; //!
  Int_t    nITSTrkls; //!
  Double_t fv0mpercentile; //!
  Double_t fv0meqpercentile; //!
  Double_t fv0apercentile; //!
  Double_t fv0cpercentile; //!
  Double_t fcl0percentile; //!
  Double_t fcl1percentile; //!
  Double_t fzvtx; //!

  Int_t multtrkcut; //!
  Double_t Qx2trkcut; //!
  Double_t Qy2trkcut; //!
  Double_t Qx3trkcut; //!
  Double_t Qy3trkcut; //!

  Double_t sumMa; 
  Double_t sumMc; 
  Double_t Qxa2; 
  Double_t Qya2; 
  Double_t Qxc2; 
  Double_t Qyc2; 
  Double_t Qxa3; 
  Double_t Qya3; 
  Double_t Qxc3; 
  Double_t Qyc3; 
   
  Double_t Qya2Cor; 
  Double_t Qxa2Cor; 
  Double_t Qya3Cor; 
  Double_t Qxa3Cor; 
  Double_t Qyc2Cor; 
  Double_t Qxc2Cor; 
  Double_t Qyc3Cor; 
  Double_t Qxc3Cor; 
  
  Double_t Qytr2Cor; 
  Double_t Qxtr2Cor; 
  Double_t Qytr3Cor; 
  Double_t Qxtr3Cor; 
  /*
  Double_t Qytr2CorESE; //!
  Double_t Qxtr2CorESE; //!
  Double_t Qytr3CorESE; //!
  Double_t Qxtr3CorESE; //!
  Double_t q2Tr; //!
  Double_t q3Tr; //!
  */
  Double_t percq2; //!
  Double_t percq3; //!
  //  Double_t percq2Tr; //!
  //  Double_t percq3Tr; //!
  
  TF1*         fResACv2; //!
  TF1*         fResATv2; //!
  TF1*         fResCTv2; //!
  TF1*         fResACv2LowESE; //!
  TF1*         fResATv2LowESE; //!
  TF1*         fResCTv2LowESE; //!
  TF1*         fResACv2HighESE; //!
  TF1*         fResATv2HighESE; //!
  TF1*         fResCTv2HighESE; //!
  TF1*         fResACv3; //!
  TF1*         fResATv3; //!
  TF1*         fResCTv3; //!
  TF1*         fResACv3LowESE; //!
  TF1*         fResATv3LowESE; //!
  TF1*         fResCTv3LowESE; //!
  TF1*         fResACv3HighESE; //!
  TF1*         fResATv3HighESE; //!
  TF1*         fResCTv3HighESE; //!
  TF1*         fResACv2vsESE[5]; //!
  TF1*         fResATv2vsESE[5]; //!
  TF1*         fResCTv2vsESE[5]; //!
  TF1*         fResACv3vsESE[5]; //!
  TF1*         fResATv3vsESE[5]; //!
  TF1*         fResCTv3vsESE[5]; //!
  
  TF1*         fLowCenCut;          //! cut low for centrality outliers
  TF1*         fHighCenCut;         //! cut high for centrality outliers
  TF1*         fV0MultOfOnCut;      //! cut V0M offline-online outliers

  TH1D*        fMultV0;             //! profile from V0 multiplicity

  TH1D*        fQx2mV0A[14];            //! <Qxn> V0A
  TH1D*        fQy2mV0A[14];            //! <Qyn> V0A
  TH1D*        fQx2sV0A[14];            //! sigma Qxn V0A
  TH1D*        fQy2sV0A[14];            //! sigma Qyn V0A
    
  TH1D*        fQx2mV0C[14];            //! <Qxn> V0C
  TH1D*        fQy2mV0C[14];            //! <Qyn> V0C
  TH1D*        fQx2sV0C[14];            //! sigma Qxn V0C
  TH1D*        fQy2sV0C[14];            //! sigma Qyn V0C
  
  TH1D*        fQx2mTrk[14];            //! <Qxn> tracklets
  TH1D*        fQy2mTrk[14];            //! <Qyn> tracklets
  TH1D*        fQx2sTrk[14];            //! sigma Qxn tracklets
  TH1D*        fQy2sTrk[14];            //! sigma Qyn tracklets
    
    
  TH1D*        fQx3mV0A[14];            //! <Qxn> V0A
  TH1D*        fQy3mV0A[14];            //! <Qyn> V0A
  TH1D*        fQx3sV0A[14];            //! sigma Qxn V0A
  TH1D*        fQy3sV0A[14];            //! sigma Qyn V0A
    
  TH1D*        fQx3mV0C[14];            //! <Qxn> V0C
  TH1D*        fQy3mV0C[14];            //! <Qyn> V0C
  TH1D*        fQx3sV0C[14];            //! sigma Qxn V0C
  TH1D*        fQy3sV0C[14];            //! sigma Qyn V0C
  
  TH1D*        fQx3mTrk[14];            //! <Qxn> tracklets
  TH1D*        fQy3mTrk[14];            //! <Qyn> tracklets
  TH1D*        fQx3sTrk[14];            //! sigma Qxn tracklets
  TH1D*        fQy3sTrk[14];            //! sigma Qyn tracklets

  
  AliAnalysisTaskSEPbPbCorrelationsJetV2(const AliAnalysisTaskSEPbPbCorrelationsJetV2&);//not implimented
  AliAnalysisTaskSEPbPbCorrelationsJetV2& operator=(const AliAnalysisTaskSEPbPbCorrelationsJetV2&);//not implimnted
  
  ClassDef(AliAnalysisTaskSEPbPbCorrelationsJetV2, 1)  // example of analysis

};

//====================================================================================================================================================

class AliBasicParticleST : public AliVParticle
{
public:
 AliBasicParticleST() : fEta(0), fPhi(0), fpT(0), fCharge(0) {}
 AliBasicParticleST(Float_t eta, Float_t phi, Float_t pt, Short_t charge) : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge) {}
  ~AliBasicParticleST() {}

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

  virtual Short_t Charge()      const { return fCharge; }
  virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
  // PID
  virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
  virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }

  virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }

  virtual void SetPhi(Double_t phi) { fPhi = phi; }

private:
  Float_t fEta;      //! eta
  Float_t fPhi;      //! phi
  Float_t fpT;       //! pT
  Short_t fCharge;   //! charge

  ClassDef(AliBasicParticleST,1) // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};

#endif
