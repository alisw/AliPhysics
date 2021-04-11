#ifndef ALIANALYSISTASKSEPBPBCORRELATIONSJETV2_H
#define ALIANALYSISTASKSEPBPBCORRELATIONSJETV2_H

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
  void SetPtOrder(Bool_t flag)       {fPtOrder = flag; }
  void SetOnlySameSign(Bool_t flag)  {fSameSign = flag; }
  void SetUseRes(Bool_t flag) { fUseRes = flag; }
  void SetAverageRes(Bool_t flag) {fAverageRes = flag; }
  void SetSubDivide(Int_t n1, Int_t n2) { fN1 = n1; fN2 = n2; }
  void SetForceCL1(Bool_t flag) {fForceCL1 = flag;} 
 
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
  void SetListContQ(TList *list) { flist_contQ = list; }
  void SetListRes(TList *list)   { flist_Res = list; }
  void SetFilterBit(Int_t iBit) {fFilterBit = iBit; }
  void SetNTPCcls(Int_t nTPC) {fTPCNcls = nTPC; }
  void SetRefMode(TString s) {fMode = s;} 

 
  void FillHistogramsV2(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
			Double_t resA2, Double_t resC2, Double_t resT2,
			Int_t index);

  void FillHistogramsdPhidEta(TObjArray *selectedArray, Int_t centrality,Double_t percentile,Int_t zvtxBin,
                                                   Double_t resA2, Double_t resC2, Double_t resT2,
                                                   Double_t resA3, Double_t resC3, Double_t resT3);
  void FillHistogramsdPhidEtaMixed(TObjArray *selectedArray, Double_t percentile,Int_t zvtxBin);

  TObjArray* CloneTrack(TObjArray*selectedTrackArray);

  TObjArray *GetAcceptedTracks(AliAODEvent *fAOD, TObjArray*tracks);
  Int_t GetCentBin();
  Int_t GetCentBin(Double_t percentile);
  Int_t GetZvtxBin(Double_t zvtx);

  Double_t CalcCorrectedPhi(Double_t phi, Double_t dPhi) const;

  void ProcessEvent(Double_t percentile, Int_t centBin, Int_t zvtxBin);

  void OpenInfoCalbration(Int_t run);
  Short_t GetVertexZ(Double_t vtxZ) const;
  Bool_t ComputeQ(AliAODEvent* aod, Double_t Zvtx);

  void CalcResolutions(Double_t percentile,
		       Double_t &resA2, Double_t &resC2, Double_t &resT2);
  
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
  Bool_t fPtOrder;
  Bool_t fUseRes; 
  Bool_t fAverageRes; 
  Bool_t fSameSign;
  Bool_t fForceCL1; 
 
  Int_t fN1; 
  Int_t fN2; 


  
  TProfile *fHistSP2A[fNMaxBinsCentrality][fNMaxBinsZvtx][3]; //!

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
  
  TProfile *fHistQxT2; //!
  TProfile *fHistQyT2; //!
  TProfile *fHistQxT3; //!
  TProfile *fHistQyT3; //!
 
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
  TString fMode;       //!

  TList *fOutputList;  //!
  TList *fOutputList1; //!
  TList *flist_contQ;  //! list_contQ
  TList *flist_Res;    //! list_Res

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
  Double_t fSPDpercentile; //!
  Double_t fzvtx; //!

  Int_t multtrkcut; //!
  Double_t Qx2trkcut; //!
  Double_t Qy2trkcut; //!

  Double_t sumMa; 
  Double_t sumMc; 
  Double_t Qxa2; 
  Double_t Qya2; 
  Double_t Qxc2; 
  Double_t Qyc2; 
   
  Double_t Qya2Cor; 
  Double_t Qxa2Cor; 
  Double_t Qyc2Cor; 
  Double_t Qxc2Cor; 
  
  Double_t Qytr2Cor; 
  Double_t Qxtr2Cor; 

  Int_t fFilterBit;
  Int_t fTPCNcls;
  
  TF1*         fResACv2; //!
  TF1*         fResATv2; //!
  TF1*         fResCTv2; //!
  
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
    
  AliAnalysisTaskSEPbPbCorrelationsJetV2(const AliAnalysisTaskSEPbPbCorrelationsJetV2&);//not implimented
  AliAnalysisTaskSEPbPbCorrelationsJetV2& operator=(const AliAnalysisTaskSEPbPbCorrelationsJetV2&);//not implimnted
  
  ClassDef(AliAnalysisTaskSEPbPbCorrelationsJetV2, 1)  // example of analysis

};

//====================================================================================================================================================
class AliBasicParticleST : public AliVParticle {
public:
  AliBasicParticleST(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
                       Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate,
                       Double_t multiplicity)
      : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fID(ID), fID1(ID1),
        fID2(ID2), fCandidate(candidate), fMultiplicity(multiplicity) {}
  virtual ~AliBasicParticleST() {}
  
  virtual Double_t Px() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Py() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Pz() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Pt() const { return fpT; }
  virtual Double_t P() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t PxPyPz(Double_t[3]) const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Xv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Yv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Zv() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Bool_t XvYvZv(Double_t[3]) const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t OneOverPt() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Phi() const { return fPhi; }
  virtual Double_t Theta() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t E() const {
    AliFatal("Not implemented"); 
    return 0;
  } 
  virtual Double_t M() const {
    AliFatal("Not implemented");
    return 0;
  } 
  virtual Double_t Eta() const { return fEta; }
  virtual Double_t Y() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Short_t Charge() const { return fCharge; }
  virtual Int_t GetLabel() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Int_t PdgCode() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual const Double_t *PID() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Short_t WhichCandidate() const { return fCandidate; }
  virtual Int_t GetID() const { return fID; }
  virtual Int_t GetIDFirstDaughter() const { return fID1; }
  virtual Int_t GetIDSecondDaughter() const { return fID2; }
  virtual Double_t Multiplicity() const { return fMultiplicity; }

  private:
  //
  Short_t fCharge;    // Charge
  Float_t fEta;       // Eta
  Float_t fPhi;       // Phi
  Float_t fpT;        // pT
  Int_t fID;          // ID
  Short_t fCandidate; // 1-pi,2-K,3-p
  Double_t fMultiplicity;
  Int_t fID1;
  Int_t fID2;
  ClassDef(AliBasicParticleST, 1);
};

#endif
