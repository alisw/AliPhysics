/*
 *****************************************************************************************/
#ifndef ALIANALYSISTASKSEPPBCORRELATIONSMCLEGOYS
#define ALIANALYSISTASKSEPPBCORRELATIONSMCLEGOYS

#include "AliAnalysisTask.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliVAD.h"
#include "AliVParticle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TString.h"
#include "AliEventCuts.h"

class TList;
class AliCFContainer;
// class AliTHn;
class AliAODEvent;
class AliEventPoolManager;
class AliVParticle;
class AliPIDResponse;
// class AliPID;
class AliAODv0;
// class AliAnalyseLeadingTrackUE;
class THnSparse;
class AliAODcascade;
class AliAODVertex;


//#ifndef ALIANALYSISTASKSEH
#include "AliAnalysisTaskSE.h"
//#endif

//---------------------------------------------------------------------------------------
class AliAnalysisTaskSEpPbCorrelationsMCLEGOYS : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEpPbCorrelationsMCLEGOYS();
  AliAnalysisTaskSEpPbCorrelationsMCLEGOYS(const char *name);
  virtual ~AliAnalysisTaskSEpPbCorrelationsMCLEGOYS();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // Analysis Correlation
  void MakeAna();

  // configration
  virtual void SetAnalysisMode(TString mode) { fAnaMode = mode; }
  virtual void SetAssociatedTrack(TString mode) { fasso = mode; }
  virtual void SetPID(Bool_t mode) { fPID = mode; }
  virtual void SetDatatype(Bool_t mode) { fDataType = mode; }
  virtual void SetCentCalib(Bool_t mode) { fcentcalib= mode; }
  virtual void SetRunType(Bool_t mode) { frun2 = mode; }
  virtual void SetFilterBit(Int_t mode) { ffilterbit = mode; }
  virtual void SetFMDcut(Bool_t mode) {fFMDcut=mode;}
  virtual void SetFMDcutpar(Int_t mode){fFMDcutmode=mode;}
  virtual void SetPtdiff(Bool_t mode){fptdiff=mode;}
  virtual void SetExtractSec(Bool_t mode){fextractsec=mode;}
  virtual void SetPtMax(Float_t mode){fPtMax=mode;}
  virtual void SetPtMin(Float_t mode){fPtMin=mode;}
  
  
  virtual void Setacceptancehole(Bool_t mode){fmakehole=mode;}
  virtual void SetAnalysisCent(TString mode) { fCentType = mode; }
  virtual void SetAnalysisCollisionType(TString mode) { fcollisiontype = mode; }
  virtual void SetMCclosure(Bool_t mode) { fMCclosure = mode; }
  virtual void SetFillCorrelation(Bool_t mode) { ffillcorrelation = mode; }
  virtual void Setmcprim(Bool_t mode){fprim=mode;}
  virtual void SetQAmode(Bool_t mode){fQA=mode;}
  
  void SetMaxNEventsInPool(Int_t events) { fPoolMaxNEvents = events; }
  void SetMinNTracksInPool(Int_t tracks) { fPoolMinNTracks = tracks; }
  void SetMinEventsToMix(Int_t events) { fMinEventsToMix = events; }

  void SetPoolPVzBinLimits(Int_t Nzvtxbins, const Double_t *ZvtxBins) {
    fNzVtxBins = Nzvtxbins;
    for (int ix = 0; ix < fNzVtxBins + 1; ix++) {
      fZvtxBins[ix] = ZvtxBins[ix];
    }
  }

  void SetPoolCentBinLimits(Int_t Ncentbins, const Double_t *CentBins) {
    fNCentBins = Ncentbins;
    for (int ix = 0; ix < fNCentBins + 1; ix++) {
      fCentBins[ix] = CentBins[ix];
    }
  }
  void DumpTObjTable(const char* note);

  
private:
  AliAnalysisTaskSEpPbCorrelationsMCLEGOYS(
      const AliAnalysisTaskSEpPbCorrelationsMCLEGOYS &det);
  AliAnalysisTaskSEpPbCorrelationsMCLEGOYS &
  operator=(const AliAnalysisTaskSEpPbCorrelationsMCLEGOYS &det);

  void DefineGeneralOutput();
  void DefineVZEROOutput();
  void DefineCorrOutput();
  void DefinedQAHistos();

  TObjArray *GetAcceptedTracksLeading(AliVEvent *faod,Bool_t leading,TObjArray*tracks);
  TObjArray *GetAcceptedTracksPID(AliVEvent *faod);
  TObjArray *GetAcceptedV0Tracks(const AliVEvent *faod);
  TObjArray *GetAcceptedCascadeTracks(AliVEvent *faod);
  TObjArray *GetAcceptedTracksAssociated(AliVEvent *faod);

  void  CalculateSP();
  TObjArray* GetFMDhitsYS(Bool_t Aside);
  Bool_t IsAcceptedDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedPhiDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedTrack(const AliAODTrack *aodTrack);
  Bool_t IsAcceptedDecayLength(const AliAODv0 *aodv0,Double_t mass,Double_t maxctau);
  Bool_t IsAcceptedV0(const AliAODv0 *aodv0);
  Bool_t IsAcceptedCascade(const AliAODcascade *casc);
  Bool_t IsAcceptedCascadeOmega(const AliAODcascade *casc);

  TObjArray* CloneTrack(TObjArray* track);
  Double_t RangePhi(Double_t DPhi);
  Double_t RangePhi_FMD(Double_t DPhi);
  Double_t RangePhi2(Double_t DPhi);
 Int_t      ConvertRunNumber(Int_t run);

/*
  void FillCorrelationTracksCentralForward(Double_t MultipOrCent, TObjArray *triggerArray,
                             TObjArray *selectedTrackArray, AliTHn *triggerHist,
                             AliTHn *associateHist, Bool_t, Float_t, Float_t,
                             Float_t, Int_t);

*/
  void FillCorrelationTracks(Double_t MultipOrCent, TObjArray *triggerArray,
                             TObjArray *selectedTrackArray, AliTHn *triggerHist, AliTHn *associateHist, Bool_t twoTrackEfficiencyCut, Float_t twoTrackEfficiencyCutValue, Float_t fTwoTrackCutMinRadius,Float_t bSign, Int_t step);
  void FillCorrelationTracksMixing(Double_t MultipOrCentMix, Double_t pvxMix,
                                   Double_t poolmax, Double_t poolmin,
                                   TObjArray *triggerArray,
                                   TObjArray *selectedV0Array,
                                   AliTHn *triggerHist, AliTHn *associateHist,
                                   Bool_t, Float_t, Float_t, Float_t, Int_t);

  inline Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1,
                             Float_t phi2, Float_t pt2, Float_t charge2,
                             Float_t radius, Float_t bSign);

  TString fcollisiontype;
  Bool_t fDataType;
  Bool_t fcentcalib;
  Bool_t frun2;
  Bool_t fQA;
  Bool_t fMCclosure;
  Bool_t fFMDcut;
  Int_t fFMDcutmode;
  Bool_t fptdiff;
  Bool_t fextractsec;
  Bool_t fmakehole;
  Bool_t fOnfly;
  TString fAnaMode;
  TString fasso;
  Bool_t fPID;

  TString fCentType;
  Bool_t ffillcorrelation;
  Bool_t fprim;
  Int_t fNEntries;
  
  Double_t lCentrality;
  Float_t bSign;
  Double_t fZVertex;

  TList *fOutputList;  // Output list
  TList *fOutputList1; // Output list
  TList *fOutputList2; // Output list

  AliPIDResponse *fPIDResponse; // PID Response
  TH2D* fhcorreffi[10];

  Int_t ffilterbit;
  Double_t fPtMin;
  Double_t fPtMax;
  Double_t fEtaMax;
  Double_t fEtaMaxExtra;
  Double_t fEtaMinExtra;
  // V0 particles
  Double_t fEtaMaxV0;
  Double_t fEtaMinV0;
  Double_t fdcaDaughtersToPrimVtx;
  Double_t fdcaBetweenDaughters;
  Double_t fRadiMin;
  Double_t fRadiMax;
  Double_t fcutcTauLam;
  Double_t fcutcTauK0;
  Double_t fcosMinK0s;
  Double_t fcosMinLambda;
  Double_t fMaxnSigmaTPCV0;
  // V0 daughter cut
  TH1D *hv0dcharge;
  Double_t fclustermin;
  Double_t fratiocluster;
  Double_t fEtaMaxDaughter;
  Double_t fEtaMinDaughter;
  THnSparseF *fHistMass_K0s;
  THnSparseF *fHistMass_K0s_MC;
  THnSparseF *fHistMass_Lambda;
  THnSparseF *fHistMass_ALambda;
  THnSparseF *fHistMass_ALambda_MC;
  THnSparseF *fHistMassXiMinus;
  THnSparseF *fHistMassXiPlus;
  THnSparseF *fHistMassOmegaMinus;
  THnSparseF *fHistMassOmegaPlus;
  TH2D *fHistMass_bumpcorr;
  THnSparseF *fHist_V0QA;
  THnSparseF *fHist_CascadeQA;
  TH2D *fHist_AP[6];
  TH2D *fHistPosNsig[6];
  TH2D *fHistNsig[6];
  TH2D *fHistNsigcorr[6];
  TH2D *fHistNegNsig[6];
  TH2D *fHistPosNsigQA[6];
  TH3D* fh3NegNsig[3];
  TH3D* fh3PosNsig[3];

  THnSparseF *fHistMass_Lambda_MC;

  //	Double_t fPtMinDaughter

  AliEventCuts fEventCuts; 
  AliAnalysisUtils* fUtils;
  AliVEvent *fEvent; //  AOD Event
  AliMCEvent* mcEvent;
  const  AliVVertex *lPrimaryBestVtx;
  Double_t tPrimaryVtxPosition[3];
  Double_t fPrimaryZVtx;

  AliVVZERO *fvzero;

  // Event Pool for mixing
  AliEventPoolManager *fPoolMgr;  //  event pool manager for Event Mixing
  AliEventPoolManager *fPoolMgr1; //  event pool manager for Event Mixing
  Double_t poolmin;
  Double_t poolmax;
  Int_t fPoolMaxNEvents;   // set maximum number of events in the pool
  Int_t fPoolMinNTracks;   // set minimum number of tracks in the pool
  Int_t fMinEventsToMix;   // set the minimum number of events want to mix
  Int_t fNzVtxBins;        // number of z vrtx bins
  Double_t fZvtxBins[100]; // [fNzVtxBinsDim]
  Int_t fNCentBins;        // number of centrality bins
  Double_t fCentBins[100]; // [fNCentBinsDim]
  // Track cuts
  Double_t fMaxnSigmaTPCTOF;

  // Global Histograms
  TH1F *fHistzvertex;
  TH1F *fHistCentrality;
  TH1F *fHistCentrality_beforecut;
  TH1F* fHistImpactparameter;
  TH2F* fHistCentzvertex;
  TH2F* mixedDist;
  TH2F* mixedDist2;
  
  
  AliTHn *fHistLeadQA;
  AliTHn *fHistPIDQA;

  AliTHn* fhistmcprim;
  AliTHn* fhistmcprimfinal;
  TH2D* fNTrackCorrMC;
  TH2D*fhmcprimvzeta;
  TH2D*   fhmcrapicent;
  TH2D* fhistmeanpt;
  

  TH1F*frefvz;
  TH1D*fhcorr[10];

  TH1D*fhmcprimpdgcode;
  TH1D*fhrefetaFMD[4];
  TH1D*fhrefphiFMD[4];

  TH2D*  fh2_FMD_acceptance_prim;
  TH2D*  fh2_FMD_eta_phi_prim;
  TH2D*  fh2_FMD_acceptance;
  TH2D*  fh2_ITS_acceptance;
  TH2F*  fh2_SPD_multcorr;
  TH2F*  fh2_SPDV0_multcorr;
  TH2F*  fh2_SPDtrack_multcorr;
  TH1F*  fhtrackletsdphi;
  TH2D*  fh2_FMD_eta_phi;
  TH1F* fHist_NeventRun;
  TH1F* fHist_V0AMultRun;
  TH1F* fHist_V0CMultRun;
  TH1F* fHist_FMDAMultRun;
  TH1F* fHist_FMDCMultRun;

  TH2D*  fhistfmdphiacc;
  AliTHn* fhistfmd;
  THnSparseF* fhistits;
  AliTHn* fhSecFMD;
  //  const TH2D& d2Ndetadphi;
  TH2F*fFMDV0;
  TH2F*fFMDV0_post;
  TH2F*fFMDV0A;
  TH2F*fFMDV0A_post;
  TH2F*fFMDV0C;
  TH2F*fFMDV0C_post;

  TH1F*fV0Amultprim;
  TH1F*fV0Amultmodi;
  TH2F*fh2_V0A;
  TH2F*fh2_V0A_all;
  TH2F*fh2_V0C;
  TH2F*fh2_V0A_comp;
  TH2F*fh2_V0A_comp_prim;
  
  TH2F *fHist_vzeromult;
  TH2F *fHist_vzeromultEqweighted;
  TH3F *fHist2dmult;
  AliTHn *fHistVZERO;

  TH1F *fHist_Stat;
  TH1F *fHist_V0Stat;

  // QA histograms
  TH2D *fHistPhiDTPCNSig;
  TH2D *fHistPhiDTOFNSig;
  TH2D *fHistPhiDTPCTOFNSig;
  THnSparseF *fHistMass_PhiMeson;
  THnSparseF *fHistMass_PhiMeson_MIX;
  THnSparseF *fHist_PhiQA;

  AliTHn *fHistTriggerTrack;
  AliTHn *fHistReconstTrack;
  AliTHn *fHistTriggerTrackMix;
  AliTHn *fHistReconstTrackMix;

  TH2D* fHistQna;
  TH2D* fHistQnc;
  TH2D* fHistCorrQna[4];
  TH2D* fHistCorrQnc[4];
  TH2D* fHistQn;
  TH2D* fHistQna_VZERO;
  TH2D* fHistQnc_VZERO;
  TH2D* fHistQn_VZERO;
  TH1D* fHistVn;
  TH1D* fHistQAQB[4];
  TH1D* fHistQAQB_VZERO[4];
  TProfile* SP_TPCATPCC;
  TProfile* SP_TPCATPCC_default;
  TProfile* SP_V0AV0C_default;
  TProfile* SP_V0ATPC_default;
  TProfile* SP_V0CTPC_default;
  TH1F* fHist_V0AV0C;
  TH1F* fHist_V0ATPC;
  TH1F* fHist_V0CTPC;
  TProfile* SP_uTPCA;
  TProfile* SP_uTPCC;
  TProfile* SP_uTPC_PP[8];
  TProfile* SP_uTPC[8];
  TProfile* SP_uTPC1[8];
  TProfile* SP_uTPC2[8];
  TProfile* SP_uTPC3[8];
  TProfile* SP_uVZEROA_PP[8];
  TProfile* SP_uVZEROA[8];
  TProfile* SP_uVZEROA1[8];
  TProfile* SP_uVZEROA2[8];
  TProfile* SP_uVZEROA3[8];
  TProfile* SP_uVZEROC_PP[8];
  TProfile* SP_uVZEROC[8];
  TProfile* SP_uVZEROC1[8];
  TProfile* SP_uVZEROC2[8];
  TProfile* SP_uVZEROC3[8];

  ClassDef(AliAnalysisTaskSEpPbCorrelationsMCLEGOYS, 2);
};
//---------------------------------------------------------------------------------------
Float_t AliAnalysisTaskSEpPbCorrelationsMCLEGOYS::GetDPhiStar(
    Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2,
    Float_t charge2, Float_t radius, Float_t bSign) {
  //
  // calculates dphistar (Copied from AliUEHistogram)
  //

  Float_t dphistar = phi1 - phi2 -
                     charge1 * bSign * TMath::ASin(0.075 * radius / pt1) +
                     charge2 * bSign * TMath::ASin(0.075 * radius / pt2);

  static const Double_t kPi = TMath::Pi();

  // circularity
  //   if (dphistar > 2 * kPi)
  //     dphistar -= 2 * kPi;
  //   if (dphistar < -2 * kPi)
  //     dphistar += 2 * kPi;

  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi * 2 - dphistar;
  if (dphistar > kPi) // might look funny but is needed
    dphistar = kPi * 2 - dphistar;

  return dphistar;
}

class AliAssociatedTrackYSLEGOMC : public AliVParticle {
public:
  AliAssociatedTrackYSLEGOMC(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
                       Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate,
                       Double_t multiplicity)
      : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fID(ID), fID1(ID1),
        fID2(ID2), fCandidate(candidate), fMultiplicity(multiplicity) {}
  virtual ~AliAssociatedTrackYSLEGOMC() {}

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
  Short_t fCharge;    // Charge
  Float_t fEta;       // Eta
  Float_t fPhi;       // Phi
  Float_t fpT;        // pT
  Int_t fID;          // ID
  Short_t fCandidate; // 1-pi,2-K,3-p
  Double_t fMultiplicity;
  Int_t fID1;
  Int_t fID2;
  ClassDef(AliAssociatedTrackYSLEGOMC, 1);
};

class AliMixTrackYSLEGOMC : public AliVParticle {
public:
  AliMixTrackYSLEGOMC(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
                Float_t px, Float_t py, Float_t pz

                )
      : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fpx(px), fpy(py),
        fpz(pz) {}
  virtual ~AliMixTrackYSLEGOMC() {}

  virtual Double_t Px() const { return fpx; }
  virtual Double_t Py() const { return fpy; }
  virtual Double_t Pz() const { return fpz; }
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
  virtual Short_t WhichCandidate() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Int_t GetID() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Int_t GetIDFirstDaughter() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Int_t GetIDSecondDaughter() const {
    AliFatal("Not implemented");
    return 0;
  }
  virtual Double_t Multiplicity() const {
    AliFatal("Not implemented");
    return 0;
  }

private:
  Short_t fCharge; // Charge
  Float_t fEta;    // Eta
  Float_t fPhi;    // Phi
  Float_t fpT;     // pT
  Float_t fpx;
  Float_t fpy;
  Float_t fpz;

  ClassDef(AliMixTrackYSLEGOMC, 1);
};
//---------------------------------------------------------------------------------------

class AliAssociatedVZEROYSLEGOMC : public AliVParticle {
public:
  AliAssociatedVZEROYSLEGOMC(Float_t multiplicity,
                       //		  Double_t multiplicity,
                       Float_t eta,
                       // Float_t phi,
                       Double_t phi, Float_t pt, Int_t ID, Short_t candidate)
      :

        fMultiplicity(multiplicity),
        fEta(eta), fPhi(phi), fpT(pt), fID(ID), fCandidate(candidate) {}
  virtual ~AliAssociatedVZEROYSLEGOMC() {}

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
  virtual Double_t Multiplicity() const { return  fMultiplicity; }
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
  virtual Short_t Charge() const {
    AliFatal("Not implemented");
    return 0;
  }

private:
  Float_t fMultiplicity; // Charge
  // Double_t fMultiplicity;
  Float_t fEta; // Eta
  // Float_t  fPhi;            // Phi
  Double_t fPhi;      // Phi
  Float_t fpT;        // pT
  Int_t fID;          // ID
  Short_t fCandidate; // 1-pi,2-K,3-p

  ClassDef(AliAssociatedVZEROYSLEGOMC, 1);
};

#endif
