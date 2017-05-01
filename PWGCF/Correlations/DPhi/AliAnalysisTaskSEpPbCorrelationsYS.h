/**************************************************************************************************                                                                                                               
 *      Leading Charged Track+V0 Correlations.(Works for Real  Data)  *                                                                                                                                           
 *                 Yuko Sekiguchi * Center for Nuclear Study(CNS) , University of Tokyo                              *                                                                                            
 *                    Email:y_sekiguchi@cns.s.u-tokyo.ac.jp *                                                                                                                                                     
 **************************************************************************************************/

#ifndef ALIANALYSISTASKSEPPBCORRELATIONS
#define ALIANALYSISTASKSEPPBCORRELATIONS

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

#ifndef ALIANALYSISTASKSEH
#include "AliAnalysisTaskSE.h"
#endif

//---------------------------------------------------------------------------------------
class AliAnalysisTaskSEpPbCorrelationsYS : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskSEpPbCorrelationsYS();
  AliAnalysisTaskSEpPbCorrelationsYS(const char *name);
  virtual ~AliAnalysisTaskSEpPbCorrelationsYS();

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
  virtual void SetRunType(Bool_t mode) { frun2 = mode; }
  virtual void SetFilterBit(Int_t mode) { ffilterbit = mode; }

  virtual void SetAnalysisCent(TString mode) { fCentType = mode; }
  virtual void SetAnalysisCollisionType(TString mode) { fcollisiontype = mode; }

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

private:
  AliAnalysisTaskSEpPbCorrelationsYS(
      const AliAnalysisTaskSEpPbCorrelationsYS &det);
  AliAnalysisTaskSEpPbCorrelationsYS &
  operator=(const AliAnalysisTaskSEpPbCorrelationsYS &det);

  void DefineGeneralOutput();
  void DefineVZEROOutput();
  void DefineCorrOutput();
  void DefinedQAHistos();

  TObjArray *GetAcceptedTracksLeading(AliAODEvent *faod,Bool_t leading);
  TObjArray *GetAcceptedTracksPID(AliAODEvent *faod);
  TObjArray *GetAcceptedV0Tracks(const AliAODEvent *faod);
  TObjArray *GetAcceptedCascadeTracks(const AliAODEvent *faod);
  TObjArray *GetAcceptedTracksAssociated(AliAODEvent *faod);

  Bool_t IsAcceptedDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedPhiDaughterTrack(const AliAODTrack *itrack);
  Bool_t IsAcceptedTrack(const AliAODTrack *aodTrack);
  Bool_t IsAcceptedDecayLength(const AliAODv0 *aodv0,Double_t mass,Double_t maxctau);
  Bool_t IsAcceptedV0(const AliAODv0 *aodv0);
  Bool_t IsAcceptedCascade(const AliAODcascade *casc);
  Bool_t IsAcceptedCascadeOmega(const AliAODcascade *casc);

  Double_t RangePhi(Double_t DPhi);
  Double_t RangePhi2(Double_t DPhi);

  void FillCorrelationTracks(Double_t MultipOrCent, TObjArray *triggerArray,
                             TObjArray *selectedTrackArray, AliTHn *triggerHist,
                             AliTHn *associateHist, Bool_t, Float_t, Float_t,
                             Float_t, Int_t);
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
  Bool_t frun2;
  Bool_t fQA;
  Bool_t fOnfly;
  TString fAnaMode;
  TString fasso;
  Bool_t fPID;

  TString fCentType;

  Double_t lCentrality;
  Float_t bSign;
  Double_t fZVertex;

  TList *fOutputList;  // Output list
  TList *fOutputList1; // Output list
  TList *fOutputList2; // Output list

  AliPIDResponse *fPIDResponse; // PID Response

  Int_t ffilterbit;
  Double_t fPtMin;
  Double_t fEtaMax;
  // V0 particles
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
  TH2D *fHistPosNsig[3];
  TH2D *fHistNsig[6];
  TH2D *fHistNsigcorr[6];
  TH2D *fHistNegNsig[3];
  TH2D *fHistPosNsigQA[6];


  THnSparseF *fHistMass_Lambda_MC;

  //	Double_t fPtMinDaughter;

  AliAODEvent *fEvent; //  AOD Event
  AliAODVertex *lPrimaryBestVtx;
  Double_t tPrimaryVtxPosition[3];
  Double_t fPrimaryZVtx;

  AliAODVZERO *fvzero;

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

  AliTHn *fHistLeadQA;
  AliTHn *fHistPIDQA;

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

  ClassDef(AliAnalysisTaskSEpPbCorrelationsYS, 2);
};
//---------------------------------------------------------------------------------------
Float_t AliAnalysisTaskSEpPbCorrelationsYS::GetDPhiStar(
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

class AliAssociatedTrackYS : public AliVParticle {
public:
  AliAssociatedTrackYS(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
                       Int_t ID, Int_t ID1, Int_t ID2, Short_t candidate,
                       Double_t multiplicity)
      : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fID(ID), fID1(ID1),
        fID2(ID2), fCandidate(candidate), fMultiplicity(multiplicity) {}
  virtual ~AliAssociatedTrackYS() {}

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
  ClassDef(AliAssociatedTrackYS, 1);
};

class AliMixTrackYS : public AliVParticle {
public:
  AliMixTrackYS(Short_t charge, Float_t eta, Float_t phi, Float_t pt,
                Float_t px, Float_t py, Float_t pz

                )
      : fCharge(charge), fEta(eta), fPhi(phi), fpT(pt), fpx(px), fpy(py),
        fpz(pz) {}
  virtual ~AliMixTrackYS() {}

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

  ClassDef(AliMixTrackYS, 1);
};
//---------------------------------------------------------------------------------------

class AliAssociatedVZEROYS : public AliVParticle {
public:
  AliAssociatedVZEROYS(Float_t multiplicity,
                       //		  Double_t multiplicity,
                       Float_t eta,
                       // Float_t phi,
                       Double_t phi, Float_t pt, Int_t ID, Short_t candidate)
      :

        fMultiplicity(multiplicity),
        fEta(eta), fPhi(phi), fpT(pt), fID(ID), fCandidate(candidate) {}
  virtual ~AliAssociatedVZEROYS() {}

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
  virtual Double_t Multiplicity() const { return fMultiplicity; }
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

  ClassDef(AliAssociatedVZEROYS, 1);
};

#endif
