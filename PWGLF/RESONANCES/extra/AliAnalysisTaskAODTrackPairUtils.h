#ifndef AliAnalysisTaskAODTrackPairUtils_cxx
#define AliAnalysisTaskAODTrackPairUtils_cxx

#include "AliAnalysisMuonUtility.h"
#include "AliMuonTrackCuts.h"
#include "AliVEventHandler.h"
#include "TFile.h"
#include "TNamed.h"
#include "TRandom1.h"
#include <iostream>

class AliAnalysisTaskAODTrackPairUtils : public TNamed {

public:
  AliAnalysisTaskAODTrackPairUtils();
  ~AliAnalysisTaskAODTrackPairUtils();

  void setPeriod(std::string period) { fPeriod = period; }

  bool setEvent(AliAODEvent *event, AliVEventHandler *handler);
  bool setMCEventInfo();
  void setMCArray(TClonesArray *array) { fMCArray = array; }
  void setMC(bool isMC) { fIsMC = isMC; }
  void setEvtSelection(bool isEvtSel) { fIsEvtSelect = isEvtSel; }

  bool isAcceptEvent();
  bool isAcceptTrackKinematics(AliAODTrack *track);

  bool isAcceptFwdMuonTrack(AliAODTrack *track);
  bool isAcceptFwdDimuon(AliAODDimuon *dimuon);

  bool isAcceptMidPrimTrackQuality(AliAODTrack *track);
  bool isAcceptMidPid(AliAODTrack *track, AliPID::EParticleType pid);

  bool isAcceptedK0sFromKpmStar(AliAODv0 *v0, double &mass);
  bool isAcceptV0TrackQuality(AliAODTrack *track);
  bool isAcceptV0Basic(AliAODv0 *v0, int charge);
  bool isAcceptV0Quality(AliAODv0 *v0, int charge);
  bool isAcceptV0Kinematics(AliAODv0 *v0);
  bool isAcceptK0s(AliAODv0 *v0);

  bool isAcceptArmenterosK0s(AliAODv0 *v0);
  bool isAcceptArmenterosK0s_Tight(AliAODv0 *v0);

  bool isAcceptTrackPairOpeningAngle(double angle) {
    if (fMaxPairCosOpeningAngleCut > angle) {
      return true;
    } else {
      return false;
    }
  }

  bool isAcceptK0sCandidateMassRange(double mass) {
    if (fMinK0sMassRange < mass && mass < fMaxK0sMassRange) {
      return true;
    } else {
      return false;
    }
  }

  bool isAcceptK0sCandidateSideBandLeft(double mass) {
    if (fMinK0sMassSideBandLeft < mass && mass < fMaxK0sMassSideBandLeft) {
      return true;
    } else {
      return false;
    }
  }

  bool isAcceptK0sCandidateSideBandRight(double mass) {
    if (fMinK0sMassSideBandRight < mass && mass < fMaxK0sMassSideBandRight) {
      return true;
    } else {
      return false;
    }
  }

  bool isAcceptRhoCandidateMassRange(double mass) {
    if (fMinRhoMassRange < mass && mass < fMaxRhoMassRange) {
      return true;
    } else {
      return false;
    }
  }
  bool isAcceptKstarCandidateMassRange(double mass) {
    if (fMinKstarMassRange < mass && mass < fMaxKstarMassRange) {
      return true;
    } else {
      return false;
    }
  }
  bool isAcceptF980CandidateMassRange(double mass) {
    if (fMinF980MassRange < mass && mass < fMaxF980MassRange) {
      return true;
    } else {
      return false;
    }
  }
  bool isAcceptF1270CandidateMassRange(double mass) {
    if (fMinF1270MassRange < mass && mass < fMaxF1270MassRange) {
      return true;
    } else {
      return false;
    }
  }

  bool isMissPidCandidateFromProtonPion(AliAODTrack *track1,
                                        AliAODTrack *track2, std::string name);
  bool isMissPidCandidateFromKaonPion(AliAODTrack *track1, AliAODTrack *track2,
                                      std::string name);
  bool isMissPidCandidateFromPionPion(AliAODTrack *track1, AliAODTrack *track2,
                                      std::string name);

  bool isSameMotherPair(AliAODTrack *track1, AliAODTrack *track2);
  bool isSameMotherPair(AliAODMCParticle *part1, AliAODMCParticle *part2);
  bool isCharmQuarkOrigin(AliAODMCParticle *particle);
  bool isBeautyQuarkOrigin(AliAODMCParticle *particle);
  bool isHeavyFlavorOrigin(AliAODMCParticle *particle);
  bool isPrimary(AliAODMCParticle *particle);

  int getMotherPdgCode(AliAODTrack *track);
  int getMotherPdgCode(AliAODMCParticle *part);
  int getMotherLabel(AliAODTrack *track);
  int getMotherLabel(AliAODMCParticle *part);

  double getTOFSigma(AliAODTrack *track1, AliPID::EParticleType pid);
  double getTPCSigma(AliAODTrack *track1, AliPID::EParticleType pid);

  bool setTrueCh();

  bool getTrueChPartInV0s(int &v0a, int &v0c) {
    v0a = fNChV0A;
    v0c = fNChV0C;
    return true;
  }

  int getNTrueChTrkInfo(int spec) {
    if (spec == 0) {
      return fNChEta05;
    } else if (spec == 1) {
      return fNChEta10;
    } else if (spec == 2) {
      return fNChEta15;
    } else {
      return fNChEta20;
    }
    return true;
  }

  bool getLeadingTrack(int &iLeading);
  void setMinLeadingTrackPt(double pt) { fMinLeadingTrackPt = pt; }

  double fMinLeadingTrackPt;

  double getSPDTrkCorr(double vtxz, int spec) {
    if (spec == 0) {
      if (!fHistSPDTrkCorrEta05) {
        return 0;
      } else {
        double delta = fHistSPDTrkCorrEta05->GetBinContent(
            fHistSPDTrkCorrEta05->GetXaxis()->FindBin(vtxz));
        return fRandom->Poisson(delta);
      }
    } else if (spec == 1) {
      if (!fHistSPDTrkCorrEta10) {
        return 0;
      } else {
        double delta = fHistSPDTrkCorrEta10->GetBinContent(
            fHistSPDTrkCorrEta10->GetXaxis()->FindBin(vtxz));
        return fRandom->Poisson(delta);
      }
    } else {
      return 0;
    }
  }

  bool isMC() { return fIsMC; }

  void setDalitzProd(bool flag) { fIsDalitzProd = flag; }
  void set2BodyProd(bool flag) { fIs2BodyProd = flag; }

  bool isDalitzProd() { return fIsDalitzProd; }
  bool is2BodyProd() { return fIs2BodyProd; }

  void setMidTrackAna(bool isMidTrack) { fIsMidTrackAna = isMidTrack; }
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Set analysis cut flags
  //////////////////////////////////////////////////////////////////////////////////////////////

  void setK0sCandidateMassRange(double min, double max) {
    fMinK0sMassRange = min;
    fMaxK0sMassRange = max;
  }
  void setRhoCandidateMassRange(double min, double max) {
    fMinRhoMassRange = min;
    fMaxRhoMassRange = max;
  }
  void setF980CandidateMassRange(double min, double max) {
    fMinF980MassRange = min;
    fMaxF980MassRange = max;
  }
  void setF1270CandidateMassRange(double min, double max) {
    fMinF1270MassRange = min;
    fMaxF1270MassRange = max;
  }
  void setKstarCandidateMassRange(double min, double max) {
    fMinKstarMassRange = min;
    fMaxKstarMassRange = max;
  }

  void setPairTargetPIDs(int pid1, int pid2) {
    if (pid1 == 11) {
      fTrackTragetPid1 = AliPID::kElectron;
    } else if (pid1 == 13) {
      fTrackTragetPid1 = AliPID::kMuon;
    } else if (pid1 == 211) {
      fTrackTragetPid1 = AliPID::kPion;
    } else if (pid1 == 321) {
      fTrackTragetPid1 = AliPID::kKaon;
    } else if (pid1 == 2212) {
      fTrackTragetPid1 = AliPID::kProton;
    } else {
      fTrackTragetPid1 = AliPID::kPion;
    }
    if (pid2 == 11) {
      fTrackTragetPid2 = AliPID::kElectron;
    } else if (pid2 == 13) {
      fTrackTragetPid2 = AliPID::kMuon;
    } else if (pid2 == 211) {
      fTrackTragetPid2 = AliPID::kPion;
    } else if (pid2 == 321) {
      fTrackTragetPid2 = AliPID::kKaon;
    } else if (pid2 == 2212) {
      fTrackTragetPid2 = AliPID::kProton;
    } else {
      fTrackTragetPid2 = AliPID::kPion;
    }
  }

  AliPID::EParticleType getPairTargetPIDs(int track) {
    if (track == 0) {
      return fTrackTragetPid1;
    } else {
      return fTrackTragetPid2;
    }
  }

  void setVertexCut(double min, double max, int min_cont) {
    fMinVertexCutZ = min;
    fMaxVertexCutZ = max;
    fMinContVtx = min_cont;
    fIsVtxZcut = true;
  }
  void setPairRapidityCut(double min, double max) {
    fMinPairRapCut = min;
    fMaxPairRapCut = max;
    fIsPairRapCut = true;
  }
  void setPairOpeningAngleCut(double val) { fMaxPairCosOpeningAngleCut = val; }
  void setPairKinematicCut(int type = 0, double min = 0.) {
    if (type == 0) {
      fIsPairPtCutForOneTrack = false;
      fIsPairPtCutForBothTracks = false;
      fMinPairPtCut = 0.;
    } else if (type == 1) {
      fIsPairPtCutForOneTrack = true;
      fIsPairPtCutForBothTracks = false;
      fMinPairPtCut = min;
    } else if (type == 2) {
      fIsPairPtCutForOneTrack = false;
      fIsPairPtCutForBothTracks = true;
      fMinPairPtCut = min;
    } else {
      std::cout << "Pair pt cut is not correct..." << std::endl;
      std::cout << "You set type = " << type
                << "  but it should be between 0 - 2" << std::endl;
    }
  }

  void setTrackKinematicCut(double min_pt, double max_pt, double min_eta,
                            double max_eta, double min_p, double max_p) {
    fMinTrackP = min_p;
    fMaxTrackP = max_p;
    fMinTrackPt = min_pt;
    fMaxTrackPt = max_pt;
    fMinTrackEta = min_eta;
    fMaxTrackEta = max_eta;
  }

  void setTrackQualities(double findable, std::string dcaxy, double dcaz,
                         double chi2tpc, double chi2its, int nclusttpc,
                         int nclustits) {
    fMinCrossRowsFindableRatio = findable;
    fMaxTrackDCAxyName = dcaxy;
    fMaxTrackDCAz = dcaz;
    fMaxReducedChi2TPC = chi2tpc;
    fMaxReducedChi2ITS = chi2its;
    fMinTrackTPCNClusts = nclusttpc;
    fMinTrackSPDNClusts = nclustits;
  }

  void setPileupRejectionCut(bool flag) { fIsPUcut = flag; }
  void setLocalBoardCut(bool flag) { fIsLBCut = flag; }

  void setPionSelectSigmaTPC(double min, double max) {
    fMinPionSigmaTPC = min;
    fMaxPionSigmaTPC = max;
  }
  void setPionSelectSigmaTOF(double min, double max) {
    fMinPionSigmaTOF = min;
    fMaxPionSigmaTOF = max;
  }
  void setKaonSelectSigmaTPC(double min, double max) {
    fMinKaonSigmaTPC = min;
    fMaxKaonSigmaTPC = max;
  }
  void setKaonSelectSigmaTOF(double min, double max) {
    fMinKaonSigmaTOF = min;
    fMaxKaonSigmaTOF = max;
  }
  void setProtonSelectSigmaTPC(double min, double max) {
    fMinProtonSigmaTPC = min;
    fMaxProtonSigmaTPC = max;
  }
  void setProtonSelectSigmaTOF(double min, double max) {
    fMinProtonSigmaTOF = min;
    fMaxProtonSigmaTOF = max;
  }
  void setElectronSelectSigmaTPC(double min, double max) {
    fMinElectronSigmaTPC = min;
    fMaxElectronSigmaTPC = max;
  }
  void setElectronSelectSigmaTOF(double min, double max) {
    fMinElectronSigmaTOF = min;
    fMaxElectronSigmaTOF = max;
  }
  void setMuonSelectSigmaTPC(double min, double max) {
    fMinMuonSigmaTPC = min;
    fMaxMuonSigmaTPC = max;
  }
  void setMuonSelectSigmaTOF(double min, double max) {
    fMinMuonSigmaTOF = min;
    fMaxMuonSigmaTOF = max;
  }
  void setArmenterosLimit(double pcm, double r0, double width) {
    fArmenterosBandWidth = width;
    fArmenterosPCM = pcm;
    fArmenterosR0 = r0;
  }
  void setV0SelectCuts(double alpha, double pangle, double v0Dca,
                       double trackDca, double min_dlength, double max_dlength,
                       double lifetime) {
    fArmenterosAlphaCutParamForPtArm = alpha;
    fMinCosPointingAngleCut = pangle;
    fMinV0DCA = v0Dca;
    fMaxTrackDCASigma = trackDca;
    fMinV0DecayLength = min_dlength;
    fMaxV0DecayLength = max_dlength;
    fMaxV0PropLifeTime = lifetime;
  }

  void setV0CutParams(double min, double max) {
    fMinV0Alpha = min;
    fMaxV0Alpha = max;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Set analysis object
  //////////////////////////////////////////////////////////////////////////////////////////////

  void setMultiEstimateMethod(std::string method) { fMultiMethod = method; }
  void setMuonTrackCut(AliMuonTrackCuts *cut) { fMuonTrackCuts = cut; }
  void setDownScalingHist(TFile *inFile) {
    fHistDsCMLL7 = (TH1F *)inFile->Get("DS_MLL7")->Clone("fHistDsCMLL7");
    fHistDsCMSL7 = (TH1F *)inFile->Get("DS_MSL7")->Clone("fHistDsCMSL7");
    fHistDsCINT7 = (TH1F *)inFile->Get("DS_INT7")->Clone("fHistDsCINT7");
  }
  void setSPDTrkCorrHist(TFile *inFile, std::string period) {
    if (inFile->Get(Form("HistDeltaSPDEta05_CINT7_%s", period.c_str()))) {
      fHistSPDTrkCorrEta05 =
          (TH1D *)inFile
              ->Get(Form("HistDeltaSPDEta05_CINT7_%s", period.c_str()))
              ->Clone("f05");
    }
    if (inFile->Get(Form("HistDeltaSPDEta10_CINT7_%s", period.c_str()))) {
      fHistSPDTrkCorrEta10 =
          (TH1D *)inFile
              ->Get(Form("HistDeltaSPDEta10_CINT7_%s", period.c_str()))
              ->Clone("f10");
    }
  }
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Get the analysis variables
  //////////////////////////////////////////////////////////////////////////////////////////////

  double getTrueVtxZ() { return fTrueVtx[2]; }
  double getVtxX() { return fVtxX; }
  double getVtxY() { return fVtxY; }
  double getVtxZ() { return fVtxZ; }
  double getCentClass(int spec) {
    if (spec == 0)
      return fCentSPDTrk;
    else if (spec == 1)
      return fCentV0M;
    else if (spec == 2)
      return fCentV0A;
    else if (spec == 3)
      return fCentV0C;
    else
      return -999;
  }
  double getVtxCont() { return fNContVtx; }
  double getCentClass() { return fCent; }
  double getPsi() { return fPsi; }
  double getDS() { return fDSfactor; }
  int getRunnumberIndex() { return fRunNumberIndex; }
  int getRunnumber() { return fRunNumber; }
  int getNSPDTrkInfo(int spec) {
    if (spec == 0)
      return fNSPDTrk05;
    else if (spec == 1)
      return fNSPDTrk10;
    else if (spec == 2)
      return fNSPDTrk15;
    else if (spec == 3)
      return fNSPDTrk20;
    else if (spec == 4)
      return fNSPDTrkAll;
    else
      return 0;
  }
  int getNCorrSPDTrkInfo(int spec) {
    if (spec == 0) {
      return fNSPDTrk05 + getSPDTrkCorr(getVtxZ(), spec);
    } else if (spec == 1) {
      return fNSPDTrk10 + getSPDTrkCorr(getVtxZ(), spec);
    } else {
      return 0;
    }
  }
  int getNSPDClustInfo(int spec) {
    if (spec == 0)
      return fNClustSPD1;
    else if (spec == 1)
      return fNClustSPD2;
    else
      return 0;
  }

  double getVzeroInfo(int spec) {
    if (spec == 0)
      return fChV0A;
    else if (spec == 1)
      return fChV0C;
    else if (spec == 2)
      return fChV0M;
    else if (spec == 3)
      return fTimeV0A - fTimeV0C;
    else if (spec == 4)
      return fTimeV0A + fTimeV0C;
    else
      return 0;
  }

  std::string getPass() { return fPass; }
  void getPeriodInfo(std::string &period, std::string &system) {
    period = fPeriod;
    system = fCollSystem;
  }
  void getTriggerInfo(bool &isCINT7, bool &isCMSL7, bool &isCMSH7,
                      bool &isCMUL7, bool &isCMLL7, bool &is0MSL, bool &is0MSH,
                      bool &is0MUL, bool &is0MLL) {
    is0MSL = fIs0MSL;
    is0MSH = fIs0MSH;
    is0MUL = fIs0MUL;
    is0MLL = fIs0MLL;

    isCINT7 = fIsCINT7;
    isCMSL7 = fIsCMSL7;
    isCMSH7 = fIsCMSH7;
    isCMUL7 = fIsCMUL7;
    isCMLL7 = fIsCMLL7;
  }
  void getTriggerInfo(bool &isCINT7, bool &isCMSL7, bool &isCMSH7,
                      bool &isCMUL7, bool &isCMLL7) {
    isCINT7 = fIsCINT7;
    isCMSL7 = fIsCMSL7;
    isCMSH7 = fIsCMSH7;
    isCMUL7 = fIsCMUL7;
    isCMLL7 = fIsCMLL7;
  }

  AliMuonTrackCuts *getMuonTrackCuts() { return fMuonTrackCuts; }

  const int fPdgCodeEta = 221;
  const int fPdgCodeRho = 113;
  const int fPdgCodeOmega = 223;
  const int fPdgCodeEtaPrime = 331;
  const int fPdgCodeKstar = 313;
  const int fPdgCodePhi = 333;

  const double fPdgLambdaMass = 1.115683;
  const double fPdgK0sMass = 0.497611;

  // private:

  void setInit();
  bool isSameRunnumber();
  bool setVtxZCentPsi();
  bool setDownScaleFactor();
  bool setPeriodInfo();
  bool setRunnumberIndex();
  bool setTriggerInfo();
  bool setSPDTrk();
  bool setSPDClust();
  bool setVZERO();

  AliAODEvent *fEvent;
  AliMultSelection *fMultSelection;
  AliMuonTrackCuts *fMuonTrackCuts;
  AliVEventHandler *fInputHandler;
  TClonesArray *fMCArray;

  int fRunNumber;
  int fRunNumberIndex;

  std::string fPeriod;
  std::string fCollSystem;
  std::string fPass;

  AliPID::EParticleType fTrackTragetPid1;
  AliPID::EParticleType fTrackTragetPid2;

  TF1 *fFuncMaxDCAxy;
  TF1 *fMinArmenterosLine;
  TF1 *fMaxArmenterosLine;
  double fArmenterosBandWidth;
  double fArmenterosPCM;
  double fArmenterosR0;
  double fArmenterosAlphaCutParamForPtArm;
  double fMinCosPointingAngleCut;
  double fMinV0DCA;
  double fMaxTrackDCASigma;
  double fMinV0DecayLength;
  double fMaxV0DecayLength;
  double fMaxV0PropLifeTime;
  double fMinV0DecayRadius;
  double fMinV0Alpha;
  double fMaxV0Alpha;

  double fMinK0sMassRange;
  double fMaxK0sMassRange;
  double fMinRhoMassRange;
  double fMaxRhoMassRange;
  double fMinKstarMassRange;
  double fMaxKstarMassRange;
  double fMinF980MassRange;
  double fMaxF980MassRange;
  double fMinF1270MassRange;
  double fMaxF1270MassRange;

  double fMinK0sMassSideBandLeft;
  double fMaxK0sMassSideBandLeft;
  double fMinK0sMassSideBandRight;
  double fMaxK0sMassSideBandRight;

  double fMinRejectMassWidthLambda;
  double fMaxRejectMassWidthLambda;

  double fMinCrossRowsFindableRatio;
  std::string fMaxTrackDCAxyName;
  double fMaxTrackDCAxy;
  double fMaxTrackDCAz;
  double fMaxReducedChi2TPC;
  double fMaxReducedChi2ITS;
  int fMinTrackTPCNClusts;
  int fMinTrackSPDNClusts;

  bool fIsMC;
  bool fIsEvtSelect;

  bool fIsVtxZcut;
  double fMaxVertexCutZ;
  double fMinVertexCutZ;
  int fNContVtx;
  int fMinContVtx;

  bool fIsPairRapCut;
  double fMinPairRapCut;
  double fMaxPairRapCut;
  double fMaxPairCosOpeningAngleCut;

  bool fIsPairPtCutForOneTrack;
  bool fIsPairPtCutForBothTracks;
  double fMinPairPtCut;

  bool fIsPUcut;
  bool fIsLBCut;

  bool fIsDalitzProd;
  bool fIs2BodyProd;

  std::string fMultiMethod;

  TH1F *fHistDsCMSL7;
  TH1F *fHistDsCINT7;
  TH1F *fHistDsCMLL7;

  double fVtxX;
  double fVtxY;
  double fVtxZ;
  double fCent;
  double fPsi;

  double fTrueVtx[3];

  double fCentSPDTrk;
  double fCentV0A;
  double fCentV0C;
  double fCentV0M;

  double fDSfactor;

  bool fIs0MSL;
  bool fIs0MSH;
  bool fIs0MUL;
  bool fIs0MLL;

  bool fIsCINT7;
  bool fIsCMSL7;
  bool fIsCMSH7;
  bool fIsCMUL7;
  bool fIsCMLL7;

  int fInput0MSH;
  int fInput0MLL;
  int fInput0MUL;
  int fInput0MSL;

  int fNSPDTrk05;
  int fNSPDTrk10;
  int fNSPDTrk15;
  int fNSPDTrk20;
  int fNSPDTrkAll;

  int fNClustSPD1;
  int fNClustSPD2;

  double fChV0A;
  double fChV0C;
  double fChV0M;
  double fTimeV0A;
  double fTimeV0C;

  int fNChV0A;
  int fNChV0C;

  int fNChEta05;
  int fNChEta10;
  int fNChEta15;
  int fNChEta20;

  AliPIDResponse *fPIDResponse;
  TRandom1 *fRandom;

  TH1D *fHistSPDTrkCorrEta05;
  TH1D *fHistSPDTrkCorrEta10;

  double fMinTrackP;
  double fMaxTrackP;
  double fMinTrackPt;
  double fMaxTrackPt;
  double fMinTrackEta;
  double fMaxTrackEta;

  double fMinPionSigmaTPC;
  double fMaxPionSigmaTPC;
  double fMinPionSigmaTOF;
  double fMaxPionSigmaTOF;

  double fMinKaonSigmaTPC;
  double fMaxKaonSigmaTPC;
  double fMinKaonSigmaTOF;
  double fMaxKaonSigmaTOF;

  double fMinProtonSigmaTPC;
  double fMaxProtonSigmaTPC;
  double fMinProtonSigmaTOF;
  double fMaxProtonSigmaTOF;

  double fMinElectronSigmaTPC;
  double fMaxElectronSigmaTPC;
  double fMinElectronSigmaTOF;
  double fMaxElectronSigmaTOF;

  double fMinMuonSigmaTPC;
  double fMaxMuonSigmaTPC;
  double fMinMuonSigmaTOF;
  double fMaxMuonSigmaTOF;

  bool fIsMidTrackAna;

  ClassDef(AliAnalysisTaskAODTrackPairUtils, 1); // example of analysis
};

#endif
