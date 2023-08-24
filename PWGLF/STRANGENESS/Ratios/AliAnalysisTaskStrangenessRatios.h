/// \class AliAnalysisTaskStrangenessRatios

#ifndef __AliAnalysisTaskStrangenessRatios__
#define __AliAnalysisTaskStrangenessRatios__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include "AliEventCuts.h"
#include "AliExternalBDT.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"

class AliPIDResponse;
class TH2F;
class TList;
class TTree;

struct MiniLambda {
  Double32_t pt;
  Double32_t ptPr;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius;      //[0,101.6,8]
  Double32_t dcaV0PV;     //[0,10.16,16]
  Double32_t dcaPiPV;     //[0,20.32,8]
  Double32_t dcaPrPV;     //[0,20.32,8]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t cosPA;       //[0.95,1,16]
  Double32_t tpcNsigmaPi; //[-5,5,8]
  Double32_t tpcNsigmaPr; //[-5,5,8]
  Double32_t itsNsigmaPr; //[-5,5,8]
  unsigned char tpcClV0Pi;
  unsigned char tpcClV0Pr;
  unsigned char centrality;
  bool matter;
  bool hasTOFhit;
  bool hasITSrefit;
};

struct MiniK0s {
  Double32_t pt;
  Double32_t ptPos;
  Double32_t ptNeg;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius;      //[0,101.6,8]
  Double32_t dcaV0PV;     //[0,10.16,16]
  Double32_t dcaPosPV;    //[0,20.32,8]
  Double32_t dcaNegPV;    //[0,20.32,8]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t cosPA;       //[0.95,1,16]
  Double32_t tpcNsigmaPos; //[-5,5,8]
  Double32_t tpcNsigmaNeg; //[-5,5,8]
  unsigned char tpcClV0Pos;
  unsigned char tpcClV0Neg;
  unsigned char centrality;
  bool hasTOFpos;
  bool hasTOFneg;
  bool hasTOFhit;
  bool hasITSrefit;
};

struct MiniLambdaMC : public MiniLambda {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  float ptMotherMC;
  float ctMotherMC;
  int pdg;
  bool isPrimary;
  bool isReconstructed;
  bool magFieldPolarity;
  unsigned char flag;
};

struct MiniK0sMC : public MiniK0s {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  bool isPrimary;
  bool isReconstructed;
  unsigned char flag;
};

struct MiniCascade {
  Double32_t pt;
  Double32_t eta;
  Double32_t mass;
  Double32_t ct;
  Double32_t radius; //[0,51.1,8]
  Double32_t radiusV0; //[0,102.0,8]
  Double32_t dcaBachPV; //[0,12.75,16]
  Double32_t dcaV0PV; //[0,10.16,16]
  Double32_t dcaV0piPV; //[0,25.5,16]
  Double32_t dcaV0prPV; //[0,12.75,16]
  Double32_t dcaV0tracks; //[0,2.54,8]
  Double32_t dcaBachV0; //[0,2.54,8]
  Double32_t cosPA; //[0.95,1,16]
  Double32_t cosPAV0; //[0.95,1,16]
  Double32_t V0invMassDelta; //[-0.01,0.01,8]
  Double32_t tpcNsigmaBach; //[-5,5,8]
  Double32_t tpcNsigmaV0Pr; //[-5,5,8]
  Double32_t tpcNsigmaV0Pi; //[-5,5,8]
  Double32_t competingMass; //[0,0.254,8]
  Double32_t bachBarCosPA;  //[0.9999, 1., 8]
  unsigned char tpcClBach;
  unsigned char tpcClV0Pr;
  unsigned char tpcClV0Pi;
  unsigned char centrality;
  bool matter;
  bool hasTOFhit;
  bool hasITSrefit;
  bool isOmega;
};

struct MiniCascadeMC : public MiniCascade {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
  unsigned char flag;
};

struct MiniCascadeTrain {
  float pt;
  float mass;
  double dcaBachPV;
  double dcaV0PV;
  double dcaV0piPV;
  double dcaV0prPV;
  double dcaV0tracks;
  double dcaBachV0;
  double cosPA;
  double cosPAV0;
  double tpcNsigmaV0Pr;
  unsigned char centrality;
};

struct MiniCascadeTrainMC : public MiniCascadeTrain {
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
  unsigned char flag;
};

struct MiniLambdaBDTOut {
  Double32_t pt;
  Double32_t ct;
  Double32_t mass;
  Double32_t bdtOutputPrompt;
  Double32_t bdtOutputNonPrompt;
  Double32_t bdtOutputBackground;
  bool matter;
  bool pileUpCheck;
};

class AliAnalysisTaskStrangenessRatios : public AliAnalysisTaskSE {
public:
  enum StatusFlag {
    kPrimary = BIT(0),
    kSecondaryFromWDXi = BIT(1),
    kSecondaryFromWDOmega = BIT(2),
    kSecondaryFromWD = BIT(3),
    kSecondaryFromMaterial = BIT(4)
  };

  AliAnalysisTaskStrangenessRatios(bool isMC = false, TString taskname = "StrangenessRatios");
  static AliAnalysisTaskStrangenessRatios* AddTask(bool isMC, TString tskname, TString suffix);
  virtual ~AliAnalysisTaskStrangenessRatios();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual Bool_t UserNotify();
  virtual void   Terminate(Option_t *) {}

  AliEventCuts  fEventCut; ///<

  void SetFillK0s(bool toogle = true) { fFillK0s = toogle; }
  void SetFillLambdas(bool toogle = true) { fFillLambdas = toogle; }
  void SetFillCascades(bool toogle = true) { fFillCascades = toogle; }
  void SetFillCascadesTrain(bool toogle = true) { fFillCascadesTrain = toogle; }
  void SetFillLambdasBDTOut(bool toogle = true) { fFillLambdasBDTOut = toogle; }
  void SetLambdaDownscaling(float scale = true) { fLambdaDownscaling = scale; }
  void SetK0sDownscaling(float scale = true) { fK0sDownscaling = scale; }

  //Setters for topological cuts
  void SetRadiusCut(float xi = 1.2, float omega = 1., float lambda = 0.5) {fCutRadius[0]=xi; fCutRadius[1]=omega; fCutRadius[2]=lambda;}
  void SetRadiusV0Cut(float cut = 3.0) {fCutRadiusV0=cut;}
  void SetDCABachToPVCut(float cut = 0.1) {fCutDCABachToPV=cut;}
  void SetDCAV0toPVCut(float cut = 0.1) {fCutDCAV0toPV=cut;}
  void SetDCAV0piToPVCut(float cut = 0.2) {fCutDCAV0piToPV=cut;}
  void SetDCAV0prToPVCut(float cut = 0.2) {fCutDCAV0prToPV=cut;}
  void SetDCAV0tracksCut(float cut = 1.0) {fCutDCAV0tracks=cut;}
  void SetCutDCABachToV0XiCut(float xi = 1.0, float omega = 0.6) {fCutDCABachToV0[0]=xi; fCutDCABachToV0[1]=omega;}
  void SetCosPACut(float cut = 0.95) {fCutCosPA=cut;}
  void SetCosPAV0Cut(float cut = 0.95) {fCutCosPAV0=cut;}
  void SetV0MassWindowCut(float cut = 0.005) {fCutV0MassWindow=cut;}
  void SetYCut(float cut = 0.5) {fCutY=cut;}
  void SetYDaughtCut(float cut = 0.8) {fCutYDaught=cut;}
  void SetNsigmaTPCCut(float cut = 4.0) {fCutNsigmaTPC=cut;}
  void SetCtCut(float cut = 4 /*cm*/) {fCutCt=cut;}
  void SetCtV0Cut(float cut = 30 /*cm*/) {fCutCtV0=cut;}
  void SetCompetingMassCut(float omega = 0.008 /*GeV/c^2*/, float xi = 0. /*GeV/c^2*/) {fCutCompetingMass[0]=xi; fCutCompetingMass[1]=omega;}
  void SetTPCcluCut(int cut = 70) {fCutTPCclu=cut;}
  void SetSaveOnlyTrueCandidates(bool cut = true) { fOnlyTrueCandidates = cut; }
  void SetSaveOnlyTrueLambdas(bool cut = true) { fOnlyTrueLambdas = cut; }
  void SetTPCRowsCut(float cut = 80.) {fCutTPCrows=cut;}
  void SetTPCRowOvFCut(float cut = 0.8) {fCutRowsOvF=cut;}
  void UseOnTheFly(bool toggle = true) { fUseOnTheFly = toggle; }
  void SetMinCentrality(int minCentrality = 0) { fMinCentrality = minCentrality; }
  void SetMaxCentrality(int maxCentrality = 90) { fMaxCentrality = maxCentrality; }
  void SetCtPreselection(double cut = 1.) { fCtPreselection = cut; }
  void SetMaxCt(double maxCt = 40.) { fMaxCt = maxCt; }
  void SetMinPt(double minPt = 0.5) { fMinPt = minPt; }
  void SetMaxPt(double maxPt = 3.5) { fMaxPt = maxPt; }
  void SetRadiusPreselection(double cut = 3.) { fRadiusPreselection = cut; }
  void SetRadiusOverflowCut(double cut = 100.) { fRadiusOverflowCut = cut; }
  void SetTpcClV0PiPreselection(int cut = 70) { fTpcClV0PiPreselection = cut; }
  void SetTpcClV0PrPreselection(int cut = 70) { fTpcClV0PrPreselection = cut; }
  void SetDCAV0piToPVOverflowCut(double cut = 20.) { fDCAV0piToPVOverflowCut = cut; }
  void SetDCAV0prToPVOverflowCut(double cut = 20.) { fDCAV0prToPVOverflowCut = cut; }
  void SetDCAV0toPVOverflowCut(double cut = 10.) { fDCAV0toPVOverflowCut = cut; }
  void SetBdtOutputBackgroundCut(double cut = 0.15) { fBdtOutputBackgroundCut = cut; }
  void SetSDDSSDclsCut(int cut = -1) { fSDDSSDclsCut = cut; }
  void SetSPDclsCut(int cut = -1) { fSPDclsCut = cut; }

  void SetUseAdditionalCuts(int toggle = true) { fUseAdditionalCuts = toggle; }

  void SetBDTPath(const char *path = "") { fBDTPath = path; }
  void SetCtBinsBDT(int nBins, double *ctBins) { fCtBinsBDT.Set(nBins+1,ctBins); }

private:
  AliAnalysisTaskStrangenessRatios (const AliAnalysisTaskStrangenessRatios &source);
  AliAnalysisTaskStrangenessRatios &operator=(const AliAnalysisTaskStrangenessRatios &source);

  void PostAllData();

  TList*          fList = nullptr;             //!<! List of the output histograms
  TTree*          fTree = nullptr;             //!<! Tree for Xis and Omegas
  TTree*          fTreeTrain = nullptr;        //!<! Tree for Xis and Omegas (training)
  TTree*          fTreeLambda = nullptr;       //!<! Tree for Lambdas
  TTree*          fTreeK0s = nullptr;          //!<! Tree for K0s
  TTree*          fTreeLambdaBDTOut = nullptr; //!<! Tree for Lambdas

  MiniCascade* fRecCascade = nullptr;          //!<! Transient fRecCascade
  MiniCascadeMC fGenCascade;
  MiniCascadeTrain* fRecCascadeTrain = nullptr;//!<! Transient fRecCascadeTrain
  MiniCascadeTrainMC fGenCascadeTrain;
  MiniLambda* fRecLambda = nullptr;          //!<! Transient fRecLambda
  MiniLambdaMC fGenLambda;
  MiniK0s* fRecK0s = nullptr;          //!<! Transient fRecK0s
  MiniK0sMC fGenK0s;
  MiniLambdaBDTOut *fRecLambdaBDTOut = nullptr;//!<! Transient fRecLambda with BDT output
  AliPIDResponse* fPID = nullptr;              //!<! ALICE PID framework
  std::vector<AliExternalBDT*> fBDT;           //!<! BDTs
  bool fMC;
  bool fOnlyTrueCandidates = false;  ///< Save only true Xi and Omegas in MC
  bool fFillK0s = false;
  bool fFillLambdas = false;
  bool fFillLambdasBDTOut = false;
  bool fFillCascades = false;
  bool fFillCascadesTrain = false;
  bool fOnlyTrueK0s = true;          ///< Save only true K0s in MC
  bool fOnlyTrueLambdas = true;      ///< Save only true Lambdas in MC
  float fLambdaDownscaling = 1.;
  float fK0sDownscaling = 1.;

  //configurable cuts
  float fCutRadius[3] = {1.2, 1.0, 3.0};
  float fCutRadiusV0 = 3.0;
  float fCutDCABachToPV = 0.1;
  float fCutDCAV0toPV = 0.1;
  float fCutDCAV0piToPV = 0.2;
  float fCutDCAV0prToPV = 0.2;
  float fCutDCAV0tracks = 1.2;
  float fCutDCABachToV0[2]{1.0, 1.0};
  float fCutCosPA = 0.95;
  float fCutCosPAV0 = 0.95;
  float fCutV0MassWindow = 0.005;
  float fCutY = 0.5;
  float fCutYDaught = 0.8;
  float fCutNsigmaTPC = 4.0;
  float fCutCt = 4;
  float fCutCtV0 = 30;
  float fCutCompetingMass[2]{0.,0.};
  float fCutBachBarCosPA = 0.99995;
  int fCutTPCclu = 70;
  float fCutTPCrows = 80.;
  float fCutRowsOvF = 0.8;
  double fCascLeastCRaws;
  double fCascLeastCRawsOvF;
  double fV0LeastCRaws;
  double fV0LeastCRawsOvF;
  double fMinCentrality = 0;
  double fMaxCentrality = 90;
  double fCtPreselection = 1.;
  double fMaxCt = 40.;
  double fMinPt = 0.5;
  double fMaxPt = 3.5;
  double fRadiusPreselection = 3;
  double fRadiusOverflowCut = 100;
  int fTpcClV0PiPreselection = 70;
  int fTpcClV0PrPreselection = 70;
  double fDCAV0piToPVOverflowCut = 20;
  double fDCAV0prToPVOverflowCut = 20;
  double fDCAV0toPVOverflowCut = 10;
  double fBdtOutputBackgroundCut = 0.15;
  int fSPDclsCut = -1;
  int fSDDSSDclsCut = -1;

  float fCosPALambda = 0.97;
  float fCutDCALambdaPrToPV = 0.08;
  float fCutDCALambdaPiToPV = 0.12;
  float fCutLambdaMass[2] = {1.09f, 1.15f};

  float fCosPAK0s = 0.97;
  float fCutDCAK0sProngToPV = 0.08;
  float fCutK0sMass[2] = {0.497611 - 0.015, 0.497611 + 0.015};

  bool fUseAdditionalCuts = false;

  bool fUseOnTheFly = false;

  TArrayD fCtBinsBDT;
  std::string fBDTPath = "";

  bool IsTopolSelected(bool isXi = true);
  bool IsTopolSelectedLambda();
  bool IsTopolSelectedK0s();
  float Eta2y(float pt, float m, float eta) const;
  int WhichBDT(double ct);
  void FindWDLambdaMother(AliAODMCParticle *track);
  bool HasTOF(AliVTrack *track);


  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskStrangenessRatios, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskStrangenessRatios__) */
