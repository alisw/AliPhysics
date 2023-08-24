/// \class AliAnalysisTaskHe3piAOD

#ifndef __AliAnalysisTaskHe3piAOD__
#define __AliAnalysisTaskHe3piAOD__

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>
#include <TString.h>
#include "AliEventCuts.h"

class AliPIDResponse;
class TH2F;
class TList;
class TTree;

struct MiniHyper
{
  Double32_t pt;
  Double32_t Rapidity;
  Double32_t m;
  Double32_t ct;
  Double32_t V0CosPA;
  Double32_t V0radius;
  Double32_t Lrec;
  Double32_t fZ;
  Double32_t TPCnSigmaPi;     //[-5,5,8]
  Double32_t TPCnSigmaHe3;    //[-5,5,8]
  Double32_t TPCmomHe3;       //[0.0,10.22,8]
  Double32_t TPCsignalHe3;    //[0.,2046.,12]
  Double32_t He3ProngPvDCAXY; //[0.0,5.10,10]
  Double32_t PiProngPvDCAXY;  //[0.0,20.46,12]
  Double32_t He3ProngPvDCA;   //[0.0,10.22,10]
  Double32_t PiProngPvDCA;    //[0.0,40.94,12]
  Double32_t ProngsDCA;       //[0.0,2.54,8]
  unsigned char NpidClustersPion;
  unsigned char NpidClustersHe3;
  unsigned char NitsClustersHe3;
  unsigned char centrality;
  unsigned char trigger;
  bool Matter;
};

struct MiniHyperMC : public MiniHyper
{
  float ptMC;
  float etaMC;
  float ctMC;
  float yMC;
  int pdg;
  bool isReconstructed;
  bool isDuplicated = false;
};

class AliAnalysisTaskHe3piAOD : public AliAnalysisTaskSE
{
public:
  enum kReducedTrigger
  {
    kINT7 = BIT(0),
    kCentral = BIT(1),
    kSemiCentral = BIT(2),
    kPositiveB = BIT(3),
    kHighMultV0 = BIT(4)
  };
  AliAnalysisTaskHe3piAOD(bool isMC = false, TString taskname = "HyperAOD");
  static AliAnalysisTaskHe3piAOD *AddTask(bool isMC, TString tskname, TString suffix);
  virtual ~AliAnalysisTaskHe3piAOD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {}

  AliEventCuts fEventCut; ///<

  void SetCustomBetheBloch(float resolution, const float bethe[5]);
  double customNsigma(double mom, double sig);
  void SaveOnlyTrueCandidates(bool toggle = true) { fOnlyTrueCandidates = toggle; }
  void UseOnTheFly(bool toggle = true) { fUseOnTheFly = toggle; }
  void SetMassRange(float min, float max)
  {
    fMassRange[0] = min;
    fMassRange[1] = max;
  }

private:
  AliAnalysisTaskHe3piAOD(const AliAnalysisTaskHe3piAOD &source);
  AliAnalysisTaskHe3piAOD &operator=(const AliAnalysisTaskHe3piAOD &source);

  void PostAllData();

  TList *fList = nullptr; //!<! List of the output histograms
  TTree *fTree = nullptr; //!<! Tree for Xis and Omegas

  MiniHyper *fRecHyper = nullptr; //!<! Transient fRecHyper
  MiniHyperMC fGenHyper;
  AliPIDResponse *fPID = nullptr; //!<! ALICE PID framework
  bool fMC;
  bool fOnlyTrueCandidates = true; ///< Save only true Hypertritons
  bool fUseOnTheFly = false;

  float fCustomBethe[5] = {-166.11733, -0.11020473, 0.10851357, 2.7018593, -0.087597824}; /// default values are for LHC18qr
  float fCustomResolution = 0.05871;                                                      /// default values are for LHC18qr
  double fMassRange[2] = {2.9, 3.1};

  float Eta2y(float pt, float m, float eta) const;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskHe3piAOD, 1);
  /// \endcond
};

#endif /* defined(__AliAnalysisTaskHe3piAOD__) */
