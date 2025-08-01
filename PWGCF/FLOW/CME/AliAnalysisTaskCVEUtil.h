/**************************************************************************
 * Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//--------------------------------------------------------------------------------
// Utilities for CVE and PID CME analysis
// Contributor: Chunzheng Wang, <chunzheng.wang@cern.ch>, Shanghai
//--------------------------------------------------------------------------------

#ifndef AliAnalysisTaskCVEUtil_h
#define AliAnalysisTaskCVEUtil_h

#include <TH1.h>
#include <TString.h>
#include <array>
#include <map>
#include <memory>

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "TF1.h"
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TProfile3D.h>

// ------------------------------------------------------------
//  ENUMS
// ------------------------------------------------------------
enum class ParticleType {
  kPosHadron, kNegHadron, kPosPion, kNegPion,
  kProton, kAntiProton, kLambda, kAntiLambda,
  COUNT           ///< keep last – number of elements
};

// ------------------------------------------------------------
//  ParticleHistos – container holding *all* histograms for one
//  particle species (both data and MC-resolved variants).
// ------------------------------------------------------------
struct MCParticleHists {
  // ---- pT -
  TH2F *h2_pt_mc              {nullptr}; //!<!
  TH2F *h2_pt_rc              {nullptr}; //!<!
  TH2F *h2_pt_rc_real         {nullptr}; //!<!

  // ---- pT dca ----
  TH3F *h3_pt_dcaXY_origin_rc       {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_origin_rc        {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_material_rc     {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_material_rc      {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_lambda_rc       {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_lambda_rc        {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_other_rc        {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_other_rc         {nullptr}; //!<!

  TH3F *h3_pt_dcaXY_origin_rc_real       {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_origin_rc_real        {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_material_rc_real     {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_material_rc_real      {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_lambda_rc_real       {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_lambda_rc_real        {nullptr}; //!<!
  TH3F *h3_pt_dcaXY_other_rc_real        {nullptr}; //!<!
  TH3F *h3_pt_dcaZ_other_rc_real         {nullptr}; //!<!

  void AddToList(TList* list) const {
    for (auto h2 : { h2_pt_mc, h2_pt_rc, h2_pt_rc_real}) {
      if (h2) list->Add(h2);
    }
    for (auto h3 : { h3_pt_dcaXY_origin_rc, h3_pt_dcaZ_origin_rc,
                     h3_pt_dcaXY_material_rc, h3_pt_dcaZ_material_rc,
                     h3_pt_dcaXY_lambda_rc, h3_pt_dcaZ_lambda_rc,
                     h3_pt_dcaXY_other_rc, h3_pt_dcaZ_other_rc}) {
      if (h3) list->Add(h3);
    }
    for (auto h3 : { h3_pt_dcaXY_origin_rc_real, h3_pt_dcaZ_origin_rc_real,
                     h3_pt_dcaXY_material_rc_real, h3_pt_dcaZ_material_rc_real,
                     h3_pt_dcaXY_lambda_rc_real, h3_pt_dcaZ_lambda_rc_real,
                     h3_pt_dcaXY_other_rc_real, h3_pt_dcaZ_other_rc_real}) {
      if (h3) list->Add(h3);
    }
  }
};

struct DataParticleHists {
  // ---- pt dca ----
  TH2F *h2_pt                        {nullptr}; //!<!
  TH3F *h3_pt_dcaXY                  {nullptr}; //!<!
  TH3F *h3_pt_dcaZ                   {nullptr}; //!<!
  TH3F *h3_pt_phi                    {nullptr}; //!<!

  void AddToList(TList* list) const {
    if (h2_pt) list->Add(h2_pt);
    for (auto h3 : { h3_pt_dcaXY, h3_pt_dcaZ, h3_pt_phi}) {
        if (h3) list->Add(h3);
    }
  }
};

using MCParticleHistsMap = std::map<ParticleType, MCParticleHists>;
using DataParticleHistsMap = std::map<ParticleType, DataParticleHists>;

// Helper to convert enum → name string (used in histogram titles)
inline const char *ParticleName(ParticleType t) {
  switch (t) {
    case ParticleType::kPosHadron:  return "poshadron";
    case ParticleType::kNegHadron:  return "neghadron";
    case ParticleType::kPosPion:    return "pospion";
    case ParticleType::kNegPion:    return "negpion";
    case ParticleType::kProton:     return "proton";
    case ParticleType::kAntiProton: return "antiproton";
    case ParticleType::kLambda:     return "lambda";
    case ParticleType::kAntiLambda: return "antilambda";
    default:                       return "unknown";
  }
}

// ------------------------------------------------------------
//  AliAnalysisTaskCVEUtil  (header-only definition)
// ------------------------------------------------------------
class AliAnalysisTaskCVEUtil : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskCVEUtil();
  AliAnalysisTaskCVEUtil(const char *name);
  virtual ~AliAnalysisTaskCVEUtil();

  //----------------------------------------------------------------
  //  Framework hooks
  //----------------------------------------------------------------
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*);

  void SetPeriod(TString period) { fPeriod = period; }
  void SetMC(bool isMC) { fIsMC = isMC; }

private:
  //----------------------------------------------------------------
  //  Internal helpers (declarations only – implementation in .cxx)
  //----------------------------------------------------------------
  bool          LoopTracks();
  bool          LoopV0s();
  bool          ProcessMCParticles();

  bool          RejectEvtTFFit(float centSPD0);
  bool          AcceptAODTrack(AliAODTrack *track);
  bool          CheckPIDofParticle(AliAODTrack *track, int pidToCheck);
  bool          GetDCA(float &dcaXY, float &dcaZ, AliAODTrack *track);

  bool          IsGoodV0(AliAODv0 *v0);
  bool          IsGoodDaughterTrack(const AliAODTrack *track);
  int           GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *nTrack);

  bool          LoadCalibHistForThisRun(); // deal with all the readin
  bool          LoadNUEGraphForThisCent(); // deal with all the readin
  float         GetNUECor(int charge, float pt);
  float         GetNUACor(int charge, float phi, float eta, float vz);

  //----------------------------------------------------------------
  //  Switches / config flags
  //----------------------------------------------------------------
  std::map<int, int>*  runNumList{nullptr};  //!<! List of run numbers
  bool        fDebug              {false};
  bool        fIsMC               {false};
  bool        fUseNarrowDCA768    {true};
  bool        fCustomProtonDCA    {true};

  //----------------------------------------------------------------
  //  Global cut parameters
  //----------------------------------------------------------------
  TString     fTrigger            {"kINT7+kSemiCentral"};
  TString     fPeriod             {"LHC18q"};

  float       fVzCut              {10.f};
  int         fFilterBit          {768};
  int         fNclsCut            {70};
  float       fChi2Max            {2.5f};
  float       fChi2Min            {0.1f};
  float       fNSigmaTPC          {3.f};
  float       fNSigmaRMS          {3.f};


  ////////////////////////
  // NUE
  ////////////////////////
  TList* fListNUE{nullptr};
  // Hadron
  TGraphErrors* gNUEPosHadron_thisCent{nullptr}; //!<!
  TGraphErrors* gNUENegHadron_thisCent{nullptr}; //!<!
  ////////////////////////
  // NUA
  ////////////////////////
  TList* fListNUA{nullptr};      //!<! read lists for NUA
  TH3F* hCorrectNUAPos{nullptr}; //!<!
  TH3F* hCorrectNUANeg{nullptr}; //!<!

  //----------------------------------------------------------------
  //  Runtime handles
  //----------------------------------------------------------------
  AliAODEvent        *fAOD           {nullptr};
  AliPIDResponse     *fPIDResponse   {nullptr};
  AliMultSelection   *fMultSel       {nullptr};
  AliMCEvent         *fMCEvent       {nullptr};
  AliEventCuts       *fEventCuts     {nullptr};

  //----------------------------------------------------------------
  //  Event-level variables (cached each entry)
  //----------------------------------------------------------------
  int fRunNum{-1};
  int fOldRunNum{-2};
  int fCentBin{-1};
  int fOldCentBin{-2};

  std::array<double,3> fVertex  {0.,0.,0.};
  float       fCent            {0.f};

  std::vector<ParticleType>    fParticles{
    ParticleType::kPosHadron,
    ParticleType::kNegHadron,
    ParticleType::kPosPion,
    ParticleType::kNegPion,
    ParticleType::kProton,
    ParticleType::kAntiProton,
    ParticleType::kLambda,
    ParticleType::kAntiLambda,
  };

  //----------------------------------------------------------------
  //  Pile-up parameterisations (unique_ptr for automatic cleanup)
  //----------------------------------------------------------------
  std::unique_ptr<TF1> fSPDCutPU{nullptr};        //!<!
  std::unique_ptr<TF1> fV0CutPU{nullptr};         //!<!
  std::unique_ptr<TF1> fCenCutLowPU{nullptr};     //!<!
  std::unique_ptr<TF1> fCenCutHighPU{nullptr};    //!<!
  std::unique_ptr<TF1> fMultCutPU{nullptr};       //!<!

  //----------------------------------------------------------------
  //  Histogram container
  //----------------------------------------------------------------
  TList                   *fOutputList {nullptr}; //!<! owner list
  TH1I                    *fHistRunNumBin{nullptr}; //!<! run number histogram
  MCParticleHistsMap      *fMCHists    {nullptr}; //!<! all mc histograms
  DataParticleHistsMap    *fDataHists  {nullptr}; //!<! all data histograms

  // TPC recenter
  float fVzBin{0.f};
  float fTPCQx{0.f}; //!<<! TPC recenter X
  float fTPCQy{0.f}; //!<<! TPC recenter Y
  float fTPCM{0.f}; //!<<! TPC recenter Z
  TProfile3D* fProfile3DQxTPCRunCentVz{nullptr}; //!<! TPC recenter profile
  TProfile3D* fProfile3DQyTPCRunCentVz{nullptr}; //!<! TPC recenter profile
  TH2F* fHistQxQyTPC{nullptr}; //!<! TPC Qx Qy
  //----------------------------------------------------------------
  //  Histogram creation helpers
  //----------------------------------------------------------------
  void CreateAllHistograms();

  //----------------------------------------------------------------
  //  Rule of five – copy disabled, move default
  //----------------------------------------------------------------
  AliAnalysisTaskCVEUtil(const AliAnalysisTaskCVEUtil &);
  AliAnalysisTaskCVEUtil &operator=(const AliAnalysisTaskCVEUtil &);

  ClassDef(AliAnalysisTaskCVEUtil, 4)
};

#endif // AliAnalysisTaskCVEUtil_h
