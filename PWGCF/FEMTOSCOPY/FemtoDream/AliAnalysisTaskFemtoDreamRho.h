/*
 * AliAnalysisTaskFemtoDreamRho.cxx
 *
 *  Created on: 19 Jul 2023
 *      Author: M. Korwieser
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_
#include "Rtypes.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include <TString.h>
#include <utility>

class AliAnalysisTaskFemtoDreamRho : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskFemtoDreamRho();
  AliAnalysisTaskFemtoDreamRho(const char *name, bool isMC, bool doMcTruth, bool doAncestors, bool doCleaning, bool doProjections, float rhoPtThreshold, float pairRapiditiySelection, bool isSameCharge, bool useNegativePairs, bool isMCTrueRhoCombBkrg, bool isMCcheckedCombs);
  virtual ~AliAnalysisTaskFemtoDreamRho();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };
  void SetPosPionCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosPionCuts = trkCuts;
  }
  void SetNegPionCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegPionCuts = trkCuts;
  }
  void SetPosPionMinvCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosPionMinvCuts = trkCuts;
  }
  void SetNegPionMinvCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegPionMinvCuts = trkCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fPosProtonCuts = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *trkCuts)
  {
    fNegProtonCuts = trkCuts;
  }
  void SetRhoCuts(AliFemtoDreamv0Cuts *rhoCuts)
  {
    fRhoCuts = rhoCuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config)
  {
    fConfig = config;
  }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; };
  void SetIsMC(bool isMC) { fIsMC = isMC; };
  void SetDoMcTruth(bool doMcTruth) { fDoMcTruth = doMcTruth; };
  void SetDoCleaning(bool doCleaning) { fDoCleaning = doCleaning; };
  void SetDoAncestors(bool doAncestors) { fDoAncestors = doAncestors; };
  void SetDoProjections(bool doProjector) { fDoProjections = doProjector; };
  void SetRhoPtThreshold(float rhoPtThreshold) { frhoPtThreshold = rhoPtThreshold; };
  void SetMinvSameCharge(bool isSameCharge) { fIsSameCharge = isSameCharge; };
  void SetUseNegativePairs(bool useNegativePairs) { fUseNegativePairs = useNegativePairs; };
  void SetMCTrueRhoCombBkgr(bool isMCTrueRhoCombBkrg) { fIsMCTrueRhoCombBkrg = isMCTrueRhoCombBkrg; };
  void SetMCcheckedCombs(bool isMCcheckedCombs) { fIsMCcheckedCombs = isMCcheckedCombs; };
  void SetPairRapiditySel(float pairRapiditiySelection) { fpairRapiditiySelection = pairRapiditiySelection; };

private:
  AliAnalysisTaskFemtoDreamRho(const AliAnalysisTaskFemtoDreamRho &);
  AliAnalysisTaskFemtoDreamRho &operator=(const AliAnalysisTaskFemtoDreamRho &);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  void SetPhiAtRadiusMCTruth(AliFemtoDreamBasePart partMC, const float bfield);
  bool CommonAncestors(const AliFemtoDreamBasePart &part1, const AliFemtoDreamBasePart &part2, AliAODEvent *Event, bool verbose);
  void FillAncestorHist2D_pTvsMinv(const AliFemtoDreamBasePart &posDaughter,
                                   const AliFemtoDreamBasePart &negDaughter,
                                   TH2F *hist2D);
  void FillAncestorHist2D_PDGvsMinv(const AliFemtoDreamBasePart &posDaughter,
                                    const AliFemtoDreamBasePart &negDaughter,
                                    TH2F *hist2D, int &PDG);
  bool AncestorIsSelected(AliFemtoDreamv0 *v0, AliFemtoDreamv0Cuts *V0Selections);
  bool CommonResonance(const AliFemtoDreamBasePart &part1, const AliFemtoDreamBasePart &part2, int &pdg_resonance, AliAODEvent *Event, bool verbose = false);
  void CalculateAlphaAndQT(const AliAODMCParticle *mcPart, const AliAODMCParticle *mcDaughterOne, const AliAODMCParticle *mcDaughterTwo, float &alpha, float &qT);
  void CalculateAlphaAndQT(const AliFemtoDreamBasePart *Part, const AliFemtoDreamBasePart *DaughterOne, const AliFemtoDreamBasePart *DaughterTwo, float &alpha, float &qT);
  void CalculateAlphaAndQT(const AliFemtoDreamBasePart Part, const TLorentzVector DaughterOne, const TLorentzVector DaughterTwo, float &alpha, float &qT);

  float RelativePairMomentum_check(const AliFemtoDreamBasePart &PartOne,
                                   const int pdg1,
                                   const AliFemtoDreamBasePart &PartTwo,
                                   const int pdg2);
  float RelativePairMomentum_check(TLorentzVector &PartOne,
                                   TLorentzVector &PartTwo);

  float CalcRapidity(const TLorentzVector &part, bool verbose = false, const float fMassForRapidity = 0.7665); // Here use the rho mass
  bool WithinPairRapidityWindow(const AliFemtoDreamBasePart &part1, const int pdg1, const AliFemtoDreamBasePart &part2,
                                const int pdg2, float pairRapiditiySelection, bool verbose = false);

  bool fIsMC;                    //
  bool fDoMcTruth;               //
  bool fDoAncestors;             //
  bool fDoProjections;           //
  float frhoPtThreshold;         //
  bool fIsSameCharge;            //
  bool fUseNegativePairs;        //
  bool fIsMCTrueRhoCombBkrg;     //
  bool fIsMCcheckedCombs;        //
  float fpairRapiditiySelection; //

  TList *fOutput;                           //!
  TList *fHistListRhoCandidatesMCAncestors; //!
  TList *fHistListRhoMCTrue;                //!
  UInt_t fTrigger;                          //
  bool fDoCleaning;                         //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamTrack *fTrackneg;            //!
  AliFemtoDreamv0 *fRhoParticle;            //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fPosPionCuts;     //
  AliFemtoDreamTrackCuts *fNegPionCuts;     //
  AliFemtoDreamTrackCuts *fPosPionMinvCuts; //
  AliFemtoDreamTrackCuts *fNegPionMinvCuts; //
  AliFemtoDreamTrackCuts *fPosProtonCuts;   //
  AliFemtoDreamTrackCuts *fNegProtonCuts;   //

  AliFemtoDreamv0Cuts *fRhoCuts; //

  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliAODTrack **fGTI;                     //!
  int fTrackBufferSize;                   //

  TH2F *fArmenterosRhoTrue;                                               //!
  TH2F *fArmenterosRhoTrue_Reconstr;                                      //!
  TH2F *fArmenterosNoCommonMother_Pos;                                    //!
  TH2F *fArmenterosNoCommonMother_Neg;                                    //!
  TH2F *fArmenterosNoCommonMother_qtDaughBoth;                            //!
  TH2F *fArmenterosNoCommonMother_alphaDaughBoth;                         //!
  TH2F *fArmenterosNoRhoTrue_Reconstr_Pos;                                //!
  TH2F *fArmenterosNoRhoTrue_Reconstr_Neg;                                //!
  TH2F *fArmenterosNoRhoTrue_Reconstr_qtDaughBoth;                        //!
  TH2F *fArmenterosNoRhoTrue_Reconstr_alphaDaughBoth;                     //!
  TH2F *fArmenterosRhoTrue_Reconstr_qtDaughBoth;                          //!
  TH2F *fArmenterosRhoTrue_Reconstr_alphaDaughBoth;                       //!
  TH2F *fHist2D_massVSpt_RhoTrue;                                         //!
  TH1F *fHist1D_pt_RhoTrue;                                               //!
  TH2F *fHist2D_pt1VSpt2_RhoTrue;                                         //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommon;                              //!
  TH2F *fHist2D_massVSpt_RhoCandidateUncommon;                            //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM;                      //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_NoResonances;         //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_kShortResonances;     //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_kStarMisidResonances; //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_rhoResonances;        //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_omegaResonances;      //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_fzeroResonances;      //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_ftwoResonances;       //!
  TH2F *fHist2D_massVSpt_RhoCandidateCommonFullInvM_otherResonances;      //!
  TH2F *fHist2D_massVSpt_RhoCandidateUncommonFullInvM;                    //!
  TH2F *fHist2D_pTvsmT_noPions;                                           //!
  TH2F *fHist2D_pTvsmT_noPrims;                                           //!
  TH2F *fHist2D_pTvsmT_noCommonMother;                                    //!
  TH2F *fHist2D_pTvsmT_noRho;                                             //!
  TH2F *fHist2D_pTvsmT_noRho_MC;                                          //!
  TH2F *fHist2D_PDGvsmT_noRho_MC;                                         //!
  TH2F *fHist2D_pTvsmT_isRho;                                             //!
  TH2F *fHist2D_pTvsmT_isRho_MC;                                          //!
  TH2F *fHist2D_PDGvsMInv_CommonAncestorResonances;                       //!

  ClassDef(AliAnalysisTaskFemtoDreamRho, 11) // Update class number
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskFemtoDreamRho_H_ */