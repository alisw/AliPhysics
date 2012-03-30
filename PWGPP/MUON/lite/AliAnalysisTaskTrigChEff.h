#ifndef ALIANALYSISTASKTRIGCHEFF_H
#define ALIANALYSISTASKTRIGCHEFF_H

/* $Id$ */ 

//
// Class for trigger chamber efficiency calculations
// and tests
//
// Author: Diego Stocco
//

#include "AliVAnalysisMuon.h"

class AliMuonTrackCuts;
class AliVParticle;
class TList;
class TObjArray;
class TString;

class AliAnalysisTaskTrigChEff : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskTrigChEff();
  AliAnalysisTaskTrigChEff(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskTrigChEff();

  void Terminate(Option_t *option);
  void FinishTaskOutput();

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);

  TList* GetEffHistoList(TString physSel, TString trigClassNames, TString centrality, TString trackSelection);

  /// Use ghost tracks in calculations
  void SetUseGhostTracks(Bool_t useGhosts = kTRUE) { fUseGhosts = useGhosts; }

 private:

  AliAnalysisTaskTrigChEff(const AliAnalysisTaskTrigChEff&);
  AliAnalysisTaskTrigChEff& operator=(const AliAnalysisTaskTrigChEff&);

  enum {
    kBendingEff,     ///< Bending plane fired
    kNonBendingEff,  ///< Non-bending plane fired
    kBothPlanesEff,  ///< Both planes fired
    kAllTracks,      ///< tracks used for calculation
    kNcounts         ///< Number of count type
  };

  enum {
    kHchamberEff,    ///< Counts per cathode histogram index
    kHslatEff,       ///< Counts per slat histogram index
    kHboardEff,      ///< Counts per board histogram index
    kHcheckBoard,    ///< Check rejected tracks per board
    kNhistoTypes     ///< Check rejected tracks per board
  };

  enum {
    kNoSelCutApt,   ///< Track matching Apt not passing selection cuts
    kMatchApt,      ///< Match All Pt
    kMatchLpt,      ///< Match Low Pt
    kMatchHpt,      ///< Match High Pt
    kNtrackSel      ///< Total number of selection types
  };
  
  enum {
    kEffFromTrack,  ///< Hit pattern from tracker track extrapolation
    kEffFromTrig,   ///< Hit pattern from trigger
    kNeffMethods    ///< Total number of efficiency methods
  };

  TString GetHistoName(Int_t itype, Int_t icount, Int_t ichamber, Int_t imatch, Int_t imethod);
  Bool_t FillEffHistoList(TString physSel, TString trigClassNames, TString centrality, TString trackSelection, TList* outList = 0x0);
  void InitLocalKeys();
 
  TObjArray* fTrackSelKeys;  ///< Selection names
  TObjArray* fCountTypeKeys; ///< Count type keys
  TObjArray* fHistoTypeKeys; ///< Base histogram name
  TObjArray* fEffMethodKeys; ///< Efficiency methods keys

  Bool_t fUseGhosts; ///< Flag to use also the trigger tracks not matching the tracker in eff. calculation
  TList*  fList;     //!<TList output object

  ClassDef(AliAnalysisTaskTrigChEff, 3); // Trigger chamber efficiencies
};

#endif
