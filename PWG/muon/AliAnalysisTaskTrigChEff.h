#ifndef ALIANALYSISTASKTRIGCHEFF_H
#define ALIANALYSISTASKTRIGCHEFF_H

/* $Id$ */ 

/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskTrigChEff
/// \brief Analysis task for trigger chamber efficiency determination
///
//  Author Diego Stocco

#include "AliAnalysisTaskSE.h"

class TList;

class AliAnalysisTaskTrigChEff : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTrigChEff();
  AliAnalysisTaskTrigChEff(const char *name);
  virtual ~AliAnalysisTaskTrigChEff();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  /// Use ghost tracks in calculations
  void SetUseGhostTracks(Bool_t useGhosts = kTRUE) { fUseGhosts = useGhosts; }

protected:
  void ResetHistos();

private:
  /// Not implemented
  AliAnalysisTaskTrigChEff(const AliAnalysisTaskTrigChEff& rhs);
  /// Not implemented
  AliAnalysisTaskTrigChEff& operator = (const AliAnalysisTaskTrigChEff& rhs);
    
  Bool_t fUseGhosts; ///< Flag to use also the trigger tracks not matching the tracker in eff. calculation

  TList*  fList; //!<TList output object

  enum {
    kNcathodes = 2,  ///< Number of cathodes
    kNchambers = 4,  ///< Number of chambers
    kNslats    = 18 ///< Number of slats
  };

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
    kHcheckBoard    ///< Check rejected tracks per board
  };
  
  Int_t GetHistoIndex(Int_t histoType, Int_t countType=-1, 
		      Int_t chamber=-1);

  ClassDef(AliAnalysisTaskTrigChEff, 1); // Trigger chamber efficiency analysis
};

#endif

