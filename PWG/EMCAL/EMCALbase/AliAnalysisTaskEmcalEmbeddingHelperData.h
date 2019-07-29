#ifndef ALIANALYSISTASKEMCALEMBEDDINGHELPERDATA_H
#define ALIANALYSISTASKEMCALEMBEDDINGHELPERDATA_H

#include <vector>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliTrackContainer.h"
#include "THistManager.h"

/**
 * \file AliAnalysisTaskEmcalEmbeddingHelperData.h
 * \brief Declaration of class AliAnalysisTaskEmcalEmbeddingHelperData
 *
 * Class inheriting from AliAnalysisTaskEmcalEmbeddingHelper which adds 
 * more functionality to embedding helper, specifically to be used
 * when embedding data to select events with tracks in a given pT range
 */

class AliAnalysisTaskEmcalEmbeddingHelperData : public AliAnalysisTaskEmcalEmbeddingHelper{
 public:

  AliAnalysisTaskEmcalEmbeddingHelperData()                          ;
  AliAnalysisTaskEmcalEmbeddingHelperData(const char *name)          ;
  virtual ~AliAnalysisTaskEmcalEmbeddingHelperData()                 ;

  void            ExecOnce() ;
  void            UserCreateOutputObjects() ;

  void            SetEmbeddedTrackContainer(AliTrackContainer *cont) {fTrackContainer = cont;}
  AliTrackContainer *GetEmbeddedTrackContainer() {return fTrackContainer;}

  void            AddPtSelection(Double_t min, Double_t max) { fPtSelection.push_back(std::make_pair(min,max));}

  // add task
  static AliAnalysisTaskEmcalEmbeddingHelperData *AddTaskEmcalEmbeddingHelperData();

  protected:
  void            RetrieveTaskPropertiesFromYAMLConfig();
  Bool_t          CheckIsEmbeddedEventSelected();

  std::vector<std::pair<Double_t,Double_t>>       fPtSelection; ///< min and max pT stored as a vector of pairs
  AliTrackContainer *fTrackContainer; ///< track container for pT selection

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalEmbeddingHelperData, 1);
  /// \endcond

};

#endif
