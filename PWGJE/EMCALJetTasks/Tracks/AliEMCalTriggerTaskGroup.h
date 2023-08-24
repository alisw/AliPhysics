/**
 * \file AliEMCalTriggerTaskGroup.h
 * \brief Container class for Analysis components with a common event selection
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERTASKGROUP_H
#define ALIEMCALTRIGGERTASKGROUP_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TNamed.h>
#include <TObjArray.h>

namespace PWGJE{
  
namespace EMCALJetTasks {

class AliEMCalTriggerAnaClassManager;
class AliEMCalTriggerBinningComponent;
class AliEMCalTriggerEventSelection;
class AliEMCalTriggerKineCuts;
class AliEMCalTriggerTracksAnalysisComponent;
class AliEMCalTriggerWeightHandler;

/**
 * \class AliEMCalTriggerTaskGroup
 * \brief Container class for analysis components with common event selection
 *
 * Group of analysis components with the same event selection, kinematical cuts, and histogram binnings.
 * Analysis components are initialised via the Initialise function, and executed, if
 * the event is selected, via the function process.
 */
class AliEMCalTriggerTaskGroup : public TNamed {
public:
  AliEMCalTriggerTaskGroup();
  AliEMCalTriggerTaskGroup(const char *name);
  AliEMCalTriggerTaskGroup(const AliEMCalTriggerTaskGroup &ref);
  AliEMCalTriggerTaskGroup &operator=(const AliEMCalTriggerTaskGroup &ref);
  virtual ~AliEMCalTriggerTaskGroup();
  virtual void Copy(TObject &other) const;

  /**
   * Set global event selection to the task group
   * \param sel The event selection handler
   */
  void SetEventSelection(const AliEMCalTriggerEventSelection *sel){ fEventSelection = sel; }

  /**
   * Set global binning handler to the task group
   * \param binning The binning handler
   */
  void SetGlobalBinning(const AliEMCalTriggerBinningComponent *const binning) { fBinning = binning; }

  void SetWeightHandler(const AliEMCalTriggerWeightHandler *handler);

  /**
   * Set common kinematical cuts for the particle selection within a group
   * \param cuts The kinematical cuts
   */
  void SetKineCuts(const AliEMCalTriggerKineCuts *cuts) { fKineCuts = cuts; }

  void AddAnalysisComponent(AliEMCalTriggerTracksAnalysisComponent * const analysis);

  TList * InitialiseAnalysisComponents(const AliEMCalTriggerAnaClassManager *mgr = NULL);
  void Process(const AliEMCalTriggerEventData * const event);

protected:
  TObjArray                                   *fAnalysisComponents;   ///< List of analysis components connected to the group
  const AliEMCalTriggerEventSelection         *fEventSelection;       ///< Common event selection for the group
  const AliEMCalTriggerBinningComponent       *fBinning;              ///< Binning handler
  const AliEMCalTriggerKineCuts               *fKineCuts;             ///< Kinematical cuts shared by the task group
  const AliEMCalTriggerWeightHandler          *fWeightHandler;        ///< Weight handler for event weighting

  ClassDef(AliEMCalTriggerTaskGroup, 1);    // Group of analysis components with common event selection
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERTASKGROUP_H */
