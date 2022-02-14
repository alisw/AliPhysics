/**
 * \file AliEMCalTriggerTracksAnalysisComponent.h
 * \brief Base class for analysis components
 *
 * Base class for analysis components. Inheriting classes have to implement the
 * functions CreateHistos and Process.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <exception>
#include <map>
#include <vector>
#include <string>
#include <TNamed.h>
#include <THistManager.h>

namespace PWGJE{ 
  
namespace EMCALJetTasks{

class AliEMCalTriggerAnaClassManager;
class AliEMCalTriggerBinningComponent;
class AliEMCalTriggerBinningDimension;
class AliEMCalTriggerEventData;
class AliEMCalTriggerKineCuts;
class AliEMCalTriggerWeightHandler;

/**
 * \class TriggerManagerNotFoundException
 * \brief Exception class for events with missing trigger configuration handler
 */
class TriggerManagerNotFoundException : public std::exception{
public:
  /**
   * Dummy constructor
   */
  TriggerManagerNotFoundException():
    fMessage("")
  {
    fMessage =  "Trigger handler not found";
  }
  /**
   * Main constructor, to be called when the exception is thrown
   * \param producer Object producing the exception
   */
  TriggerManagerNotFoundException(std::string producer):
    fMessage("")
  {
    fMessage = "Trigger handler not found for object " + producer;
  }
  /**
   * Destructor, nothing to do
   */
  virtual ~TriggerManagerNotFoundException() throw() {}

  const char *what() const throw() {
    return fMessage.c_str();
  }

private:
  std::string                     fMessage;
};

/**
 * \class AliEMCalTriggerTracksAnalysisComponent
 * \brief Base class for analysis components in the analysis of EMCAL-triggered events
 *
 * This class defines the base class of all analysis components used in the analysis of
 * EMCAL-triggered events. A new analysis component has to implement at least the functions
 *  - CreateHistos
 *  - Process
 * where the function process is abstract. The function CreateHistos contains the code for
 * the initialization of the component and is run at the start of the analysis. Process
 * implements the event loop and is called once per event.
 */
class AliEMCalTriggerTracksAnalysisComponent : public TNamed {
public:
  AliEMCalTriggerTracksAnalysisComponent();
  AliEMCalTriggerTracksAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerTracksAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data) = 0;

  /**
   * Get the list of histograms of this analysis component.
   * \return The list of histograms
   */
  THashList *GetHistList() const { return fHistos->GetListOfHistograms(); }

  /**
   * Get the common weight handler.
   * \return The weight handler
   */
  const AliEMCalTriggerWeightHandler *GetWeightHandler() const { return fWeightHandler; }

  /**
   * Set the global binning handler to this analysis component.
   * \param binning The global binning handler
   */
  void SetBinning(const AliEMCalTriggerBinningComponent * const binning) { fBinning = binning; }

  /**
   * Set the global kinematical cuts to this analysis components.
   * \param cuts The global kinematic cuts
   */
  void SetKineCuts(const AliEMCalTriggerKineCuts * const cuts) { fKineCuts = cuts; }

  /**
   * Set the global trigger class manager
   * \param classmgr The global trigger class manager
   */
  void SetTriggerClassManager(const AliEMCalTriggerAnaClassManager *classmgr) { fTriggerClassManager = classmgr; }

  /**
   * Set the global weight handler to this analysis component.
   * \param handler The global weight handler.
   */
  void SetWeightHandler(const AliEMCalTriggerWeightHandler *handler) { fWeightHandler = handler; }

  /**
   * Set the debug level for a given analysis component.
   * \param debuglevel The component debug level
   */
  void SetComponentDebugLevel(int debuglevel) { fComponentDebugLevel = debuglevel; }

protected:
  TAxis *DefineAxis(const char *name, const TBinning &binning);
  void GetMachingTriggerNames(std::vector<std::string> &triggernames) const;
  void GetAllTriggerNamesAndTitles(std::map<std::string, std::string> &triggers) const;
  void PrintTriggerNames(const std::vector<std::string> &, const std::string &componentName) const;

  THistManager                                *fHistos;                   ///< Histogram container of the analysis component
  const AliEMCalTriggerAnaClassManager        *fTriggerClassManager;      ///< Global trigger class manager
  const AliEMCalTriggerBinningComponent       *fBinning;                  ///< Global binning handler
  const AliEMCalTriggerKineCuts               *fKineCuts;                 ///< Kinematical cuts for tracks and particle selection
  const AliEMCalTriggerWeightHandler          *fWeightHandler;            ///< Event weight handler

  Int_t                                         fComponentDebugLevel;     ///< Debug level for the given analysis component

  ClassDef(AliEMCalTriggerTracksAnalysisComponent, 1)
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H */
