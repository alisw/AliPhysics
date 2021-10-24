/**
 * \file AliEMCalTriggerAnaClassManager.h
 * \brief Declaration of a management class for trigger classes
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGERANACLASSMANAGER_H
#define ALIEMCALTRIGGERANACLASSMANAGER_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <exception>
#include <TNamed.h>

class TObjArray;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerAnaTriggerClass;
class AliEMCalTriggerAnaTriggerDecision;
class AliEMCalTriggerEventData;

/**
 * \class TriggerManagerEmptyException
 * \brief Exception handling in case event selection is performed on an empty trigger manager
 */
class TriggerManagerEmptyException : public std::exception {
public:
  /**
   * Constructor
   */
  TriggerManagerEmptyException() {}
  /**
   * Destructor, nothing to do
   */
  virtual ~TriggerManagerEmptyException() throw () {}

  /**
   * Return error message
   * \return The error message
   */
  virtual const char *what() const throw() { return "Trigger manager does not contain any trigger class"; }

};

/**
 * \class AliEMCalTriggerAnaClassManager
 * \brief Manager for trigger classes
 *
 * This class manages trigger classes, meaning it serves as a
 * container and steers the event selection.
 */
class AliEMCalTriggerAnaClassManager: public TNamed {
public:
  AliEMCalTriggerAnaClassManager();
  AliEMCalTriggerAnaClassManager(const char *name);
  AliEMCalTriggerAnaClassManager(const AliEMCalTriggerAnaClassManager &ref);
  AliEMCalTriggerAnaClassManager &operator=(const AliEMCalTriggerAnaClassManager &ref);
  virtual ~AliEMCalTriggerAnaClassManager();

  void PerformEventSelection(AliEMCalTriggerEventData *trgevent);

  void AddTriggerClass(AliEMCalTriggerAnaTriggerClass *triggerclass);
  void SetTriggerDecision(AliEMCalTriggerAnaTriggerDecision *triggerdecision);

  TObjArray * GetSelectedTriggerClasses() const;
  TObjArray * GetAllTriggerClasses() const;

  bool HasMinBiasTrigger() const;

private:
  TObjArray         *fTriggerClasses;               ///< List of trigger classes
  TObjArray         *fSelected;                     ///< List of selected trigger classes

  ClassDef(AliEMCalTriggerAnaClassManager, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERANACLASSMANAGER_H */
