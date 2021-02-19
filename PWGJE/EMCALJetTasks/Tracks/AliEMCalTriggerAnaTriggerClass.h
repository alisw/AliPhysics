/**
 * \file AliEMCalTriggerAnaTriggerClass.h
 * \brief Definition of a trigger class
 *
 * The definition of a trigger class contains:
 * - a name
 * - a title
 * - a method to select events
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALTRIGGERANATRIGGERCLASS_H
#define ALIEMCALTRIGGERANATRIGGERCLASS_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <exception>
#include <string>
#include <sstream>
#include <TNamed.h>
#include <TObjArray.h>
#include <TString.h>
#include "AliEMCalTriggerAnaHelper.h"

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerAnaTriggerDecision;
class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerAnaPatternObject
 * \brief Helper class describing a trigger pattern in the trigger string
 *
 * A trigger patter has to information:
 * - The pattern string
 * - A switch whether it has to occur or is not supposed to occur in the trigger string
 *   (the later for exclusive classes)
 */
class AliEMCalTriggerAnaPatternObject: public TObject {
public:
  /**
   * Dummy constructor
   */
  AliEMCalTriggerAnaPatternObject():
    TObject(),
    fPattern(""),
    fInString(kTRUE)
  {}
  /**
   * Constructor, defining pattern and whether it has to occur in the triggerstring or is not allowed to occur there
   * \param pattern Pattern to check
   * \param inString Switch defining request for existance in the trigger string
   */
  AliEMCalTriggerAnaPatternObject(const char *pattern, bool inString):
    TObject(),
    fPattern(pattern),
    fInString(inString)
  {}
  /**
   * Destructor, nothing to do
   */
  virtual ~AliEMCalTriggerAnaPatternObject() {}

  Bool_t MatchTriggerString(const char *triggerstring) const;

private:
  TString                 fPattern;             /// Trigger pattern to check for in the trigger string
  Bool_t                  fInString;            /// Definition whether it has to be there or not

  ClassDef(AliEMCalTriggerAnaPatternObject, 1);
};

/**
 * \class AliEMCalTriggerAnaPatternContainer
 * \brief Helper class containing different trigger patterns to check in the trigger string
 */
class AliEMCalTriggerAnaPatternContainer : public TObject{
public:
  AliEMCalTriggerAnaPatternContainer():
    TObject(),
    fPatterns()
  {
    fPatterns.SetOwner(kTRUE);
  }
  /**
   * Detstructor, nothing to do
   */
  virtual ~AliEMCalTriggerAnaPatternContainer() {}

  /**
   *
   * \param pattern
   * \param isRequested
   */
  void AddPattern(const char *pattern, bool isRequested){
    fPatterns.Add(new AliEMCalTriggerAnaPatternObject(pattern, isRequested));
  }

  bool CheckTriggerString(const char * triggerstring) const;

private:
  TObjArray               fPatterns;      /// Container for patterns

  ClassDef(AliEMCalTriggerAnaPatternContainer, 1);
};

/**
 * \class AliEMCalTriggerAnaTriggerPatchTypeObject
 * \brief Wrapper class for primitive type ETATriggerType
 */
class AliEMCalTriggerAnaTriggerPatchTypeObject : public TObject {
public:
  /**
   * Dummy constructor
   */
  AliEMCalTriggerAnaTriggerPatchTypeObject():
    TObject(),
    fTriggerType(kTAUndef)
  {}
  /**
   * Main constructor
   * \param triggertype The primitive trigger type
   */
  AliEMCalTriggerAnaTriggerPatchTypeObject(ETATriggerType triggertype);
  /**
   * Destructor, nothing to do
   */
  virtual ~AliEMCalTriggerAnaTriggerPatchTypeObject() {}

  /**
   * Get the underlying trigger type
   * \return The underlying trigger type
   */
  ETATriggerType GetTriggerType() const { return fTriggerType; }

private:
  ETATriggerType                fTriggerType;       ///< Underlying trigger type

  ClassDef(AliEMCalTriggerAnaTriggerPatchTypeObject, 1);
};

/*
 * Exceptions:
 */

/**
 * \class TriggerMethodUndefinedException
 * \brief Exception class for events trigger classes which do not have a method for event selection defined
 */
class TriggerMethodUndefinedException : public std::exception{
public:
  /**
   * Constructor, creating error message
   * \param triggerclass Trigger class throwing the exception
   */
  TriggerMethodUndefinedException(std::string triggerclass):
    fMessage("")
  {
    fMessage = "Trigger method undefined for trigger class " + triggerclass;
  }
  /**
   * Destructor
   */
  virtual ~TriggerMethodUndefinedException() throw () { }

  /**
   * Returns the error message
   * \return The error message
   */
  const char *what() const throw () {
    return fMessage.c_str();
  }
private:
  std::string                     fMessage;     ///< The error message
};

/**
 * \class PatchHandlerMissingException
 * \brief Exception class for events where the trigger patch handler is not set
 */
class PatchHandlerMissingException : public std::exception{
public:
  /**
   * Constructor, creating error message
   * \param triggerclass Trigger class throwing the exception
   */
  PatchHandlerMissingException(std::string triggerclass):
    fMessage("")
  {
    fMessage = "Trigger method undefined for trigger class " + triggerclass;
  }
  /**
   * Destructor
   */
  virtual ~PatchHandlerMissingException() throw () { }

  /**
   * Returns the error message
   * \return The error message
   */
  const char *what() const throw () {
    return fMessage.c_str();
  }
private:
  std::string                     fMessage;     ///< The error message
};

/**
 * \class EventCorruptionException
 * \brief Exception class thrown when event collection is corrupted (event pointer 0 or trigger patch container not found)
 */
class EventCorruptionException : public std::exception{
public:
  /**
   * Constructor, building the error message
   * \param triggerclass Trigger class throwing the exception
   * \param reason Reason for throwing the exception
   */
  EventCorruptionException(std::string triggerclass, std::string reason):
    fMessage("")
  {
      std::stringstream messagebuilder;
      messagebuilder << "Event corrupted for trigger class " << triggerclass << ": " << reason;
      fMessage = messagebuilder.str();
  }
  /**
   * Destructor
   */
  virtual ~EventCorruptionException() throw () { }

  /**
   * Returns the error message
   * \return The error message
   */
  const char *what() const throw () {
    return fMessage.c_str();
  }

private:
  std::string                     fMessage;     ///< The error message
};

/**
 * \class AliEMCalTriggerAnaTriggerClass
 * \brief Description of a trigger class
 *
 * A trigger class is a class of events trigger under the same conditions. Triggers can be identified
 * - via the occurrence of a given trigger bit
 * - via the occurrence of a given pattern (can be complex) in the trigger string
 * - via the occurrence of a given patch type reconstructed in the EMCAL
 *
 * Consequently a trigger class has
 * - a name and a title
 * - a method how to select events
 */
class AliEMCalTriggerAnaTriggerClass: public TNamed {
public:
  AliEMCalTriggerAnaTriggerClass();
  AliEMCalTriggerAnaTriggerClass(const char *name, const char *title);
  virtual ~AliEMCalTriggerAnaTriggerClass();

  bool IsEventTriggered(const AliEMCalTriggerEventData * const) const;
  /**
   * Check whether trigger class is marked as a min. bias trigger
   * \return True if the trigger is a min. bias trigger, false otherwise
   */
  bool IsMinBiasTrigger() const { return fIsMinBiasTrigger; }

  /**
   * Add trigger bit to the class. Set the request for a trigger bit to true
   * \param triggerbit Trigger bit to check
   */
  void AddTriggerBit(UInt_t triggerbit) {
    fTriggerBits |= triggerbit;
    fDecisionFromTriggerBits = kTRUE;
  }
  /**
   * Add trigger pattern to the trigger class
   * \param triggerpattern Trigger pattern to check in the trigger string
   * \param isRequested If true the pattern has to appear in the trigger string, otherwise it is not supposed to appear there
   */
  void AddTriggerStringPattern(const char *triggerpattern, bool isRequested) {
    fTriggerStringPattern.AddPattern(triggerpattern, isRequested);
    fDecisionFromTriggerString = kTRUE;
  }
  /**
   * Add trigger patch type for selection based on trigger patchess
   * \param triggertype the patch type to select
   */
  void AddTriggerPatchType(ETATriggerType triggertype){
    fTriggerPatchTypes.Add(new AliEMCalTriggerAnaTriggerPatchTypeObject(triggertype));
  }
  /**
   * Set handling for trigger patches.
   * \param triggerhandler The trigger patch decision
   */
  void SetTriggerDecisionHandler(AliEMCalTriggerAnaTriggerDecision *triggerhandler) { fEmcalTriggerHandler = triggerhandler; }
  /**
   * Flag event as min. bias event
   * \param isMinBias If true event is flagged as a min. bias event
   */
  void SetMinBiasTrigger(Bool_t isMinBias) { fIsMinBiasTrigger = isMinBias; }

protected:
  Bool_t                                  fDecisionFromTriggerBits;             ///< Switch for using trigger bits for event selection
  Bool_t                                  fDecisionFromTriggerString;           ///< Switch for using the trigger string for event selection
  Bool_t                                  fDecisionFromTriggerPatches;          ///< Switch for using reconstructed trigger patches for event selection

  Bool_t                                  fIsMinBiasTrigger;                    ///< Flag class as min. bias trigger class
  UInt_t                                  fTriggerBits;                         ///< Trigger bits used for event selection
  AliEMCalTriggerAnaPatternContainer      fTriggerStringPattern;                ///< Trigger patterns used for event selection
  TObjArray                               fTriggerPatchTypes;                   ///< Trigger patch types used for the selection
  AliEMCalTriggerAnaTriggerDecision       *fEmcalTriggerHandler;                ///< Handling of trigger patch selection

  ClassDef(AliEMCalTriggerAnaTriggerClass, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERANATRIGGERCLASS_H */
