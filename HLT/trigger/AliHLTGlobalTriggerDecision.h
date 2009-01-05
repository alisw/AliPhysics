#ifndef ALIHLTGLOBALTRIGGERDECISION_H
#define ALIHLTGLOBALTRIGGERDECISION_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTriggerDecision.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Declaration of the AliHLTGlobalTriggerDecision class storing the global HLT decision.

#include "AliHLTTriggerDecision.h"
#include "TArrayL64.h"
#include "TObjArray.h"

class AliHLTGlobalTriggerDecision : public AliHLTTriggerDecision
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTGlobalTriggerDecision();
  
  /**
   * Constructor specifying multiple information fields.
   * \param result  The result of the global trigger decision.
   * \param triggerDomain  The trigger domain for the global trigger decision.
   * \param description  The description of (reason for) the global trigger decision.
   */
  AliHLTGlobalTriggerDecision(
      bool result, const AliHLTTriggerDomain& triggerDomain,
      const char* description = ""
    );
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTGlobalTriggerDecision();
  
  /**
   * Inherited from TObject, this prints the contents of the trigger decision.
   * \param option  Can be "short" which will print the short format or "counters"
   *    which will print only the counters or "compact" which will print only the
   *    global information but not the lists of input objects.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * Returns the number of trigger inputs that contributed to this global trigger decision.
   */
  Int_t NumberOfTriggerInputs() const { return fContributingTriggers.GetEntriesFast(); }
  
  /**
   * Returns the i'th trigger input object in fContributingTriggers.
   */
  const AliHLTTriggerDecision* TriggerInput(Int_t i) const
  {
    return static_cast<const AliHLTTriggerDecision*>( fContributingTriggers[i] );
  }
  
  /**
   * Returns the list of trigger inputs used when making the global HLT trigger decision.
   */
  const TClonesArray& TriggerInputs() const { return fContributingTriggers; }
  
  /**
   * Adds a trigger input to the list of triggers that were considered when making
   * this global trigger decision.
   * \param decision  The trigger decision object to add.
   */
  void AddTriggerInput(const AliHLTTriggerDecision& decision)
  {
    new (fContributingTriggers[fContributingTriggers.GetEntriesFast()]) AliHLTTriggerDecision(decision);
  }
  
  /**
   * Returns the number of other input objects that contributed to this global trigger decision.
   */
  Int_t NumberOfInputObjects() const { return fInputObjects.GetEntriesFast(); }
  
  /**
   * Returns the i'th input object in fInputObjects.
   */
  const TObject* InputObject(Int_t i) const { return fInputObjects[i]; }
  
  /**
   * Returns the list of other input objects used when making the global HLT trigger decision.
   */
  const TObjArray& InputObjects() const { return fInputObjects; }
  
  /**
   * Adds a input object to the list of input objects that were considered when
   * making this global trigger decision.
   * \param object  The input object to add.
   * \note  A copy of the object is made with TObject::Clone() and added.
   */
  void AddInputObject(const TObject* object)
  {
    fInputObjects.Add(object->Clone());
  }
  
  /**
   * Sets the counter array.
   */
  void SetCounters(const TArrayL64& counters) { fCounters = counters; }
  
  /**
   * Returns the event trigger counters associated with the global trigger classes.
   */
  const TArrayL64& Counters() const { return fCounters; }
  
 private:
  
  TClonesArray fContributingTriggers;  /// The list of contributing trigger decisions from all AliHLTTrigger components that were considered.
  TObjArray fInputObjects;  /// The list of other input objects.
  TArrayL64 fCounters;  /// Event trigger counters. One counter for each trigger class in the global trigger.
  
  ClassDef(AliHLTGlobalTriggerDecision, 1) // Contains the HLT global trigger decision and information contributing to the decision.
};

#endif // ALIHLTGLOBALTRIGGERDECISION_H

