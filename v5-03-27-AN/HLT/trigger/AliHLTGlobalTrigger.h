//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTGLOBALTRIGGER_H
#define ALIHLTGLOBALTRIGGER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTrigger.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Declaration of the AliHLTGlobalTrigger base class.

#include "TObject.h"
#include "AliHLTDataTypes.h"

class AliHLTTriggerDomain;
class AliHLTTriggerMenu;
class TString;
class TArrayL64;

/**
 * \class AliHLTGlobalTrigger
 * This class is an abstract class. Classes which derive from this class should
 * implement the logic for a particular trigger menu. The AliHLTGlobalTriggerComponent
 * takes the AliHLTTriggerMenu class and creates a class deriving from AliHLTGlobalTrigger
 * on the fly to implement the trigger logic for that particular trigger menu.
 */
class AliHLTGlobalTrigger
{
 public:
  
  /// Default destructor.
  virtual ~AliHLTGlobalTrigger() {}
  
  /**
   * Abstract method to fill values from a trigger menu. Specifically, the description
   * strings and domain entry values will be copied over.
   * \param  menu  The trigger menu to fill from.
   */
  virtual void FillFromMenu(const AliHLTTriggerMenu& menu) = 0;
  
  /**
   * Abstract method to indicate that a new event is being processed and the
   * internal buffers should be cleared or reset.
   */
  virtual void NewEvent() = 0;
  
  /**
   * Abstract method which should fill in the internal attributes from the given
   * object.
   * \param  object  The object to fill from.
   * \param  type  The data block type the object was found in.
   * \param  spec  The data block specification the object was found in.
   */
  virtual void Add(
      const TObject* object,
      const AliHLTComponentDataType& type,
      AliHLTUInt32_t spec
    ) = 0;
  
  /**
   * Abstract method that calculates the trigger decision
   * \param  triggerResult  The resultant decision bit of the global HLT trigger decision.
   * \param  domain  The resultant trigger domain for the global HLT result.
   * \param  description  The resultant description for the global HLT result.
   * \returns true if any of the triggers in the trigger menu were matched.
   */
  virtual bool CalculateTriggerDecision(bool& triggerResult, AliHLTTriggerDomain& domain, TString& description) = 0;
  
  /**
   * Returns the array of internal trigger counters.
   */
  virtual const TArrayL64& GetCounters() const = 0;
  
  /**
   * Sets the internal trigger counter values.
   * \param counters  The array of trigger counters to use.
   */
  virtual void SetCounters(const TArrayL64& counters) = 0;
  
  /**
   * Method for checking if the last call to one of the AliHLTGlobalTrigger
   * methods failed. This is used by the AliHLTGlobalTriggerWrapper to indicate
   * if the CINT interpreter had problems with the interpreted class.
   * \returns true if there was a problem with the method call.
   */
  virtual bool CallFailed() const { return false; }
  
  /**
   * Creates a new instance of a particular global trigger class.
   * \param name  The name of the class to create.
   * \returns the new trigger class instance which needs to be deleted by the
   *    caller with the delete operator.
   */
  static AliHLTGlobalTrigger* CreateNew(const char* name);
  
  ClassDef(AliHLTGlobalTrigger, 0) // Global HLT trigger base class which implements logic for a particular trigger menu.
};

#endif // ALIHLTGLOBALTRIGGER_H
