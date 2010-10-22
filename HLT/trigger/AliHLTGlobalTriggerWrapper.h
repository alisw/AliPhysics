//-*- Mode: C++ -*-
// $Id:  $
#ifndef ALIHLTGLOBALTRIGGERWRAPPER_H
#define ALIHLTGLOBALTRIGGERWRAPPER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTriggerWrapper.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Oct 2009
/// @brief  Declaration of the AliHLTGlobalTriggerWrapper interface class.

#include "AliHLTGlobalTrigger.h"
#include "AliHLTLogging.h"
#include "TMethodCall.h"

class TClass;

/**
 * \class AliHLTGlobalTriggerWrapper
 * This class is used to interface with an interpreted class deriving from
 * AliHLTGlobalTrigger. Technically the interpreted class just need to provide
 * the same interface, i.e. the same methods set of methods as AliHLTGlobalTrigger.
 * Direct inheritance from AliHLTGlobalTrigger and AliHLTLogging does not seam to
 * work so providing the same set of methods is enough.
 * The underlying class with which the wrapper interfaces is created and deleted
 * when the wrapper is created and deleted. Thus, the underlying object has the
 * same lifetime as the wrapper object.
 */
class AliHLTGlobalTriggerWrapper : public AliHLTGlobalTrigger, public AliHLTLogging
{
 public:
  
  /**
   * The default constructor constructs the named underlying class and interfaces with it.
   * \param classname  The name of the underlying class to interface with.
   */
  AliHLTGlobalTriggerWrapper(const char* classname);
  
  /// Default destructor will delete the underlying class.
  virtual ~AliHLTGlobalTriggerWrapper();
  
  /**
   * Fill values from a trigger menu. Specifically, the description
   * strings and domain entry values will be copied over.
   * \param  menu  The trigger menu to fill from.
   */
  virtual void FillFromMenu(const AliHLTTriggerMenu& menu);
  
  /**
   * Method to indicate that a new event is being processed and the
   * internal buffers should be reset.
   */
  virtual void NewEvent();
  
  /**
   * This method is used to fill in the internal attributes from the given object
   * which is found in the input data block list.
   * \param  object  The object to fill from.
   * \param  type  The data block type the object was found in.
   * \param  spec  The data block specification the object was found in.
   */
  virtual void Add(
      const TObject* object,
      const AliHLTComponentDataType& type,
      AliHLTUInt32_t spec
    );
  
  /**
   * Calculates the trigger decision and returns the resulting domain and description.
   * \param  triggerResult  The resultant decision bit of the global HLT trigger decision.
   * \param  domain  The resultant trigger domain for the global HLT result.
   * \param  description  The resultant description for the global HLT result.
   * \returns true if any of the triggers in the trigger menu were matched.
   */
  virtual bool CalculateTriggerDecision(bool& triggerResult, AliHLTTriggerDomain& domain, TString& description);
  
  /**
   * Returns the array of internal trigger counters.
   */
  virtual const TArrayL64& GetCounters() const;
  
  /**
   * Sets the internal trigger counter values.
   * \param counters  The array of trigger counters to use.
   */
  virtual void SetCounters(const TArrayL64& counters);
  
  /**
   * Indicates if the CINT interpreter had problems with the interpreted class
   * during a call to one of the interface methods overloaded from AliHLTGlobalTrigger.
   * \returns true if there was a problem with the method call.
   */
  virtual bool CallFailed() const { return fCallFailed; }

  /// Returns true if the wrapper object has setup the underlying object class properly.
  bool IsValid() const;
  
 private:
  
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerWrapper(const AliHLTGlobalTriggerWrapper& obj);
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerWrapper& operator = (const AliHLTGlobalTriggerWrapper& obj);
  
  TClass* fClass; /// Pointer to the object class.
  void* fObject;  /// Pointer to the object being interfaced.
  TMethodCall fFillFromMenuCall;  /// Method call object for FillFromMenu method.
  TMethodCall fNewEventCall;      /// Method call object for NewEvent method.
  TMethodCall fAddCall;           /// Method call object for Add method.
  TMethodCall fCalculateTriggerDecisionCall;  /// Method call object for CalculateTriggerDecision method.
  mutable TMethodCall fGetCountersCall;  /// Method call object for GetCounters method.
  TMethodCall fSetCountersCall;  /// Method call object for SetCounters method.
  mutable bool fCallFailed;  /// Indicates if the last method call failed or not.
  
  ClassDef(AliHLTGlobalTriggerWrapper, 0) // Wrapper class to interface with an interpreted global trigger class.
};

#endif // ALIHLTGLOBALTRIGGERWRAPPER_H
