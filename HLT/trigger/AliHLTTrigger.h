#ifndef ALIHLTTRIGGER_H
#define ALIHLTTRIGGER_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

#include "AliHLTProcessor.h"

class AliHLTTrigger : public AliHLTProcessor
{
 public:
 
  AliHLTTrigger();
  virtual ~AliHLTTrigger();

  /**
   * Returns the name of the trigger. This must be unique across the system.
   * This method is pure virtual and must be implemented by child classes.
   * @return string containing the trigger name.
   */
  virtual const char* GetTriggerName() const = 0;
  
  /**
   * Inherited from AliHLTComponent. Returns the name of the trigger by default.
   * @return string containing the trigger name as the component ID.
   */
  virtual const char* GetComponentID() { return GetTriggerName(); };

  /**
   * Get the input data types of the component.
   * This method returns kAliHLTAnyDataType by default.
   * @param list <i>[out]</i>: The list of data types to be filled.
   */
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list) const
  {
    list.push_back(kAliHLTAnyDataType);
  }

  /**
   * Returns extra output data types this trigger generates.
   * This returns an empty list by default.
   * @param list <i>[out]</i>: The list of data types to be filled.
   */
  virtual void GetOutputDataTypes(AliHLTComponentDataTypeList& list) const
  {
    list.push_back(kAliHLTDataTypeTObject);
  }

  /**
   * Get a ratio by how much the data volume is shrinked or enhanced.
   * The method returns a size proporional to the trigger name string length
   * for constBase, and 1 for inputMultiplier.
   * @param constBase        <i>[out]</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>[out]</i>: multiplication ratio
   */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
  {
    constBase = strlen(GetTriggerName()) + sizeof(TObjString) + 1;
    inputMultiplier = 1;
  }

 protected:

  /**
   * This method needs to be implemented by child classes to implement the actual
   * trigger algorithm. A possitive trigger decision is made by calling the TriggerEvent
   * method with TriggerEvent(true), or TriggerEvent(false) for a negative result
   * (no trigger).
   * If the AliHLTComponentEventData structure is needed for the current event being
   * processed then the GetEventData method can be used to fetch it.
   * Similarly, the AliHLTComponentTriggerData structure can be fetched with a call
   * to GetTriggerData.
   * @return Zero should be returned on success and a negative error code, which is
   *     compatible with the AliHLTSystem framework, on failure.
   */
  virtual int DoTrigger() = 0;
  
  /**
   * Fills the output with the trigger decision. This should be called in the DoTrigger
   * method when a trigger decision has been made.
   * @param value  The trigger decision value. True for positive trigger and false
   *     for a negative result. (true by default)
   */
  void TriggerEvent(bool value = true);
  
  /**
   * Returns the event data structure for the current event.
   * NULL is returned if this method is not called inside the DoTrigger method.
   */
  const AliHLTComponentEventData* GetEventData() const { return fEventData; }
  
  /**
   * Returns the trigger data structure for the current event.
   * NULL is returned if this method is not called inside the DoTrigger method.
   */
  AliHLTComponentTriggerData* GetTriggerData() const { return fTriggerData; }

 private:
 
  /**
   * Inherited from AliHLTComponent. This method will clear the fDecisionMade flag,
   * remember the evtData and trigData pointers, then call DoTrigger to invoke the
   * actual trigger algorithm.
   */
  virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
  
  using AliHLTProcessor::DoEvent;
 
  /**
   * Inherited from AliHLTComponent. This method is deprecated so hide it.
   * @return kAliHLTMultipleDataType is returned.
   */
  virtual AliHLTComponentDataType GetOutputDataType() { return kAliHLTMultipleDataType; };

  /**
   * Inherited from AliHLTComponent. This method is replaced with one that is
   * symmetric to GetInputDataTypes that returns void, so we make this method
   * private.
   * @param list     list to receive the output data types.
   * @return the number of elements in the list.
   */
  virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list)
  {
    GetOutputDataTypes(list);
    return list.size();
  }
  
  const AliHLTComponentEventData* fEventData; ///! Event data for the current event. Only valid inside DoTrigger.
  AliHLTComponentTriggerData* fTriggerData; ///! Trigger data for the current event. Only valid inside DoTrigger.
  bool fDecisionMade;  ///! Flag indicating if the trigger decision has been made for this trigger yet.
  int fTriggerEventResult;  ///! Result returned by PushBack method in the TriggerEvent method.
  
  ClassDef(AliHLTTrigger, 0) // Base class for HLT triggers.

};

#endif // ALIHLTTRIGGER_H

