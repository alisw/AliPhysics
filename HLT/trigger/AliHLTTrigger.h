#ifndef ALIHLTTRIGGER_H
#define ALIHLTTRIGGER_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTTrigger.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   12 Aug 2008
/// @brief  Declaration of the AliHLTTrigger base component class.

#include "AliHLTProcessor.h"
#include "AliHLTReadoutList.h"
#include "AliHLTTriggerDomain.h"

class AliHLTTriggerDecision;

/**
 * \class AliHLTTrigger
 * This is the base class from which all HLT trigger components should inherit.
 */
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
   * \note The underlying non const version of GetOutputDataTypes adds the value
   *    kAliHLTDataTypeTObject to the list.
   */
  virtual void GetOutputDataTypes(AliHLTComponentDataTypeList& /*list*/) const {}

  /**
   * Get a ratio by how much the data volume is shrunk or enhanced.
   * The method returns a size proportional to the trigger name string length
   * for constBase, and 1 for inputMultiplier.
   * @param constBase        <i>[out]</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>[out]</i>: multiplication ratio
   */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

 protected:

  /// Not implemented.
  AliHLTTrigger(const AliHLTTrigger& obj);
  /// Not implemented.
  AliHLTTrigger& operator = (const AliHLTTrigger& obj);

  /**
   * This method needs to be implemented by child classes to implement the actual
   * trigger algorithm. A positive trigger decision is made by calling the TriggerEvent
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
   * \returns zero on success and negative value on failure. The possible failure
   *    codes are:<br>
   *     -ENOSPC - If there is not enough output buffer space for the trigger decision.<br>
   *     -ENOMSG - If the trigger decision object could not be serialised.
   */
  int TriggerEvent(bool value = true);
  
  /**
   * Fills the output with the given trigger decision. This should be called in the
   * DoTrigger method when a custom trigger decision has been constructed.
   * @param value  The custom trigger decision object.
   * @param datatype  The data block type to use (set to
   *    kAliHLTDataTypeTObject|kAliHLTDataOriginOut by default).
   * @param spec  The data block specification to use (set to kAliHLTVoidDataSpec
   *    by default).
   * \returns zero on success and negative value on failure. The possible failure
   *    codes are:<br>
   *     -ENOSPC - If there is not enough output buffer space for the trigger decision.<br>
   *     -ENOMSG - If the trigger decision object could not be serialised.<br>
   *     -EINVAL - If the <i>result</i> object is NULL.
   */
  int TriggerEvent(
      AliHLTTriggerDecision* result,
      const AliHLTComponentDataType& type = kAliHLTDataTypeTObject|kAliHLTDataOriginOut,
      AliHLTUInt32_t spec = kAliHLTVoidDataSpec
    );
  
  /**
   * Method for finding out the result of the last call to TriggerEvent.
   * \returns the error code for the last call to TriggerEvent.
   */
  int GetLastTriggerEventResult() const { return fTriggerEventResult; }
  
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

  /**
   * Set a bit to 1 in the readout list which will enable that DDL for readout
   * @param ddlId     Equipment ID of DDL to readout, in decimal.
   */
  void EnableDDLBit(Int_t ddlId)
  {
    AliHLTComponent::SetDDLBit(fReadoutList, ddlId, kTRUE);
  }

  /**
   * Set a bit to 0 in the readout list which will exclude that DDL from the readout.
   * @param ddlId     Equipment ID of DDL not to readout, in decimal.
   */
  void DisableDDLBit(Int_t ddlId)
  {
    AliHLTComponent::SetDDLBit(fReadoutList, ddlId, kFALSE);
  }
  
  /**
   * Set or unset bit in the readout list.
   * @param ddlId     Equipment ID of DDL to set, in decimal.
   * @param state     kTRUE will enable readout of that DDL, kFALSE will disable readout.
   */
  void SetDDLBit(Int_t ddlId, Bool_t state)
  {
    AliHLTComponent::SetDDLBit(fReadoutList, ddlId, state);
  }
  
  /**
   * Returns the DDL readout list.
   */
  const AliHLTReadoutList& GetReadoutList() const { return fReadoutList; }

  /**
   * Returns the DDL readout list for modification by hand.
   */
  AliHLTReadoutList& GetReadoutList() { return fReadoutList; }

  /**
   * Sets the readout list object.
   * \param value  The new value to use for the readout list.
   */
  void SetReadoutList(const AliHLTReadoutList& value) { fReadoutList = value; }
  
  /**
   * Returns the trigger domain object.
   */
  const AliHLTTriggerDomain& GetTriggerDomain() const { return fTriggerDomain; }

  /**
   * Returns the trigger domain object for modification.
   */
  AliHLTTriggerDomain& GetTriggerDomain() { return fTriggerDomain; }

  /**
   * Sets the trigger domain object.
   * \param value  The new value to use for the trigger domain.
   */
  void SetTriggerDomain(const AliHLTTriggerDomain& value) { fTriggerDomain = value; }
  
  /**
   * Returns the trigger description string.
   */
  const char* GetDescription() const { return fDescription.Data(); }
  
  /**
   * Sets the trigger description string.
   * \param value  The new value to use for the description string.
   */
  void SetDescription(const char* value) { fDescription = value; }
  
  /**
   * \returns true if the trigger description, trigger domain and readout list
   *    should be cleared for each new event.
   */
  bool WillClearInfoForNewEvent() const { return fClearInfo; }
  
  /**
   * Sets the flag indicating in the trigger description, trigger domain and
   * readout list should be cleared for each new event.
   * \param value  The new value to use for the flag.
   */
  void ClearInfoForNewEvent(bool value = true) { fClearInfo = value; }

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
   * Inherited from AliHLTComponent. This method is deprecated so we hide it.
   * Rather use the const version of this method which this method calls anyway.
   * @param list     list to receive the output data types.
   */
  virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  
  /**
   * Inherited from AliHLTComponent. This method is replaced with one that is
   * symmetric to GetInputDataTypes that returns void, so we make this method
   * private. The list will always contain kAliHLTDataTypeTObject, including whatever
   * values were added by the const version of GetOutputDataTypes.
   * @param list     list to receive the output data types.
   * @return the number of elements in the list.
   */
  virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
  
  const AliHLTComponentEventData* fEventData; //! Event data for the current event. Only valid inside DoTrigger.
  AliHLTComponentTriggerData* fTriggerData; //! Trigger data for the current event. Only valid inside DoTrigger.
  bool fDecisionMade;  //! Flag indicating if the trigger decision has been made for this trigger yet.
  bool fClearInfo;  //! Flag indicating if fDescription, fReadoutList and fTriggerDomain should be cleared for each new event.
  int fTriggerEventResult;  //! Result returned by PushBack method in the TriggerEvent method.
  TString fDescription;   //! The description to use for the trigger decision.
  AliHLTReadoutList fReadoutList; //! The DDL readout list object for the current event being processed.
  AliHLTTriggerDomain fTriggerDomain; //! The trigger domain object for the current event being processed.
  
  ClassDef(AliHLTTrigger, 0) // Base class for HLT triggers.

};

#endif // ALIHLTTRIGGER_H

