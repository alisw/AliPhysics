#ifndef ALIHLTGLOBALTRIGGERCOMPONENT_H
#define ALIHLTGLOBALTRIGGERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTriggerComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Declaration of the AliHLTGlobalTriggerComponent component class.

#include "AliHLTTrigger.h"

/**
 * \class AliHLTGlobalTriggerComponent
 * This class applies the global HLT trigger to all trigger information produced
 * by components deriving from AliHLTTrigger.
 * Any information delivered by other components in data blocks that contain
 * TObjects can also be used for the trigger algorithm.
 */
class AliHLTGlobalTriggerComponent : public AliHLTTrigger
{
 public:
 
  AliHLTGlobalTriggerComponent();
  virtual ~AliHLTGlobalTriggerComponent();
  
  /**
   * Inherited from AliHLTTrigger.
   * @return string containing the global trigger name.
   */
  virtual const char* GetTriggerName() const { return "HLTGlobalTrigger"; };

  /**
   * Returns extra output data types this trigger generates.
   * This returns an kAliHLTDataTypeTObject in <i>list</i>.
   * @param list <i>[out]</i>: The list of data types to be filled.
   */
  virtual void GetOutputDataTypes(AliHLTComponentDataTypeList& list) const
  {
    list.push_back(kAliHLTDataTypeTObject);
  }

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

  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent(const AliHLTGlobalTriggerComponent& obj);
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent& operator = (const AliHLTGlobalTriggerComponent& obj);

  /**
   * Applies the global HLT trigger.
   * @return Zero is returned on success and a negative error code on failure.
   */
  virtual int DoTrigger();
  
 private:
  
  ClassDef(AliHLTGlobalTriggerComponent, 0) // Global HLT trigger component class which produces the final trigger decision and readout list.
};

#endif // ALIHLTGLOBALTRIGGERCOMPONENT_H

