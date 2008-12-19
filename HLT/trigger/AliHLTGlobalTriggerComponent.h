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

class AliHLTTriggerMenu;
class AliHLTGlobalTrigger;

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
   * Get a ratio by how much the data volume is shrunk or enhanced.
   * The method returns a size proportional to the trigger name string length
   * for constBase, and 1 for inputMultiplier.
   * @param constBase        <i>[out]</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>[out]</i>: multiplication ratio
   */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  
  /**
   * Initialise the component.
   * \param argc  The number of arguments in argv.
   * \param argv  Array of component argument strings.
   * \returns  Zero on success and negative number on failure.
   */
  virtual Int_t DoInit(int argc, const char** argv);
  
  /**
   * Cleanup the component.
   * \returns  Zero on success and negative number on failure.
   */
  virtual Int_t DoDeinit();
  
  /**
   * Spawn function creates a new object.
   * @return new class instance.
   */
  virtual AliHLTComponent* Spawn();

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
  
  /**
   * Generates the code for the global trigger to apply the given trigger menu.
   * The code will then be compiled on the fly and loaded. The name of the new
   * class is returned so that a new instance of the class can be created via:
   * \code
   *  AliHLTGlobalTrigger::CreateNew(name)
   * \endcode
   * where name is the name of the generated class as returned by this method.
   * \param  menu  The trigger menu to create the global trigger class from.
   * \param  name  The name of the generated class.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int GenerateTrigger(const AliHLTTriggerMenu* menu, TString& name);
  
  AliHLTGlobalTrigger* fTrigger;  //! Trigger object which implements the global trigger menu.
  
  ClassDef(AliHLTGlobalTriggerComponent, 0) // Global HLT trigger component class which produces the final trigger decision and readout list.
};

#endif // ALIHLTGLOBALTRIGGERCOMPONENT_H

