//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDHISTOMERGERCOMPONENT_H
#define ALIHLTTRDHISTOMERGERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *


#include "AliHLTProcessor.h"

/**
 * @class AliHLTTRDHistoMergerComponent
 * Component for adding histos from the histoComponents if those are running partition wise (SM wise) .
 * Expects all input blocks to be comparable.
 * 
 * Component ID: \b TRDHistoMerger <br>
 * Library: \b libAliHLTTRD.
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * 
 *
 * @ingroup alihlt_tpc_components
 */

class TH1;
class AliHLTTRDHistoMergerComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTRDHistoMergerComponent();
  /** destructor */
  virtual ~AliHLTTRDHistoMergerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  using AliHLTProcessor::DoEvent;
  
  int Configure(const char* arguments);

private:
  /** copy constructor prohibited */
  AliHLTTRDHistoMergerComponent(const AliHLTTRDHistoMergerComponent&);
  /** assignment operator prohibited */
  AliHLTTRDHistoMergerComponent& operator=(const AliHLTTRDHistoMergerComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 

  AliHLTUInt32_t fOutputSize;   // output size
  TH1* fHistoArr[9];            // array containing the added histos
  Bool_t fIncSM[18];            // array for telling which super module was already added

  ClassDef(AliHLTTRDHistoMergerComponent, 0);
};
#endif
