#ifndef ALIHLTV0HISTOCOMPONENT_H
#define ALIHLTV0HISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTV0HistoComponent.h
    @author Sergey Gorbunov
    @brief  Component for monitor V0 physics
*/

#include "AliHLTProcessor.h"
class TH1F;
class TH2F;

/**
 * @class AliHLTTPCV0HistoComponent
 * Component for monitor V0 physics 
 */
class AliHLTV0HistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTV0HistoComponent();
  /** destructor */
  virtual ~AliHLTV0HistoComponent();

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
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTV0HistoComponent(const AliHLTV0HistoComponent&);
  /** assignment operator prohibited */
  AliHLTV0HistoComponent& operator=(const AliHLTV0HistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
  
  TH1F *fGamma;   // Gamma inv. mass
  TH1F *fKShort;   // Ks inv. mass
  TH1F *fLambda;   // Lambda inv. mass
  TH2F *fAP;       // Armenteros-Podolanski 
  TH2F *fGammaXY;  // XY distribution of gamma convertions
  Int_t fNEvents; // n of processed events
  Int_t fNGammas; // n found total
  Int_t fNKShorts; // n found total
  Int_t fNLambdas; // n found total

  ClassDef(AliHLTV0HistoComponent, 0);

};
#endif
