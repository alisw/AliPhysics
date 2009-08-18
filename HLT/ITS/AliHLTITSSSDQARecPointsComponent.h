//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTITSSSDQARECPOINTSCOMPONENT_H
#define ALIHLTITSSSDQARECPOINTSCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        */ 
//* ALICE Experiment at CERN, All rights reserved.                         */
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTITSSSDQARecPointsComponent.h
    @author Ingrid Kielen
    @brief  Component for the SSD clusters QA
*/


#include "AliHLTProcessor.h"
#include "AliHLTITSSpacePointData.h"
#include "TClonesArray.h"
#include "AliITSRecPoint.h"

class AliHLTTPCConfMapper;

/**
 * @class AliHLTITSSSDQARecPointsComponent
 * Component for ploting charge in clusters
 * 
 * Component ID: \b ITSSSDQARecPoints <br>
 * Library: \b libAliHLTITS.
 *
 * Mandatory arguments: <br>
 * 
 * 
 * Optional arguments: <br>
 * 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTITSSSDQARecPointsComponent : public AliHLTProcessor
{
public:
  /** default constructor */
 AliHLTITSSSDQARecPointsComponent();
  /** destructor */
  virtual ~AliHLTITSSSDQARecPointsComponent();

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
  int DoInit(int argc, const char** argv);
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
 AliHLTITSSSDQARecPointsComponent(const AliHLTITSSSDQARecPointsComponent&);
  /** assignment operator prohibited */
 AliHLTITSSSDQARecPointsComponent& operator=(const AliHLTITSSSDQARecPointsComponent&);

  static const Int_t fgkSSDMODULES = 1698;      //total number of SSD modules
  static const Int_t fgkSSDLADDERSLAYER5 = 34; //ladders on layer 5
  static const Int_t fgkSSDLADDERSLAYER6 = 38; //ladders on layer 6
  static const Int_t fgkSSDMODULESPERLADDERLAYER5 = 22; //modules per ladder - layer 5
  static const Int_t fgkSSDMODULESPERLADDERLAYER6 = 25; //modules per ladder - layer 6
  static const Int_t fgkSSDMODULESLAYER5 = 748; //total number of SSD modules - layer5
  static const Int_t fgkSSDMODULESLAYER6 = 950; //total number of SSD modules - layer6
  static const Int_t fgkNumberOfPSideStrips = 768; //number of P-side strips
   
  TObjArray *fHistSSDArray; //TObjArray of the SSD QA histograms
   
  ClassDef(AliHLTITSSSDQARecPointsComponent, 0);

};
#endif
