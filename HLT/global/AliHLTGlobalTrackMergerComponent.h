// $Id$
#ifndef ALIHLTGLOBALTRACKMERGERCOMPONENT_H
#define ALIHLTGLOBALTRACKMERGERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTGlobalTrackMergerComponent.h
    @author Jacek Otwinowski
    @date   
    @brief  HLT global track merger component.
*/

#include "AliHLTProcessor.h"

class AliESDEvent;
class AliHLTGlobalTrackMerger;

/**
 * @class AliHLTGlobalTrackMergerComponent
 * The global track merger component
 *
 * @ingroup alihlt_global_components
 * @author Jacek.Otwinowski@gsi.de
 */
class AliHLTGlobalTrackMergerComponent : public AliHLTProcessor
    {
    public:
      /** standard constructor */
      AliHLTGlobalTrackMergerComponent();
      /** standard destructor */
      virtual ~AliHLTGlobalTrackMergerComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

      /** @see component interface AliHLTComponent::GetComponentID */
	const char* GetComponentID();
      
      /** @see component interface AliHLTComponent::GetInputDataTypes */
	void GetInputDataTypes(AliHLTComponentDataTypeList& list);

      /** @see component interface AliHLTComponent::GetOutputDataType */
	AliHLTComponentDataType GetOutputDataType();

      /** @see component interface AliHLTComponent::GetOutputDataSize */
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

      /** @see component interface AliHLTComponent::Spawn */
	AliHLTComponent* Spawn();

    protected:
	
       /**
       * Set the parameters
       */
	//void SetMergerParameters(Double_t maxy=2.0,Double_t maxz=3.0,
	//                         Double_t maxsnp=0.1,Double_t maxtgl=0.05, Double_t signed1Pt=0.003);
	void SetMergerParameters(Double_t maxy=1.e10,Double_t maxz=1.e10,
	                         Double_t maxsnp=1.e10,Double_t maxtgl=1.e10, Double_t signed1Pt=0.05);

	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 


      /** @see component interface AliHLTComponent::DoInit */
	int DoInit( int argc, const char** argv );

      /** @see component interface AliHLTComponent::DoDeinit */
	int DoDeinit();

      /** @see component interface @ref AliHLTProcessor::DoEvent */
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, AliHLTComponentBlockDataList& outputBlocks );

	using AliHLTProcessor::DoEvent;
	
    private:
      /** copy constructor prohibited */
      AliHLTGlobalTrackMergerComponent(const AliHLTGlobalTrackMergerComponent&);
      /** assignment operator prohibited */
      AliHLTGlobalTrackMergerComponent& operator=(const AliHLTGlobalTrackMergerComponent&);

      AliHLTGlobalTrackMerger *fGlobalTrackMerger;  //! global track merger
      AliESDEvent *fESD;                            //! AliESDEvent output from merger

      ClassDef(AliHLTGlobalTrackMergerComponent, 0)

    };
#endif
