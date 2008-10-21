// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTTPCGLOBALMERGERCOMPONENT_H
#define ALIHLTTPCGLOBALMERGERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCGlobalMergerComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  HLT TPC global merger component.
*/

#include "AliHLTProcessor.h"

class AliHLTTPCGlobalMerger;
class AliHLTTPCVertex;

/**
 * @class AliHLTTPCGlobalMergerComponent
 * The TPC global merger component
 * The component is the interface to the AliHLTGlobalMerger processing
 * class.
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCGlobalMergerComponent : public AliHLTProcessor
    {
    public:
      /** standard constructor */
      AliHLTTPCGlobalMergerComponent();
      /** standard destructor */
      virtual ~AliHLTTPCGlobalMergerComponent();

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
	void SetMergerParameters(Double_t maxy=2.0,Double_t maxz=3.0,Double_t maxkappa=0.003,
				 Double_t maxpsi=0.1,Double_t maxtgl=0.05);

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
      AliHLTTPCGlobalMergerComponent(const AliHLTTPCGlobalMergerComponent&);
      /** assignment operator prohibited */
      AliHLTTPCGlobalMergerComponent& operator=(const AliHLTTPCGlobalMergerComponent&);

      /** the global merger object */
      AliHLTTPCGlobalMerger* fGlobalMerger; //!
      /** the vertexer object */
      AliHLTTPCVertex* fVertex; //!

      struct SliceData {
	/** slice no */
	int fSlice;                                                // see above
	/** block descriptor for the vertex data block */
	const AliHLTComponentBlockData* fVertexBlock;              //! transient
	/** index of the vertex data block */
	unsigned fVertexBlockIndex;                                // see above
	/** block descriptor for the tracklet data block */
	const AliHLTComponentBlockData* fTrackletBlock;            //! transient
	/** index of the tracklet data block */
	unsigned fTrackletBlockIndex;                              // see above
      };

	ClassDef(AliHLTTPCGlobalMergerComponent, 0)

    };
#endif
