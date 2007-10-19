// @(#) $Id$

#ifndef ALIHLTTPCCLUSTERFINDERCOMPONENT_H
#define ALIHLTTPCCLUSTERFINDERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCClusterFinderComponent.h
    @author Timm Steinbeck, Matthias Richter, Jochen Thaeder
    @date   
    @brief  The TPC cluster finder component.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class AliHLTTPCClusterFinder;
class AliHLTTPCPadArray;
class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCClusterFinderComponent
 * Implementation of the cluster finder component.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * The actual cluster finding algorithm is implemented in @ref AliHLTTPCClusterFinder.
 * The component can handle unpacked and packed data of different formats via the
 * AliHLTTPCDigitReader implementations. Two components are registered, the
 * TPCClusterFinderUnpacked and the TPCClusterFinderPacked. The latter one can
 * instantiate different digit readers depending on the arguments.
 * 
 * The component has the following component arguments:
 * - rawreadermode   the mode for the @ref AliHLTTPCDigitReaderRaw, use -2 if using unsorted
 * - adc-threshold   ADC count threshold for zero suppression, if <0 the base line
 *                   calculation and subtraction is switched off
 * - pp-run          set parameters specific to a pp run; currently this switches
 *                   cluster deconvolution off for pp runs (not true for unsorted reading)
 * - unsorted        if 1 the data will be read unsorted in to a PadArray object. This should
 *                   only be done on patch level since it use a lot of memory
 * - patch           specify on which patch to resd the data unsorted
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCClusterFinderComponent : public AliHLTProcessor
    {
    public:
        /**
         * constructor 
         * @param packed    whether to use the packed or unpacked reader 
         */
	AliHLTTPCClusterFinderComponent(bool packed);
	/** destructor */
	virtual ~AliHLTTPCClusterFinderComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	
    private:
	/** copy constructor prohibited */
	AliHLTTPCClusterFinderComponent(const AliHLTTPCClusterFinderComponent&);
	/** assignment operator prohibited */
	AliHLTTPCClusterFinderComponent& operator=(const AliHLTTPCClusterFinderComponent&);
	/** the cluster finder object */
	AliHLTTPCClusterFinder* fClusterFinder;                                      //!transient
	/** the reader object for data decoding */
	AliHLTTPCDigitReader* fReader;                                               //!transient

	bool fClusterDeconv; //!transient
      float fXYClusterError; //!transient
      float fZClusterError; //!transient
      /**
       * switch to indicated the reader
       * use fPackedSwitch = true for packed inputtype "gkDDLPackedRawDataType"
       * use fPackedSwitch = false for unpacked inputtype "gkUnpackedRawDataType"
       */
      Int_t fPackedSwitch;                                                           // see above
      
      /*
       * Reads the data the new unsorted way if true
       *
       */
      Int_t fUnsorted;                                                               //!transient

      /*
       * Patch number to be read, currently given as component argument,
       * will be changed later.
       */
      Int_t fPatch;                                                                  //!transient

      /*
       * Pointer to a PadArray object containing a double array of all the pads in
       * the current patch.
       */
      AliHLTTPCPadArray * fPadArray;                                                 //!transient

      ClassDef(AliHLTTPCClusterFinderComponent, 1)

    };
#endif
