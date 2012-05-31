// $Id$
// XEmacs -*-C++-*-

#ifndef ALIHLTTPCHWCFDATAREVERTERCOMPONENT_H
#define ALIHLTTPCHWCFDATAREVERTERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCHWCFDataReverterComponent.h
/// @author Kenneth Aamodt
/// @date   
/// @brief  Component for reverting data for the HW clusterfinder
///

#include "AliHLTProcessor.h"
#include "AliHLTTPCPad.h"
#include "AliHLTDataTypes.h"

class AliHLTTPCDigitReader;
class AliHLTTPCMapping;

/**
 * @class AliHLTTPCHWCFDataReverterComponent
 * Implementation of the zero suppression component.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 *
 * The component orders the data in the format the Hardware ClusterFinder 
 * expects, and revert the 40 bit altro words.
 * 
 * The component has the following component arguments:
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCHWCFDataReverterComponent : public AliHLTProcessor
    {
    public:
        /** constructor */
	AliHLTTPCHWCFDataReverterComponent();
	/** destructor */
	virtual ~AliHLTTPCHWCFDataReverterComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

      Int_t DeInitializePadArray();
      void InitializePadArray();
      
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
        int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	using AliHLTProcessor::DoEvent;

    private:

      /** copy constructor prohibited */
      AliHLTTPCHWCFDataReverterComponent(const AliHLTTPCHWCFDataReverterComponent&);
	
      /** assignment operator prohibited */
      AliHLTTPCHWCFDataReverterComponent& operator=(const AliHLTTPCHWCFDataReverterComponent&);
	
      /** the reader object for data decoding */
      AliHLTTPCDigitReader* fDigitReader;                              //!transient

      /** Vector of pointers to pad objects */
      typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;
      
      /** 2D vector of pointers to pad objects (vector of vectors)*/
      vector<AliHLTTPCPadVector> fRowPadVector;                        //! transient
      
      /** Array containing number of pads in the different rows */
      Int_t* fNumberOfPadsInRow;                                      //! transient
      
      /** Array containing the index f forst pad on row with hwaddress > 2048 */
      Int_t* fFirstPadHigh;                                      //! transient

      /** Number of rows the patch has */
      Int_t fNumberOfRows;                                            //! transient

      /** Current patch number */
      UInt_t fCurrentPatch;                                            //! transient

      /** First row in patch */
      UInt_t fFirstRow;                                                //! transient

      /** Last row in patch */
      UInt_t fLastRow;                                                 //! transient
      
      /** Number of timebins */
      Int_t fNTimeBins;                                                //! transient

      /** Flag to check if the 2d vector is initialized */
      Bool_t fVectorInitialized;                                       //! transient

      /** pointer to mapping object */
      AliHLTTPCMapping *fMapping;

      /** Flag to check if one should interleave the data */
      Bool_t fInterleave;

      ClassDef(AliHLTTPCHWCFDataReverterComponent, 0)
};
#endif
