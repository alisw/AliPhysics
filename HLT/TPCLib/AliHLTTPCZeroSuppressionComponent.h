// XEmacs -*-C++-*-
// @(#) $Id: AliHLTTPCClusterFinderComponent.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIHLTTPCZEROSUPPRESSIONCOMPONENT_H
#define ALIHLTTPCZEROSUPPRESSIONCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCZeroSuppressionComponent.h
    @author Kenneth Aamodt
    @date   
    @brief  Component for ZeroSuppression
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"
#include "AliHLTTPCPad.h"

class AliHLTTPCDigitReader;

/**
 * @class AliHLTTPCZeroSuppressionComponent
 * Implementation of the zero suppression component.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * It reads the data pad by pad and zerosuppress the data. The output is unpacked which is 
 * sent to the clulsterfinder.
 * 
 * The component has the following component arguments:
 * - adc-threshold   ADC count threshold for zero suppression.
 *
 * - rms-threshold   RMS threshold for zero suppression.
 *          
 * - first-timebin   The first timebin for zero suppression
 *
 * - last-timebin    The last timebin for zero suppression
 *
 * - occupancy-limit Minimum number of timebins with signal
 *
 * - sort-pads Flag to switch on pad sorting(needed by the SORTED clusterfinder)
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCZeroSuppressionComponent : public AliHLTProcessor
    {
    public:
        /** constructor */
	AliHLTTPCZeroSuppressionComponent();
	/** destructor */
	virtual ~AliHLTTPCZeroSuppressionComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see @ref AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see @ref AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see @ref AliHLTComponent for description */
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
	AliHLTTPCZeroSuppressionComponent(const AliHLTTPCZeroSuppressionComponent&);
	
	/** assignment operator prohibited */
	AliHLTTPCZeroSuppressionComponent& operator=(const AliHLTTPCZeroSuppressionComponent&);
	
	/** the reader object for data decoding */
	AliHLTTPCDigitReader* fDigitReader;                                               //!transient
	
	
      typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;

      vector<AliHLTTPCPadVector> fRowPadVector;                        //! transient
      
      UInt_t* fNumberOfPadsInRow;                                      //! transient
      
      UInt_t fNumberOfRows;                                            //! transient

      UInt_t fCurrentPatch;                                            //! transient
      UInt_t fFirstRow;                                                //! transient
      UInt_t fLastRow;                                                 //! transient
      

      Int_t fStartTimeBin;                                            //! transient
      Int_t fEndTimeBin;                                              //! transient
      Int_t fNTimeBins;                                               //! transient
      Double_t fNRMSThreshold;                                           //! transient
      Int_t fSignalThreshold;                                         //! transient
      Int_t fMinimumNumberOfSignals;                                  //! transient
      UInt_t fOldRCUFormat;                                            //! transient
      Bool_t fSortPads;                                                //! transient
      Bool_t fVectorInitialized;                                       //! transient

      Int_t fValueBelowAverage;                                        //! transient
      Int_t fLeftTimeBin;                                              //! transient
      Int_t fRightTimeBin;                                             //! transient
      ClassDef(AliHLTTPCZeroSuppressionComponent, 0)  
    };
#endif
