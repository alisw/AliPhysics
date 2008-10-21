// XEmacs -*-C++-*-
// @(#) $Id: AliHLTTPCClusterFinderComponent.h 23318 2008-01-14 12:43:28Z hristov $

#ifndef ALIHLTTPCZEROSUPPRESSIONCOMPONENT_H
#define ALIHLTTPCZEROSUPPRESSIONCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCZeroSuppressionComponent.h
    @author Kenneth Aamodt
    @date   
    @brief  Component for ZeroSuppression
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCPad.h"
#include "AliHLTDataTypes.h"

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
 * @ingroup alihlt_tpc_components
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
      AliHLTTPCZeroSuppressionComponent(const AliHLTTPCZeroSuppressionComponent&);
	
      /** assignment operator prohibited */
      AliHLTTPCZeroSuppressionComponent& operator=(const AliHLTTPCZeroSuppressionComponent&);
	
      /** the reader object for data decoding */
      AliHLTTPCDigitReader* fDigitReader;                              //!transient

      /** Vector of pointers to pad objects */
      typedef vector<AliHLTTPCPad*> AliHLTTPCPadVector;
      
      /** 2D vector of pointers to pad objects (vector of vectors)*/
      vector<AliHLTTPCPadVector> fRowPadVector;                        //! transient
      
      /** Array containing number of pads in the different rows */
      Int_t* fNumberOfPadsInRow;                                      //! transient
      
      /** Number of rows the patch has */
      Int_t fNumberOfRows;                                            //! transient

      /** Current patch number */
      UInt_t fCurrentPatch;                                            //! transient

      /** First row in patch */
      UInt_t fFirstRow;                                                //! transient

      /** Last row in patch */
      UInt_t fLastRow;                                                 //! transient
      
      /** First timebin to include in zerosuppression */
      Int_t fStartTimeBin;                                             //! transient

      /** Lasr timebin to include in zerosuppression */
      Int_t fEndTimeBin;                                               //! transient

      /** Number of timebins */
      Int_t fNTimeBins;                                                //! transient

      /** Number of RMS the signal has to be larger than */
      Double_t fNRMSThreshold;                                         //! transient

      /** Signal threshold (signal has to be greater than average + this number) */
      Int_t fSignalThreshold;                                          //! transient

      /** Minimum number of signals to do zerosuppression */
      Int_t fMinimumNumberOfSignals;                                   //! transient

      /** OldRCUFormat flag */
      UInt_t fOldRCUFormat;                                            //! transient

      /** Sort pads flag */
      Bool_t fSortPads;                                                //! transient

      /** Flag to check if the 2d vector is initialized */
      Bool_t fVectorInitialized;                                       //! transient

      /** Value below average (useful for noisy pads to get the tails of a signal) */
      Int_t fValueBelowAverage;                                        //! transient

      /** Number of timebins to look left(decreasing time direction) for tails */
      Int_t fLeftTimeBin;                                              //! transient

      /** Number of timebins to look right(increasing time direction) for tails */
      Int_t fRightTimeBin;                                             //! transient

      /** Flag to switch on active pads selection */
      Bool_t fGetActivePads;                                           //! transient

      /** Flag to switch off data being sent to output */
      Bool_t fSkipSendingZSData;

      /** Flag to switch off hw list being sent to output */
      Bool_t fSendHWList;

      /** Vector of active pad hardware addresses */
      vector<AliHLTUInt16_t> fHwAddressList;                  //! transient

      ClassDef(AliHLTTPCZeroSuppressionComponent, 0)  
    };
#endif
