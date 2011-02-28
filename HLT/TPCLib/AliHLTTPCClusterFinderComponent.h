// $Id$

#ifndef ALIHLTTPCCLUSTERFINDERCOMPONENT_H
#define ALIHLTTPCCLUSTERFINDERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTPCClusterFinderComponent.h
/// @author Timm Steinbeck, Matthias Richter, Kenneth Aamodt
/// @date   
/// @brief  The TPC cluster finder component.
///

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"

class AliHLTTPCClusterFinder;
class AliHLTTPCDigitReader;
class AliTPCTransform;

/**
 * @class AliHLTTPCClusterFinderComponent
 * Implementation of the cluster finder component.
 * The component implements the interface methods of the @ref AliHLTProcessor.
 * The actual cluster finding algorithm is implemented in @ref AliHLTTPCClusterFinder.
 * Two components are registered, TPCClusterFinderUnpacked is for reading the HLT
 * internal digit data format used in the simulation. TPCClusterFinder32Bit uses
 * the AliHLTTPCDigitReader for raw data. After a phase of different decoder/raw stream
 * implementations the CF for raw data is using the default offline raw stream for
 * the 32bit RCU format AliAltroRawStreamV3, which also has a fall back to the 40bit
 * AliAltroRawStream if the old RCU format is detected.
 *
 * The clusterfinder is now using the AliTPCTransform instead of the AliHLTTPCTransform for  
 * transformations from row, pad time -> x,y,z.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCClusterFinderUnpacked and TPCClusterFinder32Bit <br>
 * Library: \b libAliHLTTPC
 * Input Data Types: @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref AliHLTTPCDefinitions::fgkClustersDataType and/or kAliHLTDataTypeHwAddr16 <br> 
 *
 *
 * Mandatory arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * Optional arguments: <br>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -deconvolute-time <br>  
 *      Turns on deconvolution in the time direction.
 * \li -deconvolute-pad <br>  
 *      Turns on deconvolution in the pad direction.
 * \li -timebins <br>  
 *      Sets the number of timebins (446 for simulated data, and 1024 for real data) Default:1024
 * \li -first-timebin <br>  
 *      First timebin taken into consideration when reading the data. 
 * \li -last-timebin <br>  
 *      Last timebin taken into consideration when reading the data. 
 * \li -sorted <br>  
 *      Switch off unsorted reading of data. Equivalent to the old argument unsorted 0.
 * \li -active-pads <br>  
 *      Switch off unsorted reading of data. Equivalent to the old argument unsorted 0.
 * \li -occupancy-limit <br>  
 *      Set the occupancy limit for the sorted clusterfinding.
 *
 *
 * Obsolete arguments: <br>
 * \li occupancy-limit <br>  
 * \li rawreadermode   <br>  
 * \li pp-run          <br>  
 * \li adc-threshold   <br>  
 * \li oldrcuformat    <br>  
 * \li unsorted        <br>  
 * \li nsigma-threshold <br>  
 *
 * <h2>Default CDB entries:</h2>
 * The component has these default CDB entries
 * \li <tt>GRP/GRP/Data</tt>.               
 * \li <tt>TPC/Calib/PadTime0</tt>.         
 * \li <tt>TPC/Calib/Parameters</tt>.       
 * \li <tt>TPC/Calib/TimeDrift</tt>.        
 * \li <tt>TPC/Calib/Temperature</tt>.      
 *
 * TODO: pad by pad gain calibration also has to be added to the clusterfinder
 *
 * And it also needs these below to avoid warnings during initialization and update of calibDB
 * \li <tt>TPC/Calib/PadGainFactor</tt>.    
 * \li <tt>TPC/Calib/TimeGain</tt>.
 * \li <tt>TPC/Calib/GainFactorDedx</tt>.
 * \li <tt>TPC/Calib/PadNoise</tt>.
 * \li <tt>TPC/Calib/Pedestals</tt>.
 * \li <tt>TPC/Calib/ClusterParam</tt>.
 * \li <tt>TPC/Calib/AltroConfig</tt>.
 * \li <tt>TPC/Calib/Pulser</tt>.
 * \li <tt>TPC/Calib/CE</tt>.
 * \li <tt>TPC/Calib/Raw</tt>.
 * \li <tt>TPC/Calib/QA</tt>.
 * \li <tt>TPC/Calib/Mapping</tt>.
 * \li <tt>TPC/Calib/Goofie</tt>.
 * \li <tt>TPC/Calib/HighVoltage</tt>.
 * \li <tt>TPC/Calib/Ref</tt>.
 *
 * These entries are used by the AliTPCTransform class to correct for T0, drift and ExB.
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCClusterFinderComponent : public AliHLTProcessor
    {
    public:
      /**
       * Defines for the cluster finder type.
       * The cluster finders can work on different formats of input data,
       * the AliHLTTPCDigitReader interface provides a transparent way to
       * read the data.
       */
      enum {
	// deprecated option for offline AliAltroRawStream
	kClusterFinderPacked,
	// Unpacked data of format AliHLTTPCUnpackedRawData */
	kClusterFinderUnpacked,
	// deprecated option for AliAltroDecoder
	kClusterFinderDecoder,
	// real data, offline altro decoder 32 bit format*/
	kClusterFinder32Bit
      };

        /**
         * constructor 
         * @param mode    input type see e.g. @ref kClusterFinderUnpacked
         */
	AliHLTTPCClusterFinderComponent(int mode);
	/** destructor */
	virtual ~AliHLTTPCClusterFinderComponent();

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
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	int Configure(const char* arguments);
	int ScanConfigurationArgument(int argc, const char** argv);
	int Reconfigure(const char* cdbEntry, const char* chainId);
	
	using AliHLTProcessor::DoEvent;

    private:
	/** standard constructor prohibited */
	AliHLTTPCClusterFinderComponent();
	/** copy constructor prohibited */
	AliHLTTPCClusterFinderComponent(const AliHLTTPCClusterFinderComponent&);
	/** assignment operator prohibited */
	AliHLTTPCClusterFinderComponent& operator=(const AliHLTTPCClusterFinderComponent&);
	/** the cluster finder object */
	AliHLTTPCClusterFinder* fClusterFinder;                                      //!transient
	/** the reader object for data decoding */
	AliHLTTPCDigitReader* fReader;                                               //!transient

	/** flag to deconvolute in time direction */
	Bool_t fDeconvTime;                                                          //!transient
        /** the object to set the time stamp */
        AliTPCTransform *fTS; //!transient

	/** flag to deconvolute in pad direction */
	Bool_t fDeconvPad;                                                           //!transient
	/** flag to switch on/off deconvolution in pad and time directions (used by sorted clusterfinding method) */
	bool fClusterDeconv;                                                         //!transient
	/** Error in xy of cluster */
	float fXYClusterError;                                                       //!transient
	/** Error in Z direction of cluster */
	float fZClusterError; //!transient
	/**
	 * switch to indicated the reader
	 * use fModeSwitch = 0 for packed inputtype "gkDDLPackedRawDataType"
	 * use fModeSwitch = 1 for unpacked inputtype "gkUnpackedRawDataType"
	 * use fModeSwitch = 2 for packed inputtype "gkDDLPackedRawDataType" with new digit reader
	 * use fModeSwitch = 3 for packed inputtype "gkDDLPackedRawDataType" with 32bit digit reader
	 */
	Int_t fModeSwitch;                                                            // see above
      
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
	 * Switch to specify if one ship out a list of active pads (pads conected to a cluster).
	 */
	Int_t fGetActivePads;                                                          //!transient

	/** First timebin taken into account when reading the data */
	Int_t fFirstTimeBin;                                                           //!transient

	/** Last timebin taken in to account when reading the data */
	Int_t fLastTimeBin;                                                            //!transient

	Bool_t fDoMC; // flag to provide MC labels
	Bool_t fReleaseMemory; // flag to release the memory after each event
	AliHLTComponentBenchmark fBenchmark; // benchmark

	ClassDef(AliHLTTPCClusterFinderComponent, 0)

};
#endif
