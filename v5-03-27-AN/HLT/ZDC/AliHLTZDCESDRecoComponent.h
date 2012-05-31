//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTZDCESDRECOCOMPONENT_H
#define ALIHLTZDCESDRECOCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTZDCESDRecoComponent.h 
    @author  Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>
    @brief   ZDC reconstruction component
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"
#include "AliESDEvent.h"

class AliZDCReconstructor;
class AliRawReaderMemory;


/**
 * @class AliHLTZDCESDRecoComponent
 * Reconstruction of ZDC data
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b ZDCESDReco <br>
 * Library: \b libAliHLTZDC.so     <br>
 * Input Data Types:  @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref kAliHLTDataTypeESDContent|kAliHLTDataOriginZDC <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting --> 
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>HLT/ConfigZDC/ZDCESDReco</tt>
 * \li -TObjString object holding a string with the configuration parameters
 *      currently empty 
 *
 * <tt>GRP/GRP/Data</tt>
 * \li -GRP object - run information
 *
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Input size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: xx Byte
 *
 * <h2>Output size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: Average : xx kByte
 *
 * <h2>Macros Tests</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>macros/makeConfigurationObjectZDCReconstruction.C</tt>
 * \li - Create configuration TObjString
 *
 * <tt>macros/HLTZDCTest.C</tt>
 * \li - Test macro for ZDC test in off-line environment
 *
 * <tt>macros/runZDCTest.sh</tt>
 * \li - Run Test macro HLTZDCTest.C
 *
 * @ingroup alihlt_zdc
 */
class AliHLTZDCESDRecoComponent : public AliHLTProcessor
{
    public:

        /** constructor */
	AliHLTZDCESDRecoComponent();
        /** destructor */
	virtual ~AliHLTZDCESDRecoComponent();

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
	        /** interface function, see @ref AliHLTComponent for description */
	int DoInit( int argc, const char** argv );
        /** interface function, see @ref AliHLTComponent for description */
	int DoDeinit();
        /** interface function, see @ref AliHLTComponent for description */
        int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
        /** interface function, see @ref AliHLTComponent for description */
        int ScanConfigurationArgument(int argc, const char** argv);
        /** interface function, see @ref AliHLTComponent for description */
        int Reconfigure(const char* cdbEntry, const char* chainId);
        /** interface function, see @ref AliHLTComponent for description */
        int ReadPreprocessorValues(const char* modules);
	
	using AliHLTProcessor::DoEvent;
		

    private:
  	/** copy constructor prohibited */
  	AliHLTZDCESDRecoComponent(const AliHLTZDCESDRecoComponent&);
  	/** assignment operator prohibited */
 	AliHLTZDCESDRecoComponent& operator=(const AliHLTZDCESDRecoComponent&);

  	/** rawreader instance */
  	AliRawReaderMemory  *fRawReader;     //! transient
      
	/** ZDC reconstructor instance */
        AliZDCReconstructor *fReconstructor; //! ZDC reconstructor

	ClassDef(AliHLTZDCESDRecoComponent, 0)

};

#endif
