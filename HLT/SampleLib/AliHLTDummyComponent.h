// $Id$

#ifndef ALIHLTDUMMYCOMPONENT_H
#define ALIHLTDUMMYCOMPONENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTDummyComponent.h
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  Declaration of a dummy component. */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

/**
 * @class AliHLTDummyComponent
 * @brief A dummy HLT processing component. 
 *
 * An implementiation of a dummy component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * <h2>General properties:</h2>
 *
 * Component ID: \b Dummy <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 *
 * Output Data Types: depending on input blocks
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -output_percentage                            <br>
 *      The fraction (%) of the input data to be copied to the output
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has no default CDB entry.
 * It does not load any configuration from the global <tt>ConfigHLT</tt>
 * folder.
 * 
 * <h2>Performance:</h2>
 * The component does not have any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not have any event data processing.
 *
 * <h2>Output size:</h2>
 * Output multiplier determined by option -output_percentage
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTDummyComponent : public AliHLTProcessor
    {
    public:
	AliHLTDummyComponent();
	virtual ~AliHLTDummyComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
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
	
	using AliHLTProcessor::DoEvent;

    private:

	/** The size of the output data produced, as a percentage of the input data's size.
	    Can be greater than 100 (%) */
	unsigned fOutputPercentage; // see above
	
	ClassDef(AliHLTDummyComponent, 0)

    };
#endif
