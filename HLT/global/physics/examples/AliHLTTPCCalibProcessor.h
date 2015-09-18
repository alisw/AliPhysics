//-*- Mode: C++ -*-
// $Id$
#ifndef ALITPCCALIBPROCESSOR_H
#define ALITPCCALIBPROCESSOR_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/// @file   AliHLTTPCCalibProcessor.h
/// @author Matthias Richter
/// @date   2012-02-01
/// @brief  A sample HLT processing component
/// @ingroup alihlt_tutorial

#include "AliHLTProcessor.h"
#include "AliHLTCalibrationProcessor.h"
class TH2F;
class TH1F;

/**
 * @class AliHLTTPCCalibProcessor
 * An example for a simple HLT processing component
 * The class features the AliHLTComponent interface for HLT processing
 * components. The interface allows to run such components in either
 * the (sequential) AliSimulation/AliReconstruction framework or the
 * parallel HLT online processing framework.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b TPCCalibProcessor                                  <br>
 * Library: \b libAliHLTTutorial.so                                    <br>
 * Input Data Types:                                                   <br>
 * Output Data Types:                                                  <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Argument scan is implemented in the function ScanConfigurationArgument().
 * see @ref alihltcomponent-initialization-ocdb.
 * Please provide specific descriptions and implementations.
 * \li -argument1     <i> parameter   </i> <br>
 *      an argument with one parameter
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -argument2                            <br>
 *      an argument without parameters
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * The component has just one default CDB entry in 
 * <tt>HLT/ConfigSample/TPCCalibProcessor</tt>.
 * \li -TObjString object holding a string with the configuration parameters
 *      explained above
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component has no output data.
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTTPCCalibProcessor : public AliHLTCalibrationProcessor {
public:
	AliHLTTPCCalibProcessor();
	virtual ~AliHLTTPCCalibProcessor();

	// AliHLTComponent interface functions
	const char* GetComponentID();
	void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
	AliHLTComponentDataType GetOutputDataType();
	void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  //	void GetOCDBObjectDescription( TMap* const targetMap);

	// Spawn function, return new class instance
	AliHLTComponent* Spawn();

protected:
	// AliHLTComponent interface functions - here they are not needed, because we are in a component that is a CalibProcessor!
  /*
	int DoInit( int argc, const char** argv );
	int DoDeinit();
	int DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
	       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
	       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
	int ScanConfigurationArgument(int argc, const char** argv);
	int Reconfigure(const char* cdbEntry, const char* chainId);
	int ReadPreprocessorValues(const char* modules);

	using AliHLTProcessor::DoEvent;
  */

  // functions of an AliHLTCalibrationProcessor

      using AliHLTCalibrationProcessor::ProcessCalibration;
      using AliHLTCalibrationProcessor::ShipDataToFXS;
      
      // Protected functions to implement AliHLTComponent's interface.
      // These functions provide initialization as well as the actual processing
      // capabilities of the component. 
      
      /** Initialize the calibration component. */
      Int_t InitCalibration();

      /** Scan commandline arguments of the calibration component. */
      Int_t ScanArgument( Int_t argc, const char** argv );

      /** DeInitialize the calibration component. */
      Int_t DeinitCalibration();

      /** Process the data in the calibration component. */
      Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData,
			  const AliHLTComponent_BlockData* blocks,
			  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
			  AliHLTUInt32_t& size,
			  vector<AliHLTComponent_BlockData>& outputBlocks);
  //Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

      /** Ship the data to the FXS at end of run or eventmodulo. */
      Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

private:
	/** copy constructor prohibited */
	AliHLTTPCCalibProcessor(const AliHLTTPCCalibProcessor&);
	/** assignment operator prohibited */
	AliHLTTPCCalibProcessor& operator=(const AliHLTTPCCalibProcessor&);

	Bool_t 	fGlobalDistr;       // flag to decide whether to store or not the global multiplicity histogram 
	Int_t   fMinTracks;         // min number of tracks to select an event
	Int_t   fMinEvents;         // min number of events to start calculating average
	Int_t   fEvents;            // number of events selected
	Float_t fMeanPt;            // sum of pt over all selected tracks
	Int_t   fTotTracks;         // number of tracks
	TH2F*   fHistoMeanPt;       //! 2D histo: mean pt vs total mult
	TH1F*   fHistoMult;         //! 1D histo: mult distribution
	TH1F*   fHistoPt;           //! 1D histo: pt distribution
	TH1F*   fHistoClusters;           //! 1D histo: cluster distribution
	TH1F*   fHistoClustersPerTrack;           //! 1D histo: cluster/track distribution


	ClassDef(AliHLTTPCCalibProcessor, 0)
};
#endif
