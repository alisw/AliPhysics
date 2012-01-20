#ifndef AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
#define AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

///
///  @file   AliHLTMUONMansoTrackerFSMComponent.h
///  @author Artur Szostak <artursz@iafrica.com>,
///          Indranil Das <indra.das@saha.ac.in>
///  @date   18 Sep 2007
///  @brief  Tracker component for the dimuon HLT using the Manso algorithm.
///
/// The tracker component performs minimal track reconstruction in stations 4 & 5.
/// It uses the Manso algorithm implemented as a finite state machine.
///

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONMansoTrackerFSMCallback.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif

class AliHLTMUONMansoTrackerFSM;
//class AliHLTMUONMansoTracksBlockWriter;
extern "C" {
struct AliHLTMUONRecHitStruct;
} // extern C


/**
 * @class AliHLTMUONMansoTrackerFSMComponent
 * @brief Dimuon HLT tracker using the Manso tracking algorithm implemented as a FSM.
 *
 * This is a tracker component for the muon spectrometer. It performs minimal track
 * reconstruction in stations 4 & 5.
 * The Manso algorithm is used, implemented as a finite state machine (FSM).
 * This makes the component fast, insensitive to residual missalignment and gives
 * reasonable pT resolution for improved event selectivity.
 * However, because of its minimalism and simplicity, it does suffer from an increased
 * fake track rate, in particular for higher multiplicity events.
 * This should not pose a problem at all for p+p or peripheral A+A events, while
 * central events should and will be triggered anyway.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b MUONMansoTrackerFSM <br>
 * Library: \b libAliHLTMUON.so  <br>
 * Input Data Types: ('RECHITS ', 'MUON'); ('TRIGRECS', 'MUON') <br>
 * Output Data Types: ('MANTRACK', 'MUON'); ('MNCANDID', 'MUON') <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None. <br>
 *
 * <h2>Optional arguments:</h2>
 * \li -zmiddle <i>number</i> <br>
 *      Sets the value for the Z coordinate for the middle of the magnetic field.
 *      The <i>number</i> should be a floating point number representing the Z
 *      coordinate in centimeters (cm).
 *      A reasponable value is around -975cm (the default). <br>
 * \li -bfieldintegral <i>number</i> <br>
 *      Sets the value for the dipole's magnetic field integral.
 *      The <i>number</i> should be a floating point number representing the integral
 *      in units of Tesla meters (T.m).
 *      The value can be negative or possitive to indicate the dipole's polarity. <br>
 * \li -a7 <i>number</i> <br>
 * \li -a8 <i>number</i> <br>
 * \li -a9 <i>number</i> <br>
 * \li -a10 <i>number</i> <br>
 * \li -b7 <i>number</i> <br>
 * \li -b8 <i>number</i> <br>
 * \li -b9 <i>number</i> <br>
 * \li -b10 <i>number</i> <br>
 *      These allow one to set the A and B region of interest parameters used by the
 *      tracker to search for correlated hits on subsequent chambers.
 *      The number indicated after -a* or -b* indicates the chamber to set the parameter for.
 *      i.e. "-a7 123" will set the A parameter for chamber 7 to the value of 123.
 *      The region of interest parameters are defined as follows: Rs = a*Rp + b,
 *      where Rs is the region of interest radius and Rp is the distance from the beam line
 *      to the centre of the region of interest (in cm).
 *      The <i>number</i> should be a floating point number in each case and the values
 *      for the B parameters should always be positive. <br>
 * \li -z7 <i>number</i> <br>
 * \li -z8 <i>number</i> <br>
 * \li -z9 <i>number</i> <br>
 * \li -z10 <i>number</i> <br>
 * \li -z11 <i>number</i> <br>
 * \li -z13 <i>number</i> <br>
 *      These allow one to set the Z coordinate used for the chamber position in the
 *      trackers internal calculations.
 *      The number after -z* indicates the chamber number to set the value for.
 *      Thus, "-z11 -123" sets chamber 11's Z coordinate position to -123 cm.
 *      The <i>number</i> should be a floating point number in each case and negative,
 *      since the spectrometer is located in the negative Z direction.
 *      The units are in centimeters (cm).
 *      The default values are the nominal chamber Z coordinates. <br>
 * \li -warn_on_unexpected_block <br>
 *      This will cause the component to generate warnings when it receives data block
 *      types it does not know how to handle. Without this option the component only
 *      generates debug messages when they are compiled in. <br>
 * \li -cdbpath <i>path</i> <br>
 *      This allows one to override the path to use for the CDB location.
 *      <i>path</i> must be a valid CDB URI. By default the HLT system framework
 *      sets the CDB path. <br>
 * \li -run <i>number</i> <br>
 *      This allows one to override the run number to use. <i>number</i> must be
 *      a positive integer number. By default the HLT system framework sets the
 *      run number. <br>
 * \li -delaysetup <br>
 *      If indicated then part of the initialisation of the component is forcefully
 *      delayed to the first event received, i.e. the Start-of-Run event. <br>
 * \li -dumponerror <br>
 *      This flag will cause the component to dump the data blocks it received if
 *      an error occurs during the processing of an event. <br>
 * \li -dumppath <i>path</i> <br>
 *      Allows one to specify the path in which to dump the received data blocks
 *      if an error occurs. <br>
 * \li -makecandidates <br>
 *      Indicates if the Manso track candidates data block should be generated.
 *      This kind of information is useful for debugging. <br>
 *
 * <h2>Standard configuration:</h2>
 * This component should normally be configured with no extra options.
 * It will automatically load the required reconstruction parameters from the CDB. <br>
 *
 * <h2>Default CDB entries:</h2>
 * The component loads the reconstruction parameters from the MUON subdirectory
 * in the CDB.
 * The magnetic field integral, A and B region of interest parameters,
 * and nominal Z coordinates are stored in a TMap under HLT/ConfigMUON/MansoTrackerFSM.
 *
 * <h2>Performance:</h2>
 * Can achieve about 3kHz processing rate for nominal event sizes containing
 * 150 tracks per event.
 *
 * <h2>Memory consumption:</h2>
 * Minimal memory consumption on the order of megabytes.
 * 5 Mbytes should be more than enough.
 *
 * <h2>Output size:</h2>
 * Output size is about equivalent to the incoming reconstructed hit and trigger
 * record input data.
 *
 * @ingroup alihlt_dimuon_component
 */
class AliHLTMUONMansoTrackerFSMComponent
	: public AliHLTMUONProcessor, public AliHLTMUONMansoTrackerFSMCallback
{
public:
	AliHLTMUONMansoTrackerFSMComponent();
	virtual ~AliHLTMUONMansoTrackerFSMComponent();

	// Public functions to implement the AliHLTProcessor interface.
	// These functions are required for the registration process.
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	
	// Inherited from AliHLTMUONMansoTrackerFSMCallback:
	virtual void RequestClusters(
			AliHLTMUONMansoTrackerFSM* tracker,
			AliHLTFloat32_t left, AliHLTFloat32_t right,
			AliHLTFloat32_t bottom, AliHLTFloat32_t top,
			AliHLTMUONChamberName chamber, const void* tag
		);
	virtual void EndOfClusterRequests(AliHLTMUONMansoTrackerFSM* tracker);
	virtual void FoundTrack(AliHLTMUONMansoTrackerFSM* tracker);
	virtual void NoTrackFound(AliHLTMUONMansoTrackerFSM* tracker);

protected:

	// Protected functions to implement the AliHLTProcessor interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 
	virtual int DoInit(int argc, const char** argv);
	virtual int Reconfigure(const char* cdbEntry, const char* componentId);
	virtual int ReadPreprocessorValues(const char* modules);
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			const AliHLTComponentBlockData* blocks,
			AliHLTComponentTriggerData& trigData,
			AliHLTUInt8_t* outputPtr,
			AliHLTUInt32_t& size,
			AliHLTComponentBlockDataList& outputBlocks
		);
	
	using AliHLTProcessor::DoEvent;

private:

	// Do not allow copying of this class.
	AliHLTMUONMansoTrackerFSMComponent(const AliHLTMUONMansoTrackerFSMComponent& /*obj*/);
	AliHLTMUONMansoTrackerFSMComponent& operator = (const AliHLTMUONMansoTrackerFSMComponent& /*obj*/);

	void Reset();
	void FreeMemory();
	
	void AddRecHits(
			AliHLTUInt32_t specification,
			const AliHLTMUONRecHitStruct* recHits,
			AliHLTUInt32_t count
		);
	
	void ResetCanLoadFlags();
	bool AtLeastOneCanLoadFlagsIsSet() const;
	
	/**
	 * Reads this component's configuration parameters from the CDB.
	 * These include the middle of the dipole Z coordinate (zmiddle), the
	 * integrated magnetic field of the dipole, Z coordinates of the chambers
	 * and the region of interest parameters used during the tracking.
	 * \note Only those parameters are loaded from CDB for which the fCanLoadxyz
	 *       flags are true.
	 * \return 0 if no errors occured and negative error code compatible with
	 *       the HLT framework on errors.
	 */
	int ReadConfigFromCDB();
	

	AliHLTMUONMansoTrackerFSM* fTracker;  //! Tracker to do the actual work.
	
	AliHLTUInt32_t fTrackCount;  //! Number of tracks currently found.
	/*AliHLTMUONMansoTracksBlockWriter*/void* fBlock;  //! The current data block we are writing.
	
	class AliRecHitBlockInfo
	{
	public:
	
		AliRecHitBlockInfo(AliHLTUInt32_t count = 0, const AliHLTMUONRecHitStruct* data = NULL) :
			fCount(count), fData(data)
		{}
	
		// Perform a shallow copy.
		AliRecHitBlockInfo(const AliRecHitBlockInfo& obj) :
			fCount(obj.fCount), fData(obj.fData)
		{}
		
		AliRecHitBlockInfo& operator = (const AliRecHitBlockInfo& obj)
		{
		  if(&obj == this) return *this;
			fCount = obj.fCount;
			fData = obj.fData;
			return *this;
		}
		
		AliHLTUInt32_t Count() const { return fCount; }
		const AliHLTMUONRecHitStruct* Data() const { return fData; }
	
	private:
		AliHLTUInt32_t fCount;  // Number of elements in fData.
		const AliHLTMUONRecHitStruct* fData; // Pointer to the array of rec hits.
	};
	
	//std::vector<AliRecHitBlockInfo> fRecHitBlock[4];  //! Arrays of rec hit block data.
	AliHLTUInt32_t fRecHitBlockArraySize;  ///< The array size of each array in fRecHitBlock.
	AliHLTUInt32_t fRecHitBlockCount[4];   ///< The number of records actually stored in fRecHitBlock[i].
	// The following are 4 dynamic arrays of AliRecHitBlockInfo structures.
	// These arrays will all have the same size = fRecHitBlockArraySize.
	// The array itself is actually allocated only once and the pointer stored in fRecHitBlock[0],
	// while the other pointers fRecHitBlock[i] {i>0} will just be set relative to fRecHitBlock[0].
	// The allocated memory is: 4 * fRecHitBlockArraySize * sizeof(AliRecHitBlockInfo).
	AliRecHitBlockInfo* fRecHitBlock[4];  //! Arrays of rec hit block data.

	bool fWarnForUnexpecedBlock;  ///< Flag indicating if we should log a warning if we got a block of an unexpected type.	
	bool fCanLoadZmiddle;  ///< Indicates if the zmiddle parameter can be loaded from CDB.
	bool fCanLoadBL;  ///< Indicates if the bfieldintegral parameter can be loaded from CDB.
	bool fCanLoadA[4];  ///< Indicates if the roi_paramA_chamber[7..10] parameter can be loaded from CDB.
	bool fCanLoadB[4];  ///< Indicates if the roi_paramB_chamber[7..10] parameter can be loaded from CDB.
	bool fCanLoadZ[6];  ///< Indicates if the chamber[7..11,13]postion parameter can be loaded from CDB.
	
	ClassDef(AliHLTMUONMansoTrackerFSMComponent, 0);  // Manso tracker component implemented as a finite state machine (FSM).
};

#endif // AliHLTMUONMANSOTRACKERFSMCOMPONENT_H
