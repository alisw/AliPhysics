//-*- Mode: C++ -*-
// $Id: $
#ifndef ALIHLTTRIGGERCOUNTERCOMPONENT_H
#define ALIHLTTRIGGERCOUNTERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// \file   AliHLTTriggerCounterComponent.h
/// \author Artur Szostak <artursz@iafrica.com>
/// \date   3 Nov 2010
/// \brief  Declaration of the AliHLTTriggerCounterComponent component class.

#include "AliHLTProcessor.h"
#include "AliHLTTriggerCounters.h"
#include "THashTable.h"
#include "TTimeStamp.h"
#include "TMap.h"

/**
 * \class AliHLTTriggerCounterComponent
 * This component is used to calculate online the trigger counters and rates for
 * output HLT and input CTP / HLT triggers. The component should be connected to
 * all the global trigger component instances in the chain and there should be
 * only one instance of AliHLTTriggerCounterComponent in the chain. If more instances
 * of AliHLTTriggerCounterComponent exist in the chain then the statistics will not
 * be calculated correctly.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b HLTTriggerCounter <br>
 * Library: \b libAliHLTTrigger.so   <br>
 * Input Data Types:  kAliHLTDataTypeGlobalTrigger | kAliHLTDataTypeTriggerDecision <br>
 * Output Data Types: kAliHLTDataTypeInputTriggerCounters | kAliHLTDataTypeOutputTriggerCounters <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None.
 *
 * <h2>Optional arguments:</h2>
 * \li -config <i>filename</i> <br>
 *      Indicates the configuration macro file to use for the initial trigger counter lists.
 * \li -publishperiod <i>period</i> <br>
 *      Sets the period between publishing of the trigger counters in seconds.
 *      The default is to publish for every event.
 * \li -skipcdb  <br>
 *      If specified then the initial counter configuration is not loaded from the CDB.
 *
 * <h2>Configuration:</h2>
 * Can only be configured with the command line arguments.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/HLTTriggerCounter - Contains the initial counter configuration.
 *
 * <h2>Performance:</h2>
 * Can run over 2kHz in HLT online system.
 *
 * <h2>Memory consumption:</h2>
 * Negligible.
 *
 * <h2>Output size:</h2>
 * Same order of magnitude as the input data size.
 *
 * \ingroup alihlt_tpc_components
 */
class AliHLTTriggerCounterComponent : public AliHLTProcessor
{
public:
	
	AliHLTTriggerCounterComponent();
	virtual ~AliHLTTriggerCounterComponent();
	
	// Methods inherited from AliHLTComponent:
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual int GetOutputDataTypes(AliHLTComponentDataTypeList& list);
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();
	virtual Int_t DoInit(int argc, const char** argv);
	virtual Int_t DoDeinit();
	
	/**
	 * Returns a reference to the initial counter configuration to use.
	 * This should be used in the configuration file passed to the component with
	 * the -config option to setup the initial counter configuration. The following
	 * is an example of such a config file.
	 * \code
	 *   void SetupConfig() {
	 *     TMap& counters = AliHLTTriggerCounterComponent::InitialCounterConfig();
	 *     counters.Add(new TObjString("TRIGGER-A"), new TObjString("Some input trigger"));
	 *     TObjString* outname = new TObjString("TRIGGER-B");
	 *     outname->SetBit(1<<14);
	 *     counters.Add(outname, new TObjString("Some output trigger"));
	 *   }
	 * \endcode
	 */
	static TMap& InitialCounterConfig()
	{
		fgInitialCounterConfig.SetOwner(kTRUE); // Make sure its the owner of all objects.
		return fgInitialCounterConfig;
	}
	
protected:
	
	// Method inherited from AliHLTProcessor:
	virtual int DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);
	using AliHLTProcessor::DoEvent;
	
private:
	
	/// Implements a ring buffer for old counter and time stamp values.
	class AliRingBuffer : public TObject
	{
	public:
		enum 
		{
			/// The maximum number of time stamp entries held in the ring buffers for rate calculations.
			/// The uncertainty in the rate should be ~ rate/sqrt(kTimeStampEntries)
			kTimeStampEntries = 100
		};
		
		/// Constructor initialises the buffer
		AliRingBuffer(const char* counterName, Double_t startTime) : TObject(), fCounterName(counterName), fPos(0)
		{
			for (size_t i = 0; i < kTimeStampEntries; ++i) {
				fCounterBuffer[i] = 0; fTimeBuffer[i] = startTime;
			}
		}
		
		/// Inherited form TObject. Returns the counter name this buffer is associated to.
		virtual const char* GetName() const { return fCounterName.Data(); }
		
		/// Inherited from TObject. Returns a hash value calculated from the counter name.
		virtual ULong_t Hash() const { return fCounterName.Hash(); }
		
		/// Returns the oldest counter value.
		ULong64_t OldestCounter() const { return fCounterBuffer[fPos]; }
		
		/// Returns the oldest time value.
		Double_t OldestTime() const { return fTimeBuffer[fPos]; }
		
		/**
		 * Increments the buffer. The new time and counter values are put into the
		 * current location to replace the existing values, then the current position
		 * is incremented (and wrapped around if necessary).
		 * \param newCounter  The new counter value to put into the buffer.
		 * \param newTime  The new time value to put into the buffer.
		 */
		void Increment(ULong64_t newCounter, Double_t newTime)
		{
			fCounterBuffer[fPos] = newCounter;
			fTimeBuffer[fPos] = newTime;
			fPos = (fPos+1) % kTimeStampEntries;
		}
		
	private:
		
		TString fCounterName;  // The name of the counter.
		UInt_t fPos;  // Current position in the buffers.
		ULong64_t fCounterBuffer[kTimeStampEntries];  // The counter elements of the buffer.
		Double_t fTimeBuffer[kTimeStampEntries];  // The time elements of the buffer.
	};
	
	// Do not allow copying of this class.
	AliHLTTriggerCounterComponent(const AliHLTTriggerCounterComponent& obj);
	AliHLTTriggerCounterComponent& operator = (const AliHLTTriggerCounterComponent& obj);
	
	/// Updates the counter rate value.
	void UpdateCounterRate(AliHLTTriggerCounters::AliCounter* counter, AliRingBuffer* timeBuf, Double_t newTime);
	
	/// Loads the initial configuration of counters from the given CDB path.
	int LoadConfigFromCDB(const char* cdbPath);
	
	/// Loads the initial configuration of counters from the given configuration file.
	int LoadConfigFromFile(const char* configFile);
	
	/**
	 * Sets the initial counter values from a TMap object containing name description pairs.
	 * \param counters  The list of counters to initialise. This TMap should contain TObjString
	 *     pairs where the key is the name of the counter and the TMap value is the description.
	 *     If the 14th bit is set in the key TObjString then the counter is considered an
	 *     output counter and an input counter otherwise.
	 */
	void SetInitialCounters(const TMap* counters);
	
	double fOutputMultiplier;    //! Multiplier value for the output size estimate.
	AliHLTTriggerCounters fInputCounters; //! Input trigger counters.
	AliHLTTriggerCounters fOutputCounters; //! Output trigger counters.
	THashTable fInputTimes;  //! Cyclic buffer storing the time stamps when the input counters were received.
	THashTable fOutputTimes;  //! Cyclic buffer storing the time stamps when the output counters were received.
	Double_t fLastPublishTime;  //! The last time the counters were pushed back to the output.
	Double_t fPublishPeriod;  //! The time between calls to push back the counters to output.
	
	static const char* fgkConfigCDBPath;  //! CDB configuration path.
	static TMap fgInitialCounterConfig;   //! Initial configuration information for the counters.
	
	ClassDef(AliHLTTriggerCounterComponent, 0)  // Component for counting HLT triggers.
};

#endif // ALIHLTTRIGGERCOUNTERCOMPONENT_H
