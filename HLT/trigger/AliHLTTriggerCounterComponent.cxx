// $Id: AliHLTTriggerCounterComponent.cxx 43071 2010-08-25 08:41:44Z richterm $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTTriggerCounterComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   3 Nov 2010
/// @brief  Implementation of the AliHLTTriggerCounterComponent class.
///
/// The AliHLTTriggerCounterComponent is used to count HLT input and output
/// triggers based on the global trigger decisions coming from the global HLT
/// trigger component.

#include "AliHLTTriggerCounterComponent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTTriggerCounters.h"
#include "AliHLTCTPData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "TObjString.h"
#include "TROOT.h"
#include <cstring>
#include <cassert>

ClassImp(AliHLTTriggerCounterComponent)

const char* AliHLTTriggerCounterComponent::fgkConfigCDBPath = "HLT/ConfigHLT/HLTTriggerCounter";
TMap AliHLTTriggerCounterComponent::fgInitialCounterConfig(TCollection::kInitHashTableCapacity, 2);


AliHLTTriggerCounterComponent::AliHLTTriggerCounterComponent() :
	AliHLTProcessor(),
	fOutputMultiplier(0),
	fInputCounters(),
	fOutputCounters(),
	fInputTimes(TCollection::kInitHashTableCapacity, 2),
	fOutputTimes(TCollection::kInitHashTableCapacity, 2),
	fLastPublishTime(-1),
	fPublishPeriod(-1)
{
	// Default constructor.
	
	fInputTimes.SetOwner(kTRUE);
	fOutputTimes.SetOwner(kTRUE);
}


AliHLTTriggerCounterComponent::~AliHLTTriggerCounterComponent()
{
	// Default destructor.
	
}


const char* AliHLTTriggerCounterComponent::GetComponentID()
{
	// Returns the component ID.
	return "HLTTriggerCounter";
}


void AliHLTTriggerCounterComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of input data types that are handled.
	list.push_back(kAliHLTDataTypeGlobalTrigger);
	list.push_back(kAliHLTDataTypeTriggerDecision);
}


AliHLTComponentDataType AliHLTTriggerCounterComponent::GetOutputDataType()
{
	// Returns kAliHLTMultipleDataType.
	return kAliHLTMultipleDataType;
}


int AliHLTTriggerCounterComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
	// Returns the list of output data block types generated.
	list.push_back(kAliHLTDataTypeInputTriggerCounters);
	list.push_back(kAliHLTDataTypeOutputTriggerCounters);
	return int(list.size());
}


void AliHLTTriggerCounterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	// Returns the buffer size requirements.
	constBase = 1024*16;
	inputMultiplier = fOutputMultiplier;
}


AliHLTComponent* AliHLTTriggerCounterComponent::Spawn()
{
	// Creates a new instance of the component.
	return new AliHLTTriggerCounterComponent;
}


Int_t AliHLTTriggerCounterComponent::DoInit(int argc, const char** argv)
{
	// Initialises the data checker component from the command line.
	
	HLTInfo("Starting HLT trigger counter component.");
	
	const char* configFileName = NULL;
	fInputCounters.Clear();
	fOutputCounters.Clear();
	fInputTimes.Clear();
	fOutputTimes.Clear();
	fLastPublishTime = -1;
	fPublishPeriod = -1;
	bool loadCDBObject = true;
	
	for (int i = 0; i < argc; ++i)
	{
		if (strcmp(argv[i], "-config") == 0)
		{
			if (configFileName != NULL)
			{
				HLTWarning("The configuration macro was already specified."
					" Will replace previous value given by -config."
				);
			}
			if (argc <= i+1)
			{
				HLTError("The configuration macro filename was not specified." );
				return -EINVAL;
			}
			configFileName = argv[i+1];
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-publishperiod") == 0)
		{
			if (fPublishPeriod != -1)
			{
				HLTWarning("The publish period was already specified."
					" Will replace previous value given by -publishperiod."
				);
			}
			if (argc <= i+1)
			{
				HLTError("Publish period value was not specified for -publishperiod.");
				return -EINVAL;
			}
			char* err = NULL;
			errno = 0;
			double tmpnum = strtod(argv[i+1], &err);
			if (err == NULL or *err != '\0')
			{
				HLTError("Cannot convert '%s' to a floating point value.", argv[i+1]);
				return false;
			}
			if (errno == ERANGE)
			{
				HLTError("The specified value '%s' is out of range.", argv[i+1]);
				return false;
			}
			fPublishPeriod = (tmpnum < 0 ? -1 : tmpnum);
			i++;
			continue;
		}
		
		if (strcmp(argv[i], "-skipcdb") == 0)
		{
			loadCDBObject = false;
			continue;
		}
		
		HLTError("Unknown option '%s'.", argv[i]);
		return -EINVAL;
	} // for loop
	
	if (configFileName != NULL)
	{
		int result = LoadConfigFromFile(configFileName);
		if (result != 0) return result;
	}
	else if (loadCDBObject)
	{
		int result = LoadConfigFromCDB(fgkConfigCDBPath);
		if (result != 0) return result;
	}
	
	return 0;
}


Int_t AliHLTTriggerCounterComponent::DoDeinit()
{
	// Cleans up the component.
	
	HLTInfo("Stopping HLT trigger counter component.");
	fInputCounters.Clear();
	fOutputCounters.Clear();
	fInputTimes.Clear();
	fOutputTimes.Clear();
	return 0;
}


int AliHLTTriggerCounterComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
	// Finds the global trigger objects and adds the triggers to the counters.
	
	fInputCounters.UpdateTimeStamp();
	fOutputCounters.UpdateTimeStamp();
	Double_t inputTime = fInputCounters.TimeStamp().AsDouble();
	Double_t outputTime = fOutputCounters.TimeStamp().AsDouble();
	
	const TObject* obj = GetFirstInputObject(kAliHLTDataTypeGlobalTrigger, "AliHLTGlobalTriggerDecision");
	while (obj != NULL)
	{
		const AliHLTGlobalTriggerDecision* decision = dynamic_cast<const AliHLTGlobalTriggerDecision*>(obj);
		if (decision != NULL)
		{
			HLTDebug("Received trigger decision object of type AliHLTGlobalTriggerDecision.");
			// Tokenise the global trigger description string which contains a
			// list of the triggers that were fired and increment the corresponding
			// counters.
			TString names = decision->Description();
			Ssiz_t from = 0;
			TString token;
			while (names.Tokenize(token, from, ","))
			{
				TObject* cntobj = fOutputCounters.FindObject(token.Data());
				if (cntobj != NULL)
				{
					AliHLTTriggerCounters::AliCounter* counter = static_cast<AliHLTTriggerCounters::AliCounter*>(cntobj);
					counter->Increment();
					UpdateCounterRate(
							counter,
							static_cast<AliRingBuffer*>( fOutputTimes.FindObject(token.Data()) ),
							outputTime
						);
				}
				else
				{
					fOutputCounters.Add(token.Data(), "New trigger output counter found during the run.", 1, 1);
					fOutputTimes.Add(new AliRingBuffer(token.Data(), outputTime));
				}
			}
			
			// Now add the CTP triggers.
			const AliHLTCTPData* ctp = dynamic_cast<const AliHLTCTPData*>( decision->InputObjects().FindObject("AliHLTCTPData") );
			if (ctp != NULL)
			{
				const TArrayL64& counters = ctp->Counters();
				for (Int_t i = 0; i < counters.GetSize(); ++i)
				{
					const char* ctpName = ctp->Name(i);
					// Check if CTP counter is initialised and skip if not.
					if (strcmp(ctpName, "AliHLTReadoutList") == 0 and counters[i] == 0) continue;
					TObject* cntobj = fInputCounters.FindObject(ctpName);
					if (cntobj != NULL)
					{
						AliHLTTriggerCounters::AliCounter* counter = static_cast<AliHLTTriggerCounters::AliCounter*>(cntobj);
						counter->Counter(counters[i]);
						UpdateCounterRate(
								counter,
								static_cast<AliRingBuffer*>( fInputTimes.FindObject(ctpName) ),
								inputTime
							);
					}
					else
					{
						fInputCounters.Add(
								ctpName,
								"New CTP trigger input counter found during the run.",
								Double_t(counters[i]),
								counters[i]
							);
						fInputTimes.Add(new AliRingBuffer(ctpName, inputTime));
					}
				}
			}
			
			// Add the list of input triggers.
			for (Int_t i = 0; i < decision->NumberOfTriggerInputs(); ++i)
			{
				const AliHLTTriggerDecision* input = decision->TriggerInput(i);
				if (input == NULL or not input->Result()) continue;
				TObject* cntobj = fInputCounters.FindObject(input->Name());
				if (cntobj != NULL)
				{
					AliHLTTriggerCounters::AliCounter* counter = static_cast<AliHLTTriggerCounters::AliCounter*>(cntobj);
					counter->Increment();
					counter->SetBit(BIT(14), true);  // mark counter as updated
					UpdateCounterRate(
							counter,
							static_cast<AliRingBuffer*>( fInputTimes.FindObject(input->Name()) ),
							inputTime
						);
				}
				else
				{
					fInputCounters.Add(input->Name(), "New trigger input counter found during the run.", 1, 1);
					fInputCounters.GetCounterN(fInputCounters.NumberOfCounters()-1).SetBit(BIT(14), true); // mark counter as updated
					fInputTimes.Add(new AliRingBuffer(input->Name(), inputTime));
				}
			}
		}
		obj = GetNextInputObject();
	}
	
	// Look for the individual trigger decision blocks to add if they are available.
	obj = GetFirstInputObject(kAliHLTDataTypeTriggerDecision, "AliHLTTriggerDecision");
	while (obj != NULL)
	{
		const AliHLTTriggerDecision* decision = dynamic_cast<const AliHLTTriggerDecision*>(obj);
		if (decision != NULL and decision->Result())
		{
			HLTDebug("Received trigger decision object of type AliHLTTriggerDecision.");
			TObject* cntobj = fInputCounters.FindObject(decision->Name());
			if (cntobj != NULL)
			{
				AliHLTTriggerCounters::AliCounter* counter = static_cast<AliHLTTriggerCounters::AliCounter*>(cntobj);
				if (not counter->TestBit(BIT(14)))  // Only update if marked as not updated.
				{
					counter->Increment();
					UpdateCounterRate(
							counter,
							static_cast<AliRingBuffer*>( fInputTimes.FindObject(decision->Name()) ),
							inputTime
						);
				}
			}
			else
			{
				fInputCounters.Add(decision->Name(), "New trigger input counter found during the run.", 1, 1);
				fInputTimes.Add(new AliRingBuffer(decision->Name(), inputTime));
			}
		}
		obj = GetNextInputObject();
	}
	// Reset bit 14 which is used temporarily to mark updated counters.
	for (UInt_t i = 0; i < fInputCounters.NumberOfCounters(); ++i)
	{
		fInputCounters.GetCounterN(i).SetBit(BIT(14), false);
	}
	
	Double_t now = TTimeStamp();
	if (fLastPublishTime == -1) fLastPublishTime = now;
	if (now - fLastPublishTime > fPublishPeriod)
	{
		fLastPublishTime = now;
		bool inputCountersNotPushed = PushBack(&fInputCounters, kAliHLTDataTypeInputTriggerCounters) != 0;
		bool outputCountersNotPushed = PushBack(&fOutputCounters, kAliHLTDataTypeOutputTriggerCounters) != 0;
		if (inputCountersNotPushed or outputCountersNotPushed)
		{
			fOutputMultiplier = (fOutputMultiplier == 0 ? 1 : fOutputMultiplier*2);
			return -ENOSPC;
		}
	}
	return 0;
}


void AliHLTTriggerCounterComponent::UpdateCounterRate(
		AliHLTTriggerCounters::AliCounter* counter, AliRingBuffer* timeBuf, Double_t newTime
	)
{
	// Updates the counter's rate value.
	
	assert(timeBuf != NULL);
	Double_t dt = newTime - timeBuf->OldestTime();
	Double_t rate = (dt != 0 ? (counter->Counter() - timeBuf->OldestCounter()) / dt : 0.);
	counter->Rate(rate);
	timeBuf->Increment(counter->Counter(), newTime);
}


int AliHLTTriggerCounterComponent::LoadConfigFromCDB(const char* cdbPath)
{
	// Loads the initial configuration of counters from the given CDB path.
	
	HLTDebug("Trying to load trigger menu from '%s'.", cdbPath);
	
	if (AliCDBManager::Instance() == NULL)
	{
		HLTError("CDB manager object not found.");
		return -EIO;
	}
	AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
	if (store == NULL)
	{
		HLTError("Could not get the the default storage for the CDB.");
		return -EIO;
	}
	Int_t version = store->GetLatestVersion(cdbPath, GetRunNo());
	if (version < 0)
	{
		HLTError("Could not find an entry in the CDB for \"%s\".", cdbPath);
		return -EIO;
	}
	Int_t subVersion = store->GetLatestSubVersion(cdbPath, GetRunNo(), version);
	AliCDBEntry* entry = AliCDBManager::Instance()->Get(cdbPath, GetRunNo(), version, subVersion);
	if (entry == NULL)
	{
		HLTError("Could not get the CDB entry for \"%s\".", cdbPath);
		return -EIO;
	}
	TObject* obj = entry->GetObject();
	if (obj == NULL)
	{
		HLTError("Configuration object for \"%s\" is missing.", cdbPath);
		return -ENOENT;
	}
	if (obj->IsA() != TMap::Class())
	{
		HLTError("Wrong type for configuration object in \"%s\". Found a %s but we expect a TMap.",
			cdbPath, obj->ClassName()
		);
		return -EPROTO;
	}
	TMap* counters = static_cast<TMap*>(obj);
	SetInitialCounters(counters);
	return 0;
}


int AliHLTTriggerCounterComponent::LoadConfigFromFile(const char* configFile)
{
	// Loads the initial configuration of counters from the given configuration file.
	
	TString cmd = ".x ";
	cmd += configFile;
	gROOT->ProcessLine(cmd);
	SetInitialCounters(&fgInitialCounterConfig);
	return 0;
}


void AliHLTTriggerCounterComponent::SetInitialCounters(const TMap* counters)
{
	// Sets the initial counter values from TMap objects containing name description pairs.
	
	Double_t now = TTimeStamp();
	TMapIter next(counters);
	TObject* key = NULL;
	while ((key = next()) != NULL)
	{
		TObject* value = counters->GetValue(key);
		if (value == NULL) continue;
		if (key->TestBit(BIT(14)))
		{
			fOutputCounters.Add(key->GetName(), value->GetName());
			fOutputTimes.Add(new AliRingBuffer(key->GetName(), now));
		}
		else
		{
			fInputCounters.Add(key->GetName(), value->GetName());
			fInputTimes.Add(new AliRingBuffer(key->GetName(), now));
		}
	}
}
