/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id: $

/**
 * \ingroup macros
 * \file DumpGlobalTrigger.C
 * \brief Macro for dumping HLT global trigger decisions.
 *
 * This macro is used to dump to screen the HLT global trigger decisions encoded
 * in AliHLTGlobalTriggerDecision objects found in the HLT raw data. The raw data
 * should be packed in a ROOT or DATE file, or in DDL data files in raw*
 * directories. An appropriate AliRawReader is created to read the data.
 *
 * The simplest way to run this macro with defaults is to run the following
 * command from a terminal shell:
 * \code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/macros/DumpGlobalTrigger.C
 * \endcode
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRawReader.h"
#include "AliHLTPluginBase.h"
#include "AliHLTOUT.h"
#include "AliHLTSystem.h"
#include "AliHLTComponent.h"
#include "AliLog.h"
#include "TSystem.h"
#include "TFile.h"
#include "Riostream.h"
#endif

/**
 * Dumps the AliHLTGlobalTriggerDecision objects found in HLT output data.
 *
 * \param dataSource  This is the path to the raw data or the ROOT/DATE file
 *     contining the raw data. (default is the current directory).
 * \param firstEvent  The event number of the first event to process. (default = 0)
 * \param lastEvent  The event number of the last event to process. If this is
 *     less than firstEvent then it is set to maximum events available
 *     automatically. (default = -1)
 * \param output  Specifies the name of a ROOT output file. This file will be
 *     generated and the objects written to it. If the value is NULL then
 *     no output is written to file. (default = NULL)
 * \param debug  Specifies if full debug messages should be printed.
 */
void DumpGlobalTrigger(
		const char* dataSource = "./",
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		const char* output = NULL,
		bool debug = false
	)
{
	if (debug)
	{
		AliLog::SetModuleDebugLevel("HLT", AliLog::kMaxType);
		AliHLTSystem* sys = AliHLTPluginBase::GetInstance();
		sys->SetGlobalLoggingLevel(kHLTLogAll);
	}
	
	gSystem->Load("libHLTrec.so");
	TFile* file = NULL;
	if (output != NULL)
	{
		file = new TFile(output, "RECREATE");
		if (file == NULL)
		{
			cerr << "ERROR: Could not create file '" << output << "'." << endl;
			return;
		}
	}
	
	// Setup the raw reader and HLTOUT handler.
	AliRawReader* rawReader = AliRawReader::Create(dataSource);
	if (rawReader == NULL)
	{
		cerr << "ERROR: Could not create raw reader for '" << dataSource << "'." << endl;
		if (file != NULL) delete file;
		return;
	}
	if (! rawReader->IsRawReaderValid())
	{
		cerr << "ERROR: Raw reader is not valid for '" << dataSource << "'." << endl;
		delete rawReader;
		if (file != NULL) delete file;
		return;
	}
	AliHLTOUT* hltout = AliHLTOUT::New(rawReader);
	if (hltout == NULL)
	{
		cerr << "ERROR: Could not create an AliHLTOUT object for '" << dataSource << "'." << endl;
		delete rawReader;
		if (file != NULL) delete file;
		return;
	}
	
	// Make sure that the lastEvent is greater than firstEvent.
	if (lastEvent < firstEvent) lastEvent = rawReader->GetNumberOfEvents();
	if (lastEvent < firstEvent) lastEvent = firstEvent;
	
	// Need to call NextEvent once here or we will start at the wrong event.
	if (! rawReader->NextEvent())
	{
		cout << "No events found in '" << dataSource << "'." << endl;
		AliHLTOUT::Delete(hltout);
		delete rawReader;
		if (file != NULL) delete file;
		return;
	}
	
	// Now step through the events.
	for (int i = 0; i < firstEvent; i++) rawReader->NextEvent();
	for (int i = firstEvent; i <= lastEvent; i++)
	{
		int result = hltout->Init();
		if (result != 0)
		{
			cerr << "ERROR: could not initialise HLTOUT." << endl;
			hltout->Reset();
			continue;
		}
		cout << "#################### Event " << i << " in " << dataSource
			<< " has event ID = " << hltout->EventId()
			<< " (0x" << hex << hltout->EventId() << dec << ")"
			<< " ####################" << endl;
		
		for (result = hltout->SelectFirstDataBlock();
		     result >= 0;
		     result = hltout->SelectNextDataBlock()
		    )
		{
			AliHLTComponentDataType dt;
			AliHLTUInt32_t spec = 0;
			hltout->GetDataBlockDescription(dt, spec);
			TObject* obj = hltout->GetDataObject();
			if (obj == NULL) continue;
			if (obj->IsA()->GetBaseClass("AliHLTGlobalTriggerDecision") != NULL)
			{
				if (dt != kAliHLTDataTypeGlobalTrigger)
				{
					cerr << "WARNING: Found an AliHLTGlobalTriggerDecision object in a data block of type '"
						<< AliHLTComponent::DataType2Text(dt).c_str()
						<< "' but expected '"
						<< AliHLTComponent::DataType2Text(kAliHLTDataTypeGlobalTrigger).c_str()
						<< "'." << endl;
				}
				if (file != NULL)
				{
					obj->Write(
						Form("HLTGlobalDecision_event_0x%llX", hltout->EventId()),
						TObject::kOverwrite
					);
				}
				obj->Print();
			}
			hltout->ReleaseDataObject(obj);
		}
		
		result = hltout->Reset();
		if (result != 0)
		{
			cerr << "ERROR: could not reset HLTOUT." << endl;
			hltout->Reset();
			continue;
		}
		rawReader->NextEvent();
	}
	
	AliHLTOUT::Delete(hltout);
	delete rawReader;
	if (file != NULL) delete file;
}
