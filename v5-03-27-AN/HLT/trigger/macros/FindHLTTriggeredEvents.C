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
 * \file FindHLTTriggeredEvents.C
 * \brief Macro for finding event numbers of HLT triggered events.
 *
 * This macro is used to generate or print a list of event numbers which had
 * a positive HLT trigger decision.
 *
 * The simplest way to run this macro with defaults (which will process any raw
 * data in the current directory) is to run the following command from a terminal
 * shell:
 * \code
 *   > aliroot -b -q $ALICE_ROOT/HLT/trigger/macros/FindHLTTriggeredEvents.C
 * \endcode
 *
 * It is also possible to use as input a file or a chunk sitting on the GRID.
 * \code 
 *   > aliroot -b -q  $ALICE_ROOT/HLT/trigger/macros/FindHLTTriggeredEvents.C'("alien:///alice/data/2010/LHC10h/000137124/raw/10000137124054.10.root")'
 *   > aliroot -b -q  $ALICE_ROOT/HLT/trigger/macros/FindHLTTriggeredEvents.C'("raw://run137124")'
 * \author Artur Szostak <artursz@iafrica.com>
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRawReader.h"
#include "AliHLTPluginBase.h"
#include "AliHLTOUT.h"
#include "AliHLTSystem.h"
#include "AliHLTComponent.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliLog.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "Riostream.h"
#endif

/**
 * Finds all the events in raw data that contains a positive HLT decision.
 *
 * \param [in] dataSource  This is the path to the raw data or the ROOT/DATE file
 *     contining the raw data.
 * \param [out] events  The output list that will be filled with event numbers of
 *     the events that were triggered by HLT for the given data source.
 * \param [in] triggerCode  The HLT trigger code to search for. These must be the same as
 *     the codes defined in the HLT trigger menu used for producing the data source. If
 *     this parameter is set to NULL then all HLT triggered events are marked. (default NULL)
 * \param [in] firstEvent  The event number of the first event to process. (default = 0)
 * \param [in] lastEvent  The event number of the last event to process. If this is
 *     less than firstEvent then the maximum events available in the data are
 *     processed automatically. (default = -1)
 * \param [in] print  Indicates if the found event information should be printed or not.
 * \returns true if there were no problems accessing the data and false otherwise.
 */
bool FindHLTTriggeredEvents(
		const char* dataSource,
		TArrayI& events,
		const char* triggerCode = NULL,
		Int_t firstEvent = 0,
		Int_t lastEvent = -1,
		bool print = true
	)
{
	gSystem->Load("libHLTrec.so");
	// FIXME: Loading the following libraries is a workaround to get rid of warnings from AliHLTOUT.
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libAliHLTUtil.so");
	gSystem->Load("libAliHLTTRD.so");
	gSystem->Load("libAliHLTMUON.so");
	gSystem->Load("libAliHLTTrigger.so");
	
	TString strfile = dataSource;
	if(strfile.BeginsWith("alien://") || strfile.BeginsWith("raw://")) TGrid::Connect("alien");

	// Setup the raw reader and HLTOUT handler.
	AliRawReader* rawReader = AliRawReader::Create(dataSource);
	if (rawReader == NULL)
	{
		cerr << "ERROR: Could not create raw reader for '" << dataSource << "'." << endl;
		return false;
	}
	if (! rawReader->IsRawReaderValid())
	{
		cerr << "ERROR: Raw reader is not valid for '" << dataSource << "'." << endl;
		delete rawReader;
		return false;
	}
	AliHLTOUT* hltout = AliHLTOUT::New(rawReader);
	if (hltout == NULL)
	{
		cerr << "ERROR: Could not create an AliHLTOUT object for '" << dataSource << "'." << endl;
		delete rawReader;
		return false;
	}
	
	// Check if we need to process all events or not.
	if (lastEvent < firstEvent) lastEvent = 0x7FFFFFFF;
	
	bool decodingOK = true;
	
	// Now step through the events.
	for (Int_t currentEvent = firstEvent, eventLoopResult = rawReader->GotoEvent(firstEvent);
	     currentEvent <= lastEvent && eventLoopResult == kTRUE;
	     ++currentEvent, eventLoopResult = rawReader->NextEvent()
	    )
	{
		int result = hltout->Init();
		if (result != 0)
		{
			cerr << "ERROR: could not initialise HLTOUT. Moving to next event." << endl;
			hltout->Reset();
			decodingOK = false;
			continue;
		}
		
		// Go through all the data blocks, looking for the HLT global trigger object.
		for (result = hltout->SelectFirstDataBlock(kAliHLTDataTypeGlobalTrigger);
		     result >= 0;
		     result = hltout->SelectNextDataBlock()
		    )
		{
			TObject* obj = hltout->GetDataObject();
			if (obj == NULL) continue;
			if (obj->IsA()->GetBaseClass("AliHLTGlobalTriggerDecision") != NULL)
			{
				AliHLTGlobalTriggerDecision* decision = (AliHLTGlobalTriggerDecision*)obj;
				bool hltTriggered = false;
				if (triggerCode != NULL)
				{
					hltTriggered = TString(decision->Description()).Contains(triggerCode);
				}
				else
				{
					hltTriggered = decision->Result();
				}
				if (hltTriggered)
				{
					if (print)
					{
						cout << "Event " << currentEvent << " in " << dataSource
							<< " with event ID = " << hltout->EventId()
							<< " (0x" << hex << hltout->EventId() << dec << ")"
							<< " was triggered by HLT with \""
							<< (triggerCode != NULL ? triggerCode : decision->Description())
							<< "\"." << endl;
					}
					events.Set(events.GetSize()+1);
					events[events.GetSize()-1] = currentEvent;
				}
			}
			hltout->ReleaseDataObject(obj);
		}
		
		result = hltout->Reset();
		if (result != 0)
		{
			cerr << "ERROR: could not reset HLTOUT properly." << endl;
			decodingOK = false;
			continue;
		}
	}
	
	AliHLTOUT::Delete(hltout);
	delete rawReader;
	return decodingOK;
}


bool FindHLTTriggeredEvents(const char* dataSource = "./")
{
	TArrayI events;
	return FindHLTTriggeredEvents(dataSource, events, NULL, 0, -1, true);
}
