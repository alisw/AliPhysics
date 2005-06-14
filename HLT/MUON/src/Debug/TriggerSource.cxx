////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Debug/TriggerSource.hpp"

#include "Debug/print.hpp"
#include <iostream>
using std::endl;
using std::cout;

namespace dHLT
{
namespace Debug
{


TriggerSource::TriggerSource()
{
	framework = NULL;
	currentevent = NULL;
	currentblock = NULL;
};

TriggerSource::~TriggerSource()
{
	Clear();
};


void TriggerSource::NewEvent(const EventID id)
{
	EventBlock newevent;
	newevent.eventid = id;
	event.push_back(newevent);
	currentevent = &event[event.size() -1];
	currentblock = NULL;
};

void TriggerSource::NewTriggerBlock()
{
	Assert( currentevent != NULL );
	currentevent->block.push_back( DataBlock() );
	currentblock = &currentevent->block[currentevent->block.size() -1];
};

void TriggerSource::AddTrigger(const TriggerRecord& t)
{
	Assert( currentblock != NULL );
	currentblock->trigger.push_back(t);
};

void TriggerSource::Clear()
{
	event.clear();
	currentevent = NULL;
	currentblock = NULL;
};


void TriggerSource::Process()
{
	Assert( framework != NULL );
	
	for (UInt i = 0; i < event.size(); i++)
	{
		for (UInt j = 0; j < event[i].block.size(); j++)
		{
			UInt count = event[i].block[j].trigger.size();
			UInt size = count * sizeof(TriggerRecord);
			TriggerRecord* triggerblock = framework->AllocateTriggerBlock(size);
			for (UInt k = 0; k < count; k++)
				triggerblock[k] = event[i].block[j].trigger[k];

			framework->ReturnTriggers( event[i].eventid, triggerblock, count );
			// Note: triggerblock is deleted by framework.
		};
		framework->EndOfTriggers(event[i].eventid);
	};
};


void TriggerSource::Dump()
{
	cout << "==================== Triggers ====================" << endl;
	for (UInt i = 0; i < event.size(); i++)
	{
		cout << "event = " << event[i].eventid << endl;
		for (UInt j = 0; j < event[i].block.size(); j++)
		{
			cout << "\tblock #" << j+1 << endl;
			for (UInt k = 0; k < event[i].block[j].trigger.size(); k++)
				cout << "\t\t" << event[i].block[j].trigger[k] << endl;
		};
	};
};


} // Debug
} // dHLT
