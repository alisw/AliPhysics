////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "BCMP/EventQueue.hpp"
#include "Utils.hpp"
#include "Error.hpp"

using namespace dHLT::System;


// TODO: replace stl::queue with internal Queue class.
//       the stl on redhat seems to be doing very strange things with memory
//       which can not be easily validated.

namespace dHLT
{
namespace BCMP
{


EventQueue::~EventQueue()
{
	// Must cleanup any memory that was not handled yet.
	while (not eventq.empty())
	{
		Record rec = eventq.front();
		if (rec.message != NULL)
			FreeMessage( (char*)rec.message );
		eventq.pop();
	};
};


void EventQueue::Push(
		const EventType eventtype, const System::Address& address,
		char* message, const UInt length
	)
{
	mutex.Lock();
	try
	{
		Record rec;
		rec.type = eventtype;
		rec.address = address;
		rec.message = message;
		rec.length = length;
		eventq.push(rec);
		queue_not_empty.Signal(false);
	}
	finally
	(
		mutex.Unlock();
	);
};


bool EventQueue::TryPop(
		EventType& eventtype, System::Address& address,
		char*& message, UInt& length
	)
{
	if ( not mutex.TryLock() ) return false;
	try
	{
		// Wait for something to be added to the queue.
		while ( eventq.empty() ) queue_not_empty.Wait(mutex);

		const Record& rec = eventq.front();
		eventtype = rec.type;
		address = rec.address;
		message = rec.message;
		length = rec.length;
		eventq.pop();
	}
	finally
	(
		mutex.Unlock();
	);
	return true;
};


bool EventQueue::Pop(
		EventType& eventtype, System::Address& address,
		char*& message, UInt& length, const UInt timeout
	)
{
	bool ready = true;
	mutex.Lock();
	try
	{
		// Wait for something to be added to the queue.
		while ( eventq.empty() and ready )
			ready = queue_not_empty.Wait(mutex, timeout);

		if (ready)
		{
			const Record& rec = eventq.front();
			eventtype = rec.type;
			address = rec.address;
			message = rec.message;
			length = rec.length;
			eventq.pop();
		};
	}
	finally
	(
		mutex.Unlock();
	);
	return ready;
};


void EventQueue::Pop(
		EventType& eventtype, System::Address& address,
		char*& message, UInt& length
	)
{
	mutex.Lock();
	try
	{
		// Wait for something to be added to the queue.
		while ( eventq.empty() ) queue_not_empty.Wait(mutex);

		const Record& rec = eventq.front();
		eventtype = rec.type;
		address = rec.address;
		message = rec.message;
		length = rec.length;
		eventq.pop();
	}
	finally
	(
		mutex.Unlock();
	);
};


} // BCMP
} // dHLT
