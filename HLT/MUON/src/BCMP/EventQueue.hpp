////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BCMP_EVENT_QUEUE_HPP
#define dHLT_BCMP_EVENT_QUEUE_HPP

#include "System/Mutex.hpp"
#include "System/MutexCondition.hpp"
#include "BCMP/Packets.hpp"

#include <queue>

namespace dHLT
{
namespace BCMP
{


enum EventType
{
	OnConnectEvent,
	OnDisconnectEvent,
	OnConnectionLostEvent,
	OnMessageEvent,
	OnErrorEvent
};


enum
{
	BAD_EVENT_TYPE         = 0x10010001,
	UNEXPECTED_EVENT       = 0x10010002
};


class UnknownEvent : public Error
{
public:

	virtual const char* Message() const throw()
	{
		return "Unknown BCMP event type.";
	};
	
	virtual Int ErrorCode() const throw()
	{
		return BAD_EVENT_TYPE;
	};
};

class UnexpectedEvent : public Error
{
public:

	virtual const char* Message() const throw()
	{
		return "Unexpected BCMP event type.";
	};
	
	virtual Int ErrorCode() const throw()
	{
		return UNEXPECTED_EVENT;
	};
};


class EventQueue
{
public:
	~EventQueue();

	// This method takes ownership of memory pointed to by message.
	// It should be allocated with AllocateMessage();
	void Push(
			const EventType eventtype, const System::Address& address,
			char* message = NULL, const UInt length = 0
		);

	// Caller must cleanup message with FreeMessage() for the following method calls:
	// TryPop and Pop.
	bool TryPop(
			EventType& eventtype, System::Address& address,
			char*& message, UInt& length
		);

	// timeout in milliseconds
	// Returns false if we timed out.
	bool Pop
			(EventType& eventtype, System::Address& address,
			char*& message, UInt& length, const UInt timeout
		);

	void Pop(
			EventType& eventtype, System::Address& address,
			char*& message, UInt& length
		);

private:

	struct Record
	{
		EventType type;
		System::Address address;
		char* message;
		UInt length;
	};

	std::queue<Record> eventq;
	System::Mutex mutex;
	System::MutexCondition queue_not_empty;
};


} // BCMP
} // dHLT

#endif // dHLT_BCMP_EVENT_QUEUE_HPP
