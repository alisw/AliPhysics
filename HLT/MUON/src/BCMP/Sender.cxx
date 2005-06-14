////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "BCMP/Sender.hpp"
#include "Utils.hpp"
#include "Error.hpp"
#include <asm/errno.h>

#include <string>
#include <iostream>
using std::cerr;
using std::endl;

using namespace dHLT::System;


#define WAIT_LIMIT          1   // in milliseconds

// Total timeout in SenderLoop = TIMEOUT_SEC + TIMEOUT_MICROSEC / 1e6 seconds
#define TIMEOUT_SEC         0
#define TIMEOUT_MICROSEC 1000

// TODO: replace stl::vector with internal class.
//       the stl on redhat seems to be doing very strange things with memory
//       which can not be easily validated.


namespace dHLT
{
namespace BCMP
{


UInt Sender::default_retry_limit = 4;

///////////////////////////////////////////////////////////////////////////////

Sender::PacketQueue::~PacketQueue()
{
	// Must handle any memory that was not handled yet.
	while (not pq.empty())
	{
		Record rec = pq.front();
		FreeMessage( (char*)rec.packet );
		pq.pop();
	};
};


void Sender::PacketQueue::Push(DataPacket* packet, const UInt size)
{
	mutex.Lock();
	try
	{
		Record rec;
		rec.packet = packet;
		rec.packetsize = size;
		pq.push(rec);
		queue_not_empty.Signal(false);
	}
	finally
	(
		mutex.Unlock();
	);
};


bool Sender::PacketQueue::Pop(DataPacket*& packet, UInt& size)
{
	bool ready = true;
	mutex.Lock();
	try
	{
		// Wait for something to be added to the queue.
		while ( pq.empty() and ready )
			ready = queue_not_empty.Wait(mutex, WAIT_LIMIT);

		if (ready) 
		{
			const Record& rec = pq.front();
			packet = rec.packet;
			size = rec.packetsize;
			pq.pop();
		};
	}
	finally
	(
		mutex.Unlock();
	);
	return ready;
};

///////////////////////////////////////////////////////////////////////////////

bool Sender::IPList::Add(const UInt ip)
{
	// Make sure the ip is not in the list already.
	std::vector<UInt>::iterator cip;
	for (cip = list.begin(); cip != list.end(); cip++)
	{
		if (*cip == ip) return false;
	};

	list.push_back(ip);
	return true;
};


bool Sender::IPList::Remove(const UInt ip)
{
	std::vector<UInt>::iterator cip;
	for (cip = list.begin(); cip != list.end(); cip++)
	{
		if (*cip == ip)
		{
			list.erase(cip);
			return true;
		};
	};
	return false;
};


bool Sender::IPList::Empty() const
{
	return list.size() == 0;
};


void Sender::IPList::Clear()
{
	list.erase(list.begin(), list.end());
};


void Sender::IPList::Copy(const IPList& rhs)
{
	list.resize(rhs.list.size());
	for (UInt i = 0; i < list.size(); i++)
		list[i] = rhs.list[i];
};


bool Sender::IPList::Contains(const UInt ip) const
{
	std::vector<UInt>::const_iterator cip;
	for (cip = list.begin(); cip != list.end(); cip++)
	{
		if (*cip == ip)	return true;
	};
	return false;
};


UInt Sender::IPList::Count() const
{
	return list.size();
};


UInt Sender::IPList::operator [] (const UInt index) const
{
	Assert( index < Count() );
	return list[index];
};


void Sender::IPList::Dump()
{
	DebugMsg(2, "hostip");
	DebugMsg(2, "-----------");
	std::vector<UInt>::iterator cip;
	for (cip = list.begin(); cip != list.end(); cip++)
	{
		DebugMsg(2, Address::AsString(*cip) );
	};
};

///////////////////////////////////////////////////////////////////////////////

Int Sender::SendThread::Execute()
{
	Int returncode = 0;
	try
	{
		Assert( sender != NULL );
		sender->SendLoop();
	}
	catch (const System::Error& e)
	{
		if (e.ErrorCode() != EINTR)
		{
			cerr << "<SendThread> Error [" << e.ErrorCode() << "]: " << e << endl;
			returncode = e.ErrorCode();
		};
	}
	catch (const Error& e)
	{
		cerr << "<SendThread> Error [" << e.ErrorCode() << "]: " << e << endl;
		returncode = e.ErrorCode();
	}
	catch (...)
	{
		cerr << "<SendThread> Unknown exception!" << endl;
		returncode = -1;
	};
	return returncode;
};


Int Sender::ReceiveThread::Execute()
{
	Int returncode = 0;
	try
	{
		Assert( sender != NULL );
		sender->ReceiveLoop();
	}
	catch (const System::Error& e)
	{
		if (e.ErrorCode() != EINTR)
		{
			cerr << "<ReceiveThread> Error [" << e.ErrorCode() << "]: " << e << endl;
			returncode = e.ErrorCode();
		};
	}
	catch (const Error& e)
	{
		cerr << "<ReceiveThread> Error [" << e.ErrorCode() << "]: " << e << endl;
		returncode = e.ErrorCode();
	}
	catch (...)
	{
		cerr << "ReceiveThread> Unknown exception!" << endl;
		returncode = -1;
	};
	return returncode;
};

///////////////////////////////////////////////////////////////////////////////

Sender::Sender(const UShort port)
{
	eventhandler = NULL;
	Init(port);
};


Sender::Sender(EventHandler* handler, const UShort port)
{
	eventhandler = handler;
	Init(port);
};


void Sender::Init(const UShort port)
{
	// Fill the packet ID stack with available packet ID values.
	for (UShort i = 0; i < 256; i++)
	{
		available_id.Push(i);
		busypackets[i].packet = NULL;
	};

	sock.AllowBroadcast();
	sock.Bind( Address::GetHostAddress(0) );
	broadcastaddr = sock.BroadcastAddress();
	broadcastaddr.Port(port);
	DebugMsg(1, "Broadcasting on: " << broadcastaddr);
	
	sendthread.sender = this;
	receivethread.sender = this;
	sendthread.Start();
	receivethread.Start();
};


Sender::~Sender()
{
	Terminate();
	sendthread.WaitForThread();
	receivethread.WaitForThread();
	
	// Cleanup memory only after the send and receive threads are finished.
	for (UShort i = 0; i < 256; i++)
	{
		if (busypackets[i].packet != NULL)
		{
			FreeMessage( (char*)busypackets[i].packet );
		};
	};
};


void Sender::Send(const char* message, const UInt length)
{
	UInt size = length + sizeof(UChar)*2;
	char* buffer = AllocateMessage(size);
	try
	{
		DataPacket* packet = (DataPacket*) buffer;
		packet->type = TypeData;
		memcpy(&packet->data, message, length);
		packetqueue.Push(packet, size);
	}
	catch (...)
	{
		// Cleanup memory on exception.
		FreeMessage(buffer);
		throw;
	};
};


bool Sender::TryHandleEvents()
{
	EventType eventtype;
	Address address;
	char* message;
	UInt length;
	bool gotmessage;
	while ( gotmessage = eventqueue.TryPop(eventtype, address, message, length) )
	{
		if (eventhandler == NULL) continue;
		switch (eventtype)
		{
		case OnConnectEvent:        eventhandler->OnConnect(address); break;
		case OnDisconnectEvent:     eventhandler->OnDisconnect(address); break;
		case OnConnectionLostEvent: eventhandler->OnConnectionLost(address); break;
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		case OnMessageEvent:        throw UnexpectedEvent();
		default:                    throw UnknownEvent();
		};
	};
	return gotmessage;
};


bool Sender::HandleEvents(const UInt timeout)
{
	EventType eventtype;
	Address address;
	char* message;
	UInt length;
	bool gotmessage;
	while ( gotmessage = eventqueue.Pop(eventtype, address, message, length, timeout) )
	{
		if (eventhandler == NULL) continue;
		switch (eventtype)
		{
		case OnConnectEvent:        eventhandler->OnConnect(address); break;
		case OnDisconnectEvent:     eventhandler->OnDisconnect(address); break;
		case OnConnectionLostEvent: eventhandler->OnConnectionLost(address); break;
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		case OnMessageEvent:        throw UnexpectedEvent();
		default:                    throw UnknownEvent();
		};
	};
	return gotmessage;
};


void Sender::HandleEvents()
{
	EventType eventtype;
	Address address;
	char* message;
	UInt length;
	do
	{
		eventqueue.Pop(eventtype, address, message, length);
		if (eventhandler == NULL) continue;
		switch (eventtype)
		{
		case OnConnectEvent:        eventhandler->OnConnect(address); break;
		case OnDisconnectEvent:     eventhandler->OnDisconnect(address); break;
		case OnConnectionLostEvent: eventhandler->OnConnectionLost(address); break;
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		case OnMessageEvent:        throw UnexpectedEvent();
		default:                    throw UnknownEvent();
		};
	}
	while (not terminate);
};


void Sender::SendLoop()
{
	SendRegisterRequest();

	terminate = false;
	bool running = true;
	bool gotpacket = false;
	DataPacket* packet;
	UInt size;
	
	// Enter processing loop.
	do
	{
		try
		{
			// Fetch the next packet to send. Will wait for one to arrive if the
			// queue is empty. gotpacket will be false if we timed out.
			gotpacket = packetqueue.Pop(packet, size);

			mutex.Lock();
			try
			{
				if (gotpacket)
				{
					bool ready = true;

					// If there are no ID values or no recievers available then wait.
					while ( (available_id.Empty() or receivers.Empty()) and ready )
						ready = queues_not_empty.Wait(mutex, WAIT_LIMIT);

					if (ready)
					{
						UChar packetid = available_id.Pop();
						packet->packetid = packetid;
						stilltoack[packetid].Copy(receivers);
						retrylimit[packetid] = default_retry_limit;
						busypackets[packetid].packet = packet;
						busypackets[packetid].size = size;

						DebugMsg(1, "send");
						sock.TargetAddress(broadcastaddr);
						sock.Send((char*)packet, size);
					};
				};
				// Make sure we only terminate when the terminate flag is set and all
				// packets have been processed.
				running = not (terminate and available_id.Full());
			}
			finally
			(
				mutex.Unlock();
			);
		}
		catch (const System::Error& e)
		{
			if (e.ErrorCode() != EINTR) throw;
		};
	}
	while (running);
	
	if (gotpacket)
	{
		// If we got a packet from the queue but promptly terminated then
		// we need to release the memory allocated, because we took ownership
		// of the packet's memory.
		FreeMessage((char*)packet);
	};
	
	canterminate = true;  // allow ReceiveLoop() to terminate.
};


void Sender::ReceiveLoop()
{
	// Enter processing loop.
	canterminate = false;
	bool running = true;
	do
	{
		try
		{
			bool socket_readable = sock.WaitUntilReadable(TIMEOUT_SEC, TIMEOUT_MICROSEC);

			mutex.Lock();
			try
			{
				if (socket_readable)
					ReadFromSocket();
				else
					HandleTimeout();

				// Make sure we only terminate when the terminate flag is set and all
				// packets have been processed.
				running = not (canterminate and available_id.Full());
			}
			finally
			(
				mutex.Unlock();
			);
		}
		catch (const System::Error& e)
		{
			if (e.ErrorCode() != EINTR) throw;
		};
	}
	while (running);
	
	SendRelease();
};


void Sender::ReadFromSocket()
{
	try
	{
		do
		{
			UInt length = sock.BytesAvailable();
			char* message = AllocateMessage(length);
			try
			{
				sock.Receive(message, length);
				ProcessMessage(message, length);
			}
			finally
			(
				// Note: message pointer gets set to NULL if it is not
				// supposed to be deleted.
				if (message != NULL) FreeMessage(message);
			);
		}
		while (sock.CanRead());
	}
	catch (const SocketError& e)
	{
		// If we get the WSAECONNRESET error that means the client aborted.
		// The previous send was unsuccessfull.
#ifdef WIN32
		if (e.ErrorCode() == WSAECONNRESET)
		{
			receivers.Remove(sock.TargetAddress().Ipv4Address());
			OnConnectionLost(sock.TargetAddress());
		}
		else
#endif // WIN32
			throw;
	};
};


void Sender::OnConnect(const Address& address)
{
	DebugMsg(1, "Sender::OnConnect( " << address << " )");
	DebugCode( receivers.Dump(); );
	eventqueue.Push(OnConnectEvent, address);
};


void Sender::OnDisconnect(const Address& address)
{
	DebugMsg(1, "Sender::OnDisconnect( " << address << " )");
	DebugCode( receivers.Dump(); )
	eventqueue.Push(OnDisconnectEvent, address);
};


void Sender::OnConnectionLost(const Address& address)
{
	DebugMsg(1, "Sender::OnConnectionLost( " << address << " )");
	eventqueue.Push(OnConnectionLostEvent, address);
};



void Sender::ProcessMessage(char*& message, const UInt length)
{
	if (length < sizeof(UChar))
	{
		SendError(1, "Malformed packet.");
		return;
	};

	PacketType type = (PacketType)(*message);
	DebugMsg(1, "Received: " << type << " packet type.");
	switch (type)
	{
	case TypeRegister:
		SendError(3, "Unexpected register message received.");
		break;

	case TypeConfirm:
		receivers.Add( sock.SourceAddress().Ipv4Address() );
		queues_not_empty.Signal(false);
		OnConnect(sock.SourceAddress());
		break;

	case TypeData:
		SendError(5, "Unexpected data message received.");
		break;

	case TypeAcknowledge:
		ProcessAcknowledge(message, length);
		break;

	case TypeRelease:
		SendError(7, "Unexpected release message received.");
		break;

	case TypeQuit:
		receivers.Remove( sock.SourceAddress().Ipv4Address() );
		OnDisconnect(sock.SourceAddress());
		break;

	case TypeError:
		ProcessError((ErrorPacket*)message, length);
		message = NULL;
		break;

	default:
		SendError(2, "Invalid protocol message type.");
	};
};


void Sender::HandleTimeout()
{
	// We timed out so resend the data packets or send register request
	// if no receivers are registered.
	// If the retry limit was reached for a particular packet then remove 
	// the destination address from the receivers list since the packet
	// didn't make it to its destination, i.e. connection lost.
	if (receivers.Empty())
	{
		SendRegisterRequest();
	}
	else
	{
		IPList lost;
		for (UShort id = 0; id < 256; id++)
		{
			if (stilltoack[id].Empty()) continue;
			
			if (retrylimit[id] == 0)
			{
				// If the retry limit was reached for packet number id then indicate 
				// a lost connection and remove it from the sending packet list.
				try
				{
					for (UInt i = 0; i < stilltoack[id].Count(); i++)
					{
						// Check if we already have this IP in the lost list.
						// If not then signal lost connection for this IP else
						// it was already done so move on.
						if (not lost.Contains(stilltoack[id][i]))
						{
							lost.Add(stilltoack[id][i]);
							Address address(stilltoack[id][i], broadcastaddr.Port());
							OnConnectionLost(address);
						};
					};
					stilltoack[id].Clear();
				}
				finally
				(
					Assert( busypackets[id].packet != NULL );
					DebugMsg(1, "Packet dropped: " << (int)id);
					FreeMessage( (char*)busypackets[id].packet );
					busypackets[id].packet = NULL;
					available_id.Push(id);
					queues_not_empty.Signal(false);
				);
			}
			else
			{
				// Resend the packet for which we did not get a reply yet.
				DebugMsg(1, "resend packet: " << (int)id << " size = " << busypackets[id].size);
				sock.TargetAddress(broadcastaddr);
				sock.Send((char*)busypackets[id].packet, busypackets[id].size);
				retrylimit[id]--;
			};
		};
	};
};


void Sender::SendRegisterRequest()
{
	RegisterPacket packet;
	packet.type = TypeRegister;
	sock.TargetAddress(broadcastaddr);
	sock.Send((char*)&packet, sizeof(RegisterPacket));
};


void Sender::SendRelease()
{
	ReleasePacket packet;
	packet.type = TypeRelease;
	sock.TargetAddress(broadcastaddr);
	sock.Send((char*)&packet, sizeof(ReleasePacket));
};


void Sender::ProcessAcknowledge(char* message, const UInt length)
{
	if (length < sizeof(UChar) * 2)
	{
		SendError(10, "Malformed acknowledge packet.");
		return;
	};

	if (not receivers.Contains(sock.SourceAddress().Ipv4Address()) )
	{
		receivers.Add( sock.SourceAddress().Ipv4Address() );
		queues_not_empty.Signal(false);
		OnConnect(sock.SourceAddress());
	};

	AcknowledgePacket* packet = (AcknowledgePacket*)message;
	if (packet->count < length - sizeof(UChar) * 2)
	{
		SendError(10, "Malformed acknowledge packet.");
		return;
	};
	
	UShort i;
	for (i = 0; i < packet->count; i++)
	{
		UChar id = packet->packetid[i];
		stilltoack[id].Remove( sock.SourceAddress().Ipv4Address() );
		DebugMsg(1, "Received ACK for packet: " << (int)id << " from " << sock.SourceAddress());
		stilltoack[id].Dump();
	};

	// Check if a packet has been delivered to all its destinations.
	// If it has, then we can clean up internal data structures.
	for (i = 0; i < 256; i++)
	{
		if (stilltoack[i].Empty() and busypackets[i].packet != NULL)
		{
			DebugMsg(1, "Done with packet: " << i);
			FreeMessage( (char*)busypackets[i].packet );
			busypackets[i].packet = NULL;
			available_id.Push(i);
			queues_not_empty.Signal(false);
		};
	};
};


void Sender::ProcessError(ErrorPacket* packet, UInt length)
{
	if (length > sizeof(UChar) + sizeof(UInt) + 1)
	{
		// Make sure there is a NULL at the end of the message string.
		((char*)packet)[length-1] = '\0';
	}
	else
	{
		char* str = "Got protocol error (The packet was malformed however!).";
		UInt newlength = strlen(str) + sizeof(UChar) + sizeof(UInt) + 1;
		packet = (ErrorPacket*) ReallocateMessage( (char*)packet, newlength );
		packet->errorcode = -1;
		strcpy(&packet->message[0], str);
	};
	
	eventqueue.Push(OnErrorEvent, sock.SourceAddress(), (char*)packet, length);
};


void Sender::SendError(const Int code, const char* message)
{
	UInt length = strlen(message) + 1 + sizeof(PacketType) + sizeof(UInt);
	char buf[1024];
	Assert( length <= sizeof(buf) );
	ErrorPacket* packet = (ErrorPacket*)&buf[0];
	packet->type = TypeError;
	packet->errorcode = code;
	strcpy(&packet->message[0], message);

	sock.TargetAddress( sock.SourceAddress() );
	sock.Send((char*)packet, length);
};


} // BCMP
} // dHLT
