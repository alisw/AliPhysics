////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "BCMP/Receiver.hpp"
#include "Utils.hpp"
#include "Error.hpp"
#include <asm/errno.h>

#include <iostream>
using std::cerr;
using std::endl;

using namespace dHLT::System;

// Total timeout in Run loop = TIMEOUT_SEC + TIMEOUT_MICROSEC / 1e6 seconds
#define TIMEOUT_SEC         1
#define TIMEOUT_MICROSEC    0

// TODO: replace stl::vector with internal class.
//       the stl on redhat seems to be doing very strange things with memory
//       which can not be easily validated.

namespace dHLT
{
namespace BCMP
{


bool Receiver::PortList::Add(const UInt ip, const UShort port)
{
	// Make sure the ip:port is not in the list already.
	std::vector<Record>::iterator rec;
	for (rec = list.begin(); rec != list.end(); rec++)
	{
		if (rec->ip == ip and rec->port == port)
			return false;
	};

	Record newrec;
	newrec.ip = ip;
	newrec.port = port;
	list.push_back(newrec);
	return true;
};


bool Receiver::PortList::Remove(const UInt ip, const UShort port)
{
	std::vector<Record>::iterator rec;
	for (rec = list.begin(); rec != list.end(); rec++)
	{
		if (rec->ip == ip and rec->port == port)
		{
			list.erase(rec);
			return true;
		};
	};
	return false;
};


bool Receiver::PortList::Contains(const UInt ip, const UShort port) const
{
	std::vector<Record>::const_iterator rec;
	for (rec = list.begin(); rec != list.end(); rec++)
	{
		if (rec->ip == ip and rec->port == port)
			return true;
	};
	return false;
};


UInt Receiver::PortList::Count() const
{
	return list.size();
};

void Receiver::PortList::Fetch(const UInt index, UInt& ip, UShort& port)
{
	Assert( index < Count() );
	ip = list[index].ip;
	port = list[index].port;
};


void Receiver::PortList::Dump()
{
	DebugMsg(2, "ip\tport");
	DebugMsg(2, "--------------------");
	std::vector<Record>::iterator rec;
	for (rec = list.begin(); rec != list.end(); rec++)
	{
		DebugMsg(2, Address::AsString(rec->ip) << "\t" << rec->port);
	};
};

///////////////////////////////////////////////////////////////////////////////

void Receiver::AckList::Add(const UInt ip, const UShort port, const UChar packetid)
{
	// Try find the ip:port in the list.
	std::vector<Record>::iterator rec;
	for (rec = list.begin(); rec != list.end(); rec++)
	{
		if (rec->IP() == ip and rec->Port() == port)
		{
			rec->AddID(packetid);
			return;
		}
	};

	Record newrec(ip, port);
	newrec.AddID(packetid);
	list.push_back(newrec);
};


bool Receiver::AckList::Empty() const
{
	return list.size() == 0;
};


void Receiver::AckList::Clear()
{
	list.erase(list.begin(), list.end());
};


UInt Receiver::AckList::Count() const
{
	return list.size();
};


const Receiver::AckList::Record& 
	Receiver::AckList::operator [] (const UInt index) const
{
	Assert( index < Count() );
	return list[index];
};

///////////////////////////////////////////////////////////////////////////////

Int Receiver::ReceiveThread::Execute()
{
	Int returncode = 0;
	try
	{
		Assert( receiver != NULL );
		receiver->Run();
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
		cerr << "<ReceiveThread> Unknown exception!" << endl;
		returncode = -1;
	};
	return returncode;
};

///////////////////////////////////////////////////////////////////////////////

Receiver::Receiver(const UShort port)
{
	DebugMsg(9, "Receiver::Receiver(" << port << ")");
	eventhandler = NULL;
	Init(port);
};


Receiver::Receiver(EventHandler* handler, const UShort port)
{
	eventhandler = handler;
	Init(port);
};


void Receiver::Init(const UShort port)
{
	sock.AllowAddressReuse();
	Address address = sock.BroadcastAddress( Address::GetHostAddress().Ipv4Address() );
	address.Port(port);
	sock.Bind(address);
	
	thread.receiver = this;
	thread.Start();
};


Receiver::~Receiver()
{
	Terminate();
	thread.WaitForThread();
};


bool Receiver::TryHandleEvents()
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
		
		case OnMessageEvent:
		{
			DataPacket* packet = (DataPacket*)message;
			eventhandler->OnMessage(&packet->data[0], length - 2 * sizeof(UChar), address);
			FreeMessage(message);
			break;
		}
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		default:                    throw UnknownEvent();
		};
	};
	return gotmessage;
};


bool Receiver::HandleEvents(const UInt timeout)
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
		
		case OnMessageEvent:
		{
			DataPacket* packet = (DataPacket*)message;
			eventhandler->OnMessage(&packet->data[0], length - 2 * sizeof(UChar), address);
			FreeMessage(message);
			break;
		}
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		default:                    throw UnknownEvent();
		};
	};
	return gotmessage;
};


void Receiver::HandleEvents()
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
		
		case OnMessageEvent:
		{
			DataPacket* packet = (DataPacket*)message;
			eventhandler->OnMessage(&packet->data[0], length - 2 * sizeof(UChar), address);
			FreeMessage(message);
			break;
		}
		
		case OnErrorEvent:
		{
			ErrorPacket* packet = (ErrorPacket*)message;
			eventhandler->OnError(&packet->message[0], packet->errorcode, address);
			FreeMessage(message);
			break;
		}
		
		default:                    throw UnknownEvent();
		};
	}
	while (not terminate);
};


void Receiver::Run()
{
	// Enter processing loop.
	terminate = false;
	do
	{
		if (not sock.WaitUntilReadable(TIMEOUT_SEC, TIMEOUT_MICROSEC)) continue;
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

			if (not acklist.Empty())
			{
				SendAcknowledges();
				acklist.Clear();
			};
		}
		catch (const SocketError& e)
		{
			// If we get the WSAECONNRESET error that means the client aborted.
			// The previous send was unsuccessfull.
#ifdef WIN32
			if (e.ErrorCode() == WSAECONNRESET)
			{
				senders.Remove(sock.TargetAddress().Ipv4Address(), sock.TargetAddress().Port());
				OnConnectionLost(sock.TargetAddress());
			}
			else
#endif // WIN32
				throw;
		}
		catch (const System::Error& e)
		{
			if (e.ErrorCode() != EINTR) throw;
		};
	}
	while (not terminate);

	// Send quit message to all senders so that we no longer receive messages.
	for (UInt i = 0; i < senders.Count(); i++)
	{
		UInt ip; UShort port;
		senders.Fetch(i, ip, port);
		Address toaddr(ip, port);
		SendQuit(toaddr);
	};
};


void Receiver::OnConnect(const Address& address)
{
	DebugMsg(1, "Connection opened by: " << address);
	DebugCode( senders.Dump(); );
	eventqueue.Push(OnConnectEvent, address);
};


void Receiver::OnReceive(const Address& from, char* message, const UInt length)
{
	DebugMsg(1, "Receiver::OnReceive(" << from << ", " << (void*)message << ", " << length << ")");
	eventqueue.Push(OnMessageEvent, from, message, length);
};


void Receiver::OnDisconnect(const Address& address)
{
	DebugMsg(1, "Receiver::OnDisconnect(" << address << ")");
	DebugCode( senders.Dump(); );
	eventqueue.Push(OnDisconnectEvent, address);
};


void Receiver::OnConnectionLost(const Address& address)
{
	DebugMsg(1, "Receiver::OnConnectionLost(" << address << ")");
	DebugCode( senders.Dump(); );
	eventqueue.Push(OnConnectionLostEvent, address);
};



void Receiver::ProcessMessage(char*& message, const UInt length)
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
		senders.Add( sock.SourceAddress().Ipv4Address(), sock.SourceAddress().Port() );
		OnConnect(sock.SourceAddress());
		SendConfirm();
		break;

	case TypeConfirm:
		SendError(4, "Unexpected confirm message received.");
		break;

	case TypeData:
		ProcessDataPacket(message, length);
		message = NULL;  // Do not delete the message yet.
		break;

	case TypeAcknowledge:
		SendError(6, "Unexpected acknowledge message received.");
		break;

	case TypeRelease:
		senders.Remove( sock.SourceAddress().Ipv4Address(), sock.SourceAddress().Port() );
		OnDisconnect(sock.SourceAddress());
		break;

	case TypeQuit:
		SendError(8, "Unexpected quit message received.");
		break;

	case TypeError:
		ProcessError((ErrorPacket*)message, length);
		message = NULL;  // Do not delete the error message yet.
		break;

	default:
		SendError(2, "Invalid protocol message type.");
	};
};


void Receiver::SendConfirm()
{
	ConfirmPacket packet;
	packet.type = TypeConfirm;
	sock.TargetAddress( sock.SourceAddress() );
	sock.Send((char*)&packet, sizeof(ConfirmPacket));
};


void Receiver::ProcessDataPacket(char* message, const UInt length)
{
	if (length < sizeof(UChar) * 2)
	{
		SendError(9, "Malformed data packet.");
		return;
	};

	if (not senders.Contains(sock.SourceAddress().Ipv4Address(), sock.SourceAddress().Port()) )
	{
		senders.Add( sock.SourceAddress().Ipv4Address(), sock.SourceAddress().Port() );
		OnConnect(sock.SourceAddress());
	};

	DataPacket* packet = (DataPacket*)message;
	acklist.Add( sock.SourceAddress().Ipv4Address(), sock.SourceAddress().Port(), packet->packetid );
	OnReceive( sock.SourceAddress(), message, length );
};


void Receiver::SendAcknowledges()
{
	union
	{
		AcknowledgePacket packet;
		char buffer[264];
	} pd;

	for (UInt i = 0; i < acklist.Count(); i++)
	{
		const AckList::Record& rec = acklist[i];
		
		pd.packet.type = TypeAcknowledge;
		pd.packet.count = rec.IDCount();
		for (UInt j = 0; j < pd.packet.count; j++)
			pd.packet.packetid[j] = rec.ID(j);

		UInt length = sizeof(UChar) * 2 + pd.packet.count;
		SendAck( (char*)&pd.packet, length );
	};
};


void Receiver::SendAck(const char* message, const UInt length)
{
	sock.TargetAddress( sock.SourceAddress() );
	sock.Send(message, length);
};


void Receiver::SendQuit(const Address& destination)
{
	QuitPacket packet;
	packet.type = TypeQuit;
	sock.TargetAddress(destination);
	sock.Send((char*)&packet, sizeof(QuitPacket));
};


void Receiver::ProcessError(ErrorPacket* packet, UInt length)
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


void Receiver::SendError(const Int code, const char* message)
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
