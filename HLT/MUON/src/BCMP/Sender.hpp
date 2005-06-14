////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BCMP_SENDER_HPP
#define dHLT_BCMP_SENDER_HPP

#include "System/Socket.hpp"
#include "System/Thread.hpp"
#include "System/Mutex.hpp"
#include "System/MutexCondition.hpp"
#include "BCMP/Packets.hpp"
#include "BCMP/EventHandler.hpp"
#include "BCMP/EventQueue.hpp"

#include <queue>

namespace dHLT
{
namespace BCMP
{


class Sender
{
public:

	Sender(const UShort port = 4900);
	Sender(EventHandler* handler, const UShort port = 4900);
	~Sender();

	void Send(const char* message, const UInt length);

	void Terminate() { terminate = true; };
	
	System::Address LocalAddress() const { return sock.LocalAddress(); };
	
	bool TryHandleEvents();
	
	// timeout in milliseconds. Returns false if we timed out.
	bool HandleEvents(const UInt timeout);

	void HandleEvents();

	EventHandler* eventhandler;

private:

	void Init(const UShort port);

	void SendLoop();
	void ReceiveLoop();

	void ReadFromSocket();
	void ProcessMessage(char*& message, const UInt length);
	void HandleTimeout();
	void SendRegisterRequest();
	void SendRelease();
	void ProcessAcknowledge(char* message, const UInt length);
	void ProcessError(ErrorPacket* message, UInt length);
	void SendError(const Int code, const char* message);

	void OnConnect(const System::Address& address);
	void OnDisconnect(const System::Address& address);
	void OnConnectionLost(const System::Address& address);
	
	class PacketQueue
	{
	public:

		~PacketQueue();
		void Push(DataPacket* packet, const UInt size);
		
		// returns false if timed out.
		bool Pop(DataPacket*& packet, UInt& size);

	private:

		struct Record
		{
			DataPacket* packet;
			UInt packetsize;
		};

		std::queue<Record> pq;
		System::Mutex mutex;
		System::MutexCondition queue_not_empty;
	};


	class IPList
	{
	public:

		// Returns true if the ip was added and not already found in the list.
		bool Add(const UInt ip);
		bool Remove(const UInt ip);
		bool Empty() const;
		void Clear();
		void Copy(const IPList& rhs);
		bool Contains(const UInt ip) const;
		UInt Count() const;
		UInt operator [] (const UInt index) const;
		void Dump();

	private:

		std::vector<UInt> list;
	};


	class PacketIDStack
	{
	public:

		PacketIDStack()
		{
			count = 0;
		};
		
		bool Empty()
		{
			return count == 0;
		};
		
		bool Full()
		{
			return count == 256;
		};

		void Push(const UChar value)
		{
			Assert( count < 256 );
			data[count++] = value;
		};

		UChar Pop()
		{
			Assert( count > 0 );
			return data[--count];
		};

	private:

		UShort count;
		UChar data[256];
	};
	
	
	class SendThread : public System::Thread
	{
	public:
		virtual Int Execute();
		Sender* sender;
	};

	class ReceiveThread : public System::Thread
	{
	public:
		virtual Int Execute();
		Sender* sender;
	};
	

	SendThread sendthread;
	ReceiveThread receivethread;

	System::Socket sock;
	System::Address broadcastaddr;
	bool terminate;
	bool canterminate;
	IPList receivers;
	PacketQueue packetqueue;
	EventQueue eventqueue;

	PacketIDStack available_id;
	IPList stilltoack[256];
	UInt retrylimit[256];

	struct PacketRecord
	{
		DataPacket* packet;
		UInt size;
	};

	PacketRecord busypackets[256];
	
	System::Mutex mutex;
	System::MutexCondition queues_not_empty;

	static UInt default_retry_limit;
};


} // BCMP
} // dHLT

#endif // dHLT_BCMP_SENDER_HPP
