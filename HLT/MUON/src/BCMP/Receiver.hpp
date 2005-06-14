////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BCMP_RECEIVER_HPP
#define dHLT_BCMP_RECEIVER_HPP

#include "System/Socket.hpp"
#include "System/Thread.hpp"
#include "BCMP/Packets.hpp"
#include "BCMP/EventQueue.hpp"
#include "BCMP/EventHandler.hpp"
#include <vector>

namespace dHLT
{
namespace BCMP
{


class Receiver
{
public:

	Receiver(const UShort port = 4900);
	Receiver(EventHandler* handler, const UShort port = 4900);
	~Receiver();

	void Terminate() { terminate = true; };
	
	System::Address LocalAddress() const
	{
		return sock.LocalAddress();
	};
	
	bool TryHandleEvents();
	
	// timeout in milliseconds. Returns false if we timed out.
	bool HandleEvents(const UInt timeout);

	void HandleEvents();
	
	EventHandler* eventhandler;

private:

	void Init(const UShort port);
	void Run();

	void ProcessMessage(char*& message, const UInt length);
	void SendConfirm();
	void ProcessDataPacket(char* message, const UInt length);
	void SendAcknowledges();
	void SendAck(const char* message, const UInt length);
	void SendQuit(const System::Address& destination);
	void ProcessError(ErrorPacket* message, UInt length);
	void SendError(const Int code, const char* message);
	
	void OnConnect(const System::Address& address);
	void OnReceive(const System::Address& from, char* packet, const UInt length);
	void OnDisconnect(const System::Address& address);
	void OnConnectionLost(const System::Address& address);


	class PortList
	{
	public:

		// Returns true if the ip:port was added and not already found in the list.
		bool Add(const UInt ip, const UShort port);

		bool Remove(const UInt ip, const UShort port);
		bool Contains(const UInt ip, const UShort port) const;
		UInt Count() const;
		void Fetch(const UInt index, UInt& ip, UShort& port);
		void Dump();

	private:

		struct Record
		{
			UInt ip;
			UShort port;
		};

		std::vector<Record> list;
	};


	class AckList
	{
	public:

		class Record
		{
		public:
			
			Record(const UInt ip = 0, const UShort port = 0)
			{
				this->ip = ip;
				this->port = port;
				count = 0;
			};


			void AddID(const UChar id)
			{
				Assert( count < 256 );
				packetid[count++] = id;
			};


			UInt IP() const
			{
				return ip;
			};


			UShort Port() const
			{
				return port;
			};

			UChar IDCount() const
			{
				return count;
			};

			UChar ID(const UChar index) const
			{
				Assert( index < IDCount() );
				return packetid[index];
			};

		private:

			UInt ip;
			UShort port;
			UChar count;
			UChar packetid[256];
		};


		void Add(const UInt ip, const UShort port, const UChar packetid);
		bool Empty() const;
		void Clear();
		UInt Count() const;
		const Record& operator [] (const UInt index) const;

	private:

		std::vector<Record> list;
	};


	class ReceiveThread : public System::Thread
	{
	public:
		virtual Int Execute();
		Receiver* receiver;
	};
	

	ReceiveThread thread;
	EventQueue eventqueue;

	System::Socket sock;
	bool terminate;
	PortList senders;
	AckList acklist;
};


} // BCMP
} // dHLT

#endif // dHLT_BCMP_RECEIVER_HPP
