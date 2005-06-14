////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Socket.hpp"
#include <windows.h>
#include <stdio.h>
#include <iostream>
#include <process.h>

#define s_addr S_un.S_addr
typedef int socklen_t;

namespace 
{

	class Initialisation
	{
	public:

		Initialisation()
		{
			using std::cerr;
			using std::endl;
			
			// Start windows socket layer:
			int major = 2;
			int minor = 0;

			WSADATA wsadata;
			WORD version = MAKEWORD(major, minor);
			int errorcode = WSAStartup(version, &wsadata);
			if ( errorcode != 0 )
			{
				cerr << "Error: " << System::Error::AsString(errorcode) << endl;
				exit(1000000);
			};
 
			// Confirm winsock version.
			if ( LOBYTE( wsadata.wVersion ) != major || HIBYTE( wsadata.wVersion ) != minor )
			{
				WSACleanup();
				cerr << "Error: Could not find required winsock version." << endl;
				exit(1000001);
			};
		};
		

		~Initialisation()
		{
			WSACleanup();
		};

	} initialisation;

};


namespace dHLT
{
namespace System
{


SocketError::SocketError() throw()
{
	errorcode = WSAGetLastError();
};

const char* SocketError::Message() const throw()
{
	return AsString(errorcode);
};


Address Address::GetHostAddress(const UShort port)
{
	char buffer[1024];
	char* hostname = (char*) &buffer[0];

	int result = gethostname(hostname, sizeof(buffer));
	if (result == SOCKET_ERROR) throw SocketError();

	struct hostent* ent = gethostbyname(hostname);
	if (ent == NULL) throw SocketError();

	UInt addr;
	if ( ent->h_length >= (short)sizeof(addr) )
		addr = *( (UInt*)ent->h_addr_list[0] );
	else
		addr = INADDR_ANY;

	return Address(htons(port), addr);
};


////////////////////////////////////////////////////////////////////////////////

bool Socket::CanBroadcast() const
{
	BOOL value;
	int valuesize = sizeof(value);
	int result = getsockopt(handle, SOL_SOCKET, SO_BROADCAST, (char*) &value, &valuesize);
	if (result == SOCKET_ERROR) throw SocketError();
	return value == TRUE;
};


void Socket::AllowBroadcast(bool value)
{
	BOOL ivalue = value ? TRUE : FALSE;
	int result = setsockopt(handle, SOL_SOCKET, SO_BROADCAST, (const char*) &ivalue, sizeof(ivalue));
	if (result == SOCKET_ERROR) throw SocketError();
};


bool CanReuseAddress() const
{
	BOOL value;
	int valuesize = sizeof(value);
	int result = getsockopt(handle, SOL_SOCKET, SO_REUSEADDR, (char*) &value, &valuesize);
	if (result == SOCKET_ERROR) throw SocketError();
	return value == TRUE;
};


void AllowAddressReuse(bool value)
{
	BOOL ivalue = value ? TRUE : FALSE;
	int result = setsockopt(handle, SOL_SOCKET, SO_REUSEADDR, (const char*) &ivalue, sizeof(ivalue));
	if (result == SOCKET_ERROR) throw SocketError();
};


Address Socket::BroadcastAddress(const UInt /*localaddress*/) const
{
	// Ignores the localaddress field for now.

	Address addr;
	unsigned long addrsize = addr.Size();

	int result = WSAIoctl(
			handle,
			SIO_GET_BROADCAST_ADDRESS,
			NULL,
			0,
			addr,
			addrsize,
			&addrsize,
			NULL,
			NULL
		);

	if (result == SOCKET_ERROR) throw SocketError();
	return addr;
};


UInt Socket::MaxMessageSize() const
{
	UInt value;
	int valuesize = sizeof(value);
	int result = getsockopt(handle, SOL_SOCKET, SO_MAX_MSG_SIZE, (char*) &value, &valuesize);
	if (result == SOCKET_ERROR) throw SocketError();
	return value;
};


UInt Socket::BytesAvailable() const
{
	unsigned long bytesavailable;
	unsigned long size = sizeof(bytesavailable);

	int result = WSAIoctl(
			handle,
			FIONREAD,
			NULL,
			0,
			bytesavailable,
			size,
			&size,
			NULL,
			NULL
		);

	if (result == SOCKET_ERROR) throw SocketError();
	return bytesavailable;
};

/*
UInt Socket::Receive(char* message, const UInt length)
{
	int sourceaddrsize = sourceaddr.Size();
	int result = recvfrom(handle, message, length, 0, sourceaddr, &sourceaddrsize);
	if (result == SOCKET_ERROR) throw SocketError();
	return (UInt)result;
};


bool Socket::CanRead() const
{
	// Initialise the socket handle set.
	fd_set set;
	FD_ZERO(&set);
	FD_SET(handle, &set);

	// Initialise the timeout so that we return immediately.
	timeval timeout;
	timeout.tv_sec = 0;
	timeout.tv_usec = 0;

	// Check if the handle is ready.
	int result = select(0, &set, NULL, NULL, &timeout);
	if (result == SOCKET_ERROR) throw SocketError();
	return FD_ISSET(handle, &set) != 0;
};


bool Socket::CanWrite() const
{
	// Initialise the socket handle set.
	fd_set set;
	FD_ZERO(&set);
	FD_SET(handle, &set);

	// Initialise the timeout so that we return immediately.
	timeval timeout;
	timeout.tv_sec = 0;
	timeout.tv_usec = 0;

	// Check if the handle is ready.
	int result = select(0, NULL, &set, NULL, &timeout);
	if (result == SOCKET_ERROR) throw SocketError();
	return FD_ISSET(handle, &set) != 0;
};


bool Socket::WaitUntilReadable(const UInt timeout_secs, const UInt timeout_usecs) const
{
	// Initialise the socket handle sets.
	fd_set set;
	FD_ZERO(&set);
	FD_SET(handle, &set);

	// Initialise the timeout structure. Set to NULL if we want to block for ever.
	timeval tv;
	timeval* timeout;
	if (timeout_usecs == 0xFFFFFFFF and timeout_secs == 0xFFFFFFFF)
		timeout = NULL;
	else
	{
		tv.tv_sec = timeout_secs;
		tv.tv_usec = timeout_usecs;
		timeout = &tv;
	};

	// Wait for something to happen on the socket.
	int result = select(0, &set, NULL, NULL, timeout);
	if (result == SOCKET_ERROR) throw SocketError();
	if (result == 0) return false;

	// Check which handle is ready to read from.
	return FD_ISSET(handle, &set) != 0;
};


bool Socket::WaitUntilWriteable(const UInt timeout_secs, const UInt timeout_usecs) const
{
	// Initialise the socket handle sets.
	fd_set set;
	FD_ZERO(&set);
	FD_SET(handle, &set);

	// Initialise the timeout structure. Set to NULL if we want to block for ever.
	timeval tv;
	timeval* timeout;
	if (timeout_usecs == 0xFFFFFFFF and timeout_secs == 0xFFFFFFFF)
		timeout = NULL;
	else
	{
		tv.tv_sec = timeout_secs;
		tv.tv_usec = timeout_usecs;
		timeout = &tv;
	};

	// Wait for something to happen on the socket.
	int result = select(0, NULL, &set, NULL, timeout);
	if (result == SOCKET_ERROR) throw SocketError();
	if (result == 0) return false;

	// Check which handle is ready to write to.
	return FD_ISSET(handle, &set) != 0;
};


struct ReadyRecord
{
	SocketStatus status;
	const Socket* sockptr;

	ReadyRecord(const SocketStatus stat = SocketStatusUnknown, const Socket* ptr = NULL)
	{
		status = stat;
		sockptr = ptr;
	};
};


void SocketList::Add(Socket& s)
{
	list.push_back(&s);
};


void SocketList::Remove(Socket& s)
{
	std::vector<Socket*>::iterator pos;
	pos = std::find(list.begin(), list.end(), &s);
	if (pos != list.end())
		list.erase(pos);
};


bool SocketList::Contains(Socket& s)
{
	std::vector<Socket*>::iterator pos;
	pos = std::find(list.begin(), list.end(), &s);
	return pos != list.end();
};


const Socket* SocketList::operator [] (UInt index) const
{
	Assert( index < list.size() );
	return list[index];
};


Socket* SocketList::operator [] (UInt index)
{
	Assert( index < list.size() );
	return list[index];
};


std::ostream& operator << (std::ostream& os, const SocketList& sl)
{
	os << "[";
	if (sl.list.size() > 0)
	{
		Int i;
		for (i = 0; i < sl.list.size() - 1; i++)
		{
			os << (UInt) sl.list[i]->Handle() << " ";
		};
		os << (UInt) sl.list[i]->Handle();
	};
	os << "]";
	return os;
};


#include <queue>

Socket* WaitFor(
		SocketStatus& status, const SocketList& sockets, const int eventflags,
		const UInt timeout_secs, const UInt timeout_usecs
	)
{
	static std::queue<ReadyRecord> readyqueue;

	if (readyqueue.empty())
	{
		UInt i;

		// Initialise the socket handle sets.
		fd_set readset, writeset, errorset;
		if (eventflags & SocketReadable)
		{
			FD_ZERO(&readset);
			for (i = 0; i < sockets.Size(); i++)
				FD_SET(sockets[i]->Handle(), &readset);
		};
		if (eventflags & SocketWriteable)
		{
			FD_ZERO(&writeset);
			for (i = 0; i < sockets.Size(); i++)
				FD_SET(sockets[i]->Handle(), &writeset);
		};
		if (eventflags & SocketFailure)
		{
			FD_ZERO(&errorset);
			for (i = 0; i < sockets.Size(); i++)
				FD_SET(sockets[i]->Handle(), &errorset);
		};

		// Initialise the timeout structure. Set to NULL if we want to block for ever.
		timeval tv;
		tv.tv_sec = timeout_secs;
		tv.tv_usec = timeout_usecs;
		timeval* timeout;
		if (timeout_usecs == 0xFFFFFFFF and timeout_secs == 0xFFFFFFFF)
			timeout = NULL;
		else
			timeout = &tv;

		fd_set* rs, *ws, *es;
		if (eventflags & SocketReadable)
			rs = &readset;
		else
			rs = NULL;

		if (eventflags & SocketWriteable)
			ws = &writeset;
		else
			ws = NULL;

		if (eventflags & SocketFailure)
			es = &errorset;
		else
			es = NULL;

		// Wait for something to happen on the sockets.
		int result = select(0, rs, ws, es, timeout);
		if (result == SOCKET_ERROR) throw SocketError();
		if (result == 0) return NULL;

		// Check which handles are ready and then add them to the ready list.
		if (eventflags & SocketReadable)
		{
			for (i = 0; i < sockets.Size(); i++)
			{
				if ( FD_ISSET(sockets[i]->Handle(), &readset) != 0 )
					readyqueue.push( ReadyRecord(CanReadSocket, sockets[i]) );
			};
		};
		if (eventflags & SocketWriteable)
		{
			for (i = 0; i < sockets.Size(); i++)
			{
				if ( FD_ISSET(sockets[i]->Handle(), &writeset) != 0 )
					readyqueue.push( ReadyRecord(CanWriteSocket, sockets[i]) );
			};
		};
		if (eventflags & SocketFailure)
		{
			for (i = 0; i < sockets.Size(); i++)
			{
				if ( FD_ISSET(sockets[i]->Handle(), &errorset) != 0 )
					readyqueue.push( ReadyRecord(SocketHasFailed, sockets[i]) );
			};
		};
	};

	ReadyRecord rec = readyqueue.front();
	readyqueue.pop();
	status = rec.status;
	return const_cast<Socket*>(rec.sockptr);
};


Socket* WaitForReadable(const SocketList& sockets, const UInt timeout_secs, const UInt timeout_usecs)
{
	static std::queue<const Socket*> readyqueue;

	if (readyqueue.empty())
	{
		UInt i;

		// Initialise the socket handle sets.
		fd_set readset;
		FD_ZERO(&readset);
		for (i = 0; i < sockets.Size(); i++)
			FD_SET(sockets[i]->Handle(), &readset);

		// Initialise the timeout structure. Set to NULL if we want to block for ever.
		timeval tv;
		tv.tv_sec = timeout_secs;
		tv.tv_usec = timeout_usecs;
		timeval* timeout;
		if (timeout_usecs == 0xFFFFFFFF and timeout_secs == 0xFFFFFFFF)
			timeout = NULL;
		else
			timeout = &tv;

		// Wait for something to happen on the sockets.
		int result = select(0, &readset, NULL, NULL, timeout);
		if (result == SOCKET_ERROR) throw SocketError();
		if (result == 0) return NULL;

		// Check which handles are ready and then add them to the ready list.
		for (i = 0; i < sockets.Size(); i++)
		{
			if ( FD_ISSET(sockets[i]->Handle(), &readset) != 0 )
				readyqueue.push(sockets[i]);
		};
	};

	Socket* sock = const_cast<Socket*>( readyqueue.front() );
	readyqueue.pop();
	return sock;
};


Socket* WaitForWriteable(const SocketList& sockets, const UInt timeout_secs, const UInt timeout_usecs)
{
	static std::queue<const Socket*> readyqueue;

	if (readyqueue.empty())
	{
		UInt i;

		// Initialise the socket handle sets.
		fd_set writeset;
		FD_ZERO(&writeset);
		for (i = 0; i < sockets.Size(); i++)
			FD_SET(sockets[i]->Handle(), &writeset);

		// Initialise the timeout structure. Set to NULL if we want to block for ever.
		timeval tv;
		tv.tv_sec = timeout_secs;
		tv.tv_usec = timeout_usecs;
		timeval* timeout;
		if (timeout_usecs == 0xFFFFFFFF and timeout_secs == 0xFFFFFFFF)
			timeout = NULL;
		else
			timeout = &tv;

		// Wait for something to happen on the sockets.
		int result = select(0, NULL, &writeset, NULL, timeout);
		if (result == SOCKET_ERROR) throw SocketError();
		if (result == 0) return NULL;

		// Check which handles are ready and then add them to the ready list.
		for (i = 0; i < sockets.Size(); i++)
		{
			if ( FD_ISSET(sockets[i]->Handle(), &writeset) != 0 )
				readyqueue.push(sockets[i]);
		};
	};

	Socket* sock = const_cast<Socket*>( readyqueue.front() );
	readyqueue.pop();
	return sock;
};
*/

} // System
} // dHLT

