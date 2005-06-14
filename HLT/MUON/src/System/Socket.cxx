////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Socket.hpp"
#include "SystemTypes.hpp"
#include "Utils.hpp"
#include <queue>


#ifdef WIN32
#	include "System/Win32/Socket.cxx"
#endif // WIN32

#ifdef LINUX
#	include "System/Linux/Socket.cxx"
#endif // LINUX


namespace dHLT
{
namespace System
{


Address::Address()
{
	data.sin_family = AF_INET;
	data.sin_port = 0;
	data.sin_addr.s_addr = INADDR_NONE;
	memset(&data.sin_zero, 0, sizeof(data.sin_zero));
};


Address::Address(const SocketAddress& address)
{
	memcpy(&data, &address, sizeof(data));
	data.sin_family = AF_INET;
};


Address::Address(const InetSocketAddress& address)
{
	memcpy(&data, &address, sizeof(data));
	data.sin_family = AF_INET;
};


Address::Address(const UInt ipv4address, const UShort port)
{
	data.sin_family = AF_INET;
	data.sin_port = htons(port);
	data.sin_addr.s_addr = htonl(ipv4address);
	memset(&data.sin_zero, 0, sizeof(data.sin_zero));
};


Address::Address(const UShort netbyte_order_port, const UInt netbyte_order_address)
{
	data.sin_family = AF_INET;
	data.sin_port = netbyte_order_port;
	data.sin_addr.s_addr = netbyte_order_address;
	memset(&data.sin_zero, 0, sizeof(data.sin_zero));
};

Address::Address(const char* string)
{
	data.sin_family = AF_INET;

	char buf[32];
	char* str = (char*) &buf;
	const char* portstr = NULL;

	// Copy the IPV4 address part to str and mark the port part
	// of the address by setting portstr to the first character
	// following the ':' character.
	for (UInt i = 0; i < sizeof(buf); i++)
	{
		buf[i] = string[i];
		if (string[i] == '\0') break;
		if (string[i] == ':')
		{
			portstr = &string[i+1];
			buf[i] = '\0';
			break;
		};
	};

	data.sin_port = htons(atoi(portstr));
	data.sin_addr.s_addr = inet_addr(str);
	memset(&data.sin_zero, 0, sizeof(data.sin_zero));
};


UInt Address::Ipv4Address() const
{
	return ntohl(data.sin_addr.s_addr);
};


void Address::Ipv4Address(const UInt value)
{
	data.sin_addr.s_addr = htonl(value);
};


UShort Address::Port() const
{
	return ntohs(data.sin_port);
};


void Address::Port(const UShort port)
{
	data.sin_port = htons(port);
};


const char* Address::AsString(const UInt ipv4address)
{
	struct in_addr addr;
	addr.s_addr = htonl(ipv4address);
	return inet_ntoa(addr);
};


const char* Address::AsStringNB(const UInt netbyte_order_address)
{
	struct in_addr addr;
	addr.s_addr = netbyte_order_address;
	return inet_ntoa(addr);
};


const char* Address::AsString() const
{
	static char buf[32];
	char* str = (char*) &buf;
	sprintf(str, "%s:%d", inet_ntoa(data.sin_addr), ntohs(data.sin_port));
	return str;
};


std::ostream& operator << (std::ostream& os, const Address& addr)
{
	os << addr.AsString();
	return os;
};

////////////////////////////////////////////////////////////////////////////////

Socket::Socket()
{
	handle = socket(PF_INET, SOCK_DGRAM, 0);
	if (handle == INVALID_SOCKET) throw SocketError();
};


Socket::Socket(const SocketHandle s)
{
	handle = s;
};


Socket::~Socket()
{
	int result = closesocket(handle);
	if (result == SOCKET_ERROR)	throw SocketError();
};


Address Socket::LocalAddress() const
{
	Address addr;
	socklen_t addrsize = addr.Size();
	int result = getsockname(handle, addr, &addrsize);
	if (result == SOCKET_ERROR) throw SocketError();
	return addr;
};


void Socket::Bind(const Address& address)
{
	int result = bind(handle, address, address.Size());
	if (result == SOCKET_ERROR) throw SocketError();
};


UInt Socket::Send(const char* message, const UInt length)
{
	int result = sendto(handle, message, length, 0, targetaddr, targetaddr.Size());
	if (result == SOCKET_ERROR) throw SocketError();
	return (UInt)result;
};


UInt Socket::Receive(char* message, const UInt length)
{
	socklen_t sourceaddrsize = sourceaddr.Size();
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
	int result = select(handle+1, &set, NULL, NULL, &timeout);
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
	int result = select(handle+1, NULL, &set, NULL, &timeout);
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
	int result = select(handle+1, &set, NULL, NULL, timeout);
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
	int result = select(handle+1, NULL, &set, NULL, timeout);
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
		UInt i;
		for (i = 0; i < sl.list.size() - 1; i++)
		{
			os << (UInt) sl.list[i]->Handle() << " ";
		};
		os << (UInt) sl.list[i]->Handle();
	};
	os << "]";
	return os;
};


Socket* WaitFor(
		SocketStatus& status, const SocketList& sockets, const int eventflags,
		const UInt timeout_secs, const UInt timeout_usecs
	)
{
	static std::queue<ReadyRecord> readyqueue;

	if (readyqueue.empty())
	{
		UInt i;
		int maxnum = 0;

		// Initialise the socket handle sets.
		fd_set readset, writeset, errorset;
		if (eventflags & SocketReadable)
		{
			FD_ZERO(&readset);
			for (i = 0; i < sockets.Size(); i++)
			{
				if (sockets[i]->Handle() > maxnum) maxnum = sockets[i]->Handle();
				FD_SET(sockets[i]->Handle(), &readset);
			};
		};
		if (eventflags & SocketWriteable)
		{
			FD_ZERO(&writeset);
			for (i = 0; i < sockets.Size(); i++)
			{
				if (sockets[i]->Handle() > maxnum) maxnum = sockets[i]->Handle();
				FD_SET(sockets[i]->Handle(), &writeset);
			};
		};
		if (eventflags & SocketFailure)
		{
			FD_ZERO(&errorset);
			for (i = 0; i < sockets.Size(); i++)
			{
				if (sockets[i]->Handle() > maxnum) maxnum = sockets[i]->Handle();
				FD_SET(sockets[i]->Handle(), &errorset);
			};
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
		int result = select(maxnum+1, rs, ws, es, timeout);
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
		int maxnum = 0;

		// Initialise the socket handle sets.
		fd_set readset;
		FD_ZERO(&readset);
		for (i = 0; i < sockets.Size(); i++)
		{
			if (sockets[i]->Handle() > maxnum) maxnum = sockets[i]->Handle();
			FD_SET(sockets[i]->Handle(), &readset);
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

		// Wait for something to happen on the sockets.
		int result = select(maxnum+1, &readset, NULL, NULL, timeout);
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
		int maxnum = 0;

		// Initialise the socket handle sets.
		fd_set writeset;
		FD_ZERO(&writeset);
		for (i = 0; i < sockets.Size(); i++)
		{
			if (sockets[i]->Handle() > maxnum) maxnum = sockets[i]->Handle();
			FD_SET(sockets[i]->Handle(), &writeset);
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

		// Wait for something to happen on the sockets.
		int result = select(maxnum+1, NULL, &writeset, NULL, timeout);
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


} // System
} // dHLT

