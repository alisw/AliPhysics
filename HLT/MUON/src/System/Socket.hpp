////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_SOCKET_HPP
#define dHLT_SYSTEM_SOCKET_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "System/SystemError.hpp"
#include "System/SystemTypes.hpp"

#include <ostream>
#include <vector>
#include <algorithm>

namespace dHLT
{
namespace System
{


class SocketError : public System::Error
{
public:
	SocketError() throw();
	virtual const char* Message() const throw();
};


class Address
{
public:

	Address();
	Address(const SocketAddress& address);
	Address(const InetSocketAddress& address);
	Address(const UInt ipv4address, const UShort port);
	Address(const UShort netbyte_order_port, const UInt netbyte_order_address);
	Address(const char* string);

	static Address GetHostAddress(const UShort port = 0);


	UInt Ipv4Address() const;
	void Ipv4Address(const UInt value);

	UShort Port() const;
	void Port(const UShort port);

	const char* AsString() const;
	static const char* AsString(const UInt ipv4address);
	static const char* AsStringNB(const UInt netbyte_order_address);

	UInt Size() const
	{
		return sizeof(data);
	};

	operator const SocketAddress* () const
	{
		return (SocketAddress*) &data;
	};

	operator SocketAddress* ()
	{
		return (SocketAddress*) &data;
	};

	operator const char* () const
	{
		return AsString();
	};


	friend std::ostream& operator << (std::ostream& os, const Address& addr);

private:

	InetSocketAddress data;
};


class Socket
{
public:

	Socket();
	Socket(const SocketHandle s);
	~Socket();

	bool CanBroadcast() const;
	void AllowBroadcast(bool value = true);
	bool CanReuseAddress() const;
	void AllowAddressReuse(bool value = true);

	UInt MaxMessageSize() const;

	Address BroadcastAddress(const UInt localaddress = INADDR_ANY) const;

	Address LocalAddress() const;

	void Bind(const Address& address);

	const Address& TargetAddress() const { return targetaddr; };
	void TargetAddress(const Address& address) { targetaddr = address; };
	UInt Send(const char* message, const UInt length);

	UInt BytesAvailable() const;

	const Address& SourceAddress() const { return sourceaddr; };
	UInt Receive(char* message, const UInt length);

	SocketHandle Handle() const { return handle; };
	operator SocketHandle () const { return handle; };

	bool CanRead() const;
	bool CanWrite() const;

	bool WaitUntilReadable(
			const UInt timeout_secs = 0xFFFFFFFF,
			const UInt timeout_usecs = 0xFFFFFFFF
		) const;

	bool WaitUntilWriteable(
			const UInt timeout_secs = 0xFFFFFFFF,
			const UInt timeout_usecs = 0xFFFFFFFF
		) const;

private:

	SocketHandle handle;
	Address targetaddr;
	Address sourceaddr;
};



class SocketList
{
public:

	void Add(Socket& s);
	void Remove(Socket& s);
	bool Contains(Socket& s);

	UInt Size() const { return list.size(); };

	const Socket* operator [] (UInt index) const;
	Socket* operator [] (UInt index);

	friend std::ostream& operator << (std::ostream& os, const SocketList& sl);

private:

	std::vector<Socket*> list;
};


enum SocketStatus
{
	SocketStatusUnknown,
	CanReadSocket,
	CanWriteSocket,
	SocketHasFailed
};

enum SocketEventTypes
{
	SocketReadable  = 0x01,
	SocketWriteable = 0x02,
	SocketFailure   = 0x04
};


Socket* WaitFor(
		SocketStatus& status,
		const SocketList& sockets,
		const int eventflags = SocketReadable | SocketWriteable | SocketFailure,
		const UInt timeout_secs = 0xFFFFFFFF,
		const UInt timeout_usecs = 0xFFFFFFFF
	);

Socket* WaitForReadable(
		const SocketList& sockets,
		const UInt timeout_secs = 0xFFFFFFFF,
		const UInt timeout_usecs = 0xFFFFFFFF
	);

Socket* WaitForWriteable(
		const SocketList& sockets,
		const UInt timeout_secs = 0xFFFFFFFF,
		const UInt timeout_usecs = 0xFFFFFFFF
	);


} // System
} // dHLT

#endif // dHLT_SYSTEM_SOCKET_HPP
