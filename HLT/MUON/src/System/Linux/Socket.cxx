////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "System/Socket.hpp"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/ioctl.h>
#include <linux/sockios.h>
#include <linux/if.h>
#include <errno.h>

#define SOCKET_ERROR -1
#define INVALID_SOCKET -1
#define closesocket close


namespace
{


/* the parameter 'buffer' should be set to NULL before calling this routine.
 */
struct ifreq* FindNetInterface(
		struct ifreq*& buffer, const dHLT::UInt ipv4address,
		const dHLT::System::SocketHandle handle
	)
{
	using dHLT::ThrowOutOfMemory;
	using dHLT::System::Address;
	using dHLT::System::SocketError;

	int structcount = 16;
	struct ifconf conf;
	
	// Allocate buffer space to contain interface descriptors.
	int buffersize = sizeof(struct ifreq) * structcount + 1;
	buffer = (struct ifreq*) malloc(buffersize);
	if (buffer == NULL) ThrowOutOfMemory();
retry:
	// Make a ioctl call to fill the buffers.
	conf.ifc_len = buffersize;
	conf.ifc_req = buffer;
	int result = ioctl(handle, SIOCGIFCONF, (char*) &conf);
	if (result == -1) throw SocketError();

	// If the ioctl call filled the whole buffer that means
	// there are more interfaces than we allocated for, so we
	// need to allocate more memory and try again.
	if (conf.ifc_len >= buffersize)
	{
		structcount = structcount * 2;
		buffersize = sizeof(struct ifreq) * structcount + 1;
		void* tmp = realloc(buffer, buffersize);
		if (tmp == NULL) ThrowOutOfMemory();
		buffer = (struct ifreq*) tmp;
		goto retry;
	};

#ifdef DEBUG
	DebugMsg(1, "================ Interfaces ==============");
	struct ifreq* lastifr_ = &buffer[ conf.ifc_len / sizeof(struct ifreq) ];
	for (struct ifreq* ifr_ = buffer; ifr_ < lastifr_ ; ifr_++)
	{
		Address address(ifr_->ifr_addr);
		DebugMsg(1, ((char*)&ifr_->ifr_name) << "\t" << address);
	};
#endif // DEBUG

	// Find the interface whose address corresponds to the sockets address.
	struct ifreq* lastifr = &buffer[ conf.ifc_len / sizeof(struct ifreq) ];
	struct ifreq* ifr;
	if (ipv4address != INADDR_ANY)
	{
		for (ifr = buffer; ifr < lastifr; ifr++)
		{
			Address ifaddress(ifr->ifr_addr);
			if ( ifaddress.Ipv4Address() == ipv4address )
				break;
		};
	}
	else
		ifr = buffer;
	// Indicate a 'device not found' error if we could not find corresponding
	// network interface for the given IPv4 address.
	if (ifr >= lastifr)
	{
		errno = ENODEV;
		throw SocketError();
	};
	
	return ifr;
};

} // end of namespace


namespace dHLT
{
namespace System
{

SocketError::SocketError() throw () : System::Error()
{
	if (errno == 0)
		errorcode = h_errno | 0x10000000;
	else
		errorcode = errno;
};

const char* SocketError::Message() const throw()
{
	if (0x10000000 & errorcode)
		return hstrerror(errorcode & ~0x10000000);
	else
		return AsString(errorcode);
};



Address Address::GetHostAddress(const UShort port)
{
	char buffer[1024];
	char* hostname = (char*) &buffer[0];

	int result = gethostname(hostname, sizeof(buffer));
	if (result == SOCKET_ERROR)
	{
		errno = 0;
		throw SocketError();
	};

	struct hostent* ent = gethostbyname(hostname);
	if (ent == NULL) throw SocketError();

	UInt addr;
	if ( ent->h_length >= (int)sizeof(addr) )
		addr = *( (UInt*)ent->h_addr_list[0] );
	else
		addr = INADDR_ANY;

	return Address(htons(port), addr);
};

////////////////////////////////////////////////////////////////////////////////

bool Socket::CanBroadcast() const
{
	bool value;
	socklen_t valuesize = sizeof(value);
	int result = getsockopt(handle, SOL_SOCKET, SO_BROADCAST, &value, &valuesize);
	if (result == SOCKET_ERROR) throw SocketError();
	return value;
};


void Socket::AllowBroadcast(bool value)
{
	int ivalue = value;
	int result = setsockopt(handle, SOL_SOCKET, SO_BROADCAST, &ivalue, sizeof(ivalue));
	if (result == SOCKET_ERROR) throw SocketError();
};


bool Socket::CanReuseAddress() const
{
	bool value;
	socklen_t valuesize = sizeof(value);
	int result = getsockopt(handle, SOL_SOCKET, SO_REUSEADDR, &value, &valuesize);
	if (result == SOCKET_ERROR) throw SocketError();
	return value;
};


void Socket::AllowAddressReuse(bool value)
{
	int ivalue = value;
	int result = setsockopt(handle, SOL_SOCKET, SO_REUSEADDR, &ivalue, sizeof(ivalue));
	if (result == SOCKET_ERROR) throw SocketError();
};


Address Socket::BroadcastAddress(const UInt interfaceaddress) const
{
	struct ifreq* buffer = NULL;
	Address broadcast_address;
	
	UInt ifaddr;
	if (interfaceaddress == INADDR_ANY)
		ifaddr = LocalAddress().Ipv4Address();
	else
		ifaddr = interfaceaddress;
	
	try
	{
		struct ifreq* ifr = FindNetInterface(buffer, ifaddr, handle);

		// We can now make a request for the broadcast address of the found interface.
		int result = ioctl(handle, SIOCGIFBRDADDR, (char*)ifr);
		if (result == -1) throw SocketError();
		broadcast_address = Address(ifr->ifr_broadaddr);
	}
	finally
	(
		if (buffer != NULL) free(buffer);
	);
	
	return broadcast_address;
};


UInt Socket::MaxMessageSize() const
{
	struct ifreq* buffer = NULL;
	int mtu;
	
	try
	{
		struct ifreq* ifr = FindNetInterface(buffer, LocalAddress().Ipv4Address(), handle);

		// We can now make a request for the MTU of the found interface.
		int result = ioctl(handle, SIOCGIFMTU, (char*)ifr);
		if (result == -1) throw SocketError();
		mtu = ifr->ifr_mtu;
	}
	finally
	(
		if (buffer != NULL) free(buffer);
	);
	
	return (UInt) mtu;
};


UInt Socket::BytesAvailable() const
{
	unsigned long bytesavailable;
	int result = ioctl(handle, FIONREAD, &bytesavailable);
	if (result == SOCKET_ERROR) throw SocketError();
	return bytesavailable;
};


} // System
} // dHLT

