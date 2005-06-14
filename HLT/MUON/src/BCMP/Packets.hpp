////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_BCMP_PACKETS_HPP
#define dHLT_BCMP_PACKETS_HPP

#include "System/Socket.hpp"
#include "BasicTypes.hpp"
#include <ostream>

namespace dHLT
{
namespace BCMP
{

// BCMP = BroadCast Message Protocol

/* This protocol tries its best to deliver message packets over UDP broadcast,
   but it might deliver duplicate packets.
 */

enum PacketType
{
	TypeRegister = 1,
	TypeConfirm = 2,
	TypeData = 3,
	TypeAcknowledge = 4,
	TypeRelease = 5,
	TypeQuit = 6,
	TypeError = 7
};

std::ostream& operator << (std::ostream& os, const PacketType type);


struct RegisterPacket
{
	UChar type;
};

struct ConfirmPacket
{
	UChar type;
};

struct DataPacket
{
	UChar type;
	UChar packetid;
	char data[1];
};

struct AcknowledgePacket
{
	UChar type;
	UChar count;
	UChar packetid[1];
};

struct ReleasePacket
{
	UChar type;
};

struct QuitPacket
{
	UChar type;
};

struct ErrorPacket
{
	UChar type;
	Int errorcode;
	char message[1];
};


char* AllocateMessage(const UInt size);
char* ReallocateMessage(char* message, const UInt newsize);
void FreeMessage(char* message);


} // BCMP
} // dHLT

#endif // dHLT_BCMP_PACKETS_HPP
