////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "BCMP/Packets.hpp"
#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "Error.hpp"
#include <stdlib.h>


#ifdef DEBUG
namespace
{
	/* Perform an assertion check at program termination.
	   Assert that we had the same number of allocations as frees.
	 */
	class Check
	{
	public:
	
		Check()
		{
			alloccount = 0;
		};

		~Check()
		{
			Assert( alloccount == 0 );
		};

		long alloccount;
	} check;

};
#endif // DEBUG


namespace dHLT
{
namespace BCMP
{


std::ostream& operator << (std::ostream& os, const PacketType type)
{
	switch (type)
	{
	case TypeRegister:    os << "Register";    break;
	case TypeConfirm:     os << "Confirm";     break;
	case TypeData:        os << "Data";        break;
	case TypeAcknowledge: os << "Acknowledge"; break;
	case TypeRelease:     os << "Release";     break;
	case TypeQuit:        os << "Quit";        break;
	case TypeError:       os << "Error";       break;
	default:              os << "unknown";
	};
	return os;
};


char* AllocateMessage(const UInt size)
{
	DebugCode( check.alloccount++; )
	char* memory = (char*)malloc(size);
	if (memory == NULL) ThrowOutOfMemory();
	DebugMsg(10, "Allocated message of " << size << " bytes: " << (void*)memory );
	return memory;
};


char* ReallocateMessage(char* message, const UInt newsize)
{
	char* memory = (char*) realloc(message, newsize);
	if (memory == NULL) ThrowOutOfMemory();
	return memory;
};


void FreeMessage(char* message)
{
	DebugMsg(10, "Free message: " << (void*)message );
	free(message);
	DebugCode( check.alloccount--; )
};


} // BCMP
} // dHLT
