////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_WIN32_SYSTEM_TYPES_HPP
#define dHLT_SYSTEM_WIN32_SYSTEM_TYPES_HPP

#include <windows.h>

namespace dHLT
{
namespace System
{


typedef HANDLE SystemThread;
typedef CRITICAL_SECTION SystemMutex;
typedef HANDLE SystemMutexCondition;
typedef SOCKET SocketHandle;
typedef struct sockaddr SocketAddress;
typedef struct sockaddr_in InetSocketAddress;
typedef HANDLE SystemFile;
typedef HANDLE SystemDirectory;


} // System
} // dHLT

#endif // dHLT_SYSTEM_WIN32_SYSTEM_TYPES_HPP
