////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_SYSTEM_LINUX_SYSTEM_TYPES_HPP
#define dHLT_SYSTEM_LINUX_SYSTEM_TYPES_HPP

#include <pthread.h>
#include <sys/socket.h>
#include <netinet/in.h>

#ifdef USE_GETDENTS_SYSCALL
#	include <linux/types.h>
#	include <linux/dirent.h>
#else // USE_GETDENTS_SYSCALL
#	include <dirent.h>
#endif // USE_GETDENTS_SYSCALL

namespace dHLT
{
namespace System
{


typedef pthread_t SystemThread;
typedef pthread_mutex_t SystemMutex;
typedef pthread_cond_t SystemMutexCondition;
typedef int SocketHandle;
typedef struct sockaddr SocketAddress;
typedef struct sockaddr_in InetSocketAddress;
typedef int SystemFile;

#ifdef USE_GETDENTS_SYSCALL
typedef int SystemDirectory;
#else // USE_GETDENTS_SYSCALL
typedef DIR* SystemDirectory;
#endif // USE_GETDENTS_SYSCALL



} // System
} // dHLT

#endif // dHLT_SYSTEM_LINUX_SYSTEM_TYPES_HPP
