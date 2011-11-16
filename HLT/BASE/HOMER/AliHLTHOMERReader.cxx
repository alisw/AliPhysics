/************************************************************************
**
**
** This file is property of and copyright by the Technical Computer
** Science Group, Kirchhoff Institute for Physics, Ruprecht-Karls-
** University, Heidelberg, Germany, 2001
** This file has been written by Timm Morten Steinbeck, 
** timm@kip.uni-heidelberg.de
**
**
** See the file license.txt for details regarding usage, modification,
** distribution and warranty.
** Important: This file is provided without any warranty, including
** fitness for any particular purpose.
**
**
** Newer versions of this file's package will be made available from 
** http://web.kip.uni-heidelberg.de/Hardwinf/L3/ 
** or the corresponding page of the Heidelberg Alice Level 3 group.
**
*************************************************************************/

/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

/** @file   AliHLTHOMERReader.cxx
    @author Timm Steinbeck
    @date   Sep 14 2007
    @brief  HLT Online Monitoring Environment including ROOT - Reader
    @note   migrated from PubSub HLT-stable-20070905.141318 (rev 2375)    */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTHOMERReader.h"
//#include <stdio.h>
#ifdef __SUNPRO_CC
#include <string.h>
#else
#include <cstring>
#endif
#include <cerrno>
#include <netdb.h>
//#include <sys/types.h>
//#include <sys/socket.h>
//#include <netinet/in.h>
//#include <netinet/tcp.h>
#include <unistd.h>
#ifndef __CYGWIN__
#include <rpc/types.h>
#endif
#include <fcntl.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#ifdef USE_ROOT
#include <Rtypes.h>
#endif


#define MOD_BIN "MOD BIN\n"
#define MOD_ASC "MOD ASC\n"
#define GET_ONE "GET ONE\n"
#define GET_ALL "GET ALL\n"

// MAXHOSTNAMELEN not defined on macosx
// 686-apple-darwin9-gcc-4.0.1
#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 64
#endif

#ifdef USE_ROOT
ClassImp(AliHLTMonitoringReader);
ClassImp(AliHLTHOMERReader);
#endif





#ifdef USE_ROOT
AliHLTHOMERReader::AliHLTHOMERReader()
  :
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// Reader implementation of the HOMER interface.
// The HLT Monitoring Environment including ROOT is
// a native interface to ship out data from the HLT chain.
// See pdf document shiped with the package
// for class documentation and tutorial.
    Init();
    }
#endif


AliHLTHOMERReader::AliHLTHOMERReader( const char* hostname, unsigned short port )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from a TCP port
    Init();
    if ( !AllocDataSources(1) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    fConnectionStatus = AddDataSource( hostname, port, fDataSources[0] );
    if ( fConnectionStatus )
	fErrorConnection = 0;
    else
	{
	fDataSourceCnt++;
	fTCPDataSourceCnt++;
	fDataSources[0].fNdx = 0;
	}
    }

AliHLTHOMERReader::AliHLTHOMERReader( unsigned int tcpCnt, const char** hostnames, const unsigned short* ports )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from multiple TCP ports
    Init();
    if ( !AllocDataSources(tcpCnt) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    for ( unsigned int n = 0; n < tcpCnt; n++, fDataSourceCnt++, fTCPDataSourceCnt++ )
	{
	fConnectionStatus = AddDataSource( hostnames[n], ports[n], fDataSources[n] );
	if ( fConnectionStatus )
	    {
	    fErrorConnection = n;
	    return;
	    }
	fDataSources[n].fNdx = n;
	}
    }

AliHLTHOMERReader::AliHLTHOMERReader( key_t shmKey, int shmSize )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from a System V shared memory segment
    Init();
    if ( !AllocDataSources(1) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    fConnectionStatus = AddDataSource( shmKey, shmSize, fDataSources[0] );
    if ( fConnectionStatus )
	fErrorConnection = 0;
    else
	{
	fDataSourceCnt++;
	fShmDataSourceCnt++;
	fDataSources[0].fNdx = 0;
	}
    }

AliHLTHOMERReader::AliHLTHOMERReader( unsigned int shmCnt, const key_t* shmKeys, const int* shmSizes )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from multiple System V shared memory segments
    Init();
    if ( !AllocDataSources(shmCnt) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    for ( unsigned int n = 0; n < shmCnt; n++, fDataSourceCnt++, fShmDataSourceCnt++ )
	{
	fConnectionStatus = AddDataSource( shmKeys[n], shmSizes[n], fDataSources[n] );
	if ( fConnectionStatus )
	    {
	    fErrorConnection = n;
	    return;
	    }
	fDataSources[n].fNdx = n;
	}
    }

AliHLTHOMERReader::AliHLTHOMERReader( unsigned int tcpCnt, const char** hostnames, const unsigned short* ports, 
			  unsigned int shmCnt, const key_t* shmKeys, const int* shmSizes )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from multiple TCP ports and multiple System V shared memory segments
    Init();
    if ( !AllocDataSources(tcpCnt+shmCnt) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    for ( unsigned int n = 0; n < tcpCnt; n++, fDataSourceCnt++, fTCPDataSourceCnt++ )
	{
	fConnectionStatus = AddDataSource( hostnames[n], ports[n], fDataSources[n] );
	if ( fConnectionStatus )
	    {
	    fErrorConnection = n;
	    return;
	    }
	fDataSources[n].fNdx = n;
	}
    for ( unsigned int n = 0; n < shmCnt; n++, fDataSourceCnt++, fShmDataSourceCnt++ )
	{
	fConnectionStatus = AddDataSource( shmKeys[n], shmSizes[n], fDataSources[tcpCnt+n] );
	if ( fConnectionStatus )
	    {
	    fErrorConnection = tcpCnt+n;
	    return;
	    }
	fDataSources[n].fNdx = n;
	}
    }

AliHLTHOMERReader::AliHLTHOMERReader( const void* pBuffer, int size )
  :
  AliHLTMonitoringReader(),
  TObject(),
  fCurrentEventType(~(homer_uint64)0),
  fCurrentEventID(~(homer_uint64)0),
  fBlockCnt(0),
  fMaxBlockCnt(0),
  fBlocks(NULL),
  fDataSourceCnt(0),
  fTCPDataSourceCnt(0),
  fShmDataSourceCnt(0),
  fDataSourceMaxCnt(0),
  fDataSources(NULL),
  fConnectionStatus(0),
  fErrorConnection(~(unsigned int)0),
  fEventRequestAdvanceTime(0)
    {
// see header file for class documentation
// For reading from a System V shared memory segment
    Init();
    if ( !AllocDataSources(1) )
	{
	fErrorConnection = 0;
	fConnectionStatus = ENOMEM;
	return;
	}
    fConnectionStatus = AddDataSource(const_cast<void*>(pBuffer), size, fDataSources[0] );
    if ( fConnectionStatus )
	fErrorConnection = 0;
    else
	{
	fDataSourceCnt++;
	fShmDataSourceCnt++;
	fDataSources[0].fNdx = 0;
	}
    }

AliHLTHOMERReader::~AliHLTHOMERReader()
    {
// see header file for class documentation
    ReleaseCurrentEvent();
    FreeDataSources();

    if (fDataSources)
      delete [] fDataSources;
    }

int  AliHLTHOMERReader::ReadNextEvent()
    {
// see header file for class documentation
// Read in the next available event
    return ReadNextEvent( false, 0 );
    }

int AliHLTHOMERReader::ReadNextEvent( unsigned long timeout )
    {
// see header file for class documentation
// Read in the next available event
    return ReadNextEvent( true, timeout );
    }

unsigned long AliHLTHOMERReader::GetBlockDataLength( unsigned long ndx ) const
    {
// see header file for class documentation
// Return the size (in bytes) of the current event's data
// block with the given block index (starting at 0).
    if ( ndx >= fBlockCnt )
	return 0;
    return fBlocks[ndx].fLength;
    }

const void* AliHLTHOMERReader::GetBlockData( unsigned long ndx ) const
    {
// see header file for class documentation
// Return a pointer to the start of the current event's data
// block with the given block index (starting at 0).
    if ( ndx >= fBlockCnt )
	return NULL;
    return fBlocks[ndx].fData;
    }

const char* AliHLTHOMERReader::GetBlockSendNodeID( unsigned long ndx ) const
    {
// see header file for class documentation
// Return IP address or hostname of node which sent the 
// current event's data block with the given block index 
// (starting at 0).
// For HOMER this is the ID of the node on which the subscriber 
// that provided this data runs/ran.
    if ( ndx >= fBlockCnt )
	return NULL;
#ifdef DEBUG
    if ( fBlocks[ndx].fSource >= fDataSourceCnt )
	{
	fprintf( stderr, "%s:%d: Internal Error: fBlocks[ndx].fSource (%u) >= fDataSourceCnt (%u)\n",
		 __FILE__, __LINE__, fBlocks[ndx].fSource, fDataSourceCnt );
	return NULL;
	}
#endif
    return fDataSources[ fBlocks[ndx].fSource ].fHostname;
    //return fBlocks[ndx].fOriginatingNodeID;
    }

homer_uint8 AliHLTHOMERReader::GetBlockByteOrder( unsigned long ndx ) const
    {
// see header file for class documentation
// Return byte order of the data stored in the 
// current event's data block with the given block index (starting at 0). 
//	   0 is unknown alignment, 
//	   1 ist little endian, 
//	   2 is big endian. */
    if ( ndx >= fBlockCnt )
	return 0;
    //return ((AliHLTRIBlockDescriptorV1*)fBlocks[ndx].fMetaData)->fType.fID;
    return *(((homer_uint8*)fBlocks[ndx].fMetaData)+kByteOrderAttribute_8b_Offset);
    }

homer_uint8 AliHLTHOMERReader::GetBlockTypeAlignment( unsigned long ndx, homer_uint8 dataType ) const
    {
// see header file for class documentation
// Return the alignment (in bytes) of the given datatype 
// in the data stored in the current event's data block
// with the given block index (starting at 0). 
// Possible values for the data type are
//	   0: homer_uint64
//	   1: homer_uint32
//	   2: uin16
//	   3: homer_uint8
//	   4: double
//	   5: float
    if ( ndx >= fBlockCnt )
	return 0;
    if ( dataType > (kFloatAlignment_8b_Offset-kAlignment_8b_StartOffset) )
	return 0;
    //return ((AliHLTRIBlockDescriptorV1*)fBlocks[ndx].fMetaData)->fType.fID;
    return *(((homer_uint8*)fBlocks[ndx].fMetaData)+kAlignment_8b_StartOffset+dataType);
    }

homer_uint64 AliHLTHOMERReader::GetBlockStatusFlags( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return 0;
    return *(((homer_uint64*)fBlocks[ndx].fMetaData)+kStatusFlags_64b_Offset);
    }

/* HOMER specific */
/* Return the type of the data in the current event's data
   block with the given block index (starting at 0). */
homer_uint64 AliHLTHOMERReader::GetBlockDataType( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return ~(homer_uint64)0;
    //return ((AliHLTRIBlockDescriptorV1*)fBlocks[ndx].fMetaData)->fType.fID;
    return *(((homer_uint64*)fBlocks[ndx].fMetaData)+kType_64b_Offset);
    }

/* Return the origin of the data in the current event's data
   block with the given block index (starting at 0). */
homer_uint32 AliHLTHOMERReader::GetBlockDataOrigin( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return ~(homer_uint32)0;
    //return (homer_uint32)( ((AliHLTRIBlockDescriptorV1*)fBlocks[ndx].fMetaData)->fSubType1.fID );
    return (homer_uint32)*(((homer_uint64*)fBlocks[ndx].fMetaData)+kSubType1_64b_Offset);
    }

/* Return a specification of the data in the current event's data
   block with the given block index (starting at 0). */
homer_uint32 AliHLTHOMERReader::GetBlockDataSpec( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return ~(homer_uint32)0;
    //return (homer_uint32)( ((AliHLTRIBlockDescriptorV1*)fBlocks[ndx].fMetaData)->fSubType2.fID );
    return (homer_uint32)*(((homer_uint64*)fBlocks[ndx].fMetaData)+kSubType2_64b_Offset);
    }

homer_uint64 AliHLTHOMERReader::GetBlockBirthSeconds( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return 0;
    return *(((homer_uint64*)fBlocks[ndx].fMetaData)+kBirth_s_64b_Offset);
    }

homer_uint64 AliHLTHOMERReader::GetBlockBirthMicroSeconds( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return 0;
    return *(((homer_uint64*)fBlocks[ndx].fMetaData)+kBirth_us_64b_Offset);
    }

/* Find the next data block in the current event with the given
   data type, origin, and specification. Returns the block's 
   index. */
unsigned long AliHLTHOMERReader::FindBlockNdx( homer_uint64 type, homer_uint32 origin, 
					 homer_uint32 spec, unsigned long startNdx ) const
    {
// see header file for class documentation
    for ( unsigned long n=startNdx; n < fBlockCnt; n++ )
	{
	if ( ( type == 0xFFFFFFFFFFFFFFFFULL || *(((homer_uint64*)fBlocks[n].fMetaData)+kType_64b_Offset)==type ) &&
	     ( origin == 0xFFFFFFFF || (homer_uint32)( *(((homer_uint64*)fBlocks[n].fMetaData)+kSubType1_64b_Offset) )==origin ) &&
	     ( spec == 0xFFFFFFFF || (homer_uint32)( *(((homer_uint64*)fBlocks[n].fMetaData)+kSubType2_64b_Offset) )==spec ) )
	    return n;
	}
    return ~(unsigned long)0;
    }

/* Find the next data block in the current event with the given
   data type, origin, and specification. Returns the block's 
   index. */
unsigned long AliHLTHOMERReader::FindBlockNdx( char type[8], char origin[4], 
					 homer_uint32 spec, unsigned long startNdx ) const
    {
// see header file for class documentation
    for ( unsigned long n=startNdx; n < fBlockCnt; n++ )
	{
	bool found1=true, found2=true;
	for ( unsigned i = 0; i < 8; i++ )
	    {
	    if ( type[i] != (char)0xFF )
		{
		found1=false;
		break;
		}
	    }
	if ( !found1 )
	    {
	    found1 = true;
	    for ( unsigned i = 0; i < 8; i++ )
		{
		//printf( "%u: Comparing type '%c' and '%c'\n", i, ((char*)(((homer_uint64*)fBlocks[n].fMetaData)+kType_64b_Offset))[i], type[i] );
		if ( ((char*)(((homer_uint64*)fBlocks[n].fMetaData)+kType_64b_Offset))[i] != type[i] )
		    {
		    found1=false;
		    break;
		    }
		}
	    }
	for ( unsigned i = 0; i < 4; i++ )
	    {
	    if ( origin[i] != (char)0xFF )
		{
		found2 = false;
		break;
		}
	    }
	if ( !found2 )
	    {
	    found2 = true;
	    for ( unsigned i = 0; i < 4; i++ )
		{
		//printf( "Comparing origin '%c' and '%c'\n", ((char*)(((homer_uint64*)fBlocks[n].fMetaData)+kSubType1_64b_Offset))[i], origin[i] );
		if ( ((char*)(((homer_uint64*)fBlocks[n].fMetaData)+kSubType1_64b_Offset))[i] != origin[i] )
		    {
		    found2=false;
		    break;
		    }
		}
	    }
	//printf( "Comparing spec '0x%08lX' and '0x%08lX'\n", (homer_uint32)( *(((homer_uint64*)fBlocks[n].fMetaData)+kSubType2_64b_Offset) ), spec );
	if ( found1 && found2 &&
	     ( spec == 0xFFFFFFFF || ((homer_uint32)( *(((homer_uint64*)fBlocks[n].fMetaData)+kSubType2_64b_Offset)) )==spec ) )
	    return n;
	}
    return ~(unsigned long)0;
    }

/* Return the ID of the node that actually produced this data block.
   This may be different from the node which sent the data to this
   monitoring object as returned by GetBlockSendNodeID. */
const char* AliHLTHOMERReader::GetBlockCreateNodeID( unsigned long ndx ) const
    {
// see header file for class documentation
    if ( ndx >= fBlockCnt )
	return NULL;
    return fBlocks[ndx].fOriginatingNodeID;
    }


void AliHLTHOMERReader::Init()
    {
// see header file for class documentation
    fCurrentEventType = ~(homer_uint64)0;
    fCurrentEventID = ~(homer_uint64)0;
    fMaxBlockCnt = fBlockCnt = 0;
    fBlocks = NULL;
	
    fDataSourceMaxCnt = fDataSourceCnt = fTCPDataSourceCnt = fShmDataSourceCnt = 0;
    fDataSources = NULL;

	
    fConnectionStatus = 0;
    fErrorConnection = ~(unsigned int)0;

    fEventRequestAdvanceTime = 0;
    }
	
bool AliHLTHOMERReader::AllocDataSources( unsigned int sourceCnt )
    {
// see header file for class documentation
    fDataSources = new DataSource[ sourceCnt ];
    if ( !fDataSources )
	return false;
    memset(fDataSources, 0, sizeof(DataSource)*sourceCnt);
    fDataSourceCnt = 0;
    fDataSourceMaxCnt = sourceCnt;
    return true;
    }

int AliHLTHOMERReader::AddDataSource( const char* hostname, unsigned short port, DataSource& source )
    {
// see header file for class documentation
    struct hostent* he;
    he = gethostbyname( hostname );
    if ( he == NULL )
	{
	//fprintf( stderr, "Unable to determine remote host address from '%s'.\n", hostname );
	return EADDRNOTAVAIL;
	}

    struct sockaddr_in remoteAddr;
    remoteAddr.sin_family = AF_INET;    // host byte order 
    remoteAddr.sin_port = htons(port);  // short, network byte order 
    remoteAddr.sin_addr = *((struct in_addr *)he->h_addr);
    memset(&(remoteAddr.sin_zero), '\0', 8);  // zero the rest of the struct

    // Create socket and connect to target program on remote node
    source.fTCPConnection = socket( AF_INET, SOCK_STREAM, 0 );
    if ( source.fTCPConnection == -1 )
	{
	return errno;
        }

    int ret;

    ret = connect( source.fTCPConnection, (struct sockaddr *)&remoteAddr, sizeof(struct sockaddr) );
    if ( ret == -1 )
	{
	ret=errno;
	close( source.fTCPConnection );
	return ret;
	} 

    ret = write( source.fTCPConnection, MOD_BIN, strlen(MOD_BIN) );
    if ( ret != (int)strlen(MOD_BIN) )
	{
	ret=errno;
	close( source.fTCPConnection );
	return ret;
	}

    unsigned hostnamelen=strlen( hostname );
    char* tmpchar = new char[ hostnamelen+1 ];
    if ( !tmpchar )
	{
	close( source.fTCPConnection );
	return ENOMEM;
	}
    strncpy( tmpchar, hostname, hostnamelen );
    tmpchar[hostnamelen]=0;
    source.fHostname = tmpchar;

    source.fType = kTCP;
    source.fTCPPort = port;
    source.fData = NULL;
    source.fDataSize = 0;
    source.fDataRead = 0;
    return 0;
    }

int AliHLTHOMERReader::AddDataSource( key_t shmKey, int shmSize, DataSource& source )
    {
// see header file for class documentation
    int ret;
    char* tmpchar = new char[ MAXHOSTNAMELEN+1 ];
    if ( !tmpchar )
	{
	return ENOMEM;
	}
    gethostname( tmpchar, MAXHOSTNAMELEN );
    tmpchar[MAXHOSTNAMELEN]=(char)0;
    source.fHostname = tmpchar;

    source.fShmID = shmget( shmKey, shmSize, 0660 );
    if ( source.fShmID == -1 )
	{
	ret = errno;
	delete [] source.fHostname;
	return ret;
	}
    
    source.fShmPtr = (void*)shmat( source.fShmID, NULL, 0 );

    if ( !source.fShmPtr )
	{
	ret = errno;
	shmctl( source.fShmID, IPC_RMID, NULL );
	delete [] source.fHostname;
	return ret;
	}

    source.fType = kShm;
    source.fShmKey = shmKey;
    source.fShmSize = shmSize;
    source.fDataSize = 0;
    source.fDataRead = 0;
    return 0;
    }

int AliHLTHOMERReader::AddDataSource( void* pBuffer, int size, DataSource& source )
    {
// see header file for class documentation
// a buffer data source is like a shm source apart from the shm attach and detach
// procedure. Furthermore, the size indicator at the beginning of the buffer is not
// cleared right before sources are read but after the reading.
    //int ret;
    if ( !pBuffer || size<=0) return EINVAL;

    char* tmpchar = new char[ MAXHOSTNAMELEN+1 ];
    if ( !tmpchar )
	{
	return ENOMEM;
	}
    gethostname( tmpchar, MAXHOSTNAMELEN );
    tmpchar[MAXHOSTNAMELEN]=(char)0;
    source.fHostname = tmpchar;

    source.fShmID = -1;
    // the data buffer does not contain a size indicator in the first 4 bytes
    // like the shm source buffer. Still we want to use the mechanism to invalidate/
    // trigger by clearing the size indicator. Take the source.fShmSize variable.
    source.fShmPtr = &source.fShmSize;
    source.fType = kBuf;
    source.fShmKey = 0;
    source.fShmSize = size;
    source.fData = pBuffer;
    source.fDataSize = 0;
    source.fDataRead = 0;
    return 0;
    }

void AliHLTHOMERReader::FreeDataSources()
    {
// see header file for class documentation
    for ( unsigned n=0; n < fDataSourceCnt; n++ )
	{
	if ( fDataSources[n].fType == kTCP )
	    FreeTCPDataSource( fDataSources[n] );
	else if ( fDataSources[n].fType == kShm )
	    FreeShmDataSource( fDataSources[n] );
	if ( fDataSources[n].fHostname )
	  delete [] fDataSources[n].fHostname;
	}
    fDataSourceCnt=0;
    }

int AliHLTHOMERReader::FreeShmDataSource( DataSource& source )
    {
// see header file for class documentation
    if ( source.fShmPtr )
      shmdt( (char*)source.fShmPtr );
//     if ( source.fShmID != -1 )
// 	shmctl( source.fShmID, IPC_RMID, NULL );
    return 0;
    }

int AliHLTHOMERReader::FreeTCPDataSource( DataSource& source )
    {
// see header file for class documentation
    if ( source.fTCPConnection )
	close( source.fTCPConnection );
    return 0;
    }

int AliHLTHOMERReader::ReadNextEvent( bool useTimeout, unsigned long timeout )
    {
// see header file for class documentation
    if ( fDataSourceCnt<=0 )
	return ENXIO;
    // Clean up currently active event.
    ReleaseCurrentEvent();
    int ret=0;
    // Trigger all configured data sources
    for ( unsigned n = 0; n<fDataSourceCnt; n++ ){
      if ( fDataSources[n].fType == kTCP )
	ret = TriggerTCPSource( fDataSources[n], useTimeout, timeout );
      else if ( fDataSources[n].fType == kShm )
	ret = TriggerShmSource( fDataSources[n], useTimeout, timeout );
      if ( ret )
	{
	  fErrorConnection = n;
	  fConnectionStatus=ret;
	  return fConnectionStatus;
	}
    }
    // Now read in data from the configured data source
    ret = ReadDataFromTCPSources( fTCPDataSourceCnt, fDataSources, useTimeout, timeout );
    
    if ( ret )
      {
	return ret;
      }
    ret = ReadDataFromShmSources( fShmDataSourceCnt, fDataSources+fTCPDataSourceCnt, useTimeout, timeout );
    if ( ret )
	{
	return ret;
	}
//     for ( unsigned n = 0; n<fDataSourceCnt; n++ )
// 	{
// 	if ( fDataSources[n].fType == kTCP )
// 	    ret = ReadDataFromTCPSource( fDataSources[n], useTimeout, timeout );
// 	else
// 	    ret = ReadDataFromShmSource( fDataSources[n], useTimeout, timeout );
// 	if ( ret )
// 	    {
// 	    fErrorConnection = n;
// 	    fConnectionStatus=ret;
// 	    return fConnectionStatus;
// 	    }
// 	}
    //Check to see that all sources contributed data for the same event
    homer_uint64 eventID;
    homer_uint64 eventType;
    if (!fDataSources[0].fData)
      {
	fErrorConnection = 0;
	fConnectionStatus=56;//ENOBUF;
	return fConnectionStatus;
      }
    eventID = GetSourceEventID( fDataSources[0] );
    eventType = GetSourceEventType( fDataSources[0] );
    for ( unsigned n = 1; n < fDataSourceCnt; n++ )
	{
	if ( !fDataSources[n].fData || GetSourceEventID( fDataSources[n] ) != eventID || GetSourceEventType( fDataSources[n] ) != eventType )
	    {
	    fErrorConnection = n;
	    fConnectionStatus=56;//EBADRQC;
	    return fConnectionStatus;
	    }
	}
    // Find all the different data blocks contained in the data from all
    // the sources.
    for ( unsigned n = 0; n < fDataSourceCnt; n++ )
	{
	ret = ParseSourceData( fDataSources[n] );
	if ( ret )
	    {
	    fErrorConnection = n;
	    fConnectionStatus=57;//EBADSLT;
	    return ret;
	    }
	}
    fCurrentEventID = eventID;
    fCurrentEventType = eventType;
    return 0;
    }

void AliHLTHOMERReader::ReleaseCurrentEvent()
    {
// see header file for class documentation
    // sources.fDataRead = 0;
    // fMaxBlockCnt
    fCurrentEventID = ~(homer_uint64)0;
    fCurrentEventType = ~(homer_uint64)0;
    for ( unsigned n = 0; n < fDataSourceCnt; n++ )
	{
	if ( fDataSources[n].fData )
	    {
	    if ( fDataSources[n].fType == kTCP )
		delete [] (homer_uint8*)fDataSources[n].fData;
	    // do not reset the data pointer for kBuf sources since this
	    // can not be set again.
	    if ( fDataSources[n].fType != kBuf )
	      fDataSources[n].fData = NULL;
	    }
	fDataSources[n].fDataSize = fDataSources[n].fDataRead = 0;
	}
    if ( fBlocks )
	{
	for ( unsigned n = 0; n < fMaxBlockCnt; n++ )
	    {
	    if ( fBlocks[n].fOriginatingNodeID )
		delete [] fBlocks[n].fOriginatingNodeID;
	    }
	delete [] fBlocks;
	fBlocks=0;
	fMaxBlockCnt = 0;
	fBlockCnt=0;
	}
    }

int AliHLTHOMERReader::TriggerTCPSource( DataSource& source, bool useTimeout, unsigned long timeoutUsec )
    {
// see header file for class documentation
    int ret=0;
    struct timeval oldSndTO, newSndTO;
    memset(&oldSndTO, 0, sizeof(oldSndTO));
    memset(&newSndTO, 0, sizeof(newSndTO));
    if ( useTimeout )
	{
	socklen_t optlen=sizeof(oldSndTO);
	ret = getsockopt( source.fTCPConnection, SOL_SOCKET, SO_SNDTIMEO, &oldSndTO, &optlen );
	if ( ret )
	    {
	    return errno;
	    }
	if ( optlen!=sizeof(oldSndTO) )
	    {
	    return ENXIO;
	    }
	newSndTO.tv_sec = timeoutUsec / 1000000;
	newSndTO.tv_usec = timeoutUsec - (newSndTO.tv_sec*1000000);
	ret = setsockopt( source.fTCPConnection, SOL_SOCKET, SO_SNDTIMEO, &newSndTO, sizeof(newSndTO) );
	if ( ret )
	    {
	    return errno;
	    }
	}
    // Send one event request
    if ( !fEventRequestAdvanceTime )
	{
	ret = write( source.fTCPConnection, GET_ONE, strlen(GET_ONE) );
	
	if ( ret != (int)strlen(GET_ONE) )
	    {
	    ret=errno;
	    setsockopt( source.fTCPConnection, SOL_SOCKET, SO_SNDTIMEO, &oldSndTO, sizeof(oldSndTO) );
	    return ret;
	    }
	}
    else
	{
	char tmpCmd[ 128 ];

	int len = snprintf( tmpCmd, 128, "FIRST ORBIT EVENT 0x%llu\n", (unsigned long long)fEventRequestAdvanceTime );
	if ( len>128 || len<0 )
	    {
	    ret=EMSGSIZE;
	    setsockopt( source.fTCPConnection, SOL_SOCKET, SO_SNDTIMEO, &oldSndTO, sizeof(oldSndTO) );
	    return ret;
	    }
	
	ret = write( source.fTCPConnection, tmpCmd, strlen(tmpCmd) );
	
	if ( ret != (int)strlen(tmpCmd) )
	    {
	    ret=errno;
	    setsockopt( source.fTCPConnection, SOL_SOCKET, SO_SNDTIMEO, &oldSndTO, sizeof(oldSndTO) );
	    return ret;
	    }
	
	}
    return 0;
    }

int AliHLTHOMERReader::TriggerShmSource( DataSource& source, bool, unsigned long ) const
    {
// see header file for class documentation
// clear the size indicator in the first 4 bytes of the buffer to request data
// from the HOMER writer.
    if ( source.fShmPtr )
	{
	*(homer_uint32*)( source.fShmPtr ) = 0;
	return 0;
	}
    else
	return EFAULT;
    }

int AliHLTHOMERReader::ReadDataFromTCPSources( unsigned sourceCnt, DataSource* sources, bool useTimeout, unsigned long timeout )
    {
// see header file for class documentation
    bool toRead = false;
    do
	{
	fd_set conns;
	FD_ZERO( &conns );
	int highestConn=0;
	toRead = false;
	unsigned firstConnection=~(unsigned)0;
	for ( unsigned long n = 0; n < sourceCnt; n++ )
	    {
	    if ( sources[n].fDataSize == 0 // size specifier not yet read
		 || sources[n].fDataRead < sources[n].fDataSize ) // Data not yet read fully
		{
		toRead = true;
		FD_SET( sources[n].fTCPConnection, &conns );
		if ( sources[n].fTCPConnection > highestConn )
		    highestConn = sources[n].fTCPConnection;
		fcntl( sources[n].fTCPConnection, F_SETFL, O_NONBLOCK );
		if ( firstConnection == ~(unsigned)0 )
		    firstConnection = n;
		}
	    else
		{
		fcntl( sources[n].fTCPConnection, F_SETFL, 0 );
		}
	    }
	if ( toRead )
	    {
	    struct timeval tv, *ptv;
	    if ( useTimeout )
		{
		tv.tv_sec = timeout / 1000000;
		tv.tv_usec = timeout - (tv.tv_sec*1000000);
		ptv = &tv;
		}
	    else
		ptv = NULL;
	    // wait until something is ready to be read
	    // either for timeout usecs or until eternity
	    int ret;
	    ret = select( highestConn+1, &conns, NULL, NULL, ptv ); 
	    if ( ret <=0 )
		{
		fErrorConnection = firstConnection;
		if ( errno )
		    fConnectionStatus = errno;
		else
		    fConnectionStatus = ETIMEDOUT;
		return fConnectionStatus;
		}
	    for ( unsigned n = 0; n < sourceCnt; n++ )
		{
		if ( FD_ISSET( sources[n].fTCPConnection, &conns ) )
		    {
		    if ( sources[n].fDataSize == 0 )
			{
			ret=read( sources[n].fTCPConnection, &(sources[n].fDataSize), sizeof(homer_uint32) );
			if ( ret != sizeof(homer_uint32) )
			    {
			    fErrorConnection = n;
			    if ( errno )
				fConnectionStatus = errno;
			    else
				fConnectionStatus = ENOMSG;
			    return fConnectionStatus;
			    }
			sources[n].fDataSize = ntohl( sources[n].fDataSize );
			sources[n].fDataRead = 0;
			sources[n].fData = new homer_uint8[ sources[n].fDataSize ];
			if ( !sources[n].fData )
			    {
			    fErrorConnection = n;
			    fConnectionStatus = ENOMEM;
			    return fConnectionStatus;
			    }
			}
		    else if ( sources[n].fData && sources[n].fDataRead < sources[n].fDataSize)
			{
			ret=read( sources[n].fTCPConnection, ((homer_uint8*)sources[n].fData)+sources[n].fDataRead, sources[n].fDataSize-sources[n].fDataRead );
			if ( ret>0 )
			    sources[n].fDataRead += ret;
			else if ( ret == 0 )
			    {
			    fErrorConnection = n;
			    fConnectionStatus = ECONNRESET;
			    return fConnectionStatus;
			    }
			else
			    {
			    fErrorConnection = n;
			    fConnectionStatus = errno;
			    return fConnectionStatus;
			    }
			}
		    else
			{
			fErrorConnection = n;
			fConnectionStatus = ENXIO;
			return fConnectionStatus;
			}
		    }
		}
	    }
	}
    while ( toRead );
    return 0;
    }

/*
int AliHLTHOMERReader::ReadDataFromTCPSources( DataSource& source, bool useTimeout, unsigned long timeout )
    {
#warning TODO If useTimeout: Set sockets to nonblocking, select + loop around GET_ONE write
    // Send one event request
    ret = write( source.fTCPConnection, GET_ONE, strlen(GET_ONE) );
    if ( ret != strlen(GET_ONE) )
	{
	return errno;
	}
    // wait for and read back size specifier
    unsigned sizeNBO;
    // The value transmitted is binary, in network byte order
    ret = read( source.fTCPConnection, &sizeNBO, sizeof(sizeNBO) );
    if ( ret != sizeof(sizeNBO) )
	{
	return errno;
	}
    // Convert back to host byte order
    source.fDataSize = ntohl( sizeNBO );
    source.fData = new homer_uint8[ source.fDataSize ];
    unsigned long dataRead=0, toRead;
    if ( !source.fData )
	{
	char buffer[1024];
	// Read in data into buffer in order not to block connection
	while ( dataRead < source.fDataSize )
	    {
	    if ( source.fDataSize-dataRead > 1024 )
		toRead = 1024;
	    else
		toRead = source.fDataSize-dataRead;
	    ret = read( source.fTCPConnection, buffer, toRead );
	    if ( ret > 0 )
		dataRead += ret;
	    else
		return errno;
	    }
	return ENOMEM;
	}
    while ( dataRead < source.fDataSize )
	{
	toRead = source.fDataSize-dataRead;
	ret = read( source.fTCPConnection, source.fData+dataRead, toRead );
	if ( ret > 0 )
	    dataRead += ret;
	else if ( ret == 0 && useTimeout )
	    {
	    struct timeval tv;
	    tv.tv_sec = timeout / 1000000;
	    tv.tv_usec = timeout - (tv.tv_sec*1000000);
	    fd_set conns;
	    FD_ZERO( &conns );
	    FD_SET( source.fTCPConnection, &conns );
	    ret = select( source.fTCPConnection+1, &conns, NULL, NULL );
	    if ( ret <=0 )
		return errno;
	    }
	else if ( ret == 0 )
	    {
	    if ( errno == EOK )
		return ECONNRESET;
	    else
		return errno;
	    }
	else
	    {
	    return errno;
	    }
	}
    return 0;
    }
*/

/*
int AliHLTHOMERReader::ReadDataFromShmSource( DataSource& source, bool useTimeout, unsigned long timeout )
    {
    
    }
*/

int AliHLTHOMERReader::ReadDataFromShmSources( unsigned sourceCnt, DataSource* sources, bool useTimeout, unsigned long timeout )
    {
// see header file for class documentation
    struct timeval tv1, tv2;
    bool found=false;
    bool all=true;
    if ( useTimeout )
	gettimeofday( &tv1, NULL );
    do
	{
	found = false;
	all = true;
	for ( unsigned n = 0; n < sourceCnt; n++ )
	    {
	    if ( !sources[n].fDataSize )
		all = false;
	    if ( sources[n].fShmPtr && *(homer_uint32*)sources[n].fShmPtr>0 && !sources[n].fDataSize )
		{
		found = true;
		sources[n].fDataSize = *(homer_uint32*)sources[n].fShmPtr;
		if (sources[n].fType==kBuf)
		  {
		    // the data buffer is already set to fData, just need to set fDataSize member
		    // to invalidate after the first reading. Subsequent calls to ReadNextEvent return 0
		    TriggerShmSource( sources[n], 0, 0 );
		  } else 
		  {
		    sources[n].fData = ((homer_uint8*)sources[n].fShmPtr)+sizeof(homer_uint32);
		  }
		}
	    }
	if ( found && useTimeout )
	    gettimeofday( &tv1, NULL );
	if ( !all && useTimeout )
	    {
	    gettimeofday( &tv2, NULL );
	    unsigned long long tdiff;
	    tdiff = tv2.tv_sec-tv1.tv_sec;
	    tdiff *= 1000000;
	    tdiff += tv2.tv_usec-tv1.tv_usec;
	    if ( tdiff > timeout )
		return ETIMEDOUT;
	    }
	if ( !all )
	    usleep( 0 );
	}
    while ( !all );
    return 0;
    }

int AliHLTHOMERReader::ParseSourceData( const DataSource& source )
    {
// see header file for class documentation
    if ( source.fData )
	{
	homer_uint8 sourceByteOrder = ((homer_uint8*)source.fData)[ kByteOrderAttribute_8b_Offset ];
	if (sourceByteOrder!=kHOMERLittleEndianByteOrder && sourceByteOrder!=kHOMERBigEndianByteOrder) return EBADMSG;
	homer_uint64 blockCnt = Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)source.fData)[ kSubType2_64b_Offset ] );
	// block count is not related to size of the data in the way the
	// following condition implies. But we can at least limit the block
	// count for the case the data is corrupted
	if (blockCnt>source.fDataSize) return EBADMSG;
	int ret=ReAllocBlocks( fMaxBlockCnt+blockCnt );
	if ( ret )
	    return ret;
	homer_uint64 descrOffset = Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)source.fData)[ kOffset_64b_Offset ] );
	for ( homer_uint64 n = 0; n < blockCnt && fBlockCnt < fMaxBlockCnt; n++, fBlockCnt++ )
	    {
	    if (descrOffset+kLength_64b_Offset>=source.fDataSize) return EBADMSG;
	    homer_uint8* descr = ((homer_uint8*)source.fData)+descrOffset;
	    unsigned descrLen = Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)descr)[ kLength_64b_Offset ] );
	    if (descrOffset+descrLen>=source.fDataSize) return EBADMSG;
	    if (Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)descr)[ kID_64b_Offset ] ) != HOMER_BLOCK_DESCRIPTOR_TYPEID) return 126/*ENOKEY*/;
	    fBlocks[fBlockCnt].fSource = source.fNdx;
	    fBlocks[fBlockCnt].fData = ((homer_uint8*)source.fData) + Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)descr)[ kOffset_64b_Offset ] );
	    fBlocks[fBlockCnt].fLength = Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)descr)[ kSize_64b_Offset ] );
	    fBlocks[fBlockCnt].fMetaData = (homer_uint64*)descr;
	    struct in_addr tmpA;
	    tmpA.s_addr = (homer_uint32)( ((homer_uint64*)descr)[ kProducerNode_64b_Offset ] );
	    char* addr = inet_ntoa( tmpA );
	    unsigned straddrlen=strlen(addr);
	    char* tmpchar = new char[ straddrlen+1 ];
	    if ( !tmpchar )
		return ENOMEM;
	    strncpy( tmpchar, addr, straddrlen );
	    tmpchar[straddrlen]=0;
	    fBlocks[fBlockCnt].fOriginatingNodeID = tmpchar;
	    descrOffset += descrLen;
	    }
	return 0;
	}
    return EFAULT;
    }
	
int AliHLTHOMERReader::ReAllocBlocks( unsigned long newCnt )
    {
// see header file for class documentation
    DataBlock* newBlocks;
    newBlocks = new DataBlock[ newCnt ];
    if ( !newBlocks )
	return ENOMEM;
    unsigned long cpCnt = (newCnt > fMaxBlockCnt) ? fMaxBlockCnt : newCnt;
    if ( fBlocks ) {
    memcpy( newBlocks, fBlocks, cpCnt*sizeof(DataBlock) );
    } else {
      fMaxBlockCnt=0;
    }
    if ( newCnt > fMaxBlockCnt )
	memset( newBlocks+fMaxBlockCnt, 0, (newCnt-fMaxBlockCnt)*sizeof(DataBlock) );
    if ( fBlocks )
	delete [] fBlocks;
    fBlocks = newBlocks;
    fMaxBlockCnt = newCnt;
    return 0;
    }

homer_uint64 AliHLTHOMERReader::GetSourceEventID( const DataSource& source )
    {
// see header file for class documentation
    homer_uint8 sourceByteOrder = ((homer_uint8*)source.fData)[ kByteOrderAttribute_8b_Offset ];
    return Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)source.fData)[ kSubType1_64b_Offset ] );
    }

homer_uint64 AliHLTHOMERReader::GetSourceEventType( const DataSource& source )
    {
// see header file for class documentation
    homer_uint8 sourceByteOrder = ((homer_uint8*)source.fData)[ kByteOrderAttribute_8b_Offset ];
    return Swap( kHOMERNativeByteOrder, sourceByteOrder, ((homer_uint64*)source.fData)[ kType_64b_Offset ] );
    }

homer_uint64 AliHLTHOMERReader::Swap( homer_uint8 destFormat, homer_uint8 sourceFormat, homer_uint64 source ) const
    {
// see header file for class documentation
    if ( destFormat == sourceFormat )
        return source;
    else
        return ((source & 0xFFULL) << 56) | 
    	((source & 0xFF00ULL) << 40) | 
    	((source & 0xFF0000ULL) << 24) | 
    	((source & 0xFF000000ULL) << 8) | 
    	((source & 0xFF00000000ULL) >> 8) | 
    	((source & 0xFF0000000000ULL) >> 24) | 
    	((source & 0xFF000000000000ULL) >>  40) | 
	((source & 0xFF00000000000000ULL) >> 56);
    }

homer_uint32 AliHLTHOMERReader::Swap( homer_uint8 destFormat, homer_uint8 sourceFormat, homer_uint32 source ) const
    {
// see header file for class documentation
    if ( destFormat == sourceFormat )
        return source;
    else
        return ((source & 0xFFUL) << 24) | 
    	((source & 0xFF00UL) << 8) | 
    	((source & 0xFF0000UL) >> 8) | 
    	((source & 0xFF000000UL) >> 24);
    }

AliHLTHOMERReader* AliHLTHOMERReaderCreateFromTCPPort(const char* hostname, unsigned short port )
    {
      // see header file for function documentation
      return new AliHLTHOMERReader(hostname, port);
    }

AliHLTHOMERReader* AliHLTHOMERReaderCreateFromTCPPorts(unsigned int tcpCnt, const char** hostnames, unsigned short* ports)
    {
      // see header file for function documentation
      return new AliHLTHOMERReader(tcpCnt, hostnames, ports);
    }

AliHLTHOMERReader* AliHLTHOMERReaderCreateFromBuffer(const void* pBuffer, int size)
    {
      // see header file for function documentation
      return new AliHLTHOMERReader(pBuffer, size);
    }

void AliHLTHOMERReaderDelete(AliHLTHOMERReader* pInstance)
    {
      // see header file for function documentation
      if (pInstance) delete pInstance;
    }

/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/
