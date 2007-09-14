// XEMacs -*-C++-*-
#ifndef _HOMER_H_
#define _HOMER_H_
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTHomerReader.h
    @author Timm Steinbeck
    @date   Sep 14 2007
    @brief  HLT Online Monitoring Environment including ROOT - Reader
    @note   migrated from PubSub HLT-stable-20070905.141318 (rev 2375)    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <limits.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include "AliHLTHOMERData.h"




class MonitoringReader
    {
    public:

	MonitoringReader() {};
	virtual ~MonitoringReader() {};
	
	/* Read in the next available event */
	virtual int ReadNextEvent() = 0;
	/* Read in the next available event, wait max. timeout microsecs. */
	virtual int ReadNextEvent( unsigned long timeout ) = 0;
	
	/* Return the type of the current event */
	virtual homer_uint64 GetEventType() const = 0;

	/* Return the ID of the current event */
	virtual homer_uint64 GetEventID() const = 0;
	
	/* Return the number of data blocks in the current event */
	virtual unsigned long GetBlockCnt() const = 0;
	
	/* Return the size (in bytes) of the current event's data
	   block with the given block index (starting at 0). */
	virtual unsigned long GetBlockDataLength( unsigned long ndx ) const = 0;
	/* Return a pointer to the start of the current event's data
	   block with the given block index (starting at 0). */
	virtual const void* GetBlockData( unsigned long ndx ) const = 0;
	/* Return IP address or hostname of node which sent the 
	   current event's data block with the given block index 
	   (starting at 0). */
	virtual const char* GetBlockSendNodeID( unsigned long ndx ) const = 0;
	/* Return byte order of the data stored in the 
	   current event's data block with the given block 
	   index (starting at 0). 
	   0 is unknown alignment, 
	   1 ist little endian, 
	   2 is big endian. */
	virtual homer_uint8 GetBlockByteOrder( unsigned long ndx ) const = 0;
	/* Return the alignment (in bytes) of the given datatype 
	   in the data stored in the current event's data block
	   with the given block index (starting at 0). 
	   Possible values for the data type are
	   0: homer_uint64
	   1: homer_uint32
	   2: uin16
	   3: homer_uint8
	   4: double
	   5: float
	*/
	virtual homer_uint8 GetBlockTypeAlignment( unsigned long ndx, homer_uint8 dataType ) const = 0;

	virtual homer_uint64 GetBlockStatusFlags( unsigned long ndx ) const = 0;

#ifdef USE_ROOT
        ClassDef(MonitoringReader,1);
#endif
    };



class HOMERReader: public MonitoringReader
    {
    public:
#ifdef USE_ROOT
	HOMERReader();
#endif

	/* Constructors & destructors, HOMER specific */
	/* For reading from a TCP port */
	HOMERReader( const char* hostname, unsigned short port );
	/* For reading from multiple TCP ports */
	HOMERReader( unsigned int tcpCnt, const char** hostnames, unsigned short* ports );
	/* For reading from a System V shared memory segment */
	HOMERReader( key_t shmKey, int shmSize );
	/* For reading from multiple System V shared memory segments */
	HOMERReader( unsigned int shmCnt, key_t* shmKey, int* shmSize );
	/* For reading from multiple TCP ports and multiple System V shared memory segments */
	HOMERReader( unsigned int tcpCnt, const char** hostnames, unsigned short* ports, 
		     unsigned int shmCnt, key_t* shmKey, int* shmSize );
	virtual ~HOMERReader();

	/* Return the status of the connection as established by one of the constructors.
	   0 means connection is ok, non-zero specifies the type of error that occured. */
	int GetConnectionStatus() const
		{
		return fConnectionStatus;
		}

	/* Return the index of the connection for which an error given by the above
	   function occured. */
	unsigned int GetErrorConnectionNdx() const
		{
		return fErrorConnection;
		}

	void SetEventRequestAdvanceTime( unsigned long time_us )
		{
		fEventRequestAdvanceTime_us = time_us;
		}

	/* Defined in MonitoringReader */
	/* Read in the next available event */
	virtual int  ReadNextEvent();
	/* Read in the next available event */
	virtual int ReadNextEvent( unsigned long timeout );

	/* Return the type of the current event */
	virtual homer_uint64 GetEventType() const
		{
		return fCurrentEventType;
		}

	/* Return the ID of the current event */
	virtual homer_uint64 GetEventID() const
		{
		return fCurrentEventID;
		}

	/* Return the number of data blocks in the current event */
	virtual unsigned long GetBlockCnt() const
		{
		return fBlockCnt;
		}

	/* Return the size (in bytes) of the current event's data
	   block with the given block index (starting at 0). */
	virtual const void* GetBlockData( unsigned long ndx ) const;
	/* Return a pointer to the start of the current event's data
	   block with the given block index (starting at 0). */
	virtual unsigned long GetBlockDataLength( unsigned long ndx ) const;
	/* Return IP address or hostname of node which sent the 
	   current event's data block with the given block index 
	   (starting at 0).
	   For HOMER this is the ID of the node on which the subscriber 
	   that provided this data runs/ran. */
	virtual const char* GetBlockSendNodeID( unsigned long ndx ) const;
	/* Return byte order of the data stored in the 
	   current event's data block with the given block 
	   index (starting at 0). 
	   0 is unknown alignment, 
	   1 ist little endian, 
	   2 is big endian. */
	virtual homer_uint8 GetBlockByteOrder( unsigned long ndx ) const;
	/* Return the alignment (in bytes) of the given datatype 
	   in the data stored in the current event's data block
	   with the given block index (starting at 0). 
	   Possible values for the data type are
	   0: homer_uint64
	   1: homer_uint32
	   2: uin16
	   3: homer_uint8
	   4: double
	   5: float
	*/
	virtual homer_uint8 GetBlockTypeAlignment( unsigned long ndx, homer_uint8 dataType ) const;

	virtual homer_uint64 GetBlockStatusFlags( unsigned long ndx ) const;

	/* HOMER specific */
	/* Return the type of the data in the current event's data
	   block with the given block index (starting at 0). */
	homer_uint64 GetBlockDataType( unsigned long ndx ) const;
	/* Return the origin of the data in the current event's data
	   block with the given block index (starting at 0). */
	homer_uint32 GetBlockDataOrigin( unsigned long ndx ) const;
	/* Return a specification of the data in the current event's data
	   block with the given block index (starting at 0). */
	homer_uint32 GetBlockDataSpec( unsigned long ndx ) const;

	/* Find the next data block in the current event with the given
	   data type, origin, and specification. Returns the block's 
	   index. */
	unsigned long FindBlockNdx( homer_uint64 type, homer_uint32 origin, 
				    homer_uint32 spec, unsigned long startNdx=0 ) const;

	/* Find the next data block in the current event with the given
	   data type, origin, and specification. Returns the block's 
	   index. */
	unsigned long FindBlockNdx( char type[8], char origin[4], 
				    homer_uint32 spec, unsigned long startNdx=0 ) const;
	
	/* Return the ID of the node that actually produced this data block.
	   This may be different from the node which sent the data to this
	   monitoring object as returned by GetBlockSendNodeID. */
	const char* GetBlockCreateNodeID( unsigned long ndx ) const;

    protected:

	enum DataSourceType { kUndef=0, kTCP, kShm };
	struct DataSource
	    {
		DataSource() { fType = kUndef; };
		DataSourceType fType;
		unsigned fNdx; // This source's index
		const char* fHostname; // Filled for both Shm and TCP
		unsigned short fTCPPort;
		key_t fShmKey;
		int fShmSize;
		int fTCPConnection; // File descriptor for the TCP connection
		int fShmID; // ID of the shared memory area
		void* fShmPtr; // Pointer to shared memory area
		void* fData; // Pointer to data read in for current event from this source
		unsigned long fDataSize; // Size of data (to be) read in for current event from this source
		unsigned long fDataRead; // Data actually read for current event
	    };

	void Init();
	
	bool AllocDataSources( unsigned int sourceCnt );
	int AddDataSource( const char* hostname, unsigned short port, DataSource& source );
	int AddDataSource( key_t shmKey, int shmSize, DataSource& source );
	void FreeDataSources();
	int FreeShmDataSource( DataSource& source );
	int FreeTCPDataSource( DataSource& source );
	int ReadNextEvent( bool useTimeout, unsigned long timeout );
	void ReleaseCurrentEvent();
	int TriggerTCPSource( DataSource& source, bool useTimeout, unsigned long timeout );
	int TriggerShmSource( DataSource& source, bool useTimeout, unsigned long timeout );
	int ReadDataFromTCPSources( unsigned sourceCnt, DataSource* sources, bool useTimeout, unsigned long timeout );
	int ReadDataFromShmSources( unsigned sourceCnt, DataSource* sources, bool useTimeout, unsigned long timeout );
	int ParseSourceData( DataSource& source );
	int ReAllocBlocks( unsigned long newCnt );
	homer_uint64 GetSourceEventID( DataSource& source );
	homer_uint64 GetSourceEventType( DataSource& source );
	homer_uint64 Swap( homer_uint8 destFormat, homer_uint8 sourceFormat, homer_uint64 source )
		{
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
	homer_uint32 Swap( homer_uint8 destFormat, homer_uint8 sourceFormat, homer_uint32 source )
		{
		if ( destFormat == sourceFormat )
		    return source;
		else
		    return ((source & 0xFFUL) << 24) | 
			((source & 0xFF00UL) << 8) | 
			((source & 0xFF0000UL) >> 8) | 
			((source & 0xFF000000UL) >> 24);
		}


	struct DataBlock
	    {
		unsigned int fSource; // Index of originating data source
		void* fData;
		unsigned long fLength;
		homer_uint64* fMetaData; // Pointer to meta data describing data itself.
		const char* fOriginatingNodeID;
	    };

	homer_uint64 fCurrentEventType;
	homer_uint64 fCurrentEventID;
	unsigned long fBlockCnt;
	unsigned long fMaxBlockCnt;
	DataBlock* fBlocks;
	
	unsigned int fDataSourceCnt;
	unsigned int fTCPDataSourceCnt;
	unsigned int fShmDataSourceCnt;
	unsigned int fDataSourceMaxCnt;
	DataSource* fDataSources;

	
	int fConnectionStatus;
	unsigned fErrorConnection;

	unsigned long fEventRequestAdvanceTime_us;

#ifdef USE_ROOT
        ClassDef(HOMERReader,2);
#endif
    };

#endif /* _HOMER_H_ */
