////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_TRACKING_RAW_DATA_MAP_HPP
#define dHLT_DDL_TRACKING_RAW_DATA_MAP_HPP

#include "Error.hpp"
#include "DDL/RawDataFormat.hpp"

namespace dHLT
{
namespace DDL
{


enum
{
	PARSE_ERROR = 0x10030001,
	BLOCK_LENGTH_TOO_SMALL = 0x10030002,
	BAD_BLOCK_LENGTH_IN_DDL_HEADER = 0x10030004
	BAD_BLOCK_LENGTH_IN_BLOCK_HEADER  = 0x10030004
	BLOCK_LENGTH_TOO_SMALL = 0x10030002,
	HEADER_LENGTH_TOO_SMALL = 0x10030002,
};


class ParseError : public Error
{
public:
	virtual const char* Message() const throw()
	{
	};
	
	virtual Int ErrorCode() const throw()
	{
	};
};


class TrackingRawDataMap
{
public:

	TrackingRawDataMap()
	{
		// Initialise all internal pointers to NULL.
		ddlheader = NULL;
		for (UInt i = 0; i < 2; i++)
		{
			blockheader[i] = NULL;
			for (UInt j = 0; j < 5; j++)
			{
				dspheader[i][j] = NULL;
				for (UInt k = 0; k < 5; k++)
				{
					patchbusheader[i][j][k] = NULL;
					data[i][j][k] = NULL;
				};
			};
		};
	};


	/* Fills the map with pointers into the DDL raw data buffer.
	   'stream' must be the pointer to the start of the DDL data buffer.
	   'size' is the total size of the buffer in bytes.
	 */
	bool Create(const void* stream, const UInt size)
	{
		ddlheader = (DDLHeader*) stream;
		if ( ddlheader->blocklength * sizeof(UInt) != size )
		{
			cerr << "Error: DDL data size indicated by DDL header does not match stream size." << endl;
			return false;
		};

		// Get a pointer marking the end of the DDL data block, and 
		// a pointer 'current' to the start of the data payload.
		UInt* current = (UInt*)(ddlheader + 1);
		UInt* end = (UInt*)ddlheader + ddlheader->blocklength;

		if (current > end)
		{
			cerr << "Error: ivalid block length in DDL header." << endl;
			return false;
		};

		// Locate all the Block headers by incrementing the current pointer
		// by the block length (in 32-bit words) each time.
		for (UInt i = 0; i < 2; i++)
		{
			blockheader[i] = (BlockHeader*)current;
			current += blockheader[i]->blocklength;
			if (current > end)
			{
				cerr << "Error: ivalid block length in block header " << i << endl;
				return false;
			};
			if ( ! LocateDSPHeaders(i, current) ) return false;
		};

		return true;
	};


	DDLHeader* GetDDLHeader() const
	{
		return ddlheader;
	};


	BlockHeader* GetBlockHeader(const UInt block) const
	{
		Assert( block < 2 );
		return blockheader[block];
	};


	DSPHeader* GetDSPHeader(const UInt block, const UInt dsp) const
	{
		Assert( block < 2 );
		Assert( dsp < 5 );
		return dspheader[block][dsp];
	};


	PatchBusHeader* GetPatchBusHeader(const UInt block, const UInt dsp, const UInt patchbus) const
	{
		Assert( block < 2 );
		Assert( dsp < 5 );
		Assert( patchbus < 5 );
		return patchbusheader[block][dsp][patchbus];
	};


	ADCData* GetADCData(const UInt block, const UInt dsp, const UInt patchbus) const
	{
		Assert( block < 2 );
		Assert( dsp < 5 );
		Assert( patchbus < 5 );
		return data[block][dsp][patchbus];
	};


private:

	bool LocateDSPHeaders(const UInt blocknumber, const UInt* end)
	{
		BlockHeader* block = blockheader[blocknumber];

		if (block->blocklength < block->datalength)
		{
			cerr << "Error: Block[" << blocknumber << "] block length is smaller than the payload length." << endl;
			return false;
		};

		UInt headerlength = (block->blocklength - block->datalength) * sizeof(UInt);
		if (headerlength < sizeof(BlockHeader))
		{
			cerr << "Error: Block[" << blocknumber << "] header length is too small." << endl;
			return false;
		};

		// Get a pointer marking the start of the Block payload.
		UInt* current = (UInt*)block + block->blocklength - block->datalength;
		if (current > end)
		{
			cerr << "Error: ivalid block length in Block[" << blocknumber << "] header." << endl;
			return false;
		};

		// Locate all the DSP headers by incrementing the current pointer
		// by the block length (in 32-bit words) each time.
		for (UInt i = 0; i < 5; i++)
		{
			dspheader[blocknumber][i] = (DSPHeader*)current;
			current += dspheader[blocknumber][i]->blocklength;
			if (current > end)
			{
				cerr << "Error: ivalid block length in DSP header " << i << endl;
				return false;
			};
			if ( ! LocatePatchBusHeaders(blocknumber, i, current) ) return false;
		};

		return true;
	};


	bool LocatePatchBusHeaders(const UInt blocknumber, const UInt dspnumber, const UInt* end)
	{
		DSPHeader* dsp = dspheader[blocknumber][dspnumber];

		if (dsp->blocklength < dsp->datalength)
		{
			cerr << "Error: DSP[" << blocknumber << "][" << dspnumber
			     << "] block length is smaller than the payload length." << endl;
			return false;
		};

		UInt headerlength = (dsp->blocklength - dsp->datalength) * sizeof(UInt);
		if (headerlength < sizeof(DSPHeader))
		{
			cerr << "Error: DSP[" << blocknumber << "][" << dspnumber
			     << "] header length is too small." << endl;
			return false;
		};

		// Get a pointer marking the start of the DSP payload.
		UInt* current = (UInt*)dsp + dsp->blocklength - dsp->datalength;
		if (current > end)
		{
			cerr << "Error: ivalid block length in DSP[" << blocknumber << "][" << dspnumber
			     << "] header." << endl;
			return false;
		};

		// Locate all the Bus Patch headers by incrementing the current pointer
		// by the block length (in 32-bit words) each time.
		for (UInt i = 0; i < 5; i++)
		{
			patchbusheader[blocknumber][dspnumber][i] = (PatchBusHeader*)current;
			current += patchbusheader[blocknumber][dspnumber][i]->blocklength;
			if (current > end)
			{
				cerr << "Error: ivalid block length in bus patch header " << i << endl;
				return false;
			}
			else if (current == end)
			{
				// Reached the end of the DSP block.
				break;
			};
			
			if (patchbusheader[blocknumber][dspnumber][i]->datalength > 0)
			{
				if (! LocateData(blocknumber, dspnumber, i, current) )
					return false;
			};
		};

		return true;
	};
	

	bool LocateData(const UInt blocknumber, const UInt dspnumber, const UInt patchbusnumber, const UInt* end)
	{
		PatchBusHeader* patchbus = patchbusheader[blocknumber][dspnumber][patchbusnumber];

		if (patchbus->blocklength < patchbus->datalength)
		{
			cerr << "Error: PatchBus[" << blocknumber << "][" << dspnumber << "][" << patchbusnumber
			     << "] block length is smaller than the payload length." << endl;
			return false;
		};

		UInt headerlength = (patchbus->blocklength - patchbus->datalength) * sizeof(UInt);
		if (headerlength < sizeof(PatchBusHeader))
		{
			cerr << "Error: PatchBus[" << blocknumber << "][" << dspnumber << "][" << patchbusnumber
			     << "] header length is too small." << endl;
			return false;
		};

		// Get a pointer marking the start of the PatchBus payload.
		UInt* current = (UInt*)patchbus + patchbus->blocklength - patchbus->datalength;
		if (current > end)
		{
			cerr << "Error: ivalid block length in PatchBus[" << blocknumber << "]["
			     << dspnumber << "][" << patchbusnumber << "] header." << endl;
			return false;
		};

		data[blocknumber][dspnumber][patchbusnumber] = (ADCData*)current;
		return true;
	};


	DDLHeader* ddlheader;
	BlockHeader* blockheader[2];
	DSPHeader* dspheader[2][5];
	PatchBusHeader* patchbusheader[2][5][5];
	ADCData* data[2][5][5];
};


}; // DDL
}; // dHLT

#endif // dHLT_DDL_TRACKING_RAW_DATA_MAP_HPP

