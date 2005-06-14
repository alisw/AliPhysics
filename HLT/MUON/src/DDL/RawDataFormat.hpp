////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_RAW_DATA_FORMAT_HPP
#define dHLT_DDL_RAW_DATA_FORMAT_HPP

#include "BasicTypes.hpp"

namespace dHLT
{
namespace DDL
{

////////////////////////////////////////////////////////////////////////////////
// Common DDL header structure.

/* Refer to:
   ALICE-INT-2002-010 v 2.0
 */
 
struct DDLHeader
{
	// Word 0:
	UInt blocklength      : 32;
	// Word 1:
	UInt L2bunchcross     : 12;
	UInt reserved1        :  4;  // Must be zero.
	UInt L1triggertype    :  8;
	UInt formatversion    :  8;
	// Word 2:
	UInt L2orbit          : 24;
	UInt reserved2        :  8;  // Must be zero.
	// Word 3:
	UInt subdetectors     : 24;
	UInt blockattribs     :  8;
	// Word 4:
	UInt TTCbunchcross    : 12;
	UInt statuserror      : 16;
	UInt reserved3        :  4;  // Must be zero.
	// Word 5:
	UInt triggerclasslow  : 32;
	// Word 6:
	UInt triggerclasshigh : 18;
	UInt reserved4        : 10;  // Must be zero.
	UInt roilow           :  4;
	// Word 7:
	UInt roihigh          : 32;
};


////////////////////////////////////////////////////////////////////////////////
// Data structures for Tracking chamber DDL streams.

/* Refer to:
   http://aliweb.cern.ch/people/tkuhr/Rawdata.html
 */

struct BlockHeader
{
	UInt blocklength;
	UInt datalength;
	UInt dspid;
	UInt triggercounter[4];
	UInt padding;    // Must be zero.
};


struct DSPHeader
{
	UInt blocklength;
	UInt datalength;
	UInt triggercounter[4];
	UInt dspid;
	UInt eventflag;   // 1 for odd 0 for even number of 32-bit words.
};


struct PatchBusHeader
{
	UInt blocklength;
	UInt datalength;
	UInt patchbusid;
	UInt counter;
};


struct ADCData
{
	struct
	{
		UInt signal    : 12;
		UInt channel   :  6;
		UInt manuid    : 11;
		UInt control0  :  1;
		UInt control1  :  1;
		UInt parity    :  1;
	};
};


////////////////////////////////////////////////////////////////////////////////
// Data structures for L0 trigger DDL streams.

/* Refer to:
   http://aliweb.cern.ch/people/tkuhr/Rawdata.html
 */
 
struct LocalData
{
	// Word 0
	UInt x1 : 16;
	UInt x2 : 16;
	// Word 1
	UInt x3 : 16;
	UInt x4 : 16;
	// Word 2
	UInt y1 : 16;
	UInt y2 : 16;
	// Word 3
	UInt y3 : 16;
	UInt y4 : 16;
	// Word 4
	UInt xpos     : 5;
	UInt xdev     : 5;
	UInt ypos     : 4;
	UInt ytrigger : 1;
	UInt decision : 4;
	UInt address  : 4;  // Board address
	UInt reserved : 9;  // Must be zero.
};


struct RegionalData
{
	// Regional header word.
	struct
	{
		UInt reserved     : 14;  // Must be zero.
		UInt version      : 8;
		UInt id           : 5;  // regional ID.
		UInt serialnumber : 5;
	};
	
	LocalData localcard[16];  // 16 * 5 words of local card data.

	// 2, 32-bit words of regional input data.
	struct
	{
		// Word 0
		UInt input0 : 4;
		UInt input1 : 4;
		UInt input2 : 4;
		UInt input3 : 4;
		UInt input4 : 4;
		UInt input5 : 4;
		UInt input6 : 4;
		UInt input7 : 4;
		// Word 1
		UInt input8 : 4;
		UInt input9 : 4;
		UInt inputA : 4;
		UInt inputB : 4;
		UInt inputC : 4;
		UInt inputD : 4;
		UInt inputE : 4;
		UInt inputF : 4;
	};
	UInt output;  // regional card output.
};


struct EnhancedHeader
{
	struct
	{
		UInt reserved     : 14;  // Must be zero.
		UInt eventtype    :  4;
		UInt serialnumber :  4;
		UInt version      :  8;
		UInt darcid       :  2;
	};
	UChar globalinput[16];
	UInt globaloutput;
};


struct TriggerData
{
	EnhancedHeader header;
	RegionalData regionalcard[8];
	UInt endmarker[2];  // 2 words indicating end of data. Each word must = 0xDEADFACE
};


}; // DDL
}; // dHLT

#endif // dHLT_DDL_RAW_DATA_FORMAT_HPP
