////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "DDL/print.hpp"

#include <stdio.h>

namespace dHLT
{
namespace DDL
{


void PrintBits(std::ostream& os, UInt x, int bitcount)
{
	for (int i = bitcount-1; i >= 0; i--)
	{
		if ( x & (1 << i) )
			os << "1";
		else
			os << "0";
	};
	os << "b";
};


void PrintHex(std::ostream& os, UInt x, UChar width)
{
	char buf[16];
	char* str = (char*)(&buf[0]);
	sprintf(str, "0x%%.%dX", width);

	char buf2[16];
	char* str2 = (char*)(&buf2[0]);
	sprintf(str2, str, x);
	os << str2;
};


std::ostream& operator << (std::ostream& os, const DDLHeader& h)
{
	os << "DDL header:" << endl;
	os << "                 block length : " << (UInt)h.blocklength << endl;
	os << "            L2 bunch crossing : " << h.L2bunchcross << endl;
	os << "  reserved word 1, bits 12-15 : "; PrintBits(os, h.reserved1, 4); os << endl;
	os << "              L1 trigger type : " << h.L1triggertype << endl;
	os << "               format version : " << h.formatversion << endl;
	os << "              L2 orbit number : " << h.L2orbit << endl;
	os << "  reserved word 2, bits 24-31 : "; PrintBits(os, h.reserved2, 8); os << endl;
	os << "  participating sub-detectors : "; PrintBits(os, h.subdetectors, 24); os << endl;
	os << "             block attributes : "; PrintBits(os, h.blockattribs, 8); os << endl;
	os << "           TTC bunch crossing : " << h.TTCbunchcross << endl;
	os << "                 Status/Error : "; PrintBits(os, h.statuserror, 16); os << endl;
	os << "  reserved word 4, bits 28-31 : "; PrintBits(os, h.reserved3, 4); os << endl;
	os << "            trigger class low : " << (UInt)h.triggerclasslow << endl;
	os << "           trigger class high : " << (UInt)h.triggerclasshigh << endl;
	os << "  reserved word 6, bits 18-27 : "; PrintBits(os, h.reserved4, 10); os << endl;
	os << "       region of interest low : " << h.roilow << endl;
	os << "      region of interest high : " << (UInt)h.roihigh << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const BlockHeader& h)
{
	os << "Block header:" << endl;
	os << "            block length : " << h.blocklength << endl;
	os << "             data length : " << h.datalength << endl;
	os << "                  DSP ID : " << h.dspid << endl;
	os << "  trigger counter word 0 : " << h.triggercounter[0] << endl;
	os << "  trigger counter word 1 : " << h.triggercounter[1] << endl;
	os << "  trigger counter word 2 : " << h.triggercounter[2] << endl;
	os << "  trigger counter word 3 : " << h.triggercounter[3] << endl;
	os << "                 padding : "; PrintHex(os, h.padding, 8); os << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const DSPHeader& h)
{
	os << "DSP header:" << endl;
	os << "            block length : " << h.blocklength << endl;
	os << "             data length : " << h.datalength << endl;
	os << "  trigger counter word 0 : " << h.triggercounter[0] << endl;
	os << "  trigger counter word 1 : " << h.triggercounter[1] << endl;
	os << "  trigger counter word 2 : " << h.triggercounter[2] << endl;
	os << "  trigger counter word 3 : " << h.triggercounter[3] << endl;
	os << "                  DSP ID : " << h.dspid << endl;
	os << "      event padding flag : "; PrintHex(os, h.eventflag, 8); os << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const PatchBusHeader& h)
{
	os << "Patch bus header:" << endl;
	os << "     block length : " << h.blocklength << endl;
	os << "      data length : " << h.datalength << endl;
	os << "     bus patch ID : " << h.patchbusid << endl;
	os << "  trigger counter : " << h.counter << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const ADCData& d)
{
	os << "[" << d.parity 
	<< "][" << d.control1
	<< "][" << d.control0
	<< "][" << d.manuid
	<< "][" << d.channel
	<< "][" << d.signal
	<< "]";
	return os;
};


std::ostream& operator << (std::ostream& os, const LocalData& d)
{
	os << "\tx1 = "; PrintBits(os, d.x1, 16); os << endl;
	os << "\tx2 = "; PrintBits(os, d.x2, 16); os << endl;
	os << "\tx3 = "; PrintBits(os, d.x3, 16); os << endl;
	os << "\tx4 = "; PrintBits(os, d.x4, 16); os << endl;
	os << "\ty1 = "; PrintBits(os, d.y1, 16); os << endl;
	os << "\ty2 = "; PrintBits(os, d.y2, 16); os << endl;
	os << "\ty3 = "; PrintBits(os, d.y3, 16); os << endl;
	os << "\ty4 = "; PrintBits(os, d.y4, 16); os << endl;

	os << "\t   x position = " << d.xpos << endl;
	os << "\t  x deviation = " << d.xdev << endl;
	os << "\t   y position = " << d.ypos << endl;
	os << "\t    y trigger : " << d.ytrigger << endl;
	os << "\t     decision : "; PrintBits(os, d.decision, 4); os << endl;
	os << "\tboard address : " << d.address << endl;
	os << "\t     reserved : "; PrintBits(os, d.reserved, 9); os << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const RegionalData& d)
{
	os << "          reserved : "; PrintBits(os, d.reserved, 14); os << endl;
	os << "  software version : " << d.version << endl;
	os << "       regional ID : " << d.id << endl;
	os << "     serial number : " << d.serialnumber << endl;

	for (UInt i = 0; i < 16; i++)
	{
		os << "=============== Local Card[" << i << "] ===============" << endl;
		os << d.localcard[i] << endl;
	};

	os << "  regional input[0] : " << d.input0 << endl;
	os << "  regional input[1] : " << d.input1 << endl;
	os << "  regional input[2] : " << d.input2 << endl;
	os << "  regional input[3] : " << d.input3 << endl;
	os << "  regional input[4] : " << d.input4 << endl;
	os << "  regional input[5] : " << d.input5 << endl;
	os << "  regional input[6] : " << d.input6 << endl;
	os << "  regional input[7] : " << d.input7 << endl;
	os << "  regional input[8] : " << d.input8 << endl;
	os << "  regional input[9] : " << d.input9 << endl;
	os << "  regional input[A] : " << d.inputA << endl;
	os << "  regional input[B] : " << d.inputB << endl;
	os << "  regional input[C] : " << d.inputC << endl;
	os << "  regional input[D] : " << d.inputD << endl;
	os << "  regional input[E] : " << d.inputE << endl;
	os << "  regional input[F] : " << d.inputF << endl;
	os << "    regional output : "; PrintHex(os, d.output, 8); os << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const EnhancedHeader& h)
{
	os << "          reserved : "; PrintBits(os, h.reserved, 14); os << endl;
	os << "        event type : " << h.eventtype << endl;
	os << "     serial number : " << h.serialnumber << endl;
	os << "  software version : " << h.version << endl;
	os << "           DARC ID : " << h.darcid << endl;

	for (UInt i = 0; i < 16; i++)
		os << "  global input [" << i << "] : " << (UInt)h.globalinput[i] << endl;
	os << "  global output : "; PrintHex(os, h.globaloutput, 8); os << endl;
	return os;
};


std::ostream& operator << (std::ostream& os, const TriggerData& d)
{
	os << d.header << endl;
	for (UInt i = 0; i < 8; i++)
	{
		os << "=============== Regional Card[" << i << "] ===============" << endl;
		os << d.regionalcard[i] << endl;
	};
	os << "End of data marker:" << endl << "\t";
	PrintHex(os, d.endmarker[0], 8);
	os << " ";
	PrintHex(os, d.endmarker[0], 8);
	return os;
};


}; // DDL
}; // dHLT
