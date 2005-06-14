////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_ADC_INPUT_INTERFACE_HPP
#define dHLT_DDL_ADC_INPUT_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"
#include "EventID.hpp"
#include "ADCStream.hpp"

namespace dHLT
{
namespace DDL
{


class ADCInputCallback
{
public:

	virtual ADCStream* AllocateADCStream(const UInt size) = 0;

	virtual void ReturnADCStream(const EventID event, ADCStream* adcstream) = 0;

	virtual void EndOfADCStreams(const EventID event) = 0;

};


} // DDL
} // dHLT

#endif // dHLT_DDL_CLUSTER_INPUT_INTERFACE_HPP
