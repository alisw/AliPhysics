////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_CONTROL_CONTROL_INTERFACE_HPP
#define dHLT_CONTROL_CONTROL_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "Error.hpp"
#include "EventID.hpp"
#include "TriggerRecord.hpp"
#include "RegionOfInterest.hpp"
#include "Cluster.hpp"
#include "Track.hpp"

namespace dHLT
{
namespace Control
{


enum StatusType
{
	UnknownStatus   = 0,
	Running         = 1,
	Suspended       = 2,
	Fault           = -1
};


class ControlInterface
{
public:

	virtual void Start() = 0;
	virtual void Stop() = 0;
	virtual StatusType Status() = 0;
	
	virtual void Reset() = 0;

};


class ControlCallback
{
public:

	virtual void ReportError(const Error& error) = 0;

};


#undef ABSTRACT_METHOD

} // Control
} // dHLT

#endif // dHLT_CONTROL_CONTROL_INTERFACE_HPP
