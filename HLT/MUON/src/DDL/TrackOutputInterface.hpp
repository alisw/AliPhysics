////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_DDL_TRACK_OUTPUT_INTERFACE_HPP
#define dHLT_DDL_TRACK_OUTPUT_INTERFACE_HPP

#include "BasicTypes.hpp"
#include "Utils.hpp"

namespace dHLT
{
namespace DDL
{


class TrackOutputInterface
{
public:

	virtual void ReturnTracks(const EventID event, const Track* newtracks, const UInt count) = 0;

	virtual void EndOfTracks(const EventID event) = 0;

};


class TrackOutputCallback
{
public:
	virtual void ReleaseTracks(const Track* tracks) = 0;
};


} // DDL
} // dHLT

#endif // dHLT_DDL_TRACK_OUTPUT_INTERFACE_HPP
