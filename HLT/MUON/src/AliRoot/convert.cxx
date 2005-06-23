////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/convert.hpp"
#include "Utils.hpp"
#include "Error.hpp"

namespace dHLT
{
namespace AliRoot
{


AliMUONHLT::Point Convert(const dHLT::Point& point)
{
	DebugMsg(5, "Convert from dHLT::Point");
	return AliMUONHLT::Point(point.x, point.y);
}


dHLT::Point Convert(const AliMUONHLT::Point& point)
{
	DebugMsg(5, "Convert from AliMUONHLT::Point");
	return dHLT::Point(point.fX, point.fY);
}


AliMUONHLT::TriggerRecord Convert(const dHLT::TriggerRecord& record, const Int_t triggernumber)
{
	DebugMsg(5, "Convert from dHLT::TriggerRecord");
	// If the trigger number is negative then set it to zero.
	if (triggernumber >= 0)
	{
		return AliMUONHLT::TriggerRecord(
					triggernumber,
					record.sign,
					record.pt,
					Convert( record.station1impact ),
					Convert( record.station2impact )
				);
	}
	else
	{
		return AliMUONHLT::TriggerRecord(
					0,
					record.sign,
					record.pt,
					Convert( record.station1impact ),
					Convert( record.station2impact )
				);
	}
}


dHLT::TriggerRecord Convert(const AliMUONHLT::TriggerRecord& record)
{
	DebugMsg(5, "Convert from AliMUONHLT::TriggerRecord");
	return dHLT::TriggerRecord(
				(dHLT::ParticleSign) record.ParticleSign(),
				record.Pt(),
				Convert( record.Station1Point() ),
				Convert( record.Station2Point() )
			);
}


AliMUONHLT::Track Convert(const dHLT::Track& track)
{
	DebugMsg(5, "Convert from dHLT::Track");
	AliMUONHLT::Track t;
	t.TriggerID( track.triggerid );
	t.ParticleSign( track.sign );
	t.P( track.p );
	t.Pt( track.pt );
	for (Int_t i = 0; i < 10; i++)
	{
		t.Hit(i) = Convert( track.point[i] );
		
		// Only convert if the ROI is valid. Otherwise the Region object
		// is filled with NaN's.
		if (track.region[i] != dHLT::INVALID_ROI)
			t.RegionOfInterest(i) = Convert( track.region[i] );
	}
	return t;
}


dHLT::Track Convert(const AliMUONHLT::Track& track)
{
	DebugMsg(5, "Convert from AliMUONHLT::Track");
	dHLT::Track t;
	t.triggerid = track.TriggerID();
	t.sign = (dHLT::ParticleSign) track.ParticleSign();
	t.p = track.P();
	t.pt = track.Pt();
	for (Int_t i = 0; i < 10; i++)
	{
		t.point[i] = Convert( track.Hit(i) );
		t.region[i] = Convert( track.RegionOfInterest(i), i );
	}
	return t;
}


AliMUONHLT::Region Convert(const dHLT::ROI region)
{
	DebugMsg(5, "Convert from dHLT::ROI");
	dHLT::RegionOfInterest roi(region);
	return AliMUONHLT::Region( roi.Left(), roi.Right(), roi.Bottom(), roi.Top() );
}


dHLT::ROI Convert(const AliMUONHLT::Region& region, const UInt_t chamber)
{
	DebugMsg(5, "Convert from AliMUONHLT::Region");
	// If the chamber number is too big then truncate it.
	if (chamber < 10)
	{
		dHLT::RegionOfInterest roi(
				region.Left(), region.Right(), region.Bottom(), region.Top(),
				(dHLT::ChamberID) chamber
			);
		return roi;
	}
	else
	{
		dHLT::RegionOfInterest roi(
				region.Left(), region.Right(), region.Bottom(), region.Top(),
				dHLT::Chamber10
			);
		return roi;
	}
}


} // AliRoot
} // dHLT
