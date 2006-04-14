////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/convert.hpp"
#include "Utils.hpp"
#include "Error.hpp"


AliHLTMUONPoint AliHLTMUONConvert(const AliHLTMUONCorePoint& point)
{
	DebugMsg(5, "Convert from AliHLTMUONCorePoint");
	return AliHLTMUONPoint(point.fX, point.fY);
}


AliHLTMUONCorePoint AliHLTMUONConvert(const AliHLTMUONPoint& point)
{
	DebugMsg(5, "Convert from AliHLTMUONPoint");
	return AliHLTMUONCorePoint(point.fX, point.fY);
}


AliHLTMUONTriggerRecord AliHLTMUONConvert(const AliHLTMUONCoreTriggerRecord& record, Int_t triggernumber)
{
	DebugMsg(5, "Convert from dHLT::TriggerRecord");
	// If the trigger number is negative then set it to zero.
	if (triggernumber >= 0)
	{
		return AliHLTMUONTriggerRecord(
					triggernumber,
					record.fSign,
					record.fPt,
					AliHLTMUONConvert( record.fStation1impact ),
					AliHLTMUONConvert( record.fStation2impact )
				);
	}
	else
	{
		return AliHLTMUONTriggerRecord(
					0,
					record.fSign,
					record.fPt,
					AliHLTMUONConvert( record.fStation1impact ),
					AliHLTMUONConvert( record.fStation2impact )
				);
	}
}


AliHLTMUONCoreTriggerRecord AliHLTMUONConvert(const AliHLTMUONTriggerRecord& record)
{
	DebugMsg(5, "Convert from AliHLTMUONTriggerRecord");
	return AliHLTMUONCoreTriggerRecord(
				(AliHLTMUONCoreParticleSign) record.ParticleSign(),
				record.Pt(),
				AliHLTMUONConvert( record.Station1Point() ),
				AliHLTMUONConvert( record.Station2Point() )
			);
}


AliHLTMUONTrack AliHLTMUONConvert(const AliHLTMUONCoreTrack& track)
{
	DebugMsg(5, "Convert from dHLT::Track");
	AliHLTMUONTrack t;
	t.TriggerID( track.fTriggerid );
	t.ParticleSign( track.fSign );
	t.P( track.fP );
	t.Pt( track.fPt );
	for (Int_t i = 0; i < 10; i++)
	{
		t.Hit(i) = AliHLTMUONConvert( track.fPoint[i] );
		
		// Only convert if the ROI is valid. Otherwise the Region object
		// is filled with NaN's.
		if (track.fRegion[i] != kInvalidROI)
			t.RegionOfInterest(i) = AliHLTMUONConvert( track.fRegion[i] );
	}
	return t;
}


AliHLTMUONCoreTrack AliHLTMUONConvert(const AliHLTMUONTrack& track)
{
	DebugMsg(5, "Convert from AliHLTMUONTrack");
	AliHLTMUONCoreTrack t;
	t.fTriggerid = track.TriggerID();
	t.fSign = (AliHLTMUONCoreParticleSign) track.ParticleSign();
	t.fP = track.P();
	t.fPt = track.Pt();
	for (Int_t i = 0; i < 10; i++)
	{
		t.fPoint[i] = AliHLTMUONConvert( track.Hit(i) );
		t.fRegion[i] = AliHLTMUONConvert( track.RegionOfInterest(i), i );
	}
	return t;
}


AliHLTMUONRegion AliHLTMUONConvert(const AliHLTMUONCoreROI region)
{
	DebugMsg(5, "Convert from AliHLTMUONCoreROI");
	AliHLTMUONCoreRegionOfInterest roi(region);
	return AliHLTMUONRegion( roi.Left(), roi.Right(), roi.Bottom(), roi.Top() );
}


AliHLTMUONCoreROI AliHLTMUONConvert(const AliHLTMUONRegion& region, UInt_t chamber)
{
	DebugMsg(5, "Convert from AliHLTMUONRegion");
	// If the chamber number is too big then truncate it.
	if (chamber < 10)
	{
		AliHLTMUONCoreRegionOfInterest roi(
				region.Left(), region.Right(), region.Bottom(), region.Top(),
				(AliHLTMUONCoreChamberID) chamber
			);
		return roi;
	}
	else
	{
		AliHLTMUONCoreRegionOfInterest roi(
				region.Left(), region.Right(), region.Bottom(), region.Top(),
				kChamber10
			);
		return roi;
	}
}

