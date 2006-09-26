////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Tracking/MansoTracker.hpp"
#include <math.h>
#include "Tracking/Calculations.hpp"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONError.h"
#include "AliHLTMUONCoreRegionOfInterest.h"


#if defined(DEBUG) || (defined(USE_ALILOG) && ! defined(LOG_NO_DEBUG))
#include <ostream>
#include "Debug/AliHLTMUONPrint.h"
namespace
{

std::ostream& operator << (std::ostream& os, AliHLTMUONCoreMansoTracker::StatesSM4 state)
{
	switch (state)
	{
	case AliHLTMUONCoreMansoTracker::kSM4Idle:          os << "kSM4Idle"; break;
	case AliHLTMUONCoreMansoTracker::kWaitChamber8:     os << "kWaitChamber8"; break;
	case AliHLTMUONCoreMansoTracker::kWaitMoreChamber8: os << "kWaitMoreChamber8"; break;
	case AliHLTMUONCoreMansoTracker::kWaitChamber7:     os << "kWaitChamber7"; break;
	case AliHLTMUONCoreMansoTracker::kWaitMoreChamber7: os << "kWaitMoreChamber7"; break;
	default:                             os << "FAULT!!"; 
	}
	return os;
}

std::ostream& operator << (std::ostream& os, AliHLTMUONCoreMansoTracker::StatesSM5 state)
{
	switch (state)
	{
	case AliHLTMUONCoreMansoTracker::kSM5Idle:           os << "kSM5Idle"; break;
	case AliHLTMUONCoreMansoTracker::kWaitChamber10:     os << "kWaitChamber10"; break;
	case AliHLTMUONCoreMansoTracker::kWaitMoreChamber10: os << "kWaitMoreChamber10"; break;
	case AliHLTMUONCoreMansoTracker::kWaitChamber9:      os << "kWaitChamber9"; break;
	case AliHLTMUONCoreMansoTracker::kWaitMoreChamber9:  os << "kWaitMoreChamber9"; break;
	case AliHLTMUONCoreMansoTracker::kSM5Done:           os << "kSM5Done"; break;
	default:                              os << "FAULT!!"; 
	}
	return os;
}

} // end of namespace
#endif // DEBUG


// Deviate from the Manso implementation by allowing a and b
// parameters per chamber and not just per station.
// The default values are derived from the work done in
//    "A first algorithm for dimuon High Level Trigger"
//    Ref ID:  ALICE-INT-2002-04 version 1.0
Float AliHLTMUONCoreMansoTracker::fgA7  = 0.016f;
Float AliHLTMUONCoreMansoTracker::fgB7  = 2.0f;
Float AliHLTMUONCoreMansoTracker::fgA8  = 0.016f;
Float AliHLTMUONCoreMansoTracker::fgB8  = 2.0f;
Float AliHLTMUONCoreMansoTracker::fgA9  = 0.020f;
Float AliHLTMUONCoreMansoTracker::fgB9  = 3.0f;
Float AliHLTMUONCoreMansoTracker::fgA10 = 0.020f;
Float AliHLTMUONCoreMansoTracker::fgB10 = 3.0f;
Float AliHLTMUONCoreMansoTracker::fgZ7  = 1274.5f;
Float AliHLTMUONCoreMansoTracker::fgZ8  = 1305.5f;
Float AliHLTMUONCoreMansoTracker::fgZ9  = 1408.6f;
Float AliHLTMUONCoreMansoTracker::fgZ10 = 1439.6f;
Float AliHLTMUONCoreMansoTracker::fgZ11 = 1603.5f;
Float AliHLTMUONCoreMansoTracker::fgZ13 = 1703.5f;


void AliHLTMUONCoreMansoTracker::AliRegionOfInterest::Create(AliHLTMUONCorePoint p, Float a, Float b)
{
// Creates a region of interest specific to the Manso algorithm from a point and
// two Manso specific parameters.

	fCentre = p;
	// Compute the radius Rp
	Float rp = (Float) sqrt( p.X() * p.X() + p.Y() * p.Y() );

	// The radius Rs for the region of interest is computed from the
	// specification given in the document:
	//   "A first algorithm for dimuon High Level Trigger"
	//   Ref ID:  ALICE-INT-2002-04 version 1.0
	//   equation:
	//     Rs = a * Rp + b
	//   given on page 3 section 4.
	fRs = a * rp + b;
}


bool AliHLTMUONCoreMansoTracker::AliRegionOfInterest::Contains(AliHLTMUONCorePoint p) const
{
	// Compute the distance between the centre of the region of interest and
	// the point p. This distance must be less than the radius of the region
	// of interest for p to be contained in the region of interest.
	register Float lx = fCentre.X() - p.X();
	register Float ly = fCentre.Y() - p.Y();
	register Float r = (Float) sqrt( lx * lx + ly * ly );
	DebugMsg(4, "\tAliRegionOfInterest::Contains : p = " << p
		<< " , centre = " << fCentre << " , r = " << r << " , Rs = " << fRs
	);
	return r <= fRs;
}


void AliHLTMUONCoreMansoTracker::AliRegionOfInterest::GetBoundaryBox(
		Float& left, Float& right, Float& bottom, Float& top
	) const
{
// Works out the smallest boundary box that will contain the region of interest.

	left = fCentre.X() - fRs;
	right = fCentre.X() + fRs;
	bottom = fCentre.Y() - fRs;
	top = fCentre.Y() + fRs;
}


AliHLTMUONCoreMansoTracker::AliVertex::AliVertex(Float x, Float y, Float z)
	: fX(x), fY(y), fZ(z)
{
// Constructor for vertex.

	fX = x;
	fY = y;
	fZ = z;
}


AliHLTMUONCoreMansoTracker::AliVertex::AliVertex(AliHLTMUONCorePoint xy, Float z)
	: fX(xy.X()), fY(xy.Y()), fZ(z)
{
// Construct vertex from a point on the XY plane and z coordinate.

	fX = xy.X();
	fY = xy.Y();
	fZ = z;
}


AliHLTMUONCoreMansoTracker::AliLine::AliLine(
        Float ax, Float ay, Float az,
        Float bx, Float by, Float bz
    ) :
	fMx(ax - bx), fMy(ay - by), fMz(az - bz),
	fCx(bx), fCy(by), fCz(bz)
{
// Construct a line defined by L = M*t + C = (A-B)*t + B
// where M and C are 3D vectors and t is a free parameter.
// A = (ax, ay, az) and B = (bx, by, bz)

	fMx = ax - bx;
	fMy = ay - by;
	fMz = az - bz;
	fCx = bx;
	fCy = by;
	fCz = bz;
}


AliHLTMUONCoreMansoTracker::AliLine::AliLine(AliVertex a, AliVertex b) :
	fMx(a.X() - b.X()), fMy(a.Y() - b.Y()), fMz(a.Z() - b.Z()),
	fCx(b.X()), fCy(b.Y()), fCz(b.Z())
{
// Contruct a line to go through two vertices a and b.

	fMx = a.X() - b.X();
	fMy = a.Y() - b.Y();
	fMz = a.Z() - b.Z();
	fCx = b.X();
	fCy = b.Y();
	fCz = b.Z();
}


AliHLTMUONCorePoint AliHLTMUONCoreMansoTracker::AliLine::FindIntersectWithXYPlain(Float z) const
{
// Find the point of intersection of the line and the XY plane at z.

	Assert( fMz != 0.0 );    // Should not have a ray perpendicular to the beam axis.
	Float t = (z - fCz) / fMz;
	Float lx = fMx*t + fCx;
	Float ly = fMy*t + fCy;

	return AliHLTMUONCorePoint(lx, ly);
}


AliHLTMUONCoreMansoTracker::AliHLTMUONCoreMansoTracker() :
	AliHLTMUONCoreTracker(),
	fSm4state(kSM4Idle), fSm5state(kSM5Idle),
	fRequestsCompleted(0), fSt4chamber(kChamber1),
	fV1(), fMc1(), fSt5z(0), fSt5data(), fSt4z(0), fSt4points(),
	fSt5rec(), fFoundPoint()
{
// Default constructor 

	fSm4state = kSM4Idle;
	fSm5state = kSM5Idle;
	fRequestsCompleted = 0;
}


void AliHLTMUONCoreMansoTracker::FindTrack(const AliHLTMUONCoreTriggerRecord& trigger)
{
// Tries to find the track from the trigger seed.

	DebugMsg(4, "SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);
	DebugMsg(1, "Processing trigger with pt = " << trigger.fPt);
	fV1 = AliVertex( trigger.fStation1impact, fgZ11 );
	AliVertex v2 = AliVertex( trigger.fStation2impact, fgZ13 );

	// Form the vector line between the above two impact points and
	// find the crossing point of the line with chamber 10 (i.e. station 5).
	fMc1.fLine = AliLine(fV1, v2);
	AliHLTMUONCorePoint p10 = fMc1.fLine.FindIntersectWithXYPlain( fgZ10 );

	// Build a region of interest for tracking station 5 (chamber 10).
	// Remember the parameters a and b are station specific.
	fMc1.fChamber = kChamber10;
	fMc1.fRoi.Create(p10, fgA10, fgB10);
	
	// Make SM5 state transition before the call to RequestClusters since
	// that method could call one of our methods again, so we need to be
	// in a consistant internal state.
	fSm5state = kWaitChamber10;

	Float left, right, bottom, top;
	fMc1.fRoi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, kChamber10, &fMc1);
}


void AliHLTMUONCoreMansoTracker::ReturnClusters(void* tag, const AliHLTMUONCoreClusterPoint* clusters, UInt count)
{
// Implementation of AliHLTMUONCoreTracker::ReturnClusters.

	Assert( count > 0 );
	Assert( clusters != NULL );
	
	AliTagData* data = (AliTagData*)tag;
	DebugMsg(4, "Got AliHLTMUONCoreMansoTracker::ReturnClusters(tag = " << tag
		<< ", chamber = " << data->fChamber
		<< ", clusters = " << clusters <<  ", count = " << count << ")"
	);
	DebugMsg(4, "SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);

	switch (data->fChamber)
	{
	case kChamber7:  ReceiveClustersChamber7(clusters, count, data); break;
	case kChamber8:  ReceiveClustersChamber8(clusters, count, data); break;
	case kChamber9:  ReceiveClustersChamber9(clusters, count); break;
	case kChamber10: ReceiveClustersChamber10(clusters, count); break;
	default:
		// Error
		DebugMsg(1, "ERROR: Got tag with an invalid value: " << data->fChamber);
	}
}


void AliHLTMUONCoreMansoTracker::EndOfClusters(void* tag)
{
// Implementation of AliHLTMUONCoreTracker::EndOfClusters.

	AliTagData* data = (AliTagData*)tag;
	DebugMsg(4, "Got AliHLTMUONCoreMansoTracker::EndOfClusters(chamber = " << data->fChamber << ")");
	DebugMsg(4, "SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);

	switch (data->fChamber)
	{
	case kChamber7:  EndOfClustersChamber7(); break;
	case kChamber8:  EndOfClustersChamber8(); break;
	case kChamber9:  EndOfClustersChamber9(); break;
	case kChamber10: EndOfClustersChamber10(); break;
	default:
		// Error
		DebugMsg(1, "ERROR: Got tag with an invalid value: " << data->fChamber);
	}
}


void AliHLTMUONCoreMansoTracker::FillTrackData(AliHLTMUONCoreTrack& track)
{
// Implementation of AliHLTMUONCoreTracker::FillTrackData

	DebugMsg(4, "FillTrack: st5 = " << fSt5rec->fClusterPoint << ", st4 = " << fFoundPoint->fClusterPoint);
	
	Float x1 = fFoundPoint->fClusterPoint.X();
	Float y1 = fFoundPoint->fClusterPoint.Y();
	Float y2 = fSt5rec->fClusterPoint.Y();
	Float momentum;
	Float pt = AliHLTMUONCoreCalculateSignedPt(x1, y1, y2, fSt4z, fSt5z, momentum);
	DebugMsg(1, "Calculated Pt = " << pt);
	DebugMsg(1, "\tusing x1 = " << x1 << " , y1 = " << y1 << " , y2 = " << y2
		<< " , z1 = " << fSt4z << " , z2 = " << fSt5z
	);

	if (pt < 0)
		track.fSign = kSignMinus;
	else if (pt > 0)
		track.fSign = kSignPlus;
	else
		track.fSign = kUnknownSign;

	track.fP = momentum;
	track.fPt = (Float) fabs(pt);
	for (UInt i = 0; i < 6; i++)
	{
		track.fPoint[i] = AliHLTMUONCorePoint(0.0, 0.0);
		track.fRegion[i] = kInvalidROI;
	}

	Float left, right, bottom, top;
	
	// Have to create the ROI numbers from the internal region of interest structures.
	fSt5rec->fTag.fRoi.GetBoundaryBox(left, right, bottom, top);
	AliHLTMUONCoreRegionOfInterest region4(left, right, bottom, top, fSt4chamber);
	fMc1.fRoi.GetBoundaryBox(left, right, bottom, top);
	AliHLTMUONCoreRegionOfInterest region5(left, right, bottom, top, fMc1.fChamber);
	
	// Depending on the chamber we received cluster points from, fill the appropriate
	// point and ROI number. This is done for station 4 then 5.
	if (fSt4chamber == kChamber8)
	{
		track.fPoint[6] = AliHLTMUONCorePoint(0.0, 0.0);
		track.fRegion[6] = kInvalidROI;
		track.fPoint[7] = fFoundPoint->fClusterPoint;
		track.fRegion[7] = region4;
	}
	else
	{
		track.fPoint[6] = fFoundPoint->fClusterPoint;
		track.fRegion[6] = region4;
		track.fPoint[7] = AliHLTMUONCorePoint(0.0, 0.0);
		track.fRegion[7] = kInvalidROI;
	}
	if (fMc1.fChamber == kChamber10)
	{
		track.fPoint[8] = AliHLTMUONCorePoint(0.0, 0.0);
		track.fRegion[8] = kInvalidROI;
		track.fPoint[9] = fSt5rec->fClusterPoint;
		track.fRegion[9] = region5;
	}
	else
	{
		track.fPoint[8] = fSt5rec->fClusterPoint;
		track.fRegion[8] = region5;
		track.fPoint[9] = AliHLTMUONCorePoint(0.0, 0.0);
		track.fRegion[9] = kInvalidROI;
	}
}


void AliHLTMUONCoreMansoTracker::Reset()
{
// Implementation of AliHLTMUONCoreTracker::Reset

	DebugMsg(4, "SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);
	fSt5data.Clear();
	fSt4points.Clear();
	fSm4state = kSM4Idle;
	fSm5state = kSM5Idle;
	fRequestsCompleted = 0;
}


// Note: In the following ReceiveClustersXXX and EndOfClustersXXX methods we make
// the state machine transitions before calls to RequestClusters, FoundTrack, 
// NoTrackFound or EndOfClusterRequests. This is important since the callback
// object will make recursive calls to the tracker's methods so we need to maintain
// a consistant internal state.
// The same would go for updating internal variables.
// In general one should only call the callback methods at the end of any of the
// following routines.

void AliHLTMUONCoreMansoTracker::ReceiveClustersChamber7(
		const AliHLTMUONCoreClusterPoint* clusters, UInt count, const AliTagData* data
	)
{
// State change method for Station 4 state machine.

	switch (fSm4state)
	{
	case kWaitChamber7:
		fSm4state = kWaitMoreChamber7;
	
	case kWaitMoreChamber7:
		for (UInt j = 0; j < count; j++)
		{
			AliHLTMUONCoreClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->fRoi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.X() << ", " << cluster.Y() << "] from chamber 7.");
				AliStation4Data* newdata = fSt4points.New();
				newdata->fClusterPoint = cluster;
				newdata->fSt5tag = data;
			}
		}
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in AliHLTMUONCoreMansoTracker::ReceiveClustersChamber7!");
	}
}


void AliHLTMUONCoreMansoTracker::ReceiveClustersChamber8(
		const AliHLTMUONCoreClusterPoint* clusters, UInt count, const AliTagData* data
	)
{
// State change method for Station 4 state machine.

	switch (fSm4state)
	{
	case kWaitChamber8:
		fSm4state = kWaitMoreChamber8;
		fSt4z = fgZ8;
		fSt4chamber = kChamber8;
	
	case kWaitMoreChamber8:
		for (UInt j = 0; j < count; j++)
		{
			AliHLTMUONCoreClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->fRoi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.X() << ", " << cluster.Y() << "] from chamber 8.");
				AliStation4Data* newdata = fSt4points.New();
				newdata->fClusterPoint = cluster;
				newdata->fSt5tag = data;
			}
		}
		break;
		
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in AliHLTMUONCoreMansoTracker::ReceiveClustersChamber8!");
	}
}


void AliHLTMUONCoreMansoTracker::ReceiveClustersChamber9(const AliHLTMUONCoreClusterPoint* clusters, UInt count)
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber9:
		fSm5state = kWaitMoreChamber9;
		fSm4state = kWaitChamber8;  // Start SM4.
	
	case kWaitMoreChamber9:
		for (UInt j = 0; j < count; j++)
		{
			AliHLTMUONCoreClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( fMc1.fRoi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.X() << ", " << cluster.Y() << "] from chamber 9.");
				AliStation5Data* data = fSt5data.New();
				data->fClusterPoint = cluster;
				ProjectToStation4(data, fgZ9);  // This adds a new request for station 4.
			}
		}
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in AliHLTMUONCoreMansoTracker::ReceiveClustersChamber9!");
	}
}


void AliHLTMUONCoreMansoTracker::ReceiveClustersChamber10(const AliHLTMUONCoreClusterPoint* clusters, UInt count)
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber10:
		fSm5state = kWaitMoreChamber10;
		fSt5z = fgZ10;
		fSm4state = kWaitChamber8;  // Start SM4.
	
	case kWaitMoreChamber10:
		for (UInt j = 0; j < count; j++)
		{
			AliHLTMUONCoreClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( fMc1.fRoi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.X() << ", " << cluster.Y() << "] from chamber 10.");
				AliStation5Data* data = fSt5data.New();
				data->fClusterPoint = cluster;
				ProjectToStation4(data, fgZ10);  // This adds a new request for station 4.
			}
		}
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in AliHLTMUONCoreMansoTracker::ReceiveClustersChamber10!");
	}
}


void AliHLTMUONCoreMansoTracker::EndOfClustersChamber7()
{
// State change method for Station 4 state machine.

	fRequestsCompleted++;  // Increment the number of requests completed for station 4.
	DebugMsg(4, "fRequestsCompleted = " << fRequestsCompleted );

	switch (fSm4state)
	{
	case kWaitChamber7:
		// If all data from station 5 is received and no data found on
		// chambers 7 or 8 then we can not find a track.
		if (fSm5state == kSM5Done) NoTrackFound();
		break;
	
	case kWaitMoreChamber7:
		if (fRequestsCompleted == fSt5data.Count() && fSm5state == kSM5Done)
			ProcessClusters();
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in AliHLTMUONCoreMansoTracker::EndOfClustersChamber7!");
	}
}


void AliHLTMUONCoreMansoTracker::EndOfClustersChamber8()
{
// State change method for Station 4 state machine.

	fRequestsCompleted++;  // Increment the number of requests completed for station 4.
	DebugMsg(4, "fRequestsCompleted = " << fRequestsCompleted );

	switch (fSm4state)
	{
	case kWaitChamber7:
		// Ignore. The requests for chamber 8 are already re-requested below.
		break;
		
	case kWaitChamber8:
		{
		fSm4state = kWaitChamber7;
		fSt4z = fgZ7;
		fSt4chamber = kChamber7;
	
		// We need to resend the requests for chamber 8, but change the request
		// to get data for chamber 7 instead:
		UInt reqlistsize = fSt5data.Count();
		DebugMsg(4, "Re-requesting clusters from chamber 7... reqlistsize = " << reqlistsize);

		Station5List::Iterator rec = fSt5data.First();
		for (UInt i = 0; i < reqlistsize; i++, rec++)
		{
			// Need to create a new st5 data block for the request.
			AliStation5Data* data = fSt5data.New();
			data->fClusterPoint = rec->fClusterPoint;
			data->fTag.fLine = rec->fTag.fLine;

			// Rebuild a region of interest for chamber 7.
			// Remember the parameters a and b are station specific.
			AliHLTMUONCorePoint p7 = data->fTag.fLine.FindIntersectWithXYPlain( fgZ7 );
			data->fTag.fChamber = kChamber7;
			data->fTag.fRoi.Create(p7, fgA7, fgB7);
			
			Float left, right, bottom, top;
			data->fTag.fRoi.GetBoundaryBox(left, right, bottom, top);
			// Make request for chamber 7 data.
			RequestClusters(left, right, bottom, top, kChamber7, &data->fTag);
		}
		}
		break;
	
	case kWaitMoreChamber8:
		if (fRequestsCompleted == fSt5data.Count() && fSm5state == kSM5Done)
			ProcessClusters();
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in AliHLTMUONCoreMansoTracker::EndOfClustersChamber8!");
	}
}


void AliHLTMUONCoreMansoTracker::EndOfClustersChamber9()
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber9:
		fSm5state = kSM5Done;
		EndOfClusterRequests();
		NoTrackFound();
		break;
		
	case kWaitMoreChamber9:
		fSm5state = kSM5Done;
		EndOfClusterRequests();
		if (fRequestsCompleted == fSt5data.Count())
			ProcessClusters();
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in AliHLTMUONCoreMansoTracker::EndOfClustersChamber9!");
	}
}


void AliHLTMUONCoreMansoTracker::EndOfClustersChamber10()
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber10:
		{
		fSm5state = kWaitChamber9;
		fSt5z = fgZ9;
		
		// No clusters found on chamber 10 so we need to make a request for
		// clusters from chamber 9:
		AliHLTMUONCorePoint p9 = fMc1.fLine.FindIntersectWithXYPlain( fgZ9 );

		// Build a region of interest for tracking station 5 (chamber 9).
		// Remember the parameters a and b are station specific.
		fMc1.fChamber = kChamber9;
		fMc1.fRoi.Create(p9, fgA9, fgB9);

		Float left, right, bottom, top;
		fMc1.fRoi.GetBoundaryBox(left, right, bottom, top);
		RequestClusters(left, right, bottom, top, kChamber9, &fMc1);
		}
		break;

	case kWaitMoreChamber10:
		fSm5state = kSM5Done;
		EndOfClusterRequests();
		if (fRequestsCompleted == fSt5data.Count())
			ProcessClusters();
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in AliHLTMUONCoreMansoTracker::EndOfClustersChamber10!");
	}
}


void AliHLTMUONCoreMansoTracker::ProjectToStation4(AliStation5Data* data, register Float station5z)
{
	// Perform chamber specific operations:
	// Since certain states of SM4 means that it is fetching for Chamber8
	// and other states are for fetching from Chamber7. We need to make
	// requests for the correct chamber.
	Assert( fSm4state == kWaitChamber8 
		|| fSm4state == kWaitMoreChamber8
		|| fSm4state == kWaitChamber7
		|| fSm4state == kWaitMoreChamber7
	);
	AliTagData* tag = &data->fTag;
	if (fSm4state == kWaitChamber8 || fSm4state == kWaitMoreChamber8)
	{
		// Form the vector line between trigger station 1 and tracking station 5,
		// and find the intersection point of the line with station 4 (chamber8).
		AliLine line51( AliVertex(data->fClusterPoint, station5z), fV1 );
		AliHLTMUONCorePoint intercept = line51.FindIntersectWithXYPlain( fgZ8 );
		tag->fLine = line51;
		
		// Build a region of interest for tracking station 4.
		tag->fChamber = kChamber8;
		tag->fRoi.Create(intercept, fgA8, fgB8);
	}
	else
	{
		// Form the vector line between trigger station 1 and tracking station 5,
		// and find the intersection point of the line with station 4 (chamber7).
		AliLine line51( AliVertex(data->fClusterPoint, station5z), fV1 );
		AliHLTMUONCorePoint intercept = line51.FindIntersectWithXYPlain( fgZ7 );
		tag->fLine = line51;
		
		// Build a region of interest for tracking station 4.
		tag->fChamber = kChamber7;
		tag->fRoi.Create(intercept, fgA7, fgB7);
	}

	// Make the request for clusters from station 4.
	Float left, right, bottom, top;
	tag->fRoi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, tag->fChamber, tag);
}


void AliHLTMUONCoreMansoTracker::ProcessClusters()
{
// Process clusters that have been received.
// This is called once all clusters have been found.

	DebugMsg(2, "ProcessClusters...");
	
	// Check if the cluster point list on station 4 is empty.
	// If it is then we have not found any tracks.
	fFoundPoint = fSt4points.First();
	if (fFoundPoint == fSt4points.End())
	{
		NoTrackFound();
		return;
	}
	
	fSt5rec = fSt5data.First();
	if (fSt5rec != fSt5data.End())
	{
		// Only look at station 5 data records that are for the found chamber number.
		// Note: either we only have chamber 8 data or we have chamber 7 data followed
		// by chamber 8 data.
		// Thus if we hit records that we are not interested in already then the list
		// contains no interesting data and we can signal no track found.
		if (fSt5rec->fTag.fChamber != fSt4chamber)
		{
			NoTrackFound();
			return;
		}
		 
		// For all combinations of cluster point pairs from station 4 and 5
		// signal a found track:
		do
		{
			DebugMsg(4, "\tfSt5rec->fTag.chamber = " << fSt5rec->fTag.fChamber
				<< " , fSt4chamber = " << fSt4chamber
			);

			for (fFoundPoint = fSt4points.First(); fFoundPoint != fSt4points.End(); fFoundPoint++)
			{
				if (fFoundPoint->fSt5tag == &fSt5rec->fTag)
					FoundTrack();
			}

			fSt5rec++;  // Get next station 5 cluster point.
		} while (fSt5rec != fSt5data.End() && fSt5rec->fTag.fChamber == fSt4chamber);
	}
	else
		NoTrackFound();
}

