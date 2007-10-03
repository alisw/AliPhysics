/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/**
 *  @file   AliHLTMUONMansoTrackerFSM.cxx
 *  @author Artur Szostak <artursz@iafrica.com>
 *  @date   
 *  @brief  Implementation of AliHLTMUONMansoTrackerFSM class.
 */

#include "AliHLTMUONMansoTrackerFSM.h"
#include "AliHLTMUONCalculations.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include <cmath>


#ifdef DEBUG
#include <ostream>
namespace
{

std::ostream& operator << (std::ostream& os, AliHLTMUONMansoTrackerFSM::StatesSM4 state)
{
	switch (state)
	{
	case AliHLTMUONMansoTrackerFSM::kSM4Idle:          os << "kSM4Idle"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitChamber8:     os << "kWaitChamber8"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitMoreChamber8: os << "kWaitMoreChamber8"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitChamber7:     os << "kWaitChamber7"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitMoreChamber7: os << "kWaitMoreChamber7"; break;
	default:                                           os << "FAULT!!"; 
	}
	return os;
}

std::ostream& operator << (std::ostream& os, AliHLTMUONMansoTrackerFSM::StatesSM5 state)
{
	switch (state)
	{
	case AliHLTMUONMansoTrackerFSM::kSM5Idle:           os << "kSM5Idle"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitChamber10:     os << "kWaitChamber10"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitMoreChamber10: os << "kWaitMoreChamber10"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitChamber9:      os << "kWaitChamber9"; break;
	case AliHLTMUONMansoTrackerFSM::kWaitMoreChamber9:  os << "kWaitMoreChamber9"; break;
	case AliHLTMUONMansoTrackerFSM::kSM5Done:           os << "kSM5Done"; break;
	default:                                            os << "FAULT!!"; 
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
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgA7  = 0.016f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgB7  = 2.0f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgA8  = 0.016f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgB8  = 2.0f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgA9  = 0.020f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgB9  = 3.0f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgA10 = 0.020f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgB10 = 3.0f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ7  = -1274.5f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ8  = -1305.5f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ9  = -1408.6f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ10 = -1439.6f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ11 = -1603.5f;
AliHLTFloat32_t AliHLTMUONMansoTrackerFSM::fgZ13 = -1703.5f;


void AliHLTMUONMansoTrackerFSM::AliRegionOfInterest::Create(
		AliHLTMUONRecHitStruct p, AliHLTFloat32_t a, AliHLTFloat32_t b
	)
{
// Creates a region of interest specific to the Manso algorithm from a point and
// two Manso specific parameters.

	fCentre = p;
	// Compute the radius Rp
	AliHLTFloat32_t rp = (AliHLTFloat32_t) sqrt( p.fX * p.fX + p.fY * p.fY );

	// The radius Rs for the region of interest is computed from the
	// specification given in the document:
	//   "A first algorithm for dimuon High Level Trigger"
	//   Ref ID:  ALICE-INT-2002-04 version 1.0
	//   equation:
	//     Rs = a * Rp + b
	//   given on page 3 section 4.
	fRs = a * rp + b;
}


bool AliHLTMUONMansoTrackerFSM::AliRegionOfInterest::Contains(AliHLTMUONRecHitStruct p) const
{
	// Compute the distance between the centre of the region of interest and
	// the point p. This distance must be less than the radius of the region
	// of interest for p to be contained in the region of interest.
	register AliHLTFloat32_t lx = fCentre.fX - p.fX;
	register AliHLTFloat32_t ly = fCentre.fY - p.fY;
	register AliHLTFloat32_t r = (AliHLTFloat32_t) sqrt( lx * lx + ly * ly );
	DebugTrace("\tAliRegionOfInterest::Contains : p = " << p
		<< " , centre = " << fCentre << " , r = " << r << " , Rs = " << fRs
	);
	return r <= fRs;
}


void AliHLTMUONMansoTrackerFSM::AliRegionOfInterest::GetBoundaryBox(
		AliHLTFloat32_t& left, AliHLTFloat32_t& right,
		AliHLTFloat32_t& bottom, AliHLTFloat32_t& top
	) const
{
// Works out the smallest boundary box that will contain the region of interest.

	left = fCentre.fX - fRs;
	right = fCentre.fX + fRs;
	bottom = fCentre.fY - fRs;
	top = fCentre.fY + fRs;
}


AliHLTMUONMansoTrackerFSM::AliVertex::AliVertex(
		AliHLTFloat32_t x, AliHLTFloat32_t y, AliHLTFloat32_t z
	)
	: fX(x), fY(y), fZ(z)
{
// Constructor for vertex.
}


AliHLTMUONMansoTrackerFSM::AliVertex::AliVertex(AliHLTMUONRecHitStruct xy, AliHLTFloat32_t z)
	: fX(xy.fX), fY(xy.fY), fZ(z)
{
// Construct vertex from a point on the XY plane and z coordinate.
}


AliHLTMUONMansoTrackerFSM::AliLine::AliLine(
		AliHLTFloat32_t ax, AliHLTFloat32_t ay, AliHLTFloat32_t az,
		AliHLTFloat32_t bx, AliHLTFloat32_t by, AliHLTFloat32_t bz
	) :
	fMx(ax - bx),
	fMy(ay - by),
	fMz(az - bz),
	fCx(bx),
	fCy(by),
	fCz(bz)
{
// Construct a line defined by L = M*t + C = (A-B)*t + B
// where M and C are 3D vectors and t is a free parameter.
// A = (ax, ay, az) and B = (bx, by, bz)
}


AliHLTMUONMansoTrackerFSM::AliLine::AliLine(AliVertex a, AliVertex b) :
	fMx(a.X() - b.X()),
	fMy(a.Y() - b.Y()),
	fMz(a.Z() - b.Z()),
	fCx(b.X()),
	fCy(b.Y()),
	fCz(b.Z())
{
// Contruct a line to go through two vertices a and b.
}


AliHLTMUONRecHitStruct AliHLTMUONMansoTrackerFSM::AliLine::FindIntersectWithXYPlain(
		AliHLTFloat32_t z
	) const
{
// Find the point of intersection of the line and the XY plane at z.

	assert( fMz != 0.0 );  // Should not have a ray perpendicular to the beam axis.
	AliHLTFloat32_t t = (z - fCz) / fMz;
	AliHLTMUONRecHitStruct p;
	p.fX = fMx*t + fCx;
	p.fY = fMy*t + fCy;
	p.fZ = z;
	return p;
}


AliHLTMUONMansoTrackerFSM::AliHLTMUONMansoTrackerFSM() :
	fCallback(NULL),
	fSm4state(kSM4Idle),
	fSm5state(kSM5Idle),
	fRequestsCompleted(0),
	fSt4chamber(kChamber1),
	fV1(),
	fMc1(),
	fSt5data(),
	fSt4points(),
	fSt5rec(),
	fFoundPoint(),
	fTriggerId(-1),
	fTrackId(0)
{
// Default constructor
}


void AliHLTMUONMansoTrackerFSM::FindTrack(const AliHLTMUONTriggerRecordStruct& trigger)
{
// Tries to find the track from the trigger seed.

	DebugTrace("SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);
	DebugTrace("Processing trigger with ID = " << trigger.fId);
	
	fTriggerId = trigger.fId;
	
	AliHLTMUONParticleSign sign;
	bool hitset[4];
	AliHLTMUONUtils::UnpackTriggerRecordFlags(trigger.fFlags, sign, hitset);	
	DebugTrace("hitset = {" << hitset[0] << ", " << hitset[1] << ", "
		<< hitset[2] << ", " << hitset[3] << "}"
	); 
	
	// Find first valid hit in the trigger record from chambers 11 and 12.
	int firstHit = -1;
	for (int i = 0; i < 2; i++)
	{
		if (hitset[i])
		{
			firstHit = i;
			break;
		}
	}
	// Find the second valid hit in the trigger record from chambers 13 and 14.
	int secondHit = -1;
	for (int i = 2; i < 4; i++)
	{
		if (hitset[i])
		{
			secondHit = i;
			break;
		}
	}
	if (firstHit == -1 or secondHit == -1)
	{
		NoTrackFound();
		return;
	}
	
	fV1 = AliVertex(
			trigger.fHit[firstHit].fX,
			trigger.fHit[firstHit].fY,
			trigger.fHit[firstHit].fZ
		);
	AliVertex v2 = AliVertex(
			trigger.fHit[secondHit].fX,
			trigger.fHit[secondHit].fY,
			trigger.fHit[secondHit].fZ
		);
	
	DebugTrace("Using fV1 = {x = " << fV1.X() << ", y = " << fV1.Y() << ", "
		<< fV1.Z() << "}, with firstHit = " << firstHit
	);
	DebugTrace("Using v2 = {x = " << v2.X() << ", y = " << v2.Y() << ", "
		<< v2.Z() << "}, with secondHit = " << secondHit
	);

	// Form the vector line between the above two impact points and
	// find the crossing point of the line with chamber 10 (i.e. station 5).
	fMc1.fLine = AliLine(fV1, v2);
	AliHLTMUONRecHitStruct p10 = fMc1.fLine.FindIntersectWithXYPlain( fgZ10 );

	// Build a region of interest for tracking station 5 (chamber 10).
	// Remember the parameters a and b are station specific.
	fMc1.fChamber = kChamber10;
	fMc1.fRoi.Create(p10, fgA10, fgB10);
	
	// Make SM5 state transition before the call to RequestClusters since
	// that method could call one of our methods again, so we need to be
	// in a consistant internal state.
	fSm5state = kWaitChamber10;

	AliHLTFloat32_t left, right, bottom, top;
	fMc1.fRoi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, kChamber10, &fMc1);
}


void AliHLTMUONMansoTrackerFSM::ReturnClusters(
		void* tag, const AliHLTMUONRecHitStruct* clusters,
		AliHLTUInt32_t count
	)
{
// Implementation of AliHLTMUONMansoTrackerFSM::ReturnClusters.

	assert( count > 0 );
	assert( clusters != NULL );
	
	AliTagData* data = (AliTagData*)tag;
	DebugTrace("Got AliHLTMUONMansoTrackerFSM::ReturnClusters(tag = " << tag
		<< ", chamber = " << data->fChamber
		<< ", clusters = " << clusters <<  ", count = " << count << ")"
	);
	DebugTrace("SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);

	switch (data->fChamber)
	{
	case kChamber7:  ReceiveClustersChamber7(clusters, count, data); break;
	case kChamber8:  ReceiveClustersChamber8(clusters, count, data); break;
	case kChamber9:  ReceiveClustersChamber9(clusters, count); break;
	case kChamber10: ReceiveClustersChamber10(clusters, count); break;
	default:
		// Error
		DebugTrace("ERROR: Got tag with an invalid value: " << data->fChamber);
	}
}


void AliHLTMUONMansoTrackerFSM::EndOfClusters(void* tag)
{
// Implementation of AliHLTMUONMansoTrackerFSM::EndOfClusters.

	AliTagData* data = (AliTagData*)tag;
	DebugTrace("Got AliHLTMUONMansoTrackerFSM::EndOfClusters(chamber = " << data->fChamber << ")");
	DebugTrace("SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);

	switch (data->fChamber)
	{
	case kChamber7:  EndOfClustersChamber7(); break;
	case kChamber8:  EndOfClustersChamber8(); break;
	case kChamber9:  EndOfClustersChamber9(); break;
	case kChamber10: EndOfClustersChamber10(); break;
	default:
		// Error
		DebugTrace("ERROR: Got tag with an invalid value: " << data->fChamber);
	}
}


bool AliHLTMUONMansoTrackerFSM::FillTrackData(AliHLTMUONMansoTrackStruct& track)
{
// Implementation of AliHLTMUONMansoTrackerFSM::FillTrackData

	DebugTrace("FillTrack: st5 = " << fSt5rec->fClusterPoint << ", st4 = " << fFoundPoint->fClusterPoint);
	
	track.fId = fTrackId;
	// Increment track Id and keep it positive.
	//TODO: handle the wrapparound better.
	if (fTrackId < 0x7FFFFFFF)
		fTrackId++;
	else
		fTrackId = 0;
	
	track.fTrigRec = fTriggerId;
	
	AliHLTFloat32_t x1 = fFoundPoint->fClusterPoint.fX;
	AliHLTFloat32_t y1 = fFoundPoint->fClusterPoint.fY;
	AliHLTFloat32_t z1 = fFoundPoint->fClusterPoint.fZ;
	AliHLTFloat32_t y2 = fSt5rec->fClusterPoint.fY;
	AliHLTFloat32_t z2 = fSt5rec->fClusterPoint.fZ;
	
	bool calculated = AliHLTMUONCalculations::ComputeMomentum(x1, y1, y2, z1, z2);
	
	track.fPx = AliHLTMUONCalculations::Px();
	track.fPy = AliHLTMUONCalculations::Py();
	track.fPz = AliHLTMUONCalculations::Pz();
	DebugTrace("Calculated Px = " << track.fPx << ", Py = " << track.fPy
		<< ", Pz = " << track.fPx
	);
	DebugTrace("\tusing x1 = " << x1 << " , y1 = " << y1 << " , y2 = " << y2
		<< " , z1 = " << z1 << " , z2 = " << z2
	);
	
	track.fChi2 = 0;
	
	bool hitset[4];
	// Depending on which chamber we found reconstructed hits, fill the hit
	// in the appropriate location. This is done for station 4 then 5.
	if (fSt4chamber == kChamber8)
	{
		track.fHit[0] = AliHLTMUONConstants::NilRecHitStruct();
		hitset[0] = false;
		track.fHit[1] = fFoundPoint->fClusterPoint;
		hitset[1] = true;
	}
	else
	{
		track.fHit[0] = fFoundPoint->fClusterPoint;
		hitset[0] = true;
		track.fHit[1] = AliHLTMUONConstants::NilRecHitStruct();
		hitset[1] = false;
	}
	if (fMc1.fChamber == kChamber10)
	{
		track.fHit[2] = AliHLTMUONConstants::NilRecHitStruct();
		hitset[2] = false;
		track.fHit[3] = fSt5rec->fClusterPoint;
		hitset[3] = true;
	}
	else
	{
		track.fHit[2] = fSt5rec->fClusterPoint;
		hitset[2] = true;
		track.fHit[3] = AliHLTMUONConstants::NilRecHitStruct();
		hitset[3] = false;
	}
	
	track.fFlags = AliHLTMUONUtils::PackMansoTrackFlags(
			AliHLTMUONCalculations::Sign(), hitset
		);
	return calculated;
}


void AliHLTMUONMansoTrackerFSM::Reset()
{
// Implementation of AliHLTMUONMansoTrackerFSM::Reset

	DebugTrace("SM5 state = " << fSm5state << " , SM4 state = " << fSm4state);
	fSt5data.Clear();
	fSt4points.Clear();
	fSm4state = kSM4Idle;
	fSm5state = kSM5Idle;
	fRequestsCompleted = 0;
	fTriggerId = -1;
}


// Note: In the following ReceiveClustersXXX and EndOfClustersXXX methods we make
// the state machine transitions before calls to RequestClusters, FoundTrack, 
// NoTrackFound or EndOfClusterRequests. This is important since the callback
// object will make recursive calls to the tracker's methods so we need to maintain
// a consistant internal state.
// The same would go for updating internal variables.
// In general one should only call the callback methods at the end of any of the
// following routines.

void AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber7(
		const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count,
		const AliTagData* data
	)
{
// State change method for Station 4 state machine.

	switch (fSm4state)
	{
	case kWaitChamber7:
		fSm4state = kWaitMoreChamber7;
	
	case kWaitMoreChamber7:
		for (AliHLTUInt32_t j = 0; j < count; j++)
		{
			AliHLTMUONRecHitStruct cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->fRoi.Contains(cluster) )
			{
				DebugTrace("Adding cluster [" << cluster.fX << ", " << cluster.fY << "] from chamber 7.");
				AliStation4Data* newdata = fSt4points.Add();
				newdata->fClusterPoint = cluster;
				newdata->fSt5tag = data;
			}
		}
		break;
	
	default:
		DebugTrace("ERROR: Unexpected state for SM4 in AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber7!");
	}
}


void AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber8(
		const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count,
		const AliTagData* data
	)
{
// State change method for Station 4 state machine.

	switch (fSm4state)
	{
	case kWaitChamber8:
		fSm4state = kWaitMoreChamber8;
		fSt4chamber = kChamber8;
	
	case kWaitMoreChamber8:
		for (AliHLTUInt32_t j = 0; j < count; j++)
		{
			AliHLTMUONRecHitStruct cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->fRoi.Contains(cluster) )
			{
				DebugTrace("Adding cluster [" << cluster.fX << ", " << cluster.fY << "] from chamber 8.");
				AliStation4Data* newdata = fSt4points.Add();
				newdata->fClusterPoint = cluster;
				newdata->fSt5tag = data;
			}
		}
		break;
		
	default:
		DebugTrace("ERROR: Unexpected state for SM4 in AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber8!");
	}
}


void AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber9(
		const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count
	)
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber9:
		fSm5state = kWaitMoreChamber9;
		fSm4state = kWaitChamber8;  // Start SM4.
	
	case kWaitMoreChamber9:
		for (AliHLTUInt32_t j = 0; j < count; j++)
		{
			AliHLTMUONRecHitStruct cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( fMc1.fRoi.Contains(cluster) )
			{
				DebugTrace("Adding cluster [" << cluster.fX << ", " << cluster.fY << "] from chamber 9.");
				AliStation5Data* data = fSt5data.Add();
				data->fClusterPoint = cluster;
				ProjectToStation4(data, fgZ9);  // This adds a new request for station 4.
			}
		}
		break;

	default:
		DebugTrace("ERROR: Unexpected state for SM5 in AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber9!");
	}
}


void AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber10(
		const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count
	)
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber10:
		fSm5state = kWaitMoreChamber10;
		fSm4state = kWaitChamber8;  // Start SM4.
	
	case kWaitMoreChamber10:
		for (AliHLTUInt32_t j = 0; j < count; j++)
		{
			AliHLTMUONRecHitStruct cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( fMc1.fRoi.Contains(cluster) )
			{
				DebugTrace("Adding cluster [" << cluster.fX << ", " << cluster.fY << "] from chamber 10.");
				AliStation5Data* data = fSt5data.Add();
				data->fClusterPoint = cluster;
				ProjectToStation4(data, fgZ10);  // This adds a new request for station 4.
			}
		}
		break;

	default:
		DebugTrace("ERROR: Unexpected state for SM5 in AliHLTMUONMansoTrackerFSM::ReceiveClustersChamber10!");
	}
}


void AliHLTMUONMansoTrackerFSM::EndOfClustersChamber7()
{
// State change method for Station 4 state machine.

	fRequestsCompleted++;  // Increment the number of requests completed for station 4.
	DebugTrace("fRequestsCompleted = " << fRequestsCompleted );

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
		DebugTrace("ERROR: Unexpected state for SM4 in AliHLTMUONMansoTrackerFSM::EndOfClustersChamber7!");
	}
}


void AliHLTMUONMansoTrackerFSM::EndOfClustersChamber8()
{
// State change method for Station 4 state machine.

	fRequestsCompleted++;  // Increment the number of requests completed for station 4.
	DebugTrace("fRequestsCompleted = " << fRequestsCompleted );

	switch (fSm4state)
	{
	case kWaitChamber7:
		// Ignore. The requests for chamber 8 are already re-requested below.
		break;
		
	case kWaitChamber8:
		{
		fSm4state = kWaitChamber7;
		fSt4chamber = kChamber7;
	
		// We need to resend the requests for chamber 8, but change the request
		// to get data for chamber 7 instead:
		AliHLTUInt32_t reqlistsize = fSt5data.Count();
		DebugTrace("Re-requesting clusters from chamber 7... reqlistsize = " << reqlistsize);

		Station5List::Iterator rec = fSt5data.First();
		for (AliHLTUInt32_t i = 0; i < reqlistsize; i++, rec++)
		{
			// Need to create a new st5 data block for the request.
			AliStation5Data* data = fSt5data.Add();
			data->fClusterPoint = rec->fClusterPoint;
			data->fTag.fLine = rec->fTag.fLine;

			// Rebuild a region of interest for chamber 7.
			// Remember the parameters a and b are station specific.
			AliHLTMUONRecHitStruct p7 = data->fTag.fLine.FindIntersectWithXYPlain( fgZ7 );
			data->fTag.fChamber = kChamber7;
			data->fTag.fRoi.Create(p7, fgA7, fgB7);
			
			AliHLTFloat32_t left, right, bottom, top;
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
		DebugTrace("ERROR: Unexpected state for SM4 in AliHLTMUONMansoTrackerFSM::EndOfClustersChamber8!");
	}
}


void AliHLTMUONMansoTrackerFSM::EndOfClustersChamber9()
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
		DebugTrace("ERROR: Unexpected state for SM5 in AliHLTMUONMansoTrackerFSM::EndOfClustersChamber9!");
	}
}


void AliHLTMUONMansoTrackerFSM::EndOfClustersChamber10()
{
// State change method for Station 5 state machine.

	switch (fSm5state)
	{
	case kWaitChamber10:
		{
		fSm5state = kWaitChamber9;
		
		// No clusters found on chamber 10 so we need to make a request for
		// clusters from chamber 9:
		AliHLTMUONRecHitStruct p9 = fMc1.fLine.FindIntersectWithXYPlain( fgZ9 );

		// Build a region of interest for tracking station 5 (chamber 9).
		// Remember the parameters a and b are station specific.
		fMc1.fChamber = kChamber9;
		fMc1.fRoi.Create(p9, fgA9, fgB9);

		AliHLTFloat32_t left, right, bottom, top;
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
		DebugTrace("ERROR: Unexpected state for SM5 in AliHLTMUONMansoTrackerFSM::EndOfClustersChamber10!");
	}
}


void AliHLTMUONMansoTrackerFSM::ProjectToStation4(
		AliStation5Data* data, register AliHLTFloat32_t station5z
	)
{
	// Perform chamber specific operations:
	// Since certain states of SM4 means that it is fetching for Chamber8
	// and other states are for fetching from Chamber7. We need to make
	// requests for the correct chamber.
	assert( fSm4state == kWaitChamber8 
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
		AliHLTMUONRecHitStruct intercept = line51.FindIntersectWithXYPlain( fgZ8 );
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
		AliHLTMUONRecHitStruct intercept = line51.FindIntersectWithXYPlain( fgZ7 );
		tag->fLine = line51;
		
		// Build a region of interest for tracking station 4.
		tag->fChamber = kChamber7;
		tag->fRoi.Create(intercept, fgA7, fgB7);
	}

	// Make the request for clusters from station 4.
	AliHLTFloat32_t left, right, bottom, top;
	tag->fRoi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, tag->fChamber, tag);
}


void AliHLTMUONMansoTrackerFSM::ProcessClusters()
{
// Process clusters that have been received.
// This is called once all clusters have been found.

	DebugTrace("ProcessClusters...");
	
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
			DebugTrace("\tfSt5rec->fTag.chamber = " << fSt5rec->fTag.fChamber
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

