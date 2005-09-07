////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "Tracking/MansoTracker.hpp"
#include <math.h>
#include "Tracking/Calculations.hpp"
#include "Utils.hpp"
#include "Error.hpp"
#include "RegionOfInterest.hpp"


#if defined(DEBUG) || (defined(USE_ALILOG) && ! defined(LOG_NO_DEBUG))
#include <ostream>
#include "Debug/print.hpp"
namespace
{

using dHLT::Tracking::MansoTracker;

std::ostream& operator << (std::ostream& os, MansoTracker::StatesSM4 state)
{
	switch (state)
	{
	case MansoTracker::SM4Idle:          os << "SM4Idle"; break;
	case MansoTracker::WaitChamber8:     os << "WaitChamber8"; break;
	case MansoTracker::WaitMoreChamber8: os << "WaitMoreChamber8"; break;
	case MansoTracker::WaitChamber7:     os << "WaitChamber7"; break;
	case MansoTracker::WaitMoreChamber7: os << "WaitMoreChamber7"; break;
	default:                             os << "FAULT!!"; 
	};
	return os;
};

std::ostream& operator << (std::ostream& os, MansoTracker::StatesSM5 state)
{
	switch (state)
	{
	case MansoTracker::SM5Idle:           os << "SM5Idle"; break;
	case MansoTracker::WaitChamber10:     os << "WaitChamber10"; break;
	case MansoTracker::WaitMoreChamber10: os << "WaitMoreChamber10"; break;
	case MansoTracker::WaitChamber9:      os << "WaitChamber9"; break;
	case MansoTracker::WaitMoreChamber9:  os << "WaitMoreChamber9"; break;
	case MansoTracker::SM5Done:           os << "SM5Done"; break;
	default:                              os << "FAULT!!"; 
	};
	return os;
};

}; // end of namespace
#endif // DEBUG


namespace dHLT
{
namespace Tracking
{


// Deviate from the Manso implementation by allowing a and b
// parameters per chamber and not just per station.
// The default values are derived from the work done in
//    "A first algorithm for dimuon High Level Trigger"
//    Ref ID:  ALICE-INT-2002-04 version 1.0
Float MansoTracker::a7  = 0.016f;
Float MansoTracker::b7  = 2.0f;
Float MansoTracker::a8  = 0.016f;
Float MansoTracker::b8  = 2.0f;
Float MansoTracker::a9  = 0.020f;
Float MansoTracker::b9  = 3.0f;
Float MansoTracker::a10 = 0.020f;
Float MansoTracker::b10 = 3.0f;
Float MansoTracker::z7  = 1251.5f;
Float MansoTracker::z8  = 1278.5f;
Float MansoTracker::z9  = 1416.5f;
Float MansoTracker::z10 = 1443.5f;
Float MansoTracker::z11 = 1603.5f;
Float MansoTracker::z13 = 1703.5f;


void MansoTracker::RegionOfInterest::Create(Point p, Float a, Float b)
{
	centre = p;
	// Compute the radius Rp
	Float Rp = (Float) sqrt( p.x * p.x + p.y * p.y );

	// The radius Rs for the region of interest is computed from the
	// specification given in the document:
	//   "A first algorithm for dimuon High Level Trigger"
	//   Ref ID:  ALICE-INT-2002-04 version 1.0
	//   equation:
	//     Rs = a * Rp + b
	//   given on page 3 section 4.
	Rs = a * Rp + b;
};


bool MansoTracker::RegionOfInterest::Contains(Point p) const
{
	// Compute the distance between the centre of the region of interest and
	// the point p. This distance must be less than the radius of the region
	// of interest for p to be contained in the region of interest.
	register Float lx = centre.x - p.x;
	register Float ly = centre.y - p.y;
	register Float r = (Float) sqrt( lx * lx + ly * ly );
	DebugMsg(4, "\tRegionOfInterest::Contains : p = " << p
		<< " , centre = " << centre << " , r = " << r << " , Rs = " << Rs
	);
	return r <= Rs;
};


void MansoTracker::RegionOfInterest::GetBoundaryBox(Float& left, Float& right, Float& bottom, Float& top)
{
	left = centre.x - Rs;
	right = centre.x + Rs;
	bottom = centre.y - Rs;
	top = centre.y + Rs;
};


MansoTracker::Vertex::Vertex(Float x, Float y, Float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
};


MansoTracker::Vertex::Vertex(Point xy, Float z)
{
	x = xy.x;
	y = xy.y;
	this->z = z;
};


MansoTracker::Line::Line(
        Float Ax, Float Ay, Float Az,
        Float Bx, Float By, Float Bz
    )
{
	Mx = Ax - Bx;
	My = Ay - By;
	Mz = Az - Bz;
	Cx = Bx;
	Cy = By;
	Cz = Bz;
};


MansoTracker::Line::Line(Vertex A, Vertex B)
{
	Mx = A.x - B.x;
	My = A.y - B.y;
	Mz = A.z - B.z;
	Cx = B.x;
	Cy = B.y;
	Cz = B.z;
};


Point MansoTracker::Line::FindIntersectWithXYPlain(Float z) const
{
	Assert( Mz != 0.0 );    // Should not have a ray perpendicular to the beam axis.
	Float t = (z - Cz) / Mz;
	Float Lx = Mx*t + Cx;
	Float Ly = My*t + Cy;

	return Point(Lx, Ly);
};


MansoTracker::MansoTracker() : Tracker()
{
	sm4state = SM4Idle;
	sm5state = SM5Idle;
	requests_completed = 0;
};


void MansoTracker::FindTrack(const TriggerRecord& trigger)
{
	DebugMsg(4, "SM5 state = " << sm5state << " , SM4 state = " << sm4state);
	DebugMsg(1, "Processing trigger with pt = " << trigger.pt);
	v1 = Vertex( trigger.station1impact, z11 );
	Vertex v2 = Vertex( trigger.station2impact, z13 );

	// Form the vector line between the above two impact points and
	// find the crossing point of the line with chamber 10 (i.e. station 5).
	mc1.line = Line(v1, v2);
	Point p10 = mc1.line.FindIntersectWithXYPlain( z10 );

	// Build a region of interest for tracking station 5 (chamber 10).
	// Remember the parameters a and b are station specific.
	mc1.chamber = Chamber10;
	mc1.roi.Create(p10, a10, b10);
	
	// Make SM5 state transition before the call to RequestClusters since
	// that method could call one of our methods again, so we need to be
	// in a consistant internal state.
	sm5state = WaitChamber10;

	Float left, right, bottom, top;
	mc1.roi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, Chamber10, &mc1);
};


void MansoTracker::ReturnClusters(void* tag, const ClusterPoint* clusters, UInt count)
{
	Assert( count > 0 );
	Assert( clusters != NULL );
	
	TagData* data = (TagData*)tag;
	DebugMsg(4, "Got MansoTracker::ReturnClusters(tag = " << tag
		<< ", chamber = " << data->chamber
		<< ", clusters = " << clusters <<  ", count = " << count << ")"
	);
	DebugMsg(4, "SM5 state = " << sm5state << " , SM4 state = " << sm4state);

	switch (data->chamber)
	{
	case Chamber7:  ReceiveClustersChamber7(clusters, count, data); break;
	case Chamber8:  ReceiveClustersChamber8(clusters, count, data); break;
	case Chamber9:  ReceiveClustersChamber9(clusters, count); break;
	case Chamber10: ReceiveClustersChamber10(clusters, count); break;
	default:
		// Error
		DebugMsg(1, "ERROR: Got tag with an invalid value: " << data->chamber);
	};
};


void MansoTracker::EndOfClusters(void* tag)
{
	TagData* data = (TagData*)tag;
	DebugMsg(4, "Got MansoTracker::EndOfClusters(chamber = " << data->chamber << ")");
	DebugMsg(4, "SM5 state = " << sm5state << " , SM4 state = " << sm4state);

	switch (data->chamber)
	{
	case Chamber7:  EndOfClustersChamber7(); break;
	case Chamber8:  EndOfClustersChamber8(); break;
	case Chamber9:  EndOfClustersChamber9(); break;
	case Chamber10: EndOfClustersChamber10(); break;
	default:
		// Error
		DebugMsg(1, "ERROR: Got tag with an invalid value: " << data->chamber);
	};
};


void MansoTracker::FillTrackData(Track& track)
{
	DebugMsg(4, "FillTrack: st5 = " << st5rec->clusterpoint << ", st4 = " << *foundpoint);
	
	Float x1 = foundpoint->x;
	Float y1 = foundpoint->y;
	Float y2 = st5rec->clusterpoint.y;
	Float momentum;
	Float pt = CalculateSignedPt(x1, y1, y2, st4z, st5z, momentum);
	DebugMsg(1, "Calculated Pt = " << pt);
	DebugMsg(1, "\tusing x1 = " << x1 << " , y1 = " << y1 << " , y2 = " << y2
		<< " , z1 = " << st4z << " , z2 = " << st5z
	);

	if (pt < 0)
		track.sign = Minus;
	else if (pt > 0)
		track.sign = Plus;
	else
		track.sign = UnknownSign;

	track.p = momentum;
	track.pt = (Float) fabs(pt);
	for (UInt i = 0; i < 6; i++)
	{
		track.point[i] = Point(0.0, 0.0);
		track.region[i] = INVALID_ROI;
	};

	Float left, right, bottom, top;
	
	// Have to create the ROI numbers from the internal region of interest structures.
	st5rec->tag.roi.GetBoundaryBox(left, right, bottom, top);
	dHLT::RegionOfInterest region4(left, right, bottom, top, st4chamber);
	mc1.roi.GetBoundaryBox(left, right, bottom, top);
	dHLT::RegionOfInterest region5(left, right, bottom, top, mc1.chamber);
	
	// Depending on the chamber we received cluster points from, fill the appropriate
	// point and ROI number. This is done for station 4 then 5.
	if (st4chamber == Chamber8)
	{
		track.point[6] = Point(0.0, 0.0);
		track.region[6] = INVALID_ROI;
		track.point[7] = *foundpoint;
		track.region[7] = region4;
	}
	else
	{
		track.point[6] = *foundpoint;
		track.region[6] = region4;
		track.point[7] = Point(0.0, 0.0);
		track.region[7] = INVALID_ROI;
	};
	if (mc1.chamber == Chamber10)
	{
		track.point[8] = Point(0.0, 0.0);
		track.region[8] = INVALID_ROI;
		track.point[9] = st5rec->clusterpoint;
		track.region[9] = region5;
	}
	else
	{
		track.point[8] = st5rec->clusterpoint;
		track.region[8] = region5;
		track.point[9] = Point(0.0, 0.0);
		track.region[9] = INVALID_ROI;
	};
};


void MansoTracker::Reset()
{
	DebugMsg(4, "SM5 state = " << sm5state << " , SM4 state = " << sm4state);
	st5data.Clear();
	st4points.Clear();
	sm4state = SM4Idle;
	sm5state = SM5Idle;
	requests_completed = 0;
};


// Note: In the following ReceiveClustersXXX and EndOfClustersXXX methods we make
// the state machine transitions before calls to RequestClusters, FoundTrack, 
// NoTrackFound or EndOfClusterRequests. This is important since the callback
// object will make recursive calls the trackers methods so we need to maintain
// a consistant internal state.
// The same would go for updating internal variables.
// In general one should only call the callback methods at the end of any of the
// following routines.

void MansoTracker::ReceiveClustersChamber7(const ClusterPoint* clusters, UInt count, const TagData* data)
{
	switch (sm4state)
	{
	case WaitChamber7:
		sm4state = WaitMoreChamber7;
	
	case WaitMoreChamber7:
		for (UInt j = 0; j < count; j++)
		{
			ClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->roi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.x << ", " << cluster.y << "] from chamber 7.");
				st4points.Add(cluster);
			};
		};
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in MansoTracker::ReceiveClustersChamber7!");
	};
};


void MansoTracker::ReceiveClustersChamber8(const ClusterPoint* clusters, UInt count, const TagData* data)
{
	switch (sm4state)
	{
	case WaitChamber8:
		sm4state = WaitMoreChamber8;
		st4z = z8;
		st4chamber = Chamber8;
	
	case WaitMoreChamber8:
		for (UInt j = 0; j < count; j++)
		{
			ClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 4.
			if ( data->roi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.x << ", " << cluster.y << "] from chamber 8.");
				st4points.Add(cluster);
			};
		};
		break;
		
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in MansoTracker::ReceiveClustersChamber8!");
	};
};


void MansoTracker::ReceiveClustersChamber9(const ClusterPoint* clusters, UInt count)
{
	switch (sm5state)
	{
	case WaitChamber9:
		sm5state = WaitMoreChamber9;
		sm4state = WaitChamber8;  // Start SM4.
	
	case WaitMoreChamber9:
		for (UInt j = 0; j < count; j++)
		{
			ClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( mc1.roi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.x << ", " << cluster.y << "] from chamber 9.");
				Station5Data* data = st5data.New();
				data->clusterpoint = cluster;
				ProjectToStation4(data, z9);  // This adds a new request for station 4.
			};
		};
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in MansoTracker::ReceiveClustersChamber9!");
	};
};


void MansoTracker::ReceiveClustersChamber10(const ClusterPoint* clusters, UInt count)
{
	switch (sm5state)
	{
	case WaitChamber10:
		sm5state = WaitMoreChamber10;
		st5z = z10;
		sm4state = WaitChamber8;  // Start SM4.
	
	case WaitMoreChamber10:
		for (UInt j = 0; j < count; j++)
		{
			ClusterPoint cluster = clusters[j];
			// Check that the cluster actually is in our region of interest on station 5.
			if ( mc1.roi.Contains(cluster) )
			{
				DebugMsg(4, "Adding cluster [" << cluster.x << ", " << cluster.y << "] from chamber 10.");
				Station5Data* data = st5data.New();
				data->clusterpoint = cluster;
				ProjectToStation4(data, z10);  // This adds a new request for station 4.
			};
		};
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in MansoTracker::ReceiveClustersChamber10!");
	};
};


void MansoTracker::EndOfClustersChamber7()
{
	requests_completed++;  // Increment the number of requests completed for station 4.
	DebugMsg(4, "requests_completed = " << requests_completed );

	switch (sm4state)
	{
	case WaitChamber7:
		// If all data from station 5 is received and no data found on
		// chambers 7 or 8 then we can not find a track.
		if (sm5state == SM5Done) NoTrackFound();
		break;
	
	case WaitMoreChamber7:
		if (requests_completed == st5data.Count() and sm5state == SM5Done)
			ProcessClusters();
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in MansoTracker::EndOfClustersChamber7!");
	};
};


void MansoTracker::EndOfClustersChamber8()
{
	requests_completed++;  // Increment the number of requests completed for station 4.
	DebugMsg(4, "requests_completed = " << requests_completed );

	switch (sm4state)
	{
	case WaitChamber7:
		// Ignore. The requests for chamber 8 are already re-requested below.
		break;
		
	case WaitChamber8:
		{
		sm4state = WaitChamber7;
		st4z = z7;
		st4chamber = Chamber7;
	
		// We need to resend the requests for chamber 8, but change the request
		// to get data for chamber 7 instead:
		UInt reqlistsize = st5data.Count();
		DebugMsg(4, "Re-requesting clusters from chamber 7... reqlistsize = " << reqlistsize);

		Station5List::Iterator rec = st5data.First();
		for (UInt i = 0; i < reqlistsize; i++, rec++)
		{
			// Need to create a new st5 data block for the request.
			Station5Data* data = st5data.New();
			data->clusterpoint = rec->clusterpoint;
			data->tag.line = rec->tag.line;

			// Rebuild a region of interest for chamber 7.
			// Remember the parameters a and b are station specific.
			Point p7 = data->tag.line.FindIntersectWithXYPlain( z7 );
			data->tag.chamber = Chamber7;
			data->tag.roi.Create(p7, a7, b7);
			
			Float left, right, bottom, top;
			data->tag.roi.GetBoundaryBox(left, right, bottom, top);
			// Make request for chamber 7 data.
			RequestClusters(left, right, bottom, top, Chamber7, &data->tag);
		};
		}
		break;
	
	case WaitMoreChamber8:
		if (requests_completed == st5data.Count() and sm5state == SM5Done)
			ProcessClusters();
		break;
	
	default:
		DebugMsg(1, "ERROR: Unexpected state for SM4 in MansoTracker::EndOfClustersChamber8!");
	};
};


void MansoTracker::EndOfClustersChamber9()
{
	switch (sm5state)
	{
	case WaitChamber9:
		sm5state = SM5Done;
		EndOfClusterRequests();
		NoTrackFound();
		break;
		
	case WaitMoreChamber9:
		sm5state = SM5Done;
		EndOfClusterRequests();
		if (requests_completed == st5data.Count())
			ProcessClusters();
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in MansoTracker::EndOfClustersChamber9!");
	};
};


void MansoTracker::EndOfClustersChamber10()
{
	switch (sm5state)
	{
	case WaitChamber10:
		{
		sm5state = WaitChamber9;
		st5z = z9;
		
		// No clusters found on chamber 10 so we need to make a request for
		// clusters from chamber 9:
		Point p9 = mc1.line.FindIntersectWithXYPlain( z9 );

		// Build a region of interest for tracking station 5 (chamber 9).
		// Remember the parameters a and b are station specific.
		mc1.chamber = Chamber9;
		mc1.roi.Create(p9, a9, b9);

		Float left, right, bottom, top;
		mc1.roi.GetBoundaryBox(left, right, bottom, top);
		RequestClusters(left, right, bottom, top, Chamber9, &mc1);
		}
		break;

	case WaitMoreChamber10:
		sm5state = SM5Done;
		EndOfClusterRequests();
		if (requests_completed == st5data.Count())
			ProcessClusters();
		break;

	default:
		DebugMsg(1, "ERROR: Unexpected state for SM5 in MansoTracker::EndOfClustersChamber10!");
	};
};


void MansoTracker::ProjectToStation4(Station5Data* data, register Float station5z)
{
	// Perform chamber specific operations:
	// Since certain states of SM4 means that it is fetching for Chamber8
	// and other states are for fetching from Chamber7. We need to make
	// requests for the correct chamber.
	Assert( sm4state == WaitChamber8 or sm4state == WaitMoreChamber8 or
		sm4state == WaitChamber7 or sm4state == WaitMoreChamber7
	);
	TagData* tag = &data->tag;
	if (sm4state == WaitChamber8 or sm4state == WaitMoreChamber8)
	{
		// Form the vector line between trigger station 1 and tracking station 5,
		// and find the intersection point of the line with station 4 (chamber8).
		Line line51( Vertex(data->clusterpoint, station5z), v1 );
		Point intercept = line51.FindIntersectWithXYPlain(z8);
		tag->line = line51;
		
		// Build a region of interest for tracking station 4.
		tag->chamber = Chamber8;
		tag->roi.Create(intercept, a8, b8);
	}
	else
	{
		// Form the vector line between trigger station 1 and tracking station 5,
		// and find the intersection point of the line with station 4 (chamber7).
		Line line51( Vertex(data->clusterpoint, station5z), v1 );
		Point intercept = line51.FindIntersectWithXYPlain(z7);
		tag->line = line51;
		
		// Build a region of interest for tracking station 4.
		tag->chamber = Chamber7;
		tag->roi.Create(intercept, a7, b7);
	};

	// Make the request for clusters from station 4.
	Float left, right, bottom, top;
	tag->roi.GetBoundaryBox(left, right, bottom, top);
	RequestClusters(left, right, bottom, top, tag->chamber, tag);
};


void MansoTracker::ProcessClusters()
{
	DebugMsg(2, "ProcessClusters...");
	
	// Check if the cluster point list on station 4 is empty.
	// If it is then we have not found any tracks.
	foundpoint = st4points.First();
	if (foundpoint == st4points.End())
	{
		NoTrackFound();
		return;
	};
	
	st5rec = st5data.First();
	if (st5rec != st5data.End())
	{
		// Only look at station 5 data records that are for the found chamber number.
		// Note: either we only have chamber 8 data or we have chamber 7 data followed
		// by chamber 8 data.
		// Thus if we hit records that we are not interested in already then the list
		// contains no interesting data and we can signal no track found.
		if (st5rec->tag.chamber != st4chamber)
		{
			NoTrackFound();
			return;
		};
		 
		// For all combinations of cluster point pairs from station 4 and 5
		// signal a found track:
		do
		{
			DebugMsg(4, "\tst5rec->tag.chamber = " << st5rec->tag.chamber
				<< " , st4chamber = " << st4chamber
			);

			for (foundpoint = st4points.First(); foundpoint != st4points.End(); foundpoint++)
				FoundTrack();

			st5rec++;  // Get next station 5 cluster point.
		} while (st5rec != st5data.End() and st5rec->tag.chamber == st4chamber);
	}
	else
		NoTrackFound();
};


} // Tracking
} // dHLT
