////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_TRACKING_MANSO_TRACKER_HPP
#define dHLT_TRACKING_MANSO_TRACKER_HPP

#include "BasicTypes.hpp"
#include "Tracking/Tracker.hpp"
#include "Buffers/List.hpp"
#include "Buffers/CountedList.hpp"


namespace dHLT
{
namespace Tracking
{


class MansoTracker : public Tracker
{
public:

	MansoTracker();
	virtual ~MansoTracker() {};


	/* Note: Reset should be called for before calling FindTrack, for the
	   second or subsequent method calls to FindTrack.
	 */
	virtual void FindTrack(const TriggerRecord& trigger);
	
	virtual void ReturnClusters(void* tag, const ClusterPoint* clusters, const UInt count);
	virtual void EndOfClusters(void* tag);
	virtual void FillTrackData(Track& track);
	virtual void Reset();


	/* Get and set methods for the a and b parameters used to build the region
	   of interests. Refer to MansoFilter() for details about a and b parameters.
	 */
	static Float GetA7()                  { return a7; };
	static void SetA7(const Float value)  { a7 = a7; };
	static Float GetA8()                  { return a8; };
	static void SetA8(const Float value)  { a8 = value; };
	static Float GetA9()                  { return a9; };
	static void SetA9(const Float value)  { a9 = value; };
	static Float GetA10()                 { return a10; };
	static void SetA10(const Float value) { a10 = value; };

	static Float GetB7()                  { return b7; };
	static void SetB7(const Float value)  { b7 = value; };
	static Float GetB8()                  { return b8; };
	static void SetB8(const Float value)  { b8 = value; };
	static Float GetB9()                  { return b9; };
	static void SetB9(const Float value)  { b9 = value; };
	static Float GetB10()                 { return b10; };
	static void SetB10(const Float value) { b10 = value; };
	
	static Float GetZ7()                  { return z7; };
	static void SetZ7(const Float value)  { z7 = value; };
	static Float GetZ8()                  { return z8; };
	static void SetZ8(const Float value)  { z8 = value; };
	static Float GetZ9()                  { return z9; };
	static void SetZ9(const Float value)  { z9 = value; };
	static Float GetZ10()                 { return z10; };
	static void SetZ10(const Float value) { z10 = value; };
	static Float GetZ11()                 { return z11; };
	static void SetZ11(const Float value) { z11 = value; };
	static Float GetZ13()                 { return z13; };
	static void SetZ13(const Float value) { z13 = value; };


protected:


	class RegionOfInterest
	{
	public:
		
		RegionOfInterest() {};

		RegionOfInterest(const Point p, const Float a, const Float b)
		{
			Create(p, a, b);
		};

		/* Creates a region of interest. In this implementation it is a
		   circular disk.

		   The point p is the intersecting point of the track with the chamber.
		   For more details and for details about the parameters a and b refer to:
		   "A first algorithm for dimuon High Level Trigger"
		   Ref ID:  ALICE-INT-2002-04 version 1.0
		   equation:
		     Rs = a * Rp + b
		   given on page 3 section 4.
		 */
		void Create(const Point p, const Float a, const Float b);

		/* Returns true if the point p is within the region of interest.
		 */
		bool Contains(const Point p) const;

		void GetBoundaryBox(Float& left, Float& right, Float& bottom, Float& top);

	private:

		Point centre;  // The centre point of the region of interest.
		Float Rs;      // The redius of the region of interest around fcentre.
	};


	class Vertex
	{
	public:

		Float x, y, z;

		Vertex(const Float x = 0.0, const Float y = 0.0, const Float z = 0.0);
		Vertex(const Point xy, const Float z);

		Point AsXYPoint() const
		{
			return Point(x, y);
		};
	};

	
	class Line
	{
	public:

		/* Creates a vector line between points A and B.
		   Ax, Ay, Az are x, y and z coordinates for space point A respectively.
		   simmilarly for B.
		 */
		Line(
			const Float Ax = 0.0, const Float Ay = 0.0, const Float Az = 0.0,
			const Float Bx = 0.0, const Float By = 0.0, const Float Bz = 0.0
		);

		/* Creates a vector line between vertices A and B.
		 */
		Line(const Vertex A, const Vertex B);

		/* Finds the intersection point with the xy plain specified by the z coordinate.
		   The z coordiante would be the distance of the n'th chamber to the interaction
		   vertex.
		 */
		Point FindIntersectWithXYPlain(const Float z) const;

	private:

		// Parameters for the vector line:  L = M*t + C
		Float Mx, My, Mz, Cx, Cy, Cz;
	};

	
	struct TagData
	{
		ChamberID chamber;     // The chamber on which the region of interest lies.
		RegionOfInterest roi;  // Region of interest on the next station.
		Line line;             // line between a cluster point and the previous station.
	};
	
	struct Station5Data
	{
		ClusterPoint clusterpoint;  // Cluster point found on station 5.
		TagData tag;
	};
	
	typedef Buffers::CountedList<Station5Data> Station5List;
	typedef Buffers::List<ClusterPoint> Station4List;
	
	
	void ReceiveClustersChamber7(const ClusterPoint* clusters, const UInt count, const TagData* data);
	void ReceiveClustersChamber8(const ClusterPoint* clusters, const UInt count, const TagData* data);
	void ReceiveClustersChamber9(const ClusterPoint* clusters, const UInt count);
	void ReceiveClustersChamber10(const ClusterPoint* clusters, const UInt count);
	void EndOfClustersChamber7();
	void EndOfClustersChamber8();
	void EndOfClustersChamber9();
	void EndOfClustersChamber10();

	void ProjectToStation4(Station5Data* data, register Float station5z);
	void ProcessClusters();

DebugCode(
public:
)
	// States for state machine 4 (SM4).
	enum StatesSM4
	{
		SM4Idle,
		WaitChamber8,
		WaitMoreChamber8,
		WaitChamber7,
		WaitMoreChamber7
	};
	
	// States for state machine 5 (SM5).
	enum StatesSM5
	{
		SM5Idle,
		WaitChamber10,
		WaitMoreChamber10,
		WaitChamber9,
		WaitMoreChamber9,
		SM5Done
	};
	
protected:
	
	StatesSM4 sm4state;  // State of SM4 used for fetching clusters on chambers 7 and 8.
	StatesSM5 sm5state;  // State of SM5 used for fetching clusters on chambers 9 and 10.
	UInt requests_completed;  // Number of requests for station 4 that have completed.
	ChamberID st4chamber;     // The chamber on station 4 that data was retreived from.
	
	Vertex v1;    // The impact (hit) vertex for trigger station 1.
	TagData mc1;  // Trigger station 1 data.

	Float st5z;   // The z coordinate to use for station 5.
	Station5List st5data;  // List of found cluster points for station 5 and their tag data.
	Float st4z;   // The z coordinate to use for station 4.
	Station4List st4points;  // The found cluster points for station 4.

	// Iterators used in the FoundTrack, FillTrackData methods.
	Station5List::Iterator st5rec;
	Station4List::Iterator foundpoint; 

private:

	static Float a7, b7;    // Parameters used to create a region of interest for the 7'th chamber.
	static Float a8, b8;    // Parameters used to create a region of interest for the 8'th chamber.
	static Float a9, b9;    // Parameters used to create a region of interest for the 9'th chamber.
	static Float a10, b10;  // Parameters used to create a region of interest for the 10'th chamber.
	static Float z7, z8, z9, z10, z11, z13;  // Z coordinates of chambers 7 to 10.

};


} // Tracking
} // dHLT

#endif // dHLT_TRACKING_MANSO_TRACKER_HPP
