////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONCOREMANSOTRACKER_H
#define ALIHLTMUONCOREMANSOTRACKER_H

#include "BasicTypes.hpp"
#include "Tracking/Tracker.hpp"
#include "Buffers/List.hpp"
#include "Buffers/CountedList.hpp"


class AliHLTMUONCoreMansoTracker : public AliHLTMUONCoreTracker
{
public:

	AliHLTMUONCoreMansoTracker();
	virtual ~AliHLTMUONCoreMansoTracker() {};


	/* Note: Reset should be called for before calling FindTrack, for the
	   second or subsequent method calls to FindTrack.
	 */
	virtual void FindTrack(const AliHLTMUONCoreTriggerRecord& trigger);
	
	virtual void ReturnClusters(void* tag, const AliHLTMUONCoreClusterPoint* clusters, UInt count);
	virtual void EndOfClusters(void* tag);
	virtual void FillTrackData(AliHLTMUONCoreTrack& track);
	virtual void Reset();


	/* Get and set methods for the a and b parameters used to build the region
	   of interests. Refer to MansoFilter() for details about a and b parameters.
	 */
	static Float GetA7()            { return fgA7; };
	static void SetA7(Float value)  { fgA7 = value; };
	static Float GetA8()            { return fgA8; };
	static void SetA8(Float value)  { fgA8 = value; };
	static Float GetA9()            { return fgA9; };
	static void SetA9(Float value)  { fgA9 = value; };
	static Float GetA10()           { return fgA10; };
	static void SetA10(Float value) { fgA10 = value; };

	static Float GetB7()            { return fgB7; };
	static void SetB7(Float value)  { fgB7 = value; };
	static Float GetB8()            { return fgB8; };
	static void SetB8(Float value)  { fgB8 = value; };
	static Float GetB9()            { return fgB9; };
	static void SetB9(Float value)  { fgB9 = value; };
	static Float GetB10()           { return fgB10; };
	static void SetB10(Float value) { fgB10 = value; };
	
	static Float GetZ7()            { return fgZ7; };
	static void SetZ7(Float value)  { fgZ7 = value; };
	static Float GetZ8()            { return fgZ8; };
	static void SetZ8(Float value)  { fgZ8 = value; };
	static Float GetZ9()            { return fgZ9; };
	static void SetZ9(Float value)  { fgZ9 = value; };
	static Float GetZ10()           { return fgZ10; };
	static void SetZ10(Float value) { fgZ10 = value; };
	static Float GetZ11()           { return fgZ11; };
	static void SetZ11(Float value) { fgZ11 = value; };
	static Float GetZ13()           { return fgZ13; };
	static void SetZ13(Float value) { fgZ13 = value; };


protected:

	class RegionOfInterest
	{
	public:
		
		RegionOfInterest() {};

		RegionOfInterest(AliHLTMUONCorePoint p, Float a, Float b)
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
		void Create(AliHLTMUONCorePoint p, Float a, Float b);

		/* Returns true if the point p is within the region of interest.
		 */
		bool Contains(AliHLTMUONCorePoint p) const;

		void GetBoundaryBox(Float& left, Float& right, Float& bottom, Float& top);

	private:

		AliHLTMUONCorePoint fCentre;  // The centre point of the region of interest.
		Float fRs;      // The redius of the region of interest around fcentre.
	};


	class Vertex
	{
	public:

		Float fX, fY, fZ;

		Vertex(Float x = 0.0, Float y = 0.0, Float z = 0.0);
		Vertex(AliHLTMUONCorePoint xy, Float z);

		AliHLTMUONCorePoint AsXYPoint() const
		{
			return AliHLTMUONCorePoint(fX, fY);
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
			Float Ax = 0.0, Float Ay = 0.0, Float Az = 0.0,
			Float Bx = 0.0, Float By = 0.0, Float Bz = 0.0
		);

		/* Creates a vector line between vertices A and B.
		 */
		Line(Vertex A, Vertex B);

		/* Finds the intersection point with the xy plain specified by the z coordinate.
		   The z coordiante would be the distance of the n'th chamber to the interaction
		   vertex.
		 */
		AliHLTMUONCorePoint FindIntersectWithXYPlain(Float z) const;

	private:

		// Parameters for the vector line:  L = M*t + C
		Float fMx, fMy, fMz, fCx, fCy, fCz;
	};

	
	struct TagData
	{
		AliHLTMUONCoreChamberID fChamber;     // The chamber on which the region of interest lies.
		RegionOfInterest fRoi;  // Region of interest on the next station.
		Line fLine;             // line between a cluster point and the previous station.
	};
	
	struct Station5Data
	{
		AliHLTMUONCoreClusterPoint fClusterPoint;  // Cluster point found on station 5.
		TagData fTag;
	};
	
	typedef AliHLTMUONCoreCountedList<Station5Data> Station5List;

	struct Station4Data
	{
		AliHLTMUONCoreClusterPoint fClusterPoint;  // Cluster point found on station 4.
		const TagData* fSt5tag;      // Corresponding station 5 tag.
	};

	typedef AliHLTMUONCoreList<Station4Data> Station4List;
	
	
	void ReceiveClustersChamber7(const AliHLTMUONCoreClusterPoint* clusters, UInt count, const TagData* data);
	void ReceiveClustersChamber8(const AliHLTMUONCoreClusterPoint* clusters, UInt count, const TagData* data);
	void ReceiveClustersChamber9(const AliHLTMUONCoreClusterPoint* clusters, UInt count);
	void ReceiveClustersChamber10(const AliHLTMUONCoreClusterPoint* clusters, UInt count);
	void EndOfClustersChamber7();
	void EndOfClustersChamber8();
	void EndOfClustersChamber9();
	void EndOfClustersChamber10();

	void ProjectToStation4(Station5Data* data, register Float station5z);
	void ProcessClusters();

#if defined(DEBUG) || (defined(USE_ALILOG) && ! defined(LOG_NO_DEBUG))
public:
#endif
	// States for state machine 4 (SM4).
	enum StatesSM4
	{
		kSM4Idle,
		kWaitChamber8,
		kWaitMoreChamber8,
		kWaitChamber7,
		kWaitMoreChamber7
	};
	
	// States for state machine 5 (SM5).
	enum StatesSM5
	{
		kSM5Idle,
		kWaitChamber10,
		kWaitMoreChamber10,
		kWaitChamber9,
		kWaitMoreChamber9,
		kSM5Done
	};
	
protected:
	
	StatesSM4 fSm4state;  // State of SM4 used for fetching clusters on chambers 7 and 8.
	StatesSM5 fSm5state;  // State of SM5 used for fetching clusters on chambers 9 and 10.
	UInt fRequestsCompleted;  // Number of requests for station 4 that have completed.
	AliHLTMUONCoreChamberID fSt4chamber;     // The chamber on station 4 that data was retreived from.
	
	Vertex fV1;    // The impact (hit) vertex for trigger station 1.
	TagData fMc1;  // Trigger station 1 data.

	Float fSt5z;   // The z coordinate to use for station 5.
	Station5List fSt5data;  // List of found cluster points for station 5 and their tag data.
	Float fSt4z;   // The z coordinate to use for station 4.
	Station4List fSt4points;  // The found cluster points for station 4.

	// Iterators used in the FoundTrack, FillTrackData methods.
	Station5List::Iterator fSt5rec;
	Station4List::Iterator fFoundPoint;

private:

	static Float fgA7, fgB7;    // Parameters used to create a region of interest for the 7'th chamber.
	static Float fgA8, fgB8;    // Parameters used to create a region of interest for the 8'th chamber.
	static Float fgA9, fgB9;    // Parameters used to create a region of interest for the 9'th chamber.
	static Float fgA10, fgB10;  // Parameters used to create a region of interest for the 10'th chamber.
	static Float fgZ7, fgZ8, fgZ9, fgZ10, fgZ11, fgZ13;  // Z coordinates of chambers 7 to 10.

};


#endif // ALIHLTMUONCOREMANSOTRACKER_H

