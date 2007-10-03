#ifndef ALIHLTMUONMANSOTRACKERFSM_H
#define ALIHLTMUONMANSOTRACKERFSM_H
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
 *  @file   AliHLTMUONMansoTrackerFSM.h
 *  @author Artur Szostak <artursz@iafrica.com>
 *  @date   
 *  @brief  The AliHLTMUONMansoTrackerFSM implements the Manso tracking
 *          algorithm as a finite state machine, which partially reconstructs
 *          tracks in the muon spectrometer.
 */

#include "AliHLTMUONDataTypes.h"
#include "AliHLTMUONList.h"
#include "AliHLTMUONCountedList.h"
#include "AliHLTMUONRecHitsBlockStruct.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONMansoTracksBlockStruct.h"
#include "AliHLTMUONMansoTrackerFSMCallback.h"
#include <cassert>


class AliHLTMUONMansoTrackerFSM
{
public:

	AliHLTMUONMansoTrackerFSM();
	virtual ~AliHLTMUONMansoTrackerFSM() {}


	/* This is the starting point for the tracking algorithm. The tracker is 
	   called at this point with the specified trigger record. It needs to figure
	   out which cluster blocks it needs and request them with calls to
	   RequestClusters.
	   Any memory allocated at this point should be released in the Reset method.
	   
	   Note: Reset should be called for before calling FindTrack, for the
	   second or subsequent method calls to FindTrack.
	 */
	virtual void FindTrack(const AliHLTMUONTriggerRecordStruct& trigger);
	
	/* When requested clusters have been found by the framework they are returned
	   to the tracker using this method.
	   This method should implement any processing of the cluster blocks. If more
	   more regions of interest are identified then appropriate request should me 
	   made using RequestClusters. The tag parameter will be the same one as was
	   passed to RequestClusters.
	 */
	virtual void ReturnClusters(
			void* tag, const AliHLTMUONRecHitStruct* clusters,
			AliHLTUInt32_t count
		);
	
	/* When no more clusters are to be expected for the request with the corresponding
	   tag value, then this method is called.
	   Any final processing can be placed in here and when the track is found then
	   the algorithm can call FoundTrack otherwise NoTrackFound to indicate end of 
	   processing.
	 */
	virtual void EndOfClusters(void* tag);
	
	/* Called to receive track information after receiving a FoundTrack call.
	   The tracker should fill the track data block with all relevant information.
	   The method will return false if the momentum could not be calculated from
	   the hits found. The hits in the track structure will however still be filled.
	 */
	virtual bool FillTrackData(AliHLTMUONMansoTrackStruct& track);
	
	/* Called when the tracker should be reset to a initial state. 
	   All extra internal allocated data structured should be released.
	 */
	virtual void Reset();
	
	/* To set the TrackerCallback callback object.
	 */
	inline void SetCallback(AliHLTMUONMansoTrackerFSMCallback* callback)
	{
		fCallback = callback;
	};


	/* Get and set methods for the a and b parameters used to build the region
	   of interests. Refer to AliRegionOfInterest for details about a and b parameters.
	 */
	static AliHLTFloat32_t GetA7()            { return fgA7; };
	static void SetA7(AliHLTFloat32_t value)  { fgA7 = value; };
	static AliHLTFloat32_t GetA8()            { return fgA8; };
	static void SetA8(AliHLTFloat32_t value)  { fgA8 = value; };
	static AliHLTFloat32_t GetA9()            { return fgA9; };
	static void SetA9(AliHLTFloat32_t value)  { fgA9 = value; };
	static AliHLTFloat32_t GetA10()           { return fgA10; };
	static void SetA10(AliHLTFloat32_t value) { fgA10 = value; };

	static AliHLTFloat32_t GetB7()            { return fgB7; };
	static void SetB7(AliHLTFloat32_t value)  { fgB7 = value; };
	static AliHLTFloat32_t GetB8()            { return fgB8; };
	static void SetB8(AliHLTFloat32_t value)  { fgB8 = value; };
	static AliHLTFloat32_t GetB9()            { return fgB9; };
	static void SetB9(AliHLTFloat32_t value)  { fgB9 = value; };
	static AliHLTFloat32_t GetB10()           { return fgB10; };
	static void SetB10(AliHLTFloat32_t value) { fgB10 = value; };
	
	static AliHLTFloat32_t GetZ7()            { return fgZ7; };
	static void SetZ7(AliHLTFloat32_t value)  { fgZ7 = value; };
	static AliHLTFloat32_t GetZ8()            { return fgZ8; };
	static void SetZ8(AliHLTFloat32_t value)  { fgZ8 = value; };
	static AliHLTFloat32_t GetZ9()            { return fgZ9; };
	static void SetZ9(AliHLTFloat32_t value)  { fgZ9 = value; };
	static AliHLTFloat32_t GetZ10()           { return fgZ10; };
	static void SetZ10(AliHLTFloat32_t value) { fgZ10 = value; };
	static AliHLTFloat32_t GetZ11()           { return fgZ11; };
	static void SetZ11(AliHLTFloat32_t value) { fgZ11 = value; };
	static AliHLTFloat32_t GetZ13()           { return fgZ13; };
	static void SetZ13(AliHLTFloat32_t value) { fgZ13 = value; };


protected:

	class AliRegionOfInterest
	{
	public:
		
		AliRegionOfInterest() : fCentre(), fRs(0.0) {};

		AliRegionOfInterest(AliHLTMUONRecHitStruct p, AliHLTFloat32_t a, AliHLTFloat32_t b)
			: fCentre(), fRs(0)
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
		void Create(AliHLTMUONRecHitStruct p, AliHLTFloat32_t a, AliHLTFloat32_t b);

		/* Returns true if the point p is within the region of interest.
		 */
		bool Contains(AliHLTMUONRecHitStruct p) const;

		void GetBoundaryBox(AliHLTFloat32_t& left, AliHLTFloat32_t& right, AliHLTFloat32_t& bottom, AliHLTFloat32_t& top) const;

	private:

		AliHLTMUONRecHitStruct fCentre;  // The centre point of the region of interest.
		AliHLTFloat32_t fRs;      // The redius of the region of interest around fcentre.
	};


	class AliVertex
	{
	public:

		AliVertex(AliHLTFloat32_t x = 0.0, AliHLTFloat32_t y = 0.0, AliHLTFloat32_t z = 0.0);
		AliVertex(AliHLTMUONRecHitStruct xy, AliHLTFloat32_t z);

		AliHLTMUONRecHitStruct AsXYPoint() const
		{
			AliHLTMUONRecHitStruct p;
			p.fX = fX;
			p.fY = fY;
			p.fZ = 0;
			return p;
		};

		// Get/set methods:
		AliHLTFloat32_t X() const { return fX; };
		AliHLTFloat32_t Y() const { return fY; };
		AliHLTFloat32_t Z() const { return fZ; };
		AliHLTFloat32_t& X() { return fX; };
		AliHLTFloat32_t& Y() { return fY; };
		AliHLTFloat32_t& Z() { return fZ; };
		void X(AliHLTFloat32_t value) { fX = value; };
		void Y(AliHLTFloat32_t value) { fY = value; };
		void Z(AliHLTFloat32_t value) { fZ = value; };

	private:

		AliHLTFloat32_t fX, fY, fZ; // 3D coordinates.
	};

	
	class AliLine
	{
	public:

		/* Creates a vector line between points A and B.
		   ax, ay, az are x, y and z coordinates for space point A respectively.
		   simmilarly for B.
		 */
		AliLine(
			AliHLTFloat32_t ax = 0.0, AliHLTFloat32_t ay = 0.0, AliHLTFloat32_t az = 0.0,
			AliHLTFloat32_t bx = 0.0, AliHLTFloat32_t by = 0.0, AliHLTFloat32_t bz = 0.0
		);

		/* Creates a vector line between vertices A and B.
		 */
		AliLine(AliVertex a, AliVertex b);

		/* Finds the intersection point with the xy plain specified by the z coordinate.
		   The z coordiante would be the distance of the n'th chamber to the interaction
		   vertex.
		 */
		AliHLTMUONRecHitStruct FindIntersectWithXYPlain(AliHLTFloat32_t z) const;

	private:

		// Parameters for the vector line:  L = M*t + C
		AliHLTFloat32_t fMx, fMy, fMz, fCx, fCy, fCz;  // line parameters.
	};

	
	struct AliTagData
	{
		AliHLTMUONChamberName fChamber;     // The chamber on which the region of interest lies.
		AliRegionOfInterest fRoi;  // Region of interest on the next station.
		AliLine fLine;             // line between a cluster point and the previous station.

		AliTagData() : fChamber(kChamber1), fRoi(), fLine() {};
	};
	
	struct AliStation5Data
	{
		AliHLTMUONRecHitStruct fClusterPoint;  // Cluster point found on station 5.
		AliTagData fTag;  // Chamber, ROI and line data for station 5.

		AliStation5Data() : fClusterPoint(), fTag() {};
	};
	
	typedef AliHLTMUONCountedList<AliStation5Data> Station5List;

	struct AliStation4Data
	{
		AliHLTMUONRecHitStruct fClusterPoint;  // Cluster point found on station 4.
		const AliTagData* fSt5tag;      // Corresponding station 5 tag.

		AliStation4Data() : fClusterPoint(), fSt5tag() {};

		AliStation4Data(const AliStation4Data& data) :
			fClusterPoint(data.fClusterPoint), fSt5tag(data.fSt5tag)
		{};

		AliStation4Data& operator = (const AliStation4Data& data)
		{
			fClusterPoint = data.fClusterPoint;
			fSt5tag = data.fSt5tag;
			return *this;
		};
	};

	typedef AliHLTMUONList<AliStation4Data> Station4List;
	
	
	void ReceiveClustersChamber7(const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count, const AliTagData* data);
	void ReceiveClustersChamber8(const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count, const AliTagData* data);
	void ReceiveClustersChamber9(const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count);
	void ReceiveClustersChamber10(const AliHLTMUONRecHitStruct* clusters, AliHLTUInt32_t count);
	void EndOfClustersChamber7();
	void EndOfClustersChamber8();
	void EndOfClustersChamber9();
	void EndOfClustersChamber10();

	void ProjectToStation4(AliStation5Data* data, register AliHLTFloat32_t station5z);
	void ProcessClusters();

#ifdef DEBUG
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

	AliHLTMUONMansoTrackerFSMCallback* fCallback;
	
	StatesSM4 fSm4state;  // State of SM4 used for fetching clusters on chambers 7 and 8.
	StatesSM5 fSm5state;  // State of SM5 used for fetching clusters on chambers 9 and 10.
	AliHLTUInt32_t fRequestsCompleted;  // Number of requests for station 4 that have completed.
	AliHLTMUONChamberName fSt4chamber;  // The chamber on station 4 that data was retreived from.
	
	AliVertex fV1;    // The impact (hit) vertex for trigger station 1.
	AliTagData fMc1;  // Trigger station 1 data.

	Station5List fSt5data;  // List of found cluster points for station 5 and their tag data.
	Station4List fSt4points;  // The found cluster points for station 4.

	// Iterators used in the FoundTrack, FillTrackData methods.
	Station5List::Iterator fSt5rec;      // current station 5 record
	Station4List::Iterator fFoundPoint;  // current found point
	AliHLTInt32_t fTriggerId;  // The current ID number of the trigger record being processed.
	AliHLTInt32_t fTrackId;   // Track ID counter for the current track.
	
	
	/* To request clusters from the boundary box specified by the 'left', 'right',
	   'top' and 'bottom' boundaries and on the given chamber use this method call.
	   Supply a tag parameter if you want the request uniquely identified. 
	   This is usefull to supply a pointer to some internal state data structure
	   to figure out where processing should continue in the ReturnClusters or
	   EndOfClusters methods.
	 */
	inline void RequestClusters(
			AliHLTFloat32_t left, AliHLTFloat32_t right, AliHLTFloat32_t bottom, AliHLTFloat32_t top,
			AliHLTMUONChamberName chamber, const void* tag = NULL
		)
	{
		assert( fCallback != NULL );
		fCallback->RequestClusters(this, left, right, bottom, top, chamber, tag);
	};

	/* When no more cluster requests will be generated by this tracker then this
	   method should be called.
	   DO NOT request more clusters after calling this method.
	 */
	inline void EndOfClusterRequests()
	{
		assert( fCallback != NULL );
		fCallback->EndOfClusterRequests(this);
	};

	/* When the tracker has found a track it should call this method to inform
	   the rest of the system. At this point all cluster blocks received with
	   ReturnClusters are to be considered released and MUST NOT be accessed.
	 */
	inline void FoundTrack()
	{
		assert( fCallback != NULL );
		fCallback->FoundTrack(this);
	};

	/* If the tracker is finished processing the trigger record but has not found 
	   a track it should call this method to inform the rest of the system.
	   At this point all cluster blocks received with ReturnClusters are to be
	   considered released and MUST NOT be accessed.
	 */
	inline void NoTrackFound()
	{
		assert( fCallback != NULL );
		fCallback->NoTrackFound(this);
	};

private:

	// Not allowed to copy this object.
	AliHLTMUONMansoTrackerFSM(const AliHLTMUONMansoTrackerFSM& tracker);
	AliHLTMUONMansoTrackerFSM& operator = (const AliHLTMUONMansoTrackerFSM& tracker);

	static AliHLTFloat32_t fgA7, fgB7;    // Parameters used to create a region of interest for the 7'th chamber.
	static AliHLTFloat32_t fgA8, fgB8;    // Parameters used to create a region of interest for the 8'th chamber.
	static AliHLTFloat32_t fgA9, fgB9;    // Parameters used to create a region of interest for the 9'th chamber.
	static AliHLTFloat32_t fgA10, fgB10;  // Parameters used to create a region of interest for the 10'th chamber.
	static AliHLTFloat32_t fgZ7, fgZ8, fgZ9, fgZ10, fgZ11, fgZ13;  // Z coordinates of chambers 7 to 10.

};


#endif // ALIHLTMUONMANSOTRACKERFSM_H
