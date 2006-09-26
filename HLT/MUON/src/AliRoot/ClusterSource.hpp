////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONClusterSource is used to extract cluster points from the data
   stored in a root file for a AliRoot simulated event.
   It is used by AliHLTMUONMicrodHLT as the input data set object.
 */

#ifndef ALIHLTMUONCLUSTERSOURCE_H
#define ALIHLTMUONCLUSTERSOURCE_H

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TClonesArray.h"
#include "AliRoot/Point.hpp"

class AliMUONDataInterface;


class AliHLTMUONClusterSource : public TObject
{
public:  // Unfortunately ROOT requires the following to be public.
	
	class AliBlockData : public TObject
	{
	public:
		AliBlockData();
		AliBlockData(Int_t chamber);
		virtual ~AliBlockData();

		Int_t& Chamber() { return fChamber; };
		TClonesArray& Clusters() { return fClusters; };

	private:

		Int_t fChamber;  // The chamber number this block of clusters came from.
		TClonesArray fClusters;  // The cluster points in this block.
		
		ClassDef(AliBlockData, 1);  // Data per block.
	};
	
	class AliEventData : public TObject
	{
	public:
		AliEventData();
		AliEventData(Int_t eventnumber);
		virtual ~AliEventData();

		Int_t& EventNumber() { return fEventNumber; };
		TClonesArray& Blocks() { return fBlocks; };

	private:

		Int_t fEventNumber;  // Event number in AliMUONDataInterface from which the clusters were taken.
		TClonesArray fBlocks;  // The list of cluster blocks for this event.
		
		ClassDef(AliEventData, 1);  // Data per event.
	};


public:

	enum AreaType
	{
		kFromWholePlane,
		kFromLeftHalfPlane,
		kFromRightHalfPlane
	};
	
	enum SourceType
	{
		kFromHits,
		kFromRawClusters
	};

	AliHLTMUONClusterSource();
	AliHLTMUONClusterSource(AliMUONDataInterface* data);
	
	virtual ~AliHLTMUONClusterSource();
	
	/* Get and set methods to specify how the FillFrom methods should fill the
	   internal data structures.
	 */
	void AreaToUse(AreaType value) { fAreaToUse = value; };
	AreaType AreaToUse() const { return fAreaToUse; };
	void DataToUse(SourceType value) { fDataToUse = value; };
	SourceType DataToUse() const { return fDataToUse; };
	void MaxBlockSize(UInt_t value) { fMaxBlockSize = value; };
	UInt_t MaxBlockSize() const { return fMaxBlockSize; };

	/* Fills the internal data structures from the specified data interface
	   for all the events found in AliMUONDataInterface.
	 */
	void FillFrom(AliMUONDataInterface* data);
	
	/* Fills the internal data structures from the specified data interface
	   for the given event.
	 */
	void FillFrom(AliMUONDataInterface* data, Int_t event);
	
	/* Fills the internal data structures from the specified data interface
	   for the given event and chamber.
	 */
	void FillFrom(AliMUONDataInterface* data, Int_t event, Int_t chamber);
	
	/* Fills the internal data structures from the specified data interface
	   for the given event, chamber and cluster number.
	   If 'newblock' is set to true then the new cluster point is added to 
	   a new block. Otherwise the point is added to the current block.
	   Note: This method ignores the fAreaToUse and fMaxBlockSize flags.
	   This is very usefull for custom cluster source filling.
	   For the case of adding data from AliMUONHit objects the 'cluster' parameter
	   becomes the track number in TreeH and not the index of the AliMUONRawCluster
	   object.
	 */
	void FillFrom(
			AliMUONDataInterface* data, 
			Int_t event, Int_t chamber, Int_t cluster,
			Bool_t newblock = kFALSE
		);
	
	/* Clears all the internal arrays.
	 */
	virtual void Clear(Option_t* option = "");
	
	// Get methods.
	TString FileName()   const { return fFilename; };
	TString FolderName() const { return fFoldername; };
	
	/* Returns the number of events stored.
	 */
	Int_t NumberOfEvents() const { return fEventList.GetEntriesFast(); };
	
	/* Fetches the specified event number stored in this AliHLTMUONClusterSource.
	   Sets the current block and cluster point to the first block and cluster
	   point respectively. If there are no blocks or clusters then these pointers
	   are set to NULL.
	   kTRUE is returned if the event was found, kFALSE otherwise.
	 */
	Bool_t GetEvent(Int_t eventnumber) const;
	
	/* Fetches the first event stored in this AliHLTMUONClusterSource.
	   Sets the current block and cluster point to the first block and cluster
	   point respectively. If there are no blocks or clusters then these pointers
	   are set to NULL.
	   kTRUE is returned if the event was found, kFALSE otherwise.
	 */
	Bool_t GetFirstEvent() const;
	
	/* Returns kTRUE if there are more events to interate over.
	 */
	Bool_t MoreEvents() const;
	
	/* Fetches the next event stored following the currently selected one.
	   kTRUE is returned if the event was found, kFALSE otherwise.
	 */
	Bool_t GetNextEvent() const;
	
	/* Returns the corresponding AliRoot event number for the current event.
	   -1 is returned if no event is selected.
	 */
	Int_t CurrentEvent() const;
	
	/* Returns the number of cluster blocks in the current event.
	   -1 is returned if no event is selected.
	 */
	Int_t NumberOfBlocks() const;
	
	/* Fetches the index'th block in the current event.
	   Sets the current cluster point to the first cluster point in the block.
	   If there are no cluster points then this pointer is set to NULL.
	   kTRUE is returned if the block was found, kFALSE otherwise.
	 */
	Bool_t GetBlock(Int_t index) const;
	
	/* Fetches the first block in the current event.
	   Sets the current cluster point to the first cluster point in the block.
	   If there are no cluster points then this pointer is set to NULL.
	   kTRUE is returned if the block was found, kFALSE otherwise.
	 */
	Bool_t GetFirstBlock() const;
	
	/* Returns kTRUE if there are more blocks to interate over.
	 */
	Bool_t MoreBlocks() const;
	
	/* Fetches the next block in the current event.
	   kTRUE is returned if the block was found, kFALSE otherwise.
	 */
	Bool_t GetNextBlock() const;
	
	/* Returns the current block intex number.
	   -1 is returned if no block is selected.
	 */
	Int_t CurrentBlock() const { return fBlockIndex; };
	
	/* Returns the chamber number of the current block.
	   -1 is returned if no block is selected.
	 */
	Int_t Chamber() const;
	
	/* Returns the number of cluster points in the current block.
	   -1 is returned if no block is selected.
	 */
	Int_t NumberOfClusters() const;
	
	/* Fetches the index'th cluster point in the current block.
	   kTRUE is returned if the point was found, kFALSE otherwise.
	 */
	const AliHLTMUONPoint* GetCluster(Int_t index) const;
	
	/* Fetches the first cluster point in the current block.
	   NULL is returned if the point was not found.
	 */
	const AliHLTMUONPoint* GetFirstCluster() const;
	
	/* Returns kTRUE if there are more cluster points to interate over.
	 */
	Bool_t MoreClusters() const;
	
	/* Fetches the next cluster point in the current block.
	   NULL is returned if the point was not found.
	 */
	const AliHLTMUONPoint* GetNextCluster() const;
	
	/* Returns the x and y coordinate of the current cluster point.
	   kFALSE is returned if there is no cluster point selected.
	 */
	Bool_t FetchCluster(Float_t& x, Float_t& y) const;
	
	/* Returns the current cluster point.
	   NULL is returned if no cluster point is selected.
	 */
	const AliHLTMUONPoint* GetCluster() const { return fCurrentCluster; };

private:

	/* Adds a new AliEventData block to the fEventList and updates the fCurrentEvent,
	   fCurrentBlock and fCurrentCluster pointers.
	 */ 
	void AddEvent(Int_t eventnumber);
	
	/* Adds a new block to the current event and updates fCurrentBlock and fCurrentCluster.
	   The chamber number is assigned to the blocks fChamber value.
	 */
	void AddBlock(Int_t chamber);
	
	/* Adds a new cluster point to the current event and block.
	   The fCurrentCluster is updated appropriately.
	 */
	void AddPoint(Float_t x, Float_t y);
	
	/* Checks if the file and folder names correspond to this AliHLTMUONClusterSource's 
	   file and folder names. kTRUE is returned if they do.
	   If the file and folder names are empty then they are assigned the names
	   as found in the data interface and kTRUE is returned.
	 */
	Bool_t FileAndFolderOk(AliMUONDataInterface* data);
	
	/* Adds the whole event from the data interface to the internal data structures.
	   It is assumed that FileAndFolderOk(data) returns true just before calling
	   this method.
	 */
	void AddEventFrom(AliMUONDataInterface* data, Int_t event);
	
	/* Adds all cluster points found on the given chamber from the specified data
	   interface. The data interface should be initialised correctly, that is the
	   event should already be selected before calling this method.
	 */
	void AddChamberFrom(AliMUONDataInterface* data, Int_t chamber);
	
	/* Adds the cluster point from the specified data interface.
	   The data interface should be initialised correctly, that is the event
	   should already be selected before calling this method.
	 */
	void AddClusterFrom(AliMUONDataInterface* data, Int_t chamber, Int_t cluster);
	
	/* Checks to see if the x and y coordinate of the cluster point are in the
	   chamber region we want to fill from.
	   kTRUE is returned if (x, y) is in the region, and kFALSE otherwise.
	 */
	Bool_t InFillRegion(Float_t x, Float_t y) const;

	/* Sets all the current pointers to NULL and indices to -1.
	 */
	void ResetAllPointers() const;
	
	/* Sets the block and trigger pointers to NULL and indices to -1.
	 */
	void ResetBlockPointers() const;
	
	/* Sets just the current cluster point pointer to NULL and index to -1.
	 */
	void ResetClusterPointers() const;
	

	// Dont allow copying.
	AliHLTMUONClusterSource(const AliHLTMUONClusterSource& /*object*/)
		: TObject(),
		  fAreaToUse(kFromWholePlane), fDataToUse(kFromRawClusters),
		  fMaxBlockSize(0xFFFFFFFF), fFilename(""), fFoldername(""),
		  fEventIndex(-1), fCurrentEvent(NULL), fBlockIndex(-1),
		  fCurrentBlock(NULL), fClusterIndex(-1), fCurrentCluster(NULL),
		  fEventList(AliHLTMUONClusterSource::AliEventData::Class())
	{}

	AliHLTMUONClusterSource& operator = (const AliHLTMUONClusterSource& /*object*/) { return *this; };

	AreaType fAreaToUse;    //! The part of the chamber to fill from.
	SourceType fDataToUse;  //! The type of raw AliRoot data to fill from.
	UInt_t fMaxBlockSize;   //! The maximum block size to create in the fill methods.

	TString fFilename;    // The file from which the cluster data was taken.
	TString fFoldername;  // The folder name from which cluster data was taken.
	
	mutable Int_t fEventIndex;             //! The index number of the currently selected event.
	mutable AliEventData* fCurrentEvent;      //! Pointer to the currently selected event.
	mutable Int_t fBlockIndex;             //! The index number of the currently selected block.
	mutable AliBlockData* fCurrentBlock;      //! Pointer to the currently selected block.
	mutable Int_t fClusterIndex;           //! The index number of the currently selected cluster point.
	mutable AliHLTMUONPoint* fCurrentCluster;        //! Pointer to the currently selected cluster point.

	TClonesArray fEventList;  // List of clusters per event.

	ClassDef(AliHLTMUONClusterSource, 1)  // The source of cluster point blocks for dHLT.
};


#endif // ALIHLTMUONCLUSTERSOURCE_H
