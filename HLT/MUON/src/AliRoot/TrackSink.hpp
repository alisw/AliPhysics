////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONTrackSink is used as the output target object by AliHLTMUONMicrodHLT.
   It is just a fancy buffer to store tracks that are found by the dHLT
   algorithm.
 */

#ifndef ALIHLTMUONTRACKSINK_H
#define ALIHLTMUONTRACKSINK_H

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TClonesArray.h"
#include "AliRoot/Track.hpp"
#include "AliRoot/TriggerSource.hpp"

#include <vector>


class AliHLTMUONTrackSink : public TObject
{
public:

	AliHLTMUONTrackSink();
	virtual ~AliHLTMUONTrackSink();
	
	/* Adds a new AliEventData block to the fEventList and updates internal pointers.
	   Cannot have duplicate event numbers so this method will display an error
	   message if one attempts to add the same event number more than once.
	 */ 
	void AddEvent(Int_t eventnumber);
	
	/* Adds a new block to the current event and updates fCurrentBlock and
	   fCurrentTrack.
	 */
	void AddBlock();
	
	/* Adds a new track to the current event and block, and returns a pointer
	   to this track object to be filled by the caller.
	   The fCurrentTrack is updated appropriately.
	   If no current block is selected then NULL is returned.
	 */
	AliHLTMUONTrack* AddTrack();
	
	/* Adds the given track to the current block.
	   If no current block is selected then an error message is displayed.
	 */
	void AddTrack(const AliHLTMUONTrack& track);
	
	/* Adds the specified track parameters as a new track.
	   The fCurrentTrack is updated appropriately.
	 */
	void AddTrack(
			Int_t triggerid, Int_t sign, Float_t momentum,
			Float_t pt, const AliHLTMUONPoint hits[10], const AliHLTMUONRegion regions[10]
		);
	
	/* Sets the internal file and folder names from the trigger source.
	 */
	void SetNames(const AliHLTMUONTriggerSource* triggersource);
	
	/* Clears all the internal arrays.
	 */
	virtual void Clear(Option_t* option = "");
	
	// Get methods for file and folder names.
	TString FileName()   const { return fFilename; }
	TString FolderName() const { return fFoldername; }
	
	/* Returns the number of events stored.
	 */
	Int_t NumberOfEvents() const { return fEventList.GetEntriesFast(); };
	
	/* Fetches the event data corresponding to the given event number.
	   Sets the current block and track to the first block and track for 
	   the newly selected event.
	   If there are no blocks or tracks then these pointers are set to NULL.
	   kTRUE is returned if the event was found. kFALSE is returned if the
	   event was not found and the internal pointers left untouched.
	 */
	Bool_t GetEvent(Int_t eventnumber) const;
	
	/* Fetches the first event stored in this AliHLTMUONTrackSink.
	   Sets the current block and track to the first block and track of the
	   first event.
	   If there are no blocks or tracks then these pointers are set to NULL.
	   kFALSE is returned if there are no events stored, kTRUE is returned
	   on success.
	 */
	Bool_t GetFirstEvent() const;
	
	/* Returns kTRUE if there are more events to iterate over,
	   kFALSE otherwise.
	 */
	Bool_t MoreEvents() const;
	
	/* Fetches the next event stored following the currently selected one.
	   Sets the current block and track pointers to the first block and track
	   in the newly selected event. These pointers are set to NULL if there
	   are no blocks or tracks for this event.
	 */
	Bool_t GetNextEvent() const;
	
	/* Returns the corresponding AliRoot event number for the current event.
	   -1 is returned if no event is selected.
	 */
	Int_t CurrentEvent() const;
	
	/* Returns the number of track blocks in the current event.
	   -1 is returned if no event is selected.
	 */
	Int_t NumberOfBlocks() const;
	
	/* Fetches the index'th block in the current event.
	   Sets the current track to the first track in the block.
	   If there are no tracks then the track pointers are reset.
	   kTRUE is returned if the block was found, kFALSE otherwise.
	 */
	Bool_t GetBlock(Int_t index) const;
	
	/* Fetches the first block in the current event.
	   Sets the current track to the first track in the block.
	   If there are no tracks then the fCurrentTracks pointer is set to NULL.
	 */
	Bool_t GetFirstBlock() const;
	
	/* Returns kTRUE if there are more blocks to iterate over.
	 */
	Bool_t MoreBlocks() const;
	
	/* Fetches the next block in the current event.
	   kTRUE is returned if the block was found, kFALSE otherwise.
	   The current track pointers are reset if no more blocks are found.
	 */
	Bool_t GetNextBlock() const;
	
	/* Returns the current block index number.
	   -1 is returned if no block was selected.
	 */
	Int_t CurrentBlock() const { return fBlockIndex; };
	
	/* Returns the number of tracks in the current block.
	   -1 is returned if no block is selected.
	 */
	Int_t NumberOfTracks() const;
	
	/* Fetches the index'th track in the current block.
	   NULL is returned if the track was not found.
	 */
	const AliHLTMUONTrack* GetTrack(Int_t index) const;
	
	/* Fetches the first track in the current block.
	   NULL is returned if the track was not found.
	 */
	const AliHLTMUONTrack* GetFirstTrack() const;
	
	/* Returns kTRUE if there are more tracks to iterate over in
	   the current block. kFALSE is returned otherwise.
	 */
	Bool_t MoreTracks() const;
	
	/* Fetches the next track in the current block.
	   NULL is returned if the track was not found.
	 */
	const AliHLTMUONTrack* GetNextTrack() const;
	
	/* Returns the currently selected track.
	   NULL is returned if there is no track selected.
	 */
	const AliHLTMUONTrack* GetTrack() const { return fCurrentTrack; };
	
	/* Returns the currently selected track index.
	   -1 is returned if no track was selected.
	 */
	Int_t CurrentTrack() const { return fTrackIndex; };


public:  // Unfortunately ROOT requires the following to be public.

	class AliEventData : public TObject
	{
	public:
		AliEventData();
		AliEventData(Int_t eventnumber);
		virtual ~AliEventData();

		Int_t& EventNumber() { return fEventNumber; }
		TClonesArray& Blocks() { return fBlocks; }

	private:
	
		Int_t fEventNumber;  // Event number for this set of track blocks.
		TClonesArray fBlocks;  // The list of track blocks for this event.
		
		ClassDef(AliEventData, 1);  // Data per event.
	};

private:

	// Do not allow copying of this object.
	AliHLTMUONTrackSink(const AliHLTMUONTrackSink& /*object*/)
		: TObject(),
		  fFilename(""), fFoldername(""), fEventIndex(-1),
		  fCurrentEvent(NULL), fBlockIndex(-1), fCurrentBlock(NULL),
		  fTrackIndex(-1), fCurrentTrack(NULL),
		  fEventList(AliHLTMUONTrackSink::AliEventData::Class())
	{}

	AliHLTMUONTrackSink& operator = (const AliHLTMUONTrackSink& /*object*/) { return *this; }

	/* Sets all the current pointers to NULL and indices to -1.
	 */
	void ResetAllPointers() const;
	
	/* Sets the block and track pointers to NULL and indices to -1.
	 */
	void ResetBlockPointers() const;
	
	/* Sets just the current track pointer to NULL and index to -1.
	 */
	void ResetTrackPointers() const;

	TString fFilename;    // The file from which the track data was taken.
	TString fFoldername;  // The folder name from which track data was taken.
	
	mutable Int_t fEventIndex;            //! The current event index number.
	mutable AliEventData* fCurrentEvent;     //! Pointer to the currently selected event.
	mutable Int_t fBlockIndex;            //! The current block index number.
	mutable TClonesArray* fCurrentBlock;  //! Pointer to the currently selected block.
	mutable Int_t fTrackIndex;            //! The current track index number.
	mutable AliHLTMUONTrack* fCurrentTrack;         //! Pointer to the currently selected track.

	TClonesArray fEventList;  // List of tracks per event.

	ClassDef(AliHLTMUONTrackSink, 1)  // The data sink for track blocks for dHLT.
};


#endif // ALIHLTMUONTRACKSINK_H
