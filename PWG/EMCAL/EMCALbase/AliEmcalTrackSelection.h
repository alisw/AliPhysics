/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIEMCALTRACKSELECTION_H_
#define ALIEMCALTRACKSELECTION_H_

#include <TObject.h>
#include <TBits.h>
#include "AliEmcalTrackSelResultPtr.h"

class TClonesArray;
class TList;
class TObjArray;
class AliVCuts;
class AliVEvent;
class AliVTrack;

namespace PWG {
namespace EMCAL {

class AliEmcalCutBase;

}
}

/**
 * @class AliEmcalManagedObject
 * @brief Smart pointer implementation for objects inheriting from TObject
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @ingroup EMCALCOREFW
 * @since March 2dn, 2017
 *
 * Simple smart pointer implementation inheriting from TObject which defines ownership
 * in an object-by-object approach. The pointer will handle deletes only in case it is
 * owner over the object.
 *
 * Attention: This is a simplified version of a smart pointer, different pointers do not
 * know from each other, which might lead into troubles once the object is copied. Can be
 * done in a much more elegant way using c++11 shared_ptr.
 */
class AliEmcalManagedObject : public TObject {
public:
  /**
   * @brief Dummy constructor
   */
  AliEmcalManagedObject();

  /**
   * @brief Default constructor, creates new managed object.
   * @param[in] object Managed object
   * @param[in] owner Flag for ownership
   */
  AliEmcalManagedObject(TObject *object, bool owner = true);

  /**
   * @brief Copy constructor
   *
   * By default the new smart pointer will not own the object
   * of the reference smart pointer.
   * @param[in] ref Reference for the copy
   */
  AliEmcalManagedObject(const AliEmcalManagedObject &ref);

  /**
   * Assignment operator
   * By default the new smart pointer will not own the object
   * of the reference smart pointer.
   * @param[in] ref Reference for assignment
   * @return Pointer after assingment
   */
  AliEmcalManagedObject &operator=(const AliEmcalManagedObject &ref);

  /**
   * Destructor, will delete the managed object in case it is the owner
   */
  virtual ~AliEmcalManagedObject() { Cleanup(); }

  /**
   * @brief Check whether a managed object is set
   * @return True if object is set, false otherwise
   */
  operator bool() const { return fManagedObject != NULL; }

  /**
   * @brief Checks whether the object is set
   * @return True if pointer is owner of the object
   */
  bool IsOwner() const { return fOwner; }

  /**
   * @brief Specifying ownership over object
   * @param[in] owner Ownership (if true then object takes over ownership)
   */
  void SetOwner(bool owner = true) { fOwner = owner; }

  /**
   * @brief Set new managed object with ownership.
   *
   * In case the pointer managed already an object it will cleanup the object if necessary
   * @param[in] managedObject New object managed by this smart pointer
   * @param[in] owenr Ownership status of the new object
   */
  void SetObject(TObject *managedObject, bool owner = true) {
    Cleanup();
    fManagedObject = managedObject;
    fOwner = owner;
  }

  /**
   * @brief Providing access to managed object
   * @return Object managed by the smart pointer
   */
  TObject *GetObject() { return fManagedObject; }

protected:

  /**
   * Cleanup managed object (remove if pointer is owner)
   */
  void Cleanup();

private:
  Bool_t            fOwner;             ///< Switch defining ownership over object
  TObject           *fManagedObject;    ///< Pointer to object handled by the smart pointer

  /// \cond CLASSIMP
  ClassDef(AliEmcalManagedObject, 1);
  /// \endcond
};

/**
 * @class AliEmcalTrackSelection
 * @brief Interface for virtual track selection
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date Jan 30, 2016
 *
 * Interface for track selection within the EMCAL framework. Enables transparent track selection
 * for ESDs and AODs by implementing a wrapper derived from this class. The following abstract
 * functions need to be implemented by inheriting classes:
 * - GetAcceptedTracks (with TClonesArray and AliVEvent as parameters)
 * - IsTrackAccepted (with AliVTrackCuts)
 * - GenerateTrackCuts
 *
 * The usage of the virtual track selection is described here: \subpage VirtualTrackSelection
 */
class AliEmcalTrackSelection : public TObject {
public:

  /**
   * @enum ETrackFilterType_t
   * @brief Pre-defined track filters
   */
  enum ETrackFilterType_t {
    kNoTrackFilter = 0,   ///< No filter (all tracks passing)
    kCustomTrackFilter,   ///< Custom (user-defined) tracks
    kHybridTracks,        ///< Hybrid tracks
    kTPCOnlyTracks,       ///< TPC-only tracks
    kITSPureTracks,       ///< ITS stand-alone tracks
		kHybridTracks2010wNoRefit,	///< Hybrid tracks using the 2010 definition including NoITSrefit tracks (ESD-only)
		kHybridTracks2010woNoRefit,	///< Hybrid tracks using the 2010 definition excluding NoITSrefit tracks (ESD-only)
		kHybridTracks2011wNoRefit,	///< Hybrid tracks using the 2011 definition including NoITSrefit tracks (ESD-only)
		kHybridTracks2011woNoRefit,	///< Hybrid tracks using the 2011 definition excluding NoITSrefit tracks (ESD-only)
		kHybridTracks2018TRD				///< Hybrid tracks using the 2018 TRD test definition
  };

  /**
   * @brief Default consturctor
   *
   * Initialising objects with NULL, sets acception mode to ALL
   */
	AliEmcalTrackSelection();

	/**
	 * @brief Copy constructor
	 *
	 * Performing a flat copy
	 * @param[in] ref
	 */
	AliEmcalTrackSelection(const AliEmcalTrackSelection &ref);

	/**
	 * @brief Assingment operator
   *
   * Makes a flat copy
	 * @param[in] ref Reference for the copy
	 * @return Result of the copy
	 */
	AliEmcalTrackSelection &operator=(const AliEmcalTrackSelection &ref);

	/**
	 * @brief Destructor

   * Deletes track and track cut arrays. In case the
   * object has ownership over the track cuts itself,
   * it also deletes those
	 */
	virtual ~AliEmcalTrackSelection();

  /**
   * @brief Select tracks from a TClonesArray of input tracks
	 * 
	 * Running track selection for all tracks in the input array 
	 * (by calling the IsTrackAccepted function) and storing the track 
	 * selection result in the output container.
	 * Attention: The output container contains the result for all
	 * tracks, not only those selected. It is the user responsibility
	 * to filter the selected tracks.
   *
   * @param[in] tracks TClonesArray of tracks (must not be null)
   * @return TObjArray of AliEmcalTrackSelResultPtr with the selection status for all tracks
   */
	TObjArray *GetAcceptedTracks(const TClonesArray * const tracks);

	/**
	 * @brief Select tracks from a virtual event
	 *
	 * Running track selection for all tracks in the input event 
	 * (by calling the IsTrackAccepted function) and storing the track 
	 * selection result in the output container.
	 * Attention: The output container contains the result for all
	 * tracks, not only those selected. It is the user responsibility
	 * to filter the selected tracks.
	 * 
	 * @param[in] event AliVEvent, via interface of virtual event (must not be null)
	 * @return TObjArray of AliEmcalTrackSelResultPtr with the selection status for all tracks
	 */
	TObjArray *GetAcceptedTracks(const AliVEvent *const event);

	/**
	 * @brief Interface for track selection code
	 *
	 * This interface has to be implemented by the child classes.
	 * Here the selection of the track takes place. The function
	 * defines whether the track provided as argument is selected.
	 *
	 * Attention: The function has to support both AliESD/AliAODtack
	 * on the one side and AliPicoTrack on the other side in order
	 * to be compatible with the old EMCAL framework.
	 *
	 * @param[in] trk Track to be checked
	 * @return True if the track is accepted, false otherwise
	 */
	virtual PWG::EMCAL::AliEmcalTrackSelResultPtr IsTrackAccepted(AliVTrack * const trk) = 0;

	/**
	 * @brief Interface for track cut generators
	 *
	 * Track cut definitions might change period-by-period. Generators
	 * are methods to initialize the track cuts once the period is known.
	 * The method has default implementations for hybrid and global tracks.
	 *
	 * As the way of implementation differs between ESDs and AODs the
	 * method has to be implemented in the child class.
	 *
	 * Attention: In order to be usable for a given dataset, the initialization
	 * for the period / dataset needs to be implemented in the child class.
	 * Not all datasets are necessarily supported.
	 *
	 * @param[in] type Pre-defined track cuts type
	 * @param[in] period Period / Dataset name for which to initialize the track cuts
	 */
	virtual void GenerateTrackCuts(ETrackFilterType_t type, const char* period = "") = 0;

	/**
	 * @brief Add new track cuts to the list of cuts
	 *
	 * Takes ownership over the cuts
	 * @param[in] cuts New cuts to add
	 */
	void AddTrackCuts(AliVCuts *cuts);

	/**
	 * @brief Add track cuts of type AliEmcalCutBase
	 * 
	 * Takes ownership over the cuts
	 * @param[in] cuts New cuts to add
	 */
	void AddTrackCuts(PWG::EMCAL::AliEmcalCutBase *cuts);

	/**
	 * @brief Add new set of track cuts to the list of cuts
	 *
	 * Takes ownership over the cuts
	 * @param[in] cuts New set of cuts to add
	 */
	void AddTrackCuts(TObjArray *cuts);

	/**
	 * @brief Get the number of cut objects assigned.
	 * @return The number of cut objects
	 */
	Int_t GetNumberOfCutObjects() const;

	/**
	 * @brief Access to track cuts at a given position
	 * @param[in] icut Cut at position in array
	 * @return The cuts (NULL for invalid positions)
	 */
	PWG::EMCAL::AliEmcalCutBase *GetTrackCuts(Int_t icut);

	/**
	 * @brief Set selection mode to any
	 *
	 * In this case tracks are accepted if any of the
	 * cuts is passed.
	 */
	void SetSelectionModeAny() { fSelectionModeAny = kTRUE ; }

	/**
	 * @brief Set selection mode to all
	 *
	 * In this case tracks are only accepted if all cuts are
	 * passed.
	 */
	void SetSelectionModeAll() { fSelectionModeAny = kFALSE; }

	/**
	 * Saving QA objects to the output list.
	 * @param[in] outputList Common output list, used to store QA objects
	 */
	virtual void SaveQAObjects(TList *outputList) {}

protected:
	TObjArray    *fListOfTracks;         ///< TObjArray with accepted tracks
	TObjArray    *fListOfCuts;           ///< List of track cut objects
	Bool_t        fSelectionModeAny;     ///< Accept track if any of the cuts is fulfilled

	/// \cond CLASSIMP

	/// \endcond
};

#endif /* ALIEMCALTRACKSELECTION_H_ */
