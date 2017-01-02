#ifndef ALIEMCALTRACKSELRESULTPTR_H
#define ALIEMCALTRACKSELRESULTPTR_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <iosfwd>

class AliVTrack;

/**
 * @class AliEmcalTrackSelResultPtr
 * @brief Structure containing the result of a given track selection step
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since July 5, 2016
 *
 * Smart pointer containing selection results. In addition to the object itself a set
 * of further information is stored.
 * - The track object itself
 * - Its selection status (true - selected, false - rejected)
 * - A flag (bitmap) determining a track category
 * - Rejection reason
 *
 * Access to the original track is done via the operators * and ->. Using the operator->
 * provides direct access to public member functions of AliVTrack. As this, the following is
 * valid C++ code:
 *
 * ~~~{.cxx}
 * AliEmcalTrackSelection sel;
 * AliVTrack testtrack;
 * AliEmcalTrackSelPtr result = sel.IsTrackAccepted(&testtrack);
 * if(result.IsGetSelectionResult()){
 *   std::cout << "Pt of the selected track: " << result->Pt() << std::endl;
 * }
 * ~~~
 *
 * This smart pointer does not take care about deleting the object as it does not own
 * the content behind. Instead the user has to take care.
 */
class AliEmcalTrackSelResultPtr : public TObject {
public:

  /**
   * Dummy constructor, nothing to do
   * Default values:
   * - Track: nullptr
   * - Selection status: false
   * - Flag: 0
   */
  AliEmcalTrackSelResultPtr();

  /**
   * Constructor, fully initializing the result with a track pointer,
   * @param[in] trk Pointer to the original track
   * @param[in] selectionStatus Status of the selection (true - selected, false - rejected)
   * @param[in] flag Flag for a track category (optional, default: 0)
   */
  AliEmcalTrackSelResultPtr(AliVTrack *trk, Bool_t selectionStatus, ULong_t flag = 0);

  /**
   * Copy constructor, copying infromation (for track pointer only the pointer itself)
   * @param[in] ref Reference for the copy
   */
  AliEmcalTrackSelResultPtr(const AliEmcalTrackSelResultPtr &ref);

  /**
   * Move constructor, implements move procedure for track selection results
   * @param[in] ref Pointer to be moved
   */
  AliEmcalTrackSelResultPtr(AliEmcalTrackSelResultPtr &&ref);

  /**
   * Assignment operator, copies information from a reference into this object. As
   * for the copy and move constructor only
   * @param[in] ref Reference for the assignment
   * @return Object after assignment
   */
  AliEmcalTrackSelResultPtr &operator=(const AliEmcalTrackSelResultPtr &ref);

  /**
   * Move assignment operator, copies information from a reference into this object. As
   * for the copy and move constructor only
   * @param[in] ref Reference for the assignment
   * @return Object after assignment
   */
  AliEmcalTrackSelResultPtr &operator=(AliEmcalTrackSelResultPtr &&ref);

  /**
   * Destructor, nothing to do
   */
  virtual ~AliEmcalTrackSelResultPtr() {}

  /**
   * Comparison operator, checking for equalness. Comparison is based
   * on the address of the track pointer.
   * @param[in] other Object to be compared to
   * @return True if the underlying track pointers are equal, false otherwise
   */
  Bool_t operator==(const AliEmcalTrackSelResultPtr &other) const;

  /**
   * Comparison operator, checking whether this object is smaller than
   * the reference object. Comparison is done based on the address of the
   * track pointer.
   * @param[in] other Object to be compared to
   * @return
   */
  Bool_t operator<(const AliEmcalTrackSelResultPtr &other) const;

  /**
   * ROOT comparison function, checking for equalness. Comparison is based
   * on the address of the track pointer.
   * @param[in] o Object to be compared to
   * @return True if the underlying track pointers are equal, false otherwise
   */
  virtual Bool_t IsEqual(const TObject *o) const;

  virtual Int_t Compare(const TObject *o) const;

  /**
   * Pointer operator, returning pointer to underlying track
   * @return Pointer to the underlying track
   */
  AliVTrack * operator*() const;

  /**
   * Member access operator, returning pointer to the
   * underlying track
   * @return Pointer to the underlying track
   */
  AliVTrack * operator->() const;

  /**
   * Boolean conversion operator. Access to the selection status
   * @return track selection status
   */
  operator bool() const { return fSelectionResult; }

  /**
   * Streaming operator, print content of the track selection result to the
   * stream
   * @param stream Stream for the printing
   * @param ref Object to be streamed
   * @return Reference to the original stream after printing the object
   */
  friend std::ostream &operator<<(std::ostream &stream, const AliEmcalTrackSelResultPtr &ref);

  /**
   * Create a stream representation of the object and put it on then stream;
   * @param[in] stream Stream used to print the object
   */
  void PrintStream(std::ostream &stream) const;

  /**
   *
   * @param track
   */
  void SetTrack(AliVTrack *track) { fTrack = track; }

  /**
   * Set the track selection status
   * @param[in] selectionResult True if the track is selected, false otherwise;
   */
  void SetSelectionResult(Bool_t selectionResult) { fSelectionResult = selectionResult; }

  /**
   * Set the track category flag
   * @param flag Flag assigned to the track
   */
  void SetFlag(ULong_t flag) { fFlag = flag; }

  /**
   * Access to underlying track
   * @return Underlying track
   */
  AliVTrack *GetTrack() const { return fTrack; }

  /**
   * Get the flag assigned to the track selection status
   * @return Track flag
   */
  ULong_t GetFlag() const { return fFlag; }

  /**
   * Get the track selection status
   * @return True if the track is selected, false otherwise
   */
  Bool_t GetSelectionResult() const { return fSelectionResult; }

protected:
  AliVTrack                               *fTrack;            ///< Pointer to selected track
  Bool_t                                   fSelectionResult;  ///< Result of the track selection (true - selected, false - rejected)
  ULong_t                                  fFlag;             ///< Selection flag (optional)

  /// \cond CLASSIMP
  ClassDef(AliEmcalTrackSelResultPtr, 1);
  /// \endcond
};

/**
 * Streaming operator, print content of the track selection result to the
 * stream
 * @param stream Stream for the printing
 * @param ref Object to be streamed
 * @return Reference to the original stream after printing the object
 */
std::ostream &operator<<(std::ostream &stream, const AliEmcalTrackSelResultPtr &ref);

#endif /* ALIEMCALTRACKSELRESULTPTR_H */
