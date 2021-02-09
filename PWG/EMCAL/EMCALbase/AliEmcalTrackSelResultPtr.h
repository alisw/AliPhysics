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
#ifndef ALIEMCALTRACKSELRESULTPTR_H
#define ALIEMCALTRACKSELRESULTPTR_H

#include <TObject.h>
#include <iosfwd>

class AliVTrack;

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalTrackSelResultUserStorage
 * @brief Helper class handling the lifetime of the user object handled by AliEmcalTrackSelResultUserPtr
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 11, 2017
 * 
 * Common storage class for multiple instances of AliEmcalTrackSelResultUserPointers.
 * Together with the data the storage keeps track of the number of pointer instances
 * connected to it via the reference count. Instances can connect and disconnect
 * by the functions Connect and Release.
 * 
 * This class is foreseen to be used only from within AliEmcalTrackSelResultUserPtr. It
 * is a suicide object which kills itself once it goes out-of-scope.
 */
class AliEmcalTrackSelResultUserStorage : public TObject {
public:
  AliEmcalTrackSelResultUserStorage();
  
  /**
   * @brief Main constructor
   * 
   * Sets the data and the reference count to 1
   * @param data 
   */
  AliEmcalTrackSelResultUserStorage(TObject *data);
  
  /**
   * @brief Connect new user pointer instance to the storage
   * 
   * Used in the copy constructor and assignment operator of 
   * AliEmcalTrackSelResultUserPtr. This function is there in
   * order to tell then storage that a new user pointer instance
   * has connected to the storage. Increases the reference count.
   * Not to be used outside of class AliEmcalTrackSelResultUserPtr
   */
  void Connect();

  /**
   * @brief Release user pointer from the storage
   * 
   * Used in the destructor of AliEmcalTrackSelUserPtr. This
   * function is there in order to tell the storage that
   * a user pointer instance has disconnected. Reduces the
   * reference count. Not to be called ountside of class 
   * AliEmcalTrackSelResultUserPtr.
   */
  void Release();

  /**
   * @brief Get the number of pointer instances connected to the storage
   * 
   * @return Number of pointer instances connected to the storage
   */
  Int_t GetReferenceCount() const { return fReferenceCount; }

  /**
   * @brief Get the user data handled by the storage
   * 
   * @return User data handled by the storage
   */
  TObject *GetData() const { return fData; }

protected:

  /**
   * @brief Destructor
   * 
   * Deletes the data associated to the storage. Only to be used
   * from within AliEmcalTrackSelResultUserPtr
   */
  virtual ~AliEmcalTrackSelResultUserStorage();

private:
  AliEmcalTrackSelResultUserStorage(const AliEmcalTrackSelResultUserStorage &);
  AliEmcalTrackSelResultUserStorage(const AliEmcalTrackSelResultUserStorage &&);
  AliEmcalTrackSelResultUserStorage &operator=(const AliEmcalTrackSelResultUserStorage &);
  AliEmcalTrackSelResultUserStorage &operator=(const AliEmcalTrackSelResultUserStorage &&);
  TObject           *fData;               ///< User data
  Int_t             fReferenceCount;      ///< Reference counter

  /// \cond CLASSIMP
  ClassDef(AliEmcalTrackSelResultUserStorage, 1);
  /// \endcond
};

/**
 * @class AliEmcalTrackSelResultUserPtr 
 * @brief Handler for user objects attached to the track selection result ptr 
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 11, 2017
 *
 * Handling the lifetime of the object attached to the AliEmcalTrackSelResultPtr
 * using reference counting technique. Pointer instances created from the current
 * pointer using the copy constructor or assignment operator connect to the 
 * common storage and disconnect once they are constructed. The last instance connected
 * to a storage deletes it.
 * 
 * This class is foreseen to be used only from within AliEmcalTrackSelResultPtr.
 */
class AliEmcalTrackSelResultUserPtr : public TObject {
public:
  /**
   * @brief Dummy constructor
   * 
   * Not handling object
   */
  AliEmcalTrackSelResultUserPtr();

  /**
   * @brief Main constructor setting the data
   * 
   * Allocating the storage for the data hanlding and initializes
   * the reference count
   * 
   * @param o Object to be handled by the user pointer
   */
  AliEmcalTrackSelResultUserPtr(TObject *o);

  /**
   * @brief Copy constructor
   * 
   * Creating a new pointer instance connected to the same storage
   * 
   * @param ref Reference for the copy
   */
  AliEmcalTrackSelResultUserPtr(const AliEmcalTrackSelResultUserPtr &ref);

  /**
   * @brief Move constructor
   * 
   * Only transferring the storage, has no consequence on the reference count
   * 
   * @param ref Object to be moved
   */
  AliEmcalTrackSelResultUserPtr(AliEmcalTrackSelResultUserPtr &&ref);

  /**
   * @brief Assignment operator
   * 
   * Connecting pointer instance to the storage connected to by ref. In
   * case the pointer handles a storage itself disconnects from the storage
   * 
   * @param ref Reference for assignment
   * @return Pointer instance after assignment
   */
  AliEmcalTrackSelResultUserPtr &operator=(const AliEmcalTrackSelResultUserPtr &ref);

  /**
   * @brief Move assignment operator
   * 
   * Connecting pointer instance to the storage connected to by ref. In
   * case the pointer handles a storage itself disconnects from the storage
   * 
   * @param ref Object to be moved
   * @return Pointer instance after move
   */
  AliEmcalTrackSelResultUserPtr &operator=(AliEmcalTrackSelResultUserPtr &&ref);

  /**
   * @brief Destructor
   * 
   * Disconnects pointer instance from the common storage. If it is the last
   * pointer instance conneted to the storage deletes the storage as well.
   */
  virtual ~AliEmcalTrackSelResultUserPtr();

  /**
   * @brief Get the object handled by the storage
   * 
   * @return User object handled by the pointer if storage exists, nullptr otherwise 
   */
  TObject *GetData() const { return fUserStorage ? fUserStorage->GetData() : nullptr; }

private:
  AliEmcalTrackSelResultUserStorage         *fUserStorage;        ///< Underlying user storage for reference counting

  /// \cond CLASSIMP
  ClassDef(AliEmcalTrackSelResultUserPtr, 1);
  /// \endcond
};

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
 * - A user-info object (of type TObject) for additional information
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
   * - User object: nullptr
   */
  AliEmcalTrackSelResultPtr();

  /**
   * Constructor, fully initializing the result with a track pointer,
   * 
   * Taking 
   * @param[in] trk Pointer to the original track
   * @param[in] selectionStatus Status of the selection (true - selected, false - rejected)
   * @param[in] userobject User information defined by the track selection (default: nullptr)
   */
  AliEmcalTrackSelResultPtr(AliVTrack *trk, Bool_t selectionStatus, TObject *userobject = nullptr);

  /**
   * Copy constructor, copying infromation (for track pointer only the pointer itself)
   * @param[in] ref Reference for the copy
   */
  AliEmcalTrackSelResultPtr(const AliEmcalTrackSelResultPtr &ref);

#if !defined(__CINT__) && !defined(__MAKECINT__)
  /**
   * Move constructor, implements move procedure for track selection results
   * @param[in] ref Pointer to be moved
   */
  AliEmcalTrackSelResultPtr(AliEmcalTrackSelResultPtr &&ref);
#endif

  /**
   * Assignment operator, copies information from a reference into this object. As
   * for the copy and move constructor only
   * @param[in] ref Reference for the assignment
   * @return Object after assignment
   */
  AliEmcalTrackSelResultPtr &operator=(const AliEmcalTrackSelResultPtr &ref);

#if !defined(__CINT__) && !defined (__MAKECINT__)
  /**
   * Move assignment operator, copies information from a reference into this object. As
   * for the copy and move constructor only
   * @param[in] ref Reference for the assignment
   * @return Object after assignment
   */
  AliEmcalTrackSelResultPtr &operator=(AliEmcalTrackSelResultPtr &&ref);
#endif

  /**
   * @brief Destructor
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
   * Create a stream representation of the object and put it on then stream;
   * @param[in] stream Stream used to print the object
   */
  void PrintStream(std::ostream &stream) const;

  /**
   * @brief Set the track object handled by the track selection
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
  void SetUserInfo(TObject *userinfo) { fUserInfo = AliEmcalTrackSelResultUserPtr(userinfo); }

  /**
   * Access to underlying track
   * @return Underlying track
   */
  AliVTrack *GetTrack() const { return fTrack; }

  /**
   * Get the flag assigned to the track selection status
   * @return User information
   */
  const TObject *GetUserInfo() const { return fUserInfo.GetData(); }

  /**
   * Get the track selection status
   * @return True if the track is selected, false otherwise
   */
  Bool_t GetSelectionResult() const { return fSelectionResult; }

protected:
  AliVTrack                               *fTrack;            ///< Pointer to selected track
  Bool_t                                   fSelectionResult;  ///< Result of the track selection (true - selected, false - rejected)
  AliEmcalTrackSelResultUserPtr            fUserInfo;         ///< Selection flag (optional)

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

/**
 * @class TestAliEmcalTrackSelResultPtr
 * @brief Unit test for class AliEmcalTrackSelResultPtr
 * @ingroup EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Dec 18, 2017
 * 
 * Unit test for AliEmcalTrackSelResultPtr, covering
 * - operator bool
 * - Copy constructor
 * - Assignment operator
 * - user info
 */
class TestAliEmcalTrackSelResultPtr : public TObject {
public:
  /**
   * @brief Constructor
   */
  TestAliEmcalTrackSelResultPtr() {}

  /**
   * @brief Destructor
   */
  virtual ~TestAliEmcalTrackSelResultPtr() {}

  /**
   * @brief Run test suite
   * 
   * Running all unit tests implemented for AliEmcalTrackSelResultPtr
   * - operator bool
   * - Assignment operator
   * - Copy constructor
   * - User information 
   * 
   * @return true All tests passed
   * @return false At least one test failed
   */
  bool RunAllTests() const;

  /**
   * @brief Test for operator bool
   * 
   * Both cases tested:
   * - track selection is true, operator must return true
   * - track selection is false, operator must return false
   * 
   * @return true Test passed
   * @return false Test failed
   */
  bool TestOperatorBool() const;

  /**
   * @brief Test copy constructor with user information
   * 
   * Test prepared with 10 instances of user objects and 10 instances without user objects
   * 
   * @return true All tests passed
   * @return false At least one test failed
   */
  bool TestCopyConstructor() const;

  /**
   * @brief Tests assignment operatator with user info
   * 
   * Test prepared with 10 instances of user objects and 10 instances without user objects
   * 
   * @return true All tests passed
   * @return false All tests failed
   */
  bool TestOperatorAssign() const;

  /**
   * @brief Test handling of user storage
   * 
   * Prepared 10 pointers with and 10 without storage. Pointer must return
   * - user object for test with user object
   * - nullptr for test without user object
   * 
   * @return true Tests passed
   * @return false At least one test failed
   */
  bool TestUserInfo() const;

protected:

  bool AssertBool(const AliEmcalTrackSelResultPtr &test, bool testvalue) const;

  bool AssertPayload(const AliEmcalTrackSelResultPtr &test, void *payload) const;

  /// \cond CLASSIMP
  ClassDef(TestAliEmcalTrackSelResultPtr, 1);
  /// \endcond
};

}

}

#endif /* ALIEMCALTRACKSELRESULTPTR_H */
