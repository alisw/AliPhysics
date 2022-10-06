/************************************************************************************
 * Copyright (C) 2021, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef __PWG__EMCAL_ALIEMCALFASTORMASKCONTAINER_H__
#define __PWG__EMCAL_ALIEMCALFASTORMASKCONTAINER_H__

#include <vector>
#include <TNamed.h>
#include <TObjArray.h>

class TBrowser;

namespace PWG {

namespace EMCAL {

/**
 * @class AliEmcalFastorMaskContainer
 * @brief Container for offline masked FastORs
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since June 1st, 2021
 * 
 * # Container for masked FastORs
 * 
 * The container structure is designed as OADB object for FastORs which
 * are tagged offline. FastORs can have two different tags:
 * - Dead: FastOR should always be ignored
 * - Bad: Warm FastORs with dedicated noise pattern, up to user analysis
 * 
 * The object has dedicated getters and setters for the channel status.
 * In order to display the masked channels, the object is made browsable.
 * As the Browse function internally relies on the trigger mapping in order
 * to convert the FastOR absolute ID into a positon in EMCAL the run version
 * (Run1 or Run2) needs to be specified when constructing the object.
 */
class AliEmcalFastorMaskContainer : public TNamed {
public:

  /**
   * @enum MaskStatus_t
   * @brief Type of FastOR masking
   */
  enum MaskStatus_t {
    kDeadFastOR = 0,            ///< FastOR is set to dead
    kBadFastOR = 1,             ///< FastOR is set to bad
    kStatusUndefined = -1       ///< Mask status not defined
  };

  /**
   * @enum RunType_t
   * @brief Type of the run
   */
  enum RunType_t {
    kRun1 = 0,                  ///< Run1 (2010-2013)
    kRun2 = 1,                  ///< Run2 (2015-2018)
    kRuntypeUndefined = -1      ///< Unknown
  };

  /**
   * @brief Dummy constructor
   */
  AliEmcalFastorMaskContainer();

  /**
   * @brief Constructor, setting also the name
   * @param name Name of the container
   */
  AliEmcalFastorMaskContainer(const char *name);

  /**
   * @brief Browse functionality: Draw map of masked FastORs
   * @param b Browser where to draw the FastORs.
   */
  virtual void Browse(TBrowser *b);

  /**
   * @brief Mark type as not folder
   * @return Always false
   */
  virtual Bool_t IsFolder() const { return false; }

  /**
   * @brief Destructor
   */
  virtual ~AliEmcalFastorMaskContainer() {}

  /**
   * @brief Clear list of masked FastORs (container empty afterwards)
   */
  virtual void Clear(Option_t * /*option*/ ="") { fMaskedFastors.Clear(); }

  /**
   * @brief Get list of all masked FastORs (bad and dead)
   * @return List of all masked FastORs
   */
  std::vector<int> GetMaskAll() const;

  /**
   * @brief Get list of dead FastORs
   * @return List of dead FastORs 
   */
  std::vector<int> GetMaskDead() const;

  /**
   * @brief Get list of bad FastORs
   * @return List of bad FastORs
   */
  std::vector<int> GetMaskBad() const;

  /**
   * @brief Get the the mask status of a dedicated FastOR
   * @param absID Absolute ID of the FastOR
   * @return Mask status of the FastOR (kUndefined if not found)
   */
  MaskStatus_t GetFastorMaskStatus(int absID) const;

  /**
   * @brief Check whether a FastOR with a given absolute ID is dead
   * @param absID Absolute ID of the FastOR
   * @return True if the FastOR is dead, false otherwise 
   */
  bool IsDeadFastor(int absID) const { return GetFastorMaskStatus(absID) == MaskStatus_t::kDeadFastOR; }

  /**
   * @brief Check whether a FastOR with a given absolute ID is bad
   * @param absID Absolute ID of the FastOR 
   * @return True if the FastOR is bad, false otherwise 
   */
  bool IsBadFastor(int absID) const { return GetFastorMaskStatus(absID) == MaskStatus_t::kBadFastOR; }

  /**
   * @brief Check whether a FastOR with a given absolute ID is bad
   * @param absID Absolute ID of the FastOR 
   * @return True if the FastOR is bad, false otherwise 
   */
  bool IsGoodFastor(int absID) const { return GetFastorMaskStatus(absID) == MaskStatus_t::kStatusUndefined; }

  /**
   * @brief Add new dead FastOR to the mask container
   * @param absID Absolute ID of the FastOR
   */
  void AddDeadFastor(int absID) { SetFastorStatus(absID, MaskStatus_t::kDeadFastOR); }

  /**
   * @brief Add new bad FastOR to the mask container
   * @param absID Absolute ID of the FastOR
   */
  void AddBadFastor(int absID) { SetFastorStatus(absID, MaskStatus_t::kBadFastOR); }

  /**
   * @brief Mask FastOR as good (remove from container if present)
   * @param absID Absolute ID of the FastOR
   */
  void SetGoodFastor(int absID) { SetFastorStatus(absID, MaskStatus_t::kStatusUndefined); }

  /**
   * @brief Add list of dead FastORs
   * @param absIDs List of absolute FastOR IDs to be masked dead
   * 
   * In case the FastORs are already present in the container only
   * the status is updated, otherwise they are added as new entry.
   */
  void AddDeadFastors(std::vector<int> absIDs);

  /**
   * @brief Add list of bad FastORs
   * @param absIDs List of absolute FastOR IDs to be masked bad
   * 
   * In case the FastORs are already present in the container only
   * the status is updated, otherwise they are added as new entry.
   */
  void AddBadFastors(std::vector<int> absIDs);

  /**
   * @brief Get the number of masked FastORs in the container
   * @return Number of FastORs  
   */
  int Size() const { return fMaskedFastors.GetEntries(); }

  /**
   * @brief Set the mask status for a given FastOR status
   * @param absID Absolute ID of the FastOR
   * @param status Status (bad/dead/undefined)
   * 
   * In case the FastOR is not yet listed, add it in case it 
   * is bad or dead, otherwise just update the status. In case
   * the status in good (undefined), remove it from the list
   */
  void SetFastorStatus(int absID, MaskStatus_t status);

  /**
   * @brief Set the run type (run1 or run2) 
   * @param runtype Run type
   */
  void SetRunType(RunType_t runtype) { fRunType = runtype; }

  /**
   * @brief Set the runtype by run number
   * @param runnumber Run number of the current run
   */
  void SetRunTypeByRunNumber(int runnumber) { if(runnumber <= 197388) fRunType = RunType_t::kRun1; else fRunType = RunType_t::kRun2; }

  /**
   * @class MaskedFastor
   * @brief Helper class storing FastOR information within the AliEmcalFastorMaskContainer 
   * 
   * The object is ROOT and C++ sortable, according to the FastOR absolute ID. As information
   * it stores the FastOR absolute ID and the mask status.
   */
  class MaskedFastor : public TObject {
  public:

    /**
     * @brief Dummy constructor
     */
    MaskedFastor() : TObject(), fFastORAbsID(-1), fMaskStatus(MaskStatus_t::kStatusUndefined) {}

    /**
     * @brief Constructor, also initializing the FastOR object
     * @param absID  Absolute ID of the FastOR
     * @param status Mask status (bad/dead)
     */
    MaskedFastor(int absID, MaskStatus_t status): TObject(), fFastORAbsID(absID), fMaskStatus(status) {}

    /**
     * @brief Destructor
     */
    virtual ~MaskedFastor() {}

    /**
     * @brief Comparison operator for equalness
     * @param other FastOR to compare to
     * @return True in case the FastORs have the same absolute ID, false otherwise 
     */
    bool operator==(const MaskedFastor &other) const { return fFastORAbsID == other.fFastORAbsID; }

    /**
     * @brief Comparison operator for smaller
     * @param other FastOR to compare to
     * @return True if the FastOR abs ID is smaller than the other, false otherwise 
     */
    bool operator<(const MaskedFastor &other) const { return fFastORAbsID < other.fFastORAbsID; }

    /**
     * @brief Mark object as ROOT-sortable
     * @return Always true
     */
    virtual Bool_t IsSortable() const { return true; }

    /**
     * @brief ROOT check for equalness
     * @param obj Object to compare to
     * @return True if the other object is a MaskedFastor and the FastOR abs IDs are equal, false otherwise 
     */
    virtual Bool_t IsEqual(const TObject *obj) const {
      const MaskedFastor *otherfastor = dynamic_cast<const MaskedFastor *>(obj);
      if(!otherfastor) return false;
      return *this == *otherfastor;
    }

    /**
     * @brief ROOT comparison function
     * @param obj Object to compare to
     * @return -1 in case this object is smaller or other is not a MaskedFastor, 0 if the FastORs are equal, 1 if this FastOR is larger than the other
     */
    virtual Int_t Compare(const TObject *obj) const {
      const MaskedFastor *otherfastor = dynamic_cast<const MaskedFastor *>(obj);
      if(!otherfastor) return -1;
      if(*this == *otherfastor) return 0;
      else if (*this < *otherfastor) return -1;
      else return 1;
    }

    /**
     * @brief Get absolute ID of the FastOR
     * @return Absolute ID of the FastOR 
     */
    int GetAbsID() const { return fFastORAbsID; }

    /**
     * @brief Get the mask status of the FastOR (dead/bad)
     * @return Mask status of the FastOR
     */
    MaskStatus_t GetMaskStatus() const { return fMaskStatus; }

    /**
     * @brief Check whether the FastOR was masked as bad
     * @return True if the FastOR is bad, false otherwise 
     */
    bool IsBad() const { return fMaskStatus == MaskStatus_t::kBadFastOR; }

    /**
     * @brief Check if the FastOR was masked as dead 
     * @return True  if the FastOR is dead, false otherwise
     */
    bool IsDead() const { return fMaskStatus == MaskStatus_t::kDeadFastOR; }

    /**
     * @brief Set the absolute ID of the FastOR
     * @param absID Absolute ID of the FastOR
     */
    void SetAbsID(int absID) { fFastORAbsID = absID; }

    /**
     * @brief Set the mask status (bad/dead) of the FastOR 
     * @param status Mask status
     */
    void SetMaskStatus(MaskStatus_t status) { fMaskStatus = status; }

    /**
     * @brief Mask FastOR as bad
     */
    void SetBad() { SetMaskStatus(MaskStatus_t::kBadFastOR); }

    /**
     * @brief Mask FastOR as dead
     */
    void SetDead() { SetMaskStatus(MaskStatus_t::kDeadFastOR); }

  private:
    int             fFastORAbsID;     ///< Absolute ID of the FastOR
    MaskStatus_t    fMaskStatus;      ///< Mask Status (dead / bad)

    ClassDef(MaskedFastor, 1);
  };

private:

  /**
   * @brief Helper function finding masked FastOR object within the container
   * @param absID 
   * @return MaskedFastor* 
   */
  MaskedFastor *FindInContainer(int absID) const;

  TObjArray fMaskedFastors;           ///< Container with masked FastORs
  RunType_t fRunType;                 ///< Period in which run was taken, needed for Browse functionality

  ClassDef(AliEmcalFastorMaskContainer, 1);
};

/**
 * @class TestAliEmcalFastorMaskContainer
 * @brief Unit test for the EMCAL Fastor mask container
 * @ingroup EMCALTESTS
 */
class TestAliEmcalFastorMaskContainer : public TObject {
public:

  /**
   * @brief Constructor
   */
  TestAliEmcalFastorMaskContainer() : TObject() {}

  /**
   * @brief Destructor
   */
  virtual ~TestAliEmcalFastorMaskContainer() {}

  /**
   * @brief Test runner for all tests
   * @return true All tests passed
   * @return false At least one test did not pass
   */
  bool RunAllTests() const;

  /**
   * @brief Run in-memory test
   * @return true Test passed
   * @return false Test failed
   * 
   * Preparing a container with example bad and dead channel list,
   * and check whether the list of bad/dead channels matches the 
   * list inserted. Container is not written to file.
   */
  bool TestInMemory() const;

  /**
   * @brief Test file IO
   * @return true Test passed
   * @return false Test failed
   * 
   * Preparing a container with example bad and dead channel list and
   * write it to file. Then the container is read again from file and 
   * compared to the input list to check wether the bad/dead channels 
   * match.
   */
  bool TestIO() const;

  /**
   * @brief Test changing the status of certain FastORs
   * @return true Test passed
   * @return false Test failed
   * 
   * Preparing a container with bad and dead channel list. Randomly
   * 5 channels are selected for which the status is changed from which
   * the status is changed from dead to bad and vice versa. Lists of 
   * dead and bad channels after the test must match with the input lists
   * including the change of the channel status
   */
  bool TestChangeStatus() const;

  /**
   * @brief Test of the Clean method
   * @return true Test passed
   * @return false Test failed
   * 
   * Preparing a container with bad and dead channel list. The clean method
   * is called then. After clean the list of bad and dead channels obtained
   * from the container must be empty.
   */
  bool TestClean() const;

private:

  /**
   * @brief Get example list with bad channels for unit tests
   * @return Example lists with bad channels
   */
  std::vector<int> GetTestBadChannels() const;

  /**
   * @brief Get example list with dead channels for unit test
   * @return Example list with dead channels 
   */
  std::vector<int> GetTestDeadChannels() const;

  /**
   * @brief Check where two vectors have the same content 
   * @param list1 First list in teh check for equalness
   * @param list2 Second list in the check for equalness
   * @return true Vectors are equal
   * @return false 
   */
  bool CheckEqual(const std::vector<int> &list1, const std::vector<int> &list2) const;

  /**
   * @brief Select randomly elements from a given inputlist
   * @param inputlist Input list
   * @param length Number of entries to select
   * @return List with entries obtained from the inputlist 
   */
  std::vector<int> GetRandomSubsample(std::vector<int> inputlist, int length) const;

  /**
   * @brief Helper function printing the content of the vector
   * @param data Vector to monitor
   * @param header Header printed above the vector content
   */
  void printVector(std::vector<int> data, const char *header) const;

  ClassDef(TestAliEmcalFastorMaskContainer, 1);
};

}

}

#endif // __PWG__EMCAL_ALIEMCALFASTORMASKCONTAINER_H__
