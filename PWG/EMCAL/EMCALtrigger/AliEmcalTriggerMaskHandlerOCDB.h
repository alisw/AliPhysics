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
#ifndef PWG_EMCAL_ALIEMCALTRIGGERMASKHANDLEROCDB_H
#define PWG_EMCAL_ALIEMCALTRIGGERMASKHANDLEROCDB_H

#include <TObject.h>
#include <exception>
#include <vector>

class AliCDBManager;
class AliEMCALGeometry;
class AliEMCALTriggerDCSConfig;

namespace PWG
{

namespace EMCAL
{

/**
 * @class AliEmcalTriggerMaskHandlerOCDB
 * @brief Handler class for trigger mask at L0 and L1 stored in the OCDB
 * @ingroup EMCALTRGFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * 
 * # Easy access to the trigger masks in the OCDB
 * 
 * For monitoring, analysis, and replay or simulation of the trigger 
 * the online trigger mask is needed to be obtained from the OCDB. For
 * most use cases the mask should be in the form of the absolute ID of
 * the FastOR channel or the position of the FastOR in a col-row space
 * covering the entire EMCAL.
 * 
 * The trigger mask in the OCDB consists of two parts:
 * 
 * - Masked FastORs in within a TRU: These channels will contribute
 *   neither to L0 nor to L0
 * - Masked TRUs in the STU region: In this case the TRU was active,
 *   so it was used in the L0 trigger, but the L1 timesums from the
 *   TRU were not used in the STU. In this case all FastORs from the
 *   TRU need to be considered as masked.
 * 
 * As masked TRUs only affect the L1 trigger there are two functions
 * to provide the trigger mask which are either corresponding to L0
 * or L1.
 * 
 * While decoding of the STU region is trivial both for run1 and run2,
 * as the indexing is linearly increasing, the TRU region is more difficult
 * to decode: The mask consists of 6 bitsets each containing 16 bits,
 * resulting in total to 96 bits where each bit corresponds to a 
 * FastOR channel in the TRU. The channel index corresponding to a 
 * certain bit number differs between run1 and run2. The mapping translating
 * the bit number into a FastOR index within the TRU is applied automatically
 * and users don't need to worry about this detail.
 */
class AliEmcalTriggerMaskHandlerOCDB : public TObject {
public:
  /**
   * @struct FastORPosition
   * @brief Row/Col position of a FastOR
   */
  struct FastORPosition {
    int column; ///< Column (#eta)
    int row;    ///< Row (#phi)

    /**
     * @brief Comparison operator for equalty
     * @param other   FastOR position to compare to
     * @return true   Postition (row and column) of the two FastORs are identical
     * @return false  Positions differ
     */
    bool operator==(const FastORPosition &other) const { return column == other.column && row == other.row; }

    /**
     * @brief Comparison operator for smaller
     * @param other  FastOR position to compare to
     * @return true  Position is smaller
     * @return false Position is larger or equal
     * 
     * Comparison done first for rows and in case the rows are the same the 
     * columns are checked
     */
    bool operator<(const FastORPosition &other) const { return (48 * row + column) < (48 * other.row + other.column); } // 48 := number of FastORs in row
  };

  /**
   * @class GeometryNotSetException
   * @brief Error condition if the EMCAL geometry handler could not be created
   */
  class GeometryNotSetException : public std::exception {
  public:

    /**
     * @brief Constructor
     */
    GeometryNotSetException() {}

    /**
     * @brief Destructor
     */
    virtual ~GeometryNotSetException() throw() {}

    /**
     * @brief Get default error message
     * @return Error message 
     */
    virtual const char *what() const throw() { return "Geometry is not initalized"; }
  };

  /**
   * @class OCDBNotInitializedException
   * @brief Error handling in case the OCDB cannot be initialized
   * 
   * Error is thrown under any circumstance the AliEMCALTriggerDCSConfig cannot
   * be loaded from the OCDB. Possible error sources:
   * - OCDB path is not initialized
   * - Run number is not set
   * - Object is not found for the run in the OCDB
   */
  class OCDBNotInitializedException : public std::exception {
  public:
    /**
     * @brief Constructor
     */
    OCDBNotInitializedException() {}

    /**
     * @brief Destructor
     */
    virtual ~OCDBNotInitializedException() throw() {}

    /**
     * @brief Get default error message
     * @return Error message 
     */
    virtual const char *what() const throw() { return "OCDB Manager was not initialized"; }
  };

  /**
   * @brief Get current instance of the mask hander.
   * @return Current mask handler
   * 
   * Creating the handler if not yet existing, otherwise
   * just returning the currnet instance.
   * 
   * Attention: The handler is relying on the AliCDBManager
   * internally. The CDB Manager is expected to be configured
   * at least with the CDB path externally. If the calls are
   * done with run number -1 (current run) then for both the
   * CDB manager and the EMCAL geometry the run must be set 
   * outside.
   */
  static AliEmcalTriggerMaskHandlerOCDB *Instance();

  /**
   * @brief Destructor, cleaning internal cache (if needed)
   */
  virtual ~AliEmcalTriggerMaskHandlerOCDB() { ClearCache(); }

  /**
   * @brief Get list of masked FastORs with absolute ID masked at L0
   * @param runnumber                   Run for which to read the masked channels (-1 if current run is requested)
   * @return                            List of FastOR absolute IDs of channels mased at L0
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * Containing only the channels which are mased at L0. In case a TRU is masked
   * in the STU region this is not included.
   */
  std::vector<int> GetMaskedFastorIndicesL0(int runnumber = -1);

  /**
   * @brief Get list of masked FastORs as position in EMCAL masked at L0
   * @param runnumber                   Run for which to read the masked channels (-1 if current run is requested)
   * @return                            List of FastOR positions in EMCAL of channels mased at L0
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * Containing only the channels which are mased at L0. In case a TRU is masked
   * in the STU region this is not included.
   */
  std::vector<PWG::EMCAL::AliEmcalTriggerMaskHandlerOCDB::FastORPosition> GetMaskedFastorPositionsL0(int runnumber = -1);

  /**
   * @brief Get list of masked FastORs with absolute ID masked at L1
   * @param runnumber                   Run for which to read the masked channels (-1 if current run is requested)
   * @return                            List of FastOR absolute IDs of channels mased at L1
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * Containing channels which are mased at L1. Those are all channels masked at L0, however
   * if a full TRU is mased in the STU region the FastORs from this TRU are masked as well. 
   */
  std::vector<int> GetMaskedFastorIndicesL1(int runnumber = -1);

  /**
   * @brief Get list of masked FastORs as position in EMCAL masked at L1
   * @param runnumber                   Run for which to read the masked channels (-1 if current run is requested)
   * @return                            List of FastOR positions in EMCAL of channels mased at L1
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * Containing channels which are mased at L1. Those are all channels masked at L0, however
   * if a full TRU is mased in the STU region the FastORs from this TRU are masked as well. 
   */
  std::vector<PWG::EMCAL::AliEmcalTriggerMaskHandlerOCDB::FastORPosition> GetMaskedFastorPositionsL1(int runnumber = -1);

  /**
   * @brief Get global indices of TRUs masked in the STU region
   * @param runnumber                   Run for which to read the masked TRUs (-1 if current run is requested)
   * @return                            List of masked TRUs in global number scheme 
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   */
  std::vector<int> GetGlobalMaskedTRUIndices(int runnumber = -1);

  /**
   * @brief Create monitoring histogram with all masked FastORs at L0
   * @param runnumber                   Run for which to obtain the monitoring histogram (-1 if current run is requested)
   * @return                            Histogram all masked FastOR positions in column row space
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * See GetMaskedFastorPositionsL0 for the definition of the masked FastORs at Level0. 
   * The monitoring histogram shows masked FastORs with index 1 and active FastORs with
   * index 0. The column-row space contains the full EMCAL no matter whether the run is from
   * run1 or run2.
   */
  TH2 *MonitorMaskedFastORsL0(int runnumber = -1);

  /**
   * @brief Create monitoring histogram with all masked FastORs at L1
   * @param runnumber                   Run for which to obtain the monitoring histogram (-1 if current run is requested)
   * @return                            Histogram all masked FastOR positions in column row space
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   * 
   * See GetMaskedFastorPositionsL1 for the definition of the masked FastORs at Level1. 
   * The monitoring histogram shows masked FastORs with index 1 and active FastORs with
   * index 0. The column-row space contains the full EMCAL no matter whether the run is from
   * run1 or run2.
   */
  TH2 *MonitorMaskedFastORsL1(int runnumber = -1);

  /**
   * @brief Conversion function between position in reg masks and channel ID for run1 setup
   * @param mask      Number of the mask
   * @param bitnumber Number of the bit in the mask
   * @param onethirdsm Mark whether the supermodule is of type 1/3
   * @return          Channel ID beloging to bit number in reg mask
   * 
   * The reg mask consists of 6 separate masks each containing 16 bits,
   * resulting in 96 bits in total. In order to determin the channel ID
   * both mask number and bit number are needed
   */
  static int GetChannelForMaskRun1(int mask, int bitnumber, bool /*onethirdsm*/);

  /**
   * @brief Conversion function between position in reg masks and channel ID for run2 setup
   * @param mask      Number of the mask
   * @param bitnumber Number of the bit in the mask
   * @param onethirdsm Mark whether the supermodule is of type 1/3
   * @return          Channel ID beloging to bit number in reg mask
   * 
   * The reg mask consists of 6 separate masks each containing 16 bits,
   * resulting in 96 bits in total. In order to determin the channel ID
   * both mask number and bit number are needed
   */
  static int GetChannelForMaskRun2(int mask, int bitnumber, bool onethirdsm);

  /**
   * @brief Draw L0 trigger mask from DCS config in cvmfs OCDB
   * @param runnumber Run number for which to draw trigger mapping
   * 
   * Creating plot with the mask channels at L0 in eta-phi space.
   * For easier reading, a frame is highlighting supermodules,
   * TRUs and FastORs. In case a canvas exists, this canvas is used
   * for plotting, otherwise a new canvas is created. 
   * 
   * Attention: Access to the OCDB on cvmfs is requested. The run
   * number must be a valid physics run. Stand-alone runs or 
   * technical runs are not guaranteed.
   */
  void drawL0MaskFromCVMFS(int runnumber) { drawMaskFromCVMFS(runnumber, false); }

  /**
   * @brief Draw L1 trigger mask from DCS config in cvmfs OCDB
   * @param runnumber Run number for which to draw trigger mapping
   * 
   * Creating plot with the mask channels at L1 in eta-phi space.
   * For easier reading, a frame is highlighting supermodules,
   * TRUs and FastORs. In case a canvas exists, this canvas is used
   * for plotting, otherwise a new canvas is created. 
   * 
   * Attention: Access to the OCDB on cvmfs is requested. The run
   * number must be a valid physics run. Stand-alone runs or 
   * technical runs are not guaranteed.
   */
  void drawL1MaskFromCVMFS(int runnumber) { drawMaskFromCVMFS(runnumber, true); }

private:
  /**
   * @brief Dummy constructor
   * 
   * Class is a singleton. Use Instance() in order to get the current mask handler
   */
  AliEmcalTriggerMaskHandlerOCDB() : TObject(), fCurrentRunnumber(-1), fCurrentConfig(nullptr) {}

  AliEmcalTriggerMaskHandlerOCDB(const AliEmcalTriggerMaskHandlerOCDB &);
  AliEmcalTriggerMaskHandlerOCDB &operator=(AliEmcalTriggerMaskHandlerOCDB &);

  /**
   * @brief Clear current cache entry for the trigger DCS config
   * 
   * Also setting cache run number to -1
   */
  void ClearCache();

  /**
   * @brief Initialize CDB manager and update cached trigger config
   * @param runnumber                   Run number for which to obtain the trigger config (-1 in case the current run should be taken)
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * 
   * In case a specific run number is requested or the run number differs
   * from the current cache run, the CDB manager is updated and the trigger 
   * DCS config for the new run is loaded into the cache.
   */
  void UpdateCache(int runnumber);

  /**
   * @brief Initializing CDB Manager
   * @param runnumber                   Run number for which to configure the CDB manager (-1 in case the current run should be taken)
   * @return                            Fully initialized CDB manager 
   * @throw OCDBNotInitializedException Failure in reading the trigger config from the OCDB
   * 
   * Instancing the CDB manager. In case the a certain run is requested
   * and the run differs from the run set in the CDB manager the CDB
   * manager is updated with the current run.
   */
  AliCDBManager *InitCDB(int runnumber);

  /**
   * @brief Get the Geometry object
   * 
   * @param runnumber 
   * @return AliEMCALGeometry* 
   * @throw GeometryNotSetException     Failure in getting the geometry for the current run
   */
  AliEMCALGeometry *GetGeometry(int runnumber);

  /**
   * @brief Creating histogram template for FastOR masks monitoring histograms
   * @param name   Name of the histogram
   * @param title  Title of the histogram
   * @param run2   Flag whether the run number is from run1 or run2
   * @return       Histogram template
   * 
   * Helper function initializing histogram template for the monitoring histograms
   * used in the fuctions MonitorMaskedFastORsL0 and MonitorMaskedFastORsL1
   */
  TH2 *PrepareHistogram(const char *name, const char *title, bool run2);

  /**
   * @brief Fill monitoring histogram with masked FastOR channels
   * @param outputhist Histogram to be filled with masked FastORs
   * @param fastors    List of masked FastORs
   * 
   * Helper function to fill the monitoring histograms (templates created in 
   * PrepareHistogram) with masked FastOR channels, used in the fuctions
   * MonitorMaskedFastORsL0 and MonitorMaskedFastORsL1
   */
  void FillMaskedFastors(TH2 *outputhist, const std::vector<AliEmcalTriggerMaskHandlerOCDB::FastORPosition> &fastors) const;

  /**
   * @brief Draw frame for run1 geometry
   * 
   * Drawing frames for supermodules, TRUs and FastORs using the run1 geometry.
   * Setup consisting of full EMCAL, TRUs are aligned in phi direction. Small 
   * lines indicate FastORs within TRU, intermediate size lines mark TRUs within 
   * supermodules and big lines supermodules.
   */
  void drawRun1Frame();

  /**
   * @brief Draw frame for run2 geometry
   * 
   * Drawing frames for supermodules, TRUs and FastORs using the run1 geometry.
   * Setup consisting of full EMCAL and DCAL, TRUs are aligned in eta direction 
   * except for the samll supermodules. Small lines indicate FastORs within TRU, 
   * intermediate size lines mark TRUs within supermodules and big lines supermodules.
   */
  void drawRun2Frame();

  /**
   * @brief Draw trigger mask at L0 or L1 for a given run number using the OCDB on cvmfs
   * @param runnumber Run number 
   * @param isLevel1 If true the mask is shown at L0, otherwise at L1
   * 
   * Helper function creating the monitoring plots of the trigger mask at L0 or at L1, 
   * used in drawL0MaskFromCVMFS and drawL1MaskFromCVMFS
   */
  void drawMaskFromCVMFS(int runnumber, bool isLevel1);


  Int_t fCurrentRunnumber;                  //!<! Current run number of cache
  AliEMCALTriggerDCSConfig *fCurrentConfig; //!<! Cache object

  static AliEmcalTriggerMaskHandlerOCDB *fgInstance; //!<! Current instance of the Mask handler (for singleton)

  ClassDef(AliEmcalTriggerMaskHandlerOCDB, 1);
};

;

}

}

#endif // PWG_EMCAL_ALIEMCALTRIGGERMASKHANDLER_H