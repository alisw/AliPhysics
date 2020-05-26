#ifndef ALIEMCALTRIGGERSTUDCSCONFIG_H
#define ALIEMCALTRIGGERSTUDCSCONFIG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "TObject.h"
#include <iosfwd>
#include <string>

class TVector2;
class TClonesArray;

//________________________________________________
/// \class AliEMCALTriggerSTUDCSConfig
/// \ingroup EMCALbase
/// \brief EMCal trigger STU DCS Config
///
///  EMCAL STU DCS parameters to be stored in OCDB
///  Add more comment
///
/// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
//________________________________________________

class AliEMCALTriggerSTUDCSConfig : public TObject 
{

public:

  //________________________________________________
  /// \class AliEMCALTriggerSTUTRUErrorCount
  /// \ingroup EMCALbase
  /// \brief EMCal trigger STU  error count handling
  ///
  ///  EMCAL STU error count error count handling
  ///  Add more comment
  ///
  /// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
  //________________________________________________

  class AliEMCALTriggerSTUTRUErrorCount : public TObject 
  {
    
  public:
    
    AliEMCALTriggerSTUTRUErrorCount(): TObject(), fTime(0), fErrorCount(0) {}
    AliEMCALTriggerSTUTRUErrorCount(Int_t time, ULong64_t errorCount) : TObject(), fTime(time), fErrorCount(errorCount) {}
    ~AliEMCALTriggerSTUTRUErrorCount() {}
    
    virtual Bool_t  IsEqual(const TObject *o) const;
    virtual Int_t   Compare(const TObject *o) const;
    
    void            SetValue(Int_t time, ULong64_t errorcount){ fTime = time; fErrorCount = errorcount; }
    Int_t           GetTime()       const { return fTime       ; }
    ULong64_t       GetErrorCount() const { return fErrorCount ; }
    
  private:
    Int_t                     fTime;        ///< Time
    ULong_t                   fErrorCount;  ///< Error count
    
    /// \cond CLASSIMP
    ClassDef(AliEMCALTriggerSTUTRUErrorCount, 1) ;
    /// \endcond
    
  };

  /**
   * @brief Default constructor
   */
  AliEMCALTriggerSTUDCSConfig();

  /**
   * @brief Destructor
   */
  virtual ~AliEMCALTriggerSTUDCSConfig();
  
  /**
   * @brief Copy constructor
   * @param cd Reference for the copy
   */
  AliEMCALTriggerSTUDCSConfig           (const AliEMCALTriggerSTUDCSConfig &cd);


  /**
   * @brief Equalty operator
   * 
   * Checks whether the two STU DCS configurations are the same. For equalty
   * all setting must be the same. Error counters are not considered.
   */
  bool operator==(const AliEMCALTriggerSTUDCSConfig &other) const;

  /**
   * @brief Streaming operator
   * 
   * Printing all STU DCS settings on the output stream, except for the
   * error counters
   */
  friend std::ostream &operator<<(std::ostream &stream, const AliEMCALTriggerSTUDCSConfig &config);

  /**
   * @brief Serialize object to JSON format
   * 
   * @return JSON-serialized TRU DCS config object 
   */
  std::string ToJSON() const;
  
  /**
   * @brief Set threshold parameter for the gamma trigger 
   * @param vzpar Index in the VZERO-dependent parabola equation (0 - A, 1 - B, 2 - C)
   * @param ithresh Threhold index (0 - high threshold, 1 - low threshold)
   * @param val Threshold value
   */
  void    SetG(Int_t vzpar, Int_t ithresh, Int_t val) { fG[vzpar][ithresh]    = val; }

  /**
   * @brief Set threshold parameter for the jet trigger
   * @param vzpar Index in the VZERO-dependent parabola equation (0 - A, 1 - B, 2 - C)
   * @param ithresh Threhold index (0 - high threshold, 1 - low threshold)
   * @param val Threshold value
   */
  void    SetJ(Int_t vzpar, Int_t ithresh, Int_t val) { fJ[vzpar][ithresh]    = val; }

  /**
   * @brief Set raw data mode
   * @param rd 1 - raw data mode, 0 - no raw data mode (stand-alone MC production)
   */
  void    SetRawData(Int_t rd)             { fGetRawData = rd; }

  /**
   * @brief Set STU region
   * @param rg STU region
   * 
   * The STU region is a bitmap of TRUs. In case a TRU is masked, the bit corresponding
   * to the TRU index is set to 0, otherwise it is set to 1 indicating that data from the 
   * TRU was used in the STU algorithms
   */
  void    SetRegion(Int_t rg)              { fRegion     = rg; }

  /**
   * @brief Set the firmware version
   * @param fv Firmware version
   */
  void    SetFw(Int_t fv)                  { fFw         = fv; }
  
  /**
   * @brief Get the scale factors calculating PHOS regions to EMCAL regions
   * @param iscale PHOS sector
   * @param val Scale factor
   */
  void    SetPHOSScale(int iscale, int val) { fPHOSScale[iscale] = val; }

  /**
   * @brief Set the TRU-STU error counters for a given timestamp
   * @param itru TRU ID (0-31: EMCAL, 32-45:DCAL)
   * @param itime Timestamp
   * @param errorcounts TRU-STU error communication errors
   */
  void    SetTRUErrorCounts(Int_t itru, Int_t itime, ULong64_t errorcounts);

  /**
   * @brief Set the patch size of the jet patch
   * @param size Patch size (0 - 8x8 patch, 1 - 16x16 patch)
   */
  void    SetPatchSize(Int_t size)         { fPatchSize =size; }

  /**
   * @brief Set median mode
   * @param mode 1 - median subtraction applied, 0 - median subtraction not applied
   */
  void    SetMedianMode(Int_t mode)        { fMedian =mode;    }
  
  /**
   * @brief Get gamma trigger threshold parameter
   * @param vzpar Index in the VZERO-dependent parabola equation (0 - A, 1 - B, 2 - C)
   * @param ithresh Threhold index (0 - high threshold, 1 - low threshold)
   * @return Threshold value for the given combination
   * 
   * In case a constant threshold is applied (all pp, p-Pb and run2 PbPb) the threshold
   * is stored in VZERO paramerer C.
   */
  Int_t   GetG(int vzpar, int ithresh) const { return fG[vzpar][ithresh];    }

  /**
   * @brief Get jet trigger threshold parameter
   * @param vzpar Index in the VZERO-dependent parabola equation (0 - A, 1 - B, 2 - C)
   * @param ithresh Threhold index (0 - high threshold, 1 - low threshold)
   * @return Threshold value for the given combination
   * 
   * In case a constant threshold is applied (all pp, p-Pb and run2 PbPb) the threshold
   * is stored in VZERO paramerer C.
   */
  Int_t   GetJ(int vzpar, int ithresh) const { return fJ[vzpar][ithresh];    }

  /**
   * @brief Get raw data mode
   * @return 1 - raw data mode, 0 - no raw data mode (stand-alone MC production)
   */
  Int_t   GetRawData()       const { return fGetRawData; }

  /**
   * @brief Get STU region size
   * @return STU region size
   * 
   * The STU region is a bitmap of TRUs. In case a TRU is masked, the bit corresponding
   * to the TRU index is set to 0, otherwise it is set to 1 indicating that data from the 
   * TRU was used in the STU algorithms
   */
  Int_t   GetRegion()        const { return fRegion;     }

  /**
   * @brief Get the Firmware version
   * @return Firmware version
   */
  Int_t   GetFw()            const { return fFw;         }
  
  /**
   * @brief Get the scale factors calculating PHOS regions to EMCAL regions
   * @param iscale PHOS sector
   * @return Scale factor
   */
  Int_t   GetPHOSScale(Int_t iscale) const { return fPHOSScale[iscale]; }

  /**
   * @brief Get the patch size of the jet patch
   * @return 0 - 8x8 patch, 1 - 16x16 patch 
   */
  Int_t   GetPatchSize()     const { return fPatchSize;  }

  /**
   * @brief Get the median mode
   * @return 1 - applied, 0 - not applied 
   */
  Int_t   GetMedianMode()    const { return fMedian;     }
  
  /**
   * @brief Get the STU segmentation
   * 
   * Segmentation is providing the size of the region and
   * subregion for the gamma/jet trigger. The vector has dimension
   * 2 with high threshold - 0, low threshold - 1
   * 
   * @param v1 Gamma subregion size
   * @param v2 Gamma number of subregions
   * @param v3 Jet subregion size
   * @param v4 Jet number of subregions
   */
  void    GetSegmentation(TVector2& v1, TVector2& v2, TVector2& v3, TVector2& v4) const;

  /**
   * @brief Get the TRU-STU error counters for a given TRU
   * @param itru TRU ID (0-31: EMCAL, 32-45:DCAL)
   * @return TRU-STU error counters
   */
  TClonesArray *GetErrorCountsForTRU(Int_t itru) const;

protected:
  AliEMCALTriggerSTUDCSConfig &operator=(const AliEMCALTriggerSTUDCSConfig &cd);

private:
	
  Int_t                   fG[3][2];                   ///< GA,B,C
  Int_t                   fJ[3][2];                   ///< JA,B,C
  Int_t                   fGetRawData;                ///< GetRawData
  Int_t                   fRegion;                    ///< Region
  Int_t                   fFw;                        ///< Firmware version
  Int_t                   fPHOSScale[4];              ///< PHOS scale factors
  Int_t                   fPatchSize;                 ///< Jet patch size: 0 for 8x8 and 2 for 16x16
  Int_t                   fMedian;                    ///< 1 in case of using EMCAL/DCAL for estimating the median
  TClonesArray            *fTRUErrorCounts[68];       ///< TRU error counts
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerSTUDCSConfig,5) ;
  /// \endcond
  
};

#endif // ALIEMCALTRIGGERSTUDCSCONFIG_H

