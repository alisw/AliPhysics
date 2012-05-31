// -*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTPCMAPPING_H
#define ALIHLTTPCMAPPING_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCMapping.h
    @author Kenneth Aamodt
    @date   
    @brief  Mapping class for TPC.
*/

#include "AliHLTLogging.h"

/**
 * @class AliHLTTPCMapping
 * This is a mapping class for the TPC. It contains the mappping for all six partitions in static arrays.
 * This ensures that the asci files containing the mapping numbers is only read once per partition.
 * The only two methods interesting for the users are GetPad(hwaddress) and GetRow(hwaddress).
 *
 * There are several possibilities to count the rows:
 * - A: absolute number: 0 to 158 over all partitions
 * - B: sectorwise: 0 to 62 for inner, 0 to 95 for outer sector
 * - C: within a partition
 *
 * This mappping class is designed to return the mapping within a partition (C), while the
 * mapping files use scheme B. The first row of each partition counted in scheme B has to
 * be subtracted.
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCMapping : public AliHLTLogging {
public:
  /** standard constructor */
  AliHLTTPCMapping(UInt_t patch);

  /** standard destructor */
  virtual ~AliHLTTPCMapping();
 
  /**
   * Create mapping for the given patch.
   */
  void InitializeMap(UInt_t patch);

  /**
   * Get the pad number belonging to hardware address.
   * @param HWAddress    The hardware address of the given pad
   * @return Pad number of given HWAddress
   */
  UInt_t GetPad(UInt_t HWAddress) const;

  /**
   * Get the row number belonging to hardware address.
   * @param HWAddress    The hardware address of the given pad you are on.
   * @return Row number of hardware address (Pad).
   */
  UInt_t GetRow(UInt_t HWAddress) const;

  /**
   * Get the HW addess (channel no) for a row/pad coordinate.
   * @param row          The row in the partition
   * @param pad          The pad in the row
   * @return hardware address (channel no).
   */
  UInt_t GetHwAddress(UInt_t row, UInt_t pad) const;

  /**
   * Read the mapping array from file.
   * @param patch           the patch (partition) to read
   * @param arraySize       size of the mapping arrays
   * @param rowArray        array of row mapping
   * @param padArray        array of pad mapping
   * @param hwaMappingSize  size of the HW address mapping
   * @param hwaArray        array of HW address mapping (backwards mapping)
   * @return kTRUE if successful
   */
  Bool_t ReadArray(UInt_t patch, UInt_t arraySize, UInt_t rowArray[], UInt_t padArray[], UInt_t hwaMappingSize, UInt_t hwaArray[]) const;
  
  /**
   * Checks if the hw address is valid
   * @param HWAddress    The hardware address of the pad
   * @return kTRUE if valid HWAddress 
   */
  Bool_t IsValidHWAddress(UInt_t HWAddress) const;

  /**
   * Method which returns the row offset of the patch. 
   * @return row offset of patch.
  */
  Int_t GetRowOffset() const {return fRowOffset;}

 private:
  /** standard constructor prohibited, pad no always required */
  AliHLTTPCMapping();
  /** copy constructor prohibited */
  AliHLTTPCMapping(const AliHLTTPCMapping&);
  /** assignment operator prohibited */
  AliHLTTPCMapping& operator=(const AliHLTTPCMapping&);

  /** the readout partition/patch */
  UInt_t fPatch;                                                     //! transient

  /** global number of patches */
  static const UChar_t fgkNofPatches=6;                              //! transient

  /** Flags to check if mapping is done for the six patches */
  static Bool_t fgMappingIsDone[fgkNofPatches];                      //! transient

  /** size of mapping arrays */
  static const UInt_t fgkMappingSize[fgkNofPatches];                 // see above

  /** array of the row mappings */
  static UInt_t* fgRowMapping[fgkNofPatches];                        //! transient

  /** array of the pad mappings */
  static UInt_t* fgPadMapping[fgkNofPatches];                        //! transient

  /** array of the HW address mappings */
  static UInt_t* fgHwaMapping[fgkNofPatches];                        //! transient

  /** size of mapping array for patch 0 */
  static const UInt_t fgkMapping0Size=3200;                          // see above
  /** size of mapping array for patch 1 */
  static const UInt_t fgkMapping1Size=3584;                          // see above
  /** size of mapping array for patch 2 */
  static const UInt_t fgkMapping2Size=3200;                          // see above
  /** size of mapping array for patch 3 */
  static const UInt_t fgkMapping3Size=3328;                          // see above
  /** size of mapping array for patch 4 */
  static const UInt_t fgkMapping4Size=3328;                          // see above
  /** size of mapping array for patch 5 */
  static const UInt_t fgkMapping5Size=3328;                          // see above

  /** name space for HW address mapping pad encoding */
  static const UInt_t fgkMappingHwaPadMask=0xff;                     //! transient
  /** bit shift for HW address mapping row encoding */
  static const UInt_t fgkMappingHwaRowMask=0x3f00;                   //! transient
  /** name space for HW address mapping row encoding */
  static const UChar_t fgkMappingHwaRowBitShift=8;                   //! transient
  /** size of row/pad to channel mapping
   * bit 0-8  encode pad
   * bit 9-11 encode row
   */
  static const UInt_t fgkMappingHwaSize=0x4000;                      // see above

  /** row mapping array for patch 0 */
  static UInt_t fgRowMapping0[fgkMapping0Size];                      // see above
  /** pad mapping array for patch 0 */
  static UInt_t fgPadMapping0[fgkMapping0Size];                      // see above
  /** hw address mapping array for patch 0 */
  static UInt_t fgHwaMapping0[fgkMappingHwaSize];                    // see above
  /** row mapping array for patch 1 */
  static UInt_t fgRowMapping1[fgkMapping1Size];                      // see above
  /** pad mapping array for patch 1 */
  static UInt_t fgPadMapping1[fgkMapping1Size];                      // see above
  /** hw address mapping array for patch 1 */
  static UInt_t fgHwaMapping1[fgkMappingHwaSize];                    // see above
  /** row mapping array for patch 2 */
  static UInt_t fgRowMapping2[fgkMapping2Size];                      // see above
  /** pad mapping array for patch 2 */
  static UInt_t fgPadMapping2[fgkMapping2Size];                      // see above
  /** hw address mapping array for patch 2 */
  static UInt_t fgHwaMapping2[fgkMappingHwaSize];                    // see above
  /** row mapping array for patch 3 */
  static UInt_t fgRowMapping3[fgkMapping3Size];                      // see above
  /** pad mapping array for patch 3 */
  static UInt_t fgPadMapping3[fgkMapping3Size];                      // see above
  /** hw addres mapping array for patch 3 */
  static UInt_t fgHwaMapping3[fgkMappingHwaSize];                    // see above
  /** row mapping array for patch 4 */
  static UInt_t fgRowMapping4[fgkMapping4Size];                      // see above
  /** pad mapping array for patch 4 */
  static UInt_t fgPadMapping4[fgkMapping4Size];                      // see above
  /** hw address mapping array for patch 4 */
  static UInt_t fgHwaMapping4[fgkMappingHwaSize];                    // see above
  /** row mapping array for patch 5 */
  static UInt_t fgRowMapping5[fgkMapping5Size];                      // see above
  /** pad mapping array for patch 5 */
  static UInt_t fgPadMapping5[fgkMapping5Size];                      // see above
  /** hw address mapping array for patch 5 */
  static UInt_t fgHwaMapping5[fgkMappingHwaSize];                    // see above

  /** current row mapping array */
  UInt_t *fCurrentRowMapping;                                        //!transient
  /** current pad mapping array */
  UInt_t *fCurrentPadMapping;                                        //!transient

  /** number of rows */
  Int_t fNofRows;                                                    // see above

  /** Maximum number of hardware addresses */
  UInt_t fMaxHWAdd;                                                  //!transient

  /** row offset according to scheme A (see class description) */
  Int_t fRowOffset;                                                  //!transient

  ClassDef(AliHLTTPCMapping, 0)
};
#endif // ALIHLTTPCMAPPING_H
