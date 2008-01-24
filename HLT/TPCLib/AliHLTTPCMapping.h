// -*- Mode: C++ -*-
// $Id: AliHLTTPCDigitReaderRaw.cxx 22518 2007-11-26 23:19:47Z richterm $

#ifndef ALIHLTTPCMAPPING_H
#define ALIHLTTPCMAPPING_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTTPCMapping.h
    @author Kenneth Aamodt
    @date   
    @brief  Mapping class for TPC.
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTLogging.h"

/**
 * @class AliHLTTPCMapping
 * This is a mapping class for the TPC. It contains the mappping for all six partitions in static arrays.
 * This ensures that the asci files containing the mapping numbers is only read once per partition.
 * The only two methods interesting for the users are GetPad(hwaddress) and GetRow(hwaddress).
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
  UInt_t GetPad(UInt_t HWAddress);

  /**
   * Get the row number belonging to hardware address.
   * @param HWAddress    The hardware address of the given pad you are on.
   * @return Row number of hardware address (Pad).
   */
  UInt_t GetRow(UInt_t HWAddress);

 private:
  /** standard constructor prohibited, pad no always required */
  AliHLTTPCMapping();
  /** copy constructor prohibited */
  AliHLTTPCMapping(const AliHLTTPCMapping&);
  /** assignment operator prohibited */
  AliHLTTPCMapping& operator=(const AliHLTTPCMapping&);
  
  //Flags to check if mapping is done for the six patches
  static Bool_t fMapping0IsDone;
  static Bool_t fMapping1IsDone;
  static Bool_t fMapping2IsDone;
  static Bool_t fMapping3IsDone;
  static Bool_t fMapping4IsDone;
  static Bool_t fMapping5IsDone;


  /** size of mapping arrays */
  static const Int_t fgkMapping0Size=3200;                         // see above
  /** size of mapping array for patch 1 */
  static const Int_t fgkMapping1Size=3584;                         // see above
  /** size of mapping array for patch 2 */
  static const Int_t fgkMapping2Size=3200;                         // see above
  /** size of mapping array for patch 3 */
  static const Int_t fgkMapping3Size=3328;                         // see above
  /** size of mapping array for patch 4 */
  static const Int_t fgkMapping4Size=3328;                         // see above
  /** size of mapping array for patch 5 */
  static const Int_t fgkMapping5Size=3328;                         // see above

  /** mapping arrays for patch 0 */
  static UInt_t fgRowMapping0[fgkMapping0Size];                      // see above
  static UInt_t fgPadMapping0[fgkMapping0Size];                      // see above
  /** mapping arrays for patch 1 */
  static UInt_t fgRowMapping1[fgkMapping1Size];                      // see above
  static UInt_t fgPadMapping1[fgkMapping1Size];                      // see above
  /** mapping arrays for patch 2 */
  static UInt_t fgRowMapping2[fgkMapping2Size];                      // see above
  static UInt_t fgPadMapping2[fgkMapping2Size];                      // see above
  /** mapping arrays for patch 3 */
  static UInt_t fgRowMapping3[fgkMapping3Size];                      // see above
  static UInt_t fgPadMapping3[fgkMapping3Size];                      // see above
  /** mapping arrays for patch 4 */
  static UInt_t fgRowMapping4[fgkMapping4Size];                      // see above
  static UInt_t fgPadMapping4[fgkMapping4Size];                      // see above
  /** mapping arrays for patch 5 */
  static UInt_t fgRowMapping5[fgkMapping5Size];                      // see above
  static UInt_t fgPadMapping5[fgkMapping5Size];                      // see above

  /** current mapping array */
  UInt_t *fCurrentRowMapping;                                     //!transient
  UInt_t *fCurrentPadMapping;                                     //!transient

  /** Number of hardware adresses */
  UInt_t fNHWAdd;                                                 // see above

  /** Maximum number of hardware addresses */
  UInt_t fMaxHWAdd;                                               // see above

  ClassDef(AliHLTTPCMapping, 0)
};
#endif // ALIHLTTPCMAPPING_H
