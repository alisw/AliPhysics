#ifndef ALIFMDALTROMAPPING_H
#define ALIFMDALTROMAPPING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDAltroMapping.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:28:11 2006
    @brief   Map HW address to detector coordinates and back again. 
*/
#ifndef ALIALTROMAPPING_H
# include <AliAltroMapping.h>
#endif
#ifndef ALIFMDUSHORTMAP_H
# include "AliFMDUShortMap.h"
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
//
// Map hardware address to detector coordinates. 
//
//    The hardware address consist of a DDL number and 12bits of ALTRO
//    addresses.  The ALTRO address are formatted as follows. 
//
//    12              7         4            0
//    |---------------|---------|------------|
//    | Board #       | ALTRO # | Channel #  |
//    +---------------+---------+------------+
//
//
//____________________________________________________________________
/** @class AliFMDAltroMapping 
    @brief Class that encodes a map to/from ALTRO hardware address to
    FMD detector coordinates.  
    
    The hardware address consist of a DDL number and 12bits of ALTRO
    addresses.  The ALTRO address are formatted as follows. 
    @verbatim 
    12              7         4            0
    |---------------|---------|------------|
    | Board #       | ALTRO # | Channel #  |
    +---------------+---------+------------+
    @endverbatim 

    @ingroup FMD_base
 */
class AliFMDAltroMapping : public AliAltroMapping
{
public:
  /** Constructor */
  AliFMDAltroMapping();
  /** Destructor */
  virtual ~AliFMDAltroMapping() {}
  /** Map a hardware address into a detector index. 
      @param ddl    Hardware DDL number 
      @param hwaddr Hardware address.  
      @param det    On return, the detector #
      @param ring   On return, the ring ID
      @param sec    On return, the sector #
      @param str    On return, the base of strip #
      @return @c true on success, false otherwise */
  Bool_t Hardware2Detector(UInt_t    ddl, UInt_t    hwaddr, 
			   UShort_t& det, Char_t&   ring, 
			   UShort_t& sec, UShort_t& str) const;
  /** Map a hardware address into a detector index. 
      @param ddl     Hardware DDL number 
      @param board   FEC number
      @param altro   ALTRO number 
      @param channel Channel number 
      @param det     On return, the detector #
      @param ring    On return, the ring ID
      @param sec     On return, the sector #
      @param str     On return, the base of strip #
      @return @c true on success, false otherwise */
  Bool_t Hardware2Detector(UInt_t    ddl,   UInt_t    board, 
			   UInt_t    altro, UInt_t    channel,
			   UShort_t& det,   Char_t&   ring, 
			   UShort_t& sec,   UShort_t& str) const;
  /** Map a detector index into a hardware address. 
      @param det     The detector #
      @param ring    The ring ID
      @param sec     The sector #
      @param str     The strip #
      @param ddl     On return, hardware DDL number 
      @param board   On return, the FEC board address (local to DDL)
      @param altro   On return, the ALTRO number (local to FEC)
      @param channel On return, the channel number (local to ALTRO)
      @return @c true on success, false otherwise */
  Bool_t Detector2Hardware(UShort_t  det,   Char_t    ring, 
			   UShort_t  sec,   UShort_t  str,
			   UInt_t&   ddl,   UInt_t&   board, 
			   UInt_t&   altro, UInt_t&   channel) const;
  /** Map a detector index into a hardware address. 
      @param det    The detector #
      @param ring   The ring ID
      @param sec    The sector #
      @param str    The strip #
      @param ddl    On return, hardware DDL number 
      @param hwaddr On return, hardware address.  
      @return @c true on success, false otherwise */
  Bool_t Detector2Hardware(UShort_t  det, Char_t    ring, 
			   UShort_t  sec, UShort_t  str,
			   UInt_t&   ddl, UInt_t&   hwaddr) const;
  /** Here to take care of a a misspelling in base class 
      @param sector Sector number
      @param str    Strip number
      @param ring   Ring ID as an integer 
      @return Hardware address */
  Int_t  GetHWAdress(Int_t sector, Int_t str, Int_t ring) const
  {
    return GetHWAddress(sector, str, ring);
  }
  /** convert a partial detector index into a hardware address
      @param sector Sector number
      @param str    Strip number
      @param ring   Ring ID as an integer 
      @return Hardware address */
  Int_t  GetHWAddress(Int_t sector, Int_t str, Int_t ring) const;
  /** Get the pad-row (or sector) corresponding to hardware address
      @param hwaddr hardware address
      @return Sector number */
  Int_t  GetPadRow(Int_t hwaddr) const;
  /** Get the pad (or strip) corresponding to hardware address
      @param hwaddr hardware address
      @return Strip number */
  Int_t  GetPad(Int_t hwaddr) const;
  /** Get the sector (or ring) corresponding to hardware address
      @param hwaddr hardware address
      @return Ring ID as an integer */
  Int_t  GetSector(Int_t hwaddr) const;
  /** Print map to standard out 
      @param option Option string (hw, or det) */
  void Print(Option_t* option="hw") const;
protected:
  /** Read map from file - not used 
      @return @c true on success */
  virtual Bool_t ReadMapping();
  /** Clear map in memory */
  virtual void   DeleteMappingArrays();
  
  ClassDef(AliFMDAltroMapping, 1) // Read raw FMD Altro data 
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
