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
  /**
   * Constructor 
   */
  AliFMDAltroMapping();
  /**
   * Destructor 
   */
  virtual ~AliFMDAltroMapping() {}
  /**
   * Return detector number corresponding to given DDL number 
   * 
   * @param ddl DDL number 
   * @return Detector number 
   */ 
  Short_t DDL2Detector(UInt_t ddl) const 
  { 
    return (ddl<=2 ? Short_t(ddl + 1) : -1); 
  }
  /**
   * Return the ring identifier corresponding to a board number 
   * 
   * @param board Board number 
   * @return Ring identifier 
   */ 
  Char_t Board2Ring(UShort_t board) const { return (board%2)?'O':'I'; }

  /**
   * Return the strip base number corresponding to a channel address 
   * 
   * @param board   Board number
   * @param altro   Altro number 
   * @param channel Channel number 
   * @param ring    On return, the ring ID 
   * @param sec     On return, the sector number 
   * @param strip   On return, the strip base offset 
   * @return @c true on success 
   */ 
  Bool_t Channel2StripBase(UShort_t  board, UShort_t  altro, 
			   UShort_t  chan,  Char_t&   ring, 
			   UShort_t& sec,   Short_t&  str) const;
  /**
   * Return the strip, sample corresponding to a timebin 
   * 
   * @param sec        Sector
   * @param timebin    Time bin 
   * @param preSamples Number of pre-samples 
   * @param sampleRate Oversampling rate 
   * @param strip      On return, the strip number in this channel
   * @param sam        On return, the sample number 
   */ 
  void Timebin2Strip(UShort_t sec,        UShort_t  timebin,
		     UShort_t preSamples, UShort_t  sampleRate, 
		     Short_t& strip,      UShort_t& sample) const;

  /**
   * Map a hardware address into a detector index. 
   * 
   * @param ddl        Hardware DDL number 
   * @param hwaddr     Hardware address.  
   * @param timebin    Timebin 
   * @param preSamples # of pre samples 
   * @param sampleRate Over sampling rate 
   * @param det        On return, the detector #
   * @param ring       On return, the ring ID
   * @param sec        On return, the sector #
   * @param str        On return, the base of strip #
   * @param sam        On return, the sample number for this strip
   * @return @c true on success, false otherwise 
   */
  Bool_t Hardware2Detector(UShort_t  ddl,        UShort_t hwaddr, 
			   UShort_t  timebin,    UShort_t preSamples, 
			   UShort_t  sampleRate,
			   UShort_t& det,        Char_t&   ring, 
			   UShort_t& sec,        Short_t&  str,
			   UShort_t& sam) const;
  /**
   * Map a hardware address into a detector index. 
   * 
   * @param ddl        Hardware DDL number 
   * @param board      FEC number
   * @param altro      ALTRO number 
   * @param channel    Channel number 
   * @param timebin    Timebin 
   * @param preSamples # of pre samples 
   * @param sampleRate Over sampling rate 
   * @param det        On return, the detector #
   * @param ring       On return, the ring ID
   * @param sec        On return, the sector #
   * @param str        On return, the base of strip #
   * @param sam        On return, the sample number for this strip
   * @return @c true on success, false otherwise 
   */
  Bool_t Hardware2Detector(UShort_t  ddl,        UShort_t  board, 
			   UShort_t  altro,      UShort_t  chan,
			   UShort_t  timebin,    UShort_t  preSamples,
			   UShort_t  sampleRate,
			   UShort_t& det,        Char_t&   ring, 
			   UShort_t& sec,        Short_t&  str,
			   UShort_t& sam) const;



  /**
   * Return DDL number corresponding to given detector number 
   * 
   * @param det Detector number 
   * @return DDL number 
   */ 
  UShort_t Detector2DDL(UShort_t det) const { return det - 1; }
  /**
   * Return board address corresponding to a sector 
   * 
   * @param ring  Ring identifier 
   * @param sec   Sector number 
   * @return The board number, or negative number in case of failure 
   */
  Short_t Sector2Board(Char_t ring, UShort_t sec) const;
    /**
   * Convert strip address to a channel address. 
   * 
   * @param ring  Ring identifier 
   * @param sec   Sector number 
   * @param str   Strip number 
   * @param board On return, contains the board number 
   * @param altro On return, contains the altro number 
   * @param chan  On return, contains the channel number 
   * @return @c true on success. 
   */
  Bool_t Strip2Channel(Char_t    ring,  UShort_t  sec,   
		       UShort_t  str,   UShort_t& board,
		       UShort_t& altro, UShort_t& chan) const;
  /**
   * Get the timebin correspoding to a strip and sample 
   * 
   * @param sec        Sector number 
   * @param str        Strip number 
   * @param sam        Sample number 
   * @param preSamples Number of pre-samples. 
   * @param sampleRate The over-sampling rate 
   * @return the timebin corresponding to the passed strip 
   */
  UShort_t Strip2Timebin(UShort_t sec, UShort_t strip, 
			 UShort_t sam, UShort_t preSamples, 
			 UShort_t sampleRate) const;
  
  /**
   * Map a detector index into a hardware address. 
   * 
   * @param det         The detector #
   * @param ring        The ring ID
   * @param sec         The sector #
   * @param str         The strip #
   * @param sam         The sample number 
   * @param preSamples  Number of pre-samples
   * @param sampleRate  The oversampling rate 
   * @param ddl         On return, hardware DDL number 
   * @param board       On return, the FEC board address (local to DDL)
   * @param altro       On return, the ALTRO number (local to FEC)
   * @param channel     On return, the channel number (local to ALTRO)
   * @param timebin     On return, the timebin number (local to ALTRO)
   * @return @c true on success, false otherwise 
   */
  Bool_t Detector2Hardware(UShort_t  det,        Char_t    ring, 
			   UShort_t  sec,        UShort_t  str,
			   UShort_t  sam, 
			   UShort_t  preSamples, UShort_t  sampleRate,
			   UShort_t& ddl,        UShort_t& board, 
			   UShort_t& altro,      UShort_t& channel, 
			   UShort_t& timebin) const;
  /**
   * Map a detector index into a hardware address. 
   * 
   * @param det         The detector #
   * @param ring        The ring ID
   * @param sec         The sector #
   * @param str         The strip #
   * @param sam         The sample number 
   * @param preSamples  Number of pre-samples
   * @param sampleRate  The oversampling rate 
   * @param ddl         On return, hardware DDL number 
   * @param hwaddr      On return, hardware address.  
   * @param timebin     On return, the timebin number (local to ALTRO)
   * @return @c true on success, false otherwise 
   */
  Bool_t Detector2Hardware(UShort_t  det,        Char_t    ring, 
			   UShort_t  sec,        UShort_t  str,
			   UShort_t  sam, 
			   UShort_t  preSamples, UShort_t  sampleRate,
			   UShort_t& ddl,        UShort_t& hwaddr, 
			   UShort_t& timebin) const;
  /**
   * Convert board, chip, channel to a hardware address 
   * 
   * @param board   Board number 
   * @param altro   Altro number 
   * @param channel Channel number 
   * @return hardware address of a channel 
   */ 
  UInt_t ChannelAddress(UShort_t board, UShort_t altro, UShort_t channel) const;
  /**
   * Convert a channel address to board, altro, channel fields 
   * 
   * @param hwaddr  Channel address
   * @param board   On return, the Board number 
   * @param altro   On return, the Altro number 
   * @param channel On return, the Channel number 
   */
  void ChannelAddress(UShort_t hwaddr, UShort_t& board, UShort_t& altro, 
		      UShort_t& channel) const;
  /**
   * convert a partial detector index into a hardware address
   * 
   * @param sector Sector number
   * @param str    Strip number
   * @param ring   Ring ID as an integer 
   * @return Hardware address 
   */
  Int_t  GetHWAddress(Int_t sector, Int_t str, Int_t ring);
  /**
   * Get the pad-row (or sector) corresponding to hardware address
   * 
   * @param hwaddr hardware address
   * @return Sector number 
   */
  Int_t  GetPadRow(Int_t hwaddr) const;
  /**
   * Get the pad (or strip) corresponding to hardware address
   * 
   * @param hwaddr hardware address
   * @return Strip number 
   */
  Int_t  GetPad(Int_t hwaddr) const;
  /**
   * Get the sector (or ring) corresponding to hardware address
   * 
   * @param hwaddr hardware address
   * @return Ring ID as an integer 
   */
  Int_t  GetSector(Int_t hwaddr) const;
  /**
   * Print map to standard out 
   * 
   * @param option Option string (hw, or det) 
   */
  void Print(Option_t* option="hw") const;
protected:
  /**
   * Read map from file - not used 
   * 
   * @return @c true on success 
   */
  virtual Bool_t ReadMapping();
  /**
   * Create the inverse mapping arrays 
   */ 
  virtual Bool_t CreateInvMapping();
  
  ClassDef(AliFMDAltroMapping, 2) // Read raw FMD Altro data 
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
