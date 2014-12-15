#ifndef ALIFMDRAWSTREAM_H
#define ALIFMDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
//
// Class to read ALTRO formated data from an AliRawReader.  
// This class is mostly here to set AliAltroRawStream::fNoAltroMapping
// to false.   Furthermore, it defines the utility function
// ReadChannel to read in a full ALTRO channel.  The data is unpacked
// into the passed array.  
//
/** @file    AliFMDRawStream.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Mar 28 12:53:26 2006
    @brief   Class to read ALTRO formated data from an AliRawReader. 
*/
#ifndef ALIALTRORAWSTREAM_H
# include <AliAltroRawStream.h>
#endif 


/** @class AliFMDRawStream 
    @brief Class to read ALTRO formated data from an AliRawReader.  
    This class is mostly here to set
    AliAltroRawStream::fNoAltroMapping to false.   Furthermore, it
    defines the utility function ReadChannel to read in a full ALTRO
    channel.  The data is unpacked into the passed array. 
 */
class AliFMDRawStream : public AliAltroRawStream 
{
public:
  /** Constructor 
      @param reader Raw reader to use */
  AliFMDRawStream(AliRawReader* reader);
  /** Destructor  */
  virtual ~AliFMDRawStream() {}

  /** Read one ALTRO channel from the raw reader
      @param ddl  On return, the DDL
      @param addr On return, the hardware address
      @param len  On return, the number of entries filled in @a data
      @param data On return, the read ADC channels.
      @return @c true on success */
  virtual Bool_t ReadChannel(UInt_t& ddl, UInt_t& addr, 
			     UInt_t& len, volatile UShort_t* data);
protected:
  
  ClassDef(AliFMDRawStream, 0) // Read raw FMD Altro data 
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
