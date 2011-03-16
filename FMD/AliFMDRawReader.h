#ifndef ALIFMDRAWREADER_H
#define ALIFMDRAWREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
// 
// Class to read ADC values from a AliRawReader object. 
// Note, that it uses an ALTRO reader, which is wrong. 
// Perhaps we need to implement it our selves
// 
/* $Id$ */
/** @file    AliFMDRawReader.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:45:23 2006
    @brief   Class to read raw data 
    @ingroup FMD_rec
*/
#ifndef ROOT_TTask
# include <TTask.h>
#endif
#include "AliFMDUShortMap.h"

//____________________________________________________________________
class AliRawReader;
class AliAltroRawStreamV3;
class TTree;
class TClonesArray;
class TArrayS;
class AliFMDCalibSampleRate;
class AliFMDCalibStripRange;

//____________________________________________________________________
/** @brief Class to read ALTRO formated raw data from an AliRawReader
    object. 
    @code 
    AliRawReader*    reader    = new AliRawReaderFile(0);
    AliFMDRawReader* fmdReader = new AliFMDRawReader(reader);
    TClonesArray*    array     = new TClonesArray("AliFMDDigit");
    fmdReader->ReadAdcs(array);
    @endcode 
    @ingroup FMD_rec
*/
class AliFMDRawReader : public TTask 
{
public:
  /** 
   * CTOR 
   *
   * @param reader Raw reader
   * @param array  Output tree 
   */
  AliFMDRawReader(AliRawReader* reader, TTree* array);
  /** 
   * DTOR 
   */
  virtual ~AliFMDRawReader() {}
  /** 
   * Read in, and store in output tree 
   *
   * @param option Not used 
   */
  virtual void   Exec(Option_t* option="");
  /**
   * Read ADC's into a TClonesArray of AliFMDDigit objects. 
   *
   * @param array       Array to read into 
   * @param summable    Create SDigits rather than digits 
   * @param pedSub      Whether to do pedestal subtraction.
   * @param noiseFactor If we do pedestal subtraction, then this is
   *        the number we use to suppress remenants of the noise.
   * 
   * @return @c true on success 
   */
  virtual Bool_t ReadAdcs(TClonesArray* array);
  /** 
   * Read SOD event into passed objects.
   * 
   * @param samplerate   The sample rate object to fill
   * @param striprange   The strip range object to fill
   * @param pulseSize    The pulse size object to fill
   * @param pulseLength  The pulse length (in events) object to fill
   * 
   * @return @c true on success
   */  
  virtual Bool_t ReadSODevent(AliFMDCalibSampleRate* samplerate, 
			      AliFMDCalibStripRange* striprange, 
			      TArrayS &pulseSize, 
			      TArrayS &pulseLength, 
			      Bool_t* detectors=0);
  /** 
   * Check of the data from DDL @a ddl is zero-suppressed
   * 
   * @param ddl DDL number (0-2)
   * 
   * @return @c true if the data from this DDL is zero-suppressed. 
   */  
  Bool_t IsZeroSuppressed(UShort_t ddl) const { return fZeroSuppress[ddl]; }
  /** 
   * The factor used to multiply the noise when making on-line
   * pedestal subtraction.
   * 
   * @param ddl DDL number (0-2)
   * 
   * @return The factor used. 
   */
  UShort_t NoiseFactor(UShort_t ddl) const { return fNoiseFactor[ddl]; }

  /** 
   * Get the next signal
   * 
   * @param det  On return, the detector
   * @param rng  On return, the ring
   * @param sec  On return, the sector
   * @param str  On return, the strip
   * @param sam  On return, the sample
   * @param rat  On return, the sample rate
   * @param adc  On return, the ADC value
   * @param zs   On return, whether zero-supp. is enabled
   * @param fac  On return, the usd noise factor
   * 
   * @return true if valid data is returned
   */
  Bool_t NextSample(UShort_t& det, Char_t&   rng, UShort_t& sec, UShort_t& str,
		    UShort_t& sam, UShort_t& rat, Short_t&  adc, 
		    Bool_t&   zs,  UShort_t& fac);
  /** 
   * Get the next signal
   * 
   * @param det  On return, the detector
   * @param rng  On return, the ring
   * @param sec  On return, the sector
   * @param str  On return, the strip
   * @param adc  On return, the ADC value
   * @param zs   On return, whether zero-supp. is enabled
   * @param fac  On return, the usd noise factor
   * 
   * @return true if valid data is returned
   */
  Bool_t NextSignal(UShort_t& det, Char_t&   rng, 
		    UShort_t& sec, UShort_t& str, 
		    Short_t&  adc, Bool_t&   zs, 
		    UShort_t& fac);
  /** 
   * Whether to keep a sample based on the rate used. 
   * 
   * @param samp Sample number 
   * @param rate Over sampling rate
   * 
   * @return Whether to keep the sample or not
   */
  static Bool_t SelectSample(UShort_t samp, UShort_t rate);
protected:
  /** 
   * Copy constructor 
   * 
   * @param o Object to construct from
   */  
  AliFMDRawReader(const AliFMDRawReader& o) 
    : TTask(o), 
      fTree(0), 
      fReader(0), 
      fData(0),
      fNbytes(0), 
      fMinStrip(0), 
      fMaxStrip(127),
      fPreSamp(14+5),
      fSeen(0)
  {}
  /** 
   * Assignment operator
   * 
   * @return Reference to this object
   */
  AliFMDRawReader& operator=(const AliFMDRawReader&) { return *this; }
  /** 
   * Process a new DDL.  Sets the internal data members fZeroSuppress, 
   * fSampleRate, and fNoiseFactor based on information in the RCU trailer. 
   * 
   * @param input Input stream
   * @param det   On return, the detector number
   * 
   * @return negative value in case of problems, the DDL number otherwise
   */
  Int_t NewDDL(AliAltroRawStreamV3& input, UShort_t& det);
  /** 
   * Processs a new channel.  Sets the internal data members
   * fMinStrip, fMaxStrip, and fPreSamp. 
   * 
   * @param input   Input stream
   * @param det     Detector number
   * @param ring    On return, the ring identifier 
   * @param sec     On return, the sector number
   * @param strbase On return, the strip base
   * 
   * @return negative value in case of problems, hardware address otherwise
   */
  Int_t NewChannel(const AliAltroRawStreamV3& input, 
		   UShort_t det, Char_t&  ring, 
		   UShort_t& sec, Short_t& strbase);

  /** 
   * Process a new timebin
   * 
   * @param input   Input stream
   * @param i       Index into bunch data
   * @param t       Time
   * @param sec     Sector number
   * @param strbase Base of strip numbers for this channel
   * @param str     On return, the strip number
   * @param samp    On return, the sample number
   * 
   * @return negative value in case of problems, ADC value otherwise
   */  
  Int_t NewSample(const AliAltroRawStreamV3& input, 
		  Int_t i, UShort_t t, UShort_t sec,
		  UShort_t  strbase, Short_t&  str, UShort_t& samp);

  /** 
   * 
   * Get the number of words 
   * 
   * @return Number of 32bit words 
   */
  ULong_t GetNwords() const {return fNbytes / 4;}
  /** 
   * Get the next 32bit word from payload
   * 
   * @param idx Which 32bit word to get
   * 
   * @return 
   */
  UInt_t Get32bitWord(Int_t idx);
  /** 
   * Get short index for a given half-ring
   * 
   * @param det   Detector number
   * @param ring  Ring identifer
   * @param board Board number
   * 
   * @return 
   */
  Int_t GetHalfringIndex(UShort_t det, Char_t ring, UShort_t board) const;
  TTree*          fTree;            //! Pointer to tree to read into 
  AliRawReader*   fReader;          //! Pointer to raw reader 
  UShort_t        fSampleRate[3];   // The sample rate (if 0,inferred from data)
  UChar_t*        fData;            // Data pointer
  ULong_t  	  fNbytes;          // Number of bytes
  Bool_t          fZeroSuppress[3]; // Zero suppression flag
  UShort_t        fNoiseFactor[3];  // Noise factor 
  UShort_t        fMinStrip;        // Current minimum strip number (0)
  UShort_t        fMaxStrip;        // Current maximum strip number (127)
  UShort_t        fPreSamp;         // Current number of pre-samples (14+5)
  AliFMDUShortMap fSeen;            // Seen strips 
  
  ClassDef(AliFMDRawReader, 0) // Read FMD raw data into a cache 
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

