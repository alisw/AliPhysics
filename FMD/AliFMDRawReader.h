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
			      TArrayS &pulseLength);
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
  AliFMDRawReader(const AliFMDRawReader& o) 
    : TTask(o), 
      fTree(0), 
      fReader(0), 
      fData(0),
      fNbytes(0), 
      fSeen(0)
  {}
  AliFMDRawReader& operator=(const AliFMDRawReader&) { return *this; }
  ULong_t GetNwords() const {return fNbytes / 4;}
  UInt_t Get32bitWord(Int_t idx);
  Int_t GetHalfringIndex(UShort_t det, Char_t ring, UShort_t board);
  TTree*          fTree;       //! Pointer to tree to read into 
  AliRawReader*   fReader;     //! Pointer to raw reader 
  UShort_t        fSampleRate[3]; // The sample rate (if 0, inferred from data)
  UChar_t*        fData; 
  ULong_t  	  fNbytes; 
  Bool_t          fZeroSuppress[3];
  UShort_t        fNoiseFactor[3];
  AliFMDUShortMap fSeen;
  
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
