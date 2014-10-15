#ifndef ALIFMDRAWWRITER_H
#define ALIFMDRAWWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/* $Id$ */
/** @file    AliFMDRawWriter.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:45:56 2006
    @brief   Class to write raw data 
*/
//____________________________________________________________________
// 
// Class to writer ADC values to a Raw File
// Uses general ALTRO class - which is wrong
// Should make it right!
//
#ifndef ROOT_TTask
# include <TTask.h>
#endif

//____________________________________________________________________
class AliFMD;
class AliAltroBuffer;
class AliAltroBufferV3;
class TArrayI;
class TArrayF;
class TClonesArray;

//____________________________________________________________________
/** @class AliFMDRawWriter
    @brief Class to write ALTRO formated raw data from an array of
    AliFMDDigit objects.
    @code 
    AliFMDRawWriter* fmdWriter = new AliFMDRawWriter(0);
    TClonesArray*    array     = fmd->DigitArray();
    fmdWriter->WriteDigits(array);
    @endcode 
    @ingroup FMD_sim
*/
class AliFMDRawWriter : public TTask 
{
public:
  /** Constructor 
      @param fmd Pointer to detector object. */
  AliFMDRawWriter(AliFMD* fmd);
  /** Destructor */
  virtual ~AliFMDRawWriter() {}

  /** Write the output from AliFMD object. 
      @param option Not used */
  virtual void Exec(Option_t* option="");
  /** Write an array of AliFMDDigit objects as raw ALTRO data. 
      @param digits Array of AliFMDDigit objects to convert to raw
      ALTRO data. */
  virtual Long_t WriteDigits(TClonesArray* digits);
  /** Do zero-suppression of channel data. 
      @param data      Contain @a nWords of valid data.  On input, it 
                       contains the full channel data.  On output it
                       will contain the zero-suppresed data. 
      @param peds      Contain @a nWords pedestals 
      @param noise     Contain @a nWords pedestal widths 
      @param threshold Zero suppression threshold. */
  void ZeroSuppress(Int_t*& data, Int_t nWords, const Float_t* peds, 
		    const Float_t* noise, UShort_t threshold) const;
protected:
  AliFMDRawWriter(const AliFMDRawWriter& o) 
    : TTask(o), 
      fFMD(o.fFMD), 
      fSampleRate(o.fSampleRate), 
      fChannelsPerAltro(o.fChannelsPerAltro), 
      fThreshold(o.fThreshold)
  {}
  AliFMDRawWriter& operator=(const AliFMDRawWriter&) { return *this; }
  void WriteRCUTrailer(AliAltroBufferV3* altro,
		       UInt_t ddl,
		       Bool_t zs,
		       UShort_t factor,
		       UShort_t rate) const;
  AliFMD*       fFMD;              //! Pointer to detector description 
  UShort_t      fSampleRate;       // The sample rate (0 -> inferred from data)
  UShort_t      fChannelsPerAltro; // Number of pre-amp. channels/adc channel 
  UShort_t      fThreshold;        // Threshold for zero-suppression

  ClassDef(AliFMDRawWriter, 0) // Write FMD raw data to a DDL file
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
