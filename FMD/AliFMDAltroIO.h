#ifndef ALIFMDALTROIO_H
#define ALIFMDALTROIO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDAltroIO.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:27:31 2006
    @brief   ALTRO Input/output
*/
#include <iosfwd>
#include <TObject.h>

//____________________________________________________________________
/** @class AliFMDAltroIO AliFMDAltroIO.h <FMD/AliFMDAltroIO.h>
    @brief Base class for ALTRO Input/Output classes.
    @ingroup FMD_base
 */
class AliFMDAltroIO  : public TObject
{
public:
  /** Type of 40 bit words (signed) */
  typedef long long W40_t;
  /** Type of 10 bit words (signed) */
  typedef Int_t W10_t;
  /** Constructor */
  AliFMDAltroIO();
  /** Destructor */
  virtual ~AliFMDAltroIO() {}
  /** Error states */
  enum {
    /** No error */
    kNoError,
    /** Bad state after open/close file */
    kBadFile,
    /** Bad bit offset specified */
    kBadBits, 
    /** Bad state after reading from file */
    kBadRead, 
    /** Bad state after writing to file */
    kBadWrite, 
    /** Bad state after seeking in file */
    kBadSeek, 
    /** Could not tell position in file */
    kBadTell, 
    /** Bad trailer 40 bit word in file */
    kBadTrailer, 
    /** Bad fill word in file */
    kBadFill
  };
  /** Trailer mask */
  static const W40_t  fgkTrailerMask;
  /** Get error string */ 
  const char* ErrorString(Int_t err)  const;
protected:
  /** I/O Buffer */
  W40_t fBuffer;
  /** Pointer into buffer */
  Int_t fIBuffer;

  /** Concatenate a 10 bit word into a 40 bit word.
      @param n Offset (0-3)
      @param w 10 bit word
      @return @a w at offset @a n in a 40 bit word on success, a
      negative error code on failure. */
  virtual W40_t ConcatW40(UShort_t n, const W10_t& w) const;
  /** Extract a 10 bit word from a 40 bit word
      @param n The number 10bit word to extract (0-3)
      @param w 40 bit word to extract from. 
      @return The 10 bit word at @a n of @a w on success, or a
      negative error code otherwise. */
  virtual W10_t ExtractW10(UShort_t n, const W40_t w) const;

  ClassDef(AliFMDAltroIO,0);
};

//____________________________________________________________________
/** @class AliFMDAltroReader AliFMDAltroIO.h <FMD/AliFMDAltroIO.h>
    @brief Class to read ALTRO formated raw data from an AliRawReader
    object. 
    @code 
    AliRawReader*    reader    = new AliRawReaderFile(0);
    AliFMDRawReader* fmdReader = new AliFMDRawReader(reader);
    TClonesArray*    array     = new TClonesArray("AliFMDDigit");
    fmdReader->ReadAdcs(array);
    @endcode 
*/
class AliFMDAltroReader : public AliFMDAltroIO
{
public:
  /** Constructor 
      @param stream Stream to read from
      @exception Int_t A negative error code in case of failure */
  AliFMDAltroReader(std::istream& stream);
  virtual ~AliFMDAltroReader() {}
  /** Read one channel from the input file.  Note, that channels are
      read from the back of the file. 
      @param board   On return, the FEC board number 
      @param chip    On return, the ALTRO chip number
      @param channel On return, the ALTRO channel number
      @param last    On return, the size of the data    
      @param data    An array to fill with the data.  note, this
      should be large enough to hold all the data (1024 is the maximum
      number of timebins that can be read, so that's a safe size). 
      @return negative error code on failure, 0 if nothing is read, or
      the number of 10 bit words read. */
  Int_t ReadChannel(UShort_t& board, UShort_t& chip, UShort_t& channel, 
		    UShort_t& last,  UShort_t* data);
  /** Read one channel from the input file.  Note, that channels are
      read from the back of the file. 
      @param hwaddr On return, the hardware address 
      @param last   On return, the size of the data    
      @param data   An array to fill with the data.  note, this
      should be large enough to hold all the data (1024 is the maximum
      number of timebins that can be read, so that's a safe size). 
      @return negative error code on failure, 0 if nothing is read, or
      the number of 10 bit words read.  */
  Int_t ReadChannel(UShort_t& hwaddr, UShort_t& last, UShort_t* data);
  /** Extract the channel trailer. 
      @param hwaddr On return, the hardware address 
      @param last   On return, the size of the data    
      @return negative error code on failure, 0 if nothing is read, or
      the number of 10 bit words read. */
  Int_t ExtractTrailer(UShort_t& hwaddr, UShort_t& last);
  /** Extract bunches from data section of a channel.  
      @param last Pointer to last meaning full data entry. 
      @param data An array to fill with the read data. 
      @return negative error code on failure, otherwise number of 10
      bit words read.  */
  Int_t ExtractBunches(UShort_t last, UShort_t* data);
  /** Extract possible fill words. 
      @param last Pointer to last meaning full data entry. 
      @return Negative error code on failure, otherwise number of fill
      words read. */
  Int_t ExtractFillWords(UShort_t last);
  /** Extract bunch information from data. 
      @param data An array to fill with the read data. 
      @return negative error code on failure, otherwise number of 10
      bit words read.  */
  Int_t ExtractBunch(UShort_t* data);
  /** Check if @a x is a valid trailer
      @param x 40 bit word to check.
      @return @c true if @a x is a valid trailer */
  Bool_t IsTrailer(W40_t x);
  /** @return @c true if we're at the beginning of the file */
  Bool_t IsBof();
protected:
  /** Input stream */
  std::istream& fInput;
  /** Current position in file */
  // std::istream::pos_type 
  UShort_t fCurrent;
  /** High water mark */
  UShort_t fBegin;
  
  /** Read a 40 bit word from the input. 
      @return negative error code on failure, current position otherwise. */
  virtual Int_t ReadW40();
  /** Get a 10 bit word from the (buffered) input. 
      @return 10 bit word on success, negative error code on failure. */
  virtual W10_t GetNextW10();
  /** Get the next 40 bit word from the (buffered) input. 
      @return The 40 bit  word, or negative error code on failure */
  virtual W40_t GetNextW40();

  ClassDef(AliFMDAltroReader,0);
};

//____________________________________________________________________
/** @class AliFMDAltroWriter AliFMDAltroIO.h <FMD/AliFMDAltroIO.h>
    @brief Class to write ALTRO formated raw data from an array of
    AliFMDDigit objects.
    @code 
    AliFMDRawWriter* fmdWriter = new AliFMDRawWriter(0);
    TClonesArray*    array     = fmd->DigitArray();
    fmdWriter->WriteDigits(array);
    @endcode 
 */
class AliFMDAltroWriter : public AliFMDAltroIO
{
public:
  /** Constructor.
      @param stream File to read from
      @exception Int_t A negative error code in case of failure */
  AliFMDAltroWriter(std::ostream& stream);
  virtual ~AliFMDAltroWriter() {}
  /** @param threshold Zero-suppresion threshold */
  void SetThreshold(UShort_t threshold) { fThreshold = threshold; }
  /** Close the output, by writing the appropriate header.  The actual
      stream should be called by the user. 
      @return  number of bytes written, or negative error code on failure  */
  Int_t Close();
  /** Flush buffered output to file (if there is any). 
      @return  0, or negative error code on failure */
  Int_t Flush();
  /** Add a signal to output.   If the signal @a adc is less then the
      current threshold, a new bunch trailer is written. 
      @param adc Signal
      @return 0 on success, or negative error code on failure */
  Int_t AddSignal(UShort_t adc);
  /** Write a channel trailer to output. 
      @param board   The FEC board number (0-31)
      @param chip    The ALTRO chip number (0-7)
      @param channel The ALTRO channel number (0-16)
      @return Number of 10 bit words written, or negative error code
      on failure  */ 
  Int_t AddChannelTrailer(UShort_t board, UShort_t chip, UShort_t channel);
  /** Write a channel trailer to output.  
      @param hwaddr Hardware address (channel address)
      @return Number of 10 bit words written, or negative error code
      on failure */
  Int_t AddChannelTrailer(UInt_t hwaddr);
protected:
  /** Add a value to output buffer. 
      @param x Value to add. 
      @return number of 10 bit words written to disk, or negative
      error code on failure  */ 
  Int_t AddToBuffer(UShort_t x);
  /** Add a bunch trailer to output. 
      @return number of 10 bit words written to disk, or negative
      error code on failure  */ 
  Int_t AddBunchTrailer();
  /** Add fill words as needed to output
      @return number of 10 bit words written to disk, or negative
      error code on failure */
  Int_t AddFillWords();
  /** Zero suppression threshold */
  UShort_t fThreshold;
  /** Current time */
  UShort_t fTime;
  /** Current bunch length */
  UShort_t fLength;
  /** Last meaning-full data */
  UShort_t fLast;
  /** High-water mark (begining of file) */
  UShort_t fBegin;
  /** High-water mark (begining of file) */
  UShort_t fHeader;
  /** Total number of bytes written */
  Long_t fTotal;
  /** output stream */
  std::ostream& fOutput;

  ClassDef(AliFMDAltroWriter,0);
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
