#ifndef ALIFMDRAWSTREAM_H
#define ALIFMDRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIALTRORAWSTREAM_H
# include <AliAltroRawStream.h>
#endif 


// TPC to FMD translations 
// 
//    TPC                FMD
//    ----------+-----------
//    pad+time  |      strip
//    row       |     sector
//    sector    |       ring
// 
class AliFMDRawStream : public AliAltroRawStream 
{
public:
  AliFMDRawStream(AliRawReader* reader, UShort_t sampleRate=0);
  virtual ~AliFMDRawStream() {}

  Short_t Sector()      const { return fRow; }
  Char_t  Ring()        const { return (fSector == 0 ? 'I' : 'O'); }
  Short_t Strip()       const { return fPad + fTime / fSampleRate; }
  Short_t Sample()      const { return fTime % fSampleRate; }
  Short_t PrevSector()  const { return fPrevRow; }
  Char_t  PrevRing()    const { return (fPrevSector == 0 ? 'I' : 'O'); }
  Short_t PrevStrip()   const { return fPrevPad + fPrevTime/fSampleRate; }
    
  Bool_t  IsNewRing()   const { return (fSector != fPrevSector); }
  Bool_t  IsNewSector() const { return (fRow != fPrevRow) || IsNewRing(); }
  Bool_t  IsNewStrip()  const { return(Strip() != PrevStrip())||IsNewSector();}
    
  Short_t Count()       const { return fSignal; }
  Short_t SampleRate()  const { return fSampleRate; }
  
  virtual Bool_t Next();
  virtual Bool_t ReadChannel(UInt_t& addr, UInt_t& len, UShort_t* data);
  virtual Bool_t DumpData();
protected:
  virtual Int_t    ReadIntoBuffer();
  virtual Int_t    ReadTrailer(UInt_t& head, UInt_t& len);
  virtual Int_t    ReadFillWords(UInt_t len);
  virtual Int_t    ReadBunch(UShort_t* data);
  virtual UShort_t Get10BitWord();
  
  UShort_t  fSampleRate;         // # of ALTRO samples per VA1_ALICE clock
  Int_t     fPrevTime;           // Last time bin
  Bool_t    fExplicitSampleRate; // True if the sample rate was set externally
  Int_t     fPos;
  Int_t     fCur;
  UChar_t*  fRead;
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
