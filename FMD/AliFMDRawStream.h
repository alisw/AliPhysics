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
private:
  UShort_t fSampleRate; // # of ALTRO samples per VA1_ALICE clock
  Int_t    fPrevTime;   // Last time bin
public:
  AliFMDRawStream(AliRawReader* reader);

  Short_t Sector()      const { return fRow; }
  Char_t  Ring()        const { return (fSector == 0 ? 'I' : 'O'); }
  Short_t Strip()       const { return fPad + fTime / fSampleRate; }
  Short_t PrevSector()  const { return fPrevRow; }
  Char_t  PrevRing()    const { return (fPrevSector == 0 ? 'I' : 'O'); }
  Short_t PrevStrip()   const { return fPrevPad + fPrevTime/fSampleRate; }

  Bool_t  IsNewRing()   const { return (fSector != fPrevSector); }
  Bool_t  IsNewSector() const { return (fRow != fPrevRow) || IsNewRing(); }
  Bool_t  IsNewStrip()  const { return(Strip() != PrevStrip())||IsNewSector();}
    
  Short_t Count()      const { return fSignal; }
  Short_t SampleRate() const { return fSampleRate; }
  
  virtual Bool_t   Next();
  
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
