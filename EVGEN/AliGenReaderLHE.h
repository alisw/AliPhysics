// -*- C++ -*-
#ifndef ALIGENREADER_LHE_H
#define ALIGENREADER_LHE_H

// Les Houches record format reader
// see, e.g.,
//  * https://arxiv.org/pdf/hep-ph/0609017.pdf
//  * https://arxiv.org/pdf/hep-ph/0109068.pdf

#include <sstream>

#include <TXMLEngine.h>
#include <TParticle.h>

#include "AliGenReader.h"

class AliGenReaderLHE : public AliGenReader
{
public:
  AliGenReaderLHE()
    : AliGenReader()
    , fConvIndicesFortranToC(kTRUE)
    , fXMLEngine()
    , fXMLDoc(NULL)
    , fXMLChild(NULL)
    , fStrStream()
    , fPosTracksBegin()
    , fParticle() {}

  virtual ~AliGenReaderLHE();

  virtual void       Init();
  virtual Int_t      NextEvent();
  virtual TParticle* NextParticle();
  virtual void       RewindEvent();

  void SetConvIndicesFortranToC(Bool_t b) { fConvIndicesFortranToC = b; }
  Bool_t GetConvIndicesFortranToC() const { return fConvIndicesFortranToC; }

protected:
private:
  AliGenReaderLHE(const AliGenReaderLHE&);
  AliGenReaderLHE& operator=(const AliGenReaderLHE&);

  Bool_t            fConvIndicesFortranToC; // if true convert mother indices from FORTRAN -> C counting convention

  TXMLEngine        fXMLEngine;      //!
  XMLDocPointer_t   fXMLDoc;         //!
  XMLNodePointer_t  fXMLChild;       //!
  std::stringstream fStrStream;      //! contains event record
  std::streampos    fPosTracksBegin; //! used in RewindEvent() method
  TParticle         fParticle;       //!

  ClassDef(AliGenReaderLHE, 2);
} ;

#endif // ALIGENREADER_LHE_H
