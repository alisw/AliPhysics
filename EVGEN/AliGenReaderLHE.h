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

protected:
private:
  AliGenReaderLHE(const AliGenReaderLHE&);
  AliGenReaderLHE& operator=(const AliGenReaderLHE&);

  TXMLEngine        fXMLEngine;      //!
  XMLDocPointer_t   fXMLDoc;         //!
  XMLNodePointer_t  fXMLChild;       //!
  std::stringstream fStrStream;      //! contains event record
  std::streampos    fPosTracksBegin; //! used in RewindEvent() method
  TParticle         fParticle;       //!

  ClassDef(AliGenReaderLHE, 1);
} ;

#endif // ALIGENREADER_LHE_H
