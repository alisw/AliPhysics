#ifndef ALIGENREADERHEPMC_H
#define ALIGENREADERHEPMC_H
// Realisations of the AliGenReader interface to be used with AliGenExFile.
// NextEvent() loops over events
// and NextParticle() loops over particles.
// This implementation reads HepMC output formats
// Author: brian.peter.thorsbro@cern.ch, brian@thorsbro.dk
// Based on AliGenReaderSL by andreas.morsch@cern.ch

#include <TClonesArray.h>

#include "AliGenReader.h"
#include "AliGenEventHeader.h"
#include "THepMCParser.h"

namespace HepMC {
  class IO_BaseClass;
  class GenEvent;
}

class TParticle;

class AliGenReaderHepMC : public AliGenReader
{
public:
  AliGenReaderHepMC();
  AliGenReaderHepMC(const AliGenReaderHepMC &reader);
  virtual ~AliGenReaderHepMC();
   AliGenEventHeader * GetGenEventHeader() const {return fGenEventHeader;};
   virtual void Init();
   virtual Int_t NextEvent();
   virtual TParticle* NextParticle();
   virtual void RewindEvent();
   AliGenReaderHepMC & operator=(const AliGenReaderHepMC & rhs);

protected:
   HepMC::IO_BaseClass * fEventsHandle;   // pointer to the HepMC file handler
   HepMC::GenEvent * fGenEvent;           // pointer to a generated event
   TClonesArray * fParticleArray;         // pointer to array containing particles of current event
   TIter * fParticleIterator;             // iterator coupled to the array
   AliGenEventHeader * fGenEventHeader;   // AliGenEventHeader

private:
   void Copy(TObject&) const;

   ClassDef(AliGenReaderHepMC, 1) //Generate particles from external file
};
#endif



