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
#include "HepMC/IO_BaseClass.h"
#include "HepMC/GenEvent.h"

class TParticle;

class AliGenReaderHepMC : public AliGenReader
{
public:
   inline AliGenReaderHepMC():fEventsHandle(0), fGenEvent(0), fParticleArray(0), fParticleIterator(0), fGenEventHeader(0) {;}
   AliGenReaderHepMC(const AliGenReaderHepMC &reader)
   :AliGenReader(reader), fEventsHandle(0), fGenEvent(0), fParticleArray(0), fParticleIterator(0), fGenEventHeader(0) {reader.Copy(*this);}
   inline virtual ~AliGenReaderHepMC(){ delete fEventsHandle; delete fGenEvent; delete fParticleArray; delete fParticleIterator;} // not deleting fGenEventHeader as it is returned out
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



