#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliGenReaderHepMC.h"
#include "AliRun.h"
#include "AliStack.h"




ClassImp(AliGenReaderHepMC)


AliGenReaderHepMC& AliGenReaderHepMC::operator=(const  AliGenReaderHepMC& rhs)
{
   // Assignment operator
   rhs.Copy(*this);
   return *this;
}

void AliGenReaderHepMC::Copy(TObject&) const
{
   //
   // Copy
   //
   Fatal("Copy","Not implemented!\n");
}

void AliGenReaderHepMC::Init()
{
   // Initialisation
   fEventsHandle = new HepMC::IO_GenEvent(fFileName, std::ios::in);
   fParticleArray = new TClonesArray("TParticle");
   fParticleIterator = new TIter(fParticleArray);



}

Int_t AliGenReaderHepMC::NextEvent()
{
   // Read the next event
   if ((fGenEvent = fEventsHandle->read_next_event())) {
      THepMCParser::ParseGenEvent2TCloneArray(fGenEvent,fParticleArray,false);
      fParticleIterator->Reset();


      // implement header... somewhere

      return fGenEvent->particles_size();
   }
   return 0;
}

TParticle* AliGenReaderHepMC::NextParticle()
{
   // Read next particle
   TParticle * particle = (TParticle*)fParticleIterator->Next();
   if (particle && particle->GetStatusCode()==1) {
      particle->SetBit(kTransportBit);
   }
   return particle;
}

void AliGenReaderHepMC::RewindEvent()
{
   fParticleIterator->Reset();
}
