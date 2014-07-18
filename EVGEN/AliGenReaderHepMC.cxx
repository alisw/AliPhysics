#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliGenReaderHepMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliGenHepMCEventHeader.h"


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
   // check if file exists, using FILE to avoid (the otherwise faster) POSIX dependencies
   if (FILE *file = fopen(fFileName,"r"))  {
      printf("File %s opened \n", fFileName);
      fclose(file);
   } else {
      printf("Couldn't open input file: %s \n", fFileName);
   }
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
      THepMCParser::HeavyIonHeader_t heavyIonHeader;
      THepMCParser::PdfHeader_t pdfHeader;
      THepMCParser::ParseGenEvent2HeaderStructs(fGenEvent,heavyIonHeader,pdfHeader,true,true);
      fGenEventHeader = new AliGenHepMCEventHeader(
            heavyIonHeader.Ncoll_hard,
            heavyIonHeader.Npart_proj,
            heavyIonHeader.Npart_targ,
            heavyIonHeader.Ncoll,
            heavyIonHeader.spectator_neutrons,
            heavyIonHeader.spectator_protons,
            heavyIonHeader.N_Nwounded_collisions,
            heavyIonHeader.Nwounded_N_collisions,
            heavyIonHeader.Nwounded_Nwounded_collisions,
            heavyIonHeader.impact_parameter,
            heavyIonHeader.event_plane_angle,
            heavyIonHeader.eccentricity,
            heavyIonHeader.sigma_inel_NN,
            pdfHeader.id1,
            pdfHeader.id2,
            pdfHeader.pdf_id1,
            pdfHeader.pdf_id2,
            pdfHeader.x1,
            pdfHeader.x2,
            pdfHeader.scalePDF,
            pdfHeader.pdf1,
            pdfHeader.pdf2
      );
      printf("Parsed event with %d particles.\n", fGenEvent->particles_size());
      return fGenEvent->particles_size();
   }
   printf("No more events in the file.\n");
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
