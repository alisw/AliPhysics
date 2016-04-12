#include <TVirtualMC.h>
#include <TDatabasePDG.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliGenReaderHepMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliGenHepMCEventHeader.h"

#include "HepMC/IO_BaseClass.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

ClassImp(AliGenReaderHepMC)

AliGenReaderHepMC::AliGenReaderHepMC():fEventsHandle(0), fGenEvent(0), fParticleArray(0), fParticleIterator(0), fGenEventHeader(0) {;}

AliGenReaderHepMC::AliGenReaderHepMC(const AliGenReaderHepMC &reader)
   :AliGenReader(reader), fEventsHandle(0), fGenEvent(0), fParticleArray(0), fParticleIterator(0), fGenEventHeader(0) {reader.Copy(*this);}


AliGenReaderHepMC& AliGenReaderHepMC::operator=(const  AliGenReaderHepMC& rhs)
{
   // Assignment operator
   rhs.Copy(*this);
   return *this;
}

AliGenReaderHepMC::~AliGenReaderHepMC(){ delete fEventsHandle; delete fGenEvent; delete fParticleArray; delete fParticleIterator;} // not deleting fGenEventHeader as it is returned out

void AliGenReaderHepMC::Copy(TObject&) const
{
   //
   // Copy
   //
   Fatal("Copy","Not implemented!\n");
}

void AliGenReaderHepMC::Init()
{
   // check if file exists
   if (access(fFileName, R_OK) != 0)
      AliError(Form("Couldn't open input file: %s", fFileName));

   // Initialisation
   fEventsHandle = new HepMC::IO_GenEvent(fFileName, std::ios::in);
   fParticleArray = new TClonesArray("TParticle");
   fParticleIterator = new TIter(fParticleArray);
}

Int_t AliGenReaderHepMC::NextEvent()
{
   // Clean memory
   if (fGenEvent) delete fGenEvent;
   // Read the next event
   if ((fGenEvent = fEventsHandle->read_next_event())) {
      THepMCParser::ParseGenEvent2TCloneArray(fGenEvent,fParticleArray,"GEV","CM",false);
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
      // propagate the event weight from HepMC to the event header
      HepMC::WeightContainer weights = fGenEvent->weights();
      if (!weights.empty())
        fGenEventHeader->SetEventWeight(weights.front());
      AliDebug(1, Form("Parsed event %d with %d particles, weight = %e", fGenEvent->event_number(), fGenEvent->particles_size(), fGenEventHeader->EventWeight()));
      return fGenEvent->particles_size();
   }
   AliError("No more events in the file.");
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
