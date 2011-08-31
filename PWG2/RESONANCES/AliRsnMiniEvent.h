#ifndef ALIRSNMINIEVENT_H
#define ALIRSNMINIEVENT_H

//
// Mini-Event
// Contains only the useful quantities computed on the event
// which can be used for event mixing, or for direct output
// when doing analysis w.r. to multiplicity or event plane, for example.
//

#include <TArrayI.h>
#include <TClonesArray.h>

class AliRsnMiniParticle;

class AliRsnMiniEvent : public TObject {
public:

   AliRsnMiniEvent() : fID(-1), fVz(0.0), fMult(0.0), fAngle(0.0), fLeading(-1), fParticles("AliRsnMiniParticle", 0) {}
   ~AliRsnMiniEvent() {fParticles.Delete();}
   
   Int_t&              ID()        {return fID;}
   Float_t&            Vz()        {return fVz;}
   Float_t&            Mult()      {return fMult;}
   Float_t&            Angle()     {return fAngle;}
   TClonesArray&       Particles() {return fParticles;}
   Bool_t              IsEmpty()   {return fParticles.IsEmpty();}
   
   TArrayI             CountParticles(Char_t charge = 0, Int_t cutID = -1);
   AliRsnMiniParticle* GetParticle(Int_t i);
   AliRsnMiniParticle* LeadingParticle();
   void                AddParticle(AliRsnMiniParticle copy);
   
private:
   
   Int_t         fID;         // ID number
   Float_t       fVz;         // z-position of vertex
   Float_t       fMult;       // multiplicity or centrality
   Float_t       fAngle;      // angle of reaction plane to main reference frame
   
   Int_t         fLeading;    // index of leading particle
   TClonesArray  fParticles;  // list of selected particles
   
   ClassDef(AliRsnMiniEvent,2)
};

#endif
