#ifndef ALIRSNMINIPARTICLE_H
#define ALIRSNMINIPARTICLE_H

//
// This object is used as lightweight temporary container
// of all information needed from any input object and
// useful for resonance analysis.
//

#include <TMath.h>
#include <TObject.h>
#include <TLorentzVector.h>

#include "AliESDtrack.h"

class AliRsnDaughter;

class AliRsnMiniParticle : public TObject {
public:

   AliRsnMiniParticle() : fIndex(-1), fCharge(0), fPDG(0), fMother(0), fMotherPDG(0), fDCA(0), fNTotSisters(0), fIsFromB(kFALSE), fIsQuarkFound(kFALSE), fCutBits(0x0), fPassesOOBPileupCut(kTRUE) {
       Int_t i = 3; while (i--) fPsim[i] = fPrec[i] = fPmother[i] = 0.0;
       fIndexDaughters[0] = fIndexDaughters[1] = fIndexDaughters[2] = -1;
       fMass[0] = fMass[1] = -1.0;
   }

   Int_t         &Index()                    {return fIndex;}
   void          SetResonance()              {fIndex=-99999;}
   Bool_t        IsResonance()               {return (fIndex==-99999);}
   Int_t         &IndexV0Pos()               {return fIndexDaughters[0];}
   Int_t         &IndexV0Neg()               {return fIndexDaughters[1];}
   Int_t         &IndexBachelor()            {return fIndexDaughters[2];}
   Char_t        &Charge()                   {return fCharge;}
   Float_t       &PsimX()                    {return fPsim[0];}
   Float_t       &PsimY()                    {return fPsim[1];}
   Float_t       &PsimZ()                    {return fPsim[2];}
   Float_t       &PrecX()                    {return fPrec[0];}
   Float_t       &PrecY()                    {return fPrec[1];}
   Float_t       &PrecZ()                    {return fPrec[2];}
   Float_t       &PmotherX()                 {return fPmother[0];}
   Float_t       &PmotherY()                 {return fPmother[1];}
   Float_t       &PmotherZ()                 {return fPmother[2];}
   Float_t       &Px(Bool_t mc)              {return (mc ? fPsim[0] : fPrec[0]);}
   Float_t       &Py(Bool_t mc)              {return (mc ? fPsim[1] : fPrec[1]);}
   Float_t       &Pz(Bool_t mc)              {return (mc ? fPsim[2] : fPrec[2]);}
   Long_t        &PDG()                      {return fPDG;}
   Long_t        PDGAbs()                    {return TMath::Abs(fPDG);}
   Double_t      Mass();                     // returns PDG mass
   Double_t      &StoredMass(Bool_t mc)      {return (mc ? fMass[0] : fMass[1]);} // store mass for resonances
   Int_t         &Mother()                   {return fMother;}
   Long_t        &MotherPDG()                {return fMotherPDG;}
   Bool_t        &IsFromB()                  {return fIsFromB;}
   Bool_t        &IsQuarkFound()             {return fIsQuarkFound;}
   UShort_t      &CutBits()                  {return fCutBits;}
   Double_t      DCA()                      {return fDCA;}
   Short_t       NTotSisters()              {return fNTotSisters;}
   Bool_t        HasCutBit(Int_t i)         {UShort_t bit = 1 << i; return ((fCutBits & bit) != 0);}
   void          SetCutBit(Int_t i)         {UShort_t bit = 1 << i; fCutBits |=   bit;}
   void          ClearCutBit(Int_t i)       {UShort_t bit = 1 << i; fCutBits &= (~bit);}
   Bool_t        &PassesOOBPileupCut()      {return fPassesOOBPileupCut;}

   void          Clear(Option_t *opt="");

   void          Set4Vector(TLorentzVector &v, Float_t mass=-1.0, Bool_t mc=kFALSE);
   void          CopyDaughter(AliRsnDaughter *daughter);

private:
    
   Bool_t    TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b);

   Int_t     fIndex;        // ID of track in its event
   Int_t     fIndexDaughters[3]; // daugher indices (0:pos, 1:neg, 2: bachelor) (if not V0/resonance then 0:ESD/AOD label, 1:(-1))
   Char_t    fCharge;       // track charge *character*: '+', '-', '0' (whatever else = undefined)
   Float_t   fPsim[3];      // MC momentum of the track
   Float_t   fPrec[3];      // reconstructed momentum of the track
   Double_t  fMass[2];      // reconstructed [0] and simulated [1] mass, only used for resonances
   Float_t   fPmother[3];   // MC momentum of the track's mother
   Long_t    fPDG;          // particle PDG code
   Int_t     fMother;       // index of mother in its container
   Long_t    fMotherPDG;    // PDG code of mother
   Double_t  fDCA;          // DCA of the particle
   Short_t   fNTotSisters;  // number of  daughters of the particle
   Bool_t    fIsFromB;	    // is the particle from B meson flag
   Bool_t    fIsQuarkFound; // is the particle from a quark flag (used to reject or accept Hijing event)
   UShort_t  fCutBits;      // list of bits used to know what cuts were passed by this track
   Bool_t    fPassesOOBPileupCut; // passes out-of-bunch pileup cut

   ClassDef(AliRsnMiniParticle, 9)
};

#endif
