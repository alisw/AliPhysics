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

class AliRsnDaughter;

class AliRsnMiniParticle : public TObject {
public:

   AliRsnMiniParticle() : fIndex(-1), fCharge(0), fPDG(0), fMother(0), fMotherPDG(0), fCutBits(0x0) {Int_t i = 3; while (i--) fPsim[i] = fPrec[i] = 0.0;}

   Int_t&         Index()                    {return fIndex;}
   Char_t&        Charge()                   {return fCharge;}
   Float_t&       PsimX()                    {return fPsim[0];}
   Float_t&       PsimY()                    {return fPsim[1];}
   Float_t&       PsimZ()                    {return fPsim[2];}
   Float_t&       PrecX()                    {return fPrec[0];}
   Float_t&       PrecY()                    {return fPrec[1];}
   Float_t&       PrecZ()                    {return fPrec[2];}
   Float_t&       Px(Bool_t mc)              {return (mc ? fPsim[0] : fPrec[0]);}
   Float_t&       Py(Bool_t mc)              {return (mc ? fPsim[1] : fPrec[1]);}
   Float_t&       Pz(Bool_t mc)              {return (mc ? fPsim[2] : fPrec[2]);}
   Short_t&       PDG()                      {return fPDG;}
   Short_t        PDGAbs()                   {return TMath::Abs(fPDG);}
   Short_t&       Mother()                   {return fMother;}
   Short_t&       MotherPDG()                {return fMotherPDG;}
   UShort_t&      CutBits()                  {return fCutBits;}
   Bool_t         HasCutBit(Int_t i)         {UShort_t bit = 1 << i; return ((fCutBits & bit) != 0);}
   void           SetCutBit(Int_t i)         {UShort_t bit = 1 << i; fCutBits |=   bit;}
   void           ClearCutBit(Int_t i)       {UShort_t bit = 1 << i; fCutBits &= (~bit);}
   
   void           Set4Vector(TLorentzVector &v, Float_t mass, Bool_t mc) {v.SetXYZM(Px(mc), Py(mc), Pz(mc), mass);}
   void           CopyDaughter(AliRsnDaughter *daughter);

private:

   Int_t     fIndex;        // ID of track in its event
   Char_t    fCharge;       // track charge *character*: '+', '-', '0' (whatever else = undefined)
   Float_t   fPsim[3];      // MC momentum of the track
   Float_t   fPrec[3];      // reconstructed momentum of the track
   Short_t   fPDG;          // particle PDG code
   Short_t   fMother;       // index of mother in its container
   Short_t   fMotherPDG;    // PDG code of mother
   UShort_t  fCutBits;      // list of bits used to know what cuts were passed by this track
   
   ClassDef(AliRsnMiniParticle,2)
};

#endif
