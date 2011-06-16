#ifndef ALIRSNMINIPAIR_H
#define ALIRSNMINIPAIR_H

//
// This object is used as lightweight temporary container
// of all information needed from any input object and 
// useful for resonance analysis.
// 

#include <TObject.h>
#include <TLorentzVector.h>

class AliRsnMiniParticle;

class AliRsnMiniPair : public TObject {
public:

   AliRsnMiniPair() : fMother(-1), fMotherPDG(0) { }

   Int_t&          Mother()    {return fMother;}
   Int_t&          MotherPDG() {return fMotherPDG;}
   void            Fill(AliRsnMiniParticle *p1, AliRsnMiniParticle *p2, Double_t m1, Double_t m2, Double_t refMass);
   void            FillRef(Double_t mass);
   void            InvertP(Bool_t first);
      
   Int_t           ID(Bool_t mc) const {if (mc) return 1; else return 0;}
   
   TLorentzVector& P1 (Bool_t mc) {return fP1 [ID(mc)];}
   TLorentzVector& P2 (Bool_t mc) {return fP2 [ID(mc)];}
   TLorentzVector& Sum(Bool_t mc) {return fSum[ID(mc)];}
   TLorentzVector& Ref(Bool_t mc) {return fRef[ID(mc)];}
   
   Double_t        Pt(Bool_t mc)             const  {return fSum[ID(mc)].Pt();}
   Double_t        Pz(Bool_t mc)             const  {return fSum[ID(mc)].Pz();}
   Double_t        Eta(Bool_t mc)            const  {return fSum[ID(mc)].Eta();}
   Double_t        InvMass(Bool_t mc)        const  {return fSum[ID(mc)].M();}
   Double_t        InvMassRes()              const;            
   Double_t        Mt(Bool_t mc)             const  {return fRef[ID(mc)].Mt();}
   Double_t        Y(Bool_t mc)              const  {return fRef[ID(mc)].Rapidity();}
   Double_t        PtRatio(Bool_t mc)        const;
   Double_t        DipAngle(Bool_t mc)       const;
   Double_t        CosThetaStar(Bool_t mc);

private:

   TLorentzVector fP1 [2];    // 1st daughter momentum
   TLorentzVector fP2 [2];    // 2nd daughter momentum
   TLorentzVector fSum[2];    // sum of momenta
   TLorentzVector fRef[2];    // same as 'fSum' but with nominal resonance mass
                           
   Int_t          fMother;    // label of mothers (when common)
   Int_t          fMotherPDG; // PDG code of mother (when common)
   
   ClassDef(AliRsnMiniPair,1)
};

inline void AliRsnMiniPair::FillRef(Double_t mass)
{
//
// Fill ref 4-vectors using the passed mass and the values in 'sum'
//

   Int_t i;
   for (i = 0; i < 2; i++) {
      fRef[i].SetXYZM(fSum[i].X(), fSum[i].Y(), fSum[i].Z(), mass);
   }
}

inline Double_t AliRsnMiniPair::InvMassRes() const
{
//
// Return invariant mass resolution
//

   if (fSum[1].M() <= 0.0) return 1E20;
   
   return (fSum[0].M() - fSum[1].M()) / fSum[1].M();
}

inline Double_t AliRsnMiniPair::PtRatio(Bool_t mc) const
{
//
// Return ratio of transverse momenta of daughters
//

   Double_t num = TMath::Abs(fP1[ID(mc)].Perp() - fP2[ID(mc)].Perp());
   Double_t den = TMath::Abs(fP1[ID(mc)].Perp() + fP2[ID(mc)].Perp());
   
   if (den <= 0.0) return 1E20;
   
   return num / den;
}

inline Double_t AliRsnMiniPair::DipAngle(Bool_t mc) const
{
//
// Opening angle in a Z-T space
//

   const TLorentzVector &p1 = fP1[ID(mc)];
   const TLorentzVector &p2 = fP2[ID(mc)];
   
   return ((p1.Perp() * p2.Perp() + p1.Z() * p2.Z()) / p1.Mag() / p2.Mag());
}

#endif
