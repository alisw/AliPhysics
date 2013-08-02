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

   Int_t          &Mother()    {return fMother;}
   Int_t          &MotherPDG() {return fMotherPDG;}
   void            Fill(AliRsnMiniParticle *p1, AliRsnMiniParticle *p2, Double_t m1, Double_t m2, Double_t refMass);
   void            FillRef(Double_t mass);
   void            InvertP(Bool_t first);

   Int_t           ID(Bool_t mc) const {if (mc) return 1; else return 0;}

   TLorentzVector &P1 (Bool_t mc) {return fP1 [ID(mc)];}
   TLorentzVector &P2 (Bool_t mc) {return fP2 [ID(mc)];}
   TLorentzVector &Sum(Bool_t mc) {return fSum[ID(mc)];}
   TLorentzVector &Ref(Bool_t mc) {return fRef[ID(mc)];}

   Double_t        Pt(Bool_t mc)             const  {return fSum[ID(mc)].Pt();}
   Double_t        Pz(Bool_t mc)             const  {return fSum[ID(mc)].Pz();}
   Double_t        Eta(Bool_t mc)            const  {return fSum[ID(mc)].Eta();}
   Double_t        InvMass(Bool_t mc)        const  {return fSum[ID(mc)].M();}
   Double_t        InvMassRes()              const;
   Double_t        InvMassDiff()             const;
   Double_t        Mt(Bool_t mc)             const  {return fRef[ID(mc)].Mt();}
   Double_t        Y(Bool_t mc)              const  {return fRef[ID(mc)].Rapidity();}
   Double_t        PtRatio(Bool_t mc)        const;
   Double_t        DipAngle(Bool_t mc)       const;
   Double_t        CosThetaStar(Bool_t mc);
   Double_t        DaughterPt(Int_t daughterId, Bool_t mc); 
   void            DaughterPxPyPz(Int_t daughterId, Bool_t mc, Double_t *pxpypz); 
      
 private:
   
   TLorentzVector fP1 [2];    // 1st daughter momentum
   TLorentzVector fP2 [2];    // 2nd daughter momentum
   TLorentzVector fSum[2];    // sum of momenta
   TLorentzVector fRef[2];    // same as 'fSum' but with nominal resonance mass
   
   Int_t          fMother;    // label of mothers (when common)
   Int_t          fMotherPDG; // PDG code of mother (when common)
   
   ClassDef(AliRsnMiniPair,1)
     };

#endif
