//Simplified TParticle class


#include "AliHBTParticle.h"

#include <TParticle.h>

ClassImp(AliHBTParticle)

//______________________________________________________________________________
AliHBTParticle::AliHBTParticle():  
                fPdgCode(0), fPx(0), fPy(0),fPz(0),fE(0), fVx(0), fVy(0),fVz(0),fVt(0)
{//empty particle
}


//______________________________________________________________________________
AliHBTParticle::AliHBTParticle(Int_t pdg, Double_t px, Double_t py, Double_t pz, Double_t etot,
                     Double_t vx, Double_t vy, Double_t vz, Double_t time):
  fPdgCode(pdg), fPx(px), fPy(py),fPz(pz),fE(etot), 
                 fVx(vx), fVy(vy),fVz(vz),fVt(time)
{
//mormal constructor
  
  if (GetPDG()) {
     fCalcMass    = GetPDG()->Mass();
  } else {
     Double_t a2 = fE*fE -fPx*fPx -fPy*fPy -fPz*fPz;
     if (a2 >= 0) fCalcMass =  TMath::Sqrt(a2);
     else         fCalcMass = -TMath::Sqrt(-a2);
  }
}

//______________________________________________________________________________
AliHBTParticle::AliHBTParticle(const TParticle &p):
   fPdgCode(p.GetPdgCode()),fCalcMass(p.GetCalcMass()),
   fPx(p.Px()),fPy(p.Py()),fPz(p.Pz()),fE(p.Energy()), 
   fVx(p.Vx()),fVy(p.Vy()),fVz(p.Vz()),fVt(p.T())
{
 //all copied in the initialization
 
}

//______________________________________________________________________________
const Char_t* AliHBTParticle::GetName() const 
{
   static char def[4] = "XXX";
   const TParticlePDG *ap = TDatabasePDG::Instance()->GetParticle(fPdgCode);
   if (ap) return ap->GetName();
   else    return def;
}


//______________________________________________________________________________


