//Simplified TParticle class
#include "AliHBTParticle.h"
#include <TParticle.h>

ClassImp(AliHBTParticle)

Int_t AliHBTParticle::fgDebug = 0;
//______________________________________________________________________________
AliHBTParticle::AliHBTParticle():  
 fPdgIdx(0), fIdxInEvent(0),fNPids(0),fPids(0x0),fPidProb(0x0),
 fCalcMass(0),fPx(0), fPy(0),fPz(0),fE(0), fVx(0), fVy(0),fVz(0),fVt(0)
{//empty particle
}
//______________________________________________________________________________

AliHBTParticle::AliHBTParticle(Int_t pdg, Int_t idx,
               Double_t px, Double_t py, Double_t pz, Double_t etot,
               Double_t vx, Double_t vy, Double_t vz, Double_t time):  
  fPdgIdx(0), fIdxInEvent(0),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time)
{
//mormal constructor
  SetPdgCode(pdg);
  if (GetPDG()) {
     fCalcMass    = GetPDG()->Mass();
  } else {
     Double_t a2 = fE*fE -fPx*fPx -fPy*fPy -fPz*fPz;
     if (a2 >= 0) fCalcMass =  TMath::Sqrt(a2);
     else         fCalcMass = -TMath::Sqrt(-a2);
  }
}
//______________________________________________________________________________

AliHBTParticle::AliHBTParticle(Int_t pdg, Float_t prob, Int_t idx, 
                               Double_t px, Double_t py, Double_t pz, Double_t etot,
                               Double_t vx, Double_t vy, Double_t vz, Double_t time):
  fPdgIdx(0), fIdxInEvent(0),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time)
{
//mormal constructor
  SetPdgCode(pdg,prob);
  if (GetPDG()) {
     fCalcMass    = GetPDG()->Mass();
  } else {
     Double_t a2 = fE*fE -fPx*fPx -fPy*fPy -fPz*fPz;
     if (a2 >= 0) fCalcMass =  TMath::Sqrt(a2);
     else         fCalcMass = -TMath::Sqrt(-a2);
  }
}
//______________________________________________________________________________
AliHBTParticle::AliHBTParticle(const AliHBTParticle& in):
   fPdgIdx(in.fPdgIdx), fIdxInEvent(in.fIdxInEvent),
   fNPids(in.fNPids),fPids(new Int_t[fNPids]),fPidProb(new Float_t[fNPids]),
   fCalcMass(in.GetCalcMass()),
   fPx(in.Px()),fPy(in.Py()),fPz(in.Pz()),fE(in.Energy()), 
   fVx(in.Vx()),fVy(in.Vy()),fVz(in.Vz()),fVt(in.T())
{
 //Copy constructor
 for(Int_t i = 0; i<fNPids; i++)
  {
    fPids[i] =  in.fPids[i];
    fPidProb[i] = in.fPidProb[i];
  }
}

//______________________________________________________________________________
AliHBTParticle::AliHBTParticle(const TParticle &p,Int_t idx):
   fPdgIdx(0), fIdxInEvent(idx),
   fNPids(0),fPids(0x0),fPidProb(0x0),
   fCalcMass(p.GetCalcMass()),
   fPx(p.Px()),fPy(p.Py()),fPz(p.Pz()),fE(p.Energy()), 
   fVx(p.Vx()),fVy(p.Vy()),fVz(p.Vz()),fVt(p.T())
{
 //all copied in the initialization
 SetPdgCode(p.GetPdgCode());
}
//______________________________________________________________________________

void AliHBTParticle::SetPdgCode(Int_t pdg,Float_t prob)
{
  SetPIDprobability(pdg,prob);
  fPdgIdx = GetPidSlot(pdg);
}

//______________________________________________________________________________
void AliHBTParticle::SetPIDprobability(Int_t pdg, Float_t prob)
{
//Sets another pdg code and corresponding probabilty
//Ids are set in decreasing order
//Check if total prbaility is not ivercoming unity is performed
//in case, warning is printed
  if (fgDebug) Info("SetPIDprobability","Setting PID %d prob %f",pdg,prob);

  Float_t totprob = 0.0;//sums up probabilities
  Int_t idx = GetPidSlot(pdg);
  Int_t i;
  if (idx > -1) 
   {
     fPidProb[idx] = prob;
     for (i = 0; i < fNPids;i++) totprob+=fPidProb[i];
     if (totprob > 1.0)
       {
         Warning("SetPIDprobability","Total probability greater than UNITY: %f",totprob);
       }
     if (fgDebug) 
      {
        Info("SetPIDprobability","Current Total probability: %f",totprob);
      }
     return;
   }
    
  Int_t currentpid = GetPdgCode();
  fNPids++;
  Float_t* aPidProbNew = new Float_t[fNPids];
  Int_t* aPidsNew = new Int_t[fNPids];
  
  for (i = 0; i < fNPids-1;i++)//find a slot
   {
     if ( fPidProb[i] > prob)
      {
        if (fgDebug>4) Info("SetPID","Copying entry %d",i);
        aPidProbNew[i] = fPidProb[i];
        aPidsNew[i] = fPids[i];
        totprob+=fPidProb[i];
      }
     else break;
   }

  if (fgDebug>4) Info("SetPID","Setting new PID on entry %d",i);
  aPidProbNew[i] = prob;
  aPidsNew[i] = pdg;
  totprob+=prob;
  

  for (Int_t j = fNPids-1; j > i ;i--)//copy rest of old araays 
   {
     if (fgDebug>4) Info("SetPID","Copying from old entry %d to new entry %d",j-1,j);
     aPidProbNew[j] = fPidProb[j-1];
     aPidsNew[j] = fPids[j-1];
     totprob+=fPidProb[j-1];
   }

  if (totprob > 1.0)
   {
     Warning("SetId","Total probability is greater than 1 !!!");
     Print();
   }
  delete [] fPidProb;
  delete [] fPids;
  
  fPidProb = aPidProbNew;
  fPids = aPidsNew;
  
  fPdgIdx = GetPidSlot(currentpid);
  if (fPdgIdx == -1) fPdgIdx = 0;
}
//______________________________________________________________________________

Float_t AliHBTParticle::GetPIDprobability(Int_t pdg)
{
  Int_t idx = GetPidSlot(pdg);
  if (idx < 0) return 0.0;//such pid was not specified for this particle
  return fPidProb[idx];
}
//______________________________________________________________________________

const Char_t* AliHBTParticle::GetName() const 
{
   static char def[4] = "XXX";
   const TParticlePDG *ap = TDatabasePDG::Instance()->GetParticle(GetPdgCode());
   if (ap) return ap->GetName();
   else    return def;
}


//______________________________________________________________________________
Int_t AliHBTParticle::GetPidSlot(Int_t pdg) const
{
 //returns position of the given PID in fPids (and fPidProb) array.
 if (fPids == 0x0) return -1;
 for (Int_t i = 0; i< fNPids; i++)
  {
   if (fPids[i] == pdg) return i;
  }
 return -1;
}

//______________________________________________________________________________
void AliHBTParticle::Print() const
{
//prints information about particle
  printf("____________________________________________________\n");
  printf("Idx: %d  PID: %d  Name",fIdxInEvent,GetPdgCode());
  TParticlePDG *pdgp = TDatabasePDG::Instance()->GetParticle(GetPdgCode());
  if (pdgp)
   {
     printf("%s Mass %f\n",pdgp->GetName(),pdgp->Mass());
   }
  else
   {
     printf("Not known\n");
   }
  
  printf("Px: %+f Py: %+f Pz: %+f E: %+f Calculated Mass: %f Vx: %+f Vy: %+f Vz: %+f T: %+f",
          Px(),Py(),Pz(),Energy(),GetCalcMass(),Vx(),Vy(),Vz(),T());

  for (Int_t i = 0; i < fNPids; i++)
   {
     printf("# %d  PID: %d  Probability %f name",i,fPids[i],fPidProb[i]);
     const TParticlePDG *ap = TDatabasePDG::Instance()->GetParticle(fPids[i]);
     if (ap)
      {
        printf("%s Mass %f\n",ap->GetName(),ap->Mass());
      }
     else
      {
        printf("Not known\n");
      }
   }
}
