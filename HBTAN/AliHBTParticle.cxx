#include "AliHBTParticle.h"
//___________________________________________________________
/////////////////////////////////////////////////////////////
//
// class AliHBTParticle
//
// Ali HBT Particle: simplified class TParticle
// Simplified in order to minimize the size of object
//  - we want to keep a lot of such a objects in memory
// Additionaly adjusted for HBT Analysies purposes
// + pointer to Track Points
// + pointer to Cluster Map(s)
//
// Piotr.Skowronski@cern.ch
//
/////////////////////////////////////////////////////////////
#include <TParticle.h>
#include "AliHBTTrackPoints.h"
#include "AliHBTClusterMap.h"

ClassImp(AliHBTParticle)

Int_t AliHBTParticle::fgDebug = 0;
//______________________________________________________________________________
AliHBTParticle::AliHBTParticle():  
 fPdgIdx(0), fIdxInEvent(0),fNPids(0),fPids(0x0),fPidProb(0x0),
 fCalcMass(0),fPx(0), fPy(0),fPz(0),fE(0), fVx(0), fVy(0),fVz(0),fVt(0),
 fTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
{//empty particle
}
//______________________________________________________________________________

AliHBTParticle::AliHBTParticle(Int_t pdg, Int_t idx,
               Double_t px, Double_t py, Double_t pz, Double_t etot,
               Double_t vx, Double_t vy, Double_t vz, Double_t time):  
  fPdgIdx(0), fIdxInEvent(idx),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time),
  fTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
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
  fPdgIdx(0), fIdxInEvent(idx),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time),
  fTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
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
   TObject(in),
   fPdgIdx(in.fPdgIdx), fIdxInEvent(in.fIdxInEvent),
   fNPids(in.fNPids),fPids(new Int_t[fNPids]),fPidProb(new Float_t[fNPids]),
   fCalcMass(in.GetCalcMass()),
   fPx(in.Px()),fPy(in.Py()),fPz(in.Pz()),fE(in.Energy()), 
   fVx(in.Vx()),fVy(in.Vy()),fVz(in.Vz()),fVt(in.T()),
   fTrackPoints(0x0), fITSTrackPoints(0x0), fClusterMap(0x0)
{
 //Copy constructor
 for(Int_t i = 0; i<fNPids; i++)
  {
    fPids[i] =  in.fPids[i];
    fPidProb[i] = in.fPidProb[i];
  }
 
 if (in.fTrackPoints)
   fTrackPoints = (AliHBTTrackPoints*)in.fTrackPoints->Clone();
 if (in.fITSTrackPoints)
   fITSTrackPoints = (AliHBTTrackPoints*)in.fITSTrackPoints->Clone();
 if (in.fClusterMap)
   fClusterMap = (AliHBTClusterMap*)in.fClusterMap->Clone();
}

//______________________________________________________________________________
AliHBTParticle::AliHBTParticle(const TParticle &p,Int_t idx):
   fPdgIdx(0), fIdxInEvent(idx),
   fNPids(0),fPids(0x0),fPidProb(0x0),
   fCalcMass(p.GetCalcMass()),
   fPx(p.Px()),fPy(p.Py()),fPz(p.Pz()),fE(p.Energy()), 
   fVx(p.Vx()),fVy(p.Vy()),fVz(p.Vz()),fVt(p.T()),
   fTrackPoints(0x0), fITSTrackPoints(0x0), fClusterMap(0x0)
{
 //all copied in the initialization
 SetPdgCode(p.GetPdgCode());
}
//______________________________________________________________________________

AliHBTParticle::~AliHBTParticle()
{
//dtor  
  delete [] fPids;
  delete [] fPidProb;
  delete fTrackPoints;
  delete fITSTrackPoints;
  delete fClusterMap;
}
//______________________________________________________________________________

AliHBTParticle& AliHBTParticle::operator=(const AliHBTParticle& in)
{
//assigment operator
  
  fNPids = in.fNPids;
  delete [] fPids;
  delete [] fPidProb;
  Int_t* fPids = new Int_t[fNPids];
  Float_t* fPidProb = new Float_t[fNPids];
  for (Int_t i = 0; i < fNPids;i++)
   {
     fPids[i]    = in.fPids[i];
     fPidProb[i] = in.fPidProb[i];
   }
   
  fPdgIdx = in.fPdgIdx; 
  fIdxInEvent = in.fIdxInEvent;  
  fCalcMass = in.GetCalcMass();
  fPx = in.Px();
  fPy = in.Py();
  fPz = in.Pz();
  fE = in.Energy(); 
  fVx = in.Vx();
  fVy = in.Vy();
  fVz = in.Vz();
  fVt = in.T();
  
  delete fTrackPoints;
  fTrackPoints = (in.fTrackPoints)?(AliHBTTrackPoints*)in.fTrackPoints->Clone():0x0;

  delete fITSTrackPoints;
  fITSTrackPoints = (in.fTrackPoints)?(AliHBTTrackPoints*)in.fITSTrackPoints->Clone():0x0;
  
  delete fClusterMap;
  fClusterMap =  (in.fClusterMap)?(AliHBTClusterMap*)in.fClusterMap->Clone():0x0;
  
  return *this;
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
  if (fgDebug > 9) Info("SetPIDprobability","Setting PID %d prob %f",pdg,prob);

  Float_t totprob = 0.0;//sums up probabilities
  Int_t idx = GetPidSlot(pdg);
  Int_t i;
  if (idx > -1) 
   {
     fPidProb[idx] = prob;
     for (i = 0; i < fNPids;i++) totprob+=fPidProb[i];
     if (totprob > (1.0+0.000001))
       {
         Warning("SetPIDprobability","Total probability greater than unity (%f)",totprob);
       }
     if (fgDebug > 9) 
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
        if (fgDebug>9) Info("SetPID","Copying entry %d",i);
        aPidProbNew[i] = fPidProb[i];
        aPidsNew[i] = fPids[i];
        totprob+=fPidProb[i];
      }
     else break;
   }

  if (fgDebug > 9) Info("SetPID","Setting new PID on entry %d",i);
  aPidProbNew[i] = prob;
  aPidsNew[i] = pdg;
  totprob+=prob;
  

  for (Int_t j = fNPids-1; j > i ;j--)//copy rest of old araays 
   {
     if (fgDebug > 9) Info("SetPID","Copying from old entry %d to new entry %d",j-1,j);
     aPidProbNew[j] = fPidProb[j-1];
     aPidsNew[j] = fPids[j-1];
     totprob+=fPidProb[j-1];
   }

  delete [] fPidProb;
  delete [] fPids;
  
  fPidProb = aPidProbNew;
  fPids = aPidsNew;
  
  fPdgIdx = GetPidSlot(currentpid);
  if (fPdgIdx == -1) fPdgIdx = 0;
  
  if (totprob > (1.0+0.000001))//place for numerical error
   {
     Warning("SetId","Total probability is greater than unity (%f)!!!",totprob);
     Print();
   }
}
//______________________________________________________________________________

Float_t AliHBTParticle::GetPIDprobability(Int_t pdg) const
{
//Returns probability that this particle is the type of pdg
  Int_t idx = GetPidSlot(pdg);
  if (idx < 0) return 0.0;//such pid was not specified for this particle
  return fPidProb[idx];
}
//______________________________________________________________________________

const Char_t* AliHBTParticle::GetName() const 
{
  //returns name of this particle 
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

Int_t AliHBTParticle::GetNthPid(Int_t idx) const
{
  //returns PID sitting on slot idx in fPids
  if ( (idx < 0) || (idx >= fNPids) )
   {
     Error("GetNthPid","Out Of Bounds");
     return 0;
   }
  return fPids[idx];
}
//______________________________________________________________________________

Float_t AliHBTParticle::GetNthPidProb(Int_t idx) const
{
  //returns PID sitting on slot idx in fPidProb
  if ( (idx < 0) || (idx >= fNPids) )
   {
     Error("GetNthPid","Out Of Bounds");
     return 0;
   }
  return fPidProb[idx];
}
//______________________________________________________________________________

void AliHBTParticle::Print() const
{
//prints information about particle
  printf("____________________________________________________\n");
  printf("Idx: %d  PID: %d  Name: ",fIdxInEvent,GetPdgCode());
  TParticlePDG *pdgp = TDatabasePDG::Instance()->GetParticle(GetPdgCode());
  if (pdgp)
   {
     printf("%s Mass: %f\n",pdgp->GetName(),pdgp->Mass());
   }
  else
   {
     printf("Not known\n");
   }
  
  printf("Px: %+f Py: %+f Pz: %+f E: %+f Calculated Mass: %f\nVx: %+f Vy: %+f Vz: %+f T: %+f\n",
          Px(),Py(),Pz(),Energy(),GetCalcMass(),Vx(),Vy(),Vz(),T());

  for (Int_t i = 0; i < fNPids; i++)
   {
     printf("# %d  PID: %d  Probability %f name ",i,fPids[i],fPidProb[i]);
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

//______________________________________________________________________________

//void AliHBTParticle::Streamer(TBuffer &b)
//{
//     // Stream all objects in the array to or from the I/O buffer.
//   UInt_t R__s, R__c;
//   Int_t i;
//   if (b.IsReading()) 
//    {
//      delete [] fPids;
//      delete [] fPidProb;
//      
//      Version_t v = b.ReadVersion(&R__s, &R__c);
//      if (v == 1)
//       {
//         AliHBTParticle::Class()->ReadBuffer(b, this);
//      }
//      else
//       {
//        TObject::Streamer(b);
//       b >> fPdgIdx;
//       b >> fIdxInEvent;
//       
//       b >> fNPids;
//       Int_t* fPids = new Int_t[fNPids];
//        Float_t* fPidProb = new Float_t[fNPids];
//        for (i = 0;i<fNPids;i++) 
//         {
//           b >> fPids[i];
//         }
//        for (i = 0;i<fNPids;i++) 
//        {
//          b >> fPidProb[i];
//         }
//        b >> fCalcMass;
//
//        b >> fPx;
//       b >> fPy;
//       b >> fPz;
//        b >> fE;
//
//       b >> fVx;
//        b >> fVy;
//        b >> fVz;
//       b >> fVt;
//       Info("Streamer","Read data");
//        Print();
//       }
//
//      b.CheckByteCount(R__s, R__c,AliHBTParticle::IsA());
//    } 
//  else 
//   {
//     R__c = b.WriteVersion(AliHBTParticle::IsA(), kTRUE);
//     TObject::Streamer(b);
//     Info("Streamer","Read data");
//     Print();
//
//     b << fPdgIdx;
//     b << fIdxInEvent;
//     b << fNPids;
//     for (i = 0;i<fNPids;i++) 
//      {
//        b << fPids[i];
//      }
//      {
//      {
//     for (i = 0;i<fNPids;i++) 
//     {
//        b << fPidProb[i];
//      }
//     b << fCalcMass;
//
//     b << fPx;
//     b << fPy;
//     b << fPz;
//     b << fE;
//
//     b << fVx;
//     b << fVy;
//     b << fVz;
//     b << fVt;
//
//    b.SetByteCount(R__c, kTRUE);
//   }
//}
