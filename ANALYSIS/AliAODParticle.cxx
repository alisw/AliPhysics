#include "AliAODParticle.h"
//___________________________________________________________
/////////////////////////////////////////////////////////////
//
// class AliAODParticle
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
#include "AliTrackPoints.h"
#include "AliClusterMap.h"

ClassImp(AliAODParticle)

//______________________________________________________________________________
AliAODParticle::AliAODParticle():  
 fPdgIdx(0), fIdxInEvent(0),fNPids(0),fPids(0x0),fPidProb(0x0),
 fCalcMass(0),fPx(0), fPy(0),fPz(0),fE(0), fVx(0), fVy(0),fVz(0),fVt(0),
 fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
{//empty particle
}
//______________________________________________________________________________

AliAODParticle::AliAODParticle(Int_t pdg, Int_t idx,
               Double_t px, Double_t py, Double_t pz, Double_t etot,
               Double_t vx, Double_t vy, Double_t vz, Double_t time):  
  fPdgIdx(0), fIdxInEvent(idx),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time),
  fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
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

AliAODParticle::AliAODParticle(Int_t pdg, Float_t prob, Int_t idx, 
                               Double_t px, Double_t py, Double_t pz, Double_t etot,
                               Double_t vx, Double_t vy, Double_t vz, Double_t time):
  fPdgIdx(0), fIdxInEvent(idx),fNPids(0),fPids(0x0),fPidProb(0x0),
  fCalcMass(0), 
  fPx(px), fPy(py),fPz(pz),fE(etot), 
  fVx(vx), fVy(vy),fVz(vz),fVt(time),
  fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
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

AliAODParticle::AliAODParticle(const AliAODParticle& in):
   AliVAODParticle(in),
   fPdgIdx(in.fPdgIdx), fIdxInEvent(in.fIdxInEvent),
   fNPids(in.fNPids),fPids(new Int_t[fNPids]),fPidProb(new Float_t[fNPids]),
   fCalcMass(in.GetCalcMass()),
   fPx(in.Px()),fPy(in.Py()),fPz(in.Pz()),fE(in.E()), 
   fVx(in.Vx()),fVy(in.Vy()),fVz(in.Vz()),fVt(in.T()),
   fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
{
 //Copy constructor
// Info("AliAODParticle(const AliAODParticle& in)","");
 for(Int_t i = 0; i<fNPids; i++)
  {
    fPids[i] =  in.fPids[i];
    fPidProb[i] = in.fPidProb[i];
  }
 
 if (in.fTPCTrackPoints)
   fTPCTrackPoints = (AliTrackPoints*)in.fTPCTrackPoints->Clone();
 if (in.fITSTrackPoints)
   fITSTrackPoints = (AliTrackPoints*)in.fITSTrackPoints->Clone();
 if (in.fClusterMap)
   fClusterMap = (AliClusterMap*)in.fClusterMap->Clone();
}
//______________________________________________________________________________

AliAODParticle::AliAODParticle(const AliVAODParticle& in):
   AliVAODParticle(in),
   fPdgIdx(0), fIdxInEvent(in.GetUID()),
   fNPids(0),fPids(0x0),fPidProb(0x0),
   fCalcMass(-1.0),
   fPx(in.Px()),fPy(in.Py()),fPz(in.Pz()),fE(in.E()), 
   fVx(in.Vx()),fVy(in.Vy()),fVz(in.Vz()),fVt(in.T()),
   fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
{
 //Copy constructor
// Info("AliAODParticle(const AliVAODParticle& in)","");
 for(Int_t i = 0; i<in.GetNumberOfPids(); i++)
  {
    SetPIDprobability(in.GetNthPid(i),in.GetNthPidProb(i));
  }
 SetPdgCode(in.GetPdgCode(),in.GetPidProb());
 
 AliTrackPoints* tpts = in.GetTPCTrackPoints();
 if (tpts)  SetTPCTrackPoints((AliTrackPoints*)tpts->Clone());
 
 tpts = in.GetITSTrackPoints();  
 if (tpts) SetITSTrackPoints((AliTrackPoints*)tpts->Clone());
   
 AliClusterMap* clmap = in.GetClusterMap();
 if (clmap) SetClusterMap((AliClusterMap*)clmap->Clone());
}
//______________________________________________________________________________

AliAODParticle::AliAODParticle(const TParticle &p,Int_t idx):
   fPdgIdx(0), fIdxInEvent(idx),
   fNPids(0),fPids(0x0),fPidProb(0x0),
   fCalcMass(p.GetCalcMass()),
   fPx(p.Px()),fPy(p.Py()),fPz(p.Pz()),fE(p.Energy()), 
   fVx(p.Vx()),fVy(p.Vy()),fVz(p.Vz()),fVt(p.T()),
   fTPCTrackPoints(0x0),fITSTrackPoints(0x0),fClusterMap(0x0)
{
 //all copied in the initialization
 SetPdgCode(p.GetPdgCode());
}
//______________________________________________________________________________

AliAODParticle::~AliAODParticle()
{
//dtor  
  delete [] fPids;
  delete [] fPidProb;
  delete fTPCTrackPoints;
  delete fITSTrackPoints;
  delete fClusterMap;
}
//______________________________________________________________________________

void AliAODParticle::Clear(Option_t*)
{
//Must be implemented in order to store this object in Clones Array
  delete [] fPids;
  delete [] fPidProb;
  delete fTPCTrackPoints;
  delete fITSTrackPoints;
  delete fClusterMap;
  
  fPids = 0x0;
  fPidProb = 0x0;
  fTPCTrackPoints = 0x0;
  fITSTrackPoints = 0x0;
  fClusterMap = 0x0;
}
//______________________________________________________________________________

AliAODParticle& AliAODParticle::operator=(const AliVAODParticle& in)
{
//operator=
//  Info("operator=(const AliVAODParticle& in)","AliAODParticle");
  
  if (&in == this) return *this;

  delete [] fPids;
  delete [] fPidProb;
  fPids = 0x0;
  fPidProb = 0x0;
  fNPids = 0;

  Int_t npids = in.GetNumberOfPids();
  for (Int_t i = 0; i < npids; i++)
   {
     SetPIDprobability(in.GetNthPid(i),in.GetNthPidProb(i));
   }
   
  SetPdgCode(in.GetPdgCode(),in.GetPidProb());

  SetUID(in.GetUID());

  fCalcMass = in.Mass();
  
  fPx = in.Px();
  fPy = in.Py();
  fPz = in.Pz();
  fE = in.E(); 
  fVx = in.Vx();
  fVy = in.Vy();
  fVz = in.Vz();
  fVt = in.T();
  
  delete fTPCTrackPoints;
  AliTrackPoints* tpts = in.GetTPCTrackPoints();
  fTPCTrackPoints = (tpts)?(AliTrackPoints*)tpts->Clone():0x0;

  delete fITSTrackPoints;
  tpts = in.GetITSTrackPoints();  
  fITSTrackPoints = (tpts)?(AliTrackPoints*)tpts->Clone():0x0;
  
  delete fClusterMap;
  AliClusterMap* incmap = in.GetClusterMap();
  fClusterMap =  (incmap)?(AliClusterMap*)incmap->Clone():0x0;
  
  return *this;
}
//______________________________________________________________________________

AliAODParticle& AliAODParticle::operator=(const AliAODParticle& in)
{
//assigment operator
//  Info("operator=(const AliAODParticle& in)","AliAODParticle");
  if (&in == this) return *this;
  fNPids = in.fNPids;
  delete [] fPids;
  delete [] fPidProb;
  fPids = new Int_t[fNPids];
  fPidProb = new Float_t[fNPids];
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
  fE = in.E(); 
  fVx = in.Vx();
  fVy = in.Vy();
  fVz = in.Vz();
  fVt = in.T();
  
  delete fTPCTrackPoints;
  fTPCTrackPoints = (in.fTPCTrackPoints)?(AliTrackPoints*)in.fTPCTrackPoints->Clone():0x0;

  delete fITSTrackPoints;
  fITSTrackPoints = (in.fITSTrackPoints)?(AliTrackPoints*)in.fITSTrackPoints->Clone():0x0;
  
  delete fClusterMap;
  fClusterMap =  (in.fClusterMap)?(AliClusterMap*)in.fClusterMap->Clone():0x0;
  
  return *this;
}
//______________________________________________________________________________

void AliAODParticle::SetPdgCode(Int_t pdg,Float_t prob)
{
//Set PDG Code
  SetPIDprobability(pdg,prob);
  fPdgIdx = GetPidSlot(pdg);
}

//______________________________________________________________________________
void AliAODParticle::SetPIDprobability(Int_t pdg, Float_t prob)
{
//Sets another pdg code and corresponding probabilty
//Ids are set in decreasing order
//Check if total probability is not overcoming unity is performed
//in case, warning is printed
  if (GetDebug() > 9) Info("SetPIDprobability","Setting PID %d prob %f",pdg,prob);

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
     if (GetDebug() > 9) 
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
        if (GetDebug()>9) Info("SetPID","Copying entry %d",i);
        aPidProbNew[i] = fPidProb[i];
        aPidsNew[i] = fPids[i];
        totprob+=fPidProb[i];
      }
     else break;
   }

  if (GetDebug() > 9) Info("SetPID","Setting new PID on entry %d",i);
  aPidProbNew[i] = prob;
  aPidsNew[i] = pdg;
  totprob+=prob;
  

  for (Int_t j = fNPids-1; j > i ;j--)//copy rest of old arays 
   {
     if (GetDebug() > 9) Info("SetPID","Copying from old entry %d to new entry %d",j-1,j);
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
  
  if (totprob > (1.0+0.000001))//space for numerical error
   {
     Warning("SetId","Total probability is greater than unity (%f)!!!",totprob);
     Print();
   }
}
//______________________________________________________________________________

Float_t AliAODParticle::GetPIDprobability(Int_t pdg) const
{
//Returns probability that this particle is the type of pdg
  Int_t idx = GetPidSlot(pdg);
  if (idx < 0) return 0.0;//such pid was not specified for this particle
  return fPidProb[idx];
}
//______________________________________________________________________________

const Char_t* AliAODParticle::GetName() const 
{
  //returns name of this particle 
   static char def[4] = "XXX";
   const TParticlePDG *ap = TDatabasePDG::Instance()->GetParticle(GetPdgCode());
   if (ap) return ap->GetName();
   else    return def;
}
//______________________________________________________________________________

Int_t AliAODParticle::GetPidSlot(Int_t pdg) const
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

Int_t AliAODParticle::GetNthPid(Int_t idx) const
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

Float_t AliAODParticle::GetNthPidProb(Int_t idx) const
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

void AliAODParticle::Print() const
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
          Px(),Py(),Pz(),E(),GetCalcMass(),Vx(),Vy(),Vz(),T());

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

//void AliAODParticle::Streamer(TBuffer &b)
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
//         AliAODParticle::Class()->ReadBuffer(b, this);
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
//      b.CheckByteCount(R__s, R__c,AliAODParticle::IsA());
//    } 
//  else 
//   {
//     R__c = b.WriteVersion(AliAODParticle::IsA(), kTRUE);
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
