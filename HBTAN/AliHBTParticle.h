#ifndef ALIHBTPARTICLE_H
#define ALIHBTPARTICLE_H
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

#include <TObject.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TDatabasePDG.h>


class TParticle;
class AliHBTTrackPoints;
class AliHBTClusterMap;

class AliHBTParticle : public TObject
{
public:
                                // ****** constructors and destructor
  AliHBTParticle();
  AliHBTParticle(const AliHBTParticle& in); 
 
  AliHBTParticle(Int_t pdg, Int_t idx, Double_t px, Double_t py, Double_t pz, Double_t etot,
                 Double_t vx, Double_t vy, Double_t vz, Double_t time);

  AliHBTParticle(Int_t pdg, Float_t prob, Int_t idx, Double_t px, Double_t py, Double_t pz, Double_t etot,
                 Double_t vx, Double_t vy, Double_t vz, Double_t time);

  AliHBTParticle(const TParticle& p,Int_t idx);

  virtual ~AliHBTParticle();
  
  AliHBTParticle& operator=(const AliHBTParticle& in); 
  
  void           SetPIDprobability(Int_t pdg, Float_t prob = 1.0);
  Float_t        GetPIDprobability(Int_t pdg) const;
  
  Int_t          GetPdgCode      () const { return (fPids)?fPids[fPdgIdx]:0;}
  Int_t          GetPid          () const { return GetPdgCode();}
  Float_t        GetPidProb      () const { return (fPidProb)?fPidProb[fPdgIdx]:0;}
  
  Int_t          GetUID          () const { return fIdxInEvent;}
  Int_t          GetNumberOfPids () const { return fNPids;}
  Int_t          GetNthPid         (Int_t idx) const;
  Float_t        GetNthPidProb     (Int_t idx) const;
      
  void           SetPdgCode(Int_t pdg, Float_t prob = 1.0);
  Double_t       GetCalcMass     () const { return fCalcMass; }
  Double_t       GetMass         ()       { return (GetPDG())?GetPDG()->Mass():-1.;}

  TParticlePDG*  GetPDG          (){return TDatabasePDG::Instance()->GetParticle(GetPdgCode());}

  Int_t          Beauty          ()  { return GetPDG()->Beauty(); }
  Int_t          Charm           ()  { return GetPDG()->Charm(); }
  Int_t          Strangeness     ()  { return GetPDG()->Strangeness();}
  void ProductionVertex(TLorentzVector &v) const { v.SetXYZT(fVx,fVy,fVz,fVt);}


  Double_t         Vx    () const { return fVx;}
  Double_t         Vy    () const { return fVy;}
  Double_t         Vz    () const { return fVz;}
  Double_t         T     () const { return fVt;}

  Double_t         Px    () const { return fPx; } //X coordinate of the momentum
  Double_t         Py    () const { return fPy; } //Y coordinate of the momentum
  Double_t         Pz    () const { return fPz; } //Z coordinate of the momentum
  Double_t         P     () const                 //momentum
    { return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz); }
  
  void Momentum(TLorentzVector &v) const { v.SetPxPyPzE(fPx,fPy,fPz,fE);}
    
  Double_t         Pt    () const  //transverse momentum
    { return TMath::Sqrt(fPx*fPx+fPy*fPy); }
  Double_t         Energy() const { return fE; }
  
                                   //Pseudo Rapidity
  Double_t         Eta   () const { if (P() != fPz) return 0.5*TMath::Log((P()+fPz)/(P()-fPz)); 
                                    else return 1.e30;}

                                   //Rapidity
  Double_t         Y     () const { if (fE  != fPz) return 0.5*TMath::Log((fE+fPz)/(fE-fPz));
                                    else return 1.e30;}

  Double_t         Phi   () const { return TMath::Pi()+TMath::ATan2(-fPy,-fPx); }

  Double_t         Theta () const { return (fPz==0)?TMath::PiOver2():TMath::ACos(fPz/P()); }

  // setters

  void           SetMomentum(Double_t px, Double_t py, Double_t pz, Double_t e)
                             {fPx=px; fPy=py; fPz=pz; fE=e;}
  void           SetMomentum(const TLorentzVector& p)
                             {SetMomentum(p.Px(),p.Py(),p.Pz(),p.Energy());}

  void           SetProductionVertex(Double_t vx, Double_t vy, Double_t vz, Double_t t)
                             {fVx=vx; fVy=vy; fVz=vz; fVt=t;}
  void           SetProductionVertex(const TLorentzVector& v)
                             {SetProductionVertex(v.X(),v.Y(),v.Z(),v.T());}
  void           SetCalcMass(Double_t mass) {fCalcMass = mass;}
  
  void           SetUID(Int_t id){fIdxInEvent = id;}
  
  const Char_t*  GetName() const; 
  void           Print() const;
  
  void           SetTrackPoints(AliHBTTrackPoints* tpts){fTrackPoints = tpts;}
  AliHBTTrackPoints* GetTrackPoints() const {return fTrackPoints;}

  void           SetITSTrackPoints(AliHBTTrackPoints* tpts){fITSTrackPoints = tpts;}
  AliHBTTrackPoints* GetITSTrackPoints() const {return fITSTrackPoints;}
  
  void           SetClusterMap(AliHBTClusterMap* cm){fClusterMap = cm;}
  AliHBTClusterMap* GetClusterMap() const {return fClusterMap;}
  
  static void    SetDebug(Int_t dbg=1){fgDebug=dbg;}
  static Int_t   GetDebug(){return fgDebug;}
  
protected:
  Int_t          GetPidSlot(Int_t pdg) const;//returns position of the given PID in fPids (and fPidProb) array.

private:
  Char_t         fPdgIdx;               // index of PDG code of the particle in fPids
  Int_t          fIdxInEvent;           // index of a particle: the same particle can appear in the event
                                        //  many times with different pid's. Idx allows to check that they are the same particles
  Int_t          fNPids;                // number of non-zero proboble Pids
  Int_t         *fPids;                 // [fNPids] Array with PIDs
  Float_t       *fPidProb;              // [fNPids] PIDs probabilities
  Double_t       fCalcMass;             // Calculated mass

  Double_t       fPx;                   // x component of momentum
  Double_t       fPy;                   // y component of momentum
  Double_t       fPz;                   // z component of momentum
  Double_t       fE;                    // Energy

  Double_t       fVx;                   // x of production vertex
  Double_t       fVy;                   // y of production vertex
  Double_t       fVz;                   // z of production vertex
  Double_t       fVt;                   // t of production vertex

  AliHBTTrackPoints* fTrackPoints;      // track positions along trajectory - used by anti-merging cut 
  AliHBTTrackPoints* fITSTrackPoints;   // track position at first pixels
  
  AliHBTClusterMap*  fClusterMap;       // bit map of cluters occupation; 1 if has cluter on given layer/padrow/...
    
  static Int_t   fgDebug; //debug printout level
  ClassDef(AliHBTParticle,4)  // TParticle vertex particle information
};

#endif
