#ifndef ROOT_GParticle
#define ROOT_GParticle
 
#include "TObject.h"
#include "TMath.h"
 
#define Keep_Bit     1
#define Children_Bit 2
#define Done_Bit     4

class GParticle : public TObject {
 
private:
  
  Int_t    fKS;            // status of particle       ( LUJETS K[1] )
  Int_t    fKF;            // Geant3 particle type
  Int_t    fParent;        // parent id 
  Int_t    fFirstChild;    // id of first child 
  Int_t    fLastChild;     // id of last  child 
  
  Float_t  fPx;            // X momentum [GeV/c]
  Float_t  fPy;            // Y momentum [GeV/c]
  Float_t  fPz;            // Z momentum [GeV/c]
  Float_t  fEnergy;        // Energy    [GeV]
  Float_t  fMass;          // Mass      [Gev/c^2]
  
  Float_t  fVx;            // X vertex  [cm]
  Float_t  fVy;            // Y vertex  [cm]
  Float_t  fVz;            // Z vertex  [cm]
  
  Float_t  fPolx;          // X component of polarisation 
  Float_t  fPoly;          // Y component of polarisation 
  Float_t  fPolz;          // Z component of polarisation 
  
  Float_t  fTime;          // time of procuction
  Float_t  fLifeTime;      // proper lifetime
  Float_t  fProcessTime;   // Cpu Time to process this track
  
  char     fOrigin[11];    // generation mechanism
  Float_t  fWgt;
  Float_t  fPT;
  Float_t  fTheta;
  Float_t  fPhi;
  Float_t  fEta;
  
public:
  GParticle() { }
  
  GParticle(const GParticle &p);
  
  GParticle(const Int_t kS, const Int_t kF, 
	    const Int_t parent, const Int_t firstchild, const Int_t lastchild, 
	    const Float_t px, const Float_t py, const Float_t pz,
	    const Float_t energy, const Float_t mass,
	    const Float_t vx, const Float_t vy, const Float_t vz, 
	    const Float_t polx, const Float_t poly, const Float_t polz,
	    const Float_t time, const Float_t lifetime, 
	    const char *Origin="Unknown",Float_t wgt=1.) :
    
    fKS(kS),
    fKF(kF),
    fParent(parent),
    fFirstChild(firstchild),
    fLastChild(lastchild),
    fPx(px),
    fPy(py),
    fPz(pz),
    fEnergy(energy),
    fMass(mass),
    fVx(vx),
    fVy(vy),
    fVz(vz),
    fPolx(polx),
    fPoly(poly),
    fPolz(polz),
    fTime(time),
    fLifeTime(lifetime),
      fProcessTime(-0.1),
      fWgt(wgt)
    {Int_t p=0; 
    while((fOrigin[p]=Origin[p]) && p++<10);fOrigin[10]=0;
    fPT    =TMath::Sqrt(fPx*fPx+fPy*fPy);
    fTheta =Float_t(TMath::ATan2(Double_t(fPT),Double_t(fPz)));
    fPhi   =Float_t(TMath::ATan2(Double_t(fPy),Double_t(fPx)))+TMath::Pi();
    fEta   =-TMath::Log(TMath::Tan(fTheta/2));
    }
  
  virtual             ~GParticle() { }
  
  inline Int_t       GetKS() const {return fKS;}
  inline Int_t       GetKF() const {return fKF;}
  inline Int_t       GetParent() const {return fParent;}
  inline Int_t       GetFirstChild() const {return fFirstChild;}
  inline Int_t       GetLastChild() const {return fLastChild;}
  
  inline Float_t     GetPx() const {return fPx;}
  inline Float_t     GetPy() const {return fPy;}
  inline Float_t     GetPz() const {return fPz;}
  
  inline Float_t     GetPolx() const {return fPolx;}
  inline Float_t     GetPoly() const {return fPoly;}
  inline Float_t     GetPolz() const {return fPolz;}
  
  inline Float_t     GetEnergy() const {return fEnergy;}
  inline Float_t     GetMass() const {return fMass;}
  inline Float_t     GetMomentum() const {return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);}
  inline void        GetOrigin(char *Origin) {strcpy(Origin,fOrigin);}
  inline Float_t     GetPT()  const {return  fPT;}
  inline Float_t     GetTheta()  const {return  fTheta;}
  inline Float_t     GetEta() const {return fEta;}
  inline Float_t     GetPhi() const {return fPhi;}
  
  inline Float_t     GetVx() const {return fVx;}
  inline Float_t     GetVy() const {return fVy;}
  inline Float_t     GetVz() const {return fVz;}
  inline Float_t     GetTime() const {return fTime;}
  inline Float_t     GetLifeTime() const {return fLifeTime;}
  inline Float_t     GetWgt() const {return fWgt;}
  inline void        SetWgt(Float_t wgt) {fWgt=wgt;}  
  inline Float_t     GetProcessTime() const {return fProcessTime;}
  virtual const Text_t     *GetName()  const;
  virtual const Text_t     *GetTitle() const;
  
  inline void        SetKS(Int_t KS) {fKS=KS;}
  inline void        SetFirstChild(const Int_t firstchild) {fFirstChild=firstchild;}
  inline void        SetLastChild(const Int_t lastchild) {fLastChild=lastchild;}
  inline void        SetParent(const Int_t parent) {fParent=parent;}
  inline void        SetProcessTime(Float_t cputime) {fProcessTime=cputime;}
  
  GParticle           &operator=(const GParticle &p);
  
  ClassDef(GParticle,1)  // MonteCarlo Particle Class
};
 
inline GParticle &GParticle::operator = (const GParticle &p)
{
  fKS                   = p.fKS;
  fKF                   = p.fKF;
  fParent               = p.fParent;
  fFirstChild           = p.fFirstChild;
  fLastChild            = p.fLastChild;
  fPx                   = p.fPx;
  fPy                   = p.fPy;
  fPz                   = p.fPz;
  fEnergy               = p.fEnergy;
  fMass                 = p.fMass;
  fVx                   = p.fVx;
  fVy                   = p.fVy;
  fVz                   = p.fVz;
  fPolx                 = p.fPolx;
  fPoly                 = p.fPoly;
  fPolz                 = p.fPolz;
  fTime                 = p.fTime;
  fLifeTime             = p.fLifeTime;
  fProcessTime          = p.fProcessTime;
  fWgt                  = p.fWgt;
  fPT                   = p.fPT;
  fTheta                = p.fTheta;
  fPhi                  = p.fPhi;
  fEta                  = p.fEta;
  strcpy(fOrigin,p.fOrigin);
  return *this;
}


inline GParticle::GParticle(const GParticle &p) { *this = p; }

#endif


