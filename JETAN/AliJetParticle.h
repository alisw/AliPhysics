#ifndef ALIJETPARTICLE_H
#define ALIJETPARTICLE_H
/* $Id$ */

//___________________________________________________________
/////////////////////////////////////////////////////////////
//                                                         //
// class AliJetParticle                                    //
//                                                         //
// loizides@ikf.uni-frankfurt.de                           //
//                                                         //
/////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>
class TParticle;

class AliJetParticle : public TObject
{
  public:
  AliJetParticle();
  AliJetParticle(const AliJetParticle& in); 
  AliJetParticle(const TParticle* p, Int_t idx, Int_t l=0);
  AliJetParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l=0);
  AliJetParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l,
                 Float_t pt, Float_t eta, Float_t phi);
  virtual ~AliJetParticle(){};

  void SetMomentum(Float_t px, Float_t py, Float_t pz, Float_t e)
                    {fPx=px; fPy=py; fPz=pz; fE=e; Calculate();}
  void SetMomentum(Float_t px, Float_t py, Float_t pz, Float_t e,
		   Float_t pt, Float_t phi, Float_t eta)
                    {fPx=px;fPy=py;fPz=pz;fE=e;fPt=pt;fEta=eta;fPhi=phi;}
  
  void SetUID(Int_t id) {fIdxInEvent = id;}
  void SetLabel(Int_t l){fLabel = l;}
  
  Float_t P()      const {return TMath::Sqrt(fPx*fPx+fPy*fPy+fPz*fPz);}
  Float_t Y()      const {if (fE  != fPz) return 0.5*TMath::Log((fE+fPz)/(fE-fPz));
                                    else return 1.e30;}
  Float_t Theta () const {return (fPz==0)?TMath::PiOver2():TMath::ACos(fPz/P());}

  Int_t GetUID()   const {return fIdxInEvent;}
  Int_t GetLabel() const {return fLabel;}

  Float_t Px()     const {return fPx;}
  Float_t Py()     const {return fPy;}
  Float_t Pz()     const {return fPz;}
  Float_t Energy() const {return fE;}

  Float_t Pt()     const {return fPt;}
  Float_t Eta()    const {return fEta;}
  Float_t Phi()    const {return fPhi;} 

  void Clear(Option_t *t="");
  void Print(Option_t *t="") const;

  protected:
  void Calculate();     //calculate values

  Float_t fPx;          // x component of momentum at vertex
  Float_t fPy;          // y component of momentum at vertex
  Float_t fPz;          // z component of momentum at vertex
  Float_t fE;           // total energy
  Int_t   fIdxInEvent;  // index of particle as appeared in complete event
  Int_t   fLabel;
  Float_t fPt;          // normally calculated 
  Float_t fEta;         // normally calculated 
  Float_t fPhi;         // normally calculated 

  ClassDef(AliJetParticle,1)  // Basic Jet Particle class
};
#endif
