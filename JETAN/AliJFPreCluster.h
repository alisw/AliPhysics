// $Id$

#ifndef ALIJFPRECLUSTERH
#define ALIJFPRECLUSTERH


#include <TClonesArray.h>

class TParticle;

class AliJFPreCluster 
{
 public:
  AliJFPreCluster();
  AliJFPreCluster(const AliJFPreCluster &copy);
  AliJFPreCluster(const TParticle *p);
  AliJFPreCluster(Float_t px, Float_t py, Float_t pz, Float_t E, const TParticle *p);
  AliJFPreCluster(Float_t px, Float_t py, Float_t pz, Float_t E, TClonesArray *parts);
  virtual ~AliJFPreCluster();

  AliJFPreCluster& operator=(const AliJFPreCluster &rhs);
  friend ostream& operator<< (ostream &o, const AliJFPreCluster &c);

  inline Float_t const GetPx()  const {return fPx ;}
  inline Float_t const GetPy()  const {return fPy ;}
  inline Float_t const GetPz()  const {return fPz ;}
  inline Float_t const GetE()   const {return fE  ;}

  inline const TClonesArray* GetParticles() const {return &fParticles;}
  TClonesArray* GetParticles() {return &fParticles;}

 protected:
  //void SetValues(Float_t px, Float_t py, Float_t pz, Float_t E=-1);

  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Float_t fE;

  TClonesArray fParticles;

  ClassDef(AliJFPreCluster,1) //AliJFPreCluster class
};

#endif /*ALIJFPRECLUSTERH*/
