// $Id$

#ifndef ALIJFCLUSTERH
#define ALIJFCLUSTERH

#include <Riostream.h>
#include <vector>
#include <TMath.h>

#include "AliJFPreCluster.h"

class AliJFCluster 
{
 public:
  AliJFCluster(Int_t n=100);
  AliJFCluster(const AliJFCluster &copy);
  AliJFCluster(AliJFPreCluster &copy);
  AliJFCluster(AliJFPreCluster *precluster);
  virtual ~AliJFCluster();

  AliJFCluster& operator=(const AliJFCluster &rhs);
  AliJFCluster& operator=(AliJFPreCluster &rhs);
  AliJFCluster& operator+=(AliJFCluster &rhs);
  AliJFCluster& operator+=(AliJFPreCluster &rhs);
  friend ostream& operator<< (ostream &o, const AliJFCluster &c);

  void CombineCluster(AliJFCluster &rhs);

  inline Int_t  const GetStatus() const {return fStatus;}
  inline Int_t  const GetNMerge() const {return fNMerge;}
  inline Bool_t const IsValid()   const {return(fStatus==1   ? kTRUE:kFALSE);}
  inline Bool_t const IsMerged()  const {return(fStatus==10  ? kTRUE:kFALSE);}
  inline Bool_t const IsJet()     const {return(fStatus==100 ? kTRUE:kFALSE);}
  inline Bool_t const IsInValid() const {return(fStatus<=0   ? kTRUE:kFALSE);}

  inline Float_t const GetPx()   const {return fPx   ;}
  inline Float_t const GetPy()   const {return fPy   ;}
  inline Float_t const GetPz()   const {return fPz   ;}
  inline Float_t const GetE()    const {return fE    ;}
  inline Float_t const GetY()    const {return fY    ;}
  inline Float_t const GetPhi()  const {return fPhi  ;}
  inline Float_t const GetPt2()  const {return fPt2  ;}
  inline Float_t const GetPt2D() const {return fPt2dD;}

  Int_t const GetNCombinedCluster() const {return fList.size();};
  vector<AliJFPreCluster*> const * GetClusterList() const {return &fList;}

  inline void MarkIsValid()     {fStatus=1  ;}
  inline void MarkIsMerged()    {fStatus=10 ;}
  inline void MarkIsJet()       {fStatus=100;}
  inline void MarkIsInValid()   {fStatus=-1 ;}
  inline void SetStatus(Int_t s){fStatus=s  ;}

  void Print();

  static void SetD(Float_t D_){D2=D_*D_;}

 protected:
  void SetValues();

  void SetValues(Float_t px, Float_t py, Float_t pz, Float_t E=-1);
  void AddValues(Float_t px, Float_t py, Float_t pz, Float_t E);

  Int_t fStatus;
  Int_t fNMerge;
  Float_t fPx;
  Float_t fPy;
  Float_t fPz;
  Float_t fE;
  Float_t fY;
  Float_t fPhi;
  Float_t fPt2;
  Float_t fPt2dD;

  vector<AliJFPreCluster*> fList; //->

  static Float_t D2; //static D*D for K_t comparison 

  ClassDef(AliJFCluster,1) //AliJFCluster class
};

#endif /*ALIJFCLUSTERH*/
