#ifndef ALITPCCLUSTERMI_H
#define ALITPCCLUSTERMI_H

//-------------------------------------------------------
//                    TPC Cluster Class
//
//   Origin: Marian Ivanov
//-------------------------------------------------------

#include "AliCluster.h"
#include "TMath.h"

//_____________________________________________________________________________
class AliTPCclusterMI : public AliCluster {
public:
  AliTPCclusterMI():AliCluster(){fQ=0; fUsed=0;}
  AliTPCclusterMI(Int_t *lab, Float_t *hit) : AliCluster(lab,hit) {fQ = (Short_t)hit[4];}
  virtual Bool_t IsSortable() const; 
  virtual Int_t Compare(const TObject* obj) const;
  inline  void Use(Int_t inc=10);
  void SetQ(Float_t q) {fQ=(Short_t)q;}
  void SetType(Char_t type) {fType=type;}
  void SetMax(Short_t max) {fMax=max;}
  Int_t IsUsed(Int_t th=10) const {return (fUsed>=th) ? 1 : 0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}
  Float_t GetMax() const {return fMax;} 
  Char_t  GetType()const {return fType;}
 
private:
  Short_t   fQ ;       //Q of cluster (in ADC counts)  
  Char_t    fType;     //type of the cluster 0 means golden 
  Short_t   fMax;      //maximal amplitude in cluster
  Char_t    fUsed;     //counter of usage  
 ClassDef(AliTPCclusterMI,1)  // Time Projection Chamber clusters
};

void AliTPCclusterMI::Use(Int_t inc) 
{ 
  if (inc>0)  fUsed+=inc; 
  else 
    fUsed=0;
}

class AliTPCclusterLMI  {
 private:
  Short_t  fCZ;       // current cluster position Z in cm - 100 mum precision
  Short_t  fCY;       // current cluster position Y in cm - 100 mum precision
  UChar_t  fSigmaZ;   // shape  Z - normalised shape - normaliziation 1 - precision 2 percent
  UChar_t  fSigmaY;   // shape  Y - normalised shape - normaliziation 1 - precision 2 percent
  UShort_t fQ;        // total charge in cluster 
  UShort_t fMax;      // charge at maximum  
  Char_t   fCType;    // type of the cluster
  Int_t    fLabel[3];    // track indexes 
 public:
  AliTPCclusterLMI(){fCZ=fCY=fSigmaZ=fSigmaY=fQ=fMax=fCType=0;}
  inline Float_t  GetZ()    {return (fCZ*0.01);}
  inline Float_t  GetY()    {return (fCY*0.01);}
  inline Float_t  GetSigmaZ() {return (fSigmaZ*0.02);}
  inline Float_t  GetSigmaY() {return (fSigmaY*0.02);}  
  inline Int_t  GetType()   {return fCType;}
  inline Int_t  GetMax()   {return fMax;}
  inline Float_t  GetQ()   {return fQ;}
  inline Int_t    GelLabel(Int_t i){return fLabel[i];}
  //
  void     SetY(Float_t y){ fCY = Short_t(TMath::Nint(y*100.));} 
  void     SetZ(Float_t z){ fCZ = Short_t(TMath::Nint(z*100.));} 
  void     SetSigmaZ(Float_t sigmaz) {fSigmaZ = UChar_t(TMath::Nint(sigmaz*50.));}
  void     SetSigmaY(Float_t sigmay) {fSigmaY = UChar_t(TMath::Nint(sigmay*50.));}
  void     SetQ(Float_t q) {fQ = UShort_t(q);}
  void     SetMax(Float_t max) {fMax = UShort_t(max);}
  void     SetType(Char_t type) {fCType = type;}
  void     SetLabels(Int_t labels[3]){fLabel[0] = labels[0];fLabel[1] = labels[1];fLabel[2] = labels[2];}
  //
 public:
  ClassDef(AliTPCclusterLMI,1)  
};



#endif


