#ifndef ALITPCCLUSTERMI_H
#define ALITPCCLUSTERMI_H

//-------------------------------------------------------
//                    TPC Cluster Class
//   Parallel tracking
//   Origin: Marian Ivanov
//-------------------------------------------------------

/* $Id$ */


#include "AliCluster.h"
#include "TMath.h"
class AliTPCclusterInfo;
//_____________________________________________________________________________
class AliTPCclusterMI : public AliCluster {
public:
  AliTPCclusterMI(Bool_t withInfo=kTRUE);
  AliTPCclusterMI(const AliTPCclusterMI & cluster);
  AliTPCclusterMI &operator = (const AliTPCclusterMI & cluster); //assignment operator
  AliTPCclusterMI(Int_t *lab, Float_t *hit);
  virtual ~AliTPCclusterMI();
  virtual Bool_t IsSortable() const; 
  virtual Int_t Compare(const TObject* obj) const;
  inline  void Use(Int_t inc=10);
  virtual Float_t GetX() const { return fX;}
  virtual void  SetX(Float_t x) { fX = x;}
  virtual Int_t GetDetector() const {return fDetector;}
  virtual Int_t GetRow() const {return fRow;}
  virtual void SetDetector(Int_t detector){fDetector = (UChar_t)(detector%256);}
  virtual void SetRow(Int_t row){fRow = (UChar_t)(row%256);}  
  virtual void SetTimeBin(Float_t timeBin){ fTimeBin= timeBin;}
  virtual void SetPad(Float_t pad){ fPad = pad;}
  //
  void SetQ(Float_t q) {fQ=(UShort_t)q;}
  void SetType(Char_t type) {fType=type;}
  void SetMax(UShort_t max) {fMax=max;}
  Int_t IsUsed(Int_t th=10) const {return (fUsed>=th) ? 1 : 0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}
  Float_t GetMax() const {return fMax;} 
  Char_t  GetType()const {return fType;}
  Float_t GetTimeBin() const { return fTimeBin;}
  Float_t GetPad() const { return fPad;}
  AliTPCclusterInfo * GetInfo() const { return fInfo;}
  void SetInfo(AliTPCclusterInfo * info) { fInfo = info;}

private:
  Float_t   fX;        //X position of cluster
  Short_t   fQ ;       //Q of cluster (in ADC counts)  
  Char_t    fType;     //type of the cluster 0 means golden 
  Short_t   fMax;      //maximal amplitude in cluster
  Char_t    fUsed;     //counter of usage  
  UChar_t   fDetector; //detector  number
  UChar_t   fRow;      //row number number
  Float_t   fTimeBin;  //time bin coordinate
  Float_t   fPad;  //pad coordinate
  AliTPCclusterInfo * fInfo;  // pointer to the cluster debug info
  ClassDef(AliTPCclusterMI,3)  // Time Projection Chamber clusters
};

void AliTPCclusterMI::Use(Int_t inc) 
{ 
  if (inc>0)  fUsed+=inc; 
  else 
    fUsed=0;
}



#endif


