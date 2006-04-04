#ifndef AliRICHCluster_h
#define AliRICHCluster_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICHParam.h"    //DigAdd(), Dig()

class AliRICHCluster :public TObject
{
public:
  enum EClusterStatus {kFormed,kCoG,kUnfolded,kEmpty=-1}; 
  AliRICHCluster():
    TObject(),
    fQdc(-1),
    fCham(-1),
    fX(-1),
    fY(-1),
    fStatus(kEmpty),
    fDigs(0) {
    // Default constructor
  }
  AliRICHCluster(Int_t c,Double_t x,Double_t y,Int_t q):
    TObject(),
    fQdc(q ),
    fCham(c) ,
    fX(x),
    fY(y),
    fStatus(kUnfolded),
    fDigs(0x0) {
    // Constructor
  }
  AliRICHCluster(const AliRICHCluster & src) :
    TObject(src),
    fQdc(src.fQdc),
    fCham(src.fCham),
    fX(src.fX),
    fY(src.fY),
    fStatus(src.fStatus),
    fDigs(src.fDigs ? new TObjArray(*src.fDigs) : 0x0) {
    // Copy constructor
  }
  AliRICHCluster & operator=(const AliRICHCluster & src) {
    // Assigment operator  
    if ( this == &src ) return *this;

    // Base class assignment
    TObject::operator=(src);

    fQdc    = src.fQdc;
    fCham   = src.fCham;
    fX      = src.fX;
    fY      = src.fY;
    fStatus = src.fStatus;
    fDigs = src.fDigs ? new TObjArray(*src.fDigs) : 0x0;
    return *this;
  }
  virtual             ~AliRICHCluster(                                     )                                                                        {DigDel();}
//framework part                   
         void          Print  (Option_t *opt=""                                  )const;                  //overloaded TObject::Print() to print cluster info
  static void          FitFunc(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);                       //fit function to be used by MINUIT
//private part  
         void          CoG      (                                         );                                                      //calculates center of gravity
         Int_t         C        (                                         )const{return fCham;                                  } //chamber number
  inline void          DigAdd   (AliRICHDigit *pDig                       );                                                      //add new digit ot the cluster
         void          DigDel   (                                         )     {if(fDigs) {delete fDigs;fDigs=0;}              } //deletes the list of digits (not digits!) 
         void          DistXY   (const TVector2 &p,Double_t &x,Double_t &y)const{ x=p.X()-fX; y=p.Y()-fY;                       } //distance in x to given point 
         AliRICHDigit* Dig      (Int_t i                                  )const{return (AliRICHDigit*)fDigs->At(i);            } //pointer to i-th digit
         TObjArray*    Digits   (                                         )const{return fDigs;                                  } //list of digits  
         TVector3      Lors2Mars()                            const{return AliRICHParam::Instance()->Lors2Mars(fCham,fX,fY);     } //cluster position in MARS
         void          Reset    (                                         )     {DigDel();fQdc=fCham=-1;fX=fY=-1;fStatus=kEmpty;} //cleans the cluster
         Int_t         Solve    (TClonesArray *pCluLst,Bool_t isUnfold    );                                                      //solve cluster: MINUIT fit or CoG
         Int_t         Size     (                                         )const{return (fDigs)?fDigs->GetEntriesFast():0;      } //number of pads in cluster
         Int_t         Q        (                                         )const{return fQdc;                                   } //cluster charge in QDC channels 
         Double_t      X        (                                         )const{return fX;                                     } //cluster x position in LRS
         Double_t      Y        (                                         )const{return fY;                                     } //cluster y position in LRS 
//test part         
  static void          Test     (Double_t x,Double_t y,Double_t e=0,Bool_t isUnfold=kTRUE);                                        //test by hit (x [cm] , y [cm] , e [GeV])
  static void          Test     (                                                        );                                        //test by predifined patterns
protected:
  Int_t         fQdc;         //QDC value
  Int_t         fCham;        //10*chamber number+sector number 
  Double_t      fX;           //local x postion, cm 
  Double_t      fY;           //local y postion, cm  
  Int_t         fStatus;      //flag to mark the quality of the cluster   
  TObjArray    *fDigs;        //! list of digits forming this cluster
  ClassDef(AliRICHCluster,3)  //RICH cluster class       
};//class AliRICHCluster
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliRICHCluster::DigAdd(AliRICHDigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster    
  if(!fDigs) {fQdc=0;fDigs = new TObjArray;}
  fDigs->Add(pDig);
  fQdc+=(Int_t)pDig->Qdc(); 
  fCham=pDig->C();
  fStatus=kFormed;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
