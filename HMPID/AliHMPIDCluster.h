#ifndef AliHMPIDCluster_h
#define AliHMPIDCluster_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHMPIDDigit.h"  //DigAdd()
#include <TObjArray.h>     //DigAdd()      
class TClonesArray;        //Solve()

class AliHMPIDCluster :public TObject
{
public:
  enum EClusterStatus {kFor,kCoG,kUnf,kEmp=-1};      //status flags    
  AliHMPIDCluster(                                   ):TObject( ),fSt(kEmp ),fCh(-1   ),fQ(-1  ),fX(-1  ),fY(-1  ),fDigs(0                                  ) {} 
  AliHMPIDCluster(Int_t c,Float_t x,Float_t y,Int_t q):TObject( ),fSt(kUnf ),fCh(c    ),fQ(q   ),fX(x   ),fY(y   ),fDigs(0                                  ) {} 
  AliHMPIDCluster(const AliHMPIDCluster &c            ):TObject(c),fSt(c.fSt),fCh(c.fCh),fQ(c.fQ),fX(c.fX),fY(c.fY),fDigs(c.fDigs ? new TObjArray(*c.fDigs):0) {}
  AliHMPIDCluster &operator=(const AliHMPIDCluster &c) {
    if(this == &c)return *this;TObject::operator=(c);            fSt=c.fSt; fCh=c.fCh; fQ=c.fQ; fX=c.fX; fY=c.fY; fDigs=c.fDigs ? new TObjArray(*c.fDigs):0; return *this;}
    
  virtual ~AliHMPIDCluster(                                     )                                                                        {DigDel();}
//framework part                   
         void          Print  (Option_t *opt=""                                  )const;                  //overloaded TObject::Print() to print cluster info
  static void          FitFunc(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);                       //fit function to be used by MINUIT
//private part  
         void          CoG      (                                         );                                                      //calculates center of gravity
         void          CorrSin  (                                         );                                                      //sinoidal correction   
         Int_t         Ch       (                                         )const{return fCh;                                    } //chamber number
  inline void          DigAdd   (AliHMPIDDigit *pDig                       );                                                      //add new digit ot the cluster
         void          DigDel   (                                         )     {if(fDigs) {delete fDigs;fDigs=0;}              } //deletes the list of digits (not digits!) 
         AliHMPIDDigit* Dig      (Int_t i                                  )const{return (AliHMPIDDigit*)fDigs->At(i);            } //pointer to i-th digit
         TObjArray*    DigLst   (                                         )const{return fDigs;                                  } //list of digits  
         void          Reset    (                                         )     {DigDel();fQ=fCh=-1;fX=fY=-1;fSt=kEmp;          } //cleans the cluster
         Int_t         Solve    (TClonesArray *pCluLst,Bool_t isUnfold    );                                                      //solve cluster: MINUIT fit or CoG
         Int_t         Size     (                                         )const{return (fDigs)?fDigs->GetEntriesFast():0;      } //number of pads in cluster
         Int_t         Q        (                                         )const{return fQ;                                     } //cluster charge in QDC channels 
         Float_t       X        (                                         )const{return fX;                                     } //cluster x position in LRS
         Float_t       Y        (                                         )const{return fY;                                     } //cluster y position in LRS 
protected:
  Int_t         fSt;          //flag to mark the quality of the cluster   
  Int_t         fCh;          //chamber number
  Int_t         fQ;           //QDC value
  Float_t       fX;           //local x postion, [cm] 
  Float_t       fY;           //local y postion, [cm]  
  TObjArray    *fDigs;        //! list of digits forming this cluster
  ClassDef(AliHMPIDCluster,5)  //HMPID cluster class       
};//class AliHMPIDCluster

typedef AliHMPIDCluster AliRICHCluster; // for backward compatibility

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::DigAdd(AliHMPIDDigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster, cluster is not owner of digits
// Arguments: pDig - pointer to digit to be added  
//   Returns: none  
  if(!fDigs) {fQ=0;fDigs = new TObjArray;}
  fDigs->Add(pDig);
  fQ+=(Int_t)pDig->Q(); 
  fCh=pDig->Ch();
  fSt=kFor;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
