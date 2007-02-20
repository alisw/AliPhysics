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
  enum EClusterStatus {kFrm,kCoG,kLo1,kUnf,kMax,kAbn,kNot,kEdg,kSi1,kNoLoc,kEmp=-1};      //status flags    
  AliHMPIDCluster():TObject( ),fCh(-1),fSi(-1),fSt(kEmp),fBox(-1),fNlocMax(-1),fMaxQpad(-1),fMaxQ(-1),fQ(-1),fQi(-1),fQe(-1),fX(-1),fXi(-1),fXe(-1),fY(-1),fYi(-1),fYe(-1),fChi2(-1),fDigs(0) {} //empty ctor
  
  
  AliHMPIDCluster           (const AliHMPIDCluster &c):TObject(c),fCh(c.fCh),fSi(c.fSi),fSt(c.fSt),fBox(c.fBox),fNlocMax(c.fNlocMax),fMaxQpad(c.fMaxQpad),fMaxQ(c.fMaxQ),
                                                                  fQ (c.fQ ),fQi(c.fQi),fQe(c.fQe),
                                                                  fX (c.fX ),fXi(c.fXi),fXe(c.fXe),
                                                                  fY (c.fY ),fYi(c.fYi),fYe(c.fYe),fChi2(c.fChi2),fDigs(0)                  {}//copy ctor
  AliHMPIDCluster &operator=(const AliHMPIDCluster &c) {if(this == &c)return *this;TObject::operator=(c);          
                                                        fSi=c.fSi;  fSt=c.fSt; fCh=c.fCh; fBox=c.fBox;fNlocMax=c.fNlocMax;fMaxQpad=c.fMaxQpad; fMaxQ=c.fMaxQ;
                                                        fQ=c.fQ; fQi=c.fQi;fQe=c.fQe; 
                                                        fX=c.fX; fXi=c.fXi;fXe=c.fXe;
                                                        fY=c.fY; fYi=c.fYi;fYe=c.fYe; fChi2=c.fChi2;fDigs=c.fDigs ? new TObjArray(*c.fDigs):0; return *this;}
    
  virtual ~AliHMPIDCluster(                                     )                                                                        {if(fDigs) delete fDigs; fDigs=0;}
//framework part                   
         void           Draw   (Option_t *opt=""                                  );                       //overloaded TObject::Print() to draw cluster in current canvas
         void           Print  (Option_t *opt=""                                  )const;                  //overloaded TObject::Print() to print cluster info
  static void           FitFunc(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);                       //fit function to be used by MINUIT
//private part  
         void           CoG      (                                         );                                                      //calculates center of gravity
         void           CorrSin  (                                         );                                                      //sinoidal correction   
         Int_t          Ch       (                                         )const{return fCh;                                    } //chamber number
  inline void           DigAdd   (AliHMPIDDigit *pDig                      );                                                      //add new digit ot the cluster
         AliHMPIDDigit* Dig      (Int_t i                                  )const{return (AliHMPIDDigit*)fDigs->At(i);           } //pointer to i-th digi 
  inline Bool_t         IsInPc   ();                                                                                        //check if is in the current PC
  inline void           Reset    (                                         );                                                      //cleans the cluster
         Int_t          Size     (                                         )const{return fSi;                                    } //returns number of pads in formed cluster 
         Int_t          Solve    (TClonesArray *pCluLst,Bool_t isUnfold    );                                                      //solve cluster: MINUIT fit or CoG
         Int_t          Status  (                                          ) const{return fSt;}                                   //Status of cluster                                  
         Double_t       Q        (                                         )const{return fQ;                                     } //cluster charge in QDC channels 
         Double_t       Qi       (                                         )const{return fQi;                                    } //cluster charge in QDC channels 
         Double_t       Qe       (                                         )const{return fQe;                                    } //cluster charge in QDC channels 
         Double_t       X        (                                         )const{return fX;                                     } //cluster x position in LRS
         Double_t       Xi       (                                         )const{return fXi;                                    } //cluster charge in QDC channels 
         Double_t       Xe       (                                         )const{return fXe;                                    } //cluster charge in QDC channels 
         Double_t       Y        (                                         )const{return fY;                                     } //cluster y position in LRS 
         Double_t       Yi       (                                         )const{return fYi;                                    } //cluster charge in QDC channels 
         Double_t       Ye       (                                         )const{return fYe;                                    } //cluster charge in QDC channels 
         Double_t       Chi2     (                                         )const{return fChi2;                                  } 
protected:
  Int_t         fCh;          //chamber number
  Int_t         fSi;          //size of the formed cluster from which this cluster deduced  
  Int_t         fSt;          //flag to mark the quality of the cluster   
  Int_t         fBox;         //box contaning this cluster  
  Int_t         fNlocMax;     //number of local maxima in formed cluster        
  Int_t         fMaxQpad;     //abs pad number of a pad with the highest charge
  Double_t      fMaxQ;        //that max charge value             
  Double_t      fQ;           //QDC value
  Double_t      fQi;          //initial Q of local maximum
  Double_t      fQe;          //error on Q 
  Double_t      fX;           //local x postion, [cm] 
  Double_t      fXi;          //initial x of local maximum, [cm] 
  Double_t      fXe;          //error on x postion, [cm] 
  Double_t      fY;           //local y postion, [cm]  
  Double_t      fYi;          //initial y of local maximum, [cm] 
  Double_t      fYe;          //error on y postion, [cm] 
  Double_t      fChi2;        //some estimatore of fit quality  
  TObjArray    *fDigs;        //! list of digits forming this cluster
  ClassDef(AliHMPIDCluster,5) //HMPID cluster class       
};//class AliHMPIDCluster

typedef AliHMPIDCluster AliRICHCluster; // for backward compatibility

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::DigAdd(AliHMPIDDigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster, cluster is not owner of digits
// Arguments: pDig - pointer to digit to be added  
//   Returns: none  
  if(!fDigs) {fSi=0;fDigs = new TObjArray;} //create list of digits in the first invocation
  fDigs->Add(pDig);
  fSt=kFrm;
  fSi++;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::Reset()
{
  if(fDigs) delete fDigs; 
  fDigs=0; 
  fQ=fCh=fSi=-1;fX=fY=-1;fSt=kEmp;      
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCluster::IsInPc()
{
  Int_t pc = ((AliHMPIDDigit*)fDigs->At(0))->Pc();
 
  
  if ( fX < AliHMPIDDigit::fMinPcX[pc] || fX > AliHMPIDDigit::fMaxPcX[pc] || 
       fY < AliHMPIDDigit::fMinPcY[pc] || fY > AliHMPIDDigit::fMaxPcY[pc] ) return kFALSE;
  
  return kTRUE;
  
}

#endif
