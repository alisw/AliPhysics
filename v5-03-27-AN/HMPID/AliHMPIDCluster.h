#ifndef AliHMPIDCluster_h
#define AliHMPIDCluster_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//
// Implementation class: AliHMPIDCluster
//   
//   class to reconstruct clester in HMPID
//   it is forseen to Solve (split) the raw cluster in 
//   several clusters (# local maxima in the raw cluster -> deconvolution
//   according to a Mathieson profile of the charge 
//
#include "AliGeomManager.h"
#include "AliCluster3D.h"
#include "AliHMPIDDigit.h"  //DigAdd()
#include <TObjArray.h>     //DigAdd()      
class TClonesArray;        //Solve()

class AliHMPIDCluster :public AliCluster3D
{
public:
  enum EClusterStatus {kFrm,kCoG,kLo1,kUnf,kMax,kNot,kEdg,kSi1,kNoLoc,kAbn,kBig,kEmp=-1};      //status flags    
      AliHMPIDCluster():AliCluster3D(),
                          fCh(-1),fSi(-1),fSt(kEmp),fBox(-1),fNlocMax(-1),fMaxQpad(-1),fMaxQ(-1),fQRaw(0),
                          fQ(0),fErrQ(-1),fXX(0),fErrX(-1),fYY(0),fErrY(-1),fChi2(-1),fDigs(0),fParam(AliHMPIDParam::Instance())
      {
      }//ctor
  

      AliHMPIDCluster(const AliHMPIDCluster &c):AliCluster3D(c),
                        fCh(c.fCh),fSi(c.fSi),fSt(c.fSt),fBox(c.fBox),fNlocMax(c.fNlocMax),fMaxQpad(c.fMaxQpad),fMaxQ(c.fMaxQ),fQRaw(c.fQRaw),
                        fQ (c.fQ ),fErrQ(c.fErrQ),
                        fXX (c.fXX ),fErrX(c.fErrX),
                        fYY (c.fYY ),fErrY(c.fErrY),fChi2(c.fChi2),fDigs(0),fParam(c.fParam)                  {}//copy ctor
   virtual ~AliHMPIDCluster();//dtor   {if(fDigs) delete fDigs; fDigs=0;}
//framework part                   
         void           Draw   (Option_t *opt=""                                  );                       //overloaded TObject::Print() to draw cluster in current canvas
         void           Print  (Option_t *opt=""                                  )const;                  //overloaded TObject::Print() to print cluster info
  static void           FitFunc(Int_t &iNpars, Double_t* deriv, Double_t &chi2, Double_t *par, Int_t iflag);//fit function to be used by MINUIT
//private part  
         Int_t          Box      (                                         )const{return fBox;                                   } //Dimension of the cluster
         void           CoG      (                                         );                                                      //calculates center of gravity
         void           CorrSin  (                                         );                                                      //sinoidal correction   
         Int_t          Ch       (                                         )const{return fCh;                                    } //chamber number
  inline void           DigAdd   (AliHMPIDDigit *pDig                      );                                                      //add new digit ot the cluster
         AliHMPIDDigit* Dig      (Int_t i                                  )const{return (AliHMPIDDigit*)fDigs->At(i);           } //pointer to i-th digi 
  inline Bool_t         IsInPc   ();                                                                                               //check if is in the current PC
  inline void           Reset    (                                         );                                                      //cleans the cluster
  void           SetClusterParams(Double_t xL,Double_t yL,Int_t iCh  );                                                            //Set AliCluster3D part
         Int_t          Size     (                                         )const{return fSi;                                    } //returns number of pads in formed cluster 
         Int_t          Solve    (TClonesArray *pCluLst,Int_t *pSigmaCut, Bool_t isUnfold);                                        //solve cluster: MINUIT fit or CoG
         Int_t          Status   (                                         ) const{return fSt;}                                    //Status of cluster                                  
         Double_t       QRaw     (                                         )const{return fQRaw;                                  } //raw cluster charge in QDC channels 
         Double_t       Q        (                                         )const{return fQ;                                     } //given cluster charge in QDC channels 
         Double_t       Qe       (                                         )const{return fErrQ;                                  } //Error in cluster charge in QDC channels 
         Double_t       X        (                                         )const{return fXX;                                    } //cluster x position in LRS
         Double_t       Xe       (                                         )const{return fErrX;                                  } //cluster charge in QDC channels 
         Double_t       Y        (                                         )const{return fYY;                                    } //cluster y position in LRS 
         Double_t       Ye       (                                         )const{return fErrY;                                  } //cluster charge in QDC channels 
         Double_t       Chi2     (                                         )const{return fChi2;                                  } //chi2 of the fit
         void           DoCorrSin(Bool_t doCorrSin                         ){fgDoCorrSin=doCorrSin;}                               // Set sinoidal correction
         void           SetX     (Double_t x                               ){fXX=x;}                                               // Setter
         void           SetY     (Double_t y                               ){fYY=y;}                                               // Setter
         void           SetQ     (Double_t q                               ){fQ=q;if(fQ>4095)fQ=4095;}                             // Setter         
         void           SetQRaw  (Double_t qRaw                            ){fQRaw=qRaw;if(fQRaw>4095)fQRaw=4095;}                 // Setter         
         void           SetSize  (Int_t size                               ){fSi=size;}                                            // Setter
         void           SetCh    (Int_t chamber                            ){fCh=chamber;}                                         // Setter
         void           SetChi2  (Double_t chi2                            ){fChi2=chi2;}                                          // Setter
         void           SetStatus(Int_t status                             ){fSt=status;}                                              // Setter
         void           FindClusterSize(Int_t i,Int_t *pSigmaCut);                                                                 //Find the clusterSize of deconvoluted clusters 
 virtual void	        Clear(const Option_t*) { delete [] fDigs; fDigs=0; delete [] fParam; fParam=0; }
         
protected:
  Int_t         fCh;          //chamber number
  Int_t         fSi;          //size of the formed cluster from which this cluster deduced
  Int_t         fSt;          //flag to mark the quality of the cluster
  Int_t         fBox;         //box contaning this cluster
  Int_t         fNlocMax;     //number of local maxima in formed cluster
  Int_t         fMaxQpad;     //abs pad number of a pad with the highest charge
  Double_t      fMaxQ;        //that max charge value
  Double_t      fQRaw;        //QDC value of the raw cluster
  Double_t      fQ;           //QDC value of the actual cluster
  Double_t      fErrQ;        //error on Q
  Double_t      fXX;           //local x postion, [cm]
  Double_t      fErrX;        //error on x postion, [cm]
  Double_t      fYY;           //local y postion, [cm]
  Double_t      fErrY;        //error on y postion, [cm]
  Double_t      fChi2;        //some estimator of the fit quality
  TObjArray    *fDigs;        //! list of digits forming this cluster
  static  Bool_t fgDoCorrSin; //flag to switch on/off correction for Sinusoidal to cluster reco
  AliHMPIDParam *fParam;      //!Pointer to AliHMPIDParam
  
private:
/*
  AliHMPIDCluster &operator=(const AliHMPIDCluster &c) {if(this == &c)return *this;AliCluster3D::operator=(c);          
                                                       fSi=c.fSi;  fSt=c.fSt; fCh=c.fCh; fBox=c.fBox;fNlocMax=c.fNlocMax;fMaxQpad=c.fMaxQpad; fMaxQ=c.fMaxQ;fQRaw=c.fQRaw;
                                                        fQ=c.fQ; fErrQ=c.fErrQ; 
                                                        fXX=c.fXX; fErrX=c.fErrX;
                                                        fYY=c.fYY; fErrY=c.fErrY; fChi2=c.fChi2;fDigs=c.fDigs ? new TObjArray(*c.fDigs):0; return *this;}
*/
//  AliHMPIDCluster(const AliHMPIDCluster& c);              //dummy copy constructor
  AliHMPIDCluster &operator=(const AliHMPIDCluster& c);   //dummy assignment operator
      
  
  ClassDef(AliHMPIDCluster,9) //HMPID cluster class
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::DigAdd(AliHMPIDDigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster, cluster is not owner of digits
// Arguments: pDig - pointer to digit to be added  
// Returns: none  
  if(!fDigs) {fSi=0;fDigs = new TObjArray;} //create list of digits in the first invocation
  fDigs->Add(pDig);
  fSt=kFrm;
  fSi++;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDCluster::Reset()
{
  //
  //
  //
  if(fDigs) delete fDigs; 
  fDigs=0; 
  fSt=kEmp;
  fQRaw=fQ=0;
  fXX=fYY=0;
  fCh=fSi=fBox=fNlocMax=fMaxQpad=-1;
  fMaxQ=fErrQ=fErrX=fErrY=fChi2=-1; //empty ctor
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDCluster::IsInPc()
{
  //Check if (X,Y) position is inside the PC limits
  //Arguments:
  //  Returns: True or False
  Int_t pc = ((AliHMPIDDigit*)fDigs->At(0))->Pc();
 
  
  if ( fXX < AliHMPIDParam::MinPcX(pc) || fXX > AliHMPIDParam::MaxPcX(pc) || 
       fYY < AliHMPIDParam::MinPcY(pc) || fYY > AliHMPIDParam::MaxPcY(pc) ) return kFALSE;
  
  return kTRUE;
  
}

#endif
