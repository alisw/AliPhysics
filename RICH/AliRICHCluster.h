#ifndef AliRICHCluster_h
#define AliRICHCluster_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>         //base class
#include <TVector2.h>        //DistTo
#include "AliRICHDigit.h"
class TMinuit;

class AliRICHCluster :public TObject
{
public:
  enum EClusterStatus {kFormed,kCoG,kUnfolded,kEmpty}; 
  AliRICHCluster()                                               :TObject(),fCFM(-1),fSize(-1),fShape(-1),fQdc(-1),fChamber(-1),fX(-1),fY(-1),fStatus(kEmpty),fDigits(0) {}  //default ctor  
  AliRICHCluster(Int_t cs,Double_t x,Double_t y,Int_t q,Int_t sm):TObject(),fCFM(-1),fSize(sm),fShape(-1),fQdc(q ),fChamber(cs),fX(x ),fY(y ),fStatus(kEmpty),fDigits(0) {}  //default ctor  
  virtual          ~AliRICHCluster()                 {}                                                  //dtor
                   // AliRICHcluster(const AliRICHcluster& clus):TObject(clus)                                                         {}  //copy ctor 
  AliRICHCluster&  operator=(const AliRICHCluster&)                 {return *this;}                                                     //copy operator
                   
         void      Print(Option_t *option="")const;                                                       //
  
  
         void       Reset()                          {DeleteDigits();fCFM=fSize=fShape=fQdc=fChamber=-1;fX=fY=-1;fStatus=kEmpty;} //cleans the cluster
         void       DeleteDigits()                   {if(fDigits) {delete fDigits;} fDigits=0;}           //deletes the list of digits  
         Int_t      Nlocmax()                   const{return fSize-10000*(fSize/10000);}                //number of local maximums
         Int_t      Size()                      const{return fSize/10000;}                              //number of digits in cluster
         Int_t      Fsize()                     const{return fSize;}                                    //
         Int_t      Shape()                     const{return fShape;}                                   //cluster shape rectangulare
         Int_t      C()                         const{return fChamber/10;}                              //chamber number
         Int_t      S()                         const{return fChamber-(fChamber/10)*10;}                //sector number
         Int_t      Fchamber()                  const{return fChamber;}                                 //
         Int_t      Q()                         const{return fQdc;}                                     //cluster charge in QDC channels 
         Double_t   X()                         const{return fX;}                                       //cluster x position in LRS
         Double_t   Y()                         const{return fY;}                                       //cluster y position in LRS 
         Int_t      Status()                    const{return fStatus;}                                  //
         void       SetStatus(Int_t status)         {fStatus=status;}                                     //
         Int_t      Nmips()                     const{return fCFM-1000000*Ncerenkovs()-1000*Nfeedbacks();} //
         Int_t      Ncerenkovs()                const{return fCFM/1000000;}                                //
         Int_t      Nfeedbacks()                const{return (fCFM-1000000*Ncerenkovs())/1000;}            //
         Bool_t     IsPureMip()                 const{return fCFM<1000;}                                   //
         Bool_t     IsPureCerenkov()            const{return Nmips()==0&&Nfeedbacks()==0;}                 //
         Bool_t     IsPureFeedback()            const{return Nmips()==0&&Ncerenkovs()==0;}                 //
         Bool_t     IsSingleMip()               const{return Nmips()==1&&Ncerenkovs()==0&&Nfeedbacks()==0;}  //
         Bool_t     IsSingleCerenkov()          const{return Nmips()==0&&Ncerenkovs()==1&&Nfeedbacks()==0;}  //
         Bool_t     IsSingleFeedback()          const{return Nmips()==0&&Ncerenkovs()==0&&Nfeedbacks()==1;}  //
         Bool_t     IsMip()                     const{return Nmips()!=0;}                                  //
         Bool_t     IsCerenkov()                const{return Ncerenkovs()!=0;}                             //
         Bool_t     IsFeedback()                const{return Nfeedbacks()!=0;}                             //
         Int_t      CombiPid()                  const{return fCFM;}                                        //
         void       CFM(Int_t c,Int_t f,Int_t m)     {fCFM=1000000*c+1000*f+m;}                            //cluster contributors
         TObjArray*    Digits()                    const{return fDigits;}                                     //  
         
  inline void          AddDigit(AliRICHDigit *pDig);                                                          //add new digit ot the cluster
         AliRICHDigit* Digit  (Int_t i                      )const{return (AliRICHDigit*)fDigits->At(i); }//get pointer to i-th digit without existence check 
         TMinuit*      Solve  (                             );                                            //calculates cluster position
         TMinuit*      Unfold (                             );                                            //decompose cluster n. loc max clusters
         void          Set    (Double_t x,Double_t y,Int_t q)     {fX=x;fY=y,fQdc=q;                     }//set some cluster properties
  static void          FitFunc(Int_t &iNpars, Double_t *, Double_t &chi2, Double_t *aPar, Int_t);         //fit function to be used by MINUIT
          
         void      CoG(Int_t iNlocmax);                                                                   //calculates center of gravity
          void      Fill(AliRICHCluster *pRaw,Double_t x,Double_t y,Double_t q,Int_t cfm)              //form new resolved cluster from raw one
                    {fCFM=cfm;fChamber=pRaw->Fchamber();fSize=pRaw->Fsize();fQdc=(Int_t)(q*pRaw->Q());fX=x;fY=y;fStatus=kUnfolded;} 
         Double_t   DistTo(TVector2 x)          const{return TMath::Sqrt((x.X()-fX)*(x.X()-fX)+(x.Y()-fY)*(x.Y()-fY));} //distance to given point 
         Double_t   DistX(TVector2 x)           const{return (x.X()-fX);} //distance in x to given point 
         Double_t   DistY(TVector2 x)           const{return (x.Y()-fY);} //distance to given point 
         void       Test(const TVector2 &x,Double_t dEloss=0);            //test cluster fuctionality by provided hit with energy in eV
protected:
  Int_t         fCFM;         //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t         fSize;        //10000*(N digits) + N maxima     
  Int_t         fShape;       //100*xdim+ydim box containing the cluster
  Int_t         fQdc;         //QDC value
  Int_t         fChamber;     //10*chamber number+sector number 
  Double_t      fX;           //local x postion 
  Double_t      fY;           //local y postion  
  Int_t         fStatus;      //flag to mark the quality of the cluster   
  TObjArray    *fDigits;      //! list of digits forming this cluster
  ClassDef(AliRICHCluster,2)  //RICH cluster class       
};//class AliRICHCluster
//__________________________________________________________________________________________________
void AliRICHCluster::AddDigit(AliRICHDigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster    
  if(!fDigits) {fQdc=fSize=0;fDigits = new TObjArray;}
  fDigits->Add(pDig);
  fQdc+=(Int_t)pDig->Qdc(); 
  fChamber=10*pDig->Chamber()+pDig->Sector();
  fSize+=10000;
  fStatus=kFormed;
}
//__________________________________________________________________________________________________
#endif
