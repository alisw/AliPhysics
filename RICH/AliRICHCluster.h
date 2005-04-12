#ifndef AliRICHCluster_h
#define AliRICHCluster_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TVector.h>
#include <TVector2.h>
#include "AliRICHDigit.h"

class AliRICHCluster :public TObject
{
public:
  enum ClusterStatus {kEdge,kShape,kSize,kRaw,kResolved,kEmpty}; 
                    AliRICHCluster():TObject(),fCFM(0),fSize(0),fShape(0),fQdc(0),fChamber(0),fX(0),fY(0),fStatus(kEmpty),fDigits(0) {}  //default ctor  
  virtual          ~AliRICHCluster()                 {AliDebug(1,"Start");/*Reset();*/}                                                  //dtor
                   // AliRICHcluster(const AliRICHcluster& clus):TObject(clus)                                                         {}  //copy ctor 
   AliRICHCluster&  operator=(const AliRICHCluster&)                 {return *this;}                                                     //copy operator
                   
  
  
         void       Reset()                          {DeleteDigits();fCFM=fSize=fShape=fQdc=fChamber=0;fX=fY=0;fStatus=kEmpty;} //cleans the cluster
         void       DeleteDigits()                   {if(fDigits) {delete fDigits;} fDigits=0;}           //deletes the list of digits  
         Int_t      Nlocals()                   const{return fSize-10000*(fSize/10000);}                //number of local maximums
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
         TObjArray* Digits()                    const{return fDigits;}                                     //  
  virtual void      Print(Option_t *option="")const;                                                   //
  inline  void      AddDigit(AliRICHDigit *pDig);                                                      //
  inline  void      CoG(Int_t nLocals);                                                                //calculates center of gravity
          void      Fill(AliRICHCluster *pRaw,Double_t x,Double_t y,Double_t q,Int_t cfm)              //form new resolved cluster from raw one
                    {fCFM=cfm;fChamber=pRaw->Fchamber();fSize=pRaw->Fsize();fQdc=(Int_t)(q*pRaw->Q());fX=x;fY=y;fStatus=kResolved;} 
         Double_t   DistTo(TVector2 x)          const{return TMath::Sqrt((x.X()-fX)*(x.X()-fX)+(x.Y()-fY)*(x.Y()-fY));} //distance to given point 
         Double_t   DistX(TVector2 x)           const{return (x.X()-fX);} //distance in x to given point 
         Double_t   DistY(TVector2 x)           const{return (x.Y()-fY);} //distance to given point 
protected:
  Int_t         fCFM;         //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t         fSize;        //10000*(how many digits belong to this cluster) + nLocalMaxima     
  Int_t         fShape;       //100*xdim+ydim box containing the cluster
  Int_t         fQdc;         //QDC value
  Int_t         fChamber;     //10*module number+sector number 
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
  if(!fDigits) {fQdc=fSize=fCFM=0;fDigits = new TObjArray;}
  fQdc+=(Int_t)pDig->Q(); fDigits->Add(pDig);
  fChamber=10*pDig->C()+pDig->S();
  fSize+=10000;
  fStatus=kRaw;
}
//__________________________________________________________________________________________________
void AliRICHCluster::CoG(Int_t nLocals)
{
// Calculates naive cluster position as a center of gravity of its digits.
  Float_t xmin=999,ymin=999,xmax=0,ymax=0;   
  fX=fY=0;
  for(Int_t iDig=0;iDig<Size();iDig++) {
    AliRICHDigit *pDig=(AliRICHDigit*)fDigits->At(iDig);
    TVector pad=pDig->Pad(); Double_t q=pDig->Q();
    TVector2 x2=AliRICHParam::Pad2Loc(pad);
    fX += x2.X()*q;fY +=x2.Y()*q;
    if(pad[0]<xmin)xmin=pad[0];if(pad[0]>xmax)xmax=pad[0];if(pad[1]<ymin)ymin=pad[1];if(pad[1]>ymax)ymax=pad[1];
   }
   fX/=fQdc;fY/=fQdc;//Center of Gravity

   TVector2 center = AliRICHParam::Pad2Loc(AliRICHParam::Loc2Pad(TVector2(fX,fY)));
   fX += AliRICHParam::CogCorr(fX-center.X());

   fShape=Int_t(100*(xmax-xmin+1)+ymax-ymin+1);//find box containing cluster
   fSize+=nLocals;
   fStatus=kRaw;
}//CoG()
//__________________________________________________________________________________________________
#endif
