#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TVector3.h>
#include <Riostream.h>
#include <AliDetector.h>
#include <AliHit.h>
#include <AliDigit.h>
#include "AliRICHDigitizer.h"
#include "AliRICHParam.h"

//__________________AliRICHhit______________________________________________________________________
class AliRICHhit : public AliHit
{
public:
           AliRICHhit():AliHit()          {fChamber=kBad;  fEloss=kBad;  fInX3.SetXYZ(0,0,0);fOutX3.SetXYZ(0,0,0);}
           AliRICHhit(Int_t c,Int_t tid,TVector3 in,TVector3 out,Double_t e):AliHit(0,tid)
                                          {fChamber=c;fInX3=in; fOutX3=out;fEloss=e; fX=out.X();fY=out.Y();fZ=out.Z();}
  virtual ~AliRICHhit()                   {;}

  Int_t    C()                       const{return fChamber;}              //chamber number 
  Int_t    Chamber()                 const{return fChamber;}              //chamber number 
  Float_t  Eloss()                   const{return fEloss;}                //energy lost by track inside amplification gap  
  TVector3 InX3()                    const{return fInX3;}                 //track position at the faceplane of the gap 
  TVector3 OutX3()                   const{return fOutX3;}                //track position at the backplane of the gap 
  Double_t Length()                  const{return (fOutX3-fInX3).Mag();}  //track length inside the amplification gap
  void     Print(Option_t *option="")const;                               //virtual
protected:
  Int_t     fChamber;                      //chamber number
  Double_t  fEloss;                        //ionisation energy lost in GAP
  TVector3  fInX3;                         //position at the entrance of the GAP   
  TVector3  fOutX3;                        //position at the exit of the GAP
  ClassDef(AliRICHhit,2)                   //RICH hit class
};//class AliRICHhit

//__________________AliRICHdigit____________________________________________________________________
class AliRICHdigit :public AliDigit
{
public:
  AliRICHdigit() {fChFbMip=fChamber=fPadX=fPadY=fTracks[0]=fTracks[1]=fTracks[2]=kBad;fQdc=kBad;}
  AliRICHdigit(Int_t c,Int_t x,Int_t y,Double_t q,Int_t cpid,Int_t tid0,Int_t tid1,Int_t tid2)  
                 {fPadX=x;fPadY=y;fQdc=q;fChamber=10*c+AliRICHParam::Sector(x,y);fChFbMip=cpid;fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}
  virtual ~AliRICHdigit() {;}  
  Int_t    Compare(const TObject *pObj) const
                 {if(Id()==((AliRICHdigit*)pObj)->Id())return 0;else if(Id()>((AliRICHdigit*)pObj)->Id())return 1;else return -1;}  //virtual      
  Bool_t   IsSortable()                 const{return kTRUE;}                              //virtual
  void     Print(Option_t *option="")   const;                                            //virtual
  Int_t    ChFbMi()                     const{return fChFbMip;}                           //particle mixture for this digit
  Int_t    C()                          const{return fChamber/10;}                        //chamber number
  Int_t    S()                          const{return fChamber-(fChamber/10)*10;}          //sector number
  Int_t    X()                          const{return fPadX;}                              //x position of the pad
  Int_t    Y()                          const{return fPadY;}                              //y postion of the pad
  Int_t    Id()                         const{return fChamber*10000000+fPadX*1000+fPadY;} //absolute id of this pad
  Double_t Q()                          const{return fQdc;}                               //charge in terms of ADC channels
  Int_t    Tid(Int_t i)                 const{return fTracks[i];}                         //track reference produced this digit
  void     AddTidOffset(Int_t offset) 
    {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;};
protected:
  Int_t    fChFbMip;  //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //10*chamber number+ sector number 
  Int_t    fPadX;     //pad number along X
  Int_t    fPadY;     //pad number along Y
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHdigit,3) //RICH digit class       
};//class AliRICHdigit

//__________________AliRICHcluster__________________________________________________________________
class AliRICHcluster :public TObject
{
public:
  enum ClusterStatus {kEdge,kShape,kSize,kRaw,kResolved};
                    AliRICHcluster()                                 {fSize=fQdc=fStatus=fChamber=fDimXY=0;fX=fY=kBad;fDigits=0;}
  virtual          ~AliRICHcluster()                                 {delete fDigits;}  
  AliRICHcluster&   operator=(const AliRICHcluster&)                 {return *this;}
         Int_t      Nlocals()                                   const{return fSize - 10000*(fSize/10000);} //
         Int_t      Size()                                      const{return fSize/10000;}                 //
         Int_t      Fsize()                                     const{return fSize;}                       //
         Int_t      DimXY()                                     const{return fDimXY;}                      //
         Int_t      C()                                         const{return fChamber/10;}                 //
         Int_t      S()                                         const{return fChamber-(fChamber/10)*10;}   //
         Int_t      Fchamber()                                  const{return fChamber;}                    //
         Int_t      Q()                                         const{return fQdc;}                        // 
         Double_t   X()                                         const{return fX;}                          //
         Double_t   Y()                                         const{return fY;}                          // 
         Int_t      Status()                                    const{return fStatus;}                     //
         void       SetStatus(Int_t status)                          {fStatus=status;}                     //
         Int_t      Nmips()                                     const{return fCombiPid-1000000*Ncerenkovs()-1000*Nfeedbacks();}  //
         Int_t      Ncerenkovs()                                const{return fCombiPid/1000000;}           //
         Int_t      Nfeedbacks()                                const{return (fCombiPid-1000000*Ncerenkovs())/1000;}      //
         Bool_t     IsPureMip()                                 const{return fCombiPid<1000;}   
         Bool_t     IsPureCerenkov()                            const{return Nmips()==0&&Nfeedbacks()==0;} //
         Bool_t     IsPureFeedback()                            const{return Nmips()==0&&Ncerenkovs()==0;} //
         Int_t      CombiPid()                                  const{return fCombiPid;} //
         void       SetCombiPid(Int_t ckov,Int_t feeds,Int_t mips)   {fCombiPid=1000000*ckov+1000*feeds+mips;}            //
         void       Fill(AliRICHcluster *pRaw,Double_t x,Double_t y, Double_t q, Int_t combipid)
             {fCombiPid=combipid;fChamber=pRaw->Fchamber();fSize=pRaw->Fsize();
              fQdc=(Int_t)(q*pRaw->Q());fX=x;fY=y;fStatus=kResolved;} //
         TObjArray* Digits()                                    const{return fDigits;}                     //  
         void       Print(Option_t *option="")const;                                                       //virtual
  inline void       AddDigit(AliRICHdigit *pDig);                                                          //
  inline void       CoG(Int_t nLocals);                                                                                 // 
         void       Reset() {fSize=fQdc=fStatus=fChamber=fDimXY=kBad;fX=fY=kBad;delete fDigits;fDigits=0;} //
protected:
  Int_t         fCombiPid;    //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t         fSize;        //10000*(how many digits belong to this cluster) + nLocalMaxima     
  Int_t         fDimXY;       //100*xdim+ydim box containing the cluster
  Int_t         fQdc;         //QDC value
  Int_t         fChamber;     //10*module number+sector number 
  Double_t      fX;           //local x postion 
  Double_t      fY;           //local y postion  
  Int_t         fStatus;      //flag to mark the quality of the cluster   
  TObjArray    *fDigits;      //! list of digits forming this cluster
  ClassDef(AliRICHcluster,2)  //RICH cluster class       
};//class AliRICHcluster
//__________________________________________________________________________________________________
void AliRICHcluster::AddDigit(AliRICHdigit *pDig)
{//    
  if(!fDigits) {fQdc=fSize=fCombiPid=0;fDigits = new TObjArray;}
  fQdc+=(Int_t)pDig->Q(); fDigits->Add(pDig);
  fChamber=10*pDig->C()+pDig->S();
  fSize+=10000;
}
//__________________________________________________________________________________________________
void AliRICHcluster::CoG(Int_t nLocals)
{//
  Int_t xmin=999,ymin=999,xmax=0,ymax=0;   
  fX=fY=0;
  for(Int_t iDig=0;iDig<Size();iDig++) {
    AliRICHdigit *pDig=(AliRICHdigit*)fDigits->At(iDig);
    Int_t padX = pDig->X();Int_t padY = pDig->Y();Double_t q=pDig->Q();
    TVector2 x2=AliRICHParam::Pad2Loc(padX,padY);
    fX += x2.X()*q;fY +=x2.Y()*q;
    if(padX<xmin)xmin=padX;if(padX>xmax)xmax=padX;if(padY<ymin)ymin=padY;if(padY>ymax)ymax=padY;
   }
   fX/=fQdc;fY/=fQdc;//Center of Gravity
   fDimXY = 100*(xmax-xmin+1)+ymax-ymin+1;//find box containing cluster
   fSize+=nLocals;
   fStatus=kRaw;
}//CoG()


class AliRICHreco: public TObject
{
public:
            AliRICHreco() {fTid=fNphotons=kBad; fThetaCherenkov=kBad;}
            AliRICHreco(Int_t tid,Double_t thetaCherenkov,Int_t nPhotons) {fTid=tid;fThetaCherenkov=thetaCherenkov;fNphotons=nPhotons;}

  virtual  ~AliRICHreco() {;}

  void    Print(Option_t *option="")const;      //virtual print
  
protected:
   Int_t fTid;                 // track Id reference
   Int_t fNphotons;            // number of photons contributed to the recontruction
   Double_t fThetaCherenkov;   // reconstructed Theta Cerenkov for a given charged track

   ClassDef(AliRICHreco,1)  //RICH reco class

};//class AliRICHreco

//__________________AliRICH_________________________________________________________________________
class AliRICHParam;
class AliRICHChamber;

class AliRICH : public AliDetector 
{
public:
            AliRICH();                                            
            AliRICH(const char *name, const char *title);
            AliRICH(const AliRICH& RICH):AliDetector(RICH) {;}  //copy ctor 
  virtual  ~AliRICH();                                            
          
  AliRICH&  operator=(const AliRICH&)                 {return *this;}
//framework part  
  virtual Int_t   IsVersion()                                            const =0;                                //virtual         
  virtual void    StepManager()                                                =0;                                //virtual
          void    Hits2SDigits();                                                                                 //virtual
  AliDigitizer*   CreateDigitizer(AliRunDigitizer* man) const {return new AliRICHDigitizer(man);}                 //virtual
          void    Print(Option_t *option)const;                                                                   //virtual
          void    SetTreeAddress();                                                                               //virtual
          void    MakeBranch(Option_t *opt=" ");                                                                  //virtual
          void    CreateMaterials();                                                                              //virtual
  virtual void    BuildGeometry();                                                                                //virtual
  virtual void    CreateGeometry();                                                                               //virtual
//private part  
          Float_t AbsoCH4(Float_t x)const;                               //calculates absorption length for methane
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)const;  //deals with Fresnel absorption
          void    GenerateFeedbacks(Int_t iChamber,Float_t eloss=0);     //generates feedback photons; eloss=0 for photon
          void    CreateChambers();                                      //creates set of chambers
  inline  void    CreateHits();                                          //create hits container as a simple list
  inline  void    CreateSDigits();                                       //create sdigits container as a simple list
  inline  void    CreateDigits();                                        //create digits container as 7 lists, one per chamber
  inline  void    CreateClusters();                                      //create clusters container  as 7 lists, one per chamber
  inline  void    CreateRecos();                                         //create recos container
//        void    ResetHits()                {AliDetector::ResetHits();}  //virtual  
          void    ResetSDigits()             {fNsdigits=0;  if(fSdigits)  fSdigits ->Clear();}                                 
          void    ResetDigits()              {if(fDigitsNew)for(int i=0;i<kNCH;i++){fDigitsNew->At(i)->Clear();fNdigitsNew[i]=0;}} //virtual
          void    ResetClusters()            {if(fClusters) for(int i=0;i<kNCH;i++){fClusters ->At(i)->Clear();fNclusters[i]=0;}}
          void    ResetRecos()               {if(fRecos) fRecos->Clear();fNrecos=0;}
  TClonesArray*   SDigits()             const{return fSdigits;}
  TClonesArray*   Digits(Int_t iC)      const{if(fDigitsNew) return (TClonesArray *)fDigitsNew->At(iC-1);else return 0;}
  TClonesArray*   Clusters(Int_t iC)    const{if(fClusters)  return (TClonesArray *)fClusters->At(iC-1);else return 0;}
  TClonesArray*   Recos()               const{return fRecos;}
  AliRICHChamber* C(Int_t iC)           const{return (AliRICHChamber*)fChambers->At(iC-1);}
  AliRICHParam*   P()                   const{return fpParam;}
          void    PrintDigits()         const{for(Int_t i=0;i<kNCH;i++) fDigitsNew->At(i)->Print();}
          void    PrintClusters()       const{for(Int_t i=0;i<kNCH;i++) fClusters->At(i)->Print();}
            
  void AddHit(Int_t chamber,Int_t tid,TVector3 iX3,TVector3 oX3,Double_t eloss=0)
       {TClonesArray &tmp=*fHits;new(tmp[fNhits++])AliRICHhit(chamber,tid,iX3,oX3,eloss);} 
  inline void AddSDigit(Int_t c,Int_t x,Int_t y,Double_t q,Int_t pid,Int_t tid); 
  void AddDigit(int c,int x,int y,int q,int cfm,int *tid)
       {TClonesArray &tmp=*((TClonesArray*)fDigitsNew->At(c-1));new(tmp[fNdigitsNew[c-1]++])AliRICHdigit(c,x,y,q,cfm,tid[0],tid[1],tid[2]);}  
  void AddCluster(AliRICHcluster &cl)                     
       {Int_t c=cl.C()-1;/*cout<<c<<endl*/;TClonesArray &tmp=*((TClonesArray*)fClusters->At(c));new(tmp[fNclusters[c]++])AliRICHcluster(cl);}
  void AddReco(Int_t tid,Double_t thetaCherenkov,Int_t nPhotons) 
       {TClonesArray &tmp=*(TClonesArray*)fRecos;new(tmp[fNrecos++])AliRICHreco(tid,thetaCherenkov,nPhotons);}  
          
protected:  
  enum                  {kCSI=6,kGAP=9};
  AliRICHParam         *fpParam;             //main RICH parametrization     
  TObjArray            *fChambers;           //list of RICH chambers
  
                                             //fHits and fDigits belong to AliDetector
  TClonesArray         *fSdigits;            //! List of sdigits  
  Int_t                 fNsdigits;           //! Current number of sdigits
  
  TObjArray            *fDigitsNew;          //! Each chamber holds it's one lists of digits
  Int_t                 fNdigitsNew[kNCH];   //! Array of current numbers of digits
  
  TObjArray            *fClusters;           //! Each chamber holds it's one lists of clusters 
  Int_t                 fNclusters[kNCH];    //! Array of current numbers of raw clusters
  
  TClonesArray         *fRecos;              //! pointer to the list of recos
  Int_t                 fNrecos;             //! number of recos

  ClassDef(AliRICH,5)                        //Main RICH class 
};//class AliRICH  

//__________________________________________________________________________________________________
void AliRICH::CreateHits()
{
  if(fHits) return;
  if(GetDebug())Info("CreateHits","creating hits container.");
  fHits=new TClonesArray("AliRICHhit",10000);   fNhits=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateSDigits()
{
  if(fSdigits) return;
  if(GetDebug())Info("CreateSDigits","creating sdigits container.");
  fSdigits=new TClonesArray("AliRICHdigit",10000); fNsdigits=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateDigits()
{
  if(fDigitsNew) return;
  if(GetDebug())Info("CreateDigits","creating digits containers.");
  fDigitsNew = new TObjArray(kNCH);  
  for(Int_t i=0;i<kNCH;i++) {fDigitsNew->AddAt(new TClonesArray("AliRICHdigit",10000), i); fNdigitsNew[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::CreateClusters()
{
  if(fClusters) return;
  if(GetDebug())Info("CreateClusters","creating clusters containers.");
  fClusters = new TObjArray(kNCH);  
  for(Int_t i=0;i<kNCH;i++) {fClusters->AddAt(new TClonesArray("AliRICHcluster",10000), i); fNclusters[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::CreateRecos()
{
  if(fRecos) return;
  if(GetDebug())Info("CreateRecos","creating recos containers.");
  fRecos = new TClonesArray("AliRICHreco",1000);fNrecos=0;  
}
//__________________________________________________________________________________________________
void AliRICH::AddSDigit(Int_t c,Int_t x,Int_t y,Double_t q,Int_t pid,Int_t tid) 
{   
  switch(pid){
    case 50000050: pid=1000000;break;//cerenkov
    case 50000051: pid=1000;   break;//feedback
    default:       pid=1;      break;//mip
  }   
  TClonesArray &tmp=*fSdigits;
  new(tmp[fNsdigits++])AliRICHdigit(c,x,y,q,pid,tid,kBad,kBad);
}//AddSDigit()   
//__________________________________________________________________________________________________
#endif//#ifndef AliRICH_h
