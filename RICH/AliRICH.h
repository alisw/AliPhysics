#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <Riostream.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TVector.h>
#include <TVector3.h>

#include "AliRICHDigitizer.h"
#include "AliRICHParam.h"
#include <AliDetector.h>
#include <AliDigit.h>
#include <AliHit.h>

//__________________AliRICHhit______________________________________________________________________
class AliRICHhit : public AliHit
{
public:
           AliRICHhit():AliHit(),fChamber(kBad),fEloss(kBad) {fInX3.SetXYZ(0,0,0);fOutX3.SetXYZ(0,0,0);}
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
  AliRICHdigit():AliDigit(),fCFM(0),fChamber(0),fPadX(0),fPadY(0),fQdc(kBad){fTracks[0]=fTracks[1]=fTracks[2]=kBad;}
  AliRICHdigit(Int_t c,TVector pad,Double_t q,Int_t cfm,Int_t tid0,Int_t tid1,Int_t tid2)  
       {fPadX=(Int_t)pad[0];fPadY=(Int_t)pad[1];fQdc=q;fChamber=10*c+AliRICHParam::Pad2Sec(pad);fCFM=cfm;fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}
  virtual ~AliRICHdigit() {;}  
  Int_t    Compare(const TObject *pObj) const
                 {if(Id()==((AliRICHdigit*)pObj)->Id())return 0;else if(Id()>((AliRICHdigit*)pObj)->Id())return 1;else return -1;}  //virtual      
  virtual Bool_t   IsSortable()                 const{return kTRUE;}                              //sort interface
  virtual void     Print(Option_t *option="")   const;                                            //virtual
  Int_t    ChFbMi()                     const{return fCFM;}                               //particle mixture for this digit
  Int_t    C()                          const{return fChamber/10;}                        //chamber number
  Int_t    S()                          const{return fChamber-(fChamber/10)*10;}          //sector number
  Int_t    X()                          const{return fPadX;}                              //x position of the pad
  Int_t    Y()                          const{return fPadY;}                              //y postion of the pad
  TVector  Pad()                        const{Float_t v[2]={fPadX,fPadY}; return TVector(2,v);}
  Int_t    Id()                         const{return fChamber*10000000+fPadX*1000+fPadY;} //absolute id of this pad
  Double_t Q()                          const{return fQdc;}                               //charge in terms of ADC channels
  void     AddTidOffset(Int_t offset) 
    {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;};
protected:
  Int_t    fCFM;  //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
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
  enum ClusterStatus {kEdge,kShape,kSize,kRaw,kResolved,kEmpty=kBad};
                    AliRICHcluster():TObject(),fCFM(0),fSize(0),fShape(0),fQdc(0),fChamber(0),fX(0),fY(0),fStatus(kEmpty),fDigits(0) {}    
  virtual          ~AliRICHcluster()                 {AliDebug(1,"Start");/*Reset();*/}  
         void       Reset()                          {DeleteDigits();fCFM=fSize=fShape=fQdc=fChamber=0;fX=fY=0;fStatus=kEmpty;} //cleans the cluster
         void       DeleteDigits()                   {if(fDigits) {delete fDigits;} fDigits=0;}           //deletes the list of digits  
  AliRICHcluster&   operator=(const AliRICHcluster&) {return *this;}                                 
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
  inline  void      AddDigit(AliRICHdigit *pDig);                                                      //
  inline  void      CoG(Int_t nLocals);                                                                //calculates center of gravity
          void      Fill(AliRICHcluster *pRaw,Double_t x,Double_t y,Double_t q,Int_t cfm)              //form new resolved cluster from raw one
                    {fCFM=cfm;fChamber=pRaw->Fchamber();fSize=pRaw->Fsize();fQdc=(Int_t)(q*pRaw->Q());fX=x;fY=y;fStatus=kResolved;} 
         Double_t   DistTo(TVector2 x)          const{return TMath::Sqrt((x.X()-fX)*(x.X()-fX)+(x.Y()-fY)*(x.Y()-fY));} //distance to given point 
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
  ClassDef(AliRICHcluster,2)  //RICH cluster class       
};//class AliRICHcluster
//__________________________________________________________________________________________________
void AliRICHcluster::AddDigit(AliRICHdigit *pDig)
{
// Adds a given digit to the list of digits belonging to this cluster    
  if(!fDigits) {fQdc=fSize=fCFM=0;fDigits = new TObjArray;}
  fQdc+=(Int_t)pDig->Q(); fDigits->Add(pDig);
  fChamber=10*pDig->C()+pDig->S();
  fSize+=10000;
  fStatus=kRaw;
}
//__________________________________________________________________________________________________
void AliRICHcluster::CoG(Int_t nLocals)
{
// Calculates naive cluster position as a center of gravity of its digits.
  Float_t xmin=999,ymin=999,xmax=0,ymax=0;   
  fX=fY=0;
  for(Int_t iDig=0;iDig<Size();iDig++) {
    AliRICHdigit *pDig=(AliRICHdigit*)fDigits->At(iDig);
    TVector pad=pDig->Pad(); Double_t q=pDig->Q();
    TVector2 x2=AliRICHParam::Pad2Loc(pad);
    fX += x2.X()*q;fY +=x2.Y()*q;
    if(pad[0]<xmin)xmin=pad[0];if(pad[0]>xmax)xmax=pad[0];if(pad[1]<ymin)ymin=pad[1];if(pad[1]>ymax)ymax=pad[1];
   }
   fX/=fQdc;fY/=fQdc;//Center of Gravity
   fShape=Int_t(100*(xmax-xmin+1)+ymax-ymin+1);//find box containing cluster
   fSize+=nLocals;
   fStatus=kRaw;
}//CoG()

//__________________AliRICH_________________________________________________________________________
class AliESD;

class AliRICH : public AliDetector 
{
public:
            AliRICH();                                            
            AliRICH(const char *name, const char *title);
            AliRICH(const AliRICH& RICH):AliDetector(RICH) {;}  //copy ctor 
  virtual  ~AliRICH();                                            
          
  AliRICH&  operator=(const AliRICH&)                 {return *this;}
//framework part  
  virtual Int_t         IsVersion()                           const =0;                                  //interface from         
  virtual void          StepManager()                               =0;                                  //interface from AliMC
  virtual void          Hits2SDigits();                                                                  //interface from AliSimulation
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* man) const {return new AliRICHDigitizer(man);}  //interface from AliSimulation
//  virtual void          Reconstruct()                         const;                                     //interface from AliReconstruction
//  virtual void          FillESD(AliESD *pESD)                 const;                                     //interface from AliReconstruction          
  virtual void          SetTreeAddress();                                                                //interface from AliLoader
  virtual void          MakeBranch(Option_t *opt=" ");                                                   //interface from AliLoader
  virtual void          CreateMaterials();                                                               //interface from AliMC
  virtual void          CreateGeometry();                                                                //interface from AliMC
  virtual void          BuildGeometry();                                                                 //interface 
//private part  
          Float_t AbsoCH4(Float_t x)const;                               //calculates absorption length for methane
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)const;  //deals with Fresnel absorption
  inline  void    CreateHits();                                          //create hits container as a simple list
  inline  void    CreateSDigits();                                       //create sdigits container as a simple list
  inline  void    CreateDigits();                                        //create digits container as 7 lists, one per chamber
  inline  void    CreateClusters();                                      //create clusters container  as 7 lists, one per chamber
//        void    ResetHits()                {AliDetector::ResetHits();}  //virtual  
          void    ResetSDigits()             {fNsdigits=0;  if(fSdigits)  fSdigits ->Clear();}                                 
          void    ResetDigits()              {if(fDigitsNew)for(int i=0;i<kNchambers;i++){fDigitsNew->At(i)->Clear();fNdigitsNew[i]=0;}} //virtual
          void    ResetClusters()            {if(fClusters) for(int i=0;i<kNchambers;i++){fClusters ->At(i)->Clear();fNclusters[i]=0;}}
  TClonesArray*   SDigits()             const{return fSdigits;}
  TClonesArray*   Digits(Int_t iC)      const{if(fDigitsNew) return (TClonesArray *)fDigitsNew->At(iC-1);else return 0;}
  TClonesArray*   Clusters(Int_t iC)    const{if(fClusters)  return (TClonesArray *)fClusters->At(iC-1);else return 0;}
  AliRICHChamber* C(Int_t iC)           const{return fpParam->C(iC);}                       //provides pointer to a given chamber
  AliRICHParam*   P()                   const{return fpParam;}                              //provides pointer to a RICH params
  AliRICH*        R()                        {return this;}                                 //provides pointer to RICH main object
  TVector         Counters()            const{return fCounters;}                            //provides a set of counters
  void            ControlPlots();                                                           //utility
  virtual void    Print(Option_t *option="")               const;                           //prints current RICH status
  void            PrintHits    (Int_t iEvent=0);                                            //utility
  void            PrintSDigits (Int_t iEvent=0);                                            //utility
  void            PrintDigits  (Int_t iEvent=0);                                            //utility
  void            PrintClusters(Int_t iEvent=0);                                            //utility
  void            PrintTracks  (Int_t iEvent=0);                                            //utility
            
  void AddHit(Int_t c,Int_t tid,TVector3 i3,TVector3 o3,Double_t eloss=0){TClonesArray &tmp=*fHits;new(tmp[fNhits++])AliRICHhit(c,tid,i3,o3,eloss);}
  inline void AddSDigit(Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid); 
  void AddDigit(int c,TVector pad,int q,int cfm,int *tid)//Add simulated digit
       {TClonesArray &tmp=*((TClonesArray*)fDigitsNew->At(c-1));new(tmp[fNdigitsNew[c-1]++])AliRICHdigit(c,pad,q,cfm,tid[0],tid[1],tid[2]);}  
  void AddDigit(Int_t c,TVector pad,Int_t q)//for real data digits
       {TClonesArray &tmp=*((TClonesArray*)fDigitsNew->At(0));new(tmp[fNdigitsNew[0]++])AliRICHdigit(c,pad,q,0,-1,-1,-1);}  
  void AddCluster(AliRICHcluster &cl)                     
       {Int_t c=cl.C()-1;TClonesArray &tmp=*((TClonesArray*)fClusters->At(c));new(tmp[fNclusters[c]++])AliRICHcluster(cl);}
  AliRICHhit* Hit(Int_t tid);           //returns pointer ot RICH hit for a given tid
protected:  
  enum                  {kAir=1,kRoha,kSiO2,kC6F14,kCH4,kCsI,kGridCu,kOpSiO2,kGap,kAl,kGlass,kCu,kW,kSteel,kPerpex,kSr90};
  AliRICHParam         *fpParam;             //main RICH parametrization     
                                             //fHits and fDigits belong to AliDetector
  TClonesArray         *fSdigits;            //! list of sdigits  
  Int_t                 fNsdigits;           //! current number of sdigits
  
  TObjArray            *fDigitsNew;          //! each chamber holds it's one lists of digits
  Int_t                 fNdigitsNew[kNchambers];   //! array of current numbers of digits
  
  TObjArray            *fClusters;           //! each chamber holds it's one lists of clusters 
  Int_t                 fNclusters[kNchambers];    //! array of current numbers of raw clusters
  
  TVector               fCounters;           //Photon history conters, explanation in StepManager() 
  ClassDef(AliRICH,7)                        //Main RICH class 
};//class AliRICH  

//__________________________________________________________________________________________________
void AliRICH::CreateHits()
{
  if(fHits) return;
  AliDebug(1,"creating hits container.");
  fHits=new TClonesArray("AliRICHhit",10000);   fNhits=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateSDigits()
{
  if(fSdigits) return;
  AliDebug(1,"creating sdigits container.");
  fSdigits=new TClonesArray("AliRICHdigit",10000); fNsdigits=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateDigits()
{
  if(fDigitsNew) return;
  AliDebug(1,"creating digits containers.");
  fDigitsNew = new TObjArray(kNchambers);  
  for(Int_t i=0;i<kNchambers;i++) {fDigitsNew->AddAt(new TClonesArray("AliRICHdigit",10000), i); fNdigitsNew[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::CreateClusters()
{
  if(fClusters) return;
  AliDebug(1,"creating clusters containers.");
  fClusters = new TObjArray(kNchambers);  
  for(Int_t i=0;i<kNchambers;i++) {fClusters->AddAt(new TClonesArray("AliRICHcluster",10000), i); fNclusters[i]=0;}
}
//__________________________________________________________________________________________________
void AliRICH::AddSDigit(Int_t c,TVector pad,Double_t q,Int_t pid,Int_t tid) 
{ 
  Int_t cfm;  
  switch(pid){
    case 50000050: cfm=1000000;break;//cerenkov
    case 50000051: cfm=1000;   break;//feedback
    default:       cfm=1;      break;//mip
  }   
  TClonesArray &tmp=*fSdigits; new(tmp[fNsdigits++])AliRICHdigit(c,pad,q,cfm,tid,kBad,kBad);
}
//__________________________________________________________________________________________________
#endif//#ifndef AliRICH_h
