#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <AliDetector.h>
#include <AliHit.h>
#include <AliDigit.h>

#include "AliRICHParam.h"

#include "AliRICHSDigit.h"

//__________________AliRICHhit______________________________________________________________________
class AliRICHhit : public AliHit
{
public:
  inline   AliRICHhit();
  inline   AliRICHhit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
  AliRICHhit(Int_t tid,TVector3 x3) :AliHit(0,tid) {fInX3=x3;}
  inline   AliRICHhit(Int_t tid,TVector3 x3in,TVector3 x3out,Double_t eloss);
  virtual ~AliRICHhit()         {;}

  Int_t   C()                   const{return fChamber;}
  Int_t   Chamber()             const{return fChamber;}
  Int_t   Pid()                 const{return fPid;}    
  Int_t   Particle()            const{return fPid;}
  Float_t Eloss()               const{return fEloss;}
  Float_t MomX()                const{return fMomX;}
  Float_t MomY()                const{return fMomY;}
  Float_t MomZ()                const{return fMomZ;}
  Float_t CerenkovAngle()       const{return fCerenkovAngle;}
  Float_t MomFreoX()            const{return fMomFreoX;}
  Float_t MomFreoY()            const{return fMomFreoY;}
  Float_t MomFreoZ()            const{return fMomFreoZ;}
  void    Print(Option_t *option="")const;      //virtual
protected:
  Int_t     fChamber;                      //chamber number
  Int_t     fPid;                          //particle code
  Double_t  fEloss;                        //ionisation energy loss in gas
  Float_t   fMomX,fMomY,fMomZ;             //momentum at photochatode entry point
  Float_t   fNPads;                        //Pads hit
  Float_t   fCerenkovAngle;                //Dummy cerenkov angle
  Float_t   fMomFreoX,fMomFreoY,fMomFreoZ; //momentum at freon entry point
  TVector3  fInX3,fOutX3;                  //3-vectors at the entrance and exit of the GAP
  ClassDef(AliRICHhit,2)                   //RICH hit class
};//class AliRICHhit

  //__________________________________________________________________________________________________
AliRICHhit::AliRICHhit()
           :AliHit() 
{//default ctor  
  fChamber=fPid=kBad;
  fEloss=kBad;
  fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=fMomFreoX=fMomFreoY=fMomFreoZ=kBad;
  fInX3.SetXYZ(0,0,0);fOutX3.SetXYZ(0,0,0);
}
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hit)
           :AliHit(shunt, track)
{//ctor
  fChamber=vol[0];
  fPid=(Int_t)hit[0];
  fX=hit[1];fY=hit[2];fZ=hit[3];
  fEloss=hit[7];
  fMomX=hit[14];fMomY=hit[15];fMomZ=hit[16];
  fCerenkovAngle=hit[18];
  fMomFreoX=hit[19];fMomFreoY=hit[20];fMomFreoZ=hit[21];
}
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t tid,TVector3 x3in,TVector3 x3out,Double_t eloss)
           :AliHit(0,tid)
{//ctor
  fX=x3out.X();fY=x3out.Y();fZ=x3out.Z();
  fInX3=x3in;
  fOutX3=x3out;
  fEloss=eloss;
}

//__________________AliRICHCerenkov_________________________________________________________________
class AliRICHCerenkov: public AliHit 
{
public:
  inline   AliRICHCerenkov();
  inline   AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
  virtual ~AliRICHCerenkov() {;}
protected:
  Int_t     fChamber;          //chamber number
  Float_t   fTheta,fPhi;       //incident theta phi angles in degrees      
  Float_t   fTlength;          //track length inside the chamber
  Float_t   fEloss;            //ionisation energy loss in gas
  Int_t     fPHfirst;          //first padhit
  Int_t     fPHlast;           //last padhit
  Int_t     fCMother;          //index of mother particle
  Float_t   fLoss;             //nature of particle loss
  Float_t   fIndex;            //index of photon
  Float_t   fProduction;       //point of production
  Float_t   fMomX,fMomY,fMomZ; //local Momentum
  Float_t   fNPads;           // Pads hit
  Float_t   fCerenkovAngle;   // Cerenkov Angle
    
  ClassDef(AliRICHCerenkov,1)  //RICH cerenkov class
};//class AliRICHCerenkov

//__________________________________________________________________________________________________
AliRICHCerenkov::AliRICHCerenkov()
{//ctor
    fChamber=kBad;
    fX=fY=fZ=fTheta=fPhi=fTlength=fEloss=kBad;
    fPHfirst=fPHlast=fCMother=kBad;
    fLoss=fIndex=fProduction=fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=kBad;
}
//__________________________________________________________________________________________________
AliRICHCerenkov::AliRICHCerenkov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
                :AliHit(shunt, track)
{//ctor
    fChamber=vol[0];
    fX=hits[1];fY=hits[2];fZ=hits[3];
    fTheta=hits[4];fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t)hits[8];fPHlast=(Int_t)hits[9];
    fCMother=Int_t(hits[10]);
    fIndex = hits[11];
    fProduction = hits[12];  
    fLoss=hits[13];
    fMomX=hits[14];fMomY=hits[15];fMomZ=hits[16];
    fNPads=hits[17];
    fCerenkovAngle=hits[18];
}

//__________________AliRICHdigit____________________________________________________________________
class AliRICHdigit :public AliDigit
{
public:
           AliRICHdigit() {fCombiPid=fChamber=fPadX=fPadY=fTracks[0]=fTracks[1]=fTracks[2]=kBad;fQdc=kBad;}
           AliRICHdigit(Int_t c,Int_t x,Int_t y,Double_t q,Int_t cpid,Int_t tid0,Int_t tid1,Int_t tid2)
           {fPadX=x;fPadY=y;fQdc=q;fChamber=10*c+AliRICHParam::Pad2Sec(x,y);fCombiPid=cpid;fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}        
  virtual ~AliRICHdigit() {;}  
  Int_t    Compare(const TObject *pObj) const                                         //virtual
           {if(Id()==((AliRICHdigit*)pObj)->Id()) return 0; else if(Id()>((AliRICHdigit*)pObj)->Id()) return 1;  else return -1;} 
           
           Bool_t   IsSortable()                                  const{return kTRUE;}//virtual
           Int_t    CombiPid()                                    const{return fCombiPid;}
           Int_t    C()                                           const{return fChamber/10;}
           Int_t    S()                                           const{return fChamber-(fChamber/10)*10;} 
           Int_t    Chamber()                                     const{return C();}
           Int_t    Sector()                                      const{return S();}
           Int_t    X()                                           const{return fPadX;}
           Int_t    Y()                                           const{return fPadY;}
           Int_t    Id()                                          const{return fChamber*10000000+fPadX*1000+fPadY;}
           Double_t Q()                                           const{return fQdc;}
           Int_t    Tid(Int_t i)                                  const{return fTracks[i];}
           void  Print(Option_t *option="")const;                 //virtual
protected:
  Int_t    fCombiPid; //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //10*chamber number+ sector number 
  Int_t    fPadX;     //pad number along X
  Int_t    fPadY;     //pad number along Y
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHdigit,2) //RICH digit class       
};//class AliRICHdigit

//__________________AliRICHcluster__________________________________________________________________
class AliRICHcluster :public TObject
{
public:
  enum ClusterStatus {kOK,kEdge,kShape,kSize,kRaw};
                    AliRICHcluster()                                 {fSize=fQdc=fStatus=fChamber=fDimXY=kBad;fX=fY=kBad;fDigits=0;}
  virtual          ~AliRICHcluster()                                 {delete fDigits;}  
  AliRICHcluster&   operator=(const AliRICHcluster&)                 {return *this;}
         Int_t      Size()                                      const{return fSize;}                       //
         Int_t      DimXY()                                     const{return fDimXY;}                      //
         Int_t      C()                                         const{return fChamber/10;}                 //
         Int_t      S()                                         const{return fChamber-(fChamber/10)*10;}   //
         Int_t      Chamber()                                   const{return C();}                         //
         Int_t      Sector()                                    const{return S();}                         //
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
         void       SetCombiPid(Int_t ckov,Int_t feeds,Int_t mips)   {fCombiPid=1000000*ckov+1000*feeds+mips;}            //
         TObjArray* Digits()                                    const{return fDigits;}                     //  
         void       Print(Option_t *option="")const;                                                       //virtual
  inline void       AddDigit(AliRICHdigit *pDig);                                                          //
  inline void       CoG();                                                                                 // 
         void       Reset() {fSize=fQdc=fStatus=fChamber=fDimXY=kBad;fX=fY=kBad;delete fDigits;fDigits=0;} //
protected:
  Int_t         fCombiPid;    //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t         fSize;        //how many digits belong to this cluster    
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
  fSize++;
}
//__________________________________________________________________________________________________
void AliRICHcluster::CoG()
{//
  Int_t xmin=999,ymin=999,xmax=0,ymax=0;   
  Double_t x,y;        
  fX=fY=0;
  for(Int_t iDig=0;iDig<Size();iDig++) {
    AliRICHdigit *pDig=(AliRICHdigit*)fDigits->At(iDig);
    Int_t padX = pDig->X();Int_t padY = pDig->Y();Double_t q=pDig->Q();
    AliRICHParam::Pad2Loc(padX,padY,x,y);
    fX += x*q;fY +=y*q;
    if(padX<xmin)xmin=padX;if(padX>xmax)xmax=padX;if(padY<ymin)ymin=padY;if(padY>ymax)ymax=padY;
   }
   fX/=fQdc;fY/=fQdc;//Center of Gravity
   fDimXY = 100*(xmax-xmin+1)+ymax-ymin+1;//find box containing cluster
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
class AliRICHSDigit;

class AliRICH : public AliDetector 
{
public:
            AliRICH();                                            
            AliRICH(const char *name, const char *title);
            AliRICH(const AliRICH& RICH):AliDetector(RICH) {;}   
  virtual  ~AliRICH();                                            
          
  AliRICH&  operator=(const AliRICH&)                 {return *this;}
  virtual Int_t   IsVersion()                                            const =0;            
          void    Hits2SDigits();                                                                                 //virtual
  AliDigitizer*   CreateDigitizer(AliRunDigitizer* manager);                                                      //virtual
          void    SDigits2Digits();                                                                               //virtual
  
  inline  void    CreateHits();    
  inline  void    CreateSDigits();  
  inline  void    CreateDigits();  
  inline  void    CreateClusters();  
  inline  void    CreateRecos();  
  void AddHit(Int_t track, Int_t *vol, Float_t *hits)                 {TClonesArray &tmp=*fHits; new(tmp[fNhits++])AliRICHhit(fIshunt,track,vol,hits);}//virtual
  void AddHit(Int_t tid,TVector3 x3)                                  {TClonesArray &tmp=*fHits;new(tmp[fNhits++])AliRICHhit(tid,x3);} 
  void AddHit(Int_t tid,TVector3 x3in,TVector3 x3out,Double_t eloss)  {TClonesArray &tmp=*fHits;new(tmp[fNhits++])AliRICHhit(tid,x3in,x3out,eloss);} 
  inline void AddSDigit(int c,int x,int y,int q,int pid,int tid); 
  void AddDigit(int c,int x,int y,int q,int cpid,int *tid){TClonesArray &tmp=*((TClonesArray*)fDigitsNew->At(c-1));new(tmp[fNdigitsNew[c-1]++])AliRICHdigit(c,x,y,q,cpid,tid[0],tid[1],tid[2]);}  
  void AddCluster(AliRICHcluster &cl)                     {TClonesArray &tmp=*((TClonesArray*)fClusters->At(cl.C()-1));new(tmp[fNclusters[cl.C()-1]++])AliRICHcluster(cl);}
  void AddReco(Int_t tid,Double_t thetaCherenkov,Int_t nPhotons) {TClonesArray &tmp=*(TClonesArray*)fRecos;new(tmp[fNrecos++])AliRICHreco(tid,thetaCherenkov,nPhotons);}  
  void ResetHits()     {AliDetector::ResetHits();fNcerenkovs=0;if(fCerenkovs)fCerenkovs->Clear();fNspecials=0;if(fSpecials)fSpecials->Clear();}  //virtual  
  void ResetSDigits()  {fNsdigits=0;  if(fSdigits)  fSdigits ->Clear();}                                 
  void ResetDigits()   {if(fDigitsNew)for(int i=0;i<kNCH;i++){fDigitsNew->At(i)->Clear();fNdigitsNew[i]=0;}}
  void ResetClusters() {if(fClusters) for(int i=0;i<kNCH;i++){fClusters ->At(i)->Clear();fNclusters[i]=0;}}
  void ResetRecos()    {if(fRecos) fRecos->Clear();fNrecos=0;}
                  //Hits provided by AliDetector
  TClonesArray*   SDigits()             const{return fSdigits;}
  TClonesArray*   Digits(Int_t iC)      const{if(fDigitsNew) return (TClonesArray *)fDigitsNew->At(iC-1);else return 0;}
  TClonesArray*   Clusters(Int_t iC)    const{if(fClusters)  return (TClonesArray *)fClusters->At(iC-1);else return 0;}
  TClonesArray*   Recos()               const{return fRecos;}
          
  AliRICHChamber* C(Int_t iC)           const{return (AliRICHChamber*)fChambers->At(iC-1);}
  AliRICHParam*   Param()               const{return fpParam;}
  AliRICHParam*   P()                   const{return fpParam;}
          void    CreateChambers();         
          void    CreateMaterials(); //virtual
  virtual void    BuildGeometry();   //virtual
  virtual void    CreateGeometry();  //virtual
          Float_t AbsoCH4(Float_t x)const;
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)const;
  
  virtual void    StepManager()=0;
          void    GenerateFeedbacks(Int_t iChamber,Float_t eloss);
          void    Print(Option_t *option)const;//virtual
          void    MakeBranch(Option_t *opt=" ");
          void    SetTreeAddress();//virtual
// OLD staff OLD staff    
  inline  void    AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
  inline  void    AddSpecialOld(Int_t *array);      
          
  inline  void    CreateCerenkovsOld();  
  inline  void    CreateSpecialsOld();   
          void    ResetSpecialsOld(){fNspecials=0; if(fSpecials) fSpecials->Clear();}   
  TClonesArray*   Specials()            const{return fSpecials;}
  TClonesArray*   Cerenkovs()           const{return fCerenkovs;}
  
  
            
  AliRICHChamber& Chamber(Int_t id)      {return *((AliRICHChamber *) (*fChambers)[id]);}  
//          Int_t DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)      {return 9999;}
    
protected:  
  enum       {kCSI=6,kGAP=9};
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

  TClonesArray         *fCerenkovs;          //! ??? List of cerenkovs
  Int_t                 fNcerenkovs;         //! ??? Current number of cerenkovs
  TClonesArray         *fSpecials;           //! ??? List of specials
  Int_t                 fNspecials;          //! ??? Current number of specials  
  Int_t fCkovNumber;                         // Number of Cerenkov photons
  Int_t fFreonProd;                          // Cerenkovs produced in freon
  Int_t fFeedbacks;                          // Number of feedback photons
    
  ClassDef(AliRICH,4)                        //Main RICH class 
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
void AliRICH::AddSDigit(int c,int x,int y,int q,int pid,int tid) 
{   
  switch(pid){
    case 50000048: pid=1000000;break;
    case 50000052: pid=1000;   break;
    default:       pid=1;      break;
  }   
  TClonesArray &tmp=*fSdigits;
  new(tmp[fNsdigits++])AliRICHdigit(c,x,y,q,pid,tid,kBad,kBad);
}//AddSDigit()   


//______OLD OLD OLD OLD_____________________________________________________________________________
void AliRICH::CreateCerenkovsOld()
{
  if(fCerenkovs) return;
  if(GetDebug())Info("CreateCerenkovs","creating cerenkovs container.");
  fCerenkovs=new TClonesArray("AliRICHCerenkov",10000);   fNcerenkovs=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateSpecialsOld()
{
  if(fSpecials) return;
  if(GetDebug())Info("CreateSpecialsOld","creating SDigits special container.");
  fSpecials=new TClonesArray("AliRICHSDigit",100000); fNspecials=0;
}
//__________________________________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{//Adds the current RICH cerenkov hit to the Cerenkovs list   
  TClonesArray &tmp=*fCerenkovs;
  new(tmp[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
}
//__________________________________________________________________________________________________
void AliRICH::AddSpecialOld(Int_t *aiSDigit)
{// Adds the current Sdigit to the RICH list of Specials
  TClonesArray &lSDigits = *fSpecials;
  new(lSDigits[fNspecials++]) AliRICHSDigit(aiSDigit);
}
#endif//#ifndef AliRICH_h
