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

#include "AliRICHConst.h"
#include "AliRICHChamber.h"

#include "AliRICHDigit.h"
#include "AliRICHSDigit.h"
#include "AliRICHRawCluster.h"
class AliRICHRecHit1D;

//__________________AliRICHhit______________________________________________________________________
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
class AliRICHhit : public AliHit
{
public:
  inline   AliRICHhit();
  inline   AliRICHhit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
  inline   AliRICHhit(Int_t track,Int_t iPID,Int_t iChamber,TLorentzVector x4,Float_t eloss);
  virtual ~AliRICHhit()         {;}
    
  Int_t   C()                   {return fChamber;}
  Int_t   Chamber()             {return fChamber;}
  Int_t   Pid()                 {return fPid;}    
  Int_t   Particle()            {return fPid;}
  Float_t Theta()               {return fTheta;}
  Float_t Phi()                 {return fPhi;}
  Float_t Tlength()             {return fTlength;}
  Float_t Eloss()               {return fEloss;}
  Float_t Loss()                {return fLoss;}
  Float_t PHfirst()             {return fPHfirst;}
  Float_t PHlast()              {return fPHlast;}
  Float_t MomX()                {return fMomX;}
  Float_t MomY()                {return fMomY;}
  Float_t MomZ()                {return fMomZ;}
  Float_t CerenkovAngle()       {return fCerenkovAngle;}
  Float_t MomFreoX()            {return fMomFreoX;}
  Float_t MomFreoY()            {return fMomFreoY;}
  Float_t MomFreoZ()            {return fMomFreoZ;}
  void    Print(Option_t *option="")const;      //virtual
protected:
  Int_t     fChamber;                      //chamber number
  Int_t     fPid;                          //particle code
  Float_t   fTheta,fPhi ;                  //incident theta phi angles in degrees
  Float_t   fTlength;                      //track length inside the chamber
  Float_t   fEloss;                        //ionisation energy loss in gas
  Float_t   fPHfirst;                      //first padhit
  Float_t   fPHlast;                       //last padhit
  Float_t   fLoss;                         // did it hit the freon?
  Float_t   fMomX,fMomY,fMomZ;             //momentum at photochatode entry point
  Float_t   fNPads;                        // Pads hit
  Float_t   fCerenkovAngle;                // Dummy cerenkov angle
  Float_t   fMomFreoX,fMomFreoY,fMomFreoZ; //momentum at freon entry point
  ClassDef(AliRICHhit,1)                   //RICH hit class
};//class AliRICHhit

  //__________________________________________________________________________________________________
AliRICHhit::AliRICHhit()
           :AliHit() 
{//default ctor  
  fChamber=fPid=kBad;
  fTheta=fPhi=fTlength=fEloss=fPHfirst=fPHlast=fLoss=kBad;
  fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=fMomFreoX=fMomFreoY=fMomFreoZ=kBad;
}
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hit):
            AliHit(shunt, track)
{//ctor
  fChamber=vol[0];
  fPid=(Int_t)hit[0];
  fX=hit[1];fY=hit[2];fZ=hit[3];
  fTheta=hit[4];fPhi=hit[5];
  fTlength=hit[6];
  fEloss=hit[7];
  fPHfirst=(Int_t)hit[8];
  fPHlast=(Int_t)hit[9];
  fLoss=hit[13];
  fMomX=hit[14];fMomY=hit[15];fMomZ=hit[16];
  fNPads=hit[17];
  fCerenkovAngle=hit[18];
  fMomFreoX=hit[19];fMomFreoY=hit[20];fMomFreoZ=hit[21];
}
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t track,Int_t iPID,Int_t iChamber,TLorentzVector x4,Float_t eloss):
            AliHit(0, track)
{//ctor
  fChamber=iChamber;
  fPid=iPID;
  fX=x4.X();fY=x4.Y();fZ=x4.Z();
  fEloss=eloss;
}

//__________________AliRICHCerenkov_________________________________________________________________
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
class AliRICHCerenkov: public AliHit 
{
public:
  inline   AliRICHCerenkov();
  inline   AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
  virtual ~AliRICHCerenkov() {;}
public:
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
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
class AliRICHdigit :public AliDigit
{
public:
           AliRICHdigit() {fPadX=fPadY=fChamber=fQdc=fTracks[0]=fTracks[1]=fTracks[2]=kBad;}
  inline   AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iQdc,Int_t iT1,Int_t iT2,Int_t iT3);
  virtual ~AliRICHdigit() {;}  
  inline   Int_t Compare(const TObject *pObj) const;//virtual
           Bool_t IsSortable()                                 const{return kTRUE;}//virtual
           Int_t C()                                           const{return fChamber;}
           Int_t X()                                           const{return fPadX;}
           Int_t Y()                                           const{return fPadY;}
           Int_t Id()                                          const{return fChamber*1000000+fPadX*1000+fPadY;}
           Int_t Qdc()                                         const{return fQdc;}
           Int_t T(Int_t i)                                    const{return fTracks[i];}
           void  Print(Option_t *option)const;      //virtual
protected:
  Int_t fChamber;  //module number 
  Int_t fPadX;     //pad number along X
  Int_t fPadY;     //pad number along Y
  Int_t fQdc;      //ADC value
  ClassDef(AliRICHdigit,1) //RICH digit class       
};//class AliRICHdigit
//__________________________________________________________________________________________________
AliRICHdigit::AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iQdc,Int_t iT0,Int_t iT1,Int_t iT2)
{
  fChamber=iC;fPadX=iX;fPadY=iY;fQdc=iQdc;
  fTracks[0]=iT0;fTracks[1]=iT1;fTracks[2]=iT2;
}
//__________________________________________________________________________________________________
Int_t AliRICHdigit::Compare(const TObject *pObj)const
{
  if(Id()==((AliRICHdigit*)pObj)->Id()) 
    return 0;
  else if(Id()>((AliRICHdigit*)pObj)->Id()) 
    return 1;
  else
    return -1;
} 
//__________________AliRICH_________________________________________________________________________
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
class AliRICHParam;
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
          void    SDigits2Digits();                                                                               //virtual
          void    Digits2Reco();                                                                                  //virtual
  
  inline  void    CreateHits();    
  inline  void    CreateSdigits();  
  inline  void    CreateDigits();  
  inline  void    CreateClusters();  
  inline  void    AddHit(Int_t track, Int_t *vol, Float_t *hits);                                                 //virtual
  inline  void    AddSdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1=kBad,Int_t iT2=kBad);
  inline  void    AddDigit (Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1=kBad,Int_t iT2=kBad);     
          void    AddCluster()                                                                                {;}     
          void    ResetHits()     {AliDetector::ResetHits();fNcerenkovs=0;if(fCerenkovs)fCerenkovs->Clear();fNspecials=0;if(fSpecials)fSpecials->Clear();}  //virtual
          void    ResetSdigits()  {fNsdigits=0;  if(fSdigits)  fSdigits ->Clear();}                                 
          void    ResetDigits()   {if(fDigitsNew)for(int i=0;i<kNCH;i++){fDigitsNew->At(i)->Clear();fNdigitsNew[i]=0;}}
          void    ResetClusters() {if(fClusters) for(int i=0;i<kNCH;i++){fClusters ->At(i)->Clear();fNclusters[i]=0;}}
                  //Hits provided by AliDetector
  TClonesArray*   Sdigits()             const{return fSdigits;}
  TClonesArray*   Digits(Int_t iC)      const{if(fDigitsNew) return (TClonesArray *)fDigitsNew->At(iC-1);else return 0;}
  TClonesArray*   Clusters(Int_t iC)    const{if(fClusters)  return (TClonesArray *)fClusters->At(iC-1);else return 0;}
          
  AliRICHChamber* C(Int_t iC)           const{return (AliRICHChamber*)fChambers->At(iC-1);}
  AliRICHParam*   Param()               const{return fpParam;}
  
  //  AliRICHhit*     FirstHit(Int_t iTrkN)      {return (AliRICHhit*)AliDetector::FirstHit(iTrkN);}                   //virtual
  //  AliRICHhit*     NextHit()                  {return (AliRICHhit*)AliDetector::NextHit();}                         //virtual 
  
          void    CreateChambers();         
          void    CreateMaterials(); //virtual
  virtual void    BuildGeometry();   //virtual
  virtual void    CreateGeometry();  //virtual
          Float_t AbsoCH4(Float_t x);
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
  
  virtual void    StepManager()=0;
          void    GenerateFeedbacks(Int_t iChamber,Float_t eloss);
          void    Print(Option_t *option)const;//virtual
          void    MakeBranch(Option_t *opt=" ");
          void    SetTreeAddress();//virtual
// OLD staff OLD staff    
  inline  void    AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
  inline  void    AddSpecialOld(Int_t *);      
  inline  void    AddDigitOld(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);          
  inline  void    AddClusterOld(Int_t iChamber, const AliRICHRawCluster& cluster);
          void    AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
          
  inline  void    CreateCerenkovsOld();  
  inline  void    CreateSpecialsOld();   
  inline  void    CreateDigitsOld();    
  inline  void    CreateRawClustersOld();  
  inline  void    CreateRecos1Old();  
          void    ResetDigitsOld()  {if(fDchambers)  for(int i=0;i<kNCH;i++){fDchambers->At(i)->Clear();fNdch[i]=0;}}   //virtual
          void    ResetSpecialsOld(){fNspecials=0; if(fSpecials) fSpecials->Clear();}   
          void    ResetRawClusters(){if(fRawClusters)for(int i=0;i<kNCH;i++){fRawClusters->At(i)->Clear();fNrawch[i]=0;}}
          void    ResetRecHits1D()  {if(fRecHits1D)  for(int i=0;i<kNCH;i++){fRecHits1D  ->At(i)->Clear();fNrechits1D[i]=0;}}
  
  TClonesArray*   DigitsOld(Int_t iC)   const{if(fDchambers) return (TClonesArray *)fDchambers->At(iC-1);else return 0;}
  TClonesArray*   ClustersOld(Int_t iC) const{if(fRawClusters)return (TClonesArray *)fRawClusters->At(iC-1);else return 0;}
  TClonesArray*   Specials()            const{return fSpecials;}
  TClonesArray*   Cerenkovs()           const{return fCerenkovs;}
  
  
            
  AliRICHChamber& Chamber(Int_t id)      {return *((AliRICHChamber *) (*fChambers)[id]);}  
  TObjArray     *Dchambers()                const{return fDchambers;}
  TObjArray     *RecHits1D()                const{return fRecHits1D;}
  Int_t         *Ndch()                          {return fNdch;}
  Int_t         *Nrechits1D()                    {return fNrechits1D;} 
  TClonesArray  *DigitsAddress(Int_t id)    const{return ((TClonesArray *) (*fDchambers)[id]);}
  TClonesArray  *RecHitsAddress1D(Int_t id) const{return ((TClonesArray *) (*fRecHits1D)[id]);}
  TClonesArray  *RawClustAddress(Int_t id)  const{return ((TClonesArray *) (*fRawClusters)[id]);}    
//          Int_t DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)      {return 9999;}
    
protected:  
  AliRICHParam         *fpParam;             //main RICH parametrization     
  TObjArray            *fChambers;           //list of RICH chambers
                //fHits and fDigits belong to AliDetector
  TClonesArray         *fSdigits;            //! List of sdigits  
  Int_t                 fNsdigits;           //Current number of sdigits
  TObjArray            *fDigitsNew;          //! Each chamber holds it's one lists of digits
  Int_t                 fNdigitsNew[kNCH];   //Array of current numbers of digits
  TObjArray            *fClusters;           //! Each chamber holds it's one lists of clusters 
  Int_t                 fNclusters[kNCH];    //Array of current numbers of raw clusters

  TClonesArray         *fCerenkovs;          //! ??? List of cerenkovs
  Int_t                 fNcerenkovs;         //??? Current number of cerenkovs
  TClonesArray         *fSpecials;           //! ??? List of specials
  Int_t                 fNspecials;          //??? Current number of specials  
  TObjArray            *fDchambers;          //! Array of lists of digits  
  Int_t                 fNdch[kNCH];         //Array of current numbers of digits
  TObjArray            *fRawClusters;        //! Array of lists of raw clusters 
  Int_t                 fNrawch[kNCH];       //Array of current numbers of raw clusters
  TObjArray            *fRecHits1D;          //!List of rec. hits
  Int_t                 fNrechits1D[kNCH];   //Array of current numbers of rec hits 1D
  Int_t fCkovNumber;                         // Number of Cerenkov photons
  Int_t fFreonProd;                          // Cerenkovs produced in freon
  Int_t fFeedbacks;                          // Number of feedback photons
    
  ClassDef(AliRICH,3)                        //Main RICH class 
};//class AliRICH  
//__________________________________________________________________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{//Adds the current hit to the RICH hits list
  TClonesArray &tmp=*fHits;
  new(tmp[fNhits++])AliRICHhit(fIshunt,track,vol,hits);
}
//__________________________________________________________________________________________________
void AliRICH::AddSdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{//Adds the current Sdigit to the RICH list of Sdigits   
  TClonesArray &tmp=*fSdigits;
  new(tmp[fNsdigits++])AliRICHdigit(iC,iX,iY,iAdc,iT0,iT1,iT2);
} 
//__________________________________________________________________________________________________
void AliRICH::AddDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{//Adds the current digit to the corresponding RICH list of digits (individual list per chamber)
  TClonesArray &tmp=*((TClonesArray*)fDigitsNew->At(iC-1));
  new(tmp[fNdigitsNew[iC-1]++]) AliRICHdigit(iC,iX,iY,iAdc,iT0,iT1,iT2);
}
//__________________________________________________________________________________________________
void AliRICH::CreateHits()
{
  if(fHits) return;
  if(GetDebug())Info("CreateHits","creating hits container.");
  fHits=new TClonesArray("AliRICHhit",10000);   fNhits=0;
}
//__________________________________________________________________________________________________
void AliRICH::CreateSdigits()
{
  if(fSdigits) return;
  if(GetDebug())Info("CreateSdigits","creating sdigits container.");
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
void AliRICH::CreateDigitsOld()
{
  if(fDchambers) return;
  if(GetDebug())Info("CreateDigitsOld","creating digits containers.");
  fDchambers = new TObjArray(kNCH);  
  for(Int_t i=0;i<kNCH;i++) fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
}
//__________________________________________________________________________________________________
void AliRICH::CreateRawClustersOld()
{
  if(fRawClusters) return;
  if(GetDebug())Info("CreateClustersOld","creating clusters containers.");
  fRawClusters = new TObjArray(kNCH);
  for(Int_t i=0; i<kNCH ;i++) fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
}
//__________________________________________________________________________________________________
void AliRICH::CreateRecos1Old()
{
  if(fRecHits1D) return;
  if(GetDebug())Info("CreateRecos1DOld","creating recos 1 containers.");
  fRecHits1D = new TObjArray(kNCH);
  for(Int_t i=0; i<kNCH ;i++)  fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
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
//__________________________________________________________________________________________________
void AliRICH::AddDigitOld(Int_t iChamber, Int_t *tracks, Int_t *charges, Int_t *digits)
{// Adds the current digit to the RICH list of S digits
   TClonesArray &ldigits = *((TClonesArray*)fDchambers->At(iChamber-1));
   new(ldigits[fNdch[iChamber-1]++]) AliRICHDigit(tracks,charges,digits);
}
//__________________________________________________________________________________________________
void AliRICH::AddClusterOld(Int_t iChamber, const AliRICHRawCluster& c)
{// Add a RICH raw cluster to the list   
  TClonesArray &tmp= *((TClonesArray*)fRawClusters->At(iChamber-1));
  new(tmp[fNrawch[iChamber-1]++]) AliRICHRawCluster(c);
}
//__________________________________________________________________________________________________
#endif//#ifndef AliRICH_h
