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

class AliRICHRawCluster;
class AliRICHRecHit1D;
class AliRICHRecHit3D;

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
    
  Int_t   Chamber()             {return fChamber;}
  Int_t   Particle()            {return fParticle;}    
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
  Float_t CerenkovAngle()       {return fMomX;}
  Float_t MomFreoX()            {return fMomX;}
  Float_t MomFreoY()            {return fMomY;}
  Float_t MomFreoZ()            {return fMomZ;}
protected:
  Int_t     fChamber;                      //chamber number
  Int_t     fParticle;                     //particle code
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
  fChamber=fParticle=kBad;
  fTheta=fPhi=fTlength=fEloss=fPHfirst=fPHlast=fLoss=kBad;
  fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=fMomFreoX=fMomFreoY=fMomFreoZ=kBad;
}//AliRICHhit::default ctor
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hit):
            AliHit(shunt, track)
{//ctor
  fChamber=vol[0];
  fParticle=(Int_t)hit[0];
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
}//AliRICHhit::ctor
//__________________________________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t track,Int_t iPID,Int_t iChamber,TLorentzVector x4,Float_t eloss):
            AliHit(0, track)
{//ctor
  fChamber=iChamber;
  fParticle=iPID;
  fX=x4.X();fY=x4.Y();fZ=x4.Z();
  fEloss=eloss;
}//AliRICHhit::ctor

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
}//AliRICHCerenkov::ctor
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
}//AliRICHCerenkov::ctor

//__________________AliRICHdigit____________________________________________________________________
//__________________________________________________________________________________________________
//__________________________________________________________________________________________________
class AliRICHdigit :public AliDigit
{
public:
           AliRICHdigit() {fPadX=fPadY=fChamber=fAdc=fTracks[0]=fTracks[1]=fTracks[2]=kBad;}
  inline   AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT1,Int_t iT2,Int_t iT3);
  virtual ~AliRICHdigit() {;}  
  Int_t C()   const{return fChamber;}
  Int_t X()   const{return fPadX;}
  Int_t Y()   const{return fPadY;}
  Int_t Pad() const{return fPad;}
  Int_t Adc() const{return fAdc;}
protected:
  Int_t fChamber;  //module number 
  Int_t fPadX;     //pad number along X
  Int_t fPadY;     //pad number along Y
  Int_t fPad;      //pad number 1000*X+Y
  Int_t fAdc;      //ADC value
  ClassDef(AliRICHdigit,1) //RICH digit class       
};//class AliRICHdigit
//__________________________________________________________________________________________________
AliRICHdigit::AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{
  fChamber=iC;fPadX=iX;fPadY=iY;fPad=1000*fPadX+fPadY;fAdc=iAdc;
  fTracks[0]=iT0;fTracks[1]=iT1;fTracks[2]=iT2;
}//AliRICHdigit::ctor  

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
  virtual Int_t  IsVersion()const =0;            
  
  inline  void    AddHit(Int_t track, Int_t *vol, Float_t *hits);//virtual
  inline  void    AddHit(Int_t track,Int_t iPID,Int_t iChamber,TLorentzVector x4,Float_t eloss);
  inline  void    AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
  inline  void    AddSDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1=kBad,Int_t iT2=kBad);
  inline  void    AddDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1=kBad,Int_t iT2=kBad);
  
  inline  void    ResetHits();    //virtual
  inline  void    ResetSDigits(); 
          void    ResetDigits();  //virtual
  
  virtual void    Hits2SDigits();   
  virtual void    SDigits2Digits();
  virtual void    Digits2Reco();   
          
  AliRICHChamber* C(Int_t i)                       const{return (AliRICHChamber*)fChambers->At(i-1);}
  AliRICHParam*   Param()                          const{return fpParam;}
  TClonesArray*   SDigits()                        const{return fSDigits;}
  TClonesArray*   Cerenkovs()                      const{return fCerenkovs;}
  
          void    CreateChambers();         
  virtual void    CreateMaterials(); //GEANT materials definition
          Float_t AbsoCH4(Float_t x);
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
  virtual void    BuildGeometry();   //TNode ROOT variant for event display
  virtual void    CreateGeometry();  //GEANT volumes tree for simulation  
  
  virtual void    StepManager()=0;
          void    GenerateFeedbacks(Float_t eloss);
            
  AliRICHChamber& Chamber(Int_t id)      {return *((AliRICHChamber *) (*fChambers)[id]);}
  
          void  AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
          void  AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
          void  AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
          void  AddRecHit3D(Int_t id, Float_t* rechit, Float_t omega, Float_t theta, Float_t phi);
          void  ResetRawClusters();
          void  ResetRecHits1D();
          void  ResetRecHits3D();
  virtual void  FindClusters(Int_t nev);
          Int_t DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)      {return 9999;}
   
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
  virtual void   SetTreeAddress();
              
  
  TObjArray     *Dchambers()                     {return fDchambers;}
  TObjArray     *RecHits3D()                const{return fRecHits3D;}
  TObjArray     *RecHits1D()                const{return fRecHits1D;}
  Int_t         *Ndch()                          {return fNdch;}
  Int_t         *Nrechits1D()                    {return fNrechits1D;} 
  Int_t         *Nrechits3D()                    {return fNrechits3D;} 
  TClonesArray  *DigitsAddress(Int_t id)         {return ((TClonesArray *) (*fDchambers)[id]);}
  TClonesArray  *RecHitsAddress1D(Int_t id) const{return ((TClonesArray *) (*fRecHits1D)[id]);}
  TClonesArray  *RecHitsAddress3D(Int_t id) const{return ((TClonesArray *) (*fRecHits3D)[id]);}
  TClonesArray  *RawClustAddress(Int_t id)  const{return ((TClonesArray *) (*fRawClusters)[id]);}    
  AliRICHSDigit* FirstPad(AliRICHhit *hit, TClonesArray *clusters);
  AliRICHSDigit* NextPad(TClonesArray *clusters);

 
  virtual void Print(Option_t *option)const;
    
protected:
  AliRICHParam         *fpParam;             //main RICH parametrization     
  TObjArray            *fChambers;           //list of RICH chambers
  Int_t                 fNsdigits;           //Current number of sdigits
  Int_t                 fNcerenkovs;         //Current number of cerenkovs
  TClonesArray         *fSDigits;            //! List of sdigits
  TObjArray            *fDchambers;          //! Array of lists of digits
  TClonesArray         *fCerenkovs;          //! List of cerenkovs
  Int_t                 fNdch[kNCH];         //Array of current numbers of digits
  TObjArray            *fRawClusters;        //!List of raw clusters
  TObjArray            *fRecHits1D;          //!List of rec. hits
  TObjArray            *fRecHits3D;          //!List of rec. hits
  Int_t                 fNrawch[kNCH];       //Array of current numbers of raw clusters
  Int_t                 fNrechits1D[kNCH];   //Array of current numbers of rec hits 1D
  Int_t                 fNrechits3D[kNCH];   //Array of current numbers of rec hits 3D 

  Int_t fCkovNumber;                         // Number of Cerenkov photons
  Int_t fFreonProd;                          // Cerenkovs produced in freon
  Int_t fFeedbacks;                          // Number of feedback photons
    
  ClassDef(AliRICH,2)                        //Main RICH class 
};//class AliRICH
  
//__________________________________________________________________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{//Adds the current hit to the RICH hits list
  TClonesArray &tmp=*fHits;
  new(tmp[fNhits++])AliRICHhit(fIshunt,track,vol,hits);
}
//__________________________________________________________________________________________________
void AliRICH::AddHit(Int_t track,Int_t iPID,Int_t iChamber,TLorentzVector x4,Float_t eloss)
{//Adds the current hit to the RICH hits list
  TClonesArray &tmp=*fHits;
  new(tmp[fNhits++])AliRICHhit(track,iPID,iChamber,x4,eloss);
}
//__________________________________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{//Adds the current RICH cerenkov hit to the Cerenkovs list   
  TClonesArray &tmp=*fCerenkovs;
  new(tmp[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
}
//__________________________________________________________________________________________________
void AliRICH::ResetHits()
{//Resets hits and cerenkovs
  AliDetector::ResetHits();
  fNcerenkovs=0;
  if(fCerenkovs)fCerenkovs->Clear();
}
//__________________________________________________________________________________________________
void AliRICH::AddSDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{//Adds the current Sdigit to the RICH list of Sdigits   
  TClonesArray &tmp=*fSDigits;
  new(tmp[fNsdigits++])AliRICHdigit(iC,iX,iY,iAdc,iT0,iT1,iT2);
} 
//__________________________________________________________________________________________________
void AliRICH::ResetSDigits()
{//Resets sdigits
  fNsdigits=0;
  if(fSDigits)fSDigits->Clear();
}
//__________________________________________________________________________________________________
void AliRICH::AddDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{//Adds the current Digit to the RICH list of digits

   TClonesArray &tmp=*fDigits;
   new(tmp[fNdigits++]) AliRICHdigit(iC,iX,iY,iAdc,iT0,iT1,iT2);
}
#endif
