#ifndef AliRICH_h
#define AliRICH_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TObjArray.h>
#include <TClonesArray.h>
#include <AliDetector.h>
#include <AliHit.h>
#include <AliDigit.h>
#include "AliRICHConst.h"
#include "AliRICHChamber.h"

static const int kNCH=7;

class TFile;

class AliRICHRawCluster;
class AliRICHRecHit1D;
class AliRICHRecHit3D;
class AliRICHClusterFinder;
class AliRICHDetect;
class AliRICHChamber;
class AliRICHCerenkov;
class AliSegmentation;
class AliRICHResponse;
class AliRICHGeometry;
class AliRICHMerger;



class AliRICHhit : public AliHit
{
public:
  inline   AliRICHhit();
  inline   AliRICHhit(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *hits);
  virtual ~AliRICHhit()         {;}
    
  Int_t   Chamber()             {return fChamber;}
  Float_t Particle()            {return fParticle;}    
  Float_t Theta()               {return fTheta;}
  Float_t Phi()                 {return fPhi;}
  Float_t Tlength()             {return fTlength;}
  Float_t Eloss()               {return fEloss;}
  Float_t Loss()                {return fLoss;}
  Float_t   PHfirst()           {return fPHfirst;}
  Float_t   PHlast()            {return fPHlast;}
  Float_t MomX()                {return fMomX;}
  Float_t MomY()                {return fMomY;}
  Float_t MomZ()                {return fMomZ;}
  Float_t CerenkovAngle()       {return fMomX;}
  Float_t MomFreoX()            {return fMomX;}
  Float_t MomFreoY()            {return fMomY;}
  Float_t MomFreoZ()            {return fMomZ;}
protected:
  Int_t     fChamber;       // Chamber number
  Float_t   fParticle;      // Geant3 particle type
  Float_t   fTheta ;        // Incident theta angle in degrees      
  Float_t   fPhi   ;        // Incident phi angle in degrees
  Float_t   fTlength;       // Track length inside the chamber
  Float_t   fEloss;         // ionisation energy loss in gas   
  Float_t   fPHfirst;       // first padhit
  Float_t   fPHlast;        // last padhit
  Float_t   fLoss;          // did it hit the freon?
  Float_t   fMomX;          // x Momentum at photochatode entry point
  Float_t   fMomY;          // y Momentum at photochatode entry point
  Float_t   fMomZ;          // z Momentum at photochatode entry point
  Float_t   fNPads;         // Pads hit
  Float_t   fCerenkovAngle; // Dummy cerenkov angle
  Float_t   fMomFreoX;      // x Momentum at freon entry point
  Float_t   fMomFreoY;      // y Momentum at freon entry point
  Float_t   fMomFreoZ;      // z Momentum at freon entry point		   
        
  ClassDef(AliRICHhit,1)  //RICH hit class
};//class AliRICHhit
//______________________________________________________________________________  
AliRICHhit::AliRICHhit()
           :AliHit() 
{//default ctor  
  fChamber=-1;
  fParticle=fTheta=fPhi=fTlength=fEloss=fPHfirst=fPHlast=fLoss=-1;
  fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=fMomFreoX=fMomFreoY=fMomFreoZ=-1;
}//AliRICHhit::default ctor
//______________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
            AliHit(shunt, track)
{//ctor
  fChamber=vol[0];
  fParticle=hits[0];
  fX=hits[1];
  fY=hits[2];
  fZ=hits[3];
  fTheta=hits[4];
  fPhi=hits[5];
  fTlength=hits[6];
  fEloss=hits[7];
  fPHfirst=(Int_t) hits[8];
  fPHlast=(Int_t) hits[9];
  fLoss=hits[13];
  fMomX=hits[14];
  fMomY=hits[15];
  fMomZ=hits[16];
  fNPads=hits[17];
  fCerenkovAngle=hits[18];
  fMomFreoX=hits[19];
  fMomFreoY=hits[20];
  fMomFreoZ=hits[21];
}//AliRICHhit::ctor
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
class AliRICHCerenkov: public AliHit 
{
public:
  inline   AliRICHCerenkov();
  inline   AliRICHCerenkov(Int_t fIshunt, Int_t track, Int_t *vol, Float_t *Cerenkovs);
  virtual ~AliRICHCerenkov() {;}
public:
  Int_t     fChamber;         // Chamber number
  Float_t   fTheta ;          // Incident theta angle in degrees      
  Float_t   fPhi   ;          // Incident phi angle in degrees
  Float_t   fTlength;         // Track length inside the chamber
  Float_t   fEloss;           // ionisation energy loss in gas
  Int_t     fPHfirst;         // first padhit
  Int_t     fPHlast;          // last padhit
  Int_t     fCMother;         // index of mother particle
  Float_t   fLoss;            // nature of particle loss
  Float_t   fIndex;           // Index of photon
  Float_t   fProduction;      // Point of production
  Float_t   fMomX;            // Local Momentum
  Float_t   fMomY;            // Local Momentum
  Float_t   fMomZ;            // Local Momentum
  Float_t   fNPads;           // Pads hit
  Float_t   fCerenkovAngle;   // Cerenkov Angle
    
  ClassDef(AliRICHCerenkov,1)  //RICH cerenkov class
};//class AliRICHCerenkov
//______________________________________________________________________________
AliRICHCerenkov::AliRICHCerenkov()
{//ctor
    fChamber=-1;
    fX=fY=fZ=fTheta=fPhi=fTlength=fEloss=-1;
    fPHfirst=fPHlast=fCMother=-1;
    fLoss=fIndex=fProduction=fMomX=fMomY=fMomZ=fNPads=fCerenkovAngle=-1;
}//AliRICHCerenkov::ctor
//______________________________________________________________________________
AliRICHCerenkov::AliRICHCerenkov(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits)
                :AliHit(shunt, track)
{//ctor
    fChamber=vol[0];
    fX=hits[1];
    fY=hits[2];
    fZ=hits[3];
    fTheta=hits[4];
    fPhi=hits[5];
    fTlength=hits[6];
    fEloss=hits[7];
    fPHfirst=(Int_t) hits[8];
    fPHlast=(Int_t) hits[9];
    fCMother=Int_t(hits[10]);
    fIndex = hits[11];
    fProduction = hits[12];  
    fLoss=hits[13];
    fMomX=hits[14];
    fMomY=hits[15];
    fMomZ=hits[16];
    fNPads=hits[17];
    fCerenkovAngle=hits[18];
}//AliRICHCerenkov::ctor
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
class AliRICHdigit :public AliDigit
{
public:
           AliRICHdigit() {fPadX=fPadY=fChamber=fAdc=fTracks[0]=fTracks[1]=fTracks[2]=-1;}
  inline   AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT1,Int_t iT2,Int_t iT3);
  virtual ~AliRICHdigit() {;}  
  Int_t C()   const{return fChamber;}
  Int_t X()   const{return fPadX;}
  Int_t Y()   const{return fPadY;}
  Int_t Adc() const{return fAdc;}
protected:
  Int_t fChamber;  //module number 
  Int_t fPadX;    //pad number along X
  Int_t fPadY;    //pad number along Y
  Int_t fAdc;     //ADC value
  ClassDef(AliRICHdigit,1) //RICH digit class       
};//class AliRICHdigit
//______________________________________________________________________________
AliRICHdigit::AliRICHdigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{
  fChamber=iC;fPadX=iX;fPadY=iY;fAdc=iAdc;
  fTracks[0]=iT0;fTracks[1]=iT1;fTracks[2]=iT2;
}//AliRICHdigit::ctor  
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
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
  AliRICHParam *Param()                               {return fpParam;}
  inline  void  AddHit(Int_t track, Int_t *vol, Float_t *hits);//virtual
  inline  void  AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs);
  inline  void  AddSDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1=-1,Int_t iT2=-1);
  inline  void  ResetHits();    //virtual
  inline  void  ResetSDigits(); //virtual
          void  Hits2SDigits(); //virtual 
          
  TClonesArray  *SDigits()                       const{return fSDigits;}
  TClonesArray  *Cerenkovs()                     const{return fCerenkovs;}
          void    CreateChambers();         
  virtual void    CreateMaterials(); //GEANT materials definition
          Float_t AbsoCH4(Float_t x);
          Float_t Fresnel(Float_t ene,Float_t pdoti, Bool_t pola);
  virtual void    BuildGeometry();   //TNode ROOT variant for event display
  virtual void    CreateGeometry();  //GEANT volumes tree for simulation  
  virtual void    StepManager()=0;
  
  AliRICHChamber* C(Int_t i)        const{return (AliRICHChamber*)fChambers->At(i);}//return pointer to requested chamber
  AliRICHChamber& Chamber(Int_t id)      {return *((AliRICHChamber *) (*fChambers)[id]);}
  TObjArray*      Chambers()        const{return fChambers;}
  
        void  AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits);
          void  AddRawCluster(Int_t id, const AliRICHRawCluster& cluster);
          void  AddRecHit1D(Int_t id, Float_t* rechit, Float_t* photons, Int_t* padsx, Int_t* padsy);
          void  AddRecHit3D(Int_t id, Float_t* rechit, Float_t omega, Float_t theta, Float_t phi);
          void  ResetDigits();  //virtual
          void  ResetRawClusters();
          void  ResetRecHits1D();
          void  ResetRecHits3D();
  virtual void  FindClusters(Int_t nev);
          Int_t Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss,Int_t id, ResponseType res);//kir ????? to be  removed
  virtual void  SDigits2Digits();
  virtual void  Digits2Reco();   
          Int_t    DistancetoPrimitive(Int_t /*px*/, Int_t /*py*/)      {return 9999;}
   
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   MakeBranchInTreeD(TTree *treeD, const char *file=0);
  virtual void   SetTreeAddress();
   
     
   
  AliRICHSDigit* FirstPad(AliRICHhit *hit, TClonesArray *clusters);
  AliRICHSDigit* NextPad(TClonesArray *clusters);
   

  void     SetGeometryModel(Int_t iC,AliRICHGeometry *pRICHGeo)                    {C(iC)->SetGeometryModel(pRICHGeo);}
  void     SetSegmentationModel(Int_t iC, AliSegmentation *pAliSeg)                {C(iC)->SetSegmentationModel(pAliSeg);}
  void     SetResponseModel(Int_t iC, AliRICHResponse *pRICHRes)                   {C(iC)->SetResponseModel(pRICHRes);}
  void     SetReconstructionModel(Int_t iC, AliRICHClusterFinder *pRICHReco)       {C(iC)->SetReconstructionModel(pRICHReco);}
  AliRICHGeometry* GetGeometryModel(Int_t iC=0)                               const{return C(iC)->GetGeometryModel();}    
  AliSegmentation* GetSegmentationModel(Int_t iC=0)                           const{return C(iC)->GetSegmentationModel();}
  AliRICHResponse* GetResponseModel(Int_t iC)                                 const{return C(iC)->GetResponseModel();}

//kir  virtual void   SetMerger(AliRICHMerger* thisMerger) {fMerger=thisMerger;}  
  
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

  void DiagnosticsFE(Int_t evNumber1=0,Int_t evNumber2=0);    // Full events
  void DiagnosticsSE(Int_t diaglevel,Int_t evNumber1=0,Int_t evNumber2=0);    // Single events
 
  virtual void Print(Option_t *option)const; // Prints debug information
    
protected:
  AliRICHParam         *fpParam;              //main RICH parametrization     
  TObjArray            *fChambers;           //! List of RICH chambers
  Int_t                 fNsdigits;           //Current number of sdigits
  Int_t                 fNcerenkovs;         //Current number of cerenkovs
  TClonesArray         *fSDigits;            //! List of sdigits
  TObjArray            *fDchambers;          //! Array of lists of digits
  TClonesArray         *fCerenkovs;          //! List of cerenkovs
  Int_t                 fNdch[kNCH];         //Array of current numbers of digits
  TObjArray            *fRawClusters;        // !List of raw clusters
  TObjArray            *fRecHits1D;          // !List of rec. hits
  TObjArray            *fRecHits3D;          // !List of rec. hits
  Int_t                 fNrawch[kNCH];       //Array of current numbers of raw clusters
  Int_t                 fNrechits1D[kNCH];   //Array of current numbers of rec hits 1D
  Int_t                 fNrechits3D[kNCH];   //Array of current numbers of rec hits 3D 

  Int_t fCkovNumber;                         // Number of Cerenkov photons
  Int_t fFreonProd;                          // Cerenkovs produced in freon
  Int_t fFeedbacks;                          // Number of feedback photons
//kir  Int_t fCkovQuarz;                          // Cerenkovs crossing quartz
//kir  Int_t fCkovGap;                            // Cerenkovs crossing gap
//kir  Int_t fCkovCsi;                            // Cerenkovs crossing csi
//kir  Int_t fLostRfreo;                          // Cerenkovs reflected in freon
//kir  Int_t fLostRquar;                          // Cerenkovs reflected in quartz
//kir  Int_t fLostAfreo;                          // Cerenkovs absorbed in freon 
//kir  Int_t fLostAquarz;                         // Cerenkovs absorbed in quartz
//kir  Int_t fLostAmeta;                          // Cerenkovs absorbed in methane
//kir  Int_t fLostCsi;                            // Cerenkovs below csi quantum efficiency 
//kir  Int_t fLostWires;                          // Cerenkovs lost in wires
//kir  Float_t fMipx;                             // x coord. of MIP
//kir  Float_t fMipy;                             // y coord. of MIP
//kir  Int_t fLostFresnel;                        // Cerenkovs lost by Fresnel reflection

//kir  Text_t *fFileName;                         //! File with background hits
//kir  AliRICHMerger *fMerger;                    //! pointer to merger
    
  ClassDef(AliRICH,2)                        //Main RICH class 
};//class AliRICH

//______________________________________________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{//Adds the current hit to the RICH hits list
  TClonesArray &tmp=*fHits;
  new(tmp[fNhits++])AliRICHhit(fIshunt,track,vol,hits);
}
//______________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{//Adds the current RICH cerenkov hit to the Cerenkovs list   
  TClonesArray &tmp=*fCerenkovs;
  new(tmp[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
}
//______________________________________________________________________________
void AliRICH::ResetHits()
{//Resets hits and cerenkovs
  AliDetector::ResetHits();
  fNcerenkovs = 0;
  if(fCerenkovs)fCerenkovs->Clear();
}
//______________________________________________________________________________
void AliRICH::AddSDigit(Int_t iC,Int_t iX,Int_t iY,Int_t iAdc,Int_t iT0,Int_t iT1,Int_t iT2)
{//Adds the current Sdigit to the RICH list of Sdigits   
  TClonesArray &tmp = *fSDigits;
  new(tmp[fNsdigits++])AliRICHdigit(iC,iX,iY,iAdc,iT0,iT1,iT2);
} 
//______________________________________________________________________________    
void AliRICH::ResetSDigits()
{//Resets sdigits
  fNsdigits=0;
  if(fSDigits)fSDigits->Clear();
}
//______________________________________________________________________________
#endif
