#ifndef AliHMPIDv3_h
#define AliHMPIDv3_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//.
//HMPID class for new geometry based on TGeo
//.
//.

#include "AliHMPID.h"             //base class 
#include "AliHMPIDDigitizer.h"    //CreateDigitizer()
#include <TGeoManager.h>   

class TGeoVolume;
class TGeoMedium;
  
class AliHMPIDv3 : public AliHMPID //TObject-TNamed-AliModule-AliDetector-AliHMPID-AliHMPIDv3
{
public:
                 AliHMPIDv3()                                   :AliHMPID(          ) {;}          //default ctor
                 AliHMPIDv3(const char *name, const char *title):AliHMPID(name,title) {;}          //named ctor
  virtual       ~AliHMPIDv3()                                                         {;}          //dtor
//framework part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
          void    AddAlignableVolumes(                               )const;                                   //from AliModule invoked from AliMC           
          void    CreateMaterials  (                                 );                                        //from AliModule invoked from AliMC
static    void    IdealPosition(Int_t iCh,TGeoHMatrix *m);                                                     //ideal position of a given chamber 
static    void    IdealPositionCradle(Int_t iCh,TGeoHMatrix *m);                                               //ideal position of a module of the cradle 
          void    CreateGeometry   (                                 );                                        //from AliModule invoked from AliMC
  AliDigitizer*   CreateDigitizer  (AliDigitizationInput *m               )const{return new AliHMPIDDigitizer(m);}  //from AliModule invoked from AliSimulation::RunDigitization()
          void    Digits2Raw       (                                 );                                        //from AliModule invoked from AliSimulation::WriteRawFiles()
  virtual void    DefineOpticalProperties(                           );                                        //from AliModule invoked from AliMC::ConstructOpGeometry() to set Cerenkov properties
         // void    InitProperties   (                                 );                                        //define the phys processes on/off (dray,eloss...)                                                 
          void    Hits2SDigits     (                                 );                                        //from AliModule invoked from AliSimulation::RunSDigitization()
          void    Init             (                                 );                                        //from AliModule invoked from AliMC::InitGeometry()
          Int_t   IsVersion        (                                 )const{return 1;                      }   //from AliModule not used
          void    Print            (const Option_t *opt=""           )const;                                   //from TObject
          Bool_t  Raw2SDigits      (AliRawReader *pRR                );                                        //from AliModule invoked from AliSimulation  
          void    StepManager      (                                 );                                        //from AliModule invoked from AliMC::Stepping()
//private part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
          void    GenFee           (Float_t qtot                     );                                        //generates feedback photons
  static  Float_t Fresnel          (Float_t geV,Float_t p, Bool_t pl );                                        //deals with Fresnel absorption on PC          
  static  void    Hit2Sdi          (TClonesArray *pH,TClonesArray *pS);                                        //hits to sdigits conversion
          Bool_t  IsLostByFresnel  (                                 );                                        //checks if the photon lost on Fresnel reflection  
          void    StepCount        (                                 );                                        //counts particles in StepManager()
          void    StepHistory      (                                 );                                        //prints history of tracking in StepManager()
  static  void    TestGeom         (                                 );                                        //tests the validity of geometry
  static  void    TestPoint        (Int_t ch,Float_t x,Float_t y     );                                        //tests the validity of geometry
  TGeoVolume*     CreateChamber    (Int_t number                     );                                        //chamber geometry building
  TGeoVolume*     CradleBaseVolume ( TGeoMedium *med,  Double_t l[7], const char *name);                       //base topology of the cradle elements
  TGeoVolume*     CreateCradle     (                                 );                                        //cradle geometry building
protected:
  enum EMedia {kAir=1,kRoha=2,kSiO2=3,kC6F14=4,kCH4=5,kCsI=6,kAl=7,kCu=8,kW=9,kNeo=10,kAr=11};                       //media ids, used for geometry creation  
  enum Ecounters {kMipEnterRad=1,kCkovNew,kCkovNewRad,kCkovNewWin,kCkovNewProxGap,kCkovNewAmpGap,kCkovEnterPc,kPhotoEle};    //counters id's
  ClassDef(AliHMPIDv3,2)                                                                //HMPID full version for simulation
};

#endif
