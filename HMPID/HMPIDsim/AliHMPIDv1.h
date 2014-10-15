#ifndef AliHMPIDv1_h
#define AliHMPIDv1_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//.
//HMPID class of "standard" geometry 
//.
#include "AliHMPID.h"             //base class 
#include "AliHMPIDDigitizer.h"    //CreateDigitizer()

class TGeoManager;

class AliHMPIDv1 : public AliHMPID //TObject-TNamed-AliModule-AliDetector-AliHMPID-AliHMPIDv0
{
public:
                 AliHMPIDv1()                                   :AliHMPID(          ),fIdRad(-1),fIdWin(-1),fIdProxGap(-1),fIdAmpGap(-1),fIdPc(-1),fIdAnod(-1),fIdCath(-1),fIdColl(-1) {;}          //default ctor
                 AliHMPIDv1(const char *name, const char *title):AliHMPID(name,title),fIdRad(-1),fIdWin(-1),fIdProxGap(-1),fIdAmpGap(-1),fIdPc(-1),fIdAnod(-1),fIdCath(-1),fIdColl(-1) {;}          //named ctor
  virtual       ~AliHMPIDv1()                                                         {;}          //dtor
//framework part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
          void    AddAlignableVolumes(                               )const;                                   //from AliModule invoked from AliMC    
          void    CreateMaterials  (                                 );                                        //from AliModule invoked from AliMC
          void    CreateGeometry   (                                 );                                        //from AliModule invoked from AliMC
  AliDigitizer*   CreateDigitizer  (AliDigitizationInput *m               )const{return new AliHMPIDDigitizer(m);}  //from AliModule invoked from AliSimulation::RunDigitization()
          void    Digits2Raw       (                                 );                                        //from AliModule invoked from AliSimulation::WriteRawFiles()
  virtual void    DefineOpticalProperties();                                                                   //from AliModule invoked from AliMC::ConstructOpGeometry() to set Cerenkov properties
          void    Hits2SDigits     (                                 );                                        //from AliModule invoked from AliSimulation::RunSDigitization()
          void    Init             (                                 );                                        //from AliModule invoked from AliMC::InitGeometry()
          Int_t   IsVersion        (                                 )const{return 1;                      }   //from AliModule not used
          void    Print            (const Option_t *opt=""           )const;                                   //from TObject
          Bool_t  Raw2SDigits      (AliRawReader *pRR                );                                        //from AliMOdule invoked from AliSimulation  
          void    StepManager      (                                 );                                        //from AliModule invoked from AliMC::Stepping()
//private part++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
          void    GenFee           (Float_t qtot                     );                                        //generates feedback photons
  static  Float_t Fresnel          (Float_t geV,Float_t p, Bool_t pl );                                        //deals with Fresnel absorption on PC          
  static  void    Hit2Sdi          (TClonesArray *pH,TClonesArray *pS); 
          Bool_t  IsLostByFresnel  (                                 );                                        //checks if the photon lost on Fresnel reflection  
          void    StepCount        (                                 );                                        //counts particles in StepManager()
          void    StepHistory      (                                 );                                        //prints history of tracking in StepManager()
protected:
  enum EMedia {kAir=1,kRoha=2,kSiO2=3,kC6F14=4,kCH4=5,kCsI=6,kAl=7,kCu=8,kW=9};                               //media ids, used for geometry creation  
  enum Ecounters {kMipEnterRad=1,kCkovNew,kCkovNewRad,kCkovNewWin,kCkovNewProxGap,kCkovNewAmpGap,kCkovEnterPc,kPhotoEle};    //counters id's
  Int_t fIdRad,fIdWin,fIdProxGap,fIdAmpGap,fIdPc,fIdAnod,fIdCath,fIdColl;              //! volumes ID's used in StepManager() and Count()
  ClassDef(AliHMPIDv1,2)                                                                //HMPID full version for simulation
};

#endif
