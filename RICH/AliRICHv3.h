#ifndef AliRICHv3_h
#define AliRICHv3_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliRICH.h"

class AliRICHv3 : public AliRICH 
{
enum EDebugBits {kDebugStart=BIT(0),kDebugParam=BIT(1),kDebugHit=BIT(2),kDebugDigit=BIT(3),kDebugReco=BIT(4)};
    
public:
    
   AliRICHv3():AliRICH()                                {} // Default ctor
   AliRICHv3(const char *pcName, const char *pcTitle);     // Named ctor 
   virtual       ~AliRICHv3()                           {} // Dtor
   
   virtual Int_t  IsVersion()                     const {return 3;}
   
   void   SetDebugStart()     {fDebugLevel+=kDebugStart;}        // Controls debug message at the entring point of methods
   void ResetDebugStart()     {fDebugLevel-=kDebugStart;}        // Controls debug message at the entring point of methods
   Bool_t  IsDebugStart()const{return fDebugLevel&kDebugStart;}  // Controls debug message at the entring point of methods
   
   void   SetDebugParam()     {fDebugLevel+=kDebugParam;}        // Controls debug printout for the parameters
   void ResetDebugParam()     {fDebugLevel-=kDebugParam;}        // Controls debug printout for the parameters
   Bool_t  IsDebugParam()const{return fDebugLevel&kDebugParam;}  // Controls debug printout for the parameters
   
   void   SetDebugHit()       {fDebugLevel+=kDebugHit;}          // Controls debug printout for hits
   void ResetDebugHit()       {fDebugLevel-=kDebugHit;}          // Controls debug printout for hits
   Bool_t  IsDebugHit()  const{return fDebugLevel&kDebugHit;}    // Controls debug printout for hits
   
   void   SetDebugDigit()     {fDebugLevel+=kDebugDigit;}        // Controls debug printout for digits
   void ResetDebugDigit()     {fDebugLevel-=kDebugDigit;}        // Controls debug printout for digits
   Bool_t  IsDebugDigit()const{return fDebugLevel&kDebugDigit;}  // Controls debug printout for digits
   
   void   SetDebugReco()      {fDebugLevel+=kDebugReco;}         // Controls debug printout for reco
   void ResetDebugReco()      {fDebugLevel-=kDebugReco;}         // Controls debug printout for reco
   Bool_t  IsDebugReco() const{return fDebugLevel&kDebugReco;}   // Controls debug printout for reco
   
   virtual void   CreateMaterials(); // Provides material definition for simulation (currently GEANT)   
   virtual void   CreateGeometry();  // Provides geometry structure for simulation (currently GEANT modules tree)
   virtual void   BuildGeometry();   // Provides geometry structure for event display (ROOT TNode tree)
   virtual void   Init();            // Makes nothing for a while 
   
private:
    ClassDef(AliRICHv3,1)  //RICH full version, configurable with azimuthal rotation	
};// class AliRICHv3
		
#endif // AliRICHv3_h
