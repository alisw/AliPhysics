#ifndef AliRICHreco_h
#define AliRICHreco_h


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */

#include <TObject.h>

class AliRICHreco : public TObject 
{
public:
// ctor & dtor:      
   inline AliRICHreco(); 
   inline AliRICHreco(Int_t iMipDigitOrClustN,Double_t dConeAngle,Double_t dPartOfRing,Double_t adPhotonAngles[],Int_t iNphotons);
   virtual ~AliRICHreco() {}
// inline methods:
       Double_t  GetConeAngle() const{return fdConeAngle;}    
          Int_t  GetNphotons()  const{return fNphotons;}
          Int_t  GetMipNumber() const{return fMipDigitOrClusterN;}
   static Int_t  GetArraySize()      {return 25;}    
protected:
   Int_t         fiArraySize;                //! Array size          
   Int_t         fMipDigitOrClusterN;        // Reference to MIP digit or MIP cluster in fDigits or fRawClusters       
   Int_t         fNphotons;                  // Actual number of photons  used in reconstruction
   Double_t      fdConeAngle;                // Cerenkov cone angle at the entrance of radiator for this track
   Double_t      fdPartOfRing;               // The fraction of the photons ring which does not touch the frame
   Double_t      fadPhotonAngles[25];        // Array of individual impact angles for fiNphotons photons
      
   ClassDef(AliRICHreco,1)                   // Reconstructed track information
};

inline AliRICHreco::AliRICHreco() 
{
   fMipDigitOrClusterN=fNphotons=0;
   fdConeAngle=fdPartOfRing=0;
   for(Int_t i=0;i<GetArraySize();i++) fadPhotonAngles[i]=0;
}//inline AliRICHreco::AliRICHreco() 

inline AliRICHreco::AliRICHreco(Int_t iMipDigitOrClusterN,Double_t dConeAngle,Double_t dPartOfRing, Double_t adPhotonAngles[], Int_t iNphotons)
{
   fMipDigitOrClusterN=iMipDigitOrClusterN;
   fdConeAngle=dConeAngle;
   fdPartOfRing=dPartOfRing;
   if(iNphotons<GetArraySize())
      fNphotons=iNphotons;
   else
      fNphotons=GetArraySize();
   
   for(Int_t i=0;i<iNphotons;i++) fadPhotonAngles[i]=adPhotonAngles[i];
}//inline AliRICHreco::AliRICHreco(Double_t dConeAngle,Double_t dPartOfRing, Double_t adPhotonAngle, Int_t iNphotons)

#endif //AliRICHreco_h
