#ifndef ALITOFPAD_H
#define ALITOFPAD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////
//  TOF Class used in reconstruction          
//  AliTOFPad  class                          
//  (see implementation file for details)                               
//
//-- Authors: Bologna-ITEP-Salerno Group
////////////////////////////////////////////////////////////////


#include "TObject.h"

//_______________________________________________________
class AliTOFPad : public TObject{

 public:
   AliTOFPad();
   AliTOFPad(Int_t sector, Int_t plate, Int_t strip, Int_t pixel);
   ~AliTOFPad(){};
   void SetGeom (Int_t sector, Int_t plate, Int_t strip, Int_t pixel);
   void SetTofChargeHit(Float_t realTime, Float_t charge, Float_t geantTime, Int_t hitnum);
   void SetTrack(Int_t   track)        {fTrack=track;}
   void SetTrackMatched(Int_t track)   {fTrackMatched=track;}
   void AddState(Int_t   state)        {fState+=state;}
   void SetRealTime (Float_t realTime) {fRealTime=realTime;}
   void SetGeantTime(Float_t geantTime){fGeantTime=geantTime;}
   void SetCharge (Float_t charge)     {fCharge=charge;}
   void SetAverageTime(Float_t averageTime) {fAverageTime=averageTime;}
   void SetHit  (Int_t   hit)          {fHit=hit;}
   Int_t   GetSector()       const {return fSector;}
   Int_t   GetPlate()        const {return fPlate;}
   Int_t   GetStrip()        const {return fStrip;}
   Int_t   GetPixel()        const {return fPixel;}
   Int_t   GetTrack()        const {return fTrack;}
   Int_t   GetTrackMatched() const {return fTrackMatched;}
   Int_t   GetState()        const {return fState;}
   Float_t GetRealTime()     const {return fRealTime;}
   Float_t GetGeantTime()    const {return fGeantTime;}
   Float_t GetCharge()       const {return fCharge;}
   Float_t GetAverageTime()  const {return fAverageTime;}
   Int_t   GetHit()          const {return fHit;}
  
private:
   Int_t   fSector, fPlate, fStrip, fPixel; // sector, plate, strip and pad number
   Int_t   fTrack;           // track number of first track fired the pixel
   Int_t   fTrackMatched;    // track number i of TrackArray[i-1] matched with the pixel
   Int_t   fState;           // =1, if the pixel is fired by the track
   Float_t fRealTime;        // real time [ns] given by the pad
   Float_t fGeantTime;       // GEANT3.21 time [ns] i.e. true time
   Float_t fCharge;          // charge related to the pad
   Float_t fAverageTime;     // average time of the pad cluster due to the Edge Effect [ns]
   Int_t   fHit;             // hit number khit of HitArray[khit-1] which fired the pixel
    
   ClassDef(AliTOFPad,1)   // TOF Class used in reconstruction
};

#endif /* ALITOFPAD_H */
