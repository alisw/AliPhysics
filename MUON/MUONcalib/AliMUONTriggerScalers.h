#ifndef ALIMUONTRIGGERSCALERS_H
#define ALIMUONTRIGGERSCALERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup calib
/// \class AliMUONTriggerScalers
/// \brief MUON trigger scalers
///
//  Author: Bogdan Vulpescu

#include <TObject.h>

class AliMUONTriggerScalers : public TObject {

 public:

   AliMUONTriggerScalers();
   virtual ~AliMUONTriggerScalers();
 
   void SetNCalibEvents(UInt_t nce)               
   { fNCalibEvents = nce; }
   void SetDeltaT(UInt_t dt)                      
   { fDeltaT = dt; }
   void SetGloScaler(UInt_t gs, Int_t i)          
   { fGloScal[i] = gs; }
   void SetLocScalerLpt(UInt_t lst, Int_t i)      
   { fLocScalLpt[i] = lst; }
   void SetLocScalerStrip(ULong64_t lss, Int_t ica, Int_t ich, Int_t ilo) 
   { fLocScalStrip[ica][ich][ilo] = lss; }
   void SetLocScalerStripOver(ULong64_t lsso, Int_t ica, Int_t ich, Int_t ilo) 
   { fLocScalStripOver[ica][ich][ilo] = lsso; }

   UInt_t GetNCalibEvents() const { 
     return fNCalibEvents; }
   UInt_t GetDeltaT() const { 
     return fDeltaT; }
   UInt_t GetGloScal(Int_t i) const { 
     return (i>=0 && i<=5) ? fGloScal[i] : 0; }
   UInt_t GetLocScalLpt(Int_t i) const { 
     return (i>=0 && i<=233) ? fLocScalLpt[i] : 0; }
   ULong64_t GetLocScalStrip(Int_t iCath, Int_t iCha, Int_t iLoc) const { 
     return (iCath>=0 && iCath<=1 && iCha>=0 && iCha <=3 && iLoc >= 0 && iLoc <= 233) ? fLocScalStrip[iCath][iCha][iLoc] : 0; }
   UInt_t GetLocScalStripOver(Int_t iCath, Int_t iCha, Int_t iLoc) const { 
     return (iCath>=0 && iCath<=1 && iCha>=0 && iCha <=3 && iLoc >= 0 && iLoc <= 233) ? fLocScalStripOver[iCath][iCha][iLoc] : 0; }

 private:

   UInt_t    fNCalibEvents;   
   UInt_t    fDeltaT;
   UInt_t    fGloScal[6];
   UInt_t    fLocScalLpt[234];
   ULong64_t fLocScalStrip[2][4][234];
   UInt_t    fLocScalStripOver[2][4][234];

   ClassDef(AliMUONTriggerScalers,1)

};

#endif
