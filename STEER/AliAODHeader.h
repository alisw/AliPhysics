#ifndef AliAODHeader_H
#define AliAODHeader_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD event base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include "AliAODVertex.h"

class AliAODHeader : public TNamed {

 public :
  AliAODHeader();
 
  AliAODHeader(Int_t nRun, UShort_t nBunchX, UInt_t nOrbit,Char_t *title="");
  AliAODHeader(Int_t nRun, 
	       UShort_t nBunchX,
	       UInt_t nOrbit,
	       Int_t refMult,
	       Int_t refMultPos,
	       Int_t refMultNeg,
	       Double_t magField,
	       Double_t cent, 
	       ULong64_t fTriggerMask,
	       UChar_t   fTriggerCluster,
	       UInt_t    fEventType,
	       Char_t *title="");
  
  virtual ~AliAODHeader();
  AliAODHeader(const AliAODHeader& evt); 
  AliAODHeader& operator=(const AliAODHeader& evt);

  Int_t     GetRunNumber()          const { return fRunNumber; }
  UShort_t  GetBunchCrossNumber()   const { return fBunchCrossNumber; }
  UInt_t    GetOrbitNumber()        const { return fOrbitNumber; }
  ULong64_t GetTriggerMask()        const { return fTriggerMask; }
  UChar_t   GetTriggerCluster()     const { return fTriggerCluster; }
  UInt_t    GetEventType()          const { return fEventType; }
  Double_t  GetMagneticField()      const { return fMagneticField; }

  Double_t  GetCentrality()         const { return fCentrality; }
  Int_t     GetRefMultiplicity()    const { return fRefMult; }
  Int_t     GetRefMultiplicityPos() const { return fRefMultPos; }
  Int_t     GetRefMultiplicityNeg() const { return fRefMultNeg; }

  void SetRunNumber(Int_t nRun)                { fRunNumber = nRun; }
  void SetBunchCrossNumber(UShort_t nBx)       { fBunchCrossNumber = nBx; }
  void SetOrbitNumber(Int_t nOr)               { fOrbitNumber = nOr; }
  void SetTriggerMask(ULong64_t trigMsk)       { fTriggerMask = trigMsk; }
  void SetTriggerCluster(UChar_t trigClus)     { fTriggerCluster = trigClus; }
  void SetEventType(UInt_t evttype)            { fEventType = evttype; }
  void SetMagneticField(Double_t magFld)       { fMagneticField = magFld; }

  void SetCentrality(Double_t cent)            { fCentrality = cent; }
  void SetRefMultiplicity(Int_t refMult)       { fRefMult = refMult; }
  void SetRefMultiplicityPos(Int_t refMultPos) { fRefMultPos = refMultPos; }
  void SetRefMultiplicityNeg(Int_t refMultNeg) { fRefMultNeg = refMultNeg; }

  void Print(Option_t* option = "") const;


 private :

   Double32_t  fMagneticField;       // Solenoid Magnetic Field in kG
   Double32_t  fCentrality;          // Centrality
   ULong64_t   fTriggerMask;         // Trigger Type (mask)
   UInt_t      fEventType;           // Type of Event
   UInt_t      fOrbitNumber;         // Orbit Number
   UShort_t    fBunchCrossNumber;    // BunchCrossingNumber
   Int_t       fRunNumber;           // Run Number
   Int_t       fRefMult;             // reference multiplicity
   Int_t       fRefMultPos;          // reference multiplicity of positive particles
   Int_t       fRefMultNeg;          // reference multiplicity of negative particles
   UChar_t     fTriggerCluster;      // Trigger cluster (mask)

   ClassDef(AliAODHeader,1);
};

#endif
