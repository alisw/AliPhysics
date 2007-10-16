#ifndef AliAODHeader_H
#define AliAODHeader_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD event header class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include "AliVHeader.h"
#include "AliAODVertex.h"

class AliAODHeader : public AliVHeader {

 public :
  AliAODHeader();
 
  AliAODHeader(Int_t nRun, UShort_t nBunchX, UInt_t nOrbit, UInt_t nPeriod, Char_t *title="");
  AliAODHeader(Int_t nRun, 
	       UShort_t nBunchX,
	       UInt_t nOrbit,
	       UInt_t nPeriod,
	       Int_t refMult,
	       Int_t refMultPos,
	       Int_t refMultNeg,
	       Double_t magField,
	       Double_t muonMagFieldScale,
	       Double_t cent,
	       Double_t n1Energy,
	       Double_t p1Energy,
	       Double_t n2Energy,
	       Double_t p2Energy,
	       Double_t emEnergy,
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
  UInt_t    GetPeriodNumber()       const { return fPeriodNumber; }
  ULong64_t GetTriggerMask()        const { return fTriggerMask; }
  UChar_t   GetTriggerCluster()     const { return fTriggerCluster; }
  UInt_t    GetEventType()          const { return fEventType; }
  Double_t  GetMagneticField()      const { return fMagneticField; }
  Double_t  GetMuonMagFieldScale()  const { return fMuonMagFieldScale; }
  
  Double_t  GetCentrality()         const { return fCentrality; }
  Double_t  GetZDCN1Energy()        const { return fZDCN1Energy; }
  Double_t  GetZDCP1Energy()        const { return fZDCP1Energy; }
  Double_t  GetZDCN2Energy()        const { return fZDCN2Energy; }
  Double_t  GetZDCP2Energy()        const { return fZDCP2Energy; }
  Double_t  GetZDCEMEnergy()        const { return fZDCEMEnergy; }
  Int_t     GetRefMultiplicity()    const { return fRefMult; }
  Int_t     GetRefMultiplicityPos() const { return fRefMultPos; }
  Int_t     GetRefMultiplicityNeg() const { return fRefMultNeg; }

  Double_t  GetQTheta(UInt_t i) const;
  UInt_t    GetNQTheta() const { return (UInt_t)fNQTheta; }
  
  void SetRunNumber(Int_t nRun)                { fRunNumber = nRun; }
  void SetBunchCrossNumber(UShort_t nBx)       { fBunchCrossNumber = nBx; }
  void SetOrbitNumber(UInt_t nOr)              { fOrbitNumber = nOr; }
  void SetPeriodNumber(UInt_t nPer)            { fPeriodNumber = nPer; }
  void SetTriggerMask(ULong64_t trigMsk)       { fTriggerMask = trigMsk; }
  void SetTriggerCluster(UChar_t trigClus)     { fTriggerCluster = trigClus; }
  void SetEventType(UInt_t evttype)            { fEventType = evttype; }
  void SetMagneticField(Double_t magFld)       { fMagneticField = magFld; }
  void SetMuonMagFieldScale(Double_t magFldScl){ fMuonMagFieldScale = magFldScl; }
  
  void SetCentrality(Double_t cent)            { fCentrality = cent; }
  void SetZDCN1Energy(Double_t n1Energy)       { fZDCN1Energy = n1Energy; }
  void SetZDCP1Energy(Double_t p1Energy)       { fZDCP1Energy = p1Energy; }
  void SetZDCN2Energy(Double_t n2Energy)       { fZDCN2Energy = n2Energy; }
  void SetZDCP2Energy(Double_t p2Energy)       { fZDCP2Energy = p2Energy; }
  void SetZDCEMEnergy(Double_t emEnergy)       { fZDCEMEnergy = emEnergy; }
  void SetRefMultiplicity(Int_t refMult)       { fRefMult = refMult; }
  void SetRefMultiplicityPos(Int_t refMultPos) { fRefMultPos = refMultPos; }
  void SetRefMultiplicityNeg(Int_t refMultNeg) { fRefMultNeg = refMultNeg; }
  
  void SetQTheta(Double_t *QTheta, UInt_t size = 5);  
  void RemoveQTheta();

  void Print(Option_t* option = "") const;
  
  
 private :
  
  Double32_t  fMagneticField;       // Solenoid Magnetic Field in kG
  Double32_t  fMuonMagFieldScale;   // magnetic field scale of muon arm magnet
  Double32_t  fCentrality;          // Centrality
  Double32_t  fZDCN1Energy;         // reconstructed energy in the neutron1 ZDC
  Double32_t  fZDCP1Energy;         // reconstructed energy in the proton1 ZDC
  Double32_t  fZDCN2Energy;         // reconstructed energy in the neutron2 ZDC
  Double32_t  fZDCP2Energy;         // reconstructed energy in the proton2 ZDC
  Double32_t  fZDCEMEnergy;         // reconstructed energy in the electromagnetic ZDC
  Int_t       fNQTheta;             // number of QTheta elements
  Double32_t *fQTheta;              // [fNQTheta] values to store Lee-Yang-Zeros
  ULong64_t   fTriggerMask;         // Trigger Type (mask)
  Int_t       fRunNumber;           // Run Number
  Int_t       fRefMult;             // reference multiplicity
  Int_t       fRefMultPos;          // reference multiplicity of positive particles
  Int_t       fRefMultNeg;          // reference multiplicity of negative particles
  UInt_t      fEventType;           // Type of Event
  UInt_t      fOrbitNumber;         // Orbit Number
  UInt_t      fPeriodNumber;        // Period Number
  UShort_t    fBunchCrossNumber;    // BunchCrossingNumber
  UChar_t     fTriggerCluster;      // Trigger cluster (mask)
  
  ClassDef(AliAODHeader,5);
};

#endif
