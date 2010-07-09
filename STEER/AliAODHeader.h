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

class TGeoHMatrix;
class TString;

class AliAODHeader : public AliVHeader {

 public :
  AliAODHeader();
 
  AliAODHeader(Int_t nRun, UShort_t nBunchX, UInt_t nOrbit, UInt_t nPeriod, const Char_t *title="");
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
	       Double_t *emEnergy,
	       ULong64_t fTriggerMask,
	       UChar_t   fTriggerCluster,
	       UInt_t    fEventType,
	       const Char_t *title="",
	       Int_t nMuons=0,
	       Int_t nDimuons=0);
  
  virtual ~AliAODHeader();
  AliAODHeader(const AliAODHeader& evt); 
  AliAODHeader& operator=(const AliAODHeader& evt);
  
  Int_t     GetRunNumber()          const { return fRunNumber; }
  UShort_t  GetBunchCrossNumber()   const { return fBunchCrossNumber; }
  UInt_t    GetOrbitNumber()        const { return fOrbitNumber; }
  UInt_t    GetPeriodNumber()       const { return fPeriodNumber; }
  ULong64_t GetTriggerMask()        const { return fTriggerMask; }
  UChar_t   GetTriggerCluster()     const { return fTriggerCluster; }
  TString   GetFiredTriggerClasses()const { return fFiredTriggers;}
  UInt_t    GetEventType()          const { return fEventType; }
  Double_t  GetMagneticField()      const { return fMagneticField; }
  Double_t  GetMuonMagFieldScale()  const { return fMuonMagFieldScale; }
  
  Double_t  GetCentrality()         const { return fCentrality; }
  Double_t  GetZDCN1Energy()        const { return fZDCN1Energy; }
  Double_t  GetZDCP1Energy()        const { return fZDCP1Energy; }
  Double_t  GetZDCN2Energy()        const { return fZDCN2Energy; }
  Double_t  GetZDCP2Energy()        const { return fZDCP2Energy; }
  Double_t  GetZDCEMEnergy(Int_t i) const { return fZDCEMEnergy[i]; }
  Int_t     GetRefMultiplicity()    const { return fRefMult; }
  Int_t     GetRefMultiplicityPos() const { return fRefMultPos; }
  Int_t     GetRefMultiplicityNeg() const { return fRefMultNeg; }
  Int_t     GetNumberOfMuons()      const { return fNMuons; }
  Int_t     GetNumberOfDimuons()    const { return fNDimuons; }

  Double_t  GetQTheta(UInt_t i) const;
  UInt_t    GetNQTheta() const { return (UInt_t)fNQTheta; }

  Double_t GetDiamondX() const {return fDiamondXY[0];}
  Double_t GetDiamondY() const {return fDiamondXY[1];}
  Double_t GetDiamondZ() const {return fDiamondZ;}
  Double_t GetSigma2DiamondX() const {return fDiamondCovXY[0];}
  Double_t GetSigma2DiamondY() const {return fDiamondCovXY[2];}
  Double_t GetSigma2DiamondZ() const {return fDiamondSig2Z;}
  void GetDiamondCovXY(Float_t cov[3]) const {
    for(Int_t i=0;i<3;i++) cov[i]=fDiamondCovXY[i]; return;
  }
  
  void SetRunNumber(Int_t nRun)                { fRunNumber = nRun; }
  void SetBunchCrossNumber(UShort_t nBx)       { fBunchCrossNumber = nBx; }
  void SetOrbitNumber(UInt_t nOr)              { fOrbitNumber = nOr; }
  void SetPeriodNumber(UInt_t nPer)            { fPeriodNumber = nPer; }
  void SetTriggerMask(ULong64_t trigMsk)       { fTriggerMask = trigMsk; }
  void SetFiredTriggerClasses(TString trig)    { fFiredTriggers = trig;}
  void SetTriggerCluster(UChar_t trigClus)     { fTriggerCluster = trigClus; }
  void SetEventType(UInt_t evttype)            { fEventType = evttype; }
  void SetMagneticField(Double_t magFld)       { fMagneticField = magFld; }
  void SetMuonMagFieldScale(Double_t magFldScl){ fMuonMagFieldScale = magFldScl; }
  
  void SetCentrality(Double_t cent)            { fCentrality = cent; }
  void SetZDCN1Energy(Double_t n1Energy)       { fZDCN1Energy = n1Energy; }
  void SetZDCP1Energy(Double_t p1Energy)       { fZDCP1Energy = p1Energy; }
  void SetZDCN2Energy(Double_t n2Energy)       { fZDCN2Energy = n2Energy; }
  void SetZDCP2Energy(Double_t p2Energy)       { fZDCP2Energy = p2Energy; }
  void SetZDCEMEnergy(Double_t emEnergy1, Double_t emEnergy2)      
  	{ fZDCEMEnergy[0] = emEnergy1; fZDCEMEnergy[1] = emEnergy2;}
  void SetRefMultiplicity(Int_t refMult)       { fRefMult = refMult; }
  void SetRefMultiplicityPos(Int_t refMultPos) { fRefMultPos = refMultPos; }
  void SetRefMultiplicityNeg(Int_t refMultNeg) { fRefMultNeg = refMultNeg; }
  void SetNumberOfMuons(Int_t nMuons) { fNMuons = nMuons; }
  void SetNumberOfDimuons(Int_t nDimuons) { fNDimuons = nDimuons; }
  
  void SetQTheta(Double_t *QTheta, UInt_t size = 5);  
  void RemoveQTheta();

  void SetDiamond(Float_t xy[2],Float_t cov[3]) { 
    for(Int_t i=0;i<3;i++) {if(i<2) fDiamondXY[i]=xy[i]; fDiamondCovXY[i]=cov[i];}
  }
  void SetDiamondZ(Float_t z, Float_t sig2z){
    fDiamondZ=z; fDiamondSig2Z=sig2z;
  }

  void Print(Option_t* option = "") const;

  void    SetPHOSMatrix(TGeoHMatrix*matrix, Int_t i) {
      if ((i >= 0) && (i < kNPHOSMatrix)) fPHOSMatrix[i] = matrix;
  }
  const TGeoHMatrix* GetPHOSMatrix(Int_t i) const {
      return ((i >= 0) && (i < kNPHOSMatrix)) ? fPHOSMatrix[i] : NULL;
  }
  
  void    SetEMCALMatrix(TGeoHMatrix*matrix, Int_t i) {
      if ((i >= 0) && (i < kNEMCALMatrix)) fEMCALMatrix[i] = matrix;
  }
  const TGeoHMatrix* GetEMCALMatrix(Int_t i) const {
      return ((i >= 0) && (i < kNEMCALMatrix)) ? fEMCALMatrix[i] : NULL;
  }
  
  enum {kNPHOSMatrix = 5};
  enum {kNEMCALMatrix = 12};
  
 private :
  
  Double32_t  fMagneticField;       // Solenoid Magnetic Field in kG
  Double32_t  fMuonMagFieldScale;   // magnetic field scale of muon arm magnet
  Double32_t  fCentrality;          // Centrality
  Double32_t  fZDCN1Energy;         // reconstructed energy in the neutron1 ZDC
  Double32_t  fZDCP1Energy;         // reconstructed energy in the proton1 ZDC
  Double32_t  fZDCN2Energy;         // reconstructed energy in the neutron2 ZDC
  Double32_t  fZDCP2Energy;         // reconstructed energy in the proton2 ZDC
  Double32_t  fZDCEMEnergy[2];      // reconstructed energy in the electromagnetic ZDCs
  Int_t       fNQTheta;             // number of QTheta elements
  Double32_t *fQTheta;              // [fNQTheta] values to store Lee-Yang-Zeros
  ULong64_t   fTriggerMask;         // Trigger Type (mask)
  TString     fFiredTriggers;       // String with fired triggers
  Int_t       fRunNumber;           // Run Number
  Int_t       fRefMult;             // reference multiplicity
  Int_t       fRefMultPos;          // reference multiplicity of positive particles
  Int_t       fRefMultNeg;          // reference multiplicity of negative particles
  Int_t       fNMuons;              // number of muons in the forward spectrometer
  Int_t       fNDimuons;            // number of dimuons in the forward spectrometer
  UInt_t      fEventType;           // Type of Event
  UInt_t      fOrbitNumber;         // Orbit Number
  UInt_t      fPeriodNumber;        // Period Number
  UShort_t    fBunchCrossNumber;    // BunchCrossingNumber
  UChar_t     fTriggerCluster;      // Trigger cluster (mask)
  Double32_t      fDiamondXY[2];    // Interaction diamond (x,y) in RUN
  Double32_t      fDiamondCovXY[3]; // Interaction diamond covariance (x,y) in RUN
  Double32_t      fDiamondZ;        // Interaction diamond (z) in RUN
  Double32_t      fDiamondSig2Z;    // Interaction diamond sigma^2 (z) in RUN
  TGeoHMatrix*    fPHOSMatrix[kNPHOSMatrix];   //PHOS module position and orientation matrices
  TGeoHMatrix*    fEMCALMatrix[kNEMCALMatrix]; //EMCAL supermodule position and orientation matrices

  ClassDef(AliAODHeader,11);
};

#endif
