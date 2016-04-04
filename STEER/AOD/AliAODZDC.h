#ifndef ALIAODZDC_H
#define ALIAODZDC_H

//-------------------------------------------------------------------------
//     Class for AOD ZDC data
//     Author: Chiara Oppedisano
//     Chiara.Oppedisano@cern.ch March 2011
//-------------------------------------------------------------------------

#include <AliVZDC.h>

class AliAODZDC : public AliVZDC {
public:

  AliAODZDC();
  AliAODZDC(const AliAODZDC& zdcAOD);
  AliAODZDC &operator=(const AliAODZDC& zdcAOD);
 
  virtual ~AliAODZDC() {};
 
  // Getters  
   
  virtual Short_t  GetZDCParticipants() const {return fZDCParticipants;}
  virtual Short_t  GetZDCPartSideA()	const {return fZDCPartSideA;}
  virtual Short_t  GetZDCPartSideC()	const {return fZDCPartSideC;}
  virtual Double_t GetImpactParameter()  const {return fImpactParameter;}
  virtual Double_t GetImpactParamSideA() const {return fImpactParamSideA;}
  virtual Double_t GetImpactParamSideC() const {return fImpactParamSideC;}

  virtual Double_t GetZNCEnergy() const {return fZNCEnergy;}
  virtual Double_t GetZPCEnergy() const {return fZPCEnergy;}
  virtual Double_t GetZNAEnergy() const {return fZNAEnergy;}
  virtual Double_t GetZPAEnergy() const {return fZPAEnergy;}
  virtual Double_t GetZEM1Energy() const {return fZEM1Energy;}
  virtual Double_t GetZEM2Energy() const {return fZEM2Energy;}
  
  virtual const Double_t *GetZNCTowerEnergy() const {return fZNCTowerEnergy;}
  virtual const Double_t *GetZNATowerEnergy() const {return fZNATowerEnergy;}
  virtual const Double_t *GetZPCTowerEnergy() const {return fZPCTowerEnergy;}
  virtual const Double_t *GetZPATowerEnergy() const {return fZPATowerEnergy;}
  virtual const Double_t *GetZNCTowerEnergyLR() const {return fZNCTowerEnergyLR;}
  virtual const Double_t *GetZNATowerEnergyLR() const {return fZNATowerEnergyLR;}
  virtual const Double_t *GetZPCTowerEnergyLR() const {return fZPCTowerEnergyLR;}
  virtual const Double_t *GetZPATowerEnergyLR() const {return fZPATowerEnergyLR;}
  
  virtual Bool_t GetZNCentroidInPbPb(Float_t beamEne, Double_t centrZNC[2], Double_t centrZNA[2]);
  virtual Bool_t GetZNCentroidInpp(Double_t centrZNC[2], Double_t centrZNA[2]);

  // Setters dealing only with the 1st stored TDC hit
  virtual Float_t GetZNCTime() const {return fZNCTDC;}
  virtual Float_t GetZNATime() const {return fZNATDC;}
  virtual Float_t GetZPCTime() const {return fZPCTDC;}
  virtual Float_t GetZPATime() const {return fZPATDC;}
  //
  virtual Float_t GetZDCTimeSum() const {return fZDCTDCSum;}
  virtual Float_t GetZDCTimeDiff() const {return fZDCTDCDifference;}

  // Jan.2016: propagating multi-hit structure of TDC hits to AODs
  virtual Float_t GetZNCTDCm(Int_t i) const {return fZNCTDCm[i];}
  virtual Float_t GetZNATDCm(Int_t i) const {return fZNATDCm[i];}
  virtual Float_t GetZPCTDCm(Int_t i) const {return fZPCTDCm[i];}
  virtual Float_t GetZPATDCm(Int_t i) const {return fZPATDCm[i];}
  //
  virtual Bool_t GetTDCSum(Float_t sum[4]); 
  virtual Bool_t GetTDCDiff(Float_t diff[4]); 
  //
  virtual Bool_t  IsZNAfired() {return fIsZNAfired;}
  virtual Bool_t  IsZNCfired() {return fIsZNCfired;}
  virtual Bool_t  IsZNANDfired() {if(IsZNAfired() &&IsZNCfired()) return kTRUE;
  				  else return kFALSE;}

  // Setters  
  
  void  SetZEM1Energy(const Double_t zem1) {fZEM1Energy = zem1;}
  void  SetZEM2Energy(const Double_t zem2) {fZEM2Energy = zem2;}
  void  SetZNCTowers(const Double_t value[5], const Double_t valueLG[5]);
  void  SetZNATowers(const Double_t value[5], const Double_t valueLG[5]);
  void  SetZPCTowers(const Double_t value[5], const Double_t valueLG[5]);
  void  SetZPATowers(const Double_t value[5], const Double_t valueLG[5]);
  
  void  SetZDCParticipants(Int_t npart, Int_t npartA, Int_t npartC) 
  	{fZDCParticipants=npart; fZDCPartSideA=npartA; fZDCPartSideC=npartC;}  
  void  SetZDCImpactParameter(Float_t b, Float_t bA, Float_t bC)
  	{fImpactParameter=b; fImpactParamSideA=bA; fImpactParamSideC=bC;}
 
  // Setters dealing only with the 1st stored TDC hit
  void  SetZDCTDCSum(Float_t tdc)  {fZDCTDCSum = tdc;}
  void  SetZDCTDCDiff(Float_t tdc) {fZDCTDCDifference = tdc;}
  //
  void  SetZNCTDC(Float_t tdc) {fZNCTDC = tdc;}
  void  SetZNATDC(Float_t tdc) {fZNATDC = tdc;}
  void  SetZPCTDC(Float_t tdc) {fZPCTDC = tdc;}
  void  SetZPATDC(Float_t tdc) {fZPATDC = tdc;}

  // Jan.2016: propagating multi-hit structure of TDC hits to AODs
  void  SetZNCTDCm(Int_t i, Float_t tdc) {fZNCTDCm[i] = tdc;}
  void  SetZNATDCm(Int_t i, Float_t tdc) {fZNATDCm[i] = tdc;}
  void  SetZPCTDCm(Int_t i, Float_t tdc) {fZPCTDCm[i] = tdc;}
  void  SetZPATDCm(Int_t i, Float_t tdc) {fZPATDCm[i] = tdc;}
  //
  void  SetZNAfired() {fIsZNAfired = kTRUE;}
  void  SetZNCfired() {fIsZNCfired = kTRUE;}
  void  SetZPAfired() {fIsZPAfired = kTRUE;}
  void  SetZPCfired() {fIsZPCfired = kTRUE;}
 
protected:
  
  mutable Double32_t   fZNCEnergy;    //!E in ZNC
  mutable Double32_t   fZNAEnergy;    //!E in ZNA
  mutable Double32_t   fZPCEnergy;    //!E in ZPC
  mutable Double32_t   fZPAEnergy;    //!E in ZPA
  Double32_t   fZEM1Energy;   	      // E in ZEM1
  Double32_t   fZEM2Energy;	      // E in ZEM2
  Double32_t   fZNCTowerEnergy[5];    // E in 5 ZNC sectors - high gain chain
  Double32_t   fZNATowerEnergy[5];    // E in 5 ZNA sectors - high gain chain
  Double32_t   fZPCTowerEnergy[5];    // E in 5 ZPC sectors - high gain chain
  Double32_t   fZPATowerEnergy[5];    // E in 5 ZPA sectors - high gain chain
  Double32_t   fZNCTowerEnergyLR[5];  // E in 5 ZNC sectors - low gain chain
  Double32_t   fZNATowerEnergyLR[5];  // E in 5 ZNA sectors - low gain chain
  Double32_t   fZPCTowerEnergyLR[5];  // E in 5 ZPC sectors - low gain chain
  Double32_t   fZPATowerEnergyLR[5];  // E in 5 ZPA sectors - low gain chain
  //
  Int_t     fZDCParticipants;	   // number of participants estimated by the ZDC (ONLY in A-A)
  Int_t     fZDCPartSideA;	   // number of participants estimated by the ZDC (ONLY in A-A)
  Int_t     fZDCPartSideC;	   // number of participants estimated by the ZDC (ONLY in A-A)
  //
  Double32_t   fImpactParameter;   // impact parameter estimated by the ZDC (ONLY in A-A)
  Double32_t   fImpactParamSideA;  // impact parameter estimated by the ZDC (ONLY in A-A)
  Double32_t   fImpactParamSideC;  // impact parameter estimated by the ZDC (ONLY in A-A)
  //
  // These data members deals only with the 1st stored TDC hit
  Float_t   fZDCTDCSum;	   	   // ZDC TDC sum in ns corrected 4 phase shift
  Float_t   fZDCTDCDifference;	   // ZDC TDC diff. in ns corrected 4 phase shift
  Float_t   fZNCTDC; 	   	   // ZNC TDC in ns corrected 4 phase shift        
  Float_t   fZNATDC;     	   // ZNA TDC in ns corrected 4 phase shift;     
  Float_t   fZPCTDC; 	   	   // ZNC TDC in ns corrected 4 phase shift        
  Float_t   fZPATDC;     	   // ZNA TDC in ns corrected 4 phase shift;     
  //
  // Jan.2016: propagating multi-hit structure of TDC hits to AODs
  Float_t   fZNCTDCm[4];	// true if ZNC TDC has at least 1 hit
  Float_t   fZNATDCm[4];	// true if ZNA TDC has at least 1 hit
  Float_t   fZPCTDCm[4];	// true if ZPC TDC has at least 1 hit
  Float_t   fZPATDCm[4];	// true if ZPA TDC has at least 1 hit
  //
  Bool_t    fIsZNAfired;	// if true ZNA is fired in the event 
  Bool_t    fIsZNCfired;	// if true ZNC is fired in the event 
  Bool_t    fIsZPAfired;	// if true ZPA is fired in the event 
  Bool_t    fIsZPCfired;	// if true ZPC is fired in the event
  

  ClassDef(AliAODZDC,4)
};

#endif
