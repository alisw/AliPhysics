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

  virtual Double_t GetZNCEnergy() const;
  virtual Double_t GetZPCEnergy() const;
  virtual Double_t GetZNAEnergy() const;
  virtual Double_t GetZPAEnergy() const;
  virtual Double_t GetZEM1Energy() const {return fZEM1Energy;}
  virtual Double_t GetZEM2Energy() const {return fZEM2Energy;}
  
  virtual const Double_t *GetZNCTowerEnergy() const {return fZNCTowerEnergy;}
  virtual const Double_t *GetZNATowerEnergy() const {return fZNATowerEnergy;}
  virtual const Double_t *GetZPCTowerEnergy() const {return fZPCTowerEnergy;}
  virtual const Double_t *GetZPATowerEnergy() const {return fZPATowerEnergy;}
  virtual const Double_t *GetZNCTowerEnergyLR() const {return fZNCTowerEnergyLR;}
  virtual const Double_t *GetZNATowerEnergyLR() const {return fZNATowerEnergyLR;}
  
  virtual Bool_t GetZNCentroidInPbPb(Float_t beamEne, Double_t centrZNC[2], Double_t centrZNA[2]);
  virtual Bool_t GetZNCentroidInpp(Double_t centrZNC[2], Double_t centrZNA[2]);

  virtual Float_t GetZNCTime() const {return fZNCTDC;}
  virtual Float_t GetZNATime() const {return fZNATDC;}

  virtual Float_t GetZDCTimeSum() const {return fZDCTDCSum;}
  virtual Float_t GetZDCTimeDiff() const {return fZDCTDCDifference;}

  // Setters  
  
  void  SetZEM1Energy(const Double_t zem1) {fZEM1Energy = zem1;}
  void  SetZEM2Energy(const Double_t zem2) {fZEM2Energy = zem2;}
  void  SetZNCTowers(const Double_t value[5], const Double_t valueLG[5]);
  void  SetZNATowers(const Double_t value[5], const Double_t valueLG[5]);
  void  SetZPCTowers(const Double_t value[5])
  	{for(Int_t i=0; i<5; i++) fZPCTowerEnergy[i] = value[i];}
  void  SetZPATowers(const Double_t value[5])
  	{for(Int_t i=0; i<5; i++) fZPATowerEnergy[i] = value[i];}
  
  void  SetZDCParticipants(Int_t npart, Int_t npartA, Int_t npartC) 
  	{fZDCParticipants=npart; fZDCPartSideA=npartA; fZDCPartSideC=npartC;}
  
  void  SetZDCImpactParameter(Float_t b, Float_t bA, Float_t bC)
  	{fImpactParameter=b; fImpactParamSideA=bA; fImpactParamSideC=bC;}
 
  void  SetZDCTDCSum(Float_t tdc)  {fZDCTDCSum = tdc;}
  void  SetZDCTDCDiff(Float_t tdc) {fZDCTDCDifference = tdc;}
  
  void  SetZNCTDC(Float_t tdc) {fZNCTDC = tdc;}
  void  SetZNATDC(Float_t tdc) {fZNATDC = tdc;}
 
 
protected:
  
  mutable Double32_t   fZNCEnergy;   	      //!E in ZNC
  mutable Double32_t   fZNAEnergy;   	      //!E in ZNA
  mutable Double32_t   fZPCEnergy;   	      //!E in ZPC
  mutable Double32_t   fZPAEnergy;   	      //!E in ZPA
  Double32_t   fZEM1Energy;   	      // E in ZEM1
  Double32_t   fZEM2Energy;	      // E in ZEM2
  Double32_t   fZNCTowerEnergy[5];    // E in 5 ZNC sectors - high gain chain
  Double32_t   fZNATowerEnergy[5];    // E in 5 ZNA sectors - high gain chain
  Double32_t   fZPCTowerEnergy[5];    // E in 5 ZPC sectors - high gain chain
  Double32_t   fZPATowerEnergy[5];    // E in 5 ZPA sectors - high gain chain
  Double32_t   fZNCTowerEnergyLR[5];  // E in 5 ZNC sectors - low gain chain
  Double32_t   fZNATowerEnergyLR[5];  // E in 5 ZNA sectors - low gain chain
  //
  Int_t     fZDCParticipants;	   // number of participants estimated by the ZDC (ONLY in A-A)
  Int_t     fZDCPartSideA;	   // number of participants estimated by the ZDC (ONLY in A-A)
  Int_t     fZDCPartSideC;	   // number of participants estimated by the ZDC (ONLY in A-A)
  //
  Double32_t   fImpactParameter;   // impact parameter estimated by the ZDC (ONLY in A-A)
  Double32_t   fImpactParamSideA;  // impact parameter estimated by the ZDC (ONLY in A-A)
  Double32_t   fImpactParamSideC;  // impact parameter estimated by the ZDC (ONLY in A-A)
  //
  Float_t   fZDCTDCSum;	   	   // ZDC TDC sum in ns corrected 4 phase shift
  Float_t   fZDCTDCDifference;	   // ZDC TDC diff. in ns corrected 4 phase shift
  Float_t   fZNCTDC; 	   	   // ZNCC TDC sum in ns corrected 4 phase shift        
  Float_t   fZNATDC;     	   // ZNA TDC diff. in ns corrected 4 phase shift;     


  ClassDef(AliAODZDC,2)
};

#endif
