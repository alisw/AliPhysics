// -*- mode: C++ -*- 
#ifndef ALIESDZDC_H
#define ALIESDZDC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                      Implementation of   Class AliESDZDC
//   This is a class that summarizes the ZDC data for the ESD   
//   Origin: Christian Klein-Boesing, CERN, Christian.Klein-Boesing@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TMath.h>


class AliESDZDC: public TObject {
public:
  AliESDZDC();
  AliESDZDC(const AliESDZDC& zdc);
  AliESDZDC& operator=(const AliESDZDC& zdc);
  virtual void Copy(TObject &obj) const;

  Double_t GetZDCN1Energy() const {return fZDCN1Energy;}
  Double_t GetZDCP1Energy() const {return fZDCP1Energy;}
  Double_t GetZDCN2Energy() const {return fZDCN2Energy;}
  Double_t GetZDCP2Energy() const {return fZDCP2Energy;}
  Double_t GetZDCEMEnergy(Int_t i) const 
  	   {if(i==0){return fZDCEMEnergy;} else if(i==1){return fZDCEMEnergy1;}
	   return 0;}
  Short_t  GetZDCParticipants() const {return fZDCParticipants;}
  Short_t  GetZDCParticipants2() const {return fZDCParticipants2;}
  const Double_t * GetZN1TowerEnergy() const {return fZN1TowerEnergy;}
  const Double_t * GetZN2TowerEnergy() const {return fZN2TowerEnergy;}
  const Double_t * GetZP1TowerEnergy() const {return fZP1TowerEnergy;}
  const Double_t * GetZP2TowerEnergy() const {return fZP2TowerEnergy;}
  const Double_t * GetZN1TowerEnergyLR() const {return fZN1TowerEnergyLR;}
  const Double_t * GetZN2TowerEnergyLR() const {return fZN2TowerEnergyLR;}
  const Double_t * GetZP1TowerEnergyLR() const {return fZP1TowerEnergyLR;}
  const Double_t * GetZP2TowerEnergyLR() const {return fZP2TowerEnergyLR;}
  //
  Double32_t * GetZNCCentroid();
  Double32_t * GetZNACentroid();
  //
  void  SetZDC(Double_t n1Energy, Double_t p1Energy, 
  		Double_t emEnergy0, Double_t emEnergy1,
	       	Double_t n2Energy, Double_t p2Energy, 
		Short_t participants, Short_t participants2) 
   	{fZDCN1Energy=n1Energy; fZDCP1Energy=p1Energy; 
    	fZDCEMEnergy=emEnergy0; fZDCEMEnergy1=emEnergy1;
    	fZDCN2Energy=n2Energy; fZDCP2Energy=p2Energy; 
	fZDCParticipants=participants; fZDCParticipants2=participants2;}
  //
  void  SetZN1TowerEnergy(Float_t tow1[5])
          {for(Int_t i=0; i<5; i++) fZN1TowerEnergy[i] = tow1[i];}
  void  SetZN2TowerEnergy(Float_t tow2[5])
          {for(Int_t i=0; i<5; i++) fZN2TowerEnergy[i] = tow2[i];}
  void  SetZP1TowerEnergy(Float_t tow1[5])
          {for(Int_t i=0; i<5; i++) fZP1TowerEnergy[i] = tow1[i];}
  void  SetZP2TowerEnergy(Float_t tow2[5])
          {for(Int_t i=0; i<5; i++) fZP2TowerEnergy[i] = tow2[i];}
  void  SetZN1TowerEnergyLR(Float_t tow1[5])
          {for(Int_t i=0; i<5; i++) fZN1TowerEnergyLR[i] = tow1[i];}
  void  SetZN2TowerEnergyLR(Float_t tow2[5])
          {for(Int_t i=0; i<5; i++) fZN2TowerEnergyLR[i] = tow2[i];}
  void  SetZP1TowerEnergyLR(Float_t tow1[5])
          {for(Int_t i=0; i<5; i++) fZP1TowerEnergyLR[i] = tow1[i];}
  void  SetZP2TowerEnergyLR(Float_t tow2[5])
          {for(Int_t i=0; i<5; i++) fZP2TowerEnergyLR[i] = tow2[i];}
  void  SetZNACentroid(Float_t centrCoord[2])
  	   {for(Int_t i=0; i<2; i++) fZNACentrCoord[i] = centrCoord[i];}
  void  SetZNCCentroid(Float_t centrCoord[2])
  	   {for(Int_t i=0; i<2; i++) fZNCCentrCoord[i] = centrCoord[i];}

  void    Reset();
  void    Print(const Option_t *opt=0) const;

private:

  Double32_t   fZDCN1Energy;  // reconstructed energy in the neutron ZDC
  Double32_t   fZDCP1Energy;  // reconstructed energy in the proton ZDC
  Double32_t   fZDCN2Energy;  // reconstructed energy in the neutron ZDC
  Double32_t   fZDCP2Energy;  // reconstructed energy in the proton ZDC
  Double32_t   fZDCEMEnergy;  // signal in the electromagnetic ZDCs
  Double32_t   fZDCEMEnergy1; // second EM signal,cannot change fZDCEMEnergy to array (not backward compatible)
  Double32_t   fZN1TowerEnergy[5];// reco E in 5 ZN1 sectors - high gain chain
  Double32_t   fZN2TowerEnergy[5];// reco E in 5 ZN2 sectors - high gain chain
  Double32_t   fZP1TowerEnergy[5];// reco E in 5 ZP1 sectors - high gain chain
  Double32_t   fZP2TowerEnergy[5];// reco E in 5 ZP2 sectors - high gain chain
  Double32_t   fZN1TowerEnergyLR[5];// reco E in 5 ZN1 sectors - low gain chain
  Double32_t   fZN2TowerEnergyLR[5];// reco E in 5 ZN2 sectors - low gain chain
  Double32_t   fZP1TowerEnergyLR[5];// reco E in 5 ZP1 sectors - low gain chain
  Double32_t   fZP2TowerEnergyLR[5];// reco E in 5 ZP2 sectors - low gain chain
  Short_t      fZDCParticipants;// number of participants estimated by the ZDC - side C
  Short_t      fZDCParticipants2;// number of participants estimated by the ZDC - side A
  Double32_t   fZNACentrCoord[2];// Coordinates of the centroid over ZNC
  Double32_t   fZNCCentrCoord[2];// Coordinates of the centroid over ZNA
  //
  ClassDef(AliESDZDC,9)
};

#endif

