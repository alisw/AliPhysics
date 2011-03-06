#ifndef ALIJETDUMMYGEO_H
#define ALIJETDUMMYGEO_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//
// Temporarily added to define part of the EMCal geometry
// necessary for the jet finder
// Magali.Estienne@cern.ch
//
#include <TObject.h>
#include <TArrayD.h>
#include <TMath.h>
#include <TVector3.h>

class TGeoMatrix;
class AliJetDummyShishKebabTrd1Module;

class AliJetDummyGeo : public TObject 
{
 public: 
  AliJetDummyGeo();
  AliJetDummyGeo(const AliJetDummyGeo& geom);
  virtual ~AliJetDummyGeo();
  static AliJetDummyGeo* GetInstance() {return new AliJetDummyGeo();}
  static AliJetDummyGeo* GetInstance(const char* /*name*/, const char* /*title*/)
    {return new AliJetDummyGeo();}
  const Char_t* GetNameOfEMCALEnvelope() const {return "XEN1";}
  Float_t GetEnvelop(Int_t index) const { return fEnvelop[index];}
  Float_t AngleFromEta(Float_t eta) const { 
    // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta) const { 
    // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }
  Int_t   GetNCells()               const  {return fNCells;}
  Float_t GetPhiModuleSize()        const  {return fPhiModuleSize;}
  Float_t GetEtaModuleSize()        const  {return fEtaModuleSize;}
  Float_t GetShellThickness()       const  {return fShellThickness;}
  Float_t GetSteelFrontThickness()  const  {return fSteelFrontThick;}
  Float_t GetLongModuleSize()       const  {return fLongModuleSize;}
  Float_t GetTrd1Angle()            const  {return fTrd1Angle;}
  Float_t Get2Trd1Dx2()             const  {return f2Trd1Dx2;}
  Int_t   GetNPhi()                 const  {return fNPhi;}
  Int_t   GetNZ()                   const  {return fNZ ;}
  Int_t   GetNumberOfSuperModules() const  {return fNumberOfSuperModules;}
  Float_t GetArm1EtaMin()           const  {return fArm1EtaMin;}
  Float_t GetArm1EtaMax()           const  {return fArm1EtaMax;}
  Float_t GetArm1PhiMin()           const  {return fArm1PhiMin;}
  Float_t GetArm1PhiMax()           const  {return fArm1PhiMax;}
  Float_t GetIPDistance()           const  {return fIPDistance;} 
  void    EtaPhiFromIndex(Int_t id, Float_t& eta, Float_t& phi) const;
  void    GetGlobal(const Double_t *loc, Double_t *glob, Int_t ind) const;
  void    GetGlobal(Int_t absId, Double_t glob[3]) const;
  void    GetGlobal(Int_t absId, TVector3 &vglob) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t loc[3]) const;
  Bool_t  RelPosCellInSModule(Int_t absId, Double_t &xr, Double_t &yr, Double_t &zr) const;
  Bool_t  CheckAbsCellId(Int_t absId) const;
  Bool_t  GetCellIndex(Int_t absId,Int_t &nSupMod,Int_t &nModule,Int_t &nIphi,Int_t &nIeta) const;
  void    GetCellPhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta, Int_t &iphi, Int_t &ieta) const;
  void    GetModulePhiEtaIndexInSModule(Int_t nSupMod, Int_t nModule,  Int_t &iphim, Int_t &ietam) const;
  Bool_t  GetAbsCellIdFromEtaPhi(Double_t eta,Double_t phi, Int_t &absId) const;
  Bool_t  SuperModuleNumberFromEtaPhi(Double_t eta, Double_t phi, Int_t &nSupMod) const;
  Int_t   GetAbsCellIdFromCellIndexes(Int_t nSupMod, Int_t iphi, Int_t ieta) const;
  void    GetModuleIndexesFromCellIndexesInSModule(Int_t nSupMod, Int_t iphi, Int_t ieta, 
						   Int_t &iphim, Int_t &ietam, Int_t &nModule) const;
  Int_t   GetAbsCellId(Int_t nSupMod, Int_t nModule, Int_t nIphi, Int_t nIeta) const;
  Int_t   GetNumberOfModuleInPhiDirection(Int_t nSupMod)  const
  {
    // inline function
    if(nSupMod>=10) return fNPhi/2;
    else            return fNPhi;
  }

  void    CreateListOfTrd1Modules();
  TList  *GetShishKebabTrd1Modules() const {return fShishKebabTrd1Modules;}
  AliJetDummyShishKebabTrd1Module *GetShishKebabModule(Int_t neta);

  Bool_t  GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const;
  void    GetTransformationForSM();
  Float_t GetSampling() const {return fSampling;}
 private:
  AliJetDummyGeo &operator=(const AliJetDummyGeo &det);
  
 protected:
  Float_t fArm1EtaMin;			// Minimum pseudorapidity position of EMCAL in Eta
  Float_t fArm1EtaMax; 			// Maximum pseudorapidity position of EMCAL in Eta
  Float_t fArm1PhiMin; 			// Minimum angular position of EMCAL in Phi (degrees)
  Float_t fArm1PhiMax;			// Maximum angular position of EMCAL in Phi (degrees)
  Int_t   fNumberOfSuperModules;        // Number of supermodules
  Float_t fSteelFrontThick;             // Thickness of the front stell face of the support box - 9-sep-04
  Float_t fLateralSteelStrip;           // 13-may-05
  Float_t fEnvelop[3];			// the GEANT TUB for the detector 
  Float_t fIPDistance;                  // Radial Distance of the inner surface of the EMCAL
  Float_t fPhiGapForSM;                 // Gap betweeen supermodules in phi direction
  Int_t   fNPhi;      			// Number of Towers in the PHI direction  
  Int_t   fNZ;		         	// Number of Towers in the Z direction
  Float_t fPhiModuleSize;               // Phi -> X
  Float_t fEtaModuleSize;               // Eta -> Y
  Int_t   fNPHIdiv;                     // number phi divizion of module
  Int_t   fNETAdiv;                     // number eta divizion of module
  Float_t fPhiTileSize;                 // Size of phi tile
  Float_t fEtaTileSize;                 // Size of eta tile
  Int_t   fNECLayers;                   // number of scintillator layers
  Float_t fECScintThick;                // cm, Thickness of the scintillators
  Float_t fECPbRadThickness;            // cm, Thickness of the Pb radiators
  Float_t fSampling;                    // Sampling factor
  Float_t fTrd1Angle;                   // angle in x-z plane (in degree) 
  Int_t   fNCellsInModule;              // number cell in module
  Int_t   fNCellsInSupMod;              // number cell in super module
  Int_t   fNCells;                      // number of cells in calo
  Float_t fLongModuleSize;              // Size of long module
  Float_t f2Trd1Dx2;                    // 2*dx2 for TRD1
  Float_t fShellThickness;	        // Total thickness in (x,y) direction
  Float_t fZLength;			// Total length in z direction
  Float_t fEtaMaxOfTRD1;                // max eta in case of TRD1 geometry (see AliEMCALShishKebabTrd1Module)
  Float_t fParSM[3];                    // SM sizes as in GEANT (TRD1)
  TArrayD fPhiBoundariesOfSM;           // phi boundaries of SM in rad; size is fNumberOfSuperModules;
  TArrayD fPhiCentersOfSM;              // phi of centers of SMl size is fNumberOfSuperModules/2
  TGeoMatrix* fMatrixOfSM[12];          //![fNumberOfSuperModules]; get from gGeoManager;
  TArrayD fCentersOfCellsEtaDir;        // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD fCentersOfCellsXDir;          // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  TArrayD fCentersOfCellsPhiDir;        // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  TArrayD fEtaCentersOfCells;           // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); 
                                        // eta depend from phi position; 
  TArrayD fPhiCentersOfCells;           // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
  TList  *fShishKebabTrd1Modules;       //  List of modules
  Int_t   fDebug;                       //  Debug flag 
  ClassDef(AliJetDummyGeo,1)
};
 
#endif
