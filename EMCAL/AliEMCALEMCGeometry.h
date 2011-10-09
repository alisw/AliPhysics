#ifndef ALIEMCALEMCGEOMETRY_H
#define ALIEMCALEMCGEOMETRY_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEMCALEMCGeometry.h 26174 2008-05-26 20:27:16Z jklay $ */

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton
// EMCAL consists of a layers of scintillator, and lead.
//                  
//*-- Author: Sahal Yacoob (LBL / UCT)
//*--   and : Yves Schutz (Subatech)
//*--   and : Aleksei Pavlinov (WSU) - shashlyk staff
//*--   and : Gustavo Conesa: Add TRU mapping. TRU parameters still not fixed.
//*--   and : Magali Estienne (Subatech): class added for new library for EMCALGeoUtils.par file

// --- ROOT system ---
#include <TMath.h>
#include <TArrayD.h>
#include <TNamed.h>
class TString ;
class TObjArray;
class Riostream;

// --- AliRoot header files ---
class AliEMCALEMCGeometry;
class AliEMCALShishKebabTrd1Module;

class AliEMCALEMCGeometry : public TNamed {
public:
  AliEMCALEMCGeometry(); // default ctor only for internal usage (singleton)
  AliEMCALEMCGeometry(const AliEMCALEMCGeometry& geom);
  // ctor only for internal usage (singleton)
  AliEMCALEMCGeometry(const Text_t* name, const Text_t* title,
                      const Text_t* mcname="", const Text_t* mctitle="");

  virtual ~AliEMCALEMCGeometry(void); 

  AliEMCALEMCGeometry & operator = (const AliEMCALEMCGeometry  & /*rvalue*/) {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented");
    return *this;
  };

  //////////
  // General
  //
  Bool_t IsInitialized(void) const { return fgInit ; }
  static const Char_t* GetDefaultGeometryName() {return fgkDefaultGeometryName;}
  void   PrintGeometry();        //*MENU*  
  
  void   Init(const Text_t* mcname="", const Text_t* mctitle=""); // initializes the parameters of EMCAL
  void   CheckAdditionalOptions();        //
  void   DefineSamplingFraction(const Text_t* mcname="", const Text_t* mctitle="");        

  //////////////////////////////////////
  // Return EMCAL geometrical parameters
  //
  
  TString GetGeoName() const {return fGeoName;}
  const Char_t* GetNameOfEMCALEnvelope() const { const Char_t* env = "XEN1"; return env ;}
  Float_t GetArm1PhiMin() const { return fArm1PhiMin ; }
  Float_t GetArm1PhiMax() const { return fArm1PhiMax ; }
  Float_t GetArm1EtaMin() const { return fArm1EtaMin;}
  Float_t GetArm1EtaMax() const { return fArm1EtaMax;}
  Float_t GetIPDistance() const { return fIPDistance;}   
  Float_t GetEnvelop(Int_t index) const { return fEnvelop[index] ; }  
  Float_t GetShellThickness() const { return fShellThickness ; }
  Float_t GetZLength() const { return fZLength ; } 
  Int_t   GetNECLayers() const {return fNECLayers ;}
  Int_t   GetNZ() const {return fNZ ;}
  Int_t   GetNEta() const {return fNZ ;}
  Int_t   GetNPhi() const {return fNPhi ;}
  Float_t GetECPbRadThick()const {return fECPbRadThickness;}
  Float_t GetECScintThick() const {return fECScintThick;}
  Float_t GetSampling() const {return fSampling ; } 
  Int_t   GetNumberOfSuperModules() const {return fNumberOfSuperModules;}
  Float_t GetfPhiGapForSuperModules() const {return fPhiGapForSM;}
  Float_t GetPhiModuleSize() const  {return fPhiModuleSize;}
  Float_t GetEtaModuleSize() const  {return fEtaModuleSize;}
  Float_t GetFrontSteelStrip() const {return fFrontSteelStrip;}
  Float_t GetLateralSteelStrip() const {return fLateralSteelStrip;}
  Float_t GetPassiveScintThick() const {return fPassiveScintThick;}
  Float_t GetPhiTileSize() const {return fPhiTileSize;}
  Float_t GetEtaTileSize() const {return fEtaTileSize;}
  Int_t   GetNPhiSuperModule() const {return fNPhiSuperModule;}
  Int_t   GetNPHIdiv() const {return fNPHIdiv ;}
  Int_t   GetNETAdiv() const {return fNETAdiv ;}
  Int_t   GetNCells()  const {return fNCells;}
  Float_t GetLongModuleSize() const {return fLongModuleSize;}
  Float_t GetTrd1Angle() const {return fTrd1Angle;}
  Float_t Get2Trd1Dx2()  const {return f2Trd1Dx2;}
  Float_t GetEtaMaxOfTRD1() const {return fEtaMaxOfTRD1;}
  Float_t GetTrd1AlFrontThick() const { return fTrd1AlFrontThick;}
  Float_t GetTrd1BondPaperThick() const {return fTrd1BondPaperThick;}
  // --
  Int_t   GetNCellsInSupMod() const {return fNCellsInSupMod;}
  Int_t   GetNCellsInModule()  const {return fNCellsInModule; }
  Int_t   GetKey110DEG()      const {return fKey110DEG;}
  Int_t   GetILOSS() const {return fILOSS;}
  Int_t   GetIHADR() const {return fIHADR;}
    // For gamma(Jet) trigger simulations
  Int_t    GetNTRU() const    {return fNTRUEta*fNTRUPhi ; }  
  Int_t    GetNTRUEta() const {return fNTRUEta ; }  
  Int_t    GetNTRUPhi() const {return fNTRUPhi ; }
  Int_t    GetNEtaSubOfTRU() const {return fNEtaSubOfTRU;}
  Int_t    GetNModulesInTRU() const {return fNModulesInTRUEta*fNModulesInTRUPhi; }
  Int_t    GetNModulesInTRUEta() const {return fNModulesInTRUEta ; }  
  Int_t    GetNModulesInTRUPhi() const {return fNModulesInTRUPhi ; }  

  // --
  Float_t GetDeltaEta() const {return (fArm1EtaMax-fArm1EtaMin)/ ((Float_t)fNZ);}
  Float_t GetDeltaPhi() const {return (fArm1PhiMax-fArm1PhiMin)/ ((Float_t)fNPhi);}
  Int_t   GetNTowers() const {return fNPhi * fNZ ;}
  //
  Double_t GetPhiCenterOfSM(Int_t nsupmod) const;
  Float_t GetSuperModulesPar(Int_t ipar) {return fParSM[ipar];}
  //
  Bool_t   GetPhiBoundariesOfSM   (Int_t nSupMod, Double_t &phiMin, Double_t &phiMax) const;
  Bool_t   GetPhiBoundariesOfSMGap(Int_t nPhiSec, Double_t &phiMin, Double_t &phiMax) const;
  //
  // Local Coordinates of SM
/*   TArrayD  GetCentersOfCellsEtaDir() const {return fCentersOfCellsEtaDir;}        // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm) */
/*   TArrayD  GetCentersOfCellsXDir()   const {return fCentersOfCellsXDir;}          // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm) */
/*   TArrayD  GetCentersOfCellsPhiDir() const {return fCentersOfCellsPhiDir;}        // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm) */
/*   // */
/*   TArrayD  GetEtaCentersOfCells() const {return fEtaCentersOfCells;}           // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position;  */
/*   TArrayD  GetPhiCentersOfCells() const {return fPhiCentersOfCells;}           // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.) */

	static int ParseString(const TString &topt, TObjArray &Opt) ; 

  ///////////////////////////////
  //Geometry data member setters
  //
  void SetNZ(Int_t nz) { fNZ= nz; 
                         printf("SetNZ: Number of modules in Z set to %d", fNZ) ; }
  void SetNPhi(Int_t nphi) { fNPhi= nphi; 
                             printf("SetNPhi: Number of modules in Phi set to %d", fNPhi) ; }
  void SetNTRUEta(Int_t ntru) {fNTRUEta = ntru;
               printf("SetNTRU: Number of TRUs per SuperModule in Etaset to %d", fNTRUEta) ;}
  void SetNTRUPhi(Int_t ntru) {fNTRUPhi = ntru;
              printf("SetNTRU: Number of TRUs per SuperModule in Phi set to %d", fNTRUPhi) ;}
  void SetSampling(Float_t samp) { fSampling = samp; 
                              printf("SetSampling: Sampling factor set to %f", fSampling) ; }

  ///////////////////
  // useful utilities
  //
  Float_t AngleFromEta(Float_t eta) const { // returns theta in radians for a given pseudorapidity
    return 2.0*TMath::ATan(TMath::Exp(-eta));
  }
  Float_t ZFromEtaR(Float_t r,Float_t eta) const { // returns z in for a given
    // pseudorapidity and r=sqrt(x*x+y*y).
    return r/TMath::Tan(AngleFromEta(eta));
  }

  //////////////////////////////////////////////////
  // Obsolete methods to be thrown out when feasible
  Float_t GetGap2Active() const {return  fGap2Active ;}
  Float_t GetSteelFrontThickness() const { return fSteelFrontThick;}
  Float_t GetTrd2AngleY()const {return fTrd2AngleY;}
  Float_t Get2Trd2Dy2()  const {return f2Trd2Dy2;}
  Float_t GetTubsR()     const {return fTubsR;}
  Float_t GetTubsTurnAngle() const {return fTubsTurnAngle;}
  //  Float_t GetIP2ECASection() const { return ( GetIPDistance() + GetAlFrontThickness() 
  //					      + GetGap2Active() ) ; }   
  //////////////////////////////////////////////////

  static Bool_t  fgInit;	        // Tells if geometry has been succesfully set up.
  static const Char_t* fgkDefaultGeometryName; // Default name of geometry

private:

  // Member data

  TString fGeoName;                     //geometry name

  TObjArray *fArrayOpts;                //! array of geometry options
  const char *fkAdditionalOpts[6];  //! some additional options for the geometry type and name
  int  fNAdditionalOpts;     //! size of additional options parameter

  Float_t fECPbRadThickness;		// cm, Thickness of the Pb radiators
  Float_t fECScintThick;		// cm, Thickness of the scintillators
  Int_t   fNECLayers;			// number of scintillator layers
  
  Float_t fArm1PhiMin; 			// Minimum angular position of EMCAL in Phi (degrees)
  Float_t fArm1PhiMax;			// Maximum angular position of EMCAL in Phi (degrees)
  Float_t fArm1EtaMin;			// Minimum pseudorapidity position of EMCAL in Eta
  Float_t fArm1EtaMax; 			// Maximum pseudorapidity position of EMCAL in Eta
  
  // Geometry Parameters
  Float_t fEnvelop[3];			// the GEANT TUB for the detector 
  Float_t fIPDistance;			// Radial Distance of the inner surface of the EMCAL
  Float_t fShellThickness;		// Total thickness in (x,y) direction
  Float_t fZLength;			// Total length in z direction
  Int_t   fNZ;				// Number of Towers in the Z direction
  Int_t   fNPhi;			// Number of Towers in the PHI direction
  Float_t fSampling;			// Sampling factor

  // Shish-kebab option - 23-aug-04 by PAI; COMPACT, TWIST, TRD1 and TRD2
  Int_t   fNumberOfSuperModules;         // default is 12 = 6 * 2 
  Float_t fFrontSteelStrip;              // 13-may-05
  Float_t fLateralSteelStrip;            // 13-may-05
  Float_t fPassiveScintThick;            // 13-may-05
  Float_t fPhiModuleSize;                // Phi -> X 
  Float_t fEtaModuleSize;                // Eta -> Y
  Float_t fPhiTileSize;                  // Size of phi tile
  Float_t fEtaTileSize;                  // Size of eta tile
  Float_t fLongModuleSize;               // Size of long module
  Int_t   fNPhiSuperModule;              // 6 - number supermodule in phi direction
  Int_t   fNPHIdiv;                      // number phi divizion of module
  Int_t   fNETAdiv;                      // number eta divizion of module
  //
  Int_t   fNCells;                       // number of cells in calo
  Int_t   fNCellsInSupMod;               // number cell in super module
  Int_t   fNCellsInModule;               // number cell in module)
  //TRU parameters
  Int_t   fNTRUEta ;                     // Number of TRUs per module in eta
  Int_t   fNTRUPhi ;                     // Number of TRUs per module in phi
  Int_t   fNModulesInTRUEta;             // Number of modules per TRU in eta 
  Int_t   fNModulesInTRUPhi;             // Number of modules per TRU in phi 
  Int_t   fNEtaSubOfTRU;                 // Number of eta (z) subregiohi

  // TRD1 options - 30-sep-04
  Float_t fTrd1Angle;                    // angle in x-z plane (in degree) 
  Float_t f2Trd1Dx2;                     // 2*dx2 for TRD1
  Float_t fPhiGapForSM;                  // Gap betweeen supermodules in phi direction
  Int_t   fKey110DEG;                    // for calculation abs cell id; 19-oct-05 
  TArrayD fPhiBoundariesOfSM;            // phi boundaries of SM in rad; size is fNumberOfSuperModules;
  TArrayD fPhiCentersOfSM;                // phi of centers of SMl size is fNumberOfSuperModules/2
  Float_t fEtaMaxOfTRD1;                 // max eta in case of TRD1 geometry (see AliEMCALShishKebabTrd1Module)
  // Oct 26,2010
  Float_t fTrd1AlFrontThick;		 // Thickness of the Al front plate  
  Float_t fTrd1BondPaperThick;		 // Thickness of the Bond Paper sheet  
  // Local Coordinates of SM
  TArrayD fCentersOfCellsEtaDir;        // size fNEta*fNETAdiv (for TRD1 only) (eta or z in SM, in cm)
  TArrayD fCentersOfCellsXDir;          // size fNEta*fNETAdiv (for TRD1 only) (       x in SM, in cm)
  TArrayD fCentersOfCellsPhiDir;        // size fNPhi*fNPHIdiv (for TRD1 only) (phi or y in SM, in cm)
  //
  TArrayD fEtaCentersOfCells;           // [fNEta*fNETAdiv*fNPhi*fNPHIdiv], positive direction (eta>0); eta depend from phi position; 
  TArrayD fPhiCentersOfCells;           // [fNPhi*fNPHIdiv] from center of SM (-10. < phi < +10.)
  // Move from AliEMCALv0 - Feb 19, 2006
  TList   *fShishKebabTrd1Modules; //! list of modules
  // Local coordinates of SM for TRD1
  Float_t fParSM[3];       // SM sizes as in GEANT (TRD1)

  Int_t   fILOSS; // Options for Geant (MIP business) - will call in AliEMCAL
  Int_t   fIHADR; // Options for Geant (MIP business) - will call in AliEMCAL

  ////////////////////////////////////////////////////////////
  //Obsolete member data that will be thrown out when feasible
  //
  Float_t fGap2Active;			// Gap between the envelop and the active material
  Float_t fSteelFrontThick;		 // Thickness of the front stell face of the support box - 9-sep-04
  // TRD2 options - 27-jan-07
  Float_t fTrd2AngleY;                   // angle in y-z plane (in degree) 
  Float_t f2Trd2Dy2;                     // 2*dy2 for TRD2
  Float_t fEmptySpace;                   // 2mm om fred drawing
  // Super module as TUBS
  Float_t fTubsR;                        // radius of tubs 
  Float_t fTubsTurnAngle;                // turn angle of tubs in degree

  ///////////////////////////////////////////////////////////

  ClassDef(AliEMCALEMCGeometry, 2) // EMCAL geometry class 
};

#endif // AliEMCALEMCGEOMETRY_H
