#ifndef ALIFMDSIMULATOR_H
#define ALIFMDSIMULATOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef TTask
# include <TTask.h>
#endif
#ifndef TLorentzVector
# include <TLorentzVector.h>
#endif
#ifndef TArrayI
# include <TArrayI.h>
#endif
class TVector3;
class AliFMD;
class AliFMDRing;
class AliFMDDetector;
class AliFMD1;
class AliFMD2;
class AliFMD3;
class TObjArray;

/** Simulation of the FMD. 
    This class builds the geometry, and processes hits in the FMD */ 
class AliFMDSimulator : public TTask  
{
public:
  AliFMDSimulator();
  /** CTOR */
  AliFMDSimulator(AliFMD* fmd, Bool_t detailed=kTRUE);
  virtual ~AliFMDSimulator() {}
  /** Initialize */
  virtual void DefineMaterials();
  /** Register */
  virtual void DefineGeometry() = 0;
  /** Deal with a hit in the FMD */
  virtual void Exec(Option_t* option="");
  virtual void EndEvent();
  /** @param use Wheher to use a divided geometry */
  virtual void UseDivided(Bool_t use=kTRUE)  { fUseDivided = use; }
  /** @param use Wheher to assemblies in the geometry definition */
  virtual void UseAssembly(Bool_t use=kTRUE) { fUseAssembly = use; }
  /** Whether to make a detailed geometry or not. 
      @param use If true, make a detailed geometry */
  virtual void SetDetailed(Bool_t use) { fDetailed = use; }
protected:  
  AliFMD*        fFMD;           //! Pointer to module 
  TLorentzVector fCurrentV;      //! Current hit postition 
  TLorentzVector fCurrentP;      //! Current hit momentum
  TArrayI        fActiveId;      //! Active volume ID's
  Int_t          fCurrentPdg;    //! Current hit particle code 
  Double_t       fCurrentDeltaE; //! Current hit energy loss
  
  Bool_t         IsActive(Int_t volId) const;
  Bool_t         VMC2FMD(Int_t copy, TLorentzVector& v,
                         UShort_t& detector, Char_t& ring,
                         UShort_t& sector, UShort_t& stripe);
  Bool_t         VMC2FMD(TLorentzVector& v, UShort_t& detector,
                         Char_t& ring, UShort_t& sector, UShort_t& strip);

  static const Char_t* fgkActiveName;	// Name of Active volumes
  static const Char_t* fgkSectorName;	// Name of Sector volumes
  static const Char_t* fgkStripName;	// Name of Strip volumes
  static const Char_t* fgkModuleName;	// Name of Module volumes
  static const Char_t* fgkPCBName;	// Name of PCB volumes
  static const Char_t* fgkLongLegName;	// Name of LongLeg volumes
  static const Char_t* fgkShortLegName;	// Name of ShortLeg volumes
  static const Char_t* fgkFrontVName;	// Name of Front volumes
  static const Char_t* fgkBackVName;	// Name of Back volumes
  static const Char_t* fgkRingName;	// Name of Ring volumes
  static const Char_t* fgkTopHCName;	// Name of TopHC volumes
  static const Char_t* fgkBotHCName;	// Name of BotHC volumes
  static const Char_t* fgkTopIHCName;	// Name of TopIHC volumes
  static const Char_t* fgkBotIHCName;	// Name of BotIHC volumes
  static const Char_t* fgkNoseName;	// Name of Nose volumes
  static const Char_t* fgkBackName;	// Name of Back volumes
  static const Char_t* fgkBeamName;	// Name of Beam volumes
  static const Char_t* fgkFlangeName;	// Name of Flange volumes
  
  enum {
    kSiId,                 // ID index of Si medium
    kAirId,                // ID index of Air medium
    kPlasticId,            // ID index of Plastic medium
    kPcbId,                // ID index of PCB medium
    kSiChipId,             // ID index of Si Chip medium
    kAlId,                 // ID index of Al medium
    kCarbonId,             // ID index of Carbon medium
    kCopperId,             // ID index of Copper Medium
    kKaptonId              // ID index of Kapton Medium
  };  

  Bool_t      fDetailed;      // Whether to make a detailed simulation 
  Bool_t      fUseDivided;    // Divided volumes
  Bool_t      fUseAssembly;   // Assembly volumes
  Int_t fSectorOff;        // Sector offset in volume tree 
  Int_t fModuleOff;        // Module offset in volume tree
  Int_t fRingOff;          // Ring offset in the volume tree 
  Int_t fDetectorOff;      // Detector offfset in the volume tree 
  TObjArray* fBad;         //! List of bad hits
  
  ClassDef(AliFMDSimulator,0) // Simulation class for the FMD
};



#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
