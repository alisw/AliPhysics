#ifndef _ALIPARTICLEYIELD_H_
#define _ALIPARTICLEYIELD_H_

// AliParticleYield
// This class implements a container for particle yields results, to
// be used e.g. in thermal model fits, with utility methods to
// read/write to ASCII files or root trees
// Author: Michele Floris, michele.floris@cern.ch

#include "TObject.h"
#include "TString.h"
#include "TCut.h"

class TClonesArray;
class TTree;
class AliParticleYield : public TObject
{
public:

  enum AliPYCSystem_t {kCSpp   = 0, 
		       kCSpPb  = 1, 
		       kCSPbPb = 2 } ;

  static const char * kSystemString[];

  enum AliPYStatusCode_t {
    kSCPublished   = 0,
    kSCPreliminary = 1,
    kSCFinalNotPublished = 2,
    kSCMayChange   = 3,
  };

  static const char * kStatusString[];

  enum AliPYMeasurementType_t { // this is a bit mask: more than one bit can be one (be careful not to set mutually exclusive ones).
    // Type of measurements
    kTypeLinearInterpolation = 0x1,
    kTypeParticleRatio       = 0x2,  // If true, this is a ratio of 2 particles where the propagation of uncertainty was properly taken into account. 
    kTypeAverageAndRefit     = 0x4,  // this means that the measurement has been obtained summing several spectra in smaller centality bins (weihgted by the width of the centrality bin) and refitting them
    kTypeExtrPionRatio       = 0x8,  // Extrapolated from a different centrality bin, assumin the ratio to pions is constant
    kTypeExtrFit             = 0x20, // Extrapolated fitting the centrality dependence vs npart

    // Type of errors
    kTypeOnlyTotError        = 0x10, // If on, only the total error is returned as "GetSystError". GetStatError should be set to 0;

    // Additional flags
    kTypeAveragePartAntiPart = 0x100, // Can only be set if isSum = 1. It indicates that it is an averrage rather than a sum TODO: add separate bits for different averages (e.g. quadratic)?           

  };
  
  AliParticleYield();
  AliParticleYield(Int_t pdg, Int_t system, Float_t sqrts, Float_t value, Float_t stat, Float_t syst, Float_t norm, Float_t ymin, Float_t ymax, Int_t status, Int_t type, TString centr, Int_t isSum = 0, TString tag = "ALICE");
  AliParticleYield(Int_t pdg, Int_t system, Float_t sqrts, Float_t value, Float_t stat, Float_t syst, Float_t normPos, Float_t normNeg, Float_t ymin, Float_t ymax, Int_t status, Int_t type, TString centr, Int_t isSum = 0, TString tag = "ALICE");
  virtual ~AliParticleYield();
  AliParticleYield(const AliParticleYield& part); 
  
  // IO
  static TClonesArray * ReadFromASCIIFile(const char * fileName, const char * separators = " \t");
  static TTree * ReadFromASCIIFileAsTree(const char * fileName, const char * separators = " \t");
  static void SaveAsASCIIFile(TClonesArray * arr, const char * filename, const char * separator = " ", Int_t colWidth = 7);  
  static void WriteThermusFile(TClonesArray * arr, const char * filename, Int_t colWidth = 10);
  static TClonesArray * GetEntriesMatchingSelection(TTree * tree, TCut selection); 
  static TTree * GetTreeFromArray(TClonesArray * arr) ;

  // Misc helpers
  Bool_t CheckTypeConsistency() const;
  Bool_t CheckForDuplicates(TClonesArray * arr) ;
  virtual void Print (Option_t * opt = "") const;
  static Float_t GetError(TString error, Float_t yield) ;
  static const char * FormatCol(const char * text, Int_t width,  const char * sep =" ") ;
  static Double_t RoundToSignificantFigures(double num, int n) ;  
  static AliParticleYield * FindParticle(TClonesArray * arr, Int_t pdg, Int_t system, Float_t sqrts, TString centrality = "", Int_t isSum = -1, Int_t status = -1, Int_t pdg2 = 0);
  static AliParticleYield * FindRatio   (TClonesArray * arr, Int_t pdg, Int_t pdg2, Int_t system, Float_t sqrts, TString centrality="", Int_t isSum = -1, Int_t status = -1) { return FindParticle(arr, pdg, system, sqrts, centrality, isSum, status, pdg2); }
  Bool_t operator==(const AliParticleYield& rhs);
  Bool_t IsTheSameMeasurement(AliParticleYield &rhs);
  static AliParticleYield * Divide(AliParticleYield * part1, AliParticleYield * part2, Double_t correlatedError = 0, Option_t * opt="");
  static AliParticleYield * Add   (AliParticleYield * part1, AliParticleYield * part2, Double_t correlatedError = 0, Option_t * opt="");
  void Scale(Float_t scale) ;

  // Getters
  TString GetCentr()           const{ return fCentr           ;}
  Int_t   GetCollisionSystem() const{ return fCollisionSystem ;}
  Int_t   GetIsSum()           const{ return fIsSum           ;}
  Int_t   GetPdgCode()         const{ return fPdgCode         ;}
  Int_t   GetPdgCode2()        const{ return fPdgCode2; }
  Float_t GetSqrtS()           const{ return fSqrtS           ;}
  Float_t GetYield()           const{ return fYield           ;}
  Float_t GetNormError()       const;
  Float_t GetNormErrorNeg()    const{ return fNormErrorNeg; }
  Float_t GetNormErrorPos()    const{ return fNormErrorPos; }

  TString GetPartName()        const{ return fPartName        ;}
  Float_t GetStatError()       const{ return fStatError       ;}
  Int_t   GetStatus()          const{ return fStatus          ;}
  Float_t GetSystError()       const{ return fSystError       ;}
  Float_t GetYMax()            const{ return fYMax            ;}
  Float_t GetYMin()            const{ return fYMin            ;}
  UInt_t  GetMeasurementType() const{ return fMeasurementType ;}
  TString GetTag()             const{ return fTag; }

  const char * GetLatexName(Int_t pdg = 0)  const;
  Float_t GetTotalError(Bool_t includeNormalization = kFALSE) const;

  Bool_t  IsTypeMeasured()     const{ CheckTypeConsistency(); return (!(fMeasurementType & kTypeLinearInterpolation) && !(fMeasurementType & kTypeExtrPionRatio));}
  Bool_t  IsTypeRatio()        const{ CheckTypeConsistency(); return (fMeasurementType & kTypeParticleRatio);}
  Bool_t  IsTypeExtrapolatedWithPionRatio() const { CheckTypeConsistency(); return (fMeasurementType & kTypeExtrPionRatio);}
  Bool_t  IsTypeLinearInterp() const{ CheckTypeConsistency(); return fMeasurementType & kTypeLinearInterpolation;}
  Bool_t  IsTypeOnlyTotErr()   const{ CheckTypeConsistency(); return fMeasurementType & kTypeOnlyTotError;       }
  Bool_t  IsTypeAverage()       const{CheckTypeConsistency(); return fMeasurementType & kTypeAveragePartAntiPart;}

  static Int_t   GetSignificantDigits()  { return fSignificantDigits; }

  // Setters
  void SetCentr           (TString var           ) { fCentr = var           ;}
  void SetCollisionSystem (AliPYCSystem_t var    ) { fCollisionSystem = var ;}
  void SetIsSum           (Int_t var             ) { fIsSum = var           ;}
  void SetPdgCode         (Int_t var             ) { fPdgCode = var ; SetPartName(var);}
  void SetPdgCode2        (Int_t var             ) { fPdgCode2 = var;        }  
  void SetSqrtS           (Float_t var           ) { fSqrtS = var           ;}
  void SetYield           (Float_t var           ) { fYield = var           ;}
  void SetNormError       (Float_t var           ) { fNormErrorPos = var    ; fNormErrorNeg=0;};
  void SetNormErrorNeg    (Float_t var           ) { fNormErrorNeg = var;}
  void SetNormErrorPos    (Float_t var           ) { fNormErrorPos = var;}
  void SetPartName        (TString var           ) { fPartName = var; SetPdgCode(var); }
  void SetStatError       (Float_t var           ) { fStatError = var       ;}
  void SetStatus          (AliPYStatusCode_t var ) { fStatus = var          ;}
  void SetSystError       (Float_t var           ) { fSystError = var       ;}
  void SetYMax            (Float_t var           ) { fYMax = var            ;}
  void SetYMin            (Float_t var           ) { fYMin = var            ;}
  void SetMeasurementType (UInt_t var            ) { fMeasurementType = var ;}
  void SetTag             (TString var           ) { fTag = var;}
  //This 2 additional setters will ensure consistency between the pdg code and the name of the particle
  void SetPartName(Int_t pdgCode);
  void SetPdgCode (TString partName);

  void SetTypeBits(UInt_t mask) { fMeasurementType |= mask; } // This switches on the bits passed. Does not affect the others! If you want to set the Whole mask, use SetMeasurementType

  static void SetSignificantDigits (Int_t var) { fSignificantDigits = var;}
  



private:

  static Bool_t   Compare2Floats(Float_t a, Float_t b) ;
  static Double_t SumErrors(AliParticleYield * part1, AliParticleYield * part2, Int_t error, Option_t * opt) ;
  void CombineMetadata(AliParticleYield *part1, AliParticleYield*part2, const char * pdgSep) ;


  Int_t   fPdgCode;         // PdgCode
  Int_t   fPdgCode2;        // The PdgCode of the second particle, only needed in case of a ratio
  TString fPartName;        // Particle name (redundant, we also have PDG code)
  Int_t   fCollisionSystem; // Collision System, see the AliPYCSystem_t enum for possible values
  Float_t fSqrtS;           // center of mass energy, in GeV
  Float_t fYield;           // The yield
  Float_t fStatError;       // StatError
  Float_t fSystError;       // SystError
  Float_t fNormErrorPos;    // Normalization error, if the error is simmetric, this is used as the symmetric error. Otherwise it is just the positive one
  Float_t fNormErrorNeg;    // Normalization error (negative)
  Float_t fYMin;            // min rapidity cut
  Float_t fYMax;            // max rapidity cut
  Int_t   fStatus;          // Status code, to determine the quality of the measurement, see AliPYStatusCode_t for possible values



  UInt_t  fMeasurementType; // Measurement Type, e.g. actually measured, interpolated from 2 centrality bins  or only total error given, etc. THIS IS A BIT MASK see AliPYMeasurementType_t for possible values and the IsType* Methods for easy access. Be carefull not to set mutually exclusive values

  TString fCentr;           // Centrality. The format is estimator 3-digits id, bin low-edge 2 digits, bin hi-edge 2 digits , e.g. V0A0005
  Int_t   fIsSum;           // A flag which indicates if the yield is for a single charge or for the sum. 0 = single charge, 1 = particle + antiparticle
  TString fTag;             // Generic text tag (to be used e.g. for the name of the experiment)


  static Int_t fSignificantDigits; // Significant Digits to be used in values and errors
  static Float_t fEpsilon; // Used for float conparisons

  ClassDef(AliParticleYield,2)
};


#endif /* _ALIPARTICLEYIELD_H_ */
