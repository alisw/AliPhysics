#ifndef ALIFMD_H
#define ALIFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */

//____________________________________________________________________
//
//  Manager class for the FMD - Base class.
//  AliFMDv1, AliFMDv0, and AliFMDAlla 
//  provides concrete implementations. 
//  This class is sooooo crowded
//
#ifndef ALIDETECTOR_H 
# include <AliDetector.h>
#endif
#ifndef ROOT_TBranch
# include <TBranch.h>
#endif
class TBranch;
class TClonesArray;
class TBrowser;
class AliDigitizer;
#ifdef USE_PRE_MOVE
class AliFMDSimulator;
#endif
class AliFMDHit;

//____________________________________________________________________
class AliFMD : public AliDetector 
{
public:
  AliFMD();
  AliFMD(const char *name, const char *title);
  AliFMD(const AliFMD& other);
  virtual ~AliFMD(); 
  AliFMD& operator=(const AliFMD& other);
  // Use old implementation
  void UseOld(Bool_t use=kTRUE) { fUseOld = use;  }
  void UseAssembly(Bool_t use=kTRUE) { fUseAssembly = use; }
  void UseDetailed(Bool_t use=kTRUE) { fDetailed = use; }
  
  // GEometry ANd Tracking (GEANT :-)
  virtual void   CreateGeometry();
  virtual void   CreateMaterials(); 
  virtual void   Init();
  virtual void   StepManager() = 0;
  virtual void   FinishEvent();
  
  // Graphics and event display 
  virtual        void   BuildGeometry();
  virtual        void   DrawDetector();
  virtual        Int_t  DistanceToPrimitive(Int_t px, Int_t py);

  // Hit and digit management 
  virtual void          MakeBranch(Option_t *opt=" ");
  virtual void          SetHitsAddressBranch(TBranch *b);
  virtual void          SetTreeAddress();
  virtual TClonesArray* SDigits() { return fSDigits; }        
  virtual void          ResetSDigits();
  virtual void          AddHit(Int_t track, Int_t *vol, Float_t *hits);
  virtual AliFMDHit*    AddHitByFields(Int_t    track, 
				       UShort_t detector, 
				       Char_t   ring, 
				       UShort_t sector, 
				       UShort_t strip, 
				       Float_t  x=0,
				       Float_t  y=0, 
				       Float_t  z=0,
				       Float_t  px=0, 
				       Float_t  py=0, 
				       Float_t  pz=0,
				       Float_t  edep=0,
				       Int_t    pdg=0,
				       Float_t  t=0, 
				       Float_t  len=0, 
				       Bool_t   stopped=kFALSE);
  virtual        void   AddDigit(Int_t *digits, Int_t* notused=0);
  virtual        void   AddDigitByFields(UShort_t detector=0, 
					 Char_t   ring='\0', 
					 UShort_t sector=0, 
					 UShort_t strip=0, 
					 UShort_t count1=0, 
					 Short_t  count2=-1, 
					 Short_t  count3=-1);
  virtual        void   AddSDigit(Int_t *digits);
  virtual        void   AddSDigitByFields(UShort_t detector=0, 
					  Char_t   ring='\0', 
					  UShort_t sector=0, 
					  UShort_t strip=0, 
					  Float_t  edep=0,
					  UShort_t count1=0, 
					  Short_t  count2=-1, 
					  Short_t  count3=-1);

  // Digitisation
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual        void   Hits2Digits();
  virtual        void   Hits2SDigits();

  // Raw data 
  virtual        void   Digits2Raw();

  // Utility
  void   Browse(TBrowser* b);
protected:
  TClonesArray*      HitsArray();
  TClonesArray*      DigitsArray();
  TClonesArray*      SDigitsArray();

  TClonesArray*      fSDigits;              // Summable digits
  Int_t              fNsdigits;             // Number of digits  
  Bool_t             fDetailed;             // Use detailed geometry
  Bool_t             fUseOld;               // Use old approx geometry
  Bool_t             fUseAssembly;          // Use divided volumes
  
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

#ifdef USE_PRE_MOVE
  AliFMDSimulator*   fSimulator;            // Simulator task
#endif
  TObjArray*         fBad;                  //! debugging - bad hits 
  
  ClassDef(AliFMD,11)     // Base class FMD entry point
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
