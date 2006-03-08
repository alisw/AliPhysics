#ifndef AliFMDInput_H
#define AliFMDInput_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
#include <TObject.h>
#ifndef ROOT_TString
# include <TString.h>
#endif
class AliRunLoader;
class AliLoader;
class AliStack;
class AliRun;
class AliFMD;
class AliFMDHit;
class TString;
class TClonesArray;
class TTree;
class TGeoManager;
class TParticle;


//___________________________________________________________________
class AliFMDInput : public TObject
{
public:
  enum ETrees {
    kHits       = 1,  // Hits
    kKinematics,      // Kinematics (from sim)
    kDigits,          // Digits
    kSDigits,         // Summable digits 
    kHeader,          // Header information 
    kRecPoints,       // Reconstructed points
    kGeometry         // Not really a tree 
  };
  AliFMDInput();
  AliFMDInput(const char* gAliceFile);
  virtual ~AliFMDInput() {}

  virtual void   AddLoad(ETrees tree)     { SETBIT(fTreeMask, tree); }
  virtual void   RemoveLoad(ETrees tree)  { CLRBIT(fTreeMask, tree); }
  virtual Int_t  NEvents() const;

  virtual Bool_t Init();
  virtual Bool_t Begin(Int_t event);
  virtual Bool_t Event() = 0;
  virtual Bool_t End();
  virtual Bool_t Finish() { return kTRUE; }
  virtual Bool_t Run();
protected:
  TString       fGAliceFile; // File name of gAlice file
  AliRunLoader* fLoader;     // Loader of FMD data 
  AliRun*       fRun;        // Run information
  AliStack*     fStack;      // Stack of particles 
  AliLoader*    fFMDLoader;  // Loader of FMD data 
  AliFMD*       fFMD;        // FMD object
  TTree*        fTreeE;      // Header tree 
  TTree*        fTreeH;      // Hits tree
  TTree*        fTreeD;      // Digit tree 
  TTree*        fTreeS;      // SDigit tree 
  TTree*        fTreeR;      // RecPoint tree
  TClonesArray* fArrayE;     // Event info array
  TClonesArray* fArrayH;     // Hit info array
  TClonesArray* fArrayD;     // Digit info array
  TClonesArray* fArrayS;     // SDigit info array
  TClonesArray* fArrayN;     // Mult (single) info array
  TClonesArray* fArrayP;     // Mult (region) info array
  TGeoManager*  fGeoManager; // Geometry manager
  Int_t         fTreeMask;   // Which tree's to load
  Bool_t        fIsInit;
  ClassDef(AliFMDInput,0)  //Hits for detector FMD
};


//____________________________________________________________________
class AliFMDHit;
class AliFMDInputHits : public AliFMDInput 
{
public:
  AliFMDInputHits(const char* file="galice.root") 
    : AliFMDInput(file) { AddLoad(kHits); }
  virtual Bool_t Event();
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle* track) = 0;
  ClassDef(AliFMDInputHits, 0);
};

//____________________________________________________________________
class AliFMDDigit;
class AliFMDInputDigits : public AliFMDInput 
{
public:
  AliFMDInputDigits(const char* file="galice.root")
    : AliFMDInput(file) { AddLoad(kDigits); }
  virtual Bool_t Event();
  virtual Bool_t ProcessDigit(AliFMDDigit* digit) = 0;
  ClassDef(AliFMDInputDigits, 0);
};

//____________________________________________________________________
class AliFMDSDigit;
class AliFMDInputSDigits : public AliFMDInput 
{
public:
  AliFMDInputSDigits(const char* file="galice.root") 
    : AliFMDInput(file) { AddLoad(kSDigits); }
  virtual Bool_t Event();
  virtual Bool_t ProcessSDigit(AliFMDSDigit* sdigit) = 0;
  ClassDef(AliFMDInputSDigits, 0);
};

//____________________________________________________________________
class AliFMDMultStrip;
class AliFMDMultRegion;
class AliFMDInputRecPoints : public AliFMDInput 
{
public:
  AliFMDInputRecPoints(const char* file="galice.root") 
    : AliFMDInput(file) { AddLoad(kRecPoints); }
  virtual Bool_t Event();
  virtual Bool_t ProcessStrip(AliFMDMultStrip* mult) = 0;
  virtual Bool_t ProcessRegion(AliFMDMultRegion* mult) = 0;
  ClassDef(AliFMDInputRecPoints, 0);
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
