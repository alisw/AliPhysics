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
class AliFMDDigit;
class AliFMDSDigit;
class AliFMDRecPoint;
class AliESD;
class AliESDFMD;
class TString;
class TClonesArray;
class TTree;
class TGeoManager;
class TParticle;
class TChain;

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
    kESD,             // Load ESD's
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
  virtual Bool_t Event();
  virtual Bool_t End();
  virtual Bool_t Finish() { return kTRUE; }
  virtual Bool_t Run();

  virtual Bool_t ProcessHits();
  virtual Bool_t ProcessDigits();
  virtual Bool_t ProcessSDigits();
  virtual Bool_t ProcessRecPoints();

  virtual Bool_t ProcessHit(AliFMDHit*, TParticle*)  { return kTRUE; }
  virtual Bool_t ProcessDigit(AliFMDDigit*)          { return kTRUE; }
  virtual Bool_t ProcessSDigit(AliFMDSDigit*)        { return kTRUE; }
  virtual Bool_t ProcessRecPoint(AliFMDRecPoint*)    { return kTRUE; }
  virtual Bool_t ProcessESD(AliESDFMD*)              { return kTRUE; }
  
protected:
  TString       fGAliceFile; // File name of gAlice file
  AliRunLoader* fLoader;     // Loader of FMD data 
  AliRun*       fRun;        // Run information
  AliStack*     fStack;      // Stack of particles 
  AliLoader*    fFMDLoader;  // Loader of FMD data 
  AliFMD*       fFMD;        // FMD object
  AliESD*       fMainESD;    // ESD Object
  AliESDFMD*    fESD;        // FMD ESD data  
  TTree*        fTreeE;      // Header tree 
  TTree*        fTreeH;      // Hits tree
  TTree*        fTreeD;      // Digit tree 
  TTree*        fTreeS;      // SDigit tree 
  TTree*        fTreeR;      // RecPoint tree
  TChain*       fChainE;     // Chain of ESD's
  TClonesArray* fArrayE;     // Event info array
  TClonesArray* fArrayH;     // Hit info array
  TClonesArray* fArrayD;     // Digit info array
  TClonesArray* fArrayS;     // SDigit info array
  TClonesArray* fArrayR;     // Mult (single) info array
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
  ClassDef(AliFMDInputHits, 0);
};

//____________________________________________________________________
class AliFMDDigit;
class AliFMDInputDigits : public AliFMDInput 
{
public:
  AliFMDInputDigits(const char* file="galice.root")
    : AliFMDInput(file) { AddLoad(kDigits); }
  ClassDef(AliFMDInputDigits, 0);
};

//____________________________________________________________________
class AliFMDSDigit;
class AliFMDInputSDigits : public AliFMDInput 
{
public:
  AliFMDInputSDigits(const char* file="galice.root") 
    : AliFMDInput(file) { AddLoad(kSDigits); }
  ClassDef(AliFMDInputSDigits, 0);
};

//____________________________________________________________________
class AliFMDRecPoint;
class AliFMDInputRecPoints : public AliFMDInput 
{
public:
  AliFMDInputRecPoints(const char* file="galice.root") 
    : AliFMDInput(file) { AddLoad(kRecPoints); }
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
