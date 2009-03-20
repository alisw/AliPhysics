#ifndef ALIFMDANAPARAMETERS_H
#define ALIFMDANAPARAMETERS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
 *
 * See cxx source for full Copyright notice                               
 */
//
//The design of this class is based on the AliFMDParameters class. Its purpose
//is to hold parameters for the analysis such as background correction and 
//fit functions.
//
//Author: Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
//
//____________________________________________________________________

#ifndef ROOT_TNamed
# include <TNamed.h>
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif
#ifndef ALIFMDUSHORTMAP_H
# include <AliFMDUShortMap.h>
#endif
#ifndef ALIFMDBOOLMAP_H
# include <AliFMDBoolMap.h>
#endif
#include "AliCDBEntry.h"

#include "TFile.h"
#include "TObjArray.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TH1F.h"
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include "AliFMDAnaCalibEnergyDistribution.h"
//____________________________________________________________________
//
//  Singleton class to handle various parameters (not geometry) of the
//  FMD
//  Should get ata fromm Conditions DB.
//

class AliFMDAnaParameters : public TNamed
{
public:
  /** Enumeration of things to initialize */ 
  enum What { 
    /** Pulser gain */ 
    kBackgroundCorrection = 0x1, // Background Correction 
    kEnergyDistributions  = 0x2 // Energy Distributions
  };
  
  /** Singleton access
      @return  single to */
  static AliFMDAnaParameters* Instance();
  
  void Init(Bool_t forceReInit=kTRUE, UInt_t what=kBackgroundCorrection|kEnergyDistributions);
  Float_t GetVtxCutZ();
  Int_t GetNvtxBins();
  Int_t GetNetaBins();
  Float_t GetEtaMin();  
  Float_t GetEtaMax();
  Float_t GetMPV(Int_t det, Char_t ring, Float_t eta);
  Float_t GetSigma(Int_t det, Char_t ring, Float_t eta);
  Float_t Get2MIPWeight(Int_t det, Char_t ring, Float_t eta);
  Float_t Get3MIPWeight(Int_t det, Char_t ring, Float_t eta);
  static const char* GetBackgroundPath() { return fgkBackgroundCorrection;}
  static const char* GetEdistPath()      { return fgkEnergyDists;}
  TH2F* GetBackgroundCorrection(Int_t det, Char_t ring, Int_t vtxbin);
  
protected:
  
  AliFMDAnaParameters();
  
  AliFMDAnaParameters(const AliFMDAnaParameters& o) 
    : TNamed(o),
      fIsInit(o.fIsInit),
      fBackground(o.fBackground),
      fEnergyDistribution(o.fEnergyDistribution)
      //  fBackgroundArray(o.fBackgroundArray), 
      //fEdistArray(o.fEdistArray)
  {}
  AliFMDAnaParameters& operator=(const AliFMDAnaParameters&) { return *this; }
  virtual ~AliFMDAnaParameters() {}
  
  static AliFMDAnaParameters* fgInstance;   // Static singleton instance
  
  AliCDBEntry* GetEntry(const char* path, Bool_t fatal=kTRUE) const ;
  void InitBackground();
  void InitEnergyDists();
  TH1F* GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta);
  TObjArray* GetBackgroundArray();
  
  TAxis* GetRefAxis();
  
  
  Bool_t fIsInit;
  //TObjArray*  fBackgroundArray;
  // TObjArray*  fEdistArray;
  AliFMDAnaCalibBackgroundCorrection* fBackground;
  AliFMDAnaCalibEnergyDistribution* fEnergyDistribution;
  static const char* fgkBackgroundCorrection;
  static const char* fgkEnergyDists;
  ClassDef(AliFMDAnaParameters,0) // Manager of parameters
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

