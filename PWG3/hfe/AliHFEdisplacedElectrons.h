/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Class for electrons from beauty study
// Counting electrons from beauty
// by DCA cuts, background subtraction 
//
// Authors:
//  Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//  Carlo Bombonati <Carlo.Bombonati@cern.ch>
// 

#ifndef ALIHFEDISPLACEDELECTRONS_H
#define ALIHFEDISPLACEDELECTRONS_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ROOT_TPDGCode
#include <TPDGCode.h>
#endif

class TChain;
class TTree;
class TFile;

class TPDGCode;

class TString;
class TList;

class AliLog;

class THnSparse;

class TObjArray;
class AliStack;
class AliMCEvent;
class AliESDEvent;
class AliVEvent;

class AliESDtrack;
class AliESDVertex;

class AliHFEdisplacedElectrons : public TObject{

 public:  

  enum{
    kPDGelectron = kElectron,
    kPDGgamma = kGamma,
    kPDGpi0 = kPi0,
    kPDGpion = kPiPlus, 
    kPDGeta = 221,
    kPDGcharm = kCharm,
    kPDGbeauty = kBottom
  };  // PDG codes to be used

  AliHFEdisplacedElectrons(); // default constructor
  AliHFEdisplacedElectrons(const AliHFEdisplacedElectrons &p); // copy constructor
  AliHFEdisplacedElectrons &operator=(const AliHFEdisplacedElectrons &ref); // assignment operator

  virtual ~AliHFEdisplacedElectrons();

  void InitAnalysis();  
  void CreateOutputs(TList* const displacedList);

  void FillMcOutput(const AliESDEvent * const fESD, AliMCEvent * const fMC, const AliMCParticle * const mctrack);
  void FillEsdOutput(const AliESDEvent * const fESDEvent, AliESDtrack * const track, AliStack *stack);
  void FillDataOutput(const AliESDEvent * const fESDEvent, AliESDtrack * const track);

  Int_t GetMCpid(AliStack* stack, Int_t label) const;

  Bool_t HasMCData() const { return TestBit(kHasMCData); };
  void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };
  void SetDebugLevel(Int_t debugLevel){ fDeDebugLevel = debugLevel; };
  void SetNitsCluster(Int_t nITScls){ fNclustersITS = nITScls;};
  void SetMinPrimVtxContrib(Int_t nContrib){fMinNprimVtxContributor = nContrib;};

  void PostAnalysis() const;


 private:

  enum{
    kHasMCData = BIT(15),             // bitset for mc data usage
    kHasESDData = BIT(16)
  };
   


  enum{
    kElePhotonConv = 0,
    kEleDirectPhotonConv = 1,
    kElePi0 = 2,
    kEleEta = 3, 
    kEleB = 4, 
    kEleC = 5, 
    kEleBC = 6,
    kEleMissID = 7,
    kEleMissIDpion = 8,
    kPion = 8
  };  // electron source index
  
  enum{
    kMcElectron = 0, 
    kEsdElectron = 1, 
    kDataElectron = 2
  };  // MC or Data

  enum{
    kNDcaMin = 42, 
    kNPtIntv = 14, 
    kNKineVar = 3
  };   // several constant to be used
 

  Int_t CheckCharm(AliStack *const stack, Int_t eleLabel);
  Bool_t IsB(Int_t pdg) const;
  Bool_t IsC(Int_t pdg) const;
 
  Int_t ElectronFromSource(AliStack *stack, Int_t eleLabel) const;
  Int_t ElePhotonDirect(AliStack *stack, Int_t label) const;
  Int_t ElectronFromCharm(AliStack *stack, Int_t eleLabel) const;
  Int_t CharmFromBeauty(AliStack *stack, Int_t charmLabel) const;
  
  Int_t GetMotherLabel(AliStack *stack, Int_t label) const;
  Float_t GetRapidity(TParticle *part) const;
  Float_t GetTrackRapidity(AliESDtrack *track) const;

  static const Float_t fgkDcaMinIntv[kNDcaMin];  // DCA cut min limit
  static const Float_t fgkDcaMinPtIntv[kNPtIntv-1]; // DCA cut min limit in different pT bins
  static const Float_t fgkPtIntv[kNPtIntv]; // all pt bins

  static const Char_t *fgkKineVar[kNKineVar];  // particle names
  static const Char_t *fgkKineVarTitle[kNKineVar];  // particle names
  
  UInt_t fDeDebugLevel;   // debug level
  Int_t fNclustersITS;  // ITS clusters
  Int_t fMinNprimVtxContributor;      // minimum number of contributors to the primary vtx

  THnSparseF *fTHnSparseDcaMcEleInfo;   //! container for MC pion part
  THnSparseF *fTHnSparseDcaEsdEleInfo;   //! container for MC electron part
  THnSparseF *fTHnSparseDcaDataEleInfo; //! container for Data electron part

  TList *fDeOutputList;  //! output container
  ClassDef(AliHFEdisplacedElectrons, 0);
};

#endif
