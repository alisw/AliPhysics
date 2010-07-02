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
#include "TH3F.h"
#include "AliFMDAnaCalibBackgroundCorrection.h"
#include "AliFMDAnaCalibEnergyDistribution.h"
#include "AliFMDAnaCalibEventSelectionEfficiency.h"
#include "AliFMDAnaCalibSharingEfficiency.h"
#include "AliPhysicsSelection.h"
#include <TVector2.h>
#include <TString.h>
//#include "AliPWG0Helper.h"
#include "AliESDEvent.h"
//class AliESDEvent;

/**
 * @ingroup FMD_ana
 */
class AliFMDAnaParameters : public TNamed
{
public:
  /** Enumeration of things to initialize */ 
  enum What { 
    /** Pulser gain */ 
    kBackgroundCorrection         = 0x1, // Background Correction 
    kEnergyDistributions          = 0x2, // Energy Distributions
    kEventSelectionEfficiency     = 0x4, // Event Selection Efficiency
    kSharingEfficiency            = 0x8  // Sharing algorithm efficiency
  };
  
  enum Trigger { kMB1 = 0, kMB2, kSPDFASTOR, kNOCTP, kEMPTY };
  
  enum Energy { k900 , k10000, k14000 , k7000, k2400, k5500};
  
  enum MagField {k0G, k5G};
  
  enum Species {kPP, kPbPb};

  /** DO NOT USE THIS - ONLY FOR IO */
  AliFMDAnaParameters();
  
  /** Singleton access
      @return  single to */
  static AliFMDAnaParameters* Instance();
  
  void Init(Bool_t forceReInit=kTRUE, UInt_t what=kBackgroundCorrection|kEnergyDistributions|kEventSelectionEfficiency|kSharingEfficiency);
  Float_t GetVtxCutZ();
  Int_t GetNvtxBins();
  Int_t GetNetaBins();
  Float_t GetEtaMin();  
  Float_t GetEtaMax();
  Int_t GetEtaBin(Float_t eta);
  Float_t GetMPV(Int_t det, Char_t ring, Float_t eta);
  Float_t GetConstant(Int_t det, Char_t ring, Float_t eta);
  Float_t GetSigma(Int_t det, Char_t ring, Float_t eta);
  Float_t Get2MIPWeight(Int_t det, Char_t ring, Float_t eta);
  Float_t Get3MIPWeight(Int_t det, Char_t ring, Float_t eta);
  //static const char* GetBackgroundPath() { return fgkBackgroundCorrection;}
  // static const char* GetEdistPath()      { return fgkEnergyDists;}
  static const char* GetBackgroundID() { return fgkBackgroundID;}
  static const char* GetEdistID()      { return fgkEnergyDistributionID;}
  static const char* GetEventSelectionEffID()      { return fgkEventSelectionEffID;}
  static const char* GetSharingEffID()      { return fgkSharingEffID;}
  TH2F* GetBackgroundCorrection(Int_t det, Char_t ring, Int_t vtxbin);
  TH1F* GetDoubleHitCorrection(Int_t det, Char_t ring);
  TH1F* GetSPDDeadCorrection(Int_t vtxbin);
  
  TH1F* GetSharingEfficiency(Int_t det, Char_t ring, Int_t vtxbin);
  TH1F* GetSharingEfficiencyTrVtx(Int_t det, Char_t ring, Int_t vtxbin);
  Float_t  GetEventSelectionEfficiency(Int_t vtxbin);
  TH2F*    GetEventSelectionEfficiency(Int_t vtxbin, Char_t ring);
  Float_t  GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec) const;
  Float_t  GetEtaFromStrip(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip, Float_t zvtx) const;
  Float_t  GetStripLength(Char_t ring, UShort_t strip)  ;
  Float_t  GetBaseStripLength(Char_t ring, UShort_t strip)  ;
  Float_t  GetMaxR(Char_t ring) const;
  Float_t  GetMinR(Char_t ring) const;
  void     SetBackgroundPath(const Char_t* bgpath) {fBackgroundPath.Form(bgpath);}
  void     SetEnergyPath(const Char_t* epath) {fEnergyPath.Form(epath);}
  void     SetEventSelectionPath(const Char_t* evpath) {fEventSelectionEffPath.Form(evpath);}
  void     SetSharingEfficiencyPath(const Char_t* sharpath) {fSharingEffPath.Form(sharpath);}
  void     SetProcessPrimary(Bool_t prim=kTRUE) {fProcessPrimary = prim;}
  void     SetProcessHits(Bool_t hits=kTRUE) {fProcessHits = hits;}
  Bool_t   GetProcessPrimary() const {return fProcessPrimary;} 
  Bool_t   GetProcessHits() const {return fProcessHits;}
  Bool_t   GetVertex(const AliESDEvent* esd, Double_t* vertexXYZ);
  void     SetTriggerDefinition(Trigger trigger) {fTrigger = trigger;}
  Trigger  GetTriggerDefinition() const {return fTrigger;}
  Bool_t   IsEventTriggered(const AliESDEvent *esd) const;
  Bool_t   IsEventTriggered(const AliESDEvent *esd, Trigger trigger) ;
  void     SetEnergy(Energy energy) {fEnergy = energy;}
  Energy   GetEnergy() {return fEnergy;}
  void     SetMagField(MagField magfield) {fMagField = magfield;}
  char*    GetPath(const char* species);
  void     SetCollisionSystem(Species collsystem) {fSpecies = collsystem;}
  Species  GetCollisionSystem() const {return fSpecies;}
  void     PrintStatus() const;
  void     Print(Option_t* /* option */) const { PrintStatus(); }
  char*    GetDndetaAnalysisName() const {return "PWG2forwardDnDeta";}
  TH1F*    GetEnergyDistribution(Int_t det, Char_t ring, Float_t eta);
  TH1F*    GetEmptyEnergyDistribution(Int_t det, Char_t ring);
  TH1F*    GetRingEnergyDistribution(Int_t det, Char_t ring);
  AliPhysicsSelection* GetPhysicsSelection() const {return fPhysicsSelection;}
  Bool_t   IsRealData() const {return fRealData; }
  void     SetRealData(Bool_t realdata) {fRealData = realdata;}
  Float_t  GetLowSPDLimit() const {return fSPDlowLimit;}
  Float_t  GetHighSPDLimit() const {return fSPDhighLimit;}
  void     SetLowSPDLimit(Float_t cut) {fSPDlowLimit = cut;}
  void     SetHighSPDLimit(Float_t cut) {fSPDhighLimit = cut;}
  void     SetCentralTriggerSelection(Bool_t selection) {fCentralSelection = selection;}
  Bool_t   SharingEffPresent() {return fSharingObjectPresent;}
  Int_t    GetFirstEtaBinToInclude(Int_t vtxbin, Int_t det, Char_t ring) ;
  Int_t    GetLastEtaBinToInclude(Int_t vtxbin, Int_t det, Char_t ring) ;
  

  void     SetNumberOfEtaBinsToCut(Int_t nbins) {fNumberOfEtaBinsToCut = nbins;}
  Int_t    GetNumberOfEtaBinsToCut() const {return fNumberOfEtaBinsToCut;}
  
protected:
  
  AliFMDAnaParameters(const AliFMDAnaParameters& o) 
    : TNamed(o),
      fIsInit(o.fIsInit),
      fBackground(o.fBackground),
      fEnergyDistribution(o.fEnergyDistribution),
      fEventSelectionEfficiency(o.fEventSelectionEfficiency),
      fSharingEfficiency(o.fSharingEfficiency),
      fCorner1(o.fCorner1),
      fCorner2(o.fCorner2),
      fEnergyPath(o.fEnergyPath),
      fBackgroundPath(o.fBackgroundPath),
      fEventSelectionEffPath(o.fEventSelectionEffPath),
      fSharingEffPath(o.fSharingEffPath),
      fProcessPrimary(o.fProcessPrimary),
      fProcessHits(o.fProcessHits),
      fTrigger(o.fTrigger),
      fEnergy(o.fEnergy),
      fMagField(o.fMagField),
      fSpecies(o.fSpecies),
      fPhysicsSelection(o.fPhysicsSelection), 
      fRealData(o.fRealData),
      fSPDlowLimit(o.fSPDlowLimit),
      fSPDhighLimit(o.fSPDhighLimit),   
      fCentralSelection(o.fCentralSelection),
      fSharingObjectPresent(o.fSharingObjectPresent),
      fNumberOfEtaBinsToCut(o.fNumberOfEtaBinsToCut),
      fEtaLowBinLimits(o.fEtaLowBinLimits),
      fEtaHighBinLimits(o.fEtaHighBinLimits)
  {}
  AliFMDAnaParameters& operator=(const AliFMDAnaParameters&) { return *this; }
  virtual ~AliFMDAnaParameters() {}
  
  static AliFMDAnaParameters* fgInstance;   // Static singleton instance
  
  //  AliCDBEntry* GetEntry(const char* path, Bool_t fatal=kTRUE) const ;
  void InitBackground();
  void InitEnergyDists();
  void InitEventSelectionEff();
  void InitSharingEff();
  
  void     FindEtaLimits();
  Int_t    GetFirstEtaBinFromMap(Int_t vtxbin, Int_t det, Char_t ring) ;
  Int_t    GetLastEtaBinFromMap(Int_t vtxbin, Int_t det, Char_t ring) ;

  
  TObjArray* GetBackgroundArray();
  
  TAxis* GetRefAxis();
  void SetCorners(Char_t ring) ;
  
  Bool_t fIsInit;                      //Have we been init ?
  //TObjArray*  fBackgroundArray;
  // TObjArray*  fEdistArray;
  AliFMDAnaCalibBackgroundCorrection*         fBackground;  //BG correction object 
  AliFMDAnaCalibEnergyDistribution*           fEnergyDistribution; //Energy dist object
  AliFMDAnaCalibEventSelectionEfficiency*     fEventSelectionEfficiency; //Event selection correction object
  AliFMDAnaCalibSharingEfficiency*            fSharingEfficiency; //Sharing efficiency correction object
  //static const char* fgkBackgroundCorrection;
  //static const char* fgkEnergyDists;
  static const char* fgkBackgroundID;           // BG cor ID
  static const char* fgkEnergyDistributionID ;  // Energy dist ID
  static const char* fgkEventSelectionEffID ;   // Event selection ID
  static const char* fgkSharingEffID ;          // Sharing efficiency ID  
  
  TVector2 fCorner1;                  //First corner of hybrid
  TVector2 fCorner2;                  //Second corner of hybrid
  TString  fEnergyPath;               //Path of energy calib
  TString  fBackgroundPath;           //Path of BG correction
  TString  fEventSelectionEffPath;    //Path of event selection eff
  TString  fSharingEffPath;           //Path of sharing eff
  Bool_t   fProcessPrimary;           //Do we process primary ?
  Bool_t   fProcessHits;              //Do we process hits ?
  Trigger  fTrigger;                  //Which trigger are we using ?
  Energy   fEnergy;                   // CM energy
  MagField fMagField;                 //Magnetic field
  Species  fSpecies;                  //PbPb or pp ?
  AliPhysicsSelection* fPhysicsSelection;  //Physics selection
  Bool_t   fRealData;                 // real or simulated
  Float_t  fSPDlowLimit ;             // low limit of SPD tracklets
  Float_t  fSPDhighLimit ;             // high limit of SPD tracklets
  Bool_t   fCentralSelection;         //if event selection is done centrally
  Bool_t   fSharingObjectPresent ;    //Do we have a sharing object ? 
  Int_t    fNumberOfEtaBinsToCut;     //Number of eta bins to remove from edge effects
  TH3F     fEtaLowBinLimits;
  TH3F     fEtaHighBinLimits;
  ClassDef(AliFMDAnaParameters,1) // Manager of parameters
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

