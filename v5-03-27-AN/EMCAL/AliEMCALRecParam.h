#ifndef ALIEMCALRECPARAM_H
#define ALIEMCALRECPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in AliEMCALClusterizerv1, AliEMCALPID,
// AliEMCALTracker and use it to configure AliEMCALRawUtils
// 
// 
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

// --- ROOT system ---

#include "AliDetectorRecoParam.h" 
#include "AliLog.h"

class AliEMCALRecParam : public AliDetectorRecoParam
{
 public:

  enum AliEMCALClusterizerFlag
  {
    kClusterizerv1  = 0,
    kClusterizerNxN = 1,
    kClusterizerv2  = 2,
    kClusterizerFW  = 3
  };
  
  AliEMCALRecParam() ;
  AliEMCALRecParam(const AliEMCALRecParam& recParam);
  AliEMCALRecParam& operator = (const AliEMCALRecParam& recParam);
  virtual ~AliEMCALRecParam() {}
  
  //Clustering (Unfolding : Cynthia)
  Float_t GetClusteringThreshold() const     {return fClusteringThreshold ;}
  Float_t GetW0                 () const     {return fW0                  ;}
  Float_t GetMinECut            () const     {return fMinECut             ;}
  Float_t GetLocMaxCut          () const     {return fLocMaxCut           ;}
  Float_t GetTimeCut            () const     {return fTimeCut             ;}
  Float_t GetTimeMin            () const     {return fTimeMin             ;}
  Float_t GetTimeMax            () const     {return fTimeMax             ;}
  Bool_t  GetUnfold             () const     {return fUnfold              ;}
  Int_t   GetNRowDiff           () const     {return fNRowDiff            ;}
  Int_t   GetNColDiff           () const     {return fNColDiff            ;}

  void SetClusteringThreshold(Float_t thrsh)     {fClusteringThreshold = thrsh;}
  void SetW0                 (Float_t w0)        {fW0        = w0         ;}
  void SetMinECut            (Float_t ecut)      {fMinECut   = ecut       ;}
  void SetLocMaxCut          (Float_t locMaxCut) {fLocMaxCut = locMaxCut  ;}
  void SetTimeCut            (Float_t t)         {fTimeCut   = t          ;}
  void SetTimeMin            (Float_t t)         {fTimeMin   = t          ;}
  void SetTimeMax            (Float_t t)         {fTimeMax   = t          ;}
  void SetUnfold             (Bool_t unfold)     {fUnfold    = unfold     ;}
  void SetNxM(Int_t rdiff, Int_t cdiff)          {fNRowDiff=rdiff; fNColDiff = cdiff; }

  //PID (Guenole)
  Double_t GetGamma(Int_t i, Int_t j) const       {return fGamma[i][j];} 
  Double_t GetGammaEnergyProb(Int_t i) const      {return fGammaEnergyProb[i];} 
  Double_t GetGamma1to10(Int_t i, Int_t j) const  {return fGamma1to10[i][j];}   // not used
  Double_t GetHadron(Int_t i, Int_t j) const      {return fHadron[i][j];}
  Double_t GetHadron1to10(Int_t i, Int_t j) const {return fHadron1to10[i][j];}   // not used
  Double_t GetHadronEnergyProb(Int_t i) const     {return fHadronEnergyProb[i];}
  Double_t GetPiZero(Int_t i, Int_t j) const      {return fPiZero[i][j];}
  Double_t GetPiZeroEnergyProb(Int_t i) const     {return fPiZeroEnergyProb[i];}
  
  void SetGamma(Int_t i, Int_t j,Double_t param )       {fGamma[i][j]=param;}
  void SetGammaEnergyProb(Int_t i, Double_t param )     {fGammaEnergyProb[i]=param;}
  void SetGamma1to10(Int_t i, Int_t j,Double_t param )  {fGamma1to10[i][j]=param;}
  void SetHadron(Int_t i, Int_t j,Double_t param )      {fHadron[i][j]=param;}
  void SetHadron1to10(Int_t i, Int_t j,Double_t param ) {fHadron1to10[i][j]=param;}
  void SetHadronEnergyProb(Int_t i,Double_t param )     {fHadronEnergyProb[i]=param;}
  void SetPiZero(Int_t i, Int_t j,Double_t param)       {fPiZero[i][j]=param;}
  void SetPiZeroEnergyProb(Int_t i,Double_t param)      {fPiZeroEnergyProb[i]=param;}
  
  //Track Matching (Alberto; Revised by Rongrong)
  /* track matching cut setters */
  void SetMthCutEta(Double_t value)        {fMthCutEta = value;}
  void SetMthCutPhi(Double_t value)        {fMthCutPhi = value;}
  void SetExtrapolateStep(Double_t value)  {fStep = value;}
  void SetTrkCutPt(Double_t value)         {fTrkCutPt = value;}
  void SetTrkCutNITS(Double_t value)       {fTrkCutNITS = value;}
  void SetTrkCutNTPC(Double_t value)       {fTrkCutNTPC = value;}
  /* track matching cut getters */
  Double_t GetMthCutEta() const         {return fMthCutEta;}
  Double_t GetMthCutPhi() const         {return fMthCutPhi;}
  Double_t GetExtrapolateStep() const   {return fStep;}
  Double_t GetTrkCutPt() const          {return fTrkCutPt;}
  Double_t GetTrkCutNITS() const        {return fTrkCutNITS;}
  Double_t GetTrkCutNTPC() const        {return fTrkCutNTPC;}
  
  //Raw signal fitting (Jenn)
  /* raw signal setters */
  void SetHighLowGainFactor(Double_t value) {fHighLowGainFactor = value;}
  void SetOrderParameter(Int_t value)       {fOrderParameter = value;}
  void SetTau(Double_t value)               {fTau = value;}
  void SetNoiseThreshold(Int_t value)       {fNoiseThreshold = value;}
  void SetNPedSamples(Int_t value)          {fNPedSamples = value;} 
  void SetRemoveBadChannels(Bool_t val)     {fRemoveBadChannels=val; }
  void SetFittingAlgorithm(Int_t val)       {fFittingAlgorithm=val; }
  void SetFALTROUsage(Bool_t val)           {fUseFALTRO=val; }
  void SetLEDFit(Bool_t val)                {fFitLEDEvents=val; }

	
  /* raw signal getters */
  Double_t GetHighLowGainFactor() const {return fHighLowGainFactor;}
  Int_t    GetOrderParameter()    const {return fOrderParameter;}
  Double_t GetTau()               const {return fTau;}
  Int_t    GetNoiseThreshold()    const {return fNoiseThreshold;}
  Int_t    GetNPedSamples()       const {return fNPedSamples;}
  Bool_t   GetRemoveBadChannels() const {return fRemoveBadChannels;}
  Int_t    GetFittingAlgorithm()  const {return fFittingAlgorithm; }
  Bool_t   UseFALTRO()            const {return fUseFALTRO; }
  Bool_t   FitLEDEvents()         const {return fFitLEDEvents; }

  //Unfolding (Adam)
  Double_t GetSSPars(Int_t i) const   {return fSSPars[i];}
  Double_t GetPar5(Int_t i) const     {return fPar5[i];}
  Double_t GetPar6(Int_t i) const     {return fPar6[i];}
  void SetSSPars(Int_t i, Double_t param )     {fSSPars[i]=param;}
  void SetPar5(Int_t i, Double_t param )       {fPar5[i]=param;}
  void SetPar6(Int_t i, Double_t param )       {fPar6[i]=param;}

  virtual void Print(Option_t * option="") const;
  
  static AliEMCALRecParam* GetDefaultParameters();
  static AliEMCALRecParam* GetLowFluxParam();
  static AliEMCALRecParam* GetHighFluxParam();
  static AliEMCALRecParam* GetCalibParam();
  static AliEMCALRecParam* GetCosmicParam();

  static const  TObjArray* GetMappings();
  
  void    SetClusterizerFlag(Short_t val) { fClusterizerFlag = val;  }
  Short_t GetClusterizerFlag() const      { return fClusterizerFlag; }
  
 private:
  //Clustering
  Float_t fClusteringThreshold ; // Minimum energy to seed a EC digit in a cluster
  Float_t fW0 ;                  // Logarithmic weight for the cluster center of gravity calculation
  Float_t fMinECut;              // Minimum energy for a digit to be a member of a cluster
  Bool_t  fUnfold;               // Flag to perform cluster unfolding
  Float_t fLocMaxCut;            // Minimum energy difference to consider local maxima in a cluster
  Float_t fTimeCut ;             // Maximum time of digits with respect to EMC cluster max.
  Float_t fTimeMin ;             // Minimum time of digits
  Float_t fTimeMax ;             // Maximum time of digits
  Short_t fClusterizerFlag ;     // Choice of the clusterizer; Default selection (v1) is zero
  Int_t   fNRowDiff;             // NxN: How many neighbors to consider along row (phi)
  Int_t   fNColDiff;             // NxN: How many neighbors to consider along col (eta)

  //PID (Guenole)
  Double_t fGamma[6][6];         // Parameter to Compute PID for photons     
  Double_t fGamma1to10[6][6];    // Parameter to Compute PID not used
  Double_t fHadron[6][6]; 	     // Parameter to Compute PID for hadrons  	 
  Double_t fHadron1to10[6][6]; 	 // Parameter to Compute PID for hadrons between 1 and 10 GeV  	 
  Double_t fHadronEnergyProb[6]; // Parameter to Compute PID for energy ponderation for hadrons  	 
  Double_t fPiZeroEnergyProb[6]; // Parameter to Compute PID for energy ponderation for Pi0  	 
  Double_t fGammaEnergyProb[6];  // Parameter to Compute PID for energy ponderation for gamma  	 
  Double_t fPiZero[6][6];        // Parameter to Compute PID for pi0  	 
  
  
  //Track-Matching (Alberto; Revised by Rongrong)
  Double_t  fMthCutEta;            // eta-difference cut for track matching
  Double_t  fMthCutPhi;            // phi-difference cut for track matching
  Double_t  fStep;                 // Extrapolate length of each step
  Double_t  fTrkCutPt;             // Minimum pT cut on tracks. Needed for Pb-Pb runs
  Double_t  fTrkCutNITS;           // Number of ITS hits for track matching
  Double_t  fTrkCutNTPC;           // Number of TPC hits for track matching
  
  //Raw signal fitting parameters (Jenn)
  Double_t fHighLowGainFactor;     // gain factor to convert between high and low gain
  Int_t    fOrderParameter;        // order parameter for raw signal fit
  Double_t fTau;                   // decay constant for raw signal fit
  Int_t    fNoiseThreshold;        // threshold to consider signal or noise
  Int_t    fNPedSamples;           // number of time samples to use in pedestal calculation
  Bool_t   fRemoveBadChannels;     // select if bad channels are removed before fitting
  Int_t    fFittingAlgorithm;      // select the fitting algorithm
  Bool_t   fUseFALTRO;             // get FALTRO (trigger) and put it on trigger digits.
  Bool_t   fFitLEDEvents;          // fit LED events or not
	
  //Shower shape parameters (Adam)
  Double_t fSSPars[8]; // Unfolding shower shape parameters
  Double_t fPar5[3];   // UF SSPar nr 5
  Double_t fPar6[3];   // UF SSPar nr 6

  static TObjArray* fgkMaps;       // ALTRO mappings for RCU0..RCUX
  
  ClassDef(AliEMCALRecParam,16)     // Reconstruction parameters
};

#endif //  ALIEMCALRECPARAM_H

