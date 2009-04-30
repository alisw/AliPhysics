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
  
  AliEMCALRecParam() ;
  AliEMCALRecParam(const AliEMCALRecParam& recParam);
  AliEMCALRecParam& operator = (const AliEMCALRecParam& recParam);
  virtual ~AliEMCALRecParam() {}
 
  //Clustering (Unfolding : Cynthia)
  Float_t GetClusteringThreshold() const     {return fClusteringThreshold;}
  Float_t GetW0                 () const     {return fW0                 ;}
  Float_t GetMinECut            () const     {return fMinECut            ;}
  Float_t GetLocMaxCut          () const     {return fLocMaxCut            ;}
  Bool_t  GetUnfold             () const     {return fUnfold            ;}
  void SetClusteringThreshold(Float_t thrsh)   {fClusteringThreshold = thrsh;}
  void SetW0                 (Float_t w0)      {fW0 = w0                    ;}
  void SetMinECut            (Float_t minEcut) {fMinECut = minEcut          ;}
  void SetLocMaxCut          (Float_t locMaxCut) {fLocMaxCut = locMaxCut    ;}
  void SetUnfold             (Bool_t unfold)     {fUnfold = unfold          ; if(fUnfold) AliWarning("Cluster Unfolding ON. Implementing only for eta=0 case!!!");}

  //PID (Guenole)
  Double_t GetGamma(Int_t i, Int_t j) const    {return fGamma[i][j];} 
  Double_t GetHadron(Int_t i, Int_t j) const    {return fHadron[i][j];}
  Double_t GetPiZero5to10(Int_t i, Int_t j) const    {return fPiZero5to10[i][j];}
  Double_t GetPiZero10to60(Int_t i, Int_t j) const    {return fPiZero10to60[i][j];}

  void SetGamma(Int_t i, Int_t j,Double_t param )   {fGamma[i][j]=param;}
  void SetHadron(Int_t i, Int_t j,Double_t param )   {fHadron[i][j]=param;}
  void SetPiZero5to10(Int_t i, Int_t j,Double_t param)   {fPiZero5to10[i][j]=param;}
  void SetPiZero10to60(Int_t i, Int_t j,Double_t param)   {fPiZero10to60[i][j]=param;}

  //Track Matching (Alberto)
  /* track matching cut setters */
  void SetTrkCutX(Double_t value)        {fTrkCutX = value;}
  void SetTrkCutY(Double_t value)        {fTrkCutY = value;}
  void SetTrkCutZ(Double_t value)        {fTrkCutZ = value;}
  void SetTrkCutR(Double_t value)        {fTrkCutR = value;}
  void SetTrkCutAlphaMin(Double_t value) {fTrkCutAlphaMin = value;}
  void SetTrkCutAlphaMax(Double_t value) {fTrkCutAlphaMax = value;}
  void SetTrkCutAngle(Double_t value)    {fTrkCutAngle = value;}
  void SetTrkCutNITS(Double_t value)        {fTrkCutNITS = value;}
  void SetTrkCutNTPC(Double_t value)        {fTrkCutNTPC = value;}
  /* track matching cut getters */
  Double_t GetTrkCutX() const        {return fTrkCutX;}
  Double_t GetTrkCutY() const        {return fTrkCutY;}
  Double_t GetTrkCutZ() const        {return fTrkCutZ;}
  Double_t GetTrkCutR() const        {return fTrkCutR;}
  Double_t GetTrkCutAlphaMin() const {return fTrkCutAlphaMin;}
  Double_t GetTrkCutAlphaMax() const {return fTrkCutAlphaMax;}
  Double_t GetTrkCutAngle() const    {return fTrkCutAngle;}
  Double_t GetTrkCutNITS() const        {return fTrkCutNITS;}
  Double_t GetTrkCutNTPC() const        {return fTrkCutNTPC;}

  //Raw signal fitting (Jenn)
  /* raw signal setters */
  void SetHighLowGainFactor(Double_t value) {fHighLowGainFactor = value;}
  void SetOrderParameter(Int_t value)       {fOrderParameter = value;}
  void SetTau(Double_t value)               {fTau = value;}
  void SetNoiseThreshold(Int_t value)       {fNoiseThreshold = value;}
  void SetNPedSamples(Int_t value)          {fNPedSamples = value;}
  /* raw signal getters */
  Double_t GetHighLowGainFactor() const {return fHighLowGainFactor;}
  Int_t    GetOrderParameter()    const {return fOrderParameter;}
  Double_t GetTau()               const {return fTau;}
  Int_t    GetNoiseThreshold()    const {return fNoiseThreshold;}
  Int_t    GetNPedSamples()       const {return fNPedSamples;}
 
  virtual void Print(Option_t * option="") const ;

  static AliEMCALRecParam* GetDefaultParameters();
  static AliEMCALRecParam* GetLowFluxParam();
  static AliEMCALRecParam* GetHighFluxParam();

  static const  TObjArray* GetMappings();

private:
  //Clustering
  Float_t fClusteringThreshold ; // minimum energy to seed a EC digit in a cluster
  Float_t fW0 ;                  // logarithmic weight for the cluster center of gravity calculation
  Float_t fMinECut;              // Minimum energy for a digit to be a member of a cluster
  Bool_t fUnfold;               // flag to perform cluster unfolding
  Float_t fLocMaxCut;            // minimum energy difference to consider local maxima in a cluster

  //PID (Guenole)
  Double_t fGamma[6][6];        // Parameter to Compute PID      
  Double_t fHadron[6][6]; 	// Parameter to Compute PID   	 
  Double_t fPiZero5to10[6][6];  // Parameter to Compute PID   	 
  Double_t fPiZero10to60[6][6]; // Parameter to Compute PID   	 

  //Track-Matching (Alberto)
  Double_t  fTrkCutX;              // X-difference cut for track matching
  Double_t  fTrkCutY;              // Y-difference cut for track matching
  Double_t  fTrkCutZ;              // Z-difference cut for track matching
  Double_t  fTrkCutR;              // cut on allowed track-cluster distance
  Double_t  fTrkCutAlphaMin;       // cut on 'alpha' parameter for track matching (min)
  Double_t  fTrkCutAlphaMax;       // cut on 'alpha' parameter for track matching (min)
  Double_t  fTrkCutAngle;          // cut on relative angle between different track points for track matching
  Double_t  fTrkCutNITS;              // Number of ITS hits for track matching
  Double_t  fTrkCutNTPC;              // Number of TPC hits for track matching
 
  //Raw signal fitting parameters (Jenn)
  Double_t fHighLowGainFactor;     //gain factor to convert between high and low gain
  Int_t    fOrderParameter;        //order parameter for raw signal fit
  Double_t fTau;                   //decay constant for raw signal fit
  Int_t    fNoiseThreshold;        //threshold to consider signal or noise
  Int_t    fNPedSamples;           //number of time samples to use in pedestal calculation

  static TObjArray* fgkMaps;       // ALTRO mappings for RCU0..RCUX

  ClassDef(AliEMCALRecParam,6)   // Reconstruction parameters

} ;

#endif //  ALIEMCALRECPARAM_H

