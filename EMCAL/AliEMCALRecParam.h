#ifndef ALIEMCALRECPARAM_H
#define ALIEMCALRECPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------------
// Container of EMCAL reconstruction parameters
// The purpose of this object is to store it to OCDB
// and retrieve it in AliEMCALClusterizerv1, AliEMCALPID and AliEMCALTracker
// Author: Yuri Kharlov
//-----------------------------------------------------------------------------

// --- ROOT system ---

#include "TObject.h" 

class AliEMCALRecParam : public TObject
{
public:
  
  AliEMCALRecParam() ;
  virtual ~AliEMCALRecParam() {}
 
 //Clustering
  Float_t GetClusteringThreshold() const     {return fClusteringThreshold;}
  Float_t GetW0                 () const     {return fW0                 ;}
  Float_t GetMinECut            () const     {return fMinECut            ;}
  void SetClusteringThreshold(Float_t thrsh)   {fClusteringThreshold = thrsh;}
  void SetW0                 (Float_t w0)      {fW0 = w0                    ;}
  void SetMinECut            (Float_t minEcut) {fMinECut = minEcut          ;}

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
  /* track matching cut getters */
  Double_t GetTrkCutX() const        {return fTrkCutX;}
  Double_t GetTrkCutY() const        {return fTrkCutY;}
  Double_t GetTrkCutZ() const        {return fTrkCutZ;}
  Double_t GetTrkCutR() const        {return fTrkCutR;}
  Double_t GetTrkCutAlphaMin() const {return fTrkCutAlphaMin;}
  Double_t GetTrkCutAlphaMax() const {return fTrkCutAlphaMax;}
  Double_t GetTrkCutAngle() const    {return fTrkCutAngle;}

  virtual void Print(Option_t * option="") const ;


private:
  //Clustering
  Float_t fClusteringThreshold ; // minimum energy to seed a EC digit in a cluster
  Float_t fW0 ;                  // logarithmic weight for the cluster center of gravity calculation
  Float_t fMinECut;              // Minimum energy for a digit to be a member of a cluster

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

  ClassDef(AliEMCALRecParam,2)   // Reconstruction parameters

} ;

#endif //  ALIEMCALRECPARAM_H

