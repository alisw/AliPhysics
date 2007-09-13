#ifndef ALIJETUNITARRAY_H
#define ALIJETUNITARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

//------------------------------------------------------------------
//  Unit used by UA1 algorithm
//  Authors: Sarah Blyth (LBL/UCT)
//           Magali Estienne (IReS) (new version for JETAN)
//------------------------------------------------------------------


#include <TObject.h>
#include "AliJetFinderTypes.h"

class AliJetUnitArray : public TObject
{
 public:
  AliJetUnitArray();
  AliJetUnitArray(Int_t absId, Int_t esdId, Float_t eta, Float_t phi, Float_t en, Float_t px, Float_t py, Float_t pz, Float_t Deta, Float_t Dphi, AliJetFinderUnitDetectorFlagType_t det, AliJetFinderUnitFlagType_t inout, AliJetFinderUnitCutFlagType_t cut, Float_t mass, Int_t clusId);
  ~AliJetUnitArray();

  // Setter
  void SetUnitEnergy(Float_t energy)  {fUnitEnergy = energy;}
  void SetUnitEta(Float_t eta)        {fUnitEta = eta;} 
  void SetUnitPhi(Float_t phi)        {fUnitPhi = phi;}
  void SetUnitDeta(Float_t deta)      {fUnitDeta = deta;} 
  void SetUnitDphi(Float_t dphi)      {fUnitDphi = dphi;}
  void SetUnitID(Int_t id)            {fUnitID = id;}
  void SetUnitTrackID(Int_t esdid)    {fUnitTrackID = esdid;}
  void SetUnitEntries(Int_t num)      {fUnitNum = num;}
  void SetUnitClusterID(Int_t id)     {fUnitClusterID = id;}
  void SetUnitFlag(AliJetFinderUnitFlagType_t flag)    
  { 
    fUnitFlag = flag; 
  } 
  void SetUnitCutFlag(AliJetFinderUnitCutFlagType_t cutFlag)    
  { 
    fUnitCutFlag = cutFlag; 
  } 
  void SetUnitSignalFlag(AliJetFinderUnitSignalFlagType_t signalFlag)    
  { 
    fUnitSignalFlag = signalFlag; 
  } 
  void SetUnitDetectorFlag(AliJetFinderUnitDetectorFlagType_t detectorflag)    
  { 
    fUnitDetectorFlag = detectorflag; 
  } 
  void SetUnitPxPyPz(Double_t *pxyz)  {fUnitPx = pxyz[0]; fUnitPy = pxyz[1]; fUnitPz = pxyz[2];}
  void SetUnitMass(Float_t mass) {fUnitMass = mass;}

  // Getter
  Float_t GetUnitEnergy() const            {return fUnitEnergy;}
  Float_t GetUnitEta() const               {return fUnitEta;}
  Float_t GetUnitPhi() const               {return fUnitPhi;}         
  Float_t GetUnitDeta() const              {return fUnitDeta;}
  Float_t GetUnitDphi() const              {return fUnitDphi;}         
  Int_t   GetUnitID() const                {return fUnitID;}
  Int_t   GetUnitTrackID() const             {return fUnitTrackID;}
  Int_t   GetUnitEntries() const           {return fUnitNum;}
  Int_t   GetUnitClusterID() const         {return fUnitClusterID;}
  Float_t GetUnitMass() const              {return fUnitMass;}
  Bool_t  GetUnitPxPyPz(Double_t* p) const {p[0]=fUnitPx; p[1]=fUnitPy; p[2]=fUnitPz; return kTRUE;}

  AliJetFinderUnitFlagType_t GetUnitFlag() const     
  {
	  return fUnitFlag;
  }
  AliJetFinderUnitCutFlagType_t GetUnitCutFlag() const     
  {
	  return fUnitCutFlag;
  }
  AliJetFinderUnitSignalFlagType_t GetUnitSignalFlag() const     
  {
	  return fUnitSignalFlag;
  }
  AliJetFinderUnitDetectorFlagType_t GetUnitDetectorFlag() const     
  {
	  return fUnitDetectorFlag;
  }


  Bool_t operator>  ( AliJetUnitArray unit1) const;
  Bool_t operator<  ( AliJetUnitArray unit1) const;
  Bool_t operator== ( AliJetUnitArray unit1) const;

 protected:
  Float_t         fUnitEnergy;                          // Energy of the unit 
  Float_t         fUnitEta;                             // Eta of the unit
  Float_t         fUnitPhi;                             // Phi of the unit
  Float_t         fUnitDeta;                            // Delta Eta of the unit
  Float_t         fUnitDphi;                            // Delta Phi of the unit
  Int_t           fUnitID;                              // ID of the unit
  Int_t           fUnitTrackID;                         // ID of the unit
  Int_t           fUnitNum;                             // number of units
  Int_t           fUnitClusterID;                       // ID of the unit
  AliJetFinderUnitFlagType_t         fUnitFlag;         // Flag of the unit
  AliJetFinderUnitCutFlagType_t      fUnitCutFlag;      // Flag of the unit
  AliJetFinderUnitSignalFlagType_t   fUnitSignalFlag;   // Flag of the unit
  AliJetFinderUnitDetectorFlagType_t fUnitDetectorFlag; // Detector flag of the unit
  Float_t         fUnitPx;                              // Px of charged track
  Float_t         fUnitPy;                              // Py of charged track
  Float_t         fUnitPz;                              // Pz of charged track
  Float_t         fUnitMass;                            // Mass of charged particle

  ClassDef(AliJetUnitArray,1)
};

#endif

