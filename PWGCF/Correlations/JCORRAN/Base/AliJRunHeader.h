/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJRunHeader.h,v 1.1 2008/02/04 13:28:47 rak Exp $
////////////////////////////////////////////////////
/*!
  \file AliJRunHeader.h
  \brief
  \author J. Rak, D.J.Kim, F.Krizek  (Jyvaskyla || HIP)
  \email: djkim@cc.jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/02/04 13:28:47 $
  */
////////////////////////////////////////////////////

#ifndef ALIJRUNHEADER_H
#define ALIJRUNHEADER_H
#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <TNamed.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <vector>
#include <map>


class AliJRunHeader : public TNamed {

//#pragma link C++ class AliJRunHeader+;

 public:
  enum { kJPP, kJPbPb, kPPb };
  enum { kJESD, kJAOD };  
  AliJRunHeader();//constructor
  ~AliJRunHeader(){;}              //destructor

  //==== General Information ====//
  Int_t  GetRunNumber()     const {return fRunNumber;}
  void   SetRunNumber(Int_t runN) { fRunNumber = runN;}

  TString GetRunType() const { return fRunType;}
  void    SetRunType( const TString type ){ fRunType = type; }

  TString  GetBeamType()const{ return fBeamType; }
  void    SetBeamType(TString t){ fBeamType = t; }

  Int_t   GetBeamTypeI()const{ return fBeamTypeI; }
  void    SetBeamTypeI(int t){ fBeamTypeI = t; }

  Bool_t  IsPP() const { return fBeamTypeI==kJPP; }
  Bool_t  IsPbPb() const { return fBeamTypeI==kJPbPb; }

  Float_t GetBeamEnergy() const { return fBeamEnergy; }
  void    SetBeamEnergy( Float_t e ){ fBeamEnergy = e; }

  Bool_t  IsMC()const{ return fIsMC; }
  void    SetIsMC(Bool_t b){ fIsMC=b; }


  //==== Production Info ====//
  Int_t   GetInputFormat() const { return fInputFormat; }
  void    SetInputFormat( Int_t v ){ fInputFormat = v; }
  Bool_t  FromESD() const { return fInputFormat==kJESD; }
  Bool_t  FromAOD() const { return fInputFormat==kJAOD; }

  Bool_t  GetWithoutSDD()const{ return fWithoutSDD; }
  void    SetWithoutSDD( bool b ){ fWithoutSDD = b; }

  void    SetStoreEventPlaneSource(Bool_t dostore ){ fStoreEventPlaneSource = dostore; }
  Bool_t  GetStoreEventPlaneSource() const { return fStoreEventPlaneSource; };

  void    SetStoreEMCalInfo(Bool_t dostore ){ fStoreEMCalInfo = dostore; }
  Bool_t  GetStoreEMCalInfo() const { return fStoreEMCalInfo; };

  UInt_t  GetStoreTPCTrackBitMask()const{ return fStoreTPCTrackBitMask; }
  void    SetStoreTPCTrackBitMask( UInt_t mask ){ fStoreTPCTrackBitMask = mask; }

  UInt_t  GetStoreGCGTrackBitMask()const{ return fStoreGCGTrackBitMask; }
  void    SetStoreGCGTrackBitMask( UInt_t mask ){ fStoreGCGTrackBitMask = mask; }

  void    SetESDInfo(const TString info ) { fESDInfo = info; }
  TString GetESDInfo() const { return fESDInfo; }

  Bool_t GetRefitESDVertexTracks() const { return fRefitESDVertexTracks; }
  void    SetRefitESDVertexTracks( Bool_t refit ){ fRefitESDVertexTracks=refit; } 

  //==== Detector Status ====//
  Short_t   GetL3MagnetFieldPolarity()  const { return fL3MagnetPolarity; }
  void      SetL3MagnetFieldPolarity(Short_t p){ fL3MagnetPolarity=p; }

  Double_t  GetL3MagnetFieldIntensity() const { return fMagneticFieldL3; }
  void      SetL3MagnetFieldIntensity(Double_t i) { fMagneticFieldL3=i; }

  Float_t   GetCurrentL3() const { return fCurrentL3; }
  void      SetCurrentL3(Float_t current){ fCurrentL3 = current; }

  Float_t   GetCurrentDip() const { return fCurrentDip; }
  void      SetCurrentDip(Float_t dip){ fCurrentDip=dip; }

  Bool_t    IsUniformBMap() const { return fUniformBMap; }
  void      SetUniformBMap(Bool_t uniform){ fUniformBMap=uniform; }

  //==== Trigger Information ====//
  TString   GetFiredTriggers()const{ return fFiredTriggers; }
  void      SetFiredTriggers( TString trgs ){ fFiredTriggers=trgs; }

  ULong64_t GetTriggerMask()const{ return fTriggerMask; }
  void      SetTriggerMask( ULong64_t mask ){ fTriggerMask=mask; }

  UChar_t   GetTriggerCluster() const { return fTriggerCluster; }
  void      SetTriggerCluster( UChar_t cluster ){ fTriggerCluster=cluster; } 

  Int_t     GetSizeOfTableJCorran() const { return fSizeOfTableJCorran; }
  void      SetSizeOfTableJCorran( Int_t size ){ fSizeOfTableJCorran=size; } 

  std::vector<TString> GetActiveTriggersJCorran() const { return fActiveTriggersJCorran; }
  void    SetActiveTriggersJCorran( std::vector<TString> strs ){ fActiveTriggersJCorran=strs; } 
  TString   GetActiveTriggersJCorran(int i) const { return fActiveTriggersJCorran[i]; }
  void      SetActiveTriggersJCorran(int i, TString str){ fActiveTriggersJCorran[i]=str; }

  std::vector<TString> GetActiveTriggersAlice() const { return fActiveTriggersAlice; }
  void      SetActiveTriggersAlice( std::vector<TString> strs ){ fActiveTriggersAlice=strs; } 
  TString   GetActiveTriggersAlice(int i) const { return fActiveTriggersAlice[i]; }
  void      SetActiveTriggersAlice(int i, TString str){ fActiveTriggersAlice[i]=str; }
  Int_t     GetActiveTriggerBitAlice(TString TriggerName);

  //==== GENERAL STUFF for Class ====//

  void PrintOut();

 protected:
  //==== General Info ====//
  Int_t       fRunNumber;        //run number 
  TString     fRunType;       // ex) LHC10h
  TString     fBeamType;        // beam type kJPP, kJPbPb, kJPPb
  Int_t       fBeamTypeI;        // beam type kJPP, kJPbPb, kJPPb
  Float_t     fBeamEnergy; // beam energy
  Bool_t      fIsMC;       // MC data or real data

  //==== Production Info ====//
  Int_t       fInputFormat; // specify the input data format (kJESD or kJAOD)
  Bool_t      fWithoutSDD;  // is SDD detector used or not
  Bool_t      fStoreEventPlaneSource; // store event plane
  Bool_t      fStoreEMCalInfo; // store event plane
  UInt_t      fStoreTPCTrackBitMask;  // TPC bit mask
  UInt_t      fStoreGCGTrackBitMask;  // GCG bit maks
  TString     fESDInfo;       // information of aliroot,  root version while esd production
  Bool_t      fRefitESDVertexTracks;   // refit to ESD vertex tracks

  //==== Detector Status ====//
  Short_t     fL3MagnetPolarity; //Polarity of magnetic filed in L3 magnet (LHC convention: + -> +Bz)
  Double32_t  fMagneticFieldL3;  //Solenoid Magnetic Field in kG   
  Float_t     fCurrentL3; // L3 current
  Float_t     fCurrentDip; // dipole current
  Bool_t      fUniformBMap; // beam uniformity
  
  //==== Trigger Information ====//
  TString     fFiredTriggers;       // String with fired triggers
  ULong64_t   fTriggerMask;         // Trigger Type (mask)
  UChar_t     fTriggerCluster;      // Trigger cluster (mask)
  Int_t           fSizeOfTableJCorran;  //size of jcorran table
  std::vector<TString> fActiveTriggersJCorran;   //array maping between trigger bit and trigger names
  std::vector<TString> fActiveTriggersAlice;   //array maping between trigger bit and trigger names
  //TBit 0 = MB 
  ClassDef(AliJRunHeader,3)

};

#endif

