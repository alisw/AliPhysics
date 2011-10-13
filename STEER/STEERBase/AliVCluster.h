#ifndef ALIVCLUSTER_H
#define ALIVCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

//-------------------------------------------------------------------------
//
//   Virtual class to access calorimeter 
//   (EMCAL, PHOS, PMD, FMD) cluster data
//   Author: Gustavo Conesa Balbastre LPSC-Grenoble
//
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TLorentzVector.h>

class AliVCluster : public TObject 
{
  
 public:
  
  AliVCluster() { ; }
  virtual ~AliVCluster() { ; }
  AliVCluster(const AliVCluster& clus);
  AliVCluster & operator=(const AliVCluster& source);
  void Clear(const Option_t*) {;}
  
  enum VClu_t {kUndef = -2, 
	       kPHOSNeutral, 
		     kPHOSCharged,
	       kEMCALClusterv1,		 
	       kPMDNeutral, 
	       kPMDCharged};
  
  enum VCluPID_t {
    kElectron = 0,
    kMuon     = 1,
    kPion     = 2,
    kKaon     = 3,
    kProton   = 4,
    kPhoton   = 5,
    kPi0      = 6,
    kNeutron  = 7,
    kKaon0    = 8,
    kEleCon   = 9,
    kUnknown  = 10,
    kCharged  = 11,//For PMD?
    kNeutral  = 12 //For PMD? 
  };
  
  //Common EMCAL/PHOS/FMD/PMD
  
  virtual void        SetID(Int_t )                 { ; }
  virtual Int_t       GetID() const                 {return 0 ; }
  
  virtual void        SetType(Char_t )              { ; }  
  virtual Char_t      GetType() const               {return kUndef ; } 
  
  virtual void        SetE(Double_t )               { ; }
  virtual Double_t       E() const                  {return 0. ; }
  
  virtual void        SetChi2(Double_t )            { ; }
  virtual Double_t       Chi2() const               {return 0. ; }
  
  virtual void        SetPositionAt(Float_t,Int_t)  { ; }
  virtual void        SetPosition(Float_t *)        { ; }
  virtual void        GetPosition(Float_t *) const  { ; }	
  
  virtual void        SetPIDAt(Float_t , Int_t)     { ; }
  virtual void        SetPID(const Float_t *)       { ; }
  virtual const Double_t *GetPID() const            { return 0 ; }
  
  //CaloClusters, PHOS/EMCAL
  
  virtual Bool_t      IsEMCAL() const               {return kFALSE ; }
  virtual Bool_t      IsPHOS()  const               {return kFALSE ; }
  
  virtual void        SetDispersion(Double_t )      { ; }
  virtual Double_t    GetDispersion() const         {return 0. ;}
  
  virtual void        SetM20(Double_t)              { ; }
  virtual Double_t    GetM20() const                {return 0. ; }
  
  virtual void        SetM02(Double_t)              { ; }
  virtual Double_t    GetM02() const                {return 0. ; }
  
  virtual void        SetNExMax(UChar_t)            { ; }
  virtual UChar_t     GetNExMax() const             {return 0 ; } 
  
  virtual void        SetTOF(Double_t)              { ; }
  virtual Double_t    GetTOF() const                {return 0. ; }
  
  virtual void        SetEmcCpvDistance(Double_t)   { ; }
  virtual Double_t    GetEmcCpvDistance() const     {return 0. ; }
  virtual void        SetTrackDistance(Double_t, Double_t ){ ; }
  virtual Double_t    GetTrackDx(void)const         {return 0. ; }
  virtual Double_t    GetTrackDz(void)const         {return 0. ; }
  
  virtual void        SetDistanceToBadChannel(Double_t) { ; }
  virtual Double_t    GetDistanceToBadChannel() const   {return 0. ; }
  
  virtual void        SetNCells(Int_t)              { ; }
  virtual Int_t       GetNCells() const             {return 0 ; }
  
  virtual UShort_t   *GetCellsAbsId()               {return 0 ; }
  virtual Double_t   *GetCellsAmplitudeFraction()   {return 0 ; }
  virtual Int_t       GetCellAbsId(Int_t) const     {return 0 ; }  
  virtual Double_t    GetCellAmplitudeFraction(Int_t) const {return 0. ; }
  
  virtual Int_t       GetLabel() const              {return -1 ;}
  virtual Int_t       GetLabelAt(UInt_t) const      {return -1 ;}
  virtual Int_t      *GetLabels() const             {return 0 ; }
  virtual UInt_t      GetNLabels() const            {return 0 ; }
  
  virtual Int_t       GetNTracksMatched() const     {return 0 ; }
  virtual TObject    *GetTrackMatched(Int_t) const  {return 0 ; }//AODCaloCluster
  virtual Int_t       GetTrackMatchedIndex() const  {return -1; }//ESDCaloCluster
  
  virtual void GetMomentum(TLorentzVector &/*tl*/, Double_t * /*v*/) { ; }
  
  ClassDef(AliVCluster,0)  //VCluster 
    };

#endif //ALIVCLUSTER_H

