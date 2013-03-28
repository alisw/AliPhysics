#ifndef ALIESDCALOCLUSTER_H
#define ALIESDCALOCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */
/* $Log $ */

//-------------------------------------------------------------------------
//                          Class AliESDCaloCluster
//   This is the class to deal with during the physics analysis of data
//
//   New container for calorimeter clusters, which are the effective 
//   "tracks" for calorimeter detectors.  Can be used by PHOS and EMCAL
//
//     J.L. Klay (LLNL)
//-------------------------------------------------------------------------

#include <AliVCluster.h>
#include "AliPID.h"
#include "TArrayS.h"
#include "TArrayI.h"
#include "AliLog.h"

class TLorentzVector;

class AliESDCaloCluster : public AliVCluster 
{
  
 public:
  
  AliESDCaloCluster();
  AliESDCaloCluster(const AliESDCaloCluster& clus);
  AliESDCaloCluster & operator=(const AliESDCaloCluster& source);
  virtual ~AliESDCaloCluster();
  virtual void Copy(TObject &) const;
  void Clear(const Option_t*);

  void  SetID(Int_t id) {fID = id;}
  Int_t GetID() const {return fID;}
  
  void   SetType(Char_t type) { fClusterType = type; }
  Char_t GetType() const {return fClusterType; }
  
  Bool_t IsEMCAL() const {if(fClusterType == kEMCALClusterv1) return kTRUE; else return kFALSE;}
  Bool_t IsPHOS()  const {if(fClusterType == kPHOSNeutral || fClusterType == kPHOSCharged) return kTRUE;
    else return kFALSE;}
  
  void GetPosition  (Float_t *x) const {
    x[0]=fGlobalPos[0]; x[1]=fGlobalPos[1]; x[2]=fGlobalPos[2];}
  void SetPosition  (Float_t *x);
  void SetPositionAt(Float_t pos, Int_t ipos) {if(ipos>=0 && ipos<3) fGlobalPos[ipos] = pos ; 
    else AliInfo(Form("Bad index for position array, i = %d\n",ipos));}
  
  void  SetE(Double_t ene) { fEnergy = ene;}
  Double_t E() const       { return fEnergy;}
  
  void     SetDispersion(Double_t disp)  { fDispersion = disp; }
  Double_t GetDispersion() const         { return fDispersion; }
  
  void  SetChi2(Double_t chi2)  { fChi2 = chi2; }
  Double_t Chi2() const         { return fChi2; }
  
  const Double_t *GetPID() const { return fPID; }
  //for(Int_t i=0; i<AliPID::kSPECIESCN; ++i) pid[i]=fPID[i];}
  void SetPID  (const Float_t *pid) ;
  void SetPIDAt(Float_t p, Int_t i) {if(i>=0 && i<AliPID::kSPECIESCN) fPID[i] = p ; 
    else AliInfo(Form("Bad index for PID array, i = %d \n",i));}
  
  void     SetM20(Double_t m20) { fM20 = m20; }
  Double_t GetM20() const       { return fM20; }
  
  void     SetM02(Double_t m02) { fM02 = m02; }
  Double_t GetM02() const       { return fM02; }
  
  void    SetNExMax(UChar_t nExMax) { fNExMax = nExMax; }
  UChar_t GetNExMax() const         { return fNExMax; }
  
  void SetEmcCpvDistance(Double_t dEmcCpv) { fEmcCpvDistance = dEmcCpv; }
  Double_t GetEmcCpvDistance() const       { return fEmcCpvDistance; }
  void SetTrackDistance(Double_t dx, Double_t dz){fTrackDx=dx; fTrackDz=dz;}
  Double_t GetTrackDx(void)const {return fTrackDx;}
  Double_t GetTrackDz(void)const {return fTrackDz;}
  
  void     SetDistanceToBadChannel(Double_t dist) {fDistToBadChannel=dist;}
  Double_t GetDistanceToBadChannel() const        {return fDistToBadChannel;}
  
  void     SetTOF(Double_t tof) { fTOF = tof; }
  Double_t GetTOF() const       { return fTOF; }
  
  void AddTracksMatched(TArrayI & array)  { 
    if(!fTracksMatched)fTracksMatched   = new TArrayI(array);
    else *fTracksMatched = array;
  }
  void AddLabels(TArrayI & array)         { 
    if(!fLabels)fLabels = new TArrayI(array) ;
    else *fLabels = array;
  }
  
  void SetLabel(Int_t *array, UInt_t size)
  {
    if(fLabels) delete fLabels ;
    fLabels = new TArrayI(size,array);
  }

  TArrayI * GetTracksMatched() const  {return  fTracksMatched;}
  TArrayI * GetLabelsArray() const    {return  fLabels;}
  Int_t   * GetLabels() const         {if (fLabels) return  fLabels->GetArray(); else return 0;}

  Int_t GetTrackMatchedIndex() const   
  {if( fTracksMatched &&  fTracksMatched->GetSize() >0)  return  fTracksMatched->At(0); 
    else return -1;} //Most likely the track associated to the cluster
  
  Int_t GetLabel() const   {
    if( fLabels &&  fLabels->GetSize() >0)  return  fLabels->At(0); 
    else return -1;} //Most likely the track associated to the cluster
  Int_t GetLabelAt(UInt_t i) const {
    if (fLabels && i < (UInt_t)fLabels->GetSize()) return fLabels->At(i);
    else return -999; }
  
  Int_t GetNTracksMatched() const { if (fTracksMatched) return  fTracksMatched->GetSize(); 
    else return -1;}
  UInt_t GetNLabels() const       { if (fLabels) return  fLabels->GetSize(); 
    else return (0);}
  
  void GetMomentum(TLorentzVector& p, Double_t * vertexPosition );
  
  void  SetNCells(Int_t n)  { fNCells = n;}
  Int_t GetNCells() const   { return fNCells;}
  
  void      SetCellsAbsId(UShort_t *array) ;
  UShort_t *GetCellsAbsId() {return  fCellsAbsId;}
  
  void        SetCellsAmplitudeFraction(Double32_t *array) ;
  Double32_t *GetCellsAmplitudeFraction() {return  fCellsAmpFraction;}
  
  Int_t GetCellAbsId(Int_t i) const {  
    if (fCellsAbsId && i >=0 && i < fNCells ) return fCellsAbsId[i];    
    else return -1;}
  
  Double_t GetCellAmplitudeFraction(Int_t i) const {  
    if (fCellsAmpFraction && i >=0 && i < fNCells ) return fCellsAmpFraction[i];    
    else return -1;}
  
 protected:
  
  TArrayI * fTracksMatched; //Index of tracks close to cluster. First entry is the most likely match.
  TArrayI * fLabels;        //list of primaries that generated the cluster, ordered in deposited energy.
  
  Int_t        fNCells ;
  UShort_t   * fCellsAbsId;       //[fNCells] array of cell absId numbers
  Double32_t * fCellsAmpFraction; //[fNCells][0.,1.,16] array with cell amplitudes fraction.
  
  Double32_t   fGlobalPos[3];     // position in global coordinate systemD
  Double32_t   fEnergy;           // energy measured by calorimeter
  Double32_t   fDispersion;       // cluster dispersion, for shape analysis
  Double32_t   fChi2;             // chi2 of cluster fi
  Double32_t   fM20;              // 2-nd moment along the main eigen axis
  Double32_t   fM02;              // 2-nd moment along the second eigen axis
  
  Double32_t   fEmcCpvDistance;   // the distance from PHOS EMC rec.point to the closest CPV rec.point
  Double32_t   fTrackDx ;         // Distance to closest track in phi
  Double32_t   fTrackDz ;         // Distance to closest track in z
  
  Double32_t   fDistToBadChannel; // Distance to nearest bad channel
  Double32_t   fPID[AliPID::kSPECIESCN]; //[0,1,8]"detector response  probabilities" (for the PID)
  Int_t        fID;                // Unique Id of the cluster
  UChar_t      fNExMax ;           // number of (Ex-)maxima before unfolding  
  Char_t       fClusterType;       // Flag for different cluster type/versions
  Double_t     fTOF;               //[0,0,12] time-of-flight
  
  
  ClassDef(AliESDCaloCluster,11)  //ESDCaloCluster 

    };

#endif 


