/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

// $Id: AliJPhoton.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJPhoton.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
  */
////////////////////////////////////////////////////

#ifndef ALIJPHOTON_H
#define ALIJPHOTON_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJBaseTrack.h"

class AliJPhoton : public AliJBaseTrack {

 public:
  enum AliJCaloType {
         kUndef = -1, 
     kPHOSCalo,  
     kEMCALCalo
  };
  
  enum AliJCaloPID{
    kElectronAli = 0,
    kMuonAli = 1,
    kPionAli = 2,
    kKaonAli = 3,
    kProtonAli = 4,
		kDeuteronAli = 5,
		kTritonAli = 6,
		kHe3Ali = 7,
		kAlpha = 8,
		kPhotonAli = 9,
    kPi0Ali = 10,
    kNeutronAli = 11,
    kKaon0Ali = 12,
    kEleConAli = 13,
    kUnknownAli = 14
  };

  enum { kIsConversion=AliJBaseTrack::kNFlag, kNFlag };


  AliJPhoton();      //default constructor
  AliJPhoton(const AliJPhoton& a); //copy constructor
  virtual ~AliJPhoton();    //destructor    

  void Clear( Option_t * /* option */ ); // cleaner  

  Double32_t  GetChi2() const {return fChi2;}
  Double32_t  GetTof() const {return fTof;}                   
  Double32_t  GetX() const {return fX;}            
  Double32_t  GetY() const {return fY;}          
  Double32_t  GetZ() const {return fZ;}
  Double32_t  GetProbPhot() const {return fProbPhot;}

  void  SetChi2(Double32_t chi2) {fChi2=chi2;}
  void  SetTof(Double32_t tof) {fTof=tof;}
  void  SetPositionX(Double32_t x) {fX=x;}
  void  SetPositionY(Double32_t y) {fY=y;}
  void  SetPositionZ(Double32_t z) {fZ=z;}
  void  SetProbPhot(Double32_t prob) {fProbPhot=prob;}

  AliJPhoton& operator=(const AliJPhoton& photon);

  //TODO
  Bool_t IsPHOS()  const {return fCaloType==kPHOSCalo ? kTRUE:kFALSE ;}
  Bool_t IsEMCAL() const {return fCaloType==kEMCALCalo ? kTRUE:kFALSE ;}
  // getters
  void GetPID(Double_t *pid) const {
    for(Int_t i=0; i<kUnknownAli+1; ++i) pid[i]=fCaloPID[i];
  }
  Double32_t  GetDistToBadChannel() const {return fDistToBadChannel;}
  Double32_t  GetDispersion()       const {return fDispersion;}
  Double32_t  GetM20()              const {return fM20;}
  Double32_t  GetM02()              const {return fM02;}
  Double32_t  GetEmcCpvDist()       const {return fEmcCpvDist;}
  Int_t       GetNCells() const   { return fNCells;}
  Int_t       GetNTracksMatched() const   { return fNTracksMatched;}
  UShort_t  *GetCellsAbsId() const {return  fCellsAbsId;}
  UShort_t   GetCellAbsId(Int_t i) const {  
         if (fCellsAbsId && i >=0 && i < fNCells ) return fCellsAbsId[i];    
        else return -1;}
  Double32_t *GetCellsAmplitudeFraction() const {return  fCellsAmpFraction;}
  Double32_t  GetCellAmplitudeFraction(Int_t i) const {  
              if (fCellsAmpFraction && i >=0 && i < fNCells ) return fCellsAmpFraction[i];    
              else return -1;}
  Int_t       GetNEMCLabel() const {return fNEMCLabel;}
  Int_t       GetEMCLabel( Int_t pos = 0 ) const { 
              if( fEMCLabel && pos >= 0 && pos < fNEMCLabel ) return fEMCLabel[pos];
              else return -1;}
//   Double32_t  *GetCellsAmplitude() const {return  fCellsAmp;}
//   Double32_t   GetCellAmplitude(Int_t i) const {  
//               if (fCellsAmp && i >=0 && i < fNCells ) return fCellsAmp[i];    
//               else return -1;}
  particleType GetParticleType();
  Int_t    GetSuperModuleID() const { return fSuperModuleId; }
  void     SetSuperModuleID(Int_t id) { fSuperModuleId = id; }
  Double32_t  GetTrackDx() const {return fTrackDx;}
  Double32_t  GetTrackDz() const {return fTrackDz;}
  Double32_t GetEMax() const {return fEMax;}
  Double32_t GetECross() const {return fECross;}
  Double32_t GetECore() const {return fECore;}
  Int_t      GetNLM() const {return fNLM;}
  Int_t  *GetCellsIndex() const {return  fCellsIndex; }
  Int_t   GetCellIndex(Int_t i) const {  
         if (fCellsIndex && i >=0 && i < fNCells ) return fCellsIndex[i];    
        else return -1;}

  //setters
  void  SetCaloType(AliJCaloType calo) {fCaloType = calo;}
  void  SetDistToBadChannel(Double32_t dist) {fDistToBadChannel = dist;}
  void  SetDispersion(Double32_t disp) {fDispersion = disp;}
  void  SetM20(Double32_t m20) {fM20 = m20;}
  void  SetM02(Double32_t m02) {fM02 = m02;}
  void  SetEmcCpvDist(Double32_t dist) {fEmcCpvDist = dist;} 
  void  SetPID(const Double32_t *pid);
  void  SetNCells(Int_t n) { fNCells = n;}
  void  SetNTracksMatched( Int_t n ) { fNTracksMatched = n;}
  void  SetCellsAbsId(const UShort_t *array);
  void  SetCellsAmplitudeFraction(const Double32_t *array);
  void  SetNEMCLabel(Int_t n) { fNEMCLabel = n;}
  void  SetEMCLabel(const Int_t *array);
  void  SetEMCLabel(const Int_t pos, const Int_t lab) {
        if( fEMCLabel && pos >= 0 && pos < fNEMCLabel )
          fEMCLabel[pos] = lab; }
//   void  SetCellsAmplitude(const Double32_t *array);
  void  SetTrackDx(Double32_t trackDx) {fTrackDx = trackDx;}
  void  SetTrackDz(Double32_t trackDz) {fTrackDz = trackDz;}
  void  SetEMax(Double32_t e) {fEMax = e; }
  void  SetECross(Double32_t e) {fECross = e; }
  void  SetECore(Double32_t e) {fECore = e; }
  void  SetNLM( Int_t n ) {fNLM = n;}
  void  SetCellIndex( const Int_t pos, const Int_t ind );
  void  SetCellsIndex(const Int_t *array);
  
  void  ClearCellsIndex();

 private:

  Double32_t  fChi2;      //chi2             
  Double32_t  fTof;       //time of flight                
  Double32_t  fX, fY, fZ; // x,y,z coordinates              
  Double32_t  fProbPhot;  //probability to be a photon
  Double32_t  fTrackDx, fTrackDz; // Distance to closest track in phi and z

 //TODO
  AliJCaloType   fCaloType;              // PHOS or EMCAL photon
  Double32_t     fCaloPID[kUnknownAli+1];          // [0.,1.,8] pointer to PID object
  Double32_t     fDistToBadChannel;      // Distance to nearest bad channel
  Double32_t     fDispersion;            // cluster dispersion, for shape analysis
  Double32_t     fM20;                   // 2-nd moment along the main eigen axis
  Double32_t     fM02;                   // 2-nd moment along the second eigen axis
  Double32_t     fEmcCpvDist;            // the distance from PHOS EMC rec.point to the closest CPV rec.point
 
  Int_t          fNCells ;                 //number of cells
  Int_t          fNTracksMatched;          // number of tracks matched
  Int_t          fSuperModuleId ;          //super module id
  UShort_t      *fCellsAbsId;            //[fNCells] array of cell absId numbers
  Double32_t    *fCellsAmpFraction;    //[fNCells][0.,1.,16] array with cell amplitudes fraction (elements are 0 if unfolding off)
//  Double32_t    *fCellsAmp;            //[fNCells] array amplitudes of cluster cellsz
  Int_t          fNEMCLabel;            // number of MC labels
  Int_t         *fEMCLabel;             //[fNEMCLabel] MC labels

  // analysis time variables
  Double32_t     fEMax;  //! maximum cell energy
  Double32_t     fECross; //! energy of cells in cross around max
  Double32_t     fECore; //! central cells energies
  Int_t          fNLM;  //! number of local maxima
  Int_t         *fCellsIndex; //! cell indices in cell list

  ClassDef(AliJPhoton,2)

};

#endif

