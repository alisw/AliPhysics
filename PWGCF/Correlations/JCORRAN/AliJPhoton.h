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
    kPhotonAli = 5,
    kPi0Ali = 6,
    kNeutronAli = 7,
    kKaon0Ali = 8,
    kEleConAli = 9,
    kUnknownAli = 10
  };



  AliJPhoton();	    //default constructor
  AliJPhoton(const AliJPhoton& a); //copy constructor
  virtual ~AliJPhoton(){		//destructor    
    if(fCellsAbsId)       delete [] fCellsAbsId; 
    if(fCellsAmpFraction) delete [] fCellsAmpFraction;
  }

  float  GetChi2() const {return fChi2;}
  float  GetTof() const {return fTof;}                   
  float  GetX() const {return fX;}            
  float  GetY() const {return fY;}          
  float  GetZ() const {return fZ;}
  float  GetProbPhot() const {return fProbPhot;}

  void  SetChi2(float chi2) {fChi2=chi2;}
  void  SetTof(float tof) {fTof=tof;}
  void  SetPositionX(float x) {fX=x;}
  void  SetPositionY(float y) {fY=y;}
  void  SetPositionZ(float z) {fZ=z;}
  void  SetProbPhot(float prob) {fProbPhot=prob;}

  AliJPhoton& operator=(const AliJPhoton& photon);

  //TODO
  Bool_t InPHOS()  const {return fCaloType==kPHOSCalo ? kTRUE:kFALSE ;}
  Bool_t InEMCAL() const {return fCaloType==kEMCALCalo ? kTRUE:kFALSE ;}
  // getters
  void GetPID(Double_t *pid) const {
    for(Int_t i=0; i<kUnknownAli+1; ++i) pid[i]=fCaloPID[i];
  }
  Double_t  GetDistToBadChannel() const {return fDistToBadChannel;}
  Double_t  GetDispersion()       const {return fDispersion;}
  Double_t  GetM20()              const {return fM20;}
  Double_t  GetM02()              const {return fM02;}
  Double_t  GetEmcCpvDist()       const {return fEmcCpvDist;}
  Double_t  GetNCells() const   { return fNCells;}
  UShort_t  *GetCellsAbsId() const {return  fCellsAbsId;}
  Int_t     GetCellAbsId(Int_t i) const {  
         if (fCellsAbsId && i >=0 && i < fNCells ) return fCellsAbsId[i];    
        else return -1;}
  Double32_t *GetCellsAmplitudeFraction() const {return  fCellsAmpFraction;}
  Double_t    GetCellAmplitudeFraction(Int_t i) const {  
              if (fCellsAmpFraction && i >=0 && i < fNCells ) return fCellsAmpFraction[i];    
              else return -1;}
  particleType GetParticleType();
  Int_t    GetSuperModuleID() const { return fSuperModuleId; }
  void     SetSuperModuleID(Int_t id) { fSuperModuleId = id; }
  Double32_t  GetTrackDx() const {return fTrackDx;}
  Double32_t  GetTrackDz() const {return fTrackDz;}

  //setters
  void  SetCaloType(AliJCaloType calo) {fCaloType = calo;}
  void  SetDistToBadChannel(Double_t dist) {fDistToBadChannel = dist;}
  void  SetDispersion(Double_t disp) {fDispersion = disp;}
  void  SetM20(Double_t m20) {fM20 = m20;}
  void  SetM02(Double_t m02) {fM02 = m02;}
  void  SetEmcCpvDist(Double_t dist) {fEmcCpvDist = dist;} 
  void  SetPID(const Double_t *pid);
  void  SetNCells(Int_t n) { fNCells = n;}
  void  SetCellsAbsId(const UShort_t *array);
  void  SetCellsAmplitudeFraction(const Double32_t *array);
  void  SetTrackDx(Double_t trackDx) {fTrackDx = trackDx;}
  void  SetTrackDz(Double_t trackDz) {fTrackDz = trackDz;}

 private:

  Double32_t  fChi2;      //chi2             
  Double32_t  fTof;       //time of flight                
  Double32_t  fX, fY, fZ; // x,y,z coordinates              
  Double32_t  fProbPhot;  //probability to be a photon
  Double32_t  fTrackDx, fTrackDz; // Distance to closest track in phi and z

 //TODO
  AliJCaloType   fCaloType;              // PHOS or EMCAL photon
  Double_t       fCaloPID[kUnknownAli+1];          // [0.,1.,8] pointer to PID object
  Double_t       fDistToBadChannel;      // Distance to nearest bad channel
  Double_t       fDispersion;            // cluster dispersion, for shape analysis
  Double_t       fM20;                   // 2-nd moment along the main eigen axis
  Double_t       fM02;                   // 2-nd moment along the second eigen axis
  Double_t       fEmcCpvDist;            // the distance from PHOS EMC rec.point to the closest CPV rec.point
 
  Int_t          fNCells ;                 //number of cells
  Int_t          fSuperModuleId ;          //super module id
  UShort_t      *fCellsAbsId;            //[fNCells] array of cell absId numbers
  Double32_t    *fCellsAmpFraction;    //[fNCells][0.,1.,16] array with cell amplitudes fraction.

  ClassDef(AliJPhoton,1)

};

#endif

