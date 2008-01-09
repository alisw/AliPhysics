#ifndef AliAODDimuon_H
#define AliAODDimuon_H

/* AliAODDimuon: a class for AODs for the MUON Arm of the ALICE Experiment
 * Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
 * INFN of Torino - Italy
 */

/* 2007/07/07 v0.00 Initial version */
/* 2007/12/06 v0.01 Introduction of AliAODEventInfo */
/* 2007/12/18 v0.02 Corrected CostCS for Like-Sign, added CostKh, CostHe and xf*/

#include "TRef.h"
#include "AliVParticle.h"
#include "TLorentzVector.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODEventInfo.h"

class AliAODDimuon: public AliVParticle {
public:
  AliAODDimuon();
  AliAODDimuon(const AliAODDimuon& dimu);
  AliAODDimuon(TObject *mu0, TObject *mu1, TObject *evpoint=0);
  ~AliAODDimuon();

  // Data members
  TRef mu[2];	// Pointers to the reconstructed muons
  TRef ei;	// Pointer to the EventInfo object
  TLorentzVector *p; //! TLorentzVector of dimuon momentum (not stored into file)

  // Methods to access kinematics
  virtual Double_t Px() const;
  virtual Double_t Py() const;
  virtual Double_t Pz() const;
  virtual Bool_t PxPyPz(Double_t* p) const { p[0]=Px(); p[1]=Py(); p[2]=Pz(); return 1;}
  virtual Double_t Pt() const;
  virtual Double_t P() const;

  virtual Double_t OneOverPt() const {return Pt()>0 ? 1./Pt() : -999999999;}
  virtual Double_t Phi() const;
  virtual Double_t Theta() const;

  virtual Double_t E() const;
  virtual Double_t M() const;
  
  virtual Double_t Eta() const;
  virtual Double_t Y() const;
  
  virtual Short_t Charge() const;

  // Dimuon vertex will be implemented when the muon track covariance matrix 
  // at vertex will be included in the ESD (and AOD)
  // It would require also the information about magnetic field when filling AOD
  virtual Double_t Xv() const {return -999999999;}
  virtual Double_t Yv() const {return -999999999;}
  virtual Double_t Zv() const {return -999999999;}
  virtual Bool_t XvYvZv(Double_t* v) const { v[0]=-999999999; v[1]=-999999999; v[2]=-999999999; return 0;}

  Double_t P();
  Double_t Phi();
  Double_t Theta();
  Double_t M();
  Double_t Eta();
  Double_t Y();

  // Added functions
  Double_t xf();     // Feynman x
  Double_t CostCS(); // Cosinus of the Collins-Soper polar decay angle
  Double_t CostHe(); // Cosinus of the Helicity polar decay angle
  Int_t AnyPt();
  Int_t LowPt();
  Int_t HighPt();
  Double_t MaxChi2Match();
  // PID
  virtual const Double_t *PID() const {return 0;} // return PID object (to be defined, still)

  // Additional getters
  AliAODTrack* GetMu(Int_t imu=0){return (imu==0||imu==1)&&(mu[imu]!=0) ? (AliAODTrack*)mu[imu].GetObject() : 0; } // Get a pointer to a muon
  AliAODTrack* Mu(Int_t imu=0){return (imu==0||imu==1)&&(mu[imu]!=0) ? (AliAODTrack*)mu[imu].GetObject() : 0; } // Get a pointer to a muon
  AliAODEventInfo* MuonHeader(){return (ei!=0) ? (AliAODEventInfo*)ei.GetObject() : 0; } // Get a pointer to the AliAODEventInfo

  // Useful constants
  Double_t MProton; //! Proton mass (not stored into file)

private:
  Int_t CheckPointers() const;
  void BookP();
  ClassDef(AliAODDimuon,1)  // AliAODDimuon track
};

#endif
