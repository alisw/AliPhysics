#ifndef ALIESDTRACK_H
#define ALIESDTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliESDtrack
//   This is the class to deal with during the physics analysis of data
//      
//         Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
/*****************************************************************************
 *  Use GetExternalParameters() and GetExternalCovariance() to access the    *
 *      track information regardless of its internal representation.         *
 * This formation is now fixed in the following way:                         *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *****************************************************************************/

#include <TBits.h>
#include <TObject.h>
#include "AliPID.h"
class AliKalmanTrack;

const Int_t kNPlane = 6;

class AliESDtrack : public TObject {
public:
  AliESDtrack();
  AliESDtrack(const AliESDtrack& track);
  virtual ~AliESDtrack();
  void SetID(Int_t id) { fID =id;}
  Int_t GetID(){ return fID;}
  void SetStatus(ULong_t flags) {fFlags|=flags;}
  void ResetStatus(ULong_t flags) {fFlags&=~flags;}
  Bool_t UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags);
  void SetImpactParameters(Float_t xy,Float_t z) {fD=xy; fZ=z;}
  void SetIntegratedLength(Double_t l) {fTrackLength=l;}
  void SetIntegratedTimes(const Double_t *times);
  void SetESDpid(const Double_t *p);
  void GetESDpid(Double_t *p) const;
  
  ULong_t GetStatus() const {return fFlags;}
  Int_t GetLabel() const {return fLabel;}
  Double_t GetAlpha() const {return fRalpha;}
  void GetExternalParameters(Double_t &x, Double_t p[5]) const;
  void GetExternalCovariance(Double_t cov[15]) const;

  Bool_t GetExternalParametersAt(Double_t x, Double_t p[5]) const;
  Bool_t GetPxPyPzAt(Double_t x, Double_t p[3]) const;
  Bool_t GetXYZAt(Double_t x, Double_t r[3]) const;

  void GetImpactParameters(Float_t &xy,Float_t &z) const {xy=fD; z=fZ;}
  Double_t GetIntegratedLength() const {return fTrackLength;}
  void GetIntegratedTimes(Double_t *times) const;
  Double_t GetMass() const;
  Double_t GetP() const;
  Bool_t GetPxPyPz(Double_t *p) const;
  Bool_t GetXYZ(Double_t *r) const;
  void GetCovariance(Double_t cov[21]) const;
  Int_t GetSign() const {return (fRp[4]>0) ? 1 : -1;} 

  void SetConstrainedTrackParams(const AliKalmanTrack *t, Double_t chi2);

  Double_t GetConstrainedAlpha() const {return fCalpha;}
  Double_t GetConstrainedChi2() const {return fCchi2;}
  void GetConstrainedExternalParameters(Double_t &x, Double_t p[5]) const;
  void GetConstrainedExternalCovariance(Double_t cov[15]) const;

  Bool_t GetConstrainedPxPyPz(Double_t *p) const;
  Bool_t GetConstrainedXYZ(Double_t *r) const;

  Bool_t GetInnerPxPyPz(Double_t *p) const;
  Bool_t GetInnerXYZ(Double_t *r) const;
  void GetInnerExternalParameters(Double_t &x, Double_t p[5]) const;//skowron
  void GetInnerExternalCovariance(Double_t cov[15]) const;//skowron
  Double_t GetInnerAlpha() const {return fIalpha;}
  
  void SetITSpid(const Double_t *p);
  void SetITSChi2MIP(const Float_t *chi2mip);
  void SetITStrack(AliKalmanTrack * track){fITStrack=track;}
  void GetITSpid(Double_t *p) const;
  Float_t GetITSsignal() const {return fITSsignal;}
  Float_t GetITSchi2() const {return fITSchi2;}
  Int_t GetITSclusters(UInt_t *idx) const;
  Int_t GetITSLabel() const {return fITSLabel;}
  Float_t GetITSFakeRatio() const {return fITSFakeRatio;}
  AliKalmanTrack * GetITStrack(){return fITStrack;}

  void SetTPCpid(const Double_t *p);
  void GetTPCpid(Double_t *p) const;
  void SetTPCPoints(Float_t points[4]){for (Int_t i=0;i<4;i++) fTPCPoints[i]=points[i];}
  void SetKinkIndexes(Int_t points[3]) {for (Int_t i=0;i<3;i++) fKinkIndexes[i] = points[i];}
  void SetV0Indexes(Int_t points[3]) {for (Int_t i=0;i<3;i++) fV0Indexes[i] = points[i];}
  Float_t GetTPCsignal() const {return fTPCsignal;}
  Float_t GetTPCchi2() const {return fTPCchi2;}
  Int_t GetTPCclusters(Int_t *idx) const;
  Int_t GetTPCLabel() const {return fTPCLabel;}
  Int_t GetKinkIndex(Int_t i) const { return fKinkIndexes[i];}
  Int_t GetV0Index(Int_t i) const { return fV0Indexes[i];}
  const TBits& GetTPCClusterMap() const {return fTPCClusterMap;}
  
  void SetTRDpid(const Double_t *p);
  void SetTRDtrack(AliKalmanTrack * track){fTRDtrack=track;}
  void SetTRDsignals(Float_t dedx, Int_t i) {fTRDsignals[i]=dedx;}
  void SetTRDTimBin(Int_t timbin, Int_t i) {fTRDTimBin[i]=timbin;}
  void GetTRDpid(Double_t *p) const;
  Float_t GetTRDsignal() const {return fTRDsignal;}
  Float_t GetTRDsignals(Int_t i) const {return fTRDsignals[i];}
  Int_t GetTRDTimBin(Int_t i) const {return fTRDTimBin[i];}
  Float_t GetTRDchi2() const {return fTRDchi2;}
  Int_t GetTRDclusters(UInt_t *idx) const;
  Int_t GetTRDncls() const {return fTRDncls;}
  void    SetTRDpid(Int_t iSpecies, Float_t p);
  Float_t GetTRDpid(Int_t iSpecies) const;
  Int_t GetTRDLabel() const {return fTRDLabel;}
  void GetTRDExternalParameters(Double_t &x, Double_t &alpha, Double_t p[5], Double_t cov[15]) const;//MI
  AliKalmanTrack * GetTRDtrack(){return fTRDtrack;}

  void SetTOFsignal(Double_t tof) {fTOFsignal=tof;}
  Float_t GetTOFsignal() const {return fTOFsignal;}
  Float_t GetTOFchi2() const {return fTOFchi2;}
  void    SetTOFpid(const Double_t *p);
  void    SetTOFLabel(const Int_t *p);
  void    GetTOFpid(Double_t *p) const;
  void    GetTOFLabel(Int_t *p) const;
  void    GetTOFInfo(Float_t *info) const;
  void    SetTOFInfo(Float_t *info);
  UInt_t  GetTOFcluster() const {return fTOFindex;}
  void  SetTOFcluster(UInt_t index) {fTOFindex=index;}
  
  void    SetRICHsignal(Double_t beta) {fRICHsignal=beta;}
  Float_t GetRICHsignal() const {return fRICHsignal;}
  void    SetRICHpid(const Double_t *p);
  void    GetRICHpid(Double_t *p) const;
  
  void SetPHOSposition(const Double_t *pos)  {
    fPHOSpos[0] = pos[0]; fPHOSpos[1]=pos[1]; fPHOSpos[2]=pos[2];
  }
  void SetPHOSsignal(Double_t ene) {fPHOSsignal = ene; }
  void SetPHOSpid(const Double_t *p);
  void GetPHOSposition(Double_t *pos) const {
    pos[0]=fPHOSpos[0]; pos[1]=fPHOSpos[1]; pos[2]=fPHOSpos[2];
  }
  Float_t GetPHOSsignal() const {return fPHOSsignal;}
  void GetPHOSpid(Double_t *p) const;  

  void SetEMCALposition(const Double_t *pos)  {
    fEMCALpos[0] = pos[0]; fEMCALpos[1]=pos[1]; fEMCALpos[2]=pos[2];
  }
  void SetEMCALsignal(Double_t ene) {fEMCALsignal = ene; }
  void SetEMCALpid(const Double_t *p);
  void GetEMCALposition(Double_t *pos) const {
    pos[0]=fEMCALpos[0]; pos[1]=fEMCALpos[1]; pos[2]=fEMCALpos[2];
  }
  Float_t GetEMCALsignal() const {return fEMCALsignal;}
  void GetEMCALpid(Double_t *p) const;  

  Bool_t IsOn(Int_t mask) const {return (fFlags&mask)>0;}
  Bool_t IsRICH()  const {return fFlags&kRICHpid;}
  Bool_t IsPHOS()  const {return fFlags&kPHOSpid;}
  Bool_t IsEMCAL() const {return fFlags&kEMCALpid;}

  virtual void Print(Option_t * opt) const ; 

  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kPHOSpid=0x10000, kRICHpid=0x20000, kEMCALpid=0x40000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
  }; 
protected:
  ULong_t   fFlags;        // Reconstruction status flags 
  Int_t     fLabel;        // Track label
  Int_t     fID;           // Unique ID of the track
  Float_t   fTrackLength;  // Track length
  Float_t   fD;            // Impact parameter in XY-plane
  Float_t   fZ;            // Impact parameter in Z 
  Float_t   fTrackTime[AliPID::kSPECIES]; // TOFs estimated by the tracking
  Float_t   fR[AliPID::kSPECIES];         // combined "detector response probability"

  Int_t     fStopVertex;          // Index of stop vertex

//Running track parameters
  Double_t fRalpha;  // track rotation angle
  Double_t fRx;      // X-coordinate of the track reference plane 
  Double_t fRp[5];   // external track parameters  
  Double_t fRc[15];  // external cov. matrix of the track parameters

//Track parameters constrained to the primary vertex
  Double_t fCalpha;   // Track rotation angle
  Double_t fCx;       // x-coordinate of the track reference plane
  Double_t fCp[5];    // external track parameters
  Double_t fCc[15];   // external cov. matrix of the track parameters
  Double_t fCchi2; //chi2 at the primary vertex

//Track parameters at the inner wall of the TPC
  Double_t fIalpha;   // Track rotation angle
  Double_t fIx;       // x-coordinate of the track reference plane
  Double_t fIp[5];    // external track parameters
  Double_t fIc[15];   // external cov. matrix of the track parameters

//Track parameters at the inner wall of the TRD 
  Double_t fTalpha;   // Track rotation angle
  Double_t fTx;       // x-coordinate of the track reference plane
  Double_t fTp[5];    // external track parameters
  Double_t fTc[15];   // external cov. matrix of the track parameters

  // ITS related track information
  Float_t fITSchi2;        // chi2 in the ITS
  Float_t fITSchi2MIP[12];     // chi2s in the ITS
  Int_t   fITSncls;        // number of clusters assigned in the ITS
  UInt_t  fITSindex[6];    //! indices of the assigned ITS clusters
  Float_t fITSsignal;      // detector's PID signal
  Float_t fITSr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fITSLabel;       // label according TPC
  Float_t fITSFakeRatio;   // ration of fake tracks
  AliKalmanTrack * fITStrack; //! OWNER: pointer to the ITS track -- currently for debug purpose
  
  // TPC related track information
  Float_t fTPCchi2;        // chi2 in the TPC
  Int_t   fTPCncls;        // number of clusters assigned in the TPC
  Int_t  fTPCindex[180];  //! indices of the assigned TPC clusters
  TBits   fTPCClusterMap;  // Map of clusters, one bit per padrow; 1 if has a cluster on given padrow
  Float_t fTPCsignal;      // detector's PID signal
  Float_t fTPCr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTPCLabel;       // label according TPC
  Float_t fTPCPoints[4];   // TPC points -first, max. dens, last and max density
  Int_t   fKinkIndexes[3]; // array of indexes of posible kink candidates 
  Int_t   fV0Indexes[3]; // array of indexes of posible kink candidates 

  // TRD related track information
  Float_t fTRDchi2;        // chi2 in the TRD
  Int_t   fTRDncls;        // number of clusters assigned in the TRD
  Int_t   fTRDncls0;       // number of clusters assigned in the TRD before first material cross
  UInt_t  fTRDindex[130];   //! indices of the assigned TRD clusters
  Float_t fTRDsignal;      // detector's PID signal
  Float_t fTRDsignals[kNPlane];  // TRD signals from all six planes
  Int_t fTRDTimBin[kNPlane];     // Time bin of Max cluster from all six planes
  Float_t fTRDr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTRDLabel;       // label according TRD
  AliKalmanTrack * fTRDtrack; //! OWNER: pointer to the TRD track -- currently for debug purpose

  // TOF related track information
  Float_t fTOFchi2;        // chi2 in the TOF
  UInt_t  fTOFindex;       // index of the assigned TOF cluster
  Float_t fTOFsignal;      // detector's PID signal
  Float_t fTOFr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTOFLabel[3];       // TOF label 
  Float_t fTOFInfo[10];       //! TOF informations

  // PHOS related track information 
  Float_t fPHOSpos[3]; // position localised by PHOS in global coordinate system
  Float_t fPHOSsignal; // energy measured by PHOS
  Float_t fPHOSr[AliPID::kSPECIESN]; // PID information from PHOS

  // EMCAL related track information 
  Float_t fEMCALpos[3]; //position localised by EMCAL in global coordinate system
  Float_t fEMCALsignal; // energy measured by EMCAL
  Float_t fEMCALr[AliPID::kSPECIESN]; // PID information from EMCAL

  // HMPID related track information
  Float_t fRICHsignal;     // detector's PID signal (beta for RICH)
  Float_t fRICHr[AliPID::kSPECIES];// "detector response probabilities" (for the PID)
  	
  ClassDef(AliESDtrack,11)  //ESDtrack 
};

#endif 

