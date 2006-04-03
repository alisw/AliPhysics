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
#include "AliExternalTrackParam.h"
#include "AliPID.h"
#include <TVector3.h>

class AliESDVertex;
class AliKalmanTrack;
class AliTrackPointArray;

const Int_t kNPlane = 6;

class AliESDtrack : public AliExternalTrackParam {
public:
  AliESDtrack();
  AliESDtrack(const AliESDtrack& track);
  virtual ~AliESDtrack();
  void MakeMiniESDtrack();
  void SetID(Int_t id) { fID =id;}
  Int_t GetID(){ return fID;}
  void SetStatus(ULong_t flags) {fFlags|=flags;}
  void ResetStatus(ULong_t flags) {fFlags&=~flags;}
  Bool_t UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags);
  void SetIntegratedLength(Double_t l) {fTrackLength=l;}
  void SetIntegratedTimes(const Double_t *times);
  void SetESDpid(const Double_t *p);
  void GetESDpid(Double_t *p) const;
  
  ULong_t GetStatus() const {return fFlags;}
  Int_t GetLabel() const {return fLabel;}
  void SetLabel(Int_t label) {fLabel = label;}

  void GetExternalParameters(Double_t &x, Double_t p[5]) const;
  void GetExternalCovariance(Double_t cov[15]) const;

  Double_t GetIntegratedLength() const {return fTrackLength;}
  void GetIntegratedTimes(Double_t *times) const;
  Double_t GetMass() const;
  TVector3 P3() const {Double_t p[3]; GetPxPyPz(p); return TVector3(p[0],p[1],p[2]);} //running track momentum
  TVector3 X3() const {Double_t x[3]; GetXYZ(x); return TVector3(x[0],x[1],x[2]);}    //running track position 


  Bool_t GetConstrainedPxPyPz(Double_t *p) const {
    if (!fCp) return kFALSE;
    return fCp->GetPxPyPz(p);
  }
  Bool_t GetConstrainedXYZ(Double_t *r) const {
    if (!fCp) return kFALSE;
    return fCp->GetXYZ(r);
  }
  Bool_t GetConstrainedExternalParameters
              (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetConstrainedExternalCovariance(Double_t cov[15]) const;
  Double_t GetConstrainedChi2() const {return fCchi2;}


  Bool_t GetInnerPxPyPz(Double_t *p) const {
    if (!fIp) return kFALSE;
    return fIp->GetPxPyPz(p);
  }
  const AliExternalTrackParam * GetInnerParam() const { return fIp;}
  Bool_t GetInnerXYZ(Double_t *r) const {
    if (!fIp) return kFALSE;
    return fIp->GetXYZ(r);
  }
  Bool_t GetInnerExternalParameters
        (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetInnerExternalCovariance(Double_t cov[15]) const;
 
  const AliExternalTrackParam * GetOuterParam() const { return fOp;}
  Bool_t GetOuterPxPyPz(Double_t *p) const {
    if (!fOp) return kFALSE;
    return fOp->GetPxPyPz(p);
  }
  Bool_t GetOuterXYZ(Double_t *r) const {
    if (!fOp) return kFALSE;
    return fOp->GetXYZ(r);
  }
  Bool_t GetOuterExternalParameters
        (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetOuterExternalCovariance(Double_t cov[15]) const;


  Int_t GetNcls(Int_t idet) const;
  Int_t GetClusters(Int_t idet, UInt_t *idx) const;
 
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
  void SetTPCPointsF(UChar_t  findable){fTPCnclsF = findable;}
  Float_t GetTPCPoints(Int_t i){return fTPCPoints[i];}
  void SetKinkIndexes(Int_t points[3]) {for (Int_t i=0;i<3;i++) fKinkIndexes[i] = points[i];}
  void SetV0Indexes(Int_t points[3]) {for (Int_t i=0;i<3;i++) fV0Indexes[i] = points[i];}
  void SetTPCsignal(Float_t signal, Float_t sigma, UChar_t npoints){ fTPCsignal = signal; fTPCsignalS = sigma; fTPCsignalN = npoints;}
  Float_t GetTPCsignal() const {return fTPCsignal;}
  Float_t GetTPCchi2() const {return fTPCchi2;}
  Int_t GetTPCclusters(Int_t *idx) const;
  Float_t GetTPCdensity(Int_t row0, Int_t row1) const;
  Int_t GetTPCLabel() const {return fTPCLabel;}
  Int_t GetKinkIndex(Int_t i) const { return fKinkIndexes[i];}
  Int_t GetV0Index(Int_t i) const { return fV0Indexes[i];}
  const TBits& GetTPCClusterMap() const {return fTPCClusterMap;}
  
  void SetTRDpid(const Double_t *p);
  void     SetTRDQuality(Float_t quality){fTRDQuality=quality;}
  Float_t  GetTRDQuality()const {return fTRDQuality;}
  void     SetTRDBudget(Float_t budget){fTRDBudget=budget;}
  Float_t  GetTRDBudget()const {return fTRDBudget;}
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


  AliKalmanTrack * GetTRDtrack(){return fTRDtrack;}

  void SetTOFsignal(Double_t tof) {fTOFsignal=tof;}
  Float_t GetTOFsignal() const {return fTOFsignal;}
  void SetTOFsignalToT(Double_t ToT) {fTOFsignalToT=ToT;}
  Float_t GetTOFsignalToT() const {return fTOFsignalToT;}
  Float_t GetTOFchi2() const {return fTOFchi2;}
  void    SetTOFpid(const Double_t *p);
  void    SetTOFLabel(const Int_t *p);
  void    GetTOFpid(Double_t *p) const;
  void    GetTOFLabel(Int_t *p) const;
  void    GetTOFInfo(Float_t *info) const;
  void    SetTOFInfo(Float_t *info);
  Int_t   GetTOFCalChannel() const {return fTOFCalChannel;}
  UInt_t  GetTOFcluster() const {return fTOFindex;}
  void  SetTOFcluster(UInt_t index) {fTOFindex=index;}
  void  SetTOFCalChannel(Int_t index) {fTOFCalChannel=index;}
  
  void    SetRICHsignal(Double_t beta) {fRICHsignal=beta;}
  Float_t GetRICHsignal() const {return fRICHsignal;}
  void    SetRICHpid(const Double_t *p);
  void    GetRICHpid(Double_t *p) const;
  void    SetRICHchi2(Double_t chi2) {fRICHchi2=chi2;}
  Float_t GetRICHchi2() const {return fRICHchi2;}
  void    SetRICHcluster(UInt_t index) {fRICHindex=index;}
  UInt_t  GetRICHcluster() const {return fRICHindex;}
  void    SetRICHnclusters(Int_t n) {fRICHncls=n;}
  Int_t   GetRICHnclusters() const {return fRICHncls;}
  void    SetRICHthetaPhi(Double_t theta, Double_t phi) {
    fRICHtheta=theta; fRICHphi=phi;
  }
  void    GetRICHthetaPhi(Double_t &theta, Double_t &phi) const {
    theta=fRICHtheta; phi=fRICHphi;
  }
  void    SetRICHdxdy(Double_t dx, Double_t dy) {
    fRICHdx=dx; fRICHdy=dy;
  }
  void    GetRICHdxdy(Double_t &dx, Double_t &dy) const {
    dx=fRICHdx; dy=fRICHdy;
  }
  
 /*  void SetPHOSposition(const Double_t *pos)  { */
/*     fPHOSpos[0] = pos[0]; fPHOSpos[1]=pos[1]; fPHOSpos[2]=pos[2]; */
/*   } */
/*   void SetPHOSsignal(Double_t ene) {fPHOSsignal = ene; } */
/*   void SetPHOSpid(const Double_t *p); */
/*   void GetPHOSposition(Double_t *pos) const { */
/*     pos[0]=fPHOSpos[0]; pos[1]=fPHOSpos[1]; pos[2]=fPHOSpos[2]; */
/*   } */
/*   Float_t GetPHOSsignal() const {return fPHOSsignal;} */
/*   void GetPHOSpid(Double_t *p) const;   */

  Bool_t IsOn(Int_t mask) const {return (fFlags&mask)>0;}
  Bool_t IsRICH()  const {return fFlags&kRICHpid;}
  Bool_t IsPHOS()  const {return fFlags&kPHOSpid;}

  void   SetTrackPointArray(AliTrackPointArray *points) { fPoints = points; }
  AliTrackPointArray *GetTrackPointArray() const { return fPoints; }

  Bool_t 
    RelateToVertex(const AliESDVertex *vtx, Double_t b, Double_t maxd);
  void GetImpactParameters(Float_t &xy,Float_t &z) const {xy=fD; z=fZ;}
  void GetImpactParameters(Float_t p[2], Float_t cov[3]) const {
    p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  }
  virtual void Print(Option_t * opt) const ; 

  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kPHOSpid=0x10000, kRICHpid=0x20000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
  }; 
protected:
  
  //AliESDtrack & operator=(const AliESDtrack & );

  ULong_t   fFlags;         // Reconstruction status flags 
  Int_t     fLabel;         // Track label
  Int_t     fID;            // Unique ID of the track
  Float_t   fTrackLength;   // Track length
  Float_t   fD;             // Impact parameter in XY plane
  Float_t   fZ;             // Impact parameter in Z
  Float_t   fCdd,fCdz,fCzz; // Covariance matrix of the impact parameters 
  Float_t   fTrackTime[AliPID::kSPECIES]; // TOFs estimated by the tracking
  Float_t   fR[AliPID::kSPECIES]; // combined "detector response probability"

  Int_t   fStopVertex;  // Index of the stop vertex

//Track parameters constrained to the primary vertex
  AliExternalTrackParam *fCp; 
  Double_t fCchi2; //chi2 at the primary vertex

//Track parameters at the inner wall of the TPC
  AliExternalTrackParam *fIp;

//Track parameters at the inner wall of the TRD 
  AliExternalTrackParam *fOp;

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
  UShort_t fTPCnclsF;      // number of findable clusters in the TPC
  Int_t  fTPCindex[180];  //! indices of the assigned TPC clusters
  TBits   fTPCClusterMap;  // Map of clusters, one bit per padrow; 1 if has a cluster on given padrow
  Float_t fTPCsignal;      // detector's PID signal
  UShort_t fTPCsignalN;      // number of points used for dEdx
  Float_t  fTPCsignalS;    // RMS of dEdx measurement
  Float_t fTPCr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTPCLabel;       // label according TPC
  Float_t fTPCPoints[4];   // TPC points -first, max. dens, last and max density
  Int_t   fKinkIndexes[3]; // array of indexes of posible kink candidates 
  Int_t   fV0Indexes[3]; // array of indexes of posible kink candidates 

  // TRD related track information
  Float_t fTRDchi2;        // chi2 in the TRD
  Int_t   fTRDncls;        // number of clusters assigned in the TRD
  Int_t   fTRDncls0;       // number of clusters assigned in the TRD before first material cross
  UInt_t  fTRDindex[180];   //! indices of the assigned TRD clusters
  Float_t fTRDsignal;      // detector's PID signal
  Float_t fTRDsignals[kNPlane];  // TRD signals from all six planes
  Int_t fTRDTimBin[kNPlane];     // Time bin of Max cluster from all six planes
  Float_t fTRDr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTRDLabel;       // label according TRD
  Float_t fTRDQuality;     //trd quality factor for TOF
  Float_t fTRDBudget;     //trd material budget
  AliKalmanTrack * fTRDtrack; //! OWNER: pointer to the TRD track -- currently for debug purpose

  // TOF related track information
  Float_t fTOFchi2;        // chi2 in the TOF
  UInt_t  fTOFindex;       // index of the assigned TOF cluster
  Int_t   fTOFCalChannel; // Channel Index of the TOF Signal 
  Float_t fTOFsignal;      // detector's PID signal
  Float_t fTOFsignalToT;   // detector's ToT signal
  Float_t fTOFr[AliPID::kSPECIES]; // "detector response probabilities" (for the PID)
  Int_t   fTOFLabel[3];       // TOF label 
  Float_t fTOFInfo[10];       //! TOF informations

  // PHOS related track information 
  //  Float_t fPHOSpos[3]; // position localised by PHOS in global coordinate system
  // Float_t fPHOSsignal; // energy measured by PHOS
  //Float_t fPHOSr[AliPID::kSPECIESN]; // PID information from PHOS

  // HMPID related track information
  Float_t fRICHchi2;       // chi2 in the RICH
  Int_t   fRICHncls;       // number of photon clusters
  UInt_t  fRICHindex;      // index of the assigned MIP cluster
  Float_t fRICHsignal;     // RICH PID signal
  Float_t fRICHr[AliPID::kSPECIES];// "detector response probabilities" (for the PID)
  Float_t fRICHtheta;      // theta of the track extrapolated to the RICH
  Float_t fRICHphi;        // phi of the track extrapolated to the RICH
  Float_t fRICHdx;         // x of the track impact minus x of the MIP
  Float_t fRICHdy;         // y of the track impact minus y of the MIP

  AliTrackPointArray *fPoints; // Array which contains the track space points in the global frame

  ClassDef(AliESDtrack,24)  //ESDtrack 
};

#endif 

