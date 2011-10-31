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
 *                                                                           *
 * The Get*Label() getters return the label of the associated MC particle.   *
 * The absolute value of this label is the index of the particle within the  *
 * MC stack. If the label is negative, this track was assigned a certain     *
 * number of clusters that did not in fact belong to this track.             *
 *****************************************************************************/

#include <TBits.h>
#include "AliExternalTrackParam.h"
#include "AliVTrack.h"
#include "AliPID.h"
#include "AliESDfriendTrack.h"
#include "AliTPCdEdxInfo.h"

class TParticle;
class AliESDVertex;
class AliESDEvent;
class AliKalmanTrack;
class AliTrackPointArray;
class TPolyMarker3D;


class AliESDtrack : public AliExternalTrackParam {
public:
  //
  enum {kNITSchi2Std=3};
  //
  AliESDtrack();
  AliESDtrack(const AliESDtrack& track);
  AliESDtrack(const AliVTrack* track);
  AliESDtrack(TParticle * part);
  virtual ~AliESDtrack();
  virtual void Copy(TObject &obj) const;
  const AliESDfriendTrack *GetFriendTrack() const {return fFriendTrack;}
  void SetFriendTrack(const AliESDfriendTrack *t) {
    delete fFriendTrack; fFriendTrack=new AliESDfriendTrack(*t);
  }
  void ReleaseESDfriendTrack() { delete fFriendTrack;  fFriendTrack=0; }
  void AddCalibObject(TObject * object);     // add calib object to the list
  TObject *  GetCalibObject(Int_t index);    // return calib objct at given position
  void MakeMiniESDtrack();
  void SetID(Short_t id) { fID =id;}
  Int_t GetID() const { return fID;}
  void SetVertexID(Char_t id) { fVertexID=id;}
  Char_t GetVertexID() const { return fVertexID;}
  void SetStatus(ULong_t flags) {fFlags|=flags;}
  void ResetStatus(ULong_t flags) {fFlags&=~flags;}
  Bool_t UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags);
  void SetIntegratedLength(Double_t l) {fTrackLength=l;}
  void SetIntegratedTimes(const Double_t *times);
  void SetESDpid(const Double_t *p);
  void GetESDpid(Double_t *p) const;
  virtual const Double_t *PID() const { return fR; }

  Bool_t IsOn(Int_t mask) const {return (fFlags&mask)>0;}
  ULong_t GetStatus() const {return fFlags;}
  Int_t GetLabel() const {return fLabel;}
  void SetLabel(Int_t label) {fLabel = label;}

  void GetExternalParameters(Double_t &x, Double_t p[5]) const;
  void GetExternalCovariance(Double_t cov[15]) const;

  Double_t GetIntegratedLength() const {return fTrackLength;}
  void GetIntegratedTimes(Double_t *times) const;
  Int_t    GetPID()  const;
  Int_t    GetTOFBunchCrossing(Double_t b=0) const;
  Double_t GetMass() const {return AliPID::ParticleMass(GetPID());}
  Double_t M() const;
  Double_t E() const;
  Double_t Y() const;

  Bool_t GetConstrainedPxPyPz(Double_t *p) const {
    if (!fCp) return kFALSE;
    return fCp->GetPxPyPz(p);
  }
  Bool_t GetConstrainedXYZ(Double_t *r) const {
    if (!fCp) return kFALSE;
    return fCp->GetXYZ(r);
  }
  const AliExternalTrackParam *GetConstrainedParam() const {return fCp;}
  Bool_t GetConstrainedExternalParameters
              (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetConstrainedExternalCovariance(Double_t cov[15]) const;
  Double_t GetConstrainedChi2() const {return fCchi2;}
  Double_t GetChi2TPCConstrainedVsGlobal(const AliESDVertex* vtx) const;
  //
  
  // global track chi2
  void SetGlobalChi2(Double_t chi2) {fGlobalChi2 = chi2;}
  Double_t GetGlobalChi2() const {return fGlobalChi2;}

  Bool_t GetInnerPxPyPz(Double_t *p) const {
    if (!fIp) return kFALSE;
    return fIp->GetPxPyPz(p);
  }
  const AliExternalTrackParam * GetInnerParam() const { return fIp;}
  const AliExternalTrackParam * GetTPCInnerParam() const {return fTPCInner;}
  Bool_t FillTPCOnlyTrack(AliESDtrack &track);
  Bool_t GetInnerXYZ(Double_t *r) const {
    if (!fIp) return kFALSE;
    return fIp->GetXYZ(r);
  }
  Bool_t GetInnerExternalParameters
        (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetInnerExternalCovariance(Double_t cov[15]) const;
 
  void SetOuterParam(const AliExternalTrackParam *p, ULong_t flags);

  void SetOuterHmpParam(const AliExternalTrackParam *p, ULong_t flags);

  const AliExternalTrackParam * GetOuterParam() const { return fOp;}

  const AliExternalTrackParam * GetOuterHmpParam() const { return fHMPIDp;}
  
  Bool_t GetOuterPxPyPz(Double_t *p) const {
    if (!fOp) return kFALSE;
    return fOp->GetPxPyPz(p);
  }
  Bool_t GetOuterHmpPxPyPz(Double_t *p) const {
    if (!fHMPIDp) return kFALSE;
    return fHMPIDp->GetPxPyPz(p);
  }
  
  Bool_t GetOuterXYZ(Double_t *r) const {
    if (!fOp) return kFALSE;
    return fOp->GetXYZ(r);
  }
    Bool_t GetOuterHmpXYZ(Double_t *r) const {
    if (!fHMPIDp) return kFALSE;
    return fHMPIDp->GetXYZ(r);
  }

  Bool_t GetOuterExternalParameters
        (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetOuterExternalCovariance(Double_t cov[15]) const;

  Bool_t GetOuterHmpExternalParameters
        (Double_t &alpha, Double_t &x, Double_t p[5]) const;
  Bool_t GetOuterHmpExternalCovariance(Double_t cov[15]) const;

  
  Int_t GetNcls(Int_t idet) const;
  Int_t GetClusters(Int_t idet, Int_t *idx) const;
 
  void    SetITSpid(const Double_t *p);
  void    GetITSpid(Double_t *p) const;

  Double_t GetITSsignal() const {return fITSsignal;}
  void    SetITSdEdxSamples(const Double_t s[4]);
  void    GetITSdEdxSamples(Double_t *s) const;

  Double_t GetITSchi2() const {return fITSchi2;}
  Double_t GetITSchi2Std(Int_t step) const {return (step>-1&&step<kNITSchi2Std) ? fITSchi2Std[step] : -1;}
  void     SetITSchi2Std(Double_t chi2, Int_t step)  { if (step>-1&&step<kNITSchi2Std) fITSchi2Std[step] = chi2;}
  Char_t   GetITSclusters(Int_t *idx) const;
  UChar_t GetITSClusterMap() const {return fITSClusterMap;}
  UChar_t GetITSSharedMap() const {return fITSSharedMap;}
  void    SetITSSharedFlag(int lr) {fITSSharedMap |= 0x1<<lr;}
  Bool_t  GetITSFakeFlag()   const {return (fITSSharedMap&BIT(7))!=0;}
  void    SetITSFakeFlag(Bool_t v=kTRUE)  {if (v) fITSSharedMap|=BIT(7); else fITSSharedMap&=~BIT(7);}  
  void    SetITSSharedMap(UChar_t map) {fITSSharedMap=map;}
  void    SetITSModuleIndex(Int_t ilayer,Int_t idx) {fITSModule[ilayer]=idx;}
  Int_t   GetITSModuleIndex(Int_t ilayer) const {return fITSModule[ilayer];}
  Bool_t  GetITSModuleIndexInfo(Int_t ilayer,Int_t &idet,Int_t &status,
				Float_t &xloc,Float_t &zloc) const;
  Int_t   GetITSLabel() const {return fITSLabel;}
  void    SetITStrack(AliKalmanTrack * track){
    if (fFriendTrack) fFriendTrack->SetITStrack(track);
  }
  AliKalmanTrack *GetITStrack(){
    return fFriendTrack!=NULL?fFriendTrack->GetITStrack():NULL;
  }
  Bool_t  HasPointOnITSLayer(Int_t i) const {return TESTBIT(fITSClusterMap,i);}
  Bool_t  HasSharedPointOnITSLayer(Int_t i) const {return TESTBIT(fITSSharedMap,i);}

  void    SetTPCpid(const Double_t *p);
  void    GetTPCpid(Double_t *p) const;
  void    SetTPCPoints(Float_t points[4]){
     for (Int_t i=0;i<4;i++) fTPCPoints[i]=points[i];
  }
  void    SetTPCPointsF(UChar_t  findable){fTPCnclsF = findable;}
  void    SetTPCPointsFIter1(UChar_t  findable){fTPCnclsFIter1 = findable;}
  UShort_t   GetTPCNcls() const { return fTPCncls;}
  UShort_t   GetTPCNclsF() const { return fTPCnclsF;}
  UShort_t   GetTPCNclsIter1() const { return fTPCnclsIter1;}
  UShort_t   GetTPCNclsFIter1() const { return fTPCnclsFIter1;}
  UShort_t   GetTPCnclsS(Int_t i0=0,Int_t i1=159) const;
  UShort_t   GetTPCncls(Int_t row0=0,Int_t row1=159) const;
  Double_t GetTPCPoints(Int_t i) const {return fTPCPoints[i];}
  void    SetKinkIndexes(Int_t points[3]) {
     for (Int_t i=0;i<3;i++) fKinkIndexes[i] = points[i];
  }
  void    SetV0Indexes(Int_t points[3]) {
     for (Int_t i=0;i<3;i++) fV0Indexes[i] = points[i];
  }
  void    SetTPCsignal(Float_t signal, Float_t sigma, UChar_t npoints){ 
     fTPCsignal = signal; fTPCsignalS = sigma; fTPCsignalN = npoints;
  }
  void    SetTPCdEdxInfo(AliTPCdEdxInfo * dEdxInfo){ fTPCdEdxInfo = dEdxInfo;}
  AliTPCdEdxInfo * GetTPCdEdxInfo(){return fTPCdEdxInfo;}
  Double_t GetTPCsignal() const {return fTPCsignal;}
  Double_t GetTPCsignalSigma() const {return fTPCsignalS;}
  UShort_t GetTPCsignalN() const {return fTPCsignalN;}
  Double_t GetTPCmomentum() const {return fIp?fIp->GetP():GetP();}
  Double_t GetTPCchi2() const {return fTPCchi2;}
  Double_t GetTPCchi2Iter1() const {return fTPCchi2Iter1;}
  UShort_t GetTPCclusters(Int_t *idx) const;
  Double_t GetTPCdensity(Int_t row0, Int_t row1) const;
  Int_t   GetTPCLabel() const {return fTPCLabel;}
  Int_t   GetKinkIndex(Int_t i) const { return fKinkIndexes[i];}
  Int_t   GetV0Index(Int_t i) const { return fV0Indexes[i];}
  const TBits& GetTPCFitMap() const {return fTPCFitMap;}
  const TBits& GetTPCClusterMap() const {return fTPCClusterMap;}
  const TBits& GetTPCSharedMap() const {return fTPCSharedMap;}
  void    SetTPCFitMap(const TBits &amap) {fTPCFitMap = amap;}
  void    SetTPCClusterMap(const TBits &amap) {fTPCClusterMap = amap;}
  void    SetTPCSharedMap(const TBits &amap) {fTPCSharedMap = amap;}
  Float_t GetTPCClusterInfo(Int_t nNeighbours=3, Int_t type=0, Int_t row0=0, Int_t row1=159) const;
  Float_t GetTPCCrossedRows() const;
  
  void    SetTRDpid(const Double_t *p);
  void    SetTRDsignal(Double_t sig) {fTRDsignal = sig;}
	  
// A.Bercuci
  void    SetTRDntracklets(UChar_t q){fTRDntracklets = q;}
  UChar_t GetTRDntracklets() const {return (fTRDntracklets>>3)&7;}
  UChar_t GetTRDntrackletsPID() const {return fTRDntracklets&7;}
  // TEMPORARY alias asked by the HFE group to allow 
  // reading of the v4-16-Release data with TRUNK related software (A.Bercuci@Apr 30th 09) 
  UChar_t GetTRDpidQuality() const {return GetTRDntrackletsPID();}
// end A.Bercuci
  
  void     SetNumberOfTRDslices(Int_t n);
  Int_t    GetNumberOfTRDslices() const;
  void     SetTRDslice(Double_t q, Int_t plane, Int_t slice);
  void     SetTRDmomentum(Double_t p, Int_t plane, Double_t *sp=0x0);
  Double_t GetTRDslice(Int_t plane, Int_t slice=-1) const;
  Double_t GetTRDmomentum(Int_t plane, Double_t *sp=0x0) const;
	
  void    SetTRDQuality(Float_t quality){fTRDQuality=quality;}
  Double_t GetTRDQuality()const {return fTRDQuality;}
  void    SetTRDBudget(Float_t budget){fTRDBudget=budget;}
  Double_t GetTRDBudget()const {return fTRDBudget;}

  void    SetTRDTimBin(Int_t timbin, Int_t i) {fTRDTimBin[i]=timbin;}
  void    GetTRDpid(Double_t *p) const;
  Double_t GetTRDsignal() const {return fTRDsignal;}

  Char_t   GetTRDTimBin(Int_t i) const {return fTRDTimBin[i];}
  Double_t GetTRDchi2() const {return fTRDchi2;}
  UChar_t   GetTRDclusters(Int_t *idx) const;
  UChar_t   GetTRDncls() const {return fTRDncls;}
  UChar_t   GetTRDncls0() const {return fTRDncls0;}
  UChar_t   GetTRDtracklets(Int_t *idx) const;
  void    SetTRDpid(Int_t iSpecies, Float_t p);
  Double_t GetTRDpid(Int_t iSpecies) const;
  Int_t   GetTRDLabel() const {return fTRDLabel;}

  void    SetTRDtrack(AliKalmanTrack * track){
    if (fFriendTrack) fFriendTrack->SetTRDtrack(track);
  }
  AliKalmanTrack *GetTRDtrack(){
    return fFriendTrack!=NULL?fFriendTrack->GetTRDtrack():NULL;
  }

  void    SetTOFsignal(Double_t tof) {fTOFsignal=tof;}
  Double_t GetTOFsignal() const {return fTOFsignal;}
  void    SetTOFsignalToT(Double_t ToT) {fTOFsignalToT=ToT;}
  Double_t GetTOFsignalToT() const {return fTOFsignalToT;}
  void    SetTOFsignalRaw(Double_t tof) {fTOFsignalRaw=tof;}
  Double_t GetTOFsignalRaw() const {return fTOFsignalRaw;}
  void    SetTOFsignalDz(Double_t dz) {fTOFsignalDz=dz;}
  Double_t GetTOFsignalDz() const {return fTOFsignalDz;}
  void    SetTOFsignalDx(Double_t dx) {fTOFsignalDx=dx;}
  Double_t GetTOFsignalDx() const {return fTOFsignalDx;}
  void     SetTOFDeltaBC(Short_t deltaBC) {fTOFdeltaBC=deltaBC;};
  Short_t  GetTOFDeltaBC() const {return fTOFdeltaBC;}
  void     SetTOFL0L1(Short_t l0l1) {fTOFl0l1=l0l1;};
  Short_t  GetTOFL0L1() const {return fTOFl0l1;}
  Double_t GetTOFchi2() const {return fTOFchi2;}
  void    SetTOFpid(const Double_t *p);
  void    SetTOFLabel(const Int_t *p);
  void    GetTOFpid(Double_t *p) const;
  void    GetTOFLabel(Int_t *p) const;
  void    GetTOFInfo(Float_t *info) const;
  void    SetTOFInfo(Float_t *info);
  Int_t   GetTOFCalChannel() const {return fTOFCalChannel;}
  Int_t   GetTOFcluster() const {return fTOFindex;}
  void    SetTOFcluster(Int_t index) {fTOFindex=index;}
  void    SetTOFCalChannel(Int_t index) {fTOFCalChannel=index;}

// HMPID methodes +++++++++++++++++++++++++++++++++ (kir)
  void    SetHMPIDsignal(Double_t theta) {fHMPIDsignal=theta;}
  Double_t GetHMPIDsignal() const {return fHMPIDsignal - (Int_t)fHMPIDsignal;}
  Double_t GetHMPIDoccupancy() const {return (Int_t)fHMPIDsignal/10.0;}
  void    SetHMPIDpid(const Double_t *p);
  void    GetHMPIDpid(Double_t *p) const;  
  void    SetHMPIDchi2(Double_t chi2) {fHMPIDchi2=chi2;}
  Double_t GetHMPIDchi2() const {return fHMPIDchi2;}
  void    SetHMPIDcluIdx(Int_t ch,Int_t idx) {fHMPIDcluIdx=ch*1000000+idx;}
  Int_t   GetHMPIDcluIdx() const {return fHMPIDcluIdx;}
  void    SetHMPIDtrk(Float_t  x, Float_t  y, Float_t  th, Float_t  ph) {
     fHMPIDtrkX=x; fHMPIDtrkY=y; fHMPIDtrkTheta=th; fHMPIDtrkPhi=ph;
  }
  void    GetHMPIDtrk(Float_t &x, Float_t &y, Float_t &th, Float_t &ph) const {
     x=fHMPIDtrkX; y=fHMPIDtrkY; th=fHMPIDtrkTheta; ph=fHMPIDtrkPhi;
  }
  void    SetHMPIDmip(Float_t  x, Float_t  y, Int_t q, Int_t nph=0) {
     fHMPIDmipX=x; fHMPIDmipY=y; fHMPIDqn=1000000*nph+q;
  }
  void    GetHMPIDmip(Float_t &x,Float_t &y,Int_t &q,Int_t &nph) const {
     x=fHMPIDmipX; y=fHMPIDmipY; q=fHMPIDqn%1000000; nph=fHMPIDqn/1000000;
  }
  Bool_t  IsHMPID() const {return fFlags&kHMPIDpid;}
  Bool_t  IsPureITSStandalone() const {return fFlags&kITSpureSA;}
  Bool_t  IsMultPrimary() const {return !(fFlags&kMultSec);}
  Bool_t  IsMultSecondary() const {return (fFlags&kMultSec);}

  Int_t GetEMCALcluster() const {return fCaloIndex;}
  void SetEMCALcluster(Int_t index) {fCaloIndex=index;}
  Bool_t IsEMCAL() const {return fFlags&kEMCALmatch;}

  Int_t GetPHOScluster() const {return fCaloIndex;}
  void SetPHOScluster(Int_t index) {fCaloIndex=index;}
  Bool_t IsPHOS() const {return fFlags&kPHOSmatch;}
  Double_t GetPHOSdx()const{return fCaloDx ;}
  Double_t GetPHOSdz()const{return fCaloDz ;}
  void SetPHOSdxdz(Double_t dx, Double_t dz){fCaloDx=dx,fCaloDz=dz;}


  void SetTrackPointArray(AliTrackPointArray *points) {
    if (fFriendTrack) fFriendTrack->SetTrackPointArray(points);
  }
  const AliTrackPointArray *GetTrackPointArray() const {
    return fFriendTrack!=NULL?fFriendTrack->GetTrackPointArray():NULL; 
  }
  Bool_t RelateToVertexTPC(const AliESDVertex *vtx, Double_t b, Double_t maxd,
                           AliExternalTrackParam *cParam=0);
  Bool_t 
  RelateToVertexTPCBxByBz(const AliESDVertex *vtx, Double_t b[3],Double_t maxd,
                           AliExternalTrackParam *cParam=0);
  void GetImpactParametersTPC(Float_t &xy,Float_t &z) const {xy=fdTPC; z=fzTPC;}
  void GetImpactParametersTPC(Float_t p[2], Float_t cov[3]) const {
    p[0]=fdTPC; p[1]=fzTPC; cov[0]=fCddTPC; cov[1]=fCdzTPC; cov[2]=fCzzTPC;
  }
  Double_t GetConstrainedChi2TPC() const {return fCchi2TPC;}

  Bool_t RelateToVertex(const AliESDVertex *vtx, Double_t b, Double_t maxd,
                        AliExternalTrackParam *cParam=0);
  Bool_t 
  RelateToVertexBxByBz(const AliESDVertex *vtx, Double_t b[3], Double_t maxd,
                        AliExternalTrackParam *cParam=0);
  void GetImpactParameters(Float_t &xy,Float_t &z) const {xy=fD; z=fZ;}
  void GetImpactParameters(Float_t p[2], Float_t cov[3]) const {
    p[0]=fD; p[1]=fZ; cov[0]=fCdd; cov[1]=fCdz; cov[2]=fCzz;
  }
  virtual void Print(Option_t * opt) const ;
  AliESDEvent* GetESDEvent() const {return fESDEvent;}
  void         SetESDEvent(AliESDEvent* evt) {fESDEvent = evt;}  
  //
  // visualization (M. Ivanov)
  //
  void FillPolymarker(TPolyMarker3D *pol, Float_t magf, Float_t minR, Float_t maxR, Float_t stepR);

  //
  // online mode Matthias.Richter@cern.ch
  // in order to optimize AliESDtrack for usage in the online HLT,
  // some functionality is disabled
  // - creation of AliESDfriendTrack
  // - set lengt of bit fields fTPCClusterMap and fTPCSharedMap to 0
  static void OnlineMode(bool mode) {fgkOnlineMode=mode;}
  static bool OnlineMode() {return fgkOnlineMode;}
protected:
  
  AliExternalTrackParam *fCp; // Track parameters constrained to the primary vertex
  AliExternalTrackParam *fIp; // Track parameters estimated at the inner wall of TPC
  AliExternalTrackParam *fTPCInner; // Track parameters estimated at the inner wall of TPC using the TPC stand-alone 
  AliExternalTrackParam *fOp; // Track parameters estimated at the point of maximal radial coordinate reached during the tracking
  AliExternalTrackParam *fHMPIDp; // Track parameters at HMPID
  AliESDfriendTrack *fFriendTrack; //! All the complementary information

  TBits    fTPCFitMap;     // Map of clusters, one bit per padrow; 1 if has a cluster on given padrow which is used in the fit
  TBits    fTPCClusterMap; // Map of clusters, one bit per padrow; 1 if has a cluster on given padrow
  TBits    fTPCSharedMap;  // Map of clusters, one bit per padrow; 1 if has a shared cluster on given padrow



  ULong_t   fFlags;          // Reconstruction status flags 
  Int_t     fID;             // Unique ID of the track
  Int_t     fLabel;          // Track label
  Int_t     fITSLabel;       // label according ITS
  Int_t     fITSModule[12];  // modules crossed by the track in the ITS 
  Int_t     fTPCLabel;       // label according TPC
  Int_t     fTRDLabel;       // label according TRD
  Int_t     fTOFLabel[3];    // TOF label 
  Int_t     fTOFCalChannel;  // Channel Index of the TOF Signal 
  Int_t     fTOFindex;       // index of the assigned TOF cluster
  Int_t     fHMPIDqn;         // 1000000*number of photon clusters + QDC
  Int_t     fHMPIDcluIdx;     // 1000000*chamber id + cluster idx of the assigned MIP cluster
  Int_t     fCaloIndex;       // index of associated EMCAL/PHOS cluster (AliESDCaloCluster)


  Int_t     fKinkIndexes[3]; // array of indexes of posible kink candidates 
  Int_t     fV0Indexes[3];   // array of indexes of posible kink candidates 

  Double32_t   fR[AliPID::kSPECIES]; //[0.,0.,8] combined "detector response probability"
  Double32_t   fITSr[AliPID::kSPECIES]; //[0.,0.,8] "detector response probabilities" (for the PID)
  Double32_t   fTPCr[AliPID::kSPECIES]; //[0.,0.,8] "detector response probabilities" (for the PID)
  Double32_t   fTRDr[AliPID::kSPECIES]; //[0.,0.,8] "detector response probabilities" (for the PID)  
  Double32_t   fTOFr[AliPID::kSPECIES]; //[0.,0.,8] "detector response probabilities" (for the PID)
  Double32_t   fHMPIDr[AliPID::kSPECIES];//[0.,0.,8] "detector response probabilities" (for the PID)

  Double32_t fHMPIDtrkTheta;//[-2*pi,2*pi,16] theta of the track extrapolated to the HMPID, LORS
  // how much of this is needed?
  Double32_t fHMPIDtrkPhi;     //[-2*pi,2*pi,16] phi of the track extrapolated to the HMPID, LORS
  Double32_t fHMPIDsignal;  // HMPID PID signal (Theta ckov, rad)

  Double32_t   fTrackTime[AliPID::kSPECIES]; // TOFs estimated by the tracking
  Double32_t   fTrackLength;   // Track length

  Double32_t   fdTPC;          // TPC-only impact parameter in XY plane
  Double32_t   fzTPC;          // TPC-only impact parameter in Z
  Double32_t   fCddTPC,fCdzTPC,fCzzTPC; // Covariance matrix of the TPC-only impact parameters 
  Double32_t   fCchi2TPC;      // [0.,0.,8] TPC-only chi2 at the primary vertex

  Double32_t   fD;             // Impact parameter in XY plane
  Double32_t   fZ;             // Impact parameter in Z
  Double32_t   fCdd,fCdz,fCzz; // Covariance matrix of the impact parameters 
  Double32_t   fCchi2;          // [0.,0.,8] chi2 at the primary vertex

  Double32_t   fITSchi2Std[kNITSchi2Std];  // [0.,0.,8] standard chi2 in the ITS (with standard errors)
  Double32_t   fITSchi2;        // [0.,0.,8] chi2 in the ITS
  Double32_t   fTPCchi2;        // [0.,0.,8] chi2 in the TPC
  Double32_t   fTPCchi2Iter1;  // [0.,0.,8] chi2 in the TPC
  Double32_t   fTRDchi2;        // [0.,0.,8] chi2 in the TRD
  Double32_t   fTOFchi2;        // [0.,0.,8] chi2 in the TOF
  Double32_t fHMPIDchi2;        // [0.,0.,8] chi2 in the HMPID

  Double32_t fGlobalChi2;       // [0.,0.,8] chi2 of the global track

  Double32_t  fITSsignal;     // [0.,0.,10] detector's PID signal
  Double32_t  fITSdEdxSamples[4]; // [0.,0.,10] ITS dE/dx samples

  Double32_t  fTPCsignal;        // [0.,0.,10] detector's PID signal
  Double32_t  fTPCsignalS;       // [0.,0.,10] RMS of dEdx measurement
  AliTPCdEdxInfo * fTPCdEdxInfo; // object containing dE/dx information for different pad regions
  Double32_t  fTPCPoints[4];     // [0.,0.,10] TPC points -first, max. dens, last and max density

  Double32_t fTRDsignal;      // detector's PID signal
  Double32_t fTRDQuality;     // trd quality factor for TOF
  Double32_t fTRDBudget;      // trd material budget

  Double32_t fTOFsignal;      // detector's PID signal
  Double32_t fTOFsignalToT;   // detector's ToT signal
  Double32_t fTOFsignalRaw;   // detector's uncorrected time signal
  Double32_t fTOFsignalDz;    // local z  of track's impact on the TOF pad 
  Double32_t fTOFsignalDx;    // local x  of track's impact on the TOF pad 
  Double32_t fTOFInfo[10];    //! TOF informations
  Short_t    fTOFdeltaBC;     // detector's Delta Bunch Crossing correction
  Short_t    fTOFl0l1;        // detector's L0L1 latency correction

  Double32_t fCaloDx ;        // [0.,0.,8] distance to calorimeter cluster in calo plain (phi direction)
  Double32_t fCaloDz ;        // [0.,0.,8] distance to calorimeter cluster in calo plain (z direction)

  Double32_t fHMPIDtrkX;       // x of the track impact, LORS 
  Double32_t fHMPIDtrkY;       // y of the track impact, LORS 
  Double32_t fHMPIDmipX;       // x of the MIP in LORS
  Double32_t fHMPIDmipY;       // y of the MIP in LORS


  UShort_t fTPCncls;       // number of clusters assigned in the TPC
  UShort_t fTPCnclsF;      // number of findable clusters in the TPC
  UShort_t fTPCsignalN;    // number of points used for dEdx
  UShort_t fTPCnclsIter1;  // number of clusters assigned in the TPC - iteration 1
  UShort_t fTPCnclsFIter1; // number of findable clusters in the TPC - iteration 1

  Char_t  fITSncls;        // number of clusters assigned in the ITS
  UChar_t fITSClusterMap;  // map of clusters, one bit per a layer
  UChar_t fITSSharedMap;   // map of shared clusters, one bit per a layer
  UChar_t fTRDncls;        // number of clusters assigned in the TRD
  UChar_t fTRDncls0;       // number of clusters assigned in the TRD before first material cross
  UChar_t fTRDntracklets;  // number of TRD tracklets used for tracking/PID

  Int_t fTRDnSlices;     // number of slices used for PID in the TRD
  Double32_t *fTRDslices;  //[fTRDnSlices] 

  Char_t  fTRDTimBin[kTRDnPlanes];   // Time bin of Max cluster from all six planes
  Char_t  fVertexID; // ID of the primary vertex this track belongs to
  AliESDEvent*   fESDEvent; //!Pointer back to event to which the track belongs
  
  mutable Float_t fCacheNCrossedRows; //! Cache for the number of crossed rows
  mutable Float_t fCacheChi2TPCConstrainedVsGlobal; //! Cache for the chi2 of constrained TPC vs global track
  mutable const AliESDVertex* fCacheChi2TPCConstrainedVsGlobalVertex; //! Vertex for which the cache is valid
  
 private:
  static bool fgkOnlineMode; //! indicate the online mode to skip some of the functionality

  AliESDtrack & operator=(const AliESDtrack & );
  ClassDef(AliESDtrack,63)  //ESDtrack 
};



#endif 

