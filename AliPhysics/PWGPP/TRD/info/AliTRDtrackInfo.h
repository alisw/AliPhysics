#ifndef ALITRDTRACKINFO_H
#define ALITRDTRACKINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrackInfo.h 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIPID_H
#include "AliPID.h"
#endif


template <typename Value> class TVectorT;
typedef class TVectorT<Double_t> TVectorD;
class TGeoHMatrix;
class AliTRDseedV1;
class AliTRDtrackV1;
class AliTrackReference;
class AliTrackPointArray;
class AliExternalTrackParam;
class AliTRDtrackInfo : public TObject{
public:
  enum{
    kNTrackRefs = 12
  };
  class AliESDinfo{
    friend class AliTRDtrackInfo;  // Friend class
  public:
    enum EESDinfo{
       kV0       = 0
      ,kElectron = 1
      ,kPion     = 2
      ,kProton   = 3
    };
    AliESDinfo();
    AliESDinfo(const AliESDinfo &esd);
    virtual ~AliESDinfo();
    AliESDinfo& operator=(const AliESDinfo &esd);
    void Delete(const Option_t *);
    Float_t     Eta() const                      { return fEta;}
    Int_t       GetId() const                    { return fId;}
    ULong_t     GetStatus() const                { return fStatus;}
    Int_t       GetKinkIndex() const             { return fKinkIndex;}
    Short_t     GetTOFbc() const                 { return fTOFbc;}
    Float_t     GetTOFbeta() const               { return fTOFbeta;}
    UShort_t    GetTPCncls() const               { return fTPCncls;}
    Float_t     GetTPCdedx() const               { return fTPCdedx;}
    AliTrackPointArray* GetTrackPointArray() const      { return fTPArray; }
    UChar_t     GetPidQuality() const            { return fTRDpidQuality;}
    Int_t       GetNSlices() const               { return fTRDnSlices;}
    Double32_t* GetSliceIter() const             { return fTRDslices;}
    const Double32_t* GetResponseIter() const    { return &fTRDr[0];}
    AliExternalTrackParam* GetOuterParam() const { return fOP;}
    AliExternalTrackParam* GetITSoutParam() const { return fITSout;}
    AliExternalTrackParam* GetTPCoutParam() const { return fTPCout;}
    const Int_t* GetV0pid() const                { return &fTRDv0pid[0];}
    Int_t       GetV0pid(Int_t i) const          { return fTRDv0pid[i];}
    Bool_t      HasV0() const                    { return TESTBIT(fSteer, kV0);}
    Bool_t      IsElectron() const               { return TESTBIT(fSteer, kElectron);}
    Bool_t      IsPion() const                   { return TESTBIT(fSteer, kPion);}
    Bool_t      IsProton() const                 { return TESTBIT(fSteer, kProton);}
    Float_t     Phi() const                      { return fPhi;}
    Float_t     Pt() const                       { return fPt;}
    void        SetElectron(Bool_t e=kTRUE)      { e?SETBIT(fSteer, kElectron):CLRBIT(fSteer, kElectron);}
    void        SetPion(Bool_t e=kTRUE)          { e?SETBIT(fSteer, kPion):CLRBIT(fSteer, kPion);}
    void        SetProton(Bool_t e=kTRUE)        { e?SETBIT(fSteer, kProton):CLRBIT(fSteer, kProton);}
    void        SetV0(Bool_t v=kTRUE)            { v?SETBIT(fSteer, kV0):CLRBIT(fSteer, kV0);}
    void        SetOuterParam(const AliExternalTrackParam *op);
    void        SetITSoutParam(const AliExternalTrackParam *op);
    void        SetTPCoutParam(const AliExternalTrackParam *op);
    void        SetTrackPointArray(const AliTrackPointArray *tps);

  protected:
    UChar_t     fSteer;                   // bit map
    Int_t       fId;                      // ESD track id
    ULong_t     fStatus;                  // ESD track status
    Int_t       fKinkIndex;               // ESD kink index
    UShort_t    fTPCncls;                 // Number of Clusters inside TPC
    Float_t     fTPCdedx;                 // dEdx in TPC
    Float_t     fTOFbeta;                 // time of flight in TOF
    Short_t     fTOFbc;                   // TOF bunch crossing index
    Double32_t  fTRDr[AliPID::kSPECIES];  // TRD radial position
    UChar_t     fTRDpidQuality;           // TRD PID quality
    Int_t       fTRDnSlices;              // number of slices used for PID
    Float_t     fPt;                      // p_t at TRD entrance
    Float_t     fPhi;                     // phi at TRD entrance
    Float_t     fEta;                     // eta at TRD entrance
    Double32_t *fTRDslices;               //[fTRDnSlices]
    AliExternalTrackParam *fOP;           // outer track param
    AliExternalTrackParam *fTPCout;       // outer TPC param
    AliExternalTrackParam *fITSout;       // outer ITS param
    AliTrackPointArray *fTPArray;         // track point array to be used for alignment
    Int_t  fTRDv0pid[AliPID::kSPECIES];   // PID from v0s

    ClassDef(AliESDinfo, 9)     // ESD info related to TRD
  };

  class AliMCinfo{
  friend class AliTRDtrackInfo;
  public:
    //typedef AliTrackReference (const* constAliTrackReference);
    AliMCinfo();
    AliMCinfo(const AliMCinfo &mc);
    virtual ~AliMCinfo();
    AliMCinfo& operator=(const AliMCinfo &mc);
    Float_t GetEta() const                                  { return fEta;}
    Int_t   GetLabel() const                                { return fLabel;}
    Int_t   GetTRDlabel() const                             { return fTRDlabel;}
    Int_t   GetNTrackRefs() const                           { return fNTrackRefs;}
    Int_t   GetPDG() const                                  { return fPDG;}
    Int_t   GetPID() const ;
    Float_t GetPhi() const                                  { return fPhi;}
    Float_t GetPt() const                                   { return fPt;}
    Bool_t GetDirections(Float_t x0, Int_t chg, TGeoHMatrix *m, Double_t dir[10], UChar_t &s) const;
    AliTrackReference const* GetTrackRef(Int_t ref=0) const { return fTrackRefs[ref];}
    static Double_t GetKalmanStep()                         { return fgKalmanStep;}
    static Bool_t IsKalmanUpdate()                          { return fgKalmanUpdate;}
    Bool_t   PropagateKalman(
        TVectorD *x, TVectorD *y, TVectorD *z,
        TVectorD *dx, TVectorD *dy, TVectorD *dz,
        TVectorD *pt, TVectorD *dpt, TVectorD *budget, TVectorD *c, Double_t mass=-1) const;
    static void SetKalmanStep(Double_t s)                   { fgKalmanStep = s;}
    static void SetKalmanUpdate(Bool_t s=kTRUE)             { fgKalmanUpdate = s;}
  protected:
    Int_t   fLabel;               // ESD label
    Int_t   fTRDlabel;            // TRD label
    Int_t   fPDG;                 // particle code
    Int_t   fNTrackRefs;    	    // number of track refs
    Float_t fEta;                 // particle eta (vertex)
    Float_t fPhi;                 // particle phi (vertex)
    Float_t fPt;                  // particle pt (vertex)
    static Double_t fgKalmanStep; // Kalman step propagation
    static Bool_t fgKalmanUpdate; // Kalman update with TRD tracklets
    AliTrackReference  *fTrackRefs[kNTrackRefs];	// track refs array

    ClassDef(AliMCinfo, 4)      // MC info related to TRD
  };

  AliTRDtrackInfo();
  AliTRDtrackInfo(const AliTRDtrackInfo &trdInfo);
  ~AliTRDtrackInfo();

//  void               Clear(const Option_t *){}
  void               Delete(const Option_t *);
  AliTRDtrackInfo&   operator=(const AliTRDtrackInfo &trdInfo);
  void               AddTrackRef(const AliTrackReference *trackRef);
  Int_t              GetTrackId() const               { return fESD.fId;}
  const AliESDinfo*  GetESDinfo() const               { return &fESD; }
  const AliMCinfo*   GetMCinfo() const                { return fMC; }
  Int_t              GetNumberOfClusters() const;
  Int_t              GetNumberOfClustersRefit() const { return fNClusters;}
  Int_t              GetNTracklets() const;
  Int_t              GetNTrackRefs() const            { return fMC ? fMC->fNTrackRefs:0;} 
  Int_t              GetLabel() const                 { return fMC ? fMC->fLabel:0; }
  Int_t              GetKinkIndex() const             { return fESD.fKinkIndex;}
  Short_t            GetTOFbc() const                 { return fESD.fTOFbc;}
  UShort_t           GetTPCncls() const               { return fESD.fTPCncls;}
  Int_t              GetPDG() const                   { return fMC ? fMC->fPDG : 0; }
  Int_t              GetPID() const                   { return fMC ? fMC->GetPID() : -1; }
  ULong_t            GetStatus() const                { return fESD.fStatus;}
  AliTRDtrackV1*     GetTrack() const                 { return fTRDtrack; }
  AliTrackReference* GetTrackRef(Int_t entry) const;
  AliTrackReference* GetTrackRef(const AliTRDseedV1* const tracklet) const;

  Bool_t             IsCurved() const                 { return TestBit(kCurv);}
  Bool_t             IsPrimary() const                { return TestBit(kPrim);}
  Bool_t             HasESDtrack() const              { return ((fTRDtrack != 0x0) ||(fESD.fOP != 0));}
  Bool_t             HasMCinfo() const                { return (Bool_t)fMC; }

  void               SetElectron(Bool_t c = kTRUE)    { fESD.SetElectron(c);}
  void               SetPion(Bool_t c = kTRUE)        { fESD.SetPion(c);}
  void               SetProton(Bool_t c = kTRUE)      { fESD.SetProton(c);}
  void               SetCurved(Bool_t curv = kTRUE)   { SetBit(kCurv, curv);}
  void               SetLabel(Int_t lab)              { if(fMC) fMC->fLabel = lab; }
  void               SetTRDlabel(Int_t lab)           { if(fMC) fMC->fTRDlabel = lab; }
  void               SetNumberOfClustersRefit(Int_t n){fNClusters = n;}
  inline void        SetMC();
  void               SetMCeta(Float_t eta)            { if(fMC) fMC->fEta = eta; }
  void               SetMCphi(Float_t phi)            { if(fMC) fMC->fPhi = phi; }
  void               SetMCpt(Float_t pt)              { if(fMC) fMC->fPt = pt; }
  void               SetPDG(Int_t pdg)                { if(fMC) fMC->fPDG = pdg; }
  void               SetPrimary(Bool_t prim = kTRUE)  {SetBit(kPrim, prim);}
  void               SetOuterParam(const AliExternalTrackParam *op)  {fESD.SetOuterParam(op);}
  void               SetITSoutParam(const AliExternalTrackParam *op) {fESD.SetITSoutParam(op);}
  void               SetTPCoutParam(const AliExternalTrackParam *op) {fESD.SetTPCoutParam(op);}
  void               SetStatus(ULong_t stat)          { fESD.fStatus = stat;}
  void               SetKinkIndex(Int_t kinkIndex)    { fESD.fKinkIndex = kinkIndex;}
  void               SetTOFbc(Int_t bc)               { fESD.fTOFbc = bc;}
  void               SetTOFbeta(Float_t b)            { fESD.fTOFbeta = b;}
  void               SetTPCncls(UShort_t TPCncls)     { fESD.fTPCncls = TPCncls;}
  void               SetTPCdedx(Float_t dedx)         { fESD.fTPCdedx = dedx;}
  void               SetTrackId(Int_t id)             { fESD.fId = id;}
  void               SetTrack(const AliTRDtrackV1 *track);
  void               SetTrackPointArray(const AliTrackPointArray *tps) {fESD.SetTrackPointArray(tps);}
  void               SetESDpidQuality(UChar_t q)      { fESD.fTRDpidQuality = q;}
  void               SetSlices(Int_t n, Double32_t *s);
  inline void        SetESDpid(Double_t *);
  void               SetESDeta(Float_t eta)           { fESD.fEta = eta;}
  void               SetESDphi(Float_t phi)           { fESD.fPhi = phi;}
  void               SetESDpt(Float_t pt)             { fESD.fPt = pt;}
  inline void        SetV0pid(Int_t *);
  void               SetV0(Bool_t v0=kTRUE)           { fESD.SetV0(v0);}
  
private:
    enum{
      kCurv = 14,
      kPrim = 15
  };
  // this 2 data members have to go to ESD header.
  Int_t              fNClusters;     	// Numer of clusters from refit
  AliTRDtrackV1      *fTRDtrack; 	    // tracklets data array
  AliMCinfo          *fMC;            // MC extract for TRD
  AliESDinfo         fESD;            // ESD extract for TRD

  ClassDef(AliTRDtrackInfo, 4)        // TRD track info
};


//________________________________________________________
inline void AliTRDtrackInfo::SetMC()
{
  if(!fMC) fMC = new AliMCinfo();
}

//________________________________________________________
inline void AliTRDtrackInfo::SetESDpid(Double_t * const r)
{ 
  for(Int_t is = AliPID::kSPECIES; is--;) fESD.fTRDr[is] = r[is];
}

//________________________________________________________
inline void AliTRDtrackInfo::SetV0pid(Int_t * const r)
{ 
  for(Int_t is = AliPID::kSPECIES; is--;) fESD.fTRDv0pid[is] = r[is];
}

#endif
