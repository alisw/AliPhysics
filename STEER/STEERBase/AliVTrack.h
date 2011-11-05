#ifndef AliVTrack_H
#define AliVTrack_H
/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//-------------------------------------------------------------------------
//     base class for ESD and AOD tracks
//     Author: A. Dainese
//-------------------------------------------------------------------------

#include <TBits.h>

#include "AliVParticle.h"

class AliVVertex;
class AliExternalTrackParam;
class AliTPCdEdxInfo;
 
class AliVTrack: public AliVParticle {

public:
  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kTOFmismatch=0x100000,
    kHMPIDout=0x10000,kHMPIDpid=0x20000,
    kEMCALmatch=0x40000,
    kPHOSmatch=0x200000,
    kTRDbackup =0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000,
    kGlobalMerge=0x08000000,
    kITSpureSA=0x10000000,
    kMultInV0 =0x2000000,    //BIT(25): assumed to be belong to V0 in multiplicity estimates
    kMultSec  =0x4000000,     //BIT(26): assumed to be secondary (due to the DCA) in multiplicity estimates
    kEmbedded =0x8000000     // BIT(27), 1<<27: Is a track that has been embedded into the event
  };
  enum {
    kTRDnPlanes = 6,
    kEMCALNoMatch = -4096,
    kTOFBCNA = -100
  };

  AliVTrack() { }
  virtual ~AliVTrack() { }
  AliVTrack(const AliVTrack& vTrack); 
  AliVTrack& operator=(const AliVTrack& vTrack);

  virtual Int_t    GetID() const = 0;
  virtual UChar_t  GetITSClusterMap() const = 0;
  virtual Float_t  GetTPCClusterInfo(Int_t /*nNeighbours*/, Int_t /*type*/, Int_t /*row0*/=0, Int_t /*row1*/=159) const {return 0.;}
  virtual AliTPCdEdxInfo * GetTPCdEdxInfo() {return 0x0;}
  virtual UShort_t GetTPCNcls() const { return 0;}
  virtual UShort_t GetTPCNclsF() const { return 0;}
  virtual Double_t GetTRDslice(Int_t /*plane*/, Int_t /*slice*/) const { return -1.; }
  virtual UChar_t  GetTRDntrackletsPID() const { return 0;}
  
  virtual Int_t GetEMCALcluster()     const {return -1;}
  virtual void SetEMCALcluster(Int_t)       {;}
  virtual Bool_t IsEMCAL()            const {return kFALSE;}

  virtual Int_t GetPHOScluster()      const {return -1;}
  virtual void SetPHOScluster(Int_t)        {;}
  virtual Bool_t IsPHOS()             const {return kFALSE;}
  
  //pid info
  virtual Double_t  GetITSsignal()       const {return 0.;}
  virtual Double_t  GetTPCsignal()       const {return 0.;}
  virtual UShort_t  GetTPCsignalN()      const {return 0 ;}
  virtual Double_t  GetTPCmomentum()     const {return 0.;}
  virtual Double_t  GetTOFsignal()       const {return 0.;}
  virtual Double_t  GetHMPIDsignal()     const {return 0.;}

  virtual void      GetIntegratedTimes(Double_t */*times*/) const { return; }
  virtual Double_t  GetTRDmomentum(Int_t /*plane*/, Double_t */*sp*/=0x0) const {return 0.;}
  virtual void      GetHMPIDpid(Double_t */*p*/) const {;}
    
  virtual ULong_t  GetStatus() const = 0;
  virtual Bool_t   GetXYZ(Double_t *p) const = 0;
  virtual Double_t GetBz() const;
  virtual void     GetBxByBz(Double_t b[3]) const;
  virtual Bool_t   GetCovarianceXYZPxPyPz(Double_t cv[21]) const = 0;
  virtual Bool_t   PropagateToDCA(const AliVVertex *vtx,Double_t b,Double_t maxd,Double_t dz[2],Double_t covar[3]) = 0;
  virtual const    AliExternalTrackParam * GetOuterParam() const { return NULL; }
  virtual Int_t    GetNcls(Int_t /*idet*/) const { return 0; }
  virtual Bool_t   GetPxPyPz(Double_t */*p*/) const { return kFALSE; }
  virtual void     SetID(Short_t /*id*/) {;}
  virtual Int_t    GetTOFBunchCrossing(Double_t = 0) const { return kTOFBCNA;}

  ClassDef(AliVTrack,1)  // base class for tracks
};

#endif
