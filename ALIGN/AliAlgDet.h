#ifndef ALIALGDET_H
#define ALIALGDET_H

#include <TNamed.h>
#include <TObjArray.h>
#include <stdio.h>
#include "AliAlgAux.h"
#include "AliESDtrack.h"
class AliAlgTrack;
class AliAlgPoint;
class AliAlgSens;
class AliAlgVol;
class AliAlgSteer;
class AliTrackPointArray;
class AliExternalTrackParam;
class TH1;

/*--------------------------------------------------------
  Base class for detector: wrapper for set of volumes
  -------------------------------------------------------*/

// Author: ruben.shahoyan@cern.ch

class AliAlgDet : public TNamed
{
 public:
  enum {kInitGeomDone=BIT(14),kInitDOFsDone=BIT(15)};
  //
  AliAlgDet();
  AliAlgDet(const char* name, const char* title="");
  virtual ~AliAlgDet();
  Int_t   GetDetID()                            const {return GetUniqueID();}
  void    SetDetID(UInt_t tp);
  //
  virtual void  CacheReferenceOCDB();
  virtual void  AcknowledgeNewRun(Int_t run);
  virtual void  UpdateL2GRecoMatrices();
  //
  Int_t   VolID2SID(Int_t vid)                  const;
  Int_t   SID2VolID(Int_t sid)                  const {return sid<GetNSensors() ? fSID2VolID[sid] : -1;} //todo
  Int_t   GetNSensors()                         const {return fSensors.GetEntriesFast();}
  Int_t   GetNVolumes()                         const {return fVolumes.GetEntriesFast();}
  Int_t   GetVolIDMin()                         const {return fVolIDMin;}
  Int_t   GetVolIDMax()                         const {return fVolIDMax;}
  Bool_t  SensorOfDetector(Int_t vid)           const {return vid>=fVolIDMin && vid<=fVolIDMax;}
  void    SetAddError(double y, double z);
  const   Double_t* GetAddError()               const {return fAddError;} 
  //
  Int_t   GetNPoints()                          const {return fNPoints;}
  //
  void        SetAlgSteer(AliAlgSteer* s)             {fAlgSteer = s;}
  AliAlgSens* GetSensor(Int_t id)               const {return (AliAlgSens*)fSensors.UncheckedAt(id);}
  AliAlgSens* GetSensorByVolId(Int_t vid)       const {int sid=VolID2SID(vid); return sid<0 ? 0:GetSensor(sid);}
  AliAlgSens* GetSensor(const char* symname)    const {return (AliAlgSens*)fSensors.FindObject(symname);}
  AliAlgVol*  GetVolume(Int_t id)               const {return (AliAlgVol*)fVolumes.UncheckedAt(id);}
  AliAlgVol*  GetVolume(const char* symname)    const {return (AliAlgVol*)fVolumes.FindObject(symname);}
  //
  Bool_t      OwnsDOFID(Int_t id)               const;
  AliAlgVol*  GetVolOfDOFID(Int_t id)           const;
  //
  virtual Int_t InitGeom();
  virtual Int_t AssignDOFs();
  virtual void  InitDOFs();
  virtual void  Terminate(TH1* hdof=0);
  virtual void  AddVolume(AliAlgVol* vol);
  virtual void  DefineVolumes();
  virtual void  DefineMatrices();
  virtual void  Print(const Option_t *opt="")    const;
  virtual Int_t ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack,Bool_t inv=kFALSE);
  virtual AliAlgPoint* TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trp);
  virtual void  UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const;
  virtual void  SetUseErrorParam(Int_t v=0);
  Int_t         GetUseErrorParam()                   const {return fUseErrorParam;}
  //
  virtual Bool_t AcceptTrack(const AliESDtrack* trc,Int_t trtype) const = 0;
  Bool_t         CheckFlags(const AliESDtrack* trc,Int_t trtype) const;
  //
  virtual AliAlgPoint* GetPointFromPool();
  virtual void ResetPool();
  //
  void      SetInitGeomDone()                             {SetBit(kInitGeomDone);}
  Bool_t    GetInitGeomDone()                       const {return TestBit(kInitGeomDone);}
  //
  void      SetInitDOFsDone()                             {SetBit(kInitDOFsDone);}
  Bool_t    GetInitDOFsDone()                       const {return TestBit(kInitDOFsDone);}
  //
  Int_t     GetNDOFs()                              const {return fNDOFs;}
  //
  void      SetDisabled(Int_t tp,Bool_t v)                {fDisabled[tp]=v;SetObligatory(tp,!v);}
  void      SetDisabled()                                 {SetDisabledColl();SetDisabledCosm();}
  void      SetDisabledColl(Bool_t v=kTRUE)               {SetDisabled(AliAlgAux::kColl,v);}
  void      SetDisabledCosm(Bool_t v=kTRUE)               {SetDisabled(AliAlgAux::kCosm,v);}
  Bool_t    IsDisabled(Int_t tp)                    const {return fDisabled[tp];}
  Bool_t    IsDisabled()                            const {return IsDisabledColl()&&IsDisabledCosm();}
  Bool_t    IsDisabledColl()                        const {return IsDisabled(AliAlgAux::kColl);}
  Bool_t    IsDisabledCosm()                        const {return IsDisabled(AliAlgAux::kCosm);}
  //
  void      SetTrackFlagSel(Int_t tp,ULong_t f)           {fTrackFlagSel[tp] = f;}
  void      SetTrackFlagSelColl(ULong_t f)                {SetTrackFlagSel(AliAlgAux::kColl,f);}
  void      SetTrackFlagSelCosm(ULong_t f)                {SetTrackFlagSel(AliAlgAux::kCosm,f);}
  ULong_t   GetTrackFlagSel(Int_t tp)               const {return fTrackFlagSel[tp];}
  ULong_t   GetTrackFlagSelColl()                   const {return GetTrackFlagSel(AliAlgAux::kColl);}
  ULong_t   GetTrackFlagSelCosm()                   const {return GetTrackFlagSel(AliAlgAux::kCosm);}
  //
  void      SetNPointsSel(Int_t tp,Int_t n)               {fNPointsSel[tp] = n;}
  void      SetNPointsSelColl(Int_t n)                    {SetNPointsSel(AliAlgAux::kColl,n);}
  void      SetNPointsSelCosm(Int_t n)                    {SetNPointsSel(AliAlgAux::kCosm,n);}
  Int_t     GetNPointsSel(Int_t tp)                 const {return fNPointsSel[tp];}
  Int_t     GetNPointsSelColl()                     const {return GetNPointsSel(AliAlgAux::kColl);}
  Int_t     GetNPointsSelCosm()                     const {return GetNPointsSel(AliAlgAux::kCosm);}
  //
  Bool_t    IsObligatory(Int_t tp)                  const {return fObligatory[tp];}
  Bool_t    IsObligatoryColl()                      const {return IsObligatory(AliAlgAux::kColl);}
  Bool_t    IsObligatoryCosm()                      const {return IsObligatory(AliAlgAux::kCosm);}
  void      SetObligatory(Int_t tp,Bool_t v=kTRUE);
  void      SetObligatoryColl(Bool_t v=kTRUE)             {SetObligatory(AliAlgAux::kColl,v);}
  void      SetObligatoryCosm(Bool_t v=kTRUE)             {SetObligatory(AliAlgAux::kCosm,v);}
  //
  void      AddAutoConstraints()                    const;
  virtual void      WritePedeInfo(FILE* parOut,const Option_t *opt="") const;
  virtual void      WriteCalibrationResults()       const;
  virtual void      WriteAlignmentResults()         const;
  //
 protected:
  void     SortSensors();
  //
  // ------- dummies ---------
  AliAlgDet(const AliAlgDet&);
  AliAlgDet& operator=(const AliAlgDet&);
  //
 protected:
  //
  Int_t     fNDOFs;                      // number of DOFs free
  Int_t     fVolIDMin;                   // min volID for this detector (for sensors only)
  Int_t     fVolIDMax;                   // max volID for this detector (for sensors only)
  Int_t     fNSensors;                   // number of sensors (i.e. volID's)
  Int_t*    fSID2VolID;                  //[fNSensors] table of conversion from VolID to sid
  Int_t     fNProcPoints;                // total number of points processed
  //
  // Track selection
  Bool_t    fDisabled[AliAlgAux::kNTrackTypes];      // detector disabled/enabled in the track
  Bool_t    fObligatory[AliAlgAux::kNTrackTypes];    // detector must be present in the track
  ULong_t   fTrackFlagSel[AliAlgAux::kNTrackTypes];  // flag for track selection
  Int_t     fNPointsSel[AliAlgAux::kNTrackTypes];    // min number of points to require                 
  //
  Int_t     fUseErrorParam;          // signal that points need to be updated using track info, 0 - no
  Double_t  fAddError[2];            // additional error increment for measurement
  TObjArray fSensors;                // all sensors of the detector
  TObjArray fVolumes;                // all volumes of the detector  
  //
  // this is transient info
  Int_t     fNPoints;                //! number of points from this detector
  Int_t     fPoolNPoints;            //! number of points in the pool
  Int_t     fPoolFreePointID;        //! id of the last free point in the pool
  TObjArray fPointsPool;             //! pool of aligment points
  //
  AliAlgSteer* fAlgSteer;            // pointer to alignment steering object
  //
  ClassDef(AliAlgDet,1);             // base class for detector global alignment
};

//_____________________________________________________
inline Bool_t AliAlgDet::CheckFlags(const AliESDtrack* trc,Int_t trtype) const 
{
  // check if flags are ok
  return (trc->GetStatus()&fTrackFlagSel[trtype]) == fTrackFlagSel[trtype];
}

#endif
