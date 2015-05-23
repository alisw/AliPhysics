/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliAlgDet.h"
#include "AliAlgSens.h"
#include "AliAlgDet.h"
#include "AliAlgSteer.h"
#include "AliAlgTrack.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliAlignObj.h"
#include "AliCDBId.h"
#include "AliExternalTrackParam.h"
#include <TString.h>
#include <TH1.h>
#include <stdio.h>

ClassImp(AliAlgDet)

using namespace AliAlgAux;

//____________________________________________
AliAlgDet::AliAlgDet()
:  fNDOFs(0)
  ,fVolIDMin(-1)
  ,fVolIDMax(-1)
  ,fNSensors(0)
  ,fSID2VolID(0)
  ,fNProcPoints(0)
  //
  ,fUseErrorParam(0)
  ,fSensors()
  ,fVolumes()
  //
  ,fNPoints(0)
  ,fPoolNPoints(0)
  ,fPoolFreePointID(0)
  ,fPointsPool()
  ,fAlgSteer(0)
{
  // def c-tor
  SetUniqueID(AliAlgSteer::kUndefined); // derived detectors must override this
  fAddError[0] = fAddError[1] = 0;
  //
  for (int i=0;i<kNTrackTypes;i++) {
    fDisabled[i] = kFALSE;
    fObligatory[i] = kFALSE;
    fTrackFlagSel[i] = 0;
    fNPointsSel[i] = 0;
  }
  //
}

//____________________________________________
AliAlgDet::~AliAlgDet()
{
  // d-tor
  fSensors.Clear(); // sensors are also attached as volumes, don't delete them here
  fVolumes.Delete(); // here all is deleted
  fPointsPool.Delete();
}


//____________________________________________
Int_t AliAlgDet::ProcessPoints(const AliESDtrack* esdTr, AliAlgTrack* algTrack, Bool_t inv)
{
  // Extract the points corresponding to this detector, recalibrate/realign them to the
  // level of the "starting point" for the alignment/calibration session.
  // If inv==true, the track propagates in direction of decreasing tracking X 
  // (i.e. upper leg of cosmic track)
  //
  const AliESDfriendTrack* trF(esdTr->GetFriendTrack());
  const AliTrackPointArray* trP(trF->GetTrackPointArray());
  //
  int np(trP->GetNPoints());
  int npSel(0);
  AliAlgPoint* apnt(0);
  for (int ip=0;ip<np;ip++) {
    if (!SensorOfDetector(trP->GetVolumeID()[ip])) continue;
    if (!(apnt=TrackPoint2AlgPoint(ip, trP))) continue;
    algTrack->AddPoint(apnt);
    if (inv) apnt->SetInvDir();
    npSel++;
    fNPoints++;
  }
  //
  return npSel;
}

//____________________________________________
AliAlgPoint* AliAlgDet::TrackPoint2AlgPoint(int pntId, const AliTrackPointArray* trpArr)
{
  // convert the pntId-th point to AliAlgPoint, detectors may override this method
  //
  // convert to detector tracking frame
  UShort_t vid = trpArr->GetVolumeID()[pntId];
  Int_t sid = VolID2SID(vid); // sensor index within the detector
  if (!sid<0) return 0;
  AliAlgSens* sens = GetSensor(sid);
  if (sens->GetSkip()) return 0;
  AliAlgPoint* pnt = GetPointFromPool();
  pnt->SetSensor(sens);
  //
  double tra[3],traId[3],loc[3],
    glo[3] = {trpArr->GetX()[pntId], trpArr->GetY()[pntId], trpArr->GetZ()[pntId]};
  const TGeoHMatrix& matL2Grec = sens->GetMatrixL2GReco(); // local to global matrix used for reconstruction
  //const TGeoHMatrix& matL2G    = sens->GetMatrixL2G();     // local to global orig matrix used as a reference 
  const TGeoHMatrix& matT2L    = sens->GetMatrixT2L();     // matrix for tracking to local frame translation
  //
  // undo reco-time alignment
  matL2Grec.MasterToLocal(glo,loc); // go to local frame using reco-time matrix 
  matT2L.MasterToLocal(loc,traId); // go to tracking frame 
  //
  sens->GetMatrixClAlg().LocalToMaster(traId,tra);   // apply alignment
  //
  if (!fUseErrorParam) {
    // convert error
    TGeoHMatrix hcov;
    Double_t hcovel[9];
    const Float_t *pntcov = trpArr->GetCov()+pntId*6; // 6 elements per error matrix
    hcovel[0] = double(pntcov[0]);
    hcovel[1] = double(pntcov[1]);
    hcovel[2] = double(pntcov[2]);
    hcovel[3] = double(pntcov[1]);
    hcovel[4] = double(pntcov[3]);
    hcovel[5] = double(pntcov[4]);
    hcovel[6] = double(pntcov[2]);
    hcovel[7] = double(pntcov[4]);
    hcovel[8] = double(pntcov[5]);
    hcov.SetRotation(hcovel);
    hcov.Multiply(&matL2Grec);                
    hcov.MultiplyLeft(&matL2Grec.Inverse());    // errors in local frame
    hcov.Multiply(&matT2L);
    hcov.MultiplyLeft(&matT2L.Inverse());       // errors in tracking frame
    //
    Double_t *hcovscl = hcov.GetRotationMatrix();
    const double *sysE = sens->GetAddError(); // additional syst error
    pnt->SetYZErrTracking(hcovscl[4]+sysE[0]*sysE[0],hcovscl[5],hcovscl[8]+sysE[1]*sysE[1]);
  }
  else { // errors will be calculated just befor using the point in the fit, using track info
    pnt->SetYZErrTracking(0,0,0);
    pnt->SetNeedUpdateFromTrack();
  }
  pnt->SetXYZTracking(tra[0],tra[1],tra[2]);
  pnt->SetAlphaSens(sens->GetAlpTracking());
  pnt->SetXSens(sens->GetXTracking());
  pnt->SetDetID(GetDetID());
  pnt->SetSID(sid);
  //
  pnt->SetContainsMeasurement();
  //
  pnt->Init();
  //
  return pnt;
  //
}

//_________________________________________________________
void AliAlgDet::AcknowledgeNewRun(Int_t run)
{
  // update parameters needed to process this run

  // detector should be able to undo alignment/calibration used during the reco
  UpdateL2GRecoMatrices();

}

//_________________________________________________________
void AliAlgDet::UpdateL2GRecoMatrices()
{
  // Update L2G matrices used for data reconstruction
  //
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBEntry* ent = man->Get(Form("%s/Align/Data",GetName()));
  const TClonesArray *algArr = (const TClonesArray*)ent->GetObject();
  //
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    // call init for root level volumes, they will take care of their children
    if (!vol->GetParent()) vol->UpdateL2GRecoMatrices(algArr,0);
  }
  //
}

//_________________________________________________________
void AliAlgDet::CacheReferenceOCDB()
{
  // if necessary, detector may fetch here some reference OCDB data
}


//_________________________________________________________
AliAlgPoint* AliAlgDet::GetPointFromPool()
{
  // fetch or create new free point from the pool.
  // detector may override this method to create its own points derived from AliAlgPoint
  //
  if (fPoolFreePointID>=fPoolNPoints) { // expand pool
    fPointsPool.AddAtAndExpand(new AliAlgPoint(), fPoolNPoints++);
  }
  //
  AliAlgPoint* pnt = (AliAlgPoint*) fPointsPool.UncheckedAt(fPoolFreePointID++);
  pnt->Clear();
  return pnt;
  //
}

//_________________________________________________________
void AliAlgDet::ResetPool()
{
  // declare pool free
  fPoolFreePointID = 0;
  fNPoints = 0;
}
 
//_________________________________________________________
void AliAlgDet::DefineVolumes()
{
  // dummy method
  AliError("This method must be implemented by specific detector");
}

//_________________________________________________________
void AliAlgDet::AddVolume(AliAlgVol* vol)
{
  // add volume
  if (GetVolume(vol->GetSymName())) {
    AliFatalF("Volume %s was already added to %s",vol->GetName(),GetName());
  }
  fVolumes.AddLast(vol);
  if (vol->IsSensor()) {
    fSensors.AddLast(vol);
    ((AliAlgSens*)vol)->SetDetector(this);
    Int_t vid = ((AliAlgSens*)vol)->GetVolID();
    if (fVolIDMin<0 || vid<fVolIDMin) fVolIDMin = vid;
    if (fVolIDMax<0 || vid>fVolIDMax) fVolIDMax = vid;
  }
  //
}

//_________________________________________________________
void AliAlgDet::DefineMatrices()
{
  // define transformation matrices. Detectors may override this method
  //
  TGeoHMatrix mtmp;
  //
  TIter next(&fVolumes);
  AliAlgVol* vol(0);
  while ( (vol=(AliAlgVol*)next()) ) {
    // modified global-local matrix
    vol->PrepareMatrixL2G();
    // ideal global-local matrix
    vol->PrepareMatrixL2GIdeal();
    //
  }
  // Now set tracking-local matrix (MUST be done after ALL L2G matrices are done!)
  // Attention: for sensor it is a real tracking matrix extracted from
  // the geometry but for container alignable volumes the tracking frame
  // is used for as the reference for the alignment parameters only,
  // see its definition in the AliAlgVol::PrepateMatrixT2L
  next.Reset();
  while ( (vol=(AliAlgVol*)next()) ) vol->PrepareMatrixT2L();
  //
}

//_________________________________________________________
void AliAlgDet::SortSensors()
{
  // build local tables for internal numbering
  fNSensors = fSensors.GetEntriesFast();
  if (!fNSensors) {
    AliWarning("No sensors defined");
    return;
  }
  fSensors.Sort();
  fSID2VolID = new Int_t[fNSensors]; // cash id's for fast binary search
  for (int i=0;i<fNSensors;i++) fSID2VolID[i] = GetSensor(i)->GetVolID();
  //
}

//_________________________________________________________
Int_t AliAlgDet::InitGeom()
{
  // define hiearchy, initialize matrices, return number of global parameters
  if (GetInitGeomDone()) return 0;
  //
  DefineVolumes();
  SortSensors();    // VolID's must be in increasing order
  DefineMatrices();
  //
  // calculate number of global parameters
  int nvol = GetNVolumes();
  fNDOFs = 0;
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    fNDOFs += vol->GetNDOFs();
  }
  //
  SetInitGeomDone();
  return fNDOFs;
}

//_________________________________________________________
Int_t AliAlgDet::AssignDOFs()
{
  // assign DOFs IDs, parameters
  //
  int gloCount0(fAlgSteer->GetNDOFs()), gloCount(fAlgSteer->GetNDOFs());
  Float_t* pars = fAlgSteer->GetGloParVal(); 
  Float_t* errs = fAlgSteer->GetGloParErr(); 
  Int_t*   labs = fAlgSteer->GetGloParLab();
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    // call init for root level volumes, they will take care of their children
    if (!vol->GetParent()) vol->AssignDOFs(gloCount,pars,errs,labs);
  }
  //
  if (fNDOFs != gloCount-gloCount0) AliFatalF("Mismatch between declared %d and initialized %d DOFs for %s",
					      fNDOFs,gloCount-gloCount0,GetName());
  
  return fNDOFs;
}

//_________________________________________________________
void AliAlgDet::InitDOFs()
{
  // initialize free parameters
  if (GetInitDOFsDone()) AliFatalF("Something is wrong, DOFs are already initialized for %s",GetName());
  //
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) GetVolume(iv)->InitDOFs();
  //
  SetInitDOFsDone();
  return;
}

//_________________________________________________________
Int_t AliAlgDet::VolID2SID(Int_t vid) const 
{
  // find SID corresponding to VolID
  int mn(0),mx(fNSensors-1);
  while (mx>=mn) {
    int md( (mx+mn)>>1 ), vids(GetSensor(md)->GetVolID());
    if (vid<vids)      mx = md-1;
    else if (vid>vids) mn = md+1;
    else return md;
  }
  return -1;
}

//____________________________________________
void AliAlgDet::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("\nDetector:%5s %5d volumes %5d sensors {VolID: %5d-%5d} Def.Sys.Err: %.4e %.4e | Stat:%d\n",
	 GetName(),GetNVolumes(),GetNSensors(),GetVolIDMin(),
	 GetVolIDMax(),fAddError[0],fAddError[1],fNProcPoints);
  //
  printf("Errors assignment: ");
  if (fUseErrorParam) printf("param %d\n",fUseErrorParam);
  else printf("from TrackPoints\n");
  //
  printf("Allowed    in Collisions: %7s | Cosmic: %7s\n",
	 IsDisabled(kColl)   ? "  NO ":" YES ",IsDisabled(kCosm)   ? "  NO ":" YES ");
  //
  printf("Obligatory in Collisions: %7s | Cosmic: %7s\n",
	 IsObligatory(kColl) ? " YES ":"  NO ",IsObligatory(kCosm) ? " YES ":"  NO ");
  //
  printf("Sel. flags in Collisions: 0x%05lx | Cosmic: 0x%05lx\n",
	 fTrackFlagSel[kColl],fTrackFlagSel[kCosm]);
  //
  printf("Min.points in Collisions: %7d | Cosmic: %7d\n",
	 fNPointsSel[kColl],fNPointsSel[kCosm]);
  //
  if (!(IsDisabledColl()&&IsDisabledCosm()) && opts.Contains("long")) 
    for (int iv=0;iv<GetNVolumes();iv++) GetVolume(iv)->Print(opt);
  //
}

//____________________________________________
void AliAlgDet::SetDetID(UInt_t tp)
{
  // assign type
  if (tp>=AliAlgSteer::kNDetectors) AliFatalF("Detector typeID %d exceeds allowed range %d:%d",
					      tp,0,AliAlgSteer::kNDetectors-1);
  SetUniqueID(tp);
}

//____________________________________________
void AliAlgDet::SetAddError(double sigy, double sigz)
{
  // add syst error to all sensors
  AliInfoF("Adding sys.error %.4e %.4e to all sensors",sigy,sigz);
  fAddError[0] = sigy;
  fAddError[1] = sigz;
  for (int isn=GetNSensors();isn--;) GetSensor(isn)->SetAddError(sigy,sigz);
  //
}

//____________________________________________
void AliAlgDet::SetUseErrorParam(Int_t v) 
{
  // set type of points error parameterization
  AliFatal("UpdatePointByTrackInfo is not implemented for this detector");
  //  
}

//____________________________________________
void AliAlgDet::UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const
{
  // update point using specific error parameterization
  AliFatal("If needed, this method has to be implemented for specific detector");
}

//____________________________________________
void AliAlgDet::SetObligatory(Int_t tp,Bool_t v)
{
  // mark detector presence obligatory in the track
  fObligatory[tp] = v;
  fAlgSteer->SetObligatoryDetector(GetDetID(),tp,v);
}

//______________________________________________________
void AliAlgDet::WritePedeInfo(FILE* parOut, const Option_t *opt) const
{
  // contribute to params and constraints template files for PEDE
  fprintf(parOut,"\n!!\t\tDetector:\t%s\tNDOFs: %d\n",GetName(),GetNDOFs());
  //
  // parameters
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {  // call for root level volumes, they will take care of their children
    AliAlgVol *vol = GetVolume(iv);
    if (!vol->GetParent()) vol->WritePedeInfo(parOut,opt);
  }
  //
}

//______________________________________________________
void AliAlgDet::WriteCalibrationResults() const
{
  // store calibration results
  WriteAlignmentResults();
  // 
  // eventually we may write other calibrations
}

//______________________________________________________
void AliAlgDet::WriteAlignmentResults() const
{
  // store updated alignment
  TClonesArray* arr = new TClonesArray("AliAlignObjParams",10);
  //
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    // call only for top level objects, they will take care of children
    if (!vol->GetParent()) vol->CreateAlignmentObjects(arr);
  }
  //
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible(fAlgSteer->GetOutCDBResponsible());
  md->SetComment(fAlgSteer->GetOutCDBResponsible());
  //
  AliCDBId id(Form("%s/Align/Data",GetName()),fAlgSteer->GetOutCDBRunMin(),fAlgSteer->GetOutCDBRunMax());
  man->Put(arr,id,md); 
  //
  delete arr;
}

//______________________________________________________
Bool_t AliAlgDet::OwnsDOFID(Int_t id) const
{
  // check if DOF ID belongs to this detector
  for (int iv=GetNVolumes();iv--;) {
    AliAlgVol* vol = GetVolume(iv); // check only top level volumes
    if (!vol->GetParent() && vol->OwnsDOFID(id)) return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________
AliAlgVol* AliAlgDet::GetVolOfDOFID(Int_t id) const
{
  // gets volume owning this DOF ID
  for (int iv=GetNVolumes();iv--;) {
    AliAlgVol* vol = GetVolume(iv);
    if (vol->GetParent()) continue; // check only top level volumes
    if ( (vol=vol->GetVolOfDOFID(id)) ) return vol;
  }
  return 0;
}

//______________________________________________________
void AliAlgDet::Terminate(TH1* h)
{
  // called at the end of processing
  if (IsDisabled()) return;
  int nvol = GetNVolumes();
  fNProcPoints = 0;
  for (int iv=0;iv<nvol;iv++) {
    AliAlgVol *vol = GetVolume(iv);
    // call init for root level volumes, they will take care of their children
    if (!vol->GetParent()) fNProcPoints += vol->FinalizeStat(h);
  }
}

//________________________________________
void AliAlgDet::AddAutoConstraints() const
{
  // adds automatic constraints
  int nvol = GetNVolumes();
  for (int iv=0;iv<nvol;iv++) {  // call for root level volumes, they will take care of their children
    AliAlgVol *vol = GetVolume(iv);
    if (!vol->GetParent()) vol->AddAutoConstraints((TObjArray*)fAlgSteer->GetConstraints());
  }
}
