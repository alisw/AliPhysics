#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliLog.h"
#include <TGeoMatrix.h>
#include <TString.h>

using namespace TMath;

ClassImp(AliITSUClusterPix)

AliITSUGeomTGeo* AliITSUClusterPix::fgGeom = 0;
UInt_t           AliITSUClusterPix::fgMode = 0;

//_____________________________________________________
AliITSUClusterPix::AliITSUClusterPix()
  : fCharge(0)
  , fNxNzN(0)
{
  // default constructor
}

//_____________________________________________________
AliITSUClusterPix::~AliITSUClusterPix()
{
  // default destructor
}

//_____________________________________________________
AliITSUClusterPix::AliITSUClusterPix(const AliITSUClusterPix& cluster) 
  :AliCluster(cluster)
  ,fCharge(cluster.fCharge)
  ,fNxNzN(cluster.fNxNzN)
{
  // copy constructor
}

//______________________________________________________________________________
AliITSUClusterPix& AliITSUClusterPix::operator=(const AliITSUClusterPix& cluster)
{
  // = op
  if(&cluster == this) return *this;
  fNxNzN = cluster.fNxNzN;
  fCharge = cluster.fCharge;
  TObject::operator=(cluster);
  AliCluster::operator=(cluster);
  return *this;
}

//______________________________________________________________________________
const TGeoHMatrix*  AliITSUClusterPix::GetTracking2LocalMatrix() const
{
  // get tracking to local matrix (sensor!!!)
  return (TGeoHMatrix*)fgGeom->GetMatrixT2L(GetVolumeId());
}

//______________________________________________________________________________
TGeoHMatrix* AliITSUClusterPix::GetMatrix(Bool_t ) const
{
  // get module matrix (sensor!)
  return (TGeoHMatrix*)fgGeom->GetMatrixSens(GetVolumeId());
}

//______________________________________________________________________________
void AliITSUClusterPix::Print(Option_t* option) const
{
  // Print cluster information.
  TString str = option; 
  str.ToLower();
  printf("Cl.in mod %5d, nx:%3d nz:%3d n:%d |Err^2:%.3e %.3e %+.3e |",GetVolumeId(),GetNx(),GetNz(),
	 GetNPix(),GetSigmaY2(),GetSigmaZ2(),GetSigmaYZ());
  printf("XYZ: (%+.4e %+.4e %+.4e ",GetX(),GetY(),GetZ());
  if      (IsFrameLoc()) printf("LOC)");
  else if (IsFrameGlo()) printf("GLO)");
  else if (IsFrameTrk()) printf("TRK)");
  if (str.Contains("glo") && !IsFrameGlo() && fgGeom) {
    Float_t g[3];
    GetGlobalXYZ(g);
    printf(" (%+.4e %+.4e %+.4e GLO)",g[0],g[1],g[2]);
  }
  printf(" MClb:");
  for (int i=0;i<3;i++) printf(" %5d",GetLabel(i));
  if (TestBit(kSplit)) printf(" Spl");
  printf("\n");
  //
}

//______________________________________________________________________________
Bool_t AliITSUClusterPix::GetGlobalXYZ(Float_t xyz[3]) const
{
  // Get the global coordinates of the cluster
  // All the needed information is taken only
  // from TGeo.
  if (IsFrameGlo()) {
    xyz[0] = GetX();
    xyz[1] = GetY();
    xyz[2] = GetZ();
  }
  //
  Double_t lxyz[3] = {0, 0, 0};
  if (IsFrameTrk()) {
    const TGeoHMatrix *mt = GetTracking2LocalMatrix();
    if (!mt) return kFALSE;
    Double_t txyz[3] = {GetX(), GetY(), GetZ()};
    mt->LocalToMaster(txyz,lxyz);
  }
  else {
    lxyz[0] = GetX(); lxyz[1] = GetY(); lxyz[2] = GetZ();
  }
  //
  TGeoHMatrix *ml = GetMatrix();
  if (!ml) return kFALSE;
  Double_t gxyz[3] = {0, 0, 0};
  ml->LocalToMaster(lxyz,gxyz);
  xyz[0] = gxyz[0]; xyz[1] = gxyz[1]; xyz[2] = gxyz[2];
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliITSUClusterPix::GetGlobalCov(Float_t cov[6]) const
{
  // Get the global covariance matrix of the cluster coordinates
  // All the needed information is taken only
  // from TGeo.
  // Note: regardless on in which frame the coordinates are, the errors are always in tracking frame
  //
  return AliCluster::GetGlobalCov(cov);
}

//______________________________________________________________________________
Bool_t AliITSUClusterPix::GetXRefPlane(Float_t &xref) const
{
  // Get the distance between the origin and the ref.plane.
  // All the needed information is taken only from TGeo.
  return AliCluster::GetXRefPlane(xref);
}

//______________________________________________________________________________
void AliITSUClusterPix::GoToFrameGlo()
{
  // convert to global frame
  if (IsFrameGlo()) return;
  double loc[3],glo[3];
  //
  if (IsFrameTrk()) {
    double curr[3]={GetX(),GetY(),GetZ()};
    GetTracking2LocalMatrix()->LocalToMaster(curr,loc);
    ResetBit(kFrameTrk);
  }
  else {
    loc[0] = GetX(); loc[1] = GetY(); loc[2] = GetZ();
    ResetBit(kFrameLoc);
  }
  GetMatrix()->LocalToMaster(loc,glo);
  SetX(glo[0]);  
  SetY(glo[1]); 
  SetZ(glo[2]);
  SetBit(kFrameGlo);
  //
}

//______________________________________________________________________________
void AliITSUClusterPix::GoToFrameLoc()
{
  // convert to local frame
  if (IsFrameLoc()) return;
  //
  double loc[3],glo[3];
  if (IsFrameTrk()) {
    double curr[3]={GetX(),GetY(),GetZ()};
    GetTracking2LocalMatrix()->LocalToMaster(curr,loc);
    ResetBit(kFrameTrk);
  }
  else {
    glo[0] = GetX(); glo[1] = GetY(); glo[2] = GetZ();
    GetMatrix()->MasterToLocal(glo,loc);
    ResetBit(kFrameLoc);
  }
  SetBit(kFrameLoc);
  SetX(loc[0]); 
  SetY(loc[1]); 
  SetZ(loc[2]);
  //
}

//______________________________________________________________________________
void AliITSUClusterPix::GetLocalXYZ(Float_t xyz[3]) const
{
  // get local coordinates
  if (IsFrameLoc()) {
    xyz[0] = GetX(); xyz[1] = 0; xyz[2] = GetZ();
    return;
  }
  double loc[3],glo[3];
  if (IsFrameTrk()) {
    double curr[3]={GetX(),GetY(),GetZ()};
    GetTracking2LocalMatrix()->LocalToMaster(curr,loc);
  }
  else {
    glo[0] = GetX(); glo[1] = GetY(); glo[2] = GetZ();
    GetMatrix()->MasterToLocal(glo,loc);
  }
  for (int i=3;i--;) xyz[i] = loc[i];
  //
}

//______________________________________________________________________________
void AliITSUClusterPix::GoToFrameTrk()
{
  // convert to tracking frame
  if (IsFrameTrk()) return;
  //
  double loc[3],trk[3];
  if (IsFrameGlo()) {
    double glo[3]={GetX(),GetY(),GetZ()};
    GetMatrix()->MasterToLocal(glo,loc);
    ResetBit(kFrameGlo);
  }
  else {
    loc[0] = GetX(); loc[1] = GetY(); loc[2] = GetZ();
    ResetBit(kFrameLoc);    
  }
  // now in local frame
  GetTracking2LocalMatrix()->MasterToLocal(loc,trk);
  SetBit(kFrameTrk);
  SetX(trk[0]);  
  SetY(trk[1]); 
  SetZ(trk[2]);
  //
}

//______________________________________________________________________________
void AliITSUClusterPix::GetTrackingXYZ(Float_t xyz[3]) const
{
  // convert to tracking frame
  if (IsFrameTrk()) {
    xyz[0] = GetX(); xyz[1] = GetY(); xyz[2] = GetZ();
    return;
  }
  //
  double loc[3],trk[3];
  if (IsFrameGlo()) {
    double glo[3]={GetX(),GetY(),GetZ()};
    GetMatrix()->MasterToLocal(glo,loc);
  }
  else {
    loc[0] = GetX(); loc[1] = GetY(); loc[2] = GetZ();
  }
  // now in local frame
  GetTracking2LocalMatrix()->MasterToLocal(loc,trk);
  for (int i=3;i--;) xyz[i] = trk[i];
  //
}

//______________________________________________________________________________
Int_t AliITSUClusterPix::Compare(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUClusterPix* px = (const AliITSUClusterPix*)obj;
  float xyz[3],xyz1[3];
  if (fgMode & kSortIdLocXZ) { // sorting in local frame
    if (GetVolumeId()==px->GetVolumeId()) {
      GetLocalXYZ(xyz);
      px->GetLocalXYZ(xyz1);
      if (xyz[0]<xyz1[0]) return -1; // sort in X
      if (xyz[0]>xyz1[0]) return  1;
      if (xyz[2]<xyz1[2]) return -1; // then in Z
      if (xyz[2]>xyz1[2]) return  1;
      return 0;
    }
    return int(GetVolumeId())-int(px->GetVolumeId());
  }
  if (fgMode & kSortIdTrkYZ) { // sorting in tracking frame
    if (GetVolumeId()==px->GetVolumeId()) {
      GetTrackingXYZ(xyz);
      px->GetTrackingXYZ(xyz1);
      if (xyz[1]<xyz1[1]) return -1; // sort in Y
      if (xyz[1]>xyz1[1]) return  1;
      if (xyz[2]<xyz1[2]) return -1; // then in Z
      if (xyz[2]>xyz1[2]) return  1;
      return 0;    
    }
    return int(GetVolumeId())-int(px->GetVolumeId());    
  }
  AliFatal(Form("Unknown modr for sorting: %d",fgMode));
  return 0;
}

//______________________________________________________________________________
Bool_t AliITSUClusterPix::IsEqual(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUClusterPix* px = (const AliITSUClusterPix*)obj;
  const Float_t kTol = 1e-5;
  float xyz[3],xyz1[3];
  if (fgMode & kSortIdLocXZ) { // sorting in local frame
    if (GetVolumeId()!=px->GetVolumeId()) return kFALSE;
    GetLocalXYZ(xyz);
    px->GetLocalXYZ(xyz1);
    return (Abs(xyz[0]-xyz1[0])<kTol && Abs(xyz[2]-xyz1[2])<kTol) ? kTRUE : kFALSE;
  }
  if (fgMode & kSortIdTrkYZ) { // sorting in tracking frame
    if (GetVolumeId()!=px->GetVolumeId()) return kFALSE;
    GetTrackingXYZ(xyz);
    px->GetTrackingXYZ(xyz1);
    return (Abs(xyz[1]-xyz1[1])<kTol && Abs(xyz[2]-xyz1[2])<kTol) ? kTRUE : kFALSE;
  }
  AliFatal(Form("Unknown modr for sorting: %d",fgMode));
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliITSUClusterPix::HasCommonTrack(const AliCluster* cl) const
{
  // check if clusters have common tracks
  int lbi,lbj;
  for (int i=0;i<3;i++) {
    if ((lbi=GetLabel(i))<0) break;
    for (int j=0;j<3;j++) {
      if ((lbj=cl->GetLabel(j))<0) break;
      if (lbi==lbj) return kTRUE;
    }
  }
  return kFALSE;
}
