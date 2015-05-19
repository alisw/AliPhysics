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

#include <stdio.h>
#include "AliAlgSens.h"
#include "AliAlgAux.h"
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliExternalTrackParam.h"
#include "AliAlgPoint.h"
#include "AliAlgDet.h"

ClassImp(AliAlgSens)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSens::AliAlgSens(const char* name,Int_t vid, Int_t iid) 
  : AliAlgVol(name,iid)
  ,fDet(0)
{
  // def c-tor
  SetVolID(vid);
  fAddError[0] = fAddError[1] = 0;
  fConstrChild = 0; // sensors don't have children
}

//_________________________________________________________
AliAlgSens::~AliAlgSens()
{
  // d-tor
}

//_________________________________________________________
void AliAlgSens::DPosTraDParGeomLOC(const double *tra, double* deriv) const
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor LOCAL frame 
  // parameters in TGeoHMatrix convention.
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  memset(delta,0,kNDOFGeom*sizeof(double));
  memset(deriv,0,kNDOFGeom*3*sizeof(double));
  //
  for (int ip=kNDOFGeom;ip--;) {
    //
    if (!IsFreeDOF(ip)) continue;
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    // variation matrix in tracking frame for variation in sensor LOCAL frame
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(tra,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(tra,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(tra,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta); 
    matMod.LocalToMaster(tra,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParGeomLOC(const double *tra, double* deriv, const AliAlgVol* parent) const
{
  // Jacobian of position in sensor tracking frame (tra) vs parent volume LOCAL frame parameters.
  // NO check of parentship is done!
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  // this is the matrix for transition from sensor to parent volume local frames: LOC=matRel*loc
  TGeoHMatrix matRel = parent->GetMatrixL2G().Inverse(); 
  matRel *= GetMatrixL2G();
  //
  memset(delta,0,kNDOFGeom*sizeof(double));
  memset(deriv,0,kNDOFGeom*3*sizeof(double));
  //
  for (int ip=kNDOFGeom;ip--;) {
    //
    if (!IsFreeDOF(ip)) continue;
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(tra,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(tra,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(tra,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodLOC(matMod, delta, matRel);
    matMod.LocalToMaster(tra,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParGeomTRA(const double *tra, double* deriv) const
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor TRACKING 
  // frame parameters in TGeoHMatrix convention, i.e. the modified parameter is
  // tra' = tau*tra
  //
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  memset(delta,0,kNDOFGeom*sizeof(double));
  memset(deriv,0,kNDOFGeom*3*sizeof(double));
  //
  for (int ip=kNDOFGeom;ip--;) {
    //
    if (!IsFreeDOF(ip)) continue;
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(tra,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(tra,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(tra,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta);
    matMod.LocalToMaster(tra,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParGeomTRA(const double *tra, double* deriv, const AliAlgVol* parent) const
{
  // Jacobian of position in sensor tracking frame (tra) vs sensor TRACKING 
  // frame parameters in TGeoHMatrix convention, i.e. the modified parameter is
  // tra' = tau*tra
  //
  // Result is stored in array deriv as linearized matrix 6x3 
  const double kDelta[kNDOFGeom]={0.1,0.1,0.1,0.5,0.5,0.5};
  double delta[kNDOFGeom],pos0[3],pos1[3],pos2[3],pos3[3];
  TGeoHMatrix matMod;
  //
  // 1st we need a matrix for transition between child and parent TRACKING frames
  // Let TRA,LOC are positions in tracking and local frame of parent, linked as LOC=T2L*TRA
  // and tra,loc are positions in tracking and local frame of child,  linked as loc=t2l*tra
  // The loc and LOC are linked as LOC=R*loc, where R = L2G^-1*l2g, with L2G and l2g 
  // local2global matrices for parent and child
  //
  // Then, TRA = T2L^-1*LOC = T2L^-1*R*loc = T2L^-1*R*t2l*tra
  // -> TRA = matRel*tra, with matRel = T2L^-1*L2G^-1 * l2g*t2l
  // Note that l2g*t2l are tracking to global matrices
  TGeoHMatrix matRel,t2gP;
  GetMatrixT2G(matRel);           // t2g matrix of child
  parent->GetMatrixT2G(t2gP);     // t2g matrix of parent
  matRel.MultiplyLeft(&t2gP.Inverse());
  //
  memset(delta,0,kNDOFGeom*sizeof(double));
  memset(deriv,0,kNDOFGeom*3*sizeof(double));
  //
  for (int ip=kNDOFGeom;ip--;) {
    //
    if (!IsFreeDOF(ip)) continue;
    //
    double var = kDelta[ip];
    delta[ip] -= var;
    GetDeltaT2LmodTRA(matMod,delta,matRel);
    matMod.LocalToMaster(tra,pos0);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta,matRel);
    matMod.LocalToMaster(tra,pos1);     // varied position in tracking frame
    //
    delta[ip] += var;
    GetDeltaT2LmodTRA(matMod,delta,matRel);
    matMod.LocalToMaster(tra,pos2);     // varied position in tracking frame
    //
    delta[ip] += 0.5*var;
    GetDeltaT2LmodTRA(matMod,delta,matRel);
    matMod.LocalToMaster(tra,pos3);     // varied position in tracking frame
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./var;
  }
  //
}

//_________________________________________________________
void AliAlgSens::DPosTraDParGeom(const double *tra, double* deriv, const AliAlgVol* parent) const
{
  // calculate point position derivatives in tracking frame of sensor
  // vs standard geometrical DOFs of its parent volume (if parent!=0) or sensor itself
  Frame_t frame = parent ? parent->GetVarFrame() : GetVarFrame();
  switch(frame) {
  case kLOC : parent ? DPosTraDParGeomLOC(tra,deriv,parent) : DPosTraDParGeomLOC(tra,deriv);
    break;
  case kTRA : parent ? DPosTraDParGeomTRA(tra,deriv,parent) : DPosTraDParGeomTRA(tra,deriv);
    break;
  default   : AliErrorF("Alignment frame %d is not implemented",parent->GetVarFrame()); 
    break;
  }
}

//__________________________________________________________________
void AliAlgSens::GetModifiedMatrixT2LmodLOC(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the sensitive module tracking2local matrix from its current T2L matrix 
  // by applying local delta of modification of LOCAL frame:
  // loc' = delta*loc = T2L'*tra = T2L'*T2L^-1*loc   ->  T2L' = delta*T2L
  Delta2Matrix(matMod, delta);
  matMod.Multiply(&GetMatrixT2L());
}


//__________________________________________________________________
void AliAlgSens::GetModifiedMatrixT2LmodTRA(TGeoHMatrix& matMod, const Double_t *delta) const
{
  // prepare the sensitive module tracking2local matrix from its current T2L matrix 
  // by applying local delta of modification of TRACKING frame:
  // loc' = T2L'*tra = T2L*delta*tra    ->  T2L' = T2L*delta
  Delta2Matrix(matMod, delta);
  matMod.MultiplyLeft(&GetMatrixT2L());
}

//__________________________________________________________________
void AliAlgSens::AddChild(AliAlgVol*)
{
  AliFatalF("Sensor volume cannot have childs: id=%d %s",GetVolID(),GetName());
}

//__________________________________________________________________
Int_t AliAlgSens::Compare(const TObject* b) const
{
  // compare VolIDs
  return GetUniqueID()<b->GetUniqueID() ? -1 : 1;
}

//__________________________________________________________________
void AliAlgSens::SetTrackingFrame()
{
  // define tracking frame of the sensor
  AliWarningF("Generic method called for %s",GetSymName());
  double tra[3]={0},loc[3],glo[3];
  const TGeoHMatrix &t2l = GetMatrixT2L();
  const double* t = t2l.GetTranslation();
  double r = TMath::Sqrt(t[0]*t[0]+t[1]*t[1]);
  // ITS defines tracking frame with origin in sensor, others at 0
  if (r>1) tra[0] = r;
  //
  t2l.LocalToMaster(tra,loc);
  GetMatrixL2GIdeal().LocalToMaster(loc,glo);
  fX = Sqrt(glo[0]*glo[0]+glo[1]*glo[1]);
  fAlp = ATan2(glo[1],glo[0]);
  AliAlgAux::BringToPiPM(fAlp);
}

//____________________________________________
void AliAlgSens::Print(const Option_t *opt) const
{
  // print info
  TString opts = opt;
  opts.ToLower();
  printf("Lev:%2d IntID:%7d %s VId:%6d X:%8.4f Alp:%+.4f | Err: %.4e %.4e | Used Points: %d\n",
	 CountParents(), GetInternalID(), GetSymName(), GetVolID(), fX, fAlp, 
	 fAddError[0],fAddError[1],fNProcPoints);
  printf("     DOFs: Tot: %d (offs: %5d) Free: %d  Geom: %d {",fNDOFs,fFirstParGloID,fNDOFFree,fNDOFGeomFree);
  for (int i=0;i<kNDOFGeom;i++) printf("%d",IsFreeDOF(i) ? 1:0); 
  printf("} in %s frame\n",fgkFrameName[fVarFrame]);
  //
  if (opts.Contains("mat")) { // print matrices
    printf("L2G ideal   : "); 
    GetMatrixL2GIdeal().Print();
    printf("L2G misalign: "); 
    GetMatrixL2G().Print();
    printf("L2G RecoTime: "); 
    GetMatrixL2GReco().Print();
    printf("T2L         : "); 
    GetMatrixT2L().Print();
  }
  //
}

//____________________________________________
void AliAlgSens::PrepareMatrixT2L()
{
  // extract from geometry T2L matrix
  const TGeoHMatrix* t2l = AliGeomManager::GetTracking2LocalMatrix(GetVolID());
  if (!t2l) {
    Print("long");
    AliFatalF("Failed to find T2L matrix for VID:%d %s",GetVolID(),GetSymName());
  }
  SetMatrixT2L(*t2l);  
  //
}

//____________________________________________
void AliAlgSens::UpdatePointByTrackInfo(AliAlgPoint* pnt, const AliExternalTrackParam* t) const
{
  fDet->UpdatePointByTrackInfo(pnt,t);
}
