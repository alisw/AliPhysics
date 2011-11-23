/************************************************************************** 
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. * 
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
 
/* $Id$    */ 
//----------------------------------------------------------------------------- 
/// \class AliITSAlignMille2Module
/// Alignment class for the ALICE ITS detector 
/// 
/// This class is used by AliITSAlignMille to build custom supermodules    
/// made of ITS sensitive modules. These supermodules are then aligned
/// 
/// Custom supermodules must have VolumeID > 14335
///
/// \author M. Lunardon  
//----------------------------------------------------------------------------- 
 
#include <TGeoManager.h> 
#include <TGeoMatrix.h> 
 
#include "AliITSAlignMille2Module.h" 
#include "AliITSgeomTGeo.h" 
#include "AliGeomManager.h" 
#include "AliAlignObjParams.h" 
#include "AliLog.h" 
#include "AliITSAlignMille2.h"
 
/// \cond CLASSIMP 
ClassImp(AliITSAlignMille2Module) 
/// \endcond 

#define CORHW_

AliAlignObjParams AliITSAlignMille2Module::fgTempAlignObj;
const Float_t AliITSAlignMille2Module::fgkDummyConstraint = 1e-2;//1.E3;
    
//-------------------------------------------------------------
AliITSAlignMille2Module::AliITSAlignMille2Module() : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fDetType(-1),
  fVolumeID(0),
  fNParTot(0),
  fNParFree(0),
  fParOffs(0),
  fNProcPoints(0),
  fParVals(0),
  fParErrs(0),
  fParCstr(0),
  fSensVolIndex(0),
  fSensVolVolumeID(0),
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL),
  fChildren(0)
{ 
  /// void constructor  
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix; 
  fSensVolIndex.Set(1);
  fSensVolVolumeID.Set(1);
  fSigmaFactor[0]=fSigmaFactor[1]=fSigmaFactor[2]=1.0;
} 

//-------------------------------------------------------------
AliITSAlignMille2Module::AliITSAlignMille2Module(Int_t index, UShort_t volid, const char* symname,
						 const TGeoHMatrix *m, Int_t nsv, const UShort_t *volidsv) : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fDetType(-1),  
  fVolumeID(0),
  fNParTot(0),
  fNParFree(0),
  fParOffs(kMaxParGeom),
  fNProcPoints(0),
  fParVals(0),
  fParErrs(0),
  fParCstr(0),
  fSensVolIndex(0),
  fSensVolVolumeID(0),  
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL),
  fChildren(0)
{ 
  /// void constructor  
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix; 
  fSigmaFactor[0]=fSigmaFactor[1]=fSigmaFactor[2]=1.0;
  for (int i=kMaxParGeom;i--;) fParOffs[i] = -1;
  if (Set(index,volid,symname,m,nsv,volidsv)) {
    AliInfo("Error in AliITSAlignMille2Module::Set() - initializing void supermodule...");
  }
  AssignDetType();
} 

//-------------------------------------------------------------
AliITSAlignMille2Module::AliITSAlignMille2Module(UShort_t volid) : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),    
  fDetType(-1),
  fVolumeID(0),
  fNParTot(0),
  fNParFree(0),
  fParOffs(kMaxParGeom),
  fNProcPoints(0),
  fParVals(0),
  fParErrs(0),
  fParCstr(0),  
  fSensVolIndex(0),
  fSensVolVolumeID(0),
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL),
  fChildren(0)
{ 
  /// simple constructor building a supermodule from a single sensitive volume 
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix;   
  // temporary align object, just use the rotation...
  fSensVolIndex.Set(1);
  fSensVolVolumeID.Set(1);
  fSigmaFactor[0]=fSigmaFactor[1]=fSigmaFactor[2]=1.0;
  for (int i=kMaxParGeom;i--;) fParOffs[i] = -1;
  //
  fIndex = GetIndexFromVolumeID(volid);  
  if (fIndex>=0 && gGeoManager) { // good sensitive module and geometry loaded
    SetName(AliGeomManager::SymName(volid));
    fVolumeID = volid;
    AddSensitiveVolume(volid);
    SetSensorsProvided(kTRUE);
    if (SensVolMatrix(volid, fMatrix))
       AliInfo("Matrix not defined");
  }
  else {
    AliInfo("Wrong VolumeID or Geometry not loaded - initializing void supermodule...");
  }
  AssignDetType();
} 


//_____________________________________________________________________________
AliITSAlignMille2Module::AliITSAlignMille2Module(const AliITSAlignMille2Module &m) :
  TNamed(m),
  fNSensVol(m.fNSensVol),
  fIndex(m.fIndex),  
  fDetType(m.fDetType),
  fVolumeID(m.fVolumeID),
  fNParTot(m.fNParTot),
  fNParFree(m.fNParFree),
  fParOffs(m.fNParTot),
  fNProcPoints(0),
  fParVals(0),
  fParErrs(0),
  fParCstr(0),  
  fSensVolIndex(m.fSensVolIndex),
  fSensVolVolumeID(m.fSensVolVolumeID),
  fMatrix(new TGeoHMatrix(*m.GetMatrix())),
  fSensVolMatrix(new TGeoHMatrix),
  fSensVolModifMatrix(new TGeoHMatrix),
  fParent(m.fParent),
  fChildren(0)
{
  // Copy constructor
  fSensVolIndex = m.fSensVolIndex;
  fSensVolVolumeID = m.fSensVolVolumeID;
  for (int i=m.fNParTot;i--;) fParOffs[i] = m.fParOffs[i];
  for (int i=3;i--;) fSigmaFactor[i] = m.fSigmaFactor[i];
  if (fNParTot) {
    fParVals = new Float_t[fNParTot];
    fParErrs = new Float_t[fNParTot];
    fParCstr = new Float_t[fNParTot];
    for (int i=fNParTot;i--;) {
      fParVals[i] = m.fParVals[i];
      fParErrs[i] = m.fParErrs[i];
      fParCstr[i] = m.fParCstr[i];
    }
  }
}

//_____________________________________________________________________________
AliITSAlignMille2Module& AliITSAlignMille2Module::operator=(const AliITSAlignMille2Module &m)  
{
  // operator =
  //
  if(this==&m) return *this;
  ((TNamed *)this)->operator=(m);
  //
  fNSensVol=m.fNSensVol;
  fIndex=m.fIndex;
  fDetType = m.fDetType;
  fVolumeID=m.fVolumeID;
  fNParTot  = m.fNParTot;
  fNParFree = m.fNParFree; 
  fNProcPoints = m.fNProcPoints; 
  delete[] fParVals; fParVals = 0;
  delete[] fParErrs; fParErrs = 0;
  delete[] fParCstr; fParCstr = 0;
  //
  if (fNParTot) {
    fParVals = new Float_t[fNParTot];
    fParErrs = new Float_t[fNParTot];
    fParCstr = new Float_t[fNParTot];
    for (int i=m.GetNParTot();i--;) {
      fParVals[i] = m.fParVals[i];
      fParErrs[i] = m.fParErrs[i];
      fParCstr[i] = m.fParCstr[i];
    }
  }
  //
  fParOffs.Set(fNParTot);
  for (int i=0;i<fNParTot;i++) fParOffs[i] = m.fParOffs[i];
  for (int i=0;i<3;i++) fSigmaFactor[i] = m.fSigmaFactor[i];
  if (fMatrix) delete fMatrix;
  fMatrix=new TGeoHMatrix(*m.GetMatrix());
  fSensVolIndex = m.fSensVolIndex;
  fSensVolVolumeID = m.fSensVolVolumeID;
  fParent = m.fParent;
  fChildren.Clear();
  for (int i=0;i<m.GetNChildren();i++) fChildren.Add(m.GetChild(i));
  return *this;
}


//-------------------------------------------------------------
AliITSAlignMille2Module::~AliITSAlignMille2Module() { 
  /// Destructor 
  delete fMatrix; 
  delete fSensVolMatrix; 
  delete fSensVolModifMatrix; 
  delete[] fParVals;
  delete[] fParErrs;
  delete[] fParCstr;
  fChildren.Clear();
} 

//-------------------------------------------------------------
Int_t AliITSAlignMille2Module::Set(Int_t index, UShort_t volid, const char* symname, 
				   const TGeoHMatrix *m, Int_t nsv, const UShort_t *volidsv) 
{
  // initialize a custom supermodule
  // index, volid, symname and matrix must be given
  // if (volidsv) add nsv sensitive volumes to the supermodules
  // return 0 if success

  if (index<2198) {
    AliInfo("Index must be >= 2198");
    return -1;
  }
  if (volid<14336) {
    AliInfo("VolumeID must be >= 14336");
    return -2;
  }
  
  if (!symname) return -3;
  for (Int_t i=0; i<2198; i++) {
    if (!strcmp(symname,AliITSgeomTGeo::GetSymName(i))) {
      AliInfo("Symname already used by a Sensitive Volume");
      return -3;
    }
  }
  
  if (!m) return -4;

  // can initialize needed stuffs
  fIndex = index;
  fVolumeID = volid;
  SetName(symname);
  //
  (*fMatrix) = (*m);
  //
  fSensVolIndex.Set(nsv);
  fSensVolVolumeID.Set(nsv);
  // add sensitive volumes
  for (Int_t i=0; i<nsv; i++) AddSensitiveVolume(volidsv[i]);

  return 0;
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::SetFreeDOF(Int_t dof,Double_t cstr)
{
  if (AliITSAlignMille2::IsZero(cstr)) fParCstr[dof] = 0;  // fixed parameter
  else if (cstr>0)                     fParCstr[dof] = fgkDummyConstraint+1.; // the parameter is free and unconstrained
  else                                 fParCstr[dof] = -cstr;                 // the parameter is free but constrained
}

//-------------------------------------------------------------
Bool_t AliITSAlignMille2Module::IsSensor(UShort_t voluid) 
{
  // Does this volid correspond to sensor ?
  AliGeomManager::ELayerID layId = AliGeomManager::VolUIDToLayerSafe(voluid);
  if (layId>0 && layId<7) {
    Int_t mId = Int_t(voluid & 0x7ff);
    if( mId>=0 && mId<AliGeomManager::LayerSize(layId)) return kTRUE;
  }
  return kFALSE;
}

//-------------------------------------------------------------
Int_t AliITSAlignMille2Module::GetIndexFromVolumeID(UShort_t voluid) {
  /// index from volume ID
  AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(voluid);
  if (lay<1|| lay>6) return -1;
  Int_t idx=Int_t(voluid)-2048*lay;
  if (idx>=AliGeomManager::LayerSize(lay)) return -1;
  for (Int_t ilay=1; ilay<lay; ilay++) 
    idx += AliGeomManager::LayerSize(ilay);
  return idx;
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::AddSensitiveVolume(UShort_t voluid)
{
  /// add a sensitive volume to this supermodule
  if (GetIndexFromVolumeID(voluid)<0) return; // bad volid
  //
  // in principle, the correct size of fSensVol... arrays was set outside but check anyway
  if (fSensVolVolumeID.GetSize()<fNSensVol+1) {
    fSensVolVolumeID.Set(fNSensVol+1);
    fSensVolIndex.Set(fNSensVol+1);
  }
  //
  fSensVolVolumeID[fNSensVol] = Short_t(voluid);
  fSensVolIndex[fNSensVol] = GetIndexFromVolumeID(voluid);
  fNSensVol++;
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::DelSensitiveVolume(Int_t at)
{
  // Suppress sensor at position "at"
  // in fact we are swapping with the last valid one 
  int lastValid = --fNSensVol;
  int tmpv = fSensVolIndex[at];
  fSensVolIndex[at] = fSensVolIndex[lastValid];
  tmpv = fSensVolVolumeID[at];
  fSensVolVolumeID[at] = fSensVolVolumeID[lastValid];
  fSensVolVolumeID[lastValid] = tmpv;
  //
}

//-------------------------------------------------------------
Bool_t AliITSAlignMille2Module::IsIn(UShort_t voluid) const 
{
  /// check if voluid is defined
  if (!voluid) return kFALSE; // only positive voluid are accepted
  for (Int_t i=0; i<fNSensVol; i++) if (UShort_t(fSensVolVolumeID[i])==voluid) return kTRUE;
  return kFALSE;
}

//-------------------------------------------------------------
Bool_t AliITSAlignMille2Module::BelongsTo(AliITSAlignMille2Module* parent) const
{
  /// check if parent contains the sensors of this volume
  if (fNSensVol<1 || fNSensVol>=parent->GetNSensitiveVolumes()) return kFALSE;
  return parent->IsIn( fSensVolVolumeID[0] );
}

//-------------------------------------------------------------
TGeoHMatrix *AliITSAlignMille2Module::GetSensitiveVolumeModifiedMatrix(UShort_t voluid, const Double_t *delta,Bool_t local)
{
  // modify the original TGeoHMatrix of the sensitive module 'voluid' according
  // with a delta transform. applied to the supermodule matrix
  // return NULL if error

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager)  return NULL;

  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  tr[0]=delta[0]; // in centimeter
  tr[1]=delta[1]; 
  tr[2]=delta[2];
  ang[0]=delta[3]; // psi   (X)  in deg
  ang[1]=delta[4]; // theta (Y)
  ang[2]=delta[5]; // phi   (Z)
  //
  fgTempAlignObj.SetRotation(ang[0],ang[1],ang[2]);
  fgTempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  TGeoHMatrix hm;
  fgTempAlignObj.GetMatrix(hm);
  //printf("\n0: delta matrix\n");hm.Print();

  // 1) start setting fSensVolModif = fSensVol
  if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
  //
  if (local) {
    // 2) set fSensVolModif = SensVolRel
    fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() );
    // 3) multiply left by delta
    fSensVolModifMatrix->MultiplyLeft( &hm );
    // 4) multiply left by fMatrix
    fSensVolModifMatrix->MultiplyLeft( fMatrix );
  }
  else fSensVolModifMatrix->MultiplyLeft( &hm );
  //
  return fSensVolModifMatrix;
}

//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeMisalignment(UShort_t voluid, const Double_t *deltalocal)
{
  // calculate misalignment of sens.vol. 'voluid' according with a displacement 'deltalocal'
  // of the mother volume. The misalignment is returned as AliAlignObjParams object

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager) return NULL;
  
  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  tr[0]=deltalocal[0]; // in centimeter
  tr[1]=deltalocal[1]; 
  tr[2]=deltalocal[2];
  ang[0]=deltalocal[3]; // psi   (X)  in deg
  ang[1]=deltalocal[4]; // theta (Y)
  ang[2]=deltalocal[5]; // phi   (Z)
  //
  fgTempAlignObj.SetRotation(ang[0],ang[1],ang[2]);
  fgTempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  //
  return GetSensitiveVolumeMisalignment(voluid,&fgTempAlignObj);
}

//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeMisalignment(UShort_t voluid, const AliAlignObjParams *a)
{
  // return the misalignment of the sens. vol. 'voluid' corresponding with 
  // a misalignment 'a' in the mother volume
  // return NULL if error

  // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
  // G'sv = Gg * Dg * Lsv,g === Gsv * Dsv
  // Gg * Dg * Gg-1 * Gsv = Gsv * Gsv-1 * Gg * Dg * Gg-1 * Gsv
  //
  // => Dsv = (Gsv-1 * Gg * Dg * Gg-1 * Gsv)
  //

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager) return NULL;

  //a->Print("");

  // prepare the Delta matrix Dg
  TGeoHMatrix dg;
  a->GetMatrix(dg);
  //dg.Print();

  // 1) start setting fSensVolModif = Gsv
  if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
  //printf("\n1: modif=orig del sensvol\n");fSensVolModifMatrix->Print();

  // 2) set fSensVolModif = Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() );
  //printf("\n2: modif=relative del sensvol\n");fSensVolModifMatrix->Print();
 
  // 3) set fSensVolModif = Dg * Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( &dg );
  //printf("\n3: modif= delta*relative\n");fSensVolModifMatrix->Print();
  
  // 4) set fSensVolModif = Gg * Dg * Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( fMatrix );
  //printf("\n4: modif=quasi,manca il Gsv-1...\n");fSensVolModifMatrix->Print();

  // 5) set fSensVolModif = Gsv-1 * Gg * Dg * Gg-1 * Gsv
  if (SensVolMatrix(voluid, &dg)) return NULL;
  fSensVolModifMatrix->MultiplyLeft( &dg.Inverse() );
  //printf("\n5: modif=finale\n");fSensVolModifMatrix->Print();
  //
  // >> RS
 // 6) mo' fSensVolModif dovrebbe essere la Dsv(loc) t.c. G'sv = Gsv*Dsv(loc)
  // per trasformarla in Dsv(loc rispetto al Gsv0, non modificato) dovrebbe essere:
  // Dsv(loc) -> Dpre * Dsv(loc) * Dpre-1
  //TGeoHMatrix dpre; // dpre = Gsv0-1*Gsv
  //if (SensVolOrigGlobalMatrix(voluid, &dg)) return NULL;
  //if (SensVolMatrix(voluid, &dpre)) return NULL;
  //dpre.MultiplyLeft( &dg.Inverse() );
  //fSensVolModifMatrix->Multiply( &dpre.Inverse() );
  //fSensVolModifMatrix->MultiplyLeft( &dpre );
  // direi che NON FUNZIONA!!!!  

  // << RS

  // reset align object (may not be needed...)
  fgTempAlignObj.SetVolUID(0);
  fgTempAlignObj.SetSymName("");
  fgTempAlignObj.SetTranslation(0,0,0);
  fgTempAlignObj.SetRotation(0,0,0);
  //
  // >> RS
#ifdef CORHW_
  // correction for SPD y-shift
  if (voluid>=2048 && voluid<4256) {
    TGeoHMatrix deltay;
    double dy[3]={0.,0.0081,0.};
    deltay.SetTranslation(dy);
    fSensVolModifMatrix->MultiplyLeft( &deltay );
    fSensVolModifMatrix->Multiply( &deltay.Inverse() );
  }
#endif
  // << RS
  if (!fgTempAlignObj.SetMatrix(*fSensVolModifMatrix)) return NULL;
  fgTempAlignObj.SetVolUID(voluid);
  fgTempAlignObj.SetSymName(AliGeomManager::SymName(voluid));
  //
  return &fgTempAlignObj;
}

// >> RS
//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeTotalMisalignment(UShort_t voluid, const Double_t *deltalocal)
{
  // calculate misalignment of sens.vol. 'voluid' according with a displacement 'deltalocal'
  // of the mother volume. The misalignment is returned as AliAlignObjParams object including
  // the (evenctual) prealignment => no merging needed

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager) return NULL;
  
  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  tr[0]=deltalocal[0]; // in centimeter
  tr[1]=deltalocal[1]; 
  tr[2]=deltalocal[2];
  ang[0]=deltalocal[3]; // psi   (X)  in deg
  ang[1]=deltalocal[4]; // theta (Y)
  ang[2]=deltalocal[5]; // phi   (Z)

  // reset align object (may not be needed...)
  fgTempAlignObj.SetVolUID(0);
  fgTempAlignObj.SetSymName("");
  fgTempAlignObj.SetRotation(ang[0],ang[1],ang[2]);
  fgTempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));

  // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
  // G'sv = Gg * Dg * Lsv,g === DGsv * Gsv 
  //
  // => Dsv = (G0sv-1 * Gg * Dg * Gg-1 * GMsv)  //
  //

  // prepare the Delta matrix Dg
  TGeoHMatrix dg;
  fgTempAlignObj.GetMatrix(dg);
  //dg.Print();

  // 1) start setting fSensVolModif = Gsv
  if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
  //printf("\n1: modif=orig del sensvol\n");fSensVolModifMatrix->Print();

  // 2) set fSensVolModif = Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() );
  //printf("\n2: modif=relative del sensvol\n");fSensVolModifMatrix->Print();
 
  // 3) set fSensVolModif = Dg * Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( &dg );
  //printf("\n3: modif= delta*relative\n");fSensVolModifMatrix->Print();
  
  // 4) set fSensVolModif = Gg * Dg * Gg-1 * Gsv
  fSensVolModifMatrix->MultiplyLeft( fMatrix );
  //printf("\n4: modif=quasi,manca il Gsv-1...\n");fSensVolModifMatrix->Print();

  // 5) set fSensVolModif = G0sv-1 * Gg * Dg * Gg-1 * Gsv
  // qui usa l'orig anziche' la prealigned...
  if (SensVolOrigGlobalMatrix(voluid, &dg)) return NULL;
  fSensVolModifMatrix->MultiplyLeft( &dg.Inverse() );
  //printf("\n5: modif=finale\n");fSensVolModifMatrix->Print();

  // reset align object (may not be needed...)
  fgTempAlignObj.SetVolUID(0);
  fgTempAlignObj.SetSymName("");
  fgTempAlignObj.SetTranslation(0,0,0);
  fgTempAlignObj.SetRotation(0,0,0);

#ifdef CORHW_
  // correction for SPD y-shift
  if (voluid>=2048 && voluid<4256) {
    TGeoHMatrix deltay;
    double dy[3]={0.,0.0081,0.};
    deltay.SetTranslation(dy);
    fSensVolModifMatrix->MultiplyLeft( &deltay );
    fSensVolModifMatrix->Multiply( &deltay.Inverse() );
  }
#endif
  if (!fgTempAlignObj.SetMatrix(*fSensVolModifMatrix)) return NULL;
  fgTempAlignObj.SetVolUID(voluid);
  fgTempAlignObj.SetSymName(AliGeomManager::SymName(voluid));

  
  //fgTempAlignObj.Print("");

  return &fgTempAlignObj;
}
//-------------------------------------------------------------

//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, const Double_t *deltalocal)
{
  // calculate misalignment of sens.vol. 'voluid' according with a displacement 'deltalocal'
  // of the mother volume. The misalignment is returned as AliAlignObjParams object

  if (!IsIn(voluid)) return NULL;
  if (!gGeoManager) return NULL;
  
  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  tr[0]=deltalocal[0]; // in centimeter
  tr[1]=deltalocal[1]; 
  tr[2]=deltalocal[2];
  ang[0]=deltalocal[3]; // psi   (X)  in deg
  ang[1]=deltalocal[4]; // theta (Y)
  ang[2]=deltalocal[5]; // phi   (Z)

  // reset align object (may not be needed...)
  fgTempAlignObj.SetTranslation(0,0,0);
  fgTempAlignObj.SetRotation(0,0,0);

  fgTempAlignObj.SetRotation(ang[0],ang[1],ang[2]);
  fgTempAlignObj.SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));

  // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
  // G'sv = Gg * Dg * Lsv,g === DGsv * Gsv 
  //
  // => DGsv = (Gg * Dg * Gg-1)
  //

  // prepare the Delta matrix Dg
  TGeoHMatrix dg;
  fgTempAlignObj.GetMatrix(dg);
  //dg.Print();

  dg.MultiplyLeft( fMatrix );
  dg.Multiply( &fMatrix->Inverse() );

  // reset align object (may not be needed...)
  fgTempAlignObj.SetTranslation(0,0,0);
  fgTempAlignObj.SetRotation(0,0,0);

  fgTempAlignObj.SetVolUID(voluid);
  fgTempAlignObj.SetSymName(AliGeomManager::SymName(voluid));

  if (!fgTempAlignObj.SetMatrix(dg)) return NULL;
  
  //fgTempAlignObj.Print("");

  return &fgTempAlignObj;
}
// << RS

//-------------------------------------------------------------
TGeoHMatrix *AliITSAlignMille2Module::GetSensitiveVolumeMatrix(UShort_t voluid)
{
  // return TGeoHMatrix of the sens.vol. 'voluid' in the current geometry
  if (SensVolMatrix(voluid,fSensVolMatrix)) return NULL;
  return fSensVolMatrix;
}

//-------------------------------------------------------------
TGeoHMatrix *AliITSAlignMille2Module::GetSensitiveVolumeOrigGlobalMatrix(UShort_t voluid)
{
  // return original ideal position (from AliGeomManager::GetOrigGlobalMatrix())
  if (SensVolOrigGlobalMatrix(voluid,fSensVolMatrix)) return NULL;
  return fSensVolMatrix;
}
//-------------------------------------------------------------
Int_t AliITSAlignMille2Module::SensVolMatrix(UShort_t volid, TGeoHMatrix *m) 
{
  // set matrix for sensitive modules (SPD corrected)
  // return 0 if success
  Double_t rot[9];
  Int_t idx=GetIndexFromVolumeID(volid);
  if (idx<0) return -1;
  if (!AliITSgeomTGeo::GetRotation(idx,rot)) return -2;
  m->SetRotation(rot);
  Double_t oLoc[3]={0,0,0};
  Double_t oGlo[3]={0,0,0};
  if (!AliITSgeomTGeo::LocalToGlobal(idx,oLoc,oGlo)) return -3;
  m->SetTranslation(oGlo);
  return 0;
}

//-------------------------------------------------------------
Int_t AliITSAlignMille2Module::SensVolOrigGlobalMatrix(UShort_t volid, TGeoHMatrix *m) 
{
  // set original global matrix for sensitive modules (SPD corrected)
  // return 0 if success
  Int_t idx=GetIndexFromVolumeID(volid);
  if (idx<0) return -1;
  TGeoHMatrix mo;
  if (!AliGeomManager::GetOrigGlobalMatrix(AliGeomManager::SymName(volid),mo)) return -1;
  (*m)=mo;
  //
#ifdef CORHW_
  // SPD y-shift by 81 mu
  if (volid<5000) { 
    Double_t oLoc[3]={0.0,0.0081,0.0};
    Double_t oGlo[3]={0,0,0};
    m->LocalToMaster(oLoc,oGlo);
    m->SetTranslation(oGlo);
  }
#endif
  return 0;
}

//-------------------------------------------------------------
UShort_t AliITSAlignMille2Module::GetVolumeIDFromSymname(const Char_t *symname) {
  /// volume ID from symname
  if (!symname) return 0;

  for (UShort_t voluid=2000; voluid<13300; voluid++) {
    Int_t modId;
    AliGeomManager::ELayerID layerId = AliGeomManager::VolUIDToLayer(voluid,modId);
    if (layerId>0 && layerId<7 && modId>=0 && modId<AliGeomManager::LayerSize(layerId)) {
      if (!strcmp(symname,AliGeomManager::SymName(layerId,modId))) return voluid;
    }
  }

  return 0;
}

//-------------------------------------------------------------
UShort_t AliITSAlignMille2Module::GetVolumeIDFromIndex(Int_t index) {
  /// volume ID from index
  if (index<0 || index>2197) return 0;
  return GetVolumeIDFromSymname(AliITSgeomTGeo::GetSymName(index));
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::Print(Option_t*) const 
{
  // print data
  //
  const char* typeName[] = {"SPD","SDD","SSD"};
  printf("*** ITS SuperModule for AliITSAlignMille ***\n");
  printf("symname  : %s (type: %s)\n",GetName(),fDetType<0 ? "N/A":typeName[fDetType]);
  printf("parent   : %s | %d children\n",fParent ? fParent->GetName() : "N/A",GetNChildren());
  printf("volumeID : %4d  | index : %4d | Geom.Params are %s\n",fVolumeID,fIndex,
	 GeomParamsGlobal() ? "Global":"Local");
  printf("Factors  : X=%.2f Y=%.2f Z=%.2f\n"
	 "DOF: %cTx:%5d| %cTy:%5d| %cTz:%5d| %cPsi:%5d| %cTheta:%5d| %cPhi:%5d|",
	 fSigmaFactor[0],fSigmaFactor[1],fSigmaFactor[2],
	 IsFreeDOF(kDOFTX) ? '+':'-',GetParOffset(kDOFTX),IsFreeDOF(kDOFTY) ? '+':'-',GetParOffset(kDOFTY),
	 IsFreeDOF(kDOFTZ) ? '+':'-',GetParOffset(kDOFTZ),IsFreeDOF(kDOFPS) ? '+':'-',GetParOffset(kDOFPS),
	 IsFreeDOF(kDOFTH) ? '+':'-',GetParOffset(kDOFTH),IsFreeDOF(kDOFPH) ? '+':'-',GetParOffset(kDOFPH));
  if (IsSDD()) {
    printf("%cT0:%5d| %cDVl:%5d| %cDVr:%5d|",IsFreeDOF(kDOFT0)?'+':'-',GetParOffset(kDOFT0),
	   IsFreeDOF(kDOFDVL)?'+':'-',GetParOffset(kDOFDVL),IsFreeDOF(kDOFDVR)?'+':'-',GetParOffset(kDOFDVR));
    if (IsVDriftLRSame()) printf("(dVL=dVR)");
  }
  printf("\n");
  fMatrix->Print();
  printf("%4d Sensitive volumes | %6d Processed Points\n",fNSensVol,fNProcPoints);
  for (Int_t i=0; i<fNSensVol; i++) printf("   voluid[%d] = %d\n",i,UShort_t(fSensVolVolumeID[i]));
}

//-------------------------------------------------------------
Bool_t AliITSAlignMille2Module::IsAlignable() const
{
  // it it alignable?
  TGeoManager* geoManager = AliGeomManager::GetGeometry();
  if (!geoManager) {
    AliInfo("Couldn't initialize geometry");
    return kFALSE;
  }
  return geoManager->GetAlignableEntry(GetName())!=0;
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetLocalMatrix(TGeoHMatrix &mat) const
{
  // return the local matrix for transformation to its parent
  mat = *fMatrix;
  if (fParent) mat.MultiplyLeft( &fParent->GetMatrix()->Inverse() );
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::AssignDetType()
{
  // assign the detector type
  TString tp = GetName();
  if      (tp.Contains("SPD",TString::kIgnoreCase)) fDetType = kSPD;
  else if (tp.Contains("SDD",TString::kIgnoreCase)) fDetType = kSDD;
  else if (tp.Contains("SSD",TString::kIgnoreCase)) fDetType = kSSD;
  else fDetType = -1;
  fNParTot = IsSDD() ? kMaxParTot:kMaxParGeom;
  fNParFree = 0;
  fParVals = new Float_t[fNParTot];
  fParErrs = new Float_t[fNParTot];  
  fParCstr = new Float_t[fNParTot];  
  if (fParOffs.GetSize()<fNParTot) fParOffs.Set(fNParTot);
  for (int i=fNParTot;i--;) {
    fParVals[i] = fParErrs[i] = 0.; 
    fParCstr[i] = 0.;
    fParOffs[i] = -1;
  }
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::EvaluateDOF()
{
  // count d.o.f.
  fNParFree = 0;
  for (int i=fNParTot;i--;) if (IsFreeDOF(i)) fNParFree++;
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetSensVolGlobalParams(UShort_t volid,Double_t *t, Double_t *r)
{
  // return global parameters of the sensor volid
  for (int i=3;i--;) t[i] = r[i] = 0.;
  if (SensVolMatrix(volid,fSensVolMatrix)) return;  
  fgTempAlignObj.SetMatrix(*fSensVolMatrix);
  fgTempAlignObj.GetPars(t,r);
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetSensVolLocalParams(UShort_t volid,Double_t *t, Double_t *r)
{
  // return parameters of the sensor volid in the current module
  for (int i=3;i--;) t[i] = r[i] = 0.;
  if (SensVolMatrix(volid,fSensVolMatrix)) return;  
  fSensVolMatrix->MultiplyLeft( &fMatrix->Inverse() );
  fgTempAlignObj.SetMatrix(*fSensVolMatrix);
  fgTempAlignObj.GetPars(t,r);
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetSensVolGlobalParams(UShort_t volid,const Double_t* loct, const Double_t* locr,Double_t *t, Double_t *r)
{
  // return global parameters of the sensor volid modified by the localDelta params
  for (int i=3;i--;) t[i] = r[i] = 0.;
  if (SensVolMatrix(volid,fSensVolMatrix)) return;  
  fgTempAlignObj.SetTranslation(loct[0],loct[1],loct[2]);
  fgTempAlignObj.SetRotation(locr[0],locr[1],locr[2]);
  //
  fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);      // obtain local delta
  fSensVolModifMatrix->MultiplyLeft( fSensVolMatrix ); // obtain global delta
  fgTempAlignObj.SetMatrix(*fSensVolModifMatrix);
  fgTempAlignObj.GetPars(t,r);                         // obtain global params
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetSensVolLocalParams(UShort_t volid,const Double_t* loct,const Double_t* locr,Double_t *t, Double_t *r)
{
  // return parameters of the sensor volid (modified by the localDelta params) in the current volume
  for (int i=3;i--;) t[i] = r[i] = 0.;
  if (SensVolMatrix(volid,fSensVolMatrix)) return;  
  fgTempAlignObj.SetTranslation(loct[0],loct[1],loct[2]);
  fgTempAlignObj.SetRotation(locr[0],locr[1],locr[2]);
  //
  fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);      // obtain local delta
  fSensVolModifMatrix->MultiplyLeft( fSensVolMatrix ); // obtain global delta
  fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() ); // obtain delta in current volume
  fgTempAlignObj.SetMatrix(*fSensVolModifMatrix);
  fgTempAlignObj.GetPars(t,r);                         // obtain params
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::SetParVals(Double_t *vl,Int_t npar)
{
  // set parameters
  for (int i=TMath::Min(npar,(Int_t)fNParTot);i--;) fParVals[i] = vl[i];
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetGeomParamsGlo(Double_t *pars)
{
  // recompute parameters from local to global frame
  //
  // is there anything to do?
  if (GeomParamsGlobal()) {for (int i=kMaxParGeom;i--;) pars[i] = fParVals[i]; return;}
  //
  // IMPORTANT: It is assumed that the parents params are defined in a same way (local or global)
  // as for the current module. Since in the mp2 the modules are stored from parents to children,
  // it is safe to call this method in loop staring from the lowest level child, i.e. from the end
  // of the modules array.
  //
  // DeltaGlobal = (ModifParents)*DeltaLocal*(ModifParents)^-1 
  //
  *fSensVolMatrix = *fMatrix;   // current global matrix
  AliITSAlignMille2Module* parent = GetParent();
  while (parent) {
    if (parent->GeomParamsGlobal()) {
      AliError("Cannot convert params to Global when the parents are already Global\n");
      for (int i=kMaxParGeom;i--;) pars[i] = 0;
      return;
    }
    fSensVolMatrix->MultiplyLeft( &parent->GetMatrix()->Inverse() ); // Local Matrix
    Float_t *parpar = parent->GetParVals();
    fgTempAlignObj.SetTranslation(parpar[0],parpar[1],parpar[2]);
    fgTempAlignObj.SetRotation(parpar[3],parpar[4],parpar[5]);
    fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);
    fSensVolMatrix->MultiplyLeft(fSensVolModifMatrix);
    fSensVolMatrix->MultiplyLeft(parent->GetMatrix());  // global matrix after parents modifications
    parent = parent->GetParent();
  }
  //
  fgTempAlignObj.SetTranslation(fParVals[0],fParVals[1],fParVals[2]);
  fgTempAlignObj.SetRotation(fParVals[3],fParVals[4],fParVals[5]);
  fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);  // local delta matrix
  fSensVolModifMatrix->Multiply( &fSensVolMatrix->Inverse() );
  fSensVolModifMatrix->MultiplyLeft( fSensVolMatrix );
  fgTempAlignObj.SetMatrix( *fSensVolModifMatrix );  // global delta matrix
  fgTempAlignObj.GetPars(pars,pars+3);
  //
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetGeomParamsLoc(Double_t *pars)
{
  // recompute parameters from global to local frame
  //
  // is there anything to do?
  if (!GeomParamsGlobal()) {for (int i=kMaxParGeom;i--;) pars[i] = fParVals[i]; return;}
  //
  // IMPORTANT: It is assumed that the parents params are defined in a same way (local or global)
  // as for the current module. Since in the mp2 the modules are stored from parents to children,
  // it is safe to call this method in loop staring from the lowest level child, i.e. from the end
  // of the modules array.
  //
  //  DeltaLocal = (DeltaParents*GlobalMat)^-1*DeltaGlobal*(DeltaParents*GlobalMat)
  //
  AliITSAlignMille2Module* parent = GetParent();
  fgTempAlignObj.SetTranslation(0.,0.,0.);
  fgTempAlignObj.SetRotation(0.,0.,0.);
  fgTempAlignObj.GetMatrix(*fSensVolMatrix); // get no-shift matrix
  //
  while (parent) { // accumulate the product of parents global modifications
    if (!parent->GeomParamsGlobal()) {
      AliError("Cannot convert params to Local when the parents are already Local\n");
      for (int i=kMaxParGeom;i--;) pars[i] = 0;
      return;
    }
    Float_t *parpar = parent->GetParVals();
    fgTempAlignObj.SetTranslation(parpar[0],parpar[1],parpar[2]);
    fgTempAlignObj.SetRotation(parpar[3],parpar[4],parpar[5]);
    fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);
    fSensVolMatrix->Multiply(fSensVolModifMatrix); 
    parent = parent->GetParent();
  }
  // global matrix after parents modifications
  fSensVolMatrix->Multiply(fMatrix);
  //
  fgTempAlignObj.SetTranslation(fParVals[0],fParVals[1],fParVals[2]);
  fgTempAlignObj.SetRotation(fParVals[3],fParVals[4],fParVals[5]);
  fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);  // global delta matrix
  fSensVolModifMatrix->MultiplyLeft( &fSensVolMatrix->Inverse() );
  fSensVolModifMatrix->Multiply( fSensVolMatrix );
  fgTempAlignObj.SetMatrix( *fSensVolModifMatrix );  // local delta matrix
  fgTempAlignObj.GetPars(pars,pars+3);
  //
}


//-------------------------------------------------------------
void AliITSAlignMille2Module::CalcDerivDPosDPar(Int_t sensVol,const Double_t* pl, Double_t *deriv)
{
  // calculate jacobian of the global position vs Parameters (dPos/dParam) 
  // for the point in the sensor sensVol
  const double kDel = 0.01;
  double pos0[3],pos1[3],pos2[3],pos3[3];
  double delta[kMaxParGeom];
  //
  for (int ip=kMaxParGeom;ip--;) delta[ip] = 0;
  //
  for (int ip=kMaxParGeom;ip--;) {
    //
    delta[ip] -= kDel;
    GetSensitiveVolumeModifiedMatrix(sensVol,delta,!GeomParamsGlobal())->LocalToMaster(pl,pos0);    
    delta[ip] += kDel/2;
    GetSensitiveVolumeModifiedMatrix(sensVol,delta,!GeomParamsGlobal())->LocalToMaster(pl,pos1);    
    delta[ip] += kDel;
    GetSensitiveVolumeModifiedMatrix(sensVol,delta,!GeomParamsGlobal())->LocalToMaster(pl,pos2);    
    delta[ip] += kDel/2;
    GetSensitiveVolumeModifiedMatrix(sensVol,delta,!GeomParamsGlobal())->LocalToMaster(pl,pos3);    
    //
    delta[ip] = 0;
    double *curd = deriv + ip*3;
    for (int i=3;i--;) curd[i] = (8.*(pos2[i]-pos1[i]) - (pos3[i]-pos0[i]))/6./kDel;
  }
  //
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::CalcDerivGloLoc(Int_t idx,Double_t *deriv)
{
  // calculate derivative of global params vs local param idx:  deriv[j] = dParGlo[j]/dParLoc[idx]
  Double_t lpar[kMaxParGeom];
  for (int i=kMaxParGeom;i--;) lpar[i] = 0.;
  //  using f(x+h),f(x-h),f(x+h/2),f(x-h2)...
  Double_t par1[kMaxParGeom]; // f(x-h)
  Double_t par2[kMaxParGeom]; // f(x-h/2)
  Double_t par3[kMaxParGeom]; // f(x+h/2)
  Double_t par4[kMaxParGeom]; // f(x+h)
  //
  const Double_t dpar = 1e-3;
  //
  // first values
  lpar[idx] -= dpar;
  GetGlobalParams(lpar,lpar+3, par1,par1+3);
  //
  // second values
  lpar[idx] += dpar/2;
  GetGlobalParams(lpar,lpar+3, par2,par2+3);
  //
  // third values
  lpar[idx] += dpar;
  GetGlobalParams(lpar,lpar+3, par3,par3+3);
  //
  // fourth values
  lpar[idx] += dpar/2;
  GetGlobalParams(lpar,lpar+3, par4,par4+3);
  //
  Double_t h2 = 1./(2.*dpar);
  for (int i=kMaxParGeom;i--;) {
    Double_t d0 = par4[i]-par1[i];
    Double_t d2 = 2.*(par3[i]-par2[i]);
    deriv[i] = h2*(4*d2 - d0)/3.;
    if (TMath::Abs(deriv[i]) < 1.0e-9) deriv[i] = 0.0;
  }
  //
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::CalcDerivLocGlo(Double_t *deriv)
{
  // calculate derivative of local params vs global params:  deriv[i][j] = dParLoc[i]/dParGlo[j]
  Double_t gpar[kMaxParGeom];
  for (int i=kMaxParGeom;i--;) gpar[i] = 0.;
  //  using f(x+h),f(x-h),f(x+h/2),f(x-h2)...
  Double_t par1[kMaxParGeom]; // f(x-h)
  Double_t par2[kMaxParGeom]; // f(x-h/2)
  Double_t par3[kMaxParGeom]; // f(x+h/2)
  Double_t par4[kMaxParGeom]; // f(x+h)
  //
  const Double_t dpar = 1e-3;
  //
  for (int ig=kMaxParGeom;ig--;) {
    // first values
    gpar[ig] -= dpar;
    GetLocalParams(gpar,gpar+3, par1,par1+3);
    //
    // second values
    gpar[ig] += dpar/2;
    GetLocalParams(gpar,gpar+3, par2,par2+3);
    //
    // third values
    gpar[ig] += dpar;
    GetLocalParams(gpar,gpar+3, par3,par3+3);
    //
    // fourth values
    gpar[ig] += dpar/2;
    GetLocalParams(gpar,gpar+3, par4,par4+3);
    //
    Double_t h2 = 1./(2.*dpar);
    for (int i=kMaxParGeom;i--;) {
      Double_t d0 = par4[i]-par1[i];
      Double_t d2 = 2.*(par3[i]-par2[i]);
      int idig = i*kMaxParGeom + ig;
      deriv[idig] = h2*(4*d2 - d0)/3.;
      if (TMath::Abs(deriv[idig]) < 1.0e-9) deriv[idig] = 0.0;
    }
  }
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2Module::CalcDerivGloLoc(Int_t sensVol,Int_t paridx,Double_t* derivative)
{
  /// calculate numerically the derivatives of global params vs local param paridx for sensor sensVol: dPglob/dPloc_paridx
  //
  Double_t lpar[kMaxParGeom];
  for (int i=kMaxParGeom;i--;) lpar[i] = 0.;
  //  using f(x+h),f(x-h),f(x+h/2),f(x-h2)...
  Double_t par1[kMaxParGeom]; // f(x-h)
  Double_t par2[kMaxParGeom]; // f(x-h/2)
  Double_t par3[kMaxParGeom]; // f(x+h/2)
  Double_t par4[kMaxParGeom]; // f(x+h)
  //
  const Double_t dpar = 1e-3;
  //
  // first values
  lpar[paridx] -= dpar;
  GetSensVolGlobalParams(sensVol,lpar,lpar+3, par1,par1+3);
  //
  // second values
  lpar[paridx] += dpar/2;
  GetSensVolGlobalParams(sensVol,lpar,lpar+3, par2,par2+3);
  //
  // third values
  lpar[paridx] += dpar;
  GetSensVolGlobalParams(sensVol,lpar,lpar+3, par3,par3+3);
  //
  // fourth values
  lpar[paridx] += dpar/2;
  GetSensVolGlobalParams(sensVol,lpar,lpar+3, par4,par4+3);
  //
  Double_t h2 = 1./(2.*dpar);
  for (int i=kMaxParGeom;i--;) {
    Double_t d0 = par4[i]-par1[i];
    Double_t d2 = 2.*(par3[i]-par2[i]);
    derivative[i] = h2*(4*d2 - d0)/3.;
    if (TMath::Abs(derivative[i]) < 1.0e-9) derivative[i] = 0.0;
  }
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2Module::CalcDerivCurLoc(Int_t sensVol,Int_t paridx,Double_t* derivative)  
{
  /// calculate numerically the derivatives of sensor params in the current volume vs sensor local param paridx
  //
  Double_t lpar[kMaxParGeom];
  for (int i=kMaxParGeom;i--;) lpar[i] = 0.;
  //  using f(x+h),f(x-h),f(x+h/2),f(x-h2)...
  Double_t par1[kMaxParGeom]; // f(x-h)
  Double_t par2[kMaxParGeom]; // f(x-h/2)
  Double_t par3[kMaxParGeom]; // f(x+h/2)
  Double_t par4[kMaxParGeom]; // f(x+h)
  //
  const Double_t dpar = 1e-3;
  //
  // first values
  lpar[paridx] -= dpar;
  GetSensVolLocalParams(sensVol,lpar,lpar+3, par1,par1+3);
  //
  // second values
  lpar[paridx] += dpar/2;
  GetSensVolLocalParams(sensVol,lpar,lpar+3, par2,par2+3);
  //
  // third values
  lpar[paridx] += dpar;
  GetSensVolLocalParams(sensVol,lpar,lpar+3, par3,par3+3);
  //
  // fourth values
  lpar[paridx] += dpar/2;
  GetSensVolLocalParams(sensVol,lpar,lpar+3, par4,par4+3);
  //
  Double_t h2 = 1./(2.*dpar);
  for (int i=kMaxParGeom;i--;) {
    Double_t d0 = par4[i]-par1[i];
    Double_t d2 = 2.*(par3[i]-par2[i]);
    derivative[i] = h2*(4*d2 - d0)/3.;
    if (TMath::Abs(derivative[i]) < 1.0e-9) derivative[i] = 0.0;
  }
  //
}


//-------------------------------------------------------------
void AliITSAlignMille2Module::GetGlobalParams(Double_t *t, Double_t *r)
{
  // global parameters of the module
  fgTempAlignObj.SetMatrix( *fMatrix );
  fgTempAlignObj.GetPars(t,r);
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetGlobalParams(const Double_t* loct, const Double_t* locr, Double_t *t, Double_t *r)
{
  // global parameters of the module after the modification by local loct,locr
  fgTempAlignObj.SetTranslation(loct[0],loct[1],loct[2]);
  fgTempAlignObj.SetRotation(locr[0],locr[1],locr[2]);
  fgTempAlignObj.GetMatrix(*fSensVolModifMatrix);  
  *fSensVolMatrix = *fMatrix;
  fSensVolMatrix->Multiply(fSensVolModifMatrix);
  fgTempAlignObj.SetMatrix(*fSensVolMatrix);
  fgTempAlignObj.GetPars(t,r);
}

//-------------------------------------------------------------
void AliITSAlignMille2Module::GetLocalParams(const Double_t* glot, const Double_t* glor, Double_t *t, Double_t *r)
{
  // obtain local delta parameters from global delta params
  fgTempAlignObj.SetTranslation(glot[0],glot[1],glot[2]);
  fgTempAlignObj.SetRotation(glor[0],glor[1],glor[2]);
  fgTempAlignObj.GetMatrix(*fSensVolMatrix);  
  fSensVolMatrix->Multiply( fMatrix );
  fSensVolMatrix->MultiplyLeft( &fMatrix->Inverse() );
  fgTempAlignObj.SetMatrix(*fSensVolMatrix);
  fgTempAlignObj.GetPars(t,r);
}
