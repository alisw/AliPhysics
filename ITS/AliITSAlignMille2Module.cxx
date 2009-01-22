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
    
//-------------------------------------------------------------
AliITSAlignMille2Module::AliITSAlignMille2Module() : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),
  fNProcPoints(0),
  fSensVolIndex(0),
  fSensVolVolumeID(0),
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL)
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
AliITSAlignMille2Module::AliITSAlignMille2Module(Int_t index,UShort_t volid,char* symname,TGeoHMatrix *m,Int_t nsv,UShort_t *volidsv) : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),
  fNProcPoints(0),
  fSensVolIndex(0),
  fSensVolVolumeID(0),  
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL)
{ 
  /// void constructor  
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix; 
  fSigmaFactor[0]=fSigmaFactor[1]=fSigmaFactor[2]=1.0;
  if (Set(index,volid,symname,m,nsv,volidsv)) {
    AliInfo("Error in AliITSAlignMille2Module::Set() - initializing void supermodule...");
  }
} 

//-------------------------------------------------------------
AliITSAlignMille2Module::AliITSAlignMille2Module(UShort_t volid) : 
  TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),
  fNProcPoints(0),
  fSensVolIndex(0),
  fSensVolVolumeID(0),
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fParent(NULL)
{ 
  /// simple constructor building a supermodule from a single sensitive volume 
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix;   
  // temporary align object, just use the rotation...
  fSensVolIndex.Set(1);
  fSensVolVolumeID.Set(1);
  fSigmaFactor[0]=fSigmaFactor[1]=fSigmaFactor[2]=1.0;
  //
  fIndex = GetIndexFromVolumeID(volid);  
  if (fIndex>=0 && gGeoManager) { // good sensitive module and geometry loaded
    SetName(AliGeomManager::SymName(volid));
    fVolumeID = volid;
    AddSensitiveVolume(volid);
    if (SensVolMatrix(volid, fMatrix))
       AliInfo("Matrix not defined");
  }
  else {
    AliInfo("Wrong VolumeID or Geometry not loaded - initializing void supermodule...");
  }
} 


//_____________________________________________________________________________
AliITSAlignMille2Module::AliITSAlignMille2Module(const AliITSAlignMille2Module &m) :
  TNamed(m),
  fNSensVol(m.fNSensVol),
  fIndex(m.fIndex),
  fVolumeID(m.fVolumeID),
  fNProcPoints(0),
  fSensVolIndex(m.fSensVolIndex),
  fSensVolVolumeID(m.fSensVolVolumeID),
  fMatrix(new TGeoHMatrix(*m.GetMatrix())),
  fSensVolMatrix(new TGeoHMatrix),
  fSensVolModifMatrix(new TGeoHMatrix),
  fParent(m.fParent)
{
  // Copy constructor
  fSensVolIndex = m.fSensVolIndex;
  fSensVolVolumeID = m.fSensVolVolumeID;
  for (int i=3;i--;) fSigmaFactor[i] = m.fSigmaFactor[i];
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
  fVolumeID=m.fVolumeID;
  for (int i=3;i--;) fSigmaFactor[i] = m.fSigmaFactor[i];
  if (fMatrix) delete fMatrix;
  fMatrix=new TGeoHMatrix(*m.GetMatrix());
  fSensVolIndex = m.fSensVolIndex;
  fSensVolVolumeID = m.fSensVolVolumeID;
  fParent = m.fParent;
  return *this;
}


//-------------------------------------------------------------
AliITSAlignMille2Module::~AliITSAlignMille2Module() { 
  /// Destructor 
  delete fMatrix; 
  delete fSensVolMatrix; 
  delete fSensVolModifMatrix; 
} 

//-------------------------------------------------------------
Int_t AliITSAlignMille2Module::Set(Int_t index, UShort_t volid, char* symname, TGeoHMatrix *m, Int_t nsv, UShort_t *volidsv) 
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
  if (fSensVolVolumeID.GetSize()<fNSensVol) {
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
TGeoHMatrix *AliITSAlignMille2Module::GetSensitiveVolumeModifiedMatrix(UShort_t voluid, Double_t *delta,Bool_t local)
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
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeMisalignment(UShort_t voluid, Double_t *deltalocal)
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
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeMisalignment(UShort_t voluid, AliAlignObjParams *a)
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
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeTotalMisalignment(UShort_t voluid, Double_t *deltalocal)
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
AliAlignObjParams *AliITSAlignMille2Module::GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, Double_t *deltalocal)
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
  if (!AliGeomManager::GetOrigGlobalMatrix(AliGeomManager::SymName(volid),mo));
  (*m)=mo;

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
  ///
  printf("*** ITS SuperModule for AliITSAlignMille ***\n");
  printf("symname  : %s\n",GetName());
  printf("parent   : %s\n",fParent ? fParent->GetName() : "N/A");
  printf("volumeID : %4d  | index : %4d\n",fVolumeID,fIndex);
  printf("Factors  : X=%.2f Y=%.2f Z=%.2f | DOF: Tx:%d Ty:%d Tz:%d Phi:%d Theta:%d Psi:%d\n",
	 fSigmaFactor[0],fSigmaFactor[1],fSigmaFactor[2],
	 IsFreeDOF(AliITSAlignMille2::kDOFTX),IsFreeDOF(AliITSAlignMille2::kDOFTY),
	 IsFreeDOF(AliITSAlignMille2::kDOFTZ),IsFreeDOF(AliITSAlignMille2::kDOFPH),
	 IsFreeDOF(AliITSAlignMille2::kDOFTH),IsFreeDOF(AliITSAlignMille2::kDOFPS));
  fMatrix->Print();
  printf("%4d Sensitive volumes | %6d Processed Points\n",fNSensVol,fNProcPoints);
  for (Int_t i=0; i<fNSensVol; i++) printf("   voluid[%d] = %d\n",i,UShort_t(fSensVolVolumeID[i]));
}

//-------------------------------------------------------------
Bool_t AliITSAlignMille2Module::IsAlignable() const
{
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
