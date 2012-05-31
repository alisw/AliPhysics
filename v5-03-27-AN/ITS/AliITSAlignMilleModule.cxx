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
/// \class AliITSAlignMilleModule
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
//#include <TGeoMatrix.h> 
 
#include "AliITSAlignMilleModule.h" 
#include "AliITSgeomTGeo.h" 
#include "AliGeomManager.h" 
#include "AliAlignObjParams.h" 
#include "AliLog.h" 
 
/// \cond CLASSIMP 
ClassImp(AliITSAlignMilleModule) 
/// \endcond 
    
//-------------------------------------------------------------
AliITSAlignMilleModule::AliITSAlignMilleModule() : TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),  
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fTempAlignObj(NULL)
{ 
  /// void constructor  
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix; 
  fTempAlignObj=new AliAlignObjParams;
  for(Int_t k=0; k<fgkSensModules; k++){
    fSensVolIndex[k] = 0;
    fSensVolVolumeID[k] = 0;
  }
} 
//-------------------------------------------------------------
AliITSAlignMilleModule::AliITSAlignMilleModule(Int_t index, UShort_t volid, char* symname, const TGeoHMatrix *m, Int_t nsv, const UShort_t *volidsv) : TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),  
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fTempAlignObj(NULL)
{ 
  /// void constructor  
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix; 
  fTempAlignObj=new AliAlignObjParams;
  if (Set(index,volid,symname,m,nsv,volidsv)) {
    AliInfo("Error in AliITSAlignMilleModule::Set() - initializing void supermodule...");

  }
  for(Int_t k=0; k<fgkSensModules; k++){
    fSensVolIndex[k] = 0;
    fSensVolVolumeID[k] = 0;
  }
}
//-------------------------------------------------------------
AliITSAlignMilleModule::AliITSAlignMilleModule(UShort_t volid) : TNamed(), 
  fNSensVol(0), 
  fIndex(-1),  
  fVolumeID(0),  
  fMatrix(NULL),
  fSensVolMatrix(NULL),
  fSensVolModifMatrix(NULL),
  fTempAlignObj(NULL)
{ 
  /// simple constructor building a supermodule from a single sensitive volume 
  fMatrix = new TGeoHMatrix; 
  fSensVolMatrix = new TGeoHMatrix; 
  fSensVolModifMatrix = new TGeoHMatrix;   
  // temporary align object, just use the rotation...
  fTempAlignObj=new AliAlignObjParams;

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
  for(Int_t k=0; k<fgkSensModules; k++){
    fSensVolIndex[k] = 0;
    fSensVolVolumeID[k] = 0;
  }
  }
} 
//-------------------------------------------------------------
AliITSAlignMilleModule::~AliITSAlignMilleModule() { 
  /// Destructor 
  delete fMatrix; 
  delete fSensVolMatrix; 
  delete fSensVolModifMatrix; 
  delete fTempAlignObj;
} 
//-------------------------------------------------------------
Int_t AliITSAlignMilleModule::Set(Int_t index, UShort_t volid, char* symname, const TGeoHMatrix * const m, Int_t nsv, const UShort_t *volidsv) 
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
  (*fMatrix) = (*m);

  // add sensitive volumes
  for (Int_t i=0; i<nsv; i++) AddSensitiveVolume(volidsv[i]);

  return 0;
}
//-------------------------------------------------------------
Int_t AliITSAlignMilleModule::GetIndexFromVolumeID(UShort_t voluid) {
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
void AliITSAlignMilleModule::AddSensitiveVolume(UShort_t voluid)
{
  /// add a sensitive volume to this supermodule
  if (GetIndexFromVolumeID(voluid)<0) return; // bad volid
  fSensVolVolumeID[fNSensVol] = voluid;
  fSensVolIndex[fNSensVol] = GetIndexFromVolumeID(voluid);
  fNSensVol++;
}
//-------------------------------------------------------------
Bool_t AliITSAlignMilleModule::IsIn(UShort_t voluid) const 
{
  /// check if voluid is defined
  if (!voluid) return kFALSE; // only positive voluid are accepted
  for (Int_t i=0; i<fNSensVol; i++) {
    if (fSensVolVolumeID[i]==voluid) return kTRUE;
  }
  return kFALSE;
}
//-------------------------------------------------------------
TGeoHMatrix *AliITSAlignMilleModule::GetSensitiveVolumeModifiedMatrix(UShort_t voluid, const Double_t * const deltalocal)
{
  // modify the original TGeoHMatrix of the sensitive module 'voluid' according
  // with a delta transform. applied to the supermodule matrix
  // return NULL if error

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
  fTempAlignObj->SetVolUID(0);
  fTempAlignObj->SetSymName("");
  fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
  fTempAlignObj->SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  TGeoHMatrix hm;
  fTempAlignObj->GetMatrix(hm);
  //printf("\n0: delta matrix\n");hm.Print();

  // 1) start setting fSensVolModif = fSensVol
  if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
  //printf("\n1: modif=orig del sensvol\n");fSensVolModifMatrix->Print();

  // 2) set fSensVolModif = SensVolRel
  fSensVolModifMatrix->MultiplyLeft( &fMatrix->Inverse() );
  //printf("\n2: modif=relative del sensvol\n");fSensVolModifMatrix->Print();
 
  // 3) multiply left by delta
  fSensVolModifMatrix->MultiplyLeft( &hm );
  //printf("\n3: modif= delta*relative\n");fSensVolModifMatrix->Print();
  
  // 4) multiply left by fMatrix
  fSensVolModifMatrix->MultiplyLeft( fMatrix );
  //printf("\n4: modif=finale\n");fSensVolModifMatrix->Print();

  return fSensVolModifMatrix;
}
//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMilleModule::GetSensitiveVolumeMisalignment(UShort_t voluid, const Double_t * const deltalocal)
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
  fTempAlignObj->SetVolUID(0);
  fTempAlignObj->SetSymName("");
  fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
  fTempAlignObj->SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  
  return GetSensitiveVolumeMisalignment(voluid,fTempAlignObj);
}
//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMilleModule::GetSensitiveVolumeMisalignment(UShort_t voluid, const AliAlignObjParams *a)
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

  // reset align object (may not be needed...)
  fTempAlignObj->SetVolUID(0);
  fTempAlignObj->SetSymName("");
  fTempAlignObj->SetTranslation(0,0,0);
  fTempAlignObj->SetRotation(0,0,0);

  // correction for SPD y-shift
  if (voluid>=2048 && voluid<4256) {
    TGeoHMatrix deltay;
    double dy[3]={0.,0.0081,0.};
    deltay.SetTranslation(dy);
    fSensVolModifMatrix->MultiplyLeft( &deltay );
    fSensVolModifMatrix->Multiply( &deltay.Inverse() );
  }
  
  if (!fTempAlignObj->SetMatrix(*fSensVolModifMatrix)) return NULL;
  fTempAlignObj->SetVolUID(voluid);
  fTempAlignObj->SetSymName(AliGeomManager::SymName(voluid));
  
  //fTempAlignObj->Print("");

  return fTempAlignObj;
}
//-------------------------------------------------------------


//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMilleModule::GetSensitiveVolumeTotalMisalignment(UShort_t voluid, const Double_t * const deltalocal)
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
  fTempAlignObj->SetVolUID(0);
  fTempAlignObj->SetSymName("");
  fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
  fTempAlignObj->SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));

  // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
  // G'sv = Gg * Dg * Lsv,g === DGsv * Gsv 
  //
  // => Dsv = (G0sv-1 * Gg * Dg * Gg-1 * GMsv)  //
  //

  // prepare the Delta matrix Dg
  TGeoHMatrix dg;
  fTempAlignObj->GetMatrix(dg);
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
  fTempAlignObj->SetVolUID(0);
  fTempAlignObj->SetSymName("");
  fTempAlignObj->SetTranslation(0,0,0);
  fTempAlignObj->SetRotation(0,0,0);

  // correction for SPD y-shift
  if (voluid>=2048 && voluid<4256) {
    TGeoHMatrix deltay;
    double dy[3]={0.,0.0081,0.};
    deltay.SetTranslation(dy);
    fSensVolModifMatrix->MultiplyLeft( &deltay );
    fSensVolModifMatrix->Multiply( &deltay.Inverse() );
  }
  if (!fTempAlignObj->SetMatrix(*fSensVolModifMatrix)) return NULL;
  fTempAlignObj->SetVolUID(voluid);
  fTempAlignObj->SetSymName(AliGeomManager::SymName(voluid));

  
  //fTempAlignObj->Print("");

  return fTempAlignObj;
}
//-------------------------------------------------------------

// //-------------------------------------------------------------
// AliAlignObjParams *AliITSAlignMilleModule::GetSensitiveVolumeTotalMisalignment(UShort_t voluid, Double_t *deltalocal)
// {
//   // calculate misalignment of sens.vol. 'voluid' according with a displacement 'deltalocal'
//   // of the mother volume. The misalignment is returned as AliAlignObjParams object including
//   // the (evenctual) prealignment => no merging needed

//   if (!IsIn(voluid)) return NULL;
//   if (!gGeoManager) return NULL;
  
//   // prepare the TGeoHMatrix
//   Double_t tr[3],ang[3];
//   tr[0]=deltalocal[0]; // in centimeter
//   tr[1]=deltalocal[1]; 
//   tr[2]=deltalocal[2];
//   ang[0]=deltalocal[3]; // psi   (X)  in deg
//   ang[1]=deltalocal[4]; // theta (Y)
//   ang[2]=deltalocal[5]; // phi   (Z)

//   // reset align object (may not be needed...)
//   fTempAlignObj->SetVolUID(0);
//   fTempAlignObj->SetSymName("");
//   fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
//   fTempAlignObj->SetTranslation(tr[0],tr[1],tr[2]);
//   AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));

//   // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
//   // G'sv = Gg * Dg * Lsv,g === DGsv * Gsv 
//   //
//   // => Dsv = (G0sv-1 * Gg * Dg * Gg-1 * GMsv)  //
//   //

//   // prepare the Delta matrix Dg
//   TGeoHMatrix dg;
//   fTempAlignObj->GetMatrix(dg);
//   //dg.Print();

//   // 1) start setting fSensVolModif = Gsv
//   if (SensVolMatrix(voluid, fSensVolModifMatrix)) return NULL;
//   //printf("\n1: modif=orig del sensvol\n");fSensVolModifMatrix->Print();

//   // 2) set dg = Gg * Dg * Gg-1 
//   dg.Multiply( &fMatrix->Inverse() );
//   dg.MultiplyLeft( fMatrix );
 
//   // 3) set dg = Gsv-1 * dg * Gsv locale nel sistema del SV preallineato
//   dg.Multiply( fSensVolModifMatrix );
//   dg.MultiplyLeft( &fSensVolModifMatrix->Inverse() );

//   // calcola la deltaPre
//   TGeoHMatrix hp;
//   if (SensVolOrigGlobalMatrix(voluid, &hp)) return NULL;
//   fSensVolModifMatrix->MultiplyLeft( &hp.Inverse() );

//   // mo' fSensVolModifMatrix e' la deltapre ideale
//   // col dg diventa la deltatot ideale
//   fSensVolModifMatrix->Multiply( &dg );
//   //fSensVolModifMatrix->MultiplyLeft( &dg );
  
//   // reset align object (may not be needed...)
//   fTempAlignObj->SetVolUID(0);
//   fTempAlignObj->SetSymName("");
//   fTempAlignObj->SetTranslation(0,0,0);
//   fTempAlignObj->SetRotation(0,0,0);

//   // correction for SPD y-shift
//   if (voluid>=2048 && voluid<4256) {
//     TGeoHMatrix deltay;
//     double dy[3]={0.,0.0081,0.};
//     deltay.SetTranslation(dy);
//     fSensVolModifMatrix->MultiplyLeft( &deltay );
//     fSensVolModifMatrix->Multiply( &deltay.Inverse() );
//   }
//   if (!fTempAlignObj->SetMatrix(*fSensVolModifMatrix)) return NULL;
//   fTempAlignObj->SetVolUID(voluid);
//   fTempAlignObj->SetSymName(AliGeomManager::SymName(voluid));

  
//   //fTempAlignObj->Print("");

//   return fTempAlignObj;
// }
// //-------------------------------------------------------------


//-------------------------------------------------------------
AliAlignObjParams *AliITSAlignMilleModule::GetSensitiveVolumeGlobalMisalignment(UShort_t voluid, const Double_t * const deltalocal)
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
  fTempAlignObj->SetTranslation(0,0,0);
  fTempAlignObj->SetRotation(0,0,0);

  fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
  fTempAlignObj->SetTranslation(tr[0],tr[1],tr[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));

  // Gsv = Gg * Gg-1 * Gsv   -> Lsv,g = Gg-1 * Gsv
  // G'sv = Gg * Dg * Lsv,g === DGsv * Gsv 
  //
  // => DGsv = (Gg * Dg * Gg-1)
  //

  // prepare the Delta matrix Dg
  TGeoHMatrix dg;
  fTempAlignObj->GetMatrix(dg);
  //dg.Print();

  dg.MultiplyLeft( fMatrix );
  dg.Multiply( &fMatrix->Inverse() );

  // reset align object (may not be needed...)
  fTempAlignObj->SetTranslation(0,0,0);
  fTempAlignObj->SetRotation(0,0,0);

  fTempAlignObj->SetVolUID(voluid);
  fTempAlignObj->SetSymName(AliGeomManager::SymName(voluid));

  if (!fTempAlignObj->SetMatrix(dg)) return NULL;
  
  //fTempAlignObj->Print("");

  return fTempAlignObj;
}
//-------------------------------------------------------------

TGeoHMatrix *AliITSAlignMilleModule::GetSensitiveVolumeMatrix(UShort_t voluid)
{
  // return TGeoHMatrix of the sens.vol. 'voluid' in the current geometry
  if (SensVolMatrix(voluid,fSensVolMatrix)) return NULL;
  return fSensVolMatrix;
}
//-------------------------------------------------------------
TGeoHMatrix *AliITSAlignMilleModule::GetSensitiveVolumeOrigGlobalMatrix(UShort_t voluid)
{
  // return original ideal position (from AliGeomManager::GetOrigGlobalMatrix())
  if (SensVolOrigGlobalMatrix(voluid,fSensVolMatrix)) return NULL;
  return fSensVolMatrix;
}
//-------------------------------------------------------------
Int_t AliITSAlignMilleModule::SensVolMatrix(UShort_t volid, TGeoHMatrix *m) 
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
Int_t AliITSAlignMilleModule::SensVolOrigGlobalMatrix(UShort_t volid, TGeoHMatrix *m) 
{
  // set original global matrix for sensitive modules (SPD corrected)
  // return 0 if success
  Int_t idx=GetIndexFromVolumeID(volid);
  if (idx<0) return -1;
  TGeoHMatrix mo;
  if (!AliGeomManager::GetOrigGlobalMatrix(AliGeomManager::SymName(volid),mo)) return -1;

  (*m)=mo;

  // SPD y-shift by 81 mu
  if (volid<5000) { 
    Double_t oLoc[3]={0.0,0.0081,0.0};
    Double_t oGlo[3]={0,0,0};
    m->LocalToMaster(oLoc,oGlo);
    m->SetTranslation(oGlo);
  }
  return 0;
}
//-------------------------------------------------------------
UShort_t AliITSAlignMilleModule::GetVolumeIDFromSymname(const Char_t *symname) {
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

UShort_t AliITSAlignMilleModule::GetVolumeIDFromIndex(Int_t index) {
  /// volume ID from index
  if (index<0 || index>2197) return 0;
  return GetVolumeIDFromSymname(AliITSgeomTGeo::GetSymName(index));
}
//-------------------------------------------------------------
void AliITSAlignMilleModule::Print(Option_t*) const 
{
  ///
  printf("*** ITS SuperModule for AliITSAlignMille ***\n");
  printf("symname  : %s\n",GetName());
  printf("volumeID : %d\n",fVolumeID);
  printf("index    : %d\n",fIndex);
  fMatrix->Print();
  printf("number of sensitive modules : %d\n",fNSensVol);
  for (Int_t i=0; i<fNSensVol; i++) printf("   voluid[%d] = %d\n",i,fSensVolVolumeID[i]);
}
//_____________________________________________________________________________
AliITSAlignMilleModule::AliITSAlignMilleModule(const AliITSAlignMilleModule &m) :
  TNamed(m),
  fNSensVol(m.fNSensVol),
  fIndex(m.fIndex),
  fVolumeID(m.fVolumeID),
  fMatrix(new TGeoHMatrix(*m.GetMatrix())),
  fSensVolMatrix(new TGeoHMatrix),
  fSensVolModifMatrix(new TGeoHMatrix),
  fTempAlignObj(new AliAlignObjParams)
{
  // Copy constructor
  for(Int_t k=0; k<fgkSensModules; k++){
    fSensVolIndex[k] = 0;
    fSensVolVolumeID[k] = 0;
  }
  for (int i=0; i<fNSensVol; i++) {
    fSensVolIndex[i]=m.fSensVolIndex[i];
    fSensVolVolumeID[i]=m.fSensVolVolumeID[i];
  }
}
//_____________________________________________________________________________
AliITSAlignMilleModule& AliITSAlignMilleModule::operator=(const AliITSAlignMilleModule &m)  
{
  // operator =
  //
  if(this==&m) return *this;
  ((TNamed *)this)->operator=(m);
  
  fNSensVol=m.fNSensVol;
  fIndex=m.fIndex;
  fVolumeID=m.fVolumeID;
  delete fMatrix;
  fMatrix=new TGeoHMatrix(*m.GetMatrix());
  for (int i=0; i<fNSensVol; i++) {
    fSensVolIndex[i]=m.fSensVolIndex[i];
    fSensVolVolumeID[i]=m.fSensVolVolumeID[i];
  }
  return *this;
}

//_____________________________________________________________________________


