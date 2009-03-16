#include "AliZDCMisAligner.h"
#include "AliGeomManager.h"
#include "TClonesArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"

ClassImp(AliZDCMisAligner)

//_______________________________________________________________________________________
AliZDCMisAligner::AliZDCMisAligner() : AliMisAligner()
{
  Printf("asdfasdfasdfasdf\n\n");

}

//_______________________________________________________________________________________
TClonesArray* AliZDCMisAligner::MakeAlObjsArray() {

  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;

  Double_t dx,dy,dz,dpsi,dtheta,dphi;
  if(TString(GetMisalType())=="ideal")
  {
    dx=0., dy=0., dz=0.;
    dpsi=0., dtheta=0., dphi=0.;
  }else if(TString(GetMisalType())=="residual" || TString(GetMisalType())=="full")
  {
    dx=0., dy=0.05, dz=0.;
    dpsi=0., dtheta=0., dphi=0.;
  }else{
    AliError(Form("\"%s\" is not a valid identifier for misalignment types. Exiting ...",GetMisalType()));
    return 0;
  }

  const char *ZDCCn="ZDC/NeutronZDC_C";
  const char *ZDCCp="ZDC/ProtonZDC_C";
  const char *ZDCAn="ZDC/NeutronZDC_A";
  const char *ZDCAp="ZDC/ProtonZDC_A";

  UShort_t iIndex=0;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  new(alobj[0]) AliAlignObjParams(ZDCCn, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[1]) AliAlignObjParams(ZDCCp, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[2]) AliAlignObjParams(ZDCAn, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
  new(alobj[3]) AliAlignObjParams(ZDCAp, volid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);

  return array;
}

//_______________________________________________________________________________________
AliCDBMetaData* AliZDCMisAligner::GetCDBMetaData() const {
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Chiara Oppedisano");

  if(TString(GetMisalType())=="ideal")
    md->SetComment("Alignment objects for ZDC ideal misalignment");
  if(TString(GetMisalType())=="residual")
    md->SetComment("Alignment objects for ZDC residual misalignment");
  if(TString(GetMisalType())=="full")
    md->SetComment("Alignment objects for ZDC full misalignment");
  return md;
}
