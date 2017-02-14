#include <AliEMCALTriggerPatchInfoV1.h>

ClassImp(AliEMCALTriggerPatchInfoV1)

AliEMCALTriggerPatchInfoV1::AliEMCALTriggerPatchInfoV1():
    AliEMCALTriggerPatchInfo(),
    fSmearedEnergyV1(0.)
{

}

Double_t AliEMCALTriggerPatchInfoV1::GetSmearedETV1() const {
  return GetET(fSmearedEnergyV1);
}

AliEMCALTriggerPatchInfoV1* AliEMCALTriggerPatchInfoV1::CreateAndInitializeV1(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, Double_t patchE, UInt_t bitmask, const TVector3& vertex, const AliEMCALGeometry* geom)
{
  AliEMCALTriggerPatchInfoV1* patch = new AliEMCALTriggerPatchInfoV1;
  patch->Initialize(col0, row0, size, adc, offlineAdc, patchE, bitmask, vertex, geom);
  return patch;
}

AliEMCALTriggerPatchInfoV1* AliEMCALTriggerPatchInfoV1::CreateAndInitializeV1(UChar_t col0, UChar_t row0, UChar_t size, UInt_t adc, UInt_t offlineAdc, UInt_t bitmask, const AliEMCALGeometry* geom)
{
  AliEMCALTriggerPatchInfoV1* patch = new AliEMCALTriggerPatchInfoV1;
  patch->Initialize(col0, row0, size, adc, offlineAdc, bitmask, geom);
  return patch;
}

