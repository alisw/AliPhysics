/**************************************************************************
 * Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved. *
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
/**
 * @file AliEmcalTriggerAlgorithm.cxx
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#include "AliEmcalTriggerAlgorithm.h"
#include "AliEmcalTriggerDataGrid.h"

#include <algorithm>


/// \cond CLASSIMP
templateClassImp(AliEmcalTriggerAlgorithm)
templateClassImp(AliEmcalJetTriggerAlgorithm)
templateClassImp(AliEmcalGammaTriggerAlgorithm)
/// \endcond

template<typename T>
AliEmcalTriggerAlgorithm<T>::AliEmcalTriggerAlgorithm():
TObject(),
fRowMin(0),
fRowMax(0),
fPatchSize(0),
fBitMask(0),
fThreshold(0)
{
}

template<typename T>
AliEmcalTriggerAlgorithm<T>::AliEmcalTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
TObject(),
fRowMin(rowmin),
fRowMax(rowmax),
fPatchSize(bitmask),
fBitMask(0),
fThreshold(0)
{
}

template<typename T>
AliEmcalTriggerAlgorithm<T>::~AliEmcalTriggerAlgorithm() {
}

template<typename T>
std::vector<AliEmcalTriggerRawPatch> AliEmcalTriggerAlgorithm<T>::FindPatches(const AliEmcalTriggerDataGrid<T> &adc) const {
  std::vector<AliEmcalTriggerRawPatch> result;
  T sumadc(0);
  for(int irow = fRowMin; irow < fRowMax - fPatchSize; irow += fPatchSize){
    for(int icol = 0; icol < adc.GetNumberOfCols() - fPatchSize; icol += fPatchSize){
      sumadc = 0;
      for(int jrow = irow; jrow < irow + fPatchSize; jrow++){
        for(int jcol = icol; jcol < icol + fPatchSize; jcol++){
          try{
            sumadc += adc(jcol, jrow);
          } catch (typename AliEmcalTriggerDataGrid<T>::OutOfBoundsException &e){

          }
        }
      }
      if(sumadc > fThreshold){
        AliEmcalTriggerRawPatch recpatch(icol, irow, fPatchSize, sumadc);
        recpatch.SetBitmask(fBitMask);
        result.push_back(recpatch);
      }
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

template<typename T>
AliEmcalJetTriggerAlgorithm<T>::AliEmcalJetTriggerAlgorithm():
AliEmcalTriggerAlgorithm<T>()
{
  this->SetPatchSize(16);
}

template<typename T>
AliEmcalJetTriggerAlgorithm<T>::AliEmcalJetTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
AliEmcalTriggerAlgorithm<T>(rowmin, rowmax, bitmask)
{
  this->SetPatchSize(16);
}

template<typename T>
AliEmcalJetTriggerAlgorithm<T>::~AliEmcalJetTriggerAlgorithm(){
}

template<typename T>
AliEmcalGammaTriggerAlgorithm<T>::AliEmcalGammaTriggerAlgorithm():
AliEmcalTriggerAlgorithm<T>()
{
  this->SetPatchSize(2);
}

template<typename T>
AliEmcalGammaTriggerAlgorithm<T>::AliEmcalGammaTriggerAlgorithm(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
AliEmcalTriggerAlgorithm<T>(rowmin, rowmax, bitmask)
{
  this->SetPatchSize(2);
}

template<typename T>
AliEmcalGammaTriggerAlgorithm<T>::~AliEmcalGammaTriggerAlgorithm(){
}

template class AliEmcalTriggerAlgorithm<int>;
template class AliEmcalTriggerAlgorithm<double>;
template class AliEmcalTriggerAlgorithm<float>;
template class AliEmcalJetTriggerAlgorithm<int>;
template class AliEmcalJetTriggerAlgorithm<double>;
template class AliEmcalJetTriggerAlgorithm<float>;
template class AliEmcalGammaTriggerAlgorithm<int>;
template class AliEmcalGammaTriggerAlgorithm<double>;
template class AliEmcalGammaTriggerAlgorithm<float>;
