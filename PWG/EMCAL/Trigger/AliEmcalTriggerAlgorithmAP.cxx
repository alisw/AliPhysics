/**
 * @file AliEmcalTriggerAlgorithmAP.cxx
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
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
#include "AliEmcalTriggerDataGridAP.h"
#include "AliEmcalTriggerAlgorithmAP.h"
#include <algorithm>


/// \cond CLASSIMP
templateClassImp(AliEmcalTriggerAlgorithmAP)
templateClassImp(AliEmcalJetTriggerAlgorithmAP)
templateClassImp(AliEmcalGammaTriggerAlgorithmAP)
/// \endcond

template<typename T>
AliEmcalTriggerAlgorithmAP<T>::AliEmcalTriggerAlgorithmAP():
TObject(),
fRowMin(0),
fRowMax(0),
fPatchSize(0),
fBitMask(0),
fThreshold(0),
fOfflineThreshold(0)
{
}

template<typename T>
AliEmcalTriggerAlgorithmAP<T>::AliEmcalTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
TObject(),
fRowMin(rowmin),
fRowMax(rowmax),
fPatchSize(bitmask),
fBitMask(0),
fThreshold(0),
fOfflineThreshold(0)
{
}

template<typename T>
AliEmcalTriggerAlgorithmAP<T>::~AliEmcalTriggerAlgorithmAP() {
}

template<typename T>
std::vector<AliEmcalTriggerRawPatchAP> AliEmcalTriggerAlgorithmAP<T>::FindPatches(const AliEmcalTriggerDataGridAP<T> &adc, const AliEmcalTriggerDataGridAP<T> &offlineAdc) const {
  std::vector<AliEmcalTriggerRawPatchAP> result;
  T sumadc(0);
  T sumofflineAdc(0);
  for(int irow = fRowMin; irow < fRowMax - fPatchSize; irow += fPatchSize){
    for(int icol = 0; icol < adc.GetNumberOfCols() - fPatchSize; icol += fPatchSize){
      sumadc = 0;
      sumofflineAdc = 0;
      for(int jrow = irow; jrow < irow + fPatchSize; jrow++){
        for(int jcol = icol; jcol < icol + fPatchSize; jcol++){
          try{
            sumadc += adc(jcol, jrow);
            sumofflineAdc += offlineAdc(jcol, jrow);
          } catch (typename AliEmcalTriggerDataGridAP<T>::OutOfBoundsException &e){

          }
        }
      }
      if(sumadc > fThreshold || sumofflineAdc > fOfflineThreshold){
        AliEmcalTriggerRawPatchAP recpatch(icol, irow, fPatchSize, sumadc, sumofflineAdc);
        recpatch.SetBitmask(fBitMask);
        result.push_back(recpatch);
      }
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

template<typename T>
AliEmcalJetTriggerAlgorithmAP<T>::AliEmcalJetTriggerAlgorithmAP():
AliEmcalTriggerAlgorithmAP<T>()
{
  this->SetPatchSize(16);
}

template<typename T>
AliEmcalJetTriggerAlgorithmAP<T>::AliEmcalJetTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
AliEmcalTriggerAlgorithmAP<T>(rowmin, rowmax, bitmask)
{
  this->SetPatchSize(16);
}

template<typename T>
AliEmcalJetTriggerAlgorithmAP<T>::~AliEmcalJetTriggerAlgorithmAP(){
}

template<typename T>
AliEmcalGammaTriggerAlgorithmAP<T>::AliEmcalGammaTriggerAlgorithmAP():
AliEmcalTriggerAlgorithmAP<T>()
{
  this->SetPatchSize(2);
}

template<typename T>
AliEmcalGammaTriggerAlgorithmAP<T>::AliEmcalGammaTriggerAlgorithmAP(Int_t rowmin, Int_t rowmax, ULong_t bitmask):
AliEmcalTriggerAlgorithmAP<T>(rowmin, rowmax, bitmask)
{
  this->SetPatchSize(2);
}

template<typename T>
AliEmcalGammaTriggerAlgorithmAP<T>::~AliEmcalGammaTriggerAlgorithmAP(){
}

template class AliEmcalTriggerAlgorithmAP<int>;
template class AliEmcalTriggerAlgorithmAP<double>;
template class AliEmcalTriggerAlgorithmAP<float>;
template class AliEmcalJetTriggerAlgorithmAP<int>;
template class AliEmcalJetTriggerAlgorithmAP<double>;
template class AliEmcalJetTriggerAlgorithmAP<float>;
template class AliEmcalGammaTriggerAlgorithmAP<int>;
template class AliEmcalGammaTriggerAlgorithmAP<double>;
template class AliEmcalGammaTriggerAlgorithmAP<float>;
