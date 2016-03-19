/**
 * @file AliEMCALTriggerAlgorithm.cxx
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
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerAlgorithm.h"
#include <algorithm>


/// \cond CLASSIMP
templateClassImp(AliEMCALTriggerAlgorithm)
/// \endcond

template<typename T>
AliEMCALTriggerAlgorithm<T>::AliEMCALTriggerAlgorithm():
  TObject(),
  fRowMin(0),
  fRowMax(0),
  fPatchSize(0),
  fSubregionSize(1),
  fBitMask(0),
  fThreshold(0),
  fOfflineThreshold(0)
{
}

template<typename T>
AliEMCALTriggerAlgorithm<T>::AliEMCALTriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask):
  TObject(),
  fRowMin(rowmin),
  fRowMax(rowmax),
  fPatchSize(0),
  fSubregionSize(1),
  fBitMask(bitmask),
  fThreshold(0),
  fOfflineThreshold(0)
{
}

template<typename T>
AliEMCALTriggerAlgorithm<T>::~AliEMCALTriggerAlgorithm() {
}

template<typename T>
std::vector<AliEMCALTriggerRawPatch> AliEMCALTriggerAlgorithm<T>::FindPatches(const AliEMCALTriggerDataGrid<T> &adc, const AliEMCALTriggerDataGrid<T> &offlineAdc) const {
  std::vector<AliEMCALTriggerRawPatch> result;
  T sumadc(0);
  T sumofflineAdc(0);

  int rowStartMax = fRowMax - (fPatchSize-1);
  int colStartMax = adc.GetNumberOfCols() - fPatchSize;

  for(int irow = fRowMin; irow <= rowStartMax; irow += fSubregionSize){
    for(int icol = 0; icol <= colStartMax; icol += fSubregionSize){
      sumadc = 0;
      sumofflineAdc = 0;
      for(int jrow = irow; jrow < irow + fPatchSize; jrow++){
        for(int jcol = icol; jcol < icol + fPatchSize; jcol++){
          try{
            sumadc += adc(jcol, jrow);
            sumofflineAdc += offlineAdc(jcol, jrow);
          } catch (typename AliEMCALTriggerDataGrid<T>::OutOfBoundsException &e){

          }
        }
      }
      if(sumadc > fThreshold || sumofflineAdc > fOfflineThreshold){
        AliEMCALTriggerRawPatch recpatch(icol, irow, fPatchSize, sumadc, sumofflineAdc);
        recpatch.SetBitmask(fBitMask);
        result.push_back(recpatch);
      }
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

template class AliEMCALTriggerAlgorithm<int>;
template class AliEMCALTriggerAlgorithm<double>;
template class AliEMCALTriggerAlgorithm<float>;
