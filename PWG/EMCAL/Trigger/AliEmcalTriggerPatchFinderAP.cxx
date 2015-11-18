/**
 * @file AliEmcalTriggerPatchFinderAP.cxx
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
#include "AliEmcalTriggerAlgorithmAP.h"
#include "AliEmcalTriggerDataGridAP.h"
#include "AliEmcalTriggerPatchFinderAP.h"
#include <vector>

/// \cond CLASSIMP
templateClassImp(AliEmcalTriggerPatchFinderAP)
/// \endcond

template<typename T>
AliEmcalTriggerPatchFinderAP<T>::AliEmcalTriggerPatchFinderAP():
  TObject(),
  fTriggerAlgorithms()
{
}

template<typename T>
AliEmcalTriggerPatchFinderAP<T>::~AliEmcalTriggerPatchFinderAP() {
  for(typename std::vector<AliEmcalTriggerAlgorithmAP<T> *>::iterator algiter = fTriggerAlgorithms.begin();
      algiter != fTriggerAlgorithms.end();
      ++algiter)
    delete *algiter;
}

template<typename T>
std::vector<AliEmcalTriggerRawPatchAP> AliEmcalTriggerPatchFinderAP<T>::FindPatches(const AliEmcalTriggerDataGridAP<T> &adc, const AliEmcalTriggerDataGridAP<T> &offlineAdc) const{
  std::vector<AliEmcalTriggerRawPatchAP> result;
  for(typename std::vector<AliEmcalTriggerAlgorithmAP<T> *>::const_iterator algiter = fTriggerAlgorithms.begin();
      algiter != fTriggerAlgorithms.end();
      ++algiter)
  {
    std::vector<AliEmcalTriggerRawPatchAP> tmp = (*algiter)->FindPatches(adc, offlineAdc);
    for(std::vector<AliEmcalTriggerRawPatchAP>::iterator patchiter = tmp.begin(); patchiter != tmp.end(); ++patchiter){
      result.push_back(*patchiter);
    }
  }
  return result;
}

template class AliEmcalTriggerPatchFinderAP<int>;
template class AliEmcalTriggerPatchFinderAP<double>;
template class AliEmcalTriggerPatchFinderAP<float>;
