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
 * @file AliEmcalTriggerPatchFinder.cxx
 * @date Oct. 23, 2015
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 */
#include "AliEMCALTriggerPatchFinder.h"
#include "AliEMCALTriggerAlgorithm.h"
#include "AliEMCALTriggerDataGrid.h"

#include <vector>

/// \cond CLASSIMP
templateClassImp(AliEMCALTriggerPatchFinder)
/// \endcond

template<typename T>
AliEMCALTriggerPatchFinder<T>::AliEMCALTriggerPatchFinder():
  TObject(),
  fTriggerAlgorithms()
{
}

template<typename T>
AliEMCALTriggerPatchFinder<T>::~AliEMCALTriggerPatchFinder() {
  for(typename std::vector<AliEMCALTriggerAlgorithm<T> *>::iterator algiter = fTriggerAlgorithms.begin();
      algiter != fTriggerAlgorithms.end();
      ++algiter)
    delete *algiter;
}

template<typename T>
std::vector<AliEMCALTriggerRawPatch> AliEMCALTriggerPatchFinder<T>::FindPatches(const AliEMCALTriggerDataGrid<T> &adc) const{
  std::vector<AliEMCALTriggerRawPatch> result;
  for(typename std::vector<AliEMCALTriggerAlgorithm<T> *>::const_iterator algiter = fTriggerAlgorithms.begin();
      algiter != fTriggerAlgorithms.end();
      ++algiter)
  {
    std::vector<AliEMCALTriggerRawPatch> tmp = (*algiter)->FindPatches(adc);
    for(std::vector<AliEMCALTriggerRawPatch>::iterator patchiter = tmp.begin(); patchiter != tmp.end(); ++patchiter){
      result.push_back(*patchiter);
    }
  }
  return result;
}

template class AliEMCALTriggerPatchFinder<int>;
template class AliEMCALTriggerPatchFinder<double>;
template class AliEMCALTriggerPatchFinder<float>;
