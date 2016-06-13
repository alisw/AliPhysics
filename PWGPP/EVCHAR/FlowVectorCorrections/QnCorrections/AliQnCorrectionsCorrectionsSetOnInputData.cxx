/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/

/// \file AliQnCorrectionsCorrectionsSetOnInputData.cxx
/// \brief Set of corrections on input data class implementation

#include "AliQnCorrectionsCorrectionsSetOnInputData.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCorrectionsSetOnInputData);
/// \endcond

/// Default constructor
AliQnCorrectionsCorrectionsSetOnInputData::AliQnCorrectionsCorrectionsSetOnInputData() : TList() {

}

/// Default destructor
AliQnCorrectionsCorrectionsSetOnInputData::~AliQnCorrectionsCorrectionsSetOnInputData() {

}

/// Adds a new correction to the set.
///
/// The correction is incorporated in its proper place according to
/// its key
void AliQnCorrectionsCorrectionsSetOnInputData::AddCorrection(AliQnCorrectionsCorrectionOnInputData *correction) {
  if (IsEmpty()) {
    AddFirst(correction);
  }
  else if (correction->Before((AliQnCorrectionsCorrectionOnInputData *) First())) {
    AddFirst(correction);
  }
  else if (((AliQnCorrectionsCorrectionOnInputData *) Last())->Before(correction)) {
    AddLast(correction);
  }
  else {
    for (Int_t ix = 0; ix < GetEntries(); ix++) {
      if (!correction->Before(At(ix))) {
        AddAt(correction, ix-1);
      }
    }
  }
}

/// Fill the global list of correction steps
/// \param correctionlist (partial) global list of corrections ordered by correction key
void AliQnCorrectionsCorrectionsSetOnInputData::FillOverallCorrectionsList(TList *correctionlist) const {
  if (!IsEmpty()) {
    if (!correctionlist->IsEmpty()) {
      for (Int_t ix = 0; ix < GetEntries(); ix++) {
        if (correctionlist->FindObject(At(ix)->GetName()) != NULL) {
          /* already in the list, skip it */
          continue;
        }
        else {
          /* not in the list, include it in its proper place */
          if (At(ix)->Before((AliQnCorrectionsCorrectionStepBase *) correctionlist->First())) {
              correctionlist->AddFirst(At(ix));
          }
          else if (((AliQnCorrectionsCorrectionStepBase *) correctionlist->Last())->Before(At(ix))) {
            correctionlist->AddLast(At(ix));
          }
          else {
            for (Int_t jx = 0; jx < correctionlist->GetEntries(); jx++) {
              if (!At(ix)->Before((AliQnCorrectionsCorrectionStepBase *) correctionlist->At(jx))) {
                correctionlist->AddAt(At(ix), jx-1);
              }
            }
          }
        }
      }
    }
    else {
      /* the passed list is empty so we include all present corrections keeping the order */
      for (Int_t ix = 0; ix < GetEntries(); ix++)
        correctionlist->Add(At(ix));
    }
  }
}


