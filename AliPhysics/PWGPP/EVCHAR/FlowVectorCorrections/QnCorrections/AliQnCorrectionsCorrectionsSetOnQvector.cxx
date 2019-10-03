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

/// \file AliQnCorrectionsCorrectionsSetOnQvector.cxx
/// \brief Set of corrections on Qn vector class implementation

#include "AliQnCorrectionsCorrectionsSetOnQvector.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCorrectionsSetOnQvector);
/// \endcond

/// Default constructor
AliQnCorrectionsCorrectionsSetOnQvector::AliQnCorrectionsCorrectionsSetOnQvector() : TList() {

}

/// Default destructor
AliQnCorrectionsCorrectionsSetOnQvector::~AliQnCorrectionsCorrectionsSetOnQvector() {

}

/// Adds a new correction to the set.
///
/// The correction is incorporated in its proper place according to
/// its key
void AliQnCorrectionsCorrectionsSetOnQvector::AddCorrection(AliQnCorrectionsCorrectionOnQvector *correction) {
  if (IsEmpty()) {
    AddFirst(correction);
  }
  else if (correction->Before((AliQnCorrectionsCorrectionOnQvector *) First())) {
    AddFirst(correction);
  }
  else if (((AliQnCorrectionsCorrectionOnQvector *) Last())->Before(correction)) {
    AddLast(correction);
  }
  else {
    for (Int_t ix = 0; ix < GetEntries(); ix++) {
      if (correction->Before(At(ix))) {
        AddAt(correction, ix);
        break;
      }
    }
  }
}

/// Fill the global list of correction steps
/// \param correctionlist (partial) global list of corrections ordered by correction key
void AliQnCorrectionsCorrectionsSetOnQvector::FillOverallCorrectionsList(TList *correctionlist) const {

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
              if (At(ix)->Before((AliQnCorrectionsCorrectionStepBase *) correctionlist->At(jx))) {
                correctionlist->AddAt(At(ix), jx);
                break;
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

/// Gets the correction on Qn vector previous to the one passed as argument
/// \param correction the correction to find the previous one
/// \return the previous correction, NULL if none
const AliQnCorrectionsCorrectionOnQvector *AliQnCorrectionsCorrectionsSetOnQvector::GetPrevious(const AliQnCorrectionsCorrectionOnQvector *correction) const {
  if (correction == NULL) return NULL;
  if (IsEmpty()) return NULL;
  if (First()->GetName() == correction->GetName()) return NULL;
  if (GetEntries() == 1) return NULL;
  for (Int_t ix = 0; ix < GetEntries() - 1; ix++) {
    if (At(ix+1)->GetName() == correction->GetName())
      return At(ix);
  }
  return NULL;
}

/// Check if a concrete correction step is bein applied on this detector configuration
/// It is not enough having the correction step configured or collecting data. To
/// get an affirmative answer the correction step must be being applied.
/// Transfer the order to each of the Qn correction steps.
/// \param step the name of the correction step
/// \return TRUE if the correction step is being applied
Bool_t AliQnCorrectionsCorrectionsSetOnQvector::IsCorrectionStepBeingApplied(const char *step) const {

  for (Int_t ix = 0; ix < GetEntries(); ix++) {
    if (TString(At(ix)->GetName()).Contains(step)) {
      return At(ix)->IsBeingApplied();
    }
  }
  return kFALSE;
}

