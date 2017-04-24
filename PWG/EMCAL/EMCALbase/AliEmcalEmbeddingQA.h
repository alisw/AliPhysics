#ifndef ALIEMCALEMBEDDINGQA_H 
#define ALIEMCALEMBEDDINGQA_H

#include <TObject.h>
class TList;

#include "THistManager.h"

/**
 * @class AliEmcalEmbeddingQA
 * @brief QA Class for EMCal Embedding Framework
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @date Apr 21, 2017
 * @ingroup EMCALFWTASKS
 *
 * This class includes QA information for the EMCal Embedding Framework.
 * It provides information such as the pythia trials, cross section, and
 * pt hard spectra corresponding to selected events. It compliments the
 * information in the main embedding task, which records that information
 * for _all_ events. It should be created as a member of the class which will
 * use embedding to record the embedding properties of events which are selected.
 * This task should always be somehow called when running the Embedding Framework!
 */

class AliEmcalEmbeddingQA : public TObject {
 public:
  AliEmcalEmbeddingQA();
  ~AliEmcalEmbeddingQA() {}

  bool IsInitialized() const { return fInitialized; }

  bool Initialize();
  bool AddQAPlotsToList(TList * list);
  void RecordEmbeddedEventProperties();

 private:
  bool fInitialized;              ///<  Notes whether the QA hists have been created
  THistManager fHistManager;      ///<  Hist manager

  /// \cond CLASSIMP
  ClassDef(AliEmcalEmbeddingQA, 1); // EMCal Embedding QA
  /// \endcond
};

#endif /* ALIEMCALEMBEDDINGQA_H */
