#ifndef ALIEMCALCORRECTIONCELLCOMBINECOLLECTIONS_H
#define ALIEMCALCORRECTIONCELLCOMBINECOLLECTIONS_H

#include "AliEmcalCorrectionComponent.h"

class AliVCaloCells;

/**
 * @class AliEmcalCorrectionCellCombineCollections
 * @ingroup EMCALCORRECTIONFW
 * @brief Combines two cell collections into one collection
 *
 * The internal cell branch should be set by the correction task! Then, the
 * external cell branch name should be set in the correction configuration.
 * The branch will be retrieved by this task.
 *
 * Note that this task must run in a seperate correction task and must be
 * configured to run after the previous ones. Otherwise, the cell branches
 * will not yet be corrected!
 *
 * With a bit of generalization, along with some small modifications in
 * AliEmcalCorrectionTask, this task could combine aribtrary collections
 * of cells together. However, this is quite unnecessary for the current
 * purposes, so it will be left for the future if it is necessary.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Dec 6, 2016
 */

class AliEmcalCorrectionCellCombineCollections : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellCombineCollections();
  virtual ~AliEmcalCorrectionCellCombineCollections();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();

  std::string GetExternalCellsBranchName()                      const { return fExternalCellsBranchName; }
  std::string GetCombinedCellsBranchName()                      const { return fCreatedCellsBranchName; }
  
  void SetExternalCellsBranchName(const char * inputName)             { fExternalCellsBranchName = inputName; }
  void SetCombinedCellsBranchName(const char * inputName)             { fCreatedCellsBranchName = inputName; }

 protected:

  void SetupCombinedCells();
  void CreateCombinedCells();
  void AddCellsToCombinedCellObject(AliVCaloCells * inputCells, int indexOffset);
  void VerifyCombinedCells(std::vector <AliVCaloCells *> inputCaloCells);
  void AddObjectToEvent(TObject *obj, AliVEvent * event, Bool_t attempt = kFALSE);

  std::string fExternalCellsBranchName; ///<  Name of the cell branch which will be copied from the external event for the combined cells.
  std::string fCreatedCellsBranchName;  ///<  Name of the cell branch which will be created for the combined cells.
  bool fVerifyCombinedCells;            ///<  True if the task should confirm that the combined cells properly copied the input cells.
  bool fInitializedCombinedCells;       //!<! True if the combined cells object has been initialized
  AliVCaloCells *fCombinedCells;        //!<! Cells combined from the input and external events.
  
  AliEmcalCorrectionCellCombineCollections(const AliEmcalCorrectionCellCombineCollections &);               // Not implemented
  AliEmcalCorrectionCellCombineCollections &operator=(const AliEmcalCorrectionCellCombineCollections &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellCombineCollections> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellCombineCollections, 2); // EMCal correction to combine cell collections
  /// \endcond
};

#endif /* AliEmcalCorrectionCellCombineCollections_H */
