#ifndef ALIEMCALCORRECTIONCLONECELLCONTAINER_H
#define ALIEMCALCORRECTIONCLONECELLCONTAINER_H

#include "AliEmcalCorrectionComponent.h"

class AliVCaloCells;

/**
 * @class AliEmcalCorrectionCellCloneContainer
 * @ingroup EMCALCORRECTIONFW
 * @brief Clones default cell list into new one to be used
 *
 * The internal cell branch should be set by the correction task! Then, the
 * external cell branch name should be set in the correction configuration.
 * The branch will be retrieved by this task.
 *
 * Note that to this task must run must choose when to execute it. 
 * One can imagine as an example an analysis with different time calibrations 
 * but the same energy and bad map as the ones applied to the default branch.  
 * For MC cross talk analysis, this component should run first. In any case, 
 * it should run before any cluster related correction.
 *    
 * Class based/copied from AliEmcalCorrectionCellCloneContainer
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, author of AliEmcalCorrectionCellCloneContainer
 * @author Gustavo Conesa Balbastre <gustavo.conesa.balbastre@cern.ch>, LPSC-Grenoble
 * @date May 1, 2020
 */

class AliEmcalCorrectionCellCloneContainer : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellCloneContainer();
  virtual ~AliEmcalCorrectionCellCloneContainer();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();

  std::string GetClonedCellsBranchName()                       const { return fClonedCellsBranchName ; }
  void SetClonedCellsBranchName(const char * cloneName)              { fClonedCellsBranchName = cloneName; }

 protected:

  void SetupClonedCells();
  void CloneCells();
  void AddCellsToClonedObject();
  void AddObjectToEvent(TObject *obj, AliVEvent * event, Bool_t attempt = kFALSE);

  std::string fClonedCellsBranchName;  ///<  Name of the cell branch which will be created for the combined cells.
  bool fInitializedClonedCells;        //!<! True if the created cells object has been initialized
  AliVCaloCells *fClonedCells;         //!<! Cells created from the input cells
  
  AliEmcalCorrectionCellCloneContainer(const AliEmcalCorrectionCellCloneContainer &);               // Not implemented
  AliEmcalCorrectionCellCloneContainer &operator=(const AliEmcalCorrectionCellCloneContainer &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellCloneContainer> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellCloneContainer, 1); // EMCal correction to combine cell collections
  /// \endcond
};

#endif /* AliEmcalCorrectionCellCloneContainer_H */
