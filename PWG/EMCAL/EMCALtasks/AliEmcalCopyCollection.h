#ifndef ALIEMCALCOPYCOLLECTION_H
#define ALIEMCALCOPYCOLLECTION_H

class AliVEvent;

#include "AliAnalysisTaskSE.h"
#include "AliEmcalCorrectionTask.h"

/**
 * Copies cell, cluster, or track collections for use in the
 * EMCal framework.
 *
 */

class AliEmcalCopyCollection : public AliAnalysisTaskSE {
 public:
  AliEmcalCopyCollection();
  AliEmcalCopyCollection(std::string name, AliEmcalCorrectionTask::InputObject_t inputObjectType, std::string collectionToCopyName, std::string newCollectionName, bool isEmbedding);
  virtual ~AliEmcalCopyCollection() {};

  // Methods from AliAnalysisTaskSE
  void UserCreateOutputObjects();
  void UserExec(Option_t * option);

  AliEmcalCorrectionTask::InputObject_t GetInputObjectType() const { return fInputObjectType; }
  std::string GetCollectionToCopyName() const { return fCollectionToCopyName; }
  std::string GetNewCollectionName() const { return fNewCollectionName; }
  bool GetIsEmbedding() const { return fIsEmbedding; }

  void SetInputObjectType(const AliEmcalCorrectionTask::InputObject_t inputObjectType) { fInputObjectType = inputObjectType; }
  void SetCollectionToCopyName(const std::string collectionToCopyName) { fCollectionToCopyName = collectionToCopyName; }
  void SetNewCollectionName(const std::string newCollectionName) { fNewCollectionName = newCollectionName; }
  void SetIsEmbedding(bool isEmbedding) { fIsEmbedding = isEmbedding; }

 protected:
  void                        CreateNewObjectBranch();
  void                        NewBranch();
  void                        CopyBranchToNewObject();
  void                        CopyClusters(TClonesArray *orig, TClonesArray *dest);

  AliEmcalCorrectionTask::InputObject_t fInputObjectType;  ///<
  std::string fCollectionToCopyName;                       ///<
  std::string fNewCollectionName;                          ///<
  bool fIsEmbedding;                                       ///<
  bool fEventInitialized;                                  //!<!
  bool fIsEsd;                                             //!<!
  AliVEvent * fEvent;                                      //!<!

  // TEMP
  std::string                 fCreatedCellBranchName;      ///< Name of created cell branch
  std::string                 fCreatedClusterBranchName;   ///< Name of created cluster branch
  std::string                 fCreatedTrackBranchName;     ///< Name of created track branch
  // ENDTEMP

  /// \cond CLASSIMP
  ClassDef(AliEmcalCopyCollection, 1); // Copy cell, cluster, or track collections
  /// \endcond
};

#endif /* ALIEMCALCOPYCOLLECTION_H */ 
