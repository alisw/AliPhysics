#ifndef ALIEMCALCOPYCOLLECTION_H
#define ALIEMCALCOPYCOLLECTION_H

class AliVEvent;

#include <AliAnalysisTaskSE.h>

#include "AliEmcalContainerUtils.h"

/**
 * \class AliEmcalCopyCollection
 * \brief Copies cell, cluster, or track collections for use in the EMCal framework.
 * \ingroup EMCALCOREFW
 *
 * Often, it is necessary to copy cells, clusters, or tracks for use within the EMCal
 * framework. For example, it is important when comparing whether two sets of corrections
 * yield the same results. More generally, it is helpful when running more than one set of
 * tasks which need the same objects with different settings.
 *
 * This task consolidates the code for copying cells, clusters, and tracks into one place.
 * Some code is drawn from AliEmcalAodTrackFilterTask, and AliAnalysisTaskEMCALClusterizeFast.
 *
 * \author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * \date Nov 3, 2016
 * 
 */

class AliEmcalCopyCollection : public AliAnalysisTaskSE {
 public:
  AliEmcalCopyCollection();
  AliEmcalCopyCollection(std::string name, AliEmcalContainerUtils::InputObject_t inputObjectType, std::string collectionToCopyName, std::string newCollectionName, bool isEmbedding);
  virtual ~AliEmcalCopyCollection() {};

  // Methods from AliAnalysisTaskSE
  void UserCreateOutputObjects();
  void UserExec(Option_t * option);

  // Get methods
  AliEmcalContainerUtils::InputObject_t GetInputObjectType() const { return fInputObjectType; }
  std::string GetCollectionToCopyName() const { return fCollectionToCopyName; }
  std::string GetNewCollectionName() const { return fNewCollectionName; }
  bool GetIsEmbedding() const { return fIsEmbedding; }

  // Set methods
  void SetInputObjectType(const AliEmcalContainerUtils::InputObject_t inputObjectType) { fInputObjectType = inputObjectType; }
  void SetCollectionToCopyName(const std::string collectionToCopyName) { fCollectionToCopyName = collectionToCopyName; }
  void SetNewCollectionName(const std::string newCollectionName) { fNewCollectionName = newCollectionName; }
  void SetIsEmbedding(bool isEmbedding) { fIsEmbedding = isEmbedding; }

  // Add Task
  static AliEmcalCopyCollection* AddTaskEmcalCopyCollection(AliEmcalContainerUtils::InputObject_t inputObjectType = AliEmcalContainerUtils::kNoDefinedInputObject,
                          const char * collectionToCopyName = "",
                          const char * newCollectionName = "",
                          bool isEmbedding = false);

 protected:
  void                        CreateNewObjectBranch();
  void                        NewBranch();
  void                        CopyBranchToNewObject();
  void                        CopyClusters(TClonesArray *orig, TClonesArray *dest);

  AliEmcalContainerUtils::InputObject_t fInputObjectType;  ///< Type of collection to copy
  std::string fCollectionToCopyName;                       ///< Name of the collection branch to copy
  std::string fNewCollectionName;                          ///< Name of the collection bracnh where it will be copied
  bool fIsEmbedding;                                       ///< Denotes whether the collection is embedded (and therefore should be copied and stored in the external event)
  bool fEventInitialized;                                  //!<! Denotes whether the event is initialized
  bool fIsEsd;                                             //!<! Denotes whether the input is ESD
  AliVEvent * fEvent;                                      //!<! Points to the event to copy from and store the new bracnh in

  /// \cond CLASSIMP
  ClassDef(AliEmcalCopyCollection, 1); // Copy cell, cluster, or track collections
  /// \endcond
};

#endif /* ALIEMCALCOPYCOLLECTION_H */ 
