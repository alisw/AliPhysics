AliEmcalCopyCollection* AddTaskEmcalCopyCollection(AliEmcalCorrectionTask::InputObject_t inputObjectType = AliEmcalCorrectionTask::kNoDefinedInputObject,
                          const char * collectionToCopyName = "",
                          const char * newCollectionName = "",
                          bool isEmbedding = false)
{  
  return AliEmcalCopyCollection::AddTaskEmcalCopyCollection(inputObjectType, collectionToCopyName, newCollectionName, isEmbedding);
}
