AliEmcalCopyCollection* AddTaskEmcalCopyCollection(AliEmcalContainerUtils::InputObject_t inputObjectType = AliEmcalContainerUtils::kNoDefinedInputObject,
                          const char * collectionToCopyName = "",
                          const char * newCollectionName = "",
                          bool isEmbedding = false)
{  
  return AliEmcalCopyCollection::AddTaskEmcalCopyCollection(inputObjectType, collectionToCopyName, newCollectionName, isEmbedding);
}
