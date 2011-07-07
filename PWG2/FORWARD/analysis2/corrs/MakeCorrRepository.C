void
MakeCorrRepository()
{
 AliForwardCorrectionManager& mgr = AliForwardCorrectionManager::Instance();
 
 UInt_t what[] = { 
   AliForwardCorrectionManager::kSecondaryMap,
   AliForwardCorrectionManager::kELossFits,
   AliForwardCorrectionManager::kVertexBias,
   AliForwardCorrectionManager::kMergingEfficiency,
   AliForwardCorrectionManager::kDoubleHit,
   0
 };
 UInt_t* ptr = what;
 while (*ptr) { 
   TString dir(gSystem->ExpandPathName(mgr.GetFileDir(*ptr)));
   if (dir.IsNull()) { 
     ptr++;
     continue;
   }

   Info("MakeCorrRepository", "Making directory %s", dir.Data());
   gSystem->MakeDirectory(dir.Data());
   ptr++;
 }
}

