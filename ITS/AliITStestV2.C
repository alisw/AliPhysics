Int_t AliITStestV2() {
   Int_t rc=0;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClustersV2.C");
   if (rc=AliITSFindClustersV2()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindTracksV2.C");
   if (rc=AliITSFindTracksV2()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSComparisonV2.C");
   if (rc=AliITSComparisonV2()) return rc;

   return rc;
}
