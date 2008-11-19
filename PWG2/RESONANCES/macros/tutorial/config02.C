AliRsnPairMgr *config02()
{
  AliRsnPairMgr  *pairMgr  = new AliRsnPairMgr("PHI->KK");

  AliRsnPairDef *defKPKM = new AliRsnPairDef(AliRsnPID::kKaon, '+', AliRsnPID::kKaon, '-');

  // =========== DEFINE PAIRS ==================
  AliRsnPair *pairKPKM_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  // =========== END DEFINE PAIRS ==============

  // =========== FUNCTIONS ==============
  AliRsnHistoDef *hdInvMass = new AliRsnHistoDef(200,0.9,1.1);
  AliRsnFunction *fcnIM = new AliRsnFunction(AliRsnFunction::kInvMass,hdInvMass/*,kTRUE*/);
  pairKPKM_R->AddFunction(fcnIM);
  // =========== END FUNCTIONS =================

  // =========== DEFINE CUTS ===================

  AliRsnCutMgr *cmMy = new AliRsnCutMgr("cmMy","My cut Manager");

  AliRsnCutSet *cutSetParticle = new AliRsnCutSet("csParticle");

  ULong_t status = AliESDtrack::kESDpid | AliESDtrack::kITSpid | AliESDtrack::kTPCpid;
//   status |= AliESDtrack::kTOFpid;
  AliRsnCut *cutStatus = new AliRsnCut("status","Status Cut",AliRsnCut::kStatus,status);

  cutSetParticle->AddCut(cutStatus);
  cutSetParticle->SetCutScheme("status");

  cmMy->SetCutSet(AliRsnCut::kParticle,cutSetParticle);

  pairKPKM_R->SetCutMgr(cmMy);
  // =========== END DEFINE CUTS ===============

  pairMgr->AddPair(pairKPKM_R);

  return pairMgr;
}
