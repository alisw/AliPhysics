AliRsnPairMgr *myRsnConfig(){
  // Create pair manager
  AliRsnPairMgr  *pairMgr  = new AliRsnPairMgr("PHI->KK");
  // Create pair definition "K+K-"
  AliRsnPairDef *defKPKM = new AliRsnPairDef(AliRsnPID::kKaon, '+', AliRsnPID::kKaon, '-');

  // =========== DEFINE PAIRS ==============
  // Create pair with Realistic PID
  AliRsnPair *pairKPKM_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  // =========== END DEFINE PAIRS ==============

  // =========== FUNCTIONS ==============
  // Define output histogram definition
  AliRsnHistoDef *hdInvMass = new AliRsnHistoDef(200,0.9,1.1);
  // Define output (we want inv.mass as output)
  AliRsnFunction *fcnIM = new AliRsnFunction(AliRsnFunction::kInvMass,hdInvMass/*,kTRUE*/);
  // and add it to the pair
  pairKPKM_R->AddFunction(fcnIM);
  // =========== END FUNCTIONS =============
  
  // Add pair to the pair manager
  pairMgr->AddPair(pairKPKM_R);

  // return pair manager
  return pairMgr;
}