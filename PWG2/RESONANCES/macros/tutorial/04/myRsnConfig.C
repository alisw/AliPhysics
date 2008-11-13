AliRsnPairMgr *myRsnConfig()
{
  AliRsnPairMgr  *pairMgr  = new AliRsnPairMgr("PHI->KK");

  AliRsnPairDef *defKPKM = new AliRsnPairDef(AliRsnPID::kKaon, '+', AliRsnPID::kKaon, '-');

  // =========== DEFINE PAIRS ==================
  AliRsnPair *pairKPKM_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  // =========== END DEFINE PAIRS ==============

  // =========== FUNCTIONS ==============
  AliRsnHistoDef *hdInvMass = new AliRsnHistoDef(200,0.9,1.1);
  AliRsnFunction *fcnIM = new AliRsnFunction(AliRsnFunction::kInvMass,hdInvMass);

  AliRsnHistoDef *hdPt = new AliRsnHistoDef(100,0.0,5.0);
  AliRsnFunction *fcnPt = new AliRsnFunction(AliRsnFunction::kPtSpectrum,hdPt,kTRUE);
  
  pairKPKM_R->AddFunction(fcnIM);
  pairKPKM_R->AddFunction(fcnPt);
  // =========== END FUNCTIONS =================

  // =========== DEFINE CUTS ===================
  ULong_t status = AliESDtrack::kESDpid | AliESDtrack::kITSpid | AliESDtrack::kTPCpid;
//   status |= AliESDtrack::kTOFpid;
  AliRsnCut *cutStatus = new AliRsnCut("cutStatus","Status Cut",AliRsnCut::kStatus,status);
  
  AliRsnCut *cutPtPart = new AliRsnCut("cutPtPair", "cutPtPair", AliRsnCut::kTransMomentum, 0.0,1.0);

  AliRsnCutSet *cutSetParticle = new AliRsnCutSet("csParticle");
  cutSetParticle->AddCut(cutStatus);
  cutSetParticle->SetCutScheme("status");
  
  AliRsnCutSet *csPtPair = new AliRsnCutSet("csPtPair");
  csPtPair->AddCut(cutPtPart);
  csPtPair->SetCutScheme("!cutPtPair");
  
  AliRsnCutMgr *cmMy = new AliRsnCutMgr("cmMy","My cut Manager");
  cmMy->SetCutSet(AliRsnCut::kParticle,cutSetParticle);
  cmMy->SetCutSet(AliRsnCut::kPair,csPtPair);
  
  pairKPKM_R->SetCutMgr(cmMy);
  // =========== END DEFINE CUTS ===============

  pairMgr->AddPair(pairKPKM_R);

  return pairMgr;
}
