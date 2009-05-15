AliRsnPairManager *RsnConfigTest(const char *name="PHI->KK")
{
  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(name);

  AliRsnPairDef *defKPKM = new AliRsnPairDef(AliAODTrack::kKaon, '+', AliAODTrack::kKaon, '-');
  AliRsnPairDef *defKPKP = new AliRsnPairDef(AliAODTrack::kKaon, '+', AliAODTrack::kKaon, '+');
  AliRsnPairDef *defKMKM = new AliRsnPairDef(AliAODTrack::kKaon, '-', AliAODTrack::kKaon, '-');

  // =========== DEFINE PAIRS ==============

  // NoPID
  AliRsnPair *pairKPKM_N  = new AliRsnPair(AliRsnPair::kNoPID, defKPKM);
  AliRsnPair *pairKPKM_NS = new AliRsnPair(AliRsnPair::kNoPID, defKPKM);
  AliRsnPair *pairKPKM_NM = new AliRsnPair(AliRsnPair::kNoPIDMix, defKPKM/*, 5*/);
  AliRsnPair *pairKPKP_N  = new AliRsnPair(AliRsnPair::kNoPID, defKPKP);
  AliRsnPair *pairKMKM_N  = new AliRsnPair(AliRsnPair::kNoPID, defKMKM);
  // end NoPID

  // REALISTIC PID
  AliRsnPair *pairKPKM_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  AliRsnPair *pairKPKP_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKP);
  AliRsnPair *pairKMKM_R = new AliRsnPair(AliRsnPair::kRealisticPID, defKMKM);
  AliRsnPair *pairKPKM_RM = new AliRsnPair(AliRsnPair::kRealisticPIDMix, defKPKM/*,5*/);
  // end REALISTIC PID

  // REALISTIC PID Signal
  AliRsnPair *pairKPKM_RS = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  AliRsnPair *pairKPKM_RS2 = new AliRsnPair(AliRsnPair::kRealisticPID, defKPKM);
  // end REALISTIC PID Signal

  // PERFECT PID
  AliRsnPair *pairKPKM_P = new AliRsnPair(AliRsnPair::kPerfectPID, defKPKM);
  // end PERFECT PID
  // =========== END DEFINE PAIRS ==============

//   // =========== CUTS ==============
  ULong_t status = AliESDtrack::kESDpid;
  // status |=| AliESDtrack::kITSpid | AliESDtrack::kTPCpid;
  // status |= AliESDtrack::kTOFpid;

  AliRsnCut *cutStatus = new AliRsnCut("statusCut_ITS_TPC",
                                       "Status Cut for ITS and TPC",
                                       AliRsnCut::kStatus,
                                       status);
  AliRsnCut *cutESDLabelEqual = new AliRsnCut("cutESDLabelEqual",
                                              "cutESDLabelEqual",
                                              AliRsnCut::kIsLabelEqual);
  AliRsnCut *cutPtESDPart = new AliRsnCut("cutPtESDPart", "cutPtESDPart", AliRsnCut::kTransMomentum, 0.2, 10000.0);
  AliRsnCut *cutEtaESDPart = new AliRsnCut("cutESDEta", "cutESDEta", AliRsnCut::kEta, -0.9,0.9);

  AliRsnCutSet *cutSetParticleESD = new AliRsnCutSet("StatusCut");
//   cutSetParticleESD->AddCut(cutStatus);
//   cutSetParticleESD->SetCutScheme("statusCut_ITS_TPC");
  cutSetParticleESD->AddCut(cutPtESDPart);
  cutSetParticleESD->AddCut(cutEtaESDPart);
//   cutSetParticleESD->AddCut(cutStatus);
//   cutSetParticleESD->SetCutScheme("cutPtESDPart&statusCut_ITS_TPC");
  cutSetParticleESD->SetCutScheme("cutPtESDPart&cutESDEta");

  AliRsnCutSet *cutSetPairESD = new AliRsnCutSet("esdLabel");
  cutSetPairESD->AddCut(cutESDLabelEqual);
  cutSetPairESD->SetCutScheme("!cutESDLabelEqual");
//
  AliRsnCut *cutPtMCPart = new AliRsnCut("cutPtMCPart", "cutPtMCPart", AliRsnCut::kTransMomentumMC, 0.2, 10000.0);
  AliRsnCut *cutEtaMCPart = new AliRsnCut("cutMcEta", "cutMcEta", AliRsnCut::kEtaMC, -0.9,0.9);
  AliRsnCutSet *cutSetParticleMC = new AliRsnCutSet("MCCut");
  cutSetParticleMC->AddCut(cutPtMCPart);
  cutSetParticleMC->AddCut(cutEtaMCPart);
  cutSetParticleMC->SetCutScheme("cutPtMCPart&cutMcEta");
//
  AliRsnCutMgr *cutMgrESD = new AliRsnCutMgr("cutMgr", "Cut Mgr");
  cutMgrESD->SetCutSet(AliRsnCut::kParticle, cutSetParticleESD);
  cutMgrESD->SetCutSet(AliRsnCut::kPair, cutSetPairESD);
//
  AliRsnCut *cutIsTruePair = new AliRsnCut("cutIsTruePair", "cutIsTruePair", AliRsnCut::kIsTruePair, 333);
//
  AliRsnCutSet *cutSetPairSignal = new AliRsnCutSet("cutSetPairSignal");
  cutSetPairSignal->AddCut(cutIsTruePair);
  cutSetPairSignal->SetCutScheme("cutIsTruePair");
//
  AliRsnCutMgr *cutMgrSignal = new AliRsnCutMgr("cutMgrESDSignal","Cut Mgr Signal");
  cutMgrSignal->SetCutSet(AliRsnCut::kParticle,cutSetParticleESD);
  cutMgrSignal->SetCutSet(AliRsnCut::kParticle,cutSetParticleMC);
  cutMgrSignal->SetCutSet(AliRsnCut::kPair,cutSetPairSignal);

  AliRsnCutMgr *cutMgrSignal2 = new AliRsnCutMgr("cutMgrESDSignalNoCut","Cut Mgr Signal");
//   cutMgrSignal2->SetCutSet(AliRsnCut::kParticle,cutSetParticleESD);
  cutMgrSignal2->SetCutSet(AliRsnCut::kPair,cutSetPairSignal);

  AliRsnCutMgr *cutMgrSignalMC = new AliRsnCutMgr("cutMgrMCSignal","Cut Mgr Signal");
//   cutMgrSignalMC->SetCutSet(AliRsnCut::kParticle,cutSetParticleMC);
  cutMgrSignalMC->SetCutSet(AliRsnCut::kPair,cutSetPairSignal);

//   pairKPKM_N->SetCutMgr(cutMgrESD);
//   // pairKPKM_NS->SetCutMgr(cutMgrESD);
//   pairKPKM_NM->SetCutMgr(cutMgrESD);
//   pairKPKP_N->SetCutMgr(cutMgrESD);
//   pairKMKM_N->SetCutMgr(cutMgrESD);
  pairKPKM_R->SetCutMgr(cutMgrESD);
//   // pairKPKM_RS->SetCutMgr(cutMgrESD);
//   pairKPKM_RM->SetCutMgr(cutMgrESD);
  pairKPKP_R->SetCutMgr(cutMgrESD);
  pairKMKM_R->SetCutMgr(cutMgrESD);
//   pairKPKM_P->SetCutMgr(cutMgrESD);
  //
// //  pairKPKM_N->SetCutMgr(cutMgrSignal);
//   pairKPKM_NS->SetCutMgr(cutMgrSignal);
// //  pairKPKM_NM->SetCutMgr(cutMgrSignal);
// //  pairKPKP_N->SetCutMgr(cutMgrSignal);
// //  pairKMKM_N->SetCutMgr(cutMgrSignal);
// //  pairKPKM_R->SetCutMgr(cutMgrSignal);
  pairKPKM_RS->SetCutMgr(cutMgrSignal);
  pairKPKM_RS2->SetCutMgr(cutMgrSignal2);
//  pairKPKM_RM->SetCutMgr(cutMgrSignal);
//  pairKPKP_R->SetCutMgr(cutMgrSignal);
//  pairKMKM_R->SetCutMgr(cutMgrSignal);
  pairKPKM_P->SetCutMgr(cutMgrSignalMC);

//   // =========== END CUTS ==============
  
  // =========== FUNCTIONS ==============

  Double_t mom[11] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};
  Double_t eta[10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
//     AliRsnHistoDef *hdInvMass = new AliRsnHistoDef(1400, 0.6, 2.0);
  AliRsnHistoDef *hdInvMass = new AliRsnHistoDef(200,0.9,1.1);
  AliRsnHistoDef *hdInvMassMC = new AliRsnHistoDef(200,0.9,1.1);

  AliRsnHistoDef *hdPt = new AliRsnHistoDef(100,0.0,5.0);
  AliRsnHistoDef *hdEta = new AliRsnHistoDef(100,-2.0,2.0);
  
  Int_t num=0;
  AliRsnFunction *fcnIM[10];

  fcnIM[num] = new AliRsnFunction(AliRsnFunction::kInvMass,hdInvMass,kTRUE);
  fcnIM[num]->SetBinningCut(AliRsnCut::kTransMomentum, 11, mom,0);
  fcnIM[num]->SetBinningCut(AliRsnCut::kEta, 10, eta,1);
  num++;

//   fcnIMMC[num] = new AliRsnFunction(AliRsnFunction::kInvMassMC,hdInvMassMC, kTRUE);
//   fcnIMMC[num]->SetBinningCut(AliRsnCut::kTransMomentum, 9, mom);
//   num++;

//   fcnIM[num] = new AliRsnFunction(AliRsnFunction::kPtSpectrum,hdPt,kTRUE);
//   fcnIM[num]->SetBinningCut(AliRsnCut::kTransMomentum, 11, mom);
//   num++;
  
  fcnIM[num] = new AliRsnFunction(AliRsnFunction::kEtaSpectrum,hdEta,kTRUE);
//   fcnIM[num]->SetBinningCut(AliRsnCut::kEta, 11, eta);
  fcnIM[num]->SetBinningCut(AliRsnCut::kTransMomentum, 11, mom,0);
//   fcnIM[num]->SetBinningCut(AliRsnCut::kEta, 10, eta,1);
  num++;

  Int_t i;
  for (i=0;i<num ;i++)
  {
    pairKPKM_N->AddFunction(fcnIM[i]);
    pairKPKM_NS->AddFunction(fcnIM[i]);
    pairKPKM_NM->AddFunction(fcnIM[i]);
    pairKPKP_N->AddFunction(fcnIM[i]);
    pairKMKM_N->AddFunction(fcnIM[i]);
    pairKPKM_R->AddFunction(fcnIM[i]);
    pairKPKM_RS->AddFunction(fcnIM[i]);
    pairKPKM_RS2->AddFunction(fcnIM[i]);
    pairKPKM_RM->AddFunction(fcnIM[i]);
    pairKPKP_R->AddFunction(fcnIM[i]);
    pairKMKM_R->AddFunction(fcnIM[i]);
    pairKPKM_P->AddFunction(fcnIM[i]);

  }
  // =========== END FUNCTIONS =============

// //   pairMgr->AddPair(pairKPKM_N);
// //   pairMgr->AddPair(pairKPKM_NS);
// //   pairMgr->AddPair(pairKPKM_NM);
// //   pairMgr->AddPair(pairKPKP_N);
// //   pairMgr->AddPair(pairKMKM_N);
   pairMgr->AddPair(pairKPKM_R);
   pairMgr->AddPair(pairKPKM_RS);
   pairMgr->AddPair(pairKPKM_RS2);
// //   pairMgr->AddPair(pairKPKM_RM);
   pairMgr->AddPair(pairKPKP_R);
   pairMgr->AddPair(pairKMKM_R);
   pairMgr->AddPair(pairKPKM_P);
  return pairMgr;

}
