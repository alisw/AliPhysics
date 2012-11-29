#ifndef __CINT__
#include <Rtypes.h>
#endif
Int_t AddRsnDaughterCutsKStarNsigma_Syst(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,Bool_t isRsnMini=kFALSE,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);
   Int_t useCommonQualityCut = AliAnalysisManager::GetGlobalInt("rsnCommonQualityCut",valid);
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);

   Bool_t usePPCut = kFALSE;

   if (isPP && (opt.Contains("usePP"))) usePPCut = kTRUE;


   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 2;

   //---------------------------------------------
   //  Define single cutsP
   //---------------------------------------------

   Printf("AddRsnDaughterCutsKStarNsigma_Syst Option : %s",opt.Data());


   // default values
   Double_t nSigmaTPC_Pi=3.0;
   Double_t nSigmaTPC_K=3.0;
   Double_t nSigmaTOF_Pi=3.0;
   Double_t nSigmaTOF_K=3.0;
   Double_t etaRange=0.8;

   //Use single track Pt Cuts
   Double_t trackPtMin = 0.15;
   Double_t trackPtMax = 1.e20;
   Bool_t useTrackPtCut = kTRUE;

   if(opt.Contains("minPt02")) trackPtMin=0.2;
   if(opt.Contains("minPt03")) trackPtMin=0.3;
   if(opt.Contains("minPt04")) trackPtMin=0.4;
   if(opt.Contains("minPt05")) trackPtMin=0.5;
   if(opt.Contains("minPt06")) trackPtMin=0.6;
   if(opt.Contains("minPt10")) trackPtMin=1.0;
   if(opt.Contains("minPt15")) trackPtMin=1.5;


   //Use min TPC cluster Cut
   Int_t minclsK,maxclsK;
   Int_t minclsPi,maxclsPi;

   Bool_t useTPC_Pi=kFALSE;
   Bool_t useTOF_Pi=kFALSE;
   Bool_t useTPC_K=kFALSE;
   Bool_t useTOF_K=kFALSE;

   if (opt.Contains("qualityonly")) {
      useTPC_Pi=kFALSE;
      useTOF_Pi=kFALSE;
      useTPC_K=kFALSE;
      useTOF_K=kFALSE;
   } else if (!opt.Contains("nsig")) {
      useTPC_Pi=kTRUE;
      useTOF_Pi=kTRUE;
      useTPC_K=kTRUE;
      useTOF_K=kTRUE;
   }

   if (opt.Contains("PiTPCnsig")) useTPC_Pi=kTRUE;
   if (opt.Contains("PiTOFnsig")) useTOF_Pi=kTRUE;
   if (opt.Contains("KTPCnsig"))  useTPC_K=kTRUE;
   if (opt.Contains("KTOFnsig"))  useTOF_K=kTRUE;

   if (opt.Contains("PiTPCnsig10")) nSigmaTPC_Pi = 1.0;
   if (opt.Contains("PiTPCnsig15")) nSigmaTPC_Pi = 1.5;
   if (opt.Contains("PiTPCnsig20")) nSigmaTPC_Pi = 2.0;
   if (opt.Contains("PiTPCnsig25")) nSigmaTPC_Pi = 2.5;
   if (opt.Contains("PiTPCnsig30")) nSigmaTPC_Pi = 3.0;

   if (opt.Contains("KTPCnsig10")) nSigmaTPC_K = 1.0;
   if (opt.Contains("KTPCnsig15")) nSigmaTPC_K = 1.5;
   if (opt.Contains("KTPCnsig20")) nSigmaTPC_K = 2.0;
   if (opt.Contains("KTPCnsig25")) nSigmaTPC_K = 2.5;
   if (opt.Contains("KTPCnsig30")) nSigmaTPC_K = 3.0;

   if (opt.Contains("PiTOFnsig10")) nSigmaTOF_Pi = 1.0;
   if (opt.Contains("PiTOFnsig15")) nSigmaTOF_Pi = 1.5;
   if (opt.Contains("PiTOFnsig20")) nSigmaTOF_Pi = 2.0;
   if (opt.Contains("PiTOFnsig25")) nSigmaTOF_Pi = 2.5;
   if (opt.Contains("PiTOFnsig30")) nSigmaTOF_Pi = 3.0;

   if (opt.Contains("KTOFnsig10")) nSigmaTOF_K = 1.0;
   if (opt.Contains("KTOFnsig15")) nSigmaTOF_K = 1.5;
   if (opt.Contains("KTOFnsig20")) nSigmaTOF_K = 2.0;
   if (opt.Contains("KTOFnsig25")) nSigmaTOF_K = 2.5;
   if (opt.Contains("KTOFnsig30")) nSigmaTOF_K = 3.0;


   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
   }

   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   // KAON SETTINGS =======================================
   TString scheme="";
   TString cutname = "K_Kstar";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsK = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutK = new AliRsnCutTrackQuality("cutQuatityK");
   if (useCommonQualityCut>=0) {
      qualityCutK->SetAODTestFilterBit(useCommonQualityCut);

   } else {
      qualityCutK->SetDefaults2010();

   }

   //No filter bit
   if(opt.Contains("NOfb")) qualityCutK->SetAODTestFilterBit(-1);

   cutsK->AddCut(qualityCutK);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCutK->GetName();


   if (useTPC_K) {
      AliRsnCutPIDNSigma *cutKTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
      cutKTPC->SinglePIDRange(nSigmaTPC_K);
      cutsK->AddCut(cutKTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTPC->GetName();
   }

   if (useTOF_K) {
      AliRsnCutPIDNSigma *cutKTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFK",AliPID::kKaon,AliRsnCutPIDNSigma::kTOF);
      cutKTOF->SinglePIDRange(nSigmaTOF_K);
      cutsK->AddCut(cutKTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutKTOF->GetName();
   }
   if (useEta) {
      AliRsnValueDaughter *valEtaK = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kEta);
      AliRsnCutValue *cutEtaK = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type2),opt.Data()),-etaRange,etaRange);
      cutEtaK->SetTargetType(AliRsnTarget::kDaughter);
      cutEtaK->SetValueObj(valEtaK);
      cutsK->AddCut(cutEtaK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEtaK->GetName();
   }
   if (usePDG) {
      AliRsnCutPID *cutPDGK = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type2),opt.Data()),type2,0.0,kTRUE);
      cutsK->AddCut(cutPDGK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPDGK->GetName();
   }



   //MinPt cut kaon
   if (useTrackPtCut) {
      Printf("Adding Pt min=%.3f max=%.3f ...",trackPtMin,trackPtMax);
      AliRsnValueDaughter *valTrackPtK = new AliRsnValueDaughter(Form("val%sTrackPt%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kPt);

      AliRsnCutValue *cutTrackPtK = new AliRsnCutValue(Form("cut%sTrackPt%s",AliPID::ParticleName(type2),opt.Data()),trackPtMin,trackPtMax);
      cutTrackPtK->SetTargetType(AliRsnTarget::kDaughter);
      cutTrackPtK->SetValueObj(valTrackPtK);
      cutsK->AddCut(cutTrackPtK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutTrackPtK->GetName();
   }

   //Ncluster cut kaon
   if(opt.Contains("tpcncl80K")) {
      Printf("***** adding 80 TPCNCL cut Kaon");
      AliRsnValueDaughter *val_tpcnclK = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclK = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),80,10000);
      cut_tpcnclK->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclK->SetValueObj(val_tpcnclK);
      cutsK->AddCut(cut_tpcnclK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclK->GetName();
   }

   if(opt.Contains("tpcncl90K")) {
      Printf("***** adding 90 TPCNCL cut Kaon");
      AliRsnValueDaughter *val_tpcnclK = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclK = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),90,10000);
      cut_tpcnclK->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclK->SetValueObj(val_tpcnclK);
      cutsK->AddCut(cut_tpcnclK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclK->GetName();
   }

   if(opt.Contains("tpcncl100K")) {
      Printf("***** adding 100 TPCNCL cut Kaon");
      AliRsnValueDaughter *val_tpcnclK = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclK = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type2),opt.Data()),100,10000);
      cut_tpcnclK->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclK->SetValueObj(val_tpcnclK);
      cutsK->AddCut(cut_tpcnclK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclK->GetName();
   }


   //Ncluster cut kaon through AliRsnCutTrackQuality
   if(opt.Contains("QTPCnclK")) {
      AliRsnCutTrackQuality *QTPCNclsCutK = new AliRsnCutTrackQuality("QTPCnclK");
      QTPCNclsCutK->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
      QTPCNclsCutK->SetAODTestFilterBit(5);//reset the filter bit cut
      QTPCNclsCutK->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually, not just the filter bit
      QTPCNclsCutK->SetPtRange(0.15,1.e20);//reset the pT cut
      QTPCNclsCutK->SetEtaRange(-0.8,0.8);//reset the eta cut

      if(opt.Contains("nclK70")) minclsK=70;
      if(opt.Contains("nclK75")) minclsK=75;
      if(opt.Contains("nclK80")) minclsK=80;
      if(opt.Contains("nclK85")) minclsK=85;
      if(opt.Contains("nclK90")) minclsK=90;
      if(opt.Contains("nclK100")) minclsK=100;

      Printf(Form("++++++++ Adding Cut: NclustersTPC Kaon >= %d",minclsK));
      QTPCNclsCutK->SetTPCminNClusters(minclsK);

      cutsK->AddCut(QTPCNclsCutK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += QTPCNclsCutK->GetName();
   }


   //pt dep dcaxy cut on kaon
   if(opt.Contains("PtDCAK")) {
      AliRsnCutTrackQuality *dcaxyCutK = new AliRsnCutTrackQuality("ptdcaK");
      dcaxyCutK->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
      dcaxyCutK->SetAODTestFilterBit(5);//reset the filter bit cut
      dcaxyCutK->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually, not just the filter bit
      dcaxyCutK->SetPtRange(0.15,1.e20);//reset the pT cut
      dcaxyCutK->SetEtaRange(-0.8,0.8);//reset the eta cut
      if(opt.Contains("DCAK7s")) {dcaxyCutK->SetDCARPtFormula("0.0182+0.0350/pt^1.01");}
      if(opt.Contains("DCAK6s")) {dcaxyCutK->SetDCARPtFormula("0.0156+0.03/pt^1.01");}
      if(opt.Contains("DCAK5s")) {dcaxyCutK->SetDCARPtFormula("0.013+0.025/pt^1.01");}
      if(opt.Contains("DCAK4s")) {dcaxyCutK->SetDCARPtFormula("0.0104+0.02/pt^1.01");}
      if(opt.Contains("DCAK3s")) {dcaxyCutK->SetDCARPtFormula("0.0078+0.015/pt^1.01");}
      if(opt.Contains("DCAK2s")) {dcaxyCutK->SetDCARPtFormula("0.0052+0.01/pt^1.01");}
      if(opt.Contains("DCAK1s")) {dcaxyCutK->SetDCARPtFormula("0.0026+0.005/pt^1.01");}
      cutsK->AddCut(dcaxyCutK);
      if (!scheme.IsNull()) scheme += "&";
      scheme += dcaxyCutK->GetName();
   }

   Printf ("CUT Scheme for KAON is '%s'",scheme.Data());
   cutsK->SetCutScheme(scheme.Data());

   // END KAON =======================================

   // Pion SETTINGS ===========================================

   scheme="";
   cutname = "Pi_Kstar";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsP = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutPi = new AliRsnCutTrackQuality("cutQuatityPi");
   if (useCommonQualityCut>=0) {
      qualityCutPi->SetAODTestFilterBit(useCommonQualityCut);

   } else {
      qualityCutPi->SetDefaults2010();
   }
   //No filter bit
   if(opt.Contains("NOfb")) qualityCutPi->SetAODTestFilterBit(-1);

   cutsP->AddCut(qualityCutPi);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCutPi->GetName();
   if (useTPC_Pi) {
      AliRsnCutPIDNSigma *cutPiTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCPi",AliPID::kPion,AliRsnCutPIDNSigma::kTPC);
      cutPiTPC->SinglePIDRange(nSigmaTPC_Pi);
      cutsP->AddCut(cutPiTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPiTPC->GetName();
   }
   if (useTOF_Pi) {
      AliRsnCutPIDNSigma *cutPiTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFPi",AliPID::kPion,AliRsnCutPIDNSigma::kTOF);
      cutPiTOF->SinglePIDRange(nSigmaTOF_Pi);
      cutsP->AddCut(cutPiTOF);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPiTOF->GetName();
   }
   if (useEta) {
      AliRsnValueDaughter *valEtaP = new AliRsnValueDaughter(Form("val%sETA%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kEta);
      AliRsnCutValue *cutEtaP = new AliRsnCutValue(Form("cut%sETA%s",AliPID::ParticleName(type1),opt.Data()),-etaRange,etaRange);
      cutEtaP->SetTargetType(AliRsnTarget::kDaughter);
      cutEtaP->SetValueObj(valEtaP);
      cutsP->AddCut(cutEtaP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutEtaP->GetName();
   }
   if (usePDG) {
      AliRsnCutPID *cutPDGP = new AliRsnCutPID(Form("cut%sPDG%s",AliPID::ParticleName(type1),opt.Data()),type1,0.0,kTRUE);
      cutsP->AddCut(cutPDGP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPDGP->GetName();
   }

   //MinPt cut pion
   if (useTrackPtCut) {
      Printf("Adding Pt min=%.3f max=%.3f ...",trackPtMin,trackPtMax);
      AliRsnValueDaughter *valTrackPtP = new AliRsnValueDaughter(Form("val%sTrackPt%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kPt);

      AliRsnCutValue *cutTrackPtP = new AliRsnCutValue(Form("cut%sTrackPt%s",AliPID::ParticleName(type1),opt.Data()),trackPtMin,trackPtMax);
      cutTrackPtP->SetTargetType(AliRsnTarget::kDaughter);
      cutTrackPtP->SetValueObj(valTrackPtP);
      cutsP->AddCut(cutTrackPtP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutTrackPtP->GetName();
   }

   //Ncluster cut pion
   if(opt.Contains("tpcncl80Pi")) {
      Printf("***** adding 80 TPCNCL cut Pion");
      AliRsnValueDaughter *val_tpcnclP = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclP = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),80,10000);
      cut_tpcnclP->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclP->SetValueObj(val_tpcnclP);
      cutsP->AddCut(cut_tpcnclP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclP->GetName();
   }

   if(opt.Contains("tpcncl90Pi")) {
      Printf("***** adding 90 TPCNCL cut Pion");
      AliRsnValueDaughter *val_tpcnclP = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclP = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),90,10000);
      cut_tpcnclP->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclP->SetValueObj(val_tpcnclP);
      cutsP->AddCut(cut_tpcnclP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclP->GetName();
   }


   if(opt.Contains("tpcncl100Pi")) {
      Printf("***** adding 100 TPCNCL cut Pion");
      AliRsnValueDaughter *val_tpcnclP = new AliRsnValueDaughter(Form("val%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kNTPCclusters);
      AliRsnCutValue *cut_tpcnclP = new AliRsnCutValue(Form("cut%s_tpcncl_%s",AliPID::ParticleName(type1),opt.Data()),100,10000);
      cut_tpcnclP->SetTargetType(AliRsnTarget::kDaughter);
      cut_tpcnclP->SetValueObj(val_tpcnclP);
      cutsP->AddCut(cut_tpcnclP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cut_tpcnclP->GetName();
   }

   //Ncluster cut on pion through AliRsnCutTrackQuality
   if(opt.Contains("QTPCnclPi")) {
      AliRsnCutTrackQuality *QTPCNclsCutPi = new AliRsnCutTrackQuality("QTPCnclPi");
      QTPCNclsCutPi->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
      QTPCNclsCutPi->SetAODTestFilterBit(5);//reset the filter bit cut
      QTPCNclsCutPi->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually, not just the filter bit
      QTPCNclsCutPi->SetPtRange(0.15,1.e20);//reset the pT cut
      QTPCNclsCutPi->SetEtaRange(-0.8,0.8);//reset the eta cut

      if(opt.Contains("nclPi70")) minclsPi=70;
      if(opt.Contains("nclPi75")) minclsPi=75;
      if(opt.Contains("nclPi80")) minclsPi=80;
      if(opt.Contains("nclPi85")) minclsPi=85;
      if(opt.Contains("nclPi90")) minclsPi=90;
      if(opt.Contains("nclPi100")) minclsPi=100;

      Printf(Form("+++++++++ Adding Cut: NclustersTPC Pion >= %d",minclsPi));
      QTPCNclsCutPi->SetTPCminNClusters(minclsPi);

      cutsP->AddCut(QTPCNclsCutPi);
      if (!scheme.IsNull()) scheme += "&";
      scheme += QTPCNclsCutPi->GetName();

   }

   //pt dep dcaxy cut pion
   if(opt.Contains("PtDCAP")) {
      AliRsnCutTrackQuality *dcaxyCutP = new AliRsnCutTrackQuality("ptdcaP6s");
      dcaxyCutP->DisableAll();//disable all cuts, filter bit, pT, eta, and DCAxy cuts will be reset later
      dcaxyCutP->SetAODTestFilterBit(5);//reset the filter bit cut
      dcaxyCutP->SetCheckOnlyFilterBit(kFALSE);//tells the cut object to check all other cuts individually, not just the filter bit
      dcaxyCutP->SetPtRange(0.15,1.e20);//reset the pT cut
      dcaxyCutP->SetEtaRange(-0.8,0.8);//reset the eta cut
      if(opt.Contains("DCAP7s")) {dcaxyCutP->SetDCARPtFormula("0.0182+0.0350/pt^1.01");}
      if(opt.Contains("DCAP6s")) {dcaxyCutP->SetDCARPtFormula("0.0156+0.03/pt^1.01");}
      if(opt.Contains("DCAP5s")) {dcaxyCutP->SetDCARPtFormula("0.013+0.025/pt^1.01");}
      if(opt.Contains("DCAP4s")) {dcaxyCutP->SetDCARPtFormula("0.0104+0.02/pt^1.01");}
      if(opt.Contains("DCAP3s")) {dcaxyCutP->SetDCARPtFormula("0.0078+0.015/pt^1.01");}
      if(opt.Contains("DCAP2s")) {dcaxyCutP->SetDCARPtFormula("0.0052+0.01/pt^1.01");}
      if(opt.Contains("DCAP1s")) {dcaxyCutP->SetDCARPtFormula("0.0026+0.005/pt^1.01");}
      cutsP->AddCut(dcaxyCutP);
      if (!scheme.IsNull()) scheme += "&";
      scheme += dcaxyCutP->GetName();
   }

   Printf ("CUT Scheme for PION is '%s'",scheme.Data());
   cutsP->SetCutScheme(scheme.Data());

   // END PION =======================================

   if (opt.Contains("mon")) {
      AddMonitorOutput(cutsP->GetMonitorOutput(),opt);
      AddMonitorOutput(cutsK->GetMonitorOutput(),opt);
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutsK);
         taskRsnMini->AddTrackCuts(cutsP);

      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
//       sel->SetLabelCheck(kFALSE);
      sel->Add(cutsP, kTRUE);
      sel->Add(cutsK, kTRUE);
   }
   return numberOfCuts;

}
