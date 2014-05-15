#ifndef __CINT__
#endif
Int_t AddRsnDaughterCutsLambdaNsigma(AliPID::EParticleType type1,AliPID::EParticleType type2,TString opt,AliRsnInputHandler *rsnIH=0,AliAnalysisTaskSE *task=0)
{

   if (!rsnIH) return 0;

   Bool_t valid = kTRUE;
   //   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t useCommonQualityCut = AliRsnTrainManager::GetGlobalInt("RsnCommonQualityCut",valid);
   TString rsnQualityCut = AliRsnTrainManager::GetGlobalStr("RsnQualityCut",valid);
   Int_t isMC = AliRsnTrainManager::GetGlobalInt("IsMC",valid);
   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);
   Int_t isMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);

   // === USER HAS TO SET CORRECT NUMBER OF CUTS SETS =====
   Int_t numberOfCuts = 2;

   //---------------------------------------------
   //  Define single cutsP
   //---------------------------------------------

   Printf("AddRsnDaughterCutsLambda Option : %s",opt.Data());

   Double_t nSigmaTPC_P=3.0;
   Double_t nSigmaTPCmin_P=3.0;
   Double_t nSigmaTPCmax_P=3.0;
   Double_t nSigmaTPC_K=3.0;
   Double_t nSigmaTPCmin_K=3.0;
   Double_t nSigmaTPCmax_K=3.0;
   Double_t nSigmaTOF_P=3.0;
   Double_t nSigmaTOFmin_P=3.0;
   Double_t nSigmaTOFmax_P=3.0;
   Double_t nSigmaTOF_K=3.0;
   Double_t nSigmaTOFmin_K=3.0;
   Double_t nSigmaTOFmax_K=3.0;

   Double_t nSigmaTPCmin_PnotPi=3.0;
   Double_t nSigmaTPCmax_PnotPi=3.0;
   Double_t nSigmaTPCmin_KnotPi=3.0;
   Double_t nSigmaTPCmax_KnotPi=3.0;
   Double_t nSigmaTOFmin_PnotPi=3.0;
   Double_t nSigmaTOFmax_PnotPi=3.0;
   Double_t nSigmaTOFmin_KnotPi=3.0;
   Double_t nSigmaTOFmax_KnotPi=3.0;

   Double_t etaRange=0.8;
   Double_t PtMin_P=0.15;
   Double_t PtMax_P=1.e10;
   Double_t PtMin_K=0.15;
   Double_t PtMax_K=1.e10;
   Double_t PMax_P=1.1;
   Double_t PMax_K=0.6;

   Bool_t useTPC_P=kFALSE;
   Bool_t asymmTPC_P=kFALSE;
   Bool_t useTPC_PnotPi=kFALSE;
   Bool_t useTOF_P=kFALSE;
   Bool_t asymmTOF_P=kFALSE;
   Bool_t useTOF_PnotPi=kFALSE;
   Bool_t rejectUnmatchedTOF_P=kTRUE;
   Bool_t useTPC_K=kFALSE;
   Bool_t asymmTPC_K=kFALSE;
   Bool_t useTPC_KnotPi=kFALSE;
   Bool_t useTOF_K=kFALSE;
   Bool_t asymmTOF_K=kFALSE;
   Bool_t useTOF_KnotPi=kFALSE;
   Bool_t rejectUnmatchedTOF_K=kTRUE;

   if (opt.Contains("qualityonly")) {
      useTPC_P=kFALSE;
      useTOF_P=kFALSE;
      useTPC_K=kFALSE;
      useTOF_K=kFALSE;
   } else if (!opt.Contains("nsig")) {
      useTPC_P=kTRUE;
      useTOF_P=kTRUE;
      useTPC_K=kTRUE;
      useTOF_K=kTRUE;
   }

   if (opt.Contains("PTPCnsig")){
     if(opt.Contains("PTPCnsigm") && opt.Contains("PTPCnsigp")){
       asymmTPC_P=kTRUE;
       if(opt.Contains("PTPCnsigm10")) nSigmaTPCmin_P=-1.;
       if(opt.Contains("PTPCnsigm15")) nSigmaTPCmin_P=-1.5;
       if(opt.Contains("PTPCnsigm20")) nSigmaTPCmin_P=-2.0;
       if(opt.Contains("PTPCnsigm25")) nSigmaTPCmin_P=-2.5;
       if(opt.Contains("PTPCnsigm30")) nSigmaTPCmin_P=-3.0;
       if(opt.Contains("PTPCnsigm40")) nSigmaTPCmin_P=-4.0;
       if(opt.Contains("PTPCnsigm50")) nSigmaTPCmin_P=-5.0;
       if(opt.Contains("PTPCnsigm1000")) nSigmaTPCmin_P=-100.0;

       if(opt.Contains("PTPCnsigp10")) nSigmaTPCmax_P=1.;
       if(opt.Contains("PTPCnsigp15")) nSigmaTPCmax_P=1.5;
       if(opt.Contains("PTPCnsigp20")) nSigmaTPCmax_P=2.0;
       if(opt.Contains("PTPCnsigp25")) nSigmaTPCmax_P=2.5;
       if(opt.Contains("PTPCnsigp30")) nSigmaTPCmax_P=3.0;
       if(opt.Contains("PTPCnsigp40")) nSigmaTPCmax_P=4.0;
       if(opt.Contains("PTPCnsigp50")) nSigmaTPCmax_P=5.0;
       if(opt.Contains("PTPCnsigp1000")) nSigmaTPCmax_P=100.0;
     }else if((opt.Contains("PTPCnsigm") && !opt.Contains("PTPCnsigp")) || (!opt.Contains("PTPCnsigm") && opt.Contains("PTPCnsigp"))){
       Printf("Unable to use asymmetric proton TPC cut.");
     }else{
       useTPC_P=kTRUE;
       if (opt.Contains("PTPCnsig10")) nSigmaTPC_P = 1.0;
       if (opt.Contains("PTPCnsig15")) nSigmaTPC_P = 1.5;
       if (opt.Contains("PTPCnsig20")) nSigmaTPC_P = 2.0;
       if (opt.Contains("PTPCnsig25")) nSigmaTPC_P = 2.5;
       if (opt.Contains("PTPCnsig30")) nSigmaTPC_P = 3.0;
       if (opt.Contains("PTPCnsig40")) nSigmaTPC_P = 4.0;
       if (opt.Contains("PTPCnsig50")) nSigmaTPC_P = 5.0;
       if (opt.Contains("PTPCnsig1000")) nSigmaTPC_P = 100.0;
     }
   }

   if(opt.Contains("PnotPiTPCnsig")){
     if(opt.Contains("PnotPiTPCnsigm") && opt.Contains("PnotPiTPCnsigp")){
       useTPC_PnotPi=kTRUE;
       if(opt.Contains("PnotPiTPCnsigm10")) nSigmaTPCmin_PnotPi=-1.;
       if(opt.Contains("PnotPiTPCnsigm15")) nSigmaTPCmin_PnotPi=-1.5;
       if(opt.Contains("PnotPiTPCnsigm20")) nSigmaTPCmin_PnotPi=-2.0;
       if(opt.Contains("PnotPiTPCnsigm25")) nSigmaTPCmin_PnotPi=-2.5;
       if(opt.Contains("PnotPiTPCnsigm30")) nSigmaTPCmin_PnotPi=-3.0;
       if(opt.Contains("PnotPiTPCnsigm40")) nSigmaTPCmin_PnotPi=-4.0;
       if(opt.Contains("PnotPiTPCnsigm50")) nSigmaTPCmin_PnotPi=-5.0;
       if(opt.Contains("PnotPiTPCnsigm1000")) nSigmaTPCmin_PnotPi=-100.0;

       if(opt.Contains("PnotPiTPCnsigp10")) nSigmaTPCmax_PnotPi=1.;
       if(opt.Contains("PnotPiTPCnsigp15")) nSigmaTPCmax_PnotPi=1.5;
       if(opt.Contains("PnotPiTPCnsigp20")) nSigmaTPCmax_PnotPi=2.0;
       if(opt.Contains("PnotPiTPCnsigp25")) nSigmaTPCmax_PnotPi=2.5;
       if(opt.Contains("PnotPiTPCnsigp30")) nSigmaTPCmax_PnotPi=3.0;
       if(opt.Contains("PnotPiTPCnsigp40")) nSigmaTPCmax_PnotPi=4.0;
       if(opt.Contains("PnotPiTPCnsigp50")) nSigmaTPCmax_PnotPi=5.0;
       if(opt.Contains("PnotPiTPCnsigp1000")) nSigmaTPCmax_PnotPi=100.0;
     }else if((opt.Contains("PnotPiTPCnsigm") && !opt.Contains("PnotPiTPCnsigp")) || (!opt.Contains("PnotPiTPCnsigm") && opt.Contains("PnotPiTPCnsigp"))){
       Printf("Unable to use TPC pion veto for protons.");
     }
   }

   if (opt.Contains("KTPCnsig")){
     if(opt.Contains("KTPCnsigm") && opt.Contains("KTPCnsigp")){
       asymmTPC_K=kTRUE;
       if(opt.Contains("KTPCnsigm10")) nSigmaTPCmin_K=-1.;
       if(opt.Contains("KTPCnsigm15")) nSigmaTPCmin_K=-1.5;
       if(opt.Contains("KTPCnsigm20")) nSigmaTPCmin_K=-2.0;
       if(opt.Contains("KTPCnsigm25")) nSigmaTPCmin_K=-2.5;
       if(opt.Contains("KTPCnsigm30")) nSigmaTPCmin_K=-3.0;
       if(opt.Contains("KTPCnsigm40")) nSigmaTPCmin_K=-4.0;
       if(opt.Contains("KTPCnsigm50")) nSigmaTPCmin_K=-5.0;
       if(opt.Contains("KTPCnsigm1000")) nSigmaTPCmin_K=-100.0;

       if(opt.Contains("KTPCnsigp10")) nSigmaTPCmax_K=1.;
       if(opt.Contains("KTPCnsigp15")) nSigmaTPCmax_K=1.5;
       if(opt.Contains("KTPCnsigp20")) nSigmaTPCmax_K=2.0;
       if(opt.Contains("KTPCnsigp25")) nSigmaTPCmax_K=2.5;
       if(opt.Contains("KTPCnsigp30")) nSigmaTPCmax_K=3.0;
       if(opt.Contains("KTPCnsigp40")) nSigmaTPCmax_K=4.0;
       if(opt.Contains("KTPCnsigp50")) nSigmaTPCmax_K=5.0;
       if(opt.Contains("KTPCnsigp1000")) nSigmaTPCmax_K=100.0;
     }else if((opt.Contains("KTPCnsigm") && !opt.Contains("KTPCnsigp")) || (!opt.Contains("KTPCnsigm") && opt.Contains("KTPCnsigp"))){
       Printf("Unable to use asymmetric kaon TPC cut.");
     }else{
       useTPC_K=kTRUE;
       if (opt.Contains("KTPCnsig10")) nSigmaTPC_K = 1.0;
       if (opt.Contains("KTPCnsig15")) nSigmaTPC_K = 1.5;
       if (opt.Contains("KTPCnsig20")) nSigmaTPC_K = 2.0;
       if (opt.Contains("KTPCnsig25")) nSigmaTPC_K = 2.5;
       if (opt.Contains("KTPCnsig30")) nSigmaTPC_K = 3.0;
       if (opt.Contains("KTPCnsig40")) nSigmaTPC_K = 4.0;
       if (opt.Contains("KTPCnsig50")) nSigmaTPC_K = 5.0;
       if (opt.Contains("KTPCnsig1000")) nSigmaTPC_K = 100.0;
     }
   }

   if(opt.Contains("KnotPiTPCnsig")){
     if(opt.Contains("KnotPiTPCnsigm") && opt.Contains("KnotPiTPCnsigp")){
       useTPC_KnotPi=kTRUE;
       if(opt.Contains("KnotPiTPCnsigm10")) nSigmaTPCmin_KnotPi=-1.;
       if(opt.Contains("KnotPiTPCnsigm15")) nSigmaTPCmin_KnotPi=-1.5;
       if(opt.Contains("KnotPiTPCnsigm20")) nSigmaTPCmin_KnotPi=-2.0;
       if(opt.Contains("KnotPiTPCnsigm25")) nSigmaTPCmin_KnotPi=-2.5;
       if(opt.Contains("KnotPiTPCnsigm30")) nSigmaTPCmin_KnotPi=-3.0;
       if(opt.Contains("KnotPiTPCnsigm40")) nSigmaTPCmin_KnotPi=-4.0;
       if(opt.Contains("KnotPiTPCnsigm50")) nSigmaTPCmin_KnotPi=-5.0;
       if(opt.Contains("KnotPiTPCnsigm1000")) nSigmaTPCmin_KnotPi=-100.0;

       if(opt.Contains("KnotPiTPCnsigp10")) nSigmaTPCmax_KnotPi=1.;
       if(opt.Contains("KnotPiTPCnsigp15")) nSigmaTPCmax_KnotPi=1.5;
       if(opt.Contains("KnotPiTPCnsigp20")) nSigmaTPCmax_KnotPi=2.0;
       if(opt.Contains("KnotPiTPCnsigp25")) nSigmaTPCmax_KnotPi=2.5;
       if(opt.Contains("KnotPiTPCnsigp30")) nSigmaTPCmax_KnotPi=3.0;
       if(opt.Contains("KnotPiTPCnsigp40")) nSigmaTPCmax_KnotPi=4.0;
       if(opt.Contains("KnotPiTPCnsigp50")) nSigmaTPCmax_KnotPi=5.0;
       if(opt.Contains("KnotPiTPCnsigp1000")) nSigmaTPCmax_KnotPi=100.0;
     }else if((opt.Contains("KnotPiTPCnsigm") && !opt.Contains("KnotPiTPCnsigp")) || (!opt.Contains("KnotPiTPCnsigm") && opt.Contains("KnotPiTPCnsigp"))){
       Printf("Unable to use TPC pion veto for kaons.");
     }
   }

   if (opt.Contains("PTOFacceptUnmatched")) rejectUnmatchedTOF_P=kFALSE;
   if (opt.Contains("PTOFnsig")){
     if(opt.Contains("PTOFnsigm") && opt.Contains("PTOFnsigp")){
       asymmTOF_P=kTRUE;
       if(opt.Contains("PTOFnsigm10")) nSigmaTOFmin_P=-1.0;
       if(opt.Contains("PTOFnsigm15")) nSigmaTOFmin_P=-1.5;
       if(opt.Contains("PTOFnsigm20")) nSigmaTOFmin_P=-2.0;
       if(opt.Contains("PTOFnsigm25")) nSigmaTOFmin_P=-2.5;
       if(opt.Contains("PTOFnsigm30")) nSigmaTOFmin_P=-3.0;
       if(opt.Contains("PTOFnsigm40")) nSigmaTOFmin_P=-4.0;
       if(opt.Contains("PTOFnsigm50")) nSigmaTOFmin_P=-5.0;
       if(opt.Contains("PTOFnsigm1000")) nSigmaTOFmin_P=-100.0;

       if(opt.Contains("PTOFnsigp10")) nSigmaTOFmax_P=1.0;
       if(opt.Contains("PTOFnsigp15")) nSigmaTOFmax_P=1.5;
       if(opt.Contains("PTOFnsigp20")) nSigmaTOFmax_P=2.0;
       if(opt.Contains("PTOFnsigp25")) nSigmaTOFmax_P=2.5;
       if(opt.Contains("PTOFnsigp30")) nSigmaTOFmax_P=3.0;
       if(opt.Contains("PTOFnsigp40")) nSigmaTOFmax_P=4.0;
       if(opt.Contains("PTOFnsigp50")) nSigmaTOFmax_P=5.0;
       if(opt.Contains("PTOFnsigp1000")) nSigmaTOFmax_P=100.0;
     }else if((opt.Contains("PTOFnsigm") && !opt.Contains("PTOFnsigp")) || (!opt.Contains("PTOFnsigm") && opt.Contains("PTOFnsigp"))){
       Printf("Unable to use asymmetric proton TOF cut.");
     }else{
       useTOF_P=kTRUE;
       if (opt.Contains("PTOFnsig10")) nSigmaTOF_P = 1.0;
       if (opt.Contains("PTOFnsig15")) nSigmaTOF_P = 1.5;
       if (opt.Contains("PTOFnsig20")) nSigmaTOF_P = 2.0;
       if (opt.Contains("PTOFnsig25")) nSigmaTOF_P = 2.5;
       if (opt.Contains("PTOFnsig30")) nSigmaTOF_P = 3.0;
       if (opt.Contains("PTOFnsig40")) nSigmaTOF_P = 4.0;
       if (opt.Contains("PTOFnsig50")) nSigmaTOF_P = 5.0;
       if (opt.Contains("PTOFnsig1000")) nSigmaTOF_P = 100.0;
     }
   }

   if(opt.Contains("PnotPiTOFnsig")){
     if(opt.Contains("PnotPiTOFnsigm") && opt.Contains("PnotPiTOFnsigp")){
       useTOF_PnotPi=kTRUE;
       if(opt.Contains("PnotPiTOFnsigm10")) nSigmaTOFmin_PnotPi=-1.;
       if(opt.Contains("PnotPiTOFnsigm15")) nSigmaTOFmin_PnotPi=-1.5;
       if(opt.Contains("PnotPiTOFnsigm20")) nSigmaTOFmin_PnotPi=-2.0;
       if(opt.Contains("PnotPiTOFnsigm25")) nSigmaTOFmin_PnotPi=-2.5;
       if(opt.Contains("PnotPiTOFnsigm30")) nSigmaTOFmin_PnotPi=-3.0;
       if(opt.Contains("PnotPiTOFnsigm40")) nSigmaTOFmin_PnotPi=-4.0;
       if(opt.Contains("PnotPiTOFnsigm50")) nSigmaTOFmin_PnotPi=-5.0;
       if(opt.Contains("PnotPiTOFnsigm1000")) nSigmaTOFmin_PnotPi=-100.0;

       if(opt.Contains("PnotPiTOFnsigp10")) nSigmaTOFmax_PnotPi=1.;
       if(opt.Contains("PnotPiTOFnsigp15")) nSigmaTOFmax_PnotPi=1.5;
       if(opt.Contains("PnotPiTOFnsigp20")) nSigmaTOFmax_PnotPi=2.0;
       if(opt.Contains("PnotPiTOFnsigp25")) nSigmaTOFmax_PnotPi=2.5;
       if(opt.Contains("PnotPiTOFnsigp30")) nSigmaTOFmax_PnotPi=3.0;
       if(opt.Contains("PnotPiTOFnsigp40")) nSigmaTOFmax_PnotPi=4.0;
       if(opt.Contains("PnotPiTOFnsigp50")) nSigmaTOFmax_PnotPi=5.0;
       if(opt.Contains("PnotPiTOFnsigp1000")) nSigmaTOFmax_PnotPi=100.0;
     }else if((opt.Contains("PnotPiTOFnsigm") && !opt.Contains("PnotPiTOFnsigp")) || (!opt.Contains("PnotPiTOFnsigm") && opt.Contains("PnotPiTOFnsigp"))){
       Printf("Unable to use TOF pion veto for protons.");
     }
   }

   if (opt.Contains("KTOFacceptUnmatched")) rejectUnmatchedTOF_K=kFALSE;
   if (opt.Contains("KTOFnsig")){
     if(opt.Contains("KTOFnsigm") && opt.Contains("KTOFnsigp")){
       asymmTOF_K=kTRUE;
       if(opt.Contains("KTOFnsigm10")) nSigmaTOFmin_K=-1.0;
       if(opt.Contains("KTOFnsigm15")) nSigmaTOFmin_K=-1.5;
       if(opt.Contains("KTOFnsigm20")) nSigmaTOFmin_K=-2.0;
       if(opt.Contains("KTOFnsigm25")) nSigmaTOFmin_K=-2.5;
       if(opt.Contains("KTOFnsigm30")) nSigmaTOFmin_K=-3.0;
       if(opt.Contains("KTOFnsigm40")) nSigmaTOFmin_K=-4.0;
       if(opt.Contains("KTOFnsigm50")) nSigmaTOFmin_K=-5.0;
       if(opt.Contains("KTOFnsigm1000")) nSigmaTOFmin_K=-100.0;

       if(opt.Contains("KTOFnsigp10")) nSigmaTOFmax_K=1.0;
       if(opt.Contains("KTOFnsigp15")) nSigmaTOFmax_K=1.5;
       if(opt.Contains("KTOFnsigp20")) nSigmaTOFmax_K=2.0;
       if(opt.Contains("KTOFnsigp25")) nSigmaTOFmax_K=2.5;
       if(opt.Contains("KTOFnsigp30")) nSigmaTOFmax_K=3.0;
       if(opt.Contains("KTOFnsigp40")) nSigmaTOFmax_K=4.0;
       if(opt.Contains("KTOFnsigp50")) nSigmaTOFmax_K=5.0;
       if(opt.Contains("KTOFnsigp1000")) nSigmaTOFmax_K=100.0;
     }else if((opt.Contains("KTOFnsigm") && !opt.Contains("KTOFnsigp")) || (!opt.Contains("KTOFnsigm") && opt.Contains("KTOFnsigp"))){
       Printf("Unable to use asymmetric kaon TOF cut.");
     }else{
       useTOF_K=kTRUE;
       if (opt.Contains("KTOFnsig10")) nSigmaTOF_K = 1.0;
       if (opt.Contains("KTOFnsig15")) nSigmaTOF_K = 1.5;
       if (opt.Contains("KTOFnsig20")) nSigmaTOF_K = 2.0;
       if (opt.Contains("KTOFnsig25")) nSigmaTOF_K = 2.5;
       if (opt.Contains("KTOFnsig30")) nSigmaTOF_K = 3.0;
       if (opt.Contains("KTOFnsig40")) nSigmaTOF_K = 4.0;
       if (opt.Contains("KTOFnsig50")) nSigmaTOF_K = 5.0;
       if (opt.Contains("KTOFnsig1000")) nSigmaTOF_K = 100.0;
     }
   }

   if(opt.Contains("KnotPiTOFnsig")){
     if(opt.Contains("KnotPiTOFnsigm") && opt.Contains("KnotPiTOFnsigp")){
       useTOF_KnotPi=kTRUE;
       if(opt.Contains("KnotPiTOFnsigm10")) nSigmaTOFmin_KnotPi=-1.;
       if(opt.Contains("KnotPiTOFnsigm15")) nSigmaTOFmin_KnotPi=-1.5;
       if(opt.Contains("KnotPiTOFnsigm20")) nSigmaTOFmin_KnotPi=-2.0;
       if(opt.Contains("KnotPiTOFnsigm25")) nSigmaTOFmin_KnotPi=-2.5;
       if(opt.Contains("KnotPiTOFnsigm30")) nSigmaTOFmin_KnotPi=-3.0;
       if(opt.Contains("KnotPiTOFnsigm40")) nSigmaTOFmin_KnotPi=-4.0;
       if(opt.Contains("KnotPiTOFnsigm50")) nSigmaTOFmin_KnotPi=-5.0;
       if(opt.Contains("KnotPiTOFnsigm1000")) nSigmaTOFmin_KnotPi=-100.0;

       if(opt.Contains("KnotPiTOFnsigp10")) nSigmaTOFmax_KnotPi=1.;
       if(opt.Contains("KnotPiTOFnsigp15")) nSigmaTOFmax_KnotPi=1.5;
       if(opt.Contains("KnotPiTOFnsigp20")) nSigmaTOFmax_KnotPi=2.0;
       if(opt.Contains("KnotPiTOFnsigp25")) nSigmaTOFmax_KnotPi=2.5;
       if(opt.Contains("KnotPiTOFnsigp30")) nSigmaTOFmax_KnotPi=3.0;
       if(opt.Contains("KnotPiTOFnsigp40")) nSigmaTOFmax_KnotPi=4.0;
       if(opt.Contains("KnotPiTOFnsigp50")) nSigmaTOFmax_KnotPi=5.0;
       if(opt.Contains("KnotPiTOFnsigp1000")) nSigmaTOFmax_KnotPi=100.0;
     }else if((opt.Contains("KnotPiTOFnsigm") && !opt.Contains("KnotPiTOFnsigp")) || (!opt.Contains("KnotPiTOFnsigm") && opt.Contains("KnotPiTOFnsigp"))){
       Printf("Unable to use TOF pion veto for kaons.");
     }
   }


   Bool_t usePDG=kFALSE;
   if (opt.Contains("pdg")) {
      Printf("Using PDG");
      usePDG = kTRUE;
   }

   Bool_t useEta = kFALSE;
   if (opt.Contains("eta")) {
     for(int j=1;j<=9;j++) if(opt.Contains(Form("eta0%i",j))) etaRange=0.1*j;

      Printf("Using ETA range (%.2f,%.2f)",-etaRange,etaRange);
      useEta = kTRUE;
   }

   Bool_t usePPt=kFALSE;
   if(opt.Contains("PPt")){
     Printf("Using Proton pT range (%.2f,%.2f)",PtMin_P,PtMax_P);
     usePPt=kTRUE;
   }

   Bool_t useKPt=kFALSE;
   if(opt.Contains("KPt")){
     Printf("Using Kaon pT range (%.2f,%.2f)",PtMin_K,PtMax_K);
     useKPt=kTRUE;
   }

   Bool_t usePMax_P=kFALSE;
   if(opt.Contains("PPMax")){
     for(int j=1;j<=9;j++) if(opt.Contains(Form("PPMax0%i",j))) PMax_P=0.1*j;
     for(int j=10;j<=30;j++) if(opt.Contains(Form("PPMax%i",j))) PMax_P=0.1*j;
     Printf("Using Proton momentum range (0,%.2f)",PMax_P);
     usePMax_P=kTRUE;
   }

   Bool_t usePMax_K=kFALSE;
   if(opt.Contains("KPMax")){
     for(int j=1;j<=9;j++) if(opt.Contains(Form("KPMax0%i",j))) PMax_K=0.1*j;
     for(int j=10;j<=30;j++) if(opt.Contains(Form("KPMax%i",j))) PMax_K=0.1*j;
     Printf("Using Kaon momentum range (0,%.2f)",PMax_K);
     usePMax_K=kTRUE;
   }

   // PROTON SETTINGS ===========================================

   TString scheme="";
   TString cutname = "p_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsP = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutP = new AliRsnCutTrackQuality("cutQualityP");
   if (useCommonQualityCut>=0) {
      qualityCutP->SetAODTestFilterBit(useCommonQualityCut);
   } else {
      qualityCutP->SetDefaults2010();
   }
   cutsP->AddCut(qualityCutP);
   if (!scheme.IsNull()) scheme += "&";
   scheme += qualityCutP->GetName();

   if (useTPC_P) {
      AliRsnCutPIDNSigma *cutPTPC = new AliRsnCutPIDNSigma("cutNSigmaTPCP",AliPID::kProton,AliRsnCutPIDNSigma::kTPC);
      cutPTPC->SinglePIDRange(nSigmaTPC_P);
      cutsP->AddCut(cutPTPC);
      if (!scheme.IsNull()) scheme += "&";
      scheme += cutPTPC->GetName();
   }

   if(asymmTPC_P){
     AliRsnCutPIDTPC* cutPTPCasymm=new AliRsnCutPIDTPC("cutNsigmaTPCP",AliPID::kProton,nSigmaTPCmin_P,nSigmaTPCmax_P);
     cutsP->AddCut(cutPTPCasymm);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPTPCasymm->GetName();
   }

   if(useTPC_PnotPi){
     AliRsnCutPIDTPC* cutPnotPiTPC=new AliRsnCutPIDTPC("cutNsigmaTPCPnotPi",AliPID::kPion,nSigmaTPCmin_PnotPi,nSigmaTPCmax_PnotPi);
     cutsP->AddCut(cutPnotPiTPC);
     if (!scheme.IsNull()) scheme += "&";
     scheme += Form("(!%s)",cutPnotPiTPC->GetName());
   }

   if (useTOF_P) {
      AliRsnCutPIDNSigma *cutPTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFP",AliPID::kProton,AliRsnCutPIDNSigma::kTOF);
      cutPTOF->SinglePIDRange(nSigmaTOF_P);
      cutsP->AddCut(cutPTOF);
      if(rejectUnmatchedTOF_P){
	if (!scheme.IsNull()) scheme += "&";
	scheme += cutPTOF->GetName();
      }else{
	AliRsnCutTOFMatch *cutPTOFMatch = new AliRsnCutTOFMatch("cutPTOFMatch");
	cutsP->AddCut(cutPTOFMatch);
	if (!scheme.IsNull()) scheme += "&";
	scheme += Form("(%s|(!%s))",cutPTOF->GetName(),cutPTOFMatch->GetName());
      }
   }

   AliRsnCutTOFMatch* cutPTOFMatch=0;
   if(asymmTOF_P){
     AliRsnCutPIDTOF* cutPTOFasymm=new AliRsnCutPIDTOF("cutNsigmaTOFP",AliPID::kProton,nSigmaTOFmin_P,nSigmaTOFmax_P);
     cutsP->AddCut(cutPTOFasymm);
     if(rejectUnmatchedTOF_P){
       if (!scheme.IsNull()) scheme += "&";
       scheme += cutPTOFasymm->GetName();
     }else{
       if(!cutPTOFMatch){
	 cutPTOFMatch = new AliRsnCutTOFMatch("cutPTOFMatch");
	 cutsP->AddCut(cutPTOFMatch);
       }
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(%s|(!%s))",cutPTOFasymm->GetName(),cutPTOFMatch->GetName());
     }
   }

   if(useTOF_PnotPi){
     AliRsnCutPIDTOF* cutPnotPiTOF=new AliRsnCutPIDTOF("cutNsigmaTOFPnotPi",AliPID::kPion,nSigmaTOFmin_PnotPi,nSigmaTOFmax_PnotPi);
     cutsP->AddCut(cutPnotPiTOF);
     if(rejectUnmatchedTOF_P){
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(!%s)",cutPnotPiTOF->GetName());
     }else{
       if(!cutPTOFMatch){
	 cutPTOFMatch = new AliRsnCutTOFMatch("cutPTOFMatch");
	 cutsP->AddCut(cutPTOFMatch);
       }
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(!(%s&%s))",cutPnotPiTOF->GetName(),cutPTOFMatch->GetName());
     }
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

   if(usePPt){
     AliRsnValueDaughter *valPtP = new AliRsnValueDaughter(Form("val%sPt%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kPt);
     AliRsnCutValue *cutPtP = new AliRsnCutValue(Form("cut%sPt%s",AliPID::ParticleName(type1),opt.Data()),PtMin_P,PtMax_P);
     cutPtP->SetTargetType(AliRsnTarget::kDaughter);
     cutPtP->SetValueObj(valPtP);
     cutsP->AddCut(cutPtP);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPtP->GetName();
   }

   if(usePMax_P){
     AliRsnValueDaughter *valPP = new AliRsnValueDaughter(Form("val%sP%s",AliPID::ParticleName(type1),opt.Data()),AliRsnValueDaughter::kP);
     AliRsnCutValue *cutPP = new AliRsnCutValue(Form("cut%sP%s",AliPID::ParticleName(type1),opt.Data()),0.,PMax_P);
     cutPP->SetTargetType(AliRsnTarget::kDaughter);
     cutPP->SetValueObj(valPP);
     cutsP->AddCut(cutPP);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPP->GetName();
   }

   Printf ("CUT Scheme for PROTON is '%s'",scheme.Data());
   cutsP->SetCutScheme(scheme.Data());

   // END PROTON =======================================

   // KAON SETTINGS =======================================
   scheme="";
   cutname = "K_Lambda";
   if (!opt.IsNull()) cutname += Form("_%s",opt.Data());
   AliRsnCutSet *cutsK = new AliRsnCutSet(cutname.Data(), AliRsnTarget::kDaughter);

   AliRsnCutTrackQuality *qualityCutK = new AliRsnCutTrackQuality("cutQualityK");
   if (useCommonQualityCut>=0) {
      qualityCutK->SetAODTestFilterBit(useCommonQualityCut);
   } else {
      qualityCutK->SetDefaults2010();
   }
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

   if(asymmTPC_K){
     AliRsnCutPIDTPC* cutKTPCasymm=new AliRsnCutPIDTPC("cutNsigmaTPCP",AliPID::kKaon,nSigmaTPCmin_K,nSigmaTPCmax_K);
     cutsK->AddCut(cutKTPCasymm);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutKTPCasymm->GetName();
   }

   if(useTPC_KnotPi){
     AliRsnCutPIDTPC* cutKnotPiTPC=new AliRsnCutPIDTPC("cutNsigmaTPCKnotPi",AliPID::kPion,nSigmaTPCmin_KnotPi,nSigmaTPCmax_KnotPi);
     cutsP->AddCut(cutKnotPiTPC);
     if (!scheme.IsNull()) scheme += "&";
     scheme += Form("(!%s)",cutKnotPiTPC->GetName());
   }

   if (useTOF_K) {
      AliRsnCutPIDNSigma *cutKTOF = new AliRsnCutPIDNSigma("cutNSigmaTOFK",AliPID::kKaon,AliRsnCutPIDNSigma::kTOF);
      cutKTOF->SinglePIDRange(nSigmaTOF_K);
      cutsK->AddCut(cutKTOF);
      if(rejectUnmatchedTOF_K){
	if (!scheme.IsNull()) scheme += "&";
	scheme += cutKTOF->GetName();
      }else{
	AliRsnCutTOFMatch *cutKTOFMatch = new AliRsnCutTOFMatch("cutKTOFMatch");
	cutsK->AddCut(cutKTOFMatch);
	if (!scheme.IsNull()) scheme += "&";
	scheme += Form("(%s|(!%s))",cutKTOF->GetName(),cutKTOFMatch->GetName());
      }
   }

   AliRsnCutTOFMatch* cutKTOFMatch=0;
   if(asymmTOF_K){
     AliRsnCutPIDTOF* cutKTOFasymm=new AliRsnCutPIDTOF("cutNsigmaTOFK",AliPID::kKaon,nSigmaTOFmin_K,nSigmaTOFmax_K);
     cutsK->AddCut(cutKTOFasymm);
     if(rejectUnmatchedTOF_K){
       if (!scheme.IsNull()) scheme += "&";
       scheme += cutKTOFasymm->GetName();
     }else{
       if(!cutKTOFMatch){
	 cutKTOFMatch = new AliRsnCutTOFMatch("cutKTOFMatch");
	 cutsK->AddCut(cutKTOFMatch);
       }
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(%s|(!%s))",cutKTOFasymm->GetName(),cutKTOFMatch->GetName());
     }
   }

   if(useTOF_KnotPi){
     AliRsnCutPIDTOF* cutKnotPiTOF=new AliRsnCutPIDTOF("cutNsigmaTOFKnotPi",AliPID::kPion,nSigmaTOFmin_KnotPi,nSigmaTOFmax_KnotPi);
     cutsP->AddCut(cutKnotPiTOF);
     if(rejectUnmatchedTOF_K){
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(!%s)",cutKnotPiTOF->GetName());
     }else{
       if(!cutKTOFMatch){
	 cutKTOFMatch = new AliRsnCutTOFMatch("cutKTOFMatch");
	 cutsK->AddCut(cutKTOFMatch);
       }
       if (!scheme.IsNull()) scheme += "&";
       scheme += Form("(!(%s&%s))",cutKnotPiTOF->GetName(),cutKTOFMatch->GetName());
     }
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

   if(useKPt){
     AliRsnValueDaughter *valPtK = new AliRsnValueDaughter(Form("val%sPt%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kPt);
     AliRsnCutValue *cutPtK = new AliRsnCutValue(Form("cut%sPt%s",AliPID::ParticleName(type2),opt.Data()),PtMin_K,PtMax_K);
     cutPtK->SetTargetType(AliRsnTarget::kDaughter);
     cutPtK->SetValueObj(valPtK);
     cutsK->AddCut(cutPtK);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPtK->GetName();
   }

   if(usePMax_K){
     AliRsnValueDaughter *valPK = new AliRsnValueDaughter(Form("val%sP%s",AliPID::ParticleName(type2),opt.Data()),AliRsnValueDaughter::kP);
     AliRsnCutValue *cutPK = new AliRsnCutValue(Form("cut%sP%s",AliPID::ParticleName(type2),opt.Data()),0.,PMax_K);
     cutPK->SetTargetType(AliRsnTarget::kDaughter);
     cutPK->SetValueObj(valPK);
     cutsK->AddCut(cutPK);
     if (!scheme.IsNull()) scheme += "&";
     scheme += cutPK->GetName();
   }

   Printf ("CUT Scheme for KAON is '%s'",scheme.Data());
   cutsK->SetCutScheme(scheme.Data());

   // END KAON =======================================

   if (opt.Contains("mon")) {
      AddMonitorOutput(cutsP->GetMonitorOutput(),opt);
      AddMonitorOutput(cutsK->GetMonitorOutput(),opt);
   }
   if (isRsnMini) {
      AliRsnMiniAnalysisTask *taskRsnMini = dynamic_cast<AliRsnMiniAnalysisTask *>(task);
      if (taskRsnMini) {
         taskRsnMini->AddTrackCuts(cutsP);
         taskRsnMini->AddTrackCuts(cutsK);
      }
   } else {
      AliRsnDaughterSelector *sel = rsnIH->GetSelector();
//       sel->SetLabelCheck(kFALSE);
      sel->Add(cutsP, kTRUE);
      sel->Add(cutsK, kTRUE);
   }
   return numberOfCuts;

}
