void AddRsnPairsKStar(AliAnalysisTaskSE *task,
                      Bool_t isMC,
                      Bool_t isMixing,
                      AliPID::EParticleType pType1,
                      Int_t listID1,
                      AliPID::EParticleType pType2,
                      Int_t listID2,
                      AliRsnCutSet *cutsEvent=0,
                      AliRsnCutSet *cutsPair=0,
                      TString suffix = "") {

   Printf("id1=%d id2=%d",listID1,listID2);

   // retrieve mass from PDG database
   Int_t         pdg  = 313;
   TDatabasePDG *db   = TDatabasePDG::Instance();
   TParticlePDG *part = db->GetParticle(pdg);
   Double_t mass = part->Mass();

   Bool_t valid;
   Int_t isRsnMini = AliAnalysisManager::GetGlobalInt("rsnUseMiniPackage",valid);
   if (isRsnMini) {
      AddPairOutputMiniKStar(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsPair,suffix);
   } else {
      // this function is common and it is located in RsnConfig.C
      // as ouptup AddPairOutputKStar from this macro will be taken
      AddPair(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsEvent,cutsPair,suffix);
   }
}
void AddPairOutputKStar(AliRsnLoopPair *pair)
{
   Bool_t valid;
   Int_t isFullOutput = AliAnalysisManager::GetGlobalInt("rsnOutputFull",valid);
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);

   // axes
   AliRsnValuePair *axisIM = new AliRsnValuePair("IM", AliRsnValuePair::kInvMass);
   axisIM     ->SetBins(900, 0.6, 1.5);

   AliRsnValuePair *axisPt = new AliRsnValuePair("PT", AliRsnValuePair::kPt);
   axisPt     ->SetBins(120, 0.0, 12.0);

   AliRsnValuePair *axisEta = new AliRsnValuePair("ETA", AliRsnValuePair::kEta);
   axisEta    ->SetBins(400, -0.5, 0.5);

   AliRsnValueEvent *axisCentrality = 0;
   if (!isPP) axisCentrality = new AliRsnValueEvent("MULTI",AliRsnValueEvent::kCentralityV0);

   // output: 2D histogram of inv. mass vs. pt
   AliRsnListOutput *outPair = 0;
   if (!isFullOutput) {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoDefault);
      outPair->AddValue(axisIM);
   } else {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoSparse);
      outPair->AddValue(axisIM);
      outPair->AddValue(axisPt);
      outPair->AddValue(axisEta);
      if (axisCentrality) axisCentrality->SetBins(20,0,100);
   }
   // add outputs to loop
   pair->AddOutput(outPair);
}

void AddPairOutputMiniKStar(AliAnalysisTaskSE *task,Bool_t isMC,Bool_t isMixing, AliPID::EParticleType pType1,Int_t listID1, AliPID::EParticleType pType2,Int_t listID2, Int_t pdgMother,Double_t massMother, AliRsnCutSet *cutsPair=0,TString suffix = "") {

   Bool_t valid;
   Int_t isFullOutput = AliAnalysisManager::GetGlobalInt("rsnOutputFull",valid);
   Int_t useMixing = AliAnalysisManager::GetGlobalInt("rsnUseMixing",valid);
   Int_t isPP = AliAnalysisManager::GetGlobalInt("rsnIsPP",valid);

   AliRsnMiniAnalysisTask *taskRsnMini =  (AliRsnMiniAnalysisTask *)task;

   if (isPP) taskRsnMini->UseMultiplicity("QUALITY");
   else {
      taskRsnMini->UseCentrality("V0M");
      Int_t multID = taskRsnMini->CreateValue(AliRsnMiniValue::kMult, kFALSE);
      AliRsnMiniOutput *outMult = taskRsnMini->CreateOutput("eventMult", "HIST", "EVENT");
      outMult->AddAxis(multID, 100, 0.0, 100.0);
      Int_t paID = taskRsnMini->CreateValue(AliRsnMiniValue::kPlaneAngle, kFALSE);
      AliRsnMiniOutput *outPa = taskRsnMini->CreateOutput("planeAngle", "HIST", "EVENT");
      outPa->AddAxis(paID, 100, 0, TMath::Pi());
   }

   /* invariant mass   */ Int_t imID   = taskRsnMini->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = taskRsnMini->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* transv. momentum */ Int_t ptID   = taskRsnMini->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = taskRsnMini->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* eta              */ Int_t etaID = taskRsnMini->CreateValue(AliRsnMiniValue::kEta, kFALSE);
   /* rapidity         */ Int_t yID = taskRsnMini->CreateValue(AliRsnMiniValue::kY, kFALSE);

   Bool_t useRapidity = kTRUE;

   Int_t nIM   = 90; Double_t minIM   = 0.6, maxIM =  1.5;
   Int_t nRes   = 200; Double_t minRes   = -0.02, maxRes =  0.02;
   Int_t nEta   = 400; Double_t minEta   = -0.5, maxEta =  0.5;
   Int_t nY   = 16; Double_t minY   = -0.8, maxY =  0.8;
   Int_t nPt   = 120; Double_t minPt   = 0.0, maxPt = 12.0;
   Int_t nCent = 100; Double_t minCent = 0.0, maxCent = 100.0;
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //

   Int_t iCutK = listID1;
   Int_t iCutPi = listID2;

   // common definitions
   TString outputType = "HIST";
   if (isFullOutput) outputType = "SPARSE";

   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
//    const Int_t num = 12;
   Bool_t  use     [12] = { 1       ,  1       , useMixing, useMixing,  1      ,  1      ,  isMC   ,   isMC   ,  isMC   ,   isMC   ,  isMC    ,   isMC     };
   Bool_t  useIM   [12] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      ,  1       ,   1        };
   TString name    [12] = {"Unlike1", "Unlike2", "Mixing1", "Mixing2", "LikePP", "LikeMM", "Trues1",  "Trues2", "Res1"  ,  "Res2"  , "Mother1",  "Mother2" };
   TString comp    [12] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  , "MOTHER" ,  "MOTHER"  };
   Char_t  charge1 [12] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     , '+'      ,  '-'       };
   Char_t  charge2 [12] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     , '-'      ,  '+'       };
   Int_t   cutID1  [12] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK     };
   Int_t   cutID2  [12] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi    };

   for (Int_t i = 0; i < 12; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = taskRsnMini->CreateOutput(Form("%s_%s", suffix.Data(), name[i].Data()), outputType.Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kPion);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(pdgMother);
      out->SetMotherMass(massMother);
      // pair cuts
      if (cutsPair) out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i])
         out->AddAxis(imID, nIM, minIM, maxIM);
      else
         out->AddAxis(resID, nRes, minRes, maxRes);

      if (isFullOutput) {
         // axis Y: transverse momentum
         out->AddAxis(ptID, nPt, minPt, maxPt);
         if (useRapidity) out->AddAxis(yID, nY, minY, maxY);
         else  out->AddAxis(etaID, nEta, minEta, maxEta);
         // axis Z: centrality
         if (!isPP) out->AddAxis(centID, nCent, minCent, maxCent);
      }
   }

   // -- Create output for MC generated ------------------------------------------------------------
   //

   if (isMC) {
      // create ouput
      AliRsnMiniOutput *outMC = taskRsnMini->CreateOutput(Form("kstar_MCGen1%s", suffix.Data()), outputType.Data(), "MOTHER");
      // selection settings
      outMC->SetDaughter(0, AliRsnDaughter::kPion);
      outMC->SetDaughter(1, AliRsnDaughter::kKaon);
      outMC->SetMotherPDG(313);
      outMC->SetMotherMass(massMother);
      // pair cuts
      if (cutsPair) outMC->SetPairCuts(cutsPair);
      // axis X: invmass
      outMC->AddAxis(imID, nIM, minIM, maxIM);
      if (isFullOutput) {
         // axis Y: transverse momentum
         outMC->AddAxis(ptID, nPt, minPt, maxPt);
         if (useRapidity) outMC->AddAxis(yID, nY, minY, maxY);
         else  outMC->AddAxis(etaID, nEta, minEta, maxEta);
         // axis Z: centrality
         if (!isPP) outMC->AddAxis(centID, nCent, minCent, maxCent);
      }
   }



   if (isMC) {
      // create ouput
      AliRsnMiniOutput *outMC1 = taskRsnMini->CreateOutput(Form("phi_MCGen2%s", suffix.Data()), outputType.Data(), "MOTHER");
      // selection settings
      outMC1->SetDaughter(0, AliRsnDaughter::kKaon);
      outMC1->SetDaughter(1, AliRsnDaughter::kPion);
      outMC1->SetMotherPDG(-313);
      outMC1->SetMotherMass(massMother);
      // pair cuts
      if (cutsPair) outMC1->SetPairCuts(cutsPair);
      // axis X: invmass
      outMC1->AddAxis(imID, nIM, minIM, maxIM);
      if (isFullOutput) {
         // axis Y: transverse momentum
         outMC1->AddAxis(ptID, nPt, minPt, maxPt);
         if (useRapidity) outMC1->AddAxis(yID, nY, minY, maxY);
         else  outMC1->AddAxis(etaID, nEta, minEta, maxEta);
         // axis Z: centrality
         if (!isPP) outMC1->AddAxis(centID, nCent, minCent, maxCent);
      }
   }


}

