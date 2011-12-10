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

   if (gRsnUseMiniPackage) {
      AddPairOutputMiniKStar(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsPair,suffix);
   } else {
      // this function is common and it is located in RsnConfig.C
      // as ouptup AddPairOutputKStar from this macro will be taken
      AddPair(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsEvent,cutsPair,suffix);
   }
}
void AddPairOutputKStar(AliRsnLoopPair *pair)
{

   // axes
   AliRsnValuePair *axisIM = new AliRsnValuePair("IM", AliRsnValuePair::kInvMass);
   AliRsnValuePair *axisPt = new AliRsnValuePair("PT", AliRsnValuePair::kPt);
   axisIM     ->SetBins(900, 0.6, 1.5);
   axisPt     ->SetBins(120, 0.0, 12.0);

   // output: 2D histogram of inv. mass vs. pt
   AliRsnListOutput *outPair = 0;
   if (!gRsnOutputFull) {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoDefault);
      outPair->AddValue(axisIM);
   } else {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoSparse);
      outPair->AddValue(axisIM);
      outPair->AddValue(axisPt);
   }
   // add outputs to loop
   pair->AddOutput(outPair);
}

void AddPairOutputMiniKStar(AliAnalysisTaskSE *task,Bool_t isMC,Bool_t isMixing, AliPID::EParticleType pType1,Int_t listID1, AliPID::EParticleType pType2,Int_t listID2, Int_t pdgMother,Double_t massMother, AliRsnCutSet *cutsPair=0,TString suffix = "") {

   AliRsnMiniAnalysisTask *taskRsnMini =  (AliRsnMiniAnalysisTask *)task;
   /* invariant mass   */ Int_t imID   = taskRsnMini->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = taskRsnMini->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* transv. momentum */ Int_t ptID   = taskRsnMini->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = taskRsnMini->CreateValue(AliRsnMiniValue::kMult, kFALSE);

   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //

   Int_t iCutK = listID1;
   Int_t iCutPi = listID2;

   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   Bool_t  use     [10] = { 1       ,  1       ,  gRsnUseMixing       ,  gRsnUseMixing       ,  1      ,  1      ,  isMC   ,   isMC   ,  isMC   ,   isMC   };
   Bool_t  useIM   [10] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      };
   TString name    [10] = {"Unlike1", "Unlike2", "Mixing1", "Mixing2", "LikePP", "LikeMM", "Trues1",  "Trues2", "Res1"  ,  "Res2"  };
   TString comp    [10] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  };
   TString output  [10] = {"HIST"   , "HIST"   , "HIST"   , "HIST"   , "HIST"  , "HIST"  , "HIST"  ,  "HIST"  , "HIST"  ,  "HIST"  };
   Char_t  charge1 [10] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     };
   Char_t  charge2 [10] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     };
   Int_t   cutID1  [10] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  };
   Int_t   cutID2  [10] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi ,  iCutPi ,  iCutPi ,   iCutPi ,  iCutPi ,   iCutPi };

   for (Int_t i = 0; i < 10; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = taskRsnMini->CreateOutput(Form("%s_%s", suffix.Data(), name[i].Data()), output[i].Data(), comp[i].Data());
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
         out->AddAxis(imID, 90, 0.6, 1.5);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);

      if (gRsnOutputFull) {
         // axis Y: transverse momentum
         out->AddAxis(ptID, 100, 0.0, 10.0);
         // axis Z: centrality
         out->AddAxis(centID, 100, 0.0, 100.0);
      }
   }

}
