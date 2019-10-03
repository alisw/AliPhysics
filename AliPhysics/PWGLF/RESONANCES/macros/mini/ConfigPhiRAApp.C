//
// *** Configuration script for phi->KK analysis with 2010 runs ***
//
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigPhiRAApp
(
   AliRsnMiniAnalysisTask *task,
   Bool_t                  isMC,
   Bool_t                  isESD,
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);

   //
   // -- Define track cuts -------------------------------------------------------------------------
   //



   gROOT->LoadMacro("AliRsnCutPhiRAA.cxx++g");
   // standard kaon cut
   AliRsnCutPhiRAA *cut = new AliRsnCutPhiRAA("cut1");
   cut->SetMode(AliRsnCutPhiRAA::k2010);
   // TPC 2 sigma pid
   AliRsnCutPIDNSigma *cutKTPC2 = new AliRsnCutPIDNSigma("cut2SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC2->SinglePIDRange(2.0);
   // TPC 3 sigma pid
   AliRsnCutPIDNSigma *cutKTPC3 = new AliRsnCutPIDNSigma("cut3SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC3->SinglePIDRange(3.0);
   // TPC 4 sigma pid
   AliRsnCutPIDNSigma *cutKTPC4 = new AliRsnCutPIDNSigma("cut4SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC4->SinglePIDRange(4.0);
   // TPC 5 sigma pid
   AliRsnCutPIDNSigma *cutKTPC5 = new AliRsnCutPIDNSigma("cut5SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC5->SinglePIDRange(5.0);
   // TPC 6 sigma pid
   AliRsnCutPIDNSigma *cutKTPC6 = new AliRsnCutPIDNSigma("cut6SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC6->SinglePIDRange(6.0);

////////////////////////// Cut Sets ///////////////////////////////////////////////

   AliRsnCutSet *cutSet = new AliRsnCutSet("set_NoPID", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutSet2 = new AliRsnCutSet("set_2sigmaTPC", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutSet3 = new AliRsnCutSet("set_3sigmaTPC", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutSet4 = new AliRsnCutSet("set_4sigmaTPC", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutSet5 = new AliRsnCutSet("set_5sigmaTPC", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutSet6 = new AliRsnCutSet("set_6sigmaTPC", AliRsnTarget::kDaughter);

////////////////////////////////////////////////////////////////////////////////////

   // no PID (only standard cuts)
   cutSet->AddCut(cut);
   cutSet->SetCutScheme(cut->GetName());

   // TPC 2 sigma cut
   cutSet2->AddCut(cut);
   cutSet2->AddCut(cutKTPC2);
   cutSet2->SetCutScheme("cut1&cut2SigmaTPCK");

   // TPC 3 sigma cut
   cutSet3->AddCut(cut);
   cutSet3->AddCut(cutKTPC3);
   cutSet3->SetCutScheme("cut1&cut3SigmaTPCK");

   // TPC 4 sigma cut
   cutSet4->AddCut(cut);
   cutSet4->AddCut(cutKTPC4);
   cutSet4->SetCutScheme("cut1&cut4SigmaTPCK");

   // TPC 5 sigma cut
   cutSet5->AddCut(cut);
   cutSet5->AddCut(cutKTPC5);
   cutSet5->SetCutScheme("cut1&cut5SigmaTPCK");

   // TPC 6 sigma cut
   cutSet6->AddCut(cut);
   cutSet6->AddCut(cutKTPC6);
   cutSet6->SetCutScheme("cut1&cut6SigmaTPCK");

//////////////////////////////////////////////////////////////////////////////
   // add to task
   Int_t icut = task->AddTrackCuts(cutSet);
   Int_t icut2 = task->AddTrackCuts(cutSet2);
   Int_t icut3 = task->AddTrackCuts(cutSet3);
   Int_t icut4 = task->AddTrackCuts(cutSet4);
   Int_t icut5 = task->AddTrackCuts(cutSet5);
   Int_t icut6 = task->AddTrackCuts(cutSet6);

   Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
   AddMonitorOutput(isMC, cutSet->GetMonitorOutput());

   //
   // -- Values ------------------------------------------------------------------------------------
   //

   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
//   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);

   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //


   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --

/////////////////// NoPID ///////////////////////////////////////////////////////////////

   Bool_t  use     [4] = { 1           ,  1           ,  1           ,  1           };
   Bool_t  useIM   [4] = { 1           ,  1           ,  1           ,  1           };
   TString name    [4] = {"UnlikeNoPID", "MixingNoPID", "LikePPNoPID", "LikeMMNoPID"};
   TString comp    [4] = {"PAIR"       , "MIX"        , "PAIR"       , "PAIR"       };
   TString output  [4] = {"HIST"       , "HIST"       , "HIST"       , "HIST"       };
   Char_t  charge1 [4] = {'+'          , '+'          , '+'          , '-'          };
   Char_t  charge2 [4] = {'-'          , '-'          , '+'          , '-'          };
   Int_t   cutID   [4] = { icut        ,  icut        ,  icut        ,  icut        };

   for (Int_t i = 0; i < 4; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID[i]);
      out->SetCutID(1, cutID[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(333);
      out->SetMotherMass(1.019455);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM)
         out->AddAxis(imID, 3000, 0.4,  7.0);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 400, 0.0, 40.0);
      // axis Z: centrality
//      out->AddAxis(centID, 100, 0.0, 100.0);
   }

////////////////////// 2 sigma cut TPC////////////////////////////////////////////////////////

   TString  name2   [4] = {"Unlike2sigmTPC", "Mixing2sigmTPC", "LikePP2sigmTPC", "LikeMM2sigmTPC"};
   Int_t    cutID2  [4] = { icut2       ,  icut2       ,  icut2       ,  icut2       };

   for (Int_t i = 0; i < 4; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name2[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID2[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(333);
      out->SetMotherMass(1.019455);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM)
         out->AddAxis(imID, 3000, 0.4,  7.0);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 400, 0.0, 40.0);
      // axis Z: centrality
//      out->AddAxis(centID, 100, 0.0, 100.0);
   }

////////////////////// 3 sigma cut TPC ////////////////////////////////////////////////////////

   TString  name3   [4] = {"Unlike3sigmTPC", "Mixing3sigmTPC", "LikePP3sigmTPC", "LikeMM3sigmTPC"};
   Int_t    cutID3  [4] = { icut3       ,  icut3       ,  icut3       ,  icut3       };

   for (Int_t i = 0; i < 4; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name3[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID3[i]);
      out->SetCutID(1, cutID3[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(333);
      out->SetMotherMass(1.019455);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM)
         out->AddAxis(imID, 3000, 0.4,  7.0);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 400, 0.0, 40.0);
      // axis Z: centrality
//      out->AddAxis(centID, 100, 0.0, 100.0);
   }

////////////////////// 4 sigma cut TPC ////////////////////////////////////////////////////////

   TString  name4   [4] = {"Unlike4sigmTPC", "Mixing4sigmTPC", "LikePP4sigmTPC", "LikeMM4sigmTPC"};
   Int_t    cutID4  [4] = { icut4       ,  icut4       ,  icut4       ,  icut4       };

   for (Int_t i = 0; i < 4; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name4[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID4[i]);
      out->SetCutID(1, cutID4[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(333);
      out->SetMotherMass(1.019455);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM)
         out->AddAxis(imID, 3000, 0.4,  7.0);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 400, 0.0, 40.0);
      // axis Z: centrality
//      out->AddAxis(centID, 100, 0.0, 100.0);
   }
////////////////////// 5s TPC  /////////////////////////////////////////////////
              TString  name5   [4] = {"Unlike5sigmTPC", "Mixing5sigmTPC", "LikePP5sigmTPC", "LikeMM5sigmTPC"};
              Int_t    cutID5  [4] = { icut5       ,  icut5       ,  icut5       ,  icut5       };

              for (Int_t i = 0; i < 4; i++) {
           	  if (!use[i]) continue;
                 // create output
                 AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name5[i].Data(), suffix), output[i].Data(), comp[i].Data());
                 // selection settings
                 out->SetCutID(0, cutID5[i]);
                 out->SetCutID(1, cutID5[i]);
                 out->SetDaughter(0, AliRsnDaughter::kKaon);
                 out->SetDaughter(1, AliRsnDaughter::kKaon);
                 out->SetCharge(0, charge1[i]);
                 out->SetCharge(1, charge2[i]);
                 out->SetMotherPDG(333);
                 out->SetMotherMass(1.019455);
                 // pair cuts
                 out->SetPairCuts(cutsPair);
                 // axis X: invmass (or resolution)
                 if (useIM)
                 out->AddAxis(imID, 3000, 0.4,  7.0);
                 else
                 out->AddAxis(resID, 200, -0.02, 0.02);
                 // axis Y: transverse momentum
                 out->AddAxis(ptID, 400, 0.0, 40.0);
                 // axis Z: centrality
                 //      out->AddAxis(centID, 100, 0.0, 100.0);
               }

////////////////////// 6s TPC  /////////////////////////////////////////////////

               TString  name6   [4] = {"Unlike6sigmTPC", "Mixing6sigmTPC", "LikePP6sigmTPC", "LikeMM6sigmTPC"};
               Int_t    cutID6  [4] = { icut6       ,  icut6       ,  icut6       ,  icut6       };

               for (Int_t i = 0; i < 4; i++) {
                 if (!use[i]) continue;
                 // create output
                 AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name6[i].Data(), suffix), output[i].Data(), comp[i].Data());
                 // selection settings
                 out->SetCutID(0, cutID6[i]);
                 out->SetCutID(1, cutID6[i]);
                 out->SetDaughter(0, AliRsnDaughter::kKaon);
                 out->SetDaughter(1, AliRsnDaughter::kKaon);
                 out->SetCharge(0, charge1[i]);
                 out->SetCharge(1, charge2[i]);
                 out->SetMotherPDG(333);
                 out->SetMotherMass(1.019455);
                 // pair cuts
                 out->SetPairCuts(cutsPair);
                 // axis X: invmass (or resolution)
                 if (useIM)
                 out->AddAxis(imID, 3000, 0.4,  7.0);
                 else
                 out->AddAxis(resID, 200, -0.02, 0.02);
                 // axis Y: transverse momentum
                 out->AddAxis(ptID, 400, 0.0, 40.0);
                 // axis Z: centrality
                 //      out->AddAxis(centID, 100, 0.0, 100.0);
               }

////////////////////// THE END! ////////////////////////////////////////////////////////


   return kTRUE;
}
