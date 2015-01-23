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
Bool_t ConfigPhiRAApPb
(
   AliRsnMiniAnalysisTask *task,
   Bool_t                  isMC,
   Bool_t                  isESD,
   const char             *suffix,
   AliRsnCutSet           *cutsPair,
   AliRsnCutSet           *cutsPair2
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
   cut->SetMode(AliRsnCutPhiRAA::k2011_1);


   // TPC 3 sigma pid
   AliRsnCutPIDNSigma *cutKTPC3 = new AliRsnCutPIDNSigma("cut3SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC3->SinglePIDRange(3.0);
   // TPC 2 sigma pid
   AliRsnCutPIDNSigma *cutKTPC2 = new AliRsnCutPIDNSigma("cut2SigmaTPCK",AliPID::kKaon,AliRsnCutPIDNSigma::kTPC);
   cutKTPC2->SinglePIDRange(2.0);
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
      AliRsnCutSet *cutSet5 = new AliRsnCutSet("set_2sigmaTPC", AliRsnTarget::kDaughter);
      AliRsnCutSet *cutSet7 = new AliRsnCutSet("set_3sigmaTPC", AliRsnTarget::kDaughter);
      AliRsnCutSet *cutSet12 = new AliRsnCutSet("set_4sigmaTPC", AliRsnTarget::kDaughter);
      AliRsnCutSet *cutSet13 = new AliRsnCutSet("set_5sigmaTPC", AliRsnTarget::kDaughter);
      AliRsnCutSet *cutSet14 = new AliRsnCutSet("set_6sigmaTPC", AliRsnTarget::kDaughter);

////////////////////////////////////////////////////////////////////////////////////
      // no PID (only standard cuts)
      cutSet->AddCut(cut);
      cutSet->SetCutScheme(cut->GetName());

      // TPC 2 sigma cut
      cutSet5->AddCut(cut);
      cutSet5->AddCut(cutKTPC2);
      cutSet5->SetCutScheme("cut1&cut2SigmaTPCK");

      // TPC 3 sigma cut
      cutSet7->AddCut(cut);
      cutSet7->AddCut(cutKTPC3);
      cutSet7->SetCutScheme("cut1&cut3SigmaTPCK");

      // TPC 4 sigma cut
      cutSet12->AddCut(cut);
      cutSet12->AddCut(cutKTPC4);
      cutSet12->SetCutScheme("cut1&cut4SigmaTPCK");

      // TPC 5 sigma cut
      cutSet13->AddCut(cut);
      cutSet13->AddCut(cutKTPC5);
      cutSet13->SetCutScheme("cut1&cut5SigmaTPCK");

      // TPC 6 sigma cut
      cutSet14->AddCut(cut);
      cutSet14->AddCut(cutKTPC6);
      cutSet14->SetCutScheme("cut1&cut6SigmaTPCK");


//////////////////////////////////////////////////////////////////////////////



   // add to task
      Int_t icut = task->AddTrackCuts(cutSet);
      Int_t icut5 = task->AddTrackCuts(cutSet5);
      Int_t icut7 = task->AddTrackCuts(cutSet7);
      Int_t icut12 = task->AddTrackCuts(cutSet12);
      Int_t icut13 = task->AddTrackCuts(cutSet13);
      Int_t icut14 = task->AddTrackCuts(cutSet14);

      Printf("======== Monitoring cut AliRsnCutSetDaughterParticle enabled");
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/mini/AddMonitorOutput.C");
      AddMonitorOutput(isMC, cutSet->GetMonitorOutput());

   //
   // -- Values ------------------------------------------------------------------------------------
   //

   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
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

/////////////////// NoPID 03 ///////////////////////////////////////////////////////////////

         Bool_t  use     [4] = { 1           ,  1           ,  1           ,  1           };
         Bool_t  useIM   [4] = { 1           ,  1           ,  1           ,  1           };
         TString name    [4] = {"UnlikeNoPID_03", "MixingNoPID_03", "LikePPNoPID_03", "LikeMMNoPID_03"};
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
               out->AddAxis(imID, 500, 0.9,  1.4);
            else
               out->AddAxis(resID, 200, -0.02, 0.02);
            // axis Y: transverse momentum
            out->AddAxis(ptID, 400, 0.0, 40.0);
            // axis Z: centrality
      //      out->AddAxis(centID, 100, 0.0, 100.0);
         }

/////////////////// NoPID 05 ///////////////////////////////////////////////////////////////

         TString name_1    [4] = {"UnlikeNoPID_05", "MixingNoPID_05", "LikePPNoPID_05", "LikeMMNoPID_05"};

         for (Int_t i = 0; i < 4; i++) {
           if (!use[i]) continue;
           // create output
           AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
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
           out->SetPairCuts(cutsPair2);
           // axis X: invmass (or resolution)
           if (useIM)
             out->AddAxis(imID, 500, 0.9,  1.4);
           else
             out->AddAxis(resID, 200, -0.02, 0.02);
           // axis Y: transverse momentum
             out->AddAxis(ptID, 400, 0.0, 40.0);
           // axis Z: centrality
           //      out->AddAxis(centID, 100, 0.0, 100.0);
         }

////////////////////// 2s TPC 03 /////////////////////////////////////////////////

        TString  name5   [4] = {"Unlike2sigmaTPC_03", "Mixing2sigmaTPC_03", "LikePP2sigmaTPC_03", "LikeMM2sigmaTPC_03"};
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
            out->AddAxis(imID, 500, 0.9,  1.4);
          else
            out->AddAxis(resID, 200, -0.02, 0.02);
          // axis Y: transverse momentum
          out->AddAxis(ptID, 400, 0.0, 40.0);
          // axis Z: centrality
          //      out->AddAxis(centID, 100, 0.0, 100.0);
        }

////////////////////// 2s TPC 05 /////////////////////////////////////////////////

        TString  name5_1   [4] = {"Unlike2sigmaTPC_05", "Mixing2sigmaTPC_05", "LikePP2sigmaTPC_05", "LikeMM2sigmaTPC_05"};

        for (Int_t i = 0; i < 4; i++) {
          if (!use[i]) continue;
          // create output
          AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name5_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
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
          out->SetPairCuts(cutsPair2);
          // axis X: invmass (or resolution)
          if (useIM)
            out->AddAxis(imID, 500, 0.9,  1.4);
          else
            out->AddAxis(resID, 200, -0.02, 0.02);
          // axis Y: transverse momentum
          out->AddAxis(ptID, 400, 0.0, 40.0);
          // axis Z: centrality
          //      out->AddAxis(centID, 100, 0.0, 100.0);
        }

////////////////////// 3s TPC  03 //////////////////////////////////////////////////////////

       TString  name7   [4] = {"Unlike3sigmTPC_03", "Mixing3sigmTPC_03", "LikePP3sigmTPC_03", "LikeMM3sigmTPC_03"};
       Int_t    cutID7  [4] = { icut7       ,  icut7       ,  icut7       ,  icut7       };

       for (Int_t i = 0; i < 4; i++) {
         if (!use[i]) continue;
         // create output
         AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name7[i].Data(), suffix), output[i].Data(), comp[i].Data());
         // selection settings
         out->SetCutID(0, cutID7[i]);
         out->SetCutID(1, cutID7[i]);
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
           out->AddAxis(imID, 500, 0.9,  1.4);
         else
           out->AddAxis(resID, 200, -0.02, 0.02);
         // axis Y: transverse momentum
         out->AddAxis(ptID, 400, 0.0, 40.0);
         // axis Z: centrality
         //      out->AddAxis(centID, 100, 0.0, 100.0);
       }

////////////////////// 3s TPC 05 //////////////////////////////////////////////////////////

              TString  name7_1   [4] = {"Unlike3sigmTPC_05", "Mixing3sigmTPC_05", "LikePP3sigmTPC_05", "LikeMM3sigmTPC_05"};

              for (Int_t i = 0; i < 4; i++) {
                if (!use[i]) continue;
                // create output
                AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name7_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
                // selection settings
                out->SetCutID(0, cutID7[i]);
                out->SetCutID(1, cutID7[i]);
                out->SetDaughter(0, AliRsnDaughter::kKaon);
                out->SetDaughter(1, AliRsnDaughter::kKaon);
                out->SetCharge(0, charge1[i]);
                out->SetCharge(1, charge2[i]);
                out->SetMotherPDG(333);
                out->SetMotherMass(1.019455);
                // pair cuts
                out->SetPairCuts(cutsPair2);
                // axis X: invmass (or resolution)
                if (useIM)
                  out->AddAxis(imID, 500, 0.9,  1.4);
                else
                  out->AddAxis(resID, 200, -0.02, 0.02);
                // axis Y: transverse momentum
                out->AddAxis(ptID, 400, 0.0, 40.0);
                // axis Z: centrality
                //      out->AddAxis(centID, 100, 0.0, 100.0);
              }


////////////////////// 4s TPC  03/////////////////////////////////////////////////

      TString  name12   [4] = {"Unlike4sigmaTPC_03", "Mixing4sigmaTPC_03", "LikePP4sigmaTPC_03", "LikeMM4sigmaTPC_03"};
      Int_t    cutID12  [4] = { icut12       ,  icut12       ,  icut12       ,  icut12       };

      for (Int_t i = 0; i < 4; i++) {
        if (!use[i]) continue;
        // create output
        AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name12[i].Data(), suffix), output[i].Data(), comp[i].Data());
        // selection settings
        out->SetCutID(0, cutID12[i]);
        out->SetCutID(1, cutID12[i]);
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
          out->AddAxis(imID, 500, 0.9,  1.4);
        else
          out->AddAxis(resID, 200, -0.02, 0.02);
        // axis Y: transverse momentum
        out->AddAxis(ptID, 400, 0.0, 40.0);
        // axis Z: centrality
        //      out->AddAxis(centID, 100, 0.0, 100.0);
      }

////////////////////// 4s TPC  05 /////////////////////////////////////////////////

            TString  name12_1   [4] = {"Unlike4sigmaTPC_05", "Mixing4sigmaTPC_05", "LikePP4sigmaTPC_05", "LikeMM4sigmaTPC_05"};

            for (Int_t i = 0; i < 4; i++) {
              if (!use[i]) continue;
              // create output
              AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name12_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
              // selection settings
              out->SetCutID(0, cutID12[i]);
              out->SetCutID(1, cutID12[i]);
              out->SetDaughter(0, AliRsnDaughter::kKaon);
              out->SetDaughter(1, AliRsnDaughter::kKaon);
              out->SetCharge(0, charge1[i]);
              out->SetCharge(1, charge2[i]);
              out->SetMotherPDG(333);
              out->SetMotherMass(1.019455);
              // pair cuts
              out->SetPairCuts(cutsPair2);
              // axis X: invmass (or resolution)
              if (useIM)
                out->AddAxis(imID, 500, 0.9,  1.4);
              else
                out->AddAxis(resID, 200, -0.02, 0.02);
              // axis Y: transverse momentum
              out->AddAxis(ptID, 400, 0.0, 40.0);
              // axis Z: centrality
              //      out->AddAxis(centID, 100, 0.0, 100.0);
            }

////////////////////// 5s TPC  03/////////////////////////////////////////////////

       TString  name13   [4] = {"Unlike5sigmaTPC_03", "Mixing5sigmaTPC_03", "LikePP5sigmaTPC_03", "LikeMM5sigmaTPC_03"};
       Int_t    cutID13  [4] = { icut13       ,  icut13       ,  icut13       ,  icut13       };

       for (Int_t i = 0; i < 4; i++) {
         if (!use[i]) continue;
         // create output
         AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name13[i].Data(), suffix), output[i].Data(), comp[i].Data());
         // selection settings
         out->SetCutID(0, cutID13[i]);
         out->SetCutID(1, cutID13[i]);
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
           out->AddAxis(imID, 500, 0.9,  1.4);
         else
           out->AddAxis(resID, 200, -0.02, 0.02);
         // axis Y: transverse momentum
         out->AddAxis(ptID, 400, 0.0, 40.0);
         // axis Z: centrality
         //      out->AddAxis(centID, 100, 0.0, 100.0);
       }

////////////////////// 5s TPC  05 /////////////////////////////////////////////////

       TString  name13_1   [4] = {"Unlike5sigmaTPC_05", "Mixing5sigmaTPC_05", "LikePP5sigmaTPC_05", "LikeMM5sigmaTPC_05"};

       for (Int_t i = 0; i < 4; i++) {
         if (!use[i]) continue;
         // create output
         AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name13_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
         // selection settings
         out->SetCutID(0, cutID13[i]);
         out->SetCutID(1, cutID13[i]);
         out->SetDaughter(0, AliRsnDaughter::kKaon);
         out->SetDaughter(1, AliRsnDaughter::kKaon);
         out->SetCharge(0, charge1[i]);
         out->SetCharge(1, charge2[i]);
         out->SetMotherPDG(333);
         out->SetMotherMass(1.019455);
         // pair cuts
         out->SetPairCuts(cutsPair2);
         // axis X: invmass (or resolution)
         if (useIM)
           out->AddAxis(imID, 500, 0.9,  1.4);
         else
           out->AddAxis(resID, 200, -0.02, 0.02);
         // axis Y: transverse momentum
         out->AddAxis(ptID, 400, 0.0, 40.0);
         // axis Z: centrality
         //      out->AddAxis(centID, 100, 0.0, 100.0);
       }


////////////////////// 6s TPC  03/////////////////////////////////////////////////

      TString  name14   [4] = {"Unlike6sigmaTPC_03", "Mixing6sigmaTPC_03", "LikePP6sigmaTPC_03", "LikeMM6sigmaTPC_03"};
      Int_t    cutID14  [4] = { icut14       ,  icut14       ,  icut14       ,  icut14       };

      for (Int_t i = 0; i < 4; i++) {
        if (!use[i]) continue;
        // create output
        AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name14[i].Data(), suffix), output[i].Data(), comp[i].Data());
        // selection settings
        out->SetCutID(0, cutID14[i]);
        out->SetCutID(1, cutID14[i]);
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
          out->AddAxis(imID, 500, 0.9,  1.4);
        else
          out->AddAxis(resID, 200, -0.02, 0.02);
        // axis Y: transverse momentum
        out->AddAxis(ptID, 400, 0.0, 40.0);
        // axis Z: centrality
        //      out->AddAxis(centID, 100, 0.0, 100.0);
      }

////////////////////// 6s TPC  05 /////////////////////////////////////////////////

       TString  name14_1   [4] = {"Unlike6sigmaTPC_05", "Mixing6sigmaTPC_05", "LikePP6sigmaTPC_05", "LikeMM6sigmaTPC_05"};

       for (Int_t i = 0; i < 4; i++) {
         if (!use[i]) continue;
         // create output
         AliRsnMiniOutput *out = task->CreateOutput(Form("phi_%s%s", name14_1[i].Data(), suffix), output[i].Data(), comp[i].Data());
         // selection settings
         out->SetCutID(0, cutID14[i]);
         out->SetCutID(1, cutID14[i]);
         out->SetDaughter(0, AliRsnDaughter::kKaon);
         out->SetDaughter(1, AliRsnDaughter::kKaon);
         out->SetCharge(0, charge1[i]);
         out->SetCharge(1, charge2[i]);
         out->SetMotherPDG(333);
         out->SetMotherMass(1.019455);
         // pair cuts
         out->SetPairCuts(cutsPair2);
         // axis X: invmass (or resolution)
         if (useIM)
           out->AddAxis(imID, 500, 0.9,  1.4);
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
