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
Bool_t ConfigPhiPbPb
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isMC, 
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   // standard kaon cut
   AliRsnCutKaonForPhi2010 *cutStd = new AliRsnCutKaonForPhi2010("cutStdPbPb", 3.0, 3.0, 0.8);
   // cut set
   AliRsnCutSet *cutSetStd = new AliRsnCutSet("setStdPbPb", AliRsnTarget::kDaughter);
   cutSetStd->AddCut(cutStd);
   cutSetStd->SetCutScheme(cutStd->GetName());
   // add to task
   Int_t icutStd = task->AddTrackCuts(cutSetStd);
   
   // TPC kaon cut
   AliRsnCutKaonForPhi2010 *cutTPC = new AliRsnCutKaonForPhi2010("cutTPCPbPb", 3.0, 3.0, 0.8);
   cutTPC->SetOnlyTPC();
   // cut set
   AliRsnCutSet *cutSetTPC = new AliRsnCutSet("setTPCPbPb", AliRsnTarget::kDaughter);
   cutSetTPC->AddCut(cutTPC);
   cutSetTPC->SetCutScheme(cutTPC->GetName());
   // add to task
   Int_t icutTPC = task->AddTrackCuts(cutSetTPC);
   
   // TOF kaon cut
   AliRsnCutKaonForPhi2010 *cutTOF = new AliRsnCutKaonForPhi2010("cutTOFPbPb", 3.0, 3.0, 0.8);
   cutTOF->SetOnlyTOF();
   // cut set
   AliRsnCutSet *cutSetTOF = new AliRsnCutSet("setTOFPbPb", AliRsnTarget::kDaughter);
   cutSetTOF->AddCut(cutTOF);
   cutSetTOF->SetCutScheme(cutTOF->GetName());
   // add to task
   Int_t icutTOF = task->AddTrackCuts(cutSetTOF);
   
   // No-PID kaon cut
   AliRsnCutKaonForPhi2010 *cutNOPID = new AliRsnCutKaonForPhi2010("cutNOPIDPbPb", 3.0, 3.0, 0.8);
   cutNOPID->SetOnlyQuality();
   // cut set
   AliRsnCutSet *cutSetNOPID = new AliRsnCutSet("setNOPIDPbPb", AliRsnTarget::kDaughter);
   cutSetNOPID->AddCut(cutNOPID);
   cutSetNOPID->SetCutScheme(cutNOPID->GetName());
   // add to task
   Int_t icutNOPID = task->AddTrackCuts(cutSetNOPID);
   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
   
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
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
   //Bool_t  use     [6] = { 1      ,  1      ,  1      ,  1      ,  isMC   ,  isMC     };
   //Bool_t  useIM   [6] = { 1      ,  1      ,  1      ,  1      ,  1      ,  0        };
   //TString name    [6] = {"Unlike", "Mixing", "LikePP", "LikeMM", "Trues" , "Res"     };
   //TString comp    [6] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE"  , "TRUE"    };
   //TString output  [6] = {"HIST"  , "HIST"  , "HIST"  , "HIST"  , "HIST"  , "HIST"    };
   //Char_t  charge1 [6] = {'+'     , '+'     , '+'     , '-'     , '+'     , '+'       };
   //Char_t  charge2 [6] = {'-'     , '-'     , '+'     , '-'     , '-'     , '-'       };
   //Int_t   cutID   [6] = { icut   ,  icut   ,  icut   ,  icut   ,  icut   ,  icut     };
   
   Bool_t  use     [12] = { 1      ,  1      ,  1       ,  1      ,  1      ,  1       ,  1      ,  1      ,  1       ,  1        ,  1        ,  1          };
   Bool_t  useIM   [12] = { 1      ,  1      ,  1       ,  1      ,  1      ,  1       ,  1      ,  1      ,  1       ,  1        ,  1        ,  1          };
   TString name    [12] = {"STDPM" , "STDMIX", "STDTRUE", "TPCPM" , "TPCMIX", "TPCTRUE", "TOFPM" , "TOFMIX", "TOFTRUE", "NOPIDPM" , "NOPIDMIX", "NOPIDTRUE" };
   TString comp    [12] = {"PAIR"  , "MIX"   , "TRUE"   , "PAIR"  , "MIX"   , "TRUE"   , "PAIR"  , "MIX"   , "TRUE"   , "PAIR"    , "MIX"     , "TRUE"      };
   TString output  [12] = {"SPARSE", "SPARSE", "SPARSE" , "SPARSE", "SPARSE", "SPARSE" , "SPARSE", "SPARSE", "SPARSE" , "SPARSE"  , "SPARSE"  , "SPARSE"    };
   Char_t  charge1 [12] = {'+'     , '+'     , '+'      , '+'     , '+'     , '+'      , '+'     , '+'     , '+'      , '+'       , '+'       , '+'         };
   Char_t  charge2 [12] = {'-'     , '-'     , '-'      , '-'     , '-'     , '-'      , '-'     , '-'     , '-'      , '-'       , '-'       , '-'         };
   Int_t   cutID   [12] = { icutStd,  icutStd,  icutStd ,  icutTPC,  icutTPC,  icutTPC ,  icutTOF,  icutTOF,  icutTOF ,  icutNOPID,  icutNOPID,  icutNOPID  };
   
   for (Int_t i = 0; i < 12; i++) {
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
      out->AddAxis(ptID, 100, 0.0, 10.0);
      // axis Z: centrality 
      out->AddAxis(centID, 100, 0.0, 100.0);
   }
   
   return kTRUE;
}
