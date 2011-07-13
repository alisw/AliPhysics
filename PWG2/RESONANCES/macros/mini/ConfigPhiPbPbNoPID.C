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
Bool_t ConfigPhiPbPbNoPID
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
   
   // BB parameterization depends on data sample (MC, data)
   Double_t bbPar[5];
   if (isMC) {
      bbPar[0] = 2.15898 / 50.0;
      bbPar[1] = 1.75295E1;
      bbPar[2] = 3.40030E-9;
      bbPar[3] = 1.96178;
      bbPar[4] = 3.91720;
   } else {
      bbPar[0] = 1.41543 / 50.0;
      bbPar[1] = 2.63394E1;
      bbPar[2] = 5.0411E-11;
      bbPar[3] = 2.12543;
      bbPar[4] = 4.88663;
   }
   
   // standard kaon cut
   AliRsnCutKaonForPhi2010 *cut = new AliRsnCutKaonForPhi2010(Form("cut%s", suffix), 3.0, 3.0, 0.8);
   
   // setup (set manually the TPC PID)
   cut->SetMode(AliRsnCutKaonForPhi2010::kQuality);
   
   // cut set
   AliRsnCutSet *cutSet = new AliRsnCutSet(Form("set%s", suffix), AliRsnTarget::kDaughter);
   cutSet->AddCut(cut);
   cutSet->SetCutScheme(cut->GetName());
   
   // add to task
   Int_t icut = task->AddTrackCuts(cutSet);
   ::Info("Config", "Cut ID = %d", icut);
   
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
   Bool_t  use     [6] = { 1      ,  0      ,  0      ,  0      ,  isMC   ,  isMC     };
   Bool_t  useIM   [6] = { 1      ,  1      ,  1      ,  1      ,  1      ,  0        };
   TString name    [6] = {"Unlike", "Mixing", "LikePP", "LikeMM", "Trues" , "Res"     };
   TString comp    [6] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE"  , "TRUE"    };
   TString output  [6] = {"SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE", "SPARSE"  };
   Char_t  charge1 [6] = {'+'     , '+'     , '+'     , '-'     , '+'     , '+'       };
   Char_t  charge2 [6] = {'-'     , '-'     , '+'     , '-'     , '-'     , '-'       };
   Int_t   cutID   [6] = { icut   ,  icut   ,  icut   ,  icut   ,  icut   ,  icut     };
   
   for (Int_t i = 0; i < 6; i++) {
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
