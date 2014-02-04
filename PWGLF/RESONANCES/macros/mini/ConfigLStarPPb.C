//
// *** Configuration script for Lambda*->P+K- analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigLStarPPb
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isMC, 
   Bool_t                  isPP, 
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   //if (strlen(suffix) > 0) suffix = "suffix_";
   
   ////////////////

   Int_t aodFilterBit = 10;

   
   ///////////////

   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   Float_t nSigmaTPC=2.0;  
   gROOT->LoadMacro("AddMonitorOutput.C");


   AliRsnCutSetDaughterParticle * cutP;
   cutP  = new AliRsnCutSetDaughterParticle("P_LStar",AliRsnCutSetDaughterParticle::kTPCpidMatchPPB2011 , AliPID::kProton, nSigmaTPC, aodFilterBit);

   AliRsnCutSetDaughterParticle * cutK;
   cutK  = new AliRsnCutSetDaughterParticle("K_LStar",AliRsnCutSetDaughterParticle::kTPCpidMatchPPB2011 , AliPID::kKaon, nSigmaTPC, aodFilterBit);
   
 
   //Int_t iCutQ = task->AddTrackCuts(cutQ);
   Int_t iCutP = task->AddTrackCuts(cutP);
   Int_t iCutK = task->AddTrackCuts(cutK);

   //AddMonitorOutput(isMC, cutQ->GetMonitorOutput(),"NoSIGN");
   AddMonitorOutput(isMC, cutP->GetMonitorOutput(),"NoSIGN");
   AddMonitorOutput(isMC, cutK->GetMonitorOutput(),"NoSIGN");
   

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
   Bool_t  use     [12] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  isMC   ,   isMC   ,  isMC   ,   isMC   , isMC , isMC };
   Bool_t  useIM   [12] = { 1       ,  1       ,  1       ,  1       ,  1      ,  1      ,  1      ,   1      ,  0      ,   0      ,  1   ,   1};
   TString name    [12] = {"Unlike1", "Unlike2", "Mixing1", "Mixing2", "LikePP", "LikeMM", "Trues1",  "Trues2", "Res1"  ,  "Res2"  ,"Mother1", "Mother2"};
   TString comp    [12] = {"PAIR"   , "PAIR"   , "MIX"    , "MIX"    , "PAIR"  , "PAIR"  , "TRUE"  ,  "TRUE"  , "TRUE"  ,  "TRUE"  , "MOTHER",  "MOTHER"};
   TString output  [12] = {"SPARSE" , "SPARSE" , "SPARSE" , "SPARSE" , "SPARSE", "SPARSE", "SPARSE",  "SPARSE", "SPARSE",  "SPARSE", "SPARSE",  "SPARSE"  };
   Char_t  charge1 [12] = {'+'      , '-'      , '+'      , '-'      , '+'     , '-'     , '+'     ,  '-'     , '+'     ,  '-'     , '+'     ,  '-'};
   Char_t  charge2 [12] = {'-'      , '+'      , '-'      , '+'      , '+'     , '-'     , '-'     ,  '+'     , '-'     ,  '+'     , '-'     ,  '+'};
   Int_t   cutID1  [12] = { iCutP   ,  iCutP   ,  iCutP   ,  iCutP   ,  iCutP  ,  iCutP  ,  iCutP  ,   iCutP  ,  iCutP  ,   iCutP  ,  iCutP  ,   iCutP  };
   Int_t   cutID2  [12] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK  ,  iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK  ,  iCutK  ,   iCutK };
   ///////


   ////////   
   for (Int_t i = 0; i < 12; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("Lstar_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kProton);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(3124);
      out->SetMotherMass(1.5195);
      // pair cuts
      out->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i]) 
         out->AddAxis(imID, 900, 1.3, 2.2);
      else
         out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 1000, 0.0, 100.0);
    // axis Z: centrality-multiplicity
    //if (!isPP)
    //  out->AddAxis(centID, 100, 0.0, 100.0);
    //else 
    //  out->AddAxis(centID, 400, 0.0, 400.0);
      
    // axis W: pseudorapidity
    //    out->AddAxis(etaID, 20, -1.0, 1.0);
    // axis J: rapidity
    // out->AddAxis(yID, 10, -0.5, 0.5);
 
   }
   
   return kTRUE;
}
