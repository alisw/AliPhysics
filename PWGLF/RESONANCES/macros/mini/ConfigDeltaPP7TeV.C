//
// *** Configuration script for delta(++), delta(--) and delta(0) analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigDeltaPP7TeV
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isMC, 
   const char             *suffix,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);


   // integrated pion cut TOF
   AliRsnCutDelta *cutPiTOF = new AliRsnCutDelta("cutPionTOF",AliPID::kPion,kFALSE);
   // cut set TOF
   AliRsnCutSet *cutSetPiTOF = new AliRsnCutSet("setPionForDeltaTOF", AliRsnTarget::kDaughter);
   cutSetPiTOF->AddCut(cutPiTOF);
   cutSetPiTOF->SetCutScheme(cutPiTOF->GetName());   
   // add to task TOF
   Int_t iCutPiTOF = task->AddTrackCuts(cutSetPiTOF);
   
   // integrated proton cut TOF
   AliRsnCutDelta *cutPTOF = new AliRsnCutDelta("cutProtonTOF",AliPID::kProton,kFALSE);
   // cut set TOF
   AliRsnCutSet *cutSetPTOF = new AliRsnCutSet("setProtonForDeltaTOF", AliRsnTarget::kDaughter);
   cutSetPTOF->AddCut(cutPTOF);
   cutSetPTOF->SetCutScheme(cutPTOF->GetName());
   // add to task TOF
   Int_t iCutPTOF = task->AddTrackCuts(cutSetPTOF);

/////////////////////////////////////////////////////SYSTEMATICS
 
   AliRsnCutDelta *cutPiTOF1 = new AliRsnCutDelta("cutPionTOF1",AliPID::kPion,kFALSE);
   cutPiTOF1->SetTPCNSigmaProton(2.0);
   cutPiTOF1->SetTPCNSigmaPion(3.0);
   cutPiTOF1->SetTOFNSigmaProton(2.0);
   cutPiTOF1->SetTOFNSigmaPion(3.0);
   AliRsnCutSet *cutSetPiTOF1 = new AliRsnCutSet("setPionForDeltaTOF1", AliRsnTarget::kDaughter);
   cutSetPiTOF1->AddCut(cutPiTOF1);
   cutSetPiTOF1->SetCutScheme(cutPiTOF1->GetName());
   Int_t iCutPiTOF1 = task->AddTrackCuts(cutSetPiTOF1);

   AliRsnCutDelta *cutPTOF1 = new AliRsnCutDelta("cutProtonTOF1",AliPID::kProton,kFALSE);
   cutPTOF1->SetTPCNSigmaProton(2.0);
   cutPTOF1->SetTPCNSigmaPion(3.0);
   cutPTOF1->SetTOFNSigmaProton(2.0);
   cutPTOF1->SetTOFNSigmaPion(3.0);
   AliRsnCutSet *cutSetPTOF1 = new AliRsnCutSet("setProtonForDeltaTOF1", AliRsnTarget::kDaughter);
   cutSetPTOF1->AddCut(cutPTOF1);
   cutSetPTOF1->SetCutScheme(cutPTOF1->GetName());
   Int_t iCutPTOF1 = task->AddTrackCuts(cutSetPTOF1);

   AliRsnCutDelta *cutPiTOF2 = new AliRsnCutDelta("cutPionTOF2",AliPID::kPion,kFALSE);
   cutPiTOF2->SetTPCNSigmaProton(3.0);
   cutPiTOF2->SetTPCNSigmaPion(2.0);
   cutPiTOF2->SetTOFNSigmaProton(3.0);
   cutPiTOF2->SetTOFNSigmaPion(2.0);
   AliRsnCutSet *cutSetPiTOF2 = new AliRsnCutSet("setPionForDeltaTOF2", AliRsnTarget::kDaughter);
   cutSetPiTOF2->AddCut(cutPiTOF2);
   cutSetPiTOF2->SetCutScheme(cutPiTOF2->GetName());
   Int_t iCutPiTOF2 = task->AddTrackCuts(cutSetPiTOF2);

   AliRsnCutDelta *cutPTOF2 = new AliRsnCutDelta("cutProtonTOF2",AliPID::kProton,kFALSE);
   cutPTOF2->SetTPCNSigmaProton(3.0);
   cutPTOF2->SetTPCNSigmaPion(2.0);
   cutPTOF2->SetTOFNSigmaProton(3.0);
   cutPTOF2->SetTOFNSigmaPion(2.0);
   AliRsnCutSet *cutSetPTOF2 = new AliRsnCutSet("setProtonForDeltaTOF2", AliRsnTarget::kDaughter);
   cutSetPTOF2->AddCut(cutPTOF2);
   cutSetPTOF2->SetCutScheme(cutPTOF2->GetName());
   Int_t iCutPTOF2 = task->AddTrackCuts(cutSetPTOF2); 

   AliRsnCutDelta *cutPiTOF3 = new AliRsnCutDelta("cutPionTOF3",AliPID::kPion,kFALSE);
   cutPiTOF3->SetTPCNSigmaProton(4.0);
   cutPiTOF3->SetTPCNSigmaPion(3.0);
   cutPiTOF3->SetTOFNSigmaProton(4.0);
   cutPiTOF3->SetTOFNSigmaPion(3.0);
   AliRsnCutSet *cutSetPiTOF3 = new AliRsnCutSet("setPionForDeltaTOF3", AliRsnTarget::kDaughter);
   cutSetPiTOF3->AddCut(cutPiTOF3);
   cutSetPiTOF3->SetCutScheme(cutPiTOF3->GetName());
   Int_t iCutPiTOF3 = task->AddTrackCuts(cutSetPiTOF3);

   AliRsnCutDelta *cutPTOF3 = new AliRsnCutDelta("cutProtonTOF3",AliPID::kProton,kFALSE);
   cutPiTOF3->SetTPCNSigmaProton(4.0);
   cutPiTOF3->SetTPCNSigmaPion(3.0);
   cutPiTOF3->SetTOFNSigmaProton(4.0);
   cutPiTOF3->SetTOFNSigmaPion(3.0);
   AliRsnCutSet *cutSetPTOF3 = new AliRsnCutSet("setProtonForDeltaTOF3", AliRsnTarget::kDaughter);
   cutSetPTOF3->AddCut(cutPTOF3);
   cutSetPTOF3->SetCutScheme(cutPTOF3->GetName());
   Int_t iCutPTOF3 = task->AddTrackCuts(cutSetPTOF3); 

   AliRsnCutDelta *cutPiTOF4 = new AliRsnCutDelta("cutPionTOF4",AliPID::kPion,kFALSE);
   cutPiTOF4->SetTPCNSigmaProton(3.0);
   cutPiTOF4->SetTPCNSigmaPion(4.0);
   cutPiTOF4->SetTOFNSigmaProton(3.0);
   cutPiTOF4->SetTOFNSigmaPion(4.0);
   AliRsnCutSet *cutSetPiTOF4 = new AliRsnCutSet("setPionForDeltaTOF4", AliRsnTarget::kDaughter);
   cutSetPiTOF4->AddCut(cutPiTOF4);
   cutSetPiTOF4->SetCutScheme(cutPiTOF4->GetName());
   Int_t iCutPiTOF4 = task->AddTrackCuts(cutSetPiTOF4);

   AliRsnCutDelta *cutPTOF4 = new AliRsnCutDelta("cutProtonTOF4",AliPID::kProton,kFALSE);
   cutPTOF4->SetTPCNSigmaProton(3.0);
   cutPTOF4->SetTPCNSigmaPion(4.0);
   cutPTOF4->SetTOFNSigmaProton(3.0);
   cutPTOF4->SetTOFNSigmaPion(4.0);
   AliRsnCutSet *cutSetPTOF4 = new AliRsnCutSet("setProtonForDeltaTOF4", AliRsnTarget::kDaughter);
   cutSetPTOF4->AddCut(cutPTOF4);
   cutSetPTOF4->SetCutScheme(cutPTOF4->GetName());
   Int_t iCutPTOF4 = task->AddTrackCuts(cutSetPTOF4);

   AliRsnCutDelta *cutPiTOF5 = new AliRsnCutDelta("cutPionTOF5",AliPID::kPion,kFALSE);
   cutPiTOF5->SetTOFMomProton(2.3);
   AliRsnCutSet *cutSetPiTOF5 = new AliRsnCutSet("setPionForDeltaTOF5", AliRsnTarget::kDaughter);
   cutSetPiTOF5->AddCut(cutPiTOF5);
   cutSetPiTOF5->SetCutScheme(cutPiTOF5->GetName());
   Int_t iCutPiTOF5 = task->AddTrackCuts(cutSetPiTOF5);

   AliRsnCutDelta *cutPTOF5 = new AliRsnCutDelta("cutProtonTOF5",AliPID::kProton,kFALSE);
   cutPTOF5->SetTOFMomProton(2.3);
   AliRsnCutSet *cutSetPTOF5 = new AliRsnCutSet("setProtonForDeltaTOF5", AliRsnTarget::kDaughter);
   cutSetPTOF5->AddCut(cutPTOF5);
   cutSetPTOF5->SetCutScheme(cutPTOF5->GetName());
   Int_t iCutPTOF5 = task->AddTrackCuts(cutSetPTOF5);

   AliRsnCutDelta *cutPiTOF6 = new AliRsnCutDelta("cutPionTOF6",AliPID::kPion,kFALSE);
   cutPiTOF6->SetTOFMomProton(2.7);
   AliRsnCutSet *cutSetPiTOF6 = new AliRsnCutSet("setPionForDeltaTOF6", AliRsnTarget::kDaughter);
   cutSetPiTOF6->AddCut(cutPiTOF6);
   cutSetPiTOF6->SetCutScheme(cutPiTOF6->GetName());
   Int_t iCutPiTOF6 = task->AddTrackCuts(cutSetPiTOF6);

   AliRsnCutDelta *cutPTOF6 = new AliRsnCutDelta("cutProtonTOF6",AliPID::kProton,kFALSE);
   cutPTOF6->SetTOFMomProton(2.7);
   AliRsnCutSet *cutSetPTOF6 = new AliRsnCutSet("setProtonForDeltaTOF6", AliRsnTarget::kDaughter);
   cutSetPTOF6->AddCut(cutPTOF6);
   cutSetPTOF6->SetCutScheme(cutPTOF6->GetName());
   Int_t iCutPTOF6 = task->AddTrackCuts(cutSetPTOF6);

/////////////////////////////////////////////////////SYSTEMATICS END


   
   //
   // -- Values ------------------------------------------------------------------------------------
   //
//CreateValue(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE); 

   /* IM resolution    */ Int_t resID  = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
   /* invariant mass   */ Int_t imID   = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
   /* transv. momentum */ Int_t ptID   = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);


   /* 1st daughter pt  */ Int_t fdpt   = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt  , kFALSE);
   /* 2nd daughter pt  */ Int_t sdpt   = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt , kFALSE);
   /* 1st daughter p   */ Int_t fdp    = task->CreateValue(AliRsnMiniValue::kFirstDaughterP   , kFALSE);
   /* 2nd daughter p   */ Int_t sdp    = task->CreateValue(AliRsnMiniValue::kSecondDaughterP  , kFALSE);
   
   /* pseudorapidity   */ Int_t etaID  = task->CreateValue(AliRsnMiniValue::kEta , kFALSE);
   /* rapidity         */ Int_t yID    = task->CreateValue(AliRsnMiniValue::kY   , kFALSE);


   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
      
   Bool_t  use      [12] = { 1       	,  1          ,  1         ,  1         ,  1         ,  1         ,  isMC       ,   isMC       ,  1          ,   1         ,  isMC       ,   isMC         };
   Bool_t  useIM    [12] = { 1       	,  1          ,  1         ,  1         ,  1         ,  1         ,  1          ,   1          ,  1          ,   1         ,  1          ,   1            };
   TString name     [12] = {"UnlikePM"	, "UnlikeMP"  , "MixingPM" , "MixingMP" , "LikePP"   , "LikeMM"   , "TruesPP"   ,  "TruesMM"   , "MixingPP"  ,  "MixingMM" , "TRUESPM"   ,  "TRUESMP"     };
   TString comp     [12] = {"PAIR"   	, "PAIR"      , "MIX"      , "MIX"      , "PAIR"     , "PAIR"     , "TRUE"      ,  "TRUE"      , "MIX"       ,  "MIX"      , "TRUE"      ,  "TRUE"        };
   TString output   [12] = {"SPARSE"   	, "SPARSE"    , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"   , "SPARSE"    ,  "SPARSE"    , "SPARSE"    ,  "SPARSE"   , "SPARSE"    ,  "SPARSE"      };
   Char_t  charge1  [12] = {'+'      	, '-'         , '+'        , '-'        , '+'        , '-'        , '+'         ,  '-'         , '+'         ,  '-'        , '+'         ,  '-'           };
   Char_t  charge2  [12] = {'-'      	, '+'         , '-'        , '+'        , '+'        , '-'        , '+'         ,  '-'         , '+'         ,  '-'        , '-'         ,  '+'           };
   Int_t   cutID1   [12] = {iCutPTOF    , iCutPTOF    , iCutPTOF   , iCutPTOF   , iCutPTOF   , iCutPTOF   , iCutPTOF    ,  iCutPTOF    , iCutPTOF    ,  iCutPTOF   , iCutPTOF    ,  iCutPTOF      };
   Int_t   cutID2   [12] = {iCutPiTOF   , iCutPiTOF   , iCutPiTOF  , iCutPiTOF  , iCutPiTOF  , iCutPiTOF  , iCutPiTOF   ,  iCutPiTOF   , iCutPiTOF   ,  iCutPiTOF  , iCutPiTOF   ,  iCutPiTOF     };
   Int_t   cutID3   [12] = {iCutPTOF1 	, iCutPTOF1   , iCutPTOF1  , iCutPTOF1  , iCutPTOF1  , iCutPTOF1  , iCutPTOF1   ,  iCutPTOF1   , iCutPTOF1   ,  iCutPTOF1  , iCutPTOF1   ,  iCutPTOF1     };
   Int_t   cutID4   [12] = {iCutPiTOF1	, iCutPiTOF1  , iCutPiTOF1 , iCutPiTOF1 , iCutPiTOF1 , iCutPiTOF1 , iCutPiTOF1  ,  iCutPiTOF1  , iCutPiTOF1  ,  iCutPiTOF1 , iCutPiTOF1  ,  iCutPiTOF1    };
   Int_t   cutID5   [12] = {iCutPTOF2 	, iCutPTOF2   , iCutPTOF2  , iCutPTOF2  , iCutPTOF2  , iCutPTOF2  , iCutPTOF2   ,  iCutPTOF2   , iCutPTOF2   ,  iCutPTOF2  , iCutPTOF2   ,  iCutPTOF2     };
   Int_t   cutID6   [12] = {iCutPiTOF2	, iCutPiTOF2  , iCutPiTOF2 , iCutPiTOF2 , iCutPiTOF2 , iCutPiTOF2 , iCutPiTOF2  ,  iCutPiTOF2  , iCutPiTOF2  ,  iCutPiTOF2 , iCutPiTOF2  ,  iCutPiTOF2    };
   Int_t   cutID7   [12] = {iCutPTOF3 	, iCutPTOF3   , iCutPTOF3  , iCutPTOF3  , iCutPTOF3  , iCutPTOF3  , iCutPTOF3   ,  iCutPTOF3   , iCutPTOF3   ,  iCutPTOF3  , iCutPTOF3   ,  iCutPTOF3     };
   Int_t   cutID8   [12] = {iCutPiTOF3	, iCutPiTOF3  , iCutPiTOF3 , iCutPiTOF3 , iCutPiTOF3 , iCutPiTOF3 , iCutPiTOF3  ,  iCutPiTOF3  , iCutPiTOF3  ,  iCutPiTOF3 , iCutPiTOF3  ,  iCutPiTOF3    };
   Int_t   cutID9   [12] = {iCutPTOF4 	, iCutPTOF4   , iCutPTOF4  , iCutPTOF4  , iCutPTOF4  , iCutPTOF4  , iCutPTOF4   ,  iCutPTOF4   , iCutPTOF4   ,  iCutPTOF4  , iCutPTOF4   ,  iCutPTOF4     };
   Int_t   cutID10  [12] = {iCutPiTOF4	, iCutPiTOF4  , iCutPiTOF4 , iCutPiTOF4 , iCutPiTOF4 , iCutPiTOF4 , iCutPiTOF4  ,  iCutPiTOF4  , iCutPiTOF4  ,  iCutPiTOF4 , iCutPiTOF4  ,  iCutPiTOF4    };
   Int_t   cutID11  [12] = {iCutPTOF5 	, iCutPTOF5   , iCutPTOF5  , iCutPTOF5  , iCutPTOF5  , iCutPTOF5  , iCutPTOF5   ,  iCutPTOF5   , iCutPTOF5   ,  iCutPTOF5  , iCutPTOF5   ,  iCutPTOF5     };
   Int_t   cutID12  [12] = {iCutPiTOF5	, iCutPiTOF5  , iCutPiTOF5 , iCutPiTOF5 , iCutPiTOF5 , iCutPiTOF5 , iCutPiTOF5  ,  iCutPiTOF5  , iCutPiTOF5  ,  iCutPiTOF5 , iCutPiTOF5  ,  iCutPiTOF5    };
   Int_t   cutID13  [12] = {iCutPTOF6 	, iCutPTOF6   , iCutPTOF6  , iCutPTOF6  , iCutPTOF6  , iCutPTOF6  , iCutPTOF6   ,  iCutPTOF6   , iCutPTOF6   ,  iCutPTOF6  , iCutPTOF6   ,  iCutPTOF6     };
   Int_t   cutID14  [12] = {iCutPiTOF6	, iCutPiTOF6  , iCutPiTOF6 , iCutPiTOF6 , iCutPiTOF6 , iCutPiTOF6 , iCutPiTOF6  ,  iCutPiTOF6  , iCutPiTOF6  ,  iCutPiTOF6 , iCutPiTOF6  ,  iCutPiTOF6    };
   Int_t   ipdg     [12] = {2114        , -2114       , 2114       , -2114      , 2224       , -2224      , 2224        ,  -2224       , 2224        ,   -2224     , 2114        ,  -2114         };
   
   

   for (Int_t i = 0; i < 12; i++) {
      if (!use[i]) continue;
 
       AliRsnMiniOutput *outTOF = task->CreateOutput(Form("deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      outTOF->SetCutID(0, cutID1[i]);
      outTOF->SetCutID(1, cutID2[i]);
      outTOF->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF->SetCharge(0, charge1[i]);
      outTOF->SetCharge(1, charge2[i]);
      outTOF->SetMotherPDG(ipdg[i]);
      outTOF->SetMotherMass(1.232);
      // pair cuts
      outTOF->SetPairCuts(cutsPair);
      // axis X: invmass (or resolution)
      if (useIM[i])
         outTOF->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      outTOF->AddAxis(ptID, 100, 0.0, 10.0);
     // axis Z: Multiplicity
      outTOF->AddAxis( centID , 120 , 0.0, 120  ); 
      
      outTOF->AddAxis( fdpt   , 100 , 0.0, 10.0 ); 
      outTOF->AddAxis( sdpt   , 100 , 0.0, 10.0 ); 
      outTOF->AddAxis( fdp    , 100 , 0.0, 10.0 ); 
      outTOF->AddAxis( sdp    , 100 , 0.0, 10.0 ); 
      outTOF->AddAxis( etaID  , 20  ,-1.0, 1.0  );
      outTOF->AddAxis( yID    , 10  ,-0.5, 0.5 ); 
 
 
 
      AliRsnMiniOutput *outTOF1 = task->CreateOutput(Form("config1deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF1->SetCutID(0, cutID3[i]);
      outTOF1->SetCutID(1, cutID4[i]);
      outTOF1->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF1->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF1->SetCharge(0, charge1[i]);
      outTOF1->SetCharge(1, charge2[i]);
      outTOF1->SetMotherPDG(ipdg[i]);
      outTOF1->SetMotherMass(1.232);
      outTOF1->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF1->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF1->AddAxis(resID, 200, -0.02, 0.02);
      outTOF1->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF1->AddAxis(centID, 120, 0.0, 120); 


      AliRsnMiniOutput *outTOF2 = task->CreateOutput(Form("config2deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF2->SetCutID(0, cutID5[i]);
      outTOF2->SetCutID(1, cutID6[i]);
      outTOF2->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF2->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF2->SetCharge(0, charge1[i]);
      outTOF2->SetCharge(1, charge2[i]);
      outTOF2->SetMotherPDG(ipdg[i]);
      outTOF2->SetMotherMass(1.232);
      outTOF2->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF2->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF2->AddAxis(resID, 200, -0.02, 0.02);
      outTOF2->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF2->AddAxis(centID, 120, 0.0, 120);


      AliRsnMiniOutput *outTOF3 = task->CreateOutput(Form("config3deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF3->SetCutID(0, cutID7[i]);
      outTOF3->SetCutID(1, cutID8[i]);
      outTOF3->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF3->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF3->SetCharge(0, charge1[i]);
      outTOF3->SetCharge(1, charge2[i]);
      outTOF3->SetMotherPDG(ipdg[i]);
      outTOF3->SetMotherMass(1.232);
      outTOF3->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF3->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF3->AddAxis(resID, 200, -0.02, 0.02);
      outTOF3->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF3->AddAxis(centID, 120, 0.0, 120);


      AliRsnMiniOutput *outTOF4 = task->CreateOutput(Form("config4deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF4->SetCutID(0, cutID9[i]);
      outTOF4->SetCutID(1, cutID10[i]);
      outTOF4->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF4->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF4->SetCharge(0, charge1[i]);
      outTOF4->SetCharge(1, charge2[i]);
      outTOF4->SetMotherPDG(ipdg[i]);
      outTOF4->SetMotherMass(1.232);
      outTOF4->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF4->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF4->AddAxis(resID, 200, -0.02, 0.02);
      outTOF4->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF4->AddAxis(centID, 120, 0.0, 120);


      AliRsnMiniOutput *outTOF5 = task->CreateOutput(Form("config5deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF5->SetCutID(0, cutID11[i]);
      outTOF5->SetCutID(1, cutID12[i]);
      outTOF5->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF5->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF5->SetCharge(0, charge1[i]);
      outTOF5->SetCharge(1, charge2[i]);
      outTOF5->SetMotherPDG(ipdg[i]);
      outTOF5->SetMotherMass(1.232);
      outTOF5->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF5->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF5->AddAxis(resID, 200, -0.02, 0.02);
      outTOF5->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF5->AddAxis(centID, 120, 0.0, 120);

      AliRsnMiniOutput *outTOF6 = task->CreateOutput(Form("config6deltatof_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      outTOF6->SetCutID(0, cutID13[i]);
      outTOF6->SetCutID(1, cutID14[i]);
      outTOF6->SetDaughter(0, AliRsnDaughter::kProton);
      outTOF6->SetDaughter(1, AliRsnDaughter::kPion);
      outTOF6->SetCharge(0, charge1[i]);
      outTOF6->SetCharge(1, charge2[i]);
      outTOF6->SetMotherPDG(ipdg[i]);
      outTOF6->SetMotherMass(1.232);
      outTOF6->SetPairCuts(cutsPair);
      if (useIM[i])
         outTOF6->AddAxis(imID, 100, 1.0, 2.0);
      else
         outTOF6->AddAxis(resID, 200, -0.02, 0.02);
      outTOF6->AddAxis(ptID, 100, 0.0, 10.0);
      outTOF6->AddAxis(centID, 120, 0.0, 120);

   }
   
   return kTRUE;
}
