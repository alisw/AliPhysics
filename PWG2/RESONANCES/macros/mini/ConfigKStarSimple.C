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
Bool_t ConfigKStarSimple
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isMC, 
   
   const char             *name,
   const char             *outType,
   const char             *computationType,
   char                    charge1,
   char                    charge2,
   Bool_t                  useIM,
   
   AliRsnCutSet           *cutsPair
)
{
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   // integrated pion cut
   AliRsnCutDaughterKStar2010PP *cutPi = new AliRsnCutDaughterKStar2010PP("cutPionForKStar", AliPID::kPion);
   // cut set
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionForKStar", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   // add to task
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   
   // integrated kaon cut
   AliRsnCutDaughterKStar2010PP *cutK = new AliRsnCutDaughterKStar2010PP("cutKaonForKStar", AliPID::kKaon);
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("setKaonForKStar", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());
   // add to task
   Int_t iCutK = task->AddTrackCuts(cutSetK);
   
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
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput(name, outType, computationType);
   
   // selection settings
   out->SetCutID(0, iCutK);
   out->SetCutID(1, iCutPi);
   
   // daughter species
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, charge1);
   out->SetCharge(1, charge2);
   
   // resonance properties
   out->SetMotherPDG(313);
   out->SetMotherMass(0.896);
   
   // pair cuts
   out->SetPairCuts(cutsPair);
   
   // axis X: invmass (or resolution)
   if (useIM) 
      out->AddAxis(imID, 90, 0.6, 1.5);
   else
      out->AddAxis(resID, 200, -0.02, 0.02);
      
   // axis Y: transverse momentum
   out->AddAxis(ptID, 100, 0.0, 10.0);
}
