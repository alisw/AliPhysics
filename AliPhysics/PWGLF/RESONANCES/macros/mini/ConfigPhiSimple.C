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
Bool_t ConfigPhiSimple
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
   
   // integrated kaon cut
   AliRsnCutKaonForPhi2010PP *cutStd = new AliRsnCutKaonForPhi2010PP("cutStdPP");
   // cut set
   AliRsnCutSet *cutSetStd = new AliRsnCutSet("kaonForPhi", AliRsnTarget::kDaughter);
   cutSetStd->AddCut(cutStd);
   cutSetStd->SetCutScheme(cutStd->GetName());
   // add to task
   Int_t icut = task->AddTrackCuts(cutSetStd);
   
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
   
   // create output object                          "HIST"   "PAIR" or "MIX" or "TRUE"
   AliRsnMiniOutput *out = task->CreateOutput(name, outType, computationType);
      
   // cut IDs
   out->SetCutID(0, icut);
   out->SetCutID(1, icut);
   
   // daughter species   
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kKaon);
   out->SetCharge(0, charge1);
   out->SetCharge(1, charge2);
   
   // resonance properties
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
   
   return kTRUE;
}
