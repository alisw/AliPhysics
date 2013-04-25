#ifndef __CINT__
#include <AliRsnListOutput.h>
#include <AliRsnValuePair.h>
#include <AliRsnValueEvent.h>
#include <AliRsnMiniAnalysisTask.h>
#endif

void AddRsnPairsPhi(AliAnalysisTaskSE *task,
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
   Int_t         pdg  = 333;
   TDatabasePDG *db   = TDatabasePDG::Instance();
   TParticlePDG *part = db->GetParticle(pdg);
   Double_t mass = part->Mass();
   Bool_t valid;

   Int_t isRsnMini = AliRsnTrainManager::GetGlobalInt("IsRsnMini",valid);

   if (isRsnMini) {
      AddPairOutputMiniPhi(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsPair,suffix);
   } else {
      // this function is common and it is located in RsnConfig.C
      // as ouptup AddPairOutputPhi from this macro will be taken
      AddPair(task,isMC,isMixing,pType1,listID1,pType2,listID2,pdg,mass,cutsEvent,cutsPair,suffix);
   }
}
void AddPairOutputPhi(AliRsnLoopPair *pair)
{
   Bool_t valid;
   Int_t isFullOutput = AliRsnTrainManager::GetGlobalInt("RsnOutputFull",valid);
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);
   Int_t useRapidity = AliRsnTrainManager::GetGlobalInt("RsnUseRapidity",valid);

   // axes
   AliRsnValuePair *axisIM = new AliRsnValuePair("IM", AliRsnValuePair::kInvMass);
   AliRsnValuePair *axisPt = new AliRsnValuePair("PT", AliRsnValuePair::kPt);
   AliRsnValuePair *axisEta = new AliRsnValuePair("ETA", AliRsnValuePair::kEta);
   AliRsnValuePair *axisRapidity = new AliRsnValuePair("Y", AliRsnValuePair::kY);

   AliRsnValueEvent *axisCentrality = 0;
   if (collisionType == 1) axisCentrality = new AliRsnValueEvent("MULTI",AliRsnValueEvent::kCentralityV0);


   axisIM     ->SetBins(300, 0.9, 1.2);
//   axisIM     ->SetBins(600, 0.9, 1.2);
//   axisIM     ->SetBins(1000, 0.9, 1.9);
   axisPt     ->SetBins(120, 0.0, 12.0);
//    axisEta    ->SetBins(400, -2.0, 2.0);
   axisEta    ->SetBins(400, -0.5, 0.5);
   axisRapidity    ->SetBins(10, -0.5, 0.5);

   if (axisCentrality) axisCentrality->SetBins(20,0,100);

   // output: 2D histogram of inv. mass vs. pt
   AliRsnListOutput *outPair = 0;
   if (!isFullOutput) {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoDefault);
      outPair->AddValue(axisIM);
   } else {
      outPair = new AliRsnListOutput("pair", AliRsnListOutput::kHistoSparse);
      outPair->AddValue(axisIM);
      outPair->AddValue(axisPt);
      if (useRapidity) outPair->AddValue(axisRapidity);
      else outPair->AddValue(axisEta);
      if (axisCentrality) outPair->AddValue(axisCentrality);
   }
   // add outputs to loop
   pair->AddOutput(outPair);
}

void AddPairOutputMiniPhi(AliAnalysisTaskSE *task, Bool_t isMC,Bool_t isMixing, AliPID::EParticleType pType1,Int_t listID1, AliPID::EParticleType pType2,Int_t listID2, Int_t pdgMother,Double_t massMother, AliRsnCutSet *cutsPair=0,TString suffix = "") {

   Bool_t valid;
   Int_t isFullOutput = AliRsnTrainManager::GetGlobalInt("RsnOutputFull",valid);
   Int_t useMixing = AliRsnTrainManager::GetGlobalInt("IsMixing",valid);
   Int_t collisionType = AliRsnTrainManager::GetGlobalInt("IsCollisionType",valid);

   Int_t useMixLike = AliRsnTrainManager::GetGlobalInt("RsnMixLike",valid);
   Int_t useRapidity = AliRsnTrainManager::GetGlobalInt("RsnUseRapidity",valid);

   AliRsnMiniAnalysisTask *taskRsnMini =  (AliRsnMiniAnalysisTask *)task;

   if (collisionType==0) taskRsnMini->UseMultiplicity("QUALITY");
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
   /* IM resolution    */ Int_t resID  = taskRsnMini->CreateValue(AliRsnMiniValue::kInvMassDiff, kTRUE);
   /* transv. momentum */ Int_t ptID   = taskRsnMini->CreateValue(AliRsnMiniValue::kPt, kFALSE);
   /* centrality       */ Int_t centID = taskRsnMini->CreateValue(AliRsnMiniValue::kMult, kFALSE);
   /* eta              */ Int_t etaID = taskRsnMini->CreateValue(AliRsnMiniValue::kEta, kFALSE);
   /* rapidity         */ Int_t yID = taskRsnMini->CreateValue(AliRsnMiniValue::kY, kFALSE);


   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --

//   Bool_t useMixLike = kFALSE;
//   useMixLike = kTRUE;

   const Int_t numAxis = 7;
   Bool_t  use     [numAxis] = { 1      ,  useMixing      ,  1      ,  1      ,  isMC  ,  useMixLike,  useMixLike};
   TString name    [numAxis] = {"Unlike", "Mixing", "LikePP", "LikeMM", "Trues", "MixingPP", "MixingMM"};
   TString comp    [numAxis] = {"PAIR"  , "MIX"   , "PAIR"  , "PAIR"  , "TRUE" , "MIX"   , "MIX"   };
   Char_t  charge1 [numAxis] = {'+'     , '+'     , '+'     , '-'     , '+'    , '+'     , '-'     };
   Char_t  charge2 [numAxis] = {'-'     , '-'     , '+'     , '-'     , '-'    , '+'     , '-'     };

   // common definitions
   TString outputType = "HIST";
   if (isFullOutput) outputType = "SPARSE";

   Int_t nIM   = 215; Double_t minIM   = 0.985, maxIM =  1.2;
//    Int_t nEta   = 400; Double_t minEta   = -2.0, maxEta =  2.0;
   Int_t nEta   = 400; Double_t minEta   = -0.5, maxEta =  0.5;
   Int_t nY   = 1; Double_t minY   = -0.5, maxY =  0.5;
//   Int_t nIM   = 1000; Double_t minIM   = 0.9, maxIM =  1.9;
   Int_t nPt   = 240; Double_t minPt   = 0.0, maxPt = 24.0;
   Int_t nCent = 20; Double_t minCent = 0.0, maxCent = 100.0;
   Int_t nRes  = 200; Double_t maxRes  = 0.01;

   // in case pp
   if (collisionType==0) {
      nIM   = 150; minIM  = 0.975; maxIM  =  1.125;
      nPt   = 160; minPt  = 0.0,   maxPt  = 8.0;
      nEta  = 30;  minEta = -1.5;  maxEta =  1.5;
      nY    = 3;   minY   = -1.5;  maxY   =  1.5;
   }


   // retrieve mass from PDG database
   Int_t         pdg  = 333;
   TDatabasePDG *db   = TDatabasePDG::Instance();
   TParticlePDG *part = db->GetParticle(pdg);


   Printf(suffix.Data());
   // create standard outputs
   for (Int_t i = 0; i < numAxis; i++) {
      if (!use[i]) continue;
      // create output
      AliRsnMiniOutput *out = taskRsnMini->CreateOutput(Form("%s_%s", suffix.Data(),name[i].Data() ), outputType.Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, listID1);
      out->SetCutID(1, listID1);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kKaon);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(pdg);
      out->SetMotherMass(part->Mass());
      // pair cuts
      if (cutsPair) out->SetPairCuts(cutsPair);
      // axis X: invmass
      out->AddAxis(imID, nIM, minIM, maxIM);

      if (isFullOutput) {
         // axis Y: transverse momentum
         out->AddAxis(ptID, nPt, minPt, maxPt);

         if (useRapidity) out->AddAxis(yID, nY, minY, maxY);
         if (collisionType==0) out->AddAxis(etaID, nEta, minEta, maxEta);
         // axis Z: centrality
         if (collisionType==1) out->AddAxis(centID, nCent, minCent, maxCent);
      }
   }

   // add output for resolution
   if (isMC) {
      AliRsnMiniOutput *outRes = taskRsnMini->CreateOutput(Form("phi_Res%s", suffix.Data()), outputType.Data(), "TRUE");
      // selection settings
      outRes->SetCutID(0, listID1);
      outRes->SetCutID(1, listID1);
      outRes->SetDaughter(0, AliRsnDaughter::kKaon);
      outRes->SetDaughter(1, AliRsnDaughter::kKaon);
      outRes->SetCharge(0, '+');
      outRes->SetCharge(1, '-');
      outRes->SetMotherPDG(pdg);
      outRes->SetMotherMass(part->Mass());
      // pair cuts
      if (cutsPair) outRes->SetPairCuts(cutsPair);
      // axis X: resolution
      outRes->AddAxis(resID, nRes, -maxRes, maxRes);

      if (isFullOutput) {
         // axis Y: transverse momentum
         outRes->AddAxis(ptID, nPt, minPt, maxPt);
         if (useRapidity) outRes->AddAxis(yID, nY, minY, maxY);
         if (collisionType==0) outRes->AddAxis(etaID, nEta, minEta, maxEta);

         // axis Z: centrality
         if (collisionType==1) outRes->AddAxis(centID, nCent, minCent, maxCent);
      }
   }

   //
   // -- Create output for MC generated ------------------------------------------------------------
   //

   if (isMC) {
      // create ouput
      AliRsnMiniOutput *outMC = taskRsnMini->CreateOutput(Form("phi_MCGen%s", suffix.Data()), outputType.Data(), "MOTHER");
      // selection settings
      outMC->SetDaughter(0, AliRsnDaughter::kKaon);
      outMC->SetDaughter(1, AliRsnDaughter::kKaon);
      outMC->SetMotherPDG(pdg);
      outMC->SetMotherMass(part->Mass());
      // pair cuts
      if (cutsPair) outMC->SetPairCuts(cutsPair);
      // axis X: invmass
      outMC->AddAxis(imID, nIM, minIM, maxIM);
      if (isFullOutput) {
         // axis Y: transverse momentum
         outMC->AddAxis(ptID, nPt, minPt, maxPt);
         if (useRapidity) outMC->AddAxis(yID, nY, minY, maxY);
         if (collisionType==0) outMC->AddAxis(etaID, nEta, minEta, maxEta);
         // axis Z: centrality
         if (collisionType==1) outMC->AddAxis(centID, nCent, minCent, maxCent);
      }
   }


}



