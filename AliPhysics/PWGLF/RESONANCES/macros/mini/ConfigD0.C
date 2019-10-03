//
// *** Configuration script for D0 analysis ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t ConfigD0
(  
   AliRsnMiniAnalysisTask *task, 
   Bool_t                  isPP,
   Bool_t                  isMC,
   Bool_t                  monitor = kTRUE,
   Bool_t      		   centortracklets = kTRUE, 
   Bool_t      		   sanityhistos = kTRUE, 
   Float_t         	   nsigmaTPCPi = 3.0,
   Float_t         	   nsigmaTPCKa = 3.0,
   Float_t         	   nsigmaTOFPi = 2.0,
   Float_t         	   nsigmaTOFKa = 2.0,
   Int_t                   aodFilterBit = 5,
   Float_t         	   trackDCAcutMax = 7.0,
   Float_t         	   trackDCAcutMin = 0.0,
   Float_t		   trackDCAZcutMax = 2.0,
   Int_t           	   NTPCcluster = 70,
   Double_t                NTPCcrratio = 0.8,
   Int_t		   minSPDclt = 0,
   Double_t                minpt = 0.15,
   Short_t     		   maxSisters = 2,
   Bool_t                  checkP = kTRUE,
   Bool_t      		   useNTPCclt = kTRUE,
   Bool_t                  minDCAcutFixed = kFALSE,
   Bool_t                  maxDCAcutFixed = kFALSE,
   Bool_t                  ptdepPIDcut = kFALSE,
   Bool_t      		   checkFeedDown = kTRUE,
   Bool_t      		   checkQuark = kTRUE,
   Bool_t      		   doCalculationInMC = kTRUE,
   UShort_t    		   originDselection = 0,
   Float_t                 mineta = -0.8,
   Float_t                 maxeta = 0.8,
   Float_t		   min_inv_mass = 0.6,
   Float_t		   max_inv_mass = 2.2,
   Int_t                   bins = 320,
   const char      	  *suffix,
   AliRsnCutSet           *cutsPairY,
   AliRsnCutSet           *cutsPair
)
{
   // manage suffix
   if (strlen(suffix) > 0) suffix = Form("_%s", suffix);
   
   TString s = ""; s+=trackDCAcutMax; s+="*(0.0015+0.0050/pt^1.01)";

   const char *formula = s;
   
   TString s2 = ""; s2+=trackDCAcutMin; s2+="*(0.0015+0.0050/pt^1.01)";

   const char *formulaMin = s2;
   
   
   // 
   // -- Define track cuts -------------------------------------------------------------------------
   //
   
   // integrated pion cut
   AliRsnCutDaughterD0 *cutPi = new AliRsnCutDaughterD0("cutPionForD0", AliPID::kPion);
   //cutPi->SetNoPID(kTRUE);
   cutPi->SetTPCPionPIDCut(nsigmaTPCPi);
   cutPi->SetTOFPionPIDCut(nsigmaTOFPi);
   cutPi->SetPtDependentPIDCut(ptdepPIDcut);
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutPi->CutQuality();
   cutQuality->SetCheckOnlyFilterBit(kFALSE);
   cutQuality->SetAODTestFilterBit(aodFilterBit);
   if(maxDCAcutFixed)cutQuality->SetDCARmax(trackDCAcutMax);	         
   else cutQuality->SetDCARPtFormula(formula);
   if(minDCAcutFixed) cutQuality->SetDCARmin(trackDCAcutMin);
   else cutQuality->SetDCARPtFormulaMin(formulaMin);
   if(useNTPCclt)cutQuality->SetTPCminNClusters(NTPCcluster);
   if(!useNTPCclt)cutQuality->SetMinNCrossedRowsTPC(NTPCcluster,kTRUE);
   if(!useNTPCclt)cutQuality->SetMinNCrossedRowsOverFindableClsTPC(NTPCcrratio,kTRUE);
   cutQuality->SetPtRange(minpt,1E20);
   cutQuality->SetEtaRange(mineta, maxeta);
   cutQuality->SetDCAZmax(trackDCAZcutMax);
   cutQuality->SetSPDminNClusters(minSPDclt);
   cutQuality->SetITSminNClusters(0);
   cutQuality->SetITSmaxChi2(1E20);
   cutQuality->SetTPCmaxChi2(1E20);
   cutQuality->SetRejectKinkDaughters();
   cutQuality->Print();
  	         

   
   // cut set
   AliRsnCutSet *cutSetPi = new AliRsnCutSet("setPionD0", AliRsnTarget::kDaughter);
   cutSetPi->AddCut(cutPi);
   cutSetPi->SetCutScheme(cutPi->GetName());
   // add to task
   Int_t iCutPi = task->AddTrackCuts(cutSetPi);
   
   
   
   // integrated kaon cut
   AliRsnCutDaughterD0 *cutK = new AliRsnCutDaughterD0("cutKaonForD0", AliPID::kKaon);
   //cutK->SetNoPID(kTRUE);
   cutK->SetTPCKaonPIDCut(nsigmaTPCKa);
   cutK->SetTOFKaonPIDCut(nsigmaTOFKa);
   cutK->SetPtDependentPIDCut(ptdepPIDcut);	
   AliRsnCutTrackQuality *cutQuality = (AliRsnCutTrackQuality*) cutK->CutQuality();
   cutQuality->SetCheckOnlyFilterBit(kFALSE);
   cutQuality->SetAODTestFilterBit(aodFilterBit);
   if(maxDCAcutFixed)cutQuality->SetDCARmax(trackDCAcutMax);	         
   else cutQuality->SetDCARPtFormula(formula);
   if(minDCAcutFixed) cutQuality->SetDCARmin(trackDCAcutMin);
   else cutQuality->SetDCARPtFormulaMin(formulaMin);
   if(useNTPCclt)cutQuality->SetTPCminNClusters(NTPCcluster);
   if(!useNTPCclt)cutQuality->SetMinNCrossedRowsTPC(NTPCcluster,kTRUE);
   if(!useNTPCclt)cutQuality->SetMinNCrossedRowsOverFindableClsTPC(NTPCcrratio,kTRUE);
   cutQuality->SetPtRange(minpt,1E20);
   cutQuality->SetEtaRange(mineta, maxeta);
   cutQuality->SetDCAZmax(trackDCAZcutMax);
   cutQuality->SetSPDminNClusters(minSPDclt);
   cutQuality->SetITSminNClusters(0);
   cutQuality->SetITSmaxChi2(1E20);
   cutQuality->SetTPCmaxChi2(1E20);
   cutQuality->SetRejectKinkDaughters();
   cutQuality->Print();
	
   
   
   // cut set
   AliRsnCutSet *cutSetK = new AliRsnCutSet("setKaonD0", AliRsnTarget::kDaughter);
   cutSetK->AddCut(cutK);
   cutSetK->SetCutScheme(cutK->GetName());
   // add to task
   Int_t iCutK = task->AddTrackCuts(cutSetK);
   
   
   
   ////////////////////////////////////////////////////////////////////////////////////////////////
  
  // -- Values ------------------------------------------------------------------------------------
  /* invariant mass     */ Int_t imID       = task->CreateValue(AliRsnMiniValue::kInvMass, kFALSE);
  /* IM resolution      */ Int_t resID      = task->CreateValue(AliRsnMiniValue::kInvMassRes, kTRUE);
  /* transv. momentum   */ Int_t ptID       = task->CreateValue(AliRsnMiniValue::kPt, kFALSE);
  /* centrality         */ Int_t centID     = task->CreateValue(AliRsnMiniValue::kMult, kFALSE);
  /* pseudorapidity     */ Int_t etaID      = task->CreateValue(AliRsnMiniValue::kEta, kFALSE);
  /* rapidity           */ Int_t yID        = task->CreateValue(AliRsnMiniValue::kY, kFALSE);
  /* dca product        */ Int_t dcapID     = task->CreateValue(AliRsnMiniValue::kDCAproduct, kFALSE);
  /* first daughter pt  */ Int_t daug1ptID  = task->CreateValue(AliRsnMiniValue::kFirstDaughterPt, kFALSE);
  /* second daughter pt */ Int_t daug2ptID  = task->CreateValue(AliRsnMiniValue::kSecondDaughterPt, kFALSE);
  /* first daughter dca */ Int_t daug1dcaID = task->CreateValue(AliRsnMiniValue::kFirstDaughterDCA, kFALSE);
  /* second daughter dca*/ Int_t daug2dcaID = task->CreateValue(AliRsnMiniValue::kSecondDaughterDCA, kFALSE);
  /* number of Sisters  */ Int_t nsistID    = task->CreateValue(AliRsnMiniValue::kNSisters, kFALSE);
  /* number of trackl.  */ Int_t trackletsID = task->CreateValue(AliRsnMiniValue::kTracklets, kFALSE);
   
   //
   // -- Create all needed outputs -----------------------------------------------------------------
   //
   
   // use an array for more compact writing, which are different on mixing and charges
   // [0] = unlike
   // [1] = mixing
   // [2] = like ++
   // [3] = like --
   
   
   if(!isMC || doCalculationInMC == kTRUE){
   
   Bool_t   use     [8] = { 1	    ,  1       ,  1	  ,  1	     ,  1        ,  1	     ,  1       ,  1	  };
   Bool_t   useIM   [8] = { 1	    ,  1       ,  1	  ,  1	     ,  1        ,  1	     ,  1       ,  1	  };
   TString  name    [8] = {"Unlike1", "Unlike2", "Mixing1", "Mixing2", "RotateK1", "RotateK2", "LikePP" , "LikeMM"};
   TString  comp    [8] = {"PAIR"   , "PAIR"   , "MIX"	  , "MIX"    , "ROTATE1" , "ROTATE1" , "PAIR"   , "PAIR"  };
   TString  output  [8] = {"SPARSE" , "SPARSE" , "SPARSE" , "SPARSE" , "SPARSE"  , "SPARSE"  , "SPARSE" , "SPARSE"};
   Char_t   charge1 [8] = {'-'	    , '+'      , '-'	  , '+'      , '-'       , '+'	     , '+'      , '-'	  };
   Char_t   charge2 [8] = {'+'	    , '-'      , '+'	  , '-'      , '+'       , '-'	     , '+'      , '-'	  };
   Int_t    cutID1  [8] = { iCutK   ,  iCutK   ,  iCutK   ,  iCutK   ,  iCutK    ,  iCutK    ,  iCutK   ,  iCutK  };
   Int_t    cutID2  [8] = { iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi  ,  iCutPi   ,  iCutPi   ,  iCutPi  ,  iCutPi };
   Int_t    ipdg    [8] = { 421     , -421     ,  421	  , -421     ,  421      , -421      ,  421     , -421	  };
   Double_t mass    [8] = { 1.86486 ,  1.86486 ,  1.86486 ,  1.86486 ,  1.86486  ,  1.86486  ,  1.86486 ,  1.86486};
   
   for (Int_t i = 0; i < 8; i++) {
      if (!use[i]) continue;
      
      // create output
      AliRsnMiniOutput *out = task->CreateOutput(Form("D0_%s%s", name[i].Data(), suffix), output[i].Data(), comp[i].Data());
      // selection settings
      out->SetCutID(0, cutID1[i]);
      out->SetCutID(1, cutID2[i]);
      out->SetDaughter(0, AliRsnDaughter::kKaon);
      out->SetDaughter(1, AliRsnDaughter::kPion);
      out->SetCharge(0, charge1[i]);
      out->SetCharge(1, charge2[i]);
      out->SetMotherPDG(ipdg[i]);
      out->SetMotherMass(mass[i]);
      // pair cuts
      out->SetPairCuts(cutsPair);

      // axis X: invmass (or resolution)
      //if (useIM[i]) 
         out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
      //else
      //   out->AddAxis(resID, 200, -0.02, 0.02);
      // axis Y: transverse momentum
      out->AddAxis(ptID, 200, 0.0, 20.0);
      
      // axiz Z: pseudorapidity
      //if(sanityhistos==kTRUE) out->AddAxis(etaID, 100, -1, 1);
         
      
      if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		   else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
      else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   }

   if(sanityhistos==kTRUE){

   // SANITY CHECK HISTOS
   
   TString mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_MixingPP", mode.Data(), "MIX");
   // selection settings
   out->SetCutID(0, iCutK);
   out->SetCutID(1, iCutPi);
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, '+');
   out->SetCharge(1, '+');
   out->SetMotherPDG(421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPair);
  
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(etaID, 100, -1, 1);


   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_MixingMM", mode.Data(), "MIX");
   // selection settings
   out->SetCutID(0, iCutK);
   out->SetCutID(1, iCutPi);
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, '-');
   out->SetCharge(1, '-');
   out->SetMotherPDG(-421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPair);
   
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(etaID, 100, -1, 1);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}


   }

   
   }
   
   if(monitor == kTRUE){
   AddMonitorOutput_PionEta(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionY(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionMinPt(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionDCA(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionTPC_PIDCut(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionTOF_PIDCut(cutSetPi->GetMonitorOutput());
   AddMonitorOutput_PionNTPC(cutSetPi->GetMonitorOutput());
   
   AddMonitorOutput_KaonEta(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonY(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonMinPt(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonDCA(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonTPC_PIDCut(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonTOF_PIDCut(cutSetK->GetMonitorOutput());
   AddMonitorOutput_KaonNTPC(cutSetK->GetMonitorOutput());
   }
   
   
   if (isMC) {
   
   // TRUE RECONSTRUCTED PAIRS
   
   TString mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_True1", mode.Data(), "TRUE");
   // selection settings
   out->SetCutID(0, iCutK);
   out->SetCutID(1, iCutPi);
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, '-');
   out->SetCharge(1, '+');
   out->SetMotherPDG(421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPair);
   out->SetMaxNSisters(maxSisters);
   out->SetCheckMomentumConservation(checkP);
   out->SetCheckFeedDown(checkFeedDown);
   out->SetRejectCandidateIfNotFromQuark(checkQuark);
   out->SetDselection(originDselection);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);
   //out->AddAxis(dcapID, 100, -0.001, 0.001);
   //out->AddAxis(nsistID, 10, 0, 5);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_True2", mode.Data(), "TRUE");
   // selection settings
   out->SetCutID(0, iCutK);
   out->SetCutID(1, iCutPi);
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetCharge(0, '+');
   out->SetCharge(1, '-');
   out->SetMotherPDG(-421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPair);
   out->SetMaxNSisters(maxSisters);
   out->SetCheckMomentumConservation(checkP);
   out->SetCheckFeedDown(checkFeedDown);
   out->SetRejectCandidateIfNotFromQuark(checkQuark);
   out->SetDselection(originDselection);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);
   //out->AddAxis(dcapID, 100, -0.001, 0.001);
   //out->AddAxis(nsistID, 10, 0, 5);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   
   // GENERATED MOTHERS
   
   TString mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_TrueMC1", mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPairY);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_TrueMC2", mode.Data(), "MOTHER");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(-421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPairY);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   
   // GENERATED MOTHERS IN THE ACCEPTANCE
   
   TString mode = "SPARSE";
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_TrueMC_InAcc1", mode.Data(), "MOTHER_IN_ACC");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPairY);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);

  if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       	       else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   // create output
   AliRsnMiniOutput *out = task->CreateOutput("D0_TrueMC_InAcc2", mode.Data(), "MOTHER_IN_ACC");
   // selection settings
   out->SetDaughter(0, AliRsnDaughter::kKaon);
   out->SetDaughter(1, AliRsnDaughter::kPion);
   out->SetMotherPDG(-421);
   out->SetMotherMass(1.86486);
   // pair cuts
   out->SetPairCuts(cutsPairY);
   // binnings
   out->AddAxis(imID, bins, min_inv_mass, max_inv_mass);
   out->AddAxis(ptID, 200, 0.0, 20.0);
   //out->AddAxis(yID, 100, -1, 1);

   if (!isPP)  {if(!centortracklets) out->AddAxis(centID, 100, 0.0, 100.0);
       		else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   else {if(!centortracklets) out->AddAxis(centID, 400, 0.0, 400.0);
	else out->AddAxis(trackletsID, 250, -0.5, 249.5);}
   
   
   }


   return kTRUE;
}


void AddMonitorOutput_PionEta(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *peta=0)
{

   // PionEta
   AliRsnValueDaughter *axisPionEta = new AliRsnValueDaughter("pion_eta", AliRsnValueDaughter::kEta);
   axisPionEta->SetBins(-1.0,1.0,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionEta = new AliRsnListOutput("Pion_Eta", AliRsnListOutput::kHistoDefault);
   outMonitorPionEta->AddValue(axisPionEta);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionEta);
   if (peta) peta->AddOutput(outMonitorPionEta);
  
}

void AddMonitorOutput_PionY(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *py=0)
{

   // PionY
   AliRsnValueDaughter *axisPionY = new AliRsnValueDaughter("pion_y", AliRsnValueDaughter::kY);
   axisPionY->SetBins(-1.0,1.0,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionY = new AliRsnListOutput("Pion_Y", AliRsnListOutput::kHistoDefault);
   outMonitorPionY->AddValue(axisPionY);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionY);
   if (py) py->AddOutput(outMonitorPionY);
  
}

void AddMonitorOutput_PionMinPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pmpt=0)
{

   // PionMinPt
   AliRsnValueDaughter *axisPionMinPt = new AliRsnValueDaughter("pion_minpt", AliRsnValueDaughter::kPt);
   axisPionMinPt->SetBins(0.0,1,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionMinPt = new AliRsnListOutput("Pion_MinPt", AliRsnListOutput::kHistoDefault);
   outMonitorPionMinPt->AddValue(axisPionMinPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionMinPt);
   if (pmpt) pmpt->AddOutput(outMonitorPionMinPt);
  
}

void AddMonitorOutput_PionDCA(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *pdca=0)
{

   // PionDCA
   AliRsnValueDaughter *axisPionDCA = new AliRsnValueDaughter("pion_dca", AliRsnValueDaughter::kDCAXY);
   axisPionDCA->SetBins(-1.0,1,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionDCA = new AliRsnListOutput("Pion_DCA", AliRsnListOutput::kHistoDefault);
   outMonitorPionDCA->AddValue(axisPionDCA);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionDCA);
   if (pdca) pdca->AddOutput(outMonitorPionDCA);
  
}

void AddMonitorOutput_PionTPC_PIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piTPCPID=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionTPCPIDCut = new AliRsnValueDaughter("pionTPCPID", AliRsnValueDaughter::kTPCnsigmaPi);
   axisPionTPCPIDCut->SetBins(0.0,5,0.1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionTPCPIDCut = new AliRsnListOutput("Pion_TPC_PID_Cut", AliRsnListOutput::kHistoDefault);
   outMonitorPionTPCPIDCut->AddValue(axisPionTPCPIDCut);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionTPCPIDCut);
   if (piTPCPID) piTPCPID->AddOutput(outMonitorPionPIDCut);
  
}

void AddMonitorOutput_PionTOF_PIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piTOFPID=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionTOFPIDCut = new AliRsnValueDaughter("pionTOFPID", AliRsnValueDaughter::kTOFnsigmaPi);
   axisPionTOFPIDCut->SetBins(0.0,5,0.1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionTOFPIDCut = new AliRsnListOutput("Pion_TOF_PID_Cut", AliRsnListOutput::kHistoDefault);
   outMonitorPionTOFPIDCut->AddValue(axisPionTOFPIDCut);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionTOFPIDCut);
   if (piTOFPID) piTOFPID->AddOutput(outMonitorPionTOFPIDCut);
  
}

void AddMonitorOutput_PionNTPC(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *piNTPC=0)
{

   // Pion PID Cut
   AliRsnValueDaughter *axisPionNTPC = new AliRsnValueDaughter("pionNTPC", AliRsnValueDaughter::kNTPCclusters);
   axisPionNTPC->SetBins(0.0,200,1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorPionNTPC = new AliRsnListOutput("Pion_NTPC", AliRsnListOutput::kHistoDefault);
   outMonitorPionNTPC->AddValue(axisPionNTPC);

   // add outputs to loop
   if (mon) mon->Add(outMonitorPionNTPC);
   if (piNTPC) piNTPC->AddOutput(outMonitorPionNTPC);
  
}

void AddMonitorOutput_KaonEta(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *keta=0)
{

   // KaonEta
   AliRsnValueDaughter *axisKaonEta = new AliRsnValueDaughter("kaon_eta", AliRsnValueDaughter::kEta);
   axisKaonEta->SetBins(-1.0,1.0,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonEta = new AliRsnListOutput("Kaon_Eta", AliRsnListOutput::kHistoDefault);
   outMonitorKaonEta->AddValue(axisKaonEta);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonEta);
   if (keta) keta->AddOutput(outMonitorKaonEta);
  
}

void AddMonitorOutput_KaonY(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *ky=0)
{

   // KaonY
   AliRsnValueDaughter *axisKaonY = new AliRsnValueDaughter("kaon_y", AliRsnValueDaughter::kY);
   axisKaonY->SetBins(-1.0,1.0,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonY = new AliRsnListOutput("Kaon_Y", AliRsnListOutput::kHistoDefault);
   outMonitorKaonY->AddValue(axisKaonY);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonY);
   if (ky) ky->AddOutput(outMonitorKaonY);
  
}

void AddMonitorOutput_KaonMinPt(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *kmpt=0)
{

   // KaonMinPt
   AliRsnValueDaughter *axisKaonMinPt = new AliRsnValueDaughter("kaon_minpt", AliRsnValueDaughter::kPt);
   axisKaonMinPt->SetBins(0.0,1,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonMinPt = new AliRsnListOutput("Kaon_MinPt", AliRsnListOutput::kHistoDefault);
   outMonitorKaonMinPt->AddValue(axisKaonMinPt);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonMinPt);
   if (kmpt) kmpt->AddOutput(outMonitorKaonMinPt);
  
}

void AddMonitorOutput_KaonDCA(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *kdca=0)
{

   // KaonDCA
   AliRsnValueDaughter *axisKaonDCA = new AliRsnValueDaughter("kaon_dca", AliRsnValueDaughter::kDCAXY);
   axisKaonDCA->SetBins(-1.0,1,0.05);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonDCA = new AliRsnListOutput("Kaon_DCA", AliRsnListOutput::kHistoDefault);
   outMonitorKaonDCA->AddValue(axisKaonDCA);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonDCA);
   if (kdca) kdca->AddOutput(outMonitorKaonDCA);
  
}

void AddMonitorOutput_KaonTPC_PIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *kTPCPID=0)
{

   // Kaon TPC PID Cut
   AliRsnValueDaughter *axisKaonTPCPIDCut = new AliRsnValueDaughter("kaonTPCPID", AliRsnValueDaughter::kTPCnsigmaK);
   axisKaonTPCPIDCut->SetBins(0.0,5,0.1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonTPCPIDCut = new AliRsnListOutput("Kaon_TPC_PID_Cut", AliRsnListOutput::kHistoDefault);
   outMonitorKaonTPCPIDCut->AddValue(axisKaonTPCPIDCut);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonTPCPIDCut);
   if (kTPCPID) kTPCPID->AddOutput(outMonitorKaonTPCPIDCut);
  
}

void AddMonitorOutput_KaonTOF_PIDCut(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *kTOFPID=0)
{

   // Kaon TOF PID Cut
   AliRsnValueDaughter *axisKaonTOFPIDCut = new AliRsnValueDaughter("kaonTOFPID", AliRsnValueDaughter::kTOFnsigmaK);
   axisKaonTOFPIDCut->SetBins(0.0,5,0.1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonTOFPIDCut = new AliRsnListOutput("Kaon_TOF_PID_Cut", AliRsnListOutput::kHistoDefault);
   outMonitorKaonTOFPIDCut->AddValue(axisKaonTOFPIDCut);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonTOFPIDCut);
   if (kTOFPID) kTOFPID->AddOutput(outMonitorKaonTOFPIDCut);
  
}

void AddMonitorOutput_KaonNTPC(TObjArray *mon=0,TString opt="",AliRsnLoopDaughter *kNTPC=0)
{

   // Kaon PID Cut
   AliRsnValueDaughter *axisKaonNTPC = new AliRsnValueDaughter("kaonNTPC", AliRsnValueDaughter::kNTPCclusters);
   axisKaonNTPC->SetBins(0.0,200,1);

   // output: 2D histogram
   AliRsnListOutput *outMonitorKaonNTPC = new AliRsnListOutput("Kaon_NTPC", AliRsnListOutput::kHistoDefault);
   outMonitorKaonNTPC->AddValue(axisKaonNTPC);

   // add outputs to loop
   if (mon) mon->Add(outMonitorKaonNTPC);
   if (kNTPC) kNTPC->AddOutput(outMonitorKaonNTPC);
  
}

