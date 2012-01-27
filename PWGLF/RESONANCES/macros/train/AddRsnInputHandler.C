//
// This macro adds the specific RSN input handler, with all single particle cuts
//
void AddRsnInputHandler(Bool_t isMC, AliMultiInputEventHandler *multi)
{
   void myError  (const char *msg) {::Error  ("AddRsnInputHandler", msg);}
   void myWarning(const char *msg) {::Warning("AddRsnInputHandler", msg);}
   void myInfo   (const char *msg) {::Info   ("AddRsnInputHandler", msg);}
   
   if (!multi) {
      ::Error("AddRsnInputHandler", "Required a WELL INITIALIZED AliMultiInputEventHandler object");
      return;
   }
   
   //---------------------------------------------
   //  Define single cuts
   //---------------------------------------------

   // Track quality for ITS standalone:
   // this cut is used to select tracks of good quality, irrespective of the PID.
   // When adding status flags, the second argument tells if each considered flag
   // must be active or not in the track status, since the ITS-SA tracks need that
   // some of them are OFF (e.g.: kTPCin)
   AliRsnCutTrackQuality *cutQualityITS = new AliRsnCutTrackQuality("cutQualityITS");
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSin    , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCin    , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSrefit , kTRUE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kTPCrefit , kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpureSA, kFALSE);
   cutQualityITS->AddStatusFlag(AliESDtrack::kITSpid   , kTRUE);
   cutQualityITS->SetPtRange(0.15, 1E+20);
   cutQualityITS->SetEtaRange(-0.8, 0.8);
   cutQualityITS->SetDCARPtFormula("0.0595+0.0182/pt^1.55");
   cutQualityITS->SetDCAZmax(2.0);
   cutQualityITS->SetSPDminNClusters(1);
   cutQualityITS->SetITSminNClusters(4);
   cutQualityITS->SetITSmaxChi2(2.0);
   cutQualityITS->SetTPCminNClusters(0);
   cutQualityITS->SetTPCmaxChi2(1E+10);
   cutQualityITS->SetRejectKinkDaughters();
      
   // Track quality for TPC+ITS:
   // works exactly like the one above, but has settings for selecting TPC+ITS tracks
   // in this case, the flags required are all necessary, so here the procedure is simpler
   AliRsnCutTrackQuality *cutQualityTPC = new AliRsnCutTrackQuality("cutQualityTPC");
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCin   , kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kTPCrefit, kTRUE);
   cutQualityTPC->AddStatusFlag(AliESDtrack::kITSrefit, kTRUE);
   cutQualityTPC->SetPtRange(0.15, 1E+20);
   cutQualityTPC->SetEtaRange(-0.8, 0.8);
   cutQualityTPC->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   cutQualityTPC->SetDCAZmax(2.0);
   cutQualityTPC->SetSPDminNClusters(1);
   cutQualityTPC->SetITSminNClusters(0);
   cutQualityTPC->SetITSmaxChi2(1E+20);
   cutQualityTPC->SetTPCminNClusters(70);
   cutQualityTPC->SetTPCmaxChi2(4.0);
   cutQualityTPC->SetRejectKinkDaughters();
   
   // ITS PID
   // -- pions: 3sigma dE/dx cut
   // -- kaons: 3sigma dE/dx cut
   //
   // NOTE:
   // --> need to know if running on MC or data (for initializing the BB)
   AliRsnCutPIDITS *cutPIDITSpion = new AliRsnCutPIDITS("cutPIDITSpion", AliPID::kPion, -3.0, 3.0);
   AliRsnCutPIDITS *cutPIDITSkaon = new AliRsnCutPIDITS("cutPIDITSkaon", AliPID::kKaon, -3.0, 3.0);
   cutPIDITSpion->SetMC(isMC);
   cutPIDITSkaon->SetMC(isMC);
   
   // TPC PID
   // -- pions: 3sigma dE/dx cut
   // -- kaons: 3(5)sigma dE/dx cut for p_TPC above(below) 350 MeV/c
   //
   // NOTE:
   // --> The initialization of the BB is different between data and MC.
   AliRsnCutPIDTPC *cutPIDTPCpion  = new AliRsnCutPIDTPC("cutPIDTPCpion" , AliPID::kPion, -3.0, 3.0);
   AliRsnCutPIDTPC *cutPIDTPCkaonL = new AliRsnCutPIDTPC("cutPIDTPCkaonL", AliPID::kKaon, -5.0, 5.0);
   AliRsnCutPIDTPC *cutPIDTPCkaonH = new AliRsnCutPIDTPC("cutPIDTPCkaonH", AliPID::kKaon, -3.0, 3.0);
   
   // assign the momentum range and tell to reject tracks outside it
   cutPIDTPCkaonL->SetMomentumRange(0.00, 0.35);
   cutPIDTPCkaonH->SetMomentumRange(0.35, 1E20);
   cutPIDTPCkaonL->SetRejectOutside(kTRUE);
   cutPIDTPCkaonH->SetRejectOutside(kTRUE);
   
   // BB parameterization depends on data sample (MC, data)
   // the momentum range is passed and tracks outside it are rejected
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
   cutPIDTPCpion ->SetBBParam(bbPar);
   cutPIDTPCkaonL->SetBBParam(bbPar);
   cutPIDTPCkaonH->SetBBParam(bbPar);
   
   // TOF PID
   // -- pions: 3sigma cut
   // -- kaons: 3sigma cut
   //
   // NOTE:
   // --> since we use also TPC, unmatched tracks are accepted
   AliRsnCutPIDTOF *cutPIDTOFpion = new AliRsnCutPIDTOF("cutPIDTOFpion", AliPID::kPion, -3.0, 3.0);
   AliRsnCutPIDTOF *cutPIDTOFkaon = new AliRsnCutPIDTOF("cutPIDTOFkaon", AliPID::kKaon, -3.0, 3.0);
   cutPIDTOFpion->SetRejectUnmatched(kFALSE);
   cutPIDTOFkaon->SetRejectUnmatched(kFALSE);
   
   //---------------------------------------------
   //  Combine cuts
   //---------------------------------------------
   
   // make several combinations of cuts:
   // ITS and TPC standard
   AliRsnCutSet *cutsQualityITS = new AliRsnCutSet("qualityITS", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsQualityTPC = new AliRsnCutSet("qualityTPC", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsPionITS    = new AliRsnCutSet("pionITS"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsPionTPC    = new AliRsnCutSet("pionTPC"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsPionAll    = new AliRsnCutSet("pionAll"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsKaonITS    = new AliRsnCutSet("kaonITS"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsKaonTPC    = new AliRsnCutSet("kaonTPC"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsKaonAll    = new AliRsnCutSet("kaonAll"   , AliRsnTarget::kDaughter);
   
   // ITS standalone: quality only
   cutsQualityITS->AddCut(cutQualityITS);
   cutsQualityITS->SetCutScheme(cutQualityITS->GetName());
   
   // TPC+ITS: quality only
   cutsQualityTPC->AddCut(cutQualityTPC);
   cutsQualityTPC->SetCutScheme(cutQualityTPC->GetName());
   
   // ITS standalone: quality and ITS PID (pions)
   cutsPionITS->AddCut(cutQualityITS);
   cutsPionITS->AddCut(cutPIDITSpion);
   cutsPionITS->SetCutScheme(Form("%s&%s", cutQualityITS->GetName(), cutPIDITSpion->GetName()));
   
   // ITS standalone: quality and ITS PID (kaons)
   cutsKaonITS->AddCut(cutQualityITS);
   cutsKaonITS->AddCut(cutPIDITSkaon);
   cutsKaonITS->SetCutScheme(Form("%s&%s", cutQualityITS->GetName(), cutPIDITSkaon->GetName()));
   
   // TPC standalone: quality and TPC + TOF PID (pions)
   cutsPionTPC->AddCut(cutQualityTPC);
   cutsPionTPC->AddCut(cutPIDTPCpion);
   cutsPionTPC->AddCut(cutPIDTOFpion);
   cutsPionTPC->SetCutScheme(Form("%s & %s & %s", cutQualityTPC->GetName(), cutPIDTPCpion->GetName(), cutPIDTOFpion->GetName()));
   
   // TPC standalone: quality and TPC + TOF PID (kaons)
   cutsKaonTPC->AddCut(cutQualityTPC);
   cutsKaonTPC->AddCut(cutPIDTPCkaonL);
   cutsKaonTPC->AddCut(cutPIDTPCkaonH);
   cutsKaonTPC->AddCut(cutPIDTOFkaon);
   cutsKaonTPC->SetCutScheme(Form("%s & (%s|%s) & %s", cutQualityTPC->GetName(), cutPIDTPCkaonL->GetName(), cutPIDTPCkaonH->GetName(), cutPIDTOFkaon->GetName()));
   
   // all tracks: both qualities and PIDs (pion)
   cutsPionAll->AddCut(cutQualityITS);
   cutsPionAll->AddCut(cutQualityTPC);
   cutsPionAll->AddCut(cutPIDITSpion);
   cutsPionAll->AddCut(cutPIDTPCpion);
   cutsPionAll->AddCut(cutPIDTOFpion);
   cutsPionAll->SetCutScheme(Form("(%s)|(%s)", cutsKaonITS->GetCutScheme().Data(), cutsKaonTPC->GetCutScheme().Data()));
   
   // all tracks: both qualities and PIDs (kaon)
   cutsKaonAll->AddCut(cutQualityITS);
   cutsKaonAll->AddCut(cutQualityTPC);
   cutsKaonAll->AddCut(cutPIDITSkaon);
   cutsKaonAll->AddCut(cutPIDTPCkaonL);
   cutsKaonAll->AddCut(cutPIDTPCkaonH);
   cutsKaonAll->AddCut(cutPIDTOFkaon);
   cutsKaonAll->SetCutScheme(Form("(%s)|(%s)", cutsKaonITS->GetCutScheme().Data(), cutsKaonTPC->GetCutScheme().Data()));
   
   // setup selector in the handler and add RSN input handler
   AliRsnInputHandler *rsnIH = new AliRsnInputHandler();
   AliRsnDaughterSelector *sel = rsnIH->GetSelector();
// sel->Add(cutsQualityITS, kTRUE);
   sel->Add(cutsQualityTPC, kTRUE);
// sel->Add(cutsPionITS   , kTRUE);
// sel->Add(cutsPionTPC   , kTRUE);
// sel->Add(cutsPionAll   , kTRUE);
// sel->Add(cutsKaonITS   , kTRUE);
   sel->Add(cutsKaonTPC   , kTRUE);
// sel->Add(cutsKaonAll   , kTRUE);
   sel->Init();
   multi->AddInputEventHandler(rsnIH);
}
