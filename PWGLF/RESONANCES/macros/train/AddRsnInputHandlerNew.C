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
   
   // PID cuts for all needed detectors:
   // use new implementation based on AliPIDResponse and its related task
   // requires that AliAnalysisTaskPIDResponse is added otherwise it will raise errors
   AliRsnCutPIDNSigma *cutPIDITSpion   = new AliRsnCutPIDNSigma("cutPIDITSpion"  , AliPID::kPion, AliRsnCutPIDNSigma::kITS, 3.0);
   AliRsnCutPIDNSigma *cutPIDTPCpion   = new AliRsnCutPIDNSigma("cutPIDTPCpion"  , AliPID::kPion, AliRsnCutPIDNSigma::kTPC, 3.0);
   AliRsnCutPIDNSigma *cutPIDTOFpion   = new AliRsnCutPIDNSigma("cutPIDITSpion"  , AliPID::kPion, AliRsnCutPIDNSigma::kTOF, 3.0);
   AliRsnCutPIDNSigma *cutPIDITSkaon   = new AliRsnCutPIDNSigma("cutPIDITSkaon"  , AliPID::kKaon, AliRsnCutPIDNSigma::kITS, 3.0);
   AliRsnCutPIDNSigma *cutPIDTPCkaonLo = new AliRsnCutPIDNSigma("cutPIDTPCkaonLo", AliPID::kKaon, AliRsnCutPIDNSigma::kTPC, 5.0);
   AliRsnCutPIDNSigma *cutPIDTPCkaonHi = new AliRsnCutPIDNSigma("cutPIDTPCkaonHi", AliPID::kKaon, AliRsnCutPIDNSigma::kTPC, 3.0);
   AliRsnCutPIDNSigma *cutPIDTOFkaon   = new AliRsnCutPIDNSigma("cutPIDITSkaon"  , AliPID::kKaon, AliRsnCutPIDNSigma::kTOF, 3.0);
   
   // ITS PID:
   // - reject unmatched tracks
   cutPIDITSpion->SetRejectUnmatched();
   cutPIDITSkaon->SetRejectUnmatched();
   
   // TPC PID:
   // - pions               --> 3 sigma cut
   // - kaons below 350 MeV --> 5 sigma cut
   // - kaons above 350 MeV --> 3 sigma cut
   // - reject unmatched tracks
   cutPIDTPCkaonLo->SetMomentumRange(0.00, 0.35);
   cutPIDTPCkaonHi->SetMomentumRange(0.35, 1E20);
   cutPIDTPCkaonLo->SetRejectOutside();
   cutPIDTPCkaonHi->SetRejectOutside();
   cutPIDTPCpion  ->SetRejectUnmatched();
   cutPIDTPCkaonLo->SetRejectUnmatched();
   cutPIDTPCkaonHi->SetRejectUnmatched();
   
   // TOF PID:
   // must specify that unmatched tracks must be accepted
   cutPIDTOFpion->SetRejectUnmatched(kFALSE);
   cutPIDTOFkaon->SetRejectUnmatched(kFALSE);
   
   //---------------------------------------------
   //  Combine cuts
   //---------------------------------------------
   
   // make several combinations of cuts:
   // ITS and TPC standard
   //AliRsnCutSet *cutsQualityITS = new AliRsnCutSet("qualityITS", AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsQualityTPC = new AliRsnCutSet("qualityTPC", AliRsnTarget::kDaughter);
   //AliRsnCutSet *cutsPionITS    = new AliRsnCutSet("pionITS"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsPionTPC    = new AliRsnCutSet("pionTPC"   , AliRsnTarget::kDaughter);
   //AliRsnCutSet *cutsKaonITS    = new AliRsnCutSet("kaonITS"   , AliRsnTarget::kDaughter);
   AliRsnCutSet *cutsKaonTPC    = new AliRsnCutSet("kaonTPC"   , AliRsnTarget::kDaughter);
   
   // ITS standalone: quality only
   //cutsQualityITS->AddCut(cutQualityITS);
   //cutsQualityITS->SetCutScheme(cutQualityITS->GetName());
   
   // TPC+ITS: quality only
   cutsQualityTPC->AddCut(cutQualityTPC);
   cutsQualityTPC->SetCutScheme(cutQualityTPC->GetName());
   
   // ITS standalone: quality and ITS PID
   //cutsPionITS->AddCut(cutQualityITS);
   //cutsPionITS->AddCut(cutPIDITSpion);
   //cutsPionITS->SetCutScheme(Form("%s&%s", cutQualityITS->GetName(), cutPIDITSpion->GetName()));
   //cutsKaonITS->AddCut(cutQualityITS);
   //cutsKaonITS->AddCut(cutPIDITSkaon);
   //cutsKaonITS->SetCutScheme(Form("%s&%s", cutQualityITS->GetName(), cutPIDITSkaon->GetName()));
   
   // TPC standalone: quality and TPC PID
   cutsPionTPC->AddCut(cutQualityTPC);
   cutsPionTPC->AddCut(cutPIDTPCpion);
   cutsPionTPC->AddCut(cutPIDTOFpion);
   cutsPionTPC->SetCutScheme(Form("%s&%s&%s", cutQualityTPC->GetName(), cutPIDTPCpion->GetName(), cutPIDTOFpion->GetName()));
   
   cutsKaonTPC->AddCut(cutQualityTPC);
   cutsKaonTPC->AddCut(cutPIDTPCkaonLo);
   cutsKaonTPC->AddCut(cutPIDTPCkaonHi);
   cutsKaonTPC->AddCut(cutPIDTOFkaon);
   cutsKaonTPC->SetCutScheme(Form("%s&(%s|%s)&%s", cutQualityTPC->GetName(), cutPIDTPCkaonLo->GetName(), cutPIDTPCkaonHi->GetName(), cutPIDTOFkaon->GetName()));
   
   // setup selector in the handler and add RSN input handler
   AliRsnInputHandler *rsnIH = new AliRsnInputHandler();
   AliRsnDaughterSelector *sel = rsnIH->GetSelector();
   //sel->Add(cutsQualityITS, kTRUE);
   sel->Add(cutsQualityTPC, kTRUE);
   //sel->Add(cutsPionITS, kTRUE);
   sel->Add(cutsPionTPC, kTRUE);
   //sel->Add(cutsKaonITS, kTRUE);
   sel->Add(cutsKaonTPC, kTRUE);
   //sel->Add(cutsKaonAll, kTRUE);
   sel->Init();
   multi->AddInputEventHandler(rsnIH);
}
