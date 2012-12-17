// #include "AliMuonForwardTrackAnalysis.h"
// #include "TDatabasePDG.h"
// #include "TGeoGlobalMagField.h"
// #include "TROOT.h"
// #include "AliMagF.h"

enum {kNoOption, kResonanceOnly, kCharmOnly, kBeautyOnly, kBackground1mu, kBackground2mu, kNoResonances};

//=============================================================================================================================================================

void AliMuonForwardTrackAnalysis(const Char_t *readDir= ".",                       // the directory with the MuonGlobalTracks.root and geometry.root files
				 Int_t option = kNoOption,                         // for resonance analysis: kResonanceOnly 
				 Double_t massMin = 0.,                            // lower limit for the cut on dimuon mass
				 Double_t massMax = 10.,                           // upper limit for the cut on dimuon mass
				 Double_t maxChi2SingleMuons = 1.5,                // upper limit for the cut on the single muon chi2
				 Double_t maxOffsetSingleMuons = 1.e9.,            // upper limit for the cut on the single muon offset w.r.t. the primary vtx
				 Bool_t correlateCutOnOffsetChi2 = kTRUE,          // if true, the cut region in the chi2-offset plane for single muons is a quadrant aorund the origin
				 Double_t maxWOffsetMuonPairsAtPrimaryVtx = 1.e9,  // upper limit for the cut on weighted offset of dimuond w.r.t. the primary vtx
				 Double_t maxWOffsetMuonPairsAtPCA = 1.e9,         // upper limit for the cut on weighted offset of dimuond w.r.t. their PCA
				 Double_t maxDistancePrimaryVtxPCA = 1.e9,         // upper limit for the cut on the distance between primary vtx and PCA
				 Double_t minPCAQuality = 0.,                      // lower limit for the cut on the PCA quality
				 Int_t triggerLevel = 1,                           // level of the trigger both muons must satisfy
				 Int_t maxNWrongClusters = 999,                    // maximum number of wrong MFT clusters for a global muon track
				 const Char_t *outDir = ".",                       // directory where the output file will be created
				 Bool_t singleMuonAnalysis = kTRUE,                // if true, the aalysis of single muons will be performed
				 Bool_t muonPairAnalysis = kTRUE,                  // if true, the aalysis of muon pairs will be performed
				 Int_t firstEvent = -1,
				 Int_t lastEvent = -1, 
				 Int_t numTag = 0,                                 // number which will tag the name of the output file
				 Double_t ptMinSingleMuons = 0.0,                  // lower limit for the cut on the single muon pt
				 Double_t trueMass = 3.097,                        // used to evaluate the pseudo proper decay length, usually for J/psi only
				 Bool_t evalDimuonVtxResolution=kFALSE,            // to be set true only if prompt dimuon sources are analyzed
				 Int_t nEventsToMix = 0,                           // if <1 or >100, mixing is not performed
				 const Char_t *tag = "noTag",                      // tag added to the output file name
				 Double_t etaMinSingleMuons = -3.6,                // lower limit for the cut on the single muon eta
				 Double_t etaMaxSingleMuons = -2.5) {              // upper limit for the cut on the single muon eta
  
  const Double_t mJpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();

  if (trueMass<0) trueMass = mJpsi;
 
  gROOT -> LoadMacro("./AliMuonForwardTrackAnalysis.cxx+");
  //  AliLog::SetClassDebugLevel("AliMuonForwardTrackPair", 1);
  //  AliLog::SetClassDebugLevel("AliMuonForwardTrackAnalysis", 2);

  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  
  AliMuonForwardTrackAnalysis *myAnalysis = new AliMuonForwardTrackAnalysis();
  myAnalysis->ReadEvents(firstEvent, lastEvent);
  myAnalysis->SetInputDir(readDir);
  myAnalysis->SetOutputDir(outDir);
  myAnalysis->SetMassRange(massMin, massMax);
  myAnalysis->SetTrueMass(trueMass);
  myAnalysis->SetSingleMuonAnalysis(singleMuonAnalysis);
  myAnalysis->SetMuonPairAnalysis(muonPairAnalysis);
  myAnalysis->SetOption(option);

  myAnalysis->SetMaxNWrongClustersMC(maxNWrongClusters);
  myAnalysis->SetMinPtSingleMuons(ptMinSingleMuons);
  //  myAnalysis->SetEtaRangeSingleMuons(etaMinSingleMuons, etaMaxSingleMuons);
  myAnalysis->SetMaxChi2SingleMuons(maxChi2SingleMuons);
  myAnalysis->SetMaxOffsetSingleMuons(maxOffsetSingleMuons);
  myAnalysis->CorrelateCutOnOffsetChi2(correlateCutOnOffsetChi2);

  myAnalysis->SetMaxWOffsetMuonPairsAtPrimaryVtx(maxWOffsetMuonPairsAtPrimaryVtx);
  myAnalysis->SetMaxWOffsetMuonPairsAtPCA(maxWOffsetMuonPairsAtPCA);
  myAnalysis->SetMaxDistancePrimaryVtxPCA(maxDistancePrimaryVtxPCA);
  myAnalysis->SetMinPCAQuality(minPCAQuality);

  myAnalysis->SetMatchTrigger(triggerLevel);

  myAnalysis->EvalDimuonVtxResolution(evalDimuonVtxResolution);    // it should be true only with prompt dimuon sources

  myAnalysis->SetNEventsToMix(nEventsToMix);

  if (myAnalysis->Init("MuonGlobalTracks.root")) {
    while (myAnalysis->LoadNextEvent()) continue;
    myAnalysis->Terminate(Form("outFiles/outFile.%d.%s.root", numTag, tag));
  }

}

//================================================================================================================================

