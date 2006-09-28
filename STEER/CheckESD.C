#if !defined( __CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TError.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TPDGCode.h>
#include <TParticle.h>

#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliESD.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliPID.h"
#else
const Int_t kXiMinus = 3312;
const Int_t kOmegaMinus = 3334;
#endif


TH1F* CreateHisto(const char* name, const char* title, 
		  Int_t nBins, Double_t xMin, Double_t xMax,
		  const char* xLabel = NULL, const char* yLabel = NULL)
{
// create a histogram

  TH1F* result = new TH1F(name, title, nBins, xMin, xMax);
  result->SetOption("E");
  if (xLabel) result->GetXaxis()->SetTitle(xLabel);
  if (yLabel) result->GetYaxis()->SetTitle(yLabel);
  result->SetMarkerStyle(kFullCircle);
  return result;
}

TH1F* CreateEffHisto(TH1F* hGen, TH1F* hRec)
{
// create an efficiency histogram

  Int_t nBins = hGen->GetNbinsX();
  TH1F* hEff = (TH1F*) hGen->Clone("hEff");
  hEff->SetTitle("");
  hEff->SetStats(kFALSE);
  hEff->SetMinimum(0.);
  hEff->SetMaximum(110.);
  hEff->GetYaxis()->SetTitle("#epsilon [%]");

  for (Int_t iBin = 0; iBin <= nBins; iBin++) {
    Double_t nGen = hGen->GetBinContent(iBin);
    Double_t nRec = hRec->GetBinContent(iBin);
    if (nGen > 0) {
      Double_t eff = nRec/nGen;
      hEff->SetBinContent(iBin, 100. * eff);
      Double_t error = sqrt(eff*(1.-eff) / nGen);
      if (error == 0) error = 0.0001;
      hEff->SetBinError(iBin, 100. * error);
    } else {
      hEff->SetBinContent(iBin, -100.);
      hEff->SetBinError(iBin, 0);
    }
  }

  return hEff;
}

Bool_t FitHisto(TH1* histo, Double_t& res, Double_t& resError)
{
// fit a gaussian to a histogram

  static TF1* fitFunc = new TF1("fitFunc", "gaus");
  fitFunc->SetLineWidth(2);
  fitFunc->SetFillStyle(0);
  Double_t maxFitRange = 2;

  if (histo->Integral() > 50) {
    Float_t mean = histo->GetMean();
    Float_t rms = histo->GetRMS();
    fitFunc->SetRange(mean - maxFitRange*rms, mean + maxFitRange*rms);
    fitFunc->SetParameters(mean, rms);
    histo->Fit(fitFunc, "QRI0");
    histo->GetFunction("fitFunc")->ResetBit(1<<9);
    res = TMath::Abs(fitFunc->GetParameter(2));
    resError = TMath::Abs(fitFunc->GetParError(2));
    return kTRUE;
  }

  return kFALSE;
}


Bool_t CheckESD(const char* gAliceFileName = "galice.root", 
		const char* esdFileName = "AliESDs.root")
{
// check the content of the ESD
 
  // check values
  Int_t    checkNGenLow = 1;

  Double_t checkEffLow = 0.5;
  Double_t checkEffSigma = 3;
  Double_t checkFakeHigh = 0.5;
  Double_t checkFakeSigma = 3;

  Double_t checkResPtInvHigh = 5;
  Double_t checkResPtInvSigma = 3;
  Double_t checkResPhiHigh = 10;
  Double_t checkResPhiSigma = 3;
  Double_t checkResThetaHigh = 10;
  Double_t checkResThetaSigma = 3;

  Double_t checkPIDEffLow = 0.5;
  Double_t checkPIDEffSigma = 3;
  Double_t checkResTOFHigh = 500;
  Double_t checkResTOFSigma = 3;

  Double_t checkPHOSNLow = 5;
  Double_t checkPHOSEnergyLow = 0.3;
  Double_t checkPHOSEnergyHigh = 1.0;
  Double_t checkEMCALNLow = 50;
  Double_t checkEMCALEnergyLow = 0.05;
  Double_t checkEMCALEnergyHigh = 1.0;

  Double_t checkMUONNLow = 1;
  Double_t checkMUONPtLow = 0.5;
  Double_t checkMUONPtHigh = 10.;

  Double_t cutPtV0 = 0.3;
  Double_t checkV0EffLow = 0.02;
  Double_t checkV0EffSigma = 3;
  Double_t cutPtCascade = 0.5;
  Double_t checkCascadeEffLow = 0.01;
  Double_t checkCascadeEffSigma = 3;

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open(gAliceFileName);
  if (!runLoader) {
    Error("CheckESD", "getting run loader from file %s failed", 
	    gAliceFileName);
    return kFALSE;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("CheckESD", "no galice object found");
    return kFALSE;
  }
  runLoader->LoadKinematics();
  runLoader->LoadHeader();

  // open the ESD file
  TFile* esdFile = TFile::Open(esdFileName);
  if (!esdFile || !esdFile->IsOpen()) {
    Error("CheckESD", "opening ESD file %s failed", esdFileName);
    return kFALSE;
  }
  AliESD* esd = new AliESD;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("CheckESD", "no ESD tree found");
    return kFALSE;
  }
  tree->SetBranchAddress("ESD", &esd);

  // efficienc and resolution histograms
  Int_t nBinsPt = 15;
  Float_t minPt = 0.1;
  Float_t maxPt = 3.1;
  TH1F* hGen = CreateHisto("hGen", "generated tracks", 
			   nBinsPt, minPt, maxPt, "p_{t} [GeV/c]", "N");
  TH1F* hRec = CreateHisto("hRec", "reconstructed tracks", 
			   nBinsPt, minPt, maxPt, "p_{t} [GeV/c]", "N");
  Int_t nGen = 0;
  Int_t nRec = 0;
  Int_t nFake = 0;

  TH1F* hResPtInv = CreateHisto("hResPtInv", "", 100, -10, 10, 
           "(p_{t,rec}^{-1}-p_{t,sim}^{-1}) / p_{t,sim}^{-1} [%]", "N");
  TH1F* hResPhi = CreateHisto("hResPhi", "", 100, -20, 20, 
			      "#phi_{rec}-#phi_{sim} [mrad]", "N");
  TH1F* hResTheta = CreateHisto("hResTheta", "", 100, -20, 20, 
				"#theta_{rec}-#theta_{sim} [mrad]", "N");

  // PID
  Int_t partCode[AliPID::kSPECIES] = 
    {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  const char* partName[AliPID::kSPECIES+1] = 
    {"electron", "muon", "pion", "kaon", "proton", "other"};
  Double_t partFrac[AliPID::kSPECIES] = 
    {0.01, 0.01, 0.85, 0.10, 0.05};
  Int_t identified[AliPID::kSPECIES+1][AliPID::kSPECIES];
  for (Int_t iGen = 0; iGen < AliPID::kSPECIES+1; iGen++) {
    for (Int_t iRec = 0; iRec < AliPID::kSPECIES; iRec++) {
      identified[iGen][iRec] = 0;
    }
  }
  Int_t nIdentified = 0;

  // dE/dx and TOF
  TH2F* hDEdxRight = new TH2F("hDEdxRight", "", 300, 0, 3, 100, 0, 400);
  hDEdxRight->SetStats(kFALSE);
  hDEdxRight->GetXaxis()->SetTitle("p [GeV/c]");
  hDEdxRight->GetYaxis()->SetTitle("dE/dx_{TPC}");
  hDEdxRight->SetMarkerStyle(kFullCircle);
  hDEdxRight->SetMarkerSize(0.4);
  TH2F* hDEdxWrong = new TH2F("hDEdxWrong", "", 300, 0, 3, 100, 0, 400);
  hDEdxWrong->SetStats(kFALSE);
  hDEdxWrong->GetXaxis()->SetTitle("p [GeV/c]");
  hDEdxWrong->GetYaxis()->SetTitle("dE/dx_{TPC}");
  hDEdxWrong->SetMarkerStyle(kFullCircle);
  hDEdxWrong->SetMarkerSize(0.4);
  hDEdxWrong->SetMarkerColor(kRed);
  TH1F* hResTOFRight = CreateHisto("hResTOFRight", "", 100, -1000, 1000, 
				   "t_{TOF}-t_{track} [ps]", "N");
  TH1F* hResTOFWrong = CreateHisto("hResTOFWrong", "", 100, -1000, 1000, 
				   "t_{TOF}-t_{track} [ps]", "N");
  hResTOFWrong->SetLineColor(kRed);

  // calorimeters
  TH1F* hEPHOS = CreateHisto("hEPHOS", "PHOS", 100, 0, 5, "E [GeV]", "N");
  TH1F* hEEMCAL = CreateHisto("hEEMCAL", "EMCAL", 100, 0, 50, "E [GeV]", "N");

  // muons
  TH1F* hPtMUON = CreateHisto("hPtMUON", "MUON", 100, 0, 20, 
			      "p_{t} [GeV/c]", "N");

  // V0s and cascades
  TH1F* hMassK0 = CreateHisto("hMassK0", "K^{0}", 100, 0.4, 0.6, 
			      "M(#pi^{+}#pi^{-}) [GeV/c^{2}]", "N");
  TH1F* hMassLambda = CreateHisto("hMassLambda", "#Lambda", 100, 1.0, 1.2, 
				  "M(p#pi^{-}) [GeV/c^{2}]", "N");
  TH1F* hMassLambdaBar = CreateHisto("hMassLambdaBar", "#bar{#Lambda}", 
				     100, 1.0, 1.2, 
				     "M(#bar{p}#pi^{+}) [GeV/c^{2}]", "N");
  Int_t nGenV0s = 0;
  Int_t nRecV0s = 0;
  TH1F* hMassXi = CreateHisto("hMassXi", "#Xi", 100, 1.2, 1.5, 
			      "M(#Lambda#pi) [GeV/c^{2}]", "N");
  TH1F* hMassOmega = CreateHisto("hMassOmega", "#Omega", 100, 1.5, 1.8, 
				 "M(#LambdaK) [GeV/c^{2}]", "N");
  Int_t nGenCascades = 0;
  Int_t nRecCascades = 0;

  // loop over events
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);

    // select simulated primary particles, V0s and cascades
    AliStack* stack = gAlice->Stack();
    Int_t nParticles = stack->GetNtrack();
    TArrayF vertex(3);
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
    TObjArray selParticles;
    TObjArray selV0s;
    TObjArray selCascades;
    for (Int_t iParticle = 0; iParticle < nParticles; iParticle++) {
      TParticle* particle = stack->Particle(iParticle);
      if (!particle) continue;
      if (particle->Pt() < 0.001) continue;
      if (TMath::Abs(particle->Eta()) > 0.9) continue;
      TVector3 dVertex(particle->Vx() - vertex[0], 
		       particle->Vy() - vertex[1],
		       particle->Vz() - vertex[2]);
      if (dVertex.Mag() > 0.0001) continue;

      switch (TMath::Abs(particle->GetPdgCode())) {
      case kElectron:
      case kMuonMinus:
      case kPiPlus:
      case kKPlus:
      case kProton: {
	if (particle->Pt() > minPt) {
	  selParticles.Add(particle);
	  nGen++;
	  hGen->Fill(particle->Pt());
	}
	break;
      }
      case kK0Short:
      case kLambda0: {
	if (particle->Pt() > cutPtV0) {
	  nGenV0s++;
	  selV0s.Add(particle);
	}
	break;
      }
      case kXiMinus:
      case kOmegaMinus: {
	if (particle->Pt() > cutPtCascade) {
	  nGenCascades++;
	  selCascades.Add(particle);
	}
	break;
      }
      default: break;
      }
    }

    // get the event summary data
    tree->GetEvent(iEvent);
    if (!esd) {
      Error("CheckESD", "no ESD object found for event %d", iEvent);
      return kFALSE;
    }

    // loop over tracks
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
      AliESDtrack* track = esd->GetTrack(iTrack);

      // select tracks of selected particles
      Int_t label = TMath::Abs(track->GetLabel());
      if (label > stack->GetNtrack()) continue;     // background
      TParticle* particle = stack->Particle(label);
      if (!selParticles.Contains(particle)) continue;
      if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) continue;
      if (track->GetConstrainedChi2() > 1e9) continue;
      selParticles.Remove(particle);   // don't count multiple tracks

      nRec++;
      hRec->Fill(particle->Pt());
      if (track->GetLabel() < 0) nFake++;

      // resolutions
      Double_t p[3];
      track->GetConstrainedPxPyPz(p);
      TVector3 pTrack(p);
      hResPtInv->Fill(100. * (1./pTrack.Pt() - 1./particle->Pt()) * 
		      particle->Pt());
      hResPhi->Fill(1000. * (pTrack.Phi() - particle->Phi()));
      hResTheta->Fill(1000. * (pTrack.Theta() - particle->Theta()));

      // PID
      if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) continue;
      Int_t iGen = 5;
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
	if (TMath::Abs(particle->GetPdgCode()) == partCode[i]) iGen = i;
      }
      Double_t probability[AliPID::kSPECIES];
      track->GetESDpid(probability);
      Double_t pMax = 0;
      Int_t iRec = 0;
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
	probability[i] *= partFrac[i];
	if (probability[i] > pMax) {
	  pMax = probability[i];
	  iRec = i;
	}
      }
      identified[iGen][iRec]++;
      if (iGen == iRec) nIdentified++;

      // dE/dx and TOF
      Double_t time[AliPID::kSPECIES];
      track->GetIntegratedTimes(time);
      if (iGen == iRec) {
	hDEdxRight->Fill(pTrack.Mag(), track->GetTPCsignal());
        if ((track->GetStatus() & AliESDtrack::kTOFpid) != 0) {
	  hResTOFRight->Fill(track->GetTOFsignal() - time[iRec]);
	}
      } else {
	hDEdxWrong->Fill(pTrack.Mag(), track->GetTPCsignal());
        if ((track->GetStatus() & AliESDtrack::kTOFpid) != 0) {
	  hResTOFWrong->Fill(track->GetTOFsignal() - time[iRec]);
	}
      }
    }

    // loop over muon tracks
    {
    for (Int_t iTrack = 0; iTrack < esd->GetNumberOfMuonTracks(); iTrack++) {
      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);
      Double_t ptInv = TMath::Abs(muonTrack->GetInverseBendingMomentum());
      if (ptInv > 0.001) {
	hPtMUON->Fill(1./ptInv);
      }
    }
    }

    // loop over V0s
    for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
      AliESDv0* v0 = esd->GetV0(iV0);
      if (v0->GetOnFlyStatus()) continue;
      v0->ChangeMassHypothesis(kK0Short);
      hMassK0->Fill(v0->GetEffMass());
      v0->ChangeMassHypothesis(kLambda0);
      hMassLambda->Fill(v0->GetEffMass());
      v0->ChangeMassHypothesis(kLambda0Bar);
      hMassLambdaBar->Fill(v0->GetEffMass());

      Int_t negLabel = TMath::Abs(esd->GetTrack(v0->GetNindex())->GetLabel());
      if (negLabel > stack->GetNtrack()) continue;     // background
      Int_t negMother = stack->Particle(negLabel)->GetMother(0);
      if (negMother < 0) continue;
      Int_t posLabel = TMath::Abs(esd->GetTrack(v0->GetPindex())->GetLabel());
      if (posLabel > stack->GetNtrack()) continue;     // background
      Int_t posMother = stack->Particle(posLabel)->GetMother(0);
      if (negMother != posMother) continue;
      TParticle* particle = stack->Particle(negMother);
      if (!selV0s.Contains(particle)) continue;
      selV0s.Remove(particle);
      nRecV0s++;
    }

    // loop over Cascades
    for (Int_t iCascade = 0; iCascade < esd->GetNumberOfCascades(); 
	 iCascade++) {
      AliESDcascade* cascade = esd->GetCascade(iCascade);
      Double_t v0q;
      cascade->ChangeMassHypothesis(v0q,kXiMinus);
      hMassXi->Fill(cascade->GetEffMass());
      cascade->ChangeMassHypothesis(v0q,kOmegaMinus);
      hMassOmega->Fill(cascade->GetEffMass());

      Int_t negLabel = TMath::Abs(esd->GetTrack(cascade->GetNindex())
				  ->GetLabel());
      if (negLabel > stack->GetNtrack()) continue;     // background
      Int_t negMother = stack->Particle(negLabel)->GetMother(0);
      if (negMother < 0) continue;
      Int_t posLabel = TMath::Abs(esd->GetTrack(cascade->GetPindex())
				  ->GetLabel());
      if (posLabel > stack->GetNtrack()) continue;     // background
      Int_t posMother = stack->Particle(posLabel)->GetMother(0);
      if (negMother != posMother) continue;
      Int_t v0Mother = stack->Particle(negMother)->GetMother(0);
      if (v0Mother < 0) continue;
      Int_t bacLabel = TMath::Abs(esd->GetTrack(cascade->GetBindex())
				  ->GetLabel());
      if (bacLabel > stack->GetNtrack()) continue;     // background
      Int_t bacMother = stack->Particle(bacLabel)->GetMother(0);
      if (v0Mother != bacMother) continue;
      TParticle* particle = stack->Particle(v0Mother);
      if (!selCascades.Contains(particle)) continue;
      selCascades.Remove(particle);
      nRecCascades++;
    }

    // loop over the PHOS clusters
    {
    Int_t firstPHOSCluster = esd->GetFirstPHOSCluster();
    Int_t lastPHOSCluster  = firstPHOSCluster + esd->GetNumberOfPHOSClusters();
    for (Int_t iCluster=firstPHOSCluster; iCluster<lastPHOSCluster; iCluster++)
      hEPHOS->Fill(esd->GetCaloCluster(iCluster)->GetClusterEnergy());
    }

    // loop over the EMCAL clusters
    {
    Int_t firstEMCALCluster = esd->GetFirstEMCALCluster();
    Int_t lastEMCALCluster  = firstEMCALCluster + esd->GetNumberOfEMCALClusters();
    for (Int_t iCluster=firstEMCALCluster; iCluster<lastEMCALCluster; iCluster++)
      hEEMCAL->Fill(esd->GetCaloCluster(iCluster)->GetClusterEnergy());
    }
  }

  // perform checks
  if (nGen < checkNGenLow) {
    Warning("CheckESD", "low number of generated particles: %d", Int_t(nGen));
  }

  TH1F* hEff = CreateEffHisto(hGen, hRec);

  Info("CheckESD", "%d out of %d tracks reconstructed including %d "
	 "fake tracks", nRec, nGen, nFake);
  if (nGen > 0) {
    // efficiency
    Double_t eff = nRec*1./nGen;
    Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGen);
    Double_t fake = nFake*1./nGen;
    Double_t fakeError = TMath::Sqrt(fake*(1.-fake) / nGen);
    Info("CheckESD", "eff = (%.1f +- %.1f) %%  fake = (%.1f +- %.1f) %%",
	 100.*eff, 100.*effError, 100.*fake, 100.*fakeError);

    if (eff < checkEffLow - checkEffSigma*effError) {
      Warning("CheckESD", "low efficiency: (%.1f +- %.1f) %%", 
	      100.*eff, 100.*effError);
    }
    if (fake > checkFakeHigh + checkFakeSigma*fakeError) {
      Warning("CheckESD", "high fake: (%.1f +- %.1f) %%", 
	      100.*fake, 100.*fakeError);
    }

    // resolutions
    Double_t res, resError;
    if (FitHisto(hResPtInv, res, resError)) {
      Info("CheckESD", "relative inverse pt resolution = (%.1f +- %.1f) %%",
	   res, resError);
      if (res > checkResPtInvHigh + checkResPtInvSigma*resError) {
	Warning("CheckESD", "bad pt resolution: (%.1f +- %.1f) %%", 
		res, resError);
      }
    }

    if (FitHisto(hResPhi, res, resError)) {
      Info("CheckESD", "phi resolution = (%.1f +- %.1f) mrad", res, resError);
      if (res > checkResPhiHigh + checkResPhiSigma*resError) {
	Warning("CheckESD", "bad phi resolution: (%.1f +- %.1f) mrad", 
		res, resError);
      }
    }

    if (FitHisto(hResTheta, res, resError)) {
      Info("CheckESD", "theta resolution = (%.1f +- %.1f) mrad", 
	   res, resError);
      if (res > checkResThetaHigh + checkResThetaSigma*resError) {
	Warning("CheckESD", "bad theta resolution: (%.1f +- %.1f) mrad", 
		res, resError);
      }
    }

    // PID
    if (nRec > 0) {
      Double_t eff = nIdentified*1./nRec;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nRec);
      Info("CheckESD", "PID eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkPIDEffLow - checkPIDEffSigma*effError) {
	Warning("CheckESD", "low PID efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }

    printf("%9s:", "gen\\rec");
    for (Int_t iRec = 0; iRec < AliPID::kSPECIES; iRec++) {
      printf("%9s", partName[iRec]);
    }
    printf("\n");
    for (Int_t iGen = 0; iGen < AliPID::kSPECIES+1; iGen++) {
      printf("%9s:", partName[iGen]);
      for (Int_t iRec = 0; iRec < AliPID::kSPECIES; iRec++) {
	printf("%9d", identified[iGen][iRec]);
      }
      printf("\n");
    }

    if (FitHisto(hResTOFRight, res, resError)) {
      Info("CheckESD", "TOF resolution = (%.1f +- %.1f) ps", res, resError);
      if (res > checkResTOFHigh + checkResTOFSigma*resError) {
	Warning("CheckESD", "bad TOF resolution: (%.1f +- %.1f) ps", 
		res, resError);
      }
    }

    // calorimeters
    if (hEPHOS->Integral() < checkPHOSNLow) {
      Warning("CheckESD", "low number of PHOS particles: %d", 
	      Int_t(hEPHOS->Integral()));
    } else {
      Double_t mean = hEPHOS->GetMean();
      if (mean < checkPHOSEnergyLow) {
	Warning("CheckESD", "low mean PHOS energy: %.1f GeV", mean);
      } else if (mean > checkPHOSEnergyHigh) {
	Warning("CheckESD", "high mean PHOS energy: %.1f GeV", mean);
      }
    }

    if (hEEMCAL->Integral() < checkEMCALNLow) {
      Warning("CheckESD", "low number of EMCAL particles: %d", 
	      Int_t(hEEMCAL->Integral()));
    } else {
      Double_t mean = hEEMCAL->GetMean();
      if (mean < checkEMCALEnergyLow) {
	Warning("CheckESD", "low mean EMCAL energy: %.1f GeV", mean);
      } else if (mean > checkEMCALEnergyHigh) {
	Warning("CheckESD", "high mean EMCAL energy: %.1f GeV", mean);
      }
    }

    // muons
    if (hPtMUON->Integral() < checkMUONNLow) {
      Warning("CheckESD", "low number of MUON particles: %d", 
	      Int_t(hPtMUON->Integral()));
    } else {
      Double_t mean = hPtMUON->GetMean();
      if (mean < checkMUONPtLow) {
	Warning("CheckESD", "low mean MUON pt: %.1f GeV/c", mean);
      } else if (mean > checkMUONPtHigh) {
	Warning("CheckESD", "high mean MUON pt: %.1f GeV/c", mean);
      }
    }

    // V0s
    if (nGenV0s > 0) {
      Double_t eff = nRecV0s*1./nGenV0s;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGenV0s);
      if (effError == 0) effError = checkV0EffLow / TMath::Sqrt(1.*nGenV0s);
      Info("CheckESD", "V0 eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkV0EffLow - checkV0EffSigma*effError) {
	Warning("CheckESD", "low V0 efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }

    // Cascades
    if (nGenCascades > 0) {
      Double_t eff = nRecCascades*1./nGenCascades;
      Double_t effError = TMath::Sqrt(eff*(1.-eff) / nGenCascades);
      if (effError == 0) effError = checkV0EffLow / 
			   TMath::Sqrt(1.*nGenCascades);
      Info("CheckESD", "Cascade eff = (%.1f +- %.1f) %%", 
	   100.*eff, 100.*effError);
      if (eff < checkCascadeEffLow - checkCascadeEffSigma*effError) {
	Warning("CheckESD", "low Cascade efficiency: (%.1f +- %.1f) %%", 
		100.*eff, 100.*effError);
      }
    }
  }

  // draw the histograms if not in batch mode
  if (!gROOT->IsBatch()) {
    new TCanvas;
    hEff->DrawCopy();
    new TCanvas;
    hResPtInv->DrawCopy("E");
    new TCanvas;
    hResPhi->DrawCopy("E");
    new TCanvas;
    hResTheta->DrawCopy("E");
    new TCanvas;
    hDEdxRight->DrawCopy();
    hDEdxWrong->DrawCopy("SAME");
    new TCanvas;
    hResTOFRight->DrawCopy("E");
    hResTOFWrong->DrawCopy("SAME");
    new TCanvas;
    hEPHOS->DrawCopy("E");
    new TCanvas;
    hEEMCAL->DrawCopy("E");
    new TCanvas;
    hPtMUON->DrawCopy("E");
    new TCanvas;
    hMassK0->DrawCopy("E");
    new TCanvas;
    hMassLambda->DrawCopy("E");
    new TCanvas;
    hMassLambdaBar->DrawCopy("E");
    new TCanvas;
    hMassXi->DrawCopy("E");
    new TCanvas;
    hMassOmega->DrawCopy("E");
  }

  // write the output histograms to a file
  TFile* outputFile = TFile::Open("check.root", "recreate");
  if (!outputFile || !outputFile->IsOpen()) {
    Error("CheckESD", "opening output file check.root failed");
    return kFALSE;
  }
  hEff->Write();
  hResPtInv->Write();
  hResPhi->Write();
  hResTheta->Write();
  hDEdxRight->Write();
  hDEdxWrong->Write();
  hResTOFRight->Write();
  hResTOFWrong->Write();
  hEPHOS->Write();
  hEEMCAL->Write();
  hPtMUON->Write();
  hMassK0->Write();
  hMassLambda->Write();
  hMassLambdaBar->Write();
  hMassXi->Write();
  hMassOmega->Write();
  outputFile->Close();
  delete outputFile;

  // clean up
  delete hGen;
  delete hRec;
  delete hEff;
  delete hResPtInv;
  delete hResPhi;
  delete hResTheta;
  delete hDEdxRight;
  delete hDEdxWrong;
  delete hResTOFRight;
  delete hResTOFWrong;
  delete hEPHOS;
  delete hEEMCAL;
  delete hPtMUON;
  delete hMassK0;
  delete hMassLambda;
  delete hMassLambdaBar;
  delete hMassXi;
  delete hMassOmega;

  delete esd;
  esdFile->Close();
  delete esdFile;

  runLoader->UnloadHeader();
  runLoader->UnloadKinematics();
  delete runLoader;

  // result of check
  Info("CheckESD", "check of ESD was successfull");
  return kTRUE;
}
