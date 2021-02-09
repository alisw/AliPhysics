#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH2I.h"

#include "THnSparse.h"

#include "AliMCParticle.h"
//#include "AliStack.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
//#include "AliTPCParamSR.h"

#include "AliTPCPIDEtaTree.h"

/*
This task determines the eta dependence of the TPC signal.
For this purpose, only tracks fulfilling some standard quality cuts are taken into account.
The obtained data can be used to derive the functional behaviour of the eta dependence.
Such a function can be plugged into this task to correct for the eta dependence and to see
if there is then really no eta dependence left.

Class written by Benjamin Hess.
Contact: bhess@cern.ch
*/

ClassImp(AliTPCPIDEtaTree)

//________________________________________________________________________
AliTPCPIDEtaTree::AliTPCPIDEtaTree()
  : AliTPCPIDBase()
  , fNumEtaCorrReqErrorsIssued(0)
  , fNumMultCorrReqErrorsIssued(0)
  , fStoreMultiplicity(kFALSE)
  , fStoreNumOfSubthresholdclusters(kFALSE)
  , fStoreNumClustersInActiveVolume(kFALSE)
  , fDoAdditionalQA(kFALSE)
  , fPtpcPionCut(1.0)
  , fPtpc(0)
  , fPt(0)
  , fDeDx(0)
  , fDeDxExpected(0)
  , fTanTheta(0)
  //, fSinAlpha(0)
  //, fY(0)
  , fPhiPrime(0)
  , fTPCsignalN(0)
  , fTPCsignalNsubthreshold(0)
  , fNumTPCClustersInActiveVolume(0)
  , fPIDtype(0)
  , fMultiplicity(0)
  , fCorrectdEdxEtaDependence(kFALSE)
  , fCorrectdEdxMultiplicityDependence(kFALSE)
  , fTree(0x0)
  , fTreePions(0x0)
  , fOutputContainer(0x0)
  , fhTOFqa(0x0)
  , fhMultiplicityQA(0x0)
{
  // default Constructor
}

//________________________________________________________________________
AliTPCPIDEtaTree::AliTPCPIDEtaTree(const char *name)
  : AliTPCPIDBase(name)
  , fNumEtaCorrReqErrorsIssued(0)
  , fNumMultCorrReqErrorsIssued(0)
  , fStoreMultiplicity(kFALSE)
  , fStoreNumOfSubthresholdclusters(kFALSE)
  , fStoreNumClustersInActiveVolume(kFALSE)
  , fDoAdditionalQA(kFALSE)
  , fPtpcPionCut(1.0)
  , fPtpc(0)
  , fPt(0)
  , fDeDx(0)
  , fDeDxExpected(0)
  , fTanTheta(0)
  //, fSinAlpha(0)
  //, fY(0)
  , fPhiPrime(0)
  , fTPCsignalN(0)
  , fTPCsignalNsubthreshold(0)
  , fNumTPCClustersInActiveVolume(0)
  , fPIDtype(0)
  , fMultiplicity(0)
  , fCorrectdEdxEtaDependence(kFALSE)
  , fCorrectdEdxMultiplicityDependence(kFALSE)
  , fTree(0x0)
  , fTreePions(0x0)
  , fOutputContainer(0x0)
  , fhTOFqa(0x0)
  , fhMultiplicityQA(0x0)
{
  // Constructor

  // Define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TObjArray::Class());
}


//________________________________________________________________________
AliTPCPIDEtaTree::~AliTPCPIDEtaTree()
{
  // dtor
  
  //delete fTree;
  //fTree = 0x0;
  
  delete fOutputContainer;
  fOutputContainer = 0x0;
}


//________________________________________________________________________
void AliTPCPIDEtaTree::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  AliTPCPIDBase::UserCreateOutputObjects();
  
  if (fDoAdditionalQA) {
    OpenFile(3);
    
    fOutputContainer = new TObjArray(2);
    fOutputContainer->SetName(GetName());
    fOutputContainer->SetOwner(kTRUE);
    
    const Int_t nBins = 6;
    // p_vert, p_TPC, eta, nSigmaTOF, beta, multiplicity
    Int_t bins[nBins] =    {   48,  48,    9,   30, 110,    4  };
    Double_t xmin[nBins] = {  0.3, 0.3, -0.9, -5.0, 0.0,   0.  };
    Double_t xmax[nBins] = {  2.7, 2.7,  0.9,  5.0, 1.1, 3200  };
    
    
    fhTOFqa = new THnSparseI("hTOFqa","", nBins, bins, xmin, xmax);
    fhTOFqa->GetAxis(0)->SetTitle("p_{Vertex} (GeV/c)");
    fhTOFqa->GetAxis(1)->SetTitle("p_{TPC_inner} (GeV/c)");
    fhTOFqa->GetAxis(2)->SetTitle("#eta");
    fhTOFqa->GetAxis(3)->SetTitle("n #sigma_{p} TOF");
    fhTOFqa->GetAxis(4)->SetTitle("#beta");
    fhTOFqa->GetAxis(5)->SetTitle("Multiplicity");
    
    fOutputContainer->Add(fhTOFqa);
    
    
    
    fhMultiplicityQA = new TH2I("hMultiplicityQA", "Multiplicity check; Contributors to primary vertex per event; Total number of tracks per Event",
                                120, 0, 6000, 400, 0, 20000);
    fOutputContainer->Add(fhMultiplicityQA);
  }
  else {
    fOutputContainer = new TObjArray(1);
    fOutputContainer->SetName(GetName());
    fOutputContainer->SetOwner(kTRUE);
  }
  
  OpenFile(1);
  
  fTree = new TTree("fTree", "Tree for analysis of #eta dependence of TPC signal");
  fTree->Branch("pTPC", &fPtpc);
  //fTree->Branch("pT", &fPt);
  fTree->Branch("dEdx", &fDeDx);
  fTree->Branch("dEdxExpected", &fDeDxExpected);
  fTree->Branch("tanTheta", &fTanTheta);
  //fTree->Branch("sinAlpha", &fSinAlpha);
  //fTree->Branch("y", &fY);
  //TODO fTree->Branch("phiPrime", &fPhiPrime);
  fTree->Branch("tpcSignalN", &fTPCsignalN);
  
  if (fStoreNumOfSubthresholdclusters)
    fTree->Branch("tpcSignalNsubthreshold", &fTPCsignalNsubthreshold);
  
  if (fStoreNumClustersInActiveVolume)
    fTree->Branch("numTPCClustersInActiveVolume", &fNumTPCClustersInActiveVolume);
  
  fTree->Branch("pidType", &fPIDtype);
  
  if (fStoreMultiplicity)  {
    fTree->Branch("fMultiplicity", &fMultiplicity);
  }
  
  OpenFile(2);
  
  fTreePions = new TTree("fTreePions", "Tree for analysis of #eta dependence of TPC signal - Pions");
  fTreePions->Branch("pTPC", &fPtpc);
  fTreePions->Branch("pT", &fPt);
  fTreePions->Branch("dEdx", &fDeDx);
  fTreePions->Branch("dEdxExpected", &fDeDxExpected);
  fTreePions->Branch("tanTheta", &fTanTheta);
  fTreePions->Branch("tpcSignalN", &fTPCsignalN);
  
  if (fStoreNumOfSubthresholdclusters)
    fTreePions->Branch("tpcSignalNsubthreshold", &fTPCsignalNsubthreshold);
  
  if (fStoreMultiplicity)  {
    fTreePions->Branch("fMultiplicity", &fMultiplicity);
  }
  
  PostData(1, fTree);
  PostData(2, fTreePions);
  PostData(3, fOutputContainer);
}


//________________________________________________________________________
void AliTPCPIDEtaTree::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  if (GetIsPileUp(fESD, kPileUpRejectionSPD))
  return;
  
  if (!fPIDResponse || !fV0KineCuts)
    return;
  
  if (!GetVertexIsOk(fESD))
    return;

  const AliVVertex* primaryVertex = fESD->GetPrimaryVertexTracks(); 
  if (!primaryVertex)
    return;
  //fMultiplicity = primaryVertex->GetNContributors();
  
  
  if(primaryVertex->GetNContributors() <= 0) 
    return;
    
  fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  
  fMultiplicity = fESD->GetNumberOfESDTracks(); 
    
  if (fDoAdditionalQA) {
    Int_t nTotTracks = fESD->GetNumberOfESDTracks();
    fhMultiplicityQA->Fill(primaryVertex->GetNContributors(), nTotTracks);
  }
  
  Double_t magField = fESD->GetMagneticField();
  
  // Fill V0 arrays with V0s
  FillV0PIDlist();
  
  //AliTPCParamSR par;
  //par.Update();
  
  
  Bool_t etaCorrAvail = fPIDResponse->UseTPCEtaCorrection();
  Bool_t multCorrAvail = fPIDResponse->UseTPCMultiplicityCorrection();
  
  static Bool_t pileupReported = kFALSE;
  if (!pileupReported) {
    if (fPIDResponse->GetTPCResponse().GetPileupCorrectionObject()) {
      Printf("INFO: Pileup correction found! The track dEdx signal will be corrected");
    }
    pileupReported = kTRUE;
  }

  // Track loop to fill a Train spectrum
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    Bool_t isPr = kFALSE;
    Bool_t isPi = kFALSE;
  
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    if (TMath::Abs(track->Eta()) > fEtaCut)
      continue;
    
    // Do not distinguish positively and negatively charged V0's
    Char_t v0tagAllCharges = TMath::Abs(GetV0tag(iTracks));
    if (v0tagAllCharges == -99) {
      AliError(Form("Problem with V0 tag list (requested V0 track for track %d from %d (list status %d))!", iTracks, fESD->GetNumberOfTracks(),
                    fV0tags != 0x0));
      continue;
    }
    
    Bool_t isV0prNotMC = (v0tagAllCharges == 16) && !fMC; // Only accept V0 protons for data, not for MC
    Bool_t isV0piNotMC = (v0tagAllCharges == 15) && !fMC; // Only accept V0 pions for data, not for MC
    
    Bool_t isV0NotMC = isV0prNotMC || isV0piNotMC;
    
    // Apply track cut
    if(!isV0NotMC && fTrackFilter && !fTrackFilter->IsSelected(track))
      continue;
    
    // Note: For V0's, the cut on ncl can be done via the tree (value is stored). One should not cut on geometry
    // for V0's since they have different topology
    if (!isV0NotMC && GetUseTPCCutMIGeo()) {
      if (!TPCCutMIGeo(track, InputEvent()))
        continue;
    }
    
    
    fPtpc = track->GetTPCmomentum();

    if (fPtpc > 5.)
      continue;
    
    fPt = track->Pt();
    fPhiPrime = GetPhiPrime(track->Phi(), magField, track->Charge());
    
    Int_t label = track->GetLabel();

    AliMCParticle* mcTrack = 0x0;
    
    if (fMC) {
      if (label < 0)
        continue; // If MC is available, reject tracks with negative ESD label
      mcTrack = dynamic_cast<AliMCParticle*>(fMC->GetTrack(TMath::Abs(label)));
      if (!mcTrack) {
        Printf("ERROR: Could not receive mcTrack with label %d for track %d", label, iTracks);
        continue;
      }
      
      /*// Only accept MC primaries
      if (!fMC->Stack()->IsPhysicalPrimary(TMath::Abs(label))) {
        continue;
      }*/
    }

    if (fMC) {
      Int_t pdg = mcTrack->PdgCode();
      
      if (TMath::Abs(pdg) == 2212) // Proton
        isPr = kTRUE;
      else if ((pdg == 111 || TMath::Abs(pdg) == 211) && fPtpc <= fPtpcPionCut) // Pion below desired momentum threshold
        isPi = kTRUE;
      else
        continue; // Only take protons and pions
      
      fPIDtype = kMCid;
      /*
      if (pdg == 111 || TMath::Abs(pdg) == 211) {//Pion
        binMC = 3;
      }
      else if (TMath::Abs(pdg) == 311 || TMath::Abs(pdg) == 321) {//Kaon
        binMC = 1;
      }
      else if (TMath::Abs(pdg) == 2212) {//Proton
        binMC = 4;
      }
      else if (TMath::Abs(pdg) == 11) {//Electron
        binMC = 0;
      }
      else if (TMath::Abs(pdg) == 13) {//Muon
        binMC = 2;
      }*/
    }
    
    const Double_t pileupCorrection = fPIDResponse->GetTPCResponse().GetPileupCorrectionValue(track);
    fDeDx = track->GetTPCsignal() - pileupCorrection;
    
    UInt_t status = track->GetStatus();
    Bool_t hasTOFout  = status&AliESDtrack::kTOFout; 
    Bool_t hasTOFtime = status&AliESDtrack::kTIME;
    Bool_t hasTOF     = kFALSE;
    if (hasTOFout && hasTOFtime)
      hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength();
    // Length check only for primaries, not for V0's!
    if (length < 350 && !isV0NotMC)
      hasTOF = kFALSE;
      
    if (!fMC) {    
      // Note: Normally, the track cuts include a cut on primaries. Therefore, if the particle is a V0, it would usually
      // NOT be found via the primary cuts. This means that the PIDtype describes more or less disjoint sets
      if (isV0prNotMC) {
        // V0
        isPr = kTRUE;
        
        if (!hasTOF) {
          fPIDtype = kV0idNoTOF;
        }
        else {
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) < 3.0)
            fPIDtype = kV0idPlusTOFaccepted;
          else 
            fPIDtype = kV0idPlusTOFrejected;
        }
      }
      else if (isV0piNotMC) {
        // V0 pion
        if (fPtpc > fPtpcPionCut || fDeDx > 140.) // Reject pions above desired momentum threshold and also reject too high dEdx
          continue;
        
        isPi = kTRUE;
        
        if (!hasTOF) {
          fPIDtype = kV0idNoTOF;
        }
        else {
          if (TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion)) < 3.0)
            fPIDtype = kV0idPlusTOFaccepted;
          else  
            fPIDtype = kV0idPlusTOFrejected;
        }
      }
      else { 
        // Non-V0
        isPr = kTRUE; // If particle is accepted, it is a proton (for pions, only V0s are used)
        
        if (fPtpc < 4.0 && //TODO was 2.7 // Do not accept non-V0s above this threshold -> High contamination!!!
            (fDeDx >= 50. / TMath::Power(fPtpc, 1.3))) {// Pattern recognition instead of TPC cut to be ~independent of old TPC expected dEdx
          if (fPtpc < 0.6) {
            fPIDtype = kTPCid;
          }
          // fPtpc >= 0.6
          else if (hasTOF && TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton)) < 3.0) {
              fPIDtype = kTPCandTOFid;
          }
          else
            continue; // Reject particle
        }
        else
          continue; // Reject particle
      }
    }
    
    if (fDoAdditionalQA) {
      Double_t tofTime = track->GetTOFsignal();//in ps
      Double_t tof = tofTime * 1E-3; // ns, average T0 fill subtracted, no info from T0detector   
      Double_t length2 = track->GetIntegratedLength();
      Double_t c = TMath::C() * 1.E-9;// m/ns
      length2 = length2 * 0.01; // in meters
      tof = tof * c;
      Double_t beta = length2 / tof;
      
      Double_t nSigmaTOF = hasTOF ? fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton) : 999;
      
      // p_vert, p_TPC, eta, nSigmaTOF, beta, multiplicity
      Double_t entry[6] = { track->P(), fPtpc, track->Eta(), nSigmaTOF, beta, (Double_t)fMultiplicity };
      fhTOFqa->Fill(entry);
    }
    
    // Prepare entry for tree (some quantities have already been set)
    // Turn eta correction off -> Important here is the pure spline value and the selection via cuts. The correction
    // can be re-done manually, if needed. But it cannot be undone!
    
    if (fCorrectdEdxEtaDependence && fNumEtaCorrReqErrorsIssued < 23 && !etaCorrAvail) {
      AliError("TPC eta correction requested, but was not activated in PID response (most likely not available)!");
      fNumEtaCorrReqErrorsIssued++;
      if (fNumEtaCorrReqErrorsIssued == 23)
        AliError("Ignoring this error from now on!");
    }
    
    if (fCorrectdEdxMultiplicityDependence && fNumMultCorrReqErrorsIssued < 23 && !multCorrAvail) {
      AliError("TPC multiplicity correction requested, but was not activated in PID response (most likely not available)!");
      fNumMultCorrReqErrorsIssued++;
      if (fNumMultCorrReqErrorsIssued == 23)
        AliError("Ignoring this error from now on!");
    }
    
    if (isPr)
      fDeDxExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, 
                                                                       fCorrectdEdxEtaDependence && etaCorrAvail,
                                                                       fCorrectdEdxMultiplicityDependence && multCorrAvail,
                                                                       kFALSE // don't do pileup correction for the expected signal, since above we correct the track dEdx
                                                                      );
    else if (isPi)
      fDeDxExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, 
                                                                       fCorrectdEdxEtaDependence && etaCorrAvail,
                                                                       fCorrectdEdxMultiplicityDependence && multCorrAvail,
                                                                       kFALSE // don't do pileup correction for the expected signal, since above we correct the track dEdx
                                                                      );
    else
      fDeDxExpected = -1.;
    fTanTheta = track->GetInnerParam()->GetTgl();
    //fSinAlpha = track->GetInnerParam()->GetSnp();
    //fY = track->GetInnerParam()->GetY();
    fTPCsignalN = track->GetTPCsignalN();
    
    if (fStoreNumOfSubthresholdclusters) {
      const AliTPCdEdxInfo *clsInfo = track->GetTPCdEdxInfo();
      
      if (clsInfo) {
        Double32_t signal[4] = {0}; // 0: IROC only; 1: OROC short padas; 2: OROC long pads; 3:OROC combined
        Char_t ncl[3] = {0};        // clusters above threshold; 0: IROC; 1: OROC short; 2: OROC long
        Char_t nrows[3] = {0};      // clusters above and below threshold; 0: IROC; 1: OROC short; 2: OROC long

        clsInfo->GetTPCSignalRegionInfo(signal, ncl, nrows);

        fTPCsignalNsubthreshold = nrows[0] + nrows[1] + nrows[2] - ncl[0] - ncl[1] - ncl[2]; 
      }
      else
        fTPCsignalNsubthreshold = 200;// Set to invalid value
    }
    
    
    if (fStoreNumClustersInActiveVolume) 
      fNumTPCClustersInActiveVolume = track->GetLengthInActiveZone(1, 1.8, 220, magField);
    else
      fNumTPCClustersInActiveVolume = 200.;//Set to invalid value
 
    
    
    /* START
    const Double_t tanTheta = track->GetInnerParam()->GetTgl();

    // Constant in formula for B in kGauss (factor 0.1 to convert B from Tesla to kGauss),
    // pT in GeV/c (factor c*1e-9 to arrive at GeV/c) and curvature in 1/cm (factor 0.01 to get from m to cm)
    const Double_t constant = TMath::C()* 1e-9 * 0.1 * 0.01;
    const Double_t curvature = magField * constant / track->Pt(); // in 1./cm
    
    //const Double_t angleIROC = TMath::ASin(TMath::Min(TMath::Abs(par.GetPadRowRadii(0,  par.GetNRow(0)  / 2.) * curvature) * 0.5, 1.));
    //const Double_t angleOROC = TMath::ASin(TMath::Min(TMath::Abs(par.GetPadRowRadii(36, par.GetNRow(36) / 2.) * curvature) * 0.5, 1.));
    
    Double_t averageddzdr = 0.;
    Int_t nParts = 0;

    for (Double_t r = 85; r < 245; r++) {
      Double_t sinPhiLocal = TMath::Abs(r*curvature*0.5);
      
      // Cut on |sin(phi)| as during reco
      if (TMath::Abs(sinPhiLocal) <= 0.95) {
        const Double_t phiLocal = TMath::ASin(sinPhiLocal);
        const Double_t tanPhiLocal = TMath::Tan(phiLocal);
        
        averageddzdr += tanTheta * TMath::Sqrt(1. + tanPhiLocal * tanPhiLocal); 
        nParts++;
      }
    }
    
    if (nParts > 0)
      averageddzdr /= nParts; 
    else {
      AliError("Problems during determination of dz/dr. Skipping track!");
      continue;
    }
    

    //printf("padRow 0 / 36: %f / %f\n", par.GetPadRowRadii(0,  par.GetNRow(0) / 2.), par.GetPadRowRadii(36, par.GetNRow(36) / 2.));
    //printf("pT: %f\nFactor/magField(kGs)/curvature^-1: %f / %f /%f\nIROC/OROC/averageOld/average: %f / %f / %f / //%f\ntanThetaGlobalFromTheta/tanTheta/dzdr: %f / %f / %f\n\n",
    //        track->Pt(), constant, magField, 1./curvature, angleIROC, angleOROC, angleAverage, averageAngle, TMath::Tan(-track->Theta() + TMath::Pi() / 2.0), tanTheta, dzdr);
    //printf("pT: %f\nFactor/magField(kGs)/curvature^-1: %f / %f /%f\ntanThetaGlobalFromTheta/tanTheta/Averageddzdr: %f / %f / %f\n\n",
    //        track->Pt(), constant, magField, 1./curvature, TMath::Tan(-track->Theta() + TMath::Pi() / 2.0), tanTheta, averageddzdr);

  
    fTanTheta = averageddzdr;
    */
    
    if (isPr)
      fTree->Fill();
    else if (isPi)
      fTreePions->Fill();
  } //track loop 

  // Post output data.
  PostData(1, fTree);
  PostData(2, fTreePions);

  if (fDoAdditionalQA)
    PostData(3, fOutputContainer);
  
  // Clear the V0 PID arrays
  ClearV0PIDlist();
}      

//________________________________________________________________________
void AliTPCPIDEtaTree::Terminate(const Option_t *)
{
  // Called once at the end of the query
}
