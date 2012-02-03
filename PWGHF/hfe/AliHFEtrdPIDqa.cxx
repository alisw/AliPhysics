/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// QA class for TRD PID
// Plot Pion Efficiency at x electron efficiency
// Calculate the threshold parametrisation and save
// them in a root file
//
// Author:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TAxis.h>
#include <TBrowser.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>
#include <TIterator.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>

#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliESDtrack.h"
#include "AliHFEtrdPIDqa.h"
#include "AliHFEtools.h"
#include "AliHFEpidTRD.h"
#include "AliLog.h"

ClassImp(AliHFEtrdPIDqa)

const Double_t AliHFEtrdPIDqa::fgkElectronEff[kNElectronEffs] = {0.7,0.75, 0.8, 0.85, 0.9, 0.95};
//_______________________________________________________________
// Definition of the common binning 
const Int_t AliHFEtrdPIDqa::fgkNBinsCommon[kQuantitiesCommon] = {
  AliPID::kSPECIES + 1,         // species
  40,                           // p-bins
  AliESDtrack::kTRDnPlanes + 1 // tracklets including 0
};
const Double_t AliHFEtrdPIDqa::fgkMinBinCommon[kQuantitiesCommon] = {
  -1,      // species
  0.1,     // p-bins:
  0        // tracklets including 0
};

const Double_t AliHFEtrdPIDqa::fgkMaxBinCommon[kQuantitiesCommon] = {
  AliPID::kSPECIES,               // species
  10.,                            // p-bins
  AliESDtrack::kTRDnPlanes + 1    // tracklets including 0
};
//_______________________________________________________________

//__________________________________________________________________
AliHFEtrdPIDqa::AliHFEtrdPIDqa():
  TNamed("trdPIDqa", ""),
  fTRDpid(NULL),
  fHistos(NULL),
  fPionEfficiencies(NULL),
  fProtonEfficiencies(NULL),
  fKaonEfficiencies(NULL),
  fThresholds(NULL),
  fShowMessage(kFALSE),
  fTotalChargeInSlice0(kFALSE),
  fCentralityBin(-1)
{
  //
  // Default Constructor
  //
}

//__________________________________________________________________
AliHFEtrdPIDqa::AliHFEtrdPIDqa(const Char_t *name):
  TNamed(name, ""),
  fTRDpid(NULL),
  fHistos(NULL),
  fPionEfficiencies(NULL),
  fProtonEfficiencies(NULL),
  fKaonEfficiencies(NULL),
  fThresholds(NULL),
  fShowMessage(kFALSE),
  fTotalChargeInSlice0(kFALSE),
  fCentralityBin(-1)
{
  //
  // Main Constructor
  //
}

//__________________________________________________________________
AliHFEtrdPIDqa::AliHFEtrdPIDqa(const AliHFEtrdPIDqa &ref):
  TNamed(ref),
  fTRDpid(NULL),
  fHistos(NULL),
  fPionEfficiencies(NULL),
  fProtonEfficiencies(NULL),
  fKaonEfficiencies(NULL),
  fThresholds(NULL),
  fShowMessage(kFALSE),
  fTotalChargeInSlice0(ref.fTotalChargeInSlice0),
  fCentralityBin(ref.fCentralityBin)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//__________________________________________________________________
AliHFEtrdPIDqa &AliHFEtrdPIDqa::operator=(const AliHFEtrdPIDqa &ref){
  //
  // Assignment operator
  //
  if(this != &ref)
    ref.Copy(*this);
  return *this;
}

//__________________________________________________________________
AliHFEtrdPIDqa::~AliHFEtrdPIDqa(){
  //
  // Destructor
  //
  if(fTRDpid) delete fTRDpid;
  if(fHistos) delete fHistos;
  if(fPionEfficiencies) delete fPionEfficiencies;
  if(fProtonEfficiencies) delete fProtonEfficiencies;
  if(fKaonEfficiencies) delete fKaonEfficiencies;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::Copy(TObject &ref) const{
  //
  // Copies content of this object into object ref
  //
  TNamed::Copy(ref);

  AliHFEtrdPIDqa &target = dynamic_cast<AliHFEtrdPIDqa &>(ref);
  target.fTRDpid = fTRDpid;
  target.fHistos = dynamic_cast<AliHFEcollection *>(fHistos->Clone());  
  target.fTotalChargeInSlice0 = fTotalChargeInSlice0;
  target.fCentralityBin = fCentralityBin;
}

//__________________________________________________________________
Long64_t AliHFEtrdPIDqa::Merge(TCollection *coll){
  //
  // Merge objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty()) return 1;
  
  AliHFEtrdPIDqa *refQA = NULL;
  TIter it(coll);
  TObject *o = NULL;
  Long64_t count = 0;
  TList listHistos;
  while((o = it())){
    refQA = dynamic_cast<AliHFEtrdPIDqa *>(o);
    if(!refQA) continue;

    listHistos.Add(refQA->fHistos);
    count++; 
  }
  fHistos->Merge(&listHistos);
  return count+1;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::Browse(TBrowser *b){
  //
  // Enable Browser functionality
  //
  if(b){
    // Add objects to the browser
    if(fHistos) b->Add(fHistos, fHistos->GetName());
    if(fPionEfficiencies) b->Add(fPionEfficiencies, "Pion Efficiencies");
    if(fProtonEfficiencies) b->Add(fProtonEfficiencies, "Proton Efficiencies");  
    if(fKaonEfficiencies) b->Add(fKaonEfficiencies, "Kaon Efficiencies");
    if(fThresholds) b->Add(fThresholds, "Thresholds");
  }
}

//__________________________________________________________________
void AliHFEtrdPIDqa::Init(){
  //
  // Initialize Object
  //
  
  fHistos = new AliHFEcollection("TRDqa", "Histos for TRD QA");

  CreateLikelihoodHistogram();
  CreateQAHistogram();
  CreatedEdxHistogram();
  CreateHistoTruncatedMean();

  fTRDpid = new AliHFEpidTRD("QAtrdPID");
}

//__________________________________________________________________
void AliHFEtrdPIDqa::CreateLikelihoodHistogram(){
  //
  // Create Histogram for TRD Likelihood Studies
  //
  Int_t nbins[kQuantitiesLike]; memcpy(nbins, fgkNBinsCommon, sizeof(Int_t) * kQuantitiesCommon);
  nbins[kElectronLike] = 100;
  nbins[kNClustersLike] = 200;
  nbins[kCentralityBin] = 12;
  Double_t binMin[kQuantitiesLike]; memcpy(binMin, fgkMinBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  Double_t binMax[kQuantitiesLike]; memcpy(binMax, fgkMaxBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMin[kElectronLike] = 0.;      
  binMin[kNClustersLike] = 0.;
  binMin[kCentralityBin] = -1.;
  binMax[kElectronLike] = 1.;
  binMax[kNClustersLike] = 200.;
  binMax[kCentralityBin] = 11.;

  fHistos->CreateTHnSparse("fLikeTRD","TRD Likelihood Studies", kQuantitiesLike, nbins, binMin, binMax);
  fHistos->BinLogAxis("fLikeTRD", kP);
}

//__________________________________________________________________
void AliHFEtrdPIDqa::CreateQAHistogram(){
  //
  // Create Histogram for Basic TRD PID QA
  //
  AliDebug(1, "Called");
  Int_t nbins[kQuantitiesQA]; memcpy(nbins, fgkNBinsCommon, sizeof(Int_t) * kQuantitiesCommon);
  nbins[kNonZeroTrackletCharge] = AliESDtrack::kTRDnPlanes + 1;
  nbins[kNClusters] = 200;
  Double_t binMin[kQuantitiesQA]; memcpy(binMin, fgkMinBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMin[kNonZeroTrackletCharge] = 0.;      
  binMin[kNClusters] = 0.;      
  Double_t binMax[kQuantitiesQA]; memcpy(binMax, fgkMaxBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMax[kNonZeroTrackletCharge] = AliESDtrack::kTRDnPlanes + 1.;
  binMax[kNClusters] = 200.;

  fHistos->CreateTHnSparse("fQAtrack","TRD QA Histogram", kQuantitiesQA, nbins, binMin, binMax);
  fHistos->BinLogAxis("fQAtrack", kP);
}

//__________________________________________________________________
void AliHFEtrdPIDqa::CreatedEdxHistogram(){
  //
  // Create QA histogram for dEdx investigations
  //
  AliDebug(1, "Called");
  Int_t nbins[kQuantitiesdEdx]; memcpy(nbins, fgkNBinsCommon, sizeof(Int_t) * kQuantitiesCommon);
  nbins[kdEdx] = 100;
  nbins[kNclusters] = 261;
  nbins[kNonZeroSlices] = 9; 
  Double_t binMin[kQuantitiesdEdx]; memcpy(binMin, fgkMinBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMin[kdEdx] = 0.;     
  binMin[kNclusters] = 0;
  binMin[kNonZeroSlices] = 0.;
  Double_t binMax[kQuantitiesdEdx]; memcpy(binMax, fgkMaxBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMax[kdEdx] = 10000.;
  binMax[kNclusters] = 260.;
  binMax[kNonZeroSlices] = 8.;

  fHistos->CreateTHnSparse("fQAdEdx","TRD summed dEdx", kQuantitiesdEdx, nbins, binMin, binMax);
  fHistos->BinLogAxis("fQAdEdx", kP);
  fHistos->Sumw2("fQAdEdx");
}

//__________________________________________________________________
void AliHFEtrdPIDqa::CreateHistoTruncatedMean(){
  //
  // Create Histogram for Basic TRD PID QA:
  //
  AliDebug(1, "Called");
  Int_t nbins[kQuantitiesTruncMean]; memcpy(nbins, fgkNBinsCommon, sizeof(Int_t) * kQuantitiesCommon);
  nbins[kTPCdEdx] = 600;
  nbins[kTRDdEdxMethod1] = 1000;
  nbins[kTRDdEdxMethod2] = 1000;
  Double_t binMin[kQuantitiesTruncMean]; memcpy(binMin, fgkMinBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMin[kTPCdEdx] = 0.;      
  binMin[kTRDdEdxMethod1] = 0.;      
  binMin[kTRDdEdxMethod2] = 0.;      
  Double_t binMax[kQuantitiesTruncMean]; memcpy(binMax, fgkMaxBinCommon, sizeof(Double_t) * kQuantitiesCommon);
  binMax[kTPCdEdx] = 600;
  binMax[kTRDdEdxMethod1] = 20000.;
  binMax[kTRDdEdxMethod2] = 20000.;

  fHistos->CreateTHnSparse("fTRDtruncMean","TRD TruncatedMean studies", kQuantitiesTruncMean, nbins, binMin, binMax);
  fHistos->BinLogAxis("fTRDtruncMean", kP);
  fHistos->CreateTH2F("fTRDslicesPions","TRD dEdx per slice for Pions", 8, 0, 8, 2000, 0, 8000);
  fHistos->CreateTH2F("fTRDslicesElectrons","TRD dEdx per slice for Electrons", 8, 0, 8, 2000, 0, 8000);
}


//__________________________________________________________________
void AliHFEtrdPIDqa::ProcessTracks(const TObjArray * const tracks, Int_t species){
  //
  // Process Electron Tracks
  //
  if(species < -1 || species >= AliPID::kSPECIES) return;
  TIter it(tracks);
  const AliVTrack *track = NULL;
  while((track = dynamic_cast<const AliVTrack *>(it()))){
    if(track) ProcessTrack(track, species);
  }
}

//__________________________________________________________________
void AliHFEtrdPIDqa::ProcessTrack(const AliVTrack * const track, Int_t species){
  //
  // Process Single Track
  //
  if(TString(track->IsA()->GetName()).CompareTo("AliESDtrack") == 0)
    ProcessTrackESD(dynamic_cast<const AliESDtrack *>(track), species);
  else if(TString(track->IsA()->GetName()).CompareTo("AliAODTrack") == 0)
    ProcessTrackAOD(dynamic_cast<const AliAODTrack *>(track), species);
}


//__________________________________________________________________
void AliHFEtrdPIDqa::ProcessTrackESD(const AliESDtrack *track, Int_t species){
  //
  // Process single ESD track
  //
  if(!track) return;
  if((track->GetStatus() & AliESDtrack::kTRDout) == 0) return;  // require TRD track
  FillTRDLikelihoods(track, species);
  FillTRDQAplots(track, species);
}

//__________________________________________________________________
void AliHFEtrdPIDqa::ProcessTrackAOD(const AliAODTrack * const track, Int_t /*species*/){
  //
  // Process single AOD track
  // AOD PID object required
  //
  if(!track) return;
  AliAODPid *trackPID = track->GetDetPid();
  if(!trackPID) return;

}

//__________________________________________________________________
void AliHFEtrdPIDqa::FillTRDLikelihoods(const AliESDtrack * const track, Int_t species){
  //
  // Fill TRD Likelihood Histogram
  //
  Double_t trdLike[AliPID::kSPECIES];
  track->GetTRDpid(trdLike);
  // Renormalize 
  Double_t norm =trdLike[AliPID::kElectron]+trdLike[AliPID::kPion];
  Double_t likeEle = norm == 0. ? 0. : trdLike[AliPID::kElectron]/norm;
  const AliExternalTrackParam *outerPars = track->GetOuterParam();

  //Int_t kQuantitiesLike;
  Double_t quantities[kQuantitiesLike]; memset(quantities, 0, sizeof(Double_t) * kQuantitiesLike);
  // we store:
  // species
  // p
  // ntracklets
  // Electron Likelihood
  quantities[kSpecies] = species;
  quantities[kP] = outerPars ? outerPars->P() : track->P();
  quantities[kNTracklets] = track->GetTRDntrackletsPID();
  quantities[kElectronLike] = likeEle;
  quantities[kNClustersLike] =  track->GetTRDncls();
  quantities[kCentralityBin] = fCentralityBin;
  fHistos->Fill("fLikeTRD", quantities);

}

//__________________________________________________________________
void AliHFEtrdPIDqa::FillTRDQAplots(const AliESDtrack * const track, Int_t species){
  //
  // Fill QA Plots containing further information
  //
  const AliExternalTrackParam *outerPars = track->GetOuterParam();

  Double_t quantitiesQA[kQuantitiesQA], quantitiesdEdx[kQuantitiesdEdx], quantitiesTruncMean[kQuantitiesTruncMean];
  // we store:
  // 1. QA
  // species
  // p
  // ntracklets
  // Non-zero tracklet charges
  // Number of clusters / full track
  // 2. dEdx
  // species
  // p
  // ntracklets
  // dEdx
  // 3. Truncated Mean
  // ...
  // TPC dEdx
  // TRD dEdx Method 1
  // TRD dEdx Method 2
  quantitiesQA[kSpecies]  = quantitiesdEdx[kSpecies] 
                          = quantitiesTruncMean[kSpecies] 
                          = species;
  quantitiesQA[kP]  = quantitiesTruncMean[kP]
                    = outerPars ? outerPars->P() : track->P();
  quantitiesQA[kNTracklets] = quantitiesdEdx[kNTracklets] 
                            = quantitiesTruncMean[kNTracklets]
                            = track->GetTRDntrackletsPID();
  quantitiesQA[kNClusters] = quantitiesdEdx[kNclusters] = track->GetTRDncls();
  

  Double_t dEdxSum = 0., qSlice = 0.;
  // remove the last slice from the histogram
  Int_t ntrackletsNonZero = 0, nSlices = track->GetNumberOfTRDslices(), nSlicesNonZero = 0;
  TString speciesname = "pions";
  Bool_t selectedForSlicemon = kFALSE;
  
  switch(species){
    case AliPID::kElectron: speciesname = "Electrons"; selectedForSlicemon = kTRUE; break;
    case AliPID::kPion: speciesname = "Pions"; selectedForSlicemon = kTRUE; break;
    default: speciesname = "undefined"; selectedForSlicemon = kFALSE; break;
  };
  AliDebug(1, Form("species %d, speciesname %s, momentum %f, selected %s", species, speciesname.Data(), track->P(), selectedForSlicemon ? "yes" : "no"));
  for(Int_t iplane = 0; iplane < AliESDtrack::kTRDnPlanes; iplane++){
    quantitiesdEdx[kP] = track->GetTRDmomentum(iplane);
    dEdxSum = 0.;
    for(Int_t islice = 0; islice < nSlices; islice++){
      if(fTotalChargeInSlice0 && islice >= 7) break;
      qSlice = track->GetTRDslice(iplane, fTotalChargeInSlice0 ? islice + 1 : islice);  // hack by mfasel: For data with the new reconstruction, slice 0 is used to store the total charge, the total number of slices is 7 instead of 8
      if(qSlice > 1e-1){
        // cut out 0 slices
        nSlicesNonZero++;
        dEdxSum += qSlice;
        // Reweighting of the slices for the truncated mean: select all pion tracks above
        // 1.5 GeV and monitor the dEdx as function of slice
        if(selectedForSlicemon && track->P() > 1.5){
          AliDebug(2, Form("Filling Histogram fTRDslices%s", speciesname.Data()));
          fHistos->Fill(Form("fTRDslices%s", speciesname.Data()), static_cast<Double_t>(islice), qSlice);
        }
      }
    }
    quantitiesdEdx[kNonZeroSlices] = nSlicesNonZero;
    quantitiesdEdx[kdEdx] = fTotalChargeInSlice0 ? track->GetTRDslice(iplane, 0) : dEdxSum; // hack by mfasel: In the new reconstruction, the total charge is stored in the first slice, in the old reconstruction it has to be calculated from the slice charges.     
    if(dEdxSum) ntrackletsNonZero++;
    // Fill dEdx histogram
    if(dEdxSum > 1e-1) fHistos->Fill("fQAdEdx", quantitiesdEdx); // Cut out 0 entries
  }
  quantitiesQA[kNonZeroTrackletCharge] = ntrackletsNonZero;
  fHistos->Fill("fQAtrack", quantitiesQA);

  quantitiesTruncMean[kTPCdEdx] = track->GetTPCsignal();
  quantitiesTruncMean[kTRDdEdxMethod1] = fTRDpid->GetTRDSignalV1(track, 0.6);
  quantitiesTruncMean[kTRDdEdxMethod2] = fTRDpid->GetTRDSignalV2(track, 0.6);
  fHistos->Fill("fTRDtruncMean", quantitiesTruncMean);
}

/////////////////////////////////////////////////////////
//
// Code for Post Processing
//
// //////////////////////////////////////////////////////

//__________________________________________________________________
void AliHFEtrdPIDqa::FinishAnalysis(Int_t fCentrality, Bool_t isGreaterEqual){
  //
  // Finish Analysis:
  // Calculate Electron Efficiency for ntracklets = 4...6
  // Calculate thresholds for ntracklets = 4...6
  //
  
  if(!fPionEfficiencies){ 
    fPionEfficiencies = new TList;
    fPionEfficiencies->SetName("pionEfficiencies");
    fPionEfficiencies->SetOwner();
  }
  if(!fProtonEfficiencies){
    fProtonEfficiencies = new TList;
    fProtonEfficiencies->SetName("protonEfficiencies");
    fProtonEfficiencies->SetOwner();
  }
  if(!fThresholds){
    fThresholds = new TList;
    fThresholds->SetName("thresholds");
    fThresholds->SetOwner();
  }

  for(Int_t itr = 4; itr <= 6; itr++){
    if(fShowMessage){
      printf("========================================\n");
      printf("Analysing %d trackltes centrality %i \n", itr, fCentrality);
      printf("========================================\n");
    }
    AnalyseNTracklets(itr, fCentrality, isGreaterEqual);
  }
}

//__________________________________________________________________
void AliHFEtrdPIDqa::StoreResults(const Char_t *filename){
  //
  // Store histos into a root file
  //
  TFile *outfile = new TFile(filename, "RECREATE");
  outfile->cd();
  fPionEfficiencies->Clone()->Write(fPionEfficiencies->GetName(), kSingleKey);
  fProtonEfficiencies->Clone()->Write(fProtonEfficiencies->GetName(), kSingleKey);
  fThresholds->Clone()->Write(fThresholds->GetName(), kSingleKey);
  outfile->Close();
  delete outfile;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::SaveThresholdParameters(const Char_t *name, Double_t lowerLimit, Double_t upperLimit, Int_t icentrality){
  //
  // Fit the threshold histograms with the given parametrisation
  // and store the TF1 in the file
  //

  if(!fThresholds){
    AliError("Threshold histograms have to be created first");
    return;
  }

    if(fShowMessage){
    printf("========================================\n");
    printf("Calculating threshold parameters\n");
    printf("========================================\n");
  }

  TList *outlist = new TList;
  outlist->SetName("thresholdTRD");
  
  TGraph *threshhist = NULL;
  TF1 *threshparam = NULL;
  TList *lHistos = NULL, *lFormulas = NULL;
  for(Int_t itracklet = 4; itracklet <= 6; itracklet++){
  
    if(fShowMessage){
      printf("-------------------------------\n");
      printf("Processing %d tracklets\n", itracklet);
      printf("-------------------------------\n");
    }

    Char_t *listname=Form("%dTracklets", itracklet);
    if(icentrality!=-1) listname=Form("%dTracklets%dCentrality", itracklet, icentrality);
    lHistos = dynamic_cast<TList *>(fThresholds->FindObject(listname));
    if(!lHistos){
      AliError(Form("Threshold histograms for the case %s not found", listname));
      continue;
    }
    lFormulas = new TList;
    lFormulas->SetName(listname);
    outlist->Add(lFormulas);
    
    for(Int_t ieff = 0; ieff <  kNElectronEffs; ieff++){
      
      if(fShowMessage){
        printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
        printf("Processing Electron Efficiency %f\n", fgkElectronEff[ieff]);
        printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
      }

      threshhist = dynamic_cast<TGraph *>(lHistos->FindObject(Form("eff%d", static_cast<Int_t>(fgkElectronEff[ieff] * 100))));
      if(!threshhist) continue;
      threshparam = MakeThresholds(threshhist, lowerLimit, upperLimit);
      threshparam->SetName(Form("thresh_%d_%d", itracklet, static_cast<Int_t>(fgkElectronEff[ieff] * 100)));
      if(icentrality!=-1) threshparam->SetName(Form("thresh_%d_%d_%d", itracklet, icentrality, static_cast<Int_t>(fgkElectronEff[ieff] * 100)));
      lFormulas->Add(threshparam);
    }
  }

  // store the output
  TFile *outfile = new TFile(name, "RECREATE");
  outfile->cd();
  outlist->Write(outlist->GetName(), kSingleKey);
  outfile->Close();
  delete outfile;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::AnalyseNTracklets(Int_t nTracklets, Int_t nCentrality, Bool_t isGreaterEqual){
  //
  // Calculate Pion Efficiency, Proton Efficiency and Kaon Efficiency at discrete
  // elPion Efficiency, Proton Efficiency and Kaon Efficiency at discrete
  // electron efficiencies
  //
  THnSparse *hLikeTRD = dynamic_cast<THnSparseF *>(fHistos->Get("fLikeTRD"));
  if(!hLikeTRD){
    AliError("Likelihood Histogram not available");
    return;
  }

  Bool_t isPbPb = kFALSE;
  if(nCentrality==-1)  isPbPb = kFALSE;
  if(nCentrality!=-1)  isPbPb = kTRUE;

  Int_t binTracklets = hLikeTRD->GetAxis(kNTracklets)->FindBin(nTracklets);
  hLikeTRD->GetAxis(kNTracklets)->SetRange(binTracklets, isGreaterEqual ? 7 : binTracklets);

  if(isPbPb){
      Int_t binCentrality = hLikeTRD->GetAxis(kCentralityBin)->FindBin(nCentrality);
      hLikeTRD->GetAxis(kCentralityBin)->SetRange(binCentrality, isGreaterEqual ? 11 : binCentrality);
      /*
       new TCanvas;
       TH2 *test = hLikeTRD->Projection(kCentralityBin, kP);
       test->Draw("colz");
       */
  }

  Int_t binElectrons = hLikeTRD->GetAxis(kSpecies)->FindBin(AliPID::kElectron); 
  AliDebug(1, Form("BinElectrons %d", binElectrons));
  Int_t binPions = hLikeTRD->GetAxis(kSpecies)->FindBin(AliPID::kPion);
  AliDebug(1, Form("BinPions %d", binPions));
  Int_t binProtons =  hLikeTRD->GetAxis(kSpecies)->FindBin(AliPID::kProton);
  AliDebug(1, Form("BinProtons %d", binProtons));
  hLikeTRD->GetAxis(kSpecies)->SetRange(binElectrons, binElectrons);
  TH2 *likeElectron = hLikeTRD->Projection(kElectronLike, kP);
  likeElectron->SetName("likeElectron");
  hLikeTRD->GetAxis(kSpecies)->SetRange(binPions, binPions);
  TH2 *likePion = hLikeTRD->Projection(kElectronLike, kP);
  likePion->SetName("likePion");
  hLikeTRD->GetAxis(kSpecies)->SetRange(binProtons, binProtons);
  TH2 *likeProton = hLikeTRD->Projection(kElectronLike, kP);
  likeProton->SetName("likeProton");
  
  // Undo ranges
  hLikeTRD->GetAxis(kSpecies)->SetRange(0, hLikeTRD->GetAxis(kSpecies)->GetNbins());
  hLikeTRD->GetAxis(kNTracklets)->SetRange(0, hLikeTRD->GetAxis(kNTracklets)->GetNbins());
  hLikeTRD->GetAxis(kNTracklets)->SetRange(0, hLikeTRD->GetAxis(kCentralityBin)->GetNbins());

  // Prepare List for output
  Char_t *listname=Form("%dTracklets", nTracklets);
  if(isPbPb) listname=Form("%dTracklets%dCentrality", nTracklets, nCentrality);


  TList *listPions = new TList; listPions->SetName(listname); listPions->SetOwner();
  TList *listProtons = new TList; listProtons->SetName(listname); listProtons->SetOwner();
  TList *listThresholds = new TList; listThresholds->SetName(listname); listThresholds->SetOwner();
  fPionEfficiencies->Add(listPions);
  fProtonEfficiencies->Add(listProtons);
  fThresholds->Add(listThresholds);

  TH1 *probsEl = NULL, *probsPi = NULL, *probsPr = NULL;
  TGraphErrors *effPi = NULL, *effPr = NULL, *thresholds = NULL;
  Double_t p = 0, dp = 0;
  Int_t threshbin = 0;
  Double_t eff, error; // value and error
  for(Int_t ieff = 0; ieff < kNElectronEffs; ieff++){
    
    if(fShowMessage){
      printf("-----------------------------------------\n");
      printf("Doing Electron Efficiency %f\n", fgkElectronEff[ieff]);
      printf("-----------------------------------------\n");
    }
    effPi = new TGraphErrors(likeElectron->GetXaxis()->GetNbins());
    effPi->SetName(Form("eff%d", static_cast<Int_t >(fgkElectronEff[ieff] * 100)));
    effPr = new TGraphErrors(likeElectron->GetXaxis()->GetNbins());
    effPr->SetName(Form("eff%d", static_cast<Int_t >(fgkElectronEff[ieff] * 100)));
    thresholds = new TGraphErrors(likeElectron->GetXaxis()->GetNbins());
    thresholds->SetName(Form("eff%d", static_cast<Int_t >(fgkElectronEff[ieff] * 100)));
  
    // Add to lists
    listPions->Add(effPi);
    listProtons->Add(effPr);
    listThresholds->Add(thresholds);

    for(Int_t imom = 1; imom <= likeElectron->GetXaxis()->GetLast(); imom++){
      p = likeElectron->GetXaxis()->GetBinCenter(imom);
      dp = likeElectron->GetXaxis()->GetBinWidth(imom)/2;

      probsEl = likeElectron->ProjectionY("el", imom, imom);
      if(!probsEl->GetEntries()) continue;
      probsPi = likePion->ProjectionY("pi", imom, imom);
      if(!probsPi->GetEntries()) continue;
      probsPr = likeProton->ProjectionY("pr", imom, imom);
      if(!probsPr->GetEntries()) continue;
      AliDebug(1, Form("Calculating Values for p = %f", p));

      // Calculate non-electronEfficiency and error
      eff = CalculateHadronEfficiency(probsPi, probsEl, fgkElectronEff[ieff], threshbin, error);
      thresholds->SetPoint(imom - 1, p, probsEl->GetXaxis()->GetBinCenter(threshbin));
      thresholds->SetPointError(imom - 1, dp, EstimateThresholdError(probsEl, threshbin));
      AliDebug(1, Form("threshold %d|%f", threshbin, probsEl->GetXaxis()->GetBinCenter(threshbin)));
      AliDebug(1, Form("Pion Efficiency %f +- %f", eff, error));
      effPi->SetPoint(imom - 1, p, eff);
      effPi->SetPointError(imom - 1, dp, error);
      eff = CalculateHadronEfficiency(probsPr, probsEl, fgkElectronEff[ieff] , threshbin, error);
      AliDebug(1, Form("Proton Efficiency %f", eff));
      effPr->SetPoint(imom - 1, p, eff);
      effPr->SetPointError(imom - 1, dp, error);
 
      // cleanup
      delete probsEl;
      delete probsPi;
      delete probsPr;
    }
  }

  // remove temporary histograms
  delete likeElectron;
  delete likePion;
  delete likeProton;
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::CalculateHadronEfficiency(const TH1 * const hadron, const TH1 *const electron, Double_t eff, Int_t &threshbin, Double_t &error){
  // 
  // Calculate non-electron efficiency
  // optionally returns sums as second parameter
  //

  TH1D eleWorking(*((const TH1D *)electron)), hadronWorking(*((const TH1D *)hadron)); // Leave the original histograms untouched and do calculation including scale on copy
  eleWorking.Scale(1./eleWorking.Integral());
  hadronWorking.Scale(1./hadronWorking.Integral());

  TArrayD sumsEl(eleWorking.GetNbinsX()), sumsHd(eleWorking.GetNbinsX());

  // calculate threshold and estimated electron efficiency the threshold was taken
  Double_t elEff = 0.;  // estimated electron efficiency at the end
  Int_t currentBin = 0, nbins = 0;
  for(Int_t ibin = eleWorking.GetXaxis()->GetLast(); ibin >= eleWorking.GetXaxis()->GetFirst(); ibin--){
    currentBin = ibin;
    nbins++;
    elEff += eleWorking.GetBinContent(ibin);
    sumsEl[eleWorking.GetXaxis()->GetLast() - ibin] = elEff;
    if(elEff >= eff){
      // we found the matching bin, break the loop
      break;
    }
  }
  threshbin = currentBin;

  Double_t hdEff = 0; 
  for(Int_t ibin = hadronWorking.GetXaxis()->GetLast(); ibin >= threshbin; ibin--) {
    hdEff += hadronWorking.GetBinContent(ibin);
    sumsHd[hadronWorking.GetXaxis()->GetLast() - ibin] = hdEff;
  }

  // search sums of electron efficiency for double counts, eliminate in electron and hadron array
  TArrayD newsumsEl(100), newsumsHd(100);
  Int_t nusable = 0;
  for(Int_t ien = 0; ien < nbins; ien++){
    if(ien==0){
      newsumsEl[0] = sumsEl[0];
      nusable++;
      continue;
    }
    Int_t index = TMath::BinarySearch(nusable, newsumsEl.GetArray(), sumsEl[ien]);
    if(TMath::Abs(sumsEl[ien] - newsumsEl[index]) < 1e-13){
      // element already counted, don't add to the new arrays
      continue; 
    }
    newsumsEl[nusable] = sumsEl[ien];
    newsumsHd[nusable] = sumsHd[ien];
    nusable++;
  }

  //printf("New array\n");
  //for(Int_t ib = 0; ib < nusable; ib++){
  //  printf("Electron Efficiency %f, Pion Efficiency %f\n", newsumsEl[ib], newsumsHd[ib]);
  //}
  //printf("Do Fit\n");

  // Calculate error
  error = 0;
  if(hadronWorking.GetEntries() > 0 && eleWorking.GetEntries() > 0 && nusable > 2){
    // Do error calculation in case the bins have enough statistics
    TGraph gevh(nusable, newsumsEl.GetArray(), newsumsHd.GetArray()); 
    TF1 evh("evh","pol2", eff-.05, eff+.05);
    gevh.Fit(&evh, "Q", "", eff-.05, eff+.05);
  
    // return the error of the pion efficiency
    if(((1.-hdEff) < 0) || ((1.- elEff) < 0)){
      AliError(" ElEffi or HdEffi > 1. Error can not be calculated. Please increase statistics!");
    }   else {
      error = TMath::Sqrt(hdEff*(1-hdEff)/hadronWorking.GetEntries()+TMath::Power(evh.Derivative(eff), 2)*elEff*(1-elEff)/eleWorking.GetEntries());
    }
    AliDebug(2, Form("Pion Effi at [%f] : [%f +/- %f], Threshold[%f]", elEff, hdEff, error, eleWorking.GetBinCenter(threshbin)));
    AliDebug(2, Form("Derivative at %4.2f : %f\n", eff, evh.Derivative(eff)));
  }

  return hdEff;
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::EstimateThresholdError(const TH1 * const electron, Int_t threshbin){
  //
  // Estimate threshold error as sqrt(N_int)/N_ent
  // where N_int are the counts integrated under the electron selection band
  // and N_ent are all entries in the histogram
  //
  Int_t nInt = 0;
  for(Int_t ibin = threshbin; ibin <= electron->GetXaxis()->GetNbins(); ibin++){
    nInt += (Int_t) electron->GetBinContent(ibin);
  }
  Int_t nEnt = (Int_t) electron->Integral();
  if(nEnt == 0) return 1.;
  return TMath::Sqrt(static_cast<Double_t>(nInt))/static_cast<Double_t>(nEnt);
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::CalculateIntegratedPionEfficiency(UInt_t nTracklets, Double_t electronEff, Double_t pmin, Double_t pmax, Int_t icentrality, Double_t *error){
  //
  // Calculate Pion Efficiency for a given electron efficiency in the specified momentum range
  //
  if(nTracklets < 4 || nTracklets > 6){
    AliError("Pion Efficiency calculation only available for 4, 5, and 6 tracklets");
    return 0.;
  }
  if(electronEff < 0.6 || electronEff > 1.){
    AliError("Pion Efficiency calculation only available in the electron efficiency range 0.6 to 1");
    return 0.;
  }
  if(pmin < 0.1 || pmin > 10 || pmax < 0.1 || pmax > 10.){
    AliError("Pion Efficiency calculation only available in the momentum range 0.1 to 10 GeV/c");
    return 0.;
  }
  if(pmax < pmin){
    AliError("pmin is expected to be >= pmax");
    return 0.;
  }

  // prerequierements fullfiled
  // prepare histos
  THnSparse *hLikeTRD = dynamic_cast<THnSparseF *>(fHistos->Get("fLikeTRD"));
  if(!hLikeTRD){
    AliError("Likelihood Histogram not available");
    return 0;
  }
  Int_t binTracklets = hLikeTRD->GetAxis(kNTracklets)->FindBin(nTracklets);
  hLikeTRD->GetAxis(kNTracklets)->SetRange(binTracklets, binTracklets);

  if(icentrality!=1){
      Int_t binCentrality = hLikeTRD->GetAxis(kCentralityBin)->FindBin(icentrality);
      hLikeTRD->GetAxis(kCentralityBin)->SetRange(binCentrality, binCentrality);
  }

  Int_t pbinMin = hLikeTRD->GetAxis(kP)->FindBin(pmax),
        pbinMax = hLikeTRD->GetAxis(kP)->FindBin(pmax);
  hLikeTRD->GetAxis(kP)->SetRange(pbinMin, pbinMax);
  Int_t binElectrons = hLikeTRD->GetAxis(kSpecies)->FindBin(AliPID::kElectron); 
  Int_t binPions = hLikeTRD->GetAxis(kSpecies)->FindBin(AliPID::kPion);
  hLikeTRD->GetAxis(kSpecies)->SetRange(binElectrons, binElectrons);
  TH1 *likeElectron = hLikeTRD->Projection(kElectronLike);
  likeElectron->Scale(1./likeElectron->Integral());
  likeElectron->SetName("likeElectron");
  hLikeTRD->GetAxis(kSpecies)->SetRange(binPions, binPions);
  TH1 *likePion = hLikeTRD->Projection(kElectronLike);
  likePion->Scale(1./likePion->Integral());
  likePion->SetName("likePion");
  
  // Undo ranges
  hLikeTRD->GetAxis(kSpecies)->SetRange(0, hLikeTRD->GetAxis(kSpecies)->GetNbins());
  hLikeTRD->GetAxis(kNTracklets)->SetRange(0, hLikeTRD->GetAxis(kNTracklets)->GetNbins());
  hLikeTRD->GetAxis(kP)->SetRange(0, hLikeTRD->GetAxis(kP)->GetNbins());

  // Do Calculation
  Int_t thresh; Double_t err;
  Double_t effpi = CalculateHadronEfficiency(likePion, likeElectron, electronEff, thresh, err);
  delete likePion; delete likeElectron;
  if(error) *error = err;
  return effpi;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::DrawTracklet(Int_t itracklet, Double_t pmin, Double_t pmax, Int_t icentrality, Bool_t doFit){
  //
  // Draw efficiencies and threshold as function of p
  //
  if(!(fPionEfficiencies && fProtonEfficiencies && fThresholds)){
    AliError("No graphs to draw available");  
    return;
  }

  Char_t *listname=Form("%dTracklets", itracklet);
  if(icentrality!=-1) listname=Form("%dTracklets%dCentrality", itracklet, icentrality);

  TList *lpions = dynamic_cast<TList *>(fPionEfficiencies->FindObject(listname));
  TList *lprotons = dynamic_cast<TList *>(fProtonEfficiencies->FindObject(listname));
  
  TList *lthresholds = dynamic_cast<TList *>(fThresholds->FindObject(listname));
  if(!(lpions && lprotons && lthresholds)){
    AliDebug(1, "Relevant lists not found. Did you forget to run FinishAnalysis()?");
    return;
  }

  TGraphErrors *pi, *pr;
  TGraph *tr;
  TLegend *leg;
  Char_t *canvasname=Form("tracklet%d", itracklet);
  if(icentrality!=-1) canvasname=Form("tracklet%dcentrality%d", itracklet, icentrality);
  TCanvas *c1 = new TCanvas(canvasname, canvasname, 1024, 768);
  c1->Divide(3,2);
  TF1 *threshfit = NULL;
  for(Int_t ieff = 0; ieff < kNElectronEffs; ieff++){
    c1->cd(ieff + 1);
    gPad->SetGrid(0,0);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.08);
    pi = dynamic_cast<TGraphErrors *>(lpions->FindObject(Form("eff%d", static_cast<Int_t>(fgkElectronEff[ieff] * 100))));
    pr = dynamic_cast<TGraphErrors *>(lprotons->FindObject(Form("eff%d", static_cast<Int_t>(fgkElectronEff[ieff] * 100))));
    tr = dynamic_cast<TGraph *>(lthresholds->FindObject(Form("eff%d", static_cast<Int_t>(fgkElectronEff[ieff] * 100))));
    if(!(pi && pr && tr)) continue;

    // Define nice plot, draw
    // Axis Title
    pi->GetXaxis()->SetTitle("p / GeV/c");
    pi->GetYaxis()->SetTitle("Efficiency");
    pr->GetXaxis()->SetTitle("p / GeV/c");
    pr->GetYaxis()->SetTitle("Efficiency");
    tr->GetXaxis()->SetTitle("p / GeV/c");
    tr->GetYaxis()->SetTitle("Efficiency");
    pi->GetYaxis()->SetTitleOffset(1.2);
    pr->GetYaxis()->SetTitleOffset(1.2);
    tr->GetYaxis()->SetTitleOffset(1.2);
    pi->GetXaxis()->SetTitleSize(0.045);
    pi->GetYaxis()->SetTitleSize(0.045);
    pr->GetXaxis()->SetTitleSize(0.045);
    pr->GetYaxis()->SetTitleSize(0.045);
    tr->GetXaxis()->SetTitleSize(0.045);
    tr->GetYaxis()->SetTitleSize(0.045);
    // Axis Range
    pi->GetYaxis()->SetRangeUser(1e-3, 1.);
    pr->GetYaxis()->SetRangeUser(1e-3, 1.);
    tr->GetYaxis()->SetRangeUser(1e-3, 1.);
    if(pmin >= 0 && pmax >= 0.){
      pi->GetXaxis()->SetRangeUser(pmin, pmax);
      pr->GetXaxis()->SetRangeUser(pmin, pmax);
      tr->GetXaxis()->SetRangeUser(pmin, pmax);
    }
    // Marker
    pi->SetMarkerColor(kRed);
    pi->SetMarkerStyle(20);
    pr->SetMarkerColor(kBlue);
    pr->SetMarkerStyle(21);
    tr->SetMarkerColor(kBlack);
    tr->SetMarkerStyle(22);
    // Title
    pi->SetTitle(Form ("%.2f Electron Efficiency", fgkElectronEff[ieff]));
    pr->SetTitle(Form ("%.2f Electron Efficiency", fgkElectronEff[ieff]));
    tr->SetTitle(Form ("%.2f Electron Efficiency", fgkElectronEff[ieff]));
    // Draw
    pi->Draw("ape"); pr->Draw("pesame"); tr->Draw("pesame");

    // Optionally do Fit
    if(doFit){
      threshfit = MakeThresholds(tr, pmin, pmax);
      threshfit->SetLineColor(kBlack);
      threshfit->Draw("same");
    }

    // Add entries to legend
    if(ieff==0){
      leg = new TLegend(0.5, 0.65, 0.89, 0.85);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(pi, "Pion Efficiency", "lp");
      leg->AddEntry(pr, "Proton Efficiency", "lp");
      leg->AddEntry(tr, "Thresholds", "lp");
      leg->Draw();
      gPad->Update();
    }
  }
  c1->cd();
}

//__________________________________________________________________
TF1 *AliHFEtrdPIDqa::MakeThresholds(TGraph *threshist, Double_t lowerLimit, Double_t upperLimit){
  //
  // Create TF1 containing the threshold parametrisation
  //

  TF1 *threshparam = new TF1("thresh", "1-[0]-[1]*x-[2]*TMath::Exp(-[3]*x)", 0.1, 10);
  threshist->Fit(threshparam, "NE", "", lowerLimit, upperLimit);
  return threshparam;
}

//__________________________________________________________________
void AliHFEtrdPIDqa::ClearLists(){
  //
  // Clear lists for particle efficiencies and thresholds
  //
  if(fPionEfficiencies){
    fPionEfficiencies->Delete();
    delete fPionEfficiencies;
    fPionEfficiencies = NULL;
  }
  if(fProtonEfficiencies){
    fProtonEfficiencies->Delete();
    delete fProtonEfficiencies;
    fProtonEfficiencies = NULL;
  }
  if(fThresholds){
    fThresholds->Delete();
    delete fThresholds;
    fThresholds = NULL;
  }
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::EvalPionEfficiency(Int_t ntls, Int_t eEff, Double_t p, Int_t icentrality){
  //
  // calculate pion efficiency
  // Arguments:
  //   Number of tracklets
  //   Electron Efficiency
  //   Momentum
  //
  TGraphErrors *measurement = GetPionEfficiency(ntls, eEff, icentrality);
  if(!measurement) return -1.;
  return measurement->Eval(p);
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::EvalProtonEfficiency(Int_t ntls, Int_t eEff, Double_t p, Int_t icentrality){
  //
  // calculate proton efficiency
  // Arguments:
  //   Number of tracklets
  //   Electron Efficiency
  //   Momentum
  //
  TGraphErrors *measurement = GetProtonEfficiency(ntls, eEff, icentrality);
  if(!measurement) return -1.;
  return measurement->Eval(p);
}

//__________________________________________________________________
Double_t AliHFEtrdPIDqa::EvalThreshold(Int_t ntls, Int_t eEff, Double_t p, Int_t icentrality){
  //
  // Get the threshold to reach a certain electron efficiency
  // Arguments:
  //   Number of tracklets
  //   Electron Efficiency
  //   Momentum
  //
  TGraph *measurement = GetThreshold(ntls, eEff, icentrality);
  if(!measurement) return -1.;
  return measurement->Eval(p);
}

//__________________________________________________________________
TGraphErrors *AliHFEtrdPIDqa::GetPionEfficiency(Int_t ntracklets, Int_t eleffpercent, Int_t icentrality){
  //
  // Get Graph with pion efficiencies
  //
  Char_t *listname=Form("%dTracklets", ntracklets);
  if(icentrality!=-1) listname=Form("%dTracklets%dCentrality", ntracklets, icentrality);
  TList *graphs = dynamic_cast<TList *>(fPionEfficiencies->FindObject(listname));
  if(!graphs) return NULL;
  return dynamic_cast<TGraphErrors *>(graphs->FindObject(Form("eff%d", eleffpercent)));
}

//__________________________________________________________________
TGraphErrors *AliHFEtrdPIDqa::GetProtonEfficiency(Int_t ntracklets, Int_t eleffpercent, Int_t icentrality){
  // 
  // Get Graph with proton efficiencies
  //
  Char_t *listname=Form("%dTracklets", ntracklets);
  if(icentrality!=-1) listname=Form("%dTracklets%dCentrality", ntracklets, icentrality);
  TList *graphs = dynamic_cast<TList *>(fProtonEfficiencies->FindObject(listname));
  if(!graphs) return NULL;
  return dynamic_cast<TGraphErrors *>(graphs->FindObject(Form("eff%d", eleffpercent)));
}

//__________________________________________________________________
TGraph *AliHFEtrdPIDqa::GetThreshold(Int_t ntracklets, Int_t eleffpercent, Int_t icentrality){
  //
  // Get Graph with threshols
  //
  Char_t *listname=Form("%dTracklets", ntracklets);
  if(icentrality!=-1) listname=Form("%dTracklets%dCentrality", ntracklets, icentrality);
  TList *graphs = dynamic_cast<TList *>(fThresholds->FindObject(listname));
  if(!graphs) return NULL;
  return dynamic_cast<TGraph *>(graphs->FindObject(Form("eff%d", eleffpercent)));
}

