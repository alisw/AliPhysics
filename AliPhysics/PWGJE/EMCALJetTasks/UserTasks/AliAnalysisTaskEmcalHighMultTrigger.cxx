// $Id$
//
// High multiplicity pp trigger analysis task.
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliLog.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliParticleContainer.h"
#include "AliVVZERO.h"

#include "AliAnalysisTaskEmcalHighMultTrigger.h"

ClassImp(AliAnalysisTaskEmcalHighMultTrigger)

//________________________________________________________________________
AliAnalysisTaskEmcalHighMultTrigger::AliAnalysisTaskEmcalHighMultTrigger() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalHighMultTrigger", kTRUE),
  fNExLP(1),
  fNAccPatches(-1),
  fMedianEnergy(0),
  fMedianEnergyExLP(0),
  fSumEnergy(0),
  fSumEnergyExLP(0),
  fTruncatedMean(0),
  fTruncateThreshold(6.),
  fHistPatchEtaPhiE(0),
  fHistEnergyMedian(0),
  fHistEnergyMedianExLP(0),
  fHistEnergySum(0),
  fHistEnergySumExLP(0),
  fHistTruncatedMean(0),
  fHistTracks(0),
  fHistTracklets(0),
  fHistV0MultSum(0),
  fHistTracksTracklets(0),
  fHistTracksV0MultSum(0),
  fHistSPDTrkClsSum(0),
  fHistSPDTrkClsSumExLP(0),
  fHistSPDTrkClsMedian(0),
  fHistSPDTrkClsMedianExLP(0),
  fHistSPDTrkClsTruncMean(0)
{
  // Default constructor.

  const Int_t nMultEst = 3;
  for(Int_t i = 0; i<nMultEst; i++) {
    fHistEnergyMedianEst[i] = 0;
    fHistEnergyMedianExLPEst[i] = 0;
    fHistEnergySumEst[i] = 0;
    fHistEnergySumExLPEst[i] = 0;
    fHistEnergySumAvgEst[i] = 0;
    fHistEnergySumAvgExLPEst[i] = 0;
    fHistTruncatedMeanEst[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalHighMultTrigger::AliAnalysisTaskEmcalHighMultTrigger(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fNExLP(1),
  fNAccPatches(-1),
  fMedianEnergy(0),
  fMedianEnergyExLP(0),
  fSumEnergy(0),
  fSumEnergyExLP(0),
  fTruncatedMean(0),
  fTruncateThreshold(6.),
  fHistPatchEtaPhiE(0),
  fHistEnergyMedian(0),
  fHistEnergyMedianExLP(0),
  fHistEnergySum(0),
  fHistEnergySumExLP(0),
  fHistTruncatedMean(0),
  fHistTracks(0),
  fHistTracklets(0),
  fHistV0MultSum(0),
  fHistTracksTracklets(0),
  fHistTracksV0MultSum(0),
  fHistSPDTrkClsSum(0),
  fHistSPDTrkClsSumExLP(0),
  fHistSPDTrkClsMedian(0),
  fHistSPDTrkClsMedianExLP(0),
  fHistSPDTrkClsTruncMean(0)
{
  // Standard constructor.

  const Int_t nMultEst = 3;
  for(Int_t i = 0; i<nMultEst; i++) {
    fHistEnergyMedianEst[i] = 0;
    fHistEnergyMedianExLPEst[i] = 0;
    fHistEnergySumEst[i] = 0;
    fHistEnergySumExLPEst[i] = 0;
    fHistEnergySumAvgEst[i] = 0;
    fHistEnergySumAvgExLPEst[i] = 0;
    fHistTruncatedMeanEst[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalHighMultTrigger::~AliAnalysisTaskEmcalHighMultTrigger()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHighMultTrigger::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Int_t nE = 500;
  Double_t minE = 0.;
  Double_t maxE = 50.;

  TString histName = "";
  TString histTitle = "";

  histName = Form("fHistPatchEtaPhiE");
  histTitle = Form("%s;#eta;#phi;E",histName.Data());
  fHistPatchEtaPhiE = new TH3F(histName.Data(),histTitle.Data(),100,-1.,1.,18*8,0.,TMath::TwoPi(),100,0.,100.);
  fOutput->Add(fHistPatchEtaPhiE);

  histName = Form("fHistEnergyMedian");
  histTitle = Form("%s;med[#it{E}]",histName.Data());
  fHistEnergyMedian = new TH1F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.);
  fOutput->Add(fHistEnergyMedian);

  histName = Form("fHistEnergyMedianExLP");
  histTitle = Form("%s;med[#it{E}]",histName.Data());
  fHistEnergyMedianExLP = new TH1F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.);
  fOutput->Add(fHistEnergyMedianExLP);

  histName = Form("fHistEnergySum");
  histTitle = Form("%s;#sum[#it{E}]",histName.Data());
  fHistEnergySum = new TH1F(histName.Data(),histTitle.Data(),nE,minE,maxE);
  fOutput->Add(fHistEnergySum);

  histName = Form("fHistEnergySumExLP");
  histTitle = Form("%s;#sum[#it{E}]",histName.Data());
  fHistEnergySumExLP = new TH1F(histName.Data(),histTitle.Data(),nE,minE,maxE);
  fOutput->Add(fHistEnergySumExLP);

  histName = Form("fHistTruncatedMean");
  histTitle = Form("%s;#sum[#it{E}]",histName.Data());
  fHistTruncatedMean = new TH1F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.);
  fOutput->Add(fHistTruncatedMean);

  histName = Form("fHistTracks");
  histTitle = Form("%s;#it{N}_{tracks}",histName.Data());
  fHistTracks = new TH1F(histName.Data(),histTitle.Data(),300,0.,300.);
  fOutput->Add(fHistTracks);

  histName = Form("fHistTracklets");
  histTitle = Form("%s;#it{N}_{tracklets}",histName.Data());
  fHistTracklets = new TH1F(histName.Data(),histTitle.Data(),300,0.,300.);
  fOutput->Add(fHistTracklets);

  histName = Form("fHistV0MultSum");
  histTitle = Form("%s;mult[V0A+V0C]",histName.Data());
  fHistV0MultSum = new TH1F(histName.Data(),histTitle.Data(),500,0.,500.);
  fOutput->Add(fHistV0MultSum);

  const Int_t nMultEst = 3;
  Int_t nBinsMultEst[nMultEst] = {300,200,500};
  Double_t multEstMax[nMultEst] = {300.,200.,500.};
  TString strMultEst[nMultEst] = {"Tracks","Tracklets","V0MultSum"};
  for(Int_t i = 0; i<nMultEst; i++) {
    histName = Form("fHistEnergyMedianEst%s",strMultEst[i].Data());
    histTitle = Form("%s;med[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergyMedianEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergyMedianEst[i]);

    histName = Form("fHistEnergyMedianExLPEst%s",strMultEst[i].Data());
    histTitle = Form("%s;med[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergyMedianExLPEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergyMedianExLPEst[i]);

    histName = Form("fHistEnergySumEst%s",strMultEst[i].Data());
    histTitle = Form("%s;#sum[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergySumEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergySumEst[i]);

    histName = Form("fHistEnergySumExLPEst%s",strMultEst[i].Data());
    histTitle = Form("%s;#sum[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergySumExLPEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergySumExLPEst[i]);

    histName = Form("fHistEnergySumAvgEst%s",strMultEst[i].Data());
    histTitle = Form("%s;#sum[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergySumAvgEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergySumAvgEst[i]);

    histName = Form("fHistEnergySumAvgExLPEst%s",strMultEst[i].Data());
    histTitle = Form("%s;#sum[#it{E}];%s",histName.Data(),strMultEst[i].Data());
    fHistEnergySumAvgExLPEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistEnergySumAvgExLPEst[i]);

    histName = Form("fHistTruncatedMeanEst%s",strMultEst[i].Data());
    histTitle = Form("%s;#LT#it{E}#GT_{trunc};%s",histName.Data(),strMultEst[i].Data());
    fHistTruncatedMeanEst[i] = new TH2F(histName.Data(),histTitle.Data(),nE,minE,maxE/4.,nBinsMultEst[i],0.,multEstMax[i]);
    fOutput->Add(fHistTruncatedMeanEst[i]);
  }

  histName = Form("fHistTracksTracklets");
  histTitle = Form("%s;%s;%s",histName.Data(),strMultEst[0].Data(),strMultEst[1].Data());
  fHistTracksTracklets = new TH2F(histName.Data(),histTitle.Data(),nBinsMultEst[0],0.,multEstMax[0],nBinsMultEst[1],0.,multEstMax[1]);
  fOutput->Add(fHistTracksTracklets);

  histName = Form("fHistTracksV0MultSum");
  histTitle = Form("%s;%s;%s",histName.Data(),strMultEst[0].Data(),strMultEst[2].Data());
  fHistTracksV0MultSum = new TH2F(histName.Data(),histTitle.Data(),nBinsMultEst[0],0.,multEstMax[0],nBinsMultEst[2],0.,multEstMax[2]);
  fOutput->Add(fHistTracksV0MultSum);

  histName = Form("fHistSPDTrkClsSum");
  histTitle = Form("%s;#it{N}_{tracklets,SPD};#it{N}_{clusters,SPD};#sum[#it{E}]",histName.Data());
  fHistSPDTrkClsSum = new TH3F(histName.Data(),histTitle.Data(),nBinsMultEst[1],0.,multEstMax[1],200,0.,1000.,nE,minE,maxE);
  fOutput->Add(fHistSPDTrkClsSum);

  histName = Form("fHistSPDTrkClsSumExLP");
  histTitle = Form("%s;#it{N}_{tracklets,SPD};#it{N}_{clusters,SPD};#sum[#it{E}]",histName.Data());
  fHistSPDTrkClsSumExLP = new TH3F(histName.Data(),histTitle.Data(),nBinsMultEst[1],0.,multEstMax[1],200,0.,1000.,nE,minE,maxE);
  fOutput->Add(fHistSPDTrkClsSumExLP);

  histName = Form("fHistSPDTrkClsMedian");
  histTitle = Form("%s;#it{N}_{tracklets,SPD};#it{N}_{clusters,SPD};med[#it{E}]",histName.Data());
  fHistSPDTrkClsMedian = new TH3F(histName.Data(),histTitle.Data(),nBinsMultEst[1],0.,multEstMax[1],200,0.,1000.,nE,minE,maxE/4.);
  fOutput->Add(fHistSPDTrkClsMedian);

  histName = Form("fHistSPDTrkClsMedianExLP");
  histTitle = Form("%s;#it{N}_{tracklets,SPD};#it{N}_{clusters,SPD};med[#it{E}]",histName.Data());
  fHistSPDTrkClsMedianExLP = new TH3F(histName.Data(),histTitle.Data(),nBinsMultEst[1],0.,multEstMax[1],200,0.,1000.,nE,minE,maxE/4.);
  fOutput->Add(fHistSPDTrkClsMedianExLP);

  histName = Form("fHistSPDTrkClsTruncMean");
  histTitle = Form("%s;#it{N}_{tracklets,SPD};#it{N}_{clusters,SPD};#LT#it{E}#GT_{trunc}",histName.Data());
  fHistSPDTrkClsTruncMean = new TH3F(histName.Data(),histTitle.Data(),nBinsMultEst[1],0.,multEstMax[1],200,0.,1000.,nE,minE,maxE/4.);
  fOutput->Add(fHistSPDTrkClsTruncMean);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHighMultTrigger::FillHistograms()
{
  // Fill histograms.

  fHistEnergyMedian->Fill(fMedianEnergy);
  fHistEnergyMedianExLP->Fill(fMedianEnergyExLP);
  fHistEnergySum->Fill(fSumEnergy);
  fHistEnergySumExLP->Fill(fSumEnergyExLP);
  fHistTruncatedMean->Fill(fTruncatedMean);

  //Multiplicity estimators
  Int_t nTracks   = -1;
  if (GetParticleContainer(0)) nTracks = GetParticleContainer(0)->GetNAcceptedParticles();
  Int_t nTracklets = InputEvent()->GetMultiplicity()->GetNumberOfTracklets();

  AliVVZERO* vV0 = InputEvent()->GetVZEROData();
  Float_t multV0A=vV0->GetMTotV0A();
  Float_t multV0C=vV0->GetMTotV0C();

  fHistTracks->Fill(nTracks);
  fHistTracklets->Fill(nTracklets);
  fHistV0MultSum->Fill(multV0A+multV0C);

  const Int_t nMultEst = 3;
  Float_t multEst[nMultEst] = {(Float_t)(nTracks),(Float_t)(nTracklets),(Float_t)(multV0A+multV0C)};
  for(Int_t i = 0; i<nMultEst; i++) {
    fHistEnergyMedianEst[i]->Fill(fMedianEnergy,multEst[i]);
    fHistEnergyMedianExLPEst[i]->Fill(fMedianEnergyExLP,multEst[i]);
    fHistEnergySumEst[i]->Fill(fSumEnergy,multEst[i]);
    fHistEnergySumExLPEst[i]->Fill(fSumEnergyExLP,multEst[i]);
    if(fNAccPatches>0)
      fHistEnergySumAvgEst[i]->Fill(fSumEnergy/((Double_t)fNAccPatches),multEst[i]);
    if(fNAccPatches>1)
      fHistEnergySumAvgExLPEst[i]->Fill(fSumEnergyExLP/((Double_t)fNAccPatches-1.),multEst[i]);
    fHistTruncatedMeanEst[i]->Fill(fTruncatedMean,multEst[i]);
  }

  fHistTracksTracklets->Fill(multEst[0],multEst[1]);
  fHistTracksV0MultSum->Fill(multEst[0],multEst[2]);

  Int_t nClustersLayer0 = InputEvent()->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = InputEvent()->GetNumberOfITSClusters(1);
  fHistSPDTrkClsSum->Fill(multEst[1],(Float_t)(nClustersLayer0+nClustersLayer1),fSumEnergy);
  fHistSPDTrkClsSumExLP->Fill(multEst[1],(Float_t)(nClustersLayer0+nClustersLayer1),fSumEnergyExLP);
  fHistSPDTrkClsMedian->Fill(multEst[1],(Float_t)(nClustersLayer0+nClustersLayer1),fMedianEnergy);
  fHistSPDTrkClsMedianExLP->Fill(multEst[1],(Float_t)(nClustersLayer0+nClustersLayer1),fMedianEnergyExLP);
  fHistSPDTrkClsTruncMean->Fill(multEst[1],(Float_t)(nClustersLayer0+nClustersLayer1),fTruncatedMean);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHighMultTrigger::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHighMultTrigger::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  if(!fTriggerPatchInfo) {
    AliDebug(11,Form("Couldn't find patch info object %s",fCaloTriggerPatchInfoName.Data()));
    return kFALSE;
  }

  Int_t nPatch = fTriggerPatchInfo->GetEntriesFast();

  //sort patches
  Double_t ptarr[999] = {0};
  Int_t indexes[999] = {0};
  Int_t iacc = 0;
  for(Int_t i = 0; i<nPatch; i++) {
    AliEMCALTriggerPatchInfo *patch = dynamic_cast<AliEMCALTriggerPatchInfo*>(fTriggerPatchInfo->At(i));
    if(!patch) continue;
    if(patch->GetPatchE()>0.) {
      ptarr[iacc] = patch->GetPatchE();
      iacc++;
    }
    fHistPatchEtaPhiE->Fill(patch->GetEtaMin(),patch->GetPhiMin(),patch->GetPatchE());
  }
  
  TMath::Sort(nPatch,ptarr,indexes);
  Double_t ptarrSort[999];
  // Double_t indexesSort[999];
  Double_t ptarrSortExLP[999];
  //  Double_t indexesSortExLP[999];
  for(Int_t i = 0; i<iacc; i++) {
    ptarrSort[i] = ptarr[indexes[i]];
    //  indexesSort[i] = indexes[i];
    if(i>=fNExLP) {
      ptarrSortExLP[i] = ptarr[indexes[i]];
      //  indexesSortExLP[i] = indexes[i];
    }
  }

  //Calculate sum energy 
  fSumEnergy = 0.;
  fSumEnergyExLP = 0.;
  for(Int_t i = 0; i<iacc; i++) {
    fSumEnergy+= ptarr[indexes[i]];
    if(i>=fNExLP)
      fSumEnergyExLP+= ptarr[indexes[i]];
  }

  //calculate median of all and excluding leading patch(es)
  fMedianEnergy = TMath::Median(iacc,ptarrSort);
  fMedianEnergyExLP = TMath::Median(iacc-fNExLP,ptarrSortExLP);
  fNAccPatches = iacc;

  //calculate truncated mean
  Double_t it = 0.;
  Double_t sum = 0.;
  for(Int_t i = 0; i<iacc; i++) {
    if(ptarr[indexes[i]]<fTruncateThreshold) {
      sum+= ptarr[indexes[i]];
      it+=1.;
    }
  }
  if(it>0)
    fTruncatedMean = sum/it;

   return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHighMultTrigger::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

