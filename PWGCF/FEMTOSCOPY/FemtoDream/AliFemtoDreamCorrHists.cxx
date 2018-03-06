/*
 * AliFemtoDreamCorrHists.cxx
 *
 *  Created on: Sep 12, 2017
 *      Author: gu74req
 */
#include <iostream>
#include <vector>
#include "AliFemtoDreamCorrHists.h"
#include "TMath.h"
ClassImp(AliFemtoDreamCorrHists)
AliFemtoDreamCorrHists::AliFemtoDreamCorrHists()
:fQA(0)
,fResults(0)
,fPairs(0)
,fPairQA(0)
,fMinimalBooking(false)
,fMomentumResolution(false)
,fSameEventDist(0)
,fSameEventMultDist(0)
,fPairCounterSE(0)
,fMixedEventDist(0)
,fMixedEventMultDist(0)
,fPairCounterME(0)
,fMomResolution(0)
,fDoMultBinning(false)
{
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(
    AliFemtoDreamCollConfig *conf,bool MinimalBooking) {
  fResults=new TList();
  fResults->SetName("Results");
  fResults->SetOwner();
  if (!fMinimalBooking) {
    fQA=new TList();
    fQA->SetName("PairQA");
    fQA->SetOwner();
  }
  fMinimalBooking=MinimalBooking;
  fMomentumResolution=conf->GetDoMomResolution();
  fDoMultBinning=conf->GetDoMultBinning();
  int multbins=conf->GetNMultBins();
  int nParticles=conf->GetNParticles();
  const int nHists=conf->GetNParticleCombinations();
  std::vector<int> NBinsHist=conf->GetNBinsHist();
  std::vector<int>::iterator itNBins=NBinsHist.begin();
  std::vector<float> kRelMin=conf->GetMinKRel();
  std::vector<float>::iterator itKMin=kRelMin.begin();
  std::vector<float> kRelMax=conf->GetMaxKRel();
  std::vector<float>::iterator itKMax=kRelMax.begin();

  if (nHists!=(int)NBinsHist.size() || nHists!=(int)kRelMin.size() ||
      nHists!=(int)kRelMax.size()) {
    //Todo: Replace by AliFatal!
    std::cout<<"something went horribly wrong!"<<std::endl;
  }
  //The way the histograms are assigned later is going to be for example for
  //4 different particle species X1,X2,X3,X4:
  //    X1  X2  X3  X4
  //X1  1   2   3   4
  //X2      5   6   7
  //X3          8   9
  //X4              10<-----Number of the Histogram=Position in input vector
  //X1 corresponds the first particle array in your vector that you hand over
  //in the AliFemtoDreamPartCollection::SetEvent Method, X2 to the second and
  //so on.
  //in case we only book the most basic things we don't need this
  if (!fMinimalBooking) {
    fPairQA=new TList*[nHists];
    fPairCounterSE=new TH2F*[nHists];
    fPairCounterME=new TH2F*[nHists];
    if (fMomentumResolution) {
      fMomResolution=new TH2F*[nHists];
    }
  }
  //we always want to do this, regardless of the booking type!
  fPairs=new TList*[nHists];
  fSameEventDist=new TH1F*[nHists];
  fMixedEventDist=new TH1F*[nHists];
  if (fDoMultBinning) {
    fSameEventMultDist=new TH2F*[nHists];
    fMixedEventMultDist=new TH2F*[nHists];
  } else {
    fSameEventMultDist=0;
    fMixedEventMultDist=0;
  }
  int Counter=0;
  for (int iPar1 = 0; iPar1 < nParticles; ++iPar1) {
    for (int iPar2 = iPar1; iPar2 < nParticles; ++iPar2) {
      fPairs[Counter]=new TList();
      TString PairFolderName=Form("Particle%d_Particle%d",iPar1,iPar2);
      fPairs[Counter]->SetName(PairFolderName.Data());
      fPairs[Counter]->SetOwner();
      fResults->Add(fPairs[Counter]);

      TString SameEventName=Form("SEDist_Particle%d_Particle%d",iPar1,iPar2);
      fSameEventDist[Counter]=new TH1F(SameEventName.Data(),
                                       SameEventName.Data(),
                                       *itNBins,*itKMin,*itKMax);
      fSameEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fSameEventDist[Counter]);

      TString MixedEventName=Form("MEDist_Particle%d_Particle%d",iPar1,iPar2);
      fMixedEventDist[Counter]=new TH1F(MixedEventName.Data(),
                                        MixedEventName.Data(),
                                        *itNBins,*itKMin,*itKMax);
      fMixedEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fMixedEventDist[Counter]);

      if (fDoMultBinning) {
        TString SameMultEventName=
            Form("SEMultDist_Particle%d_Particle%d",iPar1,iPar2);
        fSameEventMultDist[Counter]=new TH2F(SameMultEventName.Data(),
                                             SameMultEventName.Data(),
                                             *itNBins,*itKMin,*itKMax,
                                             multbins,1,multbins+1);
        fSameEventMultDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fSameEventMultDist[Counter]);

        TString MixedMultEventName=
            Form("MEMultDist_Particle%d_Particle%d",iPar1,iPar2);
        fMixedEventMultDist[Counter]=new TH2F(MixedMultEventName.Data(),
                                              MixedMultEventName.Data(),
                                              *itNBins,*itKMin,*itKMax,
                                              multbins,1,multbins+1);
        fMixedEventMultDist[Counter]->Sumw2();
        fPairs[Counter]->Add(fMixedEventMultDist[Counter]);

      }
      if (!fMinimalBooking) {
        fPairQA[Counter]=new TList();
        TString PairQAName=Form("QA_Particle%d_Particle%d",iPar1,iPar2);
        fPairQA[Counter]->SetName(PairQAName.Data());
        fPairQA[Counter]->SetOwner();
        fQA->Add(fPairQA[Counter]);
        TString PairCounterSEName=
            Form("SEPairs_Particle%d_Particle%d",iPar1,iPar2);
        fPairCounterSE[Counter]=new TH2F(
            PairCounterSEName.Data(),PairCounterSEName.Data(),20,0,20,20,0,20);
        fPairCounterSE[Counter]->Sumw2();
        fPairCounterSE[Counter]->GetXaxis()->SetTitle(Form("Particle%d",iPar1));
        fPairCounterSE[Counter]->GetYaxis()->SetTitle(Form("Particle%d",iPar2));
        fPairQA[Counter]->Add(fPairCounterSE[Counter]);

        TString PairCounterMEName=
            Form("MEPairs_Particle%d_Particle%d",iPar1,iPar2);
        fPairCounterME[Counter]=new TH2F(
            PairCounterMEName.Data(),PairCounterMEName.Data(),20,0,20,20,0,20);
        fPairCounterME[Counter]->Sumw2();
        fPairCounterME[Counter]->GetXaxis()->SetTitle(Form("Particle%d",iPar1));
        fPairCounterME[Counter]->GetYaxis()->SetTitle(Form("Particle%d",iPar2));
        fPairQA[Counter]->Add(fPairCounterME[Counter]);

        if (fMomentumResolution) {
          //take twice the number of bins we use for the CF to be sure, the range is
          //hard coded. This assumed that the input is in GeV!
//          *itNBins,*itKMin,*itKMax,
          double dNBin=(*itKMax-*itKMin)/(double)(*itNBins);
          dNBin/=2;
          int nBims=(int)(1/dNBin);
          TString MomResoName=
              Form("MomentumResolution_Particle%d_Particle%d",iPar1,iPar2);
          fMomResolution[Counter]=
              new TH2F(MomResoName.Data(),MomResoName.Data(),nBims,0,1,
                       nBims,0,1);
          fMomResolution[Counter]->Sumw2();
          fMomResolution[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolution[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolution[Counter]);
        }
      }

      ++Counter;
      ++itNBins;
      ++itKMin;
      ++itKMax;
    }
  }
}

AliFemtoDreamCorrHists::~AliFemtoDreamCorrHists() {
  if (fPairs) {
    delete [] fPairs;
    delete fPairs;
  }
  if (fSameEventDist) {
    delete [] fSameEventDist;
    delete fSameEventDist;
  }
  if (fMixedEventDist) {
    delete [] fMixedEventDist;
    delete fMixedEventDist;
  }
}

