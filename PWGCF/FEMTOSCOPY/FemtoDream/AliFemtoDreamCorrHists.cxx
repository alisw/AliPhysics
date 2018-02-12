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
,fSameEventDist(0)
,fPairCounterSE(0)
,fMixedEventDist(0)
,fPairCounterME(0)
{
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf) {

  fResults=new TList();
  fResults->SetName("Results");
  fResults->SetOwner();

  fQA=new TList();
  fQA->SetName("PairQA");
  fQA->SetOwner();

  int nParticles=conf->GetNParticles();
  const int nHists=conf->GetNParticleCombinations();
  std::vector<int> NBinsHist=conf->GetNBinsHist();
  std::vector<int>::iterator itNBins=NBinsHist.begin();
  std::vector<double> kRelMin=conf->GetMinKRel();
  std::vector<double>::iterator itKMin=kRelMin.begin();
  std::vector<double> kRelMax=conf->GetMaxKRel();
  std::vector<double>::iterator itKMax=kRelMax.begin();

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
  fPairs=new TList*[nHists];
  fPairQA=new TList*[nHists];
  fSameEventDist=new TH1F*[nHists];
  fPairCounterSE=new TH2F*[nHists];
  fMixedEventDist=new TH1F*[nHists];
  fPairCounterME=new TH2F*[nHists];

  int Counter=0;
  for (int iPar1 = 0; iPar1 < nParticles; ++iPar1) {
    for (int iPar2 = iPar1; iPar2 < nParticles; ++iPar2) {
      fPairs[Counter]=new TList();
      TString PairFolderName=Form("Particle%d_Particle%d",iPar1,iPar2);
      fPairs[Counter]->SetName(PairFolderName.Data());
      fPairs[Counter]->SetOwner();
      fResults->Add(fPairs[Counter]);

      fPairQA[Counter]=new TList();
      TString PairQAName=Form("QA_Particle%d_Particle%d",iPar1,iPar2);
      fPairQA[Counter]->SetName(PairQAName.Data());
      fPairQA[Counter]->SetOwner();
      fQA->Add(fPairQA[Counter]);

      TString SameEventName=Form("SEDist_Particle%d_Particle%d",iPar1,iPar2);
      fSameEventDist[Counter]=new TH1F(SameEventName.Data(),
                                       SameEventName.Data(),
                                       *itNBins,*itKMin,*itKMax);
      fSameEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fSameEventDist[Counter]);

      TString PairCounterSEName=
          Form("SEPairs_Particle%d_Particle%d",iPar1,iPar2);
      fPairCounterSE[Counter]=new TH2F(
          PairCounterSEName.Data(),PairCounterSEName.Data(),20,0,20,20,0,20);
      fPairCounterSE[Counter]->Sumw2();
      fPairCounterSE[Counter]->GetXaxis()->SetTitle(Form("Particle%d",iPar1));
      fPairCounterSE[Counter]->GetYaxis()->SetTitle(Form("Particle%d",iPar2));
      fPairQA[Counter]->Add(fPairCounterSE[Counter]);

      TString MixedEventName=Form("MEDist_Particle%d_Particle%d",iPar1,iPar2);
      fMixedEventDist[Counter]=new TH1F(MixedEventName.Data(),
                                        MixedEventName.Data(),
                                       *itNBins,*itKMin,*itKMax);
      fMixedEventDist[Counter]->Sumw2();
      fPairs[Counter]->Add(fMixedEventDist[Counter]);

      TString PairCounterMEName=
          Form("MEPairs_Particle%d_Particle%d",iPar1,iPar2);
      fPairCounterME[Counter]=new TH2F(
          PairCounterMEName.Data(),PairCounterMEName.Data(),20,0,20,20,0,20);
      fPairCounterME[Counter]->Sumw2();
      fPairCounterME[Counter]->GetXaxis()->SetTitle(Form("Particle%d",iPar1));
      fPairCounterME[Counter]->GetYaxis()->SetTitle(Form("Particle%d",iPar2));
      fPairQA[Counter]->Add(fPairCounterME[Counter]);

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

