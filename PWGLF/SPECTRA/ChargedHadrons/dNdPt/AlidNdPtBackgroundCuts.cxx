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

#include <iostream>
#include <TList.h>
#include <TAxis.h>
#include <THnSparse.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "TParticle.h"

#include "AlidNdPtBackgroundCuts.h"

using namespace std;

ClassImp(AlidNdPtBackgroundCuts)

//_____________________________________________________________________________
AlidNdPtBackgroundCuts::AlidNdPtBackgroundCuts(const Char_t* name,const Char_t *title) : 
AliAnalysisCuts(name, title)
, fMinEta(0)
, fMaxEta(0)
, fMinPhi(0)
, fMaxPhi(0)
, fMinPt(0)
, fMaxPt(0)
, fMaxFracSharedClust(0)
, fFillControlHisto(kFALSE)
, fControlHisto(0)
{
  // default constructor 
  
  // init data members with defaults
  Init();
}

//_____________________________________________________________________________
AlidNdPtBackgroundCuts::~AlidNdPtBackgroundCuts()  
{
  // destructor
  if(fControlHisto) delete fControlHisto;
}

//_____________________________________________________________________________
void AlidNdPtBackgroundCuts::Init()  
{
  // set default values
  SetEtaWindow();
  SetPhiWindow();
  SetPtWindow();
  SetMaxFracSharedClust();

  const Int_t ptNbins = 56;
  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};

  //etasum:dphi:dpt:eta1:eta2:pt1:fracSharedClust1:qsum
  Int_t binsControlHisto[8]=  { 201,  401,              101,  30,   30,  ptNbins, 101, 3 };
  Double_t minControlHisto[8]={-3.0, -2.*TMath::Pi(),  -10,  -1.5, -1.5, 0.,      0.,  0.}; 
  Double_t maxControlHisto[8]={ 3.0,  2.*TMath::Pi(),   10,   1.5,  1.5, 16.,     1.,  3.}; 
  
  fControlHisto = new THnSparseF("fControlHisto","etasum:dphi:dpt:eta1:eta2:pt1:fracSharedClust1:qsum",8,binsControlHisto,minControlHisto,maxControlHisto);
  fControlHisto->SetBinEdges(5,binsPt);
  fControlHisto->GetAxis(0)->SetTitle("#eta1+#eta2");
  fControlHisto->GetAxis(1)->SetTitle("#phi1-#phi2 (rad)");
  fControlHisto->GetAxis(2)->SetTitle("pt1-pt2 (GeV/c)");
  fControlHisto->GetAxis(3)->SetTitle("#eta1");
  fControlHisto->GetAxis(4)->SetTitle("#eta2");
  fControlHisto->GetAxis(5)->SetTitle("pt1 (GeV/c)");
  fControlHisto->GetAxis(6)->SetTitle("fracSharedClust1");
  fControlHisto->GetAxis(7)->SetTitle("q1+q2");
  fControlHisto->Sumw2();
}

//_____________________________________________________________________________
Bool_t AlidNdPtBackgroundCuts::IsBackgroundTrack(AliESDtrack *track1, AliESDtrack *track2)
{
  // 
  // check whether track is cosmic or splitted one
  //
  if(!track1) return kFALSE;
  if(!track2) return kFALSE;
  const AliExternalTrackParam *tpcTrack1 =  track1->GetTPCInnerParam();
  const AliExternalTrackParam *tpcTrack2 =  track2->GetTPCInnerParam();
  if(!tpcTrack1) return kFALSE;
  if(!tpcTrack2) return kFALSE;

  if( IsHistogramsOn() ) 
  {
    Float_t etasum = tpcTrack1->Eta() + tpcTrack2->Eta();
    Float_t dphi   = tpcTrack1->Phi() - tpcTrack2->Phi();
    Float_t dpt    = tpcTrack1->Pt()  - tpcTrack2->Pt();
    Float_t pt1    = tpcTrack1->Pt();
    Float_t qsum   = track1->Charge() + track2->Charge();
    if(qsum == -2) qsum = 1;

    Float_t nclust1 =  track1->GetTPCNclsIter1() ; // first tracking pass
    Float_t nclust2 =  track2->GetTPCNclsIter1() ; // first tracking pass
    Float_t fracSharedClust1 = 0.0;
    if(nclust1) fracSharedClust1 = track1->GetTPCnclsS()/Float_t(nclust1); 
  
    //Float_t dsphi = (tpcTrack1->GetSnp()-tpcTrack2->GetSnp()) / TMath::Sqrt(tpcTrack1->GetSigmaSnp2()+tpcTrack2->GetSigmaSnp2());    
    //Float_t dtanl = (tpcTrack1->GetTgl()-tpcTrack2->GetTgl()) / TMath::Sqrt(tpcTrack1->GetSigmaTgl2()+tpcTrack2->GetSigmaTgl2());
    //Float_t dsphi = tpcTrack1->GetSnp()-tpcTrack2->GetSnp();
    //Float_t dtanl = tpcTrack1->GetTgl()-tpcTrack2->GetTgl();
    //


    /*
    printf("tpcTrack1->GetSnp() %e, track1->Pt() %f, track1->Theta() %f, track1->Eta() %f, track1->Phi() %f, track1->Charge() %d  \n", tpcTrack1->GetSnp(), track1->Pt(), track1->Theta(), track1->Eta(), track1->Phi(), track1->Charge());

    printf("tpcTrack2->GetSnp() %e, track2->Pt() %f, track2->Theta() %f, track2->Eta() %f, track2->Phi() %f, track2->Charge() %d  \n", tpcTrack2->GetSnp(), track2->Pt(), track2->Theta(), track2->Eta(), track2->Phi(), track2->Charge());
    */

    Double_t vControlHisto[8] = {etasum, dphi, dpt, tpcTrack1->Eta(), tpcTrack2->Eta(), pt1, fracSharedClust1,qsum};
    if(nclust1 > 70 && nclust2 > 70)
       fControlHisto->Fill(vControlHisto);
  }

  if ( IsCosmicTrack(track1,track2) || IsSplittedTrack(track1,track2) ) return kTRUE;
  else return kFALSE;

return kFALSE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtBackgroundCuts::IsCosmicTrack(AliESDtrack *track1, AliESDtrack *track2)
{
  // 
  // check whether track is cosmic
  //
  if(!track1) return kFALSE;
  if(!track2) return kFALSE;
  const AliExternalTrackParam *tpcTrack1 =  track1->GetTPCInnerParam();
  const AliExternalTrackParam *tpcTrack2 =  track2->GetTPCInnerParam();
  if(!tpcTrack1) return kFALSE;
  if(!tpcTrack2) return kFALSE;

 /*
  Float_t etasum = tpcTrack1->Eta() + tpcTrack2->Eta();
  Float_t dphi = tpcTrack1->Phi() - tpcTrack2->Phi();
  Float_t dpt  = tpcTrack1->Pt()  - tpcTrack2->Pt();
  Float_t pt1  = tpcTrack1->Pt();
  */
  Float_t qsum = track1->Charge() + track2->Charge();

  Float_t nclust =  track1->GetTPCNclsIter1() ; // first tracking pass
  Float_t fracSharedClust = 0.0;
  if(nclust) fracSharedClust = track1->GetTPCnclsS()/Float_t(nclust); 

  if( qsum != 0) return kFALSE;

return kFALSE;
}

//_____________________________________________________________________________
Bool_t AlidNdPtBackgroundCuts::IsSplittedTrack(AliESDtrack *track1, AliESDtrack *track2)
{
  // 
  // check whether track is splitted
  //
  if(!track1) return kFALSE;
  if(!track2) return kFALSE;
  const AliExternalTrackParam *tpcTrack1 =  track1->GetTPCInnerParam();
  const AliExternalTrackParam *tpcTrack2 =  track2->GetTPCInnerParam();
  if(!tpcTrack1) return kFALSE;
  if(!tpcTrack2) return kFALSE;

  /*
  Float_t etasum = tpcTrack1->Eta() + tpcTrack2->Eta();
  Float_t dphi = tpcTrack1->Phi() - tpcTrack2->Phi();
  Float_t dpt  = tpcTrack1->Pt()  - tpcTrack2->Pt();
  Float_t pt1  = tpcTrack1->Pt();
  Float_t qsum = track1->Charge() + track2->Charge();

  Float_t nclust =  track1->GetTPCNclsIter1() ; // first tracking pass
  Float_t fracSharedClust = 0.0;
  if(nclust) fracSharedClust = track1->GetTPCnclsS()/Float_t(nclust); 
  */

return kFALSE;
}



//_____________________________________________________________________________
Long64_t AlidNdPtBackgroundCuts::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)
  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AlidNdPtBackgroundCuts* entry = dynamic_cast<AlidNdPtBackgroundCuts*>(obj);
    if (entry == 0)  
      continue; 
  
    fControlHisto->Add(entry->fControlHisto);

  count++;
  }

return count;
}
