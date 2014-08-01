
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliSpectraAODTrackCuts class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliPIDResponse.h"   
#include "AliExternalTrackParam.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODTrackCuts.h"
#include <iostream>

using namespace std;

const char * AliSpectraAODTrackCuts::kBinLabel[] ={"TrkBit",
						   "TrkCuts",
						   "TrkEta",
						   "TrkDCA",
						   "TrkP",
						   "TrkPt",
						   "TrkPtTOF",
						   "TOFMatching",
						   "kTOFout",
						   "kTIME",
						   "kTOFpid",
						   "Accepted"};


ClassImp(AliSpectraAODTrackCuts)


AliSpectraAODTrackCuts::AliSpectraAODTrackCuts(const char *name) : TNamed(name, "AOD Track Cuts"), fIsSelected(0), fTrackBits(0), fMinTPCcls(0), fRequestSPDcls(0), fEtaCutMin(0), fEtaCutMax(0), fDCACut(0), fPCut(0), fPtCut(0), fYCut(0),
  fPtCutTOFMatching(0), fHistoCuts(0), fHistoNSelectedPos(0), fHistoNSelectedNeg(0), fHistoNMatchedPos(0), fHistoNMatchedNeg(0), fHistoEtaPhiHighPt(0), fTrack(0), fPIDResponse(0)
{
  // Constructor
  fHistoCuts = new TH1I("fTrkCuts", "Track Cuts", kNTrkCuts, -0.5, kNTrkCuts - 0.5);
  for(Int_t ibin=1;ibin<=kNTrkCuts;ibin++)fHistoCuts->GetXaxis()->SetBinLabel(ibin,kBinLabel[ibin-1]);
  //standard histo
  const Double_t templBins[] = {0.20,0.30,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,5.0,6.0,7.0,8.0,9.0,10.,12.,15.};
  const Int_t nbinsTempl=26;

  fHistoNSelectedPos=new TH1F("fHistoNSelectedPos","fHistoNSelectedPos",nbinsTempl,templBins);
  fHistoNSelectedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNSelectedNeg=new TH1F("fHistoNSelectedNeg","fHistoNSelectedNeg",nbinsTempl,templBins);
  fHistoNSelectedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNMatchedPos=new TH1F("fHistoNMatchedPos","fHistoNMatchedPos",nbinsTempl,templBins);
  fHistoNMatchedPos->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoNMatchedNeg=new TH1F("fHistoNMatchedNeg","fHistoNMatchedNeg",nbinsTempl,templBins);
  fHistoNMatchedNeg->GetXaxis()->SetTitle("P_{T} (GeV / c)");
  fHistoEtaPhiHighPt=new TH2F("fHistoEtaPhiHighPt","fHistoEtaPhiHighPt",200,-1,1,400,0,7);
  fHistoEtaPhiHighPt->SetXTitle("eta");
  fHistoEtaPhiHighPt->SetYTitle("phi");
  
  fEtaCutMin = -100000.0; // default value of eta cut ~ no cut
  fEtaCutMax = 100000.0; // default value of eta cut ~ no cut
  fDCACut = 100000.0; // default value of dca cut ~ no cut
  fPCut = 100000.0; // default value of p cut ~ no cut
  fPtCut = 100000.0; // default value of pt cut ~ no cut 
  fPtCutTOFMatching=0.6; //default value fot matching with TOF
  fYCut       = 100000.0; // default value of y cut ~ no cut 
  fMinTPCcls=70; // ncls in TPC
  fRequestSPDcls=kFALSE; //request a hit in the SPD
}

//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::IsSelected(AliAODTrack * track,Bool_t FillHistStat)
{
  // Returns true if Track Cuts are selected and applied
  if (!track)
    {
      printf("ERROR: Could not receive track");
      return kFALSE;
    }
  fTrack = track;
  
  if(!CheckTrackType()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkBit);
  
  if(!CheckTrackCuts()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkCuts);
  if(!CheckEtaCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkEta);
  if(!CheckDCACut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkDCA);
  if(!CheckPCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkP);
  if(!CheckPtCut()){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kTrkPt);
  if(!CheckTOFMatching(FillHistStat)){
    return kFALSE;
  }
  if(FillHistStat)fHistoCuts->Fill(kAccepted);
  return kTRUE;
}
//_________________________________________________________

Bool_t AliSpectraAODTrackCuts::CheckTrackType()
{
  // Check track Type
  if (fTrack->TestFilterBit(fTrackBits)) return kTRUE;
  return kFALSE;
}
//_________________________________________________________

Bool_t AliSpectraAODTrackCuts::CheckTrackCuts()
{
  // Check additional track Cuts
  Bool_t PassTrackCuts=kTRUE;
  if (!fTrack->HasPointOnITSLayer(0) && !fTrack->HasPointOnITSLayer(1)  && fRequestSPDcls)PassTrackCuts=kFALSE; //FIXME 1 SPD for the moment
  if (fTrack->GetTPCNcls()<fMinTPCcls)PassTrackCuts=kFALSE;
  return PassTrackCuts;
}
//________________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckEtaCut()
{
  // Check eta cut
  if (fTrack->Eta() < fEtaCutMax && fTrack->Eta() > fEtaCutMin) return kTRUE;
  return kFALSE;
}

Bool_t AliSpectraAODTrackCuts::CheckYCut(Double_t mass) 
{
  // check if the rapidity is within the set range
  Double_t y=-1000;
  if (mass > 0.) { y = fTrack->Y(mass); }//negative mass for unidentified particles
  if (TMath::Abs(y) > fYCut || y < -998.) return kFALSE;
  return kTRUE;
}
//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckDCACut()
{
  // Check DCA cut
  if (TMath::Abs(fTrack->DCA()) < fDCACut) return kTRUE; //FIXME for newest AOD fTrack->DCA() always gives -999
  return kFALSE;
}
//________________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckPCut()
{
  // Check P cut
  if (fTrack->P() < fPCut) return kTRUE;
  return kFALSE;
}
//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckPtCut()
{
  // check Pt cut
  if (fTrack->Pt() < fPtCut) return kTRUE;
  return kFALSE;
}

//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckTOFMatching(Bool_t FillHistStat)
{
  if (fTrack->Pt() < fPtCutTOFMatching) return kTRUE;
  else{
    if(FillHistStat)fHistoCuts->Fill(kTrkPtTOF);
    if(fTrack->Charge()>0)fHistoNSelectedPos->Fill(fTrack->Pt());
    else fHistoNSelectedNeg->Fill(fTrack->Pt());
    
    // Get PID response object, if needed
    if(!fPIDResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
      fPIDResponse = inputHandler->GetPIDResponse();
    }
    if(!fPIDResponse) {
      AliFatal("Cannot get pid response");
      return 0;
    }
    
    if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fTrack)==0)return kFALSE; 
    
    //check the bits of the selected particles
    UInt_t status; 
    status=fTrack->GetStatus();
    if((status&AliAODTrack::kTOFout)&&FillHistStat)fHistoCuts->Fill(kTrTOFout);
    if((status&AliAODTrack::kTIME)&&FillHistStat)fHistoCuts->Fill(kTrTIME);
    if((status&AliAODTrack::kTOFpid)&&FillHistStat)fHistoCuts->Fill(kTrTOFpid);
    
    
    if(FillHistStat)fHistoCuts->Fill(kTOFMatching);
    if(fTrack->Charge()>0)fHistoNMatchedPos->Fill(fTrack->Pt());
    else fHistoNMatchedNeg->Fill(fTrack->Pt());
    if(fTrack->Pt()>1.5){
      //fHistoEtaPhiHighPt->Fill(fTrack->GetOuterParam()->Eta(),fTrack->GetOuterParam()->Phi());
      //Printf("AliExternalTrackParam * extpar=(AliExternalTrackParam*)fTrack->GetOuterParam();");
      //AliExternalTrackParam * extpar=(AliExternalTrackParam*)fTrack->GetOuterParam();
      fHistoEtaPhiHighPt->Fill(fTrack->Eta(),fTrack->Phi());
      //Printf("fHistoEtaPhiHighPt->Fill(extpar->Eta(),extpar->Phi());");
      //fHistoEtaPhiHighPt->Fill(extpar->Eta(),extpar->Phi());
      //delete extpar;
    }
    return kTRUE;
  }
}
//_______________________________________________________
void AliSpectraAODTrackCuts::PrintCuts() const
{
  // Print cuts
  cout << "Track Cuts" << endl;
  cout << " > TrackBit\t" << fTrackBits << endl;
  cout << " > Eta cut\t" << fEtaCutMin <<","<< fEtaCutMax << endl;
  cout << " > DCA cut\t" << fDCACut << endl;
  cout << " > P cut\t" << fPCut << endl;
  cout << " > Pt cut \t" << fPtCut << endl;
  cout << " > TPC cls \t" << fMinTPCcls << endl;
}
//_______________________________________________________
void AliSpectraAODTrackCuts::SetTrackType(UInt_t bit)
{
  // Set the type of track to be used. The argument should be the bit number. The mask is produced automatically.
  fTrackBits = (0x1 << (bit - 1));
}
//_______________________________________________________

Long64_t AliSpectraAODTrackCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODTrackCuts objects with this.
  // Returns the number of merged objects (including this).

  //  AliInfo("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;//FIXME we should only 1 collection
  TList collections_histoNSelectedPos;
  TList collections_histoNSelectedNeg;
  TList collections_histoNMatchedPos;
  TList collections_histoNMatchedNeg;
  TList collections_histoEtaPhiHighPt;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraAODTrackCuts* entry = dynamic_cast<AliSpectraAODTrackCuts*> (obj);
    if (entry == 0) 
      continue;
    
    TH1I * histo = entry->GetHistoCuts();      
    collections.Add(histo);
    TH1F * histoNSelectedPos = entry->GetHistoNSelectedPos();      
    collections_histoNSelectedPos.Add(histoNSelectedPos);
    TH1F * histoNSelectedNeg = entry->GetHistoNSelectedNeg();      
    collections_histoNSelectedNeg.Add(histoNSelectedNeg);
    TH1F * histoNMatchedPos = entry->GetHistoNMatchedPos();      
    collections_histoNMatchedPos.Add(histoNMatchedPos);
    TH1F * histoNMatchedNeg = entry->GetHistoNMatchedNeg();      
    collections_histoNMatchedNeg.Add(histoNMatchedNeg);
    TH2F * histoEtaPhiHighPt = entry->GetHistoEtaPhiHighPt();      
    collections_histoEtaPhiHighPt.Add(histoEtaPhiHighPt);
    count++;
  }
  
  fHistoCuts->Merge(&collections);
  fHistoNSelectedPos->Merge(&collections_histoNSelectedPos);
  fHistoNSelectedNeg->Merge(&collections_histoNSelectedNeg);
  fHistoNMatchedPos->Merge(&collections_histoNMatchedPos);
  fHistoNMatchedNeg->Merge(&collections_histoNMatchedNeg);
  fHistoEtaPhiHighPt->Merge(&collections_histoEtaPhiHighPt);
  
  delete iter;

  return count+1;
}

