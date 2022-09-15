// Implementation of the Intermittency analysis
// for charged particles, two dimensions (eta and phi)
//Contact: ramni.gupta@cern.ch
//Contributors: R.Gupta, S.Sharma, S.K.Malik 

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "AliAODInputHandler.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskNFMs.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisDataContainer.h"
#include "TFile.h"
class AliAnalysisTaskNFMs;
using namespace std;

ClassImp(AliAnalysisTaskNFMs)

AliAnalysisTaskNFMs::AliAnalysisTaskNFMs() 
  : AliAnalysisTaskSE(),  
  fAOD(0),
  fOutHList(0),
  fQAList(0),
  fEventCuts(0),
  fNtupleListBin1(0),
  fNtupleListBin2(0),
  fNtupleListBin3(0),
  fNtupleListBin4(0),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  minEta(-0.8),
  maxEta(0.8),
  fHistQAEta(0), 
  fHistQAPhi(0),
  fHistQAVx(0),
  fHistQAVy(0),
  fHistQAVz(0), 
  fEventCounter(0),
  fHistQACent(0),
  NoOfBins(0)
{

}
AliAnalysisTaskNFMs::AliAnalysisTaskNFMs(const char *name)
  : AliAnalysisTaskSE(name),
  fAOD(0),
  fOutHList(0),
  fQAList(0),
  fEventCuts(0),
  fNtupleListBin1(0),
  fNtupleListBin2(0),
  fNtupleListBin3(0),
  fNtupleListBin4(0),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  minEta(-0.8),
  maxEta(0.8),
  fHistQAEta(0), 
  fHistQAPhi(0),
  fHistQAVx(0),
  fHistQAVy(0),
  fHistQAVz(0), 
  fEventCounter(0),
  fHistQACent(0),
  NoOfBins(0)
{
  // Constructor
  Info("AliAnalysisTaskNFMs","Specific Constructor");
  for(Int_t j = 0; j < M; j++) 
  {
    fHEtaPhiBin1[j] = NULL;
    fHEtaPhiBin2[j] = NULL;
    fHEtaPhiBin3[j] = NULL;
    fHEtaPhiBin4[j] = NULL;

    fhMapEtaPhiBin1M[j] = NULL;
    fhMapEtaPhiBin2M[j] = NULL;
    fhMapEtaPhiBin3M[j] = NULL;
    fhMapEtaPhiBin4M[j] = NULL;

    fntpMBin1[j] = NULL;
    fntpMBin2[j] = NULL; 
    fntpMBin3[j] = NULL;
    fntpMBin4[j] = NULL;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
  DefineOutput(6, TList::Class());

}

AliAnalysisTaskNFMs::~AliAnalysisTaskNFMs()
{
  // destructor
  if (fOutHList)
  {
    delete fOutHList;
  }
}

void AliAnalysisTaskNFMs::UserCreateOutputObjects()
{
  fOutHList = new TList();
  fNtupleListBin1 = new TList();
  fNtupleListBin2 = new TList();
  fNtupleListBin3 = new TList();
  fNtupleListBin4 = new TList();

  fOutHList->SetOwner(kTRUE);
  fNtupleListBin1->SetOwner(kTRUE);
  fNtupleListBin2->SetOwner(kTRUE);
  fNtupleListBin3->SetOwner(kTRUE);
  fNtupleListBin4->SetOwner(kTRUE);
    
  fQAList = new TList();
  fQAList->SetOwner(kTRUE);
  fEventCuts.AddQAplotsToList(fQAList,kTRUE);
  PostData(6,fQAList);

  //Value of M init.:
  Int_t Mbins[M];
  
  if(!(fIsMmaxET))
  {
    for(Int_t ind = 0; ind < M; ind++) Mbins[ind] = 3 * (ind +2);
  }
  
  if(fIsMmaxET)
  {
    for(Int_t ind = 0; ind < M; ind++) Mbins[ind] = 2 * (ind +2);
  }

  //Eventcounter hist defined
  fEventCounter = new TH1D("fEventCounter","histo to keep track",10,0.5,10.5);
  fEventCounter->GetXaxis()->SetBinLabel(1,"Events before cuts");
  fEventCounter->GetXaxis()->SetBinLabel(2,"Events Accessed");
  fEventCounter->GetXaxis()->SetBinLabel(3,"Events without proper vertex");
  fEventCounter->GetXaxis()->SetBinLabel(4,"Events with a proper vertex");
  fEventCounter->GetXaxis()->SetBinLabel(5,"Events with 0-5% Centrality");
  fEventCounter->GetXaxis()->SetBinLabel(6,"Events Analyzed");
  fOutHList->Add(fEventCounter); 

  //Hists defined 
  fHistQACent = new TH1F("fHistQACent","Centrality Distribution", 10, minCent-2, maxCent+2);
  fOutHList->Add(fHistQACent);
  fHistQAEta = new TH1F("fHistQAEta", "#eta distribution", 1000, minEta, maxEta);
  fHistQAPhi = new TH1F("fHistQAPhi", "#phi distribution", 1000, 0.0, 6.5);
  fHistQAVx = new TH1F("fHistQAVx","Primary vertex distribution -x coordinate;V_{x}",1500,-0.3,0.3);
  fHistQAVy = new TH1F("fHistQAVy","Primary vertex distribution -y coordinate;V_{y}",1500,-0.3,0.3);
  fHistQAVz = new TH1F("fHistQAVz","Primary vertex distribution -z coordinate;V_{z}",1500,-20.0,20.0);
  fOutHList->Add(fHistQAEta);
  fOutHList->Add(fHistQAPhi);
  fOutHList->Add(fHistQAVx);
  fOutHList->Add(fHistQAVy);
  fOutHList->Add(fHistQAVz);

  for (int iSub = 1; iSub < 5; iSub++)
  {
    fHistPtBin[iSub-1] = new TH1F(Form("fHistPtBin%d", iSub), "p_{T}dist", 1000, 0.1, 5.5);
    fEtaBin[iSub-1] = new TH1F(Form("fHistEtaBin%d", iSub), "#eta dist",1000, minEta, maxEta);
    fPhiBin[iSub-1] = new TH1F(Form("fHistPhiBin%d", iSub), "#phi dist", 1000, 0., 6.5);
    fHistMulBin[iSub-1] = new TH1F(Form("fHistMultBin%d", iSub), "mult dist", 2000, 0.0, 9200.0);
    fOutHList->Add(fHistPtBin[iSub-1]);
    fOutHList->Add(fHistMulBin[iSub-1]);
    fOutHList->Add(fEtaBin[iSub-1]);
    fOutHList->Add(fPhiBin[iSub-1]);
  }


  for (int p = 1; p < M+1; p++) 
  {
    fHEtaPhiBin1[p-1] = new TH2D(Form("fHEtaPhipT1M%d", (p)), Form("PtCut1 (#eta) and (#phi) M%d Binning",(p)), Mbins[p-1], minEta, maxEta, Mbins[p-1], 0.0, 6.30);
    fHEtaPhiBin2[p-1] = new TH2D(Form("fHEtaPhipT2M%d", (p)), Form("PtCut2 (#eta) and (#phi) M%d Binning",(p)), Mbins[p-1], minEta, maxEta, Mbins[p-1], 0.0, 6.30);
    fHEtaPhiBin3[p-1] = new TH2D(Form("fHEtaPhipT3M%d", (p)), Form("PtCut3 (#eta) and (#phi) M%d Binning",(p)), Mbins[p-1], minEta, maxEta, Mbins[p-1], 0.0, 6.30);
    fHEtaPhiBin4[p-1] = new TH2D(Form("fHEtaPhipT4M%d", (p)), Form("PtCut4 (#eta) and (#phi) M%d Binning",(p)), Mbins[p-1], minEta, maxEta, Mbins[p-1], 0.0, 6.30);
    fOutHList->Add(fHEtaPhiBin1[p-1]);
    fOutHList->Add(fHEtaPhiBin2[p-1]);
    fOutHList->Add(fHEtaPhiBin3[p-1]);
    fOutHList->Add(fHEtaPhiBin4[p-1]);
  } 
  for (int b = 1; b < M+1; b++) 
  {
    fntpMBin1[b-1] = new TNtuple(Form("fntpMBin1%d", (b)), Form("Fmsforetaphicut_%d",(b)),"Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    fntpMBin2[b-1] = new TNtuple(Form("fntpMBin2%d", (b)), Form("Fmsforetaphicut_%d",(b)),"Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    fntpMBin3[b-1] = new TNtuple(Form("fntpMBin3%d", (b)), Form("Fmsforetaphicut_%d",(b)),"Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    fntpMBin4[b-1] = new TNtuple(Form("fntpMBin4%d", (b)), Form("Fmsforetaphicut_%d",(b)),"Mbins:Av_bincontent:Fq2e:Fq3e:Fq4e:Fq5e:Fq6e:Fq7e");
    fNtupleListBin1->Add(fntpMBin1[b-1]);
    fNtupleListBin2->Add(fntpMBin2[b-1]);
    fNtupleListBin3->Add(fntpMBin3[b-1]);
    fNtupleListBin4->Add(fntpMBin4[b-1]);
  }

  PostData(1, fOutHList);
  PostData(2, fNtupleListBin1);
  PostData(3, fNtupleListBin2);
  PostData(4, fNtupleListBin3);
  PostData(5, fNtupleListBin4);

}

//Execution - called per event:
void AliAnalysisTaskNFMs::UserExec(Option_t *) 
{
  fEventCounter->Fill(1);  

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) 
  {
    AliWarning("ERROR: AliAODEvent not available \n");
    return;
  }
  // fEventCuts.SetupLHC15o();
  if(fIsPileUpCuts) {fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);}
  
  Bool_t dummy = fEventCuts.AcceptEvent(fAOD); //for QA
  if(!dummy)return ;
  
  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  /* Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; if(!(isINT7selected)) return; */
  if(fisINT7) {Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7; if(!(isINT7selected)) return;}
  if(fiskMB) {Bool_t iskMBselected = fSelectMask& AliVEvent::kMB; if(!(iskMBselected))return ;}

  fEventCounter->Fill(2);

  AliAODVertex *vertex = fAOD->GetPrimaryVertex();
  if (vertex->GetNContributors() < 1)
  {   
    fEventCounter->Fill(3);
    return;
  }
  Float_t lvx = (Float_t) vertex->GetX();
  Float_t lvy = (Float_t) vertex->GetY();
  Float_t lvz = (Float_t) vertex->GetZ();

  if(fabs(lvx) < fVxMax) 
  {
    if(fabs(lvy) < fVyMax) 
    {
      if(fabs(lvz) < fVzMax) 
      {
        fHistQAVx->Fill(lvx);
        fHistQAVy->Fill(lvy);
        fHistQAVz->Fill(lvz);
      }
    }
  }
  else return;
  fEventCounter->Fill(4);

  if(fRun1data)
  {
    //Centrality slecetion for 2.76 TeV
    AliCentrality* centrality = fAOD->GetCentrality();
    Double_t fCentrality = centrality->GetCentralityPercentile("V0M");
    if(fCentrality < minCent || fCentrality >= maxCent) return;
    fHistQACent->Fill(fCentrality);
  }
  else
  {
    //Centrality slecetion for 5.02TeV
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    double lMultiPercentile = -1;
    lMultiPercentile = MultSelection->GetMultiplicityPercentile("V0M"); 
    if(lMultiPercentile < minCent || lMultiPercentile >= maxCent) return;
    fHistQACent->Fill(lMultiPercentile);
  }

  fEventCounter->Fill(5); 

  //counter++;

  //Reset Histos:
  for(Int_t i = 0; i < M; i++)
  {
    if(fHEtaPhiBin1[i]) fHEtaPhiBin1[i]->Reset();
    if(fHEtaPhiBin2[i]) fHEtaPhiBin2[i]->Reset();
    if(fHEtaPhiBin3[i]) fHEtaPhiBin3[i]->Reset();
    if(fHEtaPhiBin4[i]) fHEtaPhiBin4[i]->Reset();
  }
  fEventCounter->Fill(6);

  Int_t counterEtacutBin1 = 0, counterEtacutBin2 = 0, counterEtacutBin3 = 0, counterEtacutBin4 = 0;
  Int_t nTracks(fAOD->GetNumberOfTracks());
  Int_t status_write = 0;

  //Track Loop:
  for (Int_t i(0); i < nTracks; i++) 
    {
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!(track->TestFilterBit(filterBit)))
        continue;

      if(fClflag){if(!(track->GetTPCCrossedRows() < fNcls)) continue;}
      if(fCrflag){if(!(track->GetTPCNcls() < fNcrows)) continue;}

      Float_t Pt  = track->Pt();         
      Float_t Eta = track->Eta();
      Float_t Phi = track->Phi();
      if (fabs(track->Eta()) >= 0.8)
        continue;
      
      fHistQAEta->Fill(Eta);
      fHistQAPhi->Fill(Phi);

      if(( Pt >= ptbin[0] && Pt <= ptbin[1]))
      {
        counterEtacutBin1++;
        fHistPtBin[0]->Fill(Pt);  
        fEtaBin[0]->Fill(Eta);   
        fPhiBin[0]->Fill(Phi);    
        for(Int_t k =0; k < M; k++) fHEtaPhiBin1[k]->Fill(Eta,Phi);   
      } 
      if(( Pt >= ptbin[2] && Pt <= ptbin[3]))
      {
        counterEtacutBin2++;
        fHistPtBin[1]->Fill(Pt);  
        fEtaBin[1]->Fill(Eta);    
        fPhiBin[1]->Fill(Phi);    
        for(Int_t k =0; k < M; k++) fHEtaPhiBin2[k]->Fill(Eta,Phi);   
      }                          
      if(( Pt >= ptbin[4] && Pt <= ptbin[5]))
      {
        counterEtacutBin3++;
        fHistPtBin[2]->Fill(Pt);  
        fEtaBin[2]->Fill(Eta);    
        fPhiBin[2]->Fill(Phi);    
        for(Int_t k =0; k < M; k++) fHEtaPhiBin3[k]->Fill(Eta,Phi);           
      } 
      if(( Pt >= ptbin[6] && Pt <= ptbin[7]))
      {
        counterEtacutBin4++;
        fHistPtBin[3]->Fill(Pt);  
        fEtaBin[3]->Fill(Eta);    
        fPhiBin[3]->Fill(Phi);    
        for(Int_t k =0; k < M; k++) fHEtaPhiBin4[k]->Fill(Eta,Phi);                
      } 
    }
  fHistMulBin[0]->Fill(counterEtacutBin1);
  fHistMulBin[1]->Fill(counterEtacutBin2);
  fHistMulBin[2]->Fill(counterEtacutBin3);
  fHistMulBin[3]->Fill(counterEtacutBin4);

  //Calculation of the Normalized Fq Moments:
  for(Int_t fNoOfpTbins =0 ;fNoOfpTbins < 4; fNoOfpTbins++) 
  {
    for(Int_t binset = 0; binset < M; binset++)
    {   
      if(status_write == 0)
      { 
        fEventCounter->Fill(7); status_write++; 
      }

      if(fIsMmaxET)  NoOfBins = 2 * (binset+2);              
      if(!(fIsMmaxET)) NoOfBins = 3 * (binset+2);              
      Double_t MSquare = TMath::Power(NoOfBins,D);    
      Double_t SumOfbincontent = 0;                   
      Double_t FqEvent[Q];
      Double_t sumoff[Q];
      Double_t bincontent;  
      Double_t Mbin = NoOfBins;
      Int_t NofXetabins = 0;
      Int_t NofXphibins = 0;

      for(Int_t index = 0; index < Q; index++)
      { 
        FqEvent[index] = 0.0;
        sumoff[index] = 0.0;  
      }

      if(fNoOfpTbins ==0)
      {
        NofXetabins = fHEtaPhiBin1[binset]->GetNbinsX();
        NofXphibins = fHEtaPhiBin1[binset]->GetNbinsY(); 
      }
      if(fNoOfpTbins ==1)
      {
        NofXetabins = fHEtaPhiBin2[binset]->GetNbinsX();
        NofXphibins = fHEtaPhiBin2[binset]->GetNbinsY(); 
      }
      if(fNoOfpTbins ==2)
      {
        NofXetabins = fHEtaPhiBin2[binset]->GetNbinsX();
        NofXphibins = fHEtaPhiBin2[binset]->GetNbinsY(); 
      }
      if(fNoOfpTbins ==3)
      {
        NofXetabins = fHEtaPhiBin4[binset]->GetNbinsX();
        NofXphibins = fHEtaPhiBin4[binset]->GetNbinsY(); 
      }

      for(Int_t etabin = 1; etabin <= NofXetabins; etabin++) 
      {
        for(Int_t phibin = 1; phibin <= NofXphibins; phibin++) 
        { 
          bincontent = 0.0; 
          if(fNoOfpTbins ==0) bincontent = fHEtaPhiBin1[binset]->GetBinContent(etabin,phibin);
          if(fNoOfpTbins ==1) bincontent = fHEtaPhiBin2[binset]->GetBinContent(etabin,phibin);
          if(fNoOfpTbins ==2) bincontent = fHEtaPhiBin3[binset]->GetBinContent(etabin,phibin);
          if(fNoOfpTbins ==3) bincontent = fHEtaPhiBin4[binset]->GetBinContent(etabin,phibin);
          SumOfbincontent += bincontent;                          		

          for(Int_t q = 0; q < Q; q++) 
          { 
            if(bincontent  >= (q+2)) 
            {
              Double_t Fqeofbin = 0.0;
              Fqeofbin = TMath::Factorial(bincontent) / TMath::Factorial(bincontent-(q+2));            
              sumoff[q] += Fqeofbin;                                  
            }
          }
        }
      }

      Double_t Av_bincontent = SumOfbincontent / MSquare; 
      for(Int_t q = 0; q < Q; q++) {  
        if(sumoff[q] > 0.0) FqEvent[q] = sumoff[q] / (MSquare);               
      }

      Float_t Fq2e = FqEvent[0];
      Float_t Fq3e = FqEvent[1];
      Float_t Fq4e = FqEvent[2];
      Float_t Fq5e = FqEvent[3];
      Float_t Fq6e = FqEvent[4];
      Float_t Fq7e = FqEvent[5];

      if(fNoOfpTbins == 0)  fntpMBin1[binset]->Fill(Mbin, Av_bincontent, Fq2e , Fq3e, Fq4e, Fq5e, Fq6e, Fq7e); 
      if(fNoOfpTbins == 1)  fntpMBin2[binset]->Fill(Mbin, Av_bincontent, Fq2e , Fq3e, Fq4e, Fq5e, Fq6e, Fq7e);
      if(fNoOfpTbins == 2)  fntpMBin3[binset]->Fill(Mbin, Av_bincontent, Fq2e , Fq3e, Fq4e, Fq5e, Fq6e, Fq7e);
      if(fNoOfpTbins == 3)  fntpMBin4[binset]->Fill(Mbin, Av_bincontent, Fq2e , Fq3e, Fq4e, Fq5e, Fq6e, Fq7e);
    } 
  }//end of calculation loop

  PostData(1, fOutHList);
  PostData(2, fNtupleListBin1);
  PostData(3, fNtupleListBin2);
  PostData(4, fNtupleListBin3);
  PostData(5, fNtupleListBin4); 
}
//________________________________________________________________________

void AliAnalysisTaskNFMs::Terminate(Option_t *)
{
  Info("AliAnalysisTaskNFMs","Task Successfully finished");
  // terminate
  // called at the END of the analysis (when all events are processed)
}
