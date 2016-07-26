// $Id$
//
// Dcal dijet performance task
//
// Author: R. Reed

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskDcalDijetPerf.h"

ClassImp(AliAnalysisTaskDcalDijetPerf)

//________________________________________________________________________
AliAnalysisTaskDcalDijetPerf::AliAnalysisTaskDcalDijetPerf() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDcalDijetPerf", kTRUE),
  fHistTracksPt(0),
  fHistTracksEtaPhi(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistJet1(0),
  fHistJet1m(0),
  fHistJet1nm(0),
  fHistJet2(0),
  fHistJet1to2(0),
  fHistDiJet1(0),
  fHistDiJet1m(0),
  fJetsCont(0),
  fJetsCont2(0),
  fJetsCont3(0),
  fTracksCont(0),
  fCaloClustersCont(0)

{
  // Default constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistTracksEtaPhi   = new TH2*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistTracksEtaPhi[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }
    fHistJet1 = 0;
    fHistJet1m = 0;
    fHistJet1nm = 0;
    fHistJet2 = 0;
    fHistJet1to2 = 0;
    fHistDiJet1 = 0;
    fHistDiJet1m = 0;
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDcalDijetPerf::AliAnalysisTaskDcalDijetPerf(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistTracksPt(0),
  fHistTracksEtaPhi(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistJet1(0),
  fHistJet1m(0),
  fHistJet1nm(0),
  fHistJet2(0),
  fHistJet1to2(0),
  fHistDiJet1(0),
  fHistDiJet1m(0),
  fJetsCont(0),
  fJetsCont2(0),
  fJetsCont3(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Standard constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistTracksEtaPhi   = new TH2*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistTracksEtaPhi[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }
  fHistJet1 = 0;
  fHistJet1m = 0;
  fHistJet1nm = 0;
  fHistJet2 = 0;
  fHistJet1to2 = 0;
  fHistDiJet1 = 0;
  fHistDiJet1m = 0;

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDcalDijetPerf::~AliAnalysisTaskDcalDijetPerf()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskDcalDijetPerf::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont           = GetJetContainer(0);
  fJetsCont2           = GetJetContainer(1);
  fJetsCont3           = GetJetContainer(2);

  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  fTracksCont->SetClassName("AliVTrack");
  //fCaloClustersCont->SetClassName("AliVCluster");
  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fParticleCollArray.GetEntriesFast()>0) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPt[i]);
      histname = "fHistTracksEtaPhi_";
      histname += i;
        fHistTracksEtaPhi[i] = new TH2F(histname.Data(), histname.Data(), fNbins, -0.7, 0.7, fNbins, 0, TMath::Pi()*2);
      fHistTracksEtaPhi[i]->GetXaxis()->SetTitle("eta");
      fHistTracksEtaPhi[i]->GetYaxis()->SetTitle("phi");
      fOutput->Add(fHistTracksEtaPhi[i]);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    if (fJetCollArray.GetEntriesFast()>0) {
      histname = "fHistLeadingJetPt_";
      histname += i;
      fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistLeadingJetPt[i]);
      
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistJetsPhiEta[i]);
      
      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 30, 0, 3);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);

      histname = "fHistJetsPtLeadHad_";
      histname += i;
      fHistJetsPtLeadHad[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistJetsPtLeadHad[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtLeadHad[i]);
    
      if (!(GetJetContainer()->GetRhoName().IsNull())) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }

  Int_t nbins[] = {150,100,100,100,150};
  Double_t xmin[] = {  0,-0.7,             0, 0,   0};
  Double_t xmax[] = {150, 0.7,TMath::TwoPi(), 1, 150};
  fHistJet1 = new THnSparseF("Jets1Collection","Jets1Collection",5,nbins,xmin,xmax);
  fOutput->Add(fHistJet1);
  fHistJet1m = new THnSparseF("Jets1CollectionMatched","Jets1Collection",5,nbins,xmin,xmax);
  fOutput->Add(fHistJet1m);
  fHistJet1nm = new THnSparseF("Jets1CollectionNotMatched","Jets1Collection",5,nbins,xmin,xmax);
  fOutput->Add(fHistJet1nm);
  fHistJet2 = new THnSparseF("Jets2Collection","Jets2Collection",5,nbins,xmin,xmax);
  fOutput->Add(fHistJet2);
  
  Int_t nbins2[] = {150,100,100,100,150,100,100,100,100};
  Double_t xmin2[] = {0,-0.7,0,0,0,-0.7,0,0,0};
  Double_t xmax2[] = {150,0.7,6.28,1,150,0.7,6.28,1,0.2};
  fHistJet1to2 = new THnSparseF("Jets1to2Collection","Jets1to2Collection",9,nbins2,xmin2,xmax2);
    fOutput->Add(fHistJet1to2);
    
    Int_t nbins3[] = {150,100,100,100,150,100,100,100,100};
    Double_t xmin3[] = {0,-0.7,0,0,0,-0.7,0,0,0};
    Double_t xmax3[] = {150,0.7,6.28,1,150,0.7,6.28,1,1};
  fHistDiJet1 = new THnSparseF("fHistDiJet1","fHistDiJet1",9,nbins3,xmin3,xmax3);
    fOutput->Add(fHistDiJet1);
    
    Int_t nbins4[] = {150,100,100,100,150,100,100,100,100,150,100,100,100};
    Double_t xmin4[] = {0,-0.7,0,0,0,-0.7,0,0,0,0,-0.7,0,0}; //pt1 eta1 phi1 NEF1 pt2 eta2 phi2 NEF2 AJ pt3 eta3 phi3 R
    Double_t xmax4[] = {150,0.7,6.28,1,150,0.7,6.28,1,1,150,0.7,6.28,0.2};
    fHistDiJet1m = new THnSparseF("fHistDiJet1m","fHistDiJet1m",13,nbins4,xmin4,xmax4);
    fOutput->Add(fHistDiJet1m);
    
    
    
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDcalDijetPerf::FillHistograms()
{
  // Fill histograms.

  if (fTracksCont) {
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
    while(track) {
      fHistTracksPt[fCentBin]->Fill(track->Pt());
      fHistTracksEtaPhi[fCentBin]->Fill(track->Eta(),track->Phi());
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }
  
  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());

      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }

    int N1 = 0;
  if (fJetsCont) {
    fJetsCont->ResetCurrentID();
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(); 
    while(jet) {
        Float_t NEFpT = jet->Pt()*jet->NEF();
        N1++;
      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());
      Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);
         Double_t jetarray[] = {jet->Pt(),jet->Eta(),jet->Phi(),jet->NEF(),NEFpT};
        fHistJet1->Fill(jetarray);
      if (fHistJetsCorrPtArea[fCentBin]) {
	Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }
      jet = fJetsCont->GetNextAcceptJet();
    }
    
    jet = fJetsCont->GetLeadingJet();
    if(jet) fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
  }//Loop over the first collection of jets
   
    int N2 = 0;
    if(fJetsCont2){
        fJetsCont2->ResetCurrentID();
        AliEmcalJet *jet = fJetsCont2->GetNextAcceptJet();
        while(jet){
            Float_t NEFpT = jet->Pt()*jet->NEF();
            N2++;
            Double_t jetarray[] = {jet->Pt(),jet->Eta(),jet->Phi(),jet->NEF(),NEFpT};
            fHistJet2->Fill(jetarray);
            jet = fJetsCont2->GetNextAcceptJet();
        }
    } // loop over the trigger jerts.
    int N1N2 = 0;
    int N1N2m = 0;
    if (fJetsCont&&fJetsCont2) {
        fJetsCont->ResetCurrentID();
        fJetsCont2->ResetCurrentID();
        AliEmcalJet *jet1 = fJetsCont->GetNextAcceptJet();
        AliEmcalJet *jet2 = fJetsCont2->GetNextAcceptJet();
        while(jet1){
            bool ismatched = false;
            Float_t NEFpT1 = jet1->Pt()*jet1->NEF();
            Double_t jetarray1[] = {jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->NEF(),NEFpT1};
            while(jet2){
                N1N2++;
                Double_t deta = jet1->Eta()-jet2->Eta();
                Double_t dphi = RelativePhi(jet1->Phi(),jet2->Phi());
                Double_t deta2 = deta*deta;
                Double_t dphi2 = dphi*dphi;
                Double_t dR = pow(deta2+dphi2,0.5);
                Double_t jetarray[] = {jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->NEF(),jet2->Pt(),jet2->Eta(),jet2->Phi(),jet2->NEF(),dR};
                if (dR<0.2){
                    N1N2m++;
                    fHistJet1to2->Fill(jetarray);
                    ismatched = true;
                }
                jet2 = fJetsCont2->GetNextAcceptJet();
            }
            if (ismatched)
                fHistJet1m->Fill(jetarray1);
            else
                fHistJet1nm->Fill(jetarray1);
            jet1 = fJetsCont->GetNextAcceptJet();
            fJetsCont2->ResetCurrentID();
            jet2 = fJetsCont2->GetNextAcceptJet();
        }
    }

    
    if (fJetsCont&&fJetsCont3) {
        fJetsCont->ResetCurrentID();
        AliEmcalJet *jet1 = fJetsCont->GetNextAcceptJet();
        fJetsCont3->ResetCurrentID();
        AliEmcalJet *jet3 = fJetsCont3->GetNextAcceptJet();
        while(jet1){
            while(jet3){
                Double_t deta = jet1->Eta()-jet3->Eta();
                Double_t dphi = RelativePhi(jet1->Phi(),jet3->Phi());
                Double_t deta2 = deta*deta;
                Double_t dphi2 = dphi*dphi;
                Double_t Aj = (jet1->Pt()-jet3->Pt())/(jet1->Pt()+jet3->Pt());
                Double_t jetarray[] = {jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->NEF(),jet3->Pt(),jet3->Eta(),jet3->Phi(),jet3->NEF(),Aj};
                //Marta used |dphi - pi|<pi/3
                if (fabs(fabs(dphi)-TMath::Pi())< TMath::Pi()/3.0){//dijet
                    fHistDiJet1->Fill(jetarray);
                    //we have a dijet, lets see if there is also a matched trigger
                    if (fJetsCont2) {
                        fJetsCont2->ResetCurrentID();
                        AliEmcalJet *jet2 = fJetsCont2->GetNextAcceptJet();
                        while(jet2){
                            Double_t tdeta = jet1->Eta()-jet2->Eta();
                            Double_t tdphi = RelativePhi(jet1->Phi(),jet2->Phi());
                            Double_t tdeta2 = tdeta*tdeta;
                            Double_t tdphi2 = tdphi*tdphi;
                            Double_t dR = pow(tdeta2+tdphi2,0.5);
                            
                            if (dR<0.2){
                                Double_t jetarray3[] = {jet1->Pt(),jet1->Eta(),jet1->Phi(),jet1->NEF(),jet3->Pt(),jet3->Eta(),jet3->Phi(),jet3->NEF(),Aj,jet2->Pt(),jet2->Eta(),jet2->Phi(),jet2->NEF(),dR};
                                //this dijet is triggered on!
                                fHistDiJet1m->Fill(jetarray3);
                            }
                            jet2 = fJetsCont2->GetNextAcceptJet();
                        } //while jet2
                    } // if jetscont2
                }// if dijet
                jet3 = fJetsCont3->GetNextAcceptJet();
            }//while jet 3
            jet1 = fJetsCont->GetNextAcceptJet();
            fJetsCont3->ResetCurrentID();
            jet3 = fJetsCont3->GetNextAcceptJet();
        }//while jet 1
    } //if jet cont 1 and 3

  return kTRUE;
}

//________________________________________________________________________
Float_t AliAnalysisTaskDcalDijetPerf:: RelativePhi(Double_t mphi,Double_t vphi) const
{
    if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
    else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
    if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
    else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
    double dphi = mphi-vphi;
    if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
    else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
    
    return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
void AliAnalysisTaskDcalDijetPerf::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskDcalDijetPerf::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskDcalDijetPerf::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
