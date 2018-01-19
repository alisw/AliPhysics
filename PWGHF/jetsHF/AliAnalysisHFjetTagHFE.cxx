// $Id$
//
// HF jet via HFE
//
// Author: S. Sakai

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

//
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include "AliAODInputHandler.h"
//


#include "AliAnalysisManager.h"
#include "AliVCluster.h"
//#include "AliAODEvent.h"  // added
//#include "AliAODHandler.h" // added 
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
//#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"
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
#include "AliMultSelection.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliAnalysisHFjetTagHFE.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisHFjetTagHFE)

//________________________________________________________________________
AliAnalysisHFjetTagHFE::AliAnalysisHFjetTagHFE() : 
  AliAnalysisTaskEmcalJet("AliAnalysisHFjetTagHFE", kTRUE),
  fVevent(0),
  fMultSelection(0),
  ftrack(0),
  fCaloClusters(0),
  fpidResponse(0),
  fcentMim(0), 
  fcentMax(10.0), 
  idbHFEj(kFALSE),
  iHybrid(kTRUE),
  fmimSig(-1.0),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fHistMultCent(0),
  fHistZcorr(0),
  fHistCent(0),
  fHistTPCnSigma(0),
  fHistEopNsig(0),
  fHistEop(0),
  fHistEopHad(0),
  fHistJetOrg(0),
  fHistJetOrgArea(0),
  fHistJetBG(0),
  fHistJetSub(0),
  fHisteJetOrg(0),
  fHisteJetBG(0),
  fHisteJetSub(0),
  fHistIncEle(0),
  fHistIncEleInJet0(0),
  fHistIncEleInJet1(0),
  fHistHfEleMC(0),
  fHistHfEleMCreco(0),
  fHistPhoEleMC(0),
  fHistPhoEleMCreco(0),
  fHistIncjet(0),
  fHistIncjetFrac(0),
  fHistIncjetOrg(0),
  fHistIncjetBG(0),
  fHistHFjet(0), 
  fHistULSjet(0), 
  fHistLSjet(0), 
  fInvmassULS(0),
  fInvmassLS(0),
  HFjetCorr0(0),
  HFjetCorr1(0),
  HFjetParticle(0),
  fQAHistJetPhi(0),
  fQAHistTrPhiJet(0),
  fQAHistTrPhi(0),
  fQAHistNits(0),
  fQAHistEleDCAxy(0),
  fQAHistEleDCAz(0),
  fHistClustE(0),
  fHistClustEtime(0),
  fEMCClsEtaPhi(0),
  fJetsCont(0),
  fJetsContPart(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fAOD(0),
  fMCarray(0),
  fMCparticle(0),
  fMCparticleMother(0)
  //fmcData(kFALSE)
{
  // Default constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisHFjetTagHFE::AliAnalysisHFjetTagHFE(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fVevent(0),
  fMultSelection(0),
  ftrack(0),
  fCaloClusters(0),
  fpidResponse(0),
  fcentMim(0), 
  fcentMax(10.0), 
  idbHFEj(kFALSE),
  iHybrid(kTRUE),
  fmimSig(-1.0),
  fHistTracksPt(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fHistMultCent(0),
  fHistZcorr(0),
  fHistCent(0),
  fHistTPCnSigma(0),//my
  fHistEopNsig(0),
  fHistEop(0),
  fHistEopHad(0),
  fHistJetOrg(0),
  fHistJetOrgArea(0),
  fHistJetBG(0),
  fHistJetSub(0),
  fHisteJetOrg(0),
  fHisteJetBG(0),
  fHisteJetSub(0),
  fHistIncEle(0),
  fHistIncEleInJet0(0),
  fHistIncEleInJet1(0),
  fHistHfEleMC(0),
  fHistHfEleMCreco(0),
  fHistPhoEleMC(0),
  fHistPhoEleMCreco(0),
  fHistIncjet(0),
  fHistIncjetFrac(0),
  fHistIncjetOrg(0),
  fHistIncjetBG(0),
  fHistHFjet(0),
  fHistULSjet(0),
  fHistLSjet(0),
  fInvmassULS(0),
  fInvmassLS(0),//my
  HFjetCorr0(0),
  HFjetCorr1(0),
  HFjetParticle(0),
  fQAHistJetPhi(0),
  fQAHistTrPhiJet(0),
  fQAHistTrPhi(0),
  fQAHistNits(0),
  fQAHistEleDCAxy(0),
  fQAHistEleDCAz(0),
  fHistClustE(0),
  fHistClustEtime(0),
  fEMCClsEtaPhi(0),
  fJetsCont(0),
  fJetsContPart(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  //
  fAOD(0),
  fMCarray(0),
  fMCparticle(0),
  fMCparticleMother(0)
  //fmcData(kFALSE)
{
  // Standard constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisHFjetTagHFE::~AliAnalysisHFjetTagHFE()
{
  // Destructor.
  delete ftrack;
  delete fCaloClusters;
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::UserCreateOutputObjects()
{
  // Create user output.

     if(idbHFEj)cout << "+++++++ MC check ++++++++ " << fmcData <<  endl;

  for(int i=0; i<5; i++)
     {
      if(idbHFEj)if(fJetCollArray.At(i))cout << " ------- " << i << " jet array " << endl;
     }

  if(idbHFEj)cout << "+++++++ MC ++++++++ " << fmcData <<  endl;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  // reconstructed
  fJetsCont           = GetJetContainer(0);
  if(idbHFEj)cout << "+++++++ jet get ++++++++" << fJetsCont <<  endl;
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  if(fTracksCont) fTracksCont->SetClassName("AliVTrack");
  if(fCaloClustersCont) fCaloClustersCont->SetClassName("AliVCluster");

  // particle
  fJetsContPart       = GetJetContainer(1);
 
  if(idbHFEj)cout << " fJetsCont :" <<  fJetsCont << endl;
  if(idbHFEj)cout << " fJetsContPart :" <<  fJetsContPart << endl;

  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fParticleCollArray.GetEntriesFast()>0) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      //fOutput->Add(fHistTracksPt[i]);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      //fOutput->Add(fHistClustersPt[i]);
    }

    if (fJetCollArray.GetEntriesFast()>0) {
      histname = "fHistLeadingJetPt_";
      histname += i;
      fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
      //fOutput->Add(fHistLeadingJetPt[i]);
      
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      //fOutput->Add(fHistJetsPhiEta[i]);
      
      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 30, 0, 3);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      //fOutput->Add(fHistJetsPtArea[i]);

      histname = "fHistJetsPtLeadHad_";
      histname += i;
      fHistJetsPtLeadHad[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistJetsPtLeadHad[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetZaxis()->SetTitle("counts");
      //fOutput->Add(fHistJetsPtLeadHad[i]);
    
      if (!(GetJetContainer()->GetRhoName().IsNull())) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	//fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }

  histname = "fHistPtDEtaDPhiTrackClus";
  fHistPtDEtaDPhiTrackClus = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);

  histname = "fHistPtDEtaDPhiClusTrack";
  fHistPtDEtaDPhiClusTrack = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fHistClustDx = new TH1F("fHistClustDx","fHistClustDx;Dx",1000,0.,1.);
  fOutput->Add(fHistClustDx);

  fHistClustDz = new TH1F("fHistClustDz","fHistClustDz;Dz",1000,0.,1.);
  fOutput->Add(fHistClustDz);

  fHistMultCent = new TH1F("fHistMultCent","centrality distribution",100,0,100);
  fOutput->Add(fHistMultCent);

  fHistZcorr = new TH2F("fHistZcorr","Z vertex corr V0 and SPD",100,-50,50,100,-50,50);
  fOutput->Add(fHistZcorr);

  fHistCent = new TH1F("fHistCent","centrality distribution",100,0,100);
  fOutput->Add(fHistCent);

  fHistTPCnSigma = new TH2F("fHistTPCnSigma","TPC nSigma;p_{T}(GeV/c);n#sigms",100,0.,20.,250,-5.,5.);
  fOutput->Add(fHistTPCnSigma);

  fHistEopNsig = new TH2F("fHistEopNsig","E/p vs. Nsig;Nsig;E/p",200,-5,5,200,0.,4.);
  fOutput->Add(fHistEopNsig);

  fHistEop = new TH2F("fHistEop","E/p;p_{T}(GeV/c);E/p",100,0.,20.,200,0.,4.);
  fOutput->Add(fHistEop);

  fHistEopHad = new TH2F("fHistEopHad","E/p hadron ;p_{T}(GeV/c);E/p",100,0.,20.,200,0.,4.);
  fOutput->Add(fHistEopHad);

  fHistJetOrg = new TH1F("fHistJetOrg","Inclusive jet org;p_{T}",300,-100.,200.);
  fOutput->Add(fHistJetOrg);

  fHistJetOrgArea = new TH2F("fHistJetOrgArea","Inclusive jet org vs. Area;p_{T};Area",300,-100.,200.,100,0,1);
  fOutput->Add(fHistJetOrgArea);

  fHistJetBG = new TH1F("fHistJetBG","BG jet;p_{T}",300,-100.,200.);
  fOutput->Add(fHistJetBG);

  fHistJetSub = new TH1F("fHistJetSub","Sub jet;p_{T}",300,-100.,200.);
  fOutput->Add(fHistJetSub);

  fHisteJetOrg = new TH1F("fHisteJetOrg","Inclusive jet org e;p_{T}",300,-100.,200.);
  fOutput->Add(fHisteJetOrg);

  fHisteJetBG = new TH1F("fHisteJetBG","BG jet e;p_{T}",300,-100.,200.);
  fOutput->Add(fHisteJetBG);

  fHisteJetSub = new TH1F("fHisteJetSub","Sub jet e;p_{T}",300,-100.,200.);
  fOutput->Add(fHisteJetSub);

  fHistIncEle = new TH1F("fHistIncEle","Inclusive electron;p_{T}",100,0.,20.);
  fOutput->Add(fHistIncEle);

  fHistHfEleMC = new TH1F("fHistHfEleMC","HF electron;p_{T}",100,0.,20.);
  fOutput->Add(fHistHfEleMC);

  fHistHfEleMCreco = new TH1F("fHistHfEleMCreco","HF reco electron;p_{T}",100,0.,20.);
  fOutput->Add(fHistHfEleMCreco);

  fHistIncEleInJet0 = new TH1F("fHistIncEleInJet0","Inclusive electron in Jet;p_{T}",100,0.,20.);
  fOutput->Add(fHistIncEleInJet0);

  fHistIncEleInJet1 = new TH1F("fHistIncEleInJet1","Inclusive electron in Jet;p_{T}",100,0.,20.);
  fOutput->Add(fHistIncEleInJet1);

  fHistPhoEleMC = new TH1F("fHistPhoEleMC","Photonic e MC reco;p_{T}",100,0.,20.);
  fOutput->Add(fHistPhoEleMC);

  fHistPhoEleMCreco = new TH1F("fHistPhoEleMCreco","Photonic e MC reco;p_{T}",100,0.,20.);
  fOutput->Add(fHistPhoEleMCreco);

  fHistIncjet = new TH2F("fHistIncjet","Inc jet;p_{T}",20,0,20,150,0.,150.);
  fOutput->Add(fHistIncjet);

  fHistIncjetFrac = new TH2F("fHistIncjetFrac","Inc jet e frac ;p_{T}",20,0,20,150,0.,1.5);
  fOutput->Add(fHistIncjetFrac);

  fHistIncjetOrg = new TH2F("fHistIncjetOrg","Inc jet org;p_{T}",20,0,20,300,-100.,200.);
  fOutput->Add(fHistIncjetOrg);

  fHistIncjetBG = new TH2F("fHistIncjetBG","Inc BG jet;p_{T}",20,0,20,300,-100.,200.);
  fOutput->Add(fHistIncjetBG);

  fHistHFjet = new TH2F("fHistHFjet","HF jet;p_{T}",20,0,20,300,-100.,200.);
  fOutput->Add(fHistHFjet);

  fHistULSjet = new TH2F("fHistULSjet","ULS jet;p_{T}",20,0,20,300,-100.,200.);
  fOutput->Add(fHistULSjet);

  fHistLSjet = new TH2F("fHistLSjet","LS jet;p_{T}",20,0,20,300,-100.,200.);
  fOutput->Add(fHistLSjet);

  fInvmassULS = new TH2F("fInvmassULS","ULS mass;p_{T};mass",20,0,20,150,0.,0.3);
  fOutput->Add(fInvmassULS);

  fInvmassLS = new TH2F("fInvmassLS","LS mass;p_{T};mass",20,0,20,150,0.,0.3);
  fOutput->Add(fInvmassLS);

  // jet
  Int_t jetpTMax = 300;
  Int_t nBine[7] =  { 50, 50, jetpTMax, jetpTMax, jetpTMax, 100, jetpTMax};
  Double_t mimHFj[7] = {  0,   0,   0,   0,  0, 0, 0 };
  Double_t maxHFj[7] = {50, 50, (Double_t)jetpTMax, (Double_t)jetpTMax, (Double_t)jetpTMax, 1, (Double_t)jetpTMax};

  HFjetCorr0 = new THnSparseD("HFjetCorr0","HF MC Corr;p_{T}^{reco}; p_{T}^{MC}; jet_{reco}; jet_{MC}; jet_{particle}; R match; pThaed;", 7, nBine, mimHFj, maxHFj);
  HFjetCorr0->Sumw2();
  fOutput->Add(HFjetCorr0);

  HFjetCorr1 = new THnSparseD("HFjetCorr1","HF MC Corr;p_{T}^{reco}; p_{T}^{MC}; jet_{reco}; jet_{MC};  jet_{particle}; R match; pThard;", 7, nBine, mimHFj, maxHFj);
  HFjetCorr1->Sumw2();
  fOutput->Add(HFjetCorr1);

  HFjetParticle = new THnSparseD("HFjetParticle","HF particle;p_{T}^{reco}; p_{T}^{MC}; jet_{reco}; jet_{MC};  jet_{particle}; R match; pThard;", 7, nBine, mimHFj, maxHFj);
  HFjetParticle->Sumw2();
  fOutput->Add(HFjetParticle);

  // QA
  fQAHistJetPhi = new TH1F("fQAHistJetPhi","jet phi",650,0.0,6.5);
  fOutput->Add(fQAHistJetPhi);

  fQAHistTrPhiJet = new TH1F("fQAHistTrPhiJet","track phi in Jet",650,0.0,6.5);
  fOutput->Add(fQAHistTrPhiJet);

  fQAHistTrPhi = new TH1F("fQAHistTrPhi","track phi",650,0.0,6.5);
  fOutput->Add(fQAHistTrPhi);
 
  fQAHistNits = new TH1F("fQAHistNits","ITS hits",7,-0.5,6.5);
  fOutput->Add(fQAHistNits);

  fQAHistEleDCAxy = new TH2F("fQAHistEleDCAxy","pT ele check DCAxy",40,0,20,200,-10,10);
  fOutput->Add(fQAHistEleDCAxy);

  fQAHistEleDCAz = new TH2F("fQAHistEleDCAz","pT ele check DCAz",40,0,20,200,-10,10);
  fOutput->Add(fQAHistEleDCAz);

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
  fOutput->Add(fHistClustE);

  fHistClustEtime = new TH1F("fHistClustEtime", "EMCAL cluster energy distribution with time; Cluster E;counts", 500, 0.0, 50.0);
  fOutput->Add(fHistClustEtime);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",1800,-0.9,0.9,630,0,6.3);
  fOutput->Add(fEMCClsEtaPhi);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisHFjetTagHFE::FillHistograms()
{
  // Fill histograms.
  if(idbHFEj)cout << " +++ Fill histograms " << endl;



  if (fTracksCont) {
    //AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0)); 
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
    while(track) {
      fHistTracksPt[fCentBin]->Fill(track->Pt()); 
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }

  if (fCaloClustersCont) {
    //AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0); 
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster();
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());
      Double_t dx = cluster->GetTrackDx();
      Double_t dz = cluster->GetTrackDz();
      fHistClustDx->Fill(dx);
      fHistClustDz->Fill(dz);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }

  if(idbHFEj)cout << "JetsCont : " << fJetsCont << endl;
  if(idbHFEj)cout << "Rho Name : " << fJetsCont->GetRhoName() << endl;
  if(idbHFEj)cout << "Rho Param : " << fJetsCont->GetRhoParameter() << endl;

  if (fJetsCont) {
    //AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0); 
    fJetsCont->ResetCurrentID();
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(); 
    while(jet) {

       //cout << "# of jets : " << jet->GetNumberOfTracks() << endl;

      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

      Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);

      if (fHistJetsCorrPtArea[fCentBin]) {
	Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }
    
       // track
        for (unsigned j = 0; j< jet->GetNumberOfTracks(); j++)
            { 
             AliVParticle *jetcont;
             jetcont = static_cast<AliVParticle*>(jet->TrackAt(j, fTracks));
             //if(idbHFEj)cout << "+++jetcont" << jetcont << endl;
             //if(idbHFEj)cout << "jet mom = " << jetcont->Px() << " ; " << jetcont->Py() << " ; " << jetcont->Pz() << endl;
             //if(jetcont) continue;
             

             //if(idbHFEj)cout << "jet mom = " << jetcont->Px() << " ; " << jetcont->Py() << " ; " << jetcont->Pz() << endl;
            }

       //

       jet = fJetsCont->GetLeadingJet();
       if(jet) fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());

       jet = fJetsCont->GetNextAcceptJet(); 
    }
    
  }

  //CheckClusTrackMatching();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::CheckClusTrackMatching()
{
  if(idbHFEj)cout << "< --------- CheckClusTrackMatching"<<endl;  

  if(!fTracksCont || !fCaloClustersCont)
    return;

  Double_t deta = 999;
  Double_t dphi = 999;

  //Get closest cluster to track
  //AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0));
  fTracksCont->ResetCurrentID(); 
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
  while(track) {
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(clusMatch) {
	AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
	fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
      }
    }
    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  }
  
  //Get closest track to cluster
  //AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0);
  fCaloClustersCont->ResetCurrentID(); 
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
  while(cluster) {
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    fHistClustersPt[fCentBin]->Fill(nPart.Pt());
    
    //Get matched track
    AliVTrack *mt = NULL;      
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(acl) {
      if(acl->GetNTracksMatched()>1)
	mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
    }
    else {
      AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
      Int_t im = ecl->GetTrackMatchedIndex();
      if(fTracksCont && im>=0) {
	mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
      }
    }
    if(mt) {
      AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
      
      /* //debugging
	 if(mt->IsEMCAL()) {
	 Int_t emc1 = mt->GetEMCALcluster();
	 Printf("current id: %d  emc1: %d",fCaloClustersCont->GetCurrentID(),emc1);
	 AliVCluster *clm = fCaloClustersCont->GetCluster(emc1);
	 AliPicoTrack::GetEtaPhiDiff(mt, clm, dphi, deta);
	 Printf("deta: %f dphi: %f",deta,dphi);
	 }
      */
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::ExecOnce() {

  
  if(idbHFEj)cout << "<------ ExecOnce: HFtagHFE " << endl;
  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

  if(idbHFEj)cout << "<------ End:ExecOnce: HFtagHFE " << endl;

}

//________________________________________________________________________
Bool_t AliAnalysisHFjetTagHFE::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().
  
  if(idbHFEj)cout <<  endl;
  if(idbHFEj)cout << "++++++++++++++++ " << endl;
  if(idbHFEj)cout << "Run!" << endl;
  if(idbHFEj)cout << fJetsCont << endl;
  if(idbHFEj)cout << fJetsContPart << endl;
  if(idbHFEj)cout << " fJetsCont :" <<  fJetsCont->GetName()<< " ; N = "<< fJetsCont->GetNAcceptedJets() << endl;
  if(idbHFEj)cout << " fJetsContPart :" <<  fJetsContPart->GetName() << " N = " << fJetsContPart->GetNAcceptedJets() << endl;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());

  // centrality
  /*
  Double_t centrality = -1;
  AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality(); 
  centrality = fCentrality->GetCentralityPercentile("V0M");
  */

  if(idbHFEj)cout << "Run number = " << fAOD->GetRunNumber() << endl;
 
  Float_t lPercentile = 300; 
  if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  if( !fMultSelection) {
   //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
   lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
 }

  Double_t centrality = -1;
  centrality = fMultSelection->GetMultiplicityPercentile("V0M", false); 
  //if(idbHFEj)cout << "cent = " << centrality << endl; 
  //if(idbHFEj)cout << "cent cut: " << fcentMim << " ; " << fcentMax << endl; 
  //if(centrality<0.0 || centrality>fcentMax)break;
  fHistMultCent->Fill(centrality);

  // vertex
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  double Zvertex = pVtx->GetZ();  
  double Yvertex = pVtx->GetY();  
  double Xvertex = pVtx->GetX();

  const AliAODVertex* SpdVtx = fAOD->GetPrimaryVertexSPD();
  double ZvertexSPD = SpdVtx->GetZ();  

  fHistZcorr->Fill(Zvertex,ZvertexSPD);  

  double del_Z = ZvertexSPD - Zvertex;
  if(idbHFEj)cout << "Zvertex = " << Zvertex << " ; SPD vertex" << ZvertexSPD  << endl;

  //PID initialised//
  //AliPIDResponse *fpidResponse = fInputHandler->GetPIDResponse();
  fpidResponse = fInputHandler->GetPIDResponse();
   //cout << "check PID ..." << endl;

  // track
  //ftrack = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
  ftrack = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
  int ntracks = 0;
  if(ftrack)
    {
     ntracks = ftrack->GetEntries();
    }
  else
   {
     ntracks = fAOD->GetNumberOfTracks();
   }
   if(idbHFEj)cout << "ftrack = " << ftrack << " ; " << ntracks << endl;
 
   //cout << "check track ..." << endl;

  // EMCal
  //TClonesArray* fCaloClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters")); 
  TClonesArray* fCaloClusters = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
  //cout << "check cluster ..." << endl;

   // MC array
   fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
   //cout << "check fMCarray ..." << endl;
   

 // analysis

  Bool_t fcentID = kFALSE;
  if(fcentMim<-10)
    {
     fcentID = kTRUE; // pp 
    } 
  else
    {
     if(centrality>fcentMim && centrality<fcentMax)fcentID = kTRUE;
    }

  //if(TMath::Abs(Zvertex)<10.0 && TMath::Abs(del_Z)<0.1 && (centrality>fcentMim && centrality<fcentMax)) // event cuts
  if(TMath::Abs(Zvertex)<10.0 && fcentID) // event cuts
    {
     //cout << "cent cut = " << centrality << endl; 
     fHistCent->Fill(centrality);    

     // MC (particle level Jet)
     if(fmcData)
       {
        MakeParticleLevelJet();
       }
     if(idbHFEj)cout << "check fmcData ..." << endl;

  /////////////////////////////
  //EMCAL cluster information//
  /////////////////////////////
  Int_t Nclust = -999;
  Nclust = fCaloClusters->GetEntries();

  for(Int_t icl=0; icl<Nclust; icl++)
  {
    AliVCluster *clust = 0x0;
    clust = dynamic_cast<AliVCluster*>(fCaloClusters->At(icl));
    if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

    if(clust && clust->IsEMCAL())
    {
      Float_t  emcx[3]; // cluster pos
      clust->GetPosition(emcx);
      TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
      Double_t emcphi = clustpos.Phi();
      Double_t emceta = clustpos.Eta();
      if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.

      //if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
      if(emcphi > 4.53 && emcphi < 5.708) continue; //DCAL  : 260 < phi < 327

      Float_t tof = clust->GetTOF()*1e+9; // ns

      fEMCClsEtaPhi->Fill(emceta,emcphi);

      Double_t clustE = clust->E();
      fHistClustE->Fill(clustE);
      if(tof>-30 && tof<30)fHistClustEtime->Fill(clustE);
    }
  }

     

     // inclusive jet
     double rho = 0.0;
     int Ncon = 0;
     if (fJetsCont) 
        {
         //AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0);
         fJetsCont->ResetCurrentID();
         AliEmcalJet *jet = fJetsCont->GetNextAcceptJet();
         rho = fJetsCont->GetRhoVal();
         if(idbHFEj)cout << "rho = " << rho << endl; 

         while(jet) {

         // check Raw jet info
         double jetpT = jet->Pt();
         double Rho_area = fJetsCont->GetRhoVal() * jet->Area();
         double jetpTsub = jetpT - Rho_area;
         double jetEta = jet->Eta();
         double jetPhi = jet->Phi();
         int Ncont = jet->GetNumberOfConstituents();
         double Jarea = jet->Area();       

         if(idbHFEj)cout << "Ncont = " << Ncont << endl;
         //if(Ncont<2)continue;

         fQAHistJetPhi->Fill(jetPhi); // QA

         fHistJetOrgArea->Fill(jetpT,Jarea);

         //if(fabs(jetEta)<0.6 && Ncont>2 && Jarea>0.2)  // 0.2 for 0.3
         if(fabs(jetEta)<0.6 && Ncont>2)  // 0.2 for 0.3
           {
            fHistJetOrg->Fill(jetpT);
            fHistJetBG->Fill(Rho_area);
            fHistJetSub->Fill(jetpTsub);
            }

            for (unsigned j = 0; j< jet->GetNumberOfTracks(); j++) 
                 {
                  AliVParticle *jetcont;
                  jetcont = static_cast<AliVParticle*>(jet->TrackAt(j, fTracks));
                  Double_t TrPhiJet = jetcont->Phi();
                  fQAHistTrPhiJet->Fill(TrPhiJet);
                  }
 
         jet = fJetsCont->GetNextAcceptJet(); 
        }
     
     }

   if(idbHFEj)cout <<"finished check jet" << endl;

  // --- Look for kink mother for AOD//
  Int_t numberofvertices = 100;
  if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
      AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
      if(!aodvertex) continue;
      if(aodvertex->GetType()==AliAODVertex::kKink) {
        AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
        if(!mother) continue;
        Int_t idmother = mother->GetID();
        listofmotherkink[numberofmotherkink] = idmother;
        numberofmotherkink++;
    }
  } //+++

    // ----- track loop

    //Bool_t ISelectronEv = kFALSE;
    Int_t ISelectronEv = 0;

    for (Int_t itrack = 0; itrack < ntracks; itrack++) {

        //cout << "tracks = " << itrack << " ; " << ntracks << endl;

        AliVParticle* ptrack = 0x0;
        if(ftrack)ptrack = dynamic_cast<AliVTrack*>(ftrack->At(itrack));
        if(!ftrack)ptrack = fAOD->GetTrack(itrack);
        AliVTrack *track = dynamic_cast<AliVTrack*>(ptrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(track);  // to apply cuts

        if(idbHFEj)cout<< "tarck label = " << track->GetLabel() << endl;

        int MCpdg = 0;
        fMCparticle = 0x0;
        fMCparticleMother = 0x0;

        if(fmcData && track->GetLabel()!=0)
          {
	   //fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
	   fMCparticle = (AliAODMCParticle*) fMCarray->At(TMath::Abs(track->GetLabel()));
           MCpdg = fMCparticle->GetPdgCode();
          }
        if(idbHFEj)cout << "MCpdg = " << MCpdg << endl;

        Bool_t isElectron = kFALSE;
        Bool_t fFlagULS=kFALSE;
        Bool_t fFlagLS=kFALSE;
        Bool_t iMCHF = kFALSE;  // b->e + c->e
        Bool_t iMCPHO = kFALSE; // g->e + pi0->e + eta->e;
        Double_t epTarray[3];
        Double_t epTarrayMC[3]; 
        for(int i=0; i<3; i++)
           {
            epTarray[i] = 0.0;
            epTarrayMC[i] = 0.0; 
           }

        // get track information
        Double_t pt = track->Pt(); 
        Double_t px = track->Px(); 
        Double_t py = track->Py(); 
        Double_t pz = track->Pz(); 
        Double_t eta = track->Eta(); 
        Double_t phi = track->Phi(); 
        Double_t d0z0[2]={-999,-999}, cov[3];
          if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
        //cout << "DCA = " << d0z0[0] << " ; " << d0z0[1] << endl;

        if(fabs(eta)>0.6)continue;

        if(TMath::Abs(MCpdg)==11)
           {
	     if(fMCparticle->GetMother()>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
	     Int_t pdgMom = fMCparticleMother->GetPdgCode();
             iMCHF  = isHeavyFlavour(pdgMom);
            }

        fQAHistTrPhi->Fill(phi); // QA
        fQAHistNits->Fill(atrack->GetITSNcls());

        
        if(iHybrid)
          {
           if(idbHFEj)cout << "Hybrid" << endl;
           if(!atrack->IsHybridGlobalConstrainedGlobal()) continue; // AOD track level
          }
        else
          {
           if(idbHFEj)cout << "non Hybrid" << endl;
           if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; // AOD track level
          }

        //cout << "track cuts ....." << endl;

        if(TMath::Abs(MCpdg)==11 && iMCHF)
           {
            fQAHistEleDCAxy->Fill(pt,d0z0[0]);
            fQAHistEleDCAz->Fill(pt,d0z0[1]);
           }

        //if(fabs(eta)>0.6)continue;
        if(fabs(d0z0[0])>3.0)continue;
        if(fabs(d0z0[1])>3.0)continue;
        if(track->GetTPCNcls() < 80) continue;
        //if(atrack->GetITSNcls() < 2) continue;   // AOD track level
        if(atrack->GetITSNcls() < 0.9) continue;   // AOD track level
        if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) continue;    // kAny
        if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        // kink cut
        //cout << "Kink cut" << endl;
	Bool_t kinkmotherpass = kTRUE;
	for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
		if(track->GetID() == listofmotherkink[kinkmother]) {
			kinkmotherpass = kFALSE;
			continue;
		}
	      }
        if(!kinkmotherpass) continue;


        // Get TPC nSigma
        Double_t dEdx =-999, fTPCnSigma=-999;
        dEdx = track->GetTPCsignal();
        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

        if(idbHFEj)cout << "fmimSig = " << fmimSig << endl;
        //if(fTPCnSigma<fmimSig || fTPCnSigma>3)continue;  //++++++++
        fHistTPCnSigma->Fill(pt,fTPCnSigma);

               epTarray[0] = px;
               epTarray[1] = py;
               epTarray[2] = pz;
        
        // Get E/p
        Int_t EMCalIndex = -1;
        EMCalIndex = track->GetEMCALcluster();
        if(EMCalIndex < 0) continue;
        
        AliVCluster *clustMatch=0x0;
        //if(!fUseTender) clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
        //clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
        if(fCaloClusters)
           {
            clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters->At(EMCalIndex));  // updated by tender 
           }
        else
          {
            clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
          }
           

        if(clustMatch && clustMatch->IsEMCAL())
        {
            /////////////////////////////////////////////
            //Properties of tracks matched to the EMCAL//
            /////////////////////////////////////////////
            if(idbHFEj)cout << "++++++++++ match EMCal  " << endl;
      
            Bool_t fClsTypeEMC = kFALSE;
            Float_t  emcx[3]; // cluster pos
            clustMatch->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            Double_t emcphi = clustpos.Phi();
            Double_t emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
            if(!fClsTypeEMC)continue;

            if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;
            
            Double_t clustMatchE = clustMatch->E();
            Double_t m20 =clustMatch->GetM20();
         
            if(m20<0.03 || m20>0.3)continue;  // shower shape cut
     
            //EMCAL EID info
            Double_t eop = -1.0;
            if(track->P()>0)eop = clustMatchE/track->P();
            if(idbHFEj)cout << "++++++++++ eop = " << eop << " ; " << pt  << endl;
            if(pt>2.0)fHistEopNsig->Fill(fTPCnSigma,eop);

            if(fTPCnSigma<-4)fHistEopHad->Fill(pt,eop);

            if(fTPCnSigma<fmimSig || fTPCnSigma>3)continue;  // Nsigma cut
            fHistEop->Fill(pt,eop);

            if(eop>0.9 && eop<1.3 && m20<0.3 && m20>0.03)isElectron = kTRUE;  
                 
            if(isElectron)
              {

               ISelectronEv++;
               SelectPhotonicElectron(itrack, track, fFlagULS, fFlagLS); // to fo ULS, LS

               if(abs(MCpdg)==11)
                 {
	          if(fMCparticle->GetMother()>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
	          Int_t pdgMom = fMCparticleMother->GetPdgCode();
                  iMCHF  = isHeavyFlavour(pdgMom);
                  iMCPHO = isPhotonic(pdgMom);
                  epTarrayMC[0] = fMCparticle->Px();
                  epTarrayMC[1] = fMCparticle->Py();
                  epTarrayMC[2] = fMCparticle->Pz();
                 }

              }
          }

    if(!isElectron)continue;

    if(idbHFEj)cout << "electron in event" << endl;

    // data
    fHistIncEle->Fill(pt);  
    // MC
    if(iMCHF)fHistHfEleMCreco->Fill(pt);  // to do ULS, LS
    if(iMCPHO)
      {  
    cout << "iMCPHO = " << iMCPHO << endl; 
       fHistPhoEleMC->Fill(pt);
    cout << "iMCPHO = " << iMCPHO << endl; 
       if(fFlagULS && !fFlagLS)fHistPhoEleMCreco->Fill(pt);
    cout << "iMCPHO = " << iMCPHO << endl; 
      }

    if(idbHFEj)cout << " ++++++ find e in jet " << endl;
    
    // MC true
    Double_t pTeJetTrue = 0.0;
    if(fmcData && fJetsContPart)
      {
        //AliEmcalJet *jetPart = fJetsContPart->GetNextAcceptJet(0);  // full or charge ?
        fJetsContPart->ResetCurrentID();
        AliEmcalJet *jetPart = fJetsContPart->GetNextAcceptJet();  // full or charge ?
        while(jetPart) 
             {
               Bool_t iTagHFjet = tagHFjet( jetPart, epTarrayMC, 0, pt);
	       if(iTagHFjet)pTeJetTrue = jetPart->Pt();
               jetPart = fJetsContPart->GetNextAcceptJet(); 
             }
      }
    
     if(idbHFEj)cout << "pTeJetTrue = " << pTeJetTrue << endl;

    // reco
    if (fJetsCont) 
       {
        //AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(0);  // full or charge ?
        fJetsCont->ResetCurrentID();
        AliEmcalJet *jet = fJetsCont->GetNextAcceptJet();  // full or charge ?
        while(jet) 
           {
            //Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);

            double jetEta = jet->Eta();
            //double jetEtacut = 0.9-0.3; // how get R size ?
            double jetEtacut = 0.6; // how get R size ?
            Bool_t iTagHFjet = tagHFjet( jet, epTarray, 0, pt);

            if(iTagHFjet && isElectron)
                 {
                  fHisteJetOrg->Fill(jet->Pt());
                  fHistIncEleInJet0->Fill(pt);  // to do ULS, LS
                 }

            if(fabs(jetEta)<jetEtacut)  // R cut already apply ?
              { 
               //Bool_t iTagHFjet = tagHFjet( jet, epTarray, 0, pt);
               if(idbHFEj)cout << "iTagHFjet = " << iTagHFjet << endl;
	       //Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	       Float_t pTeJet = jet->Pt();
	       Float_t pTeJetBG = rho * jet->Area();
               Float_t corrPt = pTeJet - pTeJetBG;
               Float_t efrac = pt/corrPt;
            
               if(ISelectronEv==1)
                 {
                  fHisteJetOrg->Fill(pTeJet);
                  fHisteJetBG->Fill(pTeJetBG);
                  fHisteJetSub->Fill(corrPt);
                 }
               
               if(iTagHFjet && isElectron) // TPC+EMCal
                   {
                    fHistIncEleInJet1->Fill(pt);  // to do ULS, LS
                    if(idbHFEj)cout << "Fill jet tag by HFE Reco(Data)" << endl;
                    fHistIncjetOrg->Fill(pt,pTeJet); 
                    fHistIncjetBG->Fill(pt,pTeJetBG); 
                    fHistIncjet->Fill(pt,corrPt);
                    fHistIncjetFrac->Fill(pt,efrac);
                   
                    if(!fFlagULS)fHistHFjet->Fill(pt,corrPt);
                    if(fFlagULS) fHistULSjet->Fill(pt,corrPt);
                    if(fFlagLS)fHistLSjet->Fill(pt,corrPt);
               
                    if(iMCHF)
                      {
                       double HFjetVals[7];
                       HFjetVals[0]=track->Pt(); HFjetVals[1]=0.0; HFjetVals[2] = corrPt; HFjetVals[3] = pTeJet; HFjetVals[4] = pTeJetTrue; HFjetVals[5] = 0.0; HFjetVals[6] = 0.0;
                       HFjetCorr1->Fill(HFjetVals);
                      }
                   } // teg by e

             } // jet eta cut
             jet = fJetsCont->GetNextAcceptJet(); 
           }  // while jet
       } // while jet   

   }  // end of Track loop

 } // end of event selection

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.


}


Bool_t AliAnalysisHFjetTagHFE::tagHFjet(AliEmcalJet* jetC, double *epT, int MCpid, double &maxpT_e)
{
  Bool_t HFjetTag = kFALSE;

  //if(idbHFEj)cout << "tagging ... electron mom = " << epT[0] << " ; " << epT[1] << " ; " << epT[2] << endl;
  //cout << "tagHFE:jet number = " << jetC->GetNumberOfTracks() << endl; 

  for (unsigned j = 0; j< jetC->GetNumberOfTracks(); j++) 
      {
       AliVParticle *jetcont;
       jetcont = static_cast<AliVParticle*>(jetC->TrackAt(j, fTracks));
       if(!jetcont) continue;
       //cout << "tagHFE:jet mom = " << jetcont->Px() << " ; " << jetcont->Py() << " ; " << jetcont->Pz() << endl;

       double Rmom[3];
       Rmom[0] = epT[0]-jetcont->Px();
       Rmom[1] = epT[1]-jetcont->Py();
       Rmom[2] = epT[2]-jetcont->Pz();
       double Rmatch = sqrt(pow(Rmom[0],2)+pow(Rmom[1],2)+pow(Rmom[2],2));
       //if(idbHFEj)cout << "dRmom = " << Rmatch << endl;

       //if(epT[0] == jetcont->Px() && epT[1] == jetcont->Py() && epT[2] == jetcont->Pz()) // electron in jet
       if(Rmatch<1e-8) // electron in jet
         {
          HFjetTag = kTRUE;
          if(idbHFEj)cout << "electron mom = " << epT[0] << " ; " << epT[1] << " ; " << epT[2] << endl;
          if(idbHFEj)cout << "tagHFE:jet mom = " << jetcont->Px() << " ; " << jetcont->Py() << " ; " << jetcont->Pz() << endl;
          if(idbHFEj)cout << "jet tag by HFE" << endl;
         }
     
      }
 return HFjetTag;
}


void AliAnalysisHFjetTagHFE::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagConvinatElec)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    Bool_t flagPhotonicElec = kFALSE;
    Bool_t flagConvinatElec = kFALSE;
    
    Int_t ntracks = -999;
    //if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    if(ftrack)
      {
       ntracks = ftrack->GetEntries();
      }
    else
      {
       ntracks = fAOD->GetNumberOfTracks();
      }     

    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
        AliVParticle* VAssotrack = 0x0;
        //if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
        if(ftrack)
          {
           VAssotrack = dynamic_cast<AliVTrack*>(ftrack->At(jtrack)); //take tracks from Tender list
          }        
        else
          { 
           VAssotrack  = fAOD->GetTrack(jtrack);
          }
 
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }

        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

        //------reject same track
        if(jtrack==itrack) continue;

        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;

        nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //------track cuts applied
        if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(aAssotrack->GetTPCNcls() < 70) continue;
        if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        
        //-------loose cut on partner electron
        //if(ptAsso <0.3) continue;
        if(ptAsso <0.15) continue;
        if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        //-------define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        //-------Get mass
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);

        if(fFlagLS)
            if(track->Pt()>1) fInvmassLS->Fill(track->Pt(),mass);
        if(fFlagULS)
            if(track->Pt()>1) fInvmassULS->Fill(track->Pt(),mass);
        
        //if(mass<100 && fFlagULS && !flagPhotonicElec) flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
        if(mass<0.1 && fFlagULS && !flagPhotonicElec) flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
        if(mass<0.1 && fFlagLS  && !flagConvinatElec) flagConvinatElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
    fFlagConvinatElec = flagConvinatElec;
}

//Bool_t isHeavyFlavour(int Mompdg)
Bool_t AliAnalysisHFjetTagHFE::isHeavyFlavour(int Mompdg)
{
 Bool_t iCharm = kFALSE;
 Bool_t iBeauty = kFALSE;
 Bool_t iHeavy = kFALSE;

      if(abs(Mompdg)==411 || abs(Mompdg)==413 || abs(Mompdg)==421 || abs(Mompdg)==423 || abs(Mompdg)==431)iCharm = kTRUE;
      if(abs(Mompdg)==511 || abs(Mompdg)==513 || abs(Mompdg)==521 || abs(Mompdg)==523 || abs(Mompdg)==531)iBeauty = kTRUE;
 if(iCharm || iBeauty)iHeavy = kTRUE;

 return iHeavy;
 
}

//Bool_t isPhotonic(int Mompdg)
Bool_t AliAnalysisHFjetTagHFE::isPhotonic(int Mompdg)
{
 Bool_t iphotonic = kFALSE;

      if(abs(Mompdg)==22 || abs(Mompdg)==111 || abs(Mompdg)==221)iphotonic = kTRUE;

 return iphotonic;
 
}

//void AliAnalysisHFjetTagHFE::MakeParticleLevelJet(THnSparse *pJet)
void AliAnalysisHFjetTagHFE::MakeParticleLevelJet()
{  
         if(idbHFEj)cout << "Making Particle Level Jet ..." << endl;
         if(idbHFEj)cout << fJetsContPart << endl;
        //fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));

	for(Int_t iMC = 0; iMC < fMCarray->GetEntries(); iMC++)
	{
         fMCparticle = 0x0;
         fMCparticleMother = 0x0;

	fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
        if(!fMCparticle)continue;
	Int_t pdg = fMCparticle->GetPdgCode();
        if(idbHFEj)cout << "pdg = " << pdg << endl;

	Int_t iMCmom =  fMCparticle->GetMother();
        if(idbHFEj)cout << "iMCmom " << iMCmom << endl;
	if(iMCmom>0) fMCparticleMother = (AliAODMCParticle*) fMCarray->At(iMCmom);
        if(idbHFEj)cout << fMCparticleMother  << endl;
	Int_t pdgMom = 0;
        if(fMCparticleMother)pdgMom = fMCparticleMother->GetPdgCode();
        if(idbHFEj)cout << "Mom = " << pdgMom << endl;
        Double_t etaMC = fMCparticle->Eta();
 

        if(fabs(pdg)==11 && pdgMom!=0 && TMath::Abs(etaMC)<0.6)
          {
           Bool_t iMCHF = isHeavyFlavour(pdgMom);
           if(iMCHF)
             {
              double MCpTarray[3];
              MCpTarray[0]=fMCparticle->Px(); 
              MCpTarray[1]=fMCparticle->Py(); 
              MCpTarray[2]=fMCparticle->Pz(); 
              cout << MCpTarray[0] << endl;
              double MChfepT=fMCparticle->Pt(); 
              cout << MChfepT << endl;
              fHistHfEleMC->Fill(MChfepT); 

              //if (fJetsCont) 
                  {
                   //AliEmcalJet *jetPart = fJetsContPart->GetNextAcceptJet(0);  // full or charge ?
                   fJetsContPart->ResetCurrentID();
                   AliEmcalJet *jetPart = fJetsContPart->GetNextAcceptJet();  // full or charge ?
                   while(jetPart) 
                    {
                     if(idbHFEj)cout << "jetPart = " << jetPart->Pt() << endl;
                     double jetEta = jetPart->Eta();
                     //double jetEtacut = 0.9-0.3; // how get R size ?
                     double jetEtacut = 0.6; // how get R size ?
                     if(fabs(jetEta)<jetEtacut)
                        {
                         Bool_t iTagHFjet = tagHFjet( jetPart, MCpTarray, 0, MChfepT);
                         //if(idbHFEj)cout << "iTagHFjet = " << iTagHFjet << endl;
                         if(iTagHFjet)
                           {
                            //if(idbHFEj)cout << "iTagHFjetMC = " << iTagHFjet << " ; " << jetPart->Pt() << endl;
                            cout << "iTagHFjetMC = " << iTagHFjet << " ; " << jetPart->Pt() << endl;
                            double HFjetVals[7];
                            HFjetVals[0]=0.0; HFjetVals[1]=MChfepT; HFjetVals[2] = 0.0; HFjetVals[3] = 0.0; HFjetVals[4] = jetPart->Pt(); HFjetVals[5] = 0.0; HFjetVals[6] = 0.0;
                            HFjetParticle->Fill(HFjetVals); 
                            //pJet->Fill(HFjetVals); 
                          }
                        }
                        if(idbHFEj)cout << "go to next jet" << endl;
                        jetPart = fJetsContPart->GetNextAcceptJet(); 
                    }
                 }     
             }
          }
        
       }
}

//________________________________________________________________________
void AliAnalysisHFjetTagHFE::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
