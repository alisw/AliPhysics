// Jet-Hadron Correlations PID
// Event plane dependence task.
//
// Author: J Mazer

// task head include
#include "AliAnalysisTaskEmcalJetHadEPpid.h"

// general ROOT includes
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TObjArray.h>

// AliROOT includes
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliAODJet.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include <AliVEvent.h>
#include <AliVParticle.h>
#include "AliRhoParameter.h"
#include "AliEmcalParticle.h"

// Localized Rho includes
#include "AliLocalRhoParameter.h"
#include "AliAnalysisTaskLocalRho.h"

// event handler (and pico's) includes
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliESDInputHandler.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"

// PID includes
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"

// magnetic field includes
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

ClassImp(AliAnalysisTaskEmcalJetHadEPpid)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHadEPpid::AliAnalysisTaskEmcalJetHadEPpid() : 
  AliAnalysisTaskEmcalJet("correlations",kFALSE), 
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0), fTrkBias(5), fClusBias(5), fTrkEta(0.9), fJetPtcut(15.0),
  fDoEventMixing(0), fMixingTracks(50000),
  doPlotGlobalRho(kFALSE), doVariableBinning(kFALSE),
  dovarbinTHnSparse(kFALSE), useAOD(kTRUE), doPID(kTRUE),
  fLocalRhoVal(0),
  fTracksName(""), fJetsName(""),
  event(0),
  isPItpc(0), isKtpc(0), isPtpc(0), nPIDtpc(0) ,   // pid TPC
  isPIits(0), isKits(0), isPits(0), nPIDits(0) ,  // pid ITS
  isPItof(0), isKtof(0), isPtof(0), nPIDtof(0) ,  // pid TOF
  fPoolMgr(0x0),
  fPIDResponse(0x0), fTPCResponse(),
  fESD(0), fAOD(0),
  fHistTPCdEdX(0), fHistITSsignal(0), //fHistTOFsignal(0),
  fHistRhovsCent(0), fHistNjetvsCent(0), fHistCentrality(0),
  fHistZvtx(0), fHistMult(0),
  fHistJetPhi(0), fHistTrackPhi(0),
  fHistJetHaddPhiIN(0), fHistJetHaddPhiOUT(0), fHistJetHaddPhiMID(0),
  fHistJetHaddPhiBias(0), fHistJetHaddPhiINBias(0), fHistJetHaddPhiOUTBias(0), fHistJetHaddPhiMIDBias(0),
  fHistMEdPHI(0),
  fHistTrackPtallcent(0),
  fHistJetEtaPhi(0), fHistJetHEtaPhi(0),
  fHistSEphieta(0), fHistMEphieta(0),
  fHistJetHaddPHI(0),
  fHistPID(0),
  fhnPID(0x0), fhnMixedEvents(0x0), fhnJH(0x0)
{
  // Default Constructor 
  for(Int_t ilab=0; ilab<4; ilab++){
    for(Int_t ipta=0; ipta<7; ipta++){
      fHistTrackEtaPhi[ilab][ipta]=0;
    }  // end of pt-associated loop
  } // end of lab loop
    
  for(Int_t itrackpt=0; itrackpt<9; itrackpt++){
    fHistJetHadbindPhi[itrackpt]=0;
    fHistJetHadbindPhiIN[itrackpt]=0;
    fHistJetHadbindPhiMID[itrackpt]=0;
    fHistJetHadbindPhiOUT[itrackpt]=0;
  } // end of trackpt bin loop

  for(Int_t icent = 0; icent<2; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    fHistAreavsRawPt[icent]=0;
    
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){
        fHistJetH[icent][iptjet][ieta]=0;
        fHistJetHBias[icent][iptjet][ieta]=0;
        fHistJetHTT[icent][iptjet][ieta]=0;
      } // end of eta loop
    } // end of pt-jet loop
  } // end of centrality loop

  // centrality dependent histo's
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistRhovsdEP[i]            = 0;
    fHistJetEtaPhiPt[i]         = 0;
    fHistJetEtaPhiPtBias[i]     = 0;
    fHistJetPtArea[i]           = 0;
    fHistJetPtAreaBias[i]       = 0;
    fHistJetPtNcon[i]           = 0;
    fHistJetPtNconBias[i]       = 0;
    fHistJetPtNconCh[i]         = 0;
    fHistJetPtNconBiasCh[i]     = 0;
    fHistJetPtNconEm[i]         = 0;
    fHistJetPtNconBiasEm[i]     = 0;
    fHistJetHaddPhiINcent[i]    = 0;
    fHistJetHaddPhiOUTcent[i]   = 0;
    fHistJetHaddPhiMIDcent[i]   = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  
  // define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetHadEPpid::AliAnalysisTaskEmcalJetHadEPpid(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0), fTrkBias(5), fClusBias(5), fTrkEta(0.9), fJetPtcut(15.0),
  fDoEventMixing(0), fMixingTracks(50000),
  doPlotGlobalRho(kFALSE), doVariableBinning(kFALSE),
  dovarbinTHnSparse(kFALSE), useAOD(kTRUE), doPID(kTRUE),
  fLocalRhoVal(0),
  fTracksName(""), fJetsName(""),
  event(0),
  isPItpc(0), isKtpc(0), isPtpc(0), nPIDtpc(0),    // pid TPC
  isPIits(0), isKits(0), isPits(0), nPIDits(0),   // pid ITS
  isPItof(0), isKtof(0), isPtof(0), nPIDtof(0),   // pid TOF
  fPoolMgr(0x0),
  fPIDResponse(0x0), fTPCResponse(),
  fESD(0), fAOD(0),
  fHistTPCdEdX(0), fHistITSsignal(0), //fHistTOFsignal(0),
  fHistRhovsCent(0), fHistNjetvsCent(0), fHistCentrality(0),
  fHistZvtx(0), fHistMult(0),
  fHistJetPhi(0), fHistTrackPhi(0),
  fHistJetHaddPhiIN(0), fHistJetHaddPhiOUT(0), fHistJetHaddPhiMID(0),
  fHistJetHaddPhiBias(0), fHistJetHaddPhiINBias(0), fHistJetHaddPhiOUTBias(0), fHistJetHaddPhiMIDBias(0),
  fHistMEdPHI(0),
  fHistTrackPtallcent(0),
  fHistJetEtaPhi(0), fHistJetHEtaPhi(0),
  fHistSEphieta(0), fHistMEphieta(0),
  fHistJetHaddPHI(0),
  fHistPID(0),
  fhnPID(0x0), fhnMixedEvents(0x0), fhnJH(0x0)
{   
  // Default Constructor 
  for(Int_t ilab=0; ilab<4; ilab++){
    for(Int_t ipta=0; ipta<7; ipta++){
      fHistTrackEtaPhi[ilab][ipta]=0;
    }  // end of pt-associated loop
  } // end of lab loop

  for(Int_t itrackpt=0; itrackpt<9; itrackpt++){
    fHistJetHadbindPhi[itrackpt]=0;
    fHistJetHadbindPhiIN[itrackpt]=0;
    fHistJetHadbindPhiMID[itrackpt]=0;
    fHistJetHadbindPhiOUT[itrackpt]=0;
  } // end of trackpt bin loop
    
  for(Int_t icent = 0; icent<2; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    fHistAreavsRawPt[icent]=0;
    
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){
        fHistJetH[icent][iptjet][ieta]=0;
        fHistJetHBias[icent][iptjet][ieta]=0;
        fHistJetHTT[icent][iptjet][ieta]=0;
      } // end of eta loop
    } // end of pt-jet loop
  } // end of centrality loop

  // centrality dependent histo's
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistRawJetPtvsTrackPt[i]   = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
    fHistJetPtvsdEP[i]          = 0;
    fHistJetPtvsdEPBias[i]      = 0;
    fHistRhovsdEP[i]            = 0;
    fHistJetEtaPhiPt[i]         = 0;
    fHistJetEtaPhiPtBias[i]     = 0;
    fHistJetPtArea[i]           = 0;
    fHistJetPtAreaBias[i]       = 0;
    fHistJetPtNcon[i]           = 0;
    fHistJetPtNconBias[i]       = 0;
    fHistJetPtNconCh[i]         = 0;
    fHistJetPtNconBiasCh[i]     = 0;
    fHistJetPtNconEm[i]         = 0;
    fHistJetPtNconBiasEm[i]     = 0;
    fHistJetHaddPhiINcent[i]    = 0;
    fHistJetHaddPhiOUTcent[i]   = 0;
    fHistJetHaddPhiMIDcent[i]   = 0;
  }

  SetMakeGeneralHistograms(kTRUE);

  // define input and output slots here
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//_______________________________________________________________________
AliAnalysisTaskEmcalJetHadEPpid::~AliAnalysisTaskEmcalJetHadEPpid()
{
  // destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHadEPpid::UserCreateOutputObjects()
{ // This is called ONCE!
  if (!fCreateHisto) return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  OpenFile(1); // do I need the 1?

  // create histograms that arn't array
  fHistTPCdEdX = new TH2F("TPCdEdX", "TPCdEdX", 2000, 0.0, 100.0, 500, 0, 500); 
  fHistITSsignal = new TH2F("ITSsignal", "ITSsignal", 2000, 0.0, 100.0, 500, 0, 500);
//  fHistTOFsignal = new TH2F("TOFsignal", "TOFsignal", 2000, 0.0, 100.0, 500, 0, 500);
  fHistRhovsCent             = new TH2F("RhovsCent",              "RhovsCent",             100, 0.0, 100.0, 500, 0, 500);
  fHistNjetvsCent            = new TH2F("NjetvsCent",             "NjetvsCent",            100, 0.0, 100.0, 100, 0, 100);
  fHistTrackPtallcent = new TH1F("fHistTrackPtallcent", "p_{T} distribution", 1000, 0.0, 100.0);
  fHistCentrality = new TH1F("fHistCentrality","centrality",100,0,100);
  fHistZvtx = new TH1F("fHistZvertex","z vertex",60,-30,30);
  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet #eta-#phi",900,-1.8,1.8,640,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron #Delta#eta-#Delta#phi",900,-1.8,1.8,640,-1.6,4.8);
  fHistSEphieta = new TH2F("fHistSEphieta", "Single Event #phi-#eta distribution", 64,-0.5*TMath::Pi(), 1.5*TMath::Pi(), 32,-1.8,1.8);  // was 64 bins
  fHistMEphieta = new TH2F("fHistMEphieta", "Mixed Event #phi-#eta distribution", 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 32,-1.8,1.8);  // was 64 bins
  fHistJetHaddPHI = new TH1F("fHistJetHaddPHI", "Jet-Hadron #Delta#varphi", 128,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fHistJetPhi = new TH1F("fHistJetPhi", "Jet #phi Distribution", 128, -2.0*TMath::Pi(), 2.0*TMath::Pi());
  fHistTrackPhi = new TH1F("fHistTrackPhi", "Track #phi Distribution", 128, -2.0*TMath::Pi(), 2.0*TMath::Pi());
  fHistJetHaddPhiIN = new TH1F("fHistJetHaddPhiIN","Jet-Hadron #Delta#varphi IN PLANE", 128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiOUT = new TH1F("fHistJetHaddPhiOUT","Jet-Hadron #Delta#varphi OUT PLANE",128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiMID = new TH1F("fHistJetHaddPhiMID","Jet-Hadron #Delta#varphi MIDDLE of PLANE",128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiBias = new TH1F("fHistJetHaddPhiBias","Jet-Hadron #Delta#varphi with bias",128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiINBias = new TH1F("fHistJetHaddPhiINBias","Jet-Hadron #Delta#varphi IN PLANE with bias", 128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiOUTBias = new TH1F("fHistJetHaddPhiOUTBias","Jet-Hadron #Delta#varphi OUT PLANE with bias",128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  fHistJetHaddPhiMIDBias = new TH1F("fHistJetHaddPhiMIDBias","Jet-Hadron #Delta#varphi MIDDLE of PLANE with bias",128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
  
  fHistPID = new TH1F("fHistPID","PID Counter",12, -0.5, 11.5);
  SetfHistPIDcounterLabels(fHistPID);
  fOutput->Add(fHistPID);

  char name[200];

  for(Int_t ipta=0; ipta<7; ++ipta){ 
    for(Int_t ilab=0; ilab<4; ++ilab){
      snprintf(name,200, "fHistTrackEtaPhi_%i_%i", ilab,ipta);
      fHistTrackEtaPhi[ilab][ipta] = new TH2F(name,name,400,-1,1,640,0.0,2.*TMath::Pi());
///      fOutput->Add(fHistTrackEtaPhi[ilab][ipta]);
    } // end of lab loop
  } // end of pt-associated loop


  // jet hadron (binned) correlations in dPHI
  for(Int_t itrackpt=0; itrackpt<9; itrackpt++){
    snprintf(name,200,"fHistJetHadbindPhi_%i",itrackpt);
    fHistJetHadbindPhi[itrackpt] = new TH1F(name,name,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHadbindPhi[itrackpt]);

    snprintf(name,200,"fHistJetHadbindPhiIN_%i",itrackpt);
    fHistJetHadbindPhiIN[itrackpt] = new TH1F(name,name,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHadbindPhiIN[itrackpt]);

    snprintf(name,200,"fHistJetHadbindPhiMID_%i",itrackpt);
    fHistJetHadbindPhiMID[itrackpt] = new TH1F(name,name,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHadbindPhiMID[itrackpt]);

    snprintf(name,200,"fHistJetHadbindPhiOUT_%i",itrackpt);
    fHistJetHadbindPhiOUT[itrackpt] = new TH1F(name,name,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHadbindPhiOUT[itrackpt]);
  } // end of trackpt bin loop

  for(Int_t icent = 0; icent<2; ++icent){
    snprintf(name,200,"fHistJetPt_%i",icent);
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPt[icent]);

    snprintf(name,200,"fHistJetPtBias_%i",icent);
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtBias[icent]);

    snprintf(name,200,"fHistJetPtTT_%i",icent);
    fHistJetPtTT[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtTT[icent]);

    snprintf(name,200,"fHistAreavsRawPt_%i",icent);
    fHistAreavsRawPt[icent] = new TH2F(name,name,100,0,1,200,0,200);
    //fOutput->Add(fHistAreavsRawPt[icent]);

    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){
        snprintf(name,200,"fHistJetH_%i_%i_%i",icent,iptjet,ieta);
        fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
        fOutput->Add(fHistJetH[icent][iptjet][ieta]);

        snprintf(name,200,"fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);
        fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
        fOutput->Add(fHistJetHBias[icent][iptjet][ieta]);

        snprintf(name,200,"fHistJetHTT_%i_%i_%i",icent,iptjet,ieta);
        fHistJetHTT[icent][iptjet][ieta]=new TH2F(name,name,64,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
        fOutput->Add(fHistJetHTT[icent][iptjet][ieta]);
      } // end of eta loop
    } // end of pt-jet loop
  } // end of centrality loop

  // strings for titles
  TString name1;
  TString title1;

  for (Int_t i = 0;i<6;++i){
    name1 = TString(Form("JetPtvsTrackPt_%i",i));
    title1 = TString(Form("Jet p_{T} vs Leading Track p_{T} cent bin %i",i));
    fHistJetPtvsTrackPt[i] = new TH2F(name1,title1,250,-50,200,250,0,50);
    fOutput->Add(fHistJetPtvsTrackPt[i]);

    name1 = TString(Form("RawJetPtvsTrackPt_%i",i));
    title1 = TString(Form("Raw Jet p_{T} vs Leading Track p_{T} cent bin %i",i));
    fHistRawJetPtvsTrackPt[i] = new TH2F(name1,title1,250,-50,200,250,0,50);
    fOutput->Add(fHistRawJetPtvsTrackPt[i]);

    name1 = TString(Form("TrackPt_%i",i));
    title1 = TString(Form("Track p_{T} cent bin %i",i));
    fHistTrackPt[i] = new TH1F(name1,title1,1000,0,100); // up to 200?
    fOutput->Add(fHistTrackPt[i]);   

    name1 = TString(Form("EP0_%i",i));
    title1 = TString(Form("EP VZero cent bin %i",i));
    fHistEP0[i] = new TH1F(name1,title1,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0[i]);

/*
    name1 = TString(Form("EP0A_%i",i));
    title1 = TString(Form("EP VZero cent bin %i",i));
    fHistEP0A[i] = new TH1F(name1,title1,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0A[i]);

    name1 = TString(Form("EP0C_%i",i));
    title1 = TString(Form("EP VZero cent bin %i",i));
    fHistEP0C[i] = new TH1F(name1,title1,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0C[i]);

    name1 = TString(Form("EPAvsC_%i",i));
    title1 = TString(Form("EP VZero cent bin %i",i));
    fHistEPAvsC[i] = new TH2F(name1,title1,144,-TMath::Pi(),TMath::Pi(),144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEPAvsC[i]);
*/

    name1 = TString(Form("JetPtvsdEP_%i",i));
    title1 = TString(Form("Jet p_{T} vs #DeltaEP cent bin %i",i));
    fHistJetPtvsdEP[i] = new TH2F(name1,title1,250,-50,200,288,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEP[i]);

    name1 = TString(Form("JetPtvsdEPBias_%i",i));
    title1 = TString(Form("Bias Jet p_{T} vs #DeltaEP cent bin %i",i));
    fHistJetPtvsdEPBias[i] = new TH2F(name1,title1,250,-50,200,288,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistJetPtvsdEPBias[i]);

    name1 = TString(Form("RhovsdEP_%i",i));
    title1 = TString(Form("#rho vs #DeltaEP cent bin %i",i));
    fHistRhovsdEP[i] = new TH2F(name1,title1,500,0,500,288,-2*TMath::Pi(),2*TMath::Pi());
    fOutput->Add(fHistRhovsdEP[i]);

    name1 = TString(Form("JetEtaPhiPt_%i",i));
    title1 = TString(Form("Jet #eta-#phi p_{T} cent bin %i",i));
    fHistJetEtaPhiPt[i] = new TH3F(name1,title1,250,-50,200,100,-1,1,64,-3.2,3.2);
    fOutput->Add(fHistJetEtaPhiPt[i]);

    name1 = TString(Form("JetEtaPhiPtBias_%i",i));
    title1 = TString(Form("Jet #eta-#phi p_{T} Bias cent bin %i",i));
    fHistJetEtaPhiPtBias[i] = new TH3F(name1,title1,250,-50,200,100,-1,1,64,-3.2,3.2);
    fOutput->Add(fHistJetEtaPhiPtBias[i]);

    name1 = TString(Form("JetPtArea_%i",i));
    title1 = TString(Form("Jet p_{T} Area cent bin %i",i));
    fHistJetPtArea[i] = new TH2F(name1,title1,250,-50,200,100,0,1);
    fOutput->Add(fHistJetPtArea[i]);

    name1 = TString(Form("JetPtAreaBias_%i",i));
    title1 = TString(Form("Jet p_{T} Area Bias cent bin %i",i));
    fHistJetPtAreaBias[i] = new TH2F(name1,title1,250,-50,200,100,0,1);
    fOutput->Add(fHistJetPtAreaBias[i]);

    name1 = TString(Form("JetPtNcon_%i",i));
    title1 = TString(Form("Jet p_{T} Ncon cent bin %i",i));
    fHistJetPtNcon[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    //fOutput->Add(fHistJetPtNcon[i]);

    name1 = TString(Form("JetPtNconBias_%i",i));
    title1 = TString(Form("Jet p_{T} NconBias cent bin %i",i));
    fHistJetPtNconBias[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    fOutput->Add(fHistJetPtNconBias[i]);

    name1 = TString(Form("JetPtNconCh_%i",i));
    title1 = TString(Form("Jet p_{T} NconCh cent bin %i",i));
    fHistJetPtNconCh[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    //fOutput->Add(fHistJetPtNconCh[i]);

    name1 = TString(Form("JetPtNconBiasCh_%i",i));
    title1 = TString(Form("Jet p_{T} NconBiasCh cent bin %i",i));
    fHistJetPtNconBiasCh[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    //fOutput->Add(fHistJetPtNconBiasCh[i]);

    name1 = TString(Form("JetPtNconEm_%i",i));
    title1 = TString(Form("Jet p_{T} NconEm cent bin %i",i));
    fHistJetPtNconEm[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    //fOutput->Add(fHistJetPtNconEm[i]);

    name1 = TString(Form("JetPtNconBiasEm_%i",i));
    title1 = TString(Form("Jet p_{T} NconBiasEm cent bin %i",i));
    fHistJetPtNconBiasEm[i] = new TH2F(name1,title1,250,-50,200,100,0,2000);
    //fOutput->Add(fHistJetPtNconBiasEm[i]);

    name1 = TString(Form("JetHaddPhiINcent_%i",i));
    title1 = TString(Form("Jet Hadron #Delta#varphi Distribution IN PLANE cent bin %i",i));
    fHistJetHaddPhiINcent[i] = new TH1F(name1,title1,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHaddPhiINcent[i]);
  
    name1 = TString(Form("JetHaddPhiOUTcent_%i",i));
    title1 = TString(Form("Jet Hadron #Delta#varphi Distribution OUT PLANE cent bin %i",i));
    fHistJetHaddPhiOUTcent[i] = new TH1F(name1,title1,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHaddPhiOUTcent[i]);

    name1 = TString(Form("JetHaddPhiMIDcent_%i",i));
    title1 = TString(Form("Jet Hadron #Delta#varphi Distribution MIDDLE of PLANE cent bin %i",i));
    fHistJetHaddPhiMIDcent[i] = new TH1F(name1,title1,128,-0.5*TMath::Pi(), 1.5*TMath::Pi());
    fOutput->Add(fHistJetHaddPhiMIDcent[i]);
  }

  // add histo's to output list
  fOutput->Add(fHistTPCdEdX);
  fOutput->Add(fHistITSsignal);
//  fOutput->Add(fHistTOFsignal);
  fOutput->Add(fHistRhovsCent);
  fOutput->Add(fHistNjetvsCent);
  fOutput->Add(fHistTrackPtallcent);
  fOutput->Add(fHistJetEtaPhi);
  //fOutput->Add(fHistTrackEtaPhi);
  fOutput->Add(fHistJetHEtaPhi);
  fOutput->Add(fHistCentrality);
  fOutput->Add(fHistZvtx);
  fOutput->Add(fHistJetPhi);
  fOutput->Add(fHistTrackPhi);
  fOutput->Add(fHistJetHaddPhiIN);
  fOutput->Add(fHistJetHaddPhiOUT);
  fOutput->Add(fHistJetHaddPhiMID);
  fOutput->Add(fHistJetHaddPhiBias);
  fOutput->Add(fHistJetHaddPhiINBias);
  fOutput->Add(fHistJetHaddPhiOUTBias);
  fOutput->Add(fHistJetHaddPhiMIDBias);
  fOutput->Add(fHistSEphieta);
  fOutput->Add(fHistMEphieta);
  fOutput->Add(fHistJetHaddPHI);

  // set up jet-hadron sparse
  UInt_t bitcoded = 0; // bit coded, see GetDimParams() below
  UInt_t cifras = 0;
  UInt_t bitcode = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9;
  fhnJH = NewTHnSparseF("fhnJH", bitcoded);
  fhnJH->Sumw2();
///  fOutput->Add(fhnJH);

  // ****************************** PID *****************************************************
  // set up PID handler
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) {
        AliError("Input handler needed");
	return;
  }

  // PID response object
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) {
        AliError("PIDResponse object was not created");
	return;
  }

  // ****************************************************************************************
/*
  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it
  Int_t nCentralityBins  = 8;
  Double_t centralityBins[9] = {0.0, 4., 9, 15, 25, 35, 55, 100.0,500.0};  
  Double_t centralityBins[nCentralityBins+1];
  for(Int_t ic=0; ic<nCentralityBins+1; ic++){
    if(ic==nCentralityBins) centralityBins[ic]=500;
    else centralityBins[ic]=10.0*ic; 
  }
*/
  // setup for Pb-Pb collisions
  Int_t nCentralityBins  = 100;
  Double_t centralityBins[nCentralityBins];
  for(Int_t ic=0; ic<nCentralityBins; ic++){
    centralityBins[ic]=1.0*ic;
  }

  fHistMult = new TH1F("fHistMult","multiplicity",nCentralityBins,centralityBins);
  fOutput->Add(fHistMult);

  // Event Mixing
  Int_t trackDepth = fMixingTracks;
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
  Int_t nZvtxBins  = 5+1+5;
  Double_t vertexBins[] = { -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};
  Double_t* zvtxbin = vertexBins;
  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);

  // set up event mixing sparse
  if(fDoEventMixing){
     cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7;
     fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);
     fhnMixedEvents->Sumw2();
     fOutput->Add(fhnMixedEvents);
  } // end of do-eventmixing

  Double_t xlowjetPT[] = {-50,-45,-40,-35,-30,-25,-20,-18,-16,-14,-12,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400};
  Double_t xlowtrPT[] = {0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.25,2.50,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.50,4.75,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,22.0,24.0,26.0,28.0,30.0,35.0,40.0,45.0,50.0,60.0,70.0,80.0,90.0,100.0};

  // number of bins you tell histogram should be (# in array - 1) because the last bin
  // is the right-most edge of the histogram 
  // i.e. this is for PT and there are 57 numbers (bins) thus we have 56 bins in our histo
  Int_t nbinsjetPT = sizeof(xlowjetPT)/sizeof(Double_t) - 1;
  Int_t nbinstrPT = sizeof(xlowtrPT)/sizeof(Double_t) - 1;

  // set up PID sparse
  if(doPID){
    bitcode = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 |
              1<<10 | 1<<11 | 1<<12 | 1<<13 | 1<<14 | 1<<15 | 1<<16 | 1<<17 | 1<<18 | 1<<19 |
              1<<20 | 1<<21 | 1<<22;
    fhnPID = NewTHnSparseFPID("fhnPID", bitcode);
    fhnPID->Sumw2();
    fOutput->Add(fhnPID);
  } // end of do-PID

  if(dovarbinTHnSparse){
    fhnPID->GetAxis(1)->Set(nbinsjetPT, xlowjetPT);
    fhnPID->GetAxis(2)->Set(nbinstrPT, xlowtrPT);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    TH2 *h2 = dynamic_cast<TH2*>(fOutput->At(i));
    if (h2){
      h2->Sumw2();
      continue;
    }
    TH3 *h3 = dynamic_cast<TH3*>(fOutput->At(i));
    if (h3){
      h3->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  PostData(1, fOutput);
}

//_________________________________________________________________________
void AliAnalysisTaskEmcalJetHadEPpid::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

//_________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetHadEPpid::Run()
{ // Main loop called for each event
  // TEST TEST TEST TEST for OBJECTS!
  if(!fLocalRho){
    AliError(Form("Couldn't get fLocalRho object, try to get it from Event based on name\n"));
    fLocalRho = GetLocalRhoFromEvent(fLocalRhoName);
  }
  if(!fTracks){
    AliError(Form("No fTracks object!!\n"));
    return kTRUE;
  }
  if(!fJets){
    AliError(Form("No fJets object!!\n"));
    return kTRUE;
  }

  // what kind of event do we have: AOD or ESD?
  if (dynamic_cast<AliAODEvent*>(InputEvent())) useAOD = kTRUE;
  else useAOD = kFALSE;

  // if we have ESD event, set up ESD object
  if(!useAOD){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      AliError(Form("ERROR: fESD not available\n"));
      return kTRUE;
    }
  }

  // if we have AOD event, set up AOD object
  if(useAOD){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
      AliError(Form("ERROR: fAOD not available\n"));
      return kTRUE;
    }
  }

  // get centrality
  Int_t centbin = GetCentBin(fCent);
  fHistCentrality->Fill(fCent); // won't be filled in pp collision (Keep this in mind!)

  // for pp analyses we will just use the first centrality bin
  if (centbin == -1) centbin = 0;

  // get vertex information
  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  Double_t zVtx=fvertex[2];
  fHistZvtx->Fill(zVtx);

  // get z-vertex bin
  //Int_t zVbin = GetzVertexBin(zVtx);

  // apply zVtx cut
  if(fabs(zVtx)>10.0) 
    return kTRUE;

  // apply cut to event on Centrality > 90%
  if(fCent>90) return kTRUE;

  // create pointer to list of input event
  TList *list = InputEvent()->GetList();
  if(!list) {
    AliError(Form("ERROR: list not attached\n"));
    return kFALSE;
  }

  // initialize TClonesArray pointers to jets and tracks
  TClonesArray *jets = 0;
  TClonesArray *tracks = 0; 

  // get Tracks object
  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracks));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracks->GetName()));
    return kTRUE;
  } // verify existence of tracks

  // get Jets object
  jets = dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if(!jets){
    AliError(Form("Pointer to jets %s == 0", fJets->GetName()));
    return kTRUE;
  } // verify existence of jets

  // get number of jets and tracks
  const Int_t Njets = jets->GetEntries(); 
  const Int_t Ntracks = tracks->GetEntries();
  if(Ntracks<1)   return kTRUE;
  if(Njets<1)	  return kTRUE;

  fHistMult->Fill(Ntracks);  // fill multiplicity distribution

  // initialize track parameters
  Int_t iTT=-1;
  Double_t ptmax=-10;

  // loop over tracks - to get hardest track (highest pt)
  for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++){
    AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
    if (!track) {
      AliError(Form("Couldn't get VTrack track %d\n", iTracks));        
      continue;
    } // verify existence of tracks

    if(TMath::Abs(track->Eta())>0.9) continue;
    if(track->Pt()<0.15) continue;
    //iCount++;
    if(track->Pt()>ptmax){
      ptmax=track->Pt();             // max pT track
      iTT=iTracks;                   // trigger tracks
    } // check if Pt>maxpt

    fHistTrackPhi->Fill(track->Phi());
    fHistTrackPtallcent->Fill(track->Pt());
  } // end of loop over tracks

  // find tracks for event (alternative), and loop over
  for (int i = 0;i<Ntracks;i++){
    AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
    if (!track) continue;

    if(TMath::Abs(track->Eta())>0.9) continue;
    if(track->Pt()<0.15) continue;

    fHistTrackPt[centbin]->Fill(track->Pt());
  }

  // event plane histograms filled
  fHistEP0[centbin]->Fill(fEPV0);
  //fHistEP0A[centbin]->Fill(fEPV0A);
  //fHistEP0C[centbin]->Fill(fEPV0C);
  //fHistEPAvsC[centbin]->Fill(fEPV0A,fEPV0C);

  // get rho from event and fill relative histo's
  fRho = GetRhoFromEvent(fRhoName);
  fRhoVal = fRho->GetVal();
  fHistRhovsCent->Fill(fCent,fRhoVal);        // Global Rho vs Centrality
  fHistRhovsdEP[centbin]->Fill(fRhoVal,fEPV0); // Global Rho vs delta Event Plane angle

  // initialize jet parameters
  Int_t ijethi=-1;
  Double_t highestjetpt=0.0;
  Int_t passedTTcut=0;
  Int_t NjetAcc = 0;

  // loop over jets in an event - to find highest jet pT and apply some cuts
  for (Int_t ijet = 0; ijet < Njets; ijet++){
    // get our jets
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
    if (!jet) continue;

    // cuts on jets
    if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
      continue;
    if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
      continue;
    fHistAreavsRawPt[centbin]->Fill(jet->Pt(),jet->Area());
    if (jet->Area()<fAreacut)
      continue;
    //prevents 0 area jets from sneaking by when area cut == 0
    if (jet->Area()==0)
      continue;
    //exclude jets with extremely high pt tracks which are likely misreconstructed
    if(jet->MaxTrackPt()>100)
      continue;
    if(!AcceptJet(jet))
      continue;

    NjetAcc++;                     // # of accepted jets
    fHistJetPhi->Fill(jet->Phi()); // Jet Phi histogram (filled)

    // get highest Pt jet in event
    if(highestjetpt<jet->Pt()){
      ijethi=ijet;
      highestjetpt=jet->Pt();
    }
  } // end of looping over jets

  // accepted jets
  fHistNjetvsCent->Fill(fCent,NjetAcc);
  Int_t NJETAcc = 0;

  // loop over jets in event and make appropriate cuts
  for (Int_t ijet = 0; ijet < Njets; ++ijet) {
     AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ijet));
     if (!jet)  // see if we have a jet
       continue;
     if(!AcceptJet(jet))             // do we accept this jet?
       continue;

//     Double_t leadjet=0;
//     if (ijet==ijethi) leadjet=1;

     // cuts on jet
     if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
        continue;
     if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
        continue;
     if (jet->Area()==0) // make sure jet has an area
       continue;
     if (jet->Pt()<0.1) // (should probably be higher..., but makes a cut on jet pT)
       continue;
     if (jet->MaxTrackPt()>100) // elminates fake tracks
       continue;
     if (!AcceptJet(jet)) // sees if jet is accepted
       continue;
     //  jets.push_back(jet);

     // check on lead jet
     Double_t leadjet=0;
     if (ijet==ijethi) leadjet=1;

     // initialize and calculate various parameters: pt, eta, phi, rho, etc...
     Double_t jetphi = jet->Phi();      // phi of jet
     NJETAcc++;   // # accepted jets
     fLocalRhoVal = fLocalRho->GetLocalVal(jetphi, GetJetRadius(0)); //fmyJetRadius);   // get local rho value
     Double_t jeteta = jet->Eta();     // ETA of jet
//     Double_t jetptraw = jet->Pt();    // raw pT of jet
     Double_t jetPt = -500;            // initialize corr jet pt
//     jetPt = jet->Pt();
//     jetPt = jet->Pt()-jet->Area()*fRhoVal;  // corrected pT of jet from rho value
     jetPt = jet->Pt()-jet->Area()*fLocalRhoVal; // corrected pT of jet using Rho modulated for V2 and V3
     Float_t dEP = -500;                    // initialize angle between jet and event plane
     //Int_t ieta = -500;             // initialize jet eta bin
     //ieta = GetEtaBin(jeteta);         // bin of eta
     dEP = RelativeEPJET(jetphi,fEPV0); // angle betweeen jet and event plane

     // make histo's
     fHistJetPtvsTrackPt[centbin]->Fill(jetPt,jet->MaxTrackPt());
     fHistRawJetPtvsTrackPt[centbin]->Fill(jet->Pt(),jet->MaxTrackPt());
//     fHistJetPtvsdEP[centbin]->Fill(jetPt,RelativePhi((fEPV0+TMath::Pi()),jet->Phi()));
     fHistJetPtvsdEP[centbin]->Fill(jetPt, dEP);
     fHistJetEtaPhiPt[centbin]->Fill(jetPt,jet->Eta(),jet->Phi());
     fHistJetPtArea[centbin]->Fill(jetPt,jet->Area());
     fHistJetPtNcon[centbin]->Fill(jetPt,jet->GetNumberOfConstituents());
     fHistJetPtNconCh[centbin]->Fill(jetPt,jet->GetNumberOfTracks());
     fHistJetPtNconEm[centbin]->Fill(jetPt,jet->GetNumberOfClusters());

     // make histo's with BIAS applied
     if (jet->MaxTrackPt()>fTrkBias){    
       fHistJetPtvsdEPBias[centbin]->Fill(jetPt, dEP);
       fHistJetEtaPhiPtBias[centbin]->Fill(jetPt,jet->Eta(),jet->Phi());
       fHistJetPtAreaBias[centbin]->Fill(jetPt,jet->Area());
       fHistJetPtNconBias[centbin]->Fill(jetPt,jet->GetNumberOfConstituents());
       fHistJetPtNconBiasCh[centbin]->Fill(jetPt,jet->GetNumberOfTracks());
       fHistJetPtNconBiasEm[centbin]->Fill(jetPt,jet->GetNumberOfClusters());
     }

    // fill some histo's
    fHistJetPt[centbin]->Fill(jet->Pt());  // fill #jets vs pT histo
    fHistJetEtaPhi->Fill(jet->Eta(),jet->Phi());     // fill jet eta-phi distribution histo
    //fHistDeltaPtvsArea->Fill(jetPt,jet->Area());

    if(leadjet && centbin==0) fHistJetPt[centbin+1]->Fill(jet->Pt());
    if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
      fHistJetPtBias[centbin]->Fill(jet->Pt());
      if(leadjet && centbin==0) fHistJetPtBias[centbin+1]->Fill(jet->Pt());
    }  // end of MaxTrackPt>ftrkBias or maxclusterPt>fclusBias

    // do we have trigger tracks
    if(iTT>0){
      AliVTrack* TT = static_cast<AliVTrack*>(tracks->At(iTT));
      if(TMath::Abs(jet->Phi()-TT->Phi()-TMath::Pi())<0.6) passedTTcut=1;
      else passedTTcut=0;
    } // end of check on iTT > 0

    if(passedTTcut)  fHistJetPtTT[centbin]->Fill(jet->Pt());

    // cut on HIGHEST jet pt in event (15 GeV default)
    if (highestjetpt>fJetPtcut) {

      // loop over all track for an event containing a jet with a pT>fJetPtCut  (15)GeV
      for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
        //AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
        // check on even type to get object for PID info
        //char eventtype[10];

	AliVTrack *track = fAOD->GetTrack(iTracks);
	if(!track) track = fESD->GetTrack(iTracks);
//     	AliAODTrack* track = fAOD->GetTrack(iTracks);
//         //sprintf(eventtype,"fAOD");
//         if(!track) AliESDtrack* track = fESD->GetTrack(iTracks);
//         //sprintf(eventtype,"fESD");
       
        if (!track) {
          AliError(Form("Couldn't get AliESDtrack %d\n", iTracks));
          continue;
        } 

//     AliESDtrack *track = GetAcceptTrack(esdtrack);
//     if(!track) continue;

	// apply track cuts
        if(TMath::Abs(track->Eta())>fTrkEta) 
          continue;
        //fHistTrackPt->Fill(track->Pt());  // fill track pT histogram
        if (track->Pt()<0.15)             // make sure track pT > 0.15 GeV
          continue;

        // calculate and get some track parameters
        track->GetTPCsignal();            // get the TPC signal of the track
        Double_t tracketa=track->Eta();   // eta of track
        Double_t deta=tracketa-jeteta;    // dETA between track and jet
        Int_t ieta=GetEtaBin(deta);       // bin of eta
	if(ieta<0) continue;

        // dPHI between jet and hadron
        Double_t dphijh = RelativePhi(jet->Phi(), track->Phi()); // angle between jet and hadron

        // jet pt bins
        //Int_t iptjet=-1;                  // initialize jet pT bin
//        iptjet=GetpTjetBin(jetPt);        // bin of jet pT
        Int_t iptjet=GetpTjetBin(jet->Pt());
	if(iptjet<0) continue;
        Double_t dR=sqrt(deta*deta+dphijh*dphijh);                   // difference of R between jet and hadron track

        // fill some jet-hadron histo's
        fHistJetH[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());  // fill jet-hadron dPHI--track pT distribution
        fHistJetHEtaPhi->Fill(deta,dphijh);                          // fill jet-hadron  eta--phi distribution
	fHistJetHaddPHI->Fill(dphijh);

        // does our max track or cluster pass the bias?
        if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
          // set up and fill Jet-Hadron THnSparse
          Double_t triggerEntries[9] = {fCent,jetPt,track->Pt(),dR,deta,dphijh,0,leadjet,zVtx};
          fhnJH->Fill(triggerEntries);    // fill Sparse Histo with trigger entries
          
          // fill histo's
          fHistSEphieta->Fill(dphijh, deta);
          fHistJetHBias[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
          fHistJetHaddPhiBias->Fill(dphijh);
	
          // in plane and out of plane histo's
          if( dEP>0 && dEP<=(TMath::Pi()/6) ){
             // we are IN plane
             fHistJetHaddPhiINBias->Fill(dphijh);
          }else if( dEP>(TMath::Pi()/3) && dEP<=(TMath::Pi()/2) ){
             // we are OUT of PLANE
             fHistJetHaddPhiOUTBias->Fill(dphijh);
          }else if( dEP>(TMath::Pi()/6) && dEP<=(TMath::Pi()/3) ){
             // we are in middle of plane
             fHistJetHaddPhiMIDBias->Fill(dphijh);
          }
        } // end of check maxtrackpt>ftrackbias or maxclusterpt>fclustbias

        if(passedTTcut) fHistJetHTT[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());

// **************************************************************************************************************
// *********************************** PID **********************************************************************
// **************************************************************************************************************
    if(doPID){
      // data information for pid    
      Double_t eta, pt, EovP, fClsE = -1.;

      EovP = 0;
      fClsE = 0;

      // nSigma of particles in TPC, TOF, and ITS
      Double_t nSigmaPion_TPC, nSigmaProton_TPC, nSigmaKaon_TPC = -1.;
      Double_t nSigmaPion_TOF, nSigmaProton_TOF, nSigmaKaon_TOF = -1.;
      Double_t nSigmaPion_ITS, nSigmaProton_ITS, nSigmaKaon_ITS = -1.;

      // misc quantities
      Double_t dEdx, ITSsig, TOFsig  = -1.;
      Double_t charge     = -99;
      //Double_t mom=track->GetTPCmomentum();

      // get parameters of track
      charge = track->Charge();    // charge of track
      eta    = track->Eta();       // ETA of track
      pt     = track->Pt();        // pT of track

      // extra attempt 
      AliVEvent *vevent=InputEvent();
      if (!vevent||!fPIDResponse) return kTRUE; // just return, maybe put at beginning

      // get impact parameters and detector signals
      dEdx = track->GetTPCsignal();
      ITSsig = track->GetITSsignal();
      TOFsig = track->GetTOFsignal();

      // fill detector signal histograms
      fHistTPCdEdX->Fill(pt, dEdx);
      fHistITSsignal->Fill(pt, ITSsig);
//     fHistTOFsignal->Fill(pt, TOFsig);

      // TPC nSigma's
      nSigmaPion_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion);
      nSigmaKaon_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon);
      nSigmaProton_TPC = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton);

      // TOF nSigma's
      nSigmaPion_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);
      nSigmaKaon_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kKaon);
      nSigmaProton_TOF = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kProton);

      // ITS nSigma's
      nSigmaPion_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kPion);
      nSigmaKaon_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kKaon);
      nSigmaProton_ITS = fPIDResponse->NumberOfSigmasITS(track,AliPID::kProton);

      // Tests to PID pions, kaons, and protons,          (default is undentified tracks)
      //Double_t nPIDtpc, nPIDits, nPIDtof = 0;
      Double_t	nPID = -10;

      // check track has pT < 0.900 GeV  - use TPC pid
      if (pt<0.900 && dEdx>0) {
	nPIDtpc = 4;

        // PION check - TPC
        if (TMath::Abs(nSigmaPion_TPC)<2 && TMath::Abs(nSigmaKaon_TPC)>2 && TMath::Abs(nSigmaProton_TPC)>2 ){
          isPItpc = kTRUE;
          nPIDtpc = 1;
          nPID=0.5;
        }else isPItpc = kFALSE; 

        // KAON check - TPC
        if (TMath::Abs(nSigmaKaon_TPC)<2 && TMath::Abs(nSigmaPion_TPC)>3 && TMath::Abs(nSigmaProton_TPC)>2 ){
          isKtpc = kTRUE;
          nPIDtpc = 2;
          nPID=1.5;
        }else isKtpc = kFALSE;

        // PROTON check - TPC
        if (TMath::Abs(nSigmaProton_TPC)<2 && TMath::Abs(nSigmaPion_TPC)>3 && TMath::Abs(nSigmaKaon_TPC)>2 ){
          isPtpc = kTRUE;
          nPIDtpc = 3;
          nPID=2.5;
        }else isPtpc = kFALSE;
      }  // cut on track pT for TPC

      // check track has pT < 0.500 GeV - use ITS pid
      if (pt<0.500 && ITSsig>0) {
        nPIDits = 4;

        // PION check - ITS
        if (TMath::Abs(nSigmaPion_ITS)<2 && TMath::Abs(nSigmaKaon_ITS)>2 && TMath::Abs(nSigmaProton_ITS)>2 ){
          isPIits = kTRUE;
          nPIDits = 1; 
	  nPID=4.5;
        }else isPIits = kFALSE;

        // KAON check - ITS
        if (TMath::Abs(nSigmaKaon_ITS)<2 && TMath::Abs(nSigmaPion_ITS)>3 && TMath::Abs(nSigmaProton_ITS)>2 ){
          isKits = kTRUE;
          nPIDits = 2;
          nPID=5.5;
        }else isKits = kFALSE;

        // PROTON check - ITS
        if (TMath::Abs(nSigmaProton_ITS)<2 && TMath::Abs(nSigmaPion_ITS)>3 && TMath::Abs(nSigmaKaon_ITS)>2 ){
          isPits = kTRUE;
          nPIDits = 3;
	  nPID=6.5;
        }else isPits = kFALSE;
      }  // cut on track pT for ITS

      // check track has 0.900 GeV < pT < 2.500 GeV - use TOF pid
      if (pt>0.900 && pt<2.500 && TOFsig>0) {
	nPIDtof = 4;

        // PION check - TOF
        if (TMath::Abs(nSigmaPion_TOF)<2 && TMath::Abs(nSigmaKaon_TOF)>2 && TMath::Abs(nSigmaProton_TOF)>2 ){
          isPItof = kTRUE;
          nPIDtof = 1;
          nPID=8.5;
        }else isPItof = kFALSE;

        // KAON check - TOF
        if (TMath::Abs(nSigmaKaon_TOF)<2 && TMath::Abs(nSigmaPion_TOF)>3 && TMath::Abs(nSigmaProton_TOF)>2 ){
          isKtof = kTRUE;
          nPIDtof = 2;
          nPID=9.5;
        }else isKtof = kFALSE;

        // PROTON check - TOF
        if (TMath::Abs(nSigmaProton_TOF)<2 && TMath::Abs(nSigmaPion_TOF)>3 && TMath::Abs(nSigmaKaon_TOF)>2 ){
          isPtof = kTRUE;
          nPIDtof = 3;
          nPID=10.5;
        }else isPtof = kFALSE;
      }  // cut on track pT for TOF

         if (nPIDtpc == 0) nPID = -0.5;
	 if (nPIDits == 0) nPID = 3.5;
         if (nPIDtof == 0) nPID = 7.5;
	 fHistPID->Fill(nPID);

          // PID sparse getting filled 
          Double_t nontrig_tracks_Entries[23] = {fCent,jetPt,pt,charge,eta,deta,dphijh,leadjet,zVtx,EovP,fClsE,
                                      nSigmaPion_TPC, nSigmaProton_TPC, nSigmaKaon_TPC,  //nSig in TPC
                                      nSigmaPion_ITS, nSigmaProton_ITS, nSigmaKaon_ITS,  //nSig in ITS
                                      nSigmaPion_TOF, nSigmaProton_TOF, nSigmaKaon_TOF,  //nSig in TOF
				      nPIDtpc, nPIDits, nPIDtof                          //PID label for each detector
                                      }; //array for PID sparse                           
          fhnPID->Fill(nontrig_tracks_Entries);   // fill Sparse histo of PID tracks 
        } // end of doPID check

	// get track pt bin
        Int_t itrackpt = -500;              // initialize track pT bin
        itrackpt = GetpTtrackBin(track->Pt());

	// all tracks: jet hadron correlations in hadron pt bins
        fHistJetHadbindPhi[itrackpt]->Fill(dphijh);

        // in plane and out of plane jet-hadron histo's
        if( dEP>0 && dEP<=(TMath::Pi()/6) ){
           // we are IN plane
           fHistJetHaddPhiINcent[centbin]->Fill(dphijh);
           fHistJetHaddPhiIN->Fill(dphijh);
           fHistJetHadbindPhiIN[itrackpt]->Fill(dphijh);
           //fHistJetHaddPhiPtcentbinIN[itrackpt][centbin]->Fill(dphijh);
        }else if( dEP>(TMath::Pi()/3) && dEP<=(TMath::Pi()/2) ){
           // we are OUT of PLANE
           fHistJetHaddPhiOUTcent[centbin]->Fill(dphijh);
           fHistJetHaddPhiOUT->Fill(dphijh);
           fHistJetHadbindPhiOUT[itrackpt]->Fill(dphijh);
           //fHistJetHaddPhiPtcentbinOUT[itrackpt][centbin]->Fill(dphijh);
        }else if( dEP>(TMath::Pi()/6) && dEP<=(TMath::Pi()/3) ){ 
           // we are in the middle of plane
           fHistJetHaddPhiMIDcent[centbin]->Fill(dphijh);
           fHistJetHaddPhiMID->Fill(dphijh);
           fHistJetHadbindPhiMID[itrackpt]->Fill(dphijh);
        }
      } // loop over tracks found in event with highest JET pT > 10.0 GeV (change)
    } // jet pt cut
  } // jet loop

// ***************************************************************************************************************
// ******************************** Event MIXING *****************************************************************
//  TObjArray* tracksClone = CloneAndReduceTrackList(tracks); // TEST

  //Prepare to do event mixing
  if(fDoEventMixing>0){
    // event mixing

    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    // mix jets from triggered events with tracks from MB events
    // get the trigger bit
    // need to change trigger bits between different runs
    UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if (trigger==0)  return kTRUE; // return

    //Double_t Ntrks=(Double_t)Ntracks*1.0;
    //cout<<"Test.. Ntrks: "<<fPoolMgr->GetEventPool(Ntrks);

    AliEventPool* pool = fPoolMgr->GetEventPool(fCent, zVtx); // for PbPb? fcent
    //AliEventPool* pool = fPoolMgr->GetEventPool(Ntrks, zVtx); // for pp

    if (!pool){
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCent, zVtx));
      //AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", Ntrks, zVtx));
      return kTRUE;
    }

    // use only jets from EMCal-triggered events (for lhc11a use AliVEvent::kEMC1)
///    if (trigger & AliVEvent::kEMC1) {
    if (trigger & AliVEvent::kEMCEJE) {  // TEST
      //check for a trigger jet
      // fmixingtrack/10 ??
      if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 /*100*/ || pool->GetCurrentNEvents() >= 5) {
        // loop over jets (passing cuts?)
        for (Int_t ijet = 0; ijet < Njets; ijet++) {
          //Double_t leadjet=0;
          //if (ijet==ijethi) leadjet=1;

          // get jet object
          AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
          if(!AcceptJet(jet))  continue;

          // get jet stats
          Double_t jetphi = jet->Phi();  // phi of jet
          Double_t jetPt = jet->Pt();    // pT of jet

          Int_t nMix = pool->GetCurrentNEvents();  // how many particles in pool to mix

          // Fill for biased jet triggers only
          if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)) {
            // Fill mixed-event histos here  
            for (Int_t jMix=0; jMix<nMix; jMix++) {
                TObjArray* bgTracks = pool->GetEvent(jMix);
                const Int_t Nbgtrks = bgTracks->GetEntries();
                for(Int_t ibg=0; ibg<Nbgtrks; ibg++) {
                  AliPicoTrack *part = static_cast<AliPicoTrack*>(bgTracks->At(ibg));
                  if(!part) continue;

                  //Double_t DPhi = jetphi - part->Phi();
                  Double_t DEta = part->Eta()-jet->Eta();                // difference in eta
                  Double_t DPhi = RelativePhi(jet->Phi(),part->Phi());   // difference in phi
                  Float_t dEP = RelativeEPJET(jetphi,fEPV0);
		  Double_t charge = part->Charge();

                  //Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);      // difference in R
                  if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
                  if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();
                           
                  Double_t triggerEntries[8] = {fCent,jetPt,part->Pt(),DEta,DPhi,dEP,zVtx, charge}; //array for ME sparse
                  fhnMixedEvents->Fill(triggerEntries,1./nMix);   // fill Sparse histo of mixed events

                  // *********************
                  //create phi distribution from mixed events
                  //fHistMEdPHI->Fill(DPhi);
                  fHistMEphieta->Fill(DPhi,DEta, 1./nMix);

                } // end of background track loop
             } // end of filling mixed-event histo's
          } // end of check for biased jet triggers
        } // end of jet loop
      } // end of check for triggered jet
    } //end EMC triggered loop

    // use only tracks from MB events (for lhc11a use AliVEvent::kMB)
///    if (trigger & AliVEvent::kMB) {
    if (trigger & AliVEvent::kAnyINT){ // test
      // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
      TObjArray* tracksClone = CloneAndReduceTrackList(tracks);

      // update pool if jet in event or not
      pool->UpdatePool(tracksClone);
    } // check on track from MB events
  } // end of event mixing

  // print some stats on the event
  event++;
  
  //cout<<"Event #: "<<event<<endl;
  //cout<<"# of jets: "<<Njets<<"      Highest jet pt: "<<highestjetpt<<endl;
  //cout<<"# tracks: "<<Ntracks<<"      Highest track pt: "<<ptmax<<endl;
  //cout<<" =============================================== "<<endl;

  return kTRUE;  // used when the function is of type bool
}  // end of RUN

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHadEPpid::GetCentBin(Double_t cent) const
{  // Get centrality bin.
  Int_t centbin = -1;
  if (cent>=0 && cent<10)	centbin = 0;
  else if (cent>=10 && cent<20)	centbin = 1;
  else if (cent>=20 && cent<30) centbin = 2;
  else if (cent>=30 && cent<40)	centbin = 3;
  else if (cent>=40 && cent<50) centbin = 4;
  else if (cent>=50 && cent<90)	centbin = 5;
 
  return centbin;
}

//________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetHadEPpid::RelativePhi(Double_t mphi,Double_t vphi) const
{ // function to calculate relative PHI
  double dphi = mphi-vphi;

  // set dphi to operate on adjusted scale
  if(dphi<-0.5*TMath::Pi()) dphi+=2.*TMath::Pi();
  if(dphi>3./2.*TMath::Pi()) dphi-=2.*TMath::Pi();

  // test
  if( dphi < -1.*TMath::Pi()/2 || dphi > 3.*TMath::Pi()/2 )
    AliWarning(Form("%s: dPHI not in range [-0.5*Pi, 1.5*Pi]!", GetName()));

  return dphi; // dphi in [-0.5Pi, 1.5Pi]                                                                                   
}



//_________________________________________________________________________
Float_t AliAnalysisTaskEmcalJetHadEPpid:: RelativeEPJET(Double_t jetAng, Double_t EPAng) const
{ // function to calculate angle between jet and EP in the 1st quadrant (0,Pi/2)
  Double_t dphi = (EPAng - jetAng);
  
  // ran into trouble with a few dEP<-Pi so trying this...
  if( dphi<-1*TMath::Pi() ){
    dphi = dphi + 1*TMath::Pi();
  } // this assumes we are doing full jets currently 
  
  if( (dphi>0) && (dphi<1*TMath::Pi()/2) ){
    // Do nothing! we are in quadrant 1
  }else if( (dphi>1*TMath::Pi()/2) && (dphi<1*TMath::Pi()) ){
    dphi = 1*TMath::Pi() - dphi;
  }else if( (dphi<0) && (dphi>-1*TMath::Pi()/2) ){
    dphi = fabs(dphi);
  }else if( (dphi<-1*TMath::Pi()/2) && (dphi>-1*TMath::Pi()) ){
    dphi = dphi + 1*TMath::Pi();
  } 
  
  // test
  if( dphi < 0 || dphi > TMath::Pi()/2 )
    AliWarning(Form("%s: dPHI not in range [0, 0.5*Pi]!", GetName()));

  return dphi;   // dphi in [0, Pi/2]
}

//Int_t ieta=GetEtaBin(deta);
//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHadEPpid::GetEtaBin(Double_t eta) const
{
  // Get eta bin for histos.
  Int_t etabin = -1;
  if (TMath::Abs(eta)<=0.4)				etabin = 0;
  else if (TMath::Abs(eta)>0.4 && TMath::Abs(eta)<0.8)	etabin = 1;
  else if (TMath::Abs(eta)>=0.8)			etabin = 2;

  return etabin;
} // end of get-eta-bin

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHadEPpid::GetpTjetBin(Double_t pt) const
{
  // Get jet pt  bin for histos.
  Int_t ptbin = -1;
  if (pt>=15 && pt<20)		ptbin = 0;
  else if (pt>=20 && pt<25)	ptbin = 1;
  else if (pt>=25 && pt<40)	ptbin = 2;
  else if (pt>=40 && pt<60)	ptbin = 3;
  else if (pt>=60)		ptbin = 4;

  return ptbin;
} // end of get-jet-pt-bin

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHadEPpid::GetpTtrackBin(Double_t pt) const
{
  // May need to update bins for future runs... (testing locally)

  // Get track pt bin for histos.
  Int_t ptbin = -1;
  if (pt < 0.5)			ptbin = 0;
  else if (pt>=0.5 && pt<1.0)	ptbin = 1;
  else if (pt>=1.0 && pt<1.5)	ptbin = 2;
  else if (pt>=1.5 && pt<2.0)	ptbin = 3;
  else if (pt>=2.0 && pt<2.5)	ptbin = 4;
  else if (pt>=2.5 && pt<3.0)	ptbin = 5;
  else if (pt>=3.0 && pt<4.0)	ptbin = 6;
  else if (pt>=4.0 && pt<5.0)	ptbin = 7;
  else if (pt>=5.0)		ptbin = 8;

  return ptbin;
} // end of get-jet-pt-bin


//___________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHadEPpid::GetzVertexBin(Double_t zVtx) const
{
  // get z-vertex bin for histo.
  int zVbin= -1;
  if (zVtx>=-10 && zVtx<-8)	zVbin = 0;
  else if (zVtx>=-8 && zVtx<-6)	zVbin = 1;
  else if (zVtx>=-6 && zVtx<-4)	zVbin = 2;
  else if (zVtx>=-4 && zVtx<-2)	zVbin = 3; 
  else if (zVtx>=-2 && zVtx<0)	zVbin = 4;
  else if (zVtx>=0 && zVtx<2)	zVbin = 5;
  else if (zVtx>=2 && zVtx<4)	zVbin = 6;
  else if (zVtx>=4 && zVtx<6)	zVbin = 7;
  else if (zVtx>=6 && zVtx<8)	zVbin = 8;
  else if (zVtx>=8 && zVtx<10)	zVbin = 9;
  else zVbin = 10;
        
  return zVbin;
} // end of get z-vertex bin

//______________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetHadEPpid::NewTHnSparseF(const char* name, UInt_t entries)
{
   // generate new THnSparseF, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){

         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF

void AliAnalysisTaskEmcalJetHadEPpid::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "V0 centrality (%)";
      nbins = 10;
      xmin = 0.;
      xmax = 100.;
      break;

   case 1:
      label = "Corrected Jet p_{T}";
      nbins = 250;
      xmin = -50.;
      xmax = 200.;
      break;

   case 2:
      label = "Track p_{T}";
      nbins = 100;
      xmin = 0.;
      xmax = 100.;
      break;

    case 3:
      label = "Relative Eta";
      nbins = 64;
      xmin = -1.6;
      xmax = 1.6;
      break;

   case 4: 
      label = "Relative Phi";
      nbins = 144;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

  case 5:
      label = "Relative angle of Jet and Reaction Plane";
      nbins = 144;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

  case 6:
      label = "z-vertex";
      nbins = 10;
      xmin = -10;
      xmax =  10;
      break;

  case 7:
      label = "track charge";
      nbins =3;
      xmin = -1.5;
      xmax = 1.5;
      break;

  case 8: // need to update
      label = "leading track";
      nbins = 10;
      xmin = 0;
      xmax = 50;
      break;

  case 9:  // need to update
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;
   } // end of switch
} // end of getting dim-params

//_________________________________________________
// From CF event mixing code PhiCorrelations
TObjArray* AliAnalysisTaskEmcalJetHadEPpid::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliPicoTrack which uses much less memory (used for event mixing)
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  /*
  Int_t nTrax = fESD->GetNumberOfTracks();
  //cout << "nTrax " << nTrax <<endl;
  
  for (Int_t i = 0; i < nTrax; ++i) 
    {
      AliESDtrack* esdtrack = fESD->GetTrack(i);
      if (!esdtrack) 
        {
          AliError(Form("Couldn't get ESD track %d\n", i));
          continue;
        }

      AliESDtrack *particle = GetAcceptTrack(esdtrack);
      if(!particle) continue;
  
      cout << "RM Hybrid track : " << i << "  " << particle->Pt() << endl;
  */

  //cout << "nEntries " << tracks->GetEntriesFast() <<endl;
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    if(TMath::Abs(particle->Eta())>fTrkEta) continue;
    if(particle->Pt()<0.15)continue;

    Double_t trackpt=particle->Pt();   // track pT

    Int_t trklabel=-1;
    trklabel=particle->GetLabel();
    //cout << "TRACK_LABEL: " << particle->GetLabel()<<endl;

    Int_t hadbin=-1;
    if(trackpt<0.5) hadbin=0;
    else if(trackpt<1) hadbin=1;
    else if(trackpt<2) hadbin=2;
    else if(trackpt<3) hadbin=3;
    else if(trackpt<5) hadbin=4;
    else if(trackpt<8) hadbin=5;
    else if(trackpt<20) hadbin=6;

//feb10 comment out
    if(hadbin>-1 && trklabel>-1 && trklabel <3) fHistTrackEtaPhi[trklabel][hadbin]->Fill(particle->Eta(),particle->Phi());
    if(hadbin>-1) fHistTrackEtaPhi[3][hadbin]->Fill(particle->Eta(),particle->Phi());

    //if(hadbin>-1) fHistTrackEtaPhi[hadbin]->Fill(particle->Eta(),particle->Phi());  // TEST

    tracksClone->Add(new AliPicoTrack(particle->Pt(), particle->Eta(), particle->Phi(), particle->Charge(), 0, 0, 0, 0));
  } // end of looping through tracks

  return tracksClone;
}

//____________________________________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetHadEPpid::NewTHnSparseFPID(const char* name, UInt_t entries)
{
   // generate new THnSparseF PID, axes are defined in GetDimParams()
   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){

         TString label("");
         GetDimParamsPID(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }

      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF PID

//________________________________________________________________________________
void AliAnalysisTaskEmcalJetHadEPpid::GetDimParamsPID(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse
   const Double_t pi = TMath::Pi();

   switch(iEntry){

   case 0:
      label = "V0 centrality (%)";
      nbins = 10;
      xmin = 0.;
      xmax = 100.;
      break;

   case 1:
      label = "Corrected Jet p_{T}";
      nbins = 250;
      xmin = -50.;
      xmax = 200.;
      break;

   case 2:
      label = "Track p_{T}";
      nbins = 500;
      xmin = 0.;
      xmax = 50.;
      break;

   case 3:
      label = "Charge of Track";
      nbins = 3;
      xmin = -1.5;
      xmax = 1.5;
      break;

   case 4:
      label = "Track Eta";
      nbins = 64;
      xmin = -1.2;
      xmax = 1.2;
      break;

   case 5:
      label = "Relative Eta of Track and Jet";
      nbins = 64;
      xmin = -1.6;
      xmax = 1.6;
      break;

   case 6:
      label = "Relative Phi of Track and Jet";
      nbins = 64;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;

   case 7:
      label = "leading jet";
      nbins = 3;
      xmin = -.5;
      xmax = 2.5;
      break;

   case 8:
      label = "Z-vertex";
      nbins = 10;
      xmin = -10.;
      xmax = 10.;
      break;

   case 9: // don't need
      label = "E/p of track";
      nbins = 100;
      xmin = 0.;
      xmax = 1.;
      break;

   case 10: // don't need
      label = "Cluster Energy";
      nbins = 50; //250;
      xmin = 0.;
      xmax = 50.; //250.;
      break;

   case 11:
      label = "N-Sigma of pions in TPC";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 12:
      label = "N-Sigma of protons in TPC";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 13:
      label = "N-Sigma of kaons in TPC";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 14:
      label = "N-Sigma of pions in ITS";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 15:
      label = "N-Sigma of protons in ITS";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 16:
      label = "N-Sigma of kaons in ITS";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 17:
      label = "N-Sigma of pions in TOF";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 18:
      label = "N-Sigma of protons in TOF";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 19:
      label = "N-Sigma of kaons in TOF";
      nbins = 100;
      xmin = -5.;
      xmax = 5.;
      break;

   case 20:
      label = "TPC PID determination";
      nbins = 4;
      xmin = 0.;
      xmax = 4.;
      break;

   case 21:
      label = "ITS PID determination";
      nbins = 4;
      xmin = 0.;
      xmax = 4.;
      break;

   case 22:
      label = "TOF PID determination";
      nbins = 4;
      xmin = 0.;
      xmax = 4.;
      break;

   } // end of switch
} // end of get dimension parameters PID

void AliAnalysisTaskEmcalJetHadEPpid::Terminate(Option_t *)
{
  cout<<"#########################"<<endl;
  cout<<"#### DONE RUNNING!!! ####"<<endl;
  cout<<"#########################"<<endl;

} // end of terminate


//void AliAnalysisTaskEmcalJetHadEPpid::FillAnalysisSummaryHistogram() const
void AliAnalysisTaskEmcalJetHadEPpid::SetfHistPIDcounterLabels(TH1* h) const
{
    // fill the analysis summary histrogram, saves all relevant analysis settigns
    h->GetXaxis()->SetBinLabel(1, "TPC: Unidentified"); // -0.5
    h->GetXaxis()->SetBinLabel(2, "TPC: Pion"); // 0.5
    h->GetXaxis()->SetBinLabel(3, "TPC: Kaon"); // 1.5
    h->GetXaxis()->SetBinLabel(4, "TPC: Proton"); // 2.5
    h->GetXaxis()->SetBinLabel(5, "ITS: Unidentified"); // 3.5
    h->GetXaxis()->SetBinLabel(6, "ITS: Pion"); // 4.5
    h->GetXaxis()->SetBinLabel(7, "ITS: Kaon"); // 5.5
    h->GetXaxis()->SetBinLabel(8, "ITS: Proton"); // 6.5
    h->GetXaxis()->SetBinLabel(9, "TOF: Unidentified"); // 7.5
    h->GetXaxis()->SetBinLabel(10, "TOF: Pion"); // 8.5 
    h->GetXaxis()->SetBinLabel(11, "TOF: Kaon"); // 9.5
    h->GetXaxis()->SetBinLabel(12, "TOF: Proton"); // 10.5

}
