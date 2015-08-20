
// Di-Jet angular correlations class
// Author: Greeshma Koyithatta Meethaleveedu and Ragahva Varma

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "THnSparse.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"
#include "AliEventPoolManager.h"

#include "AliAnalysisTaskDiJetCorrelations.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskDiJetCorrelations)

//_____________________| Constructor
AliAnalysisTaskDiJetCorrelations::AliAnalysisTaskDiJetCorrelations():
  AliAnalysisTaskSE(),
  fSetSystemValue(kTRUE),
  fRecoOrMontecarlo(kTRUE),
  fReadMC(kFALSE),
  fSetFilterBit(kTRUE),
  fbit(272),
  farrayMC(0),
  fCentrOrMult(1),
  fMinCentrality(0),
  fMaxCentrality(100),
  fTrigger1pTLowThr(0),
  fTrigger1pTHighThr(0),
  fTrigger2pTLowThr(0),
  fTrigger2pTHighThr(0),
  fAlpha(0),
  fCutResonances(kTRUE),
  fCutConversions(kTRUE),
  ftwoTrackEfficiencyCut(kTRUE),
  fuseVarCentBins(kFALSE),
  fuseVarPtBins(kFALSE),
  fHistNEvents(0),
  fHistT1CorrTrack(0),
  fHistT2CorrTrack(0),
  fOutputQA(0),
  fOutputCorr(0),
  f3DEffCor(0),
  fPool(0x0),
  fPoolMgr(0x0),
  fMixedEvent(kTRUE),
  fMEMaxPoolEvent(1000),
  fMEMinTracks(20000),
  fMEMinEventToMix(6),
  fHistTrigDPhi(0x0),
  fControlConvResT1(0x0),
  fControlConvResT2(0x0),
  fControlConvResMT1(0x0),
  fControlConvResMT2(0x0)
{
  for ( Int_t i = 0; i < 9; i++)fHistQA[i] = NULL;
}

//_____________________| Specific Constructor
AliAnalysisTaskDiJetCorrelations::AliAnalysisTaskDiJetCorrelations(const char *name):
  AliAnalysisTaskSE(name),
  fSetSystemValue(kTRUE),
  fRecoOrMontecarlo(kTRUE),
  fReadMC(kFALSE),
  fSetFilterBit(kTRUE),
  fbit(272),
  farrayMC(0),
  fCentrOrMult(1),
  fMinCentrality(0),
  fMaxCentrality(100),
  fTrigger1pTLowThr(0),
  fTrigger1pTHighThr(0),
  fTrigger2pTLowThr(0),
  fTrigger2pTHighThr(0),
  fAlpha(0),
  fCutResonances(kTRUE),
  fCutConversions(kTRUE),
  ftwoTrackEfficiencyCut(kTRUE),
  fuseVarCentBins(kFALSE),
  fuseVarPtBins(kFALSE),
  fHistNEvents(0),
  fHistT1CorrTrack(0),
  fHistT2CorrTrack(0),
  fOutputQA(0),
  fOutputCorr(0),
  f3DEffCor(0),
  fPool(0x0),
  fPoolMgr(0x0),
  fMixedEvent(kTRUE),
  fMEMaxPoolEvent(1000),
  fMEMinTracks(2000),
  fMEMinEventToMix(6),
  fHistTrigDPhi(0x0),
  fControlConvResT1(0x0),
  fControlConvResT2(0x0),
  fControlConvResMT1(0x0),
  fControlConvResMT2(0x0)

{
  Info("AliAnalysisTaskDiJetCorrelations","Calling Constructor");
  for ( Int_t i = 0; i < 9; i++)fHistQA[i] = NULL;
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class()); // Basic output slot (..)
  DefineOutput(2, TList::Class()); // Correlations form Data and MC    
}

//___________________________________| Copy Constructor
AliAnalysisTaskDiJetCorrelations::AliAnalysisTaskDiJetCorrelations(const AliAnalysisTaskDiJetCorrelations &source):
  AliAnalysisTaskSE(source),
  fSetSystemValue(source.fSetSystemValue),
  fRecoOrMontecarlo(source.fSetSystemValue),
  fReadMC(source.fReadMC),
  fSetFilterBit(source.fSetFilterBit),
  fbit(source.fbit),
  farrayMC(source.farrayMC),
  fCentrOrMult(source.fCentrOrMult),
  fMinCentrality(source.fMinCentrality),
  fMaxCentrality(source.fMaxCentrality),
  fTrigger1pTLowThr(source.fTrigger1pTLowThr),
  fTrigger1pTHighThr(source.fTrigger1pTHighThr),
  fTrigger2pTLowThr(source.fTrigger2pTLowThr),
  fTrigger2pTHighThr(source.fTrigger2pTHighThr),
  fAlpha(source.fAlpha),
  fCutResonances(source.fCutResonances),
  fCutConversions(source.fCutConversions),
  ftwoTrackEfficiencyCut(source.ftwoTrackEfficiencyCut),
  fuseVarCentBins(source.fuseVarCentBins),
  fuseVarPtBins(source.fuseVarPtBins),
  fHistNEvents(source.fHistNEvents),
  fHistT1CorrTrack(source.fHistT1CorrTrack),
  fHistT2CorrTrack(source.fHistT2CorrTrack),
  fOutputQA(source.fOutputQA),
  fOutputCorr(source.fOutputCorr),
  f3DEffCor(source.f3DEffCor),
  fPool(source.fPool),
  fPoolMgr(source.fPoolMgr),
  fMixedEvent(source.fMixedEvent),
  fMEMaxPoolEvent(source.fMEMaxPoolEvent),
  fMEMinTracks(source.fMEMinTracks),
  fMEMinEventToMix(source.fMEMinEventToMix),
  fHistTrigDPhi(source.fHistTrigDPhi),
  fControlConvResT1(source.fControlConvResT1),
  fControlConvResT2(source.fControlConvResT2),
  fControlConvResMT1(source.fControlConvResMT1),
  fControlConvResMT2(source.fControlConvResMT2)
{
  for ( Int_t i = 0; i < 9; i++)fHistQA[i] = NULL;
}

//_____________________| Destructor
AliAnalysisTaskDiJetCorrelations::~AliAnalysisTaskDiJetCorrelations()
{
  if(fOutputQA) {delete fOutputQA; fOutputQA = 0;}
  if(fOutputCorr) {delete fOutputCorr; fOutputCorr = 0;}
  if(fHistNEvents) {delete fHistNEvents; fHistNEvents = 0;}
}

//________________________________________|  Assignment Constructor
AliAnalysisTaskDiJetCorrelations& AliAnalysisTaskDiJetCorrelations::operator=(const AliAnalysisTaskDiJetCorrelations& orig)
{
  if (&orig == this) return *this; 
  AliAnalysisTaskSE::operator=(orig); 

  fSetSystemValue= orig.fSetSystemValue;
  fRecoOrMontecarlo = orig.fRecoOrMontecarlo;
  fReadMC = orig.fReadMC;
  fSetFilterBit = orig.fSetFilterBit;
  fbit  = orig.fbit;
  farrayMC = orig.farrayMC;
  fCentrOrMult = orig.fCentrOrMult;
  fMinCentrality = orig.fMinCentrality;
  fMaxCentrality = orig.fMaxCentrality;
  fTrigger1pTLowThr = orig.fTrigger1pTLowThr;
  fTrigger1pTHighThr = orig.fTrigger1pTHighThr;
  fTrigger2pTLowThr = orig.fTrigger2pTLowThr;
  fTrigger2pTHighThr = orig.fTrigger2pTHighThr;
  fAlpha= orig.fAlpha;
  fCutResonances = orig.fCutResonances;
  fCutConversions = orig.fCutConversions;
  ftwoTrackEfficiencyCut = orig.ftwoTrackEfficiencyCut;
  fuseVarCentBins = orig.fuseVarCentBins;
  fuseVarPtBins = orig.fuseVarPtBins;
  fHistNEvents = orig.fHistNEvents;
  fHistT1CorrTrack = orig.fHistT1CorrTrack;
  fHistT2CorrTrack = orig.fHistT2CorrTrack;
  fOutputQA = orig.fOutputQA;
  fOutputCorr = orig.fOutputCorr;
  f3DEffCor = orig.f3DEffCor;
  fPool = orig.fPool;
  fPoolMgr = orig.fPoolMgr;
  fMixedEvent = orig.fMixedEvent;
  fMEMaxPoolEvent = orig.fMEMaxPoolEvent;
  fMEMinTracks = orig.fMEMinTracks;
  fMEMinEventToMix = orig.fMEMinEventToMix;
  return *this; 
}

//_____________________| UserCreate Output
void AliAnalysisTaskDiJetCorrelations::UserCreateOutputObjects()
{
  fOutputQA = new TList();
  fOutputQA->SetOwner();
  fOutputQA->SetName("BasicQAHistograms");
  
  fOutputCorr = new TList();
  fOutputCorr->SetOwner();
  fOutputCorr->SetName("CorrelationsHistograms");

  fHistNEvents = new TH1F("fHistNEvents", "number of events ", 10, -0.5, 9.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEvents analyzed");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Rejected due Null Vtx and B");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"Within choosen centrality");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"With Good Vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"With Good SPDVertex");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"nTrack all chosen");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"nTrack Trigger1");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"nTrack Trigger2");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"nEvents with T1");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"number of DiJet");
  fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutputQA->Add(fHistNEvents);
  
  fHistT1CorrTrack = new TH1F("fHistT1CorrTrack", "T1 nCorr tracks ", 5, -0.5, 4.5);
  fHistT1CorrTrack->GetXaxis()->SetBinLabel(1,"nTrack Total");
  fHistT1CorrTrack->GetXaxis()->SetBinLabel(2,"after conversion");
  fHistT1CorrTrack->GetXaxis()->SetBinLabel(3,"after K0 resonance");
  fHistT1CorrTrack->GetXaxis()->SetBinLabel(4,"after Lambda resonance");
  fHistT1CorrTrack->GetXaxis()->SetBinLabel(5,"after 2track eff");
  fHistT1CorrTrack->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistT1CorrTrack->Sumw2();
  fHistT1CorrTrack->SetMinimum(0);
  fOutputQA->Add(fHistT1CorrTrack);
  
  fHistT2CorrTrack = new TH1F("fHistT2CorrTrack", "T2 nCorr tracks ", 5, -0.5, 4.5);
  fHistT2CorrTrack->GetXaxis()->SetBinLabel(1,"nTrack Total");
  fHistT2CorrTrack->GetXaxis()->SetBinLabel(2,"after conversion");
  fHistT2CorrTrack->GetXaxis()->SetBinLabel(3,"after K0 resonance");
  fHistT2CorrTrack->GetXaxis()->SetBinLabel(4,"after Lambda resonance");
  fHistT2CorrTrack->GetXaxis()->SetBinLabel(5,"after 2track eff");
  fHistT2CorrTrack->GetXaxis()->SetNdivisions(1,kFALSE);
  fHistT2CorrTrack->Sumw2();
  fHistT2CorrTrack->SetMinimum(0);
  fOutputQA->Add(fHistT2CorrTrack);
  
  DefineHistoNames();
  
  Bool_t DefPool = DefineMixedEventPool();
  if(!DefPool){
    AliInfo("UserCreateOutput: Pool is not define properly");
    return;
  }
  
  PostData(1,fOutputQA);
  PostData(2,fOutputCorr);
}


//_____________________| User Exec Part
void  AliAnalysisTaskDiJetCorrelations::UserExec(Option_t *) 
{
  
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!aod && AODEvent() && IsStandardAOD()) {
        aod = dynamic_cast<AliAODEvent*> (AODEvent());
    } else if(!aod)  {
        printf("AliAnalysisTaskDiJetCorrelations::UserExec: AOD not found!\n");
        return;
    }
  
  fHistNEvents->Fill(0);
    
  if(!fRecoOrMontecarlo){ // MC Kine
    farrayMC = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!farrayMC){
      AliError("Array of MC particles not found");
      return;
    }
    AliAODMCHeader *mcHeader = NULL;
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskDiJetCorrelations::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
  Float_t bSign = 0;
  bSign = (aod->GetMagneticField() > 0) ? 1 : -1;
  fHistNEvents->Fill(1);
  
  AliCentrality *centralityObj = 0x0;
  if(fSetSystemValue){ // pPb, PbPb
    centralityObj = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP();
    fCentrOrMult = centralityObj->GetCentralityPercentileUnchecked("V0M");
    if(centralityObj->GetQuality()!=0) return;
    if((abs(fCentrOrMult)) < 0. || (abs(fCentrOrMult)) > 100.1)return;
  }
  else if(!fSetSystemValue){ // pp, pPb
    Double_t count = -1, mineta = -1.0, maxeta = 1.0;
    AliAODTracklets* tracklets = aod->GetTracklets();
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>mineta && eta<maxeta) count++;
    }
    fCentrOrMult = count;
  }
  fHistNEvents->Fill(2); //
  
  AliAODVertex *vtxPrm = (AliAODVertex*)aod->GetPrimaryVertex();
  Bool_t isGoodVtx=kFALSE;
  TString primTitle = vtxPrm->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtxPrm->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fHistNEvents->Fill(3);
  }
  if(!isGoodVtx) return;
  
  Float_t zVertex = vtxPrm->GetZ();
  if (TMath::Abs(zVertex) > 10) return;
  fHistQA[0]->Fill(zVertex);
  
  const AliAODVertex* vtxSPD = aod->GetPrimaryVertexSPD();
  if(vtxSPD->GetNContributors()<=0) return;
  TString vtxTyp = vtxSPD->GetTitle();
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if(vtxTyp.Contains("vertexer:Z") && (zRes > 0.25)) return;
  if(TMath::Abs(vtxSPD->GetZ() - vtxPrm->GetZ()) > 0.5) return;
  fHistNEvents->Fill(4);
   
  TObjArray* fTrackArray = new TObjArray;
  fTrackArray->SetOwner(kTRUE);
  
  Double_t refmaxpT1 = -999, refmaxpT2 = -999;
  Double_t etaMaxpT1 = -999, etaMaxpT2 = -999;
  Double_t phiMaxpT1 = -999, phiMaxpT2 = -999;
  Short_t Charge1 = -999, Charge2 = -999;
  
  TString typeData = "";
  TString SEorME = "";
  
  if(fRecoOrMontecarlo){
    if(!fReadMC){
      typeData += "Data";//data
    }else
      typeData += "MCrc"; // MC reco
  }else{
    typeData += "MCKn"; // MC Generations  
  }
  
  if(!fMixedEvent){
    SEorME += "SE";
  }else if(fMixedEvent){
    SEorME += "ME";
  }
  
  //Mixed Events..

  if(fMixedEvent){
    if(TMath::Abs(zVertex)>=10){
          AliInfo(Form("Event with Zvertex = %0.2f cm out of pool bounds, SKIPPING",zVertex));
          return;
    }
      
    fPool= fPoolMgr->GetEventPool(fCentrOrMult, zVertex);
    if(!fPool){
      AliInfo(Form("No pool found for Event: multiplicity = %f, zVtx = %f cm", fCentrOrMult, zVertex));
      return;
    }
  }
  
  for (Int_t iTracks = 0; iTracks < aod->GetNumberOfTracks(); iTracks++){
    
    AliAODTrack* fAodTracks = (AliAODTrack*)aod->GetTrack(iTracks);
    if (!fAodTracks)continue;

    if(fSetFilterBit) if (!fAodTracks->TestFilterBit(fbit)) continue;
    if(fAodTracks->Eta() < -0.8 || fAodTracks->Eta() > 0.8)continue;
    if (fAodTracks->Pt() < 0.5 || fAodTracks->Pt() > 30.)continue;
    
      
    fHistNEvents->Fill(5); 
    fHistQA[1]->Fill(fAodTracks->GetTPCClusterInfo(2,1));
    fHistQA[3]->Fill(fAodTracks->DCA());
    fHistQA[4]->Fill(fAodTracks->ZAtDCA());
    fHistQA[5]->Fill(fAodTracks->Chi2perNDF());
    fHistQA[6]->Fill(fAodTracks->Pt());
    fHistQA[7]->Fill(fAodTracks->Phi());
    fHistQA[8]->Fill(fAodTracks->Eta());
    
    if(fRecoOrMontecarlo){ // reconstruction of data and MC
      if(fReadMC){
	// is track associated to particle ? if yes + implimenting the physical primary..
	Int_t label = TMath::Abs(fAodTracks->GetLabel());
	if (label<=0){
	  AliDebug(3,"Particle not matching MC label \n");
	  continue;
	}
	
	AliAODMCParticle *mcPart  = (AliAODMCParticle*)fMCEvent->GetTrack(label);
	if (!mcPart->IsPhysicalPrimary()) continue;
	fTrackArray->Add(fAodTracks);
      }else
	fTrackArray->Add(fAodTracks); //Storing all tracks for Data
    }
    
    if(fAodTracks->Pt() >= fTrigger1pTLowThr && fAodTracks->Pt()  <= fTrigger1pTHighThr){  
      if(fAodTracks->Pt() > refmaxpT1){
	refmaxpT1 = fAodTracks->Pt();
	etaMaxpT1 = fAodTracks->Eta();
	phiMaxpT1 = fAodTracks->Phi();
	Charge1 = fAodTracks->Charge();
      }
      fHistNEvents->Fill(6);//shouldn't the line be here?
    }
  }

  
  if(!fMixedEvent)if(refmaxpT1 < fTrigger1pTLowThr || refmaxpT1 > fTrigger1pTHighThr)return;
  fHistNEvents->Fill(8);
    
  for(Int_t entryT2=0; entryT2<fTrackArray->GetEntries(); entryT2++){
        
        TObject* obj = fTrackArray->At(entryT2);
        AliAODTrack* fAodTracksT2 = (AliAODTrack*)obj;
	
        if(fAodTracksT2->Pt() >= fTrigger2pTLowThr && fAodTracksT2->Pt()  <= fTrigger2pTHighThr){
	  
	  Double_t phiT2 = fAodTracksT2->Phi();
	  Double_t TrigDPhi = phiMaxpT1-phiT2;
          
	  Double_t TrigDPhi12 =  AssignCorrectPhiRange(TrigDPhi);
	  Double_t B2BTrigDPhi = TMath::Abs(TrigDPhi12-fAlpha);
            
	  if(B2BTrigDPhi > (TMath::Pi())/8)continue;
	  fHistTrigDPhi->Fill(B2BTrigDPhi);
	    
	  if(fAodTracksT2->Pt() > refmaxpT2 ){
	    fHistNEvents->Fill(7);
	    refmaxpT2 = fAodTracksT2->Pt();
	    etaMaxpT2 = fAodTracksT2->Eta();
	    phiMaxpT2 = fAodTracksT2->Phi();
	    Charge2 = fAodTracksT2->Charge();
	  }
        }
  }
    
    if(!fMixedEvent)if(refmaxpT2 < fTrigger2pTLowThr || refmaxpT2 > fTrigger2pTHighThr)return;
    fHistNEvents->Fill(9);
    
    if (refmaxpT1 != -999 && refmaxpT2 != -999){
      
      if(TMath::Abs(zVertex)>=10)AliInfo(Form("Event with Zvertex = %.2f cm out of pool bounds, SKIPPING",zVertex));
      
      Double_t fCentZvtxpT1[3] = {fCentrOrMult, zVertex, refmaxpT1};
      if(!fMixedEvent)((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg1CentZvtxpT_%s_%s",typeData.Data(), SEorME.Data())))->Fill(fCentZvtxpT1);
      
      Double_t fCentZvtxpT2[3] = {fCentrOrMult, zVertex, refmaxpT2};
      if(!fMixedEvent)((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg2CentZvtxpT_%s_%s",typeData.Data(), SEorME.Data())))->Fill(fCentZvtxpT2);
      
      Int_t NofEventsinPool = 1, NumberOfTracksStore=0; // SE
      if(fMixedEvent){ //ME
	     Bool_t PoolQuality = ProcessMixedEventPool();
	     if(!PoolQuality){
         AliInfo("Mixed event analysis: pool is not ready");
	     return;
	    }
	     NofEventsinPool = fPool->GetCurrentNEvents();
      }
      
      TObjArray* SEMEEvtTracks = (TObjArray*)fTrackArray->Clone("SE");
      for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){
	
      if(fMixedEvent)SEMEEvtTracks = fPool->GetEvent(jMix); // replacing from pool
      if(!SEMEEvtTracks)return;
	    
	  NumberOfTracksStore = SEMEEvtTracks->GetEntriesFast();
	//cout << "Number of Tracks in this event are = " << NumberOfTracksStore << endl;
	for(int k=0; k < NumberOfTracksStore; k++){
	  
	  TObject* objSEorME = SEMEEvtTracks->At(k);
	  AliAODTrack* fAodTracksAS = (AliAODTrack*)objSEorME;
	  if(!fAodTracksAS) continue;
	  if (fAodTracksAS->Pt() > refmaxpT1)continue;
                
	  if(fMixedEvent){
	    if(fAodTracksAS->Eta() < -0.8 || fAodTracksAS->Eta() > 0.8)continue;
	    if(fAodTracksAS->Pt() < 0.5 || fAodTracksAS->Pt() > 30)continue;
	    if(fSetFilterBit) if (!fAodTracksAS->TestFilterBit(fbit)) continue;
	  }
	  
	  if(fAodTracksAS->Pt() < refmaxpT1){
	    
	    Bool_t CutForConversionResonanceTrg1 = ConversionResonanceCut(refmaxpT1, phiMaxpT1, etaMaxpT1, Charge1, fAodTracksAS,fControlConvResT1, fHistT1CorrTrack);
	    if(fCutConversions || fCutResonances)if(!CutForConversionResonanceTrg1)continue;
            
	    Bool_t CutForTwoTrackEffiTrg1 = TwoTrackEfficiencyCut(refmaxpT1, phiMaxpT1, etaMaxpT1, Charge1, fAodTracksAS, bSign);
	    if(ftwoTrackEfficiencyCut)if(!CutForTwoTrackEffiTrg1)continue;
	    fHistT1CorrTrack->Fill(4);
            
	    Double_t deltaPhi1 = AssignCorrectPhiRange(phiMaxpT1 - fAodTracksAS->Phi());
	    Double_t deltaEta1  = etaMaxpT1 - fAodTracksAS->Eta();
            
	    Double_t ptTrack1 = fAodTracksAS->Pt();
	    Double_t ptLim_Sparse1 = 0.0;
	    ptLim_Sparse1 =((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg1CentZvtxDEtaDPhi_%s_%s",typeData.Data(), SEorME.Data())))->GetAxis(4)->GetXmax();
	    if(ptTrack1 > ptLim_Sparse1)ptTrack1 = ptLim_Sparse1-0.01; //filling all pT
            
	    Double_t CentzVtxDEtaDPhiTrg1[5] = {fCentrOrMult, zVertex, deltaEta1, deltaPhi1, ptTrack1};
	    Double_t effvalueT1 = 1.0;// GetTrackbyTrackEffValue(fAodTracksAS, fCentrOrMult, f3DEffCor);
	    ((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg1CentZvtxDEtaDPhi_%s_%s",typeData.Data(), SEorME.Data())))->Fill(CentzVtxDEtaDPhiTrg1, 1/effvalueT1);
                    
	    //Double_t fBasicsTrg1[5] = {fAodTracksAS->Pt(), fAodTracksAS->Eta(), fCentrOrMult, zVertex, fAodTracksAS->Phi()};
	    //if(!fMixedEvent)((THnSparseD*)fOutputQA->FindObject(Form("ThnTrg1BasicsPlots_%s_%s",typeData.Data(), SEorME.Data())))->Fill(fBasicsTrg1);           
	  }
          
          
	  if (fAodTracksAS->Pt() < refmaxpT2){
	    
	    Bool_t CutForConversionResonanceTrg2 = ConversionResonanceCut(refmaxpT2, phiMaxpT2, etaMaxpT2, Charge2, fAodTracksAS,fControlConvResT2, fHistT2CorrTrack);
	    if(fCutConversions || fCutResonances)if(!CutForConversionResonanceTrg2)continue;
            
	    Bool_t CutForTwoTrackEffiTrg2 = TwoTrackEfficiencyCut(refmaxpT2, phiMaxpT2, etaMaxpT2, Charge2, fAodTracksAS, bSign);
	    if(ftwoTrackEfficiencyCut)if(!CutForTwoTrackEffiTrg2)continue;
	    fHistT2CorrTrack->Fill(4);
                    
	    Double_t deltaPhi2 =  AssignCorrectPhiRange(phiMaxpT2 - fAodTracksAS->Phi());
	    Double_t deltaEta2 = etaMaxpT2 - fAodTracksAS->Eta();
            
	    Double_t ptLim_Sparse2 = 0.0;
	    Double_t ptTrack2 = fAodTracksAS->Pt();
	    ptLim_Sparse2 =((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg2CentZvtxDEtaDPhi_%s_%s",typeData.Data(), SEorME.Data())))->GetAxis(4)->GetXmax();
	    if(ptTrack2 > ptLim_Sparse2) ptTrack2 = ptLim_Sparse2-0.01; //filling all the pT in last bin
            
	    Double_t CentzVtxDEtaDPhiTrg2[5] = {fCentrOrMult, zVertex, deltaEta2, deltaPhi2,ptTrack2};
	    Double_t effvalueT2 =1.0;// GetTrackbyTrackEffValue(fAodTracksAS, fCentrOrMult,f3DEffCor);
	    
	    ((THnSparseD*)fOutputCorr->FindObject(Form("ThnTrg2CentZvtxDEtaDPhi_%s_%s",typeData.Data(), SEorME.Data())))->Fill(CentzVtxDEtaDPhiTrg2,1/effvalueT2);
            
	    //Double_t fBasicsTrg2[5] = {fAodTracksAS->Pt(), fAodTracksAS->Eta(), fCentrOrMult, zVertex, fAodTracksAS->Phi()};
	    //if(!fMixedEvent)((THnSparseD*)fOutputQA->FindObject(Form("ThnTrg2BasicsPlots_%s_%s",typeData.Data(), SEorME.Data())))->Fill(fBasicsTrg2);
            
	  }
	}
      }   
    }
    
    
    TObjArray* tracksClone = NULL;//
    tracksClone = (TObjArray*)fTrackArray->Clone();
    if(fMixedEvent && tracksClone->GetEntriesFast()>0)fPool->UpdatePool(tracksClone);
    
    PostData(1, fOutputQA);
    PostData(2, fOutputCorr);    
}

//_____________________|Terminate
void AliAnalysisTaskDiJetCorrelations::Terminate(Option_t *){
  
  fOutputQA = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputQA) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
  fOutputCorr = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputCorr) {
    printf("ERROR: fOutputCorr not available\n");
    return;
  }
  return;
}

//______________________________|  Cuts of Resonance and Conversions..
Bool_t AliAnalysisTaskDiJetCorrelations::ConversionResonanceCut(Double_t refmaxpT, Double_t phiMaxpT, Double_t etaMaxpT, Double_t Charge, AliAODTrack* AodTracks, TH2F*fControlConvResT, TH1F* fHistTCorrTrack){
  
  fHistTCorrTrack->Fill(0); //
  //Conversions
  if(fCutConversions && AodTracks->Charge() * Charge < 0){
    Float_t mass = GetInvMassSquaredCheap(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.510e-3, 0.510e-3); 
    if (mass < 0.1){
      mass = GetInvMassSquared(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.510e-3, 0.510e-3);
      fControlConvResT->Fill(0.0, mass);
      if(mass < 0.04 * 0.04)	return kFALSE;
    }
  }
  fHistTCorrTrack->Fill(1); //
  
  //K0s
  if (fCutResonances && AodTracks->Charge() * Charge < 0){   
    Float_t mass = GetInvMassSquaredCheap(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.1396, 0.1396);
    const Float_t kK0smass = 0.4976;
    if (TMath::Abs(mass -kK0smass * kK0smass)  < 0.1){
      mass = GetInvMassSquared(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.1396, 0.1396);
      fControlConvResT->Fill(1, mass -kK0smass * kK0smass);
      if(mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))return kFALSE;
    }
  }
  fHistTCorrTrack->Fill(2); //
  
  //lambda
  if (fCutResonances && AodTracks->Charge() * Charge < 0){
    Float_t mass1 = GetInvMassSquaredCheap(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.1396, 0.9383);
    Float_t mass2 = GetInvMassSquaredCheap(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.9383, 0.1396);
    const Float_t kLambdaMass = 1.115;
    if (TMath::Abs(mass1 -kLambdaMass * kLambdaMass)  < 0.1){
      mass1 = GetInvMassSquared(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.1396, 0.9383);
      fControlConvResT->Fill(2, mass1 - kLambdaMass * kLambdaMass);
      if(mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))return kFALSE;  
    }
    
    if (TMath::Abs(mass2 -kLambdaMass * kLambdaMass)  < 0.1){
      mass2 = GetInvMassSquared(refmaxpT, etaMaxpT, phiMaxpT, AodTracks->Pt(), AodTracks->Eta(), AodTracks->Phi(), 0.1396, 0.9383);
      fControlConvResT->Fill(2, mass2 - kLambdaMass * kLambdaMass);
      if(mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))return kFALSE;
    }
  }
  
  fHistTCorrTrack->Fill(3); //
  return kTRUE; 
}

//_____________________|  Two track Efficiency
Bool_t AliAnalysisTaskDiJetCorrelations::TwoTrackEfficiencyCut(Double_t refmaxpT, Double_t phiMaxpT, Double_t etaMaxpT, Double_t Charge, AliAODTrack* AodTracks, Float_t bSigntmp){
  Float_t pt1  = refmaxpT;
  Float_t phi1 = phiMaxpT;
  Float_t eta1 = etaMaxpT;
  Float_t charge1 = Charge;
  
  Float_t phi2 = AodTracks->Phi();
  Float_t eta2 = AodTracks->Eta();
  Float_t pt2  = AodTracks->Pt();
  Float_t charge2 = AodTracks->Charge();
  
  Float_t deta = eta1 - eta2;
  
  if (TMath::Abs(deta) < 0.02 * 2.5 * 3){
    Float_t dphistar1 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 0.8, bSigntmp);
    Float_t dphistar2 = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, 2.5, bSigntmp);
    
    const Float_t kLimit = 0.02*3;
    Float_t dphistarminabs = 1e5;
    //Float_t dphistarmin = 1e5;
  
    if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0){
      for (Double_t rad=0.8; rad<2.51; rad+=0.01){
	Float_t dphistar = GetDPhiStar(phi1, pt1, charge1, phi2, pt2, charge2, rad, bSigntmp);
	Float_t dphistarabs = TMath::Abs(dphistar);
	if (dphistarabs < dphistarminabs){
	  //dphistarmin = dphistar;
	  dphistarminabs = dphistarabs;
	}
      }
       
      if (dphistarminabs < 0.02 && TMath::Abs(deta) < 0.02){
	return kFALSE;
      }
    }
  }
  return kTRUE;
}


//______________________________|  Nomenclature of Histograms
void AliAnalysisTaskDiJetCorrelations::DefineHistoNames(){
  
  Double_t Pi = TMath::Pi();
  //QA histograms
  fHistQA[0] = new TH1F("fHistZVtx", "Z vertex distribution", 100, -15., 15.);
  fHistQA[1] = new TH1F("fHistnTPCCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTPCClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("fHistDCAXY", "dca-XY", 100, -5., 5.);
  fHistQA[4] = new TH1F("fHistDCAZ", "dca-Z", 100, -5., 5.);
  fHistQA[5] = new TH1F("fHistChi2TPC", "Chi2 TPC", 100, 0., 10.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",1000,0.,20.);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 100, -0.5, 2*Pi+0.5);
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 100, -2, 2);
  
  for( Int_t i = 0; i < 9; i++)
    {
      fHistQA[i]->Sumw2();
      fOutputQA->Add(fHistQA[i]);
    }
  
  //1 SE Distributuons Phi, Eta for Trg1 and Trg2
  fHistTrigDPhi = new TH1F("fHistTrigDPhi", " Trig Phi Difference Same",100, 0, 0.5);
  fHistTrigDPhi->Sumw2();
  fOutputQA->Add(fHistTrigDPhi);
  
  fControlConvResT1 = new TH2F("fControlConvResT1", ";id;delta mass;T1", 3, -0.5, 2.5, 100, -0.1, 0.1);
  fControlConvResT2 = new TH2F("fControlConvResT2", ";id;delta mass;T2", 3, -0.5, 2.5, 100, -0.1, 0.1);
  fControlConvResT1->Sumw2();
  fControlConvResT2->Sumw2();
  fOutputQA->Add(fControlConvResT1);
  fOutputQA->Add(fControlConvResT2);
  
  //Thnsprase Distributuons Phi, Eta for Trg1 and Trg2
  TString nameThnTrg1CentZvtxDEtaDPhi   = "ThnTrg1CentZvtxDEtaDPhi";
  TString nameThnTrg2CentZvtxDEtaDPhi   = "ThnTrg2CentZvtxDEtaDPhi";
  TString nameThnTrg1CentZvtxpT     = "ThnTrg1CentZvtxpT" ;
  TString nameThnTrg2CentZvtxpT    = "ThnTrg2CentZvtxpT" ;
  //TString nameThnTrg1BasicsPlots    = "ThnTrg1BasicsPlots";
  //TString nameThnTrg2BasicsPlots    = "ThnTrg2BasicsPlots";
  
  if(fRecoOrMontecarlo){
    if(!fReadMC){ //data
      nameThnTrg1CentZvtxDEtaDPhi  += "_Data";
      nameThnTrg2CentZvtxDEtaDPhi  += "_Data";
      nameThnTrg1CentZvtxpT  += "_Data";
      nameThnTrg2CentZvtxpT  += "_Data";
      //nameThnTrg1BasicsPlots  += "_Data";
      //nameThnTrg2BasicsPlots  += "_Data";
    }
    else {// MC reco
      nameThnTrg1CentZvtxDEtaDPhi  += "_MCrc";
      nameThnTrg2CentZvtxDEtaDPhi  += "_MCrc";
      nameThnTrg1CentZvtxpT  += "_MCrc";
      nameThnTrg2CentZvtxpT  += "_MCrc";
      //nameThnTrg1BasicsPlots  += "_MCrc";
      //nameThnTrg2BasicsPlots  += "_MCrc";
    }
  }
  else{ // MC Generations
    nameThnTrg1CentZvtxDEtaDPhi  += "_MCKn";
    nameThnTrg2CentZvtxDEtaDPhi  += "_MCKn";
    nameThnTrg1CentZvtxpT  += "_MCKn";
    nameThnTrg2CentZvtxpT  += "_MCKn";
    //nameThnTrg1BasicsPlots  += "_MCKn";
    //nameThnTrg2BasicsPlots  += "_MCKn";
  }
  
  if(!fMixedEvent){
    nameThnTrg1CentZvtxDEtaDPhi  += "_SE";
    nameThnTrg2CentZvtxDEtaDPhi  += "_SE";
    nameThnTrg1CentZvtxpT  += "_SE";
    nameThnTrg2CentZvtxpT  += "_SE";
    //nameThnTrg1BasicsPlots  += "_SE";
    //nameThnTrg2BasicsPlots  += "_SE";
  }else if(fMixedEvent){
    nameThnTrg1CentZvtxDEtaDPhi  += "_ME";
    nameThnTrg2CentZvtxDEtaDPhi  += "_ME";
    nameThnTrg1CentZvtxpT  += "_ME";
    nameThnTrg2CentZvtxpT  += "_ME";
    }

  //Catgry :1 Trigger Particles --> Cent, Zvtx, Trigger_pT
  //_____________________________________________Trigger-1
  const Int_t pTbinTrigger1 = Int_t(fTrigger1pTHighThr - fTrigger1pTLowThr);
  Int_t   fBinsTrg1[3]   = {12,       10,        pTbinTrigger1};
  Double_t fMinTrg1[3]   = {0.0,   -10.0,    fTrigger1pTLowThr};
  Double_t fMaxTrg1[3]   = {100.0,  10.0,   fTrigger1pTHighThr};
  THnSparseD *THnTrig1CentZvtxpT = new THnSparseD(nameThnTrg1CentZvtxpT.Data(),"Cent-Zvtx-pTtr1",3, fBinsTrg1, fMinTrg1, fMaxTrg1);
  THnTrig1CentZvtxpT->Sumw2();
  fOutputCorr->Add(THnTrig1CentZvtxpT);

  //_____________________________________________Trigger-2
  const Int_t pTbinTrigger2 = Int_t(fTrigger2pTHighThr - fTrigger2pTLowThr);
  Int_t   fBinsTrg2[3]   = {12,       10,        pTbinTrigger2};
  Double_t fMinTrg2[3]   = {0.0,   -10.0,    fTrigger2pTLowThr};
  Double_t fMaxTrg2[3]   = {100.0,  10.0,   fTrigger2pTHighThr};
  THnSparseD *THnTrig2CentZvtxpT = new THnSparseD(nameThnTrg2CentZvtxpT.Data(),"Cent-Zvtx-pTtr2",3, fBinsTrg2, fMinTrg2, fMaxTrg2);
  THnTrig2CentZvtxpT->Sumw2();
  fOutputCorr->Add(THnTrig2CentZvtxpT);
  
    
  //Catgry2: Correaltions Plots for SE and ME (T1, T2)
  Int_t    fBins12[5] = {12,     10,   16,               36,    11};
  Double_t  fMin12[5] = {0.0,   -10., -1.6, -0.5*TMath::Pi(),  0.5};
  Double_t  fMax12[5] = {100.0,  10.,  1.6,  1.5*TMath::Pi(), 6.00};
  THnSparseD *THnTrig1CentZvtxDEtaDPhi   = new THnSparseD(nameThnTrg1CentZvtxDEtaDPhi.Data(),"Cent-zVtx-DEta1-DPhi1-Trk",5, fBins12, fMin12, fMax12);
  THnSparseD *THnTrig2CentZvtxDEtaDPhi   = new THnSparseD(nameThnTrg2CentZvtxDEtaDPhi.Data(),"Cent-zVtx-DEta2-DPhi2-Trk",5, fBins12, fMin12, fMax12);
  if(fuseVarCentBins){
    const Int_t nvarBinsCent = 12;
    Double_t varBinsCent[nvarBinsCent+1] = {0., 1., 2., 3., 4., 5., 7.5, 10., 20., 30., 40., 50., 100.1};
    THnTrig1CentZvtxDEtaDPhi->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
    THnTrig2CentZvtxDEtaDPhi->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
    THnTrig1CentZvtxpT->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
    THnTrig2CentZvtxpT->GetAxis(0)->Set(nvarBinsCent, varBinsCent);
  }

  //Munual pT tracks Values
  if(fuseVarPtBins){
    const Int_t nvarBinspT = 11;
    Double_t varBinspT[nvarBinspT+1] = {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
    THnTrig1CentZvtxDEtaDPhi->GetAxis(4)->Set(nvarBinspT, varBinspT);
    THnTrig2CentZvtxDEtaDPhi->GetAxis(4)->Set(nvarBinspT, varBinspT);
    }
    
  THnTrig1CentZvtxDEtaDPhi->Sumw2();
  THnTrig2CentZvtxDEtaDPhi->Sumw2();
  fOutputCorr->Add(THnTrig1CentZvtxDEtaDPhi);
  fOutputCorr->Add(THnTrig2CentZvtxDEtaDPhi);
 /*
  
  Int_t   fBinsTrk[5]   = {  20,   16,    10,    10,  36};
  Double_t fMinTrk[5]   = { 0.0, -0.8,   0.0,  -10.,  0};
  Double_t fMaxTrk[5]   = {20.0,  0.8, 100.1,   10.,  2*TMath::Pi()};
  THnSparseD *THnTrg1BasicsPlots = new THnSparseD(nameThnTrg1BasicsPlots.Data(),"pTtr1-Eta-Cent-zvtx-Phi",5, fBinsTrk, fMinTrk, fMaxTrk);
  THnSparseD *THnTrg2BasicsPlots = new THnSparseD(nameThnTrg2BasicsPlots.Data(),"pTtr2-Eta-Cent-zvtx-Phi",5, fBinsTrk, fMinTrk, fMaxTrk);
  //THnTrg1BasicsPlots->Sumw2();
  //THnTrg2BasicsPlots->Sumw2();
  //fOutputQA->Add(THnTrg1BasicsPlots);
  //fOutputQA->Add(THnTrg2BasicsPlots);*/
}


//______________________________|  Track Efficiency
Double_t AliAnalysisTaskDiJetCorrelations::GetTrackbyTrackEffValue(AliAODTrack* track, Double_t CentrOrMult, TH3F* h){
  Float_t effvalue = -999.;
        
  //if(!h)cout << "Histograms was not there" << endl;
  Int_t binX = h->GetXaxis()->FindBin(CentrOrMult);
  Int_t binY = h->GetYaxis()->FindBin(track->Eta());
  Int_t binZ = h->GetZaxis()->FindBin(track->Pt());
  effvalue = h->GetBinContent(binX, binY, binZ);
  cout <<"effi value = "<< effvalue << endl; exit(1);
  if(effvalue==0) effvalue=1.;
  return effvalue;
}


//end of file : Greeshma

