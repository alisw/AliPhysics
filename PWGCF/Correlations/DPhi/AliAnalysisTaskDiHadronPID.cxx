// ----------------------------------------------------------------------------
// This class makes di-hadron correlations, with TOF and TPC signals for
// the associated particles.
//
//   Last Update:
//     - Added a number of setters/getters.
//     - Added spectra in pt, eta and phi as a function of performed cuts.
//     - Added DCA histogram.
//     - Added pp functionality.
//     - Added the option to make a DCA cut and an ITS cut.
//     - Variable centrality.
//     - Variable maximum p_T for triggers.
//     - Removed di-hadron correlations with one PID signal.
//
// ----------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)
// Last edit: Apr 12th 2012. (v 8.00)
// ----------------------------------------------------------------------------

#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"

#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODVertex.h"
//#include "AliAODPid.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
//#include "AliTOFPIDResponse.h"

using namespace std;

#include "AliAnalysisTaskDiHadronPID.h"

ClassImp(AliAnalysisTaskDiHadronPID);

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID():
	AliAnalysisTaskSE(),
	fPIDResponse(0x0),
	fAODEvent(0x0),
	fAODHeader(0x0),
	fAODVertex(0x0),
	fAODTrack(0x0),
	fGlobalTracks(0x0),
	fCentrality(0x0),
	fVertexZ(0x0),
    fDCA(0x0),
    fDCAZoomed(0x0),
    fDCAZoomedTwice(0x0),
    fDCACut(0x0),
    fDCAZoomedCut(0x0),
    fDCAZoomedTwiceCut(0x0),
    fITSHits(0x0),
    fTrackCutsCount(0x0),
    fTrackCutsPt(0x0),
    fTrackCutsEta(0x0),
    fTrackCutsPhi(0x0),
    fEtaSpectrumTrig(0x0),
    fEtaSpectrumAssoc(0x0),
    fPhiSpectrumAssoc(0x0),
	fTPCnSigmaProton(0x0),
	fTOFnSigmaProton(0x0),
	fTPCnSigmaPion(0x0),
	fTOFnSigmaPion(0x0),
	fTPCnSigmaKaon(0x0),
	fTOFnSigmaKaon(0x0),
	fTPCSignal(0x0),
	fTOFSignal(0x0),
    fMixedEvents(0x0),
	fHistoList(0x0),
    fCalculateMixedEvents(kFALSE),
    fBeamType("PbPb"),
    fMaxEta(0.8),
    fMaxPlotEta(0.9),
    fMaxPt(10.),
    fNEtaBins(25),
    fNPhiBins(36),
    fVertexZMixedEvents(2.),
    fCentralityCutMax(0.),
    fCentralityCutMin(10.),
    fZoomed(kFALSE),
    fDoITSCut(kFALSE),
    fDoDCACut(kFALSE),
    fDemandNoMismatch(kFALSE),
    fVerbose(0),
    fPrintBufferSize(kFALSE),
	fTrigBufferIndex(0),	
	fTrigBufferSize(0),
	fTrigBufferMaxSize(1000)

{
    
	//
	// Default Constructor.
	//
	
	// Trigger buffer.
	for(Int_t ii=0; ii<25000; ii++) {
		for(Int_t jj=0; jj<4; jj++) {
			fTrigBuffer[ii][jj]=0;
		}				
	}	
	
    // The identified di-hadron correlations.
	for (Int_t ii = 0; ii < 3; ii++) {
		for (Int_t jj = 0; jj < 10; jj++) {
			fDiHadronTPCTOF[ii][jj]=0x0;
		}
	}
    
    // Track cut labels.
    for (Int_t ii=0; ii<8; ii++) {
        fTrackCutLabelNumbers[ii]=ii+1;
    }
}

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID(const char *name):
	AliAnalysisTaskSE(name),
    fPIDResponse(0x0),
    fAODEvent(0x0),
    fAODHeader(0x0),
    fAODVertex(0x0),
    fAODTrack(0x0),
    fGlobalTracks(0x0),
    fCentrality(0x0),
    fVertexZ(0x0),
    fDCA(0x0),
    fDCAZoomed(0x0),
    fDCAZoomedTwice(0x0),
    fDCACut(0x0),
    fDCAZoomedCut(0x0),
    fDCAZoomedTwiceCut(0x0),
    fITSHits(0x0),
    fTrackCutsCount(0x0),
    fTrackCutsPt(0x0),
    fTrackCutsEta(0x0),
    fTrackCutsPhi(0x0),
    fEtaSpectrumTrig(0x0),
    fEtaSpectrumAssoc(0x0),
    fPhiSpectrumAssoc(0x0),
    fTPCnSigmaProton(0x0),
    fTOFnSigmaProton(0x0),
    fTPCnSigmaPion(0x0),
    fTOFnSigmaPion(0x0),
    fTPCnSigmaKaon(0x0),
    fTOFnSigmaKaon(0x0),
    fTPCSignal(0x0),
    fTOFSignal(0x0),
    fMixedEvents(0x0),
    fHistoList(0x0),
    fCalculateMixedEvents(kFALSE),
    fBeamType("PbPb"),
    fMaxEta(0.8),
    fMaxPlotEta(0.9),
    fMaxPt(10.),
    fNEtaBins(25),
    fNPhiBins(36),
    fVertexZMixedEvents(2.),
    fCentralityCutMax(0.),
    fCentralityCutMin(10.),
    fZoomed(kFALSE),
    fDoITSCut(kFALSE),
    fDoDCACut(kFALSE),
    fDemandNoMismatch(kFALSE),
    fVerbose(0),
    fPrintBufferSize(kFALSE),
    fTrigBufferIndex(0),	
    fTrigBufferSize(0),
	fTrigBufferMaxSize(1000)

{
    
	//
	// Named Constructor.
	//
	
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
	
	// Trigger buffer.
	for(Int_t ii=0; ii<25000; ii++) {
		for(Int_t jj=0; jj<4; jj++) {
			fTrigBuffer[ii][jj]=0;
		}				
	}	
	
    // The identified di-hadron correlations.
	for (Int_t ii = 0; ii < 3; ii++) {
		for (Int_t jj = 0; jj < 10; jj++) {
			fDiHadronTPCTOF[ii][jj]=0x0;
		}
	}
    
    // Track cut labels.
    for (Int_t ii=0; ii<8; ii++) {
        fTrackCutLabelNumbers[ii]=ii+1;
    }
}

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::~AliAnalysisTaskDiHadronPID() {

	//
	// Destructor.
	//
	
    if(fGlobalTracks) {
        delete fGlobalTracks;
        fGlobalTracks=0;
    }

}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::UserCreateOutputObjects() 

{
	//
	// Output objects.
	//
    
    // Print the settings of the analysis task.
    cout<<endl;
    cout<<"-----------------------------------------------"<<endl;
    cout<<" AliAnalysisTaskDiHadronPID Settings:"<<endl;
    cout<<"-----------------------------------------------"<<endl;
    cout<<"Verbose Level: "<<fVerbose<<endl;
    cout<<"Mixed Events Calculated: "<<fCalculateMixedEvents<<endl;
    cout<<"Beam Type: "<<fBeamType<<endl;
    cout<<Form("Max eta: %3.1f",fMaxEta)<<endl;
    cout<<Form("Max eta plotted: %3.1f",fMaxPlotEta)<<endl;
    cout<<Form("Max p_T for the triggers: %3.1f GeV/c.",fMaxPt)<<endl;
    cout<<"Nbins in eta and delta eta: "<<fNEtaBins<<endl;
    cout<<"Nbins in phi and delta phi: "<<fNPhiBins<<endl;
    cout<<Form("Events mixed if vertices differ %3.1f cm.",fVertexZMixedEvents)<<endl;
    cout<<Form("Centrality between %3.1f and %3.1f percent.",fCentralityCutMax,fCentralityCutMin)<<endl;
    cout<<"Tracks are cut if less than 2 SPD hits: "<<fDoITSCut<<endl;
    cout<<"Tracks cut if DCA is too big: "<<fDoDCACut<<endl;
    cout<<"Tracks cut if there is a TPC-TOF mismatch: "<<fDemandNoMismatch<<endl;
    cout<<"Maximum number of triggers stored: "<<fTrigBufferMaxSize<<endl;
    cout<<"-----------------------------------------------"<<endl;
    cout<<endl;
    
	// Obtain a pointer to the analysis manager.
	AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
    if (!manager) {
        if (fVerbose>0) cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> ERROR: Analysis manager not found."<<endl;
        return;
    }
    if (fVerbose>1) {
        cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> Analysis manager found."<<endl;
    }
    
    // Obtain a pointer to the input handler.
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());
    if (!inputHandler) {
        if (fVerbose>0) cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> ERROR: Input handler not found."<<endl;
        return;
    }
    if (fVerbose>1) {
        cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> Input handler found."<<endl;
    }
    
	// Obtain a pointer to the PID response object.	
	fPIDResponse = inputHandler->GetPIDResponse();
    if (!fPIDResponse) {
        if (fVerbose>0) cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> ERROR: PID response object not found."<<endl;
        return;
    }
    if (fVerbose>1) {
        cout<<"AliAnalysisTaskDiHadronPID::UserCreateOutputObjects -> PID response object found."<<endl;
    }

	// Create the output of the task.	
	fHistoList = new TList();
	fHistoList->SetOwner(kTRUE); 
	
	// Ranges in dPhi, dEta, and the number of bins for di-hadron correlations.
	Int_t binsHisto[4] = {fNPhiBins,fNEtaBins};
	Double_t minHisto[4] = {-TMath::Pi()/2.,-2*fMaxPlotEta};
	Double_t maxHisto[4] = {3.*TMath::Pi()/2,2*fMaxPlotEta};
	
    // --- EVENT QA PLOTS ---
    
	fCentrality = new TH1F("fCentrality","Centrality;Centrality;N_{evt}",100,0,100);
	fHistoList->Add(fCentrality);
	fVertexZ = new TH1F("fVertexZ","Vertex Z position;z (cm);N_{evt}",30,-15,15);
	fHistoList->Add(fVertexZ);
	
    // --- TRACK QA PLOTS ---
    
    fDCA = new TH2F("fDCA","DCA positions TPC only cuts;xy (cm);z (cm)",100,-5,5,100,-5,5);
    fHistoList->Add(fDCA);
    fDCAZoomed = new TH2F("fDCAZoomed","DCA positions TPC only cuts;xy (cm);z (cm)",100,-.5,.5,100,-.5,.5);
    fHistoList->Add(fDCAZoomed);
    fDCAZoomedTwice = new TH2F("fDCAZoomedTwice","DCA positions TPC only cuts;xy (cm);z (cm)",100,-.05,.05,100,-.05,.05);
    fHistoList->Add(fDCAZoomedTwice);
    fDCACut = new TH2F("fDCACut","DCA positions after DCA cut;xy (cm);z (cm)",100,-5,5,100,-5,5);
    fHistoList->Add(fDCACut);
    fDCAZoomedCut = new TH2F("fDCAZoomedCut","DCA positions after DCA cut;xy (cm);z (cm)",100,-.5,.5,100,-.5,.5);
    fHistoList->Add(fDCAZoomedCut);
    fDCAZoomedTwiceCut = new TH2F("fDCAZoomedTwiceCut","DCA positions after DCA cut;xy (cm);z (cm)",100,-.05,.05,100,-.05,.05);
    fHistoList->Add(fDCAZoomedTwiceCut);
    
    fITSHits = new TH1F("fITSHits","ITS hits",3,0,3);
    (fITSHits->GetXaxis())->SetBinLabel(1,"No SPD hit");
    (fITSHits->GetXaxis())->SetBinLabel(2,"One SPD hit");
    (fITSHits->GetXaxis())->SetBinLabel(3,"Two SPD hits");
    fHistoList->Add(fITSHits);
    
    // Determine the bins used in the track cut histograms.
    const Int_t ncuts = 5 + fDoITSCut + fDoDCACut + fDemandNoMismatch;
    
    if (!fDoITSCut) {
        fTrackCutLabelNumbers[3]--;
        fTrackCutLabelNumbers[4]--;
        fTrackCutLabelNumbers[5]--;
        fTrackCutLabelNumbers[6]--;
        fTrackCutLabelNumbers[7]--;
    }
    
    if (!fDoDCACut) {
        fTrackCutLabelNumbers[4]--;
        fTrackCutLabelNumbers[5]--;
        fTrackCutLabelNumbers[6]--;
        fTrackCutLabelNumbers[7]--;
    }
    
    fTrackCutsCount = new TH1F("fTrackCutsCount","Track Cuts;;Count",ncuts,0,ncuts);
    (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[0],"Std. TPC Only");
    (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[1],Form("#eta < %3.1f",fMaxEta));
    if (fDoITSCut) (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[2],"1-2 SPD hits");
    if (fDoDCACut) (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[3],"DCA cut");
    (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[4],"p_{T} < 5.0 GeV/c");
    (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[5],"TPCpid (p_{T} < 5.0 GeV/c)");
    (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[6],"TOFpid (p_{T} < 5.0 GeV/c)");
    if (fDemandNoMismatch) (fTrackCutsCount->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[7],"No Mismatch (p_{T} < 5.0 GeV/c)");
    fHistoList->Add(fTrackCutsCount);
    
    fTrackCutsPt = new TH2F("fTrackCutsPt","Track Cuts vs p_{T};;p_{T};Count",ncuts,0,ncuts,40,0,fMaxPt);
    (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[0],"Std. TPC Only");
    (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[1],Form("#eta < %3.1f",fMaxEta));
    if (fDoITSCut) (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[2],"1-2 SPD hits");
    if (fDoDCACut) (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[3],"DCA cut");
    (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[4],"p_{T} < 5.0 GeV/c");
    (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[5],"TPCpid (p_{T} < 5.0 GeV/c)");
    (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[6],"TOFpid (p_{T} < 5.0 GeV/c)");
    if (fDemandNoMismatch) (fTrackCutsPt->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[7],"No Mismatch (p_{T} < 5.0 GeV/c)");
    fHistoList->Add(fTrackCutsPt);
    
    fTrackCutsEta = new TH2F("fTrackCutsEta","Track Cuts vs #eta;;#eta;Count",ncuts,0,ncuts,4*binsHisto[1],-fMaxPlotEta,fMaxPlotEta);
    (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[0],"Std. TPC Only");
    (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[1],Form("#eta < %3.1f",fMaxEta));
    if (fDoITSCut) (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[2],"1-2 SPD hits");
    if (fDoDCACut) (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[3],"DCA cut");
    (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[4],"p_{T} < 5.0 GeV/c");
    (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[5],"TPCpid (p_{T} < 5.0 GeV/c)");
    (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[6],"TOFpid (p_{T} < 5.0 GeV/c)");
    if (fDemandNoMismatch) (fTrackCutsEta->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[7],"No Mismatch (p_{T} < 5.0 GeV/c)");
    fHistoList->Add(fTrackCutsEta);
    
    fTrackCutsPhi = new TH2F("fTrackCutsPhi","Track Cuts vs #phi;;#phi;Count",ncuts,0,ncuts,72,0,2.*TMath::Pi());
    (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[0],"Std. TPC Only");
    (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[1],Form("#eta < %3.1f",fMaxEta));
    if (fDoITSCut) (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[2],"1-2 SPD hits");
    if (fDoDCACut) (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[3],"DCA cut");
    (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[4],"p_{T} < 5.0 GeV/c");
    (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[5],"TPCpid (p_{T} < 5.0 GeV/c)");
    (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[6],"TOFpid (p_{T} < 5.0 GeV/c)");
    if (fDemandNoMismatch) (fTrackCutsPhi->GetXaxis())->SetBinLabel(fTrackCutLabelNumbers[7],"No Mismatch (p_{T} < 5.0 GeV/c)");
    fHistoList->Add(fTrackCutsPhi);
    
    fEtaSpectrumTrig = new TH1F("fEtaSpectrumTrig","#eta Spectrum Triggers;#eta;Count",4*binsHisto[1],-fMaxPlotEta,fMaxPlotEta);
    fHistoList->Add(fEtaSpectrumTrig);
    fEtaSpectrumAssoc = new TH2F("fEtaSpectrumAssoc","#eta Spectrum Associateds;#eta;p_{T};Count",4*binsHisto[1],-fMaxPlotEta,fMaxPlotEta,10,0.,5.);
    fHistoList->Add(fEtaSpectrumAssoc); 
    fPhiSpectrumAssoc = new TH2F("fPhiSpectrumAssoc","#phi Spectrum Associateds;#phi;p_{T};Count",72,0,2.*TMath::Pi(),10,0.,5.);
    fHistoList->Add(fPhiSpectrumAssoc);
	
    // --- QA PLOTS PID ---
    
	fTPCnSigmaProton = new TH2F("fTPCnSigmaProton","n#sigma plot for the TPC (Protons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaProton);
	fTOFnSigmaProton = new TH2F("fTOFnSigmaProton","n#sigma plot for the TOF (Protons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaProton);
	fTPCnSigmaPion = new TH2F("fTPCnSigmaPion","n#sigma plot for the TPC (Pions);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaPion);
	fTOFnSigmaPion = new TH2F("fTOFnSigmaPion","n#sigma plot for the TOF (Pions);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaPion);
	fTPCnSigmaKaon = new TH2F("fTPCnSigmaKaon","n#sigma plot for the TPC (Kaons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaKaon);
	fTOFnSigmaKaon = new TH2F("fTOFnSigmaKaon","n#sigma plot for the TOF (Kaons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaKaon);	
	fTPCSignal = new TH3F("fTPCSignal","TPC Signal;#eta;p_{T} GeV/c;dE/dx",25,-fMaxEta,fMaxEta,100,0.,5.,150,0.,300.);
	fHistoList->Add(fTPCSignal);
	fTOFSignal = new TH3F("fTOFSignal","TOF Signal;#eta;p_{T} GeV/c;t",25,-fMaxEta,fMaxEta,100,0.,5.,150,10000.,25000.);
	fHistoList->Add(fTOFSignal);
			
    // --- DI-HADRON CORRELATIONS WITH TPC AND TOF SIGNALS ---
    
	// Di Hadron Correlations with two PID signals, for each p_T bin and particle species separately. (i.e. 30 histo's) 
	// Axes: {Dphi, Deta, TPC signal, TOF signal}
    
    // Unzoomed pictures.
	Int_t binsTPCnoZoom[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                                  {100,100,100,100,100,100,100,100,100,100},
                                  {100,100,100,100,100,100,100,100,100,100}};
		
	Double_t minTPCnoZoom[3][10] = {{   -50., -50.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30.},
                                    {	-100., -50.,-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.},
                                    {  -200.,-150.,-50.,-30.,-20.,-20.,-20.,-20.,-20.,-20.}};
	Double_t maxTPCnoZoom[3][10] = {{ 150., 100., 50., 30., 25., 25., 25., 25., 25., 25.},
                                    {  50., 100., 40., 30., 30., 30., 30., 30., 30., 30.},
                                    { 100.,  50., 50., 30., 30., 30., 30., 30., 30., 30.}};
	
	Int_t binsTOFnoZoom[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                                  {100,100,100,100,100,100,100,100,100,100},
                                  {100,100,100,100,100,100,100,100,100,100}};
	
	Double_t minTOFnoZoom[3][10] = {{-2000.,-2000.,-1000., -500., -500., -500., -500.,	-500., -500., -500.},
                                    {-2000.,-2000.,-2000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.},
                                    {-2000.,-2000.,-6000.,-3000.,-2000.,-1500.,-1500.,-1500.,-1500.,-1500.}};
	Double_t maxTOFnoZoom[3][10] = {{ 2000., 2000., 5000., 3000., 2000., 1500., 1500., 1500., 1500., 1500.},
                                    { 2000., 2000., 5000., 2000., 1500., 1500., 1500., 1500., 1500., 1500.},
                                    { 2000., 2000., 2000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.}};
    
    // Zoomed pictures.
    Int_t binsTPCZoom[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                                {100,100,100,100,100,100,100,100,100,100},
                                {100,100,100,100,100,100,100,100,100,100}};
    
	Double_t minTPCZoom[3][10] = {{   -20., -20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.},
                                  {	 -20., -20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.},
                                  {   -40., -20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.}};
	Double_t maxTPCZoom[3][10] = {{  20.,  20., 20., 20., 20., 20., 20., 20., 20., 20.},
                                  {  20.,  20., 20., 20., 20., 20., 20., 20., 20., 20.},
                                  {  40.,  20., 20., 20., 20., 30., 30., 30., 30., 30.}};
	
	Int_t binsTOFZoom[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                                {100,100,100,100,100,100,100,100,100,100},
                                {100,100,100,100,100,100,100,100,100,100}};
	
	Double_t minTOFZoom[3][10] = {{-1000.,-1000., -500., -500., -500., -400., -400., -400., -400., -400.},
                                  { -800., -800., -800., -800., -800., -600., -500., -500., -400., -400.},
                                  {-1000.,-1000.,-1000.,-1000., -800.,-1000.,-1000., -800., -700., -700.}};
	Double_t maxTOFZoom[3][10] = {{ 1000., 1000., 1000., 1000., 1000., 1000., 1000.,  900.,  800.,  700.},
                                  { 1000., 1000.,  500.,  500.,  500.,  900.,  700.,  700.,  600.,  500.},
                                  { 1000., 1000., 1000.,  500.,  500.,  500.,  500.,  500.,  500.,  500.}};
    
	TString basenameTPCTOF("fDiHadronTPCTOF");
	TString basetitle("Di-Hadron correlation");
	TString finalname, finaltitle;
	TString species[3] = {"Pion","Kaon","Proton"};
	TString ptbins[11] = {"0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"};
    
	// Recall that AliPID::kPion = 2, AliPID::kKaon = 3, AliPID::kProton = 4.
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
				
		for (Int_t iPtBin = 0; iPtBin < 10; iPtBin++) {
			
            if (fZoomed) {
                binsHisto[2] = binsTPCZoom[iSpecies][iPtBin];
                binsHisto[3] = binsTOFZoom[iSpecies][iPtBin];
                minHisto[2] = minTPCZoom[iSpecies][iPtBin];
                minHisto[3] = minTOFZoom[iSpecies][iPtBin];
                maxHisto[2] = maxTPCZoom[iSpecies][iPtBin];
                maxHisto[3] = maxTOFZoom[iSpecies][iPtBin];
            } else {
                binsHisto[2] = binsTPCnoZoom[iSpecies][iPtBin];
                binsHisto[3] = binsTOFnoZoom[iSpecies][iPtBin];
                minHisto[2] = minTPCnoZoom[iSpecies][iPtBin];
                minHisto[3] = minTOFnoZoom[iSpecies][iPtBin];
                maxHisto[2] = maxTPCnoZoom[iSpecies][iPtBin];
                maxHisto[3] = maxTOFnoZoom[iSpecies][iPtBin];
            }

            // Make the di-hadron correlations with different pid signals.
            // TODO: Rewrite this with Form().
            finaltitle = basetitle;
			(((((finaltitle += " (") += species[iSpecies]) += ") ") += ptbins[iPtBin]) += " < P_t < ") += ptbins[iPtBin+1];
            finaltitle+=";#Delta#phi;#Delta#eta;dE/dx;t (ms)";
            finalname = basenameTPCTOF;
			(((finalname += "_") += species[iSpecies]) += "_") += iPtBin;
            
            fDiHadronTPCTOF[iSpecies][iPtBin] = new THnSparseF(finalname,finaltitle,4,binsHisto,minHisto,maxHisto);
            fHistoList->Add(fDiHadronTPCTOF[iSpecies][iPtBin]);
            
		}
	}
    
    // --- MIXED EVENTS ---
    
    if (fCalculateMixedEvents) {
        fMixedEvents = new TH3F("fMixedEvents","Mixed Events;#Delta#phi;#Delta#eta;p_{T_assoc}",binsHisto[0],minHisto[0],maxHisto[0],binsHisto[1],minHisto[1],maxHisto[1],10,0.,5.);
        fMixedEvents->Sumw2();
        fHistoList->Add(fMixedEvents);
    }
    
	PostData(1, fHistoList);
	
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskDiHadronPID::ClassifyTrack(AliAODTrack *track)

{
	//
    // This function both classifies tracks, and fills the track QA histograms.
    //
    // Classifies the track in:
    //  0 -> Not Useful,
    //  1 -> Associated,
    //  2 -> Trigger.
    //
    // IDEA: later we can do this with filterbits.
    //
    // The following track cuts are applied for triggers:
    //
    // 1a) pT > 5.0 GeV/c, pT < fPtMax,
    //  2) StandardTPCOnlyTrackCuts (filterbit 7 for AOD's),
    //  3) eta < fMaxEta,
    //  4) ITS track cut, a track is only selected if it has at least one hit in the SPD,
    //     that is, in one of the first two layers of the ITS. (can be switched on and off
    //     using (Bool_t)fDoITSCut),
    //  5) DCA track cut, of all tracks with at least one SPD hit, the DCA is constrained as
    //     a function of pT, in order to remove lots of the secondaries. (can be switched on and 
    //     off by using (Bool_t)fDoDCACut),
    //
    // For associateds all the same track cuts are applied, except 1a is
    // replaced by 2b. Moreover the following cuts are applied too:
    //
    // 1b) pT < 5.0 GeV/c,
    //  6) TPCpid,
    //  7) TOFpid,
    //  8) no TPCTOF mismatch.
    //
    
    // Check if a track is supplied.
    if (!track) {
        if (fVerbose>0) cout<<"AliAnalysisTaskDiHadronPID::SelectTrack -> ERROR: No track found."<<endl;
        return 0;
    }
    
    // Get basic track information.
    Int_t classification=0;
    Double_t pt = track->Pt();
    Double_t eta = track->Eta();
    Double_t phi = track->Phi();
    	
    // 1) pT cut: First separation between triggers and associateds.
    if (pt<5.0) classification=1;
    if ((pt>5.0)&&(pt<fMaxPt)) classification = 2;
    if (!classification) return 0;
    
    // 2) StandardTPCOnlyTrackCuts.
    if (!(track->TestFilterMask(1<<7))) {
		//if (fVerbose>3) cout<<"Track Ignored: Did Not pass filterbit."<<endl;
		return 0;
	}
       
	if (fVerbose>3) {
		cout<<endl;
		cout<<"pt: "<<pt<<" eta: "<<eta<<" phi: "<<phi<<endl;
	}
	
    fTrackCutsCount->Fill(fTrackCutLabelNumbers[0]-.5);
    fTrackCutsPt->Fill(fTrackCutLabelNumbers[0]-.5,pt);
    fTrackCutsEta->Fill(fTrackCutLabelNumbers[0]-.5,eta);
    fTrackCutsPhi->Fill(fTrackCutLabelNumbers[0]-.5,phi);
    
    // 3) eta cut.
    if (TMath::Abs(eta)>fMaxEta) {
		if (fVerbose>3) cout<<"Track Ignored: Eta too large."<<endl;
		return 0;
    }
	
    fTrackCutsCount->Fill(fTrackCutLabelNumbers[1]-.5);
    fTrackCutsPt->Fill(fTrackCutLabelNumbers[1]-.5,pt);
    fTrackCutsEta->Fill(fTrackCutLabelNumbers[1]-.5,eta);
    fTrackCutsPhi->Fill(fTrackCutLabelNumbers[1]-.5,phi);
    
    // Obtaining ITS information.
    AliAODTrack* globaltrack = GetGlobalTrack(track);
    Bool_t ITSLayerHit[6];
    for (Int_t iITSLayer=0; iITSLayer<6; iITSLayer++) {
        ITSLayerHit[iITSLayer] = globaltrack->HasPointOnITSLayer(iITSLayer);
    }
    Int_t SPDHits=ITSLayerHit[0]+ITSLayerHit[1];
    
	if (fVerbose>3) cout<<"SPD hits: "<<SPDHits<<endl;
	
    // Fill the ITS hist.
    fITSHits->Fill(SPDHits+0.5);

    // 4) ITS cut.
    if (fDoITSCut) {
        if (!SPDHits) {
			if (fVerbose>3) cout<<"Track Ignored: Not enough SPD hits."<<endl;
			return 0;
		}
        fTrackCutsCount->Fill(fTrackCutLabelNumbers[2]-.5);
        fTrackCutsPt->Fill(fTrackCutLabelNumbers[2]-.5,pt);
        fTrackCutsEta->Fill(fTrackCutLabelNumbers[2]-.5,eta);
        fTrackCutsPhi->Fill(fTrackCutLabelNumbers[2]-.5,phi);
    }
        
    // Propagate the global track to the DCA.
    Double_t PosAtDCA[2] = {-999,-999};
    Double_t covar[3] = {-999,-999,-999};
    globaltrack->PropagateToDCA(fAODVertex,fAODEvent->GetMagneticField(),100.,PosAtDCA,covar);
        
    // Fill the DCA hist (before cut)
    fDCA->Fill(PosAtDCA[0],PosAtDCA[1]);
    fDCAZoomed->Fill(PosAtDCA[0],PosAtDCA[1]);
    fDCAZoomedTwice->Fill(PosAtDCA[0],PosAtDCA[1]);
    
    // 5) DCA cut (See R_AA paper).
    Double_t DCAcutvalue[2];
    DCAcutvalue[0] = 0.018 + 0.035*TMath::Power(pt,-1.01);
    DCAcutvalue[1] = 2.; 
    if (fDoDCACut) {
        if (SPDHits&&((TMath::Abs(PosAtDCA[0])>DCAcutvalue[0])||(TMath::Abs(PosAtDCA[1])>DCAcutvalue[1]))) {
			if (fVerbose>3) cout<<"Track Ignored: Enough SPD hits, but out of DCA range."<<endl;
			return 0;
		}
        fTrackCutsCount->Fill(fTrackCutLabelNumbers[3]-.5);
        fTrackCutsPt->Fill(fTrackCutLabelNumbers[3]-.5,pt);
        fTrackCutsEta->Fill(fTrackCutLabelNumbers[3]-.5,eta);
        fTrackCutsPhi->Fill(fTrackCutLabelNumbers[3]-.5,phi);
        
        // Fill the DCA hist (after cut)
        fDCACut->Fill(PosAtDCA[0],PosAtDCA[1]);
        fDCAZoomedCut->Fill(PosAtDCA[0],PosAtDCA[1]);
        fDCAZoomedTwiceCut->Fill(PosAtDCA[0],PosAtDCA[1]);
    }
    
    // Now all the common cuts have been performed. Tracks identified as a trigger will
    // be returned. Note that they will not appear in the histograms after track cut 5.
    if (classification==2) return 2;
    
    // Now we're left with only tracks with pT < 5.0 GeV/c. 
    fTrackCutsCount->Fill(fTrackCutLabelNumbers[4]-.5);
    fTrackCutsPt->Fill(fTrackCutLabelNumbers[4]-.5,pt);
    fTrackCutsEta->Fill(fTrackCutLabelNumbers[4]-.5,eta);
    fTrackCutsPhi->Fill(fTrackCutLabelNumbers[4]-.5,phi);
    
    // Obtain the status of the track.
    ULong_t status = globaltrack->GetStatus();
    
    // 6) TPCpid
    if (!((status&AliAODTrack::kTPCpid)==AliAODTrack::kTPCpid)) {
		if (fVerbose>3) cout<<"Track Ignored: No TPC pid."<<endl;
		return 0;
	}
    
    fTrackCutsCount->Fill(fTrackCutLabelNumbers[5]-.5);
    fTrackCutsPt->Fill(fTrackCutLabelNumbers[5]-.5,pt);
    fTrackCutsEta->Fill(fTrackCutLabelNumbers[5]-.5,eta);
    fTrackCutsPhi->Fill(fTrackCutLabelNumbers[5]-.5,phi);
    
    // 7) TOFpid
    if (!((status&AliAODTrack::kTOFpid)==AliAODTrack::kTOFpid)) {
		if (fVerbose>3) cout<<"Track Ignored: No TOF pid."<<endl;
		return 0;
	}
    
    fTrackCutsCount->Fill(fTrackCutLabelNumbers[6]-.5);
    fTrackCutsPt->Fill(fTrackCutLabelNumbers[6]-.5,pt);
    fTrackCutsEta->Fill(fTrackCutLabelNumbers[6]-.5,eta);
    fTrackCutsPhi->Fill(fTrackCutLabelNumbers[6]-.5,phi);
    
    // 8) TPC, TOF mismatch.
    if (fDemandNoMismatch) {
        if ((status&AliAODTrack::kTOFmismatch)==AliAODTrack::kTOFmismatch) {
			if (fVerbose>3) cout<<"Track Ignored: TOF mismatch found."<<endl;
			return 0;
		}
        fTrackCutsCount->Fill(fTrackCutLabelNumbers[7]-.5);
        fTrackCutsPt->Fill(fTrackCutLabelNumbers[7]-.5,pt);
        fTrackCutsEta->Fill(fTrackCutLabelNumbers[7]-.5,eta);
        fTrackCutsPhi->Fill(fTrackCutLabelNumbers[7]-.5,phi);
    }
    
    // All tracks which made it up to here are classified as associateds.
    
    // Fill the PID QA histograms.
    Double_t mom, nSigma;
    
    mom = globaltrack->GetTPCmomentum();
    nSigma = fPIDResponse->NumberOfSigmasTPC(globaltrack,AliPID::kProton);
    fTPCnSigmaProton->Fill(mom,nSigma);
    nSigma = fPIDResponse->NumberOfSigmasTPC(globaltrack,AliPID::kPion);
    fTPCnSigmaPion->Fill(mom,nSigma);
    nSigma = fPIDResponse->NumberOfSigmasTPC(globaltrack,AliPID::kKaon);
    fTPCnSigmaKaon->Fill(mom,nSigma);
    
    fTPCSignal->Fill(eta,pt,globaltrack->GetTPCsignal());
    
    mom =globaltrack->P();
    nSigma = fPIDResponse->NumberOfSigmasTOF(globaltrack,AliPID::kProton);
    fTOFnSigmaProton->Fill(mom,nSigma);			
    nSigma = fPIDResponse->NumberOfSigmasTOF(globaltrack,AliPID::kPion);
    fTOFnSigmaPion->Fill(mom,nSigma);	
    nSigma = fPIDResponse->NumberOfSigmasTOF(globaltrack,AliPID::kKaon);
    fTOFnSigmaKaon->Fill(mom,nSigma);
    
    fTOFSignal->Fill(eta,pt,globaltrack->GetTOFsignal());
    
    // Return associated.
    return 1;
    
} 
    
//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadronPID::SelectEvent(AliAODVertex* vertex)

{
	
	//
	// Event Selection.
	//
	
	Double_t primVtx[3];
	vertex->GetXYZ(primVtx);
	if (TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2])>10.) {
		if (fVerbose>2) cout << "AliAnalysisTaskDiHadronPID::SelectEvent -> Vertex Out of Range." << endl;
        if (fVerbose>2) cout << "AliAnalysisTaskDiHadronPID::SelectEvent -> Event not selected." << endl;
		return kFALSE;
	}
	if (fVerbose>2) cout << "AliAnalysisTaskDiHadronPID::SelectEvent -> Vertex is OK." << endl;
    
    // We also wish to make centrality cut, but only if run on Pb+Pb.
    if (fBeamType=="PbPb") {
        Double_t cent = fAODHeader->GetCentrality();
        if (cent>fCentralityCutMin||cent<fCentralityCutMax) {
            if (fVerbose>2) cout<<"AliAnalysisTaskDiHadronPID::SelectEvent -> Event did not pass centaltity cut."<<endl;
            return kFALSE;
        }
    }

    if (fVerbose>2) cout << "AliAnalysisTaskDiHadronPID::SelectEvent -> Event selected." << endl;
	return kTRUE;
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::FillGlobalTracksArray() {

	// Initialize the mapping for corresponding PID tracks. (see
	// GetGlobalTrack(AliAODTrack* track)). 
	//
	
	if (!fAODEvent) {
        if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::FillGlobalTracksArray -> ERROR: fAODEvent not set." << endl;
		return;
	}
	
	fGlobalTracks = new TObjArray();
	AliAODTrack* track = 0x0;
		
	for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++) {
		
		track = fAODEvent->GetTrack(iTrack);
		
		// I.e., if it does NOT pass the filtermask.
		if (!(track->TestFilterMask(1<<7))) {
            if (track->GetID()>-1) fGlobalTracks->AddAtAndExpand(track,track->GetID());
            //if (track->GetID()<1) cout<<"Track ID: "<<track->GetID()<<" Partner ID: "<<(-track->GetID()-1)<<endl;
		}
        
	}
	
}

//_____________________________________________________________________________
AliAODTrack* AliAnalysisTaskDiHadronPID::GetGlobalTrack(AliAODTrack* track) {
	
	//
	// Returns the "parner track" of track, which contains the pid information.
	//
	
	AliAODTrack* partner = 0x0;
	    
    partner = (AliAODTrack*)(fGlobalTracks->At(-track->GetID()-1));
	    
	if (!partner&&(fVerbose>3)) cout<<"AliAnalysisTaskDiHadronPID::GetGlobalTrack -> No Global Track Found!"<<endl;
	
	return partner;
	
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskDiHadronPID::PhiRange(Double_t DPhi)

{
	//
	// Puts the argument in the range [-pi/2,3 pi/2].
	//
	
	if (DPhi < -TMath::Pi()/2) DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	

	return DPhi;
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::UserExec(Option_t *)

{
	//
	// UserExec.
	//
	
	// Input the event.
	fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fAODEvent) {
		if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::UserExec -> ERROR: No AliAODEvent pointer could be created." << endl;
		return;
	}
	
	// Get the event header.
	fAODHeader = fAODEvent->GetHeader();
	if (!fAODHeader) {
		if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::UserExec -> ERROR: No AliAODHeader pointer could be created."<<endl;
		return;
	}
	
	// Get the event vertex.
	fAODVertex = fAODEvent->GetPrimaryVertex();
	if (!fAODVertex) {
		if (fVerbose>0) cout << "AliAnalysisTaskDiHadronPID::UserExec -> ERROR: No AliAODVertex pointer could be created." << endl;
		return;
	}

    // See if the event passes the event selection.
    if (!SelectEvent(fAODVertex)) return;
    
	// Display basic event information.
	if ((fVerbose>2)) cout << endl;
	if ((fVerbose>2)&&(fBeamType=="PbPb")) cout << "Event centrality: " << fAODHeader->GetCentrality() << endl;
	if ((fVerbose>2)) cout << "Event Vertex Z: " << fAODVertex->GetZ() << endl;
	if ((fVerbose>2)) cout << "Event tracks in AOD: " << fAODEvent->GetNumberOfTracks() << endl;
	if ((fVerbose>2)) cout << endl;

    // Filling Event QA plots.
	if (fBeamType=="PbPb") fCentrality->Fill(fAODHeader->GetCentrality());
	fVertexZ->Fill(fAODVertex->GetZ());
    
	// Fill the TObjArray which holds Global tracks.
	FillGlobalTracksArray();
	
	// Create object arrays for triggers and associateds.
	TObjArray *triggers		= new TObjArray();
	TObjArray *associateds	= new TObjArray();
		    
    // In this loop the triggers and associateds will be identified, track QA and PID QA histograms will be filled.
	for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++) {
        
        // Obtain a pointer to the track.
		fAODTrack = fAODEvent->GetTrack(iTrack);
        if (!fAODTrack&&(fVerbose>0)) {
            cout << "AliAnalysisTaskDiHadronPID::UserExec -> ERROR: Track object not found." << endl;
            continue;
        }
        
        // Find the track classification.
        Int_t tracktype = ClassifyTrack(fAODTrack);
		
        if (tracktype==0) {
			continue;
		}
		
        if (tracktype==1) {
			if (fVerbose>3) cout<<"Track added to associated buffer."<<endl;
            associateds->AddLast(fAODTrack);
            fEtaSpectrumAssoc->Fill(fAODTrack->Eta(),fAODTrack->Pt());
            fPhiSpectrumAssoc->Fill(fAODTrack->Phi(),fAODTrack->Pt());
        }
        
        if (tracktype==2) {
			if (fVerbose>3) cout<<"Track added to trigger buffer."<<endl;
            triggers->AddLast(fAODTrack);
            fEtaSpectrumTrig->Fill(fAODTrack->Eta());
        }
        
        
    }
    
    // In This Loop the di-hadron correlation will be made.
    Double_t histoFill[4];                  // {Dphi, Deta, TPC signal, TOF signal}
    AliAODTrack* currentTrigger = 0x0;
	AliAODTrack* currentAssociated = 0x0;
	AliAODTrack* currentAssociatedGlobal = 0x0;
    
    AliTPCPIDResponse& TPCPIDResponse = fPIDResponse->GetTPCResponse();
    
    for (Int_t iTrig = 0; iTrig < triggers->GetEntriesFast(); iTrig++){
        
		currentTrigger = (AliAODTrack*)(triggers->At(iTrig));
		
		for (Int_t iAssoc = 0; iAssoc < associateds->GetEntriesFast(); iAssoc++) {
            
			currentAssociated = (AliAODTrack*)(associateds->At(iAssoc));
			currentAssociatedGlobal = GetGlobalTrack(currentAssociated);
            
			Double_t pt = currentAssociated->Pt();
			histoFill[0] = PhiRange(currentTrigger->Phi() - currentAssociated->Phi());
			histoFill[1] = currentTrigger->Eta() - currentAssociated->Eta();

			// Is there a caveat here when Pt = 5.00000000?
			const Int_t ptbin = (Int_t)(2*pt);
			
			if (currentAssociatedGlobal) {
				
                // Get TPC (expected) signals.
				Double_t TPCmom = currentAssociatedGlobal->GetTPCmomentum();
                Double_t TPCsignal = currentAssociatedGlobal->GetTPCsignal();
				Double_t expectedTPCsignalPion = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kPion);
                Double_t expectedTPCsignalKaon = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kKaon);
				Double_t expectedTPCsignalProton = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kProton);

                // Get TOF (expected) signals.
                Double_t TOFsignal = currentAssociatedGlobal->GetTOFsignal();
                Double_t times[AliPID::kSPECIES];
                currentAssociatedGlobal->GetIntegratedTimes(times);
                Double_t expectedTOFsignalPion = times[AliPID::kPion];
                Double_t expectedTOFsignalKaon = times[AliPID::kKaon];
                Double_t expectedTOFsignalProton = times[AliPID::kProton]; 
            
                // Fill the histograms.
                histoFill[2] = TPCsignal - expectedTPCsignalPion;
                histoFill[3] = TOFsignal - expectedTOFsignalPion;
                fDiHadronTPCTOF[0][ptbin]->Fill(histoFill);
                
                histoFill[2] = TPCsignal - expectedTPCsignalKaon;
                histoFill[3] = TOFsignal - expectedTOFsignalKaon;
                fDiHadronTPCTOF[1][ptbin]->Fill(histoFill);
                
                histoFill[2] = TPCsignal - expectedTPCsignalProton;
                histoFill[3] = TOFsignal - expectedTOFsignalProton;
                fDiHadronTPCTOF[2][ptbin]->Fill(histoFill);
                
			}
		}
	}

    // In this loop we calculate the mixed events.
    if (fCalculateMixedEvents) {
        
        // Loop over the trigger buffer.
        if (fVerbose>3) cout << "AliAnalysisTaskDiHadronPID::UserExec -> Mixing the events with "<<fTrigBufferSize<<" triggers from the buffer." <<endl;
        if (fVerbose>3) cout << "AliAnalysisTaskDiHadronPID::UserExec -> Buffer size: "<<fTrigBufferIndex<<endl;
        
        for (Int_t iTrig=0;iTrig<fTrigBufferSize;iTrig++) {
            
            // Check if the trigger and the associated have a reconstructed
            // vertext no further than 2cm apart.
            
            // fTrigBuffer[i][0] = z
            // fTrigBuffer[i][1] = phi
            // fTrigBuffer[i][2] = eta
            // fTrigBuffer[i][3] = p_t
            
            if (TMath::Abs(fTrigBuffer[iTrig][0]-fAODVertex->GetZ())<fVertexZMixedEvents) {
                
                if (fVerbose>3) cout<<"AliAnalysisTaskDiHadronPID::UserExec -> Mixing with trigger Z: "<<fTrigBuffer[iTrig][0]<<", Pt: "<<fTrigBuffer[iTrig][3]<<endl;
                
                for (Int_t iAssoc = 0; iAssoc < associateds->GetEntriesFast(); iAssoc++) {
                    
                    currentAssociated = (AliAODTrack*)(associateds->At(iAssoc));
                    currentAssociatedGlobal = GetGlobalTrack(currentAssociated);
                    
                    if (currentAssociatedGlobal) {
                        
                        Double_t DPhi = PhiRange(fTrigBuffer[iTrig][1] - currentAssociated->Phi());
                        Double_t DEta = fTrigBuffer[iTrig][2] - currentAssociated->Eta();
                        Double_t pt = currentAssociated->Pt();
                        
                        fMixedEvents->Fill(DPhi,DEta,pt);
			
                    }
                }
            }
        }
    
        // Copy the triggers from the current event into the buffer.
        if (fAODVertex->GetZ()<10.) {
        
            if (fVerbose>3) cout<<"AliAnalysisTaskDiHadronPID::UserExec -> Copying "<<triggers->GetEntriesFast()<<" triggers with vertex z = "<<fAODVertex->GetZ()<<" to the buffer."<<endl;
        
            for (Int_t iTrig = 0; iTrig<triggers->GetEntriesFast(); iTrig++) {
            
                currentTrigger = (AliAODTrack*)(triggers->At(iTrig));
                if (fVerbose>3) cout<<"AliAnalysisTaskDiHadronPID::UserExec -> Trigger pt = "<<currentTrigger->Pt()<<endl;
            
                fTrigBuffer[fTrigBufferIndex][0] = fAODVertex->GetZ();
                fTrigBuffer[fTrigBufferIndex][1] = currentTrigger->Phi();
                fTrigBuffer[fTrigBufferIndex][2] = currentTrigger->Eta();
                fTrigBuffer[fTrigBufferIndex][3] = currentTrigger->Pt();
                fTrigBufferIndex++;
                if (fTrigBufferSize<fTrigBufferMaxSize) {fTrigBufferSize++;} // 250 triggers should be enough to get 10 times more data in mixed events.
                if (fTrigBufferIndex==fTrigBufferMaxSize) {fTrigBufferIndex=0;}
            }
        }
    }        
        
    if (fPrintBufferSize) cout<<"AliAnalysisTaskDiHadronPID::UserExec -> Trigger buffer index: "<<fTrigBufferIndex<<", and size: "<<fTrigBufferSize<<endl;
    
	delete triggers;
	delete associateds;
                          
	PostData(1,fHistoList);
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::Terminate(Option_t *)

{
	//
	// Terminate.
    //
    
}

