/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TList.h>
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "AliLog.h"

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVEvent.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskJetCoreEmcal.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetCoreEmcal);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskJetCoreEmcal::AliAnalysisTaskJetCoreEmcal() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
	fCentMin(0.),
	fCentMax(100.),
	fTTLowRef(8.),
	fTTUpRef(9.),
	fTTLowSig(20.),
	fTTUpSig(50.),
	fNRPBins(50),
	fFrac(0.8),
	fJetEtaMin(-.5),
	fJetEtaMax(.5),
	fJetHadronDeltaPhi(0.6),
	fJetContName(""),
	fFillTrackHistograms(kTRUE),
	fFillJetHistograms(kTRUE),
	fRandom(0),
	fHistEvtSelection(0x0), 
	fHJetSpec(0x0),
	fh1TrigRef(0x0),
	fh1TrigSig(0x0),
	fh2Ntriggers(0x0),
	fh2RPJetsC10(0x0),
	fh2RPJetsC20(0x0),
	fh2RPTC10(0x0),
	fh2RPTC20(0x0),
	fHJetPhiCorr(0x0),
	fhDphiPtSig(0x0),
	fhDphiPtRef(0x0)
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskJetCoreEmcal::AliAnalysisTaskJetCoreEmcal(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
	fCentMin(0.),
	fCentMax(100.),
	fTTLowRef(8),
	fTTUpRef(9.),
	fTTLowSig(20.),
	fTTUpSig(50.),
	fNRPBins(50),
	fFrac(0.8),
	fJetEtaMin(-.5),
	fJetEtaMax(.5),
	fJetHadronDeltaPhi(0.6),
	fJetContName(""),
	fFillTrackHistograms(kTRUE),
	fFillJetHistograms(kTRUE),
	fRandom(0),
	fHistEvtSelection(0x0), 
	fHJetSpec(0x0),
	fh1TrigRef(0x0),
	fh1TrigSig(0x0),
	fh2Ntriggers(0x0),
	fh2RPJetsC10(0x0),
	fh2RPJetsC20(0x0),
	fh2RPTC10(0x0),
	fh2RPTC20(0x0),
	fHJetPhiCorr(0x0),
	fhDphiPtSig(0x0),
	fhDphiPtRef(0x0)
{
  SetMakeGeneralHistograms(kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskJetCoreEmcal::~AliAnalysisTaskJetCoreEmcal()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskJetCoreEmcal::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if(fFillTrackHistograms) AllocateTrackHistograms();
  if(fFillJetHistograms) AllocateJetHistograms();
  //AllocateClusterHistograms();
  //AllocateCellHistograms();
  AllocateJetCoreHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateClusterHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The cluster containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of clusters;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNClusters";
  histtitle = TString::Format("%s;Sum of n clusters;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateCellHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname(fCaloCellsName);

  fHistManager.CreateHistoGroup(groupname);
  for (Int_t cent = 0; cent < fNcentBins; cent++) {
    histname = TString::Format("%s/histCellEnergy_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;#it{E}_{cell} (GeV);counts", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 300, 0, 150);

    histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;number of cells;events", histname.Data());
    if (fForceBeamType != kpp) {
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000);
    }
    else {
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }
  }
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

        histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
      }

      histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
      }

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
      }
    }
  }
}

void AliAnalysisTaskJetCoreEmcal::AllocateJetCoreHistograms() 
{

	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);

	// set seed
	fRandom = new TRandom3(0);

	fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
	fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
	fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
	fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");

	fh1TrigRef=new TH1D("Trig Ref","",10,0.,10);
	fh1TrigSig=new TH1D("Trig Sig","",10,0.,10);  
	fh2Ntriggers=new TH2F("# of triggers","",100,0.,100.,50,0.,50.);
	fh2RPJetsC10=new TH2F("RPJetC10","",35,0.,3.5,100,0.,100.);
	fh2RPJetsC20=new TH2F("RPJetC20","",35,0.,3.5,100,0.,100.); 
	fh2RPTC10=new TH2F("RPTriggerC10","",35,0.,3.5,50,0.,50.); 
	fh2RPTC20=new TH2F("RPTriggerC20","",35,0.,3.5,50,0.,50.);  

	fOutput->Add(fHistEvtSelection);

	fOutput->Add(fh1TrigRef);
	fOutput->Add(fh1TrigSig); 
	fOutput->Add(fh2Ntriggers);
	fOutput->Add(fh2RPJetsC10);
	fOutput->Add(fh2RPJetsC20);
	fOutput->Add(fh2RPTC10);
	fOutput->Add(fh2RPTC20);

	const Int_t dimSpec = 6;
	const Int_t nBinsSpec[dimSpec]     = {100,100, 280, 50,200, fNRPBins};
	const Double_t lowBinSpec[dimSpec] = {0,0,-80, 0,-0.5*TMath::Pi(), 0};
	const Double_t hiBinSpec[dimSpec]  = {100,1, 200, 50,1.5*TMath::Pi(),  static_cast<Double_t>(fNRPBins)};
	fHJetSpec = new THnSparseF("fHJetSpec","Recoil jet spectrum",dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);

	// comment out since I want finer binning in jet area, to make it easier
	// to change selection on jet area (Leticia used 0.8*R^2*Pi whereas 0.6 is used
	// for inclusive jets)
	//change binning in jet area
//	Double_t *xPt6 = new Double_t[7];
//	xPt6[0] = 0.;
//	xPt6[1]=0.07;
//	xPt6[2]=0.2;
//	xPt6[3]=0.4;
//	xPt6[4]=0.6;
//	xPt6[5]=0.8; 
//	xPt6[6]=1;
//	fHJetSpec->SetBinEdges(1,xPt6);
//	delete [] xPt6;

	fOutput->Add(fHJetSpec);  

	// azimuthal correlation

	fhDphiPtSig = new TH2F("hDphiPtS","recoil #Delta #phi vs jet pT signal",100,-2,5,250,-50,200);  
	fhDphiPtSig->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtSig->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 
	fhDphiPtRef = new TH2F("hDphiPtR","recoil #Delta #phi vs jet pT reference",100,-2,5,250,-50,200);  
	fhDphiPtRef->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtRef->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 

	fOutput->Add(fhDphiPtRef);  
	fOutput->Add(fhDphiPtSig);  

	// =========== Switch on Sumw2 for all histos ===========
	for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
		TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
		if (h1){
			h1->Sumw2();
			continue;
		}
		THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
		if (hn){
			hn->Sumw2();
		}	  
	}

	// add QA plots from fEventCuts
	fEventCuts.AddQAplotsToList(fOutput);

	TH1::AddDirectory(oldStatus);

	PostData(1, fOutput);
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskJetCoreEmcal::FillHistograms()
{

	fHistEvtSelection->Fill(1); // number of events before event selection
	AliVEvent *ev = InputEvent();
	if (!fEventCuts.AcceptEvent(ev)) {
		fHistEvtSelection->Fill(2);
		return kTRUE;
	}

  if(fFillJetHistograms) DoJetLoop();
  if(fFillTrackHistograms) DoTrackLoop();
  //DoClusterLoop();
  //DoCellLoop();
	DoJetCoreLoop();

  return kTRUE;
}

void AliAnalysisTaskJetCoreEmcal::DoJetCoreLoop()
{

	// Do jet core analysis and fill histograms.
	AliJetContainer *jetCont = GetJetContainer(fJetContName);
	if(!jetCont) {
		AliError(Form("jet container not found - check name %s",fJetContName.Data()));
		TIter next(&fJetCollArray);
		while ((jetCont = static_cast<AliJetContainer*>(next())))
			AliError(Form("%s",jetCont->GetName()));
		AliFatal("Exit...");
		return;
	}

	// centrality selection 
	if(fDebug) Printf("centrality: %f\n", fCent);
	if ((fCent>fCentMax) || (fCent<fCentMin)) {
		fHistEvtSelection->Fill(4);
		return;
	}
	fHistEvtSelection->Fill(0); 

	// Background
	Double_t rho = 0;
	if (jetCont->GetRhoParameter()) rho = jetCont->GetRhoVal(); 
	if(fDebug) Printf("rho = %f, rho check  = %f",rho, GetRhoVal(0));

	// Choose trigger track
	Int_t nT=0;
	TList ParticleList;
	Double_t minT=0;
	Double_t maxT=0;
	Int_t number=0;
	Double_t dice=fRandom->Uniform(0,1);
	Bool_t isSignal = kFALSE;
	if(dice>fFrac){ 
		minT=fTTLowRef;
		maxT=fTTUpRef;
	}
	if(dice<=fFrac){
		isSignal = kTRUE;
		minT=fTTLowSig;
		maxT=fTTUpSig;
	} 
	nT=SelectTrigger(&ParticleList,minT,maxT,number);
	if(fDebug) Printf("%s class ---> n triggers between %f and %f = %i, index of trigger chosen = %i",dice>fFrac?"ref.":"sig.",minT,maxT,number,nT);
	if(nT<0) return;

	if(dice>fFrac) fh1TrigRef->Fill(number);
	if(dice<=fFrac)fh1TrigSig->Fill(number);


	// particle loop - 
	for(Int_t tt=0;tt<ParticleList.GetEntries();tt++){
		// histogram 
		//if(fHardest==0||fHardest==1){if(tt!=nT) continue;}
		if(tt!=nT) continue;
		AliVParticle *partback = (AliVParticle*)ParticleList.At(tt);     
		if(!partback) continue;
		if(fDebug) Printf("trigger particle pt = %f \teta = %f \t phi = %f",partback->Pt(),partback->Eta(),partback->Phi());
		//     if(partback->Pt()<8) continue;

		Int_t injet4=0;
		Int_t injet=0; 

    fh2Ntriggers->Fill(fCent,partback->Pt());
    Double_t phiBinT = RelativePhi(partback->Phi(),fEPV0);
    if(fCent<20.) fh2RPTC20->Fill(TMath::Abs(phiBinT),partback->Pt());
    if(fCent<10.) fh2RPTC10->Fill(TMath::Abs(phiBinT),partback->Pt());

		Double_t etabig=0;
		Double_t ptbig=0;
		Double_t areabig=0;
		Double_t phibig=0.;
		//   Double_t areasmall=0;

		TString histname;
		TString groupname;
		groupname = jetCont->GetName();
		UInt_t count = 0;
		for(auto jetbig : jetCont->accepted()) {
			if (!jetbig) continue;
			count++;
			ptbig   = jetbig->Pt();
			etabig  = jetbig->Eta();
			phibig  = jetbig->Phi();
			if(ptbig==0) continue; 
			Double_t phiBin = RelativePhi(phibig,fEPV0); //relative phi between jet and ev. plane
			areabig = jetbig->Area();
			Double_t ptcorr=ptbig-rho*areabig;
			//JJJ - perhaps should change eta selection if implemented in jet container
			if((etabig<fJetEtaMin)||(etabig>fJetEtaMax)) continue; 
			if(areabig>=0.07) injet=injet+1;
			if(areabig>=0.4) injet4=injet4+1;   
			Double_t dphi=RelativePhi(partback->Phi(),phibig); 
			if(fDebug) Printf("jet properties...\n\teta = %f \t phi = %f \t pt = %f \t relativephi = %f\t area = %f\t rho = %f",etabig,phibig,ptbig,dphi,areabig,rho);

			// do azimuthal correlation analysis
			// dPhi between -0.5 < dPhi < 1.5
			Double_t dPhiShift=phibig-partback->Phi();
			if(dPhiShift>2*TMath::Pi()) dPhiShift -= 2*TMath::Pi();
			if(dPhiShift<-2*TMath::Pi()) dPhiShift += 2*TMath::Pi();
			if(dPhiShift<-0.5*TMath::Pi()) dPhiShift += 2*TMath::Pi();
			if(dPhiShift>1.5*TMath::Pi()) dPhiShift -= 2*TMath::Pi();
			if(isSignal) fhDphiPtSig->Fill(dPhiShift,ptcorr);
			else         fhDphiPtRef->Fill(dPhiShift,ptcorr);

			// selection on relative phi
			if(fJetHadronDeltaPhi>0. &&
					TMath::Abs(dphi)<TMath::Pi()-fJetHadronDeltaPhi) continue;

			if(fCent<10.) fh2RPJetsC10->Fill(TMath::Abs(phiBin), ptcorr);
			if(fCent<20.) fh2RPJetsC20->Fill(TMath::Abs(phiBin), ptcorr);

			Float_t phitt=partback->Phi();
			if(phitt<0)phitt+=TMath::Pi()*2.; 
			Int_t phiBintt = GetPhiBin(phitt-fEPV0);

			Double_t fillspec[] = {fCent,areabig,ptcorr,partback->Pt(),dPhiShift, static_cast<Double_t>(phiBintt)};
			fHJetSpec->Fill(fillspec);
		}
	}
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());
      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoTrackLoop()
{
  AliClusterContainer* clusCont = GetClusterContainer(0);

  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Pt());

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Phi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Eta());

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

        if (clusCont) {
          Int_t iCluster = track->GetEMCALcluster();
          if (iCluster >= 0) {
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            if (cluster) {
              histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
              fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
            }
          }
        }
      }
    }
    sumAcceptedTracks += count;

    histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNTracks";
  fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoClusterLoop()
{
  TString histname;
  TString groupname;
  UInt_t sumAcceptedClusters = 0;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();

    for(auto cluster : clusCont->all()) {
      if (!cluster) continue;

      if (cluster->GetIsExotic()) {
        histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->E());
      }
    }

    UInt_t count = 0;
    for(auto cluster : clusCont->accepted()) {
      if (!cluster) continue;
      count++;

      AliTLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->E());

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Phi_0_2pi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Eta());
    }
    sumAcceptedClusters += count;

    histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNClusters";
  fHistManager.FillTH1(histname, sumAcceptedClusters);
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoCellLoop()
{
  if (!fCaloCells) return;

  TString histname;

  const Short_t ncells = fCaloCells->GetNumberOfCells();

  histname = TString::Format("%s/histNCells_%d", fCaloCellsName.Data(), fCentBin);
  fHistManager.FillTH1(histname, ncells);

  histname = TString::Format("%s/histCellEnergy_%d", fCaloCellsName.Data(), fCentBin);
  for (Short_t pos = 0; pos < ncells; pos++) {
    Double_t amp   = fCaloCells->GetAmplitude(pos);
    
    fHistManager.FillTH1(histname, amp);
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskJetCoreEmcal::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskJetCoreEmcal::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskJetCoreEmcal::Terminate(Option_t *) 
{
}

Int_t  AliAnalysisTaskJetCoreEmcal::SelectTrigger(TList *list,Double_t minT,Double_t maxT,Int_t &number){

	Int_t index=-1;
	Int_t triggers[100];

	for(Int_t cr=0;cr<100;cr++) triggers[cr]=-1;

	Int_t im=0;

	TString groupname = "";
	AliParticleContainer* partCont = GetParticleContainer(0);
	groupname = partCont->GetName();
	UInt_t iCount = 0;
	for(auto part : partCont->accepted()) {
		if (!part) continue;
		list->Add(part);
		iCount++;
		if(part->Pt()>=minT && part->Pt()<maxT){
			triggers[im]=iCount-1;
			im=im+1;
		}
	}
	number=im;
	Int_t rd=0;
	if(im>0) rd=fRandom->Integer(im);
	index=triggers[rd];
	return index;
}

Double_t AliAnalysisTaskJetCoreEmcal::RelativePhi(Double_t mphi,Double_t vphi){

	// get relative DeltaPhi from -pi < DeltaPhi < pi

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}

Int_t AliAnalysisTaskJetCoreEmcal::GetPhiBin(Double_t phi)
{
    Int_t phibin=-1;
    if(!(TMath::Abs(phi)<=2*TMath::Pi())){AliError("phi w.r.t. RP out of defined range");return -1;}
    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));
    phibin=Int_t(fNRPBins*phiwrtrp/(0.5*TMath::Pi()));
    if(phibin<0||phibin>=fNRPBins){AliError("Phi Bin not defined");}
    return phibin;
}
