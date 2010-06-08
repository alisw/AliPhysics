/**************************************************************************
 * Task for Fragmentation Function in PWG4 Jet Task Force Train           *
 **************************************************************************/


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

/* $Id: */


#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TParticlePDG.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
//#include "AliAnalysisHelperJetTasks.h"
#include "AliKFVertex.h"


#include "AliAnalysisTaskFragFuncBB.h"

#include "AliStack.h"


ClassImp(AliAnalysisTaskFragFuncBB)

//________________________________________________________________________

AliAnalysisTaskFragFuncBB::AliAnalysisTaskFragFuncBB()
	: AliAnalysisTaskSE()
	,fESD(0x0)
	,fAOD(0x0)
	,fMCEvent(0x0)
	,fBranchRecJets("jets")
	,fBranchGenJets("")
        ,fLeadingRecJet(-1)
        ,fLeadingGenJet(-1)
	,fTrackTypeRec(0)
	,fTrackTypeGen(0)
	,fFilterMask(0)
	,fUseGlobalSelection(kFALSE)
	,fUseAODJetInput(kFALSE)
	,fUseAODTrackInput(kFALSE)
	,fUseAODMCInput(kFALSE)
	,fJetRadius(.4)
	,fEtaMaxPart(.9)
        ,fEtaMaxJets(.5)
	,fHistList(0x0)
        ,fh1_eventSelection(0x0)
        ,fh1_vertexNContributors(0x0)
        ,fh1_vertexZ(0x0)
	,fh1_recJets_Et(0x0)
	,fh1_genJets_Et(0x0)
	,fh2_recJets_EtaPhi(0x0)
	,fh2_genJets_EtaPhi(0x0)
	,fh1_recJetsWoC_Et(0x0)
	,fh1_genJetsWoC_Et(0x0)
	,fh2_recJetsWoC_EtaPhi(0x0)
	,fh2_genJetsWoC_EtaPhi(0x0)
	,fh1_recLJets_Et(0x0)
	,fh1_genLJets_Et(0x0)
	,fh2_recLJets_EtaPhi(0x0)
	,fh2_genLJets_EtaPhi(0x0)
	,fh2_recFF_JetEt(0x0)
	,fh2_genFF_JetEt(0x0)
	,fh2_recHumpBacked_JetEt(0x0)
	,fh2_genHumpBacked_JetEt(0x0)
	,fh2_recFF_JetWoCEt(0x0)
	,fh2_genFF_JetWoCEt(0x0)
	,fh2_recHumpBacked_JetWoCEt(0x0)
	,fh2_genHumpBacked_JetWoCEt(0x0)
	,fh2_recFF_LJetEt(0x0)
	,fh2_genFF_LJetEt(0x0)
	,fh2_recHumpBacked_LJetEt(0x0)
	,fh2_genHumpBacked_LJetEt(0x0)
	,fh1_recPart_Pt(0x0)
	,fh1_genPart_Pt(0x0)
	,fh2_recPart_EtaPhi(0x0)
	,fh2_genPart_EtaPhi(0x0)
	,fh1_recJetPart_Pt(0x0)
	,fh1_genJetPart_Pt(0x0)
	,fh2_recJetPart_EtaPhi(0x0)
	,fh2_genJetPart_EtaPhi(0x0)
        ,fh2_recJetPart_RJetPt(0x0)
        ,fh2_genJetPart_RJetPt(0x0)
	,fh1_recJetWoCPart_Pt(0x0)
	,fh1_genJetWoCPart_Pt(0x0)
	,fh2_recJetWoCPart_EtaPhi(0x0)
	,fh2_genJetWoCPart_EtaPhi(0x0)
        ,fh2_recJetWoCPart_RJetPt(0x0)
        ,fh2_genJetWoCPart_RJetPt(0x0)
	,fh1_recLJetPart_Pt(0x0)
	,fh1_genLJetPart_Pt(0x0)
	,fh2_recLJetPart_EtaPhi(0x0)
	,fh2_genLJetPart_EtaPhi(0x0)
        ,fh2_recLJetPart_RJetPt(0x0)
        ,fh2_genLJetPart_RJetPt(0x0)
        ,fh3_recPart_EtaPhiPt(0x0)

{
	//
	// Constructor
	//

}

AliAnalysisTaskFragFuncBB::AliAnalysisTaskFragFuncBB(const char *name) 
	: AliAnalysisTaskSE(name)
	,fESD(0x0)
	,fAOD(0x0)
	,fMCEvent(0x0)
	,fBranchRecJets("jets")
	,fBranchGenJets("")
        ,fLeadingRecJet(-1)
        ,fLeadingGenJet(-1)
	,fTrackTypeRec(0)
	,fTrackTypeGen(0)
	,fFilterMask(0)
	,fUseGlobalSelection(kFALSE)
	,fUseAODJetInput(kFALSE)
	,fUseAODTrackInput(kFALSE)
	,fUseAODMCInput(kFALSE)
	,fJetRadius(.4)
	,fEtaMaxPart(.9)
	,fEtaMaxJets(.5)
	,fHistList(0x0)
        ,fh1_eventSelection(0x0)
        ,fh1_vertexNContributors(0x0)
        ,fh1_vertexZ(0x0)
	,fh1_recJets_Et(0x0)
	,fh1_genJets_Et(0x0)
	,fh2_recJets_EtaPhi(0x0)
	,fh2_genJets_EtaPhi(0x0)
	,fh1_recJetsWoC_Et(0x0)
	,fh1_genJetsWoC_Et(0x0)
	,fh2_recJetsWoC_EtaPhi(0x0)
	,fh2_genJetsWoC_EtaPhi(0x0)
	,fh1_recLJets_Et(0x0)
	,fh1_genLJets_Et(0x0)
	,fh2_recLJets_EtaPhi(0x0)
	,fh2_genLJets_EtaPhi(0x0)
	,fh2_recFF_JetEt(0x0)
	,fh2_genFF_JetEt(0x0)
	,fh2_recHumpBacked_JetEt(0x0)
	,fh2_genHumpBacked_JetEt(0x0)
	,fh2_recFF_JetWoCEt(0x0)
	,fh2_genFF_JetWoCEt(0x0)
	,fh2_recHumpBacked_JetWoCEt(0x0)
	,fh2_genHumpBacked_JetWoCEt(0x0)
	,fh2_recFF_LJetEt(0x0)
	,fh2_genFF_LJetEt(0x0)
	,fh2_recHumpBacked_LJetEt(0x0)
	,fh2_genHumpBacked_LJetEt(0x0)
	,fh1_recPart_Pt(0x0)
	,fh1_genPart_Pt(0x0)
	,fh2_recPart_EtaPhi(0x0)
	,fh2_genPart_EtaPhi(0x0)
	,fh1_recJetPart_Pt(0x0)
	,fh1_genJetPart_Pt(0x0)
	,fh2_recJetPart_EtaPhi(0x0)
	,fh2_genJetPart_EtaPhi(0x0)
        ,fh2_recJetPart_RJetPt(0x0)
        ,fh2_genJetPart_RJetPt(0x0)
	,fh1_recJetWoCPart_Pt(0x0)
	,fh1_genJetWoCPart_Pt(0x0)
	,fh2_recJetWoCPart_EtaPhi(0x0)
	,fh2_genJetWoCPart_EtaPhi(0x0)
        ,fh2_recJetWoCPart_RJetPt(0x0)
        ,fh2_genJetWoCPart_RJetPt(0x0)
	,fh1_recLJetPart_Pt(0x0)
	,fh1_genLJetPart_Pt(0x0)
	,fh2_recLJetPart_EtaPhi(0x0)
	,fh2_genLJetPart_EtaPhi(0x0)
        ,fh2_recLJetPart_RJetPt(0x0)
        ,fh2_genLJetPart_RJetPt(0x0)
        ,fh3_recPart_EtaPhiPt(0x0)

{
	//
	// Constructor
	//

	DefineOutput(1,  TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskFragFuncBB::~AliAnalysisTaskFragFuncBB()
{
	  
}

//________________________________________________________________________
void AliAnalysisTaskFragFuncBB::UserCreateOutputObjects()
{

	if(fDebug > 1) Printf("AliAnalysisTaskFragFuncBB::UserCreateOutputObjects()");

	//
	// Create histograms / output container
	//
	
	// set binning limits 
	const Int_t nBinJetEt = 250;
	Double_t binLimitsJetEt[nBinJetEt+1];
	for(Int_t iEt=0; iEt<=nBinJetEt; ++iEt){
		if(iEt==0) binLimitsJetEt[iEt] = 0.0;
		else	   binLimitsJetEt[iEt] = binLimitsJetEt[iEt-1] + 1.0;
	}
	
	const Int_t nBinPartPt = 250;
	Double_t binLimitsPartPt[nBinPartPt+1];
	for(Int_t iPt=0; iPt<=nBinPartPt; ++iPt){
		if(iPt==0) binLimitsPartPt[iPt] = 0.0;
		else	   binLimitsPartPt[iPt] = binLimitsPartPt[iPt-1] + 1.0;
	}
	
	const Int_t nBinPhi = 90;
	Double_t binLimitsPhi[nBinPhi+1];
	for(Int_t iPhi=0; iPhi<=nBinPhi; ++iPhi){
		if(iPhi==0) binLimitsPhi[iPhi] = 0.; //-1.*TMath::Pi();
		else        binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + TMath::Pi()*2/nBinPhi;
	}
	
	const Int_t nBinEta = 40;
	Double_t binLimitsEta[nBinEta+1];
	for(Int_t iEta=0; iEta<=nBinEta; ++iEta){
		if(iEta==0) binLimitsEta[iEta] = -2.;
		else        binLimitsEta[iEta] = binLimitsEta[iEta-1] + 0.1; // upper limit: 2.
	}
	
	const Int_t nBinZ = 25;
	Double_t binLimitsZ[nBinZ+1];
	for(Int_t iZ=0; iZ<=nBinZ; ++iZ){
		if(iZ==0) binLimitsZ[iZ] = 0.;
		else	  binLimitsZ[iZ] = binLimitsZ[iZ-1] + 0.05;  // upper limit: 1.25
	}
	
	const Int_t nBinXi = 100;
	Double_t binLimitsXi[nBinXi+1];
	for(Int_t iXi=0; iXi<=nBinXi; ++iXi){
		if(iXi==0) binLimitsXi[iXi] = 0.;
		else	   binLimitsXi[iXi] = binLimitsXi[iXi-1] + 0.1;  // upper limit: 10.
	}
	

	OpenFile(1);
	fHistList = new TList();

	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	

	// Histograms	
        fh1_eventSelection = new TH1F("fh1_eventSelection", "Event Selection", 4, -0.5, 3.5);
        fh1_vertexNContributors = new TH1F("fh1_vertexNContributors", "Vertex N contributors", 51,-.5, 50.5);
        fh1_vertexZ = new TH1F("fh1_vertexZ", "Vertex z distribution", 200, -19.5, 19.5);

	fh1_recJets_Et = new TH1F("fh1_recJets_Et", "reconstructed jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh1_genJets_Et = new TH1F("fh1_genJets_Et", "generated jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh2_recJets_EtaPhi = new TH2F("fh2_recJets_EtaPhi", "Jet axis - rec. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genJets_EtaPhi = new TH2F("fh2_genJets_EtaPhi", "Jet axis - gen. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);

	fh1_recJetsWoC_Et = new TH1F("fh1_recJetsWoC_Et", "reconstructed jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh1_genJetsWoC_Et = new TH1F("fh1_genJetsWoC_Et", "generated jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh2_recJetsWoC_EtaPhi = new TH2F("fh2_recJetsWoC_EtaPhi", "Jet axis - rec. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genJetsWoC_EtaPhi = new TH2F("fh2_genJetsWoC_EtaPhi", "Jet axis - gen. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);

	fh1_recLJets_Et = new TH1F("fh1_recLJets_Et", "reconstructed jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh1_genLJets_Et = new TH1F("fh1_genLJets_Et", "generated jets;E_{T} [GeV]", nBinJetEt, binLimitsJetEt);
	fh2_recLJets_EtaPhi = new TH2F("fh2_recLJets_EtaPhi", "Jet axis - rec. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genLJets_EtaPhi = new TH2F("fh2_genLJets_EtaPhi", "Jet axis - gen. jets;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	
	fh2_recFF_JetEt = new TH2F("fh2_recFF_JetEt", "Fragmentation Function - reconstructed; z = p_{T} / p_{T}^{Jet}", nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_genFF_JetEt = new TH2F("fh2_genFF_JetEt", "Fragmentation Function - generated; z = p_{T} / p_{T}^{Jet}", nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_recHumpBacked_JetEt = new TH2F("fh2_recHumpBacked_JetEt", "Hump Backed Plateau - reconstructed ; #xi = log( p_{T}^{Jet} / p_{T} )", nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);
	fh2_genHumpBacked_JetEt = new TH2F("fh2_genHumpBacked_JetEt", "Hump Backed Plateau - generated ; #xi = log( p_{T}^{Jet} / p_{T} )", nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);

	fh2_recFF_JetWoCEt         = new TH2F("fh2_recFF_JetWoCEt", "Fragmentation Function - reconstructed; z = p_{T} / p_{T}^{Jet}", 
                                               nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_genFF_JetWoCEt         = new TH2F("fh2_genFF_JetWoCEt", "Fragmentation Function - generated; z = p_{T} / p_{T}^{Jet}", 
                                               nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_recHumpBacked_JetWoCEt = new TH2F("fh2_recHumpBacked_JetWoCEt", "Hump Backed Plateau - reconstructed ; #xi = log( p_{T}^{Jet} / p_{T} )", 
                                               nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);
	fh2_genHumpBacked_JetWoCEt = new TH2F("fh2_genHumpBacked_JetWoCEt", "Hump Backed Plateau - generated ; #xi = log( p_{T}^{Jet} / p_{T} )", 
                                               nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);

	fh2_recFF_LJetEt = new TH2F("fh2_recFF_LJetEt", "Fragmentation Function - reconstructed; z = p_{T} / p_{T}^{Jet}", nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_genFF_LJetEt = new TH2F("fh2_genFF_LJetEt", "Fragmentation Function - generated; z = p_{T} / p_{T}^{Jet}", nBinZ, binLimitsZ, nBinJetEt, binLimitsJetEt);
	fh2_recHumpBacked_LJetEt = new TH2F("fh2_recHumpBacked_LJetEt", "Hump Backed Plateau - reconstructed ; #xi = log( p_{T}^{Jet} / p_{T} )", nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);
	fh2_genHumpBacked_LJetEt = new TH2F("fh2_genHumpBacked_LJetEt", "Hump Backed Plateau - generated ; #xi = log( p_{T}^{Jet} / p_{T} )", nBinXi, binLimitsXi, nBinJetEt, binLimitsJetEt);
	
	fh1_recPart_Pt = new TH1F("fh1_recPart_Pt", "Particles Spectrum - reconstructed;p_{T}", nBinPartPt, binLimitsPartPt);
	fh1_genPart_Pt = new TH1F("fh1_genPart_Pt", "Particles Spectrum - generated;p_{T}", nBinPartPt, binLimitsPartPt);
	fh2_recPart_EtaPhi = new TH2F("fh2_recPart_EtaPhi", "Particles distribution - reconstructed;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genPart_EtaPhi = new TH2F("fh2_genPart_EtaPhi", "Particles distribution - generated;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);

	fh1_recJetPart_Pt = new TH1F("fh1_recJetPart_Pt", "Jet Particles Spectrum - reconstructed;p_{T}", nBinPartPt, binLimitsPartPt);
	fh1_genJetPart_Pt = new TH1F("fh1_genJetPart_Pt", "Jet Particles Spectrum - generated;p_{T}", nBinPartPt, binLimitsPartPt);
	fh2_recJetPart_EtaPhi = new TH2F("fh2_recJetPart_EtaPhi", "Jet Particles - reconstructed;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genJetPart_EtaPhi = new TH2F("fh2_genJetPart_EtaPhi", "Jet Particles - generated;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
        fh2_recJetPart_RJetPt = new TH2F("fh2_recJetPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);
        fh2_genJetPart_RJetPt = new TH2F("fh2_genJetPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);
	
        fh1_recJetWoCPart_Pt = new TH1F("fh1_recJetWoCPart_Pt", "Jet Particles Spectrum - reconstructed;p_{T}", nBinPartPt, binLimitsPartPt);
	fh1_genJetWoCPart_Pt = new TH1F("fh1_genJetWoCPart_Pt", "Jet Particles Spectrum - generated;p_{T}", nBinPartPt, binLimitsPartPt);
	fh2_recJetWoCPart_EtaPhi = new TH2F("fh2_recJetWoCPart_EtaPhi", "Jet Particles - reconstructed;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genJetWoCPart_EtaPhi = new TH2F("fh2_genJetWoCPart_EtaPhi", "Jet Particles - generated;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);

        fh2_recJetWoCPart_RJetPt = new TH2F("fh2_recJetWoCPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);
        fh2_genJetWoCPart_RJetPt = new TH2F("fh2_genJetWoCPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);
	
        fh1_recLJetPart_Pt = new TH1F("fh1_recLJetPart_Pt", "Jet Particles Spectrum - reconstructed;p_{T}", nBinPartPt, binLimitsPartPt);
	fh1_genLJetPart_Pt = new TH1F("fh1_genLJetPart_Pt", "Jet Particles Spectrum - generated;p_{T}", nBinPartPt, binLimitsPartPt);
	fh2_recLJetPart_EtaPhi = new TH2F("fh2_recLJetPart_EtaPhi", "Jet Particles - reconstructed;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);
	fh2_genLJetPart_EtaPhi = new TH2F("fh2_genLJetPart_EtaPhi", "Jet Particles - generated;#eta;#phi", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi);

        fh2_recLJetPart_RJetPt = new TH2F("fh2_recLJetPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);
        fh2_genLJetPart_RJetPt = new TH2F("fh2_genLJetPart_RJetPt", "Jet Particles - distance to jet axis;R", 101, -.05, 1.05, nBinJetEt, binLimitsJetEt);

        fh3_recPart_EtaPhiPt = new TH3F("fh3_recPart_EtaPhiPt", "Particles - reconstruced;#eta;#phi;p_{T} [GeV]", nBinEta, binLimitsEta, nBinPhi, binLimitsPhi, nBinPartPt, binLimitsPartPt);
	
	
	const Int_t saveLevel = 4;
	if(saveLevel>0){
		fHistList->Add(fh1_recJetPart_Pt);
		fHistList->Add(fh2_recJetPart_RJetPt);
		fHistList->Add(fh2_recFF_JetEt);
		fHistList->Add(fh2_recHumpBacked_JetEt);
		
		fHistList->Add(fh1_genJetPart_Pt);
		fHistList->Add(fh2_genJetPart_RJetPt);
		fHistList->Add(fh2_genFF_JetEt);
		fHistList->Add(fh2_genHumpBacked_JetEt);

		fHistList->Add(fh1_recJetWoCPart_Pt);
		fHistList->Add(fh2_recJetWoCPart_RJetPt);
		fHistList->Add(fh2_recFF_JetWoCEt);
		fHistList->Add(fh2_recHumpBacked_JetWoCEt);
		
		fHistList->Add(fh1_genJetWoCPart_Pt);
		fHistList->Add(fh2_genJetWoCPart_RJetPt);
		fHistList->Add(fh2_genFF_JetWoCEt);
		fHistList->Add(fh2_genHumpBacked_JetWoCEt);

		fHistList->Add(fh1_recLJetPart_Pt);
		fHistList->Add(fh2_recLJetPart_RJetPt);
		fHistList->Add(fh2_recFF_LJetEt);
		fHistList->Add(fh2_recHumpBacked_LJetEt);
		
		fHistList->Add(fh1_genLJetPart_Pt);
		fHistList->Add(fh2_genLJetPart_RJetPt);
		fHistList->Add(fh2_genFF_LJetEt);
		fHistList->Add(fh2_genHumpBacked_LJetEt);
	}
	if(saveLevel>1){
                fHistList->Add(fh1_eventSelection);

		fHistList->Add(fh1_recPart_Pt);
		fHistList->Add(fh1_recJets_Et);
		fHistList->Add(fh1_recJetsWoC_Et);
		fHistList->Add(fh1_recLJets_Et);
		
		fHistList->Add(fh1_genPart_Pt);
		fHistList->Add(fh1_genJets_Et);
		fHistList->Add(fh1_genJetsWoC_Et);
		fHistList->Add(fh1_genLJets_Et);
	}
	if(saveLevel>2){
	}
	if(saveLevel>3){
		fHistList->Add(fh2_recJetPart_EtaPhi);
		fHistList->Add(fh2_recJetWoCPart_EtaPhi);
		fHistList->Add(fh2_recLJetPart_EtaPhi);
		fHistList->Add(fh2_recJets_EtaPhi);
		fHistList->Add(fh2_recJetsWoC_EtaPhi);
		fHistList->Add(fh2_recLJets_EtaPhi);
		fHistList->Add(fh2_recPart_EtaPhi);
		
		fHistList->Add(fh2_genJetPart_EtaPhi);
		fHistList->Add(fh2_genJetWoCPart_EtaPhi);
		fHistList->Add(fh2_genLJetPart_EtaPhi);
		fHistList->Add(fh2_genJets_EtaPhi);
		fHistList->Add(fh2_genJetsWoC_EtaPhi);
		fHistList->Add(fh2_genLJets_EtaPhi);
		fHistList->Add(fh2_genPart_EtaPhi);

                fHistList->Add(fh3_recPart_EtaPhiPt);

                fHistList->Add(fh1_vertexNContributors);
                fHistList->Add(fh1_vertexZ);
	}

	// =========== Switch on Sumw2 for all histos ===========
	for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
		TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
		if (h1) h1->Sumw2();	
	}

	TH1::AddDirectory(oldStatus);
}

//_______________________________________________________________________
void AliAnalysisTaskFragFuncBB::Init()
{
	// Initialization
	if(fDebug > 1) Printf("AliAnalysisTaskFragFuncBB::Init()");

}


//________________________________________________________________________
void AliAnalysisTaskFragFuncBB::UserExec(Option_t *) 
{
	// Main loop
	// Called for each event
	if(fDebug > 1) Printf("AliAnalysisTaskFragFuncBB::UserExec()");
	
	/*
	if(!AliAnalysisHelperJetTasks::Selected()&&fUseGlobalSelection){
		// no selection by the service task, we continue
		if (fDebug > 1)Printf("Not selected %s:%d",(char*)__FILE__,__LINE__);
		PostData(1, fHistList);
		return;
	}
	*/

	if(fDebug > 1) Printf("Analysis event #%5d", (Int_t) fEntry);
        // Trigger selection
        
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)
               ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());
        if(inputHandler->IsEventSelected()){
		if(fDebug > 1) AliInfo(" Trigger Selection: event ACCEPTED ... ");
                fh1_eventSelection->Fill(1.);
	} else {
		if (fDebug > 1) AliInfo(" Trigger Selection: event REJECTED ... ");
                fh1_eventSelection->Fill(0.);
	        PostData(1, fHistList);
		return;
	}
  


	fESD = dynamic_cast<AliESDEvent*>(InputEvent());
	if(!fESD){
		Printf("%s:%d ESDEvent not found in the input", (char*)__FILE__,__LINE__);
	}

	fMCEvent = MCEvent();
	if(!fMCEvent){
		Printf("%s:%d MCEvent not found in the input", (char*)__FILE__,__LINE__);
	}


	if(fUseAODJetInput){
		fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
		if(!fAOD){
			Printf("%s:%d AODEvent not found in InputManager %d", (char*)__FILE__,__LINE__, fUseAODJetInput);
			return;
		}
	}
	else {
		// assume that the AOD is in the general output ...
		fAOD = dynamic_cast<AliAODEvent*>(AODEvent());
		if(!fAOD){
			Printf("%s:%d AODEvent not found in the output", (char*)__FILE__,__LINE__);
			return;
		}
	}

	//Event selection (vertex) *****************************************
         
	AliKFVertex primVtx(*(fAOD->GetPrimaryVertex()));
	Int_t nTracksPrim = primVtx.GetNContributors();
        fh1_vertexNContributors->Fill(nTracksPrim);
        fh1_vertexZ->Fill(primVtx.GetZ());
	if (fDebug > 1) AliInfo(Form(" Primary-vertex Selection: %d",nTracksPrim));
	if(!nTracksPrim){
		if (fDebug > 1) AliInfo(" Primary-vertex Selection: event REJECTED ...");
                fh1_eventSelection->Fill(2.);
	        PostData(1, fHistList);
		return;
	}
	if (fDebug > 1) AliInfo(" Primary-vertex Selection: event ACCEPTED ...");
	fh1_eventSelection->Fill(3.);
        



	//___ fetch jets __________________________________________________________________________

	// reconstructed jets
	AliAODJet recJets[fgkMaxJets];
	Int_t nRecJets = 0;

	if(fBranchRecJets.Length()>0){
		TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRecJets.Data()));
		if(aodRecJets){
			Int_t iCount = 0;
			for(Int_t ij=0; ij<aodRecJets->GetEntries(); ++ij){
				if(iCount>=fgkMaxJets){
					Printf("%s:%d max. nb. of jets %d reached, %d", (char*)__FILE__,__LINE__,fgkMaxJets, ij);
					continue;
				}
				AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ij));
				if(!tmp) continue;
				recJets[iCount] = *tmp;

                                // leading jet?
                                if(fLeadingRecJet<0 || recJets[fLeadingRecJet].Pt()<recJets[iCount].Pt()) fLeadingRecJet=iCount;

				iCount++;
			}
			nRecJets = iCount;
		}
		else {
			Printf("%s:%d no reconstructed jet array with name %s in AOD", (char*)__FILE__,__LINE__,fBranchRecJets.Data());
			if(fDebug>1)fAOD->Print();
			return;
		}
	}

	// check that the jets are sorted
	Float_t ptOld = 999999.;
	for(int ir = 0;ir < nRecJets;ir++){
		Float_t tmpPt = recJets[ir].Pt();
		if(tmpPt>ptOld){
			Printf("%s:%d Jets Not Sorted!! %d:%.3E %d%.3E",(char*)__FILE__,__LINE__,ir,tmpPt,ir-1,ptOld);
		}
		ptOld = tmpPt;
	}

	// generated jets
	AliAODJet genJets[fgkMaxJets];
	Int_t nGenJets = 0;
	//Int_t nTrials = 1;  // Trials for MC trigger
	//Double_t ptHard = -1;

	if(fBranchGenJets.Length()==0){ // if no generated jets from AOD used, get pythia jets
		if(!fMCEvent){
			Printf("%s:%d no mcEvent",(char*)__FILE__,__LINE__);
		}

		/*
		   AliGenPythiaEventHeader*  pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(fMCEvent);
		   Int_t iCount = 0;  
		   if(pythiaGenHeader){
	//nTrials = pythiaGenHeader->Trials();
	//ptHard  = pythiaGenHeader->GetPtHard();
	int iProcessType = pythiaGenHeader->ProcessType();
	// 11 f+f -> f+f
	// 12 f+barf -> f+barf
	// 13 f+barf -> g+g
	// 28 f+g -> f+g
	// 53 g+g -> f+barf
	// 68 g+g -> g+g
	if (fDebug > 10)Printf("%d iProcessType %d",__LINE__, iProcessType);
	if(fDebug>20)AliAnalysisHelperJetTasks::PrintStack(fMCEvent);

	// fetch the pythia generated jets
	for(int ip=0; ip<pythiaGenHeader->NTriggerJets(); ++ip){
	if(iCount>=fgkMaxJets){
	Printf("%s:%d max. nb. of jets %d reached, %d", (char*)__FILE__,__LINE__,fgkMaxJets, ip);
	continue;
	}
	Float_t p[4];
	pythiaGenHeader->TriggerJet(ip, p);
	genJets[iCount].SetPxPyPzE(p[0], p[1], p[2], p[3]);

	iCount++;
	}
	nGenJets = iCount;
	}
		 */		
	}
	else {
		TClonesArray *aodGenJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchGenJets.Data()));
		if(aodGenJets){
			Int_t iCount = 0;
			for(Int_t ig=0; ig<aodGenJets->GetEntries(); ++ig){
				if(iCount>=fgkMaxJets){
					Printf("%s:%d max. nb. of jets %d reached, %d", (char*)__FILE__,__LINE__,fgkMaxJets, ig);
					continue;
				}
				AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodGenJets->At(ig));
				if(!tmp) continue;
				genJets[iCount] = *tmp;

                                // check for leading jet
                                if(fLeadingGenJet<0 || genJets[fLeadingGenJet].Pt()<genJets[iCount].Pt()) fLeadingGenJet=iCount;

				iCount++;
			}
			nGenJets = iCount;
		}
		else {
			if(fDebug>0)Printf("%s:%d Generated jet branch %s not found",(char*)__FILE__,__LINE__,fBranchGenJets.Data());
			if(fDebug>1)fAOD->Print();
		}
	}


	//____ fetch particles __________________________________________________________

	TList recParticles;
	TList genParticles;

	Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
        Int_t nRecPart = 0;
        if(nT>=0) nRecPart = recParticles.GetEntries();
	if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,nRecPart);
	
        nT = GetListOfTracks(&genParticles,fTrackTypeGen);
        Int_t nGenPart = 0;
        if(nT>=0) nGenPart = genParticles.GetEntries();
	if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,nGenPart);


	//____ analysis, fill histos ___________________________________________________

	// loop over rec. particles
	for(Int_t it=0; it<nRecPart; ++it){
		AliVParticle *part = (AliVParticle*)recParticles.At(it);		
		fh1_recPart_Pt       -> Fill( part->Pt() );
		fh2_recPart_EtaPhi   -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );
                fh3_recPart_EtaPhiPt -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()), part->Pt() ); 
	}

	// loop over reconstructed jets
	for(Int_t ij=0; ij<nRecJets; ++ij){

		Double_t recJetPt  = recJets[ij].Pt();
		Double_t recJetEta = recJets[ij].Eta();
		Double_t recJetPhi = TVector2::Phi_0_2pi( recJets[ij].Phi() );

                // without cuts
                fh1_recJetsWoC_Et -> Fill( recJetPt );
                fh2_recJetsWoC_EtaPhi -> Fill( recJetEta, recJetPhi );

                // eta acceptance cut
                if(TMath::Abs(recJetEta)<fEtaMaxJets){
		   fh1_recJets_Et -> Fill( recJetPt );
		   fh2_recJets_EtaPhi -> Fill( recJetEta, recJetPhi );

                   // leading jet
                   if(fLeadingRecJet==ij){
		      fh1_recLJets_Et -> Fill( recJetPt );
		      fh2_recLJets_EtaPhi -> Fill( recJetEta, recJetPhi );
                   }
                }

		// loop over rec. particles (in loop over rec. jets)
		for(Int_t it=0; it<nRecPart; ++it){
			AliVParticle *part = (AliVParticle*)recParticles.At(it);

			if(recJets[ij].DeltaR(part)<fJetRadius){

				Double_t z = part->Pt() / recJetPt;
				Double_t xi = 0;
				if(z!=0) xi = TMath::Log(1/z);
                                 
                                // without cuts
				fh2_recFF_JetWoCEt -> Fill( z, recJetPt );
				fh2_recHumpBacked_JetWoCEt -> Fill( xi, recJetPt );

				fh1_recJetWoCPart_Pt -> Fill( part->Pt() );
				fh2_recJetWoCPart_EtaPhi -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );

                                fh2_recJetWoCPart_RJetPt -> Fill( recJets[ij].DeltaR(part), recJetPt );

                                // jet eta acceptance cut
                                if(TMath::Abs(recJetEta)<fEtaMaxJets){
				   fh2_recFF_JetEt -> Fill( z, recJetPt );
				   fh2_recHumpBacked_JetEt -> Fill( xi, recJetPt );

				   fh1_recJetPart_Pt -> Fill( part->Pt() );
				   fh2_recJetPart_EtaPhi -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );

                                   fh2_recJetPart_RJetPt -> Fill( recJets[ij].DeltaR(part), recJetPt );

                                   // leading jet
                                   if(fLeadingRecJet==ij){
				      fh2_recFF_LJetEt -> Fill( z, recJetPt );
				      fh2_recHumpBacked_LJetEt -> Fill( xi, recJetPt );

				      fh1_recLJetPart_Pt -> Fill( part->Pt() );
				      fh2_recLJetPart_EtaPhi -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );

                                      fh2_recLJetPart_RJetPt -> Fill( recJets[ij].DeltaR(part), recJetPt );
                                   }
                                }

			}
		}
	}


	// loop over gen. particles
	for(Int_t it=0; it<nGenPart; ++it){
		AliVParticle *part = (AliVParticle*)genParticles.At(it);
		fh1_genPart_Pt -> Fill( part->Pt() );
		fh2_genPart_EtaPhi -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );
	}


	// loop over gen. jets
	for(Int_t ig=0; ig<nGenJets; ++ig){
		Double_t genJetPt  = genJets[ig].Pt();
		Double_t genJetEta = genJets[ig].Eta();
		Double_t genJetPhi = TVector2::Phi_0_2pi( genJets[ig].Phi() );

		fh1_genJets_Et -> Fill( genJetPt );
		fh2_genJets_EtaPhi -> Fill( genJetEta, genJetPhi );

		// loop over gen. particles (in loop over gen. jets)
		for(Int_t it=0; it<nGenPart; ++it){
			AliVParticle *part = (AliVParticle*)genParticles.At(it);

			if(genJets[ig].DeltaR(part)<fJetRadius){
				//if(jetPart[iJType][iPType][it]) continue;  // part already in other jet-cone
				//jetPart[iJType][iPType][it] = kTRUE;

				Double_t z = part->Pt() / genJetPt;
				Double_t xi = 0;
				if(z!=0) xi = TMath::Log(1/z);

				fh2_genFF_JetEt -> Fill( z, genJetPt );
				fh2_genHumpBacked_JetEt -> Fill( xi, genJetPt );

				fh1_genJetPart_Pt -> Fill( part->Pt() );
				fh2_genJetPart_EtaPhi -> Fill( part->Eta(), TVector2::Phi_0_2pi(part->Phi()) );

                                fh2_genJetPart_RJetPt -> Fill( genJets[ig].DeltaR(part), genJetPt );
			}
		}
	}

	//Post output data.
	PostData(1, fHistList);

}	  

//________________________________________________________________________
void AliAnalysisTaskFragFuncBB::Terminate(Option_t *) 
{
	if(fDebug > 1) printf("AliAnalysisTaskFragFuncBB::Terminate() \n");
}  

//________________________________________________________________________

Int_t AliAnalysisTaskFragFuncBB::GetListOfTracks(TList *list, Int_t type){

	if(fDebug > 2) Printf("%s:%d Selecting tracks with %d", (char*)__FILE__,__LINE__,type);

        if(type==kTrackUndef) return -1;

	Int_t iCount = 0;
	if(type==kTrackAOD){
		// all rec. tracks, esd filter mask, eta range
		AliAODEvent *aod = 0;
		if(fUseAODTrackInput) aod = dynamic_cast<AliAODEvent*>(InputEvent());
		else aod = AODEvent();
		if(!aod) return iCount;

		for(Int_t it=0; it<aod->GetNumberOfTracks(); ++it){
			AliAODTrack *tr = aod->GetTrack(it);
			if((fFilterMask>0)&&!(tr->TestFilterBit(fFilterMask))) continue;
			if(TMath::Abs(tr->Eta())>fEtaMaxPart) continue;
			list->Add(tr);
			iCount++;
		}
	}
	else if (type==kTrackKineAll || type==kTrackKineCharged){
		// kine particles, all or rather charged
		if(!fMCEvent) return iCount;

		for(Int_t it=0; it<fMCEvent->GetNumberOfTracks(); ++it){
			AliMCParticle* part = (AliMCParticle*) fMCEvent->GetTrack(it);

			if(type == kTrackKineCharged){
				if(part->Particle()->GetPDG()->Charge()==0) continue;
			}

			list->Add(part);
			iCount++;
		}
	}
	else if (type==kTrackAODMCCharged || type==kTrackAODMCAll || type==kTrackAODMCChargedAcceptance) {
		// MC particles (from AOD), physical primaries, all or rather charged or rather charged within acceptance
		AliAODEvent *aod = 0;
		if(fUseAODMCInput) aod = dynamic_cast<AliAODEvent*>(InputEvent());
		else aod = AODEvent();
		if(!aod) return iCount;

		TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
		if(!tca)return iCount;

		for(int it=0; it<tca->GetEntriesFast(); ++it){
			AliAODMCParticle *part = dynamic_cast<AliAODMCParticle*>(tca->At(it));
			if(!part->IsPhysicalPrimary())continue;

			if (type==kTrackAODMCCharged || type==kTrackAODMCChargedAcceptance){
				if(part->Charge()==0) continue;
				if(type==kTrackAODMCChargedAcceptance && TMath::Abs(part->Eta())>fEtaMaxPart) continue;
			}

			list->Add(part);
			iCount++;
		}
	}

	list->Sort();
	return iCount;

}
