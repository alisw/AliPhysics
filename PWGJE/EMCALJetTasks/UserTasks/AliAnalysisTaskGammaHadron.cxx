//
// Task to estimate the number of gamma-hadron
// statistic available in the Pb+Pb run.
//
// Author: E. Epple

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"

#include "AliAnalysisTaskGammaHadron.h"

ClassImp(AliAnalysisTaskGammaHadron)

//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskGammaHadron", kTRUE),
fCellEnergyCut(0.05),
fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentMethod2(""),
fCentMethod3(""),
fDoV0QA(0),
fDoEPQA(0),
fDoLeadingObjectPosition(0),
fMaxCellsInCluster(50),
fCent2(0),
fCent3(0),
fVZERO(0),
fV0ATotMult(0),
fV0CTotMult(0),


fHistNoClus_pt(0),
fHistNoClus_ptH(0),
fHistNoClus_ptLeadH(0),
fHistNoClus_Leadpt(0),
fHistNoClus_LeadptH(0),
fHistNoClus_LeadptLeadH(0),

fHistNoClus_xEH(0),
fHistNoClus_LeadxEH(0),
fHistNoClus_xELeadH(0),
fHistNoClus_LeadxELeadH(0)

{
	// Default constructor.

	fAODfilterBits[0] = 0;
	fAODfilterBits[1] = 0;

	SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(const char *name) :
		  AliAnalysisTaskEmcalJet(name, kTRUE),
		  fCellEnergyCut(0.05),
		  fParticleLevel(kFALSE),
		  fIsMC(kFALSE),
		  fCentMethod2(""),
		  fCentMethod3(""),
		  fDoV0QA(0),
		  fDoEPQA(0),
		  fDoLeadingObjectPosition(0),
		  fMaxCellsInCluster(50),
		  fCent2(0),
		  fCent3(0),
		  fVZERO(0),
		  fV0ATotMult(0),
		  fV0CTotMult(0),


		  fHistNoClus_pt(0),
		  fHistNoClus_ptH(0),
		  fHistNoClus_ptLeadH(0),
		  fHistNoClus_Leadpt(0),
		  fHistNoClus_LeadptH(0),
		  fHistNoClus_LeadptLeadH(0),

		  fHistNoClus_xEH(0),
		  fHistNoClus_LeadxEH(0),
		  fHistNoClus_xELeadH(0),
		  fHistNoClus_LeadxELeadH(0)

{
	// Standard

	fAODfilterBits[0] = 0;
	fAODfilterBits[1] = 0;

	SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskGammaHadron::~AliAnalysisTaskGammaHadron()
{
	// Destructor
}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::AllocateHistogramArrays()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::UserCreateOutputObjects()
{
	// Create histograms

	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

	AllocateHistogramArrays();

	TString histname;


	Int_t dim = 0;
	TString title[20];
	Int_t nbins[20] = {0};
	Double_t min[20] = {0};
	Double_t max[20] = {0};

	if (fForceBeamType != AliAnalysisTaskEmcal::kpp)
	{
		title[dim] = "Centrality %";
		nbins[dim] = 101;
		min[dim] = 0;
		max[dim] = 101;
		dim++;

		if (!fCentMethod2.IsNull())
		{
			title[dim] = Form("Centrality %s %%", fCentMethod2.Data());
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = 101;
			dim++;
		}

		if (!fCentMethod3.IsNull())
		{
			title[dim] = Form("Centrality %s %%", fCentMethod3.Data());
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = 101;
			dim++;
		}

		if (fDoV0QA==1)
		{
			title[dim] = "V0A total multiplicity";
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 20000;
			dim++;

			title[dim] = "V0C total multiplicity";
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 20000;
			dim++;
		}
		else if (fDoV0QA==2)
		{
			title[dim] = "V0A+V0C total multiplicity";
			nbins[dim] = 300;
			min[dim] = 0;
			max[dim] = 30000;
			dim++;
		}

		if (!fRhoName.IsNull())
		{
			title[dim] = "#rho (GeV/c)";
			nbins[dim] = fNbins*4;
			min[dim] = 0;
			max[dim] = 400;
			dim++;
		}

		if (fDoEPQA)
		{
			title[dim] = "#psi_{RP}";
			nbins[dim] = 200;
			min[dim] = -TMath::Pi();
			max[dim] = TMath::Pi();
			dim++;
		}
	}

	if (fParticleCollArray.GetEntriesFast()>0)
	{
		title[dim] = "No. of tracks";
		if (fForceBeamType != AliAnalysisTaskEmcal::kpp)
		{
			nbins[dim] = 6000;
			min[dim] = 0;
			max[dim] = 3000;
		}
		else
		{
			nbins[dim] = 200;
			min[dim] = 0;
			max[dim] = 200;
		}
		dim++;

		title[dim] = "p_{T,track}^{leading} (GeV/c)";
		nbins[dim] = fNbins;
		min[dim] = fMinBinPt;
		max[dim] = fMaxBinPt;
		dim++;

		if (fDoLeadingObjectPosition)
		{
			title[dim] = "#eta_{track}^{leading}";
			nbins[dim] = 100;
			min[dim] = -1;
			max[dim] = 1;
			dim++;

			title[dim] = "#phi_{track}^{leading}";
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = TMath::Pi() * 2.02;
			dim++;
		}
	}

	if (fClusterCollArray.GetEntriesFast()>0)
	{
		title[dim] = "No. of clusters";
		nbins[dim] = 2000;
		min[dim] = 0;
		max[dim] = 4000-0.5;
		dim++;

		title[dim] = "E_{cluster}^{leading} (GeV)";
		nbins[dim] = fNbins;
		min[dim] = fMinBinPt;
		max[dim] = fMaxBinPt;
		dim++;

		if (fDoLeadingObjectPosition)
		{
			title[dim] = "#eta_{cluster}^{leading}";
			nbins[dim] = 100;
			min[dim] = -1;
			max[dim] = 1;
			dim++;

			title[dim] = "#phi_{cluster}^{leading}";
			nbins[dim] = 101;
			min[dim] = 0;
			max[dim] = TMath::Pi() * 2.02;
			dim++;
		}
	}

	if (!fCaloCellsName.IsNull())
	{
		title[dim] = "No. of cells";
		nbins[dim] = 3000;
		min[dim] = 0;
		max[dim] = 6000-0.5;
		dim++;
	}


	// - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
	//
	//Room for more fancyness: only one inclusive centrality bin at the moment
	//                         is there a hadron on the other side?
	//                         do I have to check for emcal or phos?

	// Create histograms
	// all clusters as a function of p_T^{Cluster}
	fHistNoClus_pt = new TH1F(Form("fHistNoClus_pt_%0d",1),Form("fHistNoClus_pt_%0d",1), 25, 0, 25);
	fHistNoClus_pt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_pt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClus_pt->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_pt);

	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_Leadpt = new TH1F(Form("fHistNoClus_Leadpt_%0d",1),Form("fHistNoClus_Leadpt_%0d",1), 25, 0, 25);
	fHistNoClus_Leadpt->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_Leadpt->GetYaxis()->SetTitle(Form("No. of lead. Clus. [counts/%0.1f GeV/c]",fHistNoClus_Leadpt->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_Leadpt);

	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_ptH = new TH1F(Form("fHistNoClus_ptH_%0d",1),Form("fHistNoClus_ptH_%0d",1), 25, 0, 25);
	fHistNoClus_ptH->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptH->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_ptH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadptH = new TH1F(Form("fHistNoClus_LeadptH_%0d",1),Form("fHistNoClus_LeadptH_%0d",1), 25, 0, 25);
	fHistNoClus_LeadptH->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_LeadptH->GetYaxis()->SetTitle(Form("No. of lead. Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_LeadptH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadptH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_ptLeadH = new TH1F(Form("fHistNoClus_ptLeadH_%0d",1),Form("fHistNoClus_ptLeadH_%0d",1), 25, 0, 25);
	fHistNoClus_ptLeadH->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptLeadH->GetYaxis()->SetTitle(Form("No. of Clus. with lead. h [counts/%0.1f GeV/c]",fHistNoClus_ptLeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_ptLeadH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadptLeadH = new TH1F(Form("fHistNoClus_LeadptLeadH_%0d",1),Form("fHistNoClus_LeadptLeadH_%0d",1), 25, 0, 25);
	fHistNoClus_LeadptLeadH->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_LeadptLeadH->GetYaxis()->SetTitle(Form("No. of lead. Clus. with lead. h [counts/%0.1f GeV/c]",fHistNoClus_LeadptLeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadptLeadH);


	//
	//   x_E
	//
	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_xEH = new TH1F(Form("fHistNoClus_xEH_%0d",1),Form("fHistNoClus_xEH_%0d",1), 20, 0, 2);
	fHistNoClus_xEH->GetXaxis()->SetTitle("x_{E} (Cluster - hadron)");
	fHistNoClus_xEH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_xEH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_xEH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadxEH = new TH1F(Form("fHistNoClus_LeadxEH_%0d",1),Form("fHistNoClus_LeadxEH_%0d",1), 20, 0, 2);
	fHistNoClus_LeadxEH->GetXaxis()->SetTitle("x_{E} (Lead. Cluster - hadron)");
	fHistNoClus_LeadxEH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_LeadxEH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadxEH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_xELeadH = new TH1F(Form("fHistNoClus_xELeadH_%0d",1),Form("fHistNoClus_xELeadH_%0d",1), 20, 0, 2);
	fHistNoClus_xELeadH->GetXaxis()->SetTitle("x_{E} (Cluster - Lead. hadron)");
	fHistNoClus_xELeadH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_xELeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_xELeadH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadxELeadH = new TH1F(Form("fHistNoClus_LeadxELeadH_%0d",1),Form("fHistNoClus_LeadxELeadH_%0d",1), 20, 0, 2);
	fHistNoClus_LeadxELeadH->GetXaxis()->SetTitle("x_{E} (Lead. Cluster - Lead. hadron)");
	fHistNoClus_LeadxELeadH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_LeadxELeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadxELeadH);


	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::ExecOnce()
{
	AliAnalysisTaskEmcalJet::ExecOnce();

	if (fDoV0QA) {
		fVZERO = InputEvent()->GetVZEROData();
		if (!fVZERO) {
			AliError("AliVVZERO not available");
		}
	}
}

//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::RetrieveEventObjects()
{
	// Retrieve event objects.

	if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
		return kFALSE;

	if (!fCentMethod2.IsNull() || !fCentMethod3.IsNull()) {
		if (fBeamType == kAA || fBeamType == kpA ) {
			AliCentrality *aliCent = InputEvent()->GetCentrality();
			if (aliCent) {
				if (!fCentMethod2.IsNull())
					fCent2 = aliCent->GetCentralityPercentile(fCentMethod2);
				if (!fCentMethod3.IsNull())
					fCent3 = aliCent->GetCentralityPercentile(fCentMethod3);
			}
		}
	}

	if (fVZERO) {
		fV0ATotMult = AliESDUtils::GetCorrV0A(fVZERO->GetMTotV0A(),fVertex[2]);
		fV0CTotMult = AliESDUtils::GetCorrV0C(fVZERO->GetMTotV0C(),fVertex[2]);
	}

	return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	// Fill histograms.

	Float_t trackSum = 0;
	Float_t clusSum = 0;
	Float_t cellSum = 0;

	Int_t ntracks = 0;
	Int_t nclusters = 0;
	Int_t ncells = 0;

	Float_t leadingClusE = 0;
	Float_t leadingClusEta = 0;
	Float_t leadingClusPhi = 0;

	Float_t leadingTrackPt = 0;
	Float_t leadingTrackEta = 0;
	Float_t leadingTrackPhi = 0;


	if (fTracks)
	{
		AliVParticle *leadingTrack = 0;

		ntracks = DoTrackLoop(trackSum, leadingTrack);
		AliDebug(2,Form("%d tracks found in the event", ntracks));

		if (leadingTrack)
		{
			leadingTrackPt = leadingTrack->Pt();
			leadingTrackEta = leadingTrack->Eta();
			leadingTrackPhi = leadingTrack->Phi();
		}
	}

	if (fCaloClusters)
	{
		AliVCluster  *leadingClus = 0;

		nclusters = DoClusterLoop(clusSum, leadingClus);
		AliDebug(2,Form("%d clusters found in the event", nclusters));

		if (leadingClus)
		{
			TLorentzVector leadingClusVect;
			leadingClus->GetMomentum(leadingClusVect, fVertex);
			leadingClusE = leadingClus->E();
			leadingClusEta = leadingClusVect.Eta();
			leadingClusPhi = leadingClusVect.Phi();
		}
	}

	if (fCaloCells)
	{
		ncells = DoCellLoop(cellSum);
		AliDebug(2,Form("%d cells found in the event", ncells));
	}

	return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::DoCellLoop(Float_t &sum)
{
	// Do cell loop.
	// Counts the number of fired cells

	AliVCaloCells *cells = InputEvent()->GetEMCALCells();

	if (!cells)
		return 0;

	const Int_t ncells = cells->GetNumberOfCells();
	Int_t nAccCells = 0;

	for (Int_t pos = 0; pos < ncells; pos++) {
		Float_t amp   = cells->GetAmplitude(pos);
		Int_t   absId = cells->GetCellNumber(pos);

		if (amp < fCellEnergyCut)
			continue;

	//	fHistCellsAbsIdEnergy[fCentBin]->Fill(absId,amp);
		nAccCells++;
		sum += amp;
	}

	return nAccCells;
}

//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::GetFcross(AliVCluster *cluster, AliVCaloCells *cells)
{
	Int_t    AbsIdseed  = -1;
	Double_t Eseed      = 0;
	for (Int_t i = 0; i < cluster->GetNCells(); i++) {
		if (cells->GetCellAmplitude(cluster->GetCellAbsId(i)) > AbsIdseed) {
			Eseed     = cells->GetCellAmplitude(cluster->GetCellAbsId(i));
			AbsIdseed = cluster->GetCellAbsId(i);
		}
	}

	if (Eseed < 1e-9)
		return 100;

	Int_t imod = -1, iphi =-1, ieta=-1,iTower = -1, iIphi = -1, iIeta = -1;
	fGeom->GetCellIndex(AbsIdseed,imod,iTower,iIphi,iIeta);
	fGeom->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,iphi,ieta);

	//Get close cells index and energy, not in corners

	Int_t absID1 = -1;
	Int_t absID2 = -1;

	if (iphi < AliEMCALGeoParams::fgkEMCALRows-1) absID1 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi+1, ieta);
	if (iphi > 0)                                 absID2 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi-1, ieta);

	// In case of cell in eta = 0 border, depending on SM shift the cross cell index

	Int_t absID3 = -1;
	Int_t absID4 = -1;

	if (ieta == AliEMCALGeoParams::fgkEMCALCols-1 && !(imod%2)) {
		absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod+1, iphi, 0);
		absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta-1);
	}
	else if (ieta == 0 && imod%2) {
		absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod,   iphi, ieta+1);
		absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod-1, iphi, AliEMCALGeoParams::fgkEMCALCols-1);
	}
	else {
		if (ieta < AliEMCALGeoParams::fgkEMCALCols-1)
			absID3 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta+1);
		if (ieta > 0)
			absID4 = fGeom->GetAbsCellIdFromCellIndexes(imod, iphi, ieta-1);
	}

	Double_t  ecell1 = cells->GetCellAmplitude(absID1);
	Double_t  ecell2 = cells->GetCellAmplitude(absID2);
	Double_t  ecell3 = cells->GetCellAmplitude(absID3);
	Double_t  ecell4 = cells->GetCellAmplitude(absID4);

	Double_t Ecross = ecell1 + ecell2 + ecell3 + ecell4;

	Double_t Fcross = 1 - Ecross/Eseed;

	return Fcross;
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::DoClusterLoop(Float_t &sum, AliVCluster* &leading)
{
	// Do cluster loop.

	AliClusterContainer* clusters = GetClusterContainer(0);
	if (!clusters) return 0;

	Int_t nAccClusters = 0;

	AliVCaloCells *cells = InputEvent()->GetEMCALCells();

	sum = 0;
	leading = 0;

	// Cluster loop

	AliVCluster* cluster = 0;
	clusters->ResetCurrentID();
	while ((cluster = clusters->GetNextAcceptCluster()))
	{
		sum += cluster->E();

		if (!leading || leading->E() < cluster->E()) leading = cluster;
		TLorentzVector nPart;
		cluster->GetMomentum(nPart, fVertex);

		//all clusters
		fHistNoClus_pt->Fill(nPart.Pt()); //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong


		//check whether there is a hadron in the opposite hemisphere of the gamma
		AliParticleContainer* tracks = GetParticleContainer(0);
		if (!tracks) std::cout<<"ELI: something wrong here"<<std::endl;
		AliVParticle *leadingTrack = 0;

		tracks->ResetCurrentID();
		AliVParticle* track = 0;
		Double_t deltaPhi=0;
		while ((track = tracks->GetNextAcceptParticle()))
		{
			if (!leadingTrack || leadingTrack->Pt() < leadingTrack->Pt()) leadingTrack = track;

			//std::cout<<"nPart.Phi(): "<<nPart.Phi()<<", track->Phi(): "<<track->Phi()<<std::endl;
			deltaPhi=fabs(nPart.Phi()-track->Phi());
			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and any hadron
				if(track->Pt()>0.5)fHistNoClus_ptH->Fill(nPart.Pt());    //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong
				if(nPart.Pt()>5 && track->Pt()>0.5)fHistNoClus_xEH->Fill(-1.0*cos(deltaPhi)*track->Pt()/nPart.Pt());
			}
		}
		//end of correlated hadron check!
		if(leadingTrack)
		{
			deltaPhi=fabs(nPart.Phi()-leadingTrack->Phi());
			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and leadig hadron
				if(leadingTrack->Pt()>0.5)fHistNoClus_ptLeadH->Fill(nPart.Pt());    //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong
				if(nPart.Pt()>5 && leadingTrack->Pt()>0.5)fHistNoClus_xELeadH->Fill(-1.0*cos(deltaPhi)*leadingTrack->Pt()/nPart.Pt());
			}
		}

		//fHistClusPhiEtaEnergy[fCentBin]->Fill(nPart.Eta(), nPart.Phi(), cluster->E());

		Double_t ep = nPart.Phi() - fEPV0;
		while (ep < 0) ep += TMath::Pi();
		while (ep >= TMath::Pi()) ep -= TMath::Pi();

		nAccClusters++;
	}
	if(leading)
	{
		TLorentzVector leadingClusVect;
		leading->GetMomentum(leadingClusVect, fVertex);

		//leading cluster
		fHistNoClus_Leadpt->Fill(leadingClusVect.Pt()); //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong

	    //check whether there is a hadron in the opposite hemisphere of the gamma
		AliParticleContainer* tracks = GetParticleContainer(0);
		if (!tracks) std::cout<<"ELI: something wrong here"<<std::endl;
		AliVParticle *leadingTrack = 0;

		tracks->ResetCurrentID();
		AliVParticle* track = 0;
		Double_t deltaPhi=0;
		while ((track = tracks->GetNextAcceptParticle()))
		{
			if (!leadingTrack || leadingTrack->Pt() < leadingTrack->Pt()) leadingTrack = track;

			deltaPhi=fabs(leadingClusVect.Phi()-track->Phi());
			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and any hadron
				if(track->Pt()>0.5)fHistNoClus_LeadptH    ->Fill(leadingClusVect.Pt());    //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong
				if(leadingClusVect.Pt()>5 && track->Pt()>0.5)fHistNoClus_LeadxEH    ->Fill(-1.0*cos(deltaPhi)*track->Pt()/leadingClusVect.Pt());
			}
		}
		if(leadingTrack)
		{
			deltaPhi=fabs(leadingClusVect.Phi()-leadingTrack->Phi());
			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and leadig hadron
				if(leadingTrack->Pt()>0.5)fHistNoClus_LeadptLeadH->Fill(leadingClusVect.Pt());    //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong
				if(leadingClusVect.Pt()>5 && leadingTrack->Pt()>0.5)fHistNoClus_LeadxELeadH->Fill(-1.0*cos(deltaPhi)*leadingTrack->Pt()/leadingClusVect.Pt());
			}
		}
	}
	//else std::cout<<"No leading cluster found"<<std::endl;

	return nAccClusters;
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::DoTrackLoop(Float_t &sum, AliVParticle* &leading)
{
	// Do track loop.
	AliParticleContainer* tracks = GetParticleContainer(0);

	if (!tracks) return 0;

	Int_t nAccTracks = 0;

	sum = 0;
	leading = 0;

	Int_t neg = 0;
	Int_t zero = 0;

	tracks->ResetCurrentID();
	AliVParticle* track = 0;
	while ((track = tracks->GetNextAcceptParticle()))
	{
		nAccTracks++;

		sum += track->P();

		if (!leading || leading->Pt() < track->Pt()) leading = track;

		if (fParticleLevel)
		{
			//fHistTrPhiEtaPt[fCentBin][0]->Fill(track->Eta(), track->Phi(), track->Pt());
		}
		else
		{
			if (track->GetLabel() == 0)
			{
				zero++;
				/*if (fHistTrPhiEtaZeroLab[fCentBin])
				{
					fHistTrPhiEtaZeroLab[fCentBin]->Fill(track->Eta(), track->Phi());
					fHistTrPtZeroLab[fCentBin]->Fill(track->Pt());
				}*/
			}

			if (track->GetLabel() < 0)
			{
				neg++;
			}

			Int_t type = 0;

			if (tracks->GetClassName() == "AliPicoTrack")
			{
				type = static_cast<AliPicoTrack*>(track)->GetTrackType();
			}
			else if (tracks->GetClassName() == "AliAODTrack")
			{
				if (fAODfilterBits[0] != 0 || fAODfilterBits[1] != 0)
				{
					type = AliPicoTrack::GetTrackType(static_cast<AliAODTrack*>(track), fAODfilterBits[0], fAODfilterBits[1]);
				}
				else
				{
					type = AliPicoTrack::GetTrackType(static_cast<AliVTrack*>(track));
				}
			}

			if (type >= 0 && type <= 3)
			{
			//	fHistTrPhiEtaPt[fCentBin][type]->Fill(track->Eta(), track->Phi(), track->Pt());
			}
			else
			{
				AliDebug(2,Form("%s: track type %d not recognized!", GetName(), type));
			}

			AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track);
			if (!vtrack) continue;

		/*	if ((vtrack->GetTrackEtaOnEMCal() == -999 || vtrack->GetTrackPhiOnEMCal() == -999) && fHistTrPhiEtaNonProp[fCentBin]) {
				fHistTrPhiEtaNonProp[fCentBin]->Fill(vtrack->Eta(), vtrack->Phi());
				fHistTrPtNonProp[fCentBin]->Fill(vtrack->Pt());
			}*/
/*
			if (fHistTrEmcPhiEta[fCentBin])
				fHistTrEmcPhiEta[fCentBin]->Fill(vtrack->GetTrackEtaOnEMCal(), vtrack->GetTrackPhiOnEMCal());
			if (fHistTrEmcPt[fCentBin])
				fHistTrEmcPt[fCentBin]->Fill(vtrack->GetTrackPtOnEMCal());
			if (fHistDeltaEtaPt[fCentBin])
				fHistDeltaEtaPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Eta() - vtrack->GetTrackEtaOnEMCal());
			if (fHistDeltaPhiPt[fCentBin])
				fHistDeltaPhiPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Phi() - vtrack->GetTrackPhiOnEMCal());
			if (fHistDeltaPtvsPt[fCentBin])
				fHistDeltaPtvsPt[fCentBin]->Fill(vtrack->Pt(), vtrack->Pt() - vtrack->GetTrackPtOnEMCal());
		*/
		}
	}

	/*if (fHistTrNegativeLabels[fCentBin]) {
		fHistTrNegativeLabels[fCentBin]->Fill(1. * neg / nAccTracks);
	}

	if (fHistTrZeroLabels[fCentBin]) {
		fHistTrZeroLabels[fCentBin]->Fill(1. * zero / nAccTracks);
	}*/

	return nAccTracks;
}


