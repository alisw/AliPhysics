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
//#include "AliCaloTrackReader.h"

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
fHistNoClus_pt_tr(0),
fHistNoClus_ptH(0),
fHistNoClus_ptH_tr(0),
fHistNoClus_ptLeadH(0),
fHistNoClus_Leadpt(0),
fHistNoClus_LeadptH(0),
fHistNoClus_LeadptLeadH(0),

fHistNoClus_xEH(0),
fHistNoClus_LeadxEH(0),
fHistNoClus_xELeadH(0),
fHistNoClus_LeadxELeadH(0),
fHistpt_assHadron(0),
fHistpt_assHadron_tr(0)

{
	// Default constructor.

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
void AliAnalysisTaskGammaHadron::UserCreateOutputObjects()
{
	// Create histograms
	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

	Int_t nbins[5] = {0};
	Double_t min[5] = {0};
	Double_t max[5] = {0};


	nbins[0] = 31;
	min[0] = 0;
	max[0] = 31;

	nbins[1] = 20;
	min[1] = 0;
	max[1] = 2;

	nbins[2] = 60;  //do 1/2 GeV bins so that you can see the 0.5 cut to set as aminimum pT to combine hadron and gamma
	min[2] = 0;
	max[2] = 30;

	// - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . - . -
	//
	//Room for more fancyness: only one inclusive centrality bin at the moment
	//                         is there a hadron on the other side?
	//                         do I have to check for emcal or phos?

	// Create histograms
	// all clusters as a function of p_T^{Cluster}
	fHistNoClus_pt = new TH1F(Form("fHistNoClus_pt_%0d",1),Form("fHistNoClus_pt_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_pt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_pt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClus_pt->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_pt);

	//trigger clusters as a function of p_T^{Cluster}
	fHistNoClus_pt_tr = new TH1F(Form("fHistNoTrClus_pt_tr_%0d",1),Form("fHistNoTrClus_pt_tr_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_pt_tr->GetXaxis()->SetTitle("p_{T}^{trig. #gamma}");
	fHistNoClus_pt_tr->GetYaxis()->SetTitle(Form("No. of trig. #gamma [counts/%0.1f GeV/c]",fHistNoClus_pt_tr->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_pt_tr);

	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_Leadpt = new TH1F(Form("fHistNoClus_Leadpt_%0d",1),Form("fHistNoClus_Leadpt_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_Leadpt->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_Leadpt->GetYaxis()->SetTitle(Form("No. of lead. Clus. [counts/%0.1f GeV/c]",fHistNoClus_Leadpt->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_Leadpt);

	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_ptH = new TH1F(Form("fHistNoClus_ptH_%0d",1),Form("fHistNoClus_ptH_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_ptH->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptH->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_ptH);

	//trigger clusters as a function of p_T^{Cluster}
	fHistNoClus_ptH_tr = new TH1F(Form("fHistNoClus_ptH_tr_%0d",1),Form("fHistNoClus_ptH_tr_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_ptH_tr->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptH_tr->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH_tr->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_ptH_tr);

	//leading clusters with h as a function of p_T^{Cluster}
	fHistNoClus_LeadptH = new TH1F(Form("fHistNoClus_LeadptH_%0d",1),Form("fHistNoClus_LeadptH_%0d",1),nbins[0], min[0], max[0]);
	fHistNoClus_LeadptH->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_LeadptH->GetYaxis()->SetTitle(Form("No. of lead. Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_LeadptH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadptH);


	//clusters with leding h as a function of p_T^{Cluster}
	fHistNoClus_ptLeadH = new TH1F(Form("fHistNoClus_ptLeadH_%0d",1),Form("fHistNoClus_ptLeadH_%0d",1),nbins[0], min[0], max[0]);
	fHistNoClus_ptLeadH->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptLeadH->GetYaxis()->SetTitle(Form("No. of Clus. with lead. h [counts/%0.1f GeV/c]",fHistNoClus_ptLeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_ptLeadH);


	//leading clusters with leadng hadron as a function of p_T^{Cluster}
	fHistNoClus_LeadptLeadH = new TH1F(Form("fHistNoClus_LeadptLeadH_%0d",1),Form("fHistNoClus_LeadptLeadH_%0d",1),nbins[0], min[0], max[0]);
	fHistNoClus_LeadptLeadH->GetXaxis()->SetTitle("p_{T}^{Leading Calo Cluster}");
	fHistNoClus_LeadptLeadH->GetYaxis()->SetTitle(Form("No. of lead. Clus. with lead. h [counts/%0.1f GeV/c]",fHistNoClus_LeadptLeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadptLeadH);


	//
	//   x_E
	//
	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_xEH = new TH1F(Form("fHistNoClus_xEH_%0d",1),Form("fHistNoClus_xEH_%0d",1), nbins[1], min[1], max[1]);
	fHistNoClus_xEH->GetXaxis()->SetTitle("x_{E} (Cluster - hadron)");
	fHistNoClus_xEH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_xEH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_xEH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadxEH = new TH1F(Form("fHistNoClus_LeadxEH_%0d",1),Form("fHistNoClus_LeadxEH_%0d",1), nbins[1], min[1], max[1]);
	fHistNoClus_LeadxEH->GetXaxis()->SetTitle("x_{E} (Lead. Cluster - hadron)");
	fHistNoClus_LeadxEH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_LeadxEH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadxEH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_xELeadH = new TH1F(Form("fHistNoClus_xELeadH_%0d",1),Form("fHistNoClus_xELeadH_%0d",1), nbins[1], min[1], max[1]);
	fHistNoClus_xELeadH->GetXaxis()->SetTitle("x_{E} (Cluster - Lead. hadron)");
	fHistNoClus_xELeadH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_xELeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_xELeadH);


	//leading clusters as a function of p_T^{Cluster}
	fHistNoClus_LeadxELeadH = new TH1F(Form("fHistNoClus_LeadxELeadH_%0d",1),Form("fHistNoClus_LeadxELeadH_%0d",1), nbins[1], min[1], max[1]);
	fHistNoClus_LeadxELeadH->GetXaxis()->SetTitle("x_{E} (Lead. Cluster - Lead. hadron)");
	fHistNoClus_LeadxELeadH->GetYaxis()->SetTitle(Form("dN/dx_{E} [counts/%0.1f]",fHistNoClus_LeadxELeadH->GetBinWidth(0)));
	fOutput->Add(fHistNoClus_LeadxELeadH);


	//Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
	//Later you can summ the histograms according to the expected statistic!

	fHistpt_assHadron    = new TH1*[nbins[0]]; //make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHistpt_assHadron_tr = new TH1*[nbins[0]]; // make a p_t histogram of the associated hadron for each p_T bin of the trigger-gamma
	//
	for(Int_t i=1; i<nbins[0]+1; i++)
	{

		change that so that the first histogram is
		also filled. the I_nt =1 was probably done for the binwidth/center function ??
		//check whether the max is the
		//same as the histogram size
		Double_t BinWidth = fHistNoClus_ptH->GetBinWidth(i);
		Double_t BinValStart = fHistNoClus_ptH->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;
		//cout<<"BinStart: "<<fHistNoClus_ptH->GetBinCenter(i)<<", width "<<BinWidth<<endl;

		fHistpt_assHadron[i] = new TH1F(Form("fHistAssHadron_pt_%0d",i),Form("fHistAssHadron_pt_%0d",i), nbins[2], min[2], max[2]);
		fHistpt_assHadron[i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
		fHistpt_assHadron[i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistpt_assHadron[i]->GetBinWidth(1)));
		fOutput->Add(fHistpt_assHadron[i]);

		fHistpt_assHadron_tr[i] = new TH1F(Form("fHistAssHadron_pt_tr_%0d",i),Form("fHistAssHadron_pt_tr_%0d",i), nbins[2], min[2], max[2]);
		fHistpt_assHadron_tr[i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{trigg. #gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
		fHistpt_assHadron_tr[i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp^{h}_{T} [counts/%0.1f GeV/c]",fHistpt_assHadron_tr[i]->GetBinWidth(1)));
		fOutput->Add(fHistpt_assHadron_tr[i]);
	}



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
Int_t AliAnalysisTaskGammaHadron::DoCellLoop(Float_t &sum)
{
	// Do cell loop.
	// Counts the number of fired cells
	AliVCaloCells *cells = InputEvent()->GetEMCALCells();

	if (!cells)
		return 0;

	const Int_t ncells = cells->GetNumberOfCells();
	Int_t nAccCells = 0;

	for (Int_t pos = 0; pos < ncells; pos++)
	{
		Float_t amp   = cells->GetAmplitude(pos);
		Int_t   absId = cells->GetCellNumber(pos);

		if (amp < fCellEnergyCut)
			continue;

		nAccCells++;
		sum += amp;
	}

	return nAccCells;
}

//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::DoClusterLoop(Float_t &sum, AliVCluster* &leading)
{
	/*
 	Trigger stuff stolen from AliAnaEMCALTriggerClusters
	Int_t  idTrig = GetReader()->GetTriggerClusterIndex();
	Bool_t exotic = GetReader()->IsExoticEvent();
	Bool_t bad    = GetReader()->IsBadCellTriggerEvent();

    TClonesArray * clusterList = 0;
 	clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetInputEvent() ->FindListObject(clusterListName));
    else if(GetReader()->GetOutputEvent())
    clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetOutputEvent()->FindListObject(clusterListName));

  AliVCluster  *  badClusTrig = 0;
  if(clusterList) badClusTrig = (AliVCluster*) clusterList->At(idTrig);
  else            badClusTrig = GetReader()->GetInputEvent()->GetCaloCluster(idTrig);


should I:
clusters->GetCluster(idTrig)

What is the difference between
>>AliClusterContainer
AliClusterContainer* clusters = GetClusterContainer(0);
cluster= clusters->GetAcceptCluster(idTrig);

and
>>clusterList
TString  clusterListName   = GetReader()->GetEMCALClusterListName();
TClonesArray * clusterList = dynamic_cast<TClonesArray*> (GetReader()->GetInputEvent() ->FindListObject(clusterListName));

cluster= (AliVCluster*) clusterList->At(idTrig);
cluster= GetReader()->GetInputEvent()->GetCaloCluster(idTrig);

*/

	// Do cluster loop.

	//get the trigger ID
	// this does not work unless I load the GammaHdronCorrelations: Int_t  idTrig = GetReader()->GetTriggerClusterIndex();
	AliClusterContainer* clusters = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;

	Int_t nAccClusters = 0;

	//AliVCaloCells *cells = InputEvent()->GetEMCALCells();

	sum = 0;
	leading = 0;
	Bool_t TRIGGER=0;
	//Cluster loop

	AliVCluster* cluster = 0;
	//need to load GammaHdronCorrelations: AliVCluster* Triggercluster = GetReader()->GetInputEvent()->GetCaloCluster(idTrig);
	clusters->ResetCurrentID();
	while ((cluster = clusters->GetNextAcceptCluster()))
	{
		TRIGGER=0;
		/*//need to load GammaHdronCorrelations:
		 * if(cluster==Triggercluster)
		{
			TRIGGER=1;
			cout<<"Found trigger particle at position: "<<idTrig<<endl;
		}
		else
		{
			cout<<"No trigger found"<<endl;
		}*/

		sum += cluster->E();

		if (!leading || leading->E() < cluster->E()) leading = cluster;
		TLorentzVector CaloClusterVec;
		cluster->GetMomentum(CaloClusterVec, fVertex);

		//all clusters
		fHistNoClus_pt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
		//all trigger clusters
		//need to load GammaHdronCorrelations: if(TRIGGER)fHistNoClus_pt_tr->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong

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

			//std::cout<<"CaloClusterVec.Phi(): "<<CaloClusterVec.Phi()<<", track->Phi(): "<<track->Phi()<<std::endl;
			deltaPhi=fabs(CaloClusterVec.Phi()-track->Phi());

			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and any hadron
				if(track->Pt()>0.5)           fHistNoClus_ptH->Fill(CaloClusterVec.Pt());    //the .pt only works for gammas (E=M) for other particle this is wrong
				//need to load GammaHdronCorrelations: if(track->Pt()>0.5 && TRIGGER)fHistNoClus_ptH_tr->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
				if(CaloClusterVec.Pt()>5 && track->Pt()>0.5)fHistNoClus_xEH->Fill(-1.0*cos(deltaPhi)*track->Pt()/CaloClusterVec.Pt());

				if(CaloClusterVec.Pt()>5)
				{
					//fill the different bins of the p_t of the associated hadron
					for(Int_t i=1;i<fHistNoClus_ptH->GetNbinsX()+1;i++)  // be careful hard coded from max[0] value -- need better variable // fHistNoClus_ptH->GetLast bin etc??
					{
						//look in each p_T bin of the gamma (fHistNoClus_ptH histogram) how the pt
						//distribution of the associated hadron looks like.
						Double_t BinWidth    = fHistNoClus_ptH->GetBinWidth(i);
						Double_t BinValStart = fHistNoClus_ptH->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;
						if(CaloClusterVec.Pt()>=BinValStart && CaloClusterVec.Pt()<BinValStart+BinWidth)
						{
							fHistpt_assHadron[i]->Fill(track->Pt());
							//if the gamma is the trigger then fill the histogram
							//need to load GammaHdronCorrelations: if(TRIGGER)fHistpt_assHadron_tr[i]->Fill(track->Pt());
						}
					}
				}
			}
		}
		//end of correlated hadron check!
		if(leadingTrack)
		{
			deltaPhi=fabs(CaloClusterVec.Phi()-leadingTrack->Phi());
			if(deltaPhi>(TMath::Pi()/2.0))
			{
				//cluster and leadig hadron
				if(leadingTrack->Pt()>0.5)fHistNoClus_ptLeadH->Fill(CaloClusterVec.Pt());    //  -- ELIANE the .pt only works for gammas (E=M) for other particl. this is wrong
				if(CaloClusterVec.Pt()>5 && leadingTrack->Pt()>0.5)fHistNoClus_xELeadH->Fill(-1.0*cos(deltaPhi)*leadingTrack->Pt()/CaloClusterVec.Pt());
			}
		}

		//fHistClusPhiEtaEnergy[fCentBin]->Fill(CaloClusterVec.Eta(), CaloClusterVec.Phi(), cluster->E());
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
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	Float_t clusSum = 0;
	Int_t nclusters = 0;
	if (fCaloClusters)
	{
		AliVCluster  *leadingClus = 0;
		nclusters = DoClusterLoop(clusSum, leadingClus);
		AliDebug(2,Form("%d clusters found in the event", nclusters));

		if (leadingClus)
		{
			TLorentzVector leadingClusVect;
			leadingClus->GetMomentum(leadingClusVect, fVertex);
		}
	}
	return kTRUE;
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
		}
		else
		{
			if (track->GetLabel() == 0)
			{
				zero++;
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
			}
			else
			{
				AliDebug(2,Form("%s: track type %d not recognized!", GetName(), type));
			}

			AliVTrack* vtrack = dynamic_cast<AliVTrack*>(track);
			if (!vtrack) continue;
		}
	}

	return nAccTracks;
}


