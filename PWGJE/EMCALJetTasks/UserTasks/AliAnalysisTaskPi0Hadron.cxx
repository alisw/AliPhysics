//
// Task to estimate the number of gamma-hadron
// statistic available in the Pb+Pb run.
//
// Author: E. Epple

#include <Riostream.h>
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

#include "AliAnalysisTaskPi0Hadron.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPi0Hadron)

//________________________________________________________________________
AliAnalysisTaskPi0Hadron::AliAnalysisTaskPi0Hadron() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskPi0Hadron", kTRUE),
fCellEnergyCut(0.05),
fMaxCellsInCluster(50),

fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentMethod2(""),
fCentMethod3(""),

fOutputList1(),
fOutputList2(),
fOutputList3(),

fHistClusPairInvarMass(0), //-<()>-//
fHistClusPairInvarMasspTSlice(0), 
fHistReadout(0), //-<()>-//
fHistClusPairInvarMasspT(0),
fHistClusPairInvarMassE(0),
fHistClusPairInvarMassPlay(0),

fHistNoClus_pt(0),
fHistNoClus_pt_tr(0),
fHistNoClus_ptH(0),
fHistNoClus_ptH_tr(0),
fHistpi0(0),
fHist_DetaDphi(),
fHistpt_assHadron(0),
fHistpt_assHadron_tr(0),
fHist_DP_gh(0)

{
	// Default constructor.
	fAODfilterBits[0] = 0;
	fAODfilterBits[1] = 0;

	fRtoD=180.0/TMath::Pi();

	SetMakeGeneralHistograms(kTRUE);  //What is this??
}
//________________________________________________________________________
AliAnalysisTaskPi0Hadron::~AliAnalysisTaskPi0Hadron()
{
	// Destructor

	//Copied from chris yaldo. Ask Salvatore about it!
	// Destructor. Clean-up the output list, but not the histograms that are put inside
	// (the list is owner and will clean-up these histograms). Protect in PROOF case.
	if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
	{
		delete fOutputList1;
	}
}
//________________________________________________________________________
void AliAnalysisTaskPi0Hadron::UserCreateOutputObjects()
{
	// Create histograms
	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

	//Define sublists/folders for a better organisation of the figures
	fOutputList1 = new TList();
	fOutputList1->SetOwner();
	fOutputList1->SetName("p_{T} distributions of the gamma");
	fOutputList2 = new TList();
	fOutputList2->SetOwner();
	fOutputList2->SetName("p_{T} distributions of the associated hadrons, for a given p_{T}^{gamma}");
	fOutputList3 = new TList();
	fOutputList3->SetOwner();
	fOutputList3->SetName("Delta phi^{g-h} for a given p_{T}^{gamma}");

	//common bins for the histograms
	Int_t nbins[5] = {0};
	Double_t min[5] = {0};
	Double_t max[5] = {0};

	//settings for p_t cluster distributon
	nbins[0] = 31;
	min[0] = 0;
	max[0] = 31;
	//settings for p_t hadron distribution
	nbins[1] = 60;  //do 1/2 GeV bins so that you can see the 0.5 cut to set as a minimum pT to combine hadron and gamma
	min[1] = 0;
	max[1] = 30;
	//settings for delta phi (g-h) distribution
	nbins[2] = 100;
	min[2] = -300;
	max[2] = 300;
	//settings for delta eta (g-h) distribution
	nbins[3] = 30;
	min[3] = -3;
	max[3] = 3;
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Create Histograms
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	fHistClusPairInvarMass = new TH1F("fHistClusPairInvarMass","Invariant Mass of Cluster Pairs;Invariant Mass of Pair;dN/dm",200, 0, 1); //-<()>-//
	fOutput->Add(fHistClusPairInvarMass); //-<()>-//
       	fHistClusPairInvarMasspTSlice = new TH2F("fHistClusPairInvarMasspTSlice","Invariant Mass of Cluster Pairs;Invariant Mass of Pair;dN/dm",200, 0, 1, 100, 0, 20); //-<()>-//
	fOutput->Add(fHistClusPairInvarMasspTSlice); //-<()>-//
	
	fHistClusPairInvarMasspT = new TH2F("fHistClusPairInvarMasspT","Invariant Mass of Cluster Pairs;Invariant Mass of Pair;p_{T} of Pi0",200, 0, 1, 100, 0, 20);
	fOutput->Add(fHistClusPairInvarMasspT); //-<()>-//
	fHistClusPairInvarMassE = new TH2F("fHistClusPairInvarMassE","Invariant Mass of Cluster Pairs;Invariant Mass of Pair;Energy of Pi0",200, 0, 1, 100, 0, 5);
	fOutput->Add(fHistClusPairInvarMassE); //-<()>-//
	fHistClusPairInvarMassPlay = new TH2F("fHistClusPairInvarMassPlay","Invariant Mass of Cluster Pairs;Invariant Mass of Pair; Angular distance",200, 0, 1, 100, 0, 6); 
	fOutput->Add(fHistClusPairInvarMassPlay); //-<()>-//
	fHistReadout = new TH1F("fHistReadout","Particle ID",1000, 0, .004); //-<()>-//
	fOutput->Add(fHistReadout); //-<()>-//
	
	//.................................
	// p_T^{Cluster} distribution under different conditions

	// all clusters
	fHistNoClus_pt = new TH1F(Form("fHistNoClus_pt_%0d",1),Form("fHistNoClus_pt_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_pt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_pt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClus_pt->GetBinWidth(0)));

	//trigger clusters as a function of p_T^{Cluster}
	fHistNoClus_pt_tr = new TH1F(Form("fHistNoTrClus_pt_tr_%0d",1),Form("fHistNoTrClus_pt_tr_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_pt_tr->GetXaxis()->SetTitle("p_{T}^{trig. #gamma}");
	fHistNoClus_pt_tr->GetYaxis()->SetTitle(Form("No. of trig. #gamma [counts/%0.1f GeV/c]",fHistNoClus_pt_tr->GetBinWidth(0)));

	//clusters p_T if there was a hadron present
	fHistNoClus_ptH = new TH1F(Form("fHistNoClus_ptH_%0d",1),Form("fHistNoClus_ptH_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_ptH->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptH->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH->GetBinWidth(0)));

	//trigger clusters p_T if there was a hadron present
	fHistNoClus_ptH_tr = new TH1F(Form("fHistNoClus_ptH_tr_%0d",1),Form("fHistNoClus_ptH_tr_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_ptH_tr->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_ptH_tr->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH_tr->GetBinWidth(0)));

	fOutputList1->Add(fHistNoClus_pt);
	fOutputList1->Add(fHistNoClus_pt_tr);
	fOutputList1->Add(fHistNoClus_ptH);
	fOutputList1->Add(fHistNoClus_ptH_tr);
	fOutput->Add(fOutputList1);

	//.................................
	// two dimensional delta eta delta phi distributions
	fHist_DetaDphi = new TH2F(Form("fHist_DetaDphi_%0d",1),Form("fHist_DetaDphi_%0d",1),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
	fHist_DetaDphi->GetXaxis()->SetTitle(Form("#Delta #Phi^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",0,100));
	fHist_DetaDphi->GetYaxis()->SetTitle(Form("#Delta #eta^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",0,100));
	fOutput->Add(fHist_DetaDphi);

	//test!!

	fHistpi0 = new TH1F(Form("fHistpi0_%0d",1),Form("fHistpi0_%0d",1), 500, 0, 0.5);
	fHistpi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistpi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistpi0);



	//.................................
	//Initiate histograms for given p_t gamma bins!
	//Later you can summ the histograms according to the expected statistic!
	Int_t NoOfDPhistos =  nbins[0];

	fHistpt_assHadron    = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHistpt_assHadron_tr = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the trigger-gamma
	fHist_DP_gh          = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma

	//
	for(Int_t i=0; i<NoOfDPhistos+1; i++)
	{
		//check whether the max is the
		//same as the histogram size
		Double_t BinWidth,BinValStart;

		if(i==0)
		{
			//these are histograms over the full p_T gamma range (no binnings)
			BinWidth    = fHistNoClus_ptH->GetBinWidth(1);
			BinValStart = fHistNoClus_ptH->GetXaxis()->GetBinCenter(1)-BinWidth/2.0;
			BinWidth    =100; //to get the full length of the histogram
		}
		else
		{
			//these are histograms for a certain p_T of the gamma
			BinWidth    = fHistNoClus_ptH->GetBinWidth(i);
			BinValStart = fHistNoClus_ptH->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;
			//cout<<"BinStart: "<<fHistNoClus_ptH->GetBinCenter(i)<<", width "<<BinWidth<<endl;
		}

		//.................................
		//Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
		fHistpt_assHadron[i] = new TH1F(Form("fHistAssHadron_pt_%0d",i),Form("fHistAssHadron_pt_%0d",i), nbins[1], min[1], max[1]);
		fHistpt_assHadron[i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
		fHistpt_assHadron[i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistpt_assHadron[i]->GetBinWidth(1)));

		fHistpt_assHadron_tr[i] = new TH1F(Form("fHistAssHadron_pt_tr_%0d",i),Form("fHistAssHadron_pt_tr_%0d",i), nbins[1], min[1], max[1]);
		fHistpt_assHadron_tr[i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{trigg. #gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
		fHistpt_assHadron_tr[i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp^{h}_{T} [counts/%0.1f GeV/c]",fHistpt_assHadron_tr[i]->GetBinWidth(1)));

		if(i!=0)fOutputList2->Add(fHistpt_assHadron[i]);
		if(i!=0)fOutputList2->Add(fHistpt_assHadron_tr[i]);

		//.................................
		//Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
		fHist_DP_gh[i] = new TH1F(Form("fHist_DP_gh_%0d",i),Form("fHist_DP_gh_%0d",i), nbins[2], min[2], max[2]);
		fHist_DP_gh[i]->GetXaxis()->SetTitle(Form("#Delta #Phi^{#gamma-h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
		fHist_DP_gh[i]->GetYaxis()->SetTitle(Form("dN^{#gamma-h}/#Delta #Phi^{#gamma-h} [counts/%0.1f^{#circ}]",fHist_DP_gh[i]->GetBinWidth(1)));

		if(i!=0)fOutputList3->Add(fHist_DP_gh[i]);
	}
	fOutput->Add(fOutputList2);
	fOutput->Add(fOutputList3);
	fOutput->Add(fHistpt_assHadron[0]);
	fOutput->Add(fHistpt_assHadron_tr[0]);
	fOutput->Add(fHist_DP_gh[0]);

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPi0Hadron::ExecOnce()
{
	AliAnalysisTaskEmcalJet::ExecOnce();
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Hadron::RetrieveEventObjects()
{
	// Retrieve event objects.
	if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
	{
		return kFALSE;
	}
	if (!fCentMethod2.IsNull() || !fCentMethod3.IsNull())
	{
		if (fBeamType == kAA || fBeamType == kpA )
		{
			AliCentrality *aliCent = InputEvent()->GetCentrality();
			if (aliCent)
			{
				if (!fCentMethod2.IsNull())
				{
					//	fCent2 = aliCent->GetCentralityPercentile(fCentMethod2);
				}
				if (!fCentMethod3.IsNull())
				{
					//	fCent3 = aliCent->GetCentralityPercentile(fCentMethod3);
				}
			}
		}
	}
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Hadron::FillHistograms()
{
	Float_t clusSum = 0;
	Int_t nclusters = 0;
	if (fCaloClusters)
	{
		AliVCluster  *leadingClus = 0;
		nclusters = DoClusterLoop(clusSum, leadingClus);
		AliDebug(2,Form("%d clusters found in the event", nclusters));
	}
	return kTRUE;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPi0Hadron::DoCellLoop(Float_t &sum)
{
	// Do cell loop.
	// Counts the number of fired cells
	AliVCaloCells *cells = InputEvent()->GetEMCALCells();

	if(!cells) return 0;

	const Int_t ncells = cells->GetNumberOfCells();
	Int_t nAccCells = 0;

	for(Int_t pos = 0; pos < ncells; pos++)
	{
		Float_t amp   = cells->GetAmplitude(pos);
		Int_t   absId = cells->GetCellNumber(pos);

		if(amp < fCellEnergyCut) continue;

		nAccCells++;
		sum += amp;
	}

	return nAccCells;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPi0Hadron::DoClusterLoop(Float_t &sum, AliVCluster* &leading)
{
	//................................
	// Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	AliParticleContainer* tracks   = GetParticleContainer(0);
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Int_t nAccClusters = 0;
	sum = 0;
	leading = 0;

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* leadingTrack = 0;
	AliVParticle* track = 0;

	Int_t Loop_counter=0;

	TLorentzVector SumCaloClusterVec; //-<()>-//	
	AliVCluster* othercluster = 0;//-<()>-//
	clusters->ResetCurrentID();
	Int_t naccept = clusters->GetNAcceptedClusters(); //<>//
	Int_t nclust = clusters->GetNClusters(); //<>//
	Int_t loop = 0;
	clusters->ResetCurrentID(); //<>//

	for( Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent-1; NoCluster1++ )
	{//first cluster

		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1);
		//cluster=(AliVCluster*) clusters->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster

		sum += cluster->E();
		if (!leading || leading->E() < cluster->E()) leading = cluster;

		TLorentzVector CaloClusterVec;
		cluster->GetMomentum(CaloClusterVec, fVertex);

		fHistNoClus_pt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong

		//................................
		// Do track loop.
		if (!tracks)  return 0;
		tracks->ResetCurrentID();

		while((track = tracks->GetNextAcceptParticle()))
		{
			if (!leadingTrack || leadingTrack->Pt() < leadingTrack->Pt()) leadingTrack = track;
			Fill_GH_Hisograms(1,CaloClusterVec,track,2,0,-360);

			if(leading)
			{
				TLorentzVector leadingClusVect;
				leading->GetMomentum(leadingClusVect, fVertex);
				/// use different identifier! Fill_GH_Hisograms(2,leadingClusVect,track,0,0,0);
			}
		}
		//end of correlated hadron check!
		nAccClusters++;

		//double cluster loop for testing
		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);

				if(!cluster2 || !AccClusterForAna(cluster2))continue; //check if the cluster is a good cluster
				TLorentzVector CaloClusterVec2;
				TLorentzVector CaloClusterVecpi0;
				cluster2->GetMomentum(CaloClusterVec2, fVertex);
				if(cluster2->E()>2 && cluster->E()>2)
				{
					CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
					fHistpi0->Fill(CaloClusterVecpi0.M());
				}
			}
		}
		Loop_counter++;


		for(Int_t i=loop;i<=nclust-1;i++)    //-<()>--
		  {
		    othercluster=clusters->GetCluster(i);
		    if((othercluster)&&(othercluster != cluster)&&(AcceptCluster(othercluster)))
		      {
			TLorentzVector OtherCaloClusterVec;
			othercluster->GetMomentum(OtherCaloClusterVec, fVertex);
			//clusters->GetMomentum(OtherCaloClusterVec, i);
			SumCaloClusterVec = OtherCaloClusterVec + CaloClusterVec;
			fHistClusPairInvarMass->Fill(SumCaloClusterVec.M());
			fHistClusPairInvarMasspTSlice->Fill(SumCaloClusterVec.M(),SumCaloClusterVec.Pt());
			fHistClusPairInvarMasspT->Fill(SumCaloClusterVec.M(),SumCaloClusterVec.Pt());
			fHistClusPairInvarMassE->Fill(SumCaloClusterVec.M(),SumCaloClusterVec.E());
			Double_t twoangsize= sqrt(pow(OtherCaloClusterVec.Theta()-CaloClusterVec.Theta(),2)+pow(OtherCaloClusterVec.Phi()-CaloClusterVec.Phi(),2));
			fHistClusPairInvarMassPlay->Fill(SumCaloClusterVec.M(),OtherCaloClusterVec.Phi()-CaloClusterVec.Phi()+3.141598265259);
		      }
		  }
		fHistReadout->Fill(CaloClusterVec.M()); //<()>//


		
	}
	//if(clusters->GetLeadingCluster("e")!=leading)cout<<"why did that happen?"<<endl;

	return nAccClusters;
}
//________________________________________________________________________
Int_t AliAnalysisTaskPi0Hadron::DoTrackLoop(Float_t &sum, AliVParticle* &leading)
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

	}

	return nAccTracks;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Hadron::AccClusterForAna(AliVCluster* cluster)
{
	//Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=0; //By default rejceted

	//double check these cuts carefully with the experts
	if(cluster->E()>0.3 && cluster->GetNCells()>1)
	{
		//Now accept the cluster as a good candidate for your analysis
		Accepted=1;
	}

	return Accepted;
}
//________________________________________________________________________
Double_t AliAnalysisTaskPi0Hadron::DeltaPhi(TLorentzVector ClusterVec,AliVParticle* TrackVec)
{
	Double_t Phi_g = ClusterVec.Phi();
	Double_t Phi_h = TrackVec->Phi();

	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();

	if(Phi_g<0)   cout<<"problem1?? phig "<<Phi_g<<endl;
	if(Phi_g>2*pi)cout<<"problem2?? phig "<<Phi_g<<endl;



	if (Phi_g < 0)         Phi_g += 2*pi;
	else if (Phi_g > 2*pi) Phi_g -= 2*pi;
	if (Phi_h < 0)         Phi_h += 2*pi;
	else if (Phi_h)        Phi_h -= 2*pi;


	dPhi = Phi_g-Phi_h;


	//if (dPhi < rangeMin)      dPhi += 2*pi;
	//else if (dPhi > rangeMax) dPhi -= 2*pi;

	return dPhi;
}
//________________________________________________________________________
void AliAnalysisTaskPi0Hadron::Fill_GH_Hisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut)
{
	//this function fill several histograms under different cut conditions.
	//it is run within a cluster{ track{}} loop to get all combinations.

	//important!! the identifier is right now useless because there is only on set of
	//histograms. Later this set can be cloned so that we can check several types of different cuts in paralell

	Double_t deltaPhi,deltaEta;
	Double_t DP_alternative = DeltaPhi(ClusterVec,TrackVec);
	deltaPhi=fRtoD*(ClusterVec.Phi()-TrackVec->Phi());
	deltaEta=ClusterVec.Eta()-TrackVec->Eta();

	if(ClusterVec.Pt()>=ClusterEcut && TrackVec->Pt()>=TrackPcut && deltaPhi>=Anglecut)
	{
		fHistNoClus_ptH->Fill(ClusterVec.Pt());    //the .pt only works for gammas (E=M) for other particle this is wrong
		fHistpt_assHadron[0]->Fill(TrackVec->Pt());
		fHist_DP_gh[0]      ->Fill(deltaPhi);
		fHist_DetaDphi      ->Fill(deltaPhi,deltaEta);

		//fill histograms for different ranges of p_t^{g}
		for(Int_t i=1;i<fHistNoClus_ptH->GetNbinsX()+1;i++)  // be careful hard coded from max[0] value -- need better variable // fHistNoClus_ptH->GetLast bin etc??
		{
			//look in each p_T bin of the gamma (fHistNoClus_ptH histogram) how the pt
			//distribution of the associated hadron looks like.
			Double_t BinWidth    = fHistNoClus_ptH->GetBinWidth(i);
			Double_t BinValStart = fHistNoClus_ptH->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;

			//right now the ranges are defined via different histogram bins of the gmma p_t histogram
			// that might be changed later.
			if(ClusterVec.Pt()>=BinValStart && ClusterVec.Pt()<BinValStart+BinWidth)
			{
				fHistpt_assHadron[i]->Fill(TrackVec->Pt());
				fHist_DP_gh[i]      ->Fill(deltaPhi);
			}
		}

	}
}
