/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

// Comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
// Simple class for the jt anlyais by Beomkyu Kim and Dongjo Kim
//===========================================================

#include <TRandom.h>
#include <TMath.h>
#include <TRegexp.h>
#include <TVector.h>
#include "AliJJet.h"
#include "AliJEfficiency.h"
#include "AliJJetJtAnalysis.h"
#include "AliJHistManager.h"
#include "TClonesArray.h"

AliJJetJtAnalysis::AliJJetJtAnalysis():
	fInputList(NULL)
	, fJetList(NULL)
	, fJetListOfList()
	, fJetBgList(NULL)
	, fJetBgListOfList()
    , fJetTriggPtBorders(NULL)
    , fJetConstPtLowLimits(NULL)
    , fJetAssocPtBorders(NULL)
    , nJetContainer(0)
	, fCard(NULL)
    , fJetFinderName(0)
    , fConeSizes(0)
	, fEfficiency(0x0)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(NULL)
    , fHMG(NULL)
    , fJetFinderBin()
    , fJetTriggerBin()
    , fTrkPtBin()
    , fTrkLimPtBin()
    , fhJetPt()
    , fhJetPtBin()
    , fhZ()
    , fhZBin()
    , fhJt()
    , fhJtBin()
    , fhJtWeightBin()
    , fhLogJtWeightBin()
    , fhJtWithPtCutWeightBinBin()
    , fhLogJtWithPtCutWeightBinBin()
    , fhJtBinLimBin()
    , fhJtWeightBinLimBin()
    , fhLogJtWeightBinLimBin()
    , fhJetBgPt()
    , fhJetBgPtBin()
    , fhBgZ()
    , fhBgZBin()
    , fhBgJt()
    , fhBgJtBin()
    , fhBgJtWeightBin()
    , fhBgLogJtWeightBin()
    , fhBgJtWithPtCutWeightBinBin()
    , fhBgLogJtWithPtCutWeightBinBin()
{

}

AliJJetJtAnalysis::AliJJetJtAnalysis( AliJCard * card ):
	fInputList(NULL)
	, fJetList(NULL)
	, fJetListOfList()
	, fJetBgList(NULL)
	, fJetBgListOfList()
    , fJetTriggPtBorders(NULL)
    , fJetConstPtLowLimits(NULL)
    , fJetAssocPtBorders(NULL)
    , nJetContainer(0)
	, fCard(card)
    , fJetFinderName(0)
    , fConeSizes(0)
	, fEfficiency(0x0)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(NULL)
    , fHMG(NULL)
    , fJetFinderBin()
    , fJetTriggerBin()
    , fTrkPtBin()
    , fTrkLimPtBin()
    , fhJetPt()
    , fhJetPtBin()
    , fhZ()
    , fhZBin()
    , fhJt()
    , fhJtBin()
    , fhJtWeightBin()
    , fhLogJtWeightBin()
    , fhJtWithPtCutWeightBinBin()
    , fhLogJtWithPtCutWeightBinBin()
    , fhJtBinLimBin()
    , fhJtWeightBinLimBin()
    , fhLogJtWeightBinLimBin()
    , fhJetBgPt()
    , fhJetBgPtBin()
    , fhBgZ()
    , fhBgZBin()
    , fhBgJt()
    , fhBgJtBin()
    , fhBgJtWeightBin()
    , fhBgLogJtWeightBin()
    , fhBgJtWithPtCutWeightBinBin()
    , fhBgLogJtWithPtCutWeightBinBin()
{

}

AliJJetJtAnalysis::AliJJetJtAnalysis(const AliJJetJtAnalysis& ap) :
	fInputList(ap.fInputList)
	, fJetList(ap.fJetList)
	, fJetListOfList(ap.fJetListOfList)
	, fJetBgList(ap.fJetBgList)
	, fJetBgListOfList(ap.fJetBgListOfList)
    , fJetTriggPtBorders(ap.fJetTriggPtBorders)
    , fJetConstPtLowLimits(ap.fJetConstPtLowLimits)
    , fJetAssocPtBorders(ap.fJetAssocPtBorders)
    , nJetContainer(ap.nJetContainer)
	, fCard(ap.fCard)
    , fJetFinderName(ap.fJetFinderName)
    , fConeSizes(ap.fConeSizes)
	, fEfficiency(ap.fEfficiency)
	, cBin(-1)
	, fcent(-999)
	, zBin(-1)
	, zVert(-999)
    , fTracks(ap.fTracks)
    , fHMG(ap.fHMG)
    , fJetFinderBin(ap.fJetFinderBin)
    , fJetTriggerBin(ap.fJetTriggerBin)
    , fTrkPtBin(ap.fTrkPtBin)
    , fTrkLimPtBin(ap.fTrkLimPtBin)
    , fhJetPt(ap.fhJetPt)
    , fhJetPtBin(ap.fhJetPtBin)
    , fhZ(ap.fhZ)
    , fhZBin(ap.fhZBin)
    , fhJt(ap.fhJt)
    , fhJtBin(ap.fhJtBin)
    , fhJtWeightBin(ap.fhJtWeightBin)
    , fhLogJtWeightBin(ap.fhLogJtWeightBin)
    , fhJtWithPtCutWeightBinBin(ap.fhJtWithPtCutWeightBinBin)
    , fhLogJtWithPtCutWeightBinBin(ap.fhLogJtWithPtCutWeightBinBin)
    , fhJtBinLimBin(ap.fhJtBinLimBin)
    , fhJtWeightBinLimBin(ap.fhJtWeightBinLimBin)
    , fhLogJtWeightBinLimBin(ap.fhLogJtWeightBinLimBin)
    , fhJetBgPt(ap.fhJetBgPt)
    , fhJetBgPtBin(ap.fhJetBgPtBin)
    , fhBgZ(ap.fhBgZ)
    , fhBgZBin(ap.fhBgZBin)
    , fhBgJt(ap.fhBgJt)
    , fhBgJtBin(ap.fhBgJtBin)
    , fhBgJtWeightBin(ap.fhBgJtWeightBin)
    , fhBgLogJtWeightBin(ap.fhBgLogJtWeightBin)
    , fhBgJtWithPtCutWeightBinBin(ap.fhBgJtWithPtCutWeightBinBin)
    , fhBgLogJtWithPtCutWeightBinBin(ap.fhBgLogJtWithPtCutWeightBinBin)
{

}

AliJJetJtAnalysis& AliJJetJtAnalysis::operator = (const AliJJetJtAnalysis& ap)
{
	// assignment operator

	this->~AliJJetJtAnalysis();
	new(this) AliJJetJtAnalysis(ap);
	return *this;
}


AliJJetJtAnalysis::~AliJJetJtAnalysis(){


    fJetFinderName.clear();
    fConeSizes.clear();
    delete fEfficiency;
    delete fHMG;
   


}



void AliJJetJtAnalysis::UserCreateOutputObjects(){
	//fJetListOfList always point one address in the whole time of this analysis.
    //Thus mustn't be cleared in it's life.     
    //fJetListOfList.Clear();



    fJetTriggPtBorders = fCard->GetVector("JetTriggPtBorders");
    fJetConstPtLowLimits = fCard->GetVector("JetConstPtLowLimits");
    fJetAssocPtBorders = fCard->GetVector("JetAssocPtBorders");

	fEfficiency = new AliJEfficiency();
	fEfficiency->SetMode( fCard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien

    TRegexp reg("R[0-9][0-9][0-9]");
    TRegexp reg2("[0-9][0-9][0-9]");

    //container name has information of cone size like **R040**
    //this cone size information will be pulled to a numerical variable
    nJetContainer = fJetFinderName.size();
	for (int i=0; i<nJetContainer; i++){
		//AliJJetJtHistos *histo = new AliJJetJtHistos(fCard);
		//histo->CreateJetJtHistos();
	//	fHistos.push_back( histo  );
        TString fullNameOfiJetContainer(fJetFinderName[i]);
        TString coneSizeName (fullNameOfiJetContainer(reg));
        TString coneSizeValue (coneSizeName(reg2));
        fConeSizes.push_back( (double) coneSizeValue.Atoi()/100.);
	}


    int NBINS=150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=150;
    double logBW = (log(LimH)-log(LimL))/NBINS;
    for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);



    fHMG = new AliJHistManager( "AliJJetJtHistManager");
    fJetFinderBin .Set("JetFinderOrder","NFin","NFin:%d", AliJBin::kSingle).SetBin(nJetContainer);
    fJetTriggerBin .Set("JetTriggerBin","JetPt","p_{T,jet}:%2.0f-%2.0f").SetBin(fCard->GetVector("JetTriggPtBorders"));
    fTrkPtBin .Set("TrkPtBin","TrkPt","p_{T,constituent}:%2.0f-%2.0f").SetBin(fCard->GetVector("JetAssocPtBorders"));
    fTrkLimPtBin .Set("TrkLimitPtBin","TrkLimitPt","p_{T,Limit}:%2.0f-%2.0f").SetBin(fCard->GetVector("JetConstPtLowLimits"));
    fhJetPt 
        << TH1D("JetPt","",NBINS, LogBinsX ) << fJetFinderBin
        <<"END";
    fhJetPtBin 
        << TH1D("JetPtBin","",NBINS, LogBinsX ) << fJetFinderBin << fJetTriggerBin
        <<"END";

    int NBINSZ=150;
    double LogBinsZ[NBINSZ+1], LimLZ=0.001, LimHZ=1.1;
    double logBWZ = (TMath::Log(LimHZ)-TMath::Log(LimLZ))/NBINSZ;
    for(int ij=0;ij<=NBINSZ;ij++) LogBinsZ[ij]=LimLZ*exp(ij*logBWZ);//

    fhZ 
        << TH1D("Z","",NBINSZ, LogBinsZ ) << fJetFinderBin
        <<"END";
    fhZBin 
        << TH1D("ZBin","",NBINSZ, LogBinsZ ) << fJetFinderBin << fJetTriggerBin
        <<"END";

    int NBINSJt=150;
    double LogBinsJt[NBINSJt+1], LimLJt=0.01, LimHJt=10;
    double logBWJt = (TMath::Log(LimHJt)-TMath::Log(LimLJt))/NBINSJt;
    for(int ij=0;ij<=NBINSJt;ij++) LogBinsJt[ij]=LimLJt*exp(ij*logBWJt);
    int NBINSJtW=150;
    double LimLJtW=TMath::Log(0.01), LimHJtW=TMath::Log(10);

    fhJt 
        << TH1D("Jt","",NBINSJt, LogBinsJt ) << fJetFinderBin
        <<"END";
    fhJtBin 
        << TH1D("JtBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhJtWeightBin 
        << TH1D("JtWeightBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhLogJtWeightBin 
        << TH1D("LogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) << fJetFinderBin << fJetTriggerBin
        <<"END";

    fhJtWithPtCutWeightBinBin 
        << TH1D("JtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";
    fhLogJtWithPtCutWeightBinBin 
        << TH1D("LogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";

    fhJtBinLimBin 
        << TH1D("JtBinLimBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhJtWeightBinLimBin 
        << TH1D("JtWeightBinLimBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";
    fhLogJtWeightBinLimBin 
        << TH1D("LogJtWeightBinLimBin","",NBINSJtW, LimLJtW, LimHJtW ) << fJetFinderBin << fJetTriggerBin << fTrkLimPtBin
        <<"END";

    fhJetBgPt 
        << TH1D("JetBgPt","",NBINS, LogBinsX ) << fJetFinderBin
        <<"END";
    fhJetBgPtBin 
        << TH1D("JetBgPtBin","",NBINS, LogBinsX ) << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhZ 
        << TH1D("Z","",NBINSZ, LogBinsZ ) << fJetFinderBin
        <<"END";
    fhZBin 
        << TH1D("ZBin","",NBINSZ, LogBinsZ ) << fJetFinderBin << fJetTriggerBin
        <<"END";


    fhBgJt 
        << TH1D("BgJt","",NBINSJt, LogBinsJt ) << fJetFinderBin
        <<"END";
    fhBgJtBin 
        << TH1D("BgJtBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhBgJtWeightBin 
        << TH1D("BgJtWeightBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin
        <<"END";
    fhBgLogJtWeightBin 
        << TH1D("BgLogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) << fJetFinderBin << fJetTriggerBin
        <<"END";

    fhBgJtWithPtCutWeightBinBin 
        << TH1D("BgJtWithPtCutWeightBinBin","",NBINSJt, LogBinsJt ) << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";
    fhBgLogJtWithPtCutWeightBinBin 
        << TH1D("BgLogJtWeightBin","",NBINSJtW, LimLJtW, LimHJtW ) << fJetFinderBin << fJetTriggerBin << fTrkPtBin
        <<"END";

    fHMG->Print();
    fHMG->WriteConfig();




}

void AliJJetJtAnalysis::ClearBeforeEvent(){
	//fJetListOfList.Clear();


}

void AliJJetJtAnalysis::UserExec(){
	for( int i=0;i<fJetListOfList.GetEntries();i++ ){
		TObjArray * Jets = (TObjArray*) fJetListOfList[i];
		if(!Jets) {
			continue;
		}
		this->FillJtHistogram(Jets,i);
	}
}

void AliJJetJtAnalysis::WriteHistograms(){


	TDirectory * cwd = gDirectory;
	//const int nJetContainer = fJetListOfList.GetEntries();


	for (int i=0; i<nJetContainer; i++){
		TDirectory *nwd = gDirectory->mkdir(fJetFinderName[i]);
		nwd->cd();
		//fHistos[i]->WriteHistograms();
		cwd->cd();
	}


}



void AliJJetJtAnalysis::FillJtHistogram( TObjArray *Jets , int iContainer)
{	





	int iBin, iptaBin=0;
    int jBin=0;
	double pT = 0;
	double conPtMax =0;

	double z; double jt;
	double pta;
	//double Y , deltaY = 0;
	//double Phi, deltaPhi;
	//double deltaR= 0;
	//cout<<"histogram filling number of jets : "<<Jets->GetEntriesFast()<<endl;

    TLorentzVector  vOrtho;
    


    int k = 0;
    double deltaR = -1;
    double deltaEta = -999;
    double deltaPhi = -999;
    double effCorrection = -1;
    double thisConeSize = fConeSizes[iContainer] ;

    // iJet loop for an event
	for (int i = 0; i<Jets->GetEntries(); i++){
        AliJJet *jet = dynamic_cast<AliJJet*>( Jets->At(i) );
		pT = jet->Pt();
		if (pT<(*fJetTriggPtBorders)[1]) continue;
		iBin = GetBin(fJetTriggPtBorders,pT); // fill jetPt histos
        if( iBin < 0 ) continue;
        cout<<"iContainer "<<iContainer<<endl;
		fhJetPt[iContainer]->Fill( pT );
		fhJetPtBin[iContainer][iBin]->Fill( pT );

		for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
			AliJBaseTrack *con = jet->GetConstituent(icon);
            if (con->Pt()>conPtMax) conPtMax = con->Pt();
        }

        for (int ii = fJetConstPtLowLimits->GetNoElements(); ii >= 1 ; ii--)   {       // could also be done using GetBin( ... )
            if (conPtMax > (*fJetConstPtLowLimits)[ii]) {               // if JetConstPtLowLimits={a,...,b} -> ConPtBinBorders={a,...,b,c}
                jBin = ii-1;                                                              // where c(>>b) is ''sufficiently'' high
                break;
            }
        }


		//iConstituent loop for the iJet
        //jt, z are calcualted and filled  
		for (int icon = 0; icon<jet->GetConstituents()->GetEntries(); icon++){
			AliJBaseTrack *constituent = jet->GetConstituent(icon);
			z = (constituent->Vect()*jet->Vect().Unit())/jet->P();
			pta = constituent->Pt();
            constituent->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
            effCorrection = 1.0/constituent->GetTrackEff();
			iptaBin = GetBin(fJetAssocPtBorders, pta);
            if( iptaBin < 0 ) continue;


			fhZ[iContainer]->Fill( z , effCorrection);
			fhZBin[iContainer][iBin]->Fill( z , effCorrection);
			jt = (constituent->Vect()-z*jet->Vect()).Mag();
			fhJt[iContainer]->Fill( jt , effCorrection);
			fhJtBin[iContainer][iBin]->Fill( jt , effCorrection);
			fhJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
			fhLogJtWeightBin[iContainer][iBin]->Fill( TMath::Log(jt), 1.0/jt * effCorrection );



			if (iptaBin < 0) continue;
			fhJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill( jt, 1.0/jt * effCorrection );
			fhLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill( TMath::Log(jt), 1.0/jt * effCorrection);

            for (int jj = 0; jj <= jBin ; jj++) {
                fhJtBinLimBin[iContainer][iBin][jj]->Fill( jt, effCorrection );
                fhJtWeightBinLimBin[iContainer][iBin][jj]->Fill( jt, 1.0/jt * effCorrection );
                fhLogJtWeightBinLimBin[iContainer][iBin][jj]->Fill( TMath::Log(jt), 1.0/jt * effCorrection );
                //histos->fHistosJT[0][0][iBin][jj]->Fill( TMath::Log(jt), 1.0/jt );
            }

		}



        vOrtho.SetVect(jet->Vect().Orthogonal());
        vOrtho.SetE(jet->E());
        //if (Log) cout<<"Before R caluation, R = "<<TMath::Sqrt(it->area()/TMath::Pi())<<endl;
        //R_area = TMath::Sqrt(it->area()/TMath::Pi())*Rs[order]/R;
        //if (Log )cout<<"Rs[order] = "<<Rs[order]<<" R = "<<R<<" Bg R area : "<<R_area<<endl;

        //Background jet (iBgJet) will be produced. This background jet is orthogonal to the iJet.  
        //If there is another jJet, then iBgJet will be consecutevely moved not to have jJet in the cone size. 
        if (Jets->GetEntries()>1){
            for (int j = 0; j<Jets->GetEntries(); j++){
                if (i == j) continue;
                AliJJet *jet2 = dynamic_cast<AliJJet*>( Jets->At(j) );
                
                if (k>15) {
                    break;
                }

                    
                deltaEta = vOrtho.Eta() - jet2->Eta();
                deltaPhi = vOrtho.Phi() - jet2->Phi();
                deltaR   = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
                if ( deltaR < thisConeSize) {
        
                    vOrtho.Rotate(TMath::Pi()/8, jet->Vect());
                    j=0;
                    k++;
                }

            }
       }

       // Filling iBgJet,  Bgjt and Bgz
       // "k<16" means that we will select a iBgJet which hasn't moved 
       // more than 16 times by the process above
       if ( k<16 ){
           pT = vOrtho.Pt(); 
           if (pT<(*fJetTriggPtBorders)[1]) continue;
           
           fhJetBgPt[iContainer]->Fill( pT );
           //bbfHistos[iContainer]->fhJetBgPtWeight->Fill( pT, 1./pT);
           iBin = GetBin(fJetTriggPtBorders, pT);
		   fhJetBgPtBin[iContainer][iBin]->Fill( pT );
            if( iBin < 0 ) continue;


		   for (int icon = 0; icon<fTracks->GetEntries(); icon++){
                AliJBaseTrack *track = dynamic_cast<AliJBaseTrack*>(fTracks->At(icon));
                if (!track) continue;
                deltaEta = vOrtho.Eta() - track->Eta();
                deltaPhi = vOrtho.Phi() - track->Phi();
                deltaR   = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
                if ( deltaR > thisConeSize) continue;

                pta = track->Pt();
                track->SetTrackEff( fEfficiency->GetCorrection( pta, 5, fcent) );
                effCorrection = 1.0/track->GetTrackEff();
			    iptaBin = GetBin(fJetAssocPtBorders, pta);
                if( iptaBin < 0 ) continue;
                 
                
                z = (track->Vect()*vOrtho.Vect().Unit())/vOrtho.P();
                fhBgZ[iContainer]->Fill( z , effCorrection);
                fhBgZBin[iContainer]->Fill( z , effCorrection);

                jt = (track->Vect()-z*vOrtho.Vect()).Mag();
                fhBgJt[iContainer]->Fill( jt , effCorrection);
                fhBgJtBin[iContainer][iBin]->Fill( jt , effCorrection);
                fhBgJtWeightBin[iContainer][iBin]->Fill( jt, 1.0/jt * effCorrection );
                fhBgLogJtWeightBin[iContainer][iBin]->Fill( TMath::Log(jt), 1.0/jt * effCorrection );

			    if (iptaBin < 0) continue;
			    fhBgJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill( jt, 1.0/jt * effCorrection );
			    fhBgLogJtWithPtCutWeightBinBin[iContainer][iBin][iptaBin]->Fill( TMath::Log(jt), 1.0/jt * effCorrection );

           } 
       }


	}
}



