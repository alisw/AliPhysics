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

/**************************************************************************
 * AliJXtHistos.h
 * This class encapsulated all histograms that the analysis provides
 *
 *
 * contact: Sami Räsänen
 *          University of Jyväskylä, Finland
 *          sami.s.rasanen@jyu.fi
 **************************************************************************/

#include  "AliJXtHistos.h"
#include  "AliJCard.h"

//______________________________________________________________________________
AliJXtHistos::AliJXtHistos(AliJCard* cardP)
      :	fUseDirectory(true),
	fTopDirectory(0x0),
	fCard(0x0),
        fmaxEtaRange(0.8),
        fhistoList(0x0),
        fhConeActivity(0x0),
        fhConeActivityIsolated(0x0),
        fhZVertRaw(0x0),
        fhCentr(0x0),
        fhiCentr(0x0),
        fhEventPerRun(0x0),
        fhVertexZTriggVtx(0x0)
{
	// constructor
	fCard=cardP;
	fmaxEtaRange = fCard->Get("EtaRange");

	fhistoList = new TList();
	fTopDirectory = gDirectory;
}

//______________________________________________________________________________
AliJXtHistos::AliJXtHistos(const AliJXtHistos& obj)
      :	fUseDirectory(obj.fUseDirectory),
	fTopDirectory(obj.fTopDirectory),
	fCard(obj.fCard),
	fmaxEtaRange(obj.fmaxEtaRange),
	fhistoList(obj.fhistoList),
	fhConeActivity(obj.fhConeActivity),
	fhConeActivityIsolated(obj.fhConeActivityIsolated),
	fhZVertRaw(obj.fhZVertRaw),
	fhCentr(obj.fhCentr),
	fhiCentr(obj.fhiCentr),
	fhEventPerRun(obj.fhEventPerRun),
	fhVertexZTriggVtx(obj.fhVertexZTriggVtx)
{
	// copy constructor
}

//______________________________________________________________________________
AliJXtHistos& AliJXtHistos::operator=(const AliJXtHistos& obj){
	// copy constructor
        this->~AliJXtHistos();
        new(this) AliJXtHistos(obj);
	return *this;
}
//______________________________________________________________________________
TDirectory * AliJXtHistos::MakeDirectory(TString name){
    // Make directory with given name and move into it
    JumpToDefaultDirectory();
    TDirectory * dir = gDirectory->GetDirectory(name);
    if( !dir ) dir = gDirectory->mkdir(name);
    dir->cd();
    return dir;
}
//______________________________________________________________________________
TDirectory * AliJXtHistos::JumpToDefaultDirectory(){
    // Move to default directory
    fTopDirectory->cd();
    return gDirectory;
}
//______________________________________________________________________________
void AliJXtHistos::FillInclusiveHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin){
    // Fill all inclusive histograms
    fhChargedEta[centBin]->Fill(eta, effCorr);
    fhChargedPhi[centBin]->Fill(phi, effCorr);
    fhChargedPt[centBin]->Fill(pT, effCorr);
    fhChargedXt[centBin]->Fill(xT, effCorr);
    fhInvariantChargedPt[centBin]->Fill(pT, pT > 0. ? effCorr/pT : 0.); // Can you do 1/pT here???
    fhInvariantChargedXt[centBin]->Fill(xT, xT > 0. ? effCorr/xT : 0.); // Can you do 1/xT here???
    // For checks, fill also without efficiency correction
    fhChargedEtaNoCorr[centBin]->Fill(eta);
    fhChargedPhiNoCorr[centBin]->Fill(phi);
    fhChargedPtNoCorr[centBin]->Fill(pT);
    fhChargedXtNoCorr[centBin]->Fill(xT);
}
//______________________________________________________________________________
void AliJXtHistos::FillIsolatedHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin){
    // Fill all inclusive histograms
    fhIsolatedChargedEta[centBin]->Fill(eta, effCorr);
    fhIsolatedChargedPhi[centBin]->Fill(phi, effCorr);
    fhIsolatedChargedPt[centBin]->Fill(pT, effCorr);
    fhIsolatedChargedXt[centBin]->Fill(xT, effCorr);
    fhInvariantIsolatedChargedPt[centBin]->Fill(pT, pT > 0. ? effCorr/pT : 0.); // Can you do 1/pT here???
    fhInvariantIsolatedChargedXt[centBin]->Fill(xT, xT > 0. ? effCorr/xT : 0.); // Can you do 1/xT here???
}
//______________________________________________________________________________
void AliJXtHistos::CreateXtHistos(){
    // all the histograms of xT analysis
    
	MakeDirectory("xT");
	TH1::SetDefaultSumw2(kTRUE);
	cout << "GetDefaultSumw2() = " << TH1::GetDefaultSumw2() << endl;
    
	// === pT binning
	int nBins=150;
	double logBinsPt[nBins+1], limL=0.1, limH=100;
	double logBW = (log(limH)-log(limL))/nBins;
	for(int ij=0;ij<=nBins;ij++) logBinsPt[ij]=limL*exp(ij*logBW);

	// === xT binning
	int nBinsXt=200;
	double logBinsXt[nBinsXt+1];
	double xTlimL = 1e-5, xTlimH = 1.0, xTlogBW = (log(xTlimH)-log(xTlimL))/nBinsXt;
	for(int ij=0;ij<=nBinsXt;ij++) logBinsXt[ij]=xTlimL*exp(ij*xTlogBW);
    
    // === create histos
	for (int hic = 0;hic < fCard->GetNoOfBins(kCentrType);hic++) {
		float b1 = fCard->GetBinBorder(kCentrType, hic);
		float b2 = fCard->GetBinBorder(kCentrType, hic + 1);
        
		fhChargedPt[hic] = new TH1D(Form("hChargedPt%d",hic), Form("C: %2.0f-%2.0f%%",b1,b2), nBins, logBinsPt );
		fhChargedPt[hic]->Sumw2();
		fhistoList->Add(fhChargedPt[hic]);
        
		fhInvariantChargedPt[hic] = new TH1D(Form("hInvariantChargedPt%d",hic), Form("C: %2.0f-%2.0f%%",b1,b2), nBins, logBinsPt );
		fhInvariantChargedPt[hic]->Sumw2();
		fhistoList->Add(fhInvariantChargedPt[hic]);
        
		fhChargedPtNoCorr[hic] = new TH1D(Form("hChargedPtNoCorr%d",hic), Form("C: %2.0f-%2.0f%%",b1,b2), nBins, logBinsPt );
		fhChargedPtNoCorr[hic]->Sumw2();
		fhistoList->Add(fhChargedPtNoCorr[hic]);
        
		fhIsolatedChargedPt[hic] = new TH1D(Form("hIsolatedChargedPt%d",hic), Form("C: %2.0f-%2.0f%%",b1,b2), nBins, logBinsPt );
		fhIsolatedChargedPt[hic]->Sumw2();
		fhistoList->Add(fhIsolatedChargedPt[hic]);
        
		fhInvariantIsolatedChargedPt[hic] = new TH1D(Form("hInvariantIsolatedChargedPt%d",hic), Form("C: %2.0f-%2.0f%%",b1,b2), nBins, logBinsPt );
		fhInvariantIsolatedChargedPt[hic]->Sumw2();
		fhistoList->Add(fhInvariantIsolatedChargedPt[hic]);
        
        fhChargedXt[hic] = new TH1D(Form("hChargedXt%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), nBinsXt, logBinsXt );
        fhChargedXt[hic]->Sumw2();
        fhistoList->Add( fhChargedXt[hic] );
        
        fhInvariantChargedXt[hic] = new TH1D(Form("hInvariantChargedXt%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), nBinsXt, logBinsXt );
        fhInvariantChargedXt[hic]->Sumw2();
        fhistoList->Add( fhInvariantChargedXt[hic] );
        
        fhChargedXtNoCorr[hic] = new TH1D(Form("hChargedXtNoCorr%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), nBinsXt, logBinsXt );
        fhChargedXtNoCorr[hic]->Sumw2();
        fhistoList->Add( fhChargedXtNoCorr[hic] );
        
        fhIsolatedChargedXt[hic] = new TH1D(Form("hIsolatedChargedXt%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), nBinsXt, logBinsXt );
        fhIsolatedChargedXt[hic]->Sumw2();
        fhistoList->Add( fhIsolatedChargedXt[hic] );
        
        fhInvariantIsolatedChargedXt[hic] = new TH1D(Form("hInvariantIsolatedChargedXt%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), nBinsXt, logBinsXt );
        fhInvariantIsolatedChargedXt[hic]->Sumw2();
        fhistoList->Add( fhInvariantIsolatedChargedXt[hic] );
        
        fhChargedEta[hic] = new TH1D(Form("hChargedEta%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),100,-1.0,1.0);
        fhChargedEta[hic]->Sumw2();
        fhistoList->Add(fhChargedEta[hic]);
        
        fhChargedEtaNoCorr[hic] = new TH1D(Form("hChargedEtaNoCorr%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),100,-1.0,1.0);
        fhChargedEtaNoCorr[hic]->Sumw2();
        fhistoList->Add(fhChargedEtaNoCorr[hic]);
        
        fhIsolatedChargedEta[hic] = new TH1D(Form("hIsolatedChargedEta%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),100,-1.0,1.0);
        fhIsolatedChargedEta[hic]->Sumw2();
        fhistoList->Add(fhIsolatedChargedEta[hic]);
        
        fhChargedPhi[hic] = new TH1D(Form("hChargedPhi%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),128,-0.2,6.4);
        fhChargedPhi[hic]->Sumw2();
        fhistoList->Add(fhChargedPhi[hic]);
        
        fhChargedPhiNoCorr[hic] = new TH1D(Form("hChargedPhiNoCorr%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),128,-0.2,6.4);
        fhChargedPhiNoCorr[hic]->Sumw2();
        fhistoList->Add(fhChargedPhiNoCorr[hic]);
        
        fhIsolatedChargedPhi[hic] = new TH1D(Form("hIsolatedChargedPhi%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2),128,-0.2,6.4);
        fhIsolatedChargedPhi[hic]->Sumw2();
        fhistoList->Add(fhIsolatedChargedPhi[hic]);
        
        fhZVert[hic] = new TH1D(Form("hZVert%d",hic),Form("C: %2.0f-%2.0f%%",b1,b2), 120, -30., 30.);
        fhistoList->Add(fhZVert[hic]);
	}

	fhZVertRaw  = new TH1D("hZVertRaw","vertex 0", 120, -30., 30.);
	fhistoList->Add(fhZVertRaw);
	fhCentr          = new TH1D("hCentr","centrality", 101, -0.5, 100.5);
	fhistoList->Add(fhCentr);
	fhiCentr         = new TH1D("hiCentr","centrality",10, -0.5, 9.5);
	fhistoList->Add(fhiCentr);
	fhEventPerRun = new TH1D("hEventPerRun","log(eve)/run",200, 0, 30.0);
	fhistoList->Add(fhEventPerRun);

	//------------------ for Abs Norm FK --------------------------------
	double   binsVertexMult[] = {0,1,2,3,4,5,10000};
	int   nBinsVertexMult  = sizeof(binsVertexMult)/sizeof(double)-1;
	double binsVertexZ[]    = {-10,-6,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,6,10};
	int   nBinsVertexZ   = sizeof(binsVertexZ)/sizeof(double)-1;
	fhVertexZTriggVtx = new TH2D("hVertexZTriggVtx","Vertex counts", nBinsVertexMult, binsVertexMult, nBinsVertexZ, binsVertexZ);
	fhistoList->Add(fhVertexZTriggVtx);

	sprintf(fhtit, "Mean activity inside cone");
	sprintf(fhname, "hActivity");
	fhConeActivity = new TProfile(fhname, fhtit, nBins, logBinsPt );
	fhistoList->Add(fhConeActivity);

	sprintf(fhtit, "Mean pion activity inside cone isolated");
	sprintf(fhname, "hActivityIsolated");
	fhConeActivityIsolated = new TProfile(fhname, fhtit, nBins, logBinsPt );
	fhistoList->Add(fhConeActivityIsolated);

	JumpToDefaultDirectory();
}
