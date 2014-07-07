void SetStyles(TH1 *histo,int marker, int color, char *name){
  //histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetName(name);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
Int_t colors[] = {0,TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};

Float_t nSecondariesTrMultScaledSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesTrMultScaledData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesClMultScaledSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesClMultScaledData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesShortTrMultScaledSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesShortTrMultScaledData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesShortClMultScaledSim[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
Float_t nSecondariesShortClMultScaledData[] = {0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0};
void WriteLatex();
Float_t secondaryCorrEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryCorrPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorEmcal[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
Float_t secondaryErrorPhos[20] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};

void PlotSecondariesCorr(TString filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root",TString datafilename = "rootFiles/LHC10hPass2/Et.ESD.realPbPb.EMCal.LHC10hPass2.Run139465.root", Bool_t effCorr = kTRUE){
    TString detector = "";
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

    if(filename.Contains("EMC")){
      detector = "Emcal";
    }
    else{
      detector = "Phos";
    }
    TString tag = "";
    if(!effCorr) tag = "NoEffCorr";

    TString outname1 = "/tmp/SecondariesVsTrMult"+detector+".eps";
    TString outname3 = "/tmp/SecondariesVsTrMultRatio"+detector+".png";
    TString outname2 = "/tmp/SecondariesVsClMult"+detector+".eps";
    ofstream myfile;
    TString textfilename = "Secondaries"+detector+tag+".dat";
    myfile.open (textfilename.Data());
    ofstream myfile2;
    TString textfilename2 = "Secondaries"+detector+tag+"Short.dat";
    myfile2.open (textfilename2.Data());

    TFile *fdata = TFile::Open(datafilename, "READ");
    TList *ldata = dynamic_cast<TList*>(fdata->Get("out1"));
    TH3F *fHistCentVsNchVsNclData = ldata->FindObject("fHistCentVsNchVsNclReco");
    fHistCentVsNchVsNclData->SetName("fHistCentVsNchVsNclRecoData");

    TFile *f = TFile::Open(filename, "READ");
    TList *l = dynamic_cast<TList*>(f->Get("out1"));
    TH2F *fHistSecondariesVsNch;
    TH2F *fHistSecondariesVsNcl;
    TH2F *fHistSecondariesOutOfAccEffCorrVsNch;
    TH2F *fHistSecondariesDetectorCoverEffCorrVsNch;
    if(effCorr){
      fHistSecondariesVsNch =(TH2F *) l->FindObject("fHistSecondariesEffCorrVsNch");
      fHistSecondariesVsNcl =(TH2F *) l->FindObject("fHistSecondariesEffCorrVsNcl");
      fHistSecondariesOutOfAccEffCorrVsNch =(TH2F *) l->FindObject("fHistSecondariesOutOfAccEffCorrVsNch");
      fHistSecondariesDetectorCoverEffCorrVsNch =(TH2F *) l->FindObject("fHistSecondariesDetectorCoverEffCorrVsNch");
    }
    else{
      fHistSecondariesVsNch =(TH2F *) l->FindObject("fHistSecondariesVsNch");
      fHistSecondariesVsNcl =(TH2F *) l->FindObject("fHistSecondariesVsNcl");
      fHistSecondariesOutOfAccEffCorrVsNch =(TH2F *) l->FindObject("fHistSecondariesOutOfAccEffCorrVsNch");
      fHistSecondariesDetectorCoverEffCorrVsNch =(TH2F *) l->FindObject("fHistSecondariesDetectorCoverEffCorrVsNch");
    }
    TH3F *fHistCentVsNchVsNcl = l->FindObject("fHistCentVsNchVsNcl");

    TCanvas *c1 = new TCanvas("c1","Secondaries vs N_{ch}",600,400);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetBorderSize(0);
    c1->SetFillColor(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->SetLogz();
    //fHistSecondariesVsNch->GetXaxis()->SetRange(1,fHistSecondariesVsNch->GetXaxis()->FindBin(10.0));
    fHistSecondariesVsNch->Draw("colz");
    TH1D *prof1 = fHistSecondariesVsNch->ProfileX("prof1");
    prof1->Draw("same");
    TCanvas *c1a = new TCanvas("c1a","In detector cover and out of acceptance",600,400);
    c1a->SetTopMargin(0.02);
    c1a->SetRightMargin(0.02);
    c1a->SetBottomMargin(0.153226);
    c1a->SetBorderSize(0);
    c1a->SetFillColor(0);
    c1a->SetFillColor(0);
    c1a->SetBorderMode(0);
    c1a->SetFrameFillColor(0);
    c1a->SetFrameBorderMode(0);
    TH1D *profNCh = fHistSecondariesVsNch->ProfileY("profNCh");
    profNCh->Draw("");
    TH1D *profNChOutOfAcc = fHistSecondariesOutOfAccEffCorrVsNch->ProfileY("profNChOutOfAcc");
    TH1D *profNChDetCover = fHistSecondariesDetectorCoverEffCorrVsNch->ProfileY("profNChDetCover");
    //profNChDetCover->Draw("same");
    //profNChOutOfAcc->Draw("same");
    cout<<"Max "<<profNChDetCover->GetMaximum()<<" "<<profNChOutOfAcc->GetMaximum()<<" "<<profNCh->GetMaximum()<<endl;
    profNChDetCover->SetLineColor(TColor::kRed);
    profNChOutOfAcc->SetLineColor(TColor::kGreen+3);
    profNCh->GetXaxis()->SetTitle("N_{Ch}");
    profNCh->GetYaxis()->SetTitle("E_{T}^{Secondaries}");
    profNCh->GetXaxis()->SetLabelSize(0.06);
    profNCh->GetXaxis()->SetTitleSize(0.06);
    profNCh->GetYaxis()->SetLabelSize(0.06);
    profNCh->GetYaxis()->SetTitleSize(0.06);
    profNCh->GetYaxis()->SetTitleOffset(0.6);
    TLatex *texDet = new TLatex(50,6,detector.Data());
    texDet->SetTextSize(0.0752688);
    texDet->Draw();
    c1a->SaveAs(outname1.Data());
    TCanvas *c1b = new TCanvas("c1b","Fraction in detector cover, out of acceptance",600,400);
    c1b->SetTopMargin(0.02);
    c1b->SetRightMargin(0.02);
    c1b->SetBottomMargin(0.153226);
    c1b->SetBorderSize(0);
    c1b->SetFillColor(0);
    c1b->SetFillColor(0);
    c1b->SetBorderMode(0);
    c1b->SetFrameFillColor(0);
    c1b->SetFrameBorderMode(0);
//     profNCh->Draw("");
//     profNChDetCover->Draw("same")
//       profNChOutOfAcc->Draw("same");
    TH1D *profNChDetCoverRatio = profNChDetCover->Clone("profNChDetCoverRatio");
    TH1D *profNChOutOfAccRatio = profNChOutOfAcc->Clone("profNChOutOfAccRatio");
    profNChDetCoverRatio->Divide(profNCh);
    profNChOutOfAccRatio->Divide(profNCh);
    profNChDetCoverRatio->Draw("l");
    profNChOutOfAccRatio->Draw("same l");
    int rebin = 4;
    profNChDetCoverRatio->Rebin(rebin);
    profNChOutOfAccRatio->Rebin(rebin);
    profNChDetCoverRatio->Scale(1.0/rebin);
    profNChOutOfAccRatio->Scale(1.0/rebin);
    profNChDetCoverRatio->SetMaximum(0.2);
    c1b->SaveAs(outname3.Data());

    TCanvas *c2 = new TCanvas("c2","Secondaries vs NCl",600,400);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.02);
    c2->SetBorderSize(0);
    c2->SetFillColor(0);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetFrameFillColor(0);
    c2->SetFrameBorderMode(0);
    c2->SetLogz();
    fHistSecondariesVsNcl->Draw("colz");
    //fHistSecondariesVsNcl->GetXaxis()->SetRange(1,fHistSecondariesVsNcl->GetXaxis()->FindBin(10.0));
    TH1D *prof1 = fHistSecondariesVsNcl->ProfileX("prof2");
    prof1->Draw("same");
    TCanvas *c2a = new TCanvas("c2a","average secondaries vs NCl",600,400);
    c2a->SetTopMargin(0.02);
    c2a->SetRightMargin(0.02);
    c2a->SetBottomMargin(0.153226);
    c2a->SetBorderSize(0);
    c2a->SetFillColor(0);
    c2a->SetFillColor(0);
    c2a->SetBorderMode(0);
    c2a->SetFrameFillColor(0);
    c2a->SetFrameBorderMode(0);
    TH1D *profNCl = fHistSecondariesVsNcl->ProfileY("profNCl");
    profNCl->Draw("same");
    profNCl->GetXaxis()->SetTitle("N_{Cl}");
    profNCl->GetYaxis()->SetTitle("E_{T}^{Secondaries}");
    profNCl->GetXaxis()->SetLabelSize(0.06);
    profNCl->GetXaxis()->SetTitleSize(0.06);
    profNCl->GetYaxis()->SetLabelSize(0.06);
    profNCl->GetYaxis()->SetTitleSize(0.06);
    profNCl->GetYaxis()->SetTitleOffset(0.6);
    texDet->Draw();
    c2a->SaveAs(outname2.Data());

    TObjArray trackmultiplicity(20);
    TObjArray clustermultiplicity(20);
    TObjArray trackmultiplicityData(20);
    TObjArray clustermultiplicityData(20);
    int nbinsChMult = fHistCentVsNchVsNcl->GetYaxis()->GetNbins();
    int nbinsClMult = fHistCentVsNchVsNcl->GetZaxis()->GetNbins();
    fHistCentVsNchVsNcl->GetXaxis()->SetTitle("Cent Bin");
    fHistCentVsNchVsNcl->GetYaxis()->SetTitle("N_{tr}");
    fHistCentVsNchVsNcl->GetZaxis()->SetTitle("N_{cl}");
    for(int cb=0;cb<20;cb++){
      //x axis is centrality
      //y axis is charged track multiplicity
      //z axis is cluster multiplicity
      fHistCentVsNchVsNcl->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicity[cb] = fHistCentVsNchVsNcl->Project3D("y");
      clustermultiplicity[cb] = fHistCentVsNchVsNcl->Project3D("z");
      SetStyles((TH1*)trackmultiplicity[cb],markers[cb],colors[cb],Form("tr%i",cb));
      SetStyles((TH1*)clustermultiplicity[cb],markers[cb],colors[cb],Form("cl%i",cb));
      fHistCentVsNchVsNclData->GetXaxis()->SetRange(cb+1,cb+1);
      trackmultiplicityData[cb] = fHistCentVsNchVsNclData->Project3D("y");
      clustermultiplicityData[cb] = fHistCentVsNchVsNclData->Project3D("z");
      SetStyles((TH1*)trackmultiplicityData[cb],markers[cb],colors[cb],Form("tr%i",cb));
      SetStyles((TH1*)clustermultiplicityData[cb],markers[cb],colors[cb],Form("cl%i",cb));
    }
    TCanvas *c3 = new TCanvas("c3","<E_{T}> Neutral Removed",600,400);
    c3->SetTopMargin(0.02);
    c3->SetRightMargin(0.02);
    c3->SetBorderSize(0);
    c3->SetFillColor(0);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetFrameFillColor(0);
    c3->SetFrameBorderMode(0);
    c3->SetLogy();
    ((TH1*)trackmultiplicity[0])->Draw();
    for(int cb=1;cb<20;cb++){
//       cout<<"cb "<<cb<<" tr "<< ((TH1*)trackmultiplicity[cb])->GetMean();
//       cout<<" cluster "<< ((TH1*)clustermultiplicity[cb])->GetMean()<<endl;
      ((TH1*)trackmultiplicity[cb])->Draw("same");
    }
    TCanvas *c4 = new TCanvas("c4","<E_{T}> Neutral Removed",600,400);
    c4->SetTopMargin(0.02);
    c4->SetRightMargin(0.02);
    c4->SetBorderSize(0);
    c4->SetFillColor(0);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetFrameFillColor(0);
    c4->SetFrameBorderMode(0);
    c4->SetLogy();
    ((TH1*)clustermultiplicity[0])->Draw();
    for(int cb=1;cb<20;cb++){
//       cout<<"cb "<<cb<<" "<< ((TH1*)trackmultiplicity[cb])->GetMean()<<endl;
      ((TH1*)clustermultiplicity[cb])->Draw("same");
    }
    cout<<" Sim cluster track data cluster track "<<endl;
      //for fewer centrality bins
    int currentShortCB = 0;
    Float_t neventsShort[] = {0,0,0,0};
    //profNCh, profNCl
    for(int cb=0;cb<19;cb++){

      int nbins = ((TH1*)clustermultiplicity[cb])->GetNbinsX();
      Float_t nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
	float neventsBinn = ((TH1*)clustermultiplicity[cb])->GetBinContent(binn);
	float mult= ((TH1*)clustermultiplicity[cb])->GetBinCenter(binn);
	float meanet = profNCl->GetBinContent(profNCl->FindBin(mult));
	nSecondariesClMultScaledSim[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nSecondariesShortClMultScaledSim[currentShortCB] += neventsBinn*meanet;
	neventsShort[0] +=neventsBinn;
      }
      if(nevents>0) nSecondariesClMultScaledSim[cb] = nSecondariesClMultScaledSim[cb]/nevents;



      nbins = ((TH1*)trackmultiplicity[cb])->GetNbinsX();
      nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
	float neventsBinn = ((TH1*)trackmultiplicity[cb])->GetBinContent(binn);
	float mult= ((TH1*)trackmultiplicity[cb])->GetBinCenter(binn);
	float meanet = profNCh->GetBinContent(profNCh->FindBin(mult));
	nSecondariesTrMultScaledSim[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nSecondariesShortTrMultScaledSim[currentShortCB] += neventsBinn*meanet;
	neventsShort[1] +=neventsBinn;
      }
      if(nevents>0) nSecondariesTrMultScaledSim[cb] = nSecondariesTrMultScaledSim[cb]/nevents;

      nbins = ((TH1*)clustermultiplicity[cb])->GetNbinsX();
      nevents = 0.0;
      float meanmult = 0.0;
      //cout<<"Data cluster scaled  = (";
      for(int binn = 1;binn<=nbins;binn++){
	float neventsBinn = ((TH1*)clustermultiplicity[cb])->GetBinContent(binn);
	float mult= ((TH1*)clustermultiplicity[cb])->GetBinCenter(binn);
	float meanet = profNCl->GetBinContent(profNCl->FindBin(mult));
	nSecondariesClMultScaledData[cb] += neventsBinn*meanet;
	meanmult += neventsBinn*mult;
	nevents +=neventsBinn;
	//if(neventsBinn>0) cout<<"("<<mult<<") "<<neventsBinn<<"*"<<meanet<<"+";
	nSecondariesShortClMultScaledData[currentShortCB] += neventsBinn*meanet;
 	neventsShort[2] +=neventsBinn;
      }
      if(nevents>0) meanmult = meanmult/nevents;
      //nSecondariesClMultScaledData[cb] = nSecondariesClMultScaledData[cb]/nevents;
      //cout<<")/"<<nevents<<endl;
      //cout<<"cb "<<cb<<" = "<<nSecondariesClMultScaledData[cb]<<"/"<<nevents;
      if(nevents>0) nSecondariesClMultScaledData[cb] = nSecondariesClMultScaledData[cb]/nevents;//<<" = ";
      //cout<<" " << nSecondariesClMultScaledData[cb];
      //cout<<" mean "<<((TH1*)clustermultiplicityData[cb])->GetMean()<<" mean2 "<<meanmult<<" mean method "<< profNCl->GetBinContent(profNCl->FindBin(meanmult));
      //cout<<" mean method "<< profNCl->GetBinContent(profNCl->FindBin(((TH1*)trackmultiplicityData[cb])->GetMean()));
      //cout<<endl;

      nbins = ((TH1*)trackmultiplicity[cb])->GetNbinsX();
      nevents = 0.0;
      for(int binn = 1;binn<=nbins;binn++){
	float neventsBinn = ((TH1*)trackmultiplicityData[cb])->GetBinContent(binn);
	float mult= ((TH1*)trackmultiplicityData[cb])->GetBinCenter(binn);
	float meanet = profNCh->GetBinContent(profNCh->FindBin(mult));
	nSecondariesTrMultScaledData[cb] += neventsBinn*meanet;
	nevents +=neventsBinn;
	nSecondariesShortTrMultScaledData[currentShortCB] += neventsBinn*meanet;
	neventsShort[3] +=neventsBinn;
      }

//       cout<<"cb "<<cb<<" = "<<nSecondariesTrMultScaledData[cb]<<"/"<<nevents;
      if(nevents>0) nSecondariesTrMultScaledData[cb] = nSecondariesTrMultScaledData[cb]/nevents;
//       cout<<" " << nSecondariesTrMultScaledData[cb];
//       cout<<" mean "<<((TH1*)trackmultiplicityData[cb])->GetMean()<<" mean method "<< profNCh->GetBinContent(profNCh->FindBin(((TH1*)trackmultiplicityData[cb])->GetMean()));
//       cout<<endl<<endl;


      float meanSim = ( nSecondariesTrMultScaledSim[cb] +  nSecondariesClMultScaledSim[cb] ) /2.0;
      float errSim =  TMath::Abs( nSecondariesTrMultScaledSim[cb] -  nSecondariesClMultScaledSim[cb] ) /2.0;
      float mean = ( nSecondariesTrMultScaledData[cb] +  nSecondariesClMultScaledData[cb] ) /2.0;
      float err =  TMath::Abs( nSecondariesTrMultScaledData[cb] -  nSecondariesClMultScaledData[cb] ) /2.0;
      myfile<<Form("%2.3f %2.3f",mean,err)<<endl;

      cout<<"cb "<<cb<<" ET from secondaries "<< nSecondariesClMultScaledSim[cb];// <<endl;
      cout<<" "<< nSecondariesTrMultScaledSim[cb];
      cout<<" "<<meanSim << "+/-"<<errSim;
      cout<<" "<< nSecondariesClMultScaledData[cb];
      cout<<" "<< nSecondariesTrMultScaledData[cb];
      cout<<" "<<mean << "+/-"<<err;
      cout<<endl;

      //return;
      if(cb<2 || (cb+1)%2==0){//For combined bins, 
	if(neventsShort[0]>0)nSecondariesShortClMultScaledSim[currentShortCB] = nSecondariesShortClMultScaledSim[currentShortCB]/neventsShort[0];
	if(neventsShort[1]>0)nSecondariesShortTrMultScaledSim[currentShortCB] = nSecondariesShortTrMultScaledSim[currentShortCB]/neventsShort[1];
	if(neventsShort[2]>0)nSecondariesShortClMultScaledData[currentShortCB] = nSecondariesShortClMultScaledData[currentShortCB]/neventsShort[2];
	if(neventsShort[3]>0)nSecondariesShortTrMultScaledData[currentShortCB] = nSecondariesShortTrMultScaledData[currentShortCB]/neventsShort[3];
	cout<<"cb Short "<<currentShortCB<<" ET from secondaries "<< nSecondariesShortClMultScaledSim[currentShortCB];// <<endl;
	cout<<" "<< nSecondariesShortTrMultScaledSim[currentShortCB];

	meanSim = ( nSecondariesShortTrMultScaledSim[currentShortCB] +  nSecondariesShortClMultScaledSim[currentShortCB] ) /2.0;
	errSim =  TMath::Abs( nSecondariesShortTrMultScaledSim[currentShortCB] -  nSecondariesShortClMultScaledSim[currentShortCB] ) /2.0;
	cout<<" "<<meanSim << "+/-"<<errSim;

	cout<<" "<< nSecondariesShortClMultScaledData[currentShortCB];
	cout<<" "<< nSecondariesShortTrMultScaledData[currentShortCB];

	for(int k=0;k<4;k++){neventsShort[k] = 0;}

	mean = ( nSecondariesShortTrMultScaledData[currentShortCB] +  nSecondariesShortClMultScaledData[currentShortCB] ) /2.0;
	err =  TMath::Abs( nSecondariesShortTrMultScaledData[currentShortCB] -  nSecondariesShortClMultScaledData[currentShortCB] ) /2.0;
	myfile2<<Form("%2.3f %2.3f",mean,err)<<endl;
	//myfile2<<mean<<" "<<err<<endl;
	cout<<" "<<mean << "+/-"<<err;
	cout<<endl;

	currentShortCB++;
	cout<<"cb "<<cb<<" incrementing short CB "<<currentShortCB<<endl<<endl;
      }

    }
    myfile.close();
    myfile2.close();

    WriteLatex();
}

void WriteLatex(){
  TString detector = "Emcal";
    string inline;
    TString secondaryInfileName = "Secondaries"+detector+".dat";
    ifstream mysecondaryfile3 (secondaryInfileName.Data());
    Float_t value = 0;
    Float_t error = 0;
    Int_t i=0;
    if (mysecondaryfile3.is_open()){
      while ( mysecondaryfile3.good() )
	{
	  getline (mysecondaryfile3,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    secondaryCorrEmcal[i] = value;
	    secondaryErrorEmcal[i] = error;
	  }
	  i++;
	}
        mysecondaryfile3.close();
    }

    detector = "Phos";
    secondaryInfileName = "Secondaries"+detector+".dat";
    ifstream mysecondaryfile4 (secondaryInfileName.Data());
    Float_t value = 0;
    Float_t error = 0;
    Int_t i=0;
    if (mysecondaryfile4.is_open()){
      while ( mysecondaryfile4.good() )
	{
	  getline (mysecondaryfile4,inline);
	  istringstream tmp(inline);
	  tmp >> value;
	  tmp >> error;
	  if(i<20){
	    secondaryCorrPhos[i] = value;
	    secondaryErrorPhos[i] = error;
	  }
	  i++;
	}
        mysecondaryfile4.close();
    }
    ofstream myfile3;
    myfile3.open ("Secondaries.tex");

    for(int i=0;i<20;i++){
      TString line = Form("%i-%i & %2.3f $\\pm$ %2.3f & %2.3f $\\pm$ %2.3f \\\\",i*5,(i+1)*5,secondaryCorrPhos[i],secondaryErrorPhos[i],secondaryCorrEmcal[i],secondaryErrorEmcal[i]);
      myfile3<<line.Data()<<endl;
    }
    myfile3.close();


}
