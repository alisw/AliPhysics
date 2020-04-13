
void ExtractScaledCov(const char* listname = "list_tree.txt", const char* outfile = "scaledcov.root", int maxfiles = -1, const char* tag="")
{
    TDatime d;
    TString ids = TString("file generated ");
    ids += d.AsString();
    ids += " inputlist=";
    ids += listname;
    ids += " outputfile=";
    ids += outfile;
    ids += " maxfiles=";
    ids += maxfiles;
    ids += " tag=";
    ids += tag;
    
    TObjString oids(ids);

    TF1* fit = new TF1("linfit","[0]+x*[1]+x*x*[2]",0,1);
    TF1* normres = new TF1("normres","x*[0]+x*x*[1]",0,1);

    //tmp    TGrid::Connect("alien://");
    //gROOT->LoadMacro("tmp_function.C+");
    gROOT->ProcessLineSync(".x tmp_function.C+");    
    gROOT->ProcessLineSync("tmp_function_ptr = normres");
      
    TChain* hpt = AliXRDPROOFtoolkit::MakeChain(listname,"highPt",0,maxfiles);   
    
    hpt->SetAlias("cut0","(esdTrack->GetTPCchi2()/esdTrack->GetTPCclusters(0)<4.0)&&(abs(esdTrack.fZ)<2.0)&&(abs(esdTrack.Eta())<0.8)&&(esdTrack.IsOn(64))&&(esdTrack.IsOn(4))&&(esdTrack.HasPointOnITSLayer(0)||esdTrack.HasPointOnITSLayer(1))&&(esdTrack->fITSchi2/esdTrack->fITSncls<36.0)&&(esdTrack->GetTPCnclsS()/esdTrack->GetTPCclusters(0)<0.4)&&(chi2TPCInnerC<36.0)&&(esdTrack->GetLengthInActiveZone(1,3,220,2)>(130.0-TMath::Power(TMath::Abs(esdTrack->GetSigned1Pt()),1.5)))&&(esdTrack->GetTPCClusterInfo(2,1)>0.85*(130.0-TMath::Power(TMath::Abs(esdTrack->GetSigned1Pt()),1.5)))&&(esdTrack->GetTPCNcls()>0.85*(130.0-TMath::Power(TMath::Abs(esdTrack->GetSigned1Pt()),1.5))) &&esdTrack->Pt()>7&&sqrt(esdTrack->GetSigma1Pt2())<0.1&&esdTrack->Pt()<50");
    
    TCanvas* c1 = new TCanvas();
    c1->cd();
    // all cuts, pt Range 7-50 gev/c
    hpt->Draw("sqrt(esdTrack->GetSigma1Pt2()):abs(esdTrack->GetSigned1Pt())","cut0","hprof");
    cout<<gPad<<endl;
    cout<<gPad->GetPrimitive("htemp")<<endl;
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    
    fit->FixParameter(2,0);
    htemp->Fit(fit,"W","",0,0.1);
    fit->ReleaseParameter(2);
    htemp->Fit(fit,"W","",0,0.1);
    htemp->Draw();
    normres->SetParameters(fit->GetParameter(1),fit->GetParameter(2));    
    
    TCanvas* c2 = new TCanvas();
    c2->cd();
    hpt->Draw("sqrt(esdTrack->GetSigma1Pt2()):abs(esdTrack->GetSigned1Pt())>>hist1","cut0","COLZ");

    TCanvas* c3 = new TCanvas();
    c3->cd();
    hpt->Draw("sqrt(esdTrack->GetSigma1Pt2())-tmp_function(abs(esdTrack->GetSigned1Pt())):abs(esdTrack->GetSigned1Pt())>>hist2","cut0","COLZ");
    
    TCanvas* c4 = new TCanvas();
    c4->cd();
    hpt->Draw("sqrt(esdTrack->GetSigma1Pt2())-tmp_function(abs(esdTrack->GetSigned1Pt()))","cut0","COLZ");
    TH1F *htemp2 = (TH1F*)gPad->GetPrimitive("htemp");
    
    // write everything to file
    TFile* fout = TFile::Open(outfile,"RECREATE");
    c1->Write("c1_FIT");
    c2->Write("c2_RESOLUTION");
    c3->Write("c3_NORMALIZED");
    c4->Write("c4_PROJECTION");
    htemp2->Write("scaledcov");
    normres->Write("normres");
    fit->Write("fitfun");   
    oids.Write("INFO");
    fout->Close();
    cout<<ids<<endl;
}