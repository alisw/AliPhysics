Double_t ApplyCorrections(const char* datafile, const char* datatask, const char* corrfile, const char* idstring ,const char* outfile, const char* gifdir = 0)
//void ApplyCorrections()
{

// tmp setting
//     const char* datafile = "/lustre/alice/train/V005.PbPb/2010-11-28_0258.4230/mergedPeriods/data/PbPb/LHC10h.pass1/mknichel_dNdPtPbPb_gt_v0_c0.root";
//     const char* datatask = "mknichel_dNdPtPbPb_gt_v0_c0";
//     const char* corrfile = "/u/mknichel/alice/dNdPt/2010-11-28/corrMatr_gt_v0_c0.root";
//     const char* idstring = "gt_v0_c0";
//     const char* outfile = "/u/mknichel/alice/dNdPt/2010-11-28/finalSpectra_v0_c0.root";
//     const char* gifdir = "/u/mknichel/alice/dNdPt/2010-11-28";
    
    // settings vor zVertex cut (event and track level)
    Double_t zVert = 10.0;
    
    // setting on eta cut (track level)
    Double_t eta = 0.8;
    
    //load required libraries
    //load required libraries    
    gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG0 -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWG2 -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -I$ALICE_ROOT/PWG4 -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG3/muon -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/ANALYSIS/Tender");
    
  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");

  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libNet");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libANALYSIS");

    
    
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWG0base.so");    
    gSystem->Load("libPWG0dep"); 
    gSystem->Load("libPWG0selectors.so");
    

    // make plots nicer
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    
    // array for all histograms to be saves
    TObjArray* Hists = new TObjArray();
    
    // open file with correction matrices
    TFile *fcorr = TFile::Open(corrfile,"READ");
    if (!fcorr) return -2;
    
    // load data
    TFile* fdata = TFile::Open(datafile,"READ");
    if (!fdata) return -1;
    TList* ldata = dynamic_cast<TList*>(fdata->Get(datatask));
    if (!ldata) return -1;
    AlidNdPtAnalysisPbPb *obj = dynamic_cast<AlidNdPtAnalysisPbPb*>(ldata->FindObject("dNdPtAnalysisPbPb"));
    if (!obj) return -1;
    
    //Event statistics
    THnSparse *fRecEventHist2 = obj->GetRecEventHist2(); //reconstructed events	
    TH2D* h2RecEvent2All = fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2All");
    fRecEventHist2->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    TH2D* h2RecEvent2 = fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2");
    
    Double_t ReconstructedEvents = h2RecEvent2->Integral();
    Double_t ReconstructedEventsAll = h2RecEvent2All->Integral();
    
    Hists->Add(h2RecEvent2);
    Hists->Add(h2RecEvent2All);
  
    printf("=== Number of events from DATA                      %lf ===\n",ReconstructedEvents);
    printf("=== Number of events from DATA (before zVertex cut) %lf ===\n",ReconstructedEventsAll);        
  
    TH1D* h1ReconstructedEvents = new TH1D("h1ReconstructedEvents","h1ReconstructedEvents",1,0,1);
    TH1D* h1ReconstructedEventsAll = new TH1D("h1ReconstructedEventsAll","h1ReconstructedEventsAll",1,0,1);
    
    h1ReconstructedEvents->Fill(0,ReconstructedEvents);
    h1ReconstructedEvents->SetEntries(ReconstructedEvents);
    h1ReconstructedEventsAll->Fill(0,ReconstructedEventsAll);
    h1ReconstructedEventsAll->SetEntries(ReconstructedEventsAll);
        
    Hists->Add(h1ReconstructedEvents);
    Hists->Add(h1ReconstructedEventsAll);

     // retrieve tracks
    THnSparse* fRecTrackHist2 = obj->GetRecTrackHist2(2); //after all cuts (2)
    fRecTrackHist2->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01); //zVertex
    fRecTrackHist2->GetAxis(2)->SetRangeUser(-eta, eta-0.01); // eta   
   
    TH1D* h1RecTrackHist2_zv = fRecTrackHist2->Projection(0)->Clone("h1RecTrackHist2_zv");
    TH1D* h1RecTrackHist2_pt = fRecTrackHist2->Projection(1)->Clone("h1RecTrackHist2_pt");
    TH1D* h1RecTrackHist2_eta = fRecTrackHist2->Projection(2)->Clone("h1RecTrackHist2_eta");
   
    Hists->Add(h1RecTrackHist2_zv);
    Hists->Add(h1RecTrackHist2_pt);
    Hists->Add(h1RecTrackHist2_eta);

    // retrieve correction matrices for tracking efficiency (note different binning!)
    TH1D* h1TrackCorr_pt  = (TH1D*)fcorr->Get("h1TrackCorr_pt");  
    TH1D* h1TrackCorr_eta = (TH1D*)fcorr->Get("h1TrackCorr_eta");
        
    // retrieve correction matrices for secondaries (note different binning!)
    TH1D* h1SecCorr_pt  = (TH1D*)fcorr->Get("h1SecCorr_pt");  
    TH1D* h1SecCorr_eta = (TH1D*)fcorr->Get("h1SecCorr_eta");

    // create corrected spectra (as clone of raw data)
    TH1D* h1Corrected_pt  = h1RecTrackHist2_pt->Clone("h1Corrected_pt");
    TH1D* h1Corrected_eta  = h1RecTrackHist2_eta->Clone("h1Corrected_eta");

    // secondaries correction for pt spectrum
    for (int i=1; i <= h1Corrected_pt->GetNbinsX() ; i++) {
        Double_t pt = h1Corrected_pt->GetBinCenter(i);
        Double_t val = h1Corrected_pt->GetBinContent(i);
        Double_t err = h1Corrected_pt->GetBinError(i);
        if (pt >= 19) { pt = 19; }  // above 20 GeV corr matr have low statistics
        Double_t secCorr    = h1SecCorr_pt->GetBinContent(h1SecCorr_pt->FindBin(pt));
        Double_t secCorrErr = h1SecCorr_pt->GetBinError(h1SecCorr_pt->FindBin(pt));
        Double_t cval = val*secCorr;
        Double_t cerr = TMath::Sqrt(val*val*secCorrErr*secCorrErr + err*err*secCorr*secCorr);
        h1Corrected_pt->SetBinContent(i,cval);
        h1Corrected_pt->SetBinError(i,cerr);
    }

    // tracking efficiency correction pt spectrum
    for (int i=1; i <= h1Corrected_pt->GetNbinsX() ;i++) {
        Double_t pt = h1Corrected_pt->GetBinCenter(i);
        Double_t val = h1Corrected_pt->GetBinContent(i);
        Double_t err = h1Corrected_pt->GetBinError(i);
        if (pt >= 19) { pt = 19; } // above 20 GeV corr matr have low statistics
        Double_t effCorr    = h1TrackCorr_pt->GetBinContent(h1TrackCorr_pt->FindBin(pt));
        Double_t effCorrErr = h1TrackCorr_pt->GetBinError(h1TrackCorr_pt->FindBin(pt));
        
        Double_t cval = val*effCorr;
        Double_t cerr = TMath::Sqrt(val*val*effCorrErr*effCorrErr + err*err*effCorr*effCorr);
        h1Corrected_pt->SetBinContent(i,cval);
        h1Corrected_pt->SetBinError(i,cerr);
    }
    
    // for eta the correction is simpler because of same binning
    h1Corrected_eta->Multiply(h1SecCorr_eta);
    h1Corrected_eta->Multiply(h1TrackCorr_eta);
    
    Hists->Add(h1Corrected_pt);
    Hists->Add(h1Corrected_eta);

    // create final spectra (as clone of corrected data)
    TH1D* dNdPt   = h1Corrected_pt->Clone("dNdPt");
    TH1D* dNdEta  = h1Corrected_eta->Clone("dNdEta");
       
    // also uncorrected spectra (as clone of raw data)
    TH1D* dNdPt_raw   = h1RecTrackHist2_pt->Clone("dNdPt_raw");
    TH1D* dNdEta_raw  = h1RecTrackHist2_eta->Clone("dNdEta_raw");
       

//TF1 *fperi  = new TF1("fperi","1.00343-0.000608425*x-6.7038e-05*x*x",5.,40.);
//TF1 *fcent  = new TF1("cent","1.01074e+00-1.98127e-03*x-1.19903e-04*x*x",5.,40.);
TString id = TString(idstring);
if ( id.Contains("c0") || id.Contains("c5") ) {
  TF1 *fun  = new TF1("cent","1.01074e+00-1.98127e-03*x-1.19903e-04*x*x",5.,40.);    
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "CENTRAL pt-resolution correction used!" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
} else {
TF1 *fun  = new TF1("fperi","1.00343-0.000608425*x-6.7038e-05*x*x",5.,40.);
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "PERIPHERAL pt-resolution correction used!" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
}

    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt->GetNbinsX() ;i++) {
        Double_t pt = dNdPt->GetBinCenter(i);
        Double_t width = dNdPt->GetBinWidth(i);
        Double_t val = dNdPt->GetBinContent(i);
        Double_t err = dNdPt->GetBinError(i);
        Double_t corrPtResol = 1.0;
        corrPtResol = fun->Eval(pt);
        if (pt < 5.) { corrPtResol = 1.0; }
        Double_t cval = (val * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        Double_t cerr = (err * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        dNdPt->SetBinContent(i,cval);
        dNdPt->SetBinError(i,cerr);
    }
    // for dndeta again simpler
    dNdEta->Scale(1,"width");
    dNdEta->Scale(1./ReconstructedEvents);
    
    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_raw->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_raw->GetBinCenter(i);
        Double_t width = dNdPt_raw->GetBinWidth(i);
        Double_t val = dNdPt_raw->GetBinContent(i);
        Double_t err = dNdPt_raw->GetBinError(i);
        Double_t cval = val/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        Double_t cerr = err/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        dNdPt_raw->SetBinContent(i,cval);
        dNdPt_raw->SetBinError(i,cerr);
    }
    // for dndeta again simpler
    dNdEta_raw->Scale(1,"width");
    dNdEta_raw->Scale(1./ReconstructedEvents);
     
    dNdEta->SetMarkerStyle(21);
    dNdPt->SetMarkerStyle(21);
    dNdPt->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdEta->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    dNdPt_raw->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdEta_raw->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    Hists->Add(dNdEta_raw);
    Hists->Add(dNdPt_raw);
    Hists->Add(dNdEta);
    Hists->Add(dNdPt);
    
    
   // plot pictures and save to gifdir
    for (i=0; i < Hists->LastIndex(); i++) {    
        TCanvas* ctmp = PlotHist(Hists->At(i),idstring);
        if (gifdir && ctmp) {
            TString gif(gifdir);
            gif += '/';
            gif += ctmp->GetName();
            gif += ".gif";
            ctmp->SaveAs(gif.Data(),"gif");     
            delete ctmp;
        }
    }  

    // save all correction matrices and control histograms to file
    if (!outfile) { return; }
    TFile *out = TFile::Open(outfile,"RECREATE");
    Hists->Write();
    out->Close();    

    return ReconstructedEvents;

}

//___________________________________________________________________________
TCanvas* PlotHist(TObject* hobj, const char* label=0)
{
    TH1* h = dynamic_cast<TH1*>(hobj);
    if (!h) return 0;
    if (h->GetDimension() > 2) return 0;
    h->SetStats(0);
    if ( TString(h->GetName()).Contains("Events")) { h->SetStats(1); } 
    TString t(label);
    if (label) t += "_";
    t += h->GetName();
    h->SetTitle(t.Data());
    TCanvas* c = new TCanvas(t.Data(),t.Data());
    if (h->GetDimension() >= 1) {
        TString xlabel(h->GetXaxis()->GetTitle());
        if (xlabel.Contains("Pt")) { c->SetLogx();  c->SetLogy();  h->GetXaxis()->SetRangeUser(0.1 , 50.); }
        if (xlabel.Contains("p_{T}")) { c->SetLogx();  c->SetLogy();  h->GetXaxis()->SetRangeUser(0.1 , 50.); }        
    }
    if (h->GetDimension() == 2) {  
        TString ylabel(h->GetYaxis()->GetTitle());
        if (ylabel.Contains("Pt")) { c->SetLogy(); h->GetYaxis()->SetRangeUser(0.1 , 50.); }
        if (ylabel.Contains("p_{T}")) { c->SetLogy(); h->GetYaxis()->SetRangeUser(0.1 , 50.); }
        h->Draw("COLZ");
    }        
    if (h->GetDimension() == 1) {
        h->Draw();
    }
    return c;

}
