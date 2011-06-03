Double_t ApplyCorrections_PbPb(const char* datafile, const char* datatask, const char* corrfile, const char* idstring ,const char* outfile, const char* gifdir = 0)
//Double_t ApplyCorrections_PbPb()
{

// tmp setting
//      const char* datafile = "/lustre/alice/train/V006.PbPb/2011-03-15_0009.5917/mergedPeriods/PbPb/2.76ATeV/LHC10h.pass1/mknichel_dNdPtPbPb_TPCITS_VZERO1.root";     
//      const char* datatask = "mknichel_dNdPtPbPb_TPCITS_VZERO";
//      const char* corrfile = "/u/mknichel/alice/dNdPt_PbPb/2011-03-15/corrMatr_LHC10h8_TPCITS.root";     
//      const char* idstring = "c0";
//      const char* outfile = "/u/mknichel/alice/dNdPt_PbPb/2011-03-15/finalSpectra_LHC10h.pass1_TPCITS.root";     
//      const char* gifdir = "/u/mknichel/alice/dNdPt_PbPb/2011-03-15/LHC10h.pass1";
 
    
Int_t c_first = 1;
Int_t c_last = 11;

TString id = TString(idstring);
if ( id.Contains("c0-5") ) { c_first = c_last = 1; }
if ( id.Contains("c5-10") ) { c_first = c_last = 2; }
if ( id.Contains("c10-20") ) { c_first = c_last = 3; }
if ( id.Contains("c20-30") ) { c_first = c_last = 4; }
if ( id.Contains("c30-40") ) { c_first = c_last = 5; }
if ( id.Contains("c40-50") ) { c_first = c_last = 6; }
if ( id.Contains("c50-60") ) { c_first = c_last = 7; }
if ( id.Contains("c60-70") ) { c_first = c_last = 8; }
if ( id.Contains("c70-80") ) { c_first = c_last = 9; }
if ( id.Contains("c80-90") ) { c_first = c_last = 10; }
if ( id.Contains("c90-100") ) { c_first = c_last = 11; }
    
    // settings vor zVertex cut (event and track level)
    Double_t zVert = 10.0;
    
    // setting on eta cut (track level)
    Double_t eta = 0.8;
    
    //load required libraries
    //load required libraries    
    gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG0 -I$ALICE_ROOT/PWG1 -I$ALICE_ROOT/PWG2 -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -I$ALICE_ROOT/PWG4 -I$ALICE_ROOT/CORRFW -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWG3/muon -I$ALICE_ROOT/JETAN -I$ALICE_ROOT/ANALYSIS/Tender");
    
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
    
    // array for all histograms to be saved
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
    fRecEventHist2->GetAxis(2)->SetRange(c_first,c_last); // select centrality    
    TH2D* h2RecEvent2All = (TH2D*) fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2All");
    fRecEventHist2->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    TH2D* h2RecEvent2 = (TH2D*) fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2");
    
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
    fRecTrackHist2->GetAxis(3)->SetRange(c_first,c_last); // select centrality    
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
        if (pt >= 50) { pt = 49.5; }  // above 50 GeV corr matr have low statistics
        Double_t secCorr    = h1SecCorr_pt->GetBinContent(h1SecCorr_pt->FindBin(pt));
        Double_t secCorrErr = h1SecCorr_pt->GetBinError(h1SecCorr_pt->FindBin(pt));
        Double_t effCorr    = h1TrackCorr_pt->GetBinContent(h1TrackCorr_pt->FindBin(pt));
        Double_t effCorrErr = h1TrackCorr_pt->GetBinError(h1TrackCorr_pt->FindBin(pt));        
        Double_t corr    = effCorr*secCorr;
        Double_t corrErr = effCorr*secCorrErr + secCorr*effCorrErr; // errors are correlated        
        Double_t cval = val*corr;
        Double_t cerr = TMath::Sqrt(val*val*corrErr*corrErr + err*err*corr*corr);        
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
    
    TH1D* dNdPt_nores   = h1Corrected_pt->Clone("dNdPt_nores");    
       

//TF1 *fperi  = new TF1("fperi","1.00343-0.000608425*x-6.7038e-05*x*x",5.,40.);
//TF1 *fcent  = new TF1("cent","1.01074e+00-1.98127e-03*x-1.19903e-04*x*x",5.,40.);
TFile* fptcorr = TFile::Open("ptcorrPbPb_150511.root");

TF1 * fun = 0;

  TF1 *fcent  = new TF1("fcent","1.01074e+00-1.98127e-03*x-1.19903e-04*x*x",5.,50.);    
  TF1 *fperi  = new TF1("fperi","1.00343-0.000608425*x-6.7038e-05*x*x",5.,50.);
Double_t downscale = 1.;
            if (c_first != c_last) {
             cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
             cout << "WARNING: pt resolution correction error!" << endl;
             cout << " (works only for single centraliy bins) " << endl;
             cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
        }
        if (c_first == 1) {
		fun = (TF1*) fptcorr->Get("ptcorr_c0");
		downscale = 0.8;
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "0-5% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	
        } else if (c_first == 2) { 
        	downscale = 0.8;
		fun =(TF1*) fptcorr->Get("ptcorr_c5");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "5-10% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	        
        } else if (c_first == 3) { 
        	downscale = 0.8;
		fun =(TF1*) fptcorr->Get("ptcorr_c10");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "10-20% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	        
        } else if (c_first == 4) { 
        	downscale = 0.9;
		fun =(TF1*) fptcorr->Get("ptcorr_c20");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "20-30% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	        
        } else if (c_first == 5) { 
		fun =(TF1*) fptcorr->Get("ptcorr_c30");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "30-40% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	        
        } else if (c_first == 6) { 
		fun =(TF1*) fptcorr->Get("ptcorr_c40");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "40-50% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	        
        } else if (c_first == 7) { 
		fun =(TF1*) fptcorr->Get("ptcorr_c50");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "50-60% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	                
	} else if (c_first == 8) { 
		fun =(TF1*) fptcorr->Get("ptcorr_c60");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;
 		cout << "60-70% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++" << endl;   	                   
        } else {
        	fun =(TF1*) fptcorr->Get("ptcorr_c70");
  		cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;  		
 		cout << "70-80% central pt-resolution correction used!" << endl;
 		cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;   
        }    
    

    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange    
    
   for (int i=1; i <= dNdPt->GetNbinsX() ;i++) {
        Double_t pt = dNdPt->GetBinCenter(i);
        Double_t width = dNdPt->GetBinWidth(i);
        Double_t val = dNdPt->GetBinContent(i);
        Double_t err = dNdPt->GetBinError(i);
        Double_t corrPtResol = 1.0;     
        corrPtResol = 1.-((1.-fun->Eval(pt))*downscale );
        if (pt < 10.) { corrPtResol = 1.0; }        
        Double_t cval = (val * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        Double_t cerr = (err * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        dNdPt->SetBinContent(i,cval);
        dNdPt->SetBinError(i,cerr);
    }
    // for dndeta again simpler
    dNdEta->Scale(1,"width");
    dNdEta->Scale(1./ReconstructedEvents);
    
    // normalization and finalization without resolution correction
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_nores->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_nores->GetBinCenter(i);
        Double_t width = dNdPt_nores->GetBinWidth(i);
        Double_t val = dNdPt_nores->GetBinContent(i);
        Double_t err = dNdPt_nores->GetBinError(i);
        Double_t cval = (val)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        Double_t cerr = (err)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * pt);
        dNdPt_nores->SetBinContent(i,cval);
        dNdPt_nores->SetBinError(i,cerr);
    }    
    
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
    dNdPt->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T}) (GeV/c)^{-2}");
    dNdEta->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    dNdPt_raw->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T}) (GeV/c)^{-2}");
    dNdEta_raw->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    Hists->Add(dNdEta_raw);
    Hists->Add(dNdPt_raw);
    Hists->Add(dNdEta);
    Hists->Add(dNdPt);
    Hists->Add(dNdPt_nores);
    
    
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
        if (xlabel.Contains("Pt")) { c->SetLogx();  c->SetLogy();  h->GetXaxis()->SetRangeUser(0.1 , 100.); }
        if (xlabel.Contains("p_{T}")) { c->SetLogx();  c->SetLogy();  h->GetXaxis()->SetRangeUser(0.1 , 100.); }        
    }
    if (h->GetDimension() == 2) {  
        TString ylabel(h->GetYaxis()->GetTitle());
        if (ylabel.Contains("Pt")) { c->SetLogy(); h->GetYaxis()->SetRangeUser(0.1 , 100.); }
        if (ylabel.Contains("p_{T}")) { c->SetLogy(); h->GetYaxis()->SetRangeUser(0.1 , 100.); }
        h->Draw("COLZ");
    }        
    if (h->GetDimension() == 1) {
        h->Draw();
    }
    return c;

}
