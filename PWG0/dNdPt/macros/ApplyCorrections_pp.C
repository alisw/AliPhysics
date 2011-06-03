Double_t ApplyCorrections_pp(const char* datafile, const char* datatask, const char* corrfile, const char* idstring ,const char* outfile, const char* gifdir = 0)
//Double_t ApplyCorrections_pp()
{

// tmp setting
// const char* mcfile = "/d/alice09/mknichel/train/V007.MC_pp/2011-05-05_2347.7147/mergedPeriods/MC_pp/7TeV/LHC10f6a/mknichel_dNdPtpp_TPCITS.root";
// const char* mctask = "mknichel_dNdPtpp_TPCITS";
// const char* idstring = "TPCITS";
// const char* gifdir = "/u/mknichel/alice/dNdPt_pp/temp";
// const char* datafile = "/d/alice09/mknichel/train/V007.MC_pp/2011-05-05_2347.7147/mergedPeriods/MC_pp/7TeV/LHC10f6a/mknichel_dNdPtpp_TPCITS.root";
// const char* outfile = "/u/mknichel/alice/dNdPt_pp/temp/FinalSpectra.root";
// const char* corrfile = "/u/mknichel/alice/dNdPt_pp/temp/corrMatr_LHC10f6a_TPCITS.root";
// const char* gifdir = "/u/mknichel/alice/dNdPt_pp/temp";

TString fname (datafile);

TString id (idstring);
TString taskname = "mknichel_dNdPtpp_" + id;
TString objname = "dNdPtAnalysis_" + id;

// tmp setting
 //     const char* datatask = "jotwinow_dNdPtAnalysis_TPCITS";
//      const char* corrfile = "/u/mknichel/alice/dNdPt_pp/2011-04-04_1251/corrMatr_LHC10f6a_TPCITS.root";
//      const char* idstring = "";
//      const char* outfile = "/u/mknichel/alice/dNdPt_pp/2011-04-04_1251/finalSpectra_LHC11a_TPCITS.root";
//      const char* gifdir = "/u/mknichel/alice/dNdPt_pp/2011-04-04_1251/plots";
    
    // settings vor zVertex cut (event and track level)
    Double_t zVert = 10.0;
    
    // setting on eta cut (track level)
    Double_t eta = 0.8;
    /*
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
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWG0base.so");    
    gSystem->Load("libPWG0dep"); 
    gSystem->Load("libPWG0selectors.so");
    */

    // make plots nicer
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    
    // array for all histograms to be saves
    TObjArray* Hists = new TObjArray();
    
    // open file with correction matrices
    TFile *fcorr = TFile::Open(corrfile,"READ");
    if (!fcorr) return -2;
    
    TH1D* dNdPt_MC      =  (TH1D*) fcorr->Get("dNdPt_MC");
    TH1D* dNdEta_MC      = (TH1D*) fcorr->Get("dNdEta_MC");
    Hists->Add(dNdPt_MC);
    Hists->Add(dNdEta_MC);
    
    // load data
    TFile* fdata = TFile::Open(datafile,"READ");
    if (!fdata) return -1;
    TList* ldata = dynamic_cast<TList*>(fdata->Get(taskname.Data()));
    if (!ldata) return -1;
    AlidNdPtAnalysis *obj = dynamic_cast<AlidNdPtAnalysis*>(ldata->FindObject(objname.Data()));
    if (!obj) return -1;
    
    //Event statistics
    obj->GetRecEventHist2()->GetAxis(0)->SetRange();
    THnSparse *fRecEventHist2 = obj->GetRecEventHist2(); //reconstructed events	
    TH1D* h1RecEventHist2_zv = (TH1D*) fRecEventHist2->Projection(0)->Clone("h1RecEventHist2_zv");  // zvertex distribution
    TH2D* h2RecEvent2All = (TH2D*) fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2All");
    fRecEventHist2->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    TH2D* h2RecEvent2 = (TH2D*) fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2");
    TH2D* h2RecEvent2Corrected = (TH2D*) fRecEventHist2->Projection(0,1)->Clone("h2RecEvent2Corrected"); //will be corrected
    
    THnSparse* fEventCount = obj->GetEventCount(); 
    Hists->Add(fEventCount);
    Double_t TriggeredEventsNoVertex = fEventCount->Projection(0)->GetBinContent(2) -  fEventCount->Projection(1)->GetBinContent(2); // all triggered events without rec. vertex
    Double_t AllTriggeredEvents = fEventCount->Projection(0)->GetBinContent(2); // all
    
    Double_t ReconstructedEvents = h2RecEvent2->Integral();
    Double_t ReconstructedEventsAll = h2RecEvent2All->Integral();
    Double_t TriggerEffInel = 0.864; //numeber from michele floris for 7TeV (I use it also for 2.76)
    TriggerEffInel = 0.9379; //mc value from mc
    if ( fname.Contains("900GeV") ) { TriggerEffInel = 0.9156; } 
    if ( fname.Contains("2.76TeV") ) { TriggerEffInel = 0.883; } 
    if ( fname.Contains("7TeV") ) { TriggerEffInel = 0.8524; } 
    
    cout << "Using Trigger_to_Inel efficiecy: " << TriggerEffInel <<endl;
    
    
    Double_t ReconstructedEventsFraction = (ReconstructedEventsAll-ReconstructedEvents) / ReconstructedEventsAll;
    Double_t SelectedEventsFraction = ReconstructedEvents / ReconstructedEventsAll;
    Double_t InelasticEventsSimple = AllTriggeredEvents/TriggerEffInel*SelectedEventsFraction;
    
    Hists->Add(h2RecEvent2);
    Hists->Add(h2RecEvent2All);
    Hists->Add(h2RecEvent2Corrected);
  
    printf("=== Number of events from DATA                      %lf ===\n",ReconstructedEvents);
    printf("=== Number of events from DATA (before zVertex cut) %lf ===\n",ReconstructedEventsAll);
    printf("=== Number of events from DATA (all triggered)      %lf ===\n",AllTriggeredEvents);
  
    TH1D* h1ReconstructedEvents = new TH1D("h1ReconstructedEvents","h1ReconstructedEvents",1,0,1);
    TH1D* h1ReconstructedEventsAll = new TH1D("h1ReconstructedEventsAll","h1ReconstructedEventsAll",1,0,1);
    
    h1ReconstructedEvents->Fill(0,ReconstructedEvents);
    h1ReconstructedEvents->SetEntries(ReconstructedEvents);
    h1ReconstructedEventsAll->Fill(0,ReconstructedEventsAll);
    h1ReconstructedEventsAll->SetEntries(ReconstructedEventsAll);
        
    Hists->Add(h1ReconstructedEvents);
    Hists->Add(h1ReconstructedEventsAll);
    
     // retrieve corrections (event level)
     TH2D* h2EventTriggerEffAll   = (TH2D*)fcorr->Get("h2EventTriggerEffAll");  
     TH2D* h2EventTriggerCorrAll  = (TH2D*)fcorr->Get("h2EventTriggerCorrAll");  
     TH2D* h2EventTriggerEff      = (TH2D*)fcorr->Get("h2EventTriggerEff");  
     TH2D* h2EventTriggerCorr     = (TH2D*)fcorr->Get("h2EventTriggerCorr");  
     TH2D* h2EventRecEffAll       = (TH2D*)fcorr->Get("h2EventRecEffAll");  
     TH2D* h2EventRecCorrAll      = (TH2D*)fcorr->Get("h2EventRecCorrAll");  
     TH2D* h2EventRecEff          = (TH2D*)fcorr->Get("h2EventRecEff");  
     TH2D* h2EventRecCorr         = (TH2D*)fcorr->Get("h2EventRecCorr");  
     TH2D* h2EventEffAll          = (TH2D*)fcorr->Get("h2EventEffAll");  
     TH2D* h2EventCorrAll         = (TH2D*)fcorr->Get("h2EventCorrAll");  
     TH2D* h2EventEff             = (TH2D*)fcorr->Get("h2EventEff");  
     TH2D* h2EventCorr            = (TH2D*)fcorr->Get("h2EventCorr"); 
     TH1D* h1TriggerEff_bin0_zv   = (TH1D*)fcorr->Get("h1TriggerEff_bin0_zv"); 
     TH1D* h1TriggerCorr_bin0_zv  = (TH1D*)fcorr->Get("h1TriggerCorr_bin0_zv"); 
     TH1D* h1Ratio_zv             = (TH1D*)fcorr->Get("h1Ratio_zv"); 
     TH1D* corr_shape_trig0_zv    = (TH1D*)fcorr->Get("corr_shape_trig0_zv"); 
     TH1D* corr_shape_notrig0_zv  = (TH1D*)fcorr->Get("corr_shape_notrig0_zv"); 
     
     
    Double_t corrTrackMatch = 1.0;
   //tracking efficiency is 2% to low in MC LHC10e13 900GeV as compared to data
   if ( fname.Contains("900GeV") ) {   
       corrTrackMatch = 1.02;
       cout << "900 GeV: correct tracking matching efficiency " <<endl;
   }

     
    // correct bin0
    //h1RecEventHist2_zv->Scale(1./h1RecEventHist2_zv->Integral());
    /*
    TH1D* h1Bin0_zv =  (TH1D*) h1RecEventHist2_zv->Clone("h1Bin0_zv");
    h1Bin0_zv->Multiply(h1Ratio_zv);
    h1Bin0_zv->Scale(1./h1Bin0_zv->Integral());
    h1Bin0_zv->Scale(TriggeredEventsNoVertex);
    //h1Bin0_zv->Multiply(h1TriggerCorr_bin0_zv);
    h1Bin0_zv->GetXaxis()->SetRangeUser(-zVert, zVert-0.01); //zVertex
    Double_t bin0EventsCorrected =  h1Bin0_zv->Integral();
    */
    // this is what is used for normalization
    Double_t InelasticEventsAll =  AllTriggeredEvents/TriggerEffInel; 
    Double_t Bin0EventsAll = TriggeredEventsNoVertex;
    Double_t UntriggeredEventsAll = InelasticEventsAll - AllTriggeredEvents;
    cout << "cross check: INEL Events 1 " << InelasticEventsAll <<endl;
    cout << "cross check: INEL Events 2 " << (ReconstructedEventsAll+Bin0EventsAll+UntriggeredEventsAll) <<endl;
    
    // correct bin0
    TH1D* h1Bin0_zv =  (TH1D*) h1RecEventHist2_zv->Clone("h1Bin0_zv");
    h1Bin0_zv->Multiply(corr_shape_trig0_zv);     //correct for shape
    h1Bin0_zv->Scale(1./h1Bin0_zv->Integral());    //normalize    
    h1Bin0_zv->Scale(Bin0EventsAll);
    Double_t Bin0Events = h1Bin0_zv->Integral(5,8);
    
    //correct untriggered
    TH1D* h1NoTrig_zv =  (TH1D*) h1RecEventHist2_zv->Clone("h1NoTrig_zv");
    h1NoTrig_zv->Multiply(corr_shape_notrig0_zv);     //correct for shape
    h1NoTrig_zv->Scale(1./h1NoTrig_zv->Integral());    //normalize    
    h1NoTrig_zv->Scale(UntriggeredEventsAll);
    Double_t UntriggeredEvents = h1NoTrig_zv->Integral(5,8);
    
    
    Double_t InelasticEvents = ReconstructedEvents + Bin0Events + UntriggeredEvents;
    
    /*
    // my way (old)
    TH1D* h1Bin0_zv =  (TH1D*) h1RecEventHist2_zv->Clone("h1Bin0_zv");
    h1Bin0_zv->Multiply(h1Ratio_zv);
    h1Bin0_zv->Multiply(h1TriggerEff_bin0_zv);
    h1Bin0_zv->Scale(1./h1Bin0_zv->Integral());
    h1Bin0_zv->Scale(TriggeredEventsNoVertex);
    h1Bin0_zv->Multiply(h1TriggerCorr_bin0_zv);
    Double_t bin0EventsCorrected = h1Bin0_zv->Integral(5,8);
     */

    /*
    //jaceks way (?)
    h1RecEventHist2_zv->Scale(1./h1RecEventHist2_zv->Integral());
    h1RecEventHist2_zv->Scale(TriggeredEventsNoVertex); 
    Double_t bin0EventsCorrected = h1RecEventHist2_zv->Integral(5,8);    
*/
     // retrieve tracks
    THnSparse* hSRecTrack; // = obj->GetRecTrackHist()->Clone("hsRecTrack"); //thnsparse to be corrected
    THnSparse* fRecTrackHist2 = obj->GetRecTrackHist(); //after all cuts (2)
    fRecTrackHist2->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01); //zVertex
    fRecTrackHist2->GetAxis(2)->SetRangeUser(-eta, eta-0.01); // eta   
    Int_t dims[4] = {0,1,2,3};
    hSRecTrack = (THnSparse*) (fRecTrackHist2->Projection(4,dims,"e")->Clone("hSRecTrack"));        
    hSRecTrackAllMult = (THnSparse*) (fRecTrackHist2->Projection(3,dims,"e")->Clone("hSRecTrackAllMult"));
    
    TH3D* h3RecTrackHist2 = fRecTrackHist2->Projection(0,1,2)->Clone("h3RecTrackHist2");
   
    TH1D* h1RecTrackHist2_zv = fRecTrackHist2->Projection(0)->Clone("h1RecTrackHist2_zv");
    TH1D* h1RecTrackHist2_pt = fRecTrackHist2->Projection(1)->Clone("h1RecTrackHist2_pt");
    TH1D* h1RecTrackHist2_eta = fRecTrackHist2->Projection(2)->Clone("h1RecTrackHist2_eta");
   
    Hists->Add(h1RecTrackHist2_zv);
    Hists->Add(h1RecTrackHist2_pt);
    Hists->Add(h1RecTrackHist2_eta);
    Hists->Add(h3RecTrackHist2);
    
    // generate corrections for events
    TH1D* h1MCGeneratedEvents     = (TH1D*)fcorr->Get("h1MCGeneratedEvents");  
    TH1D* h1MCReconstructedEvents = (TH1D*)fcorr->Get("h1MCReconstructedEvents");
    TH1D* h1MCTriggeredEvents0mult = (TH1D*)fcorr->Get("h1MCTriggeredEvents0mult");
    TH1D* h1MCTriggeredEventsAll0mult = (TH1D*)fcorr->Get("h1MCTriggeredEventsAll0mult");
    TH1D* h1MCTriggeredEventsAll = (TH1D*)fcorr->Get("h1MCTriggeredEventsAll");
    TH1D* h1MCTriggeredEvents = (TH1D*)fcorr->Get("h1MCTriggeredEvents");
        
    
    Double_t MCGeneratedEvents = h1MCGeneratedEvents->GetEntries();
    Double_t MCReconstructedEvents = h1MCReconstructedEvents->GetEntries();
    Double_t CorrEvent = MCGeneratedEvents / MCReconstructedEvents;
    
    Double_t MCTriggeredEvents = h1MCTriggeredEvents->GetEntries();
    Double_t MCTriggeredEventsAll = h1MCTriggeredEventsAll->GetEntries();
    Double_t MCTriggeredEvents0mult = h1MCTriggeredEvents0mult->GetEntries();
    Double_t MCTriggeredEventsAll0mult = h1MCTriggeredEventsAll0mult->GetEntries();        
    Double_t CorrVtxEvent0mult    = MCTriggeredEvents / (MCTriggeredEvents-MCTriggeredEvents0mult); //this is used
    Double_t CorrVtxEventAll0mult = MCTriggeredEventsAll / (MCTriggeredEventsAll-MCTriggeredEventsAll0mult); 
    
    // correct for trigger/vertex inefficiencies
    for (int xbin=1; xbin <= h2RecEvent2Corrected->GetNbinsX(); xbin++) {
        for (int ybin=1; ybin <= h2RecEvent2Corrected->GetNbinsY(); ybin++) {
            Double_t x = h2RecEvent2Corrected->GetXaxis()->GetBinCenter(xbin);
            Double_t y = h2RecEvent2Corrected->GetYaxis()->GetBinCenter(ybin);
            Int_t bin  = h2EventCorr->FindBin(x,y);            
            Double_t corr    = h2EventCorr->GetBinContent(bin);
            Double_t correrr = h2EventCorr->GetBinError(bin);
            if (corr < 0.01) { corr=1.; correrr=0.;} //bin empty in correction matrix            
            Double_t val     = h2RecEvent2Corrected->GetBinContent(xbin,ybin);
            Double_t err     = h2RecEvent2Corrected->GetBinError(xbin,ybin);
            h2RecEvent2Corrected->SetBinContent(xbin,ybin,corr*val);
            h2RecEvent2Corrected->SetBinError(xbin,ybin,TMath::Sqrt(val*val*correrr*correrr + err*err*corr*corr));
            
        }
    }
    //Double_t ReconstructedEventsCorrected = h2RecEvent2Corrected->Integral() * CorrVtxEvent0mult;
    //Double_t ReconstructedEventsCorrected = h2RecEvent2Corrected->Integral() + bin0EventsCorrected;
    printf("=== Number of events from DATA (final correction)   %lf ===\n",InelasticEvents); 
    printf("=== Number of events from DATA (corr, no bin0)      %lf ===\n",h2RecEvent2Corrected->Integral());
    printf("=== Number of events from DATA (simple corr. 7TeV)  %lf ===\n",InelasticEventsSimple);          
    printf("=== Number of events from DATA (simple correction)  %lf ===\n",ReconstructedEvents * CorrEvent);
    printf("=== Number of bin0events from DATA                  %lf ===\n",TriggeredEventsNoVertex); 
    //printf("=== Number of bin0events from DATA (corrected)      %lf ===\n",bin0EventsCorrected); 
    
    // retrieve 3D correction matrices for tracking efficiency and secondaries (note different binning!)
    TH3D* h3TrackEff   = (TH3D*)fcorr->Get("h3TrackEff");  
    TH3D* h3TrackCorr  = (TH3D*)fcorr->Get("h3TrackCorr");  
    TH3D* h3SecCont    = (TH3D*)fcorr->Get("h3SecCont");  
    TH3D* h3SecCorr    = (TH3D*)fcorr->Get("h3SecCorr");      
    
    // retrieve 3D thnsparse correction matrices for tracking efficiency and secondaries (note different binning!) --- this is used!!!!
    THnSparse* hSTrackEff = (THnSparse*) fcorr->Get("hSTrackEff");
    THnSparse* hSTrackCorr = (THnSparse*) fcorr->Get("hSTrackCorr");
    THnSparse* hSSecCont = (THnSparse*) fcorr->Get("hSSecCont");
    THnSparse* hSSecCorr = (THnSparse*) fcorr->Get("hSSecCorr");
        

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
        if (pt >= 50) { pt = 49.5; } // above 50 GeV corr matr have low statistics
        Double_t effCorr    = h1TrackCorr_pt->GetBinContent(h1TrackCorr_pt->FindBin(pt));
        Double_t effCorrErr = h1TrackCorr_pt->GetBinError(h1TrackCorr_pt->FindBin(pt));
        
        Double_t cval = corrTrackMatch * val * effCorr;
        Double_t cerr = corrTrackMatch * TMath::Sqrt(val*val*effCorrErr*effCorrErr + err*err*effCorr*effCorr);
        h1Corrected_pt->SetBinContent(i,cval);
        h1Corrected_pt->SetBinError(i,cerr);
    }
       
       
    
   // efficiency and contamination correction for thnsparse
   for (Long64_t j = 0; j < hSRecTrack->GetNbins(); j++) {
       Int_t tc[4];
       Double_t tval    = hSRecTrack->GetBinContent(j,tc);
       Double_t terr    = hSRecTrack->GetBinError(j);
       Double_t tzv     = hSRecTrack->GetAxis(0)->GetBinCenter(tc[0]);
       Double_t tpt     = hSRecTrack->GetAxis(1)->GetBinCenter(tc[1]);
       Double_t teta    = hSRecTrack->GetAxis(2)->GetBinCenter(tc[2]);
       Double_t tmultMB = hSRecTrack->GetAxis(3)->GetBinCenter(tc[3]);
       if (tzv >= 10.0) { tzv = 9.9; }
       if (teta >= 0.8) {teta = 0.79;}
       if (tzv < -10.0) { tzv = -10.; }
       if (teta < -0.8) { teta = -0.8; }              
       if (tpt >= 50.) { tpt = 49.5; } // above 50 GeV corr matr have low statistics 
       if (tpt < 0.15) { tpt = 0.175; } // also below 0.15       
       Double_t xvals[3];
       xvals[0] = tzv;
       xvals[1] = tpt;
       xvals[2] = teta;
       Double_t effCorr    = hSTrackCorr->GetBinContent(hSTrackCorr->GetBin(xvals,kFALSE));
       Double_t effCorrErr = hSTrackCorr->GetBinError(hSTrackCorr->GetBin(xvals,kFALSE));
       Double_t secCorr    = hSSecCont->GetBinContent(hSSecCont->GetBin(xvals,kFALSE));
       Double_t secCorrErr = hSSecCont->GetBinError(hSSecCont->GetBin(xvals,kFALSE));
//        Double_t effCorr    = h3TrackCorr->GetBinContent(h3TrackCorr->FindBin(tzv,tpt,teta));
//        Double_t effCorrErr = h3TrackCorr->GetBinError(h3TrackCorr->FindBin(tzv,tpt,teta));
//        Double_t secCorr    = 1. - h3SecCont->GetBinContent(h3SecCont->FindBin(tzv,tpt,teta));
//        Double_t secCorrErr = h3SecCont->GetBinError(h3SecCont->FindBin(tzv,tpt,teta));
       if (effCorr < 1.) { cout << "bin empty, use efficiency 1!" << endl; effCorr=1.; }
       if (secCorr < 0.001) { cout << "bin empty, use contamination 0!" << endl; }
       if (effCorrErr < 1e-9) { cout << "eff error empty!" << endl; }
       if (secCorrErr < 1e-9) { cout << "cont error empty!" << endl;}
       Double_t corr    = corrTrackMatch * effCorr*(1.-secCorr);
       //Double_t corrErr = corrTrackMatch * TMath::Sqrt(effCorr*secCorrErr + secCorr*effCorrErr); 
       Double_t ctval = tval*corr;
       Double_t cterr = terr*corr; // errors are not correlated        
       //Double_t cterr = TMath::Sqrt(tval*tval*corrErr*corrErr + terr*terr*corr*corr); // errors are not correlated       
       hSRecTrack->SetBinContent(j,ctval);
       hSRecTrack->SetBinError(j,cterr);
    }
    
    TH1D* effCorrPt = hSRecTrack->Projection(1)->Clone("effCorrPt");
    TH1D* effCorrPtErr2 = hSRecTrack->Projection(1)->Clone("effCorrPtErr2");
    TH1D* effCorrPtNorm = hSRecTrack->Projection(1)->Clone("effCorrPtNorm");
    effCorrPt->Reset();
    effCorrPtErr2->Reset();
//     effCorrPtNorm->Reset();
    /*
for (Long64_t j = 0; j < hSRecTrackAllMult->GetNbins(); j++) {
       Int_t tc[3];
       Double_t tval    = hSRecTrackAllMult->GetBinContent(j,tc);
       Double_t tzv     = hSRecTrackAllMult->GetAxis(0)->GetBinCenter(tc[0]);
       Double_t tpt     = hSRecTrackAllMult->GetAxis(1)->GetBinCenter(tc[1]);
       Double_t teta    = hSRecTrackAllMult->GetAxis(2)->GetBinCenter(tc[2]);
       if (tpt >= 50.) { tpt = 49.5; } // above 50 GeV corr matr have low statistics 
       Double_t effCorr    = h3TrackCorr->GetBinContent(h3TrackCorr->FindBin(tzv,tpt,teta));
       Double_t effCorrErr = h3TrackCorr->GetBinError(h3TrackCorr->FindBin(tzv,tpt,teta));
       Double_t secCorr    = 1. - h3SecCont->GetBinContent(h3SecCont->FindBin(tzv,tpt,teta));
       Double_t secCorrErr = h3SecCont->GetBinError(h3SecCont->FindBin(tzv,tpt,teta));
       if (effCorr < 0.1) { cout << "bin empty, use efficiency 1!" << endl; effCorr=1.; effCorrErr=0.; }
       if (secCorr < 0.1) { cout << "bin empty, use contamination 0!" << endl; secCorr=1.; secCorrErr=0.; }
       Double_t corr    = effCorr*secCorr;
       Double_t corrErr = effCorr*secCorrErr + secCorr*effCorrErr; // errors are correlated
       effCorrPt->Fill(tpt,tval*corr);
       effCorrPtErr2->Fill(tpt,tval*tval*corrErr*corrErr);
       effCorrPtNorm->Fill(tpt,tval);
    }
       
    effCorrPt->Divide(effCorrPtNorm);
    for (Int_t i = 4; i <= effCorrPtErr2->GetNbinsX(); i++) { effCorrPtErr2->SetBinContent(i,TMath::Sqrt(effCorrPtErr2->GetBinContent(i)));  }
    effCorrPtErr2->Divide(effCorrPtNorm);
    cout << effCorrPt->GetBinContent(4) << endl;
    cout << effCorrPtErr2->GetBinContent(4) << endl;
    cout << effCorrPtNorm->GetBinContent(4) << endl;
    
    
    
    TH1D* RecTrackPtCorrected = hSRecTrack->Projection(1)->Clone("RecTrackPtCorrected");
    RecTrackPtCorrected->Reset();    
    

    for (Int_t i = 4; i <= RecTrackPtCorrected->GetNbinsX(); i++) {    
       if (0 == effCorrPtErr2->GetBinContent(i)) continue;
       Double_t val     = effCorrPt->GetBinContent(i); 
       cout << (val*val/effCorrPtNorm->GetBinContent(i)) << " " << effCorrPtErr2->GetBinContent(i) << endl;
       Double_t err2  = (val*val/effCorrPtNorm->GetBinContent(i)) + effCorrPtErr2->GetBinContent(i);
       RecTrackPtCorrected->SetBinContent(i,val);
       RecTrackPtCorrected->SetBinError(i,TMath::Sqrt(err2));
    }
    */
    
for (Long64_t j = 0; j < hSRecTrackAllMult->GetNbins(); j++) {
       Int_t tc[3];
       Double_t tval    = hSRecTrackAllMult->GetBinContent(j,tc);
       Double_t terr    = TMath::Sqrt(tval); //hSRecTrackAllMult->GetBinError(j);
       Double_t tzv     = hSRecTrackAllMult->GetAxis(0)->GetBinCenter(tc[0]);
       Double_t tpt     = hSRecTrackAllMult->GetAxis(1)->GetBinCenter(tc[1]);
       Double_t teta    = hSRecTrackAllMult->GetAxis(2)->GetBinCenter(tc[2]);
       Double_t tmultMB = hSRecTrack->GetAxis(3)->GetBinCenter(tc[3]);
       if (tzv >= 10.0) { tzv = 9.9; }
       if (teta >= 0.8) {teta = 0.79;}
       if (tzv < -10.0) { tzv = -10.; }
       if (teta < -0.8) { teta = -0.8; }              
       if (tpt >= 50.) { tpt = 49.5; } // above 50 GeV corr matr have low statistics 
       if (tpt < 0.15) { tpt = 0.175; } // also below 0.15       
      Double_t xvals[3];
       xvals[0] = tzv;
       xvals[1] = tpt;
       xvals[2] = teta;
       Double_t effCorr    = hSTrackCorr->GetBinContent(hSTrackCorr->GetBin(xvals,kFALSE));
       Double_t effCorrErr = hSTrackCorr->GetBinError(hSTrackCorr->GetBin(xvals,kFALSE));
       Double_t secCont    = hSSecCont->GetBinContent(hSSecCont->GetBin(xvals,kFALSE));
       Double_t secContErr = hSSecCont->GetBinError(hSSecCont->GetBin(xvals,kFALSE));       
//        Double_t effCorr    = h3TrackCorr->GetBinContent(h3TrackCorr->FindBin(tzv,tpt,teta));
//        Double_t effCorrErr = h3TrackCorr->GetBinError(h3TrackCorr->FindBin(tzv,tpt,teta));
//        Double_t secCorr    = 1. - h3SecCont->GetBinContent(h3SecCont->FindBin(tzv,tpt,teta));
//        Double_t secCorrErr = h3SecCont->GetBinError(h3SecCont->FindBin(tzv,tpt,teta));       
       if (effCorr < 1.) { cout << "bin empty, use efficiency 1!" << endl; effCorr=1.; }
       if (secCont < 0.001) { cout << "bin empty, use contamination 0!" << endl; }
       if (effCorrErr < 1e-9) { cout << "eff error empty!" << endl; }
       if (secContErr < 1e-9) { cout << "cont error empty!" << endl;}
       Double_t secCorr = (1.-secCont);
       Double_t corr    = corrTrackMatch * (effCorr*secCorr);
       Double_t corrErr = corrTrackMatch *  TMath::Sqrt(effCorr*secCorrErr + secCorr*effCorrErr); 
       Double_t ctval = tval*corr;
       Double_t cterr = terr*corr; // errors are not correlated        
       //Double_t cterr = TMath::Sqrt(tval*tval*corrErr*corrErr + terr*terr*corr*corr); // errors are not correlated       
       hSRecTrackAllMult->SetBinContent(j,ctval);
       hSRecTrackAllMult->SetBinError(j,cterr);
    }
    
    // create final spectrum in multiplicity bins
    TH2D* dNdPtMult = hSRecTrack->Projection(3,1)->Clone("dNdPtMult");
    
    // create final spectra (as clone of corrected data)
    TH1D* dNdPt = hSRecTrackAllMult->Projection(1)->Clone("dNdPt");
/*
    // errors on pt spectrum from corrections
    for (int i=1; i <= dNdPt->GetNbinsX() ; i++) {
        Double_t pt = dNdPt->GetBinCenter(i);
        Double_t val = dNdPt->GetBinContent(i);
        Double_t err = dNdPt->GetBinError(i);
        if (pt >= 50) { pt = 49.5; }  // above 50 GeV corr matr have low statistics
        Double_t secCorr    = h1SecCorr_pt->GetBinContent(h1SecCorr_pt->FindBin(pt));
        Double_t secCorrErr = h1SecCorr_pt->GetBinError(h1SecCorr_pt->FindBin(pt));
        Double_t effCorr    = h1TrackCorr_pt->GetBinContent(h1TrackCorr_pt->FindBin(pt));
        Double_t effCorrErr = h1TrackCorr_pt->GetBinError(h1TrackCorr_pt->FindBin(pt));
        Double_t corr    = effCorr*secCorr;
        Double_t corrErr = effCorr*secCorrErr + secCorr*effCorrErr; // errors are correlated         
        Double_t cval = val*secCorr;
        Double_t cerr = TMath::Sqrt(val*val*corrErr*corrErr + err*err*corr*corr);
        dNdPt->SetBinError(i,cerr);
    }
    
    // errors on pt spectrum from corrections
    for (int i=1; i <= dNdPtMult->GetNbinsX() ; i++) {
    for (int k=1; k <= dNdPtMult->GetNbinsY() ; k++) {
        Double_t pt = dNdPt->GetBinCenter(i);
        Double_t val = dNdPtMult->GetBinContent(i,k);
        Double_t err = dNdPtMult->GetBinError(i,k);
        if (pt >= 50) { pt = 49.5; }  // above 50 GeV corr matr have low statistics
        Double_t secCorr    = h1SecCorr_pt->GetBinContent(h1SecCorr_pt->FindBin(pt));
        Double_t secCorrErr = h1SecCorr_pt->GetBinError(h1SecCorr_pt->FindBin(pt));
        Double_t effCorr    = h1TrackCorr_pt->GetBinContent(h1TrackCorr_pt->FindBin(pt));
        Double_t effCorrErr = h1TrackCorr_pt->GetBinError(h1TrackCorr_pt->FindBin(pt));
        Double_t corr    = effCorr*secCorr;
        Double_t corrErr = effCorr*secCorrErr + secCorr*effCorrErr; // errors are correlated         
        Double_t cval = val*secCorr;
        Double_t cerr = TMath::Sqrt(val*val*corrErr*corrErr + err*err*corr*corr);
        dNdPtMult->SetBinError(i,k,cerr);
    }
    }    
*/
    
    // for eta the correction is simpler because of same binning
    h1Corrected_eta->Multiply(h1SecCorr_eta);
    h1Corrected_eta->Multiply(h1TrackCorr_eta);
    
    Hists->Add(h1Corrected_pt);
    Hists->Add(h1Corrected_eta);


    TH1D* dNdPt_simple = (TH1D*) h1Corrected_pt->Clone("dNdPt_simple");
    TH1D* dNdPt_simple_nores  = (TH1D*) h1Corrected_pt->Clone("dNdPt_simple_nores");
    
    
    TH1D* dNdEta  = hSRecTrackAllMult->Projection(2)->Clone("dNdEta");
    TH1D* dNdEta_simple  = h1Corrected_eta->Clone("dNdEta_simple");
       
    // also uncorrected spectra (as clone of raw data)
    TH1D* dNdPt_raw   = (TH1D*) h1RecTrackHist2_pt->Clone("dNdPt_raw");
    TH1D* dNdEta_raw  = (TH1D*) h1RecTrackHist2_eta->Clone("dNdEta_raw");
    
    // also uncorrected spectra (as clone of raw data)
    TH1D* dNdPt_nores          = (TH1D*) dNdPt->Clone("dNdPt_nores");    
    //TH1D* dNdEta_nores  = h1Corrected_eta->Clone("dNdEta_nores");
    
       

//TF1 *fperi  = new TF1("fperi","1.00343-0.000608425*x-6.7038e-05*x*x",5.,50.);
//TF1 *fcent  = new TF1("cent","1.01074e+00-1.98127e-03*x-1.19903e-04*x*x",5.,50.);
TFile* fptcorr = TFile::Open("ptcorr_140511.root");

TF1 *fun = 0;
TString id = TString(idstring);
if ( fname.Contains("900GeV") ) {
  TF1 *fun  = (TF1*) fptcorr->Get("ptcorr_900");
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "900GeV pt-resolution correction used!" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++" << endl;  
} elseif ( fname.Contains("2.76TeV") ) {
  TF1 *fun  = (TF1*) fptcorr->Get("ptcorr_2760");
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "2.76TeV pt-resolution correction used!" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
} elseif ( fname.Contains("7TeV") ) {
  TF1 *fun  = (TF1*) fptcorr->Get("ptcorr_7000");
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "7TeV pt-resolution correction used!" << endl;
  cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;  
} else {
  fun = TF1("none","1.",10.,50.);
  cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "ERROR! NO pt-resolution correction used!" << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
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
        if (pt < 10.) { corrPtResol = 1.0; }
        if (corrPtResol > 1.0 ) { corrPtResol = 1.0; }        
        Double_t cval = (val * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        Double_t cerr = (err * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
//         Double_t cval = (val * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * CorrEvent* pt);
//         Double_t cerr = (err * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * ReconstructedEvents * CorrEvent* pt);        
        dNdPt->SetBinContent(i,cval);
        dNdPt->SetBinError(i,cerr);
    }
    
    // normalization and finalization without resolution correction
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_nores->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_nores->GetBinCenter(i);
        Double_t width = dNdPt_nores->GetBinWidth(i);
        Double_t val = dNdPt_nores->GetBinContent(i);
        Double_t err = dNdPt_nores->GetBinError(i);
        Double_t cval = (val)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        Double_t cerr = (err)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        dNdPt_nores->SetBinContent(i,cval);
        dNdPt_nores->SetBinError(i,cerr);
    }
    
    // normalization and finalization without resolution correction
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_simple_nores->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_simple_nores->GetBinCenter(i);
        Double_t width = dNdPt_simple_nores->GetBinWidth(i);
        Double_t val = dNdPt_simple_nores->GetBinContent(i);
        Double_t err = dNdPt_simple_nores->GetBinError(i);
        Double_t cval = (val)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        Double_t cerr = (err)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        dNdPt_simple_nores->SetBinContent(i,cval);
        dNdPt_simple_nores->SetBinError(i,cerr);
    }    
        
    // for dndeta again simpler
    dNdEta->Scale(1,"width");
    dNdEta->Scale(1./InelasticEvents);
    dNdEta_simple->Scale(1,"width");
    dNdEta_simple->Scale(1./InelasticEvents);
    
    
    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_simple->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_simple->GetBinCenter(i);
        Double_t width = dNdPt_simple->GetBinWidth(i);
        Double_t val = dNdPt_simple->GetBinContent(i);
        Double_t err = dNdPt_simple->GetBinError(i);
        Double_t corrPtResol = 1.0;
        corrPtResol = fun->Eval(pt);
        if (pt < 10.) { corrPtResol = 1.0; }
        if (corrPtResol > 1.0 ) { corrPtResol = 1.0; }
        Double_t cval = (val * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        Double_t cerr = (err * corrPtResol)/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        dNdPt_simple->SetBinContent(i,cval);
        dNdPt_simple->SetBinError(i,cerr);
    }
    /*
    TH1D* dNdPt2 = (TH1D*) dNdPt->Clone("dNdPt2");
    TH1D* dNdPt2_simple = (TH1D*) dNdPt_simple->Clone("dNdPt2_simple");
    dNdPt2->Scale(InelasticEvents / (ReconstructedEvents * CorrEvent));
    dNdPt2_simple->Scale(InelasticEvents / (ReconstructedEvents * CorrEvent));
    Hists->Add(dNdPt2);
    Hists->Add(dNdPt2_simple);    
    */
    
    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int i=1; i <= dNdPt_raw->GetNbinsX() ;i++) {
        Double_t pt = dNdPt_raw->GetBinCenter(i);
        Double_t width = dNdPt_raw->GetBinWidth(i);
        Double_t val = dNdPt_raw->GetBinContent(i);
        Double_t err = dNdPt_raw->GetBinError(i);
        Double_t cval = val/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        Double_t cerr = err/(width * 2.0 * TMath::Pi() * 1.6 * InelasticEvents * pt);
        dNdPt_raw->SetBinContent(i,cval);
        dNdPt_raw->SetBinError(i,cerr);
    }
    // for dndeta again simpler
    dNdEta_raw->Scale(1,"width");
    dNdEta_raw->Scale(1./InelasticEvents);
     
    dNdEta->SetMarkerStyle(21);
    dNdPt->SetMarkerStyle(21);
    dNdPt->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdEta->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    dNdPt_raw->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdPt_nores->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdEta_raw->SetTitle("; #eta ; 1/N_{evt} (d^{2}N_{ch})/(d#eta)");
    
    Hists->Add(dNdEta_raw);
    Hists->Add(dNdPt_raw);
    Hists->Add(dNdEta);
    Hists->Add(dNdEta_simple);
    Hists->Add(dNdPt);
    Hists->Add(dNdPt_nores);  
    Hists->Add(dNdPt_simple);  
    Hists->Add(dNdPt_simple_nores);  
    Hists->Add(hSRecTrack);
    Hists->Add(hSRecTrackAllMult);
    
    
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
