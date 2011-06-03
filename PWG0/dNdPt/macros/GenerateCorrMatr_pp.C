Double_t GenerateCorrMatr_pp(const char* mcfile, const char* mctask, const char* idstring ,const char* outfile, const char* gifdir = 0)
//Double_t GenerateCorrMatr_pp()
{
//tmp setting
// const char* mcfile = "/d/alice09/mknichel/train/V007.MC_pp/2011-05-05_2347.7147/mergedPeriods/MC_pp/7TeV/LHC10f6a/mknichel_dNdPtpp_TPCITS.root";
// const char* mctask = "mknichel_dNdPtpp_TPCITS";
// const char* idstring = "TPCITS";
// const char* outfile = "/u/mknichel/alice/dNdPt_pp/temp/corrMatr_LHC10f6a_TPCITS.root";
// const char* gifdir = "/u/mknichel/alice/dNdPt_pp/temp";

TString fname (mcfile);

TString id (idstring);
TString taskname = "mknichel_dNdPtpp_" + id;
TString objname = "dNdPtAnalysis_" + id;

    // settings vor zVertex cut (event and track level)
    Double_t zVert = 10.0;
    
    // setting on eta cut (track level)
    Double_t eta = 0.8;
    
    // strangeness scaling factor (for secondaries from strange decays)
    Double_t sscale = 2.0;  //ignored
    
    
    // Set pt binning
    const Int_t ptNbins = 31;
    Double_t bins[32] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,2.4,2.8, 3.0,  4.0, 50.0, 100.0 };
    Double_t* binsPt = new Double_t[32];
    for (int ibin=0; ibin<32; ibin++) {binsPt[ibin] = bins[ibin];}
    
 

  Int_t ptNbinsTrackEventCorr = 36;
  Int_t etaNbins = 30;
  Int_t zvNbins = 12;

  Double_t binsMultDefault[28] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,
                                     9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
				     19.5,20.5, 21.5, 22.5, 23.5, 24.5, 29.5, 149.5};

  Double_t binsPtTrackEventCorrDefault[37] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,4.0,50.0};

  Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};

  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  
    // binning for corrections (note difference for efficiency and secondaries)
    if (fname.Contains("2.76TeV") || fname.Contains("900GeV")) {
    
    const Int_t ptNbinsSecCorr  = 28;
    const Int_t ptNbinsEffCorr  = 31;
    const Int_t etaNbinsCorr = 4;
    const Int_t zvNbinsCorr  = 4;
    const Int_t multNbinsCorr = 2;
    Double_t ptBinsSecCorr[29] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 3.0, 50.0, 100.0 };    
    Double_t ptBinsEffCorr[32] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.2, 1.4, 1.6, 2.0, 2.4, 2.8, 3.0, 3.4, 4.0, 50.0, 100.0 };
    Double_t etaBinsCorr[5] = {-0.8,-0.4,0.,0.4,0.8};
    
    Double_t zvBinsCorr[5] = {-10.,-5.,0.,5.,10.};
    Double_t multBinsCorr[3] = {-0.5,0.5,500.};
    
    
    } else {
    
    const Int_t ptNbinsSecCorr  = 38;
    const Int_t ptNbinsEffCorr  = 50;
    const Int_t etaNbinsCorr = 8;
    const Int_t zvNbinsCorr  = 4;
    const Int_t multNbinsCorr = 2;
    Double_t ptBinsSecCorr[39] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,3.4,4.0,50.0,100.0};
    Double_t ptBinsEffCorr[51] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 50.0, 100.0 };
    Double_t etaBinsCorr[9] = {-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8};
    
    Double_t zvBinsCorr[5] = {-10.,-5.,0.,5.,10.};
    Double_t multBinsCorr[3] = {-0.5,0.5,500.};
    
    }
    
    
    Int_t binsTrackEffCorr[3]={zvNbinsCorr,ptNbinsEffCorr,etaNbinsCorr};
    Int_t binsTrackSecCorr[3]={zvNbinsCorr,ptNbinsSecCorr,etaNbinsCorr};
    
    THnSparse* hSTrackEffCorr = new THnSparseF("hSTrackEffCorr","Zv:pT:eta",3,binsTrackEffCorr);
    hSTrackEffCorr->SetBinEdges(0,zvBinsCorr);
    hSTrackEffCorr->SetBinEdges(1,ptBinsEffCorr);
    hSTrackEffCorr->SetBinEdges(2,etaBinsCorr);
//     hSTrackEffCorr->SetBinEdges(3,multBinsCorr);
    hSTrackEffCorr->GetAxis(0)->SetTitle("Zv (cm)");
    hSTrackEffCorr->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    hSTrackEffCorr->GetAxis(2)->SetTitle("#eta");
//     hSTrackEffCorr->GetAxis(3)->SetTitle("mult MB");
    hSTrackEffCorr->Sumw2();    
    
    THnSparse* hSTrackSecCorr = new THnSparseF("hSTrackSecCorr","Zv:pT:eta",3,binsTrackSecCorr);    
    hSTrackSecCorr->SetBinEdges(0,zvBinsCorr);
    hSTrackSecCorr->SetBinEdges(1,ptBinsSecCorr);
    hSTrackSecCorr->SetBinEdges(2,etaBinsCorr);
//     hSTrackSecCorr->SetBinEdges(3,multBinsCorr);
    hSTrackSecCorr->GetAxis(0)->SetTitle("Zv (cm)");
    hSTrackSecCorr->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    hSTrackSecCorr->GetAxis(2)->SetTitle("#eta");
//     hSTrackSecCorr->GetAxis(3)->SetTitle("mult MB");
    hSTrackSecCorr->Sumw2();
    

    // make plots nicer
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    
    // array for all correction matrices
    TObjArray* CorrMatr = new TObjArray();
    
    // array for all control histograms
    TObjArray* ContHist = new TObjArray();
    

    // load mc information
    TFile* fmc = TFile::Open(mcfile,"READ");
    if (!fmc) return -1;
    TList* lmc = dynamic_cast<TList*>(fmc->Get(taskname.Data()));
    if (!lmc) return -1;
    AlidNdPtAnalysis *obj = dynamic_cast<AlidNdPtAnalysis*>(lmc->FindObject(objname.Data()));
    if (!obj) return -1;

    //Event statistics
    THnSparse *fRecEventMatrix = obj->GetRecEventMatrix(); //all reconstructed events	
    TH2D* h2RecEventAll = fRecEventMatrix->Projection(0,1)->Clone("h2RecEventAll");
    fRecEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    TH2D* h2RecEvent = fRecEventMatrix->Projection(0,1)->Clone("h2RecEvent");
    Double_t MCReconstructedEvents = h2RecEvent->Integral();
    Double_t MCReconstructedEventsAll = h2RecEventAll->Integral();
    ContHist->Add(h2RecEvent);
    ContHist->Add(h2RecEventAll);
        
    THnSparse *fTriggerEventMatrix = obj->GetTriggerEventMatrix(); //all triggered events
    Double_t MCTriggeredEventsAll0mult =  obj->GetTriggerEventMatrix()->Projection(1)->Integral(1,1); // 0 contributers to vertex
    TH2D* h2TriggerEventAll = fTriggerEventMatrix->Projection(0,1)->Clone("h2TriggerEventAll");
    fTriggerEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    Double_t MCTriggeredEvents0mult =  obj->GetTriggerEventMatrix()->Projection(1)->Integral(1,1); // 0 contributers to vertex
    TH2D* h2TriggerEvent = fTriggerEventMatrix->Projection(0,1)->Clone("h2TriggerEvent");
    Double_t MCTriggeredEvents = h2TriggerEvent->Integral();
    Double_t MCTriggeredEventsAll = h2TriggerEventAll->Integral();
    ContHist->Add(h2TriggerEvent);
    ContHist->Add(h2TriggerEventAll);
        
    THnSparse *fGenEventMatrix = obj->GetGenEventMatrix(); //all generated events
    TH2D* h2GenEventAll = fGenEventMatrix->Projection(0,1)->Clone("h2GenEventAll");
    fGenEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    TH2D* h2GenEvent = fGenEventMatrix->Projection(0,1)->Clone("h2GenEvent");
    Double_t MCGeneratedEvents = h2GenEvent->Integral();
    Double_t MCGeneratedEventsAll = h2GenEventAll->Integral(); 
    ContHist->Add(h2RecEvent);
    ContHist->Add(h2RecEventAll);
    
    THnSparse *fRecEventHist = obj->GetRecEventHist();
    THnSparse *fRecEventHist2 = obj->GetRecEventHist2();
	
    printf("=== generated MC events for correction matrices    %lf ===\n",MCGeneratedEvents);
    printf("=== triggered MC events for correction matrices    %lf ===\n",MCTriggeredEvents);
    printf("=== recontructed MC events for correction matrices %lf ===\n",MCReconstructedEvents);
    printf("\n");
    printf("=== cut on the zVertex +- %lf ===\n",zVert);
    printf("=== generated MC events (before zVertex cut)       %lf ===\n",MCGeneratedEventsAll);
    printf("=== triggered MC events (before zVertex cut)       %lf ===\n",MCTriggeredEventsAll);
    printf("=== recontructed MC events (before zVertex cut)    %lf ===\n",MCReconstructedEventsAll);
    
    TH1D* h1MCGeneratedEvents = new TH1D("h1MCGeneratedEvents","h1MCGeneratedEvents",1,0,1);
    TH1D* h1MCTriggeredEvents = new TH1D("h1MCTriggeredEvents","h1MCTriggeredEvents",1,0,1);
    TH1D* h1MCReconstructedEvents = new TH1D("h1MCReconstructedEvents","h1MCReconstructedEvents",1,0,1);
    TH1D* h1MCGeneratedEventsAll = new TH1D("h1MCGeneratedEventsAll","h1MCGeneratedEventsAll",1,0,1);
    TH1D* h1MCTriggeredEventsAll = new TH1D("h1MCTriggeredEventsAll","h1MCTriggeredEventsAll",1,0,1);
    TH1D* h1MCReconstructedEventsAll = new TH1D("h1MCReconstructedEventsAll","h1MCReconstructedEventsAll",1,0,1);
    TH1D* h1MCTriggeredEventsAll0mult = new TH1D("h1MCTriggeredEventsAll0mult","h1MCTriggeredEventsAll0mult",1,0,1);
    TH1D* h1MCTriggeredEvents0mult = new TH1D("h1MCTriggeredEvents0mult","h1MCTriggeredEvents0mult",1,0,1);
    
    h1MCGeneratedEvents->Fill(0,MCGeneratedEvents);
    h1MCGeneratedEvents->SetEntries(MCGeneratedEvents);
    h1MCTriggeredEvents->Fill(0,MCTriggeredEvents);
    h1MCTriggeredEvents->SetEntries(MCTriggeredEvents);
    h1MCReconstructedEvents->Fill(0,MCReconstructedEvents);
    h1MCReconstructedEvents->SetEntries(MCReconstructedEvents);
    h1MCGeneratedEventsAll->Fill(0,MCGeneratedEventsAll);
    h1MCGeneratedEventsAll->SetEntries(MCGeneratedEventsAll);
    h1MCTriggeredEventsAll->Fill(0,MCTriggeredEventsAll);
    h1MCTriggeredEventsAll->SetEntries(MCTriggeredEventsAll);
    h1MCReconstructedEventsAll->Fill(0,MCReconstructedEventsAll);
    h1MCReconstructedEventsAll->SetEntries(MCReconstructedEventsAll);
    h1MCTriggeredEventsAll0mult->Fill(0,MCTriggeredEventsAll0mult);
    h1MCTriggeredEventsAll0mult->SetEntries(MCTriggeredEventsAll0mult);
    h1MCTriggeredEvents0mult->Fill(0,MCTriggeredEvents0mult);
    h1MCTriggeredEvents0mult->SetEntries(MCTriggeredEvents0mult);
    
    
    ContHist->Add(h1MCGeneratedEvents);
    ContHist->Add(h1MCTriggeredEvents);
    ContHist->Add(h1MCReconstructedEvents);
    ContHist->Add(h1MCGeneratedEventsAll);
    ContHist->Add(h1MCTriggeredEventsAll);
    ContHist->Add(h1MCReconstructedEventsAll);
    ContHist->Add(h1MCTriggeredEvents0mult);
    ContHist->Add(h1MCTriggeredEventsAll0mult);
    

    // efficienfy and correction matrices for tigger and vertex efficiency
    TH2D* h2EventTriggerEffAll  = AlidNdPtHelper::GenerateCorrMatrix(h2TriggerEventAll,h2GenEventAll,"h2EventTriggerEffAll");
    TH2D* h2EventTriggerCorrAll = AlidNdPtHelper::GenerateCorrMatrix(h2GenEventAll,h2TriggerEventAll,"h2EventTriggerCorrAll"); 
    TH2D* h2EventTriggerEff  = AlidNdPtHelper::GenerateCorrMatrix(h2TriggerEvent,h2GenEvent,"h2EventTriggerEff");
    TH2D* h2EventTriggerCorr = AlidNdPtHelper::GenerateCorrMatrix(h2GenEvent,h2TriggerEvent,"h2EventTriggerCorr"); 

    TH2D* h2EventRecEffAll  = AlidNdPtHelper::GenerateCorrMatrix(h2RecEventAll,h2TriggerEventAll,"h2EventRecEffAll");
    TH2D* h2EventRecCorrAll = AlidNdPtHelper::GenerateCorrMatrix(h2TriggerEventAll,h2RecEventAll,"h2EventRecCorrAll");
    TH2D* h2EventRecEff  = AlidNdPtHelper::GenerateCorrMatrix(h2RecEvent,h2TriggerEvent,"h2EventRecEff");
    TH2D* h2EventRecCorr = AlidNdPtHelper::GenerateCorrMatrix(h2TriggerEvent,h2RecEvent,"h2EventRecCorr");

    TH2D* h2EventEffAll  = AlidNdPtHelper::GenerateCorrMatrix(h2RecEventAll,h2GenEventAll,"h2EventEffAll");
    TH2D* h2EventCorrAll = AlidNdPtHelper::GenerateCorrMatrix(h2GenEventAll,h2RecEventAll,"h2EventCorrAll");
    TH2D* h2EventEff  = AlidNdPtHelper::GenerateCorrMatrix(h2RecEvent,h2GenEvent,"h2EventEff");
    TH2D* h2EventCorr = AlidNdPtHelper::GenerateCorrMatrix(h2GenEvent,h2RecEvent,"h2EventCorr");

    CorrMatr->Add(h2EventTriggerEffAll);
    CorrMatr->Add(h2EventTriggerCorrAll);
    CorrMatr->Add(h2EventTriggerEff);
    CorrMatr->Add(h2EventTriggerCorr);
    CorrMatr->Add(h2EventRecEffAll);
    CorrMatr->Add(h2EventRecCorrAll);
    CorrMatr->Add(h2EventRecEff);
    CorrMatr->Add(h2EventRecCorr);
    CorrMatr->Add(h2EventEffAll);
    CorrMatr->Add(h2EventCorrAll);
    CorrMatr->Add(h2EventEff);
    CorrMatr->Add(h2EventCorr);
    
    THnSparse* hSMultCorrelation = obj->GetEventMultCorrelationMatrix()->Clone("hSMultCorrelation");
    CorrMatr->Add(hSMultCorrelation);
    TH2D* h2MultCorrelation = hSMultCorrelation->Projection(0,1)->Clone("h2MultCorrelation");
    CorrMatr->Add(h2MultCorrelation);

    // all recontructed
    THnSparse *fRecTrackMatrix = obj->GetRecTrackMatrix();
    fRecTrackMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fRecTrackMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    TH3D* h3RecTrack = fRecTrackMatrix->Projection(0,1,2)->Clone("h3RecTrack");
    TH2D* h2RecTrack_zv_pt  = fRecTrackMatrix->Projection(0,1)->Clone("h2RecTrack_zv_pt");
    TH2D* h2RecTrack_zv_eta = fRecTrackMatrix->Projection(0,2)->Clone("h2RecTrack_zv_eta");
    TH2D* h2RecTrack_pt_eta = fRecTrackMatrix->Projection(1,2)->Clone("h2RecTrack_pt_eta");
    TH1D* h1RecTrack_zv  = fRecTrackMatrix->Projection(0)->Clone("h1RecTrack_zv");
    TH1D* h1RecTrack_pt  = fRecTrackMatrix->Projection(1)->Clone("h1RecTrack_pt");
    TH1D* h1RecTrack_eta = fRecTrackMatrix->Projection(2)->Clone("h1RecTrack_eta");
    Double_t MCReconstructedTracks = h3RecTrack->Integral();    

    ContHist->Add(h3RecTrack);
    ContHist->Add(h2RecTrack_zv_pt);
    ContHist->Add(h2RecTrack_zv_eta);
    ContHist->Add(h2RecTrack_pt_eta);
    ContHist->Add(h1RecTrack_zv);
    ContHist->Add(h1RecTrack_pt);
    ContHist->Add(h1RecTrack_eta);

     // recontructed primary tracks
    THnSparse *fRecPrimTrackMatrix = obj->GetRecPrimTrackMatrix();
    THnSparse *fRecTrackMatrixScaled = fRecPrimTrackMatrix->Clone("fRecTrackMatrixScaled"); //used later for secondaries scaling
    fRecPrimTrackMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fRecPrimTrackMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    TH3D* h3RecPrimTrack = fRecPrimTrackMatrix->Projection(0,1,2)->Clone("h3RecPrimTrack");
    TH2D* h2RecPrimTrack_zv_pt  = fRecPrimTrackMatrix->Projection(0,1)->Clone("h2RecPrimTrack_zv_pt");
    TH2D* h2RecPrimTrack_zv_eta = fRecPrimTrackMatrix->Projection(0,2)->Clone("h2RecPrimTrack_zv_eta");
    TH2D* h2RecPrimTrack_pt_eta = fRecPrimTrackMatrix->Projection(1,2)->Clone("h2RecPrimTrack_pt_eta");
    TH1D* h1RecPrimTrack_zv  = fRecPrimTrackMatrix->Projection(0)->Clone("h1RecPrimTrack_zv");
    TH1D* h1RecPrimTrack_pt  = fRecPrimTrackMatrix->Projection(1)->Clone("h1RecPrimTrack_pt");
    TH1D* h1RecPrimTrack_eta = fRecPrimTrackMatrix->Projection(2)->Clone("h1RecPrimTrack_eta");
    Double_t MCReconstructedPrimTracks = h3RecPrimTrack->Integral();

    ContHist->Add(h3RecPrimTrack);
    ContHist->Add(h2RecPrimTrack_zv_pt);
    ContHist->Add(h2RecPrimTrack_zv_eta);
    ContHist->Add(h2RecPrimTrack_pt_eta);
    ContHist->Add(h1RecPrimTrack_zv);
    ContHist->Add(h1RecPrimTrack_pt);
    ContHist->Add(h1RecPrimTrack_eta);
    
    // recontructed secondary tracks
    THnSparse *fRecSecTrackMatrix = obj->GetRecSecTrackMatrix();
    THnSparse *fRecSecTrackMatrixScaled = fRecSecTrackMatrix->Clone("fRecSecTrackMatrixScaled"); //used later for secondaries scaling
    fRecSecTrackMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fRecSecTrackMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    TH3D* h3RecSecTrack = fRecSecTrackMatrix->Projection(0,1,2)->Clone("h3RecSecTrack");
    TH2D* h2RecSecTrack_zv_pt  = fRecSecTrackMatrix->Projection(0,1)->Clone("h2RecSecTrack_zv_pt");
    TH2D* h2RecSecTrack_zv_eta = fRecSecTrackMatrix->Projection(0,2)->Clone("h2RecSecTrack_zv_eta");
    TH2D* h2RecSecTrack_pt_eta = fRecSecTrackMatrix->Projection(1,2)->Clone("h2RecSecTrack_pt_eta");
    TH1D* h1RecSecTrack_zv  = fRecSecTrackMatrix->Projection(0)->Clone("h1RecSecTrack_zv");
    TH1D* h1RecSecTrack_pt  = fRecSecTrackMatrix->Projection(1)->Clone("h1RecSecTrack_pt");
    TH1D* h1RecSecTrack_eta = fRecSecTrackMatrix->Projection(2)->Clone("h1RecSecTrack_eta");
    Double_t MCReconstructedSecTracks = h3RecSecTrack->Integral();

    ContHist->Add(h3RecSecTrack);
    ContHist->Add(h2RecSecTrack_zv_pt);
    ContHist->Add(h2RecSecTrack_zv_eta);
    ContHist->Add(h2RecSecTrack_pt_eta);
    ContHist->Add(h1RecSecTrack_zv);
    ContHist->Add(h1RecSecTrack_pt);
    ContHist->Add(h1RecSecTrack_eta);
    
    // generated primary tracks
    THnSparse *fGenPrimTrackMatrix = obj->GetGenPrimTrackMatrix();
    fGenPrimTrackMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fGenPrimTrackMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    TH3D* h3GenPrimTrack = fGenPrimTrackMatrix->Projection(0,1,2)->Clone("h3GenPrimTrack");
    TH2D* h2GenPrimTrack_zv_Pt  = fGenPrimTrackMatrix->Projection(0,1)->Clone("h2GenPrimTrack_zv_pt");
    TH2D* h2GenPrimTrack_zv_eta = fGenPrimTrackMatrix->Projection(0,2)->Clone("h2GenPrimTrack_zv_eta");
    TH2D* h2GenPrimTrack_pt_eta = fGenPrimTrackMatrix->Projection(1,2)->Clone("h2GenPrimTrack_pt_eta");
    TH1D* h1GenPrimTrack_zv  = fGenPrimTrackMatrix->Projection(0)->Clone("h1GenPrimTrack_zv");
    TH1D* h1GenPrimTrack_pt  = fGenPrimTrackMatrix->Projection(1)->Clone("h1GenPrimTrack_pt");
    TH1D* h1GenPrimTrack_eta = fGenPrimTrackMatrix->Projection(2)->Clone("h1GenPrimTrack_eta");
    Double_t MCGeneratedPrimTracks = h3GenPrimTrack->Integral();

    ContHist->Add(h3GenPrimTrack);
    ContHist->Add(h2GenPrimTrack_zv_pt);
    ContHist->Add(h2GenPrimTrack_zv_eta);
    ContHist->Add(h2GenPrimTrack_pt_eta);
    ContHist->Add(h1GenPrimTrack_zv);
    ContHist->Add(h1GenPrimTrack_pt);
    ContHist->Add(h1GenPrimTrack_eta);
    printf("\n");
    printf("==============================================================\n");    
    printf("=== recontructed MC tracks              %lf ===\n",MCReconstructedTracks);
    printf("=== recontructed MC secondary tracks    %lf ===\n",MCReconstructedSecTracks);
    printf("=== recontructed MC primary tracks      %lf ===\n",MCReconstructedPrimTracks);
    printf("=== generated MC primary track          %lf ===\n",MCGeneratedPrimTracks);
    printf("==============================================================\n");    
    printf("\n");
    
    // mc truth histogram (for self-consistency check)
    TH1D* dNdPt_MC   = (TH1D*) h1GenPrimTrack_pt->Clone("dNdPt_MC");
    TH1D* dNdEta_MC   =(TH1D*) h1GenPrimTrack_eta->Clone("dNdEta_MC");
    // normalization and finalization
    // 1/N_evt 1/(2 pi pt) 1/width 1/etarange
   for (int ii=1; ii <= dNdPt_MC->GetNbinsX() ;ii++) {
        Double_t pt = dNdPt_MC->GetBinCenter(ii);
        Double_t width = dNdPt_MC->GetBinWidth(ii);
        Double_t val = dNdPt_MC->GetBinContent(ii);
        Double_t err = dNdPt_MC->GetBinError(ii);        
        Double_t cval = (val)/(width * 2.0 * TMath::Pi() * 1.6 * MCGeneratedEvents * pt);
        Double_t cerr = (err)/(width * 2.0 * TMath::Pi() * 1.6 * MCGeneratedEvents * pt);
        dNdPt_MC->SetBinContent(ii,cval);
        dNdPt_MC->SetBinError(ii,cerr);
    }    
    dNdPt_MC->SetMarkerStyle(21);
    dNdPt_MC->SetTitle("; p_{T} (GeV/c) ; 1/N_{evt} 1/(2#pi p_{T}) (d^{2}N_{ch})/(d#eta dp_{T})^{-2}");
    dNdEta_MC->Scale(1,"width");
    dNdEta_MC->Scale(1./MCGeneratedEvents);        
    ContHist->Add(dNdPt_MC);
    ContHist->Add(dNdEta_MC);
	

   // Rebin for corrections (pt)
   cout << "rebinning...." << endl;
   
   TH1D* h1RecPrimTrack_pt_Rebin = h1RecPrimTrack_pt->Rebin(ptNbins,"h1RecPrimTrack_pt_Rebin",binsPt);
   TH1D* h1GenPrimTrack_pt_Rebin = h1GenPrimTrack_pt->Rebin(ptNbins,"h1GenPrimTrack_pt_Rebin",binsPt);
   ContHist->Add(h1RecPrimTrack_pt_Rebin);
   ContHist->Add(h1GenPrimTrack_pt_Rebin);

	
//    THnSparse *fSparseTriggerTrackEvent = obj->GetTriggerTrackEventMatrix();//Tracks from triggered events
//    THnSparse *fSparseVtxTrackEvent = obj->GetRecTrackEventMatrix();//Tracks from events with rec. vtx
//    THnSparse *fSparseGenTrackEvent = obj->GetGenTrackEventMatrix();//generated TrackEvent matrix

    // tracking efficiencies + corrections  
   TH2D* h2TrackEff_zv_pt   = AlidNdPtHelper::GenerateCorrMatrix(h2RecPrimTrack_zv_pt,h2GenPrimTrack_zv_pt,"h2TrackEff_zv_pt");
   TH2D* h2TrackCorr_zv_pt  = AlidNdPtHelper::GenerateCorrMatrix(h2GenPrimTrack_zv_pt,h2RecPrimTrack_zv_pt,"h2TrackCorr_zv_pt");
   TH2D* h2TrackEff_zv_eta  = AlidNdPtHelper::GenerateCorrMatrix(h2RecPrimTrack_zv_eta,h2GenPrimTrack_zv_eta,"h2TrackEff_zv_eta");
   TH2D* h2TrackCorr_zv_eta = AlidNdPtHelper::GenerateCorrMatrix(h2GenPrimTrack_zv_eta,h2RecPrimTrack_zv_eta,"h2TrackCorr_zv_eta");
   TH2D* h2TrackEff_pt_eta  = AlidNdPtHelper::GenerateCorrMatrix(h2RecPrimTrack_pt_eta,h2GenPrimTrack_pt_eta,"h2TrackEff_pt_eta");
   TH2D* h2TrackCorr_pt_eta = AlidNdPtHelper::GenerateCorrMatrix(h2GenPrimTrack_pt_eta,h2RecPrimTrack_pt_eta,"h2TrackCorr_pt_eta");
  
    
   TH1D* h1TrackEff_zv   = AlidNdPtHelper::GenerateCorrMatrix(h1RecPrimTrack_zv,h1GenPrimTrack_zv,"h1TrackEff_zv");
   TH1D* h1TrackCorr_zv  = AlidNdPtHelper::GenerateCorrMatrix(h1GenPrimTrack_zv,h1RecPrimTrack_zv,"h1TrackCorr_zv");
   TH1D* h1TrackEff_pt   = AlidNdPtHelper::GenerateCorrMatrix(h1RecPrimTrack_pt_Rebin,h1GenPrimTrack_pt_Rebin,"h1TrackEff_pt");
   TH1D* h1TrackCorr_pt  = AlidNdPtHelper::GenerateCorrMatrix(h1GenPrimTrack_pt_Rebin,h1RecPrimTrack_pt_Rebin,"h1TrackCorr_pt");
   TH1D* h1TrackEff_eta  = AlidNdPtHelper::GenerateCorrMatrix(h1RecPrimTrack_eta,h1GenPrimTrack_eta,"h1TrackEff_eta");
   TH1D* h1TrackCorr_eta = AlidNdPtHelper::GenerateCorrMatrix(h1GenPrimTrack_eta,h1RecPrimTrack_eta,"h1TrackCorr_eta");
   CorrMatr->Add(h2TrackEff_zv_pt);
   CorrMatr->Add(h2TrackCorr_zv_pt);
   CorrMatr->Add(h2TrackEff_zv_eta);
   CorrMatr->Add(h2TrackCorr_zv_eta);
   CorrMatr->Add(h2TrackEff_pt_eta);
   CorrMatr->Add(h2TrackCorr_pt_eta);
   CorrMatr->Add(h1TrackEff_zv);
   CorrMatr->Add(h1TrackCorr_zv);
   CorrMatr->Add(h1TrackEff_pt);
   CorrMatr->Add(h1TrackCorr_pt);
   CorrMatr->Add(h1TrackEff_eta);
   CorrMatr->Add(h1TrackCorr_eta);

   // scale the secondaries before calculating correction matrices
   for (Long64_t i = 0; i < fRecSecTrackMatrixScaled->GetNbins(); i++) {
       Int_t c[3];
       Double_t val = fRecSecTrackMatrixScaled->GetBinContent(i,c);
       Double_t err = fRecSecTrackMatrixScaled->GetBinError(i);
       Double_t pt =  fRecSecTrackMatrixScaled->GetAxis(1)->GetBinCenter(c[1]);
       Double_t scale = GetStrangenessCorrFactor(pt,sscale);
//        Double_t scale = AlidNdPtHelper::GetStrangenessCorrFactor(pt);
       fRecSecTrackMatrixScaled->SetBinContent(c,val*scale);
       fRecSecTrackMatrixScaled->SetBinError(c,err*scale);
    }
    
    // for correct determination of secondaries contamination, also the total total tracks have to be scaled
    // this is done by taking primaries and adding the scaled secondaries
    fRecTrackMatrixScaled->Add(fRecSecTrackMatrixScaled);

    fRecSecTrackMatrixScaled->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fRecSecTrackMatrixScaled->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    
    TH3D* h3RecSecTrackScaled = fRecSecTrackMatrixScaled->Projection(0,1,2)->Clone("h3RecSecTrackScaled");
    TH2D* h2RecSecTrackScaled_zv_pt  = fRecSecTrackMatrixScaled->Projection(0,1)->Clone("h2RecSecTrackScaled_zv_pt");
    TH2D* h2RecSecTrackScaled_zv_eta = fRecSecTrackMatrixScaled->Projection(0,2)->Clone("h2RecSecTrackScaled_zv_eta");
    TH2D* h2RecSecTrackScaled_pt_eta = fRecSecTrackMatrixScaled->Projection(1,2)->Clone("h2RecSecTrackScaled_pt_eta");
    TH1D* h1RecSecTrackScaled_zv  = fRecSecTrackMatrixScaled->Projection(0)->Clone("h1RecSecTrackScaled_zv");
    TH1D* h1RecSecTrackScaled_pt  = fRecSecTrackMatrixScaled->Projection(1)->Clone("h1RecSecTrackScaled_pt");
    TH1D* h1RecSecTrackScaled_eta = fRecSecTrackMatrixScaled->Projection(2)->Clone("h1RecSecTrackScaled_eta");

    ContHist->Add(h3RecSecTrackScaled);
    ContHist->Add(h2RecSecTrackScaled_zv_pt);
    ContHist->Add(h2RecSecTrackScaled_zv_eta);
    ContHist->Add(h2RecSecTrackScaled_pt_eta);
    ContHist->Add(h1RecSecTrackScaled_zv);
    ContHist->Add(h1RecSecTrackScaled_pt);
    ContHist->Add(h1RecSecTrackScaled_eta);    
    
    fRecTrackMatrixScaled->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    fRecTrackMatrixScaled->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta
    
    TH3D* h3RecTrackScaled = fRecTrackMatrixScaled->Projection(0,1,2)->Clone("h3RecTrackScaled");
    TH2D* h2RecTrackScaled_zv_pt  = fRecTrackMatrixScaled->Projection(0,1)->Clone("h2RecTrackScaled_zv_pt");
    TH2D* h2RecTrackScaled_zv_eta = fRecTrackMatrixScaled->Projection(0,2)->Clone("h2RecTrackScaled_zv_eta");
    TH2D* h2RecTrackScaled_pt_eta = fRecTrackMatrixScaled->Projection(1,2)->Clone("h2RecTrackScaled_pt_eta");
    TH1D* h1RecTrackScaled_zv  = fRecTrackMatrixScaled->Projection(0)->Clone("h1RecTrackScaled_zv");
    TH1D* h1RecTrackScaled_pt  = fRecTrackMatrixScaled->Projection(1)->Clone("h1RecTrackScaled_pt");
    TH1D* h1RecTrackScaled_eta = fRecTrackMatrixScaled->Projection(2)->Clone("h1RecTrackScaled_eta");

    ContHist->Add(h3RecTrackScaled);
    ContHist->Add(h2RecTrackScaled_zv_pt);
    ContHist->Add(h2RecTrackScaled_zv_eta);
    ContHist->Add(h2RecTrackScaled_pt_eta);
    ContHist->Add(h1RecTrackScaled_zv);
    ContHist->Add(h1RecTrackScaled_pt);
    ContHist->Add(h1RecTrackScaled_eta);
    
   // Rebin for corrections (pt)
   TH1D* h1RecTrackScaled_pt_Rebin = h1RecTrackScaled_pt->Rebin(ptNbins,"h1RecTrackScaled_pt_Rebin",binsPt);
   TH1D* h1RecSecTrackScaled_pt_Rebin = h1RecSecTrackScaled_pt->Rebin(ptNbins,"h1RecSecTrackScaled_pt_Rebin",binsPt);
   ContHist->Add(h1RecTrackScaled_pt_Rebin);
   ContHist->Add(h1RecSecTrackScaled_pt_Rebin);
   
    // create histograms for secondaries contamination and correction
    
    TH2D* h2SecCont_zv_pt  = AlidNdPtHelper::GenerateCorrMatrix(h2RecSecTrackScaled_zv_pt,h2RecTrackScaled_zv_pt,"h2SecCont_zv_pt");
    TH2D* h2SecCorr_zv_pt  = AlidNdPtHelper::GenerateContCorrMatrix(h2RecSecTrackScaled_zv_pt,h2RecTrackScaled_zv_pt,"h2SecCorr_zv_pt");
    TH2D* h2SecCont_zv_eta = AlidNdPtHelper::GenerateCorrMatrix(h2RecSecTrackScaled_zv_eta,h2RecTrackScaled_zv_eta,"h2SecCont_zv_eta");
    TH2D* h2SecCorr_zv_eta = AlidNdPtHelper::GenerateContCorrMatrix(h2RecSecTrackScaled_zv_eta,h2RecTrackScaled_zv_eta,"h2SecCorr_zv_eta");
    TH2D* h2SecCont_pt_eta = AlidNdPtHelper::GenerateCorrMatrix(h2RecSecTrackScaled_pt_eta,h2RecTrackScaled_pt_eta,"h2SecCont_pt_eta");
    TH2D* h2SecCorr_pt_eta = AlidNdPtHelper::GenerateContCorrMatrix(h2RecSecTrackScaled_pt_eta,h2RecTrackScaled_pt_eta,"h2SecCorr_pt_eta");
    TH1D* h1SecCont_zv = AlidNdPtHelper::GenerateCorrMatrix(h1RecSecTrackScaled_zv,h1RecTrackScaled_zv,"h1SecCont_zv");
    TH1D* h1SecCorr_zv = AlidNdPtHelper::GenerateContCorrMatrix(h1RecSecTrackScaled_zv,h1RecTrackScaled_zv,"h1SecCorr_zv");
    TH1D* h1SecCont_pt = AlidNdPtHelper::GenerateCorrMatrix(h1RecSecTrackScaled_pt_Rebin,h1RecTrackScaled_pt_Rebin,"h1SecCont_pt");
    TH1D* h1SecCorr_pt = AlidNdPtHelper::GenerateContCorrMatrix(h1RecSecTrackScaled_pt_Rebin,h1RecTrackScaled_pt_Rebin,"h1SecCorr_pt");
    TH1D* h1SecCont_eta = AlidNdPtHelper::GenerateCorrMatrix(h1RecSecTrackScaled_eta,h1RecTrackScaled_eta,"h1SecCont_eta");
    TH1D* h1SecCorr_eta = AlidNdPtHelper::GenerateContCorrMatrix(h1RecSecTrackScaled_eta,h1RecTrackScaled_eta,"h1SecCorr_eta");

    CorrMatr->Add(h2SecCont_zv_pt);
    CorrMatr->Add(h2SecCorr_zv_pt);
    CorrMatr->Add(h2SecCont_zv_eta);
    CorrMatr->Add(h2SecCorr_zv_eta);
    CorrMatr->Add(h2SecCont_pt_eta);
    CorrMatr->Add(h2SecCorr_pt_eta);
    CorrMatr->Add(h1SecCont_zv);
    CorrMatr->Add(h1SecCorr_zv);
    CorrMatr->Add(h1SecCont_pt);
    CorrMatr->Add(h1SecCorr_pt);
    CorrMatr->Add(h1SecCont_eta);
    CorrMatr->Add(h1SecCorr_eta);
    
    //vertex distribution for events with no vertex
    fGenEventMatrix->GetAxis(0)->SetRange();
    fGenEventMatrix->GetAxis(1)->SetRange(1,1);
    TH1D* h1GenEventMatrix_bin0_zv = (TH1D*) fGenEventMatrix->Projection(0)->Clone("h1GenEventMatrix_bin0_zv");
    ContHist->Add(h1GenEventMatrix_bin0_zv);
    
    fTriggerEventMatrix->GetAxis(0)->SetRange();
    fTriggerEventMatrix->GetAxis(1)->SetRange(1,1);
    TH1D* h1TriggerEventMatrix_bin0_zv = (TH1D*) fTriggerEventMatrix->Projection(0)->Clone("h1TriggerEventMatrix_bin0_zv");
    ContHist->Add(h1TriggerEventMatrix_bin0_zv);


    fRecEventMatrix->GetAxis(0)->SetRange();
    fRecEventMatrix->GetAxis(1)->SetRange();
    TH1D* h1RecEventMatrix_zv = (TH1D*) fRecEventMatrix->Projection(0)->Clone("h1RecEventMatrix_zv");
    h1RecEventMatrix_zv->Scale(1./h1RecEventMatrix_zv->Integral());
    //ContHist->Add(h1TriggerEventMatrix_zv);
   
    TH1D* h1TriggerEff_bin0_zv   = AlidNdPtHelper::GenerateCorrMatrix(h1TriggerEventMatrix_bin0_zv,h1GenEventMatrix_bin0_zv,"h1TriggerEff_bin0_zv");
    TH1D* h1TriggerCorr_bin0_zv  = AlidNdPtHelper::GenerateCorrMatrix(h1GenEventMatrix_bin0_zv,h1TriggerEventMatrix_bin0_zv,"h1TriggerCorr_bin0_zv");
    CorrMatr->Add(h1TriggerEff_bin0_zv);
    CorrMatr->Add(h1TriggerCorr_bin0_zv);
    
    TH1D* h1Ratio_zv = (TH1D*) h1GenEventMatrix_bin0_zv->Clone("h1Ratio_zv");
    h1Ratio_zv->Scale(1./h1Ratio_zv->Integral());
    h1Ratio_zv->Divide(h1RecEventMatrix_zv);
    CorrMatr->Add(h1Ratio_zv);
    
    
     // aded to correct sperately for zvert shape of triggered and untriggered bin0 events
     fRecEventMatrix->GetAxis(0)->SetRange();
     fRecEventMatrix->GetAxis(1)->SetRange();
     fTriggerEventMatrix->GetAxis(1)->SetRangeUser(0,0);
     fGenEventMatrix->GetAxis(1)->SetRangeUser(0,0);
     
     TH1D* rec_zv = fRecEventMatrix->Projection(0)->Clone("rec_zv");
     TH1D* trig0_zv = fTriggerEventMatrix->Projection(0)->Clone("trig0_zv");
     TH1D* corr_shape_trig0_zv = AlidNdPtHelper::GenerateCorrMatrix(trig0_zv,rec_zv,"corr_shape_trig0_zv");
     CorrMatr->Add(corr_shape_trig0_zv);     
     
     TH1D* gen0_zv = fGenEventMatrix->Projection(0)->Clone("gen0_zv");
     TH1D* gen0notrig_zv = gen0_zv->Clone("gen0notrig_zv");
     gen0notrig_zv->Add(trig0_zv,-1.);
     TH1D* corr_shape_notrig0_zv = AlidNdPtHelper::GenerateCorrMatrix(gen0notrig_zv,rec_zv,"corr_shape_notrig0_zv");
     CorrMatr->Add(corr_shape_notrig0_zv); 
         
         
   // efficienfy and correction for THnSparse
    THnSparse* hSRecSecTrackMatrixScaled = AlidNdPtHelper::RebinTHnSparse(fRecSecTrackMatrixScaled,hSTrackSecCorr,"hSRecSecTrackMatrixScaled");
    THnSparse* hSRecTrackMatrixScaled = AlidNdPtHelper::RebinTHnSparse(fRecTrackMatrixScaled,hSTrackSecCorr,"hSRecTrackMatrixScaled");
    THnSparse* hSRecPrimTrackMatrix = AlidNdPtHelper::RebinTHnSparse(fRecPrimTrackMatrix,hSTrackEffCorr,"hSRecPrimTrackMatrix");
    THnSparse* hSGenPrimTrackMatrix = AlidNdPtHelper::RebinTHnSparse(fGenPrimTrackMatrix,hSTrackEffCorr,"hSGenPrimTrackMatrix");
    ContHist->Add(hSRecSecTrackMatrixScaled);
    ContHist->Add(hSRecTrackMatrixScaled);
    ContHist->Add(hSRecPrimTrackMatrix);
    ContHist->Add(hSGenPrimTrackMatrix);
    
    THnSparse* hSTrackEff  = AlidNdPtHelper::GenerateCorrMatrix(hSRecPrimTrackMatrix,hSGenPrimTrackMatrix,"hSTrackEff");
    THnSparse* hSTrackCorr = AlidNdPtHelper::GenerateCorrMatrix(hSGenPrimTrackMatrix,hSRecPrimTrackMatrix,"hSTrackCorr");
    THnSparse* hSSecCont   = AlidNdPtHelper::GenerateCorrMatrix(hSRecSecTrackMatrixScaled,hSRecTrackMatrixScaled,"hSSecCont");
    THnSparse* hSSecCorr   = AlidNdPtHelper::GenerateContCorrMatrix(hSRecSecTrackMatrixScaled,hSRecTrackMatrixScaled,"hSSecCorr");
    CorrMatr->Add(hSTrackEff);
    CorrMatr->Add(hSTrackCorr);
    CorrMatr->Add(hSSecCont);
    CorrMatr->Add(hSSecCorr);
    
    // create th3 from thnsparse, used for corrections
    TH3D* h3RecPrimTrack_Rebin = hSRecPrimTrackMatrix->Projection(0,1,2)->Clone("h3RecPrimTrack_Rebin");
    TH3D* h3GenPrimTrack_Rebin = hSGenPrimTrackMatrix->Projection(0,1,2)->Clone("h3GenPrimTrack_Rebin");
    TH3D* h3RecSecTrackScaled_Rebin = hSRecSecTrackMatrixScaled->Projection(0,1,2)->Clone("h3RecSecTrackScaled_Rebin");
    TH3D* h3RecTrackScaled_Rebin = hSRecTrackMatrixScaled->Projection(0,1,2)->Clone("h3RecTrackScaled_Rebin");
   
    TH3D* h3TrackEff  = AlidNdPtHelper::GenerateCorrMatrix(h3RecPrimTrack_Rebin,h3GenPrimTrack_Rebin,"h3TrackEff");
    TH3D* h3TrackCorr = AlidNdPtHelper::GenerateCorrMatrix(h3GenPrimTrack_Rebin,h3RecPrimTrack_Rebin,"h3TrackCorr");    
    TH3D* h3SecCont = AlidNdPtHelper::GenerateCorrMatrix(h3RecSecTrackScaled_Rebin,h3RecTrackScaled_Rebin,"h3SecCont");
    TH3D* h3SecCorr = AlidNdPtHelper::GenerateContCorrMatrix(h3RecSecTrackScaled_Rebin,h3RecTrackScaled_Rebin,"h3SecCorr");
    CorrMatr->Add(h3TrackEff);
    CorrMatr->Add(h3TrackCorr);
    CorrMatr->Add(h3SecCont);
    CorrMatr->Add(h3SecCorr);
 
    // check for trigger/vertex bias on track level
    THnSparse* hSRecTrackEventMatrix = (THnSparse*) obj->GetRecTrackEventMatrix()->Clone("hSRecTrackEventMatrix");
    THnSparse* hSGenTrackEventMatrix = (THnSparse*) obj->GetGenTrackEventMatrix()->Clone("hSGenTrackEventMatrix");
    THnSparse* hSTriggerTrackEventMatrix = (THnSparse*) obj->GetTriggerTrackEventMatrix()->Clone("hSTriggerTrackEventMatrix");
    
    hSRecTrackEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    hSRecTrackEventMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta    
    hSGenTrackEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    hSGenTrackEventMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta    
    hSTriggerTrackEventMatrix->GetAxis(0)->SetRangeUser(-zVert, zVert-0.01);//zVer
    hSTriggerTrackEventMatrix->GetAxis(2)->SetRangeUser(-eta, eta-0.01);//eta    
    
    
    TH1D* h1TriggerBiasCorr_zv = AlidNdPtHelper::GenerateCorrMatrix(hSTriggerTrackEventMatrix->Projection(0),hSGenTrackEventMatrix->Projection(0),"h1TriggerBiasCorr_zv");
    TH1D* h1TriggerBiasCorr_pt = AlidNdPtHelper::GenerateCorrMatrix(hSTriggerTrackEventMatrix->Projection(1),hSGenTrackEventMatrix->Projection(1),"h1TriggerBiasCorr_pt");
    TH1D* h1TriggerBiasCorr_eta = AlidNdPtHelper::GenerateCorrMatrix(hSTriggerTrackEventMatrix->Projection(2),hSGenTrackEventMatrix->Projection(2),"h1TriggerBiasCorr_eta");
    
    TH1D* h1VertexBiasCorr_zv  = AlidNdPtHelper::GenerateCorrMatrix(hSRecTrackEventMatrix->Projection(0),hSTriggerTrackEventMatrix->Projection(0),"h1VertexBiasCorr_zv");
    TH1D* h1VertexBiasCorr_pt  = AlidNdPtHelper::GenerateCorrMatrix(hSRecTrackEventMatrix->Projection(1),hSTriggerTrackEventMatrix->Projection(1),"h1VertexBiasCorr_pt");   
    TH1D* h1VertexBiasCorr_eta  = AlidNdPtHelper::GenerateCorrMatrix(hSRecTrackEventMatrix->Projection(2),hSTriggerTrackEventMatrix->Projection(2),"h1VertexBiasCorr_eta");    
    
   ContHist->Add(h1TriggerBiasCorr_zv);    
   ContHist->Add(h1TriggerBiasCorr_pt);    
   ContHist->Add(h1TriggerBiasCorr_eta);    
   ContHist->Add(h1VertexBiasCorr_zv);    
   ContHist->Add(h1VertexBiasCorr_pt);    
   ContHist->Add(h1VertexBiasCorr_eta);    

 
    
    // plot pictures and save to gifdir
    for (i=0; i < CorrMatr->LastIndex(); i++) {    
        TCanvas* ctmp = PlotHist(CorrMatr->At(i),idstring);
        if (gifdir && ctmp) {
            TString gif(gifdir);
            gif += '/';
            gif += ctmp->GetName();
            gif += ".gif";
            ctmp->SaveAs(gif.Data(),"gif");     
            delete ctmp;
        }
    }
    for (i=0; i < ContHist->LastIndex(); i++) {    
        TCanvas* ctmp = PlotHist(ContHist->At(i),idstring);
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
    CorrMatr->Write();
    ContHist->Write();
    out->Close();
    
    return MCReconstructedEvents;

}


//_____________________________________________________________________________
Double_t GetStrangenessCorrFactor(Double_t pt, Double_t s)
{
    // data driven correction factor for secondaries (PbPb)

    if (pt <= 0.17) return 1.0;
    if (pt <= 0.4) return GetLinearInterpolationValue(0.17,1.0,0.4,1.07, pt);
    if (pt <= 0.6) return GetLinearInterpolationValue(0.4,1.07,0.6,1.25, pt);
    if (pt <= 1.2) return GetLinearInterpolationValue(0.6,1.25,1.2,1.5,  pt);
    return 1.5;
}


//___________________________________________________________________________
Double_t GetLinearInterpolationValue(const Double_t x1,const  Double_t y1,const  Double_t x2,const  Double_t y2, const Double_t pt)
{
    //
    // linear interpolation
    //
    return ((y2-y1)/(x2-x1))*pt+(y2-(((y2-y1)/(x2-x1))*x2)); 
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
        if (xlabel.Contains("Pt")) { c->SetLogx();  h->GetXaxis()->SetRangeUser(0.1 , 100.); }
    }
    if (h->GetDimension() == 1) {
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

Int_t CheckLoadLibrary(const char* library)
{
  // checks if a library is already loaded, if not loads the library

  if (strlen(gSystem->GetLibraries(Form("%s.so", library), "", kFALSE)) > 0)
    return 1;

  return gSystem->Load(library);
}
