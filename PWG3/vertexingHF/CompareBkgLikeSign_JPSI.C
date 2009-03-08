void CompareBkgLikeSign_JPSI(const char *aodFileName="AliAOD.root",
			     const char *aodHFFileName="AliAOD.VertexingHF.root") 
{
  //
  // Example macro to read D0->Kpi, JPSI->ee, Like Sign 
  // candidates from AOD (having the standard AOD + 
  // a friend heavy-flavour AOD) and apply cuts
  // C. Di Giglio
  //
  Bool_t useParFiles=kFALSE;
  gROOT->LoadMacro("$ALICE_ROOT/PWG3/vertexingHF/LoadLibraries.C");
  LoadLibraries(useParFiles);

  //create histograms
  TH2F *hCPtaVSd0d0 = new TH2F("hCPtaVSd0d0","Like Sign pairs correlation plot",1000,-50000,50000,1000,-1,1);
  hCPtaVSd0d0->SetXTitle("Product of impact parameters [#mu m^{2}]");
  hCPtaVSd0d0->SetYTitle("Cosine of pointing angle");
  TH1F *hMass = new TH1F("hMass","Like Sign mass plot",100,0.,3.2);
  hMass->SetXTitle("Invariant mass [GeV]");
  hMass->SetYTitle("Entries");
  hMass->Sumw2();
  TH1F *hMassJPSI = new TH1F("hMassJPSI","JPSI mass plot",100,0.,3.2);
  hMassJPSI->SetXTitle("Invariant mass [GeV]");
  hMassJPSI->SetYTitle("Entries");
  hMassJPSI->Sumw2();
  TH1F *hSecVtxZ = new TH1F("hSecVtxZ","LikeSign decay vertex z",1000,-10,10);
  hSecVtxZ->SetXTitle("z of decay vertex [cm]");
  hSecVtxZ->SetYTitle("Entries");
  // Cosine theta star
  TH1F *hCtsJPSI = new TH1F("hCtsJPSI","JPSI cosine of decay angle",50,-1.,1.);
  hCtsJPSI->SetXTitle("cos #theta^{*}");
  hCtsJPSI->SetYTitle("Entries");
  hCtsJPSI->Sumw2();
  TH1F *hCtsLikeSign = new TH1F("hCtsLikeSign","LikeSign cosine of decay angle",50,-1.,1.);
  hCtsLikeSign->SetXTitle("cos #theta^{*}");
  hCtsLikeSign->SetYTitle("Entries");
  hCtsLikeSign->Sumw2();
  TH1F *hCtsLikeSignPos = new TH1F("hCtsLikeSignPos","LikeSign cosine of decay angle for ++ pairs",50,-1.,1.);
  hCtsLikeSignPos->SetXTitle("cos #theta^{*}");
  hCtsLikeSignPos->SetYTitle("Entries");
  hCtsLikeSignPos->Sumw2();
  TH1F *hCtsLikeSignNeg = new TH1F("hCtsLikeSignNeg","LikeSign cosine of decay angle for -- pairs",50,-1.,1.);
  hCtsLikeSignNeg->SetXTitle("cos #theta^{*}");
  hCtsLikeSignNeg->SetYTitle("Entries");
  hCtsLikeSignNeg->Sumw2();
  // Cosine of pointing angle
  TH1F *hCPtaJPSI = new TH1F("hCPtaJPSI","JPSI cosine of pointing angle distribution",100,-1.,1.);
  hCPtaJPSI->SetXTitle("cos #theta_{point}");
  hCPtaJPSI->SetYTitle("Entries");
  hCPtaJPSI->Sumw2();
  TH1F *hCPtaLikeSign = new TH1F("hCPtaLikeSign","LikeSign cosine of pointing angle distribution",100,-1.,1.);
  hCPtaLikeSign->SetXTitle("cos #theta_{point}");
  hCPtaLikeSign->SetYTitle("Entries");
  hCPtaLikeSign->Sumw2();
  // d0 x d0
  TH1F *hd0d0JPSI = new TH1F("hd0d0JPSI","JPSI Product of the impact parameters",100,-100000,100000);
  hd0d0JPSI->SetXTitle("d_{0}^{k} #times d_{0}^{#pi} [#mu m^{2}]");
  hd0d0JPSI->SetYTitle("Entries");
  hd0d0JPSI->Sumw2();
  TH1F *hd0d0LikeSign = new TH1F("hd0d0LikeSign","LikeSign Product of the impact parameters",100,-100000,100000);
  hd0d0LikeSign->SetXTitle("d_{0}^{+/-} #times d_{0}^{+/-} [#mu m^{2}]");
  hd0d0LikeSign->SetYTitle("Entries");
  hd0d0LikeSign->Sumw2();
  // DCA
  TH1F *hDCAJPSI = new TH1F("hDCAJPSI","JPSI distance of closest approach",100,0.,5.);
  hDCAJPSI->SetXTitle("dca [10^{2}#mu m]");
  hDCAJPSI->SetYTitle("Entries*10^{2}");
  hDCAJPSI->Sumw2();
  TH1F *hDCALikeSign = new TH1F("hDCALikeSign","LikeSign distance of closest approach",100,0.,5.);
  hDCALikeSign->SetXTitle("dca [10^{2}#mu m]");
  hDCALikeSign->SetYTitle("Entries*10^{2}");
  hDCALikeSign->Sumw2();

  Int_t nTotHF=0,nTotLikeSign=0;
  //Double_t normalizationFactorAvg=0;
  Double_t weightLS=1.;
  Double_t nTotPosPairs=0;
  Double_t nTotNegPairs=0;
  //Double_t nTotPosPairsAvg=0;
  //Double_t nTotNegPairsAvg=0;//comment in case of one directory processing
  Int_t totEvents=0;

  // open input file and get the TTree
  TFile inFile(aodFileName,"READ");
  if (!inFile.IsOpen()) return;

  TTree *aodTree = (TTree*)inFile.Get("aodTree");
  aodTree->AddFriend("aodTree",aodHFFileName);

  AliAODEvent *aod = new AliAODEvent();

  aod->ReadFromTree(aodTree);

  // load heavy flavour vertices
  TClonesArray *arrayVerticesHF = 
    (TClonesArray*)aod->GetList()->FindObject("VerticesHF"); 
    
  // load like sign candidates
  TClonesArray *arrayLikeSign =
    (TClonesArray*)aod->GetList()->FindObject("LikeSign2Prong");

  // load B->JPSI->e+e- candidates
  TClonesArray *arrayBtoJpsiToEle =
    (TClonesArray*)aod->GetList()->FindObject("JPSItoEle");

  Bool_t normalize = kTRUE;

  Double_t cuts[9]=
  // cuts[0] = inv. mass half width [GeV]
  // cuts[1] = dca [cm]
  // cuts[2] = cosThetaStar (negative electron)
  // cuts[3] = pTP [GeV/c]
  // cuts[4] = pTN [GeV/c]
  // cuts[5] = d0P [cm]   upper limit!
  // cuts[6] = d0N [cm]  upper limit!
  // cuts[7] = d0d0 [cm^2]
  // cuts[8] = cosThetaPoint
                     {1000.,
		      100000.,
		      1.1,
		      0.,
		      0.,
		      100000.,
		      100000.,
		      100000000.,
		      -1.1}; 
  //Double_t cuts[9]= {1000., 100000., 0.8, 0., 0., 100000., 100000., 80000., 0.8};

  Int_t nTotBtoJpsiToEle=0;
  Int_t nTotD0toKpi=0;
  AliAODVertex *vtx1=0;

  // loop over events
  Int_t nEvents = aodTree->GetEntries();
  for (Int_t nEv = 0; nEv < nEvents; nEv++) {
    cout<<"\n------------ Event: "<<nEv<<" ------------------"<<endl;

    // read event
    aodTree->GetEvent(nEv);

    //print event info
    aod->GetHeader()->Print();

    // primary vertex
    vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    vtx1->Print();

    // make trkIDtoEntry register (temporary)
    Int_t trkIDtoEntry[100000];
    for(Int_t it=0;it<aod->GetNumberOfTracks();it++) {
      AliAODTrack *track = aod->GetTrack(it);
      trkIDtoEntry[track->GetID()]=it;
    }
    
    // loop over Like sign candidates
    Int_t nLikeSign = arrayLikeSign->GetEntriesFast();
    nTotLikeSign += nLikeSign;
    cout<<"Number of LikeSign pairs: "<<nLikeSign<<endl;
    Int_t nPosPairs=0,nNegPairs=0;
    
    for (Int_t iLikeSign = 0; iLikeSign < nLikeSign; iLikeSign++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayLikeSign->UncheckedAt(iLikeSign);
      Bool_t unsetvtx=kFALSE;
      if(!d->GetOwnPrimaryVtx()) {
	d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
	unsetvtx=kTRUE;
      }
      Int_t okBtoJPSIls=0;
      if(d->SelectBtoJPSI(cuts,okBtoJPSIls)) {
	hMass->Fill(d->InvMassJPSIee());
	hCPtaVSd0d0->Fill(1e8*d->Prodd0d0(),d->CosPointingAngle());
 	hSecVtxZ->Fill(d->GetSecVtxZ());
        hCPtaLikeSign->Fill(d->CosPointingAngle());
        hd0d0LikeSign->Fill(1e8*d->Prodd0d0());
        hCtsLikeSign->Fill(d->CosThetaStarJPSI());
	hDCALikeSign->Fill(100*d->GetDCA());
        AliAODTrack *trk0 = (AliAODTrack*)d->GetDaughter(0);
 	if(!trk0) {
	  trk0=aod->GetTrack(trkIDtoEntry[d->GetProngID(0)]);
	  cout<<"references to standard AOD not available"<<endl;
	}
        if((trk0->Charge())==1) {
	  nPosPairs++;
          hCtsLikeSignPos->Fill(d->CosThetaStarJPSI()); 
	} else {
	  nNegPairs++;
          hCtsLikeSignNeg->Fill(d->CosThetaStarJPSI());
	}
      }
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
    }
    cout<<"\n------------ N. of positive pairs in Event: "<<nEv<<" ----- "<< nPosPairs<<endl;
    cout<<"\n------------ N. of negative pairs in Event: "<<nEv<<" ----- "<< nNegPairs<<endl;
    
    nTotPosPairs += nPosPairs;
    nTotNegPairs += nNegPairs;

    // loop over JPSI candidates
    Int_t nBtoJpsiToEle = arrayBtoJpsiToEle->GetEntriesFast();
    nTotBtoJpsiToEle += nBtoJpsiToEle;
    cout<<"Number of B->J/Psi->e+e-: "<<nBtoJpsiToEle<<endl;

    for (Int_t iBtoJpsiToEle = 0; iBtoJpsiToEle < nBtoJpsiToEle; iBtoJpsiToEle++) {
      AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)arrayBtoJpsiToEle->UncheckedAt(iBtoJpsiToEle);
      Bool_t unsetvtx=kFALSE;
      if(!d->GetOwnPrimaryVtx()) {
        d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
        unsetvtx=kTRUE;
      }
      Int_t okBtoJPSI=0;
      if(d->SelectBtoJPSI(cuts,okBtoJPSI)) {
        hMassJPSI->Fill(d->InvMassJPSIee());
        hCtsJPSI->Fill(d->CosThetaStarJPSI());
        hd0d0JPSI->Fill(1e8*d->Prodd0d0());
        hCPtaJPSI->Fill(d->CosPointingAngle());
        hDCAJPSI->Fill(100*d->GetDCA());
      }
      if(unsetvtx) d->UnsetOwnPrimaryVtx();
    }

    Int_t nJPSItoEle = arrayBtoJpsiToEle->GetEntriesFast();
    nTotBtoJpsiToEle += nJPSItoEle;
    cout<<"Number of JPSI->ee: "<<nJPSItoEle<<endl;

  }

  //nTotPosPairsAvg = nTotPosPairs/nEvents;
  //nTotNegPairsAvg = nTotNegPairs/nEvents;
  //normalizationFactorAvg = 2.*TMath::Sqrt(nTotPosPairsAvg*nTotNegPairsAvg);
  Double_t normalizationFactor = 2.*TMath::Sqrt(nTotPosPairs*nTotNegPairs);

  //cout << "\n------------ Normalization factor averaged over events: -------- " << normalizationFactorAvg << endl;

  if(normalize) weightLS = normalizationFactor;

  cout <<"Applied weight to like-sign spectrum ------- " << weightLS << endl;
  printf("\n Total HF vertices: %d\n",nTotHF);
  printf("\n Total LikeSign pairs: %d\n",nTotLikeSign);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,1000);
  c1->Divide(2,3);
  c1->cd(1);
  hCPtaLikeSign->SetMarkerStyle(20);
  hCPtaLikeSign->SetMarkerSize(0.7);
  hCPtaLikeSign->SetMarkerColor(4);
  hCPtaLikeSign->Scale((1/weightLS)*hCPtaLikeSign->GetEntries());
  hCPtaLikeSign->Draw();
  hCPtaJPSI->SetMarkerStyle(20);
  hCPtaJPSI->SetMarkerSize(0.7);
  hCPtaJPSI->SetMarkerColor(2);
  hCPtaJPSI->Draw("same");
  c1->cd(2)->SetLogy();
  hd0d0LikeSign->SetMarkerStyle(20);
  hd0d0LikeSign->SetMarkerSize(0.7);
  hd0d0LikeSign->SetMarkerColor(4);
  hd0d0LikeSign->Scale((1/weightLS)*hd0d0LikeSign->GetEntries());
  hd0d0LikeSign->Draw();
  hd0d0JPSI->SetMarkerStyle(20);
  hd0d0JPSI->SetMarkerSize(0.7);
  hd0d0JPSI->SetMarkerColor(2);
  hd0d0JPSI->Draw("same");
  c1->cd(3)->SetLogy();
  hCtsLikeSign->SetMarkerStyle(20);
  hCtsLikeSign->SetMarkerSize(0.7);
  hCtsLikeSign->SetMarkerColor(4); // all LS couples
  hCtsLikeSign->Scale((1/weightLS)*hCtsLikeSign->GetEntries());
  hCtsLikeSignPos->SetMarkerStyle(20);
  hCtsLikeSignPos->SetMarkerSize(0.6);
  hCtsLikeSignPos->SetMarkerColor(3); // ++ 
  hCtsLikeSignPos->Scale((1/weightLS)*hCtsLikeSignPos->GetEntries());
  hCtsLikeSignNeg->SetMarkerStyle(21);
  hCtsLikeSignNeg->SetMarkerSize(0.6);
  hCtsLikeSignNeg->SetMarkerColor(7); // --  
  hCtsLikeSignNeg->Scale((1/weightLS)*hCtsLikeSignNeg->GetEntries());
  hCtsLikeSign->Draw();
  hCtsLikeSignPos->Draw("same");
  hCtsLikeSignNeg->Draw("same");
  hCtsJPSI->SetMarkerStyle(20);
  hCtsJPSI->SetMarkerSize(0.7);
  hCtsJPSI->SetMarkerColor(2); // JPSI is red
  hCtsJPSI->Draw("same");
  c1->cd(4);
  hDCALikeSign->SetMarkerStyle(20);
  hDCALikeSign->SetMarkerSize(0.7);
  hDCALikeSign->SetMarkerColor(4);
  hDCALikeSign->Scale((1/weightLS)*hDCALikeSign->GetEntries());
  hDCALikeSign->Draw();
  hDCAJPSI->SetMarkerStyle(20);
  hDCAJPSI->SetMarkerSize(0.7);
  hDCAJPSI->SetMarkerColor(2);
  hDCAJPSI->Draw("same");
  c1->cd(5)->SetLogy();
  hMass->SetMarkerStyle(20);
  hMass->SetMarkerSize(0.7);
  hMass->SetMarkerColor(4);
  hMass->Scale((1/weightLS)*hMass->GetEntries());
  hMass->Draw();
  hMassJPSI->SetMarkerStyle(20);
  hMassJPSI->SetMarkerSize(0.7);
  hMassJPSI->SetMarkerColor(2);
  hMassJPSI->Draw("same");
  c1->cd(6)->SetLogy();
  hMassJPSI->Draw();

  return;
}
