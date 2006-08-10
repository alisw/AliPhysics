#ifndef __CINT__
  #include <iostream.h>
  #include "AliTRDtracker.h"
  #include "AliTRDcluster.h" 
  #include "AliTRDv1.h"
  #include "AliTRDgeometry.h"    
  #include "AliTRDparameter.h"    
  #include "alles.h"  
  #include "AliTRDmcTrack.h"
  #include "AliTRDtrack.h"
  
  #include "TFile.h"
  #include "TParticle.h"
  #include "TStopwatch.h"
#endif


void AliTRDseedsAnalysis() {

  const Int_t nPrimaries = (Int_t) 3*86030/4;
  //  const Int_t nPrimaries = 500;

  Bool_t page[6]; for(Int_t i=0; i<6; i++) page[i]=kTRUE;

  const Int_t nPtSteps = 30;
  const Float_t maxPt = 3.;
  const Float_t minPt = 0.;
  
  const Int_t nEtaSteps = 40;
  const Float_t maxEta = 1.;
  const Float_t minEta = -1.;
  
  // page[0]
  TH1F *hNcl=new TH1F("hNcl","Seeds No of Clusters",255,-9.5,500.5); 
  hNcl->SetFillColor(4);
  hNcl->SetXTitle("No of clusters"); 
  hNcl->SetYTitle("counts"); 
  TH1F *hNtb=new TH1F("hNtb","Seeds No of TB with clusters",160,-9.5,150.5); 
  hNtb->SetFillColor(4);
  hNtb->SetXTitle("No of timebins with clusters"); 
  hNtb->SetYTitle("counts"); 
  TH1F *hNep=new TH1F("hNep","Seeds end point",160, -9.5, 150.5); 
  hNep->SetFillColor(4);
  hNep->SetXTitle("outermost timebin with cluster"); 
  hNep->SetYTitle("counts"); 
  TH1F *hNmg=new TH1F("hNmg","Seeds max gap",160, -9.5, 150.5); 
  hNmg->SetFillColor(4);
  hNmg->SetXTitle("max gap (tb)"); 
  hNmg->SetYTitle("counts"); 

  // page[1]
  Float_t fMin = -5.5, fMax = 134.5;
  Int_t   iChan = 140; 
  TH2F *h2ep = new TH2F("h2ep","MC vs RT (end point)",iChan,fMin,fMax,iChan,fMin,fMax); 
  h2ep->SetMarkerColor(4);
  h2ep->SetMarkerSize(2);
  TH2F *h2ntb = new TH2F("h2ntb","MC vs RT (TB with Clusters)",iChan,fMin,fMax,iChan,fMin,fMax); 
  h2ntb->SetMarkerColor(4);
  h2ntb->SetMarkerSize(2);
  TH2F *h2mg = new TH2F("h2mg","MC vs RT (max gap)",iChan/2,fMin,fMax/2,iChan/2,fMin,fMax/2); 
  h2mg->SetMarkerColor(4);
  h2mg->SetMarkerSize(2);

  // page[2]
  TH1F *hPt_all=new TH1F("hPt_all","Seeds Pt",nPtSteps,minPt,maxPt);
  hPt_all->SetLineColor(4);
  hPt_all->SetXTitle("Pt (GeV/c)"); 
  hPt_all->SetYTitle("counts"); 

  TH1F *hPt_short=new TH1F("hPt_short","Short seeds Pt",nPtSteps,minPt,maxPt);
  hPt_short->SetLineColor(8);
  TH1F *hPt_long=new TH1F("hPt_long","Long seeds Pt",nPtSteps,minPt,maxPt);
  hPt_long->SetLineColor(2);
  
  TH1F *hrtPt_short=new TH1F("hrtPt_short","RT from short seeds",nPtSteps,minPt,maxPt);
  hrtPt_short->SetLineColor(8);
  TH1F *hrtPt_long=new TH1F("hrtPt_long","RT from long seeds",nPtSteps,minPt,maxPt);
  hrtPt_long->SetLineColor(2);

  TH1F *hEta_all=new TH1F("hEta_all","Seeds Eta",nEtaSteps,minEta,maxEta);
  hEta_all->SetLineColor(4);
  hEta_all->SetXTitle("Eta"); 
  hEta_all->SetYTitle("counts"); 
  TH1F *hEta_short=new TH1F("hEta_short","Short seeds Eta",nEtaSteps,minEta,maxEta);
  hEta_short->SetLineColor(8);
  TH1F *hEta_long=new TH1F("hEta_long","Long seeds Eta",nEtaSteps,minEta,maxEta);
  hEta_long->SetLineColor(2);
  
  TH1F *hrtEta_short=new TH1F("hrtEta_short","RT from short seeds",nEtaSteps,minEta,maxEta);
  hrtEta_short->SetLineColor(8);
  TH1F *hrtEta_long=new TH1F("hrtEta_long","RT from long seeds",nEtaSteps,minEta,maxEta);
  hrtPt_long->SetLineColor(2);

  // page[3]
  TH1F *hEff_long = new TH1F("hEff_long","Efficiency vs Pt",nPtSteps,minPt,maxPt); 
  hEff_long->SetLineColor(2); hEff_long->SetLineWidth(2);  
  hEff_long->SetXTitle("Pt, GeV/c"); 
  hEff_long->SetYTitle("Efficiency"); 
  
  TH1F *hEff_short = new TH1F("hEff_short","Efficiency short",nPtSteps,minPt,maxPt); 
  hEff_short->SetLineColor(8); hEff_short->SetLineWidth(2);  
  hEff_short->SetXTitle("Pt, GeV/c"); 
  hEff_short->SetYTitle("Efficiency"); 
  
  TH1F *hEff_long_eta = new TH1F("hEff_long_eta","Efficiency vs Eta",nEtaSteps,minEta,maxEta); 
  hEff_long_eta->SetLineColor(2); hEff_long_eta->SetLineWidth(2);  
  hEff_long_eta->SetXTitle("Pt, GeV/c"); 
  hEff_long_eta->SetYTitle("Efficiency"); 
  
  TH1F *hEff_short_eta = new TH1F("hEff_short_eta","Efficiency short",nEtaSteps,minEta,maxEta); 
  hEff_short_eta->SetLineColor(8); hEff_short_eta->SetLineWidth(2);  
  hEff_short_eta->SetXTitle("Pt, GeV/c"); 
  hEff_short_eta->SetYTitle("Efficiency"); 
  
  // page[4]
  TH1F *hY=new TH1F("hY","Y resolution",50,-0.5,0.5);hY->SetLineColor(4);
  hY->SetLineWidth(2);
  hY->SetXTitle("(cm)");
  TH1F *hZ=new TH1F("hZ","Z resolution",80,-4,4);hZ->SetLineColor(4);
  hZ->SetLineWidth(2);
  hZ->SetXTitle("(cm)");
  TH1F *hPhi=new TH1F("hPhi","PHI resolution",100,-20.,20.); 
  hPhi->SetFillColor(4);
  hPhi->SetXTitle("(mrad)");
  TH1F *hLambda=new TH1F("hLambda","Lambda resolution",100,-100,100);
  hLambda->SetFillColor(17);
  hLambda->SetXTitle("(mrad)");
  TH1F *hPt=new TH1F("hPt","Relative Pt resolution",30,-10.,10.);  
  
  // page[5]
  TH1F *hY_n=new TH1F("hY_n","Normalized Y resolution",50,-0.5,0.5); hY_n->SetLineColor(4);
  hY_n->SetLineWidth(2);
  hY_n->SetXTitle("  ");
  TH1F *hZ_n=new TH1F("hZ_n","Normalized Z resolution",50,-0.5,0.5); hZ_n->SetLineColor(2);
  hZ_n->SetLineWidth(2);
  hZ_n->SetXTitle("   ");
  TH1F *hLambda_n=new TH1F("hLambda_n","Normalized Lambda resolution",50,-0.5,0.5);
  hLambda_n->SetFillColor(17);
  TH1F *hPt_n=new TH1F("hPt_n","Normalized Pt resolution",50,-0.5,0.5); 
  

  Int_t nEvent = 0;

  // Load Seeds saved as MC tracks 
  TObjArray mctarray(2000);
  TFile *mctf=TFile::Open("AliTRDtrackableSeeds.root");
  if (!mctf->IsOpen()) {cerr<<"Can't open AliTRDtrackableSeeds.root !\n"; return;}
  TTree *mctracktree=(TTree*)mctf->Get("MCtracks");
  TBranch *mctbranch=mctracktree->GetBranch("MCtracks");
  Int_t const nMCtracks = (Int_t) mctracktree->GetEntries();
  cerr<<"Found "<<nMCtracks<<" entries in the MC tracks tree"<<endl;
  for (Int_t i=0; i<nMCtracks; i++) {
    AliTRDmcTrack *ioMCtrack=new AliTRDmcTrack();
    mctbranch->SetAddress(&ioMCtrack);
    mctracktree->GetEvent(i);
    mctarray.AddLast(ioMCtrack);
  }
  mctf->Close();                 

  TFile *tf=TFile::Open("AliTRDtracks.root");

  if (!tf->IsOpen()) {cerr<<"Can't open AliTRDtracks.root !\n"; return;}
  TObjArray rtarray(2000);

  char   tname[100];
  sprintf(tname,"TRDb_%d",nEvent);     
  TTree *tracktree=(TTree*)tf->Get(tname);

  TBranch *tbranch=tracktree->GetBranch("tracks");

  Int_t nRecTracks = (Int_t) tracktree->GetEntries();
  cerr<<"Found "<<nRecTracks<<" entries in the track tree"<<endl;

  for (Int_t i=0; i<nRecTracks; i++) {
    AliTRDtrack *iotrack=new AliTRDtrack();
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(i);
    rtarray.AddLast(iotrack);
  }
  tf->Close();                 

  //  return;

  AliTRDmcTrack *mct;
  AliTRDtrack *rt;
  Int_t rtIndex[nMCtracks];
  for(Int_t j = 0; j < nMCtracks; j++) {
    mct = (AliTRDmcTrack*) mctarray.UncheckedAt(j);
    rtIndex[j] = 0;
    for (Int_t i=0; i<nRecTracks; i++) {
      rt = (AliTRDtrack*) rtarray.UncheckedAt(i);
      Int_t label = rt->GetSeedLabel();
      if(mct->GetTrackIndex() == label) rtIndex[j] = i;
    }
  } 

  // Load clusters

  TFile *geofile =TFile::Open("AliTRDclusters.root");   
  AliTRDtracker *Tracker = new AliTRDtracker(geofile);
  Tracker->SetEventNumber(nEvent);

  AliTRDgeometry *fGeom   = (AliTRDgeometry*) geofile->Get("TRDgeometry"); 
  AliTRDparameter *fPar   = (AliTRDparameter*) geofile->Get("TRDparameter"); 
  
  Char_t *alifile = "AliTRDclusters.root";
  TObjArray carray(2000);
  TObjArray *ClustersArray = &carray;
  Tracker->ReadClusters(ClustersArray,alifile);   


  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
  alifile = "galice.root";
  TFile *gafl = (TFile*) gROOT->GetListOfFiles()->FindObject(alifile);
  if (!gafl) {
    cout << "Open the ALIROOT-file " << alifile << endl;
    gafl = new TFile(alifile);
  }
  else {
    cout << alifile << " is already open" << endl;
  }
  gAlice = (AliRun*) gafl->Get("gAlice");
  if (gAlice)
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");

  AliTRDv1       *TRD = (AliTRDv1*) gAlice->GetDetector("TRD");    
  const Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) return;

  TParticle *p;
  Bool_t electrons[300000];
  for(Int_t i = 0; i < 300000; i++) electrons[i] = kFALSE;

  Bool_t mark_electrons = kFALSE;
  if(mark_electrons) {
    printf("mark electrons\n");

    for(Int_t i = nPrimaries; i < nparticles; i++) {
      p = gAlice->Particle(i);
      if(p->GetMass() > 0.01) continue;
      if(p->GetMass() < 0.00001) continue;
      electrons[i] = kTRUE;
    }
  }

  AliTRDcluster *cl = 0;
  Int_t nw, label, index, ti[3];
  Int_t mct_tbwc, rt_tbwc, mct_max_gap, rt_max_gap, mct_sector, rt_sector; 
  Int_t mct_max_tb, rt_max_tb, mct_min_tb, rt_min_tb;

  Int_t det, plane, ltb, gtb, gap, max_gap, sector;
  Double_t Pt, Px, Py, Pz;
  Double_t mcPt, mcPx, mcPy, mcPz;
  Double_t x,y,z, mcX;
  Int_t rtClusters, rtCorrect;

  Int_t nToFind_long, nFound_long, nToFind_short, nFound_short;

  printf("\n");

  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region
  Double_t dx = (Double_t) fPar->GetTimeBinSize();   

  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = 0;
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx);

  tbDrift = TMath::Min(tbDrift,maxDrift);
  tbAmp = TMath::Min(tbAmp,maxAmp);

  const Int_t nPlanes = fGeom->Nplan();
  const Int_t tbpp = Tracker->GetTimeBinsPerPlane();
  const Int_t nTB = tbpp * nPlanes;
  Int_t mask[nTB][2];

  nToFind_long = 0; nFound_long = 0;
  nToFind_short = 0; nFound_short = 0;

  for(Int_t i=0; i < nMCtracks; i++) {
    mct = (AliTRDmcTrack *) mctarray.UncheckedAt(i);
    label = mct->GetTrackIndex();

    // Waveform of the MC track
    for(nw = 0; nw < nTB; nw++) { mask[nw][0] = -1; mask[nw][1] = -1; }
    Int_t mct_ncl = mct->GetNumberOfClusters();

    for(ltb = 0; ltb < kMAX_TB; ltb++) {
      for(plane = 0; plane < nPlanes; plane++) {
	for(Int_t n = 0; n < 2; n++) {
	  index = mct->GetClusterIndex(ltb,plane,n);
	  if(index < 0) continue; 
	  cl = (AliTRDcluster *) ClustersArray->UncheckedAt(index);

	  for(nw = 0; nw < 3; nw++) ti[nw] = cl->GetLabel(nw);

	  if((ti[0] != label) && (ti[1] != label) &&  (ti[2] != label)) {
	    printf("mc wrong track label: %d, %d, %d != %d \n",
		   ti[0], ti[1], ti[2], label);
	  } 
	  det=cl->GetDetector();
	  if(fGeom->GetPlane(det) != plane) 
	    printf("cluster plane = %d != %d expected plane\n", 
		   fGeom->GetPlane(det), plane);
	  if(cl->GetLocalTimeBin() != ltb) 
	    printf("cluster ltb = %d != %d expected ltb\n", 
		   cl->GetLocalTimeBin(), ltb);
	  gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	  mask[gtb][n] = index;
	}
      }
    }
    
    for(plane = 0; plane < nPlanes; plane++) {
      for(ltb = tbDrift-1; ltb >= -tbAmp; ltb--) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb][0] > -1) break;
      }
      if(ltb >= -tbAmp) break;
    }  
    if((plane == nPlanes) && (ltb == -tbAmp-1)) {
      // printf("warning: for track %d min tb is not found and set to %d!\n",
      //	 label, nTB-1);
      mct_min_tb = nTB-1;
      // for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
    }
    else {
      mct_min_tb = Tracker->GetGlobalTimeBin(0,plane,ltb);
    }

    for(plane = nPlanes-1 ; plane>=0; plane--) {
      for(ltb = -tbAmp; ltb < tbDrift; ltb++) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb][0] > -1) break;
      }
      if(ltb < tbDrift) break;
    }  
    if((plane == -1) && (ltb == tbDrift)) {
      //      printf("warning: for track %d max tb is not found and set to 0!\n",label);
      //      for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
      mct_max_tb = 0;
      //      mcX = Tracker->GetX(0,0,tbDrift-1);
    }
    else {
      mct_max_tb = Tracker-> GetGlobalTimeBin(0,plane,ltb);
      //      mcX = Tracker->GetX(0,plane,ltb);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(mask[gtb][0]);
      det=cl->GetDetector();
      sector = fGeom->GetSector(det);  
      mct_sector = sector;
    }

    mct_tbwc = 0;
    mct_max_gap = 0;
    gap = -1;
    
    for(nw = mct_min_tb; nw < mct_max_tb+1; nw++) {
      gap++;
      if(mask[nw][0] > -1) {
	mct_tbwc++;
	if(gap > mct_max_gap) mct_max_gap = gap;
	gap = 0;
      }
    }

    //  Waveform of the reconstructed track
    if(rtIndex[i] >= 0) rt = (AliTRDtrack *) rtarray.UncheckedAt(rtIndex[i]);

    for(nw = 0; nw < nTB; nw++) { mask[nw][0] = -1; mask[nw][1] = -1; }
    Int_t rt_ncl = rt->GetNumberOfClusters();

    for(Int_t n = 0; n < rt_ncl; n++) {
      index = rt->GetClusterIndex(n);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(index);
      
      for(nw = 0; nw < 3; nw++) ti[nw] = cl->GetLabel(nw);
      
      if((ti[0] != label) && (ti[1] != label) &&  (ti[2] != label)) {
	// printf("rt wrong track label: %d, %d, %d != %d \n", ti[0], ti[1], ti[2], label);
	continue;
      } 
      
      det=cl->GetDetector();
      plane = fGeom->GetPlane(det); 
      ltb = cl->GetLocalTimeBin(); 
      gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
      mask[gtb][0] = index;
    }

    for(plane = 0; plane < nPlanes; plane++) {
      for(ltb = tbDrift-1; ltb >= -tbAmp; ltb--) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb][0] > -1) break;
      }
      if(ltb >= -tbAmp) break;
    }  
    if((plane == nPlanes) && (ltb == -tbAmp-1)) {
      // printf("warning: for track %d min tb is not found and set to %d!\n",
      //	     label, nTB-1);
      rt_min_tb = nTB-1;
      // for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
    }
    else {
      rt_min_tb = Tracker->GetGlobalTimeBin(0,plane,ltb);
    }

    for(plane = nPlanes-1 ; plane>=0; plane--) {
      for(ltb = -tbAmp; ltb < tbDrift; ltb++) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb][0] > -1) break;
      }
      if(ltb < tbDrift) break;
    }  
    if((plane == -1) && (ltb == tbDrift)) {
      // printf("warning: for track %d max tb is not found and set to 0!\n",label);
      //      for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
      rt_max_tb = 0;
      //      mcX = Tracker->GetX(0,0,tbDrift-1);
    }
    else {
      rt_max_tb = Tracker-> GetGlobalTimeBin(0,plane,ltb);
      //      mcX = Tracker->GetX(0,plane,ltb);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(mask[gtb][0]);
      det=cl->GetDetector();
      sector = fGeom->GetSector(det);  
      rt_sector = sector;
    }

    rt_tbwc = 0;
    rt_max_gap = 0;
    gap = -1;
    
    for(nw = rt_min_tb; nw < rt_max_tb+1; nw++) {
      gap++;
      if(mask[nw][0] > -1) {
	rt_tbwc++;
	if(gap > rt_max_gap) rt_max_gap = gap;
	gap = 0;
      }
    }

    // Fill the histoes

    if(page[0]) {
      hNcl->Fill((Float_t) mct_ncl);
      hNtb->Fill((Float_t) mct_tbwc);
      hNep->Fill((Float_t) mct_max_tb);
      hNmg->Fill((Float_t) mct_max_gap);
    }
    if(page[1]) {
      h2ep->Fill((Float_t) rt_max_tb, (Float_t) mct_max_tb);
      h2ntb->Fill((Float_t) rt_tbwc, (Float_t) mct_tbwc);
      h2mg->Fill((Float_t) rt_max_gap, (Float_t) mct_max_gap);
    }
    if(page[2]) {
      p = gAlice->Particle(label);
      hPt_all->Fill(p->Pt());
      hEta_all->Fill(p->Eta());
      if(mct_max_tb > 60) {
	nToFind_long++;
	hPt_long->Fill(p->Pt());
	hEta_long->Fill(p->Eta());
	if(((mct_max_tb - rt_max_tb) < 10) && 
	   (((Float_t) rt_tbwc) / ((Float_t) mct_tbwc) > 0.7)) {
	  nFound_long++;
	  hrtPt_long->Fill(p->Pt());
	  hrtEta_long->Fill(p->Eta());
	}	  
      }
      if((mct_max_tb < 60) && (mct_max_tb > 10)) {
	nToFind_short++;
	hPt_short->Fill(p->Pt());
	hEta_short->Fill(p->Eta());
	if(((mct_max_tb - rt_max_tb) < 10) && 
	   (((Float_t) rt_tbwc) / ((Float_t) mct_tbwc) > 0.7)) {
	  nFound_short++;
	  hrtPt_short->Fill(p->Pt());
	  hrtEta_short->Fill(p->Eta());
	}	  
      }
    }
    if(page[4] && page[5]) {
      if((mct_tbwc > 50) && (rt_tbwc > 50)) {
	if(rt->GetSeedLabel() != label) {
	  printf("mct vs rt index mismatch: %d != %d \n",
		 label, rt->GetSeedLabel());
	  return;
	}
	Pt = TMath::Abs(rt->GetPt()); 
	Double_t cc = TMath::Abs(rt->GetSigmaC2()); 
	mct->GetPxPyPzXYZ(mcPx,mcPy,mcPz,x,y,z,-1);
	rt->GetPxPyPz(Px,Py,Pz);      

	printf("\n\ntrack %d \n", label);
	printf("rt Px, Py, Pz: %f, %f, %f \n", Px, Py, Pz); 
	printf("mc Px, Py, Pz: %f, %f, %f \n", mcPx, mcPy, mcPz); 
    
	mcPt = TMath::Sqrt(mcPx * mcPx + mcPy * mcPy);
	if(mcPt < 0.0001) mcPt = 0.0001;
    
	Float_t lamg=TMath::ATan2(mcPz,mcPt);
	Float_t lam =TMath::ATan2(Pz,Pt);
	if(TMath::Abs(mcPt) < 0.0001) printf("attempt to divide by mcPt = %f\n",mcPt);
	else hPt_n->Fill((0.3*0.4/100*(1/Pt - 1/mcPt))/TMath::Sqrt(cc));
    
	Float_t phig=TMath::ATan2(mcPy,mcPx);
	Float_t phi =TMath::ATan2(Py,Px);    
	//	if(!(rt->PropagateTo(x))) continue;
	//    if((mcPt > Pt_min) && (mcPt < Pt_max)) {
	hLambda->Fill((lam - lamg)*1000.);
	hLambda_n->Fill((lam - lamg)/TMath::Sqrt(rt->GetSigmaTgl2()));
	if(TMath::Abs(mcPt) < 0.0001) printf("attempt to divide by mcPt = %f\n",mcPt);
	else hPt->Fill((1/Pt - 1/mcPt)/(1/mcPt)*100.); 
	hPhi->Fill((phi - phig)*1000.);
	hY->Fill(rt->GetY() - y);
	hZ->Fill(rt->GetZ() - z);
	hY_n->Fill((rt->GetY() - y)/TMath::Sqrt(rt->GetSigmaY2()));
	hZ_n->Fill((rt->GetZ() - z)/TMath::Sqrt(rt->GetSigmaZ2()));
      }
    }  
  }

  // plot pictures

  if(page[0]) {
    TCanvas *c0=new TCanvas("c0","",0,0,700,850);
    gStyle->SetOptStat(111110);
    c0->Divide(2,2);
    c0->cd(1); gPad->SetLogy(); hNcl->Draw();
    c0->cd(2); hNtb->Draw();
    c0->cd(3); hNep->Draw();
    c0->cd(4); hNmg->Draw();
    c0->Print("c0.ps","ps"); 
  }
  if(page[1]) {
    TCanvas *c1=new TCanvas("c1","",0,0,700,850);
    gStyle->SetOptStat(0);
    c1->Divide(2,2);
    c1->cd(1); h2ep->Draw();
    c1->cd(2); h2ntb->Draw();
    c1->cd(3); h2mg->Draw();
    //    c1->cd(4); hNmg->Draw();
    c1->Print("c1.ps","ps"); 
  }
  if(page[2]) {
    TCanvas *c2=new TCanvas("c2","",0,0,700,850);
    gStyle->SetOptStat(0);
    c2->Divide(2,2);
    c2->cd(1); hPt_all->Draw(); hPt_long->Draw("same"); hPt_short->Draw("same"); 
    c2->cd(2); hEta_all->Draw(); hEta_long->Draw("same"); hEta_short->Draw("same"); 
    //    c2->cd(3); h2mg->Draw();
    //    c2->cd(4); hNmg->Draw();
    c2->Print("c2.ps","ps"); 
  }
  if(page[3]) {
    TCanvas *c3=new TCanvas("c3","",0,0,700,850);
    gStyle->SetOptStat(0);
    c3->Divide(1,2);
    c3->cd(1);
    hrtPt_long->Sumw2(); hPt_long->Sumw2(); 
    hEff_long->Divide(hrtPt_long,hPt_long,1,1.,"b");
    hEff_long->SetMaximum(1.4);
    hEff_long->SetYTitle("Matching Efficiency");
    hEff_long->SetXTitle("Pt (GeV/c)");
    hEff_long->Draw();          
    hrtPt_short->Sumw2(); hPt_short->Sumw2(); 
    hEff_short->Divide(hrtPt_short,hPt_short,1,1.,"b");
    hEff_short->SetMaximum(1.4);
    hEff_short->SetYTitle("Matching Efficiency");
    hEff_short->SetXTitle("Pt (GeV/c)");
    hEff_short->Draw("same");          

    c3->cd(2);
    hrtEta_long->Sumw2(); hEta_long->Sumw2(); 
    hEff_long_eta->Divide(hrtEta_long,hEta_long,1,1.,"b");
    hEff_long_eta->SetMaximum(1.4);
    hEff_long_eta->SetYTitle("Matching Efficiency");
    hEff_long_eta->SetXTitle("Eta");
    hEff_long_eta->Draw();          
    hrtEta_short->Sumw2(); hEta_short->Sumw2(); 
    hEff_short_eta->Divide(hrtEta_short,hEta_short,1,1.,"b");
    hEff_short_eta->SetMaximum(1.4);
    hEff_short_eta->SetYTitle("Matching Efficiency");
    hEff_short_eta->SetXTitle("Eta");
    hEff_short_eta->Draw("same");          
    c3->Print("c3.ps","ps"); 
  }
  if(page[4]) {
    TCanvas *c4=new TCanvas("c4","",0,0,700,850);
    gStyle->SetOptStat(111110);
    c4->Divide(2,3);
    c4->cd(1); hY->Draw();
    c4->cd(2); hZ->Draw();
    c4->cd(3); hPhi->Draw();
    c4->cd(4); hLambda->Draw();
    c4->cd(5); hPt->Draw();
    c4->Print("c4.ps","ps"); 
  }
  if(page[5]) {
    TCanvas *c5=new TCanvas("c5","",0,0,700,850);
    gStyle->SetOptStat(111110);
    c5->Divide(2,3);
    c5->cd(1); hY_n->Draw();
    c5->cd(2); hZ_n->Draw();
    c5->cd(3); hLambda_n->Draw();
    c5->cd(4); hPt_n->Draw();
    c5->Print("c5.ps","ps"); 
  }
  printf("Efficiency (long) = %d/%d = %f\n",nFound_long,nToFind_long,
	 ((Float_t) nFound_long / (Float_t) nToFind_long));
  printf("Efficiency (shrt) = %d/%d = %f\n",nFound_short,nToFind_short,
	 ((Float_t) nFound_short / (Float_t) nToFind_short));
}


