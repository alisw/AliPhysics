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

void AliTRDbackTrackAnalysis() {

  //  const Int_t nPrimaries = 84210/400;
  const Int_t nPrimaries = 84210/16;


  Float_t Pt_min = 0.;
  Float_t Pt_max = 20.;

  TH1F *hp=new TH1F("hp","PHI resolution",100,-20.,20.); hp->SetFillColor(4);
  TH1F *hl=new TH1F("hl","LAMBDA resolution",100,-100,100);hl->SetFillColor(4);
  TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.);

  TH1F *hpta=new TH1F("hpta","norm. Pt resolution",50,-5.,5.);
  hpta->SetFillColor(17);
  hpta->SetXTitle("(%)");

  TH1F *hla=new TH1F("hla","norm. Lambda resolution",50,-5.,5.);
  hla->SetFillColor(17);
  TH1F *hya=new TH1F("hya","norm. Y resolution",50,-5.,5.);
  hya->SetFillColor(17);
  TH1F *hza=new TH1F("hza","norm. Z resolution",50,-5.,5.);
  hza->SetFillColor(17);

  hpt->SetFillColor(2);             

  TH1F *hy=new TH1F("hy","Y resolution",50,-0.5,0.5);hy->SetLineColor(4);
  hy->SetLineWidth(2);
  hy->SetXTitle("(cm)");

  TH1F *hz=new TH1F("hz","Z resolution",80,-4,4);hz->SetLineColor(2);
  hz->SetLineWidth(2);
  hz->SetXTitle("(cm)");

  TH1F *hx=new TH1F("hx","X(out)",150,250.,400.); hx->SetFillColor(4);

  const Int_t nPtSteps = 30;
  const Float_t maxPt = 3.;

  TH1F *hgood=new TH1F("hgood","long TRD tracks, TPC seeds",nPtSteps,0,maxPt);
  hgood->SetYTitle("Counts");
  hgood->SetXTitle("Pt (GeV/c)");
  hgood->SetLineColor(3);

  TH1F *hGoodAndSeed = new TH1F("GoodAndSeed","TPC seed and good",nPtSteps,0,maxPt);
  hGoodAndSeed->SetLineColor(2);

  TH2F *hRTvsMC = new TH2F("RTvsMC","RTvsMC",100,-4.5,95.5,100,-4.5,95.5); 
  hRTvsMC->SetMarkerColor(4);
  hRTvsMC->SetMarkerSize(2);

  TH2F *hXvsMCX = new TH2F("XvsMCX","XvsMCX",150,250,400,150,250,400); 
  hXvsMCX->SetMarkerColor(4);
  hXvsMCX->SetMarkerSize(2);

  TH2F *h2seed = new TH2F("2dGood","TPC seeds",60,0.,60,50,0.,3.); 
  h2seed->SetMarkerColor(4);
  h2seed->SetMarkerSize(2);

  TH2F *h2lost = new TH2F("2dSeedAndGood","SeedButNotGood",60,0.,60.,50,0.,3.);  
  h2lost->SetMarkerColor(2);
  h2lost->SetMarkerSize(2);


  TH1F *hseed=new TH1F("seed","TPC seeds",nPtSteps,0,maxPt);
  hseed->SetLineColor(4);


  TH1F *hfound=new TH1F("hfound","Found tracks",nPtSteps,0,maxPt);  
  hfound->SetYTitle("Counts");
  hfound->SetXTitle("Pt (GeV/c)");

  TH1F *heff=new TH1F("heff","Matching Efficiency",nPtSteps,0,maxPt); // efficiency, good tracks
  heff->SetLineColor(4); heff->SetLineWidth(2);  
  heff->SetXTitle("Pt, GeV/c"); 
  heff->SetYTitle("Efficiency"); 

  TH1F *hSeedEff = new TH1F("hSeedEff","TPC Efficiency",nPtSteps,0,maxPt); 
  hSeedEff->SetLineColor(4); hSeedEff->SetLineWidth(2);  
  hSeedEff->SetXTitle("Pt, GeV/c"); 
  hSeedEff->SetYTitle("Efficiency"); 


  TH1F *hol=new TH1F("hol","Overlap fraction",105,-2.5,102.5); 
  hol->SetLineColor(4); hol->SetLineWidth(2);  
  hol->SetXTitle("Fraction,(%)"); 
  hol->SetYTitle("Counts"); 

  TH1F *hend=new TH1F("end","missing tail",80,-10.5,69.5);

  TH1F *hFraction=new TH1F("fraction","Fraction of found clusters",110,0,1.1);
  TH1F *hCorrect=new TH1F("correct","Fraction of correct clusters",110,0,1.1);


  Int_t nEvent = 0;
  const Int_t maxIndex = nPrimaries;
  Bool_t seedLabel[maxIndex];
  Int_t mcIndex[maxIndex];
  Int_t rtIndex[maxIndex];

  for(Int_t i = 0; i < maxIndex; i++) {
    seedLabel[i] = kFALSE;
    mcIndex[i] = -1;
    rtIndex[i] = -1;
  }
  
  // mark available seeds from TPC
  Int_t nSeeds = 0;
  Int_t nPrimarySeeds = 0;
  
  printf("marking found seeds from TPC\n"); 
  TDirectory *savedir=gDirectory;  

  TFile *in=TFile::Open("AliTPCBackTracks.root");
  if (!in->IsOpen()) {
    cerr<<"can't open file AliTPCBackTracks.root  !\n"; return;
  }      
 
  char   tname[100];
  sprintf(tname,"seedsTPCtoTRD_%d",nEvent);
  TTree *seedTree=(TTree*)in->Get(tname);
  if (!seedTree) {
     cerr<<"AliTRDtracker::PropagateBack(): ";
     cerr<<"can't get a tree with seeds from TPC !\n";
  }

  AliTPCtrack *seed=new AliTPCtrack;
  seedTree->SetBranchAddress("tracks",&seed);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     Int_t lbl = seed->GetLabel();
     if(lbl < 0) { 
       printf("negative seed label %d \n",lbl); 
       continue;
     }
     if(lbl >= maxIndex) continue;
     seedLabel[lbl] = kTRUE;
     if(lbl < nPrimaries) nPrimarySeeds++;
     nSeeds++;
  }
  delete seed;
  delete seedTree;                                

  printf("Found %d seeds from primaries among overall %d seeds \n",
	 nPrimarySeeds, nSeeds);

  savedir->cd(); 
  // done with marking TPC seeds


  TFile *tf=TFile::Open("AliTRDtracks.root");

  if (!tf->IsOpen()) {cerr<<"Can't open AliTRDtracks.root !\n"; return;}
  TObjArray tarray(2000);

  sprintf(tname,"TRDb_%d",nEvent);     
  TTree *tracktree=(TTree*)tf->Get(tname);

  TBranch *tbranch=tracktree->GetBranch("tracks");

  Int_t nRecTracks = (Int_t) tracktree->GetEntries();
  cerr<<"Found "<<nRecTracks<<" entries in the track tree"<<endl;

  for (Int_t i=0; i<nRecTracks; i++) {
    AliTRDtrack *iotrack=new AliTRDtrack();
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(i);
    tarray.AddLast(iotrack);
    Int_t trackLabel = iotrack->GetLabel();

    //    printf("rt with %d clusters and label %d \n",
    //	   iotrack->GetNumberOfClusters(), trackLabel);

    if(trackLabel < 0) continue;
    if(trackLabel >= maxIndex) continue;
    rtIndex[trackLabel] = i;
  }
  tf->Close();                 

  //  return;

  // Load MC tracks 
  TFile *mctf=TFile::Open("AliTRDmcTracks.root");
  if (!mctf->IsOpen()) {cerr<<"Can't open AliTRDmcTracks.root !\n"; return;}
  TObjArray mctarray(2000);
  TTree *mctracktree=(TTree*)mctf->Get("MCtracks");
  TBranch *mctbranch=mctracktree->GetBranch("MCtracks");
  Int_t nMCtracks = (Int_t) mctracktree->GetEntries();
  cerr<<"Found "<<nMCtracks<<" entries in the MC tracks tree"<<endl;
  for (Int_t i=0; i<nMCtracks; i++) {
    AliTRDmcTrack *ioMCtrack=new AliTRDmcTrack;
    mctbranch->SetAddress(&ioMCtrack);
    mctracktree->GetEvent(i);
    mctarray.AddLast(ioMCtrack);
    Int_t mcLabel = ioMCtrack->GetTrackIndex();
    if(mcLabel < 0) {printf("negative mc label detected!\n"); continue;}
    if(mcLabel >= maxIndex) continue;
    mcIndex[mcLabel] = i;
  }
  mctf->Close();                 


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

  // Get AliRun object from file or create it if not on file
  gAlice = (AliRun*) gafl->Get("gAlice");
  if (gAlice)
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");


  // Define the objects
  AliTRDv1       *TRD = (AliTRDv1*) gAlice->GetDetector("TRD");    

  // Import the Trees for the event nEvent in the file
  const Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) return;

  TParticle *p;
  Bool_t electrons[300000] = { kFALSE };

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
  Int_t nw, label, index, ti[3], tbwc;
  Int_t det, plane, ltb, gtb, gap, max_gap, sector, mc_sector, min_tb, max_tb;
  Double_t Pt, Px, Py, Pz;
  Double_t mcPt, mcPx, mcPy, mcPz;
  Double_t x,y,z, mcX;
  Int_t rtClusters, rtCorrect;

  printf("\n");

  AliTRDmcTrack *mct = 0;
  AliTRDtrack *rt = 0;

  Double_t dxAmp = (Double_t) fGeom->CamHght();   // Amplification region
  Double_t dxDrift = (Double_t) fGeom->CdrHght(); // Drift region
  Double_t dx = (Double_t) fPar->GetTimeBinSize();   

  Int_t tbAmp = fPar->GetTimeBefore();
  Int_t maxAmp = (Int_t) ((dxAmp+0.000001)/dx);
  Int_t tbDrift = fPar->GetTimeMax();
  Int_t maxDrift = (Int_t) ((dxDrift+0.000001)/dx);

  tbDrift = TMath::Min(tbDrift,maxDrift);
  tbAmp = TMath::Min(tbAmp,maxAmp);

  const Int_t nPlanes = fGeom->Nplan();
  const Int_t tbpp = Tracker->GetTimeBinsPerPlane();
  const Int_t nTB = tbpp * nPlanes;
  Int_t mask[nTB];

  for(Int_t i=0; i < maxIndex; i++) {

    rt = 0; mct = 0;
    if(rtIndex[i] >= 0) rt = (AliTRDtrack *) tarray.UncheckedAt(rtIndex[i]);
    if(mcIndex[i] >= 0) mct = (AliTRDmcTrack *) mctarray.UncheckedAt(mcIndex[i]);

    if(!mct) continue;
    label = mct->GetTrackIndex();
    
    //    if(TMath::Abs(mct->GetMass()-0.136) > 0.01) continue;

    Int_t ncl = mct->GetNumberOfClusters();

    // check how many time bins have a cluster

    for(nw = 0; nw < nTB; nw++) mask[nw] = -1;

    for(Int_t j = 0; j < ncl; j++) {
      index = mct->GetClusterIndex(j);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(index);

      for(nw = 0; nw < 3; nw++) ti[nw] = cl->GetLabel(nw);

      if((ti[0] != label) && (ti[1] != label) &&  (ti[2] != label)) {
	printf("wrong track label: %d, %d, %d != %d \n",
	       ti[0], ti[1], ti[2], label);
      } 
      det=cl->GetDetector();
      plane = fGeom->GetPlane(det); 
      ltb = cl->GetLocalTimeBin();
      gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
      mask[gtb] = index;
    }




    for(plane = 0; plane < nPlanes; plane++) {
      for(ltb = tbDrift-1; ltb >= -tbAmp; ltb--) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb] > -1) break;
      }
      if(ltb >= -tbAmp) break;
    }  
    if((plane == nPlanes) && (ltb == -tbAmp-1)) {
      printf("warning: for track %d min tb is not found and set to %d!\n",
	     label, nTB-1);
      min_tb = nTB-1;
      //      for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
    }
    else {
      min_tb = Tracker->GetGlobalTimeBin(0,plane,ltb);
    }

    for(plane = nPlanes-1 ; plane>=0; plane--) {
      for(ltb = -tbAmp; ltb < tbDrift; ltb++) {
	gtb = Tracker->GetGlobalTimeBin(0,plane,ltb);
	if(mask[gtb] > -1) break;
      }
      if(ltb < tbDrift) break;
    }  
    if((plane == -1) && (ltb == tbDrift)) {
      printf("warning: for track %d max tb is not found and set to 0!\n",label);
      //      for(Int_t tb = 0; tb<nTB; tb++) printf("gtb %d; cl index %d\n",tb,mask[tb]);
      max_tb = 0;
      mcX = Tracker->GetX(0,0,tbDrift-1);
    }
    else {
      max_tb = Tracker-> GetGlobalTimeBin(0,plane,ltb);
      mcX = Tracker->GetX(0,plane,ltb);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(mask[gtb]);
      det=cl->GetDetector();
      sector = fGeom->GetSector(det);  
      mc_sector = sector;
    }


    tbwc = 0;
    max_gap = 0;
    gap = -1;
    
    for(nw = min_tb; nw < max_tb+1; nw++) {
      gap++;
      if(mask[nw] > -1) {
	tbwc++;
	if(gap > max_gap) max_gap = gap;
	gap = 0;
      }
    }

    if(tbwc < ((Int_t) (nTB * Tracker->GetMinClustersInTrack()))) continue;
    if(gap > Tracker->GetMaxGap()) continue; 
    p = gAlice->Particle(label);
    
    printf("good track %d with min_tb %d, max_tb %d, gap %d\n",
	   label,min_tb,max_tb,gap);

    hgood->Fill(p->Pt());

    if(!rt) continue;

    if(rt->GetLabel() != label) {
      printf("mct vs rt index mismatch: %d != %d \n",
	     label, rt->GetLabel());
      return;
    }

    if(!seedLabel[i]) continue;

    hGoodAndSeed->Fill(p->Pt());

    hXvsMCX->Fill(mcX,rt->GetX());
    
    rtClusters = rt->GetNumberOfClusters();

    // find number of tb with correct cluster
    rtCorrect = 0;
    
    for(Int_t j = 0; j < rtClusters; j++) {
      index = rt->GetClusterIndex(j);
      cl = (AliTRDcluster *) ClustersArray->UncheckedAt(index);

      for(nw = 0; nw < 3; nw++) ti[nw] = cl->GetLabel(nw);

      if((ti[0] != label)&&(ti[1] != label)&&(ti[2] != label)) continue; 
      rtCorrect++;
    }
    
    Float_t  foundFraction = ((Float_t) rtCorrect) / (tbwc + 0.00001);
    Float_t  correctFraction = ((Float_t) rtCorrect) / ((Float_t) rtClusters);

    if(foundFraction > 1) printf("fraction = %f/%f \n",
				 (Float_t) rtCorrect, (Float_t) tbwc);

    
    if(foundFraction <= 1) hFraction->Fill(foundFraction);
    hCorrect->Fill(correctFraction);

    hRTvsMC->Fill((Float_t) tbwc, (Float_t) rtClusters);
    
    if((foundFraction < 0.7) || (correctFraction < 0.7)) {
      printf("not found track %d with FrctnFound %f and FrctnCorrect %f\n",
	     label, foundFraction, correctFraction);
      continue;
    }

    hfound->Fill(p->Pt());

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
    else hpta->Fill((0.3*0.4/100*(1/Pt - 1/mcPt))/TMath::Sqrt(cc));
    
    Float_t phig=TMath::ATan2(mcPy,mcPx);
    Float_t phi =TMath::ATan2(Py,Px);
    

    if(!(rt->PropagateTo(x))) continue;

    
    if((mcPt > Pt_min) && (mcPt < Pt_max)) {
      hl->Fill((lam - lamg)*1000.);
      hla->Fill((lam - lamg)/TMath::Sqrt(rt->GetSigmaTgl2()));
      if(TMath::Abs(mcPt) < 0.0001) printf("attempt to divide by mcPt = %f\n",mcPt);
      else hpt->Fill((1/Pt - 1/mcPt)/(1/mcPt)*100.); 
      hp->Fill((phi - phig)*1000.);
      hy->Fill(rt->GetY() - y);
      hz->Fill(rt->GetZ() - z);
      hya->Fill((rt->GetY() - y)/TMath::Sqrt(rt->GetSigmaY2()));
      hza->Fill((rt->GetZ() - z)/TMath::Sqrt(rt->GetSigmaZ2()));
      hx->Fill((Float_t) x);
    }  
  }


  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1); 


  /*  
  TCanvas *c1=new TCanvas("c1","",0,0,700,850);

  //  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat(111110);
  gStyle->SetOptFit(1); 
    
  TPad *p1=new TPad("p1","",0,0.3,.5,.6); p1->Draw();
  p1->cd(); p1->SetFillColor(10); p1->SetFrameFillColor(10);
  hgood->Draw(); hGoodAndSeed->Draw("same"); c1->cd();  
  
  TPad *p2=new TPad("p2","",0.5,.3,1,.6); p2->Draw();
  p2->cd(); p2->SetFillColor(10); p2->SetFrameFillColor(10);
  hFraction->SetYTitle("Counts");
  hFraction->SetXTitle("Fraction");
  hFraction->Draw(); c1->cd();      
  
  TPad *p3=new TPad("p3","",0,0,0.5,0.3); p3->Draw();
  p3->cd(); p3->SetFillColor(10); p3->SetFrameFillColor(10);
  hfound->Draw(); c1->cd();   
    
  TPad *p4=new TPad("p4","",0.5,0,1.,0.3); p4->Draw();
  p4->cd(); p4->SetFillColor(10); p4->SetFrameFillColor(10);
  hCorrect->SetYTitle("Counts");
  hCorrect->SetXTitle("Fraction");
  hCorrect->Draw(); c1->cd();

  TPad *p5=new TPad("p5","",0,0.6,1,1); p5->Draw(); p5->cd();
  p5->SetFillColor(10); p5->SetFrameFillColor(10);
  hgood->Sumw2(); hGoodAndSeed->Sumw2(); // hfake->Sumw2();
  hSeedEff->Divide(hGoodAndSeed,hgood,1,1.,"b");
  //  hf->Divide(hfake,hgood,1,1.,"b");
  hSeedEff->SetMaximum(1.4);
  hSeedEff->SetYTitle("TPC efficiency");
  hSeedEff->SetXTitle("Pt (GeV/c)");
  hSeedEff->Draw();          

  TLine *line1 = new TLine(0,1.0,maxPt,1.0); line1->SetLineStyle(4);
  line1->Draw("same");
  TLine *line2 = new TLine(0,0.9,maxPt,0.9); line2->SetLineStyle(4);
  line2->Draw("same");    
  */  

  TCanvas *c2=new TCanvas("c2","",0,0,700,850);
    
  TPad *p12=new TPad("p12","",0,0.3,.5,.6); p12->Draw();
  p12->cd(); p12->SetFillColor(42); p12->SetFrameFillColor(10);
  hp->SetFillColor(4);  hp->SetXTitle("(mrad)"); hp->Fit("gaus"); c2->cd();  
  
  TPad *p22=new TPad("p22","",0.5,.3,1,.6); p22->Draw();
  p22->cd(); p22->SetFillColor(42); p22->SetFrameFillColor(10);
  hl->SetXTitle("(mrad)"); hl->Fit("gaus"); c2->cd();      
  
  TPad *p32=new TPad("p32","",0,0,0.5,0.3); p32->Draw();
  p32->cd(); p32->SetFillColor(42); p32->SetFrameFillColor(10);
  hpt->SetXTitle("(%)"); hpt->Fit("gaus"); c2->cd();   
    
  TPad *p42=new TPad("p42","",0.5,0,1.,0.3); p42->Draw();
  p42->cd(); p42->SetFillColor(42); p42->SetFrameFillColor(10);
  hgood->Draw(); hGoodAndSeed->Draw("same"); c2->cd();  

  TPad *p52=new TPad("p52","",0,0.6,1,1); p52->Draw(); p52->cd();
  p52->SetFillColor(41); p52->SetFrameFillColor(10);
  hfound->Sumw2();
  heff->Divide(hfound,hGoodAndSeed,1,1.,"b");
  //  hf->Divide(hfake,hgood,1,1.,"b");
  heff->SetMaximum(1.4);
  heff->SetXTitle("Pt (GeV/c)");
  heff->Draw();          

  TLine *line12 = new TLine(0,1.0,maxPt,1.0); line12->SetLineStyle(4);
  line12->Draw("same");
  TLine *line22 = new TLine(0,0.9,maxPt,0.9); line22->SetLineStyle(4);
  line22->Draw("same");    

  c2->Print("matching.ps","ps"); 


  /*  
  TCanvas *cxyz = new TCanvas("cxyz","",50,50,750,900);
  cxyz->Divide(2,2);
  cxyz->cd(1); hx->Draw();
  cxyz->cd(2); hy->Draw();
  cxyz->cd(3); hz->Draw();
  cxyz->cd(4); hXvsMCX->Draw();
  */

  TCanvas *cs = new TCanvas("cs","",0,0,700,850);
  cs->Divide(2,3);

  cs->cd(1); hy->Fit("gaus");
  cs->cd(2); hz->Fit("gaus");
  cs->cd(3); hpta->Fit("gaus");
  cs->cd(4); hla->Fit("gaus");
  cs->cd(5); hya->Fit("gaus");
  cs->cd(6); hza->Fit("gaus");

  cs->Print("resolution.ps","ps"); 

  /*
  TCanvas *cvs = new TCanvas("cvs","",0,0,700,850);
  cvs->cd(); hRTvsMC->Draw();
  */

}


