void ITSsddanalysis (Int_t evNumber1=0,Int_t evNumber2=0)
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L anal.C   //this loads the macro in memory
//     Root > anal();     //by default process first event   
//     Root > anal(2);    //process third event
//Begin_Html
/*
<img src="gif/anal.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////
    
// Dynamically link some shared libs

  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  TString *str = new TString("galice.root");
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(str->Data());
  if (!file) file = new TFile(str->Data(),"UPDATE");
  
  // Get AliRun object from file or create it if not on file
  //   if (!gAlice) {
  gAlice = (AliRun*)file->Get("gAlice");
  if (gAlice) printf("AliRun object found on file\n");
  if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  //}
  
  TH2F *local = new TH2F("local","time vs anode difference",101,-1.01,1.01,210,-210.,210.);
  TH2F *local1 = new TH2F("local1","time vs anode difference",101,-5.01,5.01,210,-2100.,2100.);
  TH2F *dtah = new TH2F("dtah","anode difference vs drift time (hits)",210,-21000.,21000.,101,-10000.01,10000.01);
  TH2F *dtap = new TH2F("dtap","anode difference vs drift time (points)",21000,-21000.,210.,101,-10000.01,10000.01);
  TH1F *dclutim = new TH1F("dclutim","cluster time difference (dan < 1)",200,-2000.,2000.);
  TH1F *dfkctim = new TH1F("dfkctim","cluster time difference (dan > 5)",200,-2000.,2000.);
  TH2F *xy3 = new TH2F("xy3","y vs x",100,-200.02,200.02,100,-200.02,200.02);
  TH2F *xz3 = new TH2F("xz3","x vs z",100,-200.02,200.02,100,-200.02,200.02);
  TH2F *yz3 = new TH2F("yz3","y vs z",100,-200.02,200.02,100,-200.02,200.02);
  TH2F *xy4 = new TH2F("xy4","y vs x",100,-200.02,200.02,100,-200.02,200.02);
  TH2F *xz4 = new TH2F("xz4","x vs z",100,-200.02,200.02,100,-200.02,200.02);
  TH2F *yz4 = new TH2F("yz4","y vs z",100,-200.02,200.02,100,-200.02,200.02);
  TH1F *chd = new TH1F("chargediff","Charge difference (Gen-Rec)",100,-1000.,1000.);
  TH1F *chr = new TH1F("chargeratio","Charge ratio (Gen/Rec)",100,0.,0.1);
  TH2F *chp = new TH2F("chp","Point Charge vs time",28,0.,7000.,300,0.,3000.);
  TH2F *chh = new TH2F("chh","Hit Charge vs time",28,0.,7000.,80,0.,40.);
  TH2F *dadth = new TH2F("dadth","danode vs dtime (hits)",560,-14000.,14000.,601,-300.5,300.5);
  TH2F *dadt = new TH2F("dadt","danode vs dtime (points)",280,-7000.,7000.,601,-300.5,300.5);
  TH2F *aa = new TH2F("aa","anode hit vs point",204,-0.5,50.5.,204,-.5,50.5);
  TH2F *at = new TH2F("at","anode difference vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *tt = new TH2F("tt","time coord. difference (um) vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *at1 = new TH2F("at1","(1a) anode difference vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *tt1 = new TH2F("tt1","(1a) time coord. difference (um) vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *at2 = new TH2F("at2","(2a) anode difference vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *tt2 = new TH2F("tt2","(2a) time coord. difference (um) vs drift path (mm)",18,0.,36.,100,-200.,200.);
  TH2F *asigma = new TH2F("asigma","anode sigma vs drift path (mm)",18,0.,36.,300,0.,300.);
  TH2F *tsigma = new TH2F("tsigma","tau. vs drift path (mm)",18,0.,36.,200,0.,100.);
  TH2F *asigma2 = new TH2F("asigma2","2 anode sigma vs drift path (mm)",18,0.,36.,150,0.,300.);
  
  TH1F *dtrp = new TH1F("dtrp","Double track separation (microns)",40,0.,2000.);
  TH1F *dtrpall = new TH1F("dtrpall","Double track separation (mm)",100,0.,100.);
  TH1F *dtrh = new TH1F("dtrh","Double track separation (microns)",40,0.,2000.);
  TH1F *dtrhall = new TH1F("dtrhall","Double track separation (mm)",100,0.,100.);

  TH1F *p = new TH1F("p","Momentum ",100,0.,20.);
  
  TH1F *effh = new TH1F("effh","Hit Multiplicity vs drift path (mm)",18,0.,36.);
  TH1F *effp = new TH1F("effp","Point Multiplicity vs drift path (mm)",18,0.,36.);
  
  TH2F *anodes = new TH2F("nanodes","Anode Multiplicity vs drift time",28,0.,7000.,5,0.5,5.5);
  TH2F *andtsm = new TH2F("nand_ntsm","Anode Mult vs Time Mult",15,0.5,15.5,5,0.5,5.5);
  TH2F *tsampl = new TH2F("nsampls","Sample Multiplicity vs drift time",28,0.,7000.,15,0.5,15.5);
  TH2F *ntotal = new TH2F("ntotal","Cluster Multiplicity vs drift time",28,0.,7000.,60,0.5,60.5);
  TH2F *clmult = new TH2F("clmult","Anode Multiplicity vs Total Multiplicity",50,0.5,50.5,5,0.5,5.5);
  TH2F *amplit1 = new TH2F("amplit1","Point Amplitude vs drift path (mm)",28,0.,7000.,64,0.5,1024.5);
  TH2F *amplit = new TH2F("amplit","Point Amplitude vs drift path (mm)",28,0.,7000.,60,0.5,600.5);
  TH1F *hitpnt = new TH1F("hitpnt","Hit-Point Multiplicity",21,-10.5,10.5);
  
  TH1F *nmatch = new TH1F("nmatch","Number of matched points",5,-0.5.,4.5);
  TH1F *rec_vs_time = new TH1F("rec_vs_time","Point Rec. vs drift path",36,0.,36.);
  TH1F *hit_vs_time = new TH1F("hit_vs_time","Hit vs drift path",36,0.,36.);
  TH1F *rec_vs_time3 = new TH1F("rec_vs_time3","Point Rec. vs drift path",36,0.,36.);
  TH1F *hit_vs_time3 = new TH1F("hit_vs_time3","Hit vs drift path",36,0.,36.);
  TH1F *rec_vs_time4 = new TH1F("rec_vs_time4","Point Rec. vs drift path",36,0.,36.);
  TH1F *hit_vs_time4 = new TH1F("hit_vs_time4","Hit vs drift path",36,0.,36.);
  TH1F *rec_vs_time1 = new TH1F("rec_vs_time1","Point Rec. vs drift path",36,0.,36.);
  TH1F *hit_vs_time1 = new TH1F("hit_vs_time1","Hit vs drift path",36,0.,36.);
  TH1F *fake_vs_time = new TH1F("fake_vs_time","fake points vs drift path",36,0.,36.);
  
  TH1F *noihist = new TH1F("noisehist","noise",80,10.,30.);
  
  TH1F *occupancy3 = new TH1F("occupancy3","Occupancy vs Detector Number, Layer 3",20,0.5,20.5);
  TH1F *occupancy4 = new TH1F("occupancy4","Occupancy vs Detector Number, Layer 4",20,0.5,20.5);
  
  TH2F *pntmap3 = new TH2F("pntmap3","Point map Layer 3",20,0.5,20.5,10,0.5,10.5);
  TH2F *hitmap3 = new TH2F("hitmap3","Hit map Layer 3",20,0.5,20.5,10,0.5,10.5);
  TH2F *map3 = new TH2F("map3","Hit/Point map Layer 3",20,0.5,20.5,10,0.5,10.5);
  TH2F *pntmap4 = new TH2F("pntmap4","Point map Layer 4",30,0.5,30.5,10,0.5,10.5);
  TH2F *hitmap4 = new TH2F("hitmap4","Hit map Layer 4",30,0.5,30.5,10,0.5,10.5);
  TH2F *map4 = new TH2F("map4","Hit/Point map Layer 4",30,0.5,30.5,10,0.5,10.5);
  TH2F *xz = new TH2F("xz","X vs Z",50,-5,5.,50,-5.,5.);
  TH2F *and_tim = new TH2F("and_tim","Tim vs Anode",30,-100,356.,30,-8000.,8000.);
  TH2F *pand_ptim = new TH2F("pand_ptim","Tim vs Anode",30,-100,356.,30,-8000.,8000.);

  //Int_t nanodes = 256;
  //TH2F *mappa3hit[14][6][2];
  //TH2F *mappa4hit[22][8][2];
  //TH2F *mappa3pnt[14][6][2];
  //TH2F *mappa4pnt[22][8][2];
  /*
    for(Int_t i=0;i<22;i++) {
    for(Int_t j=0;j<8;j++) {
    for(Int_t k=0;k<2;k++) {
    TString *hname = new TString("hitmap_");
    TString *cname = new TString("pntmap_");
    Char_t lad[2];
    sprintf(lad,"%d",i+1);
    hname->Append(lad);
    hname->Append("_");
    cname->Append(lad);
    cname->Append("_");
    Char_t det[2];
    sprintf(det,"%d",j+1);
    hname->Append(det);
    hname->Append("_");
    cname->Append(det);
    cname->Append("_");
    Char_t wng[2];
    sprintf(wng,"%d",k+1);
    hname->Append(wng);
    cname->Append(wng);
    //mappa4hit[i][j][k] = new TH2F(hname->Data(),hname->Data(),nanodes,0.5,nanodes+0.5,256,0.5,256.5);
    //mappa4pnt[i][j][k] = new TH2F(cname->Data(),cname->Data(),nanodes,0.5,nanodes+0.5,256,0.5,256.5);
    if(i<14 && j<6) {
    mappa3hit[i][j][k] = new TH2F(hname->Data(),hname->Data(),nanodes,0.5,nanodes+0.5,256,0.5,256.5);
    mappa3pnt[i][j][k] = new TH2F(cname->Data(),cname->Data(),nanodes,0.5,nanodes+0.5,256,0.5,256.5);
    }
    }
    }
    }
  */
 
  AliITS *ITS  = (AliITS*) gAlice->GetModule("ITS");
  if (!ITS) { cout << "no ITS" << endl; return; }
    
  Int_t nparticles = gAlice->GetEvent(0);

  Int_t cp[8]={0,0,0,0,0,0,0,0};
  
  AliITSDetType *iDetType=ITS->DetType(1);
  
  AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
  if (!res1) {
    res1=new AliITSresponseSDD();
    ITS->SetResponseModel(1,res1);
  }
  res1->SetZeroSupp("2D");  // 1D
  res1->SetParamOptions("same","same");
  //res1->SetFilenames(" ","$(ALICE_ROOT)/ITS/base.dat","$(ALICE_ROOT)/ITS/2D.dat ");
  res1->SetCompressParam(cp);
  res1->SetDriftSpeed(7.3);
  Float_t vdrift = res1->DriftSpeed();
  
  AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
  AliITSgeom *aliitsgeo = ITS->GetITSgeom();
  
  Int_t cp[8]={0,0,0,0,0,0,0,0};
    
  Int_t dum = 0;
  Float_t apitch = seg1->Dpz(dum);
  Float_t tstep = seg1->Dpx(dum);
  Float_t maxand = seg1->Npz()/2.;
  //cout << "anodes: " << maxand << ", tstep: " << tstep << endl;
  
  Float_t n,b;
  res1->GetNoiseParam(n,b);
  printf("SDD: noise baseline %f %f zs option %s data type %s\n",n,b,res1->ZeroSuppOption(),res1->DataType());
  printf("SDD: DriftSpeed %f TopValue %f\n",res1->DriftSpeed(),res1->DynamicRange());
  Float_t dif0,dif1;
  res1->DiffCoeff(dif0,dif1);
  printf("SDD: dif0 %f dif1 %f\n",dif0,dif1);
  
  AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
  ITS->SetSimulationModel(1,sim1);

  //
  //   Loop over events
  //
  
  Int_t Nh=0;
  Int_t Nh1=0;
  for (int nev=0; nev<= evNumber2; nev++) {
    if(nev>0) {
      nparticles = gAlice->GetEvent(nev);
      gAlice->SetEvent(nev);
    }
    cout << "nparticles  " <<nparticles<<endl;
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return;

    // Reset Pointers
    AliITShit *itsHit;
    AliITSRecPoint *itsPnt = 0;
    AliITSRawClusterSDD *itsClu = 0;

    // Reset Event Counters
    Int_t nGoodTotalHits = 0;
    Int_t nGoodTotalPnts = 0;
        
    // Get Hit, Cluster & Recpoints Tree Pointers
    
    TTree *TH = gAlice->TreeH();
    Int_t nenthit=TH->GetEntries();
    printf("Found %d entries in the Hit tree (must be one per track per event!)\n",nenthit);
    
    ITS->GetTreeC(nev);
    TTree *TC=ITS->TreeC();
    Int_t nentclu=TC->GetEntries();
    printf("Found %d entries in the Hit tree (must be one per module per event!)\n",nentclu);
    
    TTree *TR = gAlice->TreeR();
    Int_t nentrec=TR->GetEntries();
    printf("Found %d entries in the Rec tree (must be one per module per event!)\n",nentrec);
    
    // Get Pointers to Clusters & Recpoints TClonesArrays
    
    Int_t iDet = 1;  // 1 = SDD

    TClonesArray *ITSclu  = ITS->ClustersAddress(iDet);
    TClonesArray *ITSrec  = ITS->RecPoints(); 
    
    // check recpoints
    
    Int_t nbytes = 0;
    Int_t totpoints = 0;
    Int_t totclust = 0;
    
    // check hits
    
    Int_t nmodules=0;    
    ITS->InitModules(-1,nmodules); 
    ITS->FillModules(nev,0,nmodules,"","");
    
    TObjArray *fITSmodules = ITS->GetModules();
    
    Int_t first = aliitsgeo->GetStartDet(iDet);
    Int_t last = aliitsgeo->GetLastDet(iDet);
    printf("det type %d first, last %d %d \n",iDet,first,last);
    
    for (Int_t mod=0; mod<last-first+1; mod++) {
      cout << "Module: " << mod+1 << endl;
      TTree *TR = gAlice->TreeR();
      Int_t nentrec=TR->GetEntries();
      TClonesArray *ITSrec  = ITS->RecPoints(); 
      
      ITS->ResetClusters();
      TC->GetEvent(mod);
      ITS->ResetRecPoints();
      nbytes += TR->GetEvent(mod);
      
      Int_t nrecp = ITSrec->GetEntries();
      totpoints += nrecp;
      //if (nrecp) printf("Found %d rec points for module %d\n",nrecp,mod);
      if (!nrecp) continue;
      
      Int_t nrecc = ITSclu->GetEntries();
      totclust += nrecc;
      //if (nrecc) printf("Found %d clusters for module %d\n",nrecc,mod);
      
      Int_t nrecp = ITSrec->GetEntries();
      Int_t startSDD = aliitsgeo->GetStartSDD();
      Int_t *flagP = new Int_t [nrecp];     
      memset( flagP, 0, sizeof(Int_t)*nrecp ); 

      //printf("point loop\n");

      Int_t nGoodPnts = 0;
      for (Int_t pnt=0;pnt<nrecp;pnt++) {
	itsPnt  = (AliITSRecPoint*)ITSrec->At(pnt);
	if(!itsPnt) continue;
	itsClu  = (AliITSRawClusterSDD*)ITSclu->At(pnt);
	if(!itsClu) continue;
        //itsClu->PrintInfo();	
	nGoodPnts++;
	nGoodTotalPnts++;
	
	Int_t pntlayer;
	Int_t pntladder;
	Int_t pntdetector;
	aliitsgeo->GetModuleId(mod+first,pntlayer,pntladder,pntdetector);
	Int_t pntmult = itsClu->Samples();
	Int_t pntands = itsClu->Anodes();
	Float_t pnttime = itsClu->T();
	Float_t pntanod = itsClu->A();
	Float_t pntchrg = itsClu->Q();
	Float_t pntampl = itsClu->PeakAmpl();
	Float_t pntpath = pnttime*vdrift/1000.;

	Float_t wy = 0.;
	if(itsClu->Anodes() != 0.) {
	  wy = pntmult/((Float_t) pntands);
	}
	clmult->Fill((Float_t)pntmult,(Float_t) pntands);
	ntotal->Fill(pnttime,(Float_t)pntmult);
	tsampl->Fill(pnttime,wy);
	amplit->Fill(pnttime,pntampl);
	amplit1->Fill(pnttime,pntampl);
	
	//  Detector  occupancy
	
	if(pntlayer == 3) {
	  occupancy3->Fill((Float_t) pntdetector,(Float_t) pntmult);
	  //mappa3pnt[pntladder-1][pntdetector-1][0]->Fill(pntanod,pnttime);
	}
	if(pntlayer == 4) {
	  occupancy4->Fill((Float_t) pntdetector,(Float_t) pntmult);
	  //mappa4pnt[pntladder-1][pntdetector-1][0]->Fill(pntanod,pnttime);
	}
	
	//  Point  Efficiency  vs  time.
	
	effp->Fill(pntpath);
	anodes->Fill(pnttime,pntands);
	andtsm->Fill(wy,pntands);
	chp->Fill(pnttime,pntchrg);
      }
      
      //printf("hit loop\n");

      Float_t sddLength = seg1->Dx();
      Float_t sddWidth = seg1->Dz();
      Float_t driftSpeed=res1->DriftSpeed();    
      
      Int_t nGoodHits = 0;

      AliITSmodule *Mod = (AliITSmodule *)fITSmodules->At(mod+first);
      Int_t nhits = Mod->GetNhits();
      for (Int_t hit=0;hit<nhits;hit++) {
	itsHit   = (AliITShit*)Mod->GetHit(hit);
	
	Float_t avx = 0.;
	Float_t avy = 0.;
	Float_t avz = 0.;
	Int_t ifl = 0;
	Float_t DepEnergy = 100000.*itsHit->GetIonization();
	AliITShit *itsHit1 = 0;
	if(DepEnergy == 0.) { 
	  hit++;
	  if(hit == nhits) break;
	  itsHit1 = (AliITShit*) Mod->GetHit(hit);
	  avx = itsHit1->GetXG();
	  avy = itsHit1->GetYG();
	  avz = itsHit1->GetZG();
	  ifl = 1;
	}
	avx += itsHit->GetXG();
	avy += itsHit->GetYG();
	avz += itsHit->GetZG();
	if(DepEnergy == 0.) { 
	  avx /= 2.;
	  avy /= 2.;
	  avz /= 2.;
	}
	if(ifl == 0) continue;
	
	Float_t px; Float_t py; Float_t pz;
	itsHit->GetMomentumG(px,py,pz);
	Float_t ptot = TMath::Sqrt(px*px+py*py+pz*pz);
	p->Fill(ptot*100);
	if(ptot < 0.05) continue;
	
	Int_t Layer = itsHit->GetLayer();
	Int_t Ladder = itsHit->GetLadder();
	Int_t Det = itsHit->GetDetector();
	
	Float_t And;
	Float_t Tim;
	Float_t x = itsHit->GetXL();
	Float_t z = itsHit->GetZL();
	xz->Fill(z,x);
	seg1->GetPadTxz(x,z);
	And = z;
	Tim = x*tstep;
	and_tim->Fill(And,Tim);

	Float_t And1;
	Float_t Tim1;
	Float_t x1;
	Float_t z1;
	if(itsHit1) {
	  x1 = itsHit1->GetXL();
	  z1 = itsHit1->GetZL();
	  xz->Fill(z1,x1);
	  seg1->GetPadTxz(x1,z1);
	  And1 = z1;
	  Tim1 = x1*tstep;
	  and_tim->Fill(And1,Tim1);
	}
	Float_t DepEnergy = 100000.*itsHit->GetIonization();
	if(DepEnergy == 0.) DepEnergy = 100000.*itsHit1->GetIonization(); 
	if(DepEnergy < 5.) continue;	

	if(itsHit1) {
	  Tim += Tim1;
	  Tim /= 2.;
	  And += And1;
	  And /= 2.;
	}
	if(And < 0. || And > maxand) { cout << "And: " << And << endl; continue; }
	Float_t path = TMath::Abs(Tim)*vdrift/1000.;
	hit_vs_time->Fill(path);
	if(Layer==3) hit_vs_time3->Fill(path);
	if(Layer==4) hit_vs_time4->Fill(path);
	
	nGoodHits++;
	nGoodTotalHits++;
	
	//if(Layer == 3) mappa3hit[Ladder-1][Det-1][0]->Fill(And,Tim);
	//if(Layer == 4) mappa4hit[Ladder-1][Det-1][0]->Fill(And,Tim);
	
	effh->Fill(path);
	Float_t ww = DepEnergy;
	chh->Fill(TMath::Abs(Tim),ww);
	Int_t inmatches = 0;
	
	Float_t diffmin = 100000.;
	Int_t pntmin = -1;

	//printf("point loop\n");

	for (Int_t pnt=0;pnt<nrecp;pnt++) {
	  itsPnt  = (AliITSRecPoint*)ITSrec->At(pnt);
	  if(!itsPnt) continue;
	  itsClu  = (AliITSRawClusterSDD*)ITSclu->At(pnt);
	  if(!itsClu) continue;
	  
	  Int_t LayerP;
	  Int_t LadderP;
	  Int_t DetP;
	  aliitsgeo->GetModuleId(mod+first,LayerP,LadderP,DetP);
	  Int_t LayerH =  itsHit->GetLayer();
	  Int_t LadderH =  itsHit->GetLadder();
	  Int_t DetH =  itsHit->GetDetector();
	  if(LayerH != LayerP) continue;
	  if(LadderH != LadderP) continue;
	  if(DetH != DetP) continue;

	  Float_t Pand = (Float_t) itsClu->A();
	  if(Pand < 0 || Pand > maxand) { cout << "Pand: " << Pand << endl; continue; }
	  Float_t Ptim = (Float_t) itsClu->T();
	  Float_t Pwng = (Float_t) itsClu->W();
	  if(Pwng == 1) Ptim *= -1.;
          pand_ptim->Fill(Pand,Ptim);
	  Float_t adiff = And-Pand;
	  Float_t tdiff = Tim-Ptim;  
	  if(And < 0) {
	    printf("tim %f\n",Tim);
	    printf("and %f\n",And);
	  }
	  if(Pwng == 1) tdiff *=-1.;
	  local1->Fill(adiff,tdiff);
	  
	  if(TMath::Abs(adiff) >= 1) continue;
	  if(TMath::Abs(tdiff) >= 100) continue;

	  Float_t apdiff = adiff*apitch;
	  Float_t tpdiff = tdiff*vdrift;
	 
          Float_t diff = TMath::Sqrt( apdiff*apdiff+tpdiff*tpdiff ); 
	  if( diff < diffmin ){
	    diffmin = diff;
	    pntmin = pnt;
	  }
	  	  
          if(TMath::Abs(adiff) < 1. && TMath::Abs(tdiff) < 100.) {
            inmatches++;
	  }
	  if( pntmin > -1 ) {
	    if( flagP[pntmin] == 1) continue;
	    flagP[pntmin] = 1;
	    itsClu  = (AliITSRawClusterSDD*)ITSclu->At( pntmin );
	    Float_t Pand = (Float_t) itsClu->A();
	    Float_t Ptim = (Float_t) itsClu->T();
	    Float_t Pwng = (Float_t) itsClu->W();
	    Float_t sigma = itsClu->Asigma();
	    Float_t tau = itsClu->Tsigma();
	    Int_t pntands = itsClu->Anodes();
	    if(Pwng == 1) Ptim *= -1.;
	    Float_t adiff = And-Pand;
	    Float_t tdiff = Tim-Ptim;
	    if(Pwng == 1) tdiff *=-1.;
	    local->Fill(adiff,tdiff);
	    
	    Float_t dpath = Ptim*vdrift/1000.;
	    Float_t dpathh = Tim*vdrift/1000.;
	    Float_t adpath = TMath::Abs(dpath);
	    Float_t adpathh = TMath::Abs(dpathh);
	    Float_t apdiff = adiff*apitch;
	    Float_t tpdiff = tdiff*vdrift;
	    aa->Fill(Pand,And);
	    
	    Int_t pntands = itsClu->Anodes();
	    if(pntands == 1) {
	      at1->Fill(adpath,apdiff);
	      tt1->Fill(adpath,tpdiff);
	    } 
	    if(pntands == 2) {
	      at2->Fill(adpath,apdiff);
	      tt2->Fill(adpath,tpdiff);
	    } 
	    
	    at->Fill(adpathh,apdiff);
	    tt->Fill(adpathh,tpdiff);
	    asigma->Fill(adpathh,sigma);
	    tsigma->Fill(adpathh,tau);
	    if(pntands == 2) asigma2->Fill(adpathh,sigma);
	    
	    Float_t *lP = new Float_t[3];
	    lP[0] = itsPnt->GetX();
	    lP[1] = 0.;
	    lP[2] = itsPnt->GetZ();
	    Float_t *gP = new Float_t[3];
	    aliitsgeo->LtoG(LayerH,LadderH,DetH,lP,gP);
	    Float_t dx = avx - gP[0];
	    Float_t dy = avy - gP[1];
	    Float_t dz = avz - gP[2];
	    delete lP;
	    delete gP;
	    
	    Float_t pntchrg = itsClu->Q();
	    Float_t dq = DepEnergy/0.122 - pntchrg;
	    Float_t rq = 0;
	    if(pntchrg != 0) rq = DepEnergy/0.122/((Float_t) pntchrg);
	    if(LayerH == 3) {
	      xy3->Fill(dx,dy);
	      xz3->Fill(dz,dx);
	      yz3->Fill(dz,dy);
	    } else if(LayerH == 4) {
	      xy4->Fill(dx,dy);
	      xz4->Fill(dz,dx);
	      yz4->Fill(dz,dy);
	    }
	    chd->Fill(dq);
	    if(rq != 0.) chr->Fill(rq); 
	    
	    rec_vs_time->Fill(adpathh);
	    if(Layer==3) rec_vs_time3->Fill(adpathh);
	    if(Layer==4) rec_vs_time4->Fill(adpathh);
	  }
	}
	nmatch->Fill(inmatches);
      }  // loop hits
      
      if(nGoodHits != nGoodPnts) {
	printf("module: %d",mod+1);
	printf(", nGoodHits: %d",nGoodHits);
	printf(", nGoodPnts: %d\n",nGoodPnts);
      }
      Float_t nHP = (Float_t) nGoodHits-nGoodPnts;
      hitpnt->Fill(nHP);
      
      Int_t www = 0.;
      if(nGoodHits != 0) www = nGoodPnts/((Float_t) nGoodHits);
      if(Layer == 3) {
	pntmap3->Fill(Ladder,Det,(Float_t) nGoodPnts);
	hitmap3->Fill(Ladder,Det,(Float_t) nGoodHits);
	map3->Fill(Ladder,Det,www);
      }
      if(Layer == 4) {
	pntmap4->Fill(Ladder,Det,(Float_t) nGoodPnts);
	hitmap4->Fill(Ladder,Det,(Float_t) nGoodHits);
	map4->Fill(Ladder,Det,www);
      }
      
      //printf("double hit loop\n");
      
      Stat_t wh = 1.;
      if(nGoodHits > 1) 
	wh /= (((Float_t) nGoodHits)*((Float_t) nGoodHits)-1)/2.;
      else
	wh = 0.;
      
      Int_t *flag = new Int_t[nhits];
      Int_t nGoodHitsOK = 0;
      for (Int_t hit=0;hit<nhits;hit++) {
	flag[hit] = 0;
	itsHit   = (AliITShit*)Mod->GetHit(hit);
	Float_t avx = 0.;
	Float_t avy = 0.;
	Float_t avz = 0.;
	Int_t ifl = 0;
	Float_t DepEnergy = 100000.*itsHit->GetIonization();
	AliITShit *itsHit1 = 0;
	if(DepEnergy == 0.) { 
	  hit++;
	  flag[hit] = 0;
	  if(hit == nhits) break;
	  itsHit1 = (AliITShit*) Mod->GetHit(hit);
	  avx = itsHit1->GetXG();
	  avy = itsHit1->GetYG();
	  avz = itsHit1->GetZG();
	  ifl = 1;
	}
	avx += itsHit->GetXG();
	avy += itsHit->GetYG();
	avz += itsHit->GetZG();
	if(DepEnergy == 0.) { 
	  avx /= 2.;
	  avy /= 2.;
	  avz /= 2.;
	}
	if(DepEnergy < 5. && DepEnergy > 0.) continue;
	if(ifl == 0) continue;

	Float_t px; Float_t py; Float_t pz;
	itsHit->GetMomentumG(px,py,pz);
	Float_t ptot = TMath::Sqrt(px*px+py*py+pz*pz);
	if(ptot < 0.05) continue;
	
	for (Int_t hit1=hit+1;hit1<nhits;hit1++) {
	  itsHit2   = (AliITShit*)Mod->GetHit(hit1);
	  
	  Float_t avx2 = 0.;
	  Float_t avy2 = 0.;
	  Float_t avz2 = 0.;
	  Int_t ifl2 = 0;
	  Float_t DepEnergy2 = 100000.*itsHit2->GetIonization();
	  AliITShit *itsHit3 = 0;
	  if(DepEnergy2 == 0.) { 
	    hit1++;
	    itsHit3 = (AliITShit*) Mod->GetHit(hit1);
	    avx2 = itsHit3->GetXG();
	    avy2 = itsHit3->GetYG();
	    avz2 = itsHit3->GetZG();
	    ifl2 = 1;
	  }
	  avx2 += itsHit2->GetXG();
	  avy2 += itsHit2->GetYG();
	  avz2 += itsHit2->GetZG();
	  if(DepEnergy2 == 0.) { 
	    avx2 /= 2.;
	    avy2 /= 2.;
	    avz2 /= 2.;
	  }
	  if(DepEnergy2 < 5. && DepEnergy2 > 0.) continue;
	  if(itsHit->GetLayer() != itsHit2->GetLayer()) continue;
	  if(itsHit->GetLadder() != itsHit2->GetLadder()) continue;
	  if(itsHit->GetDetector() != itsHit2->GetDetector()) continue;
	  if(ifl2 == 0) continue;
	  
	  Float_t px1; Float_t py1; Float_t pz1;
	  itsHit2->GetMomentumG(px1,py1,pz1);
	  Float_t ptot1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1);
	  if(ptot1 < 0.05) continue;
	  
	  Float_t And;
	  Float_t Tim;
	  Float_t x = itsHit->GetXL();
	  Float_t z = itsHit->GetZL();
	  seg1->GetPadTxz(x,z);
	  And = z;
	  Tim = x*tstep;
	  if(And < 0 || And > maxand) continue;
	  Float_t And2;
	  Float_t Tim2;
	  Float_t x2 = itsHit2->GetXL();
	  Float_t z2 = itsHit2->GetZL();
	  seg1->GetPadTxz(x2,z2);
	  And2 = z2;
	  Tim2 = x2*tstep;
	  if(And2 < 0 || And2 > maxand) continue;
	  Float_t da = apitch*(And-And2);
          Float_t dt = vdrift*(Tim-Tim2);
	  Float_t danh = And-And2;
	  Float_t dtmh = Tim-Tim2;
	  Float_t dist = TMath::Sqrt(da*da+dt*dt);
	  if(dt < 1000.) {
	    Float_t wx = dt*clock/(1000.*vdrift);
	    Float_t wy = da/apitch;
	    dtah->Fill(wx,wy);
	  }
	  if(dist<20.) { cout << "skip hit " << hit1 << endl; flag[hit] = 1; }
	 
	  if(ifl == 1 && ifl2 == 1) {
	    if(dist>10)dtrh->Fill(dist,wh);
	    dtrhall->Fill(dist/1000.,wh);
	    dadth->Fill(dtmh,danh);
	  }
	}  // end cluster loop
	
	Float_t path = TMath::Abs(Tim)*vdrift/1000.;
	if(flag[hit] == 0) { nGoodHitsOK++; hit_vs_time1->Fill(path); }
      } // end hit loop
      delete [] flag;
      printf("nGoodHits: %d",nGoodHits);
      printf(", nGoodHitsOK: %d\n",nGoodHitsOK);

      //printf("cluster loop\n");
 
      AliITSRawClusterSDD *itsCluFake = 0;
      Int_t nGoodPntsOK = 0;
      for( int ip=0; ip<nrecp; ip++) {
	itsClu  = (AliITSRawClusterSDD*)ITSclu->At(ip);
	if(!itsClu) continue;
	Float_t Ptim = (Float_t) itsClu->T();
	Float_t dpath = Ptim*vdrift/1000.;
	Float_t adpath = TMath::Abs(dpath);
	if( flagP[ip] == 1) { nGoodPntsOK++; rec_vs_time1->Fill(dpath); }
	else {
	  cout << "ip: " << ip << endl;
	  itsCluFake  = (AliITSRawClusterSDD*)ITSclu->At( ip );
	  if(!itsCluFake) continue;
	  cout << "pointer: " << itsCluFake << endl;
	  itsCluFake->PrintInfo();
	  Float_t Ptim = (Float_t) itsCluFake->T();
	  Float_t dpath = Ptim*vdrift/1000.;
	  fake_vs_time->Fill(dpath);  
	}
      }
      
      Stat_t wp = 1.;
      if(nGoodPntsOK > 1) 
	wp /= (((Float_t) nGoodPntsOK)*((Float_t) nGoodPntsOK)-1)/2.;
      else
	wp = 0.;
      
      //printf("double cluster loop\n");
      for (Int_t pnt=0;pnt<nrecp;pnt++) {
        if( flagP[pnt] == 0) continue;
	itsPnt  = (AliITSRecPoint*)ITSrec->At(pnt);
	if(!itsPnt) continue;
	itsClu  = (AliITSRawClusterSDD*)ITSclu->At(pnt);
	if(!itsClu) continue;
	AliITSRecPoint *itsPnt1;
	AliITSRawClusterSDD *itsClu1;
	for (Int_t pnt1=pnt+1;pnt1<nrecp;pnt1++) {
	  if( flagP[pnt1] == 0) continue;
	  itsPnt1 = (AliITSRecPoint*)ITSrec->At(pnt1);
	  if(!itsPnt1) continue;
	  itsClu1 = (AliITSRawClusterSDD*)ITSclu->At(pnt1);
	  if(!itsClu1) continue;
	  
	  Float_t dan = itsClu->A()-itsClu1->A();
	  Float_t Pwng = (Float_t) itsClu->W();
	  Float_t dt1 = itsClu->T();
	  if(Pwng == 1) dt1 *= -1.;	  
	  Float_t dt2 = itsClu1->T();
	  Float_t Pwng = (Float_t) itsClu1->W();
	  if(Pwng == 1) dt2 *= -1.;	  
	  Float_t dtm = dt1-dt2;
	  dadt->Fill(dtm,dan);
	  Float_t dap = apitch*(itsClu->A()-itsClu1->A());
          Float_t dtp = vdrift*(dt1-dt2);
	  Float_t distp = TMath::Sqrt(dap*dap+dtp*dtp);
	  if(TMath::Abs(dan) < 1.) dclutim->Fill(dtp);
	  if(TMath::Abs(dan) > 10.) dfkctim->Fill(dtp); 
	  if(dtp < 1000.) {
	    Float_t wx = dtp*clock;
	    Float_t wy = dap/apitch;
	    dtap->Fill(wx,wy);
	  }
	  dtrp->Fill(distp,wp);
	  dtrpall->Fill(distp/1000.,wp);
	} // end loop points
      } // end loop points
      
      delete [] flagP;
    } // end loop modules

    ITS->ClearModules();
    gAlice->CleanDetectors();
    
  } // end loop events

  cout << "open output file" << endl;
  TFile fhistos("SDD_histos_test.root","RECREATE");
  local->Write();
  local1->Write();
  
  aa->Write();
  at->Write();
  tt->Write();
  at1->Write();
  tt1->Write();
  at2->Write();
  tt2->Write();
  asigma->Write();
  tsigma->Write();
  asigma2->Write();
  dadt->Write();
  dadth->Write();
  dclutim->Write();
  dfkctim->Write();
  xy3->Write();
  xz3->Write();
  yz3->Write();
  xy4->Write();
  xz4->Write();
  yz4->Write();
  chd->Write();
  chr->Write();
  chh->Write();
  chp->Write();
  dtrp->Write();
  dtrpall->Write();
  dtah->Write();
  dtap->Write();
  dtrh->Write();
  dtrhall->Write();
  effh->Write();
  effp->Write();
  rec_vs_time->Write();
  hit_vs_time->Write();
  hit_vs_time1->Write();
  rec_vs_time1->Write();
  rec_vs_time3->Write();
  hit_vs_time3->Write();
  rec_vs_time4->Write();
  hit_vs_time4->Write();
  fake_vs_time->Write();
  p->Write();
  nmatch->Write();
  anodes->Write();
  tsampl->Write();
  ntotal->Write();
  clmult->Write();
  andtsm->Write();
  amplit->Write();
  amplit1->Write();
  hitpnt->Write();
  noihist->Write();
  occupancy3->Write();
  occupancy4->Write();
  map3->Write();
  hitmap3->Write();
  pntmap3->Write();
  map4->Write();
  hitmap4->Write();
  pntmap4->Write();
  /*
    for(Int_t i=0;i<22;i++) {
    for(Int_t j=0;j<8;j++) {
    for(Int_t k=0;k<2;k++) {
    //mappa4hit[i][j][k]->Write();
    //mappa4pnt[i][j][k]->Write();
    if(i<14 && j<6) {
    //mappa3hit[i][j][k]->Write();
    //mappa3pnt[i][j][k]->Write();
    }
    }
    }
    }
  */
  xz->Write();
  and_tim->Write();
  pand_ptim->Write();
  fhistos.Close();
  file->Close();
}






