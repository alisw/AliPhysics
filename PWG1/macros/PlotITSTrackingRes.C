Bool_t gboxbin=kFALSE;
Int_t  gMarkerColor=2;

void PlotITSTrackingRes(TString filename="ITS.Performance.root",
			Int_t pdgcode=211,
			Bool_t useAbsPdg=kTRUE,
			Bool_t box=kTRUE,
			Int_t minITSpoints=6,
			Bool_t nofakes=kTRUE,
			Bool_t askITSrefit=kTRUE,
			Int_t minTPCcls=1) 
{
  //
  // Plot ITS tracking resolutions from ITS.Performance.root
  // A. Dainese
  // 


  //Open File
  if(gSystem->AccessPathName(filename.Data())) {
    printf("file not found!\n");  
    return;
  }
  TFile *file= TFile::Open(filename.Data());
  cout<<"Opening file "<<filename.Data()<<endl;   
    
  TList *list = (TList*)file->Get("cOutputITS");
  TNtuple *ntTracks = (TNtuple*)list->FindObject("fNtupleESDTracks");

      
  //Getting and Addressing  NTuples
    
  Float_t pt,ptmes,eta,phi,pdg,d0True,d0TV,d0zTV,sigmad0zTV,d0All,d0Oth,sigmad0Oth,ITSflag,sigmad0TV,sigmad0All;
  
  ntTracks->SetBranchAddress("pt",&ptmes);
  ntTracks->SetBranchAddress("ptMC",&pt);
  ntTracks->SetBranchAddress("eta",&eta);
  ntTracks->SetBranchAddress("phi",&phi);
  ntTracks->SetBranchAddress("pdgMC",&pdg);
  ntTracks->SetBranchAddress("d0MC",&d0True);
  ntTracks->SetBranchAddress("d0MCv",&d0TV);
  ntTracks->SetBranchAddress("z0MCv",&d0zTV);
  ntTracks->SetBranchAddress("sigmad0MCv",&sigmad0TV);
  ntTracks->SetBranchAddress("sigmaz0MCv",&sigmad0zTV);
  ntTracks->SetBranchAddress("d0",&d0All);
  ntTracks->SetBranchAddress("sigmad0",&sigmad0All);
  ntTracks->SetBranchAddress("ITSflag",&ITSflag);

    
  //make STANDARD Pt BINNING
  gboxbin=box;
  const Int_t nPtBins=10;
  Double_t ptbinning[2*nPtBins];
  const Int_t numint=10;
  Double_t ptmin[numint],lenghtbin[numint],stepbin[numint],interval[numint],steplenght;
  Int_t nbinning[numint];
  Int_t runningbin=0;

  for(Int_t i=0;i<nPtBins;i++) nbinning[i]=1;
  
  ptbinning[0]=0.2;
  ptbinning[1]=0.35;
  ptbinning[2]=0.4;
  ptbinning[3]=0.6;
  ptbinning[4]=0.7;
  ptbinning[5]=0.85;
  ptbinning[6]=0.85;
  ptbinning[7]=1.2;
  ptbinning[8]=1.7;
  ptbinning[9]=2.3;
  ptbinning[10]=3.6;
  ptbinning[11]=4.4;
  ptbinning[12]=6.5;
  ptbinning[13]=7.5;
  ptbinning[14]=9.5;
  ptbinning[15]=10.5;
  ptbinning[16]=19.5;
  ptbinning[17]=20.5;
  ptbinning[18]=29.5;
  ptbinning[19]=30.5;
 
  interval[0]=600;
  interval[1]=300;
  interval[2]=200;
  interval[3]=150;
  interval[4]=100;
  interval[5]=100;
  interval[6]=100;
  interval[7]=70;
  interval[8]=50;
  interval[9]=50;

  
  const Int_t dipnbin=10;
  const Double_t lenghtbindip=0.1;//0.03;
  const Double_t stepbindip=0.0;
  const Double_t dipinterval=200;
  
  
  //  Double_t radice=TMath::Sqrt(2.)/2.; // |eta|<0.9
  Double_t radice=0.; // all eta
 
  Int_t ntracks,entries,ntotaltracks=0,nfaketr=0;
  Double_t bin,sinteta,rms;
  Double_t ptbin[nPtBins],dipbin[dipnbin],multptbin[nPtBins],multdipbin[dipnbin],ex[nPtBins],exdip[dipnbin];
  Double_t sigmaResolutionTV[dipnbin][nPtBins],sigmaResolutionPtTV[nPtBins],errsigmResPtTV[nPtBins],meanResolutionPtTV[nPtBins],errmeanResPtTV[nPtBins],sigmaResolutionDipTV[dipnbin],errsigmResDipTV[dipnbin];
  Double_t binpt,bindip,sigmaResolutionOV[dipnbin][nPtBins],sigmaResolutionPtOV[nPtBins],errsigmResPtOV[nPtBins],sigmaResolutionDipOV[dipnbin],errsigmResDipOV[dipnbin];
  Double_t sigmaResolutionRV[dipnbin][nPtBins],sigmaResolutionPtRV[nPtBins],errsigmResPtRV[nPtBins],sigmaResolutionDipRV[dipnbin],errsigmResDipRV[dipnbin];
  Double_t mediumsigmad0rphiPtTV[nPtBins],errmediumsigmad0rphiPtTV[nPtBins],mediumsigmad0rphiPtRV[nPtBins],errmediumsigmad0rphiPtRV[nPtBins],mediumsigmad0rphiPtOV[nPtBins],errmediumsigmad0rphiPtOV[nPtBins];
  Double_t sigmaPullPtTV[nPtBins],errsigmPullPtTV[nPtBins],sigmaPullPtRV[nPtBins],errsigmPullPtRV[nPtBins],sigmaPullPtOV[nPtBins],errsigmPullPtOV[nPtBins];
  Double_t sigmad0zResPtTV[nPtBins],errsigmd0zResPtTV[nPtBins];
  // Double_t sigmad0zResPtRV[nPtBins],errsigmd0zResPtRV[nPtBins],sigmad0zResPtOV[nPtBins],errsigmd0zResPtOV[nPtBins];
  Double_t sigmad0zPullPtTV[nPtBins],errsigmd0zPullPtTV[nPtBins];
  Double_t sigmaPtResTV[nPtBins],errsigmPtResTV[nPtBins],meanPtResTV[nPtBins],errmeanPtResTV[nPtBins],sigmaCurvResTV[nPtBins],errsigmCurvResTV[nPtBins];
  // Double_t sigmad0zPullPtRV[nPtBins],errsigmd0zPullPtRV[nPtBins],sigmad0zPullPtOV[nPtBins],errsigmd0zPullPtOV[nPtBins];
  
  Int_t ITSnpoints;
  
  Float_t hptbin[nPtBins+1],hdipbin[dipnbin+1];
  
  TF1 *gauss = new TF1("gauss","gaus",-10,10); 

  TGraphErrors *gr0TV;
  TGraphErrors *gr0OV;
  TGraphErrors *gr0RV; 
  TGraphErrors *gr1TV;
  TGraphErrors *gr1OV;
  TGraphErrors *gr1RV;

  TGraphErrors *grPullTV;
  TGraphErrors *grPullOV;
  TGraphErrors *grPullRV; 

  TGraphErrors *grsigmad0rphiTV;
  TGraphErrors *grsigmad0rphiRV;
  TGraphErrors *grsigmad0rphiOV;

  TGraphErrors *grd0zResTV;
  //TGraphErrors *grd0zRes0OV;
  //TGraphErrors *grd0zRes0RV; 
  
  TGraphErrors *grd0zPullTV;
  //TGraphErrors *grd0zPull0OV;
  //TGraphErrors *grd0zPull0RV; 
  TGraphErrors *grPtResTV;
  TGraphErrors *grPtResMeanTV;
  TGraphErrors *grCurvResTV;

  TString numerobin,numerobin2,numerobinSigma,numerobinSigma2,numerobinPull,numerobinPull2,numerobind0zRes,numerobind0zPull,numerobind0zRes2,numerobind0zPull2,numerobinCurvRes,numerobinCurvRes2,numerobinPtRes,numerobinPtRes2,numerobin3,titlegraph1,titlegraph2,str=" ";
  

  
  //DECLARING AND CONSTRUCTING HISTOGRAMS
   
  TH1F **hFitResolutionTV=new TH1F*[dipnbin*nPtBins];
  TH1F **hFitResolutionOV=new TH1F*[dipnbin*nPtBins];
  TH1F **hFitResolutionRV=new TH1F*[dipnbin*nPtBins];

  TH1F **hFitResolutionPtTV=new TH1F*[nPtBins];
  TH1F **hFitResolutionPtRV=new TH1F*[nPtBins];
  TH1F **hFitResolutionPtOV=new TH1F*[nPtBins];

  TH1F **hFitsigmad0rphiPtTV=new TH1F*[nPtBins];
  TH1F **hFitsigmad0rphiPtRV=new TH1F*[nPtBins];
  TH1F **hFitsigmad0rphiPtOV=new TH1F*[nPtBins];

  TH1F **hFitPullPtTV=new TH1F*[nPtBins];
  TH1F **hFitPullPtRV=new TH1F*[nPtBins];
  TH1F **hFitPullPtOV=new TH1F*[nPtBins];

  TH1F **hFitd0zResPtTV=new TH1F*[nPtBins];
  //  TH1F **hFitd0zResPtRV=new TH1F*[nPtBins];
  //TH1F **hFitd0zResPtOV=new TH1F*[nPtBins];

  TH1F **hFitd0zPullPtTV=new TH1F*[nPtBins];
  //TH1F **hFitd0zPullPtRV=new TH1F*[nPtBins];
  //TH1F **hFitd0zPullPtOV=new TH1F*[nPtBins];

  TH1F **hFitResolutionDipTV=new TH1F*[dipnbin];
  TH1F **hFitResolutionDipRV=new TH1F*[dipnbin];
  TH1F **hFitResolutionDipOV=new TH1F*[dipnbin];
  
  TH1F **hFitPtResTV=new TH1F*[nPtBins];
  TH1F **hFitCurvResTV=new TH1F*[nPtBins];
 
  TH2F *hFitResPtDipTV;
  TH2F *hFitResPtDipRV;
  TH2F *hFitResPtDipOV;

  
  Int_t incycle=0;
  
  for(Int_t v=0;v<dipnbin;v++) {
    incycle=0;
    for (Int_t nint=0;nint<numint;nint++) {
      for (Int_t k=0;k<nbinning[nint];k++) { 
	numerobin ="TranvDipResolution";
	numerobin+=v;
	numerobin.Append("_");
	numerobin+=(k+incycle);	
	numerobin.Append("bin");
	
	numerobin2=numerobin;
	numerobin.Append("TV");
	hFitResolutionTV[v*nPtBins+k+incycle]=new TH1F(numerobin,numerobin,3*100,-2*interval[nint],2*interval[nint]);
	//      numerobin.Replace(24,25,"RV");
	//	cout<<numerobin<<endl;
	
	numerobin=numerobin2;
	numerobin2.Append("RV");
	hFitResolutionRV[v*nPtBins+k+incycle]=new TH1F(numerobin2,numerobin2,3*100,-2*interval[nint],2*interval[nint]);
	//	cout<<numerobin2<<endl;
	
	numerobin.Append("OV");
	hFitResolutionOV[v*nPtBins+k+incycle]=new TH1F(numerobin,numerobin,3*100,-2*interval[nint],2*interval[nint]);
	//	cout<<numerobin<<endl;
	
	//	cout<<(v*nPtBins+k+incycle)<<endl;
      }
      incycle+=nbinning[nint];
    }
  }

  for(Int_t v=0;v<dipnbin;v++) {
    numerobin="DipResolution";
    numerobin+=v;
    numerobin.Append("bin");
    
    numerobin2=numerobin;
    numerobin.Append("TV");
    hFitResolutionDipTV[v]=new TH1F(numerobin,numerobin,3*100,-3*dipinterval,3*dipinterval);

    numerobin=numerobin2;
    numerobin2.Append("RV");
    hFitResolutionDipRV[v]=new TH1F(numerobin2,numerobin2,3*100,-3*dipinterval,3*dipinterval);

    numerobin.Append("OV");
    hFitResolutionDipOV[v]=new TH1F(numerobin,numerobin,3*100,-3*dipinterval,3*dipinterval);

    dipbin[v]=0;
    multdipbin[v]=0;    
  }
 
  incycle=0;
  
  for (Int_t nint=0;nint<numint;nint++) {   
    for(Int_t k=0;k<nbinning[nint];k++) { 
      numerobin = "d0PtResolution";
      numerobinSigma="d0rphiPtsigma";
      numerobinPull="d0PtPull";
      numerobind0zRes="d0zPtResolution";
      numerobind0zPull="d0zPtPull";
      numerobinPtRes="PtResolution";
      numerobinCurvRes="CurvRes";


      numerobin+=k+incycle;
      numerobinSigma+=k+incycle;
      numerobinPull+=k+incycle;
      numerobind0zRes+=k+incycle;
      numerobind0zPull+=k+incycle;
      numerobinPtRes+=k+incycle;
      numerobinCurvRes+=k+incycle;

      numerobin.Append("bin");
      numerobinSigma.Append("bin");
      numerobinPull.Append("bin");
      numerobind0zRes.Append("bin");
      numerobind0zPull.Append("bin");
      numerobinPtRes.Append("bin");
      numerobinCurvRes.Append("bin");

	
      numerobin2=numerobin;
      numerobinSigma2=numerobinSigma;
      numerobinPull2=numerobinPull;
      numerobind0zRes2=numerobind0zRes;
      numerobind0zPull2=numerobind0zPull;
      numerobinPtRes2=numerobinPtRes;
      numerobinCurvRes2=numerobinCurvRes;

      numerobin.Append("TV");
      hFitResolutionPtTV[k+incycle]=new TH1F(numerobin,numerobin,100,-2*interval[nint],2*interval[nint]);
      hFitsigmad0rphiPtTV[k+incycle]=new TH1F(numerobinSigma,numerobinSigma,3*100,0.,4*interval[nint]);
      numerobinPull.Append("TV");
      hFitPullPtTV[k+incycle]=new TH1F(numerobinPull,numerobinPull,100,-4.,4.);
      numerobind0zRes.Append("TV");
      hFitd0zResPtTV[k+incycle]=new TH1F(numerobind0zRes,numerobind0zRes,3*100,-3*interval[nint],3*interval[nint]);
      numerobind0zPull.Append("TV");
      hFitd0zPullPtTV[k+incycle]=new TH1F(numerobind0zPull,numerobind0zPull,100,-4.,4.);
      numerobinPtRes.Append("TV");
      hFitPtResTV[k+incycle]=new TH1F(numerobinPtRes,numerobinPtRes,300,-1.,1.);
      numerobinCurvRes.Append("TV");
      hFitCurvResTV[k+incycle]=new TH1F(numerobinCurvRes,numerobinCurvRes,300,-1.,1.);

      numerobin=numerobin2;
      numerobinSigma=numerobinSigma2;
      numerobin.Append("RV");
      numerobinSigma.Append("RV");
      hFitResolutionPtRV[k+incycle]=new TH1F(numerobin,numerobin,3*100,-2*interval[nint],2*interval[nint]);
      hFitsigmad0rphiPtRV[k+incycle]=new TH1F(numerobinSigma,numerobinSigma,3*100,0.,4*interval[nint]);
      numerobinPull=numerobinPull2;
      numerobinPull.Append("RV");
      hFitPullPtRV[k+incycle]=new TH1F(numerobinPull,numerobinPull,100,-4,4);

      /*      numerobind0zRes=numerobind0zRes2;
      numerobind0zRes.Append("RV");
      hFitd0zResPtRV[k+incycle]=new TH1F(numerobind0zRes,numerobind0zRes,3*100,-3*interval[nint],3*interval[nint]);
      numerobind0zPull=numerobind0zPull2;
      numerobind0zPull.Append("RV");
      hFitd0zPullPtRV[k+incycle]=new TH1F(numerobind0zPull,numerobind0zPull,100,-4,4);
      */

      numerobin=numerobin2;
      numerobin.Append("OV");
      hFitResolutionPtOV[k+incycle]=new TH1F(numerobin,numerobin,3*100,-2*interval[nint],2*interval[nint]);

      numerobinSigma=numerobinSigma2;
      numerobinSigma.Append("OV");
      hFitsigmad0rphiPtOV[k+incycle]=new TH1F(numerobinSigma,numerobinSigma,3*100,0.,4.*interval[nint]);
      numerobinPull=numerobinPull2;
      numerobinPull.Append("OV");
      hFitPullPtOV[k+incycle]=new TH1F(numerobinPull,numerobinPull,100,-4,4);

      /* numerobind0zRes=numerobind0zRes2;
      numerobind0zRes.Append("OV");
      hFitd0zResPtOV[k+incycle]=new TH1F(numerobind0zRes,numerobind0zRes,3*100,-3*interval[nint],3*interval[nint]);
      numerobind0zPull=numerobind0zPull2;
      numerobind0zPull.Append("OV");
      hFitd0zPullPtOV[k+incycle]=new TH1F(numerobind0zPull,numerobind0zPull,100,-4,4);
      */
      
      ptbin[k+incycle]=0;
      multptbin[k+incycle]=0;
      hptbin[k+incycle]=0;
      //      hptbin[(k+incycle)+1]=0;
      //hptbin[(k+incycle)]=0;
    }
    incycle+=nbinning[nint];
  }
    
  //  binpt=lenghtbin+stepbin;

  bindip=lenghtbindip+stepbindip;
  Int_t np=0;  
  Int_t kbox=0;
  if(gboxbin) kbox=1;



  // ------------ Loop on Tracks ----------
  ntracks=ntTracks->GetEntries();
  cout<<"Number of Tracks: "<<ntracks<<endl;  
  for (Int_t j=0;j<ntracks;j++) {
    if(j%5000==0) printf("Reading track %d\n",j);
    ntTracks->GetEvent(j);    
    d0True*=1.e4;
    d0TV*=1.e4;
    d0zTV*=1.e4;
    sigmad0zTV*=1.e4;
    d0All*=1.e4;
    d0Oth*1.e4;
    sigmad0Oth*=1.e4;
    sigmad0TV*=1.e4;
    sigmad0All*=1.e4;

    if(TMath::Abs(eta)>0.9) continue;
    //if(TMath::Abs(eta)<0.9 || TMath::Abs(eta)>1.4) continue;
    Float_t theta = 2.*TMath::ATan(TMath::Exp(-eta));
    //if(TMath::Abs(TMath::Abs(phi)-0.5*TMath::Pi())>0.5 || TMath::Abs(theta-0.5*TMath::Pi())>0.5) continue;
    Bool_t isFake=kFALSE;
    if(ITSflag<0) isFake=kTRUE;
    //printf("  %d\n",ITSflag);
    ITSflag=TMath::Abs(ITSflag); 
    //printf("%d\n",ITSflag);
    Int_t ITSflagorig=ITSflag;
    Int_t nTPCcls =       (Int_t)(ITSflag/100000);
    ITSflag -= nTPCcls*100000;
    Int_t nITSclspart =   (Int_t)(ITSflag/10000);
    ITSflag -= nITSclspart*10000;
    //printf("%d\n",ITSflag);
    Int_t nITSclsassign = (Int_t)(ITSflag/1000);
    ITSflag -= nITSclsassign*1000;
    //printf("%d\n",ITSflag);
    //printf("%d\n",ITSflag);
    Int_t nITSsel=ITSnCluster(ITSflag);


    if(nTPCcls<=minTPCcls) continue;

    if(useAbsPdg) pdg=TMath::Abs(pdg);
    if(pdg!=pdgcode) continue;


    
    if(TMath::Abs(nITSclsassign)<minITSpoints) continue;

    if(askITSrefit && sigmad0TV<0.) continue;   

    if(TMath::Abs(d0True)>0.5) continue;//skipping secondaries tracks if asked
    //consider only "part" (if "all" is set to 1 all particles are considered)

    if(isFake) nfaketr++;
    if(nofakes && isFake) continue;        // reject fakes
    np++;

    //hPtResolutionOV->Fill(pt-ptmes);
    //hPtResolutionRV->Fill(ptRV-ptmesRV);
    //hPzResolutionOV->Fill(pz-pzmes);
    //hPzResolutionRV->Fill(pzRV-pzmesRV);
    
    sinteta=TMath::Sin(2*TMath::ATan(TMath::Exp(-eta)));

    for(Int_t v=0;v<dipnbin;v++) {
      if (v*bindip+radice<sinteta&&sinteta<=v*bindip+radice+lenghtbindip) { //printf("Prova2 %d \n",v);
	
	for (Int_t k=0;k<nPtBins;k=k+1) {
	  
	  hFitResolutionDipTV[v]->Fill(d0TV-d0True);
	  
	  if(sigmad0All>0.) {
	    hFitResolutionDipRV[v]->Fill(d0All-d0True);
	  }
	  //if(sigmad0Oth>0.) {
	  //  hFitResolutionDipOV[v]->Fill(d0Oth-d0True);
	  //}
	  
	  //Transverse Momentum bin
	  if(ptbinning[(1+kbox)*k]<pt&&pt<=ptbinning[(1+kbox)*k+1]) { //printf("Prova3 %d \n",k/2)
	    if(ptmes!=0.&&pt!=0.) {
	      hFitPtResTV[k]->Fill(ptmes/pt-1.);
	      hFitCurvResTV[k]->Fill(1/ptmes-1/pt);
	    }
	    hFitResolutionTV[v*nPtBins+k]->Fill(d0TV-d0True);
	    hFitResolutionPtTV[k]->Fill(d0TV-d0True);
	    //to calculate resolution on d0z we should have d0zTrue but we haven't. So we ask for the track being primary and then d0zTV is equal to the residual
	    if(TMath::Abs(d0True)<0.1) {
	      hFitd0zResPtTV[k]->Fill(d0zTV);
	      if(sigmad0zTV>0)hFitd0zPullPtTV[k]->Fill(d0zTV/sigmad0zTV);
	    }
	    if(sigmad0TV>0) {
	      hFitPullPtTV[k]->Fill((d0TV-d0True)/sigmad0TV);
	      hFitsigmad0rphiPtTV[k]->Fill(sigmad0TV);
	    }
	    if(sigmad0All>=0.) {
	      hFitResolutionRV[v*nPtBins+k]->Fill(d0All-d0True);	
	      hFitResolutionPtRV[k]->Fill(d0All-d0True);
	      if(sigmad0All>0) {
		hFitPullPtRV[k]->Fill((d0All-d0True)/sigmad0All);
		hFitsigmad0rphiPtRV[k]->Fill(sigmad0All);
	      }
	    }
	    /*if(sigmad0Oth>=0.) {
	      hFitResolutionOV[v*nPtBins+k]->Fill(d0Oth-d0True);
	      hFitResolutionPtOV[k]->Fill(d0Oth-d0True);
	      if(sigmad0Oth>0) {
	      hFitPullPtOV[k]->Fill((d0Oth-d0True)/sigmad0Oth);
	      hFitsigmad0rphiPtOV[k]->Fill(sigmad0Oth);
	      }
	      }*/
	    
	    ptbin[k]+=pt;
	    multptbin[k]++;  
	    dipbin[v]+=sinteta;
	    multdipbin[v]++;
	  }
	}
      }
    }
  }

  printf("Loop on ntuple finished \n");

  /*
    for (Int_t k=0;k<nPtBins;k++) {
    for (Int_t j=0;j<dipnbin;j++) {
    hFitResolutionPtTV[k]->Add(hFitResolutionTV[j*nPtBins+k],1.);
    hFitResolutionPtRV[k]->Add(hFitResolutionRV[j*nPtBins+k],1.);
    hFitResolutionPtOV[k]->Add(hFitResolutionOV[j*nPtBins+k],1.);
    
    hFitResolutionDipTV[j]->Add(hFitResolutionTV[j*nPtBins+k],1.);
    hFitResolutionDipRV[j]->Add(hFitResolutionRV[j*nPtBins+k],1.);
    hFitResolutionDipOV[j]->Add(hFitResolutionOV[j*nPtBins+k],1.);
    }
    }
  */
  //Loop and Fits
  gauss->SetRange(-dipinterval,dipinterval);
  
  // pt and dip angle binning array for TGraphs: it reproduces the mean pt and angle of the intervals
  for (Int_t v=0;v<dipnbin;v++) {
    if (multdipbin[v]!=0) {
      exdip[v]=0;
      dipbin[v]=(dipbin[v]/multdipbin[v]);
      
      hFitResolutionDipTV[v]->Fit("gauss","N,R");
      sigmaResolutionDipTV[v]=gauss->GetParameter(2);
      errsigmResDipTV[v]=gauss->GetParError(2);
      
      hFitResolutionDipRV[v]->Fit("gauss","N,R");
      sigmaResolutionDipRV[v]=gauss->GetParameter(2);
      errsigmResDipRV[v]=gauss->GetParError(2);
      
      hFitResolutionDipOV[v]->Fit("gauss","N,R");
      sigmaResolutionDipOV[v]=gauss->GetParameter(2);
      errsigmResDipOV[v]=gauss->GetParError(2);
      
      hdipbin[v]=(Float_t)radice+v*bindip;
    } else {
      dipbin[v]=radice+v*bindip+lenghtbindip/2.;
      
      exdip[v]=0;
      sigmaResolutionDipTV[v]=0;
      errsigmResDipTV[v]=0;
      sigmaResolutionDipRV[v]=0;
      errsigmResDipRV[v]=0;
      sigmaResolutionDipOV[v]=0;
      errsigmResDipTV[v]=0;
      
      hdipbin[v]=(Float_t)radice+v*bindip;
    }
  }
  
  hdipbin[dipnbin]=(Float_t)hdipbin[dipnbin-1]+bindip;
  
  for (Int_t j=0;j<nPtBins;j++) {
    if (multptbin[j]!=0) {
      ptbin[j]=ptbin[j]/multptbin[j];
    } else {
      ptbin[j]=(ptbinning[j]+ptbinning[j+1])/2.;
    }
    hptbin[j]=(Float_t)ptbinning[j];
  }
  hptbin[nPtBins]=ptbinning[nPtBins];
  
  incycle=0; 
  
  for (Int_t nint=0;nint<numint;nint++) {	
    gauss->SetRange(-interval[nint],interval[nint]);
    for (Int_t k=0;k<nbinning[nint];k++) {
      ex[k+incycle]=0;
      rms=hFitPtResTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitPtResTV[k+incycle]->Fit("gauss","N,R");
      sigmaPtResTV[k+incycle]=gauss->GetParameter(2);
      errsigmPtResTV[k+incycle]=gauss->GetParError(2);
      meanPtResTV[k+incycle]=gauss->GetParameter(1);
      errmeanPtResTV[k+incycle]=gauss->GetParError(1);

      rms=hFitCurvResTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitCurvResTV[k+incycle]->Fit("gauss","N,R");
      sigmaCurvResTV[k+incycle]=gauss->GetParameter(2);
      errsigmCurvResTV[k+incycle]=gauss->GetParError(2);

      rms=hFitResolutionPtTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitResolutionPtTV[k+incycle]->Fit("gauss","N,R");
      sigmaResolutionPtTV[k+incycle]=gauss->GetParameter(2);
      errsigmResPtTV[k+incycle]=gauss->GetParError(2);
      meanResolutionPtTV[k+incycle]=gauss->GetParameter(1);
      errmeanResPtTV[k+incycle]=gauss->GetParError(1);
      
      rms=hFitResolutionPtRV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitResolutionPtRV[k+incycle]->Fit("gauss","N,R");
      sigmaResolutionPtRV[k+incycle]=gauss->GetParameter(2);
      errsigmResPtRV[k+incycle]=gauss->GetParError(2);
      
      rms=hFitResolutionPtOV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitResolutionPtOV[k+incycle]->Fit("gauss","N,R");
      sigmaResolutionPtOV[k+incycle]=gauss->GetParameter(2); 
      errsigmResPtOV[k+incycle]=gauss->GetParError(2);

      /*
      rms=hFitsigmad0rphiPtTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitsigmad0rphiPtTV[k+incycle]->Fit("gauss","N,R");
      mediumsigmad0rphiPtTV[k+incycle]=gauss->GetParameter(2);
      errmediumsigmad0rphiPtTV[k+incycle]=gauss->GetParError(2);*/
      mediumsigmad0rphiPtTV[k+incycle]=hFitsigmad0rphiPtTV[k+incycle]->GetMean();
      errmediumsigmad0rphiPtTV[k+incycle]=0;
      
      /*
      rms=hFitsigmad0rphiPtRV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitsigmad0rphiPtRV[k+incycle]->Fit("gauss","N,R");
      mediumsigmad0rphiPtRV[k+incycle]=gauss->GetParameter(2);
      errmediumsigmad0rphiPtRV[k+incycle]=gauss->GetParError(2);*/
      mediumsigmad0rphiPtRV[k+incycle]=hFitsigmad0rphiPtRV[k+incycle]->GetMean();
      errmediumsigmad0rphiPtRV[k+incycle]=0;
      
      /*  
      rms=hFitsigmad0rphiPtOV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitsigmad0rphiPtOV[k+incycle]->Fit("gauss","N,R");
      mediumsigmad0rphiPtOV[k+incycle]=gauss->GetParameter(2); 
      errmediumsigmad0rphiPtOV[k+incycle]=gauss->GetParError(2);*/
      mediumsigmad0rphiPtOV[k+incycle]=hFitsigmad0rphiPtOV[k+incycle]->GetMean();
      errmediumsigmad0rphiPtOV[k+incycle]=0;

      rms=hFitPullPtTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitPullPtTV[k+incycle]->Fit("gauss","N,R");
      sigmaPullPtTV[k+incycle]=gauss->GetParameter(2);
      errsigmPullPtTV[k+incycle]=gauss->GetParError(2);
      
      rms=hFitPullPtRV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitPullPtRV[k+incycle]->Fit("gauss","N,R");
      sigmaPullPtRV[k+incycle]=gauss->GetParameter(2);
      errsigmPullPtRV[k+incycle]=gauss->GetParError(2);
      
      rms=hFitPullPtOV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitPullPtOV[k+incycle]->Fit("gauss","N,R");
      sigmaPullPtOV[k+incycle]=gauss->GetParameter(2); 
      errsigmPullPtOV[k+incycle]=gauss->GetParError(2);

      rms=hFitd0zResPtTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitd0zResPtTV[k+incycle]->Fit("gauss","N,R");
      sigmad0zResPtTV[k+incycle]=gauss->GetParameter(2);
      errsigmd0zResPtTV[k+incycle]=gauss->GetParError(2);

      rms=hFitd0zPullPtTV[k+incycle]->GetRMS();
      gauss->SetRange(-3*rms,3*rms);
      hFitd0zPullPtTV[k+incycle]->Fit("gauss","N,R");
      sigmad0zPullPtTV[k+incycle]=gauss->GetParameter(2);
      errsigmd0zPullPtTV[k+incycle]=gauss->GetParError(2);      
    }
    incycle+=nbinning[nint];
  }
  
  // hptbin[nPtBins]=hptbin[nPtBins-1]+lenghtbin;
  
  
  hFitResPtDipTV=new TH2F("hFitResPtDipTV","The 3D histogram for Resolution in Pt*sin#theta space",nPtBins,hptbin,dipnbin,hdipbin);
  hFitResPtDipRV=new TH2F("hFitResPtDipRV","The 3D histogram for Resolution in Pt*sin#theta space",nPtBins,hptbin,dipnbin,hdipbin);
  hFitResPtDipOV=new TH2F("hFitResPtDipOV","The 3D histogram for Resolution in Pt*sin#theta space",nPtBins,hptbin,dipnbin,hdipbin);
  
  
  for (Int_t v=0;v<dipnbin;v++) { 
    incycle=0;
    for (Int_t nint=0;nint<numint;nint++) {	
      gauss->SetRange(-interval[nint],interval[nint]);
      for (Int_t k=0;k<nbinning[nint];k++) {
	//gauss->SetRange(-interval,interval);
	if(hFitResolutionOV[v*nPtBins+k+incycle]->GetEntries()>1.) {
	  // cout<<"v "<<v<<"ptbin"<<k+incycle<<endl;
	  cout<<hFitResolutionTV[v*nPtBins+k+incycle]->GetEntries()<<endl;
	  hFitResolutionTV[v*nPtBins+k+incycle]->Fit("gauss","N,R");
	  sigmaResolutionTV[v][k+incycle]=gauss->GetParameter(2);
	  hFitResPtDipTV->SetBinContent((k+incycle+1),v+1,sigmaResolutionTV[v][k+incycle]);
	  
	  hFitResolutionRV[v*nPtBins+k+incycle]->Fit("gauss","N,R");
	  sigmaResolutionRV[v][k+incycle]=gauss->GetParameter(2);
	  hFitResPtDipRV->SetBinContent((k+incycle+1),v+1,sigmaResolutionRV[v][k+incycle]);
	  
	  hFitResolutionOV[v*nPtBins+k+incycle]->Fit("gauss","N,R");
	  sigmaResolutionOV[v][k+incycle]=gauss->GetParameter(2);
	  hFitResPtDipOV->SetBinContent((k+incycle+1),v+1,sigmaResolutionOV[v][k+incycle]);
	} else {
	  sigmaResolutionTV[v][k+incycle]=0;
	  sigmaResolutionRV[v][k+incycle]=0;
	  sigmaResolutionOV[v][k+incycle]=0;
	  hFitResPtDipTV->SetBinContent((k+incycle+1),v+1,sigmaResolutionTV[v][k+incycle]);
	  hFitResPtDipRV->SetBinContent((k+incycle+1),v+1,sigmaResolutionRV[v][k+incycle]);
	  hFitResPtDipOV->SetBinContent((k+incycle+1),v+1,sigmaResolutionOV[v][k+incycle]);
	}
      }
      incycle+=nbinning[nint];
    }
  }  
  
  printf("###########  Number of good tracks: %d  ###########\n",np);
  printf("###########  Number of fake tracks: %d  ###########\n",nfaketr);


  //-------DRAWING  TGRAPH--------  
  
  const int cWidth = 1400;
  const int cHeight = 800;
  
  TCanvas* c0 = new TCanvas("c0","c0",cWidth,cHeight);
  c0->Divide(3,1);
  c0->cd(1);
  gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();    
  gr0TV = new TGraphErrors(nPtBins,ptbin,sigmaResolutionPtTV,ex,errsigmResPtTV);
  gr0TV->SetName("mygraphd0TV");
  gr0TV->SetLineColor(1);
  gr0TV->SetLineWidth(1);
  gr0TV->SetMarkerColor(gMarkerColor);
  gr0TV->SetMarkerStyle(21);
  gr0TV->SetMarkerSize(1);
  titlegraph1="d_{0}(r#phi) Resolution  (MC vertex)";
  gr0TV->SetTitle(titlegraph1);
  gr0TV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr0TV->GetYaxis()->SetTitle("#sigma [#mu m]");
  gr0TV->GetXaxis()->SetTitleSize(0.05);
  gr0TV->GetYaxis()->SetTitleSize(0.05);
  gr0TV->GetXaxis()->SetLabelSize(0.05);
  gr0TV->GetYaxis()->SetLabelSize(0.05);
  gr0TV->Draw("AP");

  c0->cd(3);
  gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();    
  TGraphErrors *gr0meanTV = new TGraphErrors(nPtBins,ptbin,meanResolutionPtTV,ex,errmeanResPtTV);
  gr0meanTV->SetName("mygraphmeand0TV");
  gr0meanTV->SetLineColor(1);
  gr0meanTV->SetLineWidth(1);
  gr0meanTV->SetMarkerColor(gMarkerColor);
  gr0meanTV->SetMarkerStyle(21);
  gr0meanTV->SetMarkerSize(1);
  titlegraph1="Mean of d_{0}(r#phi) residuals  (MC vertex)";
  gr0meanTV->SetTitle(titlegraph1);
  gr0meanTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr0meanTV->GetYaxis()->SetTitle("mean of residuals [#mu m]");
  gr0meanTV->GetXaxis()->SetTitleSize(0.05);
  gr0meanTV->GetYaxis()->SetTitleSize(0.05);
  gr0meanTV->GetXaxis()->SetLabelSize(0.05);
  gr0meanTV->GetYaxis()->SetLabelSize(0.05);
  gr0meanTV->Draw("AP");
  
  /*
    TCanvas* c1 = new TCanvas("c1","Resolution analysis dip Angle",cWidth,cHeight);
    c1->Divide(3,1);
    c1->cd(1);
    gr1TV = new TGraphErrors(dipnbin,dipbin,sigmaResolutionDipTV,exdip,errsigmResDipTV);
    gr1TV->SetName("mygraph1TV");
    gr1TV->SetLineColor(1);
    gr1TV->SetLineWidth(1);
    gr1TV->SetMarkerColor(gMarkerColor);
    gr1TV->SetMarkerStyle(21);
    gr1TV->SetMarkerSize(1);
    titlegraph1="d_{0}(z) Resolution for ";
    titlegraph1.Append(partforgraph);
    titlegraph1.Append(" using True Vertex");
    gr1TV->SetTitle(titlegraph1);
    gr1TV->GetXaxis()->SetTitle("sin#theta ");
    gr1TV->GetYaxis()->SetTitle("#sigma [#mu m]");
    gr1TV->GetXaxis()->SetTitleSize(0.05);
    gr1TV->GetYaxis()->SetTitleSize(0.05);
    gr1TV->GetXaxis()->SetLabelSize(0.05);
    gr1TV->GetYaxis()->SetLabelSize(0.05);
    gr1TV->GetXaxis()->SetRangeUser(0.7,1);
    gr1TV->GetYaxis()->SetRangeUser(0,100);
    gr1TV->Draw("AP");
    c1->Update();
        
    c1->cd(2);
    // gPad->SetLogy();
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    gr1RV = new TGraphErrors(dipnbin,dipbin,sigmaResolutionDipRV,exdip,errsigmResDipRV);
    gr1RV->SetName("mygraph1RV");
    gr1RV->SetLineColor(1);
    gr1RV->SetLineWidth(1);
    gr1RV->SetMarkerColor(gMarkerColor);
    gr1RV->SetMarkerStyle(21);
    gr1RV->SetMarkerSize(1);
    titlegraph1="d_{0}(z) Resolution for ";
    titlegraph1.Append(partforgraph);
    titlegraph1.Append(" using Reconstructed Vertex");
    gr1RV->SetTitle(titlegraph1);
    gr1RV->GetXaxis()->SetTitle("sin#theta ");
    gr1RV->GetYaxis()->SetTitle("#sigma [#mu m]");
    gr1RV->GetXaxis()->SetTitleSize(0.05);
    gr1RV->GetYaxis()->SetTitleSize(0.05);
    gr1RV->GetXaxis()->SetLabelSize(0.05);
    gr1RV->GetYaxis()->SetLabelSize(0.05);
    gr1RV->GetXaxis()->SetRangeUser(0.7,1);
    gr1RV->GetYaxis()->SetRangeUser(0,100);
    gr1RV->Draw("AP");
    c1->Update();


  
    c1->cd(3);
    // gPad->SetLogy();
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    gr1OV = new TGraphErrors(dipnbin,dipbin,sigmaResolutionDipOV,exdip,errsigmResDipOV);
    gr1OV->SetName("mygraph1OV");
    gr1OV->SetLineColor(1);
    gr1OV->SetLineWidth(1);
    gr1OV->SetMarkerColor(gMarkerColor);
    gr1OV->SetMarkerStyle(21);
    gr1OV->SetMarkerSize(1);
    titlegraph1="d_{0}(z) Resolution for ";
    titlegraph1.Append(partforgraph);
    titlegraph1.Append(" using Vertex on the Fly");
    gr1OV->SetTitle(titlegraph1);
    gr1OV->GetXaxis()->SetTitle("sin#theta ");
    gr1OV->GetYaxis()->SetTitle("#sigma [#mu m]");
    gr1OV->GetXaxis()->SetTitleSize(0.05);
    gr1OV->GetYaxis()->SetTitleSize(0.05);
    gr1OV->GetXaxis()->SetLabelSize(0.05);
    gr1OV->GetYaxis()->SetLabelSize(0.05);
    gr1OV->GetXaxis()->SetRangeUser(0.7,1);
    gr1OV->GetYaxis()->SetRangeUser(0,100);
    gr1OV->Draw("AP");
    c1->Update();
    */    

    TCanvas* c2 = new TCanvas("c2","c2",cWidth,cHeight);
    c2->Divide(3,1);
    c2->cd(1);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();    
    grPullTV = new TGraphErrors(nPtBins,ptbin,sigmaPullPtTV,ex,errsigmPullPtTV);
    grPullTV->SetName("mygraphPullTV");
    grPullTV->SetMinimum(0);
    grPullTV->SetMaximum(1.5);
    grPullTV->SetLineColor(1);
    grPullTV->SetLineWidth(1);
    grPullTV->SetMarkerColor(gMarkerColor);
    grPullTV->SetMarkerStyle(21);
    grPullTV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Pull  (MC vertex)";
    grPullTV->SetTitle(titlegraph1);
    grPullTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grPullTV->GetYaxis()->SetTitle("");
    grPullTV->GetXaxis()->SetTitleSize(0.05);
    grPullTV->GetYaxis()->SetTitleSize(0.05);
    grPullTV->GetXaxis()->SetLabelSize(0.05);
    grPullTV->GetYaxis()->SetLabelSize(0.05);
    grPullTV->Draw("AP");
    c0->cd(2);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grPullTV->Draw("AP");

    c2->cd(2);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grPullRV = new TGraphErrors(nPtBins,ptbin,sigmaPullPtRV,ex,errsigmPullPtRV);
    grPullRV->SetName("mygraphPullRV");
    grPullRV->SetMinimum(0);
    grPullRV->SetMaximum(1.2);
    grPullRV->SetLineColor(1);
    grPullRV->SetLineWidth(1);
    grPullRV->SetMarkerColor(gMarkerColor);
    grPullRV->SetMarkerStyle(21);
    grPullRV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Pull for using Reconstructed Vertex";
    grPullRV->SetTitle(titlegraph1);
    grPullRV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grPullRV->GetYaxis()->SetTitle("");
    grPullRV->GetXaxis()->SetTitleSize(0.05);
    grPullRV->GetYaxis()->SetTitleSize(0.05);
    grPullRV->GetXaxis()->SetLabelSize(0.05);
    grPullRV->GetYaxis()->SetLabelSize(0.05);
    grPullRV->Draw("AP");

    c2->cd(3);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grPullOV = new TGraphErrors(nPtBins,ptbin,sigmaPullPtOV,ex,errsigmPullPtOV);
    grPullOV->SetName("mygraphPullOV");
    grPullOV->SetMinimum(0);
    grPullOV->SetMaximum(1.2);
    grPullOV->SetLineColor(1);
    grPullOV->SetLineWidth(1);
    grPullOV->SetMarkerColor(gMarkerColor);
    grPullOV->SetMarkerStyle(21);
    grPullOV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Pull using Vertex on the Fly";
    grPullOV->SetTitle(titlegraph1);
    grPullOV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grPullOV->GetYaxis()->SetTitle("");
    grPullOV->GetXaxis()->SetTitleSize(0.05);
    grPullOV->GetYaxis()->SetTitleSize(0.05);
    grPullOV->GetXaxis()->SetLabelSize(0.05);
    grPullOV->GetYaxis()->SetLabelSize(0.05);
    grPullOV->Draw("AP");


    TCanvas* c3 = new TCanvas("c3","c3",cWidth,cHeight);
    c3->Divide(2,1);   
    c3->cd(1);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grd0zResTV = new TGraphErrors(nPtBins,ptbin,sigmad0zResPtTV,ex,errsigmd0zResPtTV);
    grd0zResTV->SetName("mygraphd0zResTV");
    grd0zResTV->SetLineColor(1);
    grd0zResTV->SetLineWidth(1);
    grd0zResTV->SetMarkerColor(gMarkerColor);
    grd0zResTV->SetMarkerStyle(21);
    grd0zResTV->SetMarkerSize(1);
    titlegraph1="d_{0}(z) Resolution  (MC vertex)";
    grd0zResTV->SetTitle(titlegraph1);
    grd0zResTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grd0zResTV->GetYaxis()->SetTitle("#sigma [#mu m]");
    grd0zResTV->GetXaxis()->SetTitleSize(0.05);
    grd0zResTV->GetYaxis()->SetTitleSize(0.05);
    grd0zResTV->GetXaxis()->SetLabelSize(0.05);
    grd0zResTV->GetYaxis()->SetLabelSize(0.05);
    grd0zResTV->Draw("AP");
    c3->cd(2);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grd0zPullTV = new TGraphErrors(nPtBins,ptbin,sigmad0zPullPtTV,ex,errsigmd0zPullPtTV);
    grd0zPullTV->SetName("mygraphd0zPullTV");
    grd0zPullTV->SetMinimum(0);
    grd0zPullTV->SetMaximum(1.5);
    grd0zPullTV->SetLineColor(1);
    grd0zPullTV->SetLineWidth(1);
    grd0zPullTV->SetMarkerColor(gMarkerColor);
    grd0zPullTV->SetMarkerStyle(21);
    grd0zPullTV->SetMarkerSize(1);
    titlegraph1="d_{0}(z) Pull  (MC vertex)";
    grd0zPullTV->SetTitle(titlegraph1);
    grd0zPullTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grd0zPullTV->GetYaxis()->SetTitle("");
    grd0zPullTV->GetXaxis()->SetTitleSize(0.05);
    grd0zPullTV->GetYaxis()->SetTitleSize(0.05);
    grd0zPullTV->GetXaxis()->SetLabelSize(0.05);
    grd0zPullTV->GetYaxis()->SetLabelSize(0.05);
    grd0zPullTV->Draw("AP");

    TCanvas* c4 = new TCanvas("c4","c4",cWidth,cHeight);
    c4->Divide(3,1);
    c4->cd(1);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grPtResTV = new TGraphErrors(nPtBins,ptbin,sigmaPtResTV,ex,errsigmPtResTV);
    grPtResTV->SetName("mygraphPtResTV");
    grPtResTV->SetLineColor(1);
    grPtResTV->SetLineWidth(1);
    grPtResTV->SetMarkerColor(gMarkerColor);
    grPtResTV->SetMarkerStyle(21);
    grPtResTV->SetMarkerSize(1);
    titlegraph1="p_{t}mes/p_{t}true - 1 ";
    grPtResTV->SetTitle(titlegraph1);
    grPtResTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grPtResTV->GetYaxis()->SetTitle("#sigma(p_{t}mes/p_{t}true- 1 ) ");
    grPtResTV->GetXaxis()->SetTitleSize(0.05);
    grPtResTV->GetYaxis()->SetTitleSize(0.05);
    grPtResTV->GetXaxis()->SetLabelSize(0.05);
    grPtResTV->GetYaxis()->SetLabelSize(0.05);
    grPtResTV->GetYaxis()->SetRangeUser(0,0.025);
    grPtResTV->Draw("AP");
    c4->cd(2);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grPtResMeanTV = new TGraphErrors(nPtBins,ptbin,meanPtResTV,ex,errmeanPtResTV);
    grPtResMeanTV->SetName("mygraphPtResMeanTV");
    grPtResMeanTV->SetLineColor(1);
    grPtResMeanTV->SetLineWidth(1);
    grPtResMeanTV->SetMarkerColor(gMarkerColor);
    grPtResMeanTV->SetMarkerStyle(21);
    grPtResMeanTV->SetMarkerSize(1);
    titlegraph1="p_{t}mes/p_{t}true - 1 ";
    grPtResMeanTV->SetTitle(titlegraph1);
    grPtResMeanTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grPtResMeanTV->GetYaxis()->SetTitle("mean(p_{t}mes/p_{t}true - 1) ");
    grPtResMeanTV->GetXaxis()->SetTitleSize(0.05);
    grPtResMeanTV->GetYaxis()->SetTitleSize(0.05);
    grPtResMeanTV->GetXaxis()->SetLabelSize(0.05);
    grPtResMeanTV->GetYaxis()->SetLabelSize(0.05);
    grPtResMeanTV->GetYaxis()->SetRangeUser(0,0.025);
    grPtResMeanTV->Draw("AP");
    c4->cd(3);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grCurvResTV = new TGraphErrors(nPtBins,ptbin,sigmaCurvResTV,ex,errsigmCurvResTV);
    grCurvResTV->SetName("mygraphCurvResTV");
    grCurvResTV->SetLineColor(1);
    grCurvResTV->SetLineWidth(1);
    grCurvResTV->SetMarkerColor(gMarkerColor);
    grCurvResTV->SetMarkerStyle(21);
    grCurvResTV->SetMarkerSize(1);
    titlegraph1="1/p_{t} Resolution ";
    grCurvResTV->SetTitle(titlegraph1);
    grCurvResTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grCurvResTV->GetYaxis()->SetTitle("#sigma(1/p_{t}) [1/(GeV/c)]");
    grCurvResTV->GetXaxis()->SetTitleSize(0.05);
    grCurvResTV->GetYaxis()->SetTitleSize(0.05);
    grCurvResTV->GetXaxis()->SetLabelSize(0.05);
    grCurvResTV->GetYaxis()->SetLabelSize(0.05);
    grCurvResTV->Draw("AP");


    /*
    TCanvas* c5 = new TCanvas("c5","Resolution analysis 2",cWidth,cHeight);
    c5->Divide(3,1);
    c5->cd(1);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();    
    grsigmad0rphiTV = new TGraphErrors(nPtBins,ptbin,mediumsigmad0rphiPtTV,ex,errmediumsigmad0rphiPtTV);
    grsigmad0rphiTV->SetName("mygraphsigmad0rphiTV");
    grsigmad0rphiTV->SetLineColor(1);
    grsigmad0rphiTV->SetLineWidth(1);
    grsigmad0rphiTV->SetMarkerColor(gMarkerColor);
    grsigmad0rphiTV->SetMarkerStyle(21);
    grsigmad0rphiTV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Resolution Estimated from Cov. Matr (RelateToVtx) for ";
    titlegraph1.Append(" using TV");
    grsigmad0rphiTV->SetTitle(titlegraph1);
    grsigmad0rphiTV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grsigmad0rphiTV->GetYaxis()->SetTitle("#sigma [#mu m]");
    grsigmad0rphiTV->GetXaxis()->SetTitleSize(0.05);
    grsigmad0rphiTV->GetYaxis()->SetTitleSize(0.05);
    grsigmad0rphiTV->GetXaxis()->SetLabelSize(0.05);
    grsigmad0rphiTV->GetYaxis()->SetLabelSize(0.05);
    grsigmad0rphiTV->Draw("AP");
    c5->cd(2);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grsigmad0rphiRV = new TGraphErrors(nPtBins,ptbin,mediumsigmad0rphiPtRV,ex,errmediumsigmad0rphiPtRV);
    grsigmad0rphiRV->SetName("mygraphsigmad0rphiRV");
    grsigmad0rphiRV->SetLineColor(1);
    grsigmad0rphiRV->SetLineWidth(1);
    grsigmad0rphiRV->SetMarkerColor(gMarkerColor);
    grsigmad0rphiRV->SetMarkerStyle(21);
    grsigmad0rphiRV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Resolution Estimated from Cov. Matr (RelateToVtx) for ";
    titlegraph1.Append(" using RVtx");
    grsigmad0rphiRV->SetTitle(titlegraph1);
    grsigmad0rphiRV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grsigmad0rphiRV->GetYaxis()->SetTitle("#sigma [#mu m]");
    grsigmad0rphiRV->GetXaxis()->SetTitleSize(0.05);
    grsigmad0rphiRV->GetYaxis()->SetTitleSize(0.05);
    grsigmad0rphiRV->GetXaxis()->SetLabelSize(0.05);
    grsigmad0rphiRV->GetYaxis()->SetLabelSize(0.05);
    grsigmad0rphiRV->Draw("AP");
    c5->cd(3);
    gPad->SetLogx(); gPad->SetGridx(); gPad->SetGridy();
    grsigmad0rphiOV = new TGraphErrors(nPtBins,ptbin,mediumsigmad0rphiPtOV,ex,errmediumsigmad0rphiPtOV);
    grsigmad0rphiOV->SetName("mygraphsigmad0rphiOV");
    grsigmad0rphiOV->SetLineColor(1);
    grsigmad0rphiOV->SetLineWidth(1);
    grsigmad0rphiOV->SetMarkerColor(gMarkerColor);
    grsigmad0rphiOV->SetMarkerStyle(21);
    grsigmad0rphiOV->SetMarkerSize(1);
    titlegraph1="d_{0}(r#phi) Resolution Estimated from Cov. Matr (RelateToVtx) for ";
    titlegraph1.Append(" using VtxOth");
    grsigmad0rphiOV->SetTitle(titlegraph1);
    grsigmad0rphiOV->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
    grsigmad0rphiOV->GetYaxis()->SetTitle("#sigma [#mu m]");
    grsigmad0rphiOV->GetXaxis()->SetTitleSize(0.05);
    grsigmad0rphiOV->GetYaxis()->SetTitleSize(0.05);
    grsigmad0rphiOV->GetXaxis()->SetLabelSize(0.05);
    grsigmad0rphiOV->GetYaxis()->SetLabelSize(0.05);
    grsigmad0rphiOV->Draw("AP");
    */



    
    //---------CLOSE EVRYTHING AND SAVE------------
   TFile *outfile = new TFile("ResolutionsAnalysis.root","recreate");
    outfile->cd(); 
    for (Int_t i=0;i<nPtBins;i++) {
      hFitResolutionPtTV[i]->Write();
      hFitResolutionPtRV[i]->Write();
      hFitResolutionPtOV[i]->Write();
      delete hFitResolutionPtTV[i];//=0;
      delete	hFitResolutionPtRV[i];//=0;
      delete	hFitResolutionPtOV[i];
      hFitsigmad0rphiPtTV[i]->Write();
      hFitsigmad0rphiPtRV[i]->Write();
      hFitsigmad0rphiPtOV[i]->Write();
      delete hFitsigmad0rphiPtTV[i];//=0;
      delete	hFitsigmad0rphiPtRV[i];//=0;
      delete	hFitsigmad0rphiPtOV[i];
      hFitPullPtTV[i]->Write();
      hFitPullPtRV[i]->Write();
      hFitPullPtOV[i]->Write();
      delete hFitPullPtTV[i];//=0;
      delete	hFitPullPtRV[i];//=0;
      delete	hFitPullPtOV[i];
      hFitd0zResPtTV[i]->Write();
      hFitd0zPullPtTV[i]->Write();
      delete hFitd0zResPtTV[i];//=0;
      delete hFitd0zPullPtTV[i];//=0;
      hFitPtResTV[i]->Write();
      hFitCurvResTV[i]->Write();
      delete hFitPtResTV[i];//=0;
      delete hFitCurvResTV[i];//=0;
      for(Int_t j=0;j<dipnbin;j++) {
	hFitResolutionTV[j*nPtBins+i]->Write();
	hFitResolutionRV[j*nPtBins+i]->Write();
	hFitResolutionOV[j*nPtBins+i]->Write();
	delete  hFitResolutionTV[j*nPtBins+i];//=0;
	delete hFitResolutionRV[j*nPtBins+i];//=0;
	delete hFitResolutionOV[j*nPtBins+i];
      }
    }
    for (Int_t i=0;i<dipnbin;i++) {
      hFitResolutionDipTV[i]->Write();
      hFitResolutionDipRV[i]->Write();
      hFitResolutionDipOV[i]->Write();
      delete hFitResolutionDipTV[i];//=0;
      delete hFitResolutionDipRV[i];//=0;
      delete hFitResolutionDipOV[i];//=0;
    }
  
    //    delete hFitResolutionPtTV;
    hFitResPtDipTV->Write();
    hFitResPtDipRV->Write();
    hFitResPtDipOV->Write();
    delete hFitResPtDipOV;
    delete hFitResPtDipTV;
    delete hFitResPtDipRV;

    gr0TV->Write();
    gr0meanTV->Write();
    //gr0RV->Write();
    //gr0OV->Write();
    //gr1TV->Write();
    //gr1RV->Write();
    //gr1OV->Write();

    //grsigmad0rphiTV->Write();
    //grsigmad0rphiRV->Write();
    //grsigmad0rphiOV->Write();
  
    grd0zResTV->Write();
    grd0zPullTV->Write();
    grPtResTV->Write();
    grPtResMeanTV->Write();
    grCurvResTV->Write();
    

    grPullTV->Write();
    grPullRV->Write();
    //grPullOV->Write();
    
    outfile->Close();
    
    return; 
}
//-----------------------------------------------------------------------------
Int_t ITSnCluster(Int_t ITSf) {

  Int_t nTPCclusters = (Int_t)(ITSf/1000);
  ITSf -= nTPCclusters*1000;
  Int_t sign,aITSf,n0,n3,aux,ITSsel;//be careful!:working with integer!
  n0=0;
  n3=0;
  aITSf=TMath::Abs(ITSf);
  if(aITSf==0)sign=0;
  else if (aITSf==ITSf)sign=1;
  else if(aITSf==-ITSf)sign=-1;
  
  if (aITSf/100<1)n0++;
  n3=aITSf/300;
  aux=(aITSf-(aITSf/100)*100);
  if(aux/10<1)n0++;
  n3+=aux/30;
  aux=(aux-(aux/10)*10);
  if (aux==0)n0++;
  n3+=aux/3;
  ITSsel=3+n3-n0;
  
  if(ITSsel>6) {
    printf("Wrong ITSflag assignment! \n");
    return 99;
  }
  return ITSsel*sign;
}
//-----------------------------------------------------------------------------
Bool_t kITSrefit(Double_t sigmd0TV) {
  //return TRUE if track was kITSrefit, FALSE otherwise (see AliTrackProperties.C)
  Bool_t kITSref=kTRUE;
  if (sigmd0TV<0.)kITSref=kFALSE;
  return kITSref;
}
//-----------------------------------------------------------------------------
void PlotResolutions() {

  TCanvas *c1a = new TCanvas("c1a","d0 resolution",0,0,800,800);
  c1a->SetLogx();
  c1a->SetGridx();
  c1a->SetGridy();
  TCanvas *c1b = new TCanvas("c1b","d0 mean",0,0,800,800);
  c1b->SetLogx();
  c1b->SetGridx();
  c1b->SetGridy();
  TCanvas *c1c = new TCanvas("c1c","d0 pull",0,0,800,800);
  c1c->SetLogx();
  c1c->SetGridx();
  c1c->SetGridy();

  TCanvas *c2a = new TCanvas("c2a","z0 resolution",0,0,800,800);
  c2a->SetLogx();
  c2a->SetGridx();
  c2a->SetGridy();
  TCanvas *c2c = new TCanvas("c2c","z0 pull",0,0,800,800);
  c2c->SetLogx();
  c2c->SetGridx();
  c2c->SetGridy();

  TCanvas *c3a = new TCanvas("c3a","pt resolution",0,0,800,800);
  c3a->SetLogx();
  c3a->SetGridx();
  c3a->SetGridy();
  TCanvas *c3b = new TCanvas("c3b","pt mean",0,0,800,800);
  c3b->SetLogx();
  c3b->SetGridx();
  c3b->SetGridy();

  TCanvas *c4 = new TCanvas("c4","d0 residuals",0,0,800,800);
  c4->Divide(3,3);
  c4_1->SetLogy();
  c4_2->SetLogy();
  c4_3->SetLogy();
  c4_4->SetLogy();
  c4_5->SetLogy();
  c4_6->SetLogy();
  c4_7->SetLogy();
  c4_8->SetLogy();
  c4_9->SetLogy();

  TCanvas *c5 = new TCanvas("c5","pt residuals",0,0,800,800);
  c5->Divide(3,3);
  c5_1->SetLogy();
  c5_2->SetLogy();
  c5_3->SetLogy();
  c5_4->SetLogy();
  c5_5->SetLogy();
  c5_6->SetLogy();
  c5_7->SetLogy();
  c5_8->SetLogy();
  c5_9->SetLogy();


  TLegend *leg1=new TLegend(0.5,0.5,0.9,0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  TLegend *leg2=new TLegend(0.5,0.5,0.9,0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  TGraph *mark20=new TGraph();
  mark20->SetMarkerColor(4);
  mark20->SetLineColor(4);
  mark20->SetMarkerStyle(20);
  leg1->AddEntry(mark20,"pions","p");
  leg2->AddEntry(mark20,"pions","l");
  TGraph *mark21=new TGraph();
  mark21->SetMarkerColor(2);
  mark21->SetLineColor(2);
  mark21->SetMarkerStyle(21);
  leg1->AddEntry(mark21,"electrons","p");
  leg2->AddEntry(mark21,"electrons","l");
  TGraph *mark25=new TGraph();
  mark25->SetMarkerColor(1);
  mark25->SetLineColor(1);
  mark25->SetMarkerStyle(25);
  leg1->AddEntry(mark25,"electrons (no brem)","p");
  leg2->AddEntry(mark25,"electrons (no brem)","l");
  TGraph *mark22=new TGraph();
  mark22->SetMarkerColor(6);
  mark22->SetLineColor(6);
  mark22->SetMarkerStyle(22);
  leg1->AddEntry(mark22,"kaons","p");
  leg2->AddEntry(mark22,"kaons","l");
  TGraph *mark23=new TGraph();
  mark23->SetMarkerColor(8);
  mark23->SetLineColor(8);
  mark23->SetMarkerStyle(23);
  leg1->AddEntry(mark23,"protons","p");
  leg2->AddEntry(mark23,"protons","l");

  // pions
  TFile *f1=new TFile("boxPiv416Release_zero/ResolutionsAnalysis_Pi.root");
  c1a->cd();
  mygraphd0TV->Draw("a,p");
  mygraphd0TV->SetMarkerColor(4);
  mygraphd0TV->SetMarkerStyle(20);
  leg1->Draw();
  c1b->cd();
  mygraphmeand0TV->Draw("a,p");
  mygraphmeand0TV->SetMarkerColor(4);
  mygraphmeand0TV->SetMarkerStyle(20);
  leg1->Draw();
  c1c->cd();
  mygraphPullTV->Draw("a,p");
  mygraphPullTV->SetMarkerColor(4);
  mygraphPullTV->SetMarkerStyle(20);
  leg1->Draw();
  c2a->cd();
  mygraphd0zResTV->Draw("a,p");
  mygraphd0zResTV->SetMarkerColor(4);
  mygraphd0zResTV->SetMarkerStyle(20);
  leg1->Draw();
  c2c->cd();
  mygraphd0zPullTV->Draw("a,p");
  mygraphd0zPullTV->SetMarkerColor(4);
  mygraphd0zPullTV->SetMarkerStyle(20);
  leg1->Draw();
  c3a->cd();
  mygraphPtResTV->Draw("a,p");
  mygraphPtResTV->SetMarkerColor(4);
  mygraphPtResTV->SetMarkerStyle(20);
  leg1->Draw();
  c3b->cd();
  mygraphPtResMeanTV->Draw("a,p");
  mygraphPtResMeanTV->SetMarkerColor(4);
  mygraphPtResMeanTV->SetMarkerStyle(20);
  leg1->Draw();
  c4->cd(1);
  d0PtResolution0binTV->Draw();
  d0PtResolution0binTV->SetLineColor(4);
  c4->cd(2);
  d0PtResolution1binTV->Draw();
  d0PtResolution1binTV->SetLineColor(4);
  c4->cd(3);
  d0PtResolution2binTV->Draw();
  d0PtResolution2binTV->SetLineColor(4);
  c4->cd(4);
  d0PtResolution3binTV->Draw();
  d0PtResolution3binTV->SetLineColor(4);
  c4->cd(5);
  d0PtResolution4binTV->Draw();
  d0PtResolution4binTV->SetLineColor(4);
  c4->cd(6);
  d0PtResolution5binTV->Draw();
  d0PtResolution5binTV->SetLineColor(4);
  c4->cd(7);
  d0PtResolution6binTV->Draw();
  d0PtResolution6binTV->SetLineColor(4);
  c4->cd(8);
  d0PtResolution7binTV->Draw();
  d0PtResolution7binTV->SetLineColor(4);
  c4->cd(9);
  d0PtResolution8binTV->Draw();
  d0PtResolution8binTV->SetLineColor(4);
  leg2->Draw();
  c5->cd(1);
  PtResolution0binTV->Draw();
  PtResolution0binTV->SetLineColor(4);
  c5->cd(2);
  PtResolution1binTV->Draw();
  PtResolution1binTV->SetLineColor(4);
  c5->cd(3);
  PtResolution2binTV->Draw();
  PtResolution2binTV->SetLineColor(4);
  c5->cd(4);
  PtResolution3binTV->Draw();
  PtResolution3binTV->SetLineColor(4);
  c5->cd(5);
  PtResolution4binTV->Draw();
  PtResolution4binTV->SetLineColor(4);
  c5->cd(6);
  PtResolution5binTV->Draw();
  PtResolution5binTV->SetLineColor(4);
  c5->cd(7);
  PtResolution6binTV->Draw();
  PtResolution6binTV->SetLineColor(4);
  c5->cd(8);
  PtResolution7binTV->Draw();
  PtResolution7binTV->SetLineColor(4);
  c5->cd(9);
  PtResolution8binTV->Draw();
  PtResolution8binTV->SetLineColor(4);
  leg2->Draw();



  // kaons
  TFile *f3=new TFile("boxKav416Release_zero/ResolutionsAnalysis_Ka.root");
  c1a->cd();
  mygraphd0TV->Draw("p");
  mygraphd0TV->SetMarkerColor(6);
  mygraphd0TV->SetMarkerStyle(22);
  c1b->cd();
  mygraphmeand0TV->Draw("p");
  mygraphmeand0TV->SetMarkerColor(6);
  mygraphmeand0TV->SetMarkerStyle(22);
  c1c->cd();
  mygraphPullTV->Draw("p");
  mygraphPullTV->SetMarkerColor(6);
  mygraphPullTV->SetMarkerStyle(22);
  c2a->cd();
  mygraphd0zResTV->Draw("p");
  mygraphd0zResTV->SetMarkerColor(6);
  mygraphd0zResTV->SetMarkerStyle(22);
  c2c->cd();
  mygraphd0zPullTV->Draw("p");
  mygraphd0zPullTV->SetMarkerColor(6);
  mygraphd0zPullTV->SetMarkerStyle(22);
  c3a->cd();
  mygraphPtResTV->Draw("p");
  mygraphPtResTV->SetMarkerColor(6);
  mygraphPtResTV->SetMarkerStyle(22);
  c3b->cd();
  mygraphPtResMeanTV->Draw("p");
  mygraphPtResMeanTV->SetMarkerColor(6);
  mygraphPtResMeanTV->SetMarkerStyle(22);
  c4->cd(1);
  d0PtResolution0binTV->Draw("same");
  d0PtResolution0binTV->SetLineColor(6);
  c4->cd(2);
  d0PtResolution1binTV->Draw("same");
  d0PtResolution1binTV->SetLineColor(6);
  c4->cd(3);
  d0PtResolution2binTV->Draw("same");
  d0PtResolution2binTV->SetLineColor(6);
  c4->cd(4);
  d0PtResolution3binTV->Draw("same");
  d0PtResolution3binTV->SetLineColor(6);
  c4->cd(5);
  d0PtResolution4binTV->Draw("same");
  d0PtResolution4binTV->SetLineColor(6);
  c4->cd(6);
  d0PtResolution5binTV->Draw("same");
  d0PtResolution5binTV->SetLineColor(6);
  c4->cd(7);
  d0PtResolution6binTV->Draw("same");
  d0PtResolution6binTV->SetLineColor(6);
  c4->cd(8);
  d0PtResolution7binTV->Draw("same");
  d0PtResolution7binTV->SetLineColor(6);
  c4->cd(9);
  d0PtResolution8binTV->Draw("same");
  d0PtResolution8binTV->SetLineColor(6);
  c5->cd(1);
  PtResolution0binTV->Draw("same");
  PtResolution0binTV->SetLineColor(6);
  c5->cd(2);
  PtResolution1binTV->Draw("same");
  PtResolution1binTV->SetLineColor(6);
  c5->cd(3);
  PtResolution2binTV->Draw("same");
  PtResolution2binTV->SetLineColor(6);
  c5->cd(4);
  PtResolution3binTV->Draw("same");
  PtResolution3binTV->SetLineColor(6);
  c5->cd(5);
  PtResolution4binTV->Draw("same");
  PtResolution4binTV->SetLineColor(6);
  c5->cd(6);
  PtResolution5binTV->Draw("same");
  PtResolution5binTV->SetLineColor(6);
  c5->cd(7);
  PtResolution6binTV->Draw("same");
  PtResolution6binTV->SetLineColor(6);
  c5->cd(8);
  PtResolution7binTV->Draw("same");
  PtResolution7binTV->SetLineColor(6);
  c5->cd(9);
  PtResolution8binTV->Draw("same");
  PtResolution8binTV->SetLineColor(6);


  // protons
  TFile *f4=new TFile("boxPrv416Release_zero/ResolutionsAnalysis_Pr.root");
  c1a->cd();
  mygraphd0TV->Draw("p");
  mygraphd0TV->SetMarkerColor(8);
  mygraphd0TV->SetMarkerStyle(23);
  c1b->cd();
  mygraphmeand0TV->Draw("p");
  mygraphmeand0TV->SetMarkerColor(8);
  mygraphmeand0TV->SetMarkerStyle(23);
  c1c->cd();
  mygraphPullTV->Draw("p");
  mygraphPullTV->SetMarkerColor(8);
  mygraphPullTV->SetMarkerStyle(23);
  c2a->cd();
  mygraphd0zResTV->Draw("p");
  mygraphd0zResTV->SetMarkerColor(8);
  mygraphd0zResTV->SetMarkerStyle(23);
  c2c->cd();
  mygraphd0zPullTV->Draw("p");
  mygraphd0zPullTV->SetMarkerColor(8);
  mygraphd0zPullTV->SetMarkerStyle(23);
  c3a->cd();
  mygraphPtResTV->Draw("p");
  mygraphPtResTV->SetMarkerColor(8);
  mygraphPtResTV->SetMarkerStyle(23);
  c3b->cd();
  mygraphPtResMeanTV->Draw("p");
  mygraphPtResMeanTV->SetMarkerColor(8);
  mygraphPtResMeanTV->SetMarkerStyle(23);
  c4->cd(1);
  d0PtResolution0binTV->Draw("same");
  d0PtResolution0binTV->SetLineColor(8);
  c4->cd(2);
  d0PtResolution1binTV->Draw("same");
  d0PtResolution1binTV->SetLineColor(8);
  c4->cd(3);
  d0PtResolution2binTV->Draw("same");
  d0PtResolution2binTV->SetLineColor(8);
  c4->cd(4);
  d0PtResolution3binTV->Draw("same");
  d0PtResolution3binTV->SetLineColor(8);
  c4->cd(5);
  d0PtResolution4binTV->Draw("same");
  d0PtResolution4binTV->SetLineColor(8);
  c4->cd(6);
  d0PtResolution5binTV->Draw("same");
  d0PtResolution5binTV->SetLineColor(8);
  c4->cd(7);
  d0PtResolution6binTV->Draw("same");
  d0PtResolution6binTV->SetLineColor(8);
  c4->cd(8);
  d0PtResolution7binTV->Draw("same");
  d0PtResolution7binTV->SetLineColor(8);
  c4->cd(9);
  d0PtResolution8binTV->Draw("same");
  d0PtResolution8binTV->SetLineColor(8);
  c5->cd(1);
  PtResolution0binTV->Draw("same");
  PtResolution0binTV->SetLineColor(8);
  c5->cd(2);
  PtResolution1binTV->Draw("same");
  PtResolution1binTV->SetLineColor(8);
  c5->cd(3);
  PtResolution2binTV->Draw("same");
  PtResolution2binTV->SetLineColor(8);
  c5->cd(4);
  PtResolution3binTV->Draw("same");
  PtResolution3binTV->SetLineColor(8);
  c5->cd(5);
  PtResolution4binTV->Draw("same");
  PtResolution4binTV->SetLineColor(8);
  c5->cd(6);
  PtResolution5binTV->Draw("same");
  PtResolution5binTV->SetLineColor(8);
  c5->cd(7);
  PtResolution6binTV->Draw("same");
  PtResolution6binTV->SetLineColor(8);
  c5->cd(8);
  PtResolution7binTV->Draw("same");
  PtResolution7binTV->SetLineColor(8);
  c5->cd(9);
  PtResolution8binTV->Draw("same");
  PtResolution8binTV->SetLineColor(8);


  // electrons
  TFile *f2=new TFile("boxElv416Release_zero/ResolutionsAnalysis_El.root");
  c1a->cd();
  mygraphd0TV->Draw("p");
  mygraphd0TV->SetMarkerColor(2);
  mygraphd0TV->SetMarkerStyle(21);
  c1b->cd();
  mygraphmeand0TV->Draw("p");
  mygraphmeand0TV->SetMarkerColor(2);
  mygraphmeand0TV->SetMarkerStyle(21);
  c1c->cd();
  mygraphPullTV->Draw("p");
  mygraphPullTV->SetMarkerColor(2);
  mygraphPullTV->SetMarkerStyle(21);
  c2a->cd();
  mygraphd0zResTV->Draw("p");
  mygraphd0zResTV->SetMarkerColor(2);
  mygraphd0zResTV->SetMarkerStyle(21);
  c2c->cd();
  mygraphd0zPullTV->Draw("p");
  mygraphd0zPullTV->SetMarkerColor(2);
  mygraphd0zPullTV->SetMarkerStyle(21);
  c3a->cd();
  mygraphPtResTV->Draw("p");
  mygraphPtResTV->SetMarkerColor(2);
  mygraphPtResTV->SetMarkerStyle(21);
  c3b->cd();
  mygraphPtResMeanTV->Draw("p");
  mygraphPtResMeanTV->SetMarkerColor(2);
  mygraphPtResMeanTV->SetMarkerStyle(21);
  c4->cd(1);
  d0PtResolution0binTV->Draw("same");
  d0PtResolution0binTV->SetLineColor(2);
  c4->cd(2);
  d0PtResolution1binTV->Draw("same");
  d0PtResolution1binTV->SetLineColor(2);
  c4->cd(3);
  d0PtResolution2binTV->Draw("same");
  d0PtResolution2binTV->SetLineColor(2);
  c4->cd(4);
  d0PtResolution3binTV->Draw("same");
  d0PtResolution3binTV->SetLineColor(2);
  c4->cd(5);
  d0PtResolution4binTV->Draw("same");
  d0PtResolution4binTV->SetLineColor(2);
  c4->cd(6);
  d0PtResolution5binTV->Draw("same");
  d0PtResolution5binTV->SetLineColor(2);
  c4->cd(7);
  d0PtResolution6binTV->Draw("same");
  d0PtResolution6binTV->SetLineColor(2);
  c4->cd(8);
  d0PtResolution7binTV->Draw("same");
  d0PtResolution7binTV->SetLineColor(2);
  c4->cd(9);
  d0PtResolution8binTV->Draw("same");
  d0PtResolution8binTV->SetLineColor(2);
  c5->cd(1);
  PtResolution0binTV->Draw("same");
  PtResolution0binTV->SetLineColor(2);
  c5->cd(2);
  PtResolution1binTV->Draw("same");
  PtResolution1binTV->SetLineColor(2);
  c5->cd(3);
  PtResolution2binTV->Draw("same");
  PtResolution2binTV->SetLineColor(2);
  c5->cd(4);
  PtResolution3binTV->Draw("same");
  PtResolution3binTV->SetLineColor(2);
  c5->cd(5);
  PtResolution4binTV->Draw("same");
  PtResolution4binTV->SetLineColor(2);
  c5->cd(6);
  PtResolution5binTV->Draw("same");
  PtResolution5binTV->SetLineColor(2);
  c5->cd(7);
  PtResolution6binTV->Draw("same");
  PtResolution6binTV->SetLineColor(2);
  c5->cd(8);
  PtResolution7binTV->Draw("same");
  PtResolution7binTV->SetLineColor(2);
  c5->cd(9);
  PtResolution8binTV->Draw("same");
  PtResolution8binTV->SetLineColor(2);

  // electrons, no bremsstrahlung
  TFile *f2=new TFile("boxElv416Release_zero/ResolutionsAnalysis_ElNoBrem.root");
  c1a->cd();
  mygraphd0TV->Draw("p");
  mygraphd0TV->SetMarkerColor(1);
  mygraphd0TV->SetMarkerStyle(25);
  c1b->cd();
  mygraphmeand0TV->Draw("p");
  mygraphmeand0TV->SetMarkerColor(1);
  mygraphmeand0TV->SetMarkerStyle(25);
  c1c->cd();
  mygraphPullTV->Draw("p");
  mygraphPullTV->SetMarkerColor(1);
  mygraphPullTV->SetMarkerStyle(25);
  c2a->cd();
  mygraphd0zResTV->Draw("p");
  mygraphd0zResTV->SetMarkerColor(1);
  mygraphd0zResTV->SetMarkerStyle(25);
  c2c->cd();
  mygraphd0zPullTV->Draw("p");
  mygraphd0zPullTV->SetMarkerColor(1);
  mygraphd0zPullTV->SetMarkerStyle(25);
  c3a->cd();
  mygraphPtResTV->Draw("p");
  mygraphPtResTV->SetMarkerColor(1);
  mygraphPtResTV->SetMarkerStyle(25);
  c3b->cd();
  mygraphPtResMeanTV->Draw("p");
  mygraphPtResMeanTV->SetMarkerColor(1);
  mygraphPtResMeanTV->SetMarkerStyle(25);
  c4->cd(1);
  d0PtResolution0binTV->Draw("same");
  d0PtResolution0binTV->SetLineColor(1);
  c4->cd(2);
  d0PtResolution1binTV->Draw("same");
  d0PtResolution1binTV->SetLineColor(1);
  c4->cd(3);
  d0PtResolution2binTV->Draw("same");
  d0PtResolution2binTV->SetLineColor(1);
  c4->cd(4);
  d0PtResolution3binTV->Draw("same");
  d0PtResolution3binTV->SetLineColor(1);
  c4->cd(5);
  d0PtResolution4binTV->Draw("same");
  d0PtResolution4binTV->SetLineColor(1);
  c4->cd(6);
  d0PtResolution5binTV->Draw("same");
  d0PtResolution5binTV->SetLineColor(1);
  c4->cd(7);
  d0PtResolution6binTV->Draw("same");
  d0PtResolution6binTV->SetLineColor(1);
  c4->cd(8);
  d0PtResolution7binTV->Draw("same");
  d0PtResolution7binTV->SetLineColor(1);
  c4->cd(9);
  d0PtResolution8binTV->Draw("same");
  d0PtResolution8binTV->SetLineColor(1);
  c5->cd(1);
  PtResolution0binTV->Draw("same");
  PtResolution0binTV->SetLineColor(1);
  c5->cd(2);
  PtResolution1binTV->Draw("same");
  PtResolution1binTV->SetLineColor(1);
  c5->cd(3);
  PtResolution2binTV->Draw("same");
  PtResolution2binTV->SetLineColor(1);
  c5->cd(4);
  PtResolution3binTV->Draw("same");
  PtResolution3binTV->SetLineColor(1);
  c5->cd(5);
  PtResolution4binTV->Draw("same");
  PtResolution4binTV->SetLineColor(1);
  c5->cd(6);
  PtResolution5binTV->Draw("same");
  PtResolution5binTV->SetLineColor(1);
  c5->cd(7);
  PtResolution6binTV->Draw("same");
  PtResolution6binTV->SetLineColor(1);
  c5->cd(8);
  PtResolution7binTV->Draw("same");
  PtResolution7binTV->SetLineColor(1);
  c5->cd(9);
  PtResolution8binTV->Draw("same");
  PtResolution8binTV->SetLineColor(1);


  return;
}
