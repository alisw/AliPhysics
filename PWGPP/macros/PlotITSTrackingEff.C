void PlotITSTrackingEff(TString filename="ITS.Performance.root",
			TString part="Pi",
			TString partforgraph="pions",
			Int_t pdgcode=211,
			Bool_t useAbsPdg=kTRUE,
			Bool_t nofakes=kTRUE,
			Bool_t askITSrefit=kTRUE,
			Int_t minTPCcls=1) 
{
  //
  // Macro to plot ITS tracking efficiency from ITS.Performance.root
  // A. Dainese
  //  

  //Open File
  if(gSystem->AccessPathName(filename.Data())) {
    printf("file not found!\n");  
    return;
  }
  TFile *file= TFile::Open(filename.Data());
  cout<<"Open File "<<filename<<endl;   

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

  
  TString titlegraph1,titlefile;  
  Int_t  ntracks,nITSsel,nfakes=0;
  const Int_t nPtBins=10;  
  Double_t meanpt[nPtBins];
  ///////////////////////////////////////////////////////////////////////////
  Double_t nprimgen[nPtBins]; 
  Float_t nfiles=48.; if(!useAbsPdg) nfiles*=0.5; 
  for(Int_t j=0; j<nPtBins; j++) nprimgen[j]=nfiles*50.;// *18.;
  // the factor 50 is for the cut |eta|<0.9
  // the factor 18 is for the cut 0.9<|eta|<1.4
  ///////////////////////////////////////////////////////////////////////////
  Double_t n6ITSclspart[nPtBins],n5ITSclspart[nPtBins],n4ITSclspart[nPtBins],n3ITSclspart[nPtBins],n2ITSclspart[nPtBins],n1ITSclspart[nPtBins];
  Double_t n6ITSclsparttoTPC[nPtBins],n5ITSclsparttoTPC[nPtBins],n4ITSclsparttoTPC[nPtBins],n3ITSclsparttoTPC[nPtBins],n2ITSclsparttoTPC[nPtBins];
  Double_t nSPD[nPtBins],n2ITS[nPtBins],n3ITS[nPtBins],n4ITS[nPtBins],n5ITS[nPtBins],n6ITS[nPtBins],n456ITS[nPtBins],ptbinning[2*nPtBins],ntottr[nPtBins];
  Double_t nTPCnorm6[nPtBins],n1ITSnorm6[nPtBins],n2ITSnorm6[nPtBins],n3ITSnorm6[nPtBins],n4ITSnorm6[nPtBins],n5ITSnorm6[nPtBins],n6ITSnorm6[nPtBins],n456ITSnorm6[nPtBins],n4ITSnorm5[nPtBins],n5ITSnorm5[nPtBins],n3ITSnorm4[nPtBins],n4ITSnorm4[nPtBins],n2ITSnorm3[nPtBins],n3ITSnorm3[nPtBins],n1ITSnorm2[nPtBins],n2ITSnorm2[nPtBins];
  Double_t nSPDtoGen[nPtBins],n2ITStoGen[nPtBins],n3ITStoGen[nPtBins],n4ITStoGen[nPtBins],n5ITStoGen[nPtBins],n6ITStoGen[nPtBins],n456ITStoGen[nPtBins],nTPCtoGen[nPtBins];
  Double_t nTPCto6cls[nPtBins],n1ITSto6cls[nPtBins],n2ITSto6cls[nPtBins],n3ITSto6cls[nPtBins],n4ITSto6cls[nPtBins],n5ITSto6cls[nPtBins],n4ITSto5cls[nPtBins],n5ITSto5cls[nPtBins],n3ITSto4cls[nPtBins],n4ITSto4cls[nPtBins],n2ITSto3cls[nPtBins],n3ITSto3cls[nPtBins],n1ITSto2cls[nPtBins],n2ITSto2cls[nPtBins],n6ITSto6cls[nPtBins],n456ITSto6cls[nPtBins];
  Double_t n330[nPtBins],n331[nPtBins],n332[nPtBins],n313[nPtBins],n323[nPtBins],n133[nPtBins],n233[nPtBins],nOTHERS[nPtBins],nSSD[nPtBins],nSDD[nPtBins];

  Double_t sigmnSPD[nPtBins],sigmn4ITS[nPtBins],sigmn5ITS[nPtBins],sigmn6ITS[nPtBins],sigmn456ITS[nPtBins],ptbinning[2*nPtBins],ntottr[nPtBins];
  Double_t sigmnSPDtoGen[nPtBins],sigmn2ITStoGen[nPtBins],sigmn3ITStoGen[nPtBins],sigmn4ITStoGen[nPtBins],sigmn5ITStoGen[nPtBins],sigmn6ITStoGen[nPtBins],sigmn456ITStoGen[nPtBins],sigmnTPCtoGen[nPtBins];
  Double_t sigmnTPCto6cls[nPtBins],sigmn1ITSto6cls[nPtBins],sigmn2ITSto6cls[nPtBins],sigmn3ITSto6cls[nPtBins],sigmn4ITSto6cls[nPtBins],sigmn5ITSto6cls[nPtBins],sigmn6ITSto6cls[nPtBins],sigmn456ITSto6cls[nPtBins],sigmn4ITSto5cls[nPtBins],sigmn5ITSto5cls[nPtBins],sigmn3ITSto4cls[nPtBins],sigmn4ITSto4cls[nPtBins],sigmn2ITSto3cls[nPtBins],sigmn3ITSto3cls[nPtBins],sigmn1ITSto2cls[nPtBins],sigmn2ITSto2cls[nPtBins];
  Double_t sigmn330[nPtBins],sigmn331[nPtBins],sigmn332[nPtBins],sigmn313[nPtBins],sigmn323[nPtBins],sigmn133[nPtBins],sigmn233[nPtBins],sigmnOTHERS[nPtBins],sigmnSSD[nPtBins],sigmnSDD[nPtBins];
  Double_t sigmn6ITSclsparttoTPC[nPtBins],sigmn5ITSclsparttoTPC[nPtBins],sigmn4ITSclsparttoTPC[nPtBins],sigmn3ITSclsparttoTPC[nPtBins],sigmn2ITSclsparttoTPC[nPtBins];

  Double_t errx[nPtBins];
  
  TGraphErrors *gr4,*gr5,*gr6,*gr456,*grSPDtoGen,*gr2toGen,*gr3toGen,*gr4toGen,*gr5toGen,*gr6toGen,*gr456toGen,*grTPCtoGen,*grTPCto6cls,*gr1to6cls,*gr2to6cls,*gr3to6cls,*gr4to6cls,*gr5to6cls,*gr6to6cls,*gr456to6cls,*gr4to5cls,*gr5to5cls,*gr3to4cls,*gr4to4cls,*gr2to3cls,*gr3to3cls,*gr1to2cls,*gr2to2cls,*gr330,*gr331,*gr332,*gr313,*gr323,*gr133,*gr233,*grOTHERS,*grSSD,*grSPD,*gr6clspart,*gr5clspart,*gr4clspart,*gr2ITSclsparttoTPC,*gr3ITSclsparttoTPC,*gr4ITSclsparttoTPC,*gr5ITSclsparttoTPC,*gr6ITSclsparttoTPC;
  
  for(Int_t j=0;j<nPtBins;j++) {
    n1ITSclspart[j]=0.;
    n2ITSclspart[j]=0.;
    n3ITSclspart[j]=0.;
    n4ITSclspart[j]=0.;
    n5ITSclspart[j]=0.;
    n6ITSclspart[j]=0.;
    nSPD[j]=0.;
    n2ITS[j]=0.;
    n3ITS[j]=0.;
    n4ITS[j]=0.;
    n5ITS[j]=0.;
    n6ITS[j]=0.;
    n456ITS[j]=0.;
    n1ITSnorm6[j]=0.;
    n2ITSnorm6[j]=0.;
    n3ITSnorm6[j]=0.;
    n4ITSnorm6[j]=0.;
    n5ITSnorm6[j]=0.;
    n6ITSnorm6[j]=0.;
    n456ITSnorm6[j]=0.;
    n4ITSnorm5[j]=0.;
    n5ITSnorm5[j]=0.;
    n3ITSnorm4[j]=0.;
    n4ITSnorm4[j]=0.;
    n2ITSnorm3[j]=0.;
    n3ITSnorm3[j]=0.;
    n1ITSnorm2[j]=0.;
    n2ITSnorm2[j]=0.;
    nSPDtoGen[j]=0.;
    n2ITStoGen[j]=0.;
    n3ITStoGen[j]=0.;
    n4ITStoGen[j]=0.;
    n5ITStoGen[j]=0.;
    n6ITStoGen[j]=0.;
    n456ITStoGen[j]=0.;
    nTPCtoGen[j]=0.;
    n1ITSto6cls[j]=0.;
    n2ITSto6cls[j]=0.;
    n3ITSto6cls[j]=0.;
    n4ITSto6cls[j]=0.;
    n5ITSto6cls[j]=0.;
    n6ITSto6cls[j]=0.;
    n4ITSto5cls[j]=0.;
    n5ITSto5cls[j]=0.;
    n3ITSto4cls[j]=0.;
    n4ITSto4cls[j]=0.;
    n2ITSto3cls[j]=0.;
    n3ITSto3cls[j]=0.;
    n1ITSto2cls[j]=0.;
    n2ITSto2cls[j]=0.;
    n456ITSto6cls[j]=0.;
    nTPCto6cls[j]=0.;
    n330[j]=0;
    n331[j]=0;
    n332[j]=0;
    n313[j]=0;
    n323[j]=0;
    n133[j]=0;
    n233[j]=0;
    nOTHERS[j]=0;
    nSDD[j]=0.;
    nSSD[j]=0.;
    errx[j]=0.;
    ntottr[j]=0;
    meanpt[j]=0;
  }

  ptbinning[0]=0.15;
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
  ptbinning[14]=8.5;
  ptbinning[15]=11.5;
  ptbinning[16]=18.5;
  ptbinning[17]=21.5;
  ptbinning[18]=27.5;
  ptbinning[19]=32.5;
  

  
  

  // ------------ Starting Loop on Tracks ----------
  ntracks=ntTracks->GetEntries();
  printf("number of Tracks %d \n",ntracks);
  
  for (Int_t j=0;j<ntracks;j++) {
    if(j%10000==0) printf("Reading track %d\n",j);
    ntTracks->GetEvent(j); 
    if(TMath::Abs(d0True)>1.) continue; // only primaries
    if(TMath::Abs(eta)>0.9) continue;   // only tracks with |etapart|<0.9
    //if(TMath::Abs(eta)<0.9 || TMath::Abs(eta)>1.4) continue;   // only tracks with 0.9<|etapart|<1.4
    Float_t theta = 2.*TMath::ATan(TMath::Exp(-eta));
    //if(TMath::Abs(TMath::Abs(phi)-0.5*TMath::Pi())>0.5 || TMath::Abs(theta-0.5*TMath::Pi())>0.5) continue;
    //if(phi<(.5+0.25)*TMath::Pi() || phi>(.5+0.45)*TMath::Pi()) continue;
    if(ITSflag<0) nfakes++;
    //printf("  %d\n",ITSflag);
    if(nofakes && ITSflag<0) continue;        // reject fakes
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
    nITSsel=ITSnCluster(ITSflag);


    if(nTPCcls<=minTPCcls) continue;

    if(useAbsPdg) pdg=TMath::Abs(pdg);
    if(pdg!=pdgcode) continue;

    for (Int_t k=0; k<nPtBins; k++) {
      if(ptbinning[2*k]<pt&&pt<ptbinning[2*k+1]) {
	//if(nITSclspart==3) {nITSclspart=5;}else{nITSclspart=-1;}
	if(ptmes>0.) {
	  ntottr[k]+=1.;
	  meanpt[k]+=pt;
	  if(nITSclspart==6) n6ITSclspart[k]+=1.;
	  if(nITSclspart==5) n5ITSclspart[k]+=1.;
	  if(nITSclspart==4) n4ITSclspart[k]+=1.;
	  if(nITSclspart==3) n3ITSclspart[k]+=1.;
	  if(nITSclspart==2) n2ITSclspart[k]+=1.;
	}	
	if(askITSrefit && sigmad0TV<0.) continue;   
	if(nITSclsassign!=nITSsel) cout<<" ERROR  "<<ITSflag<<"  "<<nITSclsassign<<"  "<<nITSsel<<"   "<<ITSflagorig<<endl; 
	
	if(nITSsel==2) n2ITS[k]+=1.;
	if(nITSsel==3) n3ITS[k]+=1.;
	if(nITSsel==4) n4ITS[k]+=1.;
	if(nITSsel==5) n5ITS[k]+=1.;
	if(nITSsel==6) n6ITS[k]+=1.;
	
	if(((Int_t)ITSflag)%10==3) nSPD[k]+=1.; // both SPD
	
	if(nITSclspart==6) {
	  if(nITSsel==1) n1ITSnorm6[k]+=1.;
	  if(nITSsel==2) n2ITSnorm6[k]+=1.;
	  if(nITSsel==3) n3ITSnorm6[k]+=1.;
	  if(nITSsel==4) n4ITSnorm6[k]+=1.;
	  if(nITSsel==5) n5ITSnorm6[k]+=1.;
	  if(nITSsel==6) n6ITSnorm6[k]+=1.;
	}
	if(nITSclspart==5) {
	  if(nITSsel==4) n4ITSnorm5[k]+=1.;
	  if(nITSsel==5) n5ITSnorm5[k]+=1.;
	}
	if(nITSclspart==4) {
	  if(nITSsel==3) n3ITSnorm4[k]+=1.;
	  if(nITSsel==4) n4ITSnorm4[k]+=1.;
	}
	if(nITSclspart==3) {
	  if(nITSsel==2) n2ITSnorm3[k]+=1.;
	  if(nITSsel==3) n3ITSnorm3[k]+=1.;
	}
	if(nITSclspart==2) {
	  if(nITSsel==1) n1ITSnorm2[k]+=1.;
	  if(nITSsel==2) n2ITSnorm2[k]+=1.;
	}

	if(k==nPtBins-1 && nITSclspart==5 && nITSsel==4) printf("nClustersPart %d nClustersTrack %d Map %d\n",nITSclspart,nITSsel,ITSflag);
	
	//5 points array
	if(nITSclspart==5) {
	  if(ITSflag==332) {
	    n332[k]+=1.;
	  } else if(ITSflag==331) {
	    n331[k]+=1.;
	  } else if(ITSflag==323) {
	    n323[k]+=1.;
	  } else if(ITSflag==313) {
	    n313[k]+=1.;
	  } else if(ITSflag==233) {
	    n233[k]+=1.;
	  } else if(ITSflag==133) {
	    n133[k]+=1.;
	  }
	}
	if(ITSflag==330) {
	  n330[k]+=1.;
	} else if(ITSflag!=333) {
	  nOTHERS[k]+=1.;
	}
	
      }
    }

  } //end loop on tracks

  printf("%d  %d  %d  %d  %d\n",n2ITS[9],n3ITS[9],n4ITS[9],n5ITS[9],n6ITS[9]);

  for(Int_t k=0;k<nPtBins;k++) {
    if(n6ITSclspart[k]<1) n6ITSclspart[k]=1.;
    if(n5ITSclspart[k]<1) n5ITSclspart[k]=1.;
    if(n4ITSclspart[k]<1) n4ITSclspart[k]=1.;
    if(n3ITSclspart[k]<1) n3ITSclspart[k]=1.;
    if(n2ITSclspart[k]<1) n2ITSclspart[k]=1.;
    if(n1ITSclspart[k]<1) n1ITSclspart[k]=1.;

    //cout<<n1ITSclspart[k]<<" "<<n2ITSclspart[k]<<" "<<n3ITSclspart[k]<<" "<<n4ITSclspart[k]<<" "<<n5ITSclspart[k]<<" "<<n6ITSclspart[k]<<endl;

    if(ntottr[k]!=0) {
      meanpt[k]=meanpt[k]/ntottr[k];

      nSPDtoGen[k]=nSPD[k]/nprimgen[k];
      n2ITStoGen[k]=n2ITS[k]/nprimgen[k];
      n3ITStoGen[k]=n3ITS[k]/nprimgen[k];
      n4ITStoGen[k]=n4ITS[k]/nprimgen[k];
      n5ITStoGen[k]=n5ITS[k]/nprimgen[k];
      n6ITStoGen[k]=n6ITS[k]/nprimgen[k];
      nTPCtoGen[k]=ntottr[k]/nprimgen[k];
      n456ITStoGen[k]=n4ITStoGen[k]+n5ITStoGen[k]+n6ITStoGen[k];
      
      n1ITSto6cls[k]=n1ITSnorm6[k]/n6ITSclspart[k];
      n2ITSto6cls[k]=n2ITSnorm6[k]/n6ITSclspart[k];
      n3ITSto6cls[k]=n3ITSnorm6[k]/n6ITSclspart[k];
      n4ITSto6cls[k]=n4ITSnorm6[k]/n6ITSclspart[k];
      n5ITSto6cls[k]=n5ITSnorm6[k]/n6ITSclspart[k];
      n6ITSto6cls[k]=n6ITSnorm6[k]/n6ITSclspart[k];
      n4ITSto5cls[k]=n4ITSnorm5[k]/n5ITSclspart[k];
      n5ITSto5cls[k]=n5ITSnorm5[k]/n5ITSclspart[k];
      n3ITSto4cls[k]=n3ITSnorm4[k]/n4ITSclspart[k];
      n4ITSto4cls[k]=n4ITSnorm4[k]/n4ITSclspart[k];
      n2ITSto3cls[k]=n2ITSnorm3[k]/n3ITSclspart[k];
      n3ITSto3cls[k]=n3ITSnorm3[k]/n3ITSclspart[k];
      n1ITSto2cls[k]=n1ITSnorm2[k]/n2ITSclspart[k];
      n2ITSto2cls[k]=n2ITSnorm2[k]/n2ITSclspart[k];

      nTPCto6cls[k]=ntottr[k]/n6ITSclspart[k];
      n456ITSto6cls[k]=n4ITSto6cls[k]+n5ITSto6cls[k]+n6ITSto6cls[k];
      //cout<<meanpt[k]<<"  "<<n6ITSnorm6[k]<<"  "<<n6ITSclspart[k]<<endl;

      n6ITSclsparttoTPC[k]=n6ITSclspart[k]/ntottr[k];
      n5ITSclsparttoTPC[k]=n5ITSclspart[k]/ntottr[k];
      n4ITSclsparttoTPC[k]=n4ITSclspart[k]/ntottr[k];
      n3ITSclsparttoTPC[k]=n3ITSclspart[k]/ntottr[k];
      n2ITSclsparttoTPC[k]=n2ITSclspart[k]/ntottr[k];
      
      n4ITS[k]=n4ITS[k]/ntottr[k];
      n5ITS[k]=n5ITS[k]/ntottr[k];
      n6ITS[k]=n6ITS[k]/ntottr[k];

      n330[k]=n330[k]/ntottr[k];
      n331[k]=n331[k]/n5ITSclspart[k];
      n332[k]=n332[k]/n5ITSclspart[k];
      n313[k]=n313[k]/n5ITSclspart[k];
      n323[k]=n323[k]/n5ITSclspart[k];
      n133[k]=n133[k]/n5ITSclspart[k];
      n233[k]=n233[k]/n5ITSclspart[k];
      nOTHERS[k]=nOTHERS[k]/ntottr[k];
      n456ITS[k]=n4ITS[k]+n5ITS[k]+n6ITS[k];
      nSSD[k]=n133[k]+n233[k];
      nSDD[k]=n323[k]+n313[k];

      
      sigmn4ITS[k]=TMath::Sqrt(n4ITS[k]*(1-n4ITS[k])/ntottr[k]);
      sigmn5ITS[k]=TMath::Sqrt(n5ITS[k]*(1-n5ITS[k])/ntottr[k]);
      sigmn6ITS[k]=TMath::Sqrt(n6ITS[k]*(1-n6ITS[k])/ntottr[k]);
      sigmn456ITS[k]=TMath::Sqrt(n456ITS[k]*(1-n456ITS[k])/ntottr[k]);

      sigmn2ITSclsparttoTPC[k]=TMath::Sqrt(n2ITSclsparttoTPC[k]*(1-n2ITSclsparttoTPC[k])/ntottr[k]);
      sigmn3ITSclsparttoTPC[k]=TMath::Sqrt(n3ITSclsparttoTPC[k]*(1-n3ITSclsparttoTPC[k])/ntottr[k]);
      sigmn4ITSclsparttoTPC[k]=TMath::Sqrt(n4ITSclsparttoTPC[k]*(1-n4ITSclsparttoTPC[k])/ntottr[k]);
      sigmn5ITSclsparttoTPC[k]=TMath::Sqrt(n5ITSclsparttoTPC[k]*(1-n5ITSclsparttoTPC[k])/ntottr[k]);
      sigmn6ITSclsparttoTPC[k]=TMath::Sqrt(n6ITSclsparttoTPC[k]*(1-n6ITSclsparttoTPC[k])/ntottr[k]);

      sigmnSPDtoGen[k]=TMath::Sqrt(nSPDtoGen[k]*(1-nSPDtoGen[k])/nprimgen[k]);
      sigmn2ITStoGen[k]=TMath::Sqrt(n2ITStoGen[k]*(1-n2ITStoGen[k])/nprimgen[k]);
      sigmn3ITStoGen[k]=TMath::Sqrt(n3ITStoGen[k]*(1-n3ITStoGen[k])/nprimgen[k]);
      sigmn4ITStoGen[k]=TMath::Sqrt(n4ITStoGen[k]*(1-n4ITStoGen[k])/nprimgen[k]);
      sigmn5ITStoGen[k]=TMath::Sqrt(n5ITStoGen[k]*(1-n5ITStoGen[k])/nprimgen[k]);
      sigmn6ITStoGen[k]=TMath::Sqrt(n6ITStoGen[k]*(1-n6ITStoGen[k])/nprimgen[k]);
      sigmn456ITStoGen[k]=TMath::Sqrt(n456ITStoGen[k]*(1-n456ITStoGen[k])/nprimgen[k]);
      sigmnTPCtoGen[k]=TMath::Sqrt(nTPCtoGen[k]*(1-nTPCtoGen[k])/nprimgen[k]);
      sigmnTPCto6cls[k]=TMath::Sqrt(nTPCto6cls[k]*(1-nTPCto6cls[k])/n6ITSclspart[k]);
      
      sigmn1ITSto6cls[k]=TMath::Sqrt(n1ITSto6cls[k]*(1-n1ITSto6cls[k])/n6ITSclspart[k]);
      sigmn2ITSto6cls[k]=TMath::Sqrt(n2ITSto6cls[k]*(1-n2ITSto6cls[k])/n6ITSclspart[k]);
      sigmn3ITSto6cls[k]=TMath::Sqrt(n3ITSto6cls[k]*(1-n3ITSto6cls[k])/n6ITSclspart[k]);
      sigmn4ITSto6cls[k]=TMath::Sqrt(n4ITSto6cls[k]*(1-n4ITSto6cls[k])/n6ITSclspart[k]);
      sigmn5ITSto6cls[k]=TMath::Sqrt(n5ITSto6cls[k]*(1-n5ITSto6cls[k])/n6ITSclspart[k]);
      sigmn6ITSto6cls[k]=TMath::Sqrt(n6ITSto6cls[k]*(1-n6ITSto6cls[k])/n6ITSclspart[k]);
      sigmn4ITSto5cls[k]=TMath::Sqrt(n4ITSto5cls[k]*(1-n4ITSto5cls[k])/n5ITSclspart[k]);
      sigmn5ITSto5cls[k]=TMath::Sqrt(n5ITSto5cls[k]*(1-n5ITSto5cls[k])/n5ITSclspart[k]);
      sigmn3ITSto4cls[k]=TMath::Sqrt(n3ITSto4cls[k]*(1-n3ITSto4cls[k])/n4ITSclspart[k]);
      sigmn4ITSto4cls[k]=TMath::Sqrt(n4ITSto4cls[k]*(1-n4ITSto4cls[k])/n4ITSclspart[k]);
      sigmn2ITSto3cls[k]=TMath::Sqrt(n2ITSto3cls[k]*(1-n2ITSto3cls[k])/n3ITSclspart[k]);
      sigmn3ITSto3cls[k]=TMath::Sqrt(n3ITSto3cls[k]*(1-n3ITSto3cls[k])/n3ITSclspart[k]);
      sigmn1ITSto2cls[k]=TMath::Sqrt(n1ITSto2cls[k]*(1-n1ITSto2cls[k])/n2ITSclspart[k]);
      sigmn2ITSto2cls[k]=TMath::Sqrt(n2ITSto2cls[k]*(1-n2ITSto2cls[k])/n2ITSclspart[k]);
      sigmn456ITSto6cls[k]=TMath::Sqrt(n456ITSto6cls[k]*(1-n456ITSto6cls[k])/n6ITSclspart[k]);
      
      sigmn330[k]=TMath::Sqrt(n330[k]*(1-n330[k])/ntottr[k]);
      sigmn331[k]=TMath::Sqrt(n331[k]*(1-n331[k])/ntottr[k]);
      sigmn332[k]=TMath::Sqrt(n332[k]*(1-n332[k])/ntottr[k]);
      sigmn313[k]=TMath::Sqrt(n313[k]*(1-n313[k])/ntottr[k]);
      sigmn323[k]=TMath::Sqrt(n323[k]*(1-n323[k])/ntottr[k]);
      sigmn133[k]=TMath::Sqrt(n133[k]*(1-n133[k])/ntottr[k]);
      sigmn233[k]=TMath::Sqrt(n233[k]*(1-n233[k])/ntottr[k]);

     
      sigmnOTHERS[k]=TMath::Sqrt(nOTHERS[k]*(1-nOTHERS[k])/ntottr[k]);
      sigmnSSD[k]=TMath::Sqrt((n133[k]+n233[k])*(1-n133[k]-n233[k])/ntottr[k]);
      sigmnSDD[k]=TMath::Sqrt((n323[k]+n313[k])*(1-n323[k]-n313[k])/ntottr[k]);
    } else {
      meanpt[k]=(ptbinning[2*k+1]+ptbinning[2*k])/2.;
    }
    
  } // end loop on pt bins

  //-------DRAWING  TGRAPH--------  
  
  const int cWidth = 500;
  const int cHeight = (int)(500*(29.7/21.));
  
  TCanvas* ceff0 = new TCanvas("ceff0","ceff0",cWidth,cHeight);
  ceff0->SetLogx();
  ceff0->SetGridy();

  gr4 = new TGraphErrors(nPtBins,meanpt,n4ITS,errx,sigmn4ITS);
  gr4->SetName("mygraph04ITS");
  gr4->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr4->SetLineColor(1);
  gr4->SetLineWidth(1);
  gr4->SetMarkerColor(2);
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerSize(1);
  gr4->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4->GetYaxis()->SetTitle("ITS efficiency (normalized to TPC trks)");
  gr4->SetMaximum(1.2);
  gr4->SetMinimum(0.);
  gr4->Draw("AP");
  
  gr5 = new TGraphErrors(nPtBins,meanpt,n5ITS,errx,sigmn5ITS);
  gr5->SetName("mygraph05ITS");
  gr5->SetLineColor(1);
  gr5->SetLineWidth(1);
  gr5->SetMarkerColor(8);
  gr5->SetMarkerStyle(21);
  gr5->SetMarkerSize(1.);
  gr5->SetTitle(titlegraph1);
  gr5->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr5->GetYaxis()->SetTitle("");
  gr5->Draw("P");

  gr6 = new TGraphErrors(nPtBins,meanpt,n6ITS,errx,sigmn6ITS);
  gr6->SetName("mygraph06ITS");
  gr6->SetLineColor(1);
  gr6->SetLineWidth(1);
  gr6->SetMarkerColor(4);
  gr6->SetMarkerStyle(21);
  gr6->SetMarkerSize(1.);
  gr6->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr6->GetYaxis()->SetTitle("");
  gr6->Draw("P");

  gr456 = new TGraphErrors(nPtBins,meanpt,n456ITS,errx,sigmn456ITS);
  gr456->SetName("mygraph456ITS");
  gr456->SetLineColor(1);
  gr456->SetLineWidth(1);
  gr456->SetMarkerColor(1);
  gr456->SetMarkerStyle(21);
  gr456->SetMarkerSize(1.);
  gr456->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr456->GetYaxis()->SetTitle("");
  gr456->Draw("P");

  TLegend *l0 = new TLegend(0.5,0.5,0.9,0.9);
  l0->SetFillStyle(0);
  l0->SetBorderSize(0);
  l0->AddEntry(gr4,"kTPCin + 4 ITS cls","p");
  l0->AddEntry(gr5,"kTPCin + 5 ITS cls","p");
  l0->AddEntry(gr6,"kTPCin + 6 ITS cls","p");
  l0->AddEntry(gr456,"kTPCin + 4or5or6 ITS cls","p");
  l0->Draw();

  TCanvas* ceff1 = new TCanvas("ceff1","ceff1",cWidth,cHeight);
  ceff1->SetLogx();
  ceff1->SetGridy();
  
  gr4toGen = new TGraphErrors(nPtBins,meanpt,n4ITStoGen,errx,sigmn4ITStoGen);
  gr4toGen->SetName("mygraph04ITStoGen");
  gr4toGen->SetTitle("Physical tracking efficiency (B=0.5T, no misal.)");
  gr4toGen->SetLineColor(1);
  gr4toGen->SetLineWidth(1);
  gr4toGen->SetMarkerColor(2);
  gr4toGen->SetMarkerStyle(24);
  gr4toGen->SetMarkerSize(1.);
  gr4toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4toGen->GetYaxis()->SetTitle("Physical efficiency (normalized to generated)");
  gr4toGen->SetMaximum(1.2);
  gr4toGen->SetMinimum(0.);
  gr4toGen->Draw("AP");

  gr2toGen = new TGraphErrors(nPtBins,meanpt,n2ITStoGen,errx,sigmn2ITStoGen);
  gr2toGen->SetName("mygraph02ITStoGen");
  gr2toGen->SetLineColor(1);
  gr2toGen->SetLineWidth(1);
  gr2toGen->SetMarkerColor(1);
  gr2toGen->SetMarkerStyle(20);
  gr2toGen->SetMarkerSize(1.);
  gr2toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr2toGen->GetYaxis()->SetTitle("");
  gr2toGen->Draw("P");

  gr3toGen = new TGraphErrors(nPtBins,meanpt,n3ITStoGen,errx,sigmn3ITStoGen);
  gr3toGen->SetName("mygraph03ITStoGen");
  gr3toGen->SetLineColor(1);
  gr3toGen->SetLineWidth(1);
  gr3toGen->SetMarkerColor(2);
  gr3toGen->SetMarkerStyle(21);
  gr3toGen->SetMarkerSize(1.);
  gr3toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr3toGen->GetYaxis()->SetTitle("");
  gr3toGen->Draw("P");

  gr5toGen = new TGraphErrors(nPtBins,meanpt,n5ITStoGen,errx,sigmn5ITStoGen);
  gr5toGen->SetName("mygraph05ITStoGen");
  gr5toGen->SetLineColor(1);
  gr5toGen->SetLineWidth(1);
  gr5toGen->SetMarkerColor(8);
  gr5toGen->SetMarkerStyle(24);
  gr5toGen->SetMarkerSize(1.);
  gr5toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr5toGen->GetYaxis()->SetTitle("");
  gr5toGen->Draw("P");

  gr6toGen = new TGraphErrors(nPtBins,meanpt,n6ITStoGen,errx,sigmn6ITStoGen);
  gr6toGen->SetName("mygraph06ITStoGen");
  gr6toGen->SetLineColor(1);
  gr6toGen->SetLineWidth(1);
  gr6toGen->SetMarkerColor(4);
  gr6toGen->SetMarkerStyle(24);
  gr6toGen->SetMarkerSize(1.);
  gr6toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr6toGen->GetYaxis()->SetTitle("");
  gr6toGen->Draw("P");

  gr456toGen = new TGraphErrors(nPtBins,meanpt,n456ITStoGen,errx,sigmn456ITStoGen);
  gr456toGen->SetName("mygraph0456ITStoGen");
  gr456toGen->SetLineColor(1);
  gr456toGen->SetLineWidth(1);
  gr456toGen->SetMarkerColor(1);
  gr456toGen->SetMarkerStyle(24);
  gr456toGen->SetMarkerSize(1.);
  gr456toGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr456toGen->GetYaxis()->SetTitle("");
  //gr456toGen->Draw("P");

  
  grTPCtoGen = new TGraphErrors(nPtBins,meanpt,nTPCtoGen,errx,sigmnTPCtoGen);
  grTPCtoGen->SetName("mygraph0TPCtoGen");
  grTPCtoGen->SetLineColor(1);
  grTPCtoGen->SetLineWidth(1);
  grTPCtoGen->SetMarkerColor(6);
  grTPCtoGen->SetMarkerStyle(20);
  grTPCtoGen->SetMarkerSize(1.);
  grTPCtoGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  grTPCtoGen->GetYaxis()->SetTitle("");
  grTPCtoGen->Draw("P");

  grSPDtoGen = new TGraphErrors(nPtBins,meanpt,nSPDtoGen,errx,sigmnSPDtoGen);
  grSPDtoGen->SetName("mygraph0SPDtoGen");
  grSPDtoGen->SetLineColor(1);
  grSPDtoGen->SetLineWidth(1);
  grSPDtoGen->SetMarkerColor(9);
  grSPDtoGen->SetMarkerStyle(22);
  grSPDtoGen->SetMarkerSize(1.3);
  grSPDtoGen->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  grSPDtoGen->GetYaxis()->SetTitle("");
  grSPDtoGen->Draw("P");
  

  TLegend *l1 = new TLegend(0.5,0.5,0.9,0.9);
  l1->SetFillStyle(0);
  l1->SetBorderSize(0);
  l1->AddEntry(grTPCtoGen,"kTPCin","p");
  l1->AddEntry(grSPDtoGen,"kTPCin + 2 SPD cls","p");
  l1->AddEntry(gr2toGen,"kTPCin + 2 ITS cls","p");
  l1->AddEntry(gr3toGen,"kTPCin + 3 ITS cls","p");
  l1->AddEntry(gr4toGen,"kTPCin + 4 ITS cls","p");
  l1->AddEntry(gr5toGen,"kTPCin + 5 ITS cls","p");
  l1->AddEntry(gr6toGen,"kTPCin + 6 ITS cls","p");
  //l1->AddEntry(gr456toGen,"kTPCin + 4or5or6 ITS cls","p");
  l1->Draw();
  


  TCanvas* ceff2_6 = new TCanvas("ceff2_6","ceff2_6",cWidth,cHeight);
  ceff2_6->SetLogx();
  ceff2_6->SetGridy();
  
  gr4to6cls = new TGraphErrors(nPtBins,meanpt,n4ITSto6cls,errx,sigmn4ITSto6cls);
  gr4to6cls->SetName("mygraph04ITSto6cls");
  gr4to6cls->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr4to6cls->SetLineColor(1);
  gr4to6cls->SetLineWidth(1);
  gr4to6cls->SetMarkerColor(2);
  gr4to6cls->SetMarkerStyle(24);
  gr4to6cls->SetMarkerSize(1.);
  gr4to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4to6cls->GetYaxis()->SetTitle("ITS efficiency (normalized to kTPCin with 6 ITS cls)");
  gr4to6cls->SetMaximum(1.2);
  gr4to6cls->SetMinimum(0.);
  gr4to6cls->Draw("AP");

  gr5to6cls = new TGraphErrors(nPtBins,meanpt,n5ITSto6cls,errx,sigmn5ITSto6cls);
  gr5to6cls->SetName("mygraph05ITSto6cls");
  gr5to6cls->SetLineColor(1);
  gr5to6cls->SetLineWidth(1);
  gr5to6cls->SetMarkerColor(8);
  gr5to6cls->SetMarkerStyle(25);
  gr5to6cls->SetMarkerSize(1.);
  gr5to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr5to6cls->GetYaxis()->SetTitle("");
  gr5to6cls->Draw("P");

  gr6to6cls = new TGraphErrors(nPtBins,meanpt,n6ITSto6cls,errx,sigmn6ITSto6cls);
  gr6to6cls->SetName("mygraph06ITSto6cls");
  gr6to6cls->SetLineColor(1);
  gr6to6cls->SetLineWidth(1);
  gr6to6cls->SetMarkerColor(4);
  gr6to6cls->SetMarkerStyle(26);
  gr6to6cls->SetMarkerSize(1.);
  gr6to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr6to6cls->GetYaxis()->SetTitle("");
  gr6to6cls->Draw("P");

  gr3to6cls = new TGraphErrors(nPtBins,meanpt,n3ITSto6cls,errx,sigmn3ITSto6cls);
  gr3to6cls->SetName("mygraph03ITSto6cls");
  gr3to6cls->SetLineColor(1);
  gr3to6cls->SetLineWidth(1);
  gr3to6cls->SetMarkerColor(9);
  gr3to6cls->SetMarkerStyle(27);
  gr3to6cls->SetMarkerSize(1.);
  gr3to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr3to6cls->GetYaxis()->SetTitle("");
  gr3to6cls->Draw("P");

  gr2to6cls = new TGraphErrors(nPtBins,meanpt,n2ITSto6cls,errx,sigmn2ITSto6cls);
  gr2to6cls->SetName("mygraph02ITSto6cls");
  gr2to6cls->SetLineColor(1);
  gr2to6cls->SetLineWidth(1);
  gr2to6cls->SetMarkerColor(1);
  gr2to6cls->SetMarkerStyle(28);
  gr2to6cls->SetMarkerSize(1.);
  gr2to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr2to6cls->GetYaxis()->SetTitle("");
  gr2to6cls->Draw("P");

  gr1to6cls = new TGraphErrors(nPtBins,meanpt,n1ITSto6cls,errx,sigmn1ITSto6cls);
  gr1to6cls->SetName("mygraph01ITSto6cls");
  gr1to6cls->SetLineColor(1);
  gr1to6cls->SetLineWidth(1);
  gr1to6cls->SetMarkerColor(6);
  gr1to6cls->SetMarkerStyle(29);
  gr1to6cls->SetMarkerSize(1.);
  gr1to6cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr1to6cls->GetYaxis()->SetTitle("");
  gr1to6cls->Draw("P");


  TLegend *l2_6 = new TLegend(0.5,0.5,0.9,0.9);
  l2_6->SetFillStyle(0);
  l2_6->SetBorderSize(0);
  l2_6->AddEntry(gr6to6cls,"kTPCin + 6 ITS cls","p");
  l2_6->AddEntry(gr5to6cls,"kTPCin + 5 ITS cls","p");
  l2_6->AddEntry(gr4to6cls,"kTPCin + 4 ITS cls","p");
  l2_6->AddEntry(gr3to6cls,"kTPCin + 3 ITS cls","p");
  l2_6->AddEntry(gr2to6cls,"kTPCin + 2 ITS cls","p");
  l2_6->AddEntry(gr1to6cls,"kTPCin + 1 ITS cls","p");
  l2_6->Draw();

  TCanvas* ceff2_5 = new TCanvas("ceff2_5","ceff2_5",cWidth,cHeight);
  ceff2_5->SetLogx();
  ceff2_5->SetGridy();
  
  gr4to5cls = new TGraphErrors(nPtBins,meanpt,n4ITSto5cls,errx,sigmn4ITSto5cls);
  gr4to5cls->SetName("mygraph04ITSto5cls");
  gr4to5cls->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr4to5cls->SetLineColor(1);
  gr4to5cls->SetLineWidth(1);
  gr4to5cls->SetMarkerColor(2);
  gr4to5cls->SetMarkerStyle(24);
  gr4to5cls->SetMarkerSize(1.);
  gr4to5cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4to5cls->GetYaxis()->SetTitle("ITS efficiency (normalized to kTPCin with 5 ITS cls)");
  gr4to5cls->SetMaximum(1.2);
  gr4to5cls->SetMinimum(0.);
  gr4to5cls->Draw("AP");

  gr5to5cls = new TGraphErrors(nPtBins,meanpt,n5ITSto5cls,errx,sigmn5ITSto5cls);
  gr5to5cls->SetName("mygraph05ITSto5cls");
  gr5to5cls->SetLineColor(1);
  gr5to5cls->SetLineWidth(1);
  gr5to5cls->SetMarkerColor(1);
  gr5to5cls->SetMarkerStyle(25);
  gr5to5cls->SetMarkerSize(1.);
  gr5to5cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr5to5cls->GetYaxis()->SetTitle("");
  gr5to5cls->Draw("P");

  TLegend *l2_5 = new TLegend(0.5,0.5,0.9,0.9);
  l2_5->SetFillStyle(0);
  l2_5->SetBorderSize(0);
  l2_5->AddEntry(gr5to5cls,"kTPCin + 5 ITS cls","p");
  l2_5->AddEntry(gr4to5cls,"kTPCin + 4 ITS cls","p");
  l2_5->Draw();

  TCanvas* ceff2_4 = new TCanvas("ceff2_4","ceff2_4",cWidth,cHeight);
  ceff2_4->SetLogx();
  ceff2_4->SetGridy();
  
  gr4to4cls = new TGraphErrors(nPtBins,meanpt,n4ITSto4cls,errx,sigmn4ITSto4cls);
  gr4to4cls->SetName("mygraph04ITSto4cls");
  gr4to4cls->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr4to4cls->SetLineColor(1);
  gr4to4cls->SetLineWidth(1);
  gr4to4cls->SetMarkerColor(2);
  gr4to4cls->SetMarkerStyle(24);
  gr4to4cls->SetMarkerSize(1.);
  gr4to4cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4to4cls->GetYaxis()->SetTitle("ITS efficiency (normalized to kTPCin with 4 ITS cls)");
  gr4to4cls->SetMaximum(1.2);
  gr4to4cls->SetMinimum(0.);
  gr4to4cls->Draw("AP");

  gr3to4cls = new TGraphErrors(nPtBins,meanpt,n3ITSto4cls,errx,sigmn3ITSto4cls);
  gr3to4cls->SetName("mygraph03ITSto4cls");
  gr3to4cls->SetLineColor(1);
  gr3to4cls->SetLineWidth(1);
  gr3to4cls->SetMarkerColor(1);
  gr3to4cls->SetMarkerStyle(25);
  gr3to4cls->SetMarkerSize(1.);
  gr3to4cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr3to4cls->GetYaxis()->SetTitle("");
  gr3to4cls->Draw("P");

  TLegend *l2_4 = new TLegend(0.5,0.5,0.9,0.9);
  l2_4->SetFillStyle(0);
  l2_4->SetBorderSize(0);
  l2_4->AddEntry(gr4to4cls,"kTPCin + 4 ITS cls","p");
  l2_4->AddEntry(gr3to4cls,"kTPCin + 3 ITS cls","p");
  l2_4->Draw();

  TCanvas* ceff2_3 = new TCanvas("ceff2_3","ceff2_3",cWidth,cHeight);
  ceff2_3->SetLogx();
  ceff2_3->SetGridy();
  
  gr3to3cls = new TGraphErrors(nPtBins,meanpt,n3ITSto3cls,errx,sigmn3ITSto3cls);
  gr3to3cls->SetName("mygraph03ITSto3cls");
  gr3to3cls->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr3to3cls->SetLineColor(1);
  gr3to3cls->SetLineWidth(1);
  gr3to3cls->SetMarkerColor(2);
  gr3to3cls->SetMarkerStyle(24);
  gr3to3cls->SetMarkerSize(1.);
  gr3to3cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr3to3cls->GetYaxis()->SetTitle("ITS efficiency (normalized to kTPCin with 3 ITS cls)");
  gr3to3cls->SetMaximum(1.2);
  gr3to3cls->SetMinimum(0.);
  gr3to3cls->Draw("AP");

  gr2to3cls = new TGraphErrors(nPtBins,meanpt,n2ITSto3cls,errx,sigmn2ITSto3cls);
  gr2to3cls->SetName("mygraph02ITSto3cls");
  gr2to3cls->SetLineColor(1);
  gr2to3cls->SetLineWidth(1);
  gr2to3cls->SetMarkerColor(1);
  gr2to3cls->SetMarkerStyle(25);
  gr2to3cls->SetMarkerSize(1.);
  gr2to3cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr2to3cls->GetYaxis()->SetTitle("");
  gr2to3cls->Draw("P");

  TLegend *l2_3 = new TLegend(0.5,0.5,0.9,0.9);
  l2_3->SetFillStyle(0);
  l2_3->SetBorderSize(0);
  l2_3->AddEntry(gr3to3cls,"kTPCin + 3 ITS cls","p");
  l2_3->AddEntry(gr2to3cls,"kTPCin + 2 ITS cls","p");
  l2_3->Draw();

  TCanvas* ceff2_2 = new TCanvas("ceff2_2","ceff2_2",cWidth,cHeight);
  ceff2_2->SetLogx();
  ceff2_2->SetGridy();
  
  gr1to2cls = new TGraphErrors(nPtBins,meanpt,n1ITSto2cls,errx,sigmn1ITSto2cls);
  gr1to2cls->SetName("mygraph01ITSto2cls");
  gr1to2cls->SetTitle("ITS tracking efficiency (B=0.5T, no misal.)");
  gr1to2cls->SetLineColor(1);
  gr1to2cls->SetLineWidth(1);
  gr1to2cls->SetMarkerColor(2);
  gr1to2cls->SetMarkerStyle(24);
  gr1to2cls->SetMarkerSize(1.);
  gr1to2cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr1to2cls->GetYaxis()->SetTitle("ITS efficiency (normalized to kTPCin with 2 ITS cls)");
  gr1to2cls->SetMaximum(1.2);
  gr1to2cls->SetMinimum(0.);
  gr1to2cls->Draw("AP");

  gr2to2cls = new TGraphErrors(nPtBins,meanpt,n2ITSto2cls,errx,sigmn2ITSto2cls);
  gr2to2cls->SetName("mygraph02ITSto2cls");
  gr2to2cls->SetLineColor(1);
  gr2to2cls->SetLineWidth(1);
  gr2to2cls->SetMarkerColor(1);
  gr2to2cls->SetMarkerStyle(25);
  gr2to2cls->SetMarkerSize(1.);
  gr2to2cls->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr2to2cls->GetYaxis()->SetTitle("");
  gr2to2cls->Draw("P");

  TLegend *l2_2 = new TLegend(0.5,0.5,0.9,0.9);
  l2_2->SetFillStyle(0);
  l2_2->SetBorderSize(0);
  l2_2->AddEntry(gr2to2cls,"kTPCin + 2 ITS cls","p");
  l2_2->AddEntry(gr1to2cls,"kTPCin + 1 ITS cls","p");
  l2_2->Draw();


  TCanvas* ceff3 = new TCanvas("ceff3","ceff3",cWidth,cHeight);
  ceff3->SetLogx();
  ceff3->SetGridy();
  
  gr4ITSclsparttoTPC = new TGraphErrors(nPtBins,meanpt,n4ITSclsparttoTPC,errx,sigmn4ITSclsparttoTPC);
  gr4ITSclsparttoTPC->SetName("mygraph04ITSclsparttoTPC");
  gr4ITSclsparttoTPC->SetTitle("Number of produced clusters (B=0.5T, no misal.)");
  gr4ITSclsparttoTPC->SetLineColor(1);
  gr4ITSclsparttoTPC->SetLineWidth(1);
  gr4ITSclsparttoTPC->SetMarkerColor(2);
  gr4ITSclsparttoTPC->SetMarkerStyle(22);
  gr4ITSclsparttoTPC->SetMarkerSize(1.);
  gr4ITSclsparttoTPC->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr4ITSclsparttoTPC->GetYaxis()->SetTitle("Fraction of primary kTPCin tracks producing N clusters");
  gr4ITSclsparttoTPC->SetMaximum(1.2);
  gr4ITSclsparttoTPC->SetMinimum(0.);
  gr4ITSclsparttoTPC->Draw("AP");

  gr5ITSclsparttoTPC = new TGraphErrors(nPtBins,meanpt,n5ITSclsparttoTPC,errx,sigmn5ITSclsparttoTPC);
  gr5ITSclsparttoTPC->SetName("mygraph05ITSclsparttoTPC");
  gr5ITSclsparttoTPC->SetLineColor(1);
  gr5ITSclsparttoTPC->SetLineWidth(1);
  gr5ITSclsparttoTPC->SetMarkerColor(8);
  gr5ITSclsparttoTPC->SetMarkerStyle(22);
  gr5ITSclsparttoTPC->SetMarkerSize(1.);
  gr5ITSclsparttoTPC->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr5ITSclsparttoTPC->GetYaxis()->SetTitle("");
  gr5ITSclsparttoTPC->Draw("P");

  gr6ITSclsparttoTPC = new TGraphErrors(nPtBins,meanpt,n6ITSclsparttoTPC,errx,sigmn6ITSclsparttoTPC);
  gr6ITSclsparttoTPC->SetName("mygraph06ITSITSclsparttoTPC");
  gr6ITSclsparttoTPC->SetLineColor(1);
  gr6ITSclsparttoTPC->SetLineWidth(1);
  gr6ITSclsparttoTPC->SetMarkerColor(4);
  gr6ITSclsparttoTPC->SetMarkerStyle(22);
  gr6ITSclsparttoTPC->SetMarkerSize(1.);
  gr6ITSclsparttoTPC->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr6ITSclsparttoTPC->GetYaxis()->SetTitle("");
  gr6ITSclsparttoTPC->Draw("P");

  gr3ITSclsparttoTPC = new TGraphErrors(nPtBins,meanpt,n3ITSclsparttoTPC,errx,sigmn3ITSclsparttoTPC);
  gr3ITSclsparttoTPC->SetName("mygraph03ITSITSclsparttoTPC");
  gr3ITSclsparttoTPC->SetLineColor(1);
  gr3ITSclsparttoTPC->SetLineWidth(1);
  gr3ITSclsparttoTPC->SetMarkerColor(6);
  gr3ITSclsparttoTPC->SetMarkerStyle(22);
  gr3ITSclsparttoTPC->SetMarkerSize(1.);
  gr3ITSclsparttoTPC->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr3ITSclsparttoTPC->GetYaxis()->SetTitle("");
  gr3ITSclsparttoTPC->Draw("P");

  gr2ITSclsparttoTPC = new TGraphErrors(nPtBins,meanpt,n2ITSclsparttoTPC,errx,sigmn2ITSclsparttoTPC);
  gr2ITSclsparttoTPC->SetName("mygraph06ITSITSclsparttoTPC");
  gr2ITSclsparttoTPC->SetLineColor(1);
  gr2ITSclsparttoTPC->SetLineWidth(1);
  gr2ITSclsparttoTPC->SetMarkerColor(1);
  gr2ITSclsparttoTPC->SetMarkerStyle(22);
  gr2ITSclsparttoTPC->SetMarkerSize(1.);
  gr2ITSclsparttoTPC->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr2ITSclsparttoTPC->GetYaxis()->SetTitle("");
  gr2ITSclsparttoTPC->Draw("P");

  TLegend *l3 = new TLegend(0.5,0.5,0.9,0.9);
  l3->SetFillStyle(0);
  l3->SetBorderSize(0);
  l3->AddEntry(gr2ITSclsparttoTPC,"2 ITS cls","p");
  l3->AddEntry(gr3ITSclsparttoTPC,"3 ITS cls","p");
  l3->AddEntry(gr4ITSclsparttoTPC,"4 ITS cls","p");
  l3->AddEntry(gr5ITSclsparttoTPC,"5 ITS cls","p");
  l3->AddEntry(gr6ITSclsparttoTPC,"6 ITS cls","p");
  l3->Draw();
  

  ////////////--------------5 POINTS SCENARIO
  
  TCanvas* fitting1 = new TCanvas("fitting1","5 Points Scenario tracks analysis",cWidth,cHeight);
  fitting1->Divide(3,2);
  fitting1->cd(1);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr133 = new TGraphErrors(nPtBins,meanpt,n133,errx,sigmn133);
  gr133->SetName("mygraph133");
  gr133->SetLineColor(1);
  gr133->SetLineWidth(1);
  gr133->SetMarkerColor(4);
  gr133->SetMarkerStyle(21);
  gr133->SetMarkerSize(.5);
  titlegraph1="NTracks without SSD2 point for ";
  //titlegraph1="d_{0}(r#phi) Resolution for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using True Vertex");
  gr133->SetTitle(titlegraph1);
  gr133->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr133->GetYaxis()->SetTitle("");
  gr133->Draw("AP");
  
  fitting1->cd(4);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr233 = new TGraphErrors(nPtBins,meanpt,n233,errx,sigmn233);
  gr233->SetName("mygraph233");
  gr233->SetLineColor(1);
  gr233->SetLineWidth(1);
  gr233->SetMarkerColor(4);
  gr233->SetMarkerStyle(21);
  gr233->SetMarkerSize(.5);
  //titlegraph1="d_{0}(r#phi) Resolution for  ";
  titlegraph1="NTracks without SSD1 point for ";
  titlegraph1.Append(partforgraph);
  // titlegraph1.Append(" using Reconstructed Vertex");
  gr233->SetTitle(titlegraph1);
  gr233->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr233->GetYaxis()->SetTitle("");
  gr233->Draw("AP");
  
  fitting1->cd(2);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr313 = new TGraphErrors(nPtBins,meanpt,n313,errx,sigmn313);
  gr313->SetName("mygraph313");
  gr313->SetLineColor(1);
  gr313->SetLineWidth(1);
  gr313->SetMarkerColor(4);
  gr313->SetMarkerStyle(21);
  gr313->SetMarkerSize(.5);
  titlegraph1="NTracks without SDD2 point for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using Vertex on the Fly");
  gr313->SetTitle(titlegraph1);
  gr313->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr313->GetYaxis()->SetTitle("");
  gr313->Draw("AP");

  fitting1->cd(5);
  // gPad->SetLogy();
  gPad->SetLogx();
 gPad->SetGridx();
  gPad->SetGridy();  
  gr323 = new TGraphErrors(nPtBins,meanpt,n323,errx,sigmn323);
  gr323->SetName("mygraph323");
  gr323->SetLineColor(1);
  gr323->SetLineWidth(1);
  gr323->SetMarkerColor(4);
  gr323->SetMarkerStyle(21);
  gr323->SetMarkerSize(.5);
  titlegraph1="NTracks without SDD1 point for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using Vertex on the Fly");
  gr323->SetTitle(titlegraph1);
  gr323->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr323->GetYaxis()->SetTitle("");
  gr323->Draw("AP");

  fitting1->cd(3);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr331 = new TGraphErrors(nPtBins,meanpt,n331,errx,sigmn331);
  gr331->SetName("mygraph331");
  gr331->SetLineColor(1);
  gr331->SetLineWidth(1);
  gr331->SetMarkerColor(4);
  gr331->SetMarkerStyle(21);
  gr331->SetMarkerSize(.5);
  titlegraph1="NTracks without SPD2 point for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using Vertex on the Fly");
  gr331->SetTitle(titlegraph1);
  gr331->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr331->GetYaxis()->SetTitle("");
  gr331->Draw("AP");

  fitting1->cd(6);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr332 = new TGraphErrors(nPtBins,meanpt,n332,errx,sigmn332);
  gr332->SetName("mygraph332");
  gr332->SetLineColor(1);
  gr332->SetLineWidth(1);
  gr332->SetMarkerColor(4);
  gr332->SetMarkerStyle(21);
  gr332->SetMarkerSize(.5);
  titlegraph1="NTracks without SPD1 point for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using Vertex on the Fly");
  gr332->SetTitle(titlegraph1);
  gr332->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr332->GetYaxis()->SetTitle("");
  gr332->Draw("AP");


  TCanvas* fitting2 = new TCanvas("fitting2","SPD Scenario tracks analysis",cWidth,cHeight);
  fitting2->Divide(3,2);
  
  fitting2->cd(1);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr331->Draw("AP");

  fitting2->cd(2);
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr332->Draw("AP");

  fitting2->cd(3);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  gr330 = new TGraphErrors(nPtBins,meanpt,n330,errx,sigmn330);
  gr330->SetName("mygraph330");
  gr330->SetLineColor(1);
  gr330->SetLineWidth(1);
  gr330->SetMarkerColor(4);
  gr330->SetMarkerStyle(21);
  gr330->SetMarkerSize(.5);
  titlegraph1="NTracks without 2 SPD points for ";
  //titlegraph1="d_{0}(r#phi) Resolution for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using True Vertex");
  gr330->SetTitle(titlegraph1);
  gr330->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  gr330->GetYaxis()->SetTitle("");
  gr330->Draw("AP");

  fitting2->cd(4);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  grSSD = new TGraphErrors(nPtBins,meanpt,nSSD,errx,sigmnSSD);
  grSSD->SetName("mygraphSSD");
  grSSD->SetLineColor(1);
  grSSD->SetLineWidth(1);
  grSSD->SetMarkerColor(4);
  grSSD->SetMarkerStyle(21);
  grSSD->SetMarkerSize(.5);
  titlegraph1="NTracks without one SSD point for ";
  //titlegraph1="d_{0}(r#phi) Resolution for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using True Vertex");
  grSSD->SetTitle(titlegraph1);
  grSSD->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  grSSD->GetYaxis()->SetTitle("");
  grSSD->Draw("AP");

  fitting2->cd(5);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  grSDD = new TGraphErrors(nPtBins,meanpt,nSDD,errx,sigmnSDD);
  grSDD->SetName("mygraphSDD");
  grSDD->SetLineColor(1);
  grSDD->SetLineWidth(1);
  grSDD->SetMarkerColor(4);
  grSDD->SetMarkerStyle(21);
  grSDD->SetMarkerSize(.5);
  titlegraph1="NTracks without one SDD point for ";
  //titlegraph1="d_{0}(r#phi) Resolution for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using True Vertex");
  grSDD->SetTitle(titlegraph1);
  grSDD->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  grSDD->GetYaxis()->SetTitle("");
  grSDD->Draw("AP");

  fitting2->cd(6);
  // gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridx();
  gPad->SetGridy();  
  grOTHERS = new TGraphErrors(nPtBins,meanpt,nOTHERS,errx,sigmnOTHERS);
  grOTHERS->SetName("mygraphOTHERS");
  grOTHERS->SetLineColor(1);
  grOTHERS->SetLineWidth(1);
  grOTHERS->SetMarkerColor(4);
  grOTHERS->SetMarkerStyle(21);
  grOTHERS->SetMarkerSize(.5);
  titlegraph1="NTracks with 4 points (2 SPD+ ?)for ";
  //titlegraph1="d_{0}(r#phi) Resolution for ";
  titlegraph1.Append(partforgraph);
  //    titlegraph1.Append(" using True Vertex");
  grOTHERS->SetTitle(titlegraph1);
  grOTHERS->GetXaxis()->SetTitle("p_{t} [GeV/c] ");
  grOTHERS->GetYaxis()->SetTitle("");
  grOTHERS->Draw("AP");

  ceff0->cd();
  gr330->Draw("p");
  gr331->Draw("p");
  gr332->Draw("p");
  grSSD->Draw("p");
  grSDD->Draw("p");
  //  gr313->Draw("p");
  // gr323->Draw("p");
  //gr133->Draw("p");
  //gr233->Draw("p");
  

  //---------CLOSE EVRYTHING AND SAVE------------
  TFile *outfile = new TFile("EfficienciesAnalysis.root","recreate");
  gr4->Write();
  gr5->Write();
  gr6->Write();
  gr456->Write();
  gr4toGen->Write();
  gr5toGen->Write();
  gr6toGen->Write();
  gr456toGen->Write();
  grTPCtoGen->Write();
  gr4to6cls->Write();
  gr5to6cls->Write();
  gr6to6cls->Write();
  //gr456to6cls->Write();
  //gr330->Write();
  //gr331->Write();
  //gr332->Write();
  //gr323->Write();
  //gr313->Write();
  //gr133->Write();
  //gr233->Write();
  //grOTHERS->Write();
  //grSPD->Write();
  //grSSD->Write();
  //grSDD->Write();
  outfile->Close();


  printf("Number of fakes: %d\n",nfakes);
  return;
}
//--------------------------------------------------------------------------
Int_t ITSnCluster(Int_t ITSf) 
{
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
    cout<<"Wrong ITSflag assignment!"<<endl;
    return 99;
  }
  
  return ITSsel*sign;
}
//----------------------------------------------------------------------------
