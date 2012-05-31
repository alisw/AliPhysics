void tzero15tracks(const char* datafile, Float_t globalTimeRes=1.2e-10)
{
  
  //
  // Calculate the t0 distribution 
  // selecting 15 tracks
  // (primary particles reaching TOF)
  // NB: see the description of the analogous AliTOFT0 class
  // NB II: apply this macro to a realistic Pb-Pb event (e.g. Hijing event)
  // Use case
  // - start root
  // root [0] .L tzero15tracksopt.C
  // // exec the macro with the default time resolution 120 ps
  // root [1] tzero15tracksopt("hij25evCentralTOFv4T0Set10.04T.root")

  // // exec the macro with 150 ps global time resolution
  // // NB: time resolution has to be given in [s]
  // root [1] tzero15tracksopt("hij25evCentralTOFv4T0Set10.04T.root", 1.5e-10))


  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  Int_t ngood;
  Int_t nmisidentified=0;
  Int_t nmisidentified0=0;
  Int_t nmisidentified1=0;
  Int_t nmisidentified2=0;
  Int_t nmisidentified3=0;
  Int_t nmisidentified4=0;
  Int_t nmisidentified5=0;
  Int_t nmisidentified6=0;
  Int_t nmisidentified7=0;
  Int_t nmisidentified8=0;
  Int_t nmisidentified9=0;
  Int_t nmisidentified10=0;
  Int_t nmisidentified11=0;
  Int_t nmisidentified12=0;
  Int_t nmisidentified13=0;
  Int_t nmisidentified14=0;
  Int_t nbytes = 0;
  Int_t ipartold = -1;
  Int_t dblehit = 0;
  Int_t j,hit,ipart;
  Int_t nhits;
  Int_t selected=0;
  Int_t istop=0;
  const Float_t timeresolution=globalTimeRes; // in [s] 
  Float_t timeresolutioninns=timeresolution*(1.e+9); // convert in [ns]
  cout << "Global Time Resolution " << timeresolutioninns << " ns" << endl; //

  // output (histos!) filename
  char outFileName[100];
  strcpy(outFileName,"ht015tr");
  strcat(outFileName,datafile);

  gBenchmark->Start("TOFT0");
  
  TH1F* htzerobest= new TH1F("htzerobest","T0 for best assignment",200,-1.,1.);
  TH1F* hchibest  = new TH1F("hchibest","ChiSquare Min Distribution",80,0.,40.);
  TH1F* hchibestconflevel  = new TH1F("hchibestconflevel","ChiSquare Min Confidence Level",10,0.,1.);
  
  Float_t t0best=999.;
  Int_t assparticle[15]={3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  Int_t truparticle[15]={3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
  Float_t timeofflight[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t momentum[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t timezero[15];
  Float_t weightedtimezero[15];
  Float_t beta[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t sqMomError[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t sqTrackError[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t massarray[3]={0.13957,0.493677,0.9382723};
  Float_t chisquare=999.;
  Float_t tracktoflen[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  TRandom* rnd = new TRandom();
  TParticle *particle;
  AliTOFhitT0  *tofHit;
  
  TFile *file = TFile::Open(datafile,"old");
  
  // Get AliRun object from file or create it if not on file
  gAlice = (AliRun*)file->Get("gAlice");
  
  // Import the Kine and Hits Trees for the event evNumber in the file

  Int_t nallevents = (Int_t) gAlice->TreeE()->GetEntries();
  
  for (Int_t nev=0; nev<nallevents;nev++) {
    
    Int_t nparticles = gAlice->GetEvent(nev);
    if (nparticles <= 0) return;
    //     ngood=0;
    // Get pointers to Alice detectors and Hits containers
    
    AliDetector *TOF  = gAlice->GetDetector("TOF");
    TObjArray *Particles = gAlice->Particles();
    
    if (TOF) TClonesArray *TOFhits   = TOF->Hits();
    
    TTree *TH = gAlice->TreeH();
    Int_t ntracks = TH->GetEntries();
    
    
    // Start loop on primary tracks in the hits containers
    Int_t lasttrack=-1;
    Int_t nset=0;
    for (Int_t track=lasttrack+1; track<ntracks;track++) {
      
      if(nset>=3) break; // check on the number of set analyzed
      
      gAlice->ResetHits();
      nbytes += TH->GetEvent(track);
      
      // =======>Histogram TOF
      if (TOF) {
	nhits = TOFhits->GetEntriesFast();
	
	for (hit=0;hit<nhits;hit++) {
	  
	  tofHit   = (AliTOFhitT0*)TOFhits->UncheckedAt(hit);
	  ipart    = tofHit->GetTrack();
	  
	  // check to discard the case when the same particle is selected more than one
	  // time 
	  
	  if (ipart != ipartold){
	    
	    particle = (TParticle*)gAlice->Particle(ipart);
	    
	    Float_t idealtime=tofHit->GetTof();
	    //	   Float_t time=idealtime;
	    Float_t time   = rnd->Gaus(idealtime, timeresolution);
	    Float_t toflen=tofHit->GetLen();
	    toflen=toflen/100.; // toflen given in m because GEANT gives cm!
	    Int_t pdg   = particle->GetPdgCode();
	    Int_t abspdg   =TMath::Abs(pdg);
	    Float_t idealmom  = particle->P();
	    Float_t momres=idealmom*0.025; // 2.5% res token into account for all momenta
	    //Float_t mom =idealmom+2.*(rnd->Rndm()-0.5)*momres;
	    Float_t mom =rnd->Gaus(idealmom,momres);
	    Bool_t isgoodpart=(abspdg==211 || abspdg==2212 || abspdg==321);
	    
	    time*=1.e+9; // tof given in nanoseconds	   
	    if (particle->GetFirstMother() < 0 && isgoodpart && mom<=1.75 && mom>=1.25){
	      
	      selected+=1;
	      istop=selected;
	      if(istop>15) break;
	      Int_t index=selected-1;
	      timeofflight[index]=time;
	      tracktoflen[index]=toflen;
	      momentum[index]=mom;
	      switch (abspdg) {
	      case 211:
		truparticle[index]=0;
		break ;
	      case 321:
		truparticle[index]=1;
		break ;
	      case 2212: 
		truparticle[index]=2;
		break ;
	      }
	      
	    }
	    ipartold = ipart;
	    
	    if(istop==15){ // start analysis on current set
	      nset+=1;
	      lasttrack=track;
	      istop=0;
	      selected=0;
	      cout << "starting t0 calculation for current set" << endl;
	      for (Int_t i1=0; i1<3;i1++) {
		beta[0]=momentum[0]/TMath::Sqrt(massarray[i1]*massarray[i1]+momentum[0]*momentum[0]);
		for (Int_t i2=0; i2<3;i2++) { 
		  beta[1]=momentum[1]/TMath::Sqrt(massarray[i2]*massarray[i2]+momentum[1]*momentum[1]);
		  for (Int_t i3=0; i3<3;i3++) {
		    beta[2]=momentum[2]/TMath::Sqrt(massarray[i3]*massarray[i3]+momentum[2]*momentum[2]);
		    for (Int_t i4=0; i4<3;i4++) {
		      beta[3]=momentum[3]/TMath::Sqrt(massarray[i4]*massarray[i4]+momentum[3]*momentum[3]);
		      for (Int_t i5=0; i5<3;i5++) {
			beta[4]=momentum[4]/TMath::Sqrt(massarray[i5]*massarray[i5]+momentum[4]*momentum[4]);
			for (Int_t i6=0; i6<3;i6++) {
			  beta[5]=momentum[5]/TMath::Sqrt(massarray[i6]*massarray[i6]+momentum[5]*momentum[5]);
			  for (Int_t i7=0; i7<3;i7++) { 
			    beta[6]=momentum[6]/TMath::Sqrt(massarray[i7]*massarray[i7]+momentum[6]*momentum[6]);
			    for (Int_t i8=0; i8<3;i8++) {
			      beta[7]=momentum[7]/TMath::Sqrt(massarray[i8]*massarray[i8]+momentum[7]*momentum[7]);
			      for (Int_t i9=0; i9<3;i9++) {
				beta[8]=momentum[8]/TMath::Sqrt(massarray[i9]*massarray[i9]+momentum[8]*momentum[8]);
				for (Int_t i10=0; i10<3;i10++) { 	
				  beta[9]=momentum[9]/TMath::Sqrt(massarray[i10]*massarray[i10]+momentum[9]*momentum[9]);
				  for (Int_t i11=0; i11<3;i11++) {
				    beta[10]=momentum[10]/TMath::Sqrt(massarray[i11]*massarray[i11]+momentum[10]*momentum[10]);
				    for (Int_t i12=0; i12<3;i12++) { 
				      beta[11]=momentum[11]/TMath::Sqrt(massarray[i12]*massarray[i12]+momentum[11]*momentum[11]);
				      for (Int_t i13=0; i13<3;i13++) {
					beta[12]=momentum[12]/TMath::Sqrt(massarray[i13]*massarray[i13]+momentum[12]*momentum[12]);
					for (Int_t i14=0; i14<3;i14++) {
					  beta[13]=momentum[13]/TMath::Sqrt(massarray[i14]*massarray[i14]+momentum[13]*momentum[13]);
					  for (Int_t i15=0; i15<3;i15++) {
					    beta[14]=momentum[14]/TMath::Sqrt(massarray[i15]*massarray[i15]+momentum[14]*momentum[14]);
				    
					    Float_t meantzero=0.;
					    Float_t sumAllweights=0.;
					    for (Int_t itz=0; itz<15;itz++) {
					      sqMomError[itz]=((1.-beta[itz]*beta[itz])*0.025)*((1.-beta[itz]*beta[itz])*0.025)*(tracktoflen[itz]/(0.299792*beta[itz]))*(tracktoflen[itz]/(0.299792*beta[itz])); // this gives the square of the momentum error in nanoseconds
					      sqTrackError[itz]=(timeresolutioninns*timeresolutioninns+sqMomError[itz]); // total error for the current track
					      sumAllweights+=1./sqTrackError[itz];
					      timezero[itz]=(tracktoflen[itz]/(beta[itz]*0.299792))-timeofflight[itz];
					      weightedtimezero[itz]=((tracktoflen[itz]/(beta[itz]*0.299792))-timeofflight[itz])/sqTrackError[itz];// weighted time zero for current track
					      meantzero+=weightedtimezero[itz];
					    }
					    
					    meantzero=meantzero/sumAllweights; // it is given in [ns]
					    
					    Float_t dummychisquare=0.;
					    // calculate the chisquare for the current assignment
					    for (Int_t itz=0; itz<15;itz++) {
					      dummychisquare+=(timezero[itz]-meantzero)*(timezero[itz]-meantzero)/sqTrackError[itz];
					    }
					    
					    
					    if(dummychisquare<=chisquare){
					      assparticle[0]=i1;
					      assparticle[1]=i2;
					      assparticle[2]=i3;
					      assparticle[3]=i4;
					      assparticle[4]=i5;
					      assparticle[5]=i6;
					      assparticle[6]=i7;
					      assparticle[7]=i8;
					      assparticle[8]=i9;
					      assparticle[9]=i10;
					      assparticle[10]=i11;
					      assparticle[11]=i12;
					      assparticle[12]=i13;
					      assparticle[13]=i14;
					      assparticle[14]=i15;					      
					      chisquare=dummychisquare;
					      t0best=meantzero;
					    }
					  }
					}
				      }
				    }
				  }					  
				}
			      }
			    }
			  }
			}	
			
		      }
		    }
		  }
		}
	      }	
	      
	      cout << "true" << truparticle[0] << truparticle[1] << truparticle[2] << truparticle[3] << truparticle[4] << truparticle[5] << truparticle[6] << truparticle[7] << truparticle[8] << truparticle[9] <<truparticle[10] << truparticle[11] << truparticle[1
2] << truparticle[13] << truparticle[14] <<endl;
	      cout << "best" << assparticle[0] << assparticle[1] << assparticle[2] << assparticle[3] << assparticle[4] << assparticle[5] << assparticle[6] << assparticle[7] << assparticle[8] << assparticle[9] <<assparticle[10] << assparticle[11] << assparticle[1
2] << assparticle[13] << assparticle[14] << endl;
	      if(truparticle[0]==assparticle[0] && truparticle[1]==assparticle[1] && truparticle[2]==assparticle[2]  && truparticle[3]==assparticle[3] && truparticle[4]==assparticle[4]&& truparticle[5]==assparticle[5] && truparticle[6]==assparticle[6] && trupart
icle[7]==assparticle[7]  && truparticle[8]==assparticle[8] && truparticle[9]==assparticle[9] &&truparticle[10]==assparticle[10] && truparticle[11]==assparticle[11] && truparticle[12]==assparticle[12]  && truparticle[13]==assparticle[13] && truparticle[14]
==assparticle[14]) ngood+=1;
	      if(truparticle[0]!=assparticle[0]) nmisidentified0+=1;
	      if(truparticle[1]!=assparticle[1]) nmisidentified1+=1;
	      if(truparticle[2]!=assparticle[2]) nmisidentified2+=1;
	      if(truparticle[3]!=assparticle[3]) nmisidentified3+=1;
	      if(truparticle[4]!=assparticle[4]) nmisidentified4+=1;
	      if(truparticle[5]!=assparticle[5]) nmisidentified5+=1;
	      if(truparticle[6]!=assparticle[6]) nmisidentified6+=1;
	      if(truparticle[7]!=assparticle[7]) nmisidentified7+=1;
	      if(truparticle[8]!=assparticle[8]) nmisidentified8+=1;
	      if(truparticle[9]!=assparticle[9]) nmisidentified9+=1;
	      if(truparticle[10]!=assparticle[10]) nmisidentified10+=1;
	      if(truparticle[11]!=assparticle[11]) nmisidentified11+=1;
	      if(truparticle[12]!=assparticle[12]) nmisidentified12+=1;
	      if(truparticle[13]!=assparticle[13]) nmisidentified13+=1;
	      if(truparticle[14]!=assparticle[14]) nmisidentified14+=1;
	      cout << "chisquare for current set" << chisquare << endl;
	      htzerobest->Fill(t0best);
	      hchibest->Fill(chisquare);
	      Double_t dblechisquare=(Double_t)chisquare;
	      Float_t confLevel=(Float_t)TMath::Prob(dblechisquare,14);
	      cout << "conf level " << confLevel << endl;
	      hchibestconflevel->Fill(confLevel);
	      chisquare=999.;
	      t0best=999.;
	      
	    } // end for the current set. close if(istop==15)
	    
	    
	  } // end condition on ipartold
	  
	  
	}// end loop on hits for the current track
      }// end TOF
      if(istop>=15) break;
    }// end loop on Tracks
    cout << "ngood for current event " << ngood << endl;
    
    
  }//end loop on events
  
  nmisidentified=(nmisidentified0+nmisidentified1+nmisidentified2+nmisidentified3+nmisidentified4+nmisidentified5+nmisidentified6+nmisidentified7+nmisidentified8+nmisidentified9+nmisidentified10+nmisidentified11+nmisidentified12+nmisidentified13+nmisident
ified14);
  cout << "total misidentified " << nmisidentified << endl;
  TFile *houtfile = new TFile(outFileName,"recreate");
  houtfile->cd();
  
  htzerobest->Write();
  hchibest->Write();
  hchibestconflevel->Write();
  houtfile->Close();  

  gBenchmark->Stop("TOFT0");
  cout << "15 tracks:" << endl ;
  cout << "   took " << gBenchmark->GetCpuTime("TOFT0") << " seconds in order to calculate T0 " 
       <<  gBenchmark->GetCpuTime("TOFT0")/nallevents << " seconds per event " << endl ;
  cout << endl ;

}


