void MUONTestAbso (Int_t evNumber1=0,Int_t evNumber2=0) 
{
  // Useful to test the absorber correction in the reconstruction procedure
  // (mean energy loss + Branson correction)
  // The AliMUONTrackParam class from the reconstruction is directly checked
  // in this macro using the GEANT particle momentum downstream of the absorber.
  // Histograms are saved in file MUONTestAbso.root.
  // Use DrawTestAbso.C to display control histograms.

    
  // Dynamically link some shared libs                    
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Connect the Root Galice file containing Geometry, Kine and Hits

  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");

    
  if (!file) {
    printf("\n Creating galice.root \n");
    file = new TFile("galice.root");
  } else {
    printf("\n galice.root found in file list");
  }
    
// Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)(file->Get("gAlice"));
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) {
      printf("\n Create new gAlice object");
      gAlice = new AliRun("gAlice","Alice test program");
    }
  }

    
  //  Create some histograms

  TH1F *hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2) around Upsilon", 120, 8., 11.);
  TH1F *hDeltaP1 = new TH1F("hDeltaP1", " delta P for muons r < rLimit", 100, -10., 10.);
  TH1F *hDeltaAngle1 = new TH1F("hDeltaAngle1", " delta angle for muons r < rLimit", 100, 0., 0.5);
  //
  TH1F *hDeltaP2 = new TH1F("hDeltaP2", " delta P for muons r > rLimit", 100, -10., 10.);
  TH1F *hDeltaAngle2 = new TH1F("hDeltaAngle2", " delta angle for muons r > rLimit", 100, 0., 0.5);
  
  // Define constant parameters
  Float_t massMuon = 0.105658389;  // muon mass
  Float_t massUpsilon = 9.46037;   // Upsilon mass
  Double_t PI = 3.141592654;
  Double_t zEndAbso = 503;         // z position of the end of the absorber 
  Double_t rLimit = zEndAbso * TMath::Tan(3*PI/180); // 3 deg. limit (different absorber composition)
  Double_t printLevel = 0;

  if (printLevel >= 1)
  cout <<" **** z end Absorber : "<< zEndAbso <<" rLimit : "<< rLimit << endl;

  Double_t pX[2],pY[2],pZ[2],pTot[2],theta[2],radius[2]; // extrapolated parameters of particles from chamber 1 to vertex 
   
  for (Int_t nev=0; nev<= evNumber2; nev++) {  // loop over events
    if (printLevel >= 0) {
      cout << "  "<<endl;
      cout << " **** Event # " << nev <<endl;
    }
    
    Int_t nparticles = gAlice->GetEvent(nev);
    if (printLevel >= 1) {
      cout << " Total number of nparticles  " << nparticles <<endl;
    }
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return;
     
    Double_t PxGen[2],PyGen[2],PzGen[2],PTotGen[2],ThetaGen[2]; // parameters of generated  particles at the vertex
     
    for (int iPart = 0; iPart < nparticles  ; iPart++) {  // loop over particles
      
      TParticle *particle;
      particle = gAlice->Particle(iPart);    

      if ( particle->GetPdgCode() == kMuonPlus ) {
	PxGen[0] = particle->Px();
	PyGen[0] = particle->Py();
	PzGen[0] = particle->Pz();
	PTotGen[0] = TMath::Sqrt(PxGen[0]*PxGen[0]+PyGen[0]*PyGen[0]+PzGen[0]*PzGen[0]);
	ThetaGen[0] = TMath::ATan2(TMath::Sqrt(PxGen[0]*PxGen[0]+PyGen[0]*PyGen[0]),PzGen[0])*180/PI;
	if (printLevel >= 1) {
	  cout << " Particle id : "<<particle->GetPdgCode()<<" px,py,pz,theta : "<<PxGen[0]<<" "<<PyGen[0]<<" "<<PzGen[0]<<" "<<ThetaGen[0]<<endl;
	}
      }
      if ( particle->GetPdgCode() == kMuonMinus ) {
	PxGen[1] = particle->Px();
	PyGen[1] = particle->Py();
	PzGen[1] = particle->Pz();
	PTotGen[1] = TMath::Sqrt(PxGen[1]*PxGen[1]+PyGen[1]*PyGen[1]+PzGen[1]*PzGen[1]);
	ThetaGen[1] = TMath::ATan2(TMath::Sqrt(PxGen[1]*PxGen[1]+PyGen[1]*PyGen[1]),PzGen[1])*180/PI;
	if (printLevel >= 1)
	  cout << " Particle #: "<<particle->GetPdgCode()<<" px,py,pz,theta : "<<PxGen[1]<<" "<<PyGen[1]<<" "<<PzGen[1]<<" "<<ThetaGen[1]<<endl;
      }
    }


    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
     
     
    TTree *TH = gAlice->TreeH();
    Int_t ntracks =(Int_t) TH->GetEntries();
     
    for (Int_t i = 0; i < 2; i++){  
      pX[i] = pY[i] = pZ[i] = pTot[i] = theta[i] = radius[i] = 0;
    }
    
    for (Int_t track=0; track<ntracks;track++) { // loop over tracks
      gAlice->ResetHits();
      Int_t nbytes += TH->GetEvent(track);
      if (MUON)  {
	for(AliMUONHit* mHit=(AliMUONHit*)MUON->FirstHit(-1); 
	    mHit;
	    mHit=(AliMUONHit*)MUON->NextHit())  // loop over GEANT muon hits
	  {
	    Int_t   nch   = mHit->Chamber();  // chamber number
	     
	    if (nch != 1) continue;  
	    
	    Int_t ipart = mHit->Particle();
	    Double_t x     = mHit->X();        // x-pos of hit in chamber #1
	    Double_t y     = mHit->Y();        // y-pos of hit in chamber #1
	    Double_t z     = mHit->Z();        // z-pos of hit in chamber #1
	     
	     
	    if (ipart == kMuonPlus || ipart == kMuonMinus) {
	      Double_t px0=mHit->Px();
	      Double_t py0=mHit->Py();		   
	      Double_t pz0=mHit->Pz();
	      Double_t nonBendingSlope=px0/pz0;
	      Double_t bendingSlope=py0/pz0;

	      // create an AliMUONTrackParam object in chamber #1
	      AliMUONTrackParam *trackParam = new AliMUONTrackParam();
	      Double_t bendingMometum = TMath::Sqrt(pz0*pz0+px0*px0);
	      Double_t inverseBendingMomentum = 1/bendingMometum;
	       
	      trackParam->SetBendingCoor(y);
	      trackParam->SetNonBendingCoor(x);
	      trackParam->SetInverseBendingMomentum(inverseBendingMomentum);
	      trackParam->SetBendingSlope(bendingSlope);
	      trackParam->SetNonBendingSlope(nonBendingSlope);
	      trackParam->SetZ(z);	
	       
	      if (printLevel >= 1) {
		cout <<" **** before extrap " <<endl;
		cout << " px, py, pz = " <<px0<<" "<<py0<<" "<<pz0<<endl;
		trackParam->Dump();
	      }
	      trackParam->ExtrapToVertex(); // extrapolate track parameters through the absorber
	      if (printLevel >= 1) {
		cout <<" **** after extrap " <<endl;
		trackParam->Dump();
	      }

	      Int_t iMuon;
	      if (ipart == kMuonPlus) iMuon = 0;
	      else iMuon = 1;
	       
	      // calculate track parameters extrapolated to vertex.
	      bendingSlope = trackParam->GetBendingSlope();
	      nonBendingSlope = trackParam->GetNonBendingSlope();
	      Double_t pYZ = 1/TMath::Abs(trackParam->GetInverseBendingMomentum());
	      pZ[iMuon] = pYZ/TMath::Sqrt(1+bendingSlope*bendingSlope);
	      pX[iMuon] = pZ[iMuon] * nonBendingSlope;
	      pY[iMuon] = pZ[iMuon] * bendingSlope;
	      pTot[iMuon] = TMath::Sqrt(pYZ*pYZ+pX[iMuon]*pX[iMuon]);
	      theta[iMuon] = TMath::ATan2(TMath::Sqrt(pX[iMuon]*pX[iMuon]+pY[iMuon]*pY[iMuon]),pZ[iMuon])*180/PI;
	      radius[iMuon] = TMath::Sqrt(x*x+y*y);
	      
	      if (printLevel >= 1)
		cout << " px, py, pz, theta, radius = " <<pX[iMuon]<<" "<<pY[iMuon]<<" "<<pZ[iMuon]<<" "<<theta[iMuon]<<" "<<radius[iMuon]<<endl;
	      delete trackParam;
	      
	    } // if MuonPlus or MuonMinus
	  } // loop over MUON hits 
      } // if MUON 
    } // loop over tracks	   

    if (pX[0] != 0 && pX[1] != 0) { // if track 1 & 2 exist 
      Float_t massResonance = MuPlusMuMinusMass(pX[0],pY[0],pZ[0],pX[1],pY[1],pZ[1],massMuon);
      hInvMassRes->Fill(massResonance);
      for (Int_t i = 0; i<2; i++){
	Double_t cosAngle = pX[i]*PxGen[i]+pY[i]*PyGen[i]+pZ[i]*PzGen[i];
	cosAngle = cosAngle/(pTot[i]*PTotGen[i]);
	Double_t deltaAngle = TMath::ACos(cosAngle)*180/PI;
	if (radius[i] < rLimit) {
	  hDeltaP1->Fill(pTot[i]-PTotGen[i]);
	  hDeltaAngle1->Fill(deltaAngle);
	}else{
	  hDeltaP2->Fill(pTot[i]-PTotGen[i]);
	  hDeltaAngle2->Fill(deltaAngle);
	}
      }
    }
  } // loop over event

  // save histograms in MUONTestAbso.root
  TFile *histoFile = new TFile("MUONTestAbso.root", "RECREATE");
  hInvMassRes->Write();
  hDeltaP1->Write();
  hDeltaAngle1->Write();
  hDeltaP2->Write();
  hDeltaAngle2->Write();
  histoFile->Close();
  
}


Float_t MuPlusMuMinusMass(Float_t Px1, Float_t Py1, Float_t Pz1, Float_t Px2, Float_t Py2, Float_t Pz2, Float_t MuonMass)
{
  // return invariant mass for particle 1 & 2 whose masses are equal to MuonMass

  Float_t e1 = TMath::Sqrt(MuonMass * MuonMass + Px1 * Px1 + Py1 * Py1 + Pz1 * Pz1);
  Float_t e2 = TMath::Sqrt(MuonMass * MuonMass + Px2 * Px2 + Py2 * Py2 + Pz2 * Pz2);
  return (TMath::Sqrt(2.0 * (MuonMass * MuonMass + e1 * e2 - Px1 * Px2 - Py1 * Py2 - Pz1 * Pz2)));
}
