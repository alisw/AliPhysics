void MUONTestAbso (Int_t evNumber1=0,Int_t evNumber2=0,Int_t testSingle=0) 
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

  TH1F *hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2)", 1200, 0., 12.);
//   TH1F *hInvMassRes = new TH1F("hInvMassRes", "Mu+Mu- invariant mass (GeV/c2)", 100, 9., 10.5);
  TH1F *hDeltaP1 = new TH1F("hDeltaP1", " delta P for muons theta < 3 deg", 100, -10., 10.);
  TH1F *hDeltaAngle1 = new TH1F("hDeltaAngle1", " delta angle for muons theta < 3 deg", 100, 0., 0.5);
  //
  TH1F *hDeltaP2 = new TH1F("hDeltaP2", " delta P for muons theta > 3 deg", 100, -10., 10.);
  TH1F *hDeltaAngle2 = new TH1F("hDeltaAngle2", " delta angle for muons theta > 3 deg", 100, 0., 0.5);

  //
  TH1F *hRap = new TH1F("hRap"," Muon Rapidity gen.",20,2.5,4);
  //

  Int_t nBinProf=19;
  Float_t xProfMin=5;
  Float_t xProfMax=290;
  char titleHist[50];
  char numberHist[10];
  TH1F *DeltaP1[100],*DeltaP2[100];
  TH1F *DeltaAngleX1[100],*DeltaAngleX2[100];
  TH1F *DeltaAngleY1[100],*DeltaAngleY2[100];
  
  TH1F *hP1 = new TH1F("hP1"," P  theta < 3 deg ",nBinProf,xProfMin,xProfMax);
  TProfile *hProfDeltaPvsP1 = new TProfile("hProfDeltaPvsP1"," delta P  vs  P  theta < 3 deg ",nBinProf,xProfMin,xProfMax,-50,50,"s");
  TH2F *h2DeltaPvsP1 = new TH2F("h2DeltaPvsP1"," delta P  vs  P  theta < 3 deg",nBinProf,xProfMin,xProfMax,50,-50,50);
   TH1F *hSigmaPvsP1 = new TH1F("hSigmaPvsP1"," deltaP/P vs P  theta < 3 deg",nBinProf,xProfMin,xProfMax);
   for (Int_t iHist=0; iHist < nBinProf; iHist++){
      sprintf(titleHist,"%s","hDelta1P");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaP1[iHist] = new TH1F(titleHist," deltaP  theta < 3 deg ",400,-50,50);

      sprintf(titleHist,"%s","hDelta1AngleX");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaAngleX1[iHist] = new TH1F(titleHist," delta Angle X theta < 3 deg ",400,-0.05,0.05);
    
      sprintf(titleHist,"%s","hDelta1AngleY");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaAngleY1[iHist] = new TH1F(titleHist," delta Angle Y theta < 3 deg ",400,-0.05,0.05);
   }   
   TH1F *hSigGausPvsP1 = new TH1F("hSigGausPvsP1"," deltaP/P vs P gauss theta < 3 deg",nBinProf,xProfMin,xProfMax);
   TH1F *hSigGausAngleXvsP1 = new TH1F("hSigGausAngleXvsP1"," delta theta X vs P gauss theta < 3 deg",nBinProf,xProfMin,xProfMax);
   TH1F *hSigGausAngleYvsP1 = new TH1F("hSigGausAngleYvsP1"," delta theta Y vs P gauss theta < 3 deg",nBinProf,xProfMin,xProfMax);
 
  TH1F *hP2 = new TH1F("hP2"," P  theta > 3 deg ",nBinProf,xProfMin,xProfMax);
  TProfile *hProfDeltaPvsP2 = new TProfile("hProfDeltaPvsP2"," delta P  vs  P  theta > 3 deg",nBinProf,xProfMin,xProfMax,-50,50,"s");
  TH2F *h2DeltaPvsP2 = new TH2F("h2DeltaPvsP2"," delta P  vs  P  theta > 3 deg",nBinProf,xProfMin,xProfMax,50,-50,50);
   TH1F *hSigmaPvsP2 = new TH1F("hSigmaPvsP2"," deltaP/P vs P theta > 3 deg",nBinProf,xProfMin,xProfMax);
   for (Int_t iHist=0; iHist < nBinProf; iHist++){
      sprintf(titleHist,"%s","hDelta2P");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaP2[iHist] = new TH1F(titleHist," deltaP  theta > 3 deg ",400,-50,50);
    
      sprintf(titleHist,"%s","hDelta2AngleX");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaAngleX2[iHist] = new TH1F(titleHist," delta Angle X theta > 3 deg ",400,-0.05,0.05);

      sprintf(titleHist,"%s","hDelta2AngleY");
      sprintf(numberHist,"%d",iHist+1);
      strcat(titleHist,numberHist);
      DeltaAngleY2[iHist] = new TH1F(titleHist," delta Angle Y theta > 3 deg ",400,-0.05,0.05);
   }   
   TH1F *hSigGausPvsP2 = new TH1F("hSigGausPvsP2"," deltaP/P vs P gauss theta > 3 deg",nBinProf,xProfMin,xProfMax);
   
  TH1F *hSigGausAngleXvsP2 = new TH1F("hSigGausAngleXvsP2"," delta theta X vs P gauss theta > 3 deg",nBinProf,xProfMin,xProfMax);
   TH1F *hSigGausAngleYvsP2 = new TH1F("hSigGausAngleYvsP2"," delta theta Y vs P gauss theta > 3 deg",nBinProf,xProfMin,xProfMax);

   TProfile *hDeltaPxvsPhi = new TProfile("hDeltaPxvsPhi"," delta Px  vs  Phi",45,-180,180,-1,1,"s");
  TProfile *hDeltaPyvsPhi = new TProfile("hDeltaPyvsPhi"," delta Py  vs  Phi",45,-180,180,-1,1,"s");
  TProfile *hDeltaPhivsPz = new TProfile("hDeltaPhivsPz"," delta phi  vs pZ",25,0,100,-4,4);

  // Define constant parameters
  Float_t massMuon = 0.105658389;  // muon mass
  Float_t massUpsilon = 9.46037;   // Upsilon mass
  Double_t zEndAbso = 503;         // z position of the end of the absorber 
  Double_t rLimit = zEndAbso * TMath::Tan(3*TMath::Pi()/180); // 3 deg. limit (different absorber composition)
  Double_t printLevel = 0;

  if (printLevel >= 1)
  cout <<" **** z end Absorber : "<< zEndAbso <<" rLimit : "<< rLimit << endl;

  Double_t pX[2],pY[2],pZ[2],pTot[2],theta[2],radius[2]; // extrapolated parameters of particles from chamber 1 to vertex 
   
  for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {  // loop over events
    if (printLevel >= 0 && nev%100 == 0) {
      cout << " **** Event # " << nev <<endl;
    }
    
    Int_t nparticles = gAlice->GetEvent(nev);
    if (printLevel >= 1) {
      cout << " Total number of nparticles  " << nparticles <<endl;
    }
    if (nev < evNumber1) continue;
    if (nparticles <= 0) return;
     
    Double_t PxGen[2],PyGen[2],PzGen[2],PTotGen[2],ThetaGen[2],RapGen[2]; // parameters of generated  particles at the vertex
    ThetaGen[0] = ThetaGen[1] = 0;
    RapGen[0] = RapGen[1] = 0;

    for (int iPart = 0; iPart < nparticles  ; iPart++) {  // loop over particles
      
      TParticle *particle;
      particle = gAlice->Particle(iPart);    

      if ( particle->GetPdgCode() == kMuonPlus ) {
	PxGen[0] = particle->Px();
	PyGen[0] = particle->Py();
	PzGen[0] = particle->Pz();
	PTotGen[0] = TMath::Sqrt(PxGen[0]*PxGen[0]+PyGen[0]*PyGen[0]+PzGen[0]*PzGen[0]);
	ThetaGen[0] = TMath::ATan2(TMath::Sqrt(PxGen[0]*PxGen[0]+PyGen[0]*PyGen[0]),PzGen[0]);
	RapGen[0] = rapParticle(PxGen[0],PyGen[0],PzGen[0],massMuon);
	if (printLevel >= 1) {
	  cout << " Particle id : "<<particle->GetPdgCode()<<" px,py,pz,theta : "<<PxGen[0]<<" "<<PyGen[0]<<" "<<PzGen[0]<<" "<<ThetaGen[0]<<endl;
	}
      }
      if ( particle->GetPdgCode() == kMuonMinus ) {
	PxGen[1] = particle->Px();
	PyGen[1] = particle->Py();
	PzGen[1] = particle->Pz();
	PTotGen[1] = TMath::Sqrt(PxGen[1]*PxGen[1]+PyGen[1]*PyGen[1]+PzGen[1]*PzGen[1]);
	ThetaGen[1] = TMath::ATan2(TMath::Sqrt(PxGen[1]*PxGen[1]+PyGen[1]*PyGen[1]),PzGen[1]);
	RapGen[1] = rapParticle(PxGen[1],PyGen[1],PzGen[1],massMuon);
	if (printLevel >= 1)
	  cout << " Particle #: "<<particle->GetPdgCode()<<" px,py,pz,theta : "<<PxGen[1]<<" "<<PyGen[1]<<" "<<PzGen[1]<<" "<<ThetaGen[1]<<endl;
      }
    } // end loop over particles


    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
     
     
    TTree *TH = gAlice->TreeH();
    Int_t ntracks =(Int_t) TH->GetEntries();

    Bool_t hitChamberMUPlus[10],hitChamberMUMinus[10];
    Bool_t hitStationMUPlus[5],hitStationMUMinus[5];
    for (Int_t iCh=0; iCh < 10; iCh++){
       hitChamberMUPlus[iCh] = kFALSE;
       hitChamberMUMinus[iCh] = kFALSE;
    }
    for (Int_t iSt=0; iSt < 5; iSt++){
       hitStationMUPlus[iSt] = kFALSE;
       hitStationMUMinus[iSt] = kFALSE;
    }

    Int_t nHitChamber=0;
    
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
	     
	    
	    Int_t ipart = mHit->Particle();
	    Double_t x     = mHit->X();        // x-pos of hit in chamber #1
	    Double_t y     = mHit->Y();        // y-pos of hit in chamber #1
	    Double_t z     = mHit->Z();        // z-pos of hit in chamber #1

	    Int_t iSt = (nch-1)/2; 
	    if ((nch-1) < 10) {
	      if (ipart == kMuonPlus){
		hitChamberMUPlus[nch-1] = kTRUE;
		hitStationMUPlus[iSt] = kTRUE;
	      }
	      if (ipart == kMuonMinus){
		hitChamberMUMinus[nch-1] = kTRUE;
		hitStationMUMinus[iSt] = kTRUE;
	      }
	      nHitChamber++;
	    }
	     
	    if (nch != 1) continue;  
	     
	    if (ipart == kMuonPlus || ipart == kMuonMinus) {
	      Double_t px0=mHit->Px();
	      Double_t py0=mHit->Py();		   
	      Double_t pz0=mHit->Pz();
	      Double_t nonBendingSlope=px0/pz0;
	      Double_t bendingSlope=py0/pz0;

	      // create an AliMUONTrackParam object in chamber #1
	      AliMUONTrackParam *trackParam = new AliMUONTrackParam();
	      Int_t sign = -TMath::Sign(1.0,ipart);
	      Double_t bendingMometum = sign* TMath::Sqrt(pz0*pz0+py0*py0); // Bug px -> py !!!
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
	      if (!testSingle)
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
	      theta[iMuon] = TMath::ATan2(TMath::Sqrt(pX[iMuon]*pX[iMuon]+pY[iMuon]*pY[iMuon]),pZ[iMuon]);
	      radius[iMuon] = TMath::Sqrt(x*x+y*y);
	      
	      if (printLevel >= 1)
		cout << " px, py, pz, theta, radius = " <<pX[iMuon]<<" "<<pY[iMuon]<<" "<<pZ[iMuon]<<" "<<theta[iMuon]<<" "<<radius[iMuon]<<endl;
	      delete trackParam;
	      
	    } // if MuonPlus or MuonMinus
	  } // loop over MUON hits 
      } // if MUON 
    } // loop over tracks	


    Bool_t goodMUPlus = kTRUE;
    Bool_t goodMUMinus = kTRUE;
    for (Int_t iCh=0; iCh < 1; iCh++) { // !!!!  10 -> 1
      if (!hitChamberMUPlus[iCh]) {
	goodMUPlus = kFALSE;
	printf(" evt number %d no muon+ in chamber %d \n",nev,iCh+1);
      }
      if (!hitChamberMUMinus[iCh]) {
	goodMUMinus = kFALSE;
	printf(" evt number %d no muon- in chamber %d \n",nev,iCh+1);
      }
    }

//     for (Int_t iSt=0; iSt < 5; iSt++) { 
//       if (!hitStationMUPlus[iSt]) {
// 	goodMUPlus = kFALSE;
// // 	printf(" evt number %d no muon+ in chamber %d \n",nev,iCh+1);
//       }
//       if (!hitStationMUMinus[iSt]) {
// 	goodMUMinus = kFALSE;
// // 	printf(" evt number %d no muon- in chamber %d \n",nev,iCh+1);
//       }
//     }

    if (pX[0] != 0 && pX[1] != 0) { // if track 1 & 2 exist 
      Float_t massResonance = MuPlusMuMinusMass(pX[0],pY[0],pZ[0],pX[1],pY[1],pZ[1],massMuon);
      hInvMassRes->Fill(massResonance);
    }

   
    for (Int_t i = 0; i<2; i++){  // !!! 1 -> 2
      if (i == 0 && !goodMUPlus) continue; 
      if (i == 1 && !goodMUMinus) continue;
      if (pX[i] == 0) continue;

      hRap->Fill(RapGen[i]);
      
      Double_t cosAngle = pX[i]*PxGen[i]+pY[i]*PyGen[i]+pZ[i]*PzGen[i];
      cosAngle = cosAngle/(pTot[i]*PTotGen[i]);
      Double_t deltaAngle = TMath::ACos(cosAngle)*180/TMath::Pi();
      Double_t deltaPNorm = 0;
      if (testSingle) {
	deltaPNorm = (PTotGen[i]-pTot[i])* TMath::Cos(ThetaGen[i]);
      } else {
	deltaPNorm = PTotGen[i]-pTot[i];	  
      }

      Float_t thetaGX = TMath::ATan2(PxGen[i],PzGen[i]);
      Float_t thetaRX = TMath::ATan2(pX[i],pZ[i]);
      Float_t deltaAngleX = thetaRX - thetaGX; 

      Float_t thetaGY = TMath::ATan2(PyGen[i],PzGen[i]);
      Float_t thetaRY = TMath::ATan2(pY[i],pZ[i]);
      Float_t deltaAngleY = thetaRY - thetaGY; 
      
      Float_t phiG = TMath::ATan2(PyGen[i],PxGen[i])*180/TMath::Pi();
      Float_t phiR = TMath::ATan2(pY[i],pX[i])*180/TMath::Pi();
      
      hDeltaPxvsPhi->Fill(phiG,PxGen[i]-pX[i]);
      hDeltaPyvsPhi->Fill(phiG,PyGen[i]-pY[i]);
      hDeltaPhivsPz->Fill(PzGen[i],phiR-phiG);

      Int_t iHist=0;
      if (ThetaGen[i] < (3*TMath::Pi()/180)) {
	hDeltaP1->Fill(pTot[i]-PTotGen[i]);
	hDeltaAngle1->Fill(deltaAngle);
	hP1->Fill(PTotGen[i]);
	hProfDeltaPvsP1->Fill(PTotGen[i],deltaPNorm);
	h2DeltaPvsP1->Fill(PTotGen[i],deltaPNorm);
	iHist = (PTotGen[i]-xProfMin)*nBinProf/(xProfMax-xProfMin);
	if (PTotGen[i] > xProfMin && PTotGen[i] < xProfMax  ) {
	  DeltaP1[iHist]->Fill(deltaPNorm);
	  DeltaAngleX1[iHist]->Fill(deltaAngleX);
	  DeltaAngleY1[iHist]->Fill(deltaAngleY);
	}
      } else {
	hDeltaP2->Fill(pTot[i]-PTotGen[i]);
	hDeltaAngle2->Fill(deltaAngle);
	hP2->Fill(PTotGen[i]);
	hProfDeltaPvsP2->Fill(PTotGen[i],deltaPNorm);
	h2DeltaPvsP2->Fill(PTotGen[i],deltaPNorm);
	iHist = (PTotGen[i]-xProfMin)*nBinProf/(xProfMax-xProfMin);
	if (PTotGen[i] > xProfMin && PTotGen[i] < xProfMax ) {
	  DeltaP2[iHist]->Fill(deltaPNorm);	
	  DeltaAngleX2[iHist]->Fill(deltaAngleX);
	  DeltaAngleY2[iHist]->Fill(deltaAngleY);
	}
      }
    }
  } // loop over event
  

  Float_t sigmaP = 0;
  Float_t xBin = 0;
  for (Int_t iBin = 1; iBin <= nBinProf; iBin++){
    sigmaP = hProfDeltaPvsP1->GetBinError(iBin);
    xBin = hProfDeltaPvsP1->GetBinCenter(iBin);
    hSigmaPvsP1->SetBinContent(iBin,sigmaP/xBin); 
    sigmaP = hProfDeltaPvsP2->GetBinError(iBin);
    xBin = hProfDeltaPvsP2->GetBinCenter(iBin);
    hSigmaPvsP2->SetBinContent(iBin,sigmaP/xBin); 
  }

  TF1  *g1= new TF1("g1","gaus",-25,25) ; 
  TF1  *g2= new TF1("g2","gaus",-25,25) ; 
  Float_t sigmaPGaus;
  Float_t xBinGaus;
  for (Int_t iHist=0; iHist < nBinProf; iHist++){
    DeltaP1[iHist]->Fit("g1","RQ");
    sigmaPGaus =  g1->GetParameter(2);
    printf(" ** 1 ** iHist= %d , sigmaPGaus = %f \n",iHist,sigmaPGaus);
    xBinGaus = hSigGausPvsP1->GetBinCenter(iHist+1);
    hSigGausPvsP1->SetBinContent(iHist+1,sigmaPGaus/xBinGaus);
    
    DeltaP2[iHist]->Fit("g2","RQ");
    sigmaPGaus =  g2->GetParameter(2);
    printf(" ** 2 ** iHist= %d , sigmaGaus = %f \n",iHist,sigmaPGaus);
    xBinGaus = hSigGausPvsP2->GetBinCenter(iHist+1);
    hSigGausPvsP2->SetBinContent(iHist+1,sigmaPGaus/xBinGaus); 
  }   


  // Angles 
  TF1  *g3= new TF1("g3","gaus") ; 
  TF1  *g4= new TF1("g4","gaus") ; 
  TF1  *g5= new TF1("g5","gaus") ; 
  TF1  *g6= new TF1("g6","gaus") ; 
  Float_t sigmaAngleGaus,limitGaus,errorSigma;
  Float_t nSigFit = 3;
  for (Int_t iHist=0; iHist < nBinProf; iHist++){
    limitGaus = nSigFit * (DeltaAngleX1[iHist]->GetRMS());
    g3->SetRange(-limitGaus,limitGaus);
    DeltaAngleX1[iHist]->Fit("g3","RQ");
    sigmaAngleGaus =  g3->GetParameter(2);
    errorSigma = g3->GetParError(2);
    printf(" ** 1 ** iHist= %d , sigmaAngleGaus X = %f \n",iHist,sigmaAngleGaus);
    hSigGausAngleXvsP1->SetBinContent(iHist+1,sigmaAngleGaus);
    hSigGausAngleXvsP1->SetBinError(iHist+1,errorSigma);

    limitGaus = nSigFit * (DeltaAngleY1[iHist]->GetRMS());
    g4->SetRange(-limitGaus,limitGaus);
    DeltaAngleY1[iHist]->Fit("g4","RQ");
    sigmaAngleGaus =  g4->GetParameter(2);
    errorSigma = g4->GetParError(2);
    printf(" ** 1 ** iHist= %d , sigmaAngleGaus Y = %f \n",iHist,sigmaAngleGaus);
    hSigGausAngleYvsP1->SetBinContent(iHist+1,sigmaAngleGaus);
    hSigGausAngleYvsP1->SetBinError(iHist+1,errorSigma);
    
    limitGaus = nSigFit * (DeltaAngleX2[iHist]->GetRMS());
    g5->SetRange(-limitGaus,limitGaus);
    DeltaAngleX2[iHist]->Fit("g5","RQ");
    sigmaAngleGaus =  g5->GetParameter(2);
    errorSigma = g5->GetParError(2);
    printf(" ** 1 ** iHist= %d , sigmaAngleGaus X = %f \n",iHist,sigmaAngleGaus);
    hSigGausAngleXvsP2->SetBinContent(iHist+1,sigmaAngleGaus);
    hSigGausAngleXvsP2->SetBinError(iHist+1,errorSigma);

    limitGaus = nSigFit * (DeltaAngleY2[iHist]->GetRMS());
    g6->SetRange(-limitGaus,limitGaus);
    DeltaAngleY2[iHist]->Fit("g6","RQ");
    sigmaAngleGaus =  g6->GetParameter(2);
    errorSigma = g6->GetParError(2);
    printf(" ** 1 ** iHist= %d , sigmaAngleGaus Y = %f \n",iHist,sigmaAngleGaus);
    hSigGausAngleYvsP2->SetBinContent(iHist+1,sigmaAngleGaus);
    hSigGausAngleYvsP2->SetBinError(iHist+1,errorSigma);
    
  }  
  // save histograms in MUONTestAbso.root
  TFile *histoFile = new TFile("MUONTestAbso.root", "RECREATE");

  hInvMassRes->Write();
  hRap->Write();

  hDeltaP1->Write();
  hDeltaAngle1->Write();
  hP1->Write();
  hProfDeltaPvsP1->Write();
  h2DeltaPvsP1->Write();
  hSigmaPvsP1->Write();
  hSigGausPvsP1->Write();
  hSigGausAngleXvsP1->Write();
  hSigGausAngleYvsP1->Write();

  hDeltaP2->Write();
  hDeltaAngle2->Write();
  hP2->Write();
  hProfDeltaPvsP2->Write();
  h2DeltaPvsP2->Write();
  hSigmaPvsP2->Write();
  hSigGausPvsP2->Write();
  hSigGausAngleXvsP2->Write();
  hSigGausAngleYvsP2->Write();

  for (Int_t iHist=0; iHist < nBinProf; iHist++){
     DeltaP1[iHist]->Write();
     DeltaP2[iHist]->Write();
     DeltaAngleX1[iHist]->Write();
     DeltaAngleY1[iHist]->Write();
     DeltaAngleX2[iHist]->Write();
     DeltaAngleY2[iHist]->Write();
  }

  hDeltaPxvsPhi->Write();
  hDeltaPyvsPhi->Write();
  hDeltaPhivsPz->Write();

  histoFile->Close();
  
}


Float_t MuPlusMuMinusMass(Float_t Px1, Float_t Py1, Float_t Pz1, Float_t Px2, Float_t Py2, Float_t Pz2, Float_t MuonMass)
{
  // return invariant mass for particle 1 & 2 whose masses are equal to MuonMass

  Float_t e1 = TMath::Sqrt(MuonMass * MuonMass + Px1 * Px1 + Py1 * Py1 + Pz1 * Pz1);
  Float_t e2 = TMath::Sqrt(MuonMass * MuonMass + Px2 * Px2 + Py2 * Py2 + Pz2 * Pz2);
  return (TMath::Sqrt(2.0 * (MuonMass * MuonMass + e1 * e2 - Px1 * Px2 - Py1 * Py2 - Pz1 * Pz2)));
}

Float_t rapParticle(Float_t Px, Float_t Py, Float_t Pz, Float_t Mass)
{
  // return rapidity for particle 

  // Particle energy
  Float_t Energy = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz+Mass*Mass);
  // Rapidity
  Float_t Rapidity = 0.5*TMath::Log((Energy+Pz)/(Energy-Pz));
  return Rapidity;
}
