void Fragment(Int_t nev=0, Int_t debug=0)
{
 gROOT->Reset();
 
 Float_t b;
 Int_t nspectp, nspectn, nspectpfree, nspectnfree; 
 Int_t zz[100], nn[100], nalpha, NFrag, ztot, ntot, ndeu;
 void spectator(Float_t, Int_t*, Int_t*);
 
 TH2F *hsp = new TH2F("hsp","Spectator protons vs b",50,0.,20.,90,0.,90.);
 hsp -> SetXTitle("b (fm)");
 hsp -> SetYTitle("Num. of spectator protons");
 
 TH2F *hsn = new TH2F("hsn","Spectator neutrons vs b",50,0.,20.,130,0.,130.);
 hsn -> SetXTitle("b (fm)");
 hsn -> SetYTitle("Num. of spectator neutrons");
 
 TH2F *hFragp = new TH2F("hFragp","Free spectator protons vs b",50,0.,20.,60,0.,60.);
 hFragp -> SetXTitle("b (fm)");
 hFragp -> SetYTitle("Num. of free spectator protons");
 
 TH2F *hFragn = new TH2F("hFragn","Free spectator neutrons vs b",50,0.,20.,80,0.,80.);
 hFragn -> SetXTitle("b (fm)");
 hFragn -> SetYTitle("Num. of free spectator neutrons");

 TF1 *funb = new  TF1("funb","x",0,20);
 for(Int_t ievent=0; ievent<nev; ievent++){  
   if(ievent%1000==0){printf("Processing event %d\n",ievent);}
   b= Float_t(funb->GetRandom()); 
   spectator(b, &nspectp, &nspectn);
   hsp -> Fill(b,nspectp);
   hsn -> Fill(b,nspectn);
   AliZDCFragment *gallio = new AliZDCFragment(b);
   for(Int_t j=0; j<=99; j++){
      zz[j] =0;
      nn[j] =0;
   }
//
// Generation of fragments
   gallio->GenerateIMF();
   nalpha = gallio->GetNalpha();
   NFrag = gallio->GetFragmentNum();
//
   // Attach neutrons
   ztot = gallio->GetZtot();
   ntot = gallio->GetNtot();
   gallio->AttachNeutrons();
   //
   nspectpfree = (nspectp-ztot-2*nalpha);
   nspectnfree = (nspectn-ntot-2*nalpha);
   
   // Removing deuterons
   ndeu = (Int_t) (nspectnfree*gallio->DeuteronNumber());
   nspectpfree -= ndeu;
   nspectnfree -= ndeu;
   //
   hFragp -> Fill(b, nspectpfree);
   hFragn -> Fill(b, nspectnfree);
//
// Print
   if(debug == 1){
     printf("\n	b = %f",b);
     printf("	#spect p = %d, #spect n = %d\n",nspectp,nspectn); 
     printf("	#spect p free = %d, #spect n free = %d\n",nspectpfree,nspectnfree); 
     printf("	#fragments = %f ",NFrag);
     /*for(Int_t i=0; i<NFrag; i++){
   	printf("\n	ZZ[%d] = %d, NN[%d] = %d\n",i,zz[i],i,nn[i]);
     }*/
     printf("	NAlpha = %d, Ztot = %d, Ntot = %d\n\n",nalpha,ztot,ntot);
   }
   delete gallio;
  } //Event loop
  
  TProfile *profsp = hsp->ProfileX("profsp",-1,-1,"s");
  TProfile *profsn = hsn->ProfileX("profsn",-1,-1,"s");
  TProfile *profFragp = hFragp->ProfileX("profFragp",-1,-1,"s");
  TProfile *profFragn = hFragn->ProfileX("profFragn",-1,-1,"s");

   
  //***********************************************************
  // #### ROOT initialization
  gROOT->Reset();
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  //
  TCanvas *c1 = new TCanvas("c1","Fragmentation I",0,0,700,700);
  c1->cd();
  c1->SetFillColor(kAzure+7);
  c1->Divide(2,2);
  c1->cd(1);
  hsp -> Draw("colz");   
  c1->cd(2);
  hsn -> Draw("colz");    
  c1->cd(3);
  hFragp -> Draw("colz");   
  c1->cd(4);
  hFragn -> Draw("colz"); 
  //
  c1->Print("Fragment1.gif");   
  //
  TCanvas *c2 = new TCanvas("c2","Fragmentation II",300,10,700,700);
  c2->cd();
  c2->SetFillColor(kAzure+10);
  c2->Divide(2,2);
  c2->cd(1);
  profsp ->SetLineColor(kRed);  profsp ->SetLineWidth(2);
  profsp -> Draw("colz");   
  c2->cd(2);
  profsn ->SetLineColor(kRed);  profsn ->SetLineWidth(2);
  profsn -> Draw("colz");    
  c2->cd(3);
  profFragp -> Draw("colz");   
  profFragp ->SetLineColor(kAzure);  profFragp ->SetLineWidth(2);
  c2->cd(4);
  profFragn -> Draw("colz");    
  profFragn ->SetLineColor(kAzure);  profFragn ->SetLineWidth(2);
  //
  c1->Print("Fragment2.gif");   

}


void spectator(Float_t b, Int_t* NSpectp, Int_t* NSpectn)
{
   Float_t SppvsB[6] = {3.633,-1.518,1.360,-.4899e-1,-.2398e-2,.1066e-3};
   Float_t SpnvsB[6] = {5.639,-1.685,1.803,-.3129e-1,-.6618e-2,.2352e-3};
   Float_t Sigmap[4] = {.5668,-.2200e-1,.3657e-3,-.2201e-5};
   Float_t Sigman[4] = {.4185,-.9798e-2,.1052e-3,-.4238e-6};
   
   Float_t rnsp = SppvsB[0]+SppvsB[1]*b+SppvsB[2]*(b*b)+SppvsB[3]*(b*b*b)+
   	          SppvsB[4]*(b*b*b*b)+SppvsB[5]*(b*b*b*b*b);
   Float_t rnsn = SpnvsB[0]+SpnvsB[1]*b+SpnvsB[2]*(b*b)+SpnvsB[3]*(b*b*b)+
   	          SpnvsB[4]*(b*b*b*b)+SpnvsB[5]*(b*b*b*b*b);
   Float_t snsp = Sigmap[0]+Sigmap[1]*rnsp+Sigmap[2]*(rnsp*rnsp)+Sigmap[3]*
   		  (rnsp*rnsp*rnsp); 
   Float_t snsn = Sigman[0]+Sigman[1]*rnsn+Sigman[2]*(rnsn*rnsn)+Sigman[3]*
   		  (rnsn*rnsn*rnsn); 
   
   snsp = snsp*rnsp;
   snsn = snsn*rnsn;
   
   Float_t xgaup = gRandom->Gaus(0.0,1.0);
   snsp = snsp*xgaup;
   Float_t xgaun = gRandom->Gaus(0.0,1.0);
   snsn = snsn*xgaun;
   rnsp=rnsp+snsp;
   rnsn=rnsn+snsn;
   
   *NSpectp = Int_t(rnsp);
   *NSpectn = Int_t(rnsn);

}
