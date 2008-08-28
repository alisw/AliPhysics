void Fragment(Int_t nev=0, Int_t debug=0)
{
 gROOT->Reset();
 
 Float_t b;
 Int_t nspectp, nspectn, nspectpfree, nspectnfree; 
 Int_t zz[100], nn[100], nalpha, ztot, ntot;
 void spectator(Float_t, Int_t*, Int_t*);
 
 TH2F *hsp = new TH2F("hsp","Spectator protons vs b",100,0.,20.,100,0.,150.);
 hsp -> SetXTitle("b (fm)");
 hsp -> SetYTitle("Num. of spectator protons");
 
 TH2F *hsn = new TH2F("hsn","Spectator neutrons vs b",100,0.,20.,100,0.,150.);
 hsn -> SetXTitle("b (fm)");
 hsn -> SetYTitle("Num. of spectator neutrons");
 
 TH2F *hFragp = new TH2F("hFragp","Free spectator protons vs b",100,0.,20.,100,0.,150.);
 hFragp -> SetXTitle("b (fm)");
 hFragp -> SetYTitle("Num. of free spectator protons");
 
 TH2F *hFragn = new TH2F("hFragn","Free spectator neutrons vs b",100,0.,20.,100,0.,150.);
 hFragn -> SetXTitle("b (fm)");
 hFragn -> SetYTitle("Num. of free spectator neutrons");

 TF1 *funb = new  TF1("funb","x",0,15);
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
   Int_t NFrag = gallio->GetFragmentNum();
//
// Attach neutrons
   ztot=0;
   ntot=0;
   gallio->AttachNeutrons();
   nspectpfree = nspectp-ztot-2*nalpha;
   nspectnfree = nspectn-ntot-2*nalpha;
   hFragp -> Fill(b,nspectpfree);
   hFragn -> Fill(b,nspectnfree);
//
// Print
   if(debug == 1){
     printf("\n	b = %f\n",b);
     printf("\n	#spect p = %d, #spect n = %d\n",nspectp,nspectn); 
     printf("\n	#spect p free = %d, #spect n free = %d\n",nspectpfree,nspectnfree); 
     printf("\n	#fragments = %f\n",NFrag);
     for(Int_t i=0; i<NFrag; i++){
   	printf("\n	ZZ[%d] = %d, NN[%d] = %d\n",i,zz[i],i,nn[i]);
     }
     printf("\n	NAlpha = %d, Ztot = %d, Ntot = %d\n\n",nalpha,ztot,ntot);
   }
   delete gallio;
  } //Event loop

   TCanvas *c1 = new TCanvas("c1","Fragmentation",0,10,580,700);
   c1->cd();
   c1->SetFillColor(29);
   c1->Divide(2,2);
   c1->cd(1);
   hsp -> Draw("box");   
   c1->cd(2);
   hsn -> Draw("box");    
   c1->cd(3);
   hFragp -> Draw("box");   
   c1->cd(4);
   hFragn -> Draw("box");    

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
