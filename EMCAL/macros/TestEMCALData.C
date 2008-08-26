


void TestEMCALData() {

  TH1F* hE = new TH1F("hEmcalHits",    "Hits energy distribution in EMCAL",       100, 0., 10.) ;
  hE->Sumw2() ;
  TH1I* hM = new TH1I("hEmcalHitsMul", "Hits multiplicity distribution in EMCAL", 500, 0., 10000) ;
  hM->Sumw2() ;

  TH1I * dA = new TH1I("hEmcalDigits",    "Digits amplitude distribution in EMCAL",    500, 0, 5000) ;
  dA->Sumw2() ;
  TH1I * dM = new TH1I("hEmcalDigitsMul", "Digits multiplicity distribution in EMCAL", 500, 0, 1000) ;
  dM->Sumw2() ;

  TH1F * sE = new TH1F("hEmcalSDigits",    "SDigits energy distribution in EMCAL",    200, 0, 1000) ;
  sE->Sumw2() ;
  TH1I * sM = new TH1I("hEmcalSDigitsMul", "SDigits multiplicity distribution in EMCAL", 500, 0, 1000) ;
  sM->Sumw2() ;

  TH1F * cE = new TH1F("hEmcalRecPoints",    "RecPoints energy distribution in EMCAL",    100, 0, 20) ;
  cE->Sumw2() ;
  TH1I * cM = new TH1I("hEmcalRecPointsMul", "RecPoints multiplicity distribution in EMCAL", 500, 0, 1000) ;
  cM->Sumw2() ;

  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  if (rl == 0x0)
    cout<<"Can not instatiate the Run Loader"<<endl;

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));

   //Load Hits
   rl->LoadHits("EMCAL");
   //Load Digits
   rl->LoadDigits("EMCAL");
   //Load SDigits
   rl->LoadSDigits("EMCAL");
   //Load RecPoints
   rl->LoadRecPoints("EMCAL");

   
   //one way to do it that works for all three
   AliEMCALHit* hit;
   AliEMCALDigit* dig;
   AliEMCALRecPoint* rp;

   rl->GetEvent(0);
   
   //Fill array of hits                                                                        
   TClonesArray *hits = emcalLoader->Hits();
   //Fill array of digits                                                                        
   TClonesArray *digits = emcalLoader->Digits();
   //Fill array of clusters                                                                        
   TClonesArray *clusters = emcalLoader->RecPoints();

   //Get hits from the list                                                                    
   hM->Fill(hits->GetEntries());
   for(Int_t ihit = 0; ihit< hits->GetEntries();ihit++){
     //cout<<">> idig "<<idig<<endl;                                                             
     hit = static_cast<AliEMCALHit *>(hits->At(ihit)) ;

     if(hit != 0){
       hE->Fill(hit->GetEnergy());
       cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy() << endl;
     }
   }

   //Get digits from the list                                                                    
   dM->Fill(digits->GetEntries());
   for(Int_t idig = 0; idig< digits->GetEntries();idig++){
     //cout<<">> idig "<<idig<<endl;                                                             
     dig = static_cast<AliEMCALDigit *>(digits->At(idig)) ;

     if(dig != 0){
       dA->Fill(dig->GetAmp());
       cout << "Digit info " << dig->GetId() << " " << dig->GetAmp() << endl;
     }
   }


   //Get clusters from the list                                                                    
   cM->Fill(clusters->GetEntries());
   for(Int_t iclu = 0; iclu< clusters->GetEntries();iclu++){
     //cout<<">> idig "<<idig<<endl;                                                             
     rp = static_cast<AliEMCALRecPoint *>(clusters->At(iclu)) ;

     if(rp != 0){
       cE->Fill(rp->GetEnergy());
       cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
     }
   }


   /*
   //another way to do it
   TTree *hTree = rl->GetTreeH("EMCAL",false);
   TTree *dTree = rl->GetTreeD("EMCAL",false);
   TTree *sTree = rl->GetTreeS("EMCAL",false);
   TTree *cTree = rl->GetTreeR("EMCAL",false);
   
   TObjArray *harr=NULL;
   TBranch *hbranch=hTree->GetBranch("EMCAL");
   hbranch->SetAddress(&harr);
   
   TObjArray *darr=NULL;
   TBranch *dbranch=dTree->GetBranch("EMCAL");
   dbranch->SetAddress(&darr);
   
   TObjArray *sarr=NULL;
   TBranch *sbranch=sTree->GetBranch("EMCAL");
   sbranch->SetAddress(&sarr);
   
   TObjArray *carr=NULL;
   TBranch *cbranch=cTree->GetBranch("EMCALECARP");
   cbranch->SetAddress(&carr);
   

   if(hbranch->GetEvent(0)) {
     for(Int_t ih = 0; ih < harr->GetEntriesFast(); ih++) {
       hM->Fill(harr->GetEntriesFast());
       AliEMCALHit* hit =(AliEMCALHit*)harr->UncheckedAt(ih);
       if(hit != 0){
	 hE->Fill(hit->GetEnergy());
	 cout << "Hit info " << hit->GetId() << " " << hit->GetEnergy()*10.5 << endl;
       }
     }
   }

   if(dbranch->GetEvent(0)) {
     for(Int_t id = 0; id < darr->GetEntriesFast(); id++) {
       dM->Fill(darr->GetEntriesFast());
       AliEMCALDigit* dig =(AliEMCALDigit*)darr->UncheckedAt(id);
       if(dig != 0){
	 dA->Fill(dig->GetAmp());
	 cout << "Digit info " << dig->GetId() << " " << dig->GetAmp() << endl;
       }
     }
   }

   if(sbranch->GetEvent(0)) {
     for(Int_t id = 0; id < sarr->GetEntriesFast(); id++) {
       sM->Fill(sarr->GetEntriesFast());
       AliEMCALDigit* sdig =(AliEMCALDigit*)sarr->UncheckedAt(id);
       if(sdig != 0){
	 sE->Fill(sdig->GetAmp()/1.e+6);
	 cout << "SDigit info " << sdig->GetId() << " " << sdig->GetAmp()/1.e+6 << endl;
       }
     }
   }

   if(cbranch->GetEvent(0)) {
     for(Int_t ic = 0; ic < carr->GetEntriesFast(); ic++) {
       cE->Fill(carr->GetEntriesFast());
       AliEMCALRecPoint* rp =(AliEMCALRecPoint*)carr->UncheckedAt(ic);
       if(rp != 0){
	 cE->Fill(rp->GetEnergy());
	 cout << "RecPoint info " << rp->GetAbsId() << " " << rp->GetEnergy() << endl;
       }
     }
   }
   
   */
       
   TCanvas *chits = new TCanvas("chits","Hits",20,20,800,400);
   chits->Divide(2,1);
   chits->cd(1);
   hE->Draw();
   chits->cd(2);
   hM->Draw();
       
   TCanvas *cdig = new TCanvas("cdig","Digits",20,40,800,400);
   cdig->Divide(2,1);
   cdig->cd(1);
   dA->Draw();
   cdig->cd(2);
   dM->Draw();

   /*
   TCanvas *csdig = new TCanvas("csdig","SDigits",20,60,800,400);
   csdig->Divide(2,1);
   csdig->cd(1);
   sE->Draw();
   csdig->cd(2);
   sM->Draw();
   */

   TCanvas *cclu = new TCanvas("cclu","Clusters",20,60,800,400);
   cclu->Divide(2,1);
   cclu->cd(1);
   cE->Draw();
   cclu->cd(2);
   cM->Draw();

}
