void MUONcombi (Int_t evNumber=0) 
{
    const Float_t runWeight = 4.e8;
    
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice = 0;
   }

   
//
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   TClonesArray *PartArray = gAlice->Particles();
   TParticle *Part;
   

// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
//
//
   TH1F *dmass  = new TH1F("dmass","Dimuon-Mass Distribution"
			  ,100,0.,5.);
   
   TH1F *dmassc = new TH1F("dmassc","Dimuon-Mass Distribution"
			   ,200,0.,10.);
   TH1F *dmassd = new TH1F("dmassd","Dimuon-Mass Distribution"
			   ,200,0.,10.);

   TH1F *pt     = new TH1F("pt","pT-single"
			   ,50,0.,10.);
   TH1F *hCont[10];
   
//
//  Generator Loop
//

   Int_t i=0;
   
   AliGenCocktailEntry *Entry, *e1, *e2;
   AliGenCocktail*  Cocktail = (AliGenCocktail*) gAlice->Generator();
//                   Single Generator
   for (Entry=Cocktail->FirstGenerator();
	Entry;
	Entry=Cocktail->NextGenerator()
       ) {
       Entry->PrintInfo();
     }
//
// Initialize Combinator 
//
   AliDimuCombinator* Combinator = new AliDimuCombinator();
   Combinator->SetEtaCut(2.5, 4.);
   Combinator->SetPtMin(1.0);   
   
   i=0;
   
//
// Single Muon  Loop  
//
   Combinator->ResetRange();
   
   

   for(Muon=Combinator->FirstMuon();
       Muon; 
       Muon=Combinator->NextMuon()) {
//
       Int_t chfirst = Muon->GetFirstDaughter();
       Int_t chlast  = Muon->GetLastDaughter();
       Int_t parent  = Muon->GetFirstMother();
       Float_t ptm   = Muon->Pt();
       Float_t eta   = Muon->Eta();
       
//       printf("\n Particle %d Parent %d first child %d last child %d",
//	      i,parent, chfirst, chlast);
//       printf("\n Particle pt, eta: %f , %f ", ptm, eta);
       i++;       
   }
//
// Di-Muon Loop

   Float_t pt1,pt2;
   TParticle* Muon1;
   TParticle* Muon2;

   Combinator->ResetRange();

//
//   Dimuon Loop controlled by Generator Loop
//
   Float_t  sig = 0;
   Int_t  icont = 0;
   char name[30];   
   for (Cocktail->FirstGeneratorPair(e1,e2);
	(e1&&e2);
	Cocktail->NextGeneratorPair(e1,e2)
       ){

       sprintf(name, "%s-%s", e1->GetName(), e2->GetName());
       hCont[icont] = new TH1F(name,"Dimuon-Mass Distribution",100,0.,5.);
       printf("\n ----------------------------------------------------");
       e1->PrintInfo();
       e2->PrintInfo();
       printf("\n ----------------------------------------------------");
       Combinator->SetFirstRange (e1->GetFirst(), e1->GetLast());
       Combinator->SetSecondRange(e2->GetFirst(), e2->GetLast());
       Combinator->SetRate(e1->Rate(), e2->Rate());
       for (Combinator->FirstMuonPairSelected(Muon1,Muon2);
	    (Muon1 && Muon2);
	    Combinator->NextMuonPairSelected(Muon1,Muon2))
       {
	   pt1 = Muon1->Pt();
	   pt2 = Muon2->Pt();       
	   Float_t mass = Combinator->Mass(Muon1, Muon2);
	   Float_t wgt  = runWeight*Combinator->Weight(Muon1, Muon2)*
	       Combinator->DecayProbability(Muon1)*
	       Combinator->DecayProbability(Muon2);
	   pt->Fill(pt1, wgt);
	   pt->Fill(pt2, wgt);
	   Float_t smeared_mass = mass;
	   Combinator->SmearGauss(0.02*mass, smeared_mass);
	   
	   if (TMath::Min(pt1,pt2) > -0.5*TMath::Max(pt1,pt2)+2.) {
	       if (Combinator->Correlated(Muon1, Muon2)) {
		   dmassc->Fill(smeared_mass, wgt);
		   sig += wgt;
	       } else {
// account for the fact that we sum like-sign and unlike-sign		   
		   wgt *= 0.5;
	       }
	       dmass->Fill(smeared_mass, wgt);
	       hCont[icont]->Fill(smeared_mass, wgt);
	   }
       } // Dimuon Loop
       icont++;
   }// Generator Loop
   printf("\n Signal %e \n \n", sig);
   
//
//Create a canvas, set the view range, show histograms
//
   TCanvas *c1 = new TCanvas("c1","Dimuon Plots",400,10,600,700);
   TPad* pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
   TPad* pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
   TPad* pad13 = new TPad("pad13"," ",0.01,0.01,0.49,0.49);
   TPad* pad14 = new TPad("pad14"," ",0.51,0.01,0.99,0.49);
   pad11->SetFillColor(11);
   pad12->SetFillColor(11);
   pad13->SetFillColor(11);
   pad14->SetFillColor(11);
   pad11->Draw();
   pad12->Draw();
   pad13->Draw();
   pad14->Draw();
   
   pad11->cd();
   pt->SetFillColor(42);
   pt->SetXTitle("pT (GeV)");
   pt->Draw();

   pad12->cd();
   pad12->SetLogy(0);
   dmass->SetFillColor(42);
   dmass->SetXTitle("m (GeV)");
   dmass->Draw();

   pad13->cd();
   pad13->SetLogy(0);
   dmassc->SetFillColor(42);
   dmassc->SetXTitle("m (GeV)");
   dmassc->Draw();

   pad14->cd();
   pad14->SetLogy(0);
   dmassd->SetFillColor(42);
   dmassd->SetXTitle("m (GeV)");
   dmassd->Draw();


}










