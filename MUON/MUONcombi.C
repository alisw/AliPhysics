void MUONcombi (Int_t evNumber=0) 
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
	gSystem->Load("$ALITOP/cern.so/lib/libpdfDUMMY.so");
	gSystem->Load("$ALITOP/cern.so/lib/libPythia.so");
	gSystem->Load("$ROOTSYS/lib/libEG.so");       
	gSystem->Load("$ROOTSYS/lib/libEGPythia.so");       
	gSystem->Load("libGeant3Dummy.so");    //a dummy version of Geant3
	gSystem->Load("PHOS/libPHOSdummy.so"); //the standard Alice classes 
	gSystem->Load("libgalice.so");         // the standard Alice classes 
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
//   TClonesArray &Partarray = *(gAlice->Particles());
   TClonesArray *PartArray = gAlice->Particles();
   TParticle *Part;
   

// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
//
//
   TH1F *dmass = new TH1F("dmass","Dimuon-Mass Distribution"
			 ,50,0.,10.);
   TH1F *pt    = new TH1F("pt","pT-single"
			 ,50,0.,10.);
//
//  Generator Loop
//
   
   
   AliGenCocktailEntry *Entry, *e1, *e2;
   AliGenCocktail*  Cocktail = (AliGenCocktail*) gAlice->Generator();
//                   Single Generator
   for (Entry=Cocktail->FirstGenerator();
	Entry;
	Entry=Cocktail->NextGenerator()
       ) {
       Entry->PrintInfo();
     }
//                   Pairs of Generators

//
// Initialize Combinator 
   DimuonCombinator Combinator = DimuonCombinator(PartArray);
   Combinator->SetEtaCut(-4, 4);
   Combinator->SetPtMin(0.5);   

   Int_t i=0;
//
// Single Muon  Loop  
//
/*

   for(Muon=Combinator.FirstMuon();
       Muon; 
       Muon=Combinator.NextMuon()) {
//
       Int_t chfirst= Muon->GetFirstChild();
       Int_t chlast = Muon->GetLastChild();
       Int_t parent = Muon->GetParent();
       Float_t ptm  = Muon->GetPT();
       Float_t eta  = Muon->GetEta();
       
       printf("\n Particle %d Parent %d first child %d last child %d",
	      i,parent, chfirst, chlast);
       printf("\n Particle pt, eta: %f , %f ",pt,eta);
       i++;       
   }
*/
//
// Di-Muon Loop

   Float_t pt1,pt2;
   TParticle* Muon1, *Muon2;
/*
   Combinator->ResetRange();
   for (Combinator->FirstMuonPairSelected(Muon1,Muon2);
	(Muon1 && Muon2);
	Combinator->NextMuonPairSelected(Muon1,Muon2))
   {
       pt1=Muon1->GetPT();
       pt2=Muon2->GetPT();       
       Float_t mass=Combinator->Mass(Muon1, Muon2);
       Float_t wgt =Combinator->Weight(Muon1, Muon2);
       pt->Fill(pt1, wgt);
       pt->Fill(pt2, wgt);
       Float_t smeared_mass=mass;
       Combinator->SmearGauss(0.05*mass, smeared_mass);
       dmass->Fill(smeared_mass , wgt);
    }
*/
//
//   Dimuon Loop controlled by Generator Loop
//
   for (Cocktail->FirstGeneratorPair(e1,e2);
	(e1&&e2);
	Cocktail->NextGeneratorPair(e1,e2)
       ){
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
	   pt1=Muon1->GetPT();
	   pt2=Muon2->GetPT();       
	   Float_t mass=Combinator->Mass(Muon1, Muon2);
	   Float_t wgt =Combinator->Weight(Muon1, Muon2);
	   pt->Fill(pt1, wgt);
	   pt->Fill(pt2, wgt);
	   Float_t smeared_mass=mass;
	   Combinator->SmearGauss(0.05*mass, smeared_mass);
	   dmass->Fill(smeared_mass , wgt);
       } // Dimuon Loop
   }// Generator Loop

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
   dmass->SetFillColor(42);
   dmass->SetXTitle("m (GeV)");
   dmass->Draw();
}










