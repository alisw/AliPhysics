// option 1 : idres=114    if j/psi
//            idres=116    if upsilon

// option 2 : ireadgeant=1 if geant hits
//            ireadgeant=0 if space points

// option 3 : ibgr=1 if upsilon+background (option 3 usefull only for geant hits) 
//            ibgr=0 if no background

void MUONtestzaza (Int_t evNumber1=0, Int_t evNumber2=99, Int_t idres=116, Int_t ireadgeant=1, Int_t ibgr=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

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
    } else {	printf("\n galice.root found in file list");
    }
    //    file->ls();
    
// Get AliRun object from file or create it if not on file
    if (!gAlice) {
	gAlice = (AliRun*)(file->Get("gAlice"));
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) {
	    printf("\n create new gAlice object");
	    gAlice = new AliRun("gAlice","Alice test program");
	}
    }
    
    Double_t seff  = 0.95;
    Double_t sb0   = 0.7;
    Double_t sbl3  = 0.2;

    Int_t ifit   = 0;
    Int_t idebug = 1;

    AliMUON *MUON  = (AliMUON*) gAlice->GetModule("MUON");
    
    Int_t nparticles = gAlice->GetEvent(evNumber1);
    if (nparticles <= 0) return;
    MUON->Init(seff,sb0,sbl3);
//   Loop over events 
//
    Int_t inev_bgd=0;
    Int_t nev_bgd=4;
    
    for (Int_t nev= evNumber1; nev<= evNumber2; nev++) 
    {
	printf("nev=%d\n",nev);
	if (nev != evNumber1) Int_t nparticles = gAlice->GetEvent(nev);
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	MUON->FinishEvent();
	if (ireadgeant==1 && ibgr==1) {
	  if (inev_bgd==nev_bgd) inev_bgd=0;
	  MUON->Reconst(ifit,idebug,inev_bgd,nev,idres,ireadgeant,"Add","galice_bgr.root");
	  inev_bgd++;
	}
	else
	  MUON->Reconst(ifit,idebug,inev_bgd,nev,idres,ireadgeant,"rien1","rien2");
	MUON->FinishEvent();


	
    } // event loop 
    
    MUON->Close();
}



