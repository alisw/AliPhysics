#include "AliGenFLUKAsource.h"
#include "AliGenMUONlib.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliConst.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
 ClassImp(AliGenFLUKAsource)
     AliGenFLUKAsource::AliGenFLUKAsource()
	 :AliGenerator(-1)
{
    //
    fName="FLUKA";
    fTitle="FLUKA Boundary Source";
    // Read in all particle types by default
    fIkine=6;
    // Set maximum admitted age of particles to 1.e-05 by default 
    fAgeMax=1.e-05;
    // Do not add weight
    fAddWeight=1.;
    // Shift the z-coordinate of the impact point by 4.5 cm only if it reads 
    // from  specific boundary source to the chamber (fZshift=4.5;),else there 
    // is no need to shift as it reads boundary source for the whole volume of 
    // the Muon Arm; the default value corresponds to boundary source for the
    // whole volume of the MUON Arm 
    fZshift=0;
    // Set the default file 
    fFileName="flukasource.root";

    fTreeFluka=0;
//
//  Read all particles
    fNpart=-1;
    
}

AliGenFLUKAsource::AliGenFLUKAsource(Int_t npart)
    :AliGenerator(npart)
{
    //
    fName="FLUKA";
    fTitle="FLUKA Boundary Source";
    // Read in all particle types by default
    fIkine=6;
    // Set maximum admitted age of particles to 1.e-05 by default 
    fAgeMax=1.e-05;
    // Do not add weight
    fAddWeight=1.;
    // Shift the z-coordinate of the impact point by 4.5 cm only if it reads 
    // from  specific boundary source to the chamber (fZshift=4.5;),else there 
    // is no need to shift as it reads boundary source for the whole volume of 
    // the Muon Arm; the default value corresponds to boundary source for the
    // whole volume of the MUON Arm 
    fZshift=0;
    // Set the default file 
    fFileName="flukasource.root";

    fTreeFluka=0;

}

//____________________________________________________________
AliGenFLUKAsource::~AliGenFLUKAsource()
{
    delete fTreeFluka;
}

//____________________________________________________________
void AliGenFLUKAsource::FlukaInit() 
{
//
// reset the existing file environment and open a new root file if
// the pointer to the Fluka tree is null
    
    TFile *File=0;
    if (fTreeFluka==0) {
        if (!File) {
	    File = new TFile(fFileName);
	    File->cd();
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
        }
// get the tree address in the Fluka boundary source file
	fTreeFluka = (TTree*)gDirectory->Get("h1");
    } else {
        File = fTreeFluka->GetCurrentFile();
        File->cd();
    }

    TTree *h2=fTreeFluka;
//Set branch addresses
    h2->SetBranchAddress("Ip",&Ip);
    h2->SetBranchAddress("Ipp",&Ipp);
    h2->SetBranchAddress("Xi",&Xi);
    h2->SetBranchAddress("Yi",&Yi);
    h2->SetBranchAddress("Zi",&Zi);
    h2->SetBranchAddress("Px",&Px);
    h2->SetBranchAddress("Py",&Py);
    h2->SetBranchAddress("Pz",&Pz);
    h2->SetBranchAddress("Ekin",&Ekin);
    h2->SetBranchAddress("Zv",&Zv);
    h2->SetBranchAddress("Rv",&Rv);
    h2->SetBranchAddress("Itra",&Itra);
    h2->SetBranchAddress("Igas",&Igas);
    h2->SetBranchAddress("Wgt",&Wgt);
    h2->SetBranchAddress("Etag",&Etag);
    h2->SetBranchAddress("Ptg",&Ptg);
    h2->SetBranchAddress("Age",&Age);
}

//____________________________________________________________
void AliGenFLUKAsource::Generate()
{

  AliMC* pMC = AliMC::GetMC();

  const Int_t ifluge[28]={kProton, kProtonBar, kElectron, kPositron,
			  kNuE, kNuEBar, kGamma, kNeutron, kNeutronBar,
			  kMuonPlus, kMuonMinus, kK0Long , kPiPlus, kPiMinus,
			  kKPlus, kKMinus, kLambda0, kLambda0Bar, kK0Short,
			  kSigmaMinus, kSigmaPlus, kSigma0, kPi0, kK0, kK0Bar,
			  0,kNuMu,kNuMuBar};
  Float_t polar[3]= {0,0,0};
  //
  Float_t origin[3];
  Float_t p[3];
  Float_t prwn;
  Float_t wgt, fwgt;
  Float_t phi;
  char name[100];
  Float_t amass, charge, tlife;
  Int_t itrtyp;
  Int_t iwgt;
  Int_t i, j, part, nt;
  static Int_t irwn;
  //
  Float_t random[2];
  //
  FlukaInit();
  TTree *h2=fTreeFluka;
  Int_t nentries = (Int_t) h2->GetEntries();
  if (fNpart == -1) fNpart=Int_t(nentries*fFrac);
  
  // loop over number of particles
  Int_t nb=0;
  for (i=0; i<fNpart;i++) {
    Int_t ev=pMC->CurrentEvent();
    Int_t entry=fNpart*(ev-1)+i; 
    nb = (Int_t)h2->GetEvent(entry); 
    if (irwn > nentries) {
      printf("No more entries in the FLUKA boundary source file\n");
      TFile *File=0;
      // Get AliRun object or create it 
      if (!gAlice) {
        gAlice = (AliRun*)File->Get("gAlice");
        if (gAlice) printf("AliRun object found on file\n");
        if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
      }
      TTree *fAli=gAlice->TreeK();
      if (fAli) File =fAli->GetCurrentFile();
      File->cd();
      printf("Generate - I'm out \n");
      return;
    }   
    if (Ip > 28 || Ip < 0) {
      irwn++;
      continue;
    }
 
    if ((Ip != fIkine && fIkine != 6 && fIkine != 9) || Age > fAgeMax){
       irwn++;
       continue;
    } else if (fIkine == 9) {
       if (Ip == 7 || Ip == 8 || Age > fAgeMax) { 
          irwn++;
          continue;
       }
    }
    

    irwn++;
    printf("\n Particle type: %f \n \n ", Ip);
    
    if (Ip ==7){
       prwn=Ekin;
       part=1;
    } else if (Ip == 8) {
       prwn=sqrt(Ekin*Ekin + 2.*0.93956563);
       part=13;
    } else {
       part=ifluge[int(Ip)-1];
       pMC->Gfpart(part, name, itrtyp,  
		   amass, charge, tlife); 
       prwn=sqrt(Ekin*Ekin + 2.*amass);
    }
    origin[0]=Xi;
    origin[1]=Yi;
    origin[2]=Zi;

    p[0]=Px*prwn;
    p[1]=Py*prwn;
    p[2]=Pz*prwn;

    //handle particle weight correctly
    wgt = (part == 13) ? Wgt*fAddWeight : Wgt;
    iwgt=Int_t(wgt);
    fwgt=wgt-Float_t(iwgt);
    pMC->Rndm(random,2);
    if (random[0] < fwgt) iwgt++;
    if (part==1 && iwgt>100) iwgt=100;
    Int_t nstack=0;
    for (j=0; j<iwgt; j++) {
	gAlice->SetTrack(1,-1,part,p,origin,polar,0,"Primary",nt);
	pMC->Rndm(random,2);
	phi=2*random[1]*TMath::Pi();
	Float_t pn1=p[0]*TMath::Sin(phi) - p[1]*TMath::Cos(phi);
	Float_t pn2=p[0]*TMath::Cos(phi) + p[1]*TMath::Sin(phi);
	p[0]=pn1;
	p[1]=pn2;
	Float_t on1=origin[0]*TMath::Sin(phi)-origin[1]*TMath::Cos(phi);
	Float_t on2=origin[0]*TMath::Cos(phi)+origin[1]*TMath::Sin(phi);
	origin[0]=on1;
	origin[1]=on2;
	nstack++;
    }
    if (nstack == 0) continue;
  }
 
    TFile *File=0;
// Get AliRun object or create it 
    if (!gAlice) {
	gAlice = (AliRun*)File->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    TTree *fAli=gAlice->TreeK();
    if (fAli) File =fAli->GetCurrentFile();
    File->cd();
}









