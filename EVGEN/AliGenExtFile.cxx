#include "AliGenExtFile.h"
#include "AliGenMUONlib.h"
#include "AliMC.h"
#include "AliRun.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <stdlib.h>
 ClassImp(AliGenExtFile)
     AliGenExtFile::AliGenExtFile()
	 :AliGenerator(-1)
{
    //
    fName="ExtFile";
    fTitle="Primaries from ext. File";
    fFileName="dtujet93.root";
    fTreeNtuple=0;
    fNcurrent=0;
//
//  Read all particles
    fNpart=-1;
}

AliGenExtFile::AliGenExtFile(Int_t npart)
    :AliGenerator(npart)
{
    //
    fName="ExtFile";
    fTitle="Primaries from ext. File";
    fFileName="dtujet93.root";
    fTreeNtuple=0;
    fNcurrent=0;
}

//____________________________________________________________
AliGenExtFile::~AliGenExtFile()
{
    delete fTreeNtuple;
}

//____________________________________________________________
void AliGenExtFile::NtupleInit() 
{
//
// reset the existing file environment and open a new root file if
// the pointer to the Fluka tree is null
    
    TFile *File=0;
    if (fTreeNtuple==0) {
        if (!File) {
	    File = new TFile(fFileName);
	    File->cd();
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
        }
// get the tree address in the Fluka boundary source file
	fTreeNtuple = (TTree*)gDirectory->Get("h888");
    } else {
        File = fTreeNtuple->GetCurrentFile();
        File->cd();
    }

    TTree *h2=fTreeNtuple;
//Set branch addresses
//Set branch addresses
    h2->SetBranchAddress("Nihead",&Nihead);
    h2->SetBranchAddress("Ihead",Ihead);
    h2->SetBranchAddress("Nrhead",&Nrhead);
    h2->SetBranchAddress("Rhead",Rhead);
    h2->SetBranchAddress("Idpart",&Idpart);
    h2->SetBranchAddress("Theta",&Theta);
    h2->SetBranchAddress("Phi",&Phi);
    h2->SetBranchAddress("P",&P);
    h2->SetBranchAddress("E",&E);
}


//____________________________________________________________
void AliGenExtFile::Generate()
{

  Float_t polar[3]= {0,0,0};
  //
  Float_t origin[3]={0,0,0};
  Float_t p[3];
  Float_t random[6];
  Float_t prwn;
  char name[100];
  Float_t amass, charge, tlife;
  Int_t itrtyp;
  Int_t i, j, nt, Ntracks=0;
  //
  NtupleInit();
  TTree *h2=fTreeNtuple;
  Int_t nentries = (Int_t) h2->GetEntries();
  // loop over number of particles
  Int_t nb = (Int_t)h2->GetEvent(fNcurrent);
  Int_t i5=Ihead[4];
  Int_t i6=Ihead[5];

  for (j=0;j<3;j++) origin[j]=fOrigin[j];
  if(fVertexSmear==perEvent) {
    gMC->Rndm(random,6);
    for (j=0;j<3;j++) {
	origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }

  if (fNcurrent >= nentries) {
      printf("\n No more entries !!! !\n");
      return;
  }
  
	  
  if (i5==0) {
      printf("\n This should never happen !\n");
  } else {
      printf("\n Next event contains %d tracks! \n", i6);
      Ntracks=i6;
  }
  for (i=0; i<Ntracks; i++) {

      gMC->Gfpart(Idpart, name, itrtyp,amass, charge, tlife); 
      prwn=sqrt((E+amass)*(E-amass));

      Theta *= TMath::Pi()/180.;
      Phi    = (Phi-180)*TMath::Pi()/180.;      
      if(Theta<fThetaMin || Theta>fThetaMax ||
	 Phi<fPhiMin || Phi>fPhiMax         ||
	 prwn<fPMin || prwn>fPMax)          
      {
	  ;
      } else {
	  p[0]=prwn*TMath::Sin(Theta)*TMath::Cos(Phi);
	  p[1]=prwn*TMath::Sin(Theta)*TMath::Sin(Phi);      
	  p[2]=prwn*TMath::Cos(Theta);
	  
	  if(fVertexSmear==perTrack) {
	      gMC->Rndm(random,6);
	      for (j=0;j<3;j++) {
		  origin[j]=fOrigin[j]
		      +fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	      }
	  }
	  gAlice->SetTrack(fTrackIt,-1,Idpart,p,origin,polar,0,"Primary",nt);
      }
      fNcurrent++;
      nb = (Int_t)h2->GetEvent(fNcurrent); 
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









