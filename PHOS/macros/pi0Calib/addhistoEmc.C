#include "TFile.h"
#include "TList.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1.h"

void addhistoEmc(const char* list="list.txt", char *newname=0)
{
  // Add histograms from a set of files listed in file list
  // and write the result into a new file.

  char SumName[20];

  if (!newname) {
    sprintf(SumName,"Sum_All_Emc.root");
    newname=SumName;
    printf("\n  === Default output file name is %s ===\n",newname);
  }

  TFile *f=NULL;
  TH1F* hst[5][64][56];
  TH1F* hadd=NULL;
  char hnam[80]; char htit[80];
  TList* hlist = 0;

  char fpath[80];
  int ifscanf=0;
  Int_t ifile=0;

  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
        sprintf(hnam,"%d_%d_%d",iMod,iX,iZ);
        sprintf(htit,"Two-gamma inv. mass for mod %d, cell (%d,%d)",iMod,iX,iZ);
        hst[iMod][iX][iZ] = new TH1F(hnam,htit,100,0.,300.);
      }
    }
  }
  
  TH1F*  hmgg = new TH1F("hmgg","2-cluster invariant mass",100,0.,300.);
  
  FILE* fd = fopen(list,"r");
  while( ifscanf = fscanf(fd,"%s",fpath) != EOF) {
    f=new TFile(fpath);
    hlist = (TList*)f->Get("histos");
    
    for(Int_t iList=0; iList<hlist->GetEntries(); iList++) {
      hadd = (TH1F*)hlist->At(iList);
      const char* str = hadd->GetName();
      int md,X,Z;
      if (sscanf(str,"%d_%d_%d",&md,&X,&Z)) {
	hst[md][X][Z]->Add(hadd);
	//printf("Added hst[%d][%d][%d]\n",md,X,Z);
      }
      else {
	printf("Trying to add histogram %s to hmgg.\n",hadd->GetName());
	hmgg->Add(hadd);
      }
    }
    
    printf("Deleting list..\n");
    hlist->Delete();
    printf("OK!\n");
    ifile++;
   
    printf("File %s processed.\n",fpath);
    if(f) delete f;
  }
  
  printf("%d processed.\n",ifile);

  TFile outfile(newname,"recreate");

  for(Int_t iMod=0; iMod<5; iMod++) {
    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	if(hst[iMod][iX][iZ]->GetEntries()) hst[iMod][iX][iZ]->Print();
        hst[iMod][iX][iZ]->Write();
      }
    }
  }
  
  hmgg->Write();
  outfile.Close();
}
