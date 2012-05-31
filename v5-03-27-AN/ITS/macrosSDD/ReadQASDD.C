#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TFileMerger.h>
#include <TAlienFile.h>
//#include <TExec.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <Riostream.h>
#include <TObjArray.h>
#include <TClass.h>
#endif
void ReadQASDD(Int_t runNb = 101498,Int_t year=2009,Char_t period[10]="LHC09c",Char_t pass[8]="pass1",Int_t maxfiles=300,Char_t filetosearch[50]="Merged.QA.Data.root",Char_t initfileout[50]="File.QA")
{

  //****************** Connection to alien *****************************************
  gSystem->Load("libNetx.so") ; 
  gSystem->Load("libRAliEn.so"); 
  TGrid::Connect("alien://",0,0,"t"); 
  //TGrid *gGrid = TGrid::Connect("alien"); 
  if(!gGrid||!gGrid->IsConnected()) {
    printf("gGrid not found! exit macro\n");
    return;
  }

  TFileMerger merger;
  merger.SetFastMethod(kTRUE);
  Char_t fileName[100];
  Char_t directory[100];
  sprintf(fileName,"%s.%i.%s.%s.Run.%i.root",initfileout,year,period,pass,runNb);
  merger.OutputFile(fileName);//metto il nome del file QA
  //sprintf(directory,"local://%s",gSystem->pwd());
  Char_t path[200];

  sprintf(path,"/alice/data/%04i/%s/%09i/ESDs/%s/%02i%09i*.*",year,period,runNb,pass,year-2000,runNb);
  printf("path %s\n",path);

  TGridResult *gr = gGrid->Query(path,filetosearch);
  if (gr->GetEntries() < 1) {
    printf("In this run there are not QA files: Exit macro\n");
    return;
  }
  else{printf("================>%d files found\n", gr->GetEntries());}

  Int_t mergedFiles = 0;
  Int_t nFiles = gr->GetEntries();
  if(nFiles>maxfiles) nFiles=maxfiles; 
  for (Int_t i = 3; i <nFiles ; i++) { 
    printf("File %i/%i\n",i+1,nFiles); 
    sprintf(directory,"%s",gr->GetKey(i,"turl"));
    printf("%s\n\n", directory);
       if(i==0) 
      {
	TFile *checkfile=TFile::Open(directory);
	if(checkfile->GetKey("ITS")==0x0){
	  printf("file: %s \n Run %d, In this run ITS QA has not been executed.-- Exit macro\n",directory, runNb);
	  break;
	  //mergedFiles=-1;
	}
	checkfile->Close();
	delete checkfile;
	checkfile=NULL;
      }
       if (merger.AddFile(directory))
	 mergedFiles++;
  }
// ------------------------------in this section we create a file that will contain the number of chunks to normalize later
  if(mergedFiles>0){
    FILE * pChunkNumber;
    pChunkNumber = fopen ("ChunkNumber.txt","w+");
    fprintf (pChunkNumber, "%f\n", mergedFiles);
    fclose (pChunkNumber);
    
    //printf("Add done\n");
    if(merger.Merge()==kTRUE){
      printf("merged %d files\n", mergedFiles);
      printf("output written on %s\n", fileName);
      printf("Merge done!\n");
    }
    else{printf("no files merged\n");return;}
  } 
  else{printf("no files merged\n");return;}
}
