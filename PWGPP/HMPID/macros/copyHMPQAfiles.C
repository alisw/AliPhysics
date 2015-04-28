#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"

//______________________________________________________________________________
const char *workdir = "$HOME/QA_hmpid";


Int_t copyHMPQAfiles(char * dataType, Int_t year, char *period, char *pass, char *filename){
//  available variables:
//  dataType     e.g. data or sim
//  year         e.g. 2011
//  period       e.g. LHC13g
//  runNumber    e.g. 169123
//  pass         e.g. cpass1,pass1,passMC

  ifstream pFile;

  pFile.open( Form("%s_%i_%s_%s.data",dataType,year,period,pass) );

  if(!pFile){cout<<"DataSet: "<<Form("%s_%i_%s_%s.data",dataType,year,period,pass)<<" File not present. Check file name!!!"<<endl; return 0;}

  //Create output directories
  gSystem->Exec(Form("mkdir -p $HOME/QA_hmpid/%s/%i/%s/%s",dataType,year,period,pass));

  //download QA histos:
  Int_t number = 0;
  pFile>>number;
  cout<<number<<endl;
  //download QA histos:

  for(Int_t ir=0;ir<number;ir++){
    Int_t run=0;
    pFile>>run;
    cout<<run<<endl;


    // download QA files
    if(!gSystem->IsFileInIncludePath( Form("%s/%s/%i/%s/%s/QAresults_%i.root",workdir, dataType, year, period, pass, run) )) {

      printf("File: QAresults_%i.root does NOT exists ... DOWNLOADING ...\n",run);

      if(strcmp(dataType,"data")==0)
        gSystem->Exec(Form("alien_cp alien:///alice/data/%i/%s/000%i/%s/%s    $HOME/QA_hmpid/data/%i/%s/%s/QAresults_%i.root",year,period,run,pass,filename,year,period,pass,run ));

      if(strcmp(dataType,"sim")==0)
	gSystem->Exec(Form("alien_cp alien:///alice/sim/%i/%s/%i/%s    $HOME/QA_hmpid/QA_%s_%s/QAresults_%i.root",year,period,run,filename,period,pass,number ));

    } else printf("File: QAresults_%i.root exists\n",number);

  }

  return 1;

}
