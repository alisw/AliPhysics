//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/EMCALBadChannels.root.
*/

#include "TGeoManager.h"

/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new BC maps *
 *******************************************************************/
void UpdateEMCAL_OADB_RunByRunTimeCalib(const char *fileNameOADBAli="$ALICE_DATA/OADB/EMCAL/EMCALTimeL1PhaseCalib.root")
{
    gSystem->Load("libOADB");  
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libEMCALrec");

    const char *fileNameOADB                ="EMCALTimeL1PhaseCalib_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // LHC17pq - 20180310
//      updateFile(fileNameOADB,"TimeCalibL1LHC17o","LHC17opq/LHC17o_mcp1_L1phases.root","LHC17opq/runlist_17o.txt","pass1");
//      updateFile(fileNameOADB,"TimeCalibL1LHC17p","LHC17opq/LHC17p_mcp1_L1phases.root","LHC17opq/runlist_17p.txt","pass1");
//      updateFile(fileNameOADB,"TimeCalibL1LHC17q","LHC17opq/LHC17q_mcp1_L1phases.root","LHC17opq/runlist_17q.txt","pass1");
    // LHC17j - 20180320
//      updateFile(fileNameOADB,"TimeCalibL1LHC17j","LHC17j/LHC17j_mcp1_L1phases.root","LHC17j/runlist_17j.txt","pass1");


    // the final output will be sorted by runnumber
    sortOutput("EMCALTimeL1PhaseCalib_temp.root");
    rebuildTimeContainer("EMCALTimeL1PhaseCalib.root");

}

/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString filename, TString runlist, TString passname){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    gSystem->Load("libOADB");  
  AliOADBContainer *con    = new AliOADBContainer("");
  con->InitFromFile(fileNameOADB, "AliEMCALTimeL1PhaseCalib");

  //list with run numbers to update
  ifstream fList;
  fList.open(runlist.Data());

  Int_t runNumber  = 0;
  TString string;
  Int_t nRuns=0;
  Int_t nSM = 20;

  Bool_t firstTime=kFALSE;

  //open file with corrections and check
  TFile *referenceFile = TFile::Open(filename.Data());//<<---change it here
  if(referenceFile==0x0){
    AliFatal("No reference file with L1 phases run by run");
    return;
  } 

  TH1C *tmpRefRun=NULL;

  if (fList.good()) 
  {

    while( string.ReadLine(fList, kFALSE) ) 
    {
      sscanf(string.Data(), "%d",&runNumber);
      
      if     (runNumber < 140000) nSM = 4;
      else if(runNumber < 200000) nSM = 10;
      else nSM = 20;

      if(runNumber>200000){//L1 phase is only in LHC15 periods

    tmpRefRun = (TH1C*)referenceFile->Get(Form("h%d",runNumber));
    if(tmpRefRun==0x0) continue;

    //get run (period) array
    TObjArray *arrayPeriod=NULL;
    arrayPeriod=(TObjArray*)con->GetObject(runNumber);
    if(arrayPeriod==0x0){//not exist in OADB need to be created
      arrayPeriod=new TObjArray(1);
      arrayPeriod->SetName(Form("%d",runNumber));
      firstTime=kTRUE;
    }
    //cout<<"arrayPeriod (not null)"<<arrayPeriod<<endl;

    //create pass array
    TObjArray *arrayPass=new TObjArray(1);
    arrayPass->SetName(passname.Data());


    //add histogram to pass array
    //arrayPass->Add(h);
    arrayPass->Add(tmpRefRun);
    
    //add pass array to period
    arrayPeriod->Add(arrayPass);

    //When updating object that has already been created: for instance, adding pass2,3 etc.
    //Just get the object and add new array. Append of runnumber is already done in this case.

    if(firstTime){
      //add arrayperiod to main oadb
      con->AddDefaultObject((TObject*)arrayPeriod);
      //Establishing run number with the correct objects
      con->AppendObject(arrayPeriod,runNumber,runNumber);
      firstTime=kFALSE;
    }
    
    nRuns++;
      }
    }
  }
  fList.close();
  printf(" *** nRuns ***  %d\n",nRuns);
  con->WriteToFile("EMCALTimeL1PhaseCalib_temp.root");   
  //con->WriteToFile("EMCALTimeL1PhaseCalib.root");   
  printf(" written to file\n");

  tmpRefRun->Delete();

  referenceFile->Close();    
}


/*******************************************************************
 *  NOTE: Sorting function to sort the final OADB file             *
 *                  by ascending runnumber                         *
 *******************************************************************/
void sortOutput(const char *fileNameOADB=""){

    TFile *f                                    = TFile::Open(fileNameOADB);
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALTimeL1PhaseCalib");
    con->SetName("Old");

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALTimeL1PhaseCalib");
    // First entry needs to be added before sorting loop
//     con2->AddDefaultObject(con->GetObjectByIndex(0));
//     con2->AppendObject(con->GetObjectByIndex(0),con->LowerLimit(0),con->UpperLimit(0));
    Int_t lowestRunIndex = con->GetIndexForRun(235716);
//     con2->AddDefaultObject(con->GetObjectByIndex(lowestRunIndex));
    con2->AppendObject(con->GetObjectByIndex(lowestRunIndex),con->LowerLimit(lowestRunIndex),con->UpperLimit(lowestRunIndex));
    // sorting magic happens here
    for(int i=1;i<con->GetNumberOfEntries();i++){
        largerthan                              = con2->UpperLimit(con2->GetNumberOfEntries()-1);
        currentvalue                            = -1;
        indexAdd                                = 0;
        for(int j=0;j<con->GetNumberOfEntries();j++){
            if(con->UpperLimit(j)<=largerthan) 
                continue;
            else if(currentvalue < 0){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
            else if(con->UpperLimit(j)<currentvalue){
                currentvalue                    = con->UpperLimit(j);
                indexAdd                        = j;
            }
        }
//         con2->AddDefaultObject(con->GetObjectByIndex(indexAdd));
        con2->AppendObject(con->GetObjectByIndex(indexAdd),con->LowerLimit(indexAdd),con->UpperLimit(indexAdd));
    }

    printf("\n\n");
    Int_t nentries2                             = con2->GetNumberOfEntries();
    for(int i=0;i<nentries2;i++){
        printf("\n Entry2 --> %d/%d -->",i,nentries2);
        printf("%d -- %d --> obj = %p , %s", con2->LowerLimit(i),con2->UpperLimit(i),con2->GetObjectByIndex(i),con2->GetObjectByIndex(i)->GetName());
    }
    printf("\n\n");

    con2->WriteToFile("EMCALTimeL1PhaseCalib.root");
    gSystem->Exec(Form("rm %s",fileNameOADB));

}


TObjArray *CreatePeriodContainer(TObjArray *inputcont){
  TObjArray *newcont = new TObjArray(inputcont->GetEntries());
  newcont->SetName(inputcont->GetName());
  for(int i = 0; i < inputcont->GetEntries(); i++){
    newcont->AddAt(inputcont->At(i)->Clone(), i);
  }
  return newcont;
}

/*******************************************************************
 *  NOTE: Function required to fix OADB ownership                  *
 *                                                                 *
 *******************************************************************/
void rebuildTimeContainer(const char *fileNameOADB=""){
  TFile *reader = TFile::Open(fileNameOADB);
  AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALTimeL1PhaseCalib"));
  delete reader;

  AliOADBContainer *newcont = new AliOADBContainer("AliEMCALTimeL1PhaseCalib");
  for(int irun = 0; irun < cont->GetNumberOfEntries(); irun++){
    newcont->AppendObject(CreatePeriodContainer(static_cast<TObjArray *>(cont->GetObjArray()->At(irun))), cont->LowerLimit(irun), cont->UpperLimit(irun));
  }

  newcont->WriteToFile("EMCALTimeL1PhaseCalibOADBfix.root");

  TFile *reader = TFile::Open("EMCALTimeL1PhaseCalibOADBfix.root", "READ");
    AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALTimeL1PhaseCalib"));
    delete reader;
    delete cont;
}
