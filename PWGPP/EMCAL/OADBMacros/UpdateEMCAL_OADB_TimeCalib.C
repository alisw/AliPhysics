//!  EMCal bad channel map creation macro
/*!
    With this macro, the 2D histogram bad channel maps can be added to OADB/EMCAL/EMCALBadChannels.root.
*/

#include "TGeoManager.h"

/*******************************************************************
 *  NOTE: Main function which needs to be adjusted for new BC maps *
 *******************************************************************/
void UpdateEMCAL_OADB_TimeCalib(const char *fileNameOADBAli="$ALICE_DATA/OADB/EMCAL/EMCALTimeCalib.root")
{
    gSystem->Load("libOADB");  
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libEMCALrec");

    const char *fileNameOADB                ="EMCALTimeCalib_temp.root";
    gSystem->Exec(Form("cp %s %s",fileNameOADBAli,fileNameOADB));

    // LHC17pq - 20180310
//     updateFile(fileNameOADB,"TimeCalibLHC17pq","LHC17opq/LHC17pqmerge_finalcalib.root","pass1",282000,282441,0);
//     updateFile(fileNameOADB,"TimeCalibLHC17o","LHC17opq/LHC17o_finalcalib.root","pass1",280282,281999,0);

    // LHC17j - 20180320
//     updateFile(fileNameOADB,"TimeCalibLHC17j","LHC17j/LHC17j_finalcalib.root","pass1",274591,274671,0);


    // the final output will be sorted by runnumber
    sortOutput(fileNameOADB);
    rebuildTimeContainer("EMCALTimeCalib.root");

}

/*******************************************************************
 *  NOTE: Update function that removes existing BC maps from old   *
 *          file and adds new BC maps at the end of the file       *
 *  NOTE: The final file is saved and directly renamed as new      *
 *      input file for the next loop (this reduces total file size)*
 *******************************************************************/
void updateFile(const char *fileNameOADB,TString arrName, TString fileNameTimeCalib, TString passname ,Int_t runMin, Int_t runMax, Int_t updateExistingMap = 0){
    printf("\n\nAdding %s maps to OADB file:\n",arrName.Data());

    gSystem->Load("libOADB");
    AliOADBContainer *con   = new AliOADBContainer("");
    con->InitFromFile(fileNameOADB, "AliEMCALTimeCalib");
    con->SetName("Old"); 

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALTimeCalib");

    Int_t runNumber  = runMin;
    Bool_t firstTime=kFALSE;

    //open file with corrections and check
    TFile *referenceFile = TFile::Open(fileNameTimeCalib.Data());//<<---change it here
    if(referenceFile==0x0){
        AliFatal("No reference file with time calibration");
        return;
    }

    //get run (period) array
    TObjArray *arrayPeriod=NULL;
    arrayPeriod=(TObjArray*)con->GetObject(runNumber);
    if(arrayPeriod==0x0){
        arrayPeriod=new TObjArray(1);
        arrayPeriod->SetName(Form("%s",arrName.Data()));
        firstTime=kTRUE;
    }
    for(int i=0;i<con->GetNumberOfEntries();i++){
        if (runMin >= con->LowerLimit(i) && runMin <= con->UpperLimit(i)){
            printf("\n!!! Found object #%d for runrange %d--%d as low run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),runMin);
        } else if (runMax >= con->LowerLimit(i) && runMax <= con->UpperLimit(i)){
            printf("\n!!! Found object #%d for runrange %d--%d as high run limit %d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),runMax);
        } else if (runMin <= con->LowerLimit(i) && runMax >= con->UpperLimit(i)){
            printf("\n!!! Found object #%d for runrange %d--%d as full run range %d--%d is contained\n", con->GetObjectByIndex(i), con->LowerLimit(i),con->UpperLimit(i),runMin,runMax);
        }else{
            con2->AddDefaultObject(con->GetObjectByIndex(i));
            con2->AppendObject(con->GetObjectByIndex(i),con->LowerLimit(i),con->UpperLimit(i));
        }
    }

    //create pass array
    TObjArray *arrayPass=new TObjArray(8);
    arrayPass->SetName(passname.Data());

    //add histograms to pass array
    TH1D *tmpRefRun=NULL;
    for(Int_t iBC=0;iBC<4;iBC++) {//high gain
        tmpRefRun = (TH1D*)referenceFile->Get(Form("hAllTimeAvBC%d",iBC));
        if(tmpRefRun==0x0) continue;
        arrayPass->Add(tmpRefRun);
    }
    if(runNumber>200000){//Low and high gain calibration startsfrom LHC15 periods
        for(Int_t iBC=0;iBC<4;iBC++) {//low gain
        tmpRefRun = (TH1D*)referenceFile->Get(Form("hAllTimeAvLGBC%d",iBC));
        if(tmpRefRun==0x0) continue;
        arrayPass->Add(tmpRefRun);
        }
    }

    //add pass array to period
    arrayPeriod->Add(arrayPass);

    if(firstTime){
        //add arrayperiod to main oadb
        con2->AddDefaultObject((TObject*)arrayPeriod);
        //Establishing run number with the correct objects
        con2->AppendObject(arrayPeriod,runMin,runMax);
        firstTime=kFALSE;
    }

    con2->WriteToFile("EMCALTimeCalib_temp.root");
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
    AliOADBContainer *con                       =(AliOADBContainer*)f->Get("AliEMCALTimeCalib");
    con->SetName("Old");

    Int_t indexAdd                              = 0;
    Int_t largerthan                            = 0;
    Int_t currentvalue                          = 0;

    AliOADBContainer *con2                      = new AliOADBContainer("AliEMCALTimeCalib");
    // First entry needs to be added before sorting loop
//     con2->AddDefaultObject(con->GetObjectByIndex(0));
//     con2->AppendObject(con->GetObjectByIndex(0),con->LowerLimit(0),con->UpperLimit(0));
    Int_t lowestRunIndex = con->GetIndexForRun(114737);
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

    con2->WriteToFile("EMCALTimeCalib.root");
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
  AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALTimeCalib"));
  delete reader;

  AliOADBContainer *newcont = new AliOADBContainer("AliEMCALTimeCalib");
  for(int irun = 0; irun < cont->GetNumberOfEntries(); irun++){
    newcont->AppendObject(CreatePeriodContainer(static_cast<TObjArray *>(cont->GetObjArray()->At(irun))), cont->LowerLimit(irun), cont->UpperLimit(irun));
  }

  newcont->WriteToFile("EMCALTimeCalibOADBfix.root");

  TFile *reader = TFile::Open("EMCALTimeCalibOADBfix.root", "READ");
    AliOADBContainer *cont = static_cast<AliOADBContainer *>(reader->Get("AliEMCALTimeCalib"));
    delete reader;
    delete cont;
}
