//This a modification in a code sent by Marco Bregant, which, originally, created the OADB for misalignment matrices

//In this macro, the histograms with Recalibraton Factors are loaded and some TObjarrays are filled with these histograms.
// At the end, a OADB container is created receiving these arrays.
// This UpdateEMCAL_OADB_Recalib updates the information of a original OADB file and writes the output to BetaRecalib.root

//#include "AliOADBContainer.h"
//#include "TObjArray.h"
//#include "TFile.h"
//#include "Riostream.h"
//#include "TSystem.h"

 const int kNbSMEMCAL=10;
 const int kNbSMEMCALthird=2;
 const int kNbSMDCAL=6;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;


void UpdateEMCAL_OADB_Recalib_Run2
(const char *fileNameOADB="EMCALRecalib.root",
 const char *fileNameRecalibFactors="multiplyPi0CalibrationFactors_TextToHisto_Final.root")
{
  
  gSystem->Load("libOADB");  
  
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); //Updating the original OADB file, output will be written into BetaRecalib.root 
  
  //----------------------------------------------------------------------------
  // Move Run2 calibration that includes EMCal 2012-13 and no DCal to periods >= LHC15n
  // duplicate pass2 to pass2 and pass3
  //----------------------------------------------------------------------------

  TObjArray *array15 = (TObjArray*)con->GetObject(210000,"LHC15");
  
  con->UpdateObject(con->GetIndexForRun(210000), array15, 244340,999999);
  
  TObjArray *pass1=(TObjArray*)array15->FindObject("pass1");
  TObjArray *recal=(TObjArray*)pass1  ->FindObject("Recalib"); // needed later for first Run2 period
  
  TObjArray *array15_pass2 = new TObjArray(kNbSMtot);
  TObjArray *array15_pass3 = new TObjArray(kNbSMtot);
  TObjArray *array15_pass4 = new TObjArray(kNbSMtot);
  array15_pass2->SetName("pass2");
  array15_pass3->SetName("pass3");
  array15_pass4->SetName("pass4");
  array15_pass2->Add(*&recal);
  array15_pass3->Add(*&recal);
  array15_pass4->Add(*&recal);
  array15->Add(*&array15_pass2);
  array15->Add(*&array15_pass3);
  array15->Add(*&array15_pass4);
  
  //----------------------------------------------------------------------------
  // Create new calibration, including DCal for periods LHC15a to LHC15m
  // Recover the EMCal calibration parameters from 2012-2013, add them to the new array
  //----------------------------------------------------------------------------
  // **** Loading the root files with DCal Recalibration Factors:
  TFile* fRecalibFactors=new TFile(fileNameRecalibFactors);
  
  TObjArray *array15fm = new TObjArray(kNbSMtot);
  array15fm->SetName("LHC15fm");
  TObjArray *array15fm_pass1 = new TObjArray(kNbSMtot);
  TObjArray *array15fm_pass2 = new TObjArray(kNbSMtot);
  TObjArray *array15fm_pass3 = new TObjArray(kNbSMtot);
  array15fm_pass1->SetName("pass1");
  array15fm_pass2->SetName("pass2");
  array15fm_pass3->SetName("pass3");
  
  TObjArray *arrayRecalibFactors = new TObjArray(kNbSMtot);
  arrayRecalibFactors->SetName("Recalib");
  
  // Filling The objects above with the EMCALRecalFactors_SM Histos for DCal:
  char name[30];  
  for (Int_t iSM=10;iSM<kNbSMtot;iSM++)
  {
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
    cout<<"New SM "<<iSM<<"; Recalib : "<<name<<endl;
    arrayRecalibFactors->Add(fRecalibFactors->Get(name));
  } //iSM
    //fRecalibFactors->Close();
  
  // Add EMCal Run1 recalibration factors, factors obtained from previous update:
  for (Int_t iSM=0;iSM<10;iSM++)
  {
    arrayRecalibFactors->Add(recal->At(iSM));
    cout<<"Old SM  "<<iSM<<"; Recalib : "<<(recal->At(iSM))->GetName()<<endl;
  } //iSM
  
  array15fm_pass1->Add(*&arrayRecalibFactors);
  array15fm_pass2->Add(*&arrayRecalibFactors);
  array15fm_pass3->Add(*&arrayRecalibFactors);

  array15fm->Add(*&array15fm_pass1);
  array15fm->Add(*&array15fm_pass2);
  array15fm->Add(*&array15fm_pass3);
  
  //con->AddDefaultObject((TObject*) &array15fm);
  //con->AppendObject((TObject*) &array15fm,209122,244284);
  con->AddDefaultObject((TObject*)array15fm);
  con->AppendObject((TObject*)array15fm,209122,244284);
    
  //----------------------------------------------------------------------------

  
  con->WriteToFile("BetaRecalib.root");
  
  //test(195935); // If someone wants to test container
  
}






void test(int runnumber=195345){
//
// let's read back the file
AliOADBContainer *cont=new AliOADBContainer("");
cont->InitFromFile("BetaRecalib.root", "AliEMCALRecalib");
// 
cout<<"_________--------------- dump ---------------------___________"<<endl;
cont->Dump();
cout<<"_________--------------- list ---------------------___________"<<endl;
//cont0.List();
cout<<"cont->GetDefaultList()->Print()"<<endl;
cont->GetDefaultList()->Print();

TObjArray *recal=cont->GetObject(runnumber); //GetObject(int runnumber)
recal->ls();

TObjArray *recalpass=recal->FindObject("pass1");

if(!recalpass){
  cout<<" norecalpass"<<endl;
  return;
}

TObjArray *recalib=recalpass->FindObject("Recalib");

if(!recalib){
  cout<<" no recalib found"<<endl;
  return;
}

TH2F *h2=(TH2F*)recalib->FindObject("EMCALRecalFactors_SM0");
if(!h2){
  return;
cout<<" no histo found"<<endl;
}
h2->DrawCopy("colz");
cout<<"That's all folks!"<<endl;

  
}






