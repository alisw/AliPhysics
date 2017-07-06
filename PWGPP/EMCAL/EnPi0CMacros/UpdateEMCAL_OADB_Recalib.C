
 const int kNbSMEMCAL=10;
 const int kNbSMEMCALthird=2;
 const int kNbSMDCAL=6;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;

///
/// \file UpdateEMCAL_OADB_Recalib.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Update OADB file with energy recalibration factors.
///
/// The histograms with energy recalibraton Factors are loaded and some TObjarrays 
/// are filled with these histograms. At the end, a OADB container is created 
/// receiving these arrays.
/// This UpdateEMCAL_OADB_Recalib updates the information of a original OADB file and writes the output to BetaRecalib.root///
///
/// Similar macro can be found in $ALICE_ROOT/EMCAL/macros/OADB
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS ???
/// \author Marcel Figueredo, <marcel.figueredo@cern.ch>, Sao Paulo
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Update OADB Container for EMCal energy recalibration factors
/// from external file. 
///
/// \param fileNameOADB: OADB file name and path
/// \param fileNameRecalibFactors: name and path of input file with new factors
///
void UpdateEMCAL_OADB_Recalib
(
 const char *fileNameOADB="EMCALRecalib_input4test.root",
 const char *fileNameRecalibFactors="/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data_Nov2016finalCalib/multiplyPi0CalibrationFactors_TextToHisto_Final.root"
 )
{

gSystem->Load("libOADB");  

AliOADBContainer *con	= new AliOADBContainer("");
con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); //Updating the original OADB file, output will be written into BetaRecalib.root 

// **** Loading the root files with Recalibration Factors:
TFile* fRecalibFactors=new TFile(fileNameRecalibFactors);

TObjArray *arrayRecalibFactors = new TObjArray(kNbSMtot);
arrayRecalibFactors->SetName("Recalib");

char name[30];


// Filling The objects above with the EMCALRecalFactors_SM Histos:
for (Int_t iSM=0;iSM<kNbSMtot;iSM++){
    cout<<"SM "<< iSM<<endl;
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
    cout<<"Recalib : "<<name<<endl;
    arrayRecalibFactors->Add(fRecalibFactors->Get(name));
    } //iSM


//********************************************************************
// Setting Periods
//**** Adding pass object to period Object ****/
//When updating object that has already been created. For instance, adding pass2,3 etc.
//Just get the object and add new array. Append of runnumber is already done in this case.

/*TObjArray *array13b = (TObjArray*)con->GetObject(195345,"LHC13b");
TObjArray *array13c = (TObjArray*)con->GetObject(195529,"LHC13c");          
TObjArray *array13d = (TObjArray*)con->GetObject(195681,"LHC13d");
TObjArray *array13e = (TObjArray*)con->GetObject(195935,"LHC13e");
TObjArray *array13f = (TObjArray*)con->GetObject(196433,"LHC13f");


TObjArray *array13b_pass3 = new TObjArray(kNbSMtot); array13b_pass3->SetName("pass3");  array13b_pass3->Add(*&arrayRecalibFactors);
TObjArray *array13c_pass2 = new TObjArray(kNbSMtot); array13c_pass2->SetName("pass2");  array13c_pass2->Add(*&arrayRecalibFactors);
TObjArray *array13d_pass2 = new TObjArray(kNbSMtot); array13d_pass2->SetName("pass2");  array13d_pass2->Add(*&arrayRecalibFactors);
TObjArray *array13e_pass2 = new TObjArray(kNbSMtot); array13e_pass2->SetName("pass2");  array13e_pass2->Add(*&arrayRecalibFactors);
TObjArray *array13f_pass2 = new TObjArray(kNbSMtot); array13f_pass2->SetName("pass2");  array13f_pass2->Add(*&arrayRecalibFactors);

array13b->Add(*&array13b_pass3);
array13c->Add(*&array13c_pass2);
array13d->Add(*&array13d_pass2);
array13e->Add(*&array13e_pass2);
array13f->Add(*&array13f_pass2);*/

 TObjArray *array15_pass1 = new TObjArray(kNbSMtot);
 TObjArray *array15_pass2 = new TObjArray(kNbSMtot);
 TObjArray *array15_pass3 = new TObjArray(kNbSMtot);
 array15_pass1->SetName("pass1");
 array15_pass2->SetName("pass2");
 array15_pass3->SetName("pass3");
 array15_pass1->Add(*&arrayRecalibFactors);
 array15_pass2->Add(*&arrayRecalibFactors);
 array15_pass3->Add(*&arrayRecalibFactors);

con->WriteToFile("BetaRecalib.root");

//test(195935); // If someone wants to test container

}





///
/// Test what was updated, let's read back the file
///
/// \param runnumber: reference run number
///
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






