//This a modification in a code sent by Marco Bregant, which, originally, created the OADB for misalignment matrices

//In this macro, the histograms with Recalibraton Factors are loaded and some TObjarrays are filled with these histograms.
// At the end, a OADB container is created receiving these arrays.
// This UpdateEMCAL_OADB_Recalib updates the information of a original OADB file and writes the output to BetaRecalib.root


void UpdateEMCAL_OADB_Recalib(const char *fileNameOADB="$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root",
			      const char *fileName12="RecalDB/RecalibrationFactors2012_10SM_iter8.root")
{

gSystem->Load("libOADB");  

Bool_t is2012=1;
Bool_t is2013=1;

AliOADBContainer *con	= new AliOADBContainer("");
con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); //Updating the original OADB file, output will be written into BetaRecalib.root 

// **** Loading the root files with Recalibration Factors:
TFile* f12=new TFile(fileName12);

TObjArray *array12_13 = new TObjArray(10); // 2012--2013 ---> Same R.F. for both 2012 and 1013 pass1.
array12_13->SetName("Recalib");

char name[30];


// Filling The objects above with the EMCALRecalFactors_SM Histos:
for (Int_t mod=0;mod<10;mod++){
    cout<<"SM "<< mod<<endl;
    // Recalib Objects for 2012: Still using 10 SM's
    sprintf(name,"EMCALRecalFactors_SM%d",mod);
    cout<<"Recalib2012 and 2013:"<<name<<endl;
    array12_13->Add(f12->Get(name));
    } //mod

// So far, SM11 and SM12 receive 1. 
TH2F *h0  = (TH2F*)f12->Get("EMCALRecalFactors_SM0");
TH2F *h10 = (TH2F*)h0->Clone("EMCALRecalFactors_SM10");
TH2F *h11 = (TH2F*)h0->Clone("EMCALRecalFactors_SM11");
h10->SetName("EMCALRecalFactors_SM10");
h10->SetTitle("EMCALRecalFactors_SM10");
h11->SetName("EMCALRecalFactors_SM11");
h11->SetTitle("EMCALRecalFactors_SM11");
int nbinsx = h10->GetNbinsX();
int nbinsy = h10->GetNbinsY();
for(int i=0;i<nbinsx;i++)
  for(int j=0;j<nbinsy;j++){
      h10->SetBinContent(i,j,1.);
      h11->SetBinContent(i,j,1.);
  }
if(is2012||is2013) array12_13->Add(h10);
if(is2012||is2013) array12_13->Add(h11);

//********************************************************************
// ******* 2012 -- 2013 ******************
// Setting Periods
//**** Adding pass object to period Object ****/
//When updating object that has already been created. For instance, adding pass2,3 etc.
//Just get the object and add new array. Append of runnumber is already done in this case.

TObjArray *array13b = (TObjArray*)con->GetObject(195345,"LHC13b");
TObjArray *array13c = (TObjArray*)con->GetObject(195529,"LHC13c");          
TObjArray *array13d = (TObjArray*)con->GetObject(195681,"LHC13d");
TObjArray *array13e = (TObjArray*)con->GetObject(195935,"LHC13e");
TObjArray *array13f = (TObjArray*)con->GetObject(196433,"LHC13f");


TObjArray *array13b_pass3 = new TObjArray(10); array13b_pass3->SetName("pass3");  array13b_pass3->Add(*&array12_13);
TObjArray *array13c_pass2 = new TObjArray(10); array13c_pass2->SetName("pass2");  array13c_pass2->Add(*&array12_13);
TObjArray *array13d_pass2 = new TObjArray(10); array13d_pass2->SetName("pass2");  array13d_pass2->Add(*&array12_13);
TObjArray *array13e_pass2 = new TObjArray(10); array13e_pass2->SetName("pass2");  array13e_pass2->Add(*&array12_13);
TObjArray *array13f_pass2 = new TObjArray(10); array13f_pass2->SetName("pass2");  array13f_pass2->Add(*&array12_13);

array13b->Add(*&array13b_pass3);
array13c->Add(*&array13c_pass2);
array13d->Add(*&array13d_pass2);
array13e->Add(*&array13e_pass2);
array13f->Add(*&array13f_pass2);

con->WriteToFile("BetaRecalib.root");

test(195935); // If someone wants to test container

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
