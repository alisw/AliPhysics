//This a modification in a code sent by Marco Bregant, which, originally, created the OADB for misalignment matrices

//In this macro, the histograms with Recalibraton Factors are loaded and some TObjarrays are filled with these histograms.
// At the end, a OADB container is created receiving these arrays.

void CreateEMCAL_OADB_Recalib(const char *fileName10s_d="RecalDB/summer_december_2010/RecalibrationFactors.root",const char *fileName10d="RecalDB/december2010/RecalibrationFactors.root",
			      const char *fileName10s="RecalDB/summer2010/RecalibrationFactors.root", const char *fileName11v1="RecalDB/2011_v1/RecalibrationFactors.root")
{

gSystem->Load("libOADB");  
  
Bool_t Recalib2010=kTRUE;
Bool_t Recalib2011=kTRUE;


// **** Loading the root files with Recalibration Factors:
TFile* f10s_d=new TFile(fileName10s_d);
TFile* f10d=new TFile(fileName10d);
TFile* f10s=new TFile(fileName10s);
TFile* f11v1=new TFile(fileName11v1);

// ********* Arrays with the objects *********************
TObjArray array10s_d(4); //summer_december_2010/
TObjArray array10dec(4); //december2010
TObjArray array10s(4); //summer2010
TObjArray array11v1(10); // 2011_v1

array10s_d.SetName("Recalib");
array10dec.SetName("Recalib");
array10s.SetName("Recalib");
array11v1.SetName("Recalib");

char name[30];
// Filling The objects above with the EMCALRecalFactors_SM Histos:
for (Int_t mod=0;mod<10;mod++){

    cout<<"SM "<< mod<<endl;

    // Recalib Objects for 2010:
    if (Recalib2010&&mod<4)  { //if 2010 4 SMs
	    sprintf(name,"EMCALRecalFactors_SM%d",mod);
	    cout<<"Recalib2010:"<<name<<endl;

	    //Summer_December:
	    array10s_d.Add(f10s_d->Get(name));

	    //December:
	    array10dec.Add(f10d->Get(name)); 

	    //Summer:
	    array10s.Add(&array10s_d); 
	    }

    // Recalib Objects for 2011:
    if (Recalib2011)  { //if 2011 10 SMs
	    sprintf(name,"EMCALRecalFactors_SM%d",mod);
	    cout<<"Recalib2011:"<<name<<endl;
	    array11v1.Add(f11v1->Get(name));
	    }
	    
    } //mod
//********************************************************************

// ************** Establishing different configuration according to the pass ***********
// ************ Latest Tender Information : (September 3rd 2011) AliRootTrunk revision 51405 ***********************

// ******* 2010 ******************
TObjArray array10b_pass2(4);//runRC >= 114737 && runRC <= 117223 
array10b_pass2.SetName("pass2");

TObjArray array10c_pass2(4);//runRC >= 118503 && runRC <= 121040 
array10c_pass2.SetName("pass2");

TObjArray array10c_pass3(4);//runRC >= 118503 && runRC <= 121040 
array10c_pass3.SetName("pass3");

TObjArray array10d_pass1(4);//runRC >= 122195 && runRC <= 126437
array10d_pass1.SetName("pass1");

TObjArray array10e_pass1(4);//runRC >= 127712 && runRC <= 130850 
array10e_pass1.SetName("pass1");

TObjArray array10f_pass1a(4);//runRC >= 133004 && runRC <  134657
array10f_pass1a.SetName("pass1");

//December:
TObjArray array10d_pass2(4); //runRC >= 122195 && runRC <= 126437 
array10d_pass2.SetName("pass2");

TObjArray array10f_pass1b(4);//runRC >= 134657 && runRC <= 135029 
array10f_pass1b.SetName("pass1");

TObjArray array10g_pass1(4); //runRC >= 135654 && runRC <= 136377
array10g_pass1.SetName("pass1");

TObjArray array10h_pass1a(4); //Until Christmas: runRC >= 136851 && runRC < 137231  
array10h_pass1a.SetName("pass1");
	
//Summer:
TObjArray array10h_pass1b(4); //runRC >= 137231 && runRC <= 139517
array10h_pass1b.SetName("pass1");

//Summer_December:
array10b_pass2.Add(&array10s_d);//runRC >= 114737 && runRC <= 117223 
array10c_pass2.Add(&array10s_d);//runRC >= 118503 && runRC <= 121040 
array10c_pass3.Add(&array10s_d);//runRC >= 118503 && runRC <= 121040 
array10d_pass1.Add(&array10s_d);//runRC >= 122195 && runRC <= 126437
array10e_pass1.Add(&array10s_d);//runRC >= 127712 && runRC <= 130850 
array10f_pass1a.Add(&array10s_d);//runRC >= 133004 && runRC <  134657

//December:
array10d_pass2.Add(&array10dec); //runRC >= 122195 && runRC <= 126437 
array10f_pass1b.Add(&array10dec);//runRC >= 134657 && runRC <= 135029 
array10g_pass1.Add(&array10dec); //runRC >= 135654 && runRC <= 136377
array10h_pass1a.Add(&array10dec); //Until Christmas: runRC >= 136851 && runRC < 137231  
	
//Summer:
array10h_pass1b.Add(&array10s); //runRC >= 137231 && runRC <= 139517

//*********** 2011 ***************************************
TObjArray array11a_pass1(10); // LHC11a pass1
array11a_pass1.SetName("pass1");

TObjArray array11a_pass2(10); // LHC11a pass2
array11a_pass2.SetName("pass2");

TObjArray array11c_pass1(10); // LHC11c pass1
array11c_pass1.SetName("pass1");

array11a_pass1.Add(&array11v1);
array11a_pass2.Add(&array11v1);
array11c_pass1.Add(&array11v1);
// *************************************************************************************

// ******** Establishing the objects for the specific Periods ************************

// 2010
TObjArray array10b(4);
array10b.SetName("LHC10b");

TObjArray array10c(4);
array10c.SetName("LHC10c");

TObjArray array10d(4);
array10d.SetName("LHC10d");

TObjArray array10e(4);
array10e.SetName("LHC10e");

TObjArray array10fa(4);
array10fa.SetName("LHC10fa");

TObjArray array10fb(4);
array10fb.SetName("LHC10fb");

TObjArray array10g(4);
array10g.SetName("LHC10g");

TObjArray array10ha(4);
array10ha.SetName("LHC10ha");

TObjArray array10hb(4);
array10hb.SetName("LHC10hb");


array10b.Add(&array10b_pass2);

array10c.Add(&array10c_pass2);
array10c.Add(&array10c_pass3);

array10d.Add(&array10d_pass1);
array10d.Add(&array10d_pass2);

array10e.Add(&array10e_pass1);

array10fa.Add(&array10f_pass1a);
array10fb.Add(&array10f_pass1b);

array10g.Add(&array10g_pass1); 

array10ha.Add(&array10h_pass1a);
array10hb.Add(&array10h_pass1b);

// ************

//2011
TObjArray array11a(10);
array11a.SetName("LHC11a");

TObjArray array11c(10);
array11c.SetName("LHC11c");

array11a.Add(&array11a_pass1);
array11a.Add(&array11a_pass2);
array11c.Add(&array11c_pass1);
// *****************


// Creating Container 
AliOADBContainer* con = new AliOADBContainer("AliEMCALRecalib");

con->AddDefaultObject((TObject*) &array10b);
con->AddDefaultObject((TObject*) &array10c);
con->AddDefaultObject((TObject*) &array10d);
con->AddDefaultObject((TObject*) &array10e);
con->AddDefaultObject((TObject*) &array10fa);
con->AddDefaultObject((TObject*) &array10fb);
con->AddDefaultObject((TObject*) &array10g);
con->AddDefaultObject((TObject*) &array10ha);
con->AddDefaultObject((TObject*) &array10hb);
con->AddDefaultObject((TObject*) &array11a);
con->AddDefaultObject((TObject*) &array11c);

// Appending objects to their specific Run number
con->AppendObject(&array10b,114737,117223);
con->AppendObject(&array10c,118503,121040);
con->AppendObject(&array10d,122195,126437);
con->AppendObject(&array10e,127712,130850);
con->AppendObject(&array10fa,133004,134656);
con->AppendObject(&array10fb,134657,135029);
con->AppendObject(&array10g,135654,136377);
con->AppendObject(&array10ha,136851,137230);
con->AppendObject(&array10hb,137231,139517);
con->AppendObject(&array11a,144871,146860);
con->AppendObject(&array11c,151636,155384);
con->WriteToFile("BetaRecalib.root");

test(); // If someone wants to test container
}

void test(){
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

TObjArray *recal=cont->GetObject(118503); //GetObject(int runnumber)
recal->ls();

TObjArray *recalpass=recal->FindObject("pass2");

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
