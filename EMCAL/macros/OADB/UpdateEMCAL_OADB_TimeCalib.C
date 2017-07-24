///
/// \file UpdateEMCAL_OADB_TimeCalib.C
/// \ingroup EMCAL_OADB
/// \brief Update OADB file with EMCal Time calibration factors for different BCs
///
/// Histograms with time calibration offsets are loaded and some TObjArrays are filled with these histograms.
/// This UpdateEMCAL_OADB_TimeCalib.C updates the information of a original OADB file EMCALTimeCalib.root
/// and writes the output to EMCALTimeCalib_update1.root
/// Remember to change EMCALTimeCalib_updateX.root to EMCALTimeCalib.root when update done
///
/// ******* Example to Update OADB Container 
/// ******* for EMCal Time calibration
/// ******* for different BCs
///// update LHC15ijno and LHC16hfg
///  .x UpdateEMCAL_OADB_TimeCalib.C("OADBfiles_v0_notModified/EMCALTimeCalib.root","CorrectionFilesTime/Reference_LHC15j_calib_final.root",236967,238622,"TimeCalib15j","pass0","update")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update.root","CorrectionFilesTime/Reference_LHC15i_calib.root",235709,236849,"TimeCalib15i","pass0","update2")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update2.root","CorrectionFilesTime/Reference_LHC15o_final.root",244917,246994,"TimeCalib15o","pass1","update3")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update3.root","CorrectionFilesTime/Reference_LHC15n_final.root",244340,244628,"TimeCalib15n","pass1","update4")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update4.root","CorrectionFilesTime/Reference_LHC15i_calib_mcp2_final.root",235709,236849,"TimeCalib15i","pass1","update5")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update5.root","CorrectionFilesTime/Reference_LHC15j_calib_mcp2_final.root",236967,238622,"TimeCalib15j","pass1","update6")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update6.root","CorrectionFilesTime/Reference_LHC16h_mcp1_novdM_calib_final.root",254381,255467,"TimeCalib16h","pass1","update7")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update7.root","CorrectionFilesTime/Reference_LHC16f_final.root",253658,253961,"TimeCalib16f","pass1","update8")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update8.root","CorrectionFilesTime/Reference_LHC16g_final.root",254124,254332,"TimeCalib16g","pass1","update9")
///
///// update LHC10b-f periods pass4
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update9.root","CorrectionFilesTime/Reference_LHC10b_calib.root",114737,117223,"TimeCalib10b","pass4","update10")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update10.root","CorrectionFilesTime/Reference_LHC10c_calib.root",118359,121040,"TimeCalib10c","pass4","update11")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update11.root","CorrectionFilesTime/Reference_LHC10d_calib.root",122195,126437,"TimeCalib10d","pass4","update12")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update12.root","CorrectionFilesTime/Reference_LHC10e_calib.root",127102,130850,"TimeCalib10e","pass4","update13")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update13.root","CorrectionFilesTime/Reference_LHC10f_calib.root",133004,135031,"TimeCalib10f","pass4","update14")
///
///// update LHC16klqrst only pass1
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update14.root","CorrectionFilesTime/Reference_LHC16k_final.root",256504,258574,"TimeCalib16k","pass1","update15")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update15.root","CorrectionFilesTime/Reference_LHC16l_final.root",258883,260187,"TimeCalib16l","pass1","update16")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update16.root","CorrectionFilesTime/Reference_LHC16q_final.root",265015,265525,"TimeCalib16q","pass1","update17")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update17.root","CorrectionFilesTime/Reference_LHC16r_final.root",265630,266318,"TimeCalib16r","pass1","update18")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update18.root","CorrectionFilesTime/Reference_LHC16s_final.root",266405,267131,"TimeCalib16s","pass1","update19")
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update19.root","CorrectionFilesTime/Reference_LHC16q_final.root",267161,267166,"TimeCalib16t","pass1","update20")
/// \author Adam Matyja, <adam.tomasz.matyja@cern.ch>, IFJ PAN Cracow
///



///
/// update the root file
//
/// \param fileNameOADB: TString, input OADB file
/// \param filename: TString, input file with constants for given period
/// \param runMin: Int_t, minimum run of validity for calib. constants
/// \param runMax: Int_t, maximum run of validity for calib. constants
/// \param periodname: TString, period name
/// \param passname: TString, pass name
/// \param outputPrefix: TString, output prefix
///
void UpdateEMCAL_OADB_TimeCalib(TString fileNameOADB="$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root",
				TString filename="CorrectionFilesTime/Reference_LHC15j_calib_final.root",
				Int_t runMin, Int_t runMax, TString periodname="TimeCalib15j",
				TString passname="pass1", TString outputPrefix="update1")
{
  gSystem->Load("libOADB");  
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB.Data(), "AliEMCALTimeCalib");

  Int_t runNumber  = runMin;
  Bool_t firstTime=kFALSE;

  //open file with corrections and check
  TFile *referenceFile = TFile::Open(filename.Data());//<<---change it here
  if(referenceFile==0x0){
    AliFatal("No reference file with time calibration");
    return;
  } 

  //get run (period) array
  TObjArray *arrayPeriod=NULL;
  //cout<<"arrayPeriod (null) "<<arrayPeriod<<endl;
  arrayPeriod=(TObjArray*)con->GetObject(runNumber);
  //cout<<"arrayPeriod (still null) "<<arrayPeriod<<endl;
  if(arrayPeriod==0x0){//not exist in OADB need to be created
    arrayPeriod=new TObjArray(1);
    arrayPeriod->SetName(Form("%s",periodname.Data()));
    firstTime=kTRUE;
  }
  //cout<<"arrayPeriod (not null)"<<arrayPeriod<<endl;
  
  //create pass array
  TObjArray *arrayPass=new TObjArray(8);
  arrayPass->SetName(passname.Data());
  
  //add histograms to pass array
  TH1D *tmpRefRun=NULL;

  if(runNumber>200000){//Low and high gain calibration starts from LHC15 periods
    for(Int_t iBC=0;iBC<4;iBC++) {//high gain
      tmpRefRun = (TH1D*)referenceFile->Get(Form("hAllTimeAvBC%d",iBC));
      if(tmpRefRun==0x0) continue;
      arrayPass->Add(tmpRefRun);
    }

    for(Int_t iBC=0;iBC<4;iBC++) {//low gain
      tmpRefRun = (TH1D*)referenceFile->Get(Form("hAllTimeAvLGBC%d",iBC));
      if(tmpRefRun==0x0) continue;
      arrayPass->Add(tmpRefRun);
    }
  } else {//for Run1 calibration
    for(Int_t iBC=0;iBC<4;iBC++) {//high gain gistograms are stored in low gain histos in new 2015 calibration class
      tmpRefRun = (TH1D*)referenceFile->Get(Form("hAllTimeAvLGBC%d",iBC));
      tmpRefRun->SetNameTitle(Form("hAllTimeAvBC%d",iBC),Form("hAllTimeAvBC%d",iBC));
      if(tmpRefRun==0x0) continue;
      arrayPass->Add(tmpRefRun);
    }

  }
  
  //add pass array to period
  arrayPeriod->Add(arrayPass);
  
  //When updating object that has already been created: for instance, adding pass2,3 etc.
  //Just get the object and add new array. Append of runnumber is already done in this case.
  
  if(firstTime){
    //add arrayperiod to main oadb
    con->AddDefaultObject((TObject*)arrayPeriod);
    //Establishing run number with the correct objects
    con->AppendObject(arrayPeriod,runMin,runMax);
    firstTime=kFALSE;
  }

  con->WriteToFile(Form("EMCALTimeCalib_%s.root",outputPrefix.Data()));   
  printf(" written to file\n");
  
  tmpRefRun->Delete();
  
  referenceFile->Close();    
}

/// \param runnumber: Int_t, run number
/// \param passname: TString, pass name
///
void Read(Int_t runnumber=236967,TString passname="pass1" ){
  //
  // let's read back the file
  //
  AliOADBContainer *cont=new AliOADBContainer("");
  cont->InitFromFile(Form("EMCALTimeCalib_update.root"), "AliEMCALTimeCalib");
  // 
  cout<<"_________--------------- dump ---------------------___________"<<endl;
  cont->Dump();
  cout<<"_________--------------- list ---------------------___________"<<endl;
  //cont0.List();
  cout<<"cont->GetDefaultList()->Print()"<<endl;
  cont->GetDefaultList()->Print();

  TObjArray *recal=cont->GetObject(runnumber); //GetObject(int runnumber)
  recal->ls();

  TObjArray *recalpass=(TObjArray *)recal->FindObject(passname.Data());
  if(!recalpass){
    cout<<" no pass "<<passname.Data() <<endl;
    return;
  }

  TH1D *h=(TH1D *)recalpass->FindObject(Form("hAllTimeAvBC0"));
  if (h) {
    printf("runNumber %d found\n", runNumber);
  }
  else {
    printf("runNumber %d not found - Error - do not use this run - not calibrated\n", runNumber);
    return;
  }
  h->Print(); // tmp debug  

  // Read parameter file line-by-line  
  // Get number of lines first
  
//  Int_t nSM = 20;
//  
//  for(Int_t iSM = 0; iSM < nSM; iSM++)
//  {
//    printf("SM %d, content %d\n",iSM,h->GetBinContent(iSM));
//  }
  
  h->Draw();
}
