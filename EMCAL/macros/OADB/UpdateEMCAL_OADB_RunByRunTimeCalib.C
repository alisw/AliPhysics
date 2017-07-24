/// \file UpdateEMCAL_OADB_RunByRunTimeCalib.C
/// \ingroup EMCAL_OADB
/// \brief  Update OADB Container for EMCal Run by Run calibration dependent on Time L1 phase 
///
/// Histograms with L1Phases are loaded and some TObjArrays are filled with these histograms.
/// This UpdateEMCAL_OADB_RunByRunTimeCalib.C updates the information of a original OADB file 
/// and writes the output to EMCALTimeL1PhaseCalib_update.root
///
/// ******* Example to Update OADB Container 
/// ******* for EMCal Run by Run calibration dependent 
/// ******* on Time L1 phase 
/// to run it:
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib.root","CorrectionFiles/ReferenceSM_LHC15i_calib_v2.root","CorrectionFiles/runlist_LHC15i_calib.txt","pass0","EMCALTimeL1PhaseCalib_update.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update.root","CorrectionFiles/ReferenceSM_LHC15j_muon_calo_pass1_v1.root","CorrectionFiles/runlist_LHC15j.txt","pass0","EMCALTimeL1PhaseCalib_update2.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update2.root","CorrectionFiles/ReferenceSM_LHC15n_v1.root","CorrectionFiles/runlist_LHC15n.txt","pass1","EMCALTimeL1PhaseCalib_update3.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update3.root","CorrectionFiles/ReferenceSM_LHC15o_muon_calo_pass1_v1.root","CorrectionFiles/runlist_LHC15o_mcp1.txt","pass1","EMCALTimeL1PhaseCalib_update4.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update4.root","CorrectionFiles/ReferenceSM_LHC15i_calib_muon_calo_pass2_v1.root","CorrectionFiles/runlist_LHC15i_calib.txt","pass1","EMCALTimeL1PhaseCalib_update5.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update5.root","CorrectionFiles/ReferenceSM_LHC15j_calib_muon_calo_pass2_v1.root","CorrectionFiles/runlist_LHC15j_calib.txt","pass1","EMCALTimeL1PhaseCalib_update6.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update6.root","CorrectionFiles/ReferenceSM_LHC15j_mcp2_v1.root","CorrectionFiles/runlist_LHC15j.txt","pass1","EMCALTimeL1PhaseCalib_update7.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update7.root","CorrectionFiles/ReferenceSM_LHC16h_mcp1_pass1_v1.root","CorrectionFiles/runlist_LHC16h.txt","pass1","EMCALTimeL1PhaseCalib_update8.root")
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update8.root","CorrectionFiles/ReferenceSM_LHC16f_mcp1_v1.root","CorrectionFiles/runlist_LHC16f.txt","pass1","EMCALTimeL1PhaseCalib_update9.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update9.root","CorrectionFiles/ReferenceSM_LHC16g_mcp1_v1.root","CorrectionFiles/runlist_LHC16g.txt","pass1","EMCALTimeL1PhaseCalib_update10.root")  
///
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update10.root","CorrectionFiles/ReferenceSM_LHC16k_pass1_v2.root","CorrectionFiles/runlist_LHC16k.txt","pass1","EMCALTimeL1PhaseCalib_update11.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update11.root","CorrectionFiles/ReferenceSM_LHC16l_pass1_v2.root","CorrectionFiles/runlist_LHC16l.txt","pass1","EMCALTimeL1PhaseCalib_update12.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update12.root","CorrectionFiles/ReferenceSM_LHC16q_mcp1_v1.root","CorrectionFiles/runlist_LHC16q.txt","pass1","EMCALTimeL1PhaseCalib_update13.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update13.root","CorrectionFiles/ReferenceSM_LHC16r_mcp2_v2.root","CorrectionFiles/runlist_LHC16r.txt","pass1","EMCALTimeL1PhaseCalib_update14.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update14.root","CorrectionFiles/ReferenceSM_LHC16s_mcp1_v1.root","CorrectionFiles/runlist_LHC16s.txt","pass1","EMCALTimeL1PhaseCalib_update15.root")  
///.x UpdateEMCAL_OADB_RunByRunTimeCalib.C("EMCALTimeL1PhaseCalib_update15.root","CorrectionFiles/ReferenceSM_LHC16t_mcp1_v1.root","CorrectionFiles/runlist_LHC16t.txt","pass1","EMCALTimeL1PhaseCalib_update16.root")  
///
/// \author Adam Matyja, <adam.tomasz.matyja@cern.ch>, IFJ PAN Cracow
///

///
/// update the root file
//
/// \param fileNameOADB: TString, input OADB file
/// \param filename: TString, input file with constants for given period
/// \param runlist: TString, list file with runs for calibration
/// \param passname: TString, pass name
/// \param output: TString, output filename
void UpdateEMCAL_OADB_RunByRunTimeCalib(TString fileNameOADB="$ALICE_ROOT/OADB/EMCAL/EMCALTimeL1PhaseCalib.root",
					TString filename="CorrectionFiles/ReferenceSM_LHC15j_calib_v2.root",
					TString runlist="CorrectionFiles/runlist_LHC15j_calib.txt",
					TString passname="pass1", TString output="EMCALTimeL1PhaseCalib_update.root")
{
  gSystem->Load("libOADB");  
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB.Data(), "AliEMCALTimeL1PhaseCalib");

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
	//cout<<"arrayPeriod (null) "<<arrayPeriod<<endl;
	//arrayPeriod=(TObjArray*)con->GetObject(runNumber,Form("%d",runNumber));
	arrayPeriod=(TObjArray*)con->GetObject(runNumber);
	//cout<<"arrayPeriod (still null) "<<arrayPeriod<<endl;
	if(arrayPeriod==0x0){//not exist in OADB need to be created
	  arrayPeriod=new TObjArray(1);
	  arrayPeriod->SetName(Form("%d",runNumber));
	  firstTime=kTRUE;
	}
	//cout<<"arrayPeriod (not null)"<<arrayPeriod<<endl;

	//create pass array
	TObjArray *arrayPass=new TObjArray(1);
	arrayPass->SetName(passname.Data());

//	// Init the histogram
//	TH1C *h = new TH1C(Form("hh%d",runNumber),"",nSM-1,0,nSM-1);
//	for(Int_t ism = 0; ism < nSM; ism++) {
//	  Short_t recalFactor = tmpRefRun->GetBinContent(ism);
//	  h->SetBinContent(ism,(Short_t)(recalFactor));
//	}

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
  //con->WriteToFile("EMCALTimeL1PhaseCalib_update2.root");   
  con->WriteToFile(output.Data());   
  //con->WriteToFile("EMCALTimeL1PhaseCalib.root");   
  printf(" written to file %s\n",output.Data());

  tmpRefRun->Delete();

  referenceFile->Close();    
}

/// \param runnumber: Int_t, run number
/// \param passname: TString, pass name
///
void Read(int runnumber=236967,TString passname="pass1" ){
  //
  // let's read back the file
  //
  AliOADBContainer *cont=new AliOADBContainer("");
  cont->InitFromFile("EMCALTimeL1PhaseCalib_update.root", "AliEMCALTimeL1PhaseCalib");
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

  TH1C *h=(TH1C *)recalpass->FindObject(Form("h%d",runNumber));
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
  
  Int_t nSM = 20;
  
  for(Int_t iSM = 0; iSM < nSM; iSM++)
  {
    printf("SM %d, content %d\n",iSM,h->GetBinContent(iSM));
  }
  
  h->Draw();
}
