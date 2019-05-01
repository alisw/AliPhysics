/// \file CreateEMCAL_OADB_RunByRunTimeCalib.C
/// \ingroup EMCAL_OADB
/// \brief Create OADB Container for EMCal Run by Run calibration dependent on Time L1 phase 
///
/// Histograms with L1Phases are loaded and some TObjArrays are filled with these histograms.
/// This UpdateEMCAL_OADB_RunByRunTimeCalib.C updates the information of a original OADB file 
/// and writes the output to EMCALTimeL1PhaseCalib_update.root
///
/// ******* Example to Create or Read OADB Container 
/// ******* for EMCal Run by Run calibration dependent 
/// ******* on Time L1 phase 
/// \author Adam Matyja, <adam.tomasz.matyja@cern.ch>, IFJ PAN Cracow
///

/// Create or read EMCAL OADB file with L1 phases
///
/// \param opt: Int_t, 1 - create, 0 - read
/// \param runNumber: Int_t, run number
/// \param filename: TString, input file with constants for runs
/// \param runlist: TString, list file with runs for calibration
/// \param passname: TString, pass name
///
void CreateEMCAL_OADB_RunByRunTimeCalib(Int_t opt = 1, Int_t runNumber = 236967,
					TString filename="CorrectionFiles/ReferenceSM_LHC15j_calib_v2.root",
					TString runlist="CorrectionFiles/runlist_LHC15j_calib.txt",
					TString passname="pass1")
{
  //to create run it:
  //.x CreateEMCAL_OADB_RunByRunTimeCalib.C(1,0,"CorrectionFiles/ReferenceSM_LHC15j_calib_v2.root","CorrectionFiles/runlist_LHC15j_calib.txt","pass0")
  //
  // to read run it:
  //.x CreateEMCAL_OADB_RunByRunTimeCalib.C(0,236967,"EMCALTimeL1PhaseCalib.root","","pass0")

  if(opt == 0) Read(runNumber,passname,filename);
  if(opt == 1) Create(filename,runlist,passname);
}

/// Create EMCAL OADB file with L1 phases
///
/// \param filename: TString, input file with constants for runs
/// \param runlist: TString, list file with runs for calibration
/// \param passname: TString, pass name
///
void Create(TString filename="CorrectionFiles/ReferenceSM_LHC15j_calib_v2.root",
	    TString runlist="CorrectionFiles/runlist_LHC15j_calib.txt",
	    TString passname="pass1")
{
  //Create OADB container for Time L1 phase calibration parameters

  gSystem->Load("libOADB");            
  AliOADBContainer* con = new AliOADBContainer("AliEMCALTimeL1PhaseCalib");
      


  // Get the list of run numbers to be added to the OADB, parameters provided usually in a 
  // root file per period
  
  ifstream fList;
  //fList.open("CorrectionFiles/runlist_LHC15j_calib.txt");//<<---change it here
  fList.open(runlist.Data());//<<---change it here
  
  Int_t runNumber  = 0;
  TString string;
  Int_t nRuns=0;
  Int_t nSM = 20;

  //open file with corrections and check
  //TFile *referenceFile = TFile::Open("CorrectionFiles/ReferenceSM_LHC15j_calib_v2.root");//<<---change it here
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
	// Access class that contains methods to read the content of 
	// the calibration file per run
//	AliEMCALCalibTimeDepCorrection  *corr =  new AliEMCALCalibTimeDepCorrection();
//	corr->ReadRootInfo(Form("CorrectionFiles/Run%d_Correction.root",runNumber));

	tmpRefRun = (TH1C*)referenceFile->Get(Form("h%d",runNumber));
	if(tmpRefRun==0x0) continue;

	//create period array
	TObjArray *arrayPeriod=new TObjArray(1);
	//arrayPeriod->SetName("LHC15j");//<<-----change here
	arrayPeriod->SetName(Form("%d",runNumber));//<<-----change here
	//create pass array
	TObjArray *arrayPass=new TObjArray(1);
	//arrayPass->SetName("muon_calo_pass1");//<<-----change here
	arrayPass->SetName(passname.Data());//<<-----change here

	// Init the histogram
//	TH1C *h = new TH1C(Form("hh%d",runNumber),"",nSM-1,0,nSM-1);
//	
//	for(Int_t ism = 0; ism < nSM; ism++) {
//	  Short_t recalFactor = tmpRefRun->GetBinContent(ism);
//	  h->SetBinContent(ism,(Short_t)(recalFactor));
//	}

	//add histogram to pass array
	//	arrayPass->Add(h);
	arrayPass->Add(tmpRefRun);

	//add pass array to period
	arrayPeriod->Add(arrayPass);

	//add arrayperiod to main oadb
	//	con->AddDefaultObject(h);
	//con->AddDefaultObject((TObject*)arrayPass);
	con->AddDefaultObject((TObject*)arrayPeriod);
      
	//Establishing run number with the correct objects
	//con->AppendObject(h,runNumber,runNumber);
	//con->AppendObject(arrayPass,runNumber,runNumber);
	con->AppendObject(arrayPeriod,runNumber,runNumber);
	
	//	tmpRefRun->Delete();
	
	nRuns++;
      }
    }
  }
  fList.close();
  printf(" *** nRuns ***  %d\n",nRuns);
  con->WriteToFile("EMCALTimeL1PhaseCalib.root");   
  printf(" written to file\n");

  tmpRefRun->Delete();

  referenceFile->Close();    
}

/// Read EMCAL OADB file with L1 phases
///
/// \param runNumber: Int_t, run number
/// \param passname: TString, pass name
/// \param filename: TString, input file with constants for runs
///
void Read(Int_t runNumber = 236967,TString passname="pass1",TString filename="EMCALTimeL1PhaseCalib.root")
{
  
  gSystem->Load("libOADB");            
  
  AliOADBContainer *cont=new AliOADBContainer("");
  //  cont->InitFromFile("$ALICE_ROOT/OADB/EMCAL/EMCALTimeL1PhaseCalib.root", "AliEMCALTimeL1PhaseCalib");
  cont->InitFromFile(filename.Data(), "AliEMCALTimeL1PhaseCalib");

  cout<<"_________--------------- dump ---------------------___________"<<endl;
  cont->Dump();
  
  cout<<"cont->GetDefaultList()->Print()"<<endl;
  cont->GetDefaultList()->Print();
  
  //TH1C *h=cont->GetObject(runNumber); //GetObject(int runnumber)
  TObjArray *recal=(TObjArray *)cont->GetObject(runNumber); //GetObject(int runnumber)
  recal->ls();

  //TObjArray *recalpass=(TObjArray *)recal->FindObject("muon_calo_pass1");
  TObjArray *recalpass=(TObjArray *)recal->FindObject(passname.Data());
  if(!recalpass){
    cout<<" no pass"<<passname.Data() <<endl;
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


