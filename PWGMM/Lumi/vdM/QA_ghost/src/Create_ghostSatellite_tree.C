#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//Global scan time variables
UInt_t gScan0StartX = 0;
UInt_t gScan0EndX   = 0;
UInt_t gScan0StartY = 0;
UInt_t gScan0EndY   = 0;
UInt_t gScan1StartX = 0;
UInt_t gScan1EndX   = 0;
UInt_t gScan1StartY = 0;
UInt_t gScan1EndY   = 0;

// Checks the variable names from the input data file
Bool_t CheckVariableName(TString name, Int_t &variable)
{
	Bool_t useVariable = kFALSE;
	TObjArray *s = name.Tokenize(":");
	TString temp = ((TObjString *)s->At(2))->GetString();
	//cout << "variable: " << temp.Data() << endl;

	variable = -1;
	//OBS: CompareTo returns 0 for identical text
	if (!temp.CompareTo("Q_ghost_bunched"))     { useVariable=kTRUE;  variable = 1;}
	if (!temp.CompareTo("Q_sat_bunched"))       { useVariable=kTRUE;  variable = 2;}
	if (!temp.CompareTo("sigma_ghost_bunched")) { useVariable=kTRUE;  variable = 3;}
	if (!temp.CompareTo("sigma_sat_bunched"))   { useVariable=kTRUE;  variable = 4;}
	if (!temp.CompareTo("sigma_ghost_bunched_stat")) { useVariable=kTRUE;  variable = 5;}
	if (!temp.CompareTo("sigma_ghost_bunched_syst")) { useVariable=kTRUE;  variable = 6;}
	if (!temp.CompareTo("sigma_sat_bunched_stat"))   { useVariable=kTRUE;  variable = 7;}
	if (!temp.CompareTo("sigma_sat_bunched_syst"))   { useVariable=kTRUE;  variable = 8;}

	if (useVariable) cout << "Use variable " << temp.Data() << endl;
	return useVariable;
}

// Returns the comma separated values from a given string
void GetVariables(TString name, UInt_t &timestamp, Double_t &value)
{
	TObjArray *s = name.Tokenize(",");
	TString tmpTime = ((TObjString *)s->At(0))->GetString();
	TString tmpValue = ((TObjString *)s->At(1))->GetString();

	//convert strings into correct variable type
	TDatime t(tmpTime.Data());
	timestamp = t.Convert();
	value = tmpValue.Atof();

	//cout << "variables " << timestamp << " and " << t.Convert() << endl;
}

// Set the start and end times of the scans
void SetScanTimes(){
  //set set branches
  Double_t time;
  g_vdm_Tree->ResetBranchAddresses();
  g_vdm_Tree->SetBranchAddress("time",&time);

  g_vdm_Tree->GetEntry(g_Idx_Start_Scan_x[0]);
  gScan0StartX = (UInt_t) time;
  cout << "Start X0: " << gScan0StartX << endl;

  g_vdm_Tree->GetEntry(g_Idx_End_Scan_x[0]);
  gScan0EndX = (UInt_t) time;
  cout << "End X0: " << gScan0EndX << endl;

  g_vdm_Tree->GetEntry(g_Idx_Start_Scan_y[0]);
  gScan0StartY = (UInt_t) time;
  cout << "Start Y0: " << gScan0StartY << endl;

  g_vdm_Tree->GetEntry(g_Idx_End_Scan_y[0]);
  gScan0EndY = (UInt_t) time;
  cout << "End Y0: " << gScan0EndY << endl;

  g_vdm_Tree->GetEntry(g_Idx_Start_Scan_x[1]);
  gScan1StartX = (UInt_t) time;
  cout << "Start X1: " << gScan1StartX << endl;

  g_vdm_Tree->GetEntry(g_Idx_End_Scan_x[1]);
  gScan1EndX = (UInt_t) time;
  cout << "End X1: " << gScan1EndX << endl;

  g_vdm_Tree->GetEntry(g_Idx_Start_Scan_y[1]);
  gScan1StartY = (UInt_t) time;
  cout << "Start Y1: " << gScan1StartY << endl;

  g_vdm_Tree->GetEntry(g_Idx_End_Scan_y[1]);
  gScan1EndY = (UInt_t) time;
  cout << "End Y1: " << gScan1EndY << endl;
}

//Check if given timestamp falls inside scan period
Bool_t InsideScanPeriod(UInt_t time){
  Bool_t isInside = kFALSE;
  if (time>=gScan0StartX && time<=gScan0EndX) isInside = kTRUE;
  if (time>=gScan0StartY && time<=gScan0EndY) isInside = kTRUE;
  if (time>=gScan1StartX && time<=gScan1EndX) isInside = kTRUE;
  if (time>=gScan1StartY && time<=gScan1EndY) isInside = kTRUE;
//  if (time>=gScan0StartX && time<=gScan1EndY) isInside = kTRUE;
  return isInside;

}

//-------------------------------------------------------
// Get the ghost and satellite charge from LHC-files
// The files include measurements for every five minuits the whole fill
// Must select only measurements falling inside scan
//-------------------------------------------------------

void Get_ghostSatellite_tree(Int_t Fill, Int_t opt)
// opt: 0 => beam 1 device 1
// opt: 1 => beam 1 device 2
// opt: 2 => beam 2
{

  // create names for ghost satellite files from martino
  // ** WARNING: Directory structure and nameing convention
  //    valid for fill 4937 ... the corresponding structure
  //    has to be added for other fills ... **
  char *data_file_name = new char[kg_string_size];
  if (opt == 0)
    sprintf(data_file_name,"../Fill-%d/ghostCharge/%d_B1.SYS_A.CSV",
	    Fill,Fill);
  else if (opt == 1)
    sprintf(data_file_name,"../Fill-%d/ghostCharge/%d_B1.SYS_B.CSV",
	    Fill,Fill);
  else
    sprintf(data_file_name,"../Fill-%d/ghostCharge/%d_B2.SYS_A.CSV",
      Fill,Fill);

  // open data file
  ifstream data_file;
  data_file.open(data_file_name);

  // now prepare names of output
  char *ghost_file_name = new char[kg_string_size];
  if (opt == 0)
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B1_A.root",
	    Fill);
  else if (opt == 1)
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B1_B.root",
	    Fill);
  else
    sprintf(ghost_file_name,"../Fill-%d/ghostCharge/ghostCharge_B2_A.root",
	    Fill);


  //Get information from data file
  UInt_t Q_ghost_time_array[1000];
  UInt_t Q_sat_time_array[1000];
  UInt_t sigma_ghost_time_array[1000];
  UInt_t sigma_sat_time_array[1000];
  Double_t Q_ghost_value_array[1000];
  Double_t Q_sat_value_array[1000];
  Double_t sigma_ghost_value_array[1000];
  Double_t sigma_ghost_stat_value_array[1000];
  Double_t sigma_ghost_syst_value_array[1000];
  Double_t sigma_sat_value_array[1000];
  Double_t sigma_sat_stat_value_array[1000];
  Double_t sigma_sat_syst_value_array[1000];
  Int_t ighost    = 0;
  Int_t ighostErr = 0;
  Int_t isat      = 0;
  Int_t isatErr   = 0;
  Int_t ighostErrStat = 0;
  Int_t ighostErrSyst = 0;
  Int_t isatErrStat = 0;
  Int_t isatErrSyst = 0;

  Bool_t fillInfo = kFALSE;
  TString line;
  Int_t variable = -1;
  while (!line.ReadLine(data_file).eof()) {
    UInt_t time = 0;
    Double_t value = 0;
    Bool_t isInsideScan = kFALSE;
    //initialize variables
    if (line.Contains("VARIABLE")){
      fillInfo = CheckVariableName(line,variable);
    }
    else if (line.Contains("DataType")|| line.Contains("Unit") || line.Contains("Length") ||line.Contains("Description") ||line.Contains("Timestamp") )
      continue;
    else if (fillInfo) {
      GetVariables(line,time,value);
      isInsideScan = InsideScanPeriod(time);
      if (isInsideScan) {
      //  cout << time << endl;
        switch (variable) {
          case 1:
              Q_ghost_value_array[ighost] = value;
              Q_ghost_time_array[ighost] = time;
                cout << line.Data() << endl;
                cout << "timestamp" << time << " value " << value << endl;
              ighost++;
          break;
          case 2:
              Q_sat_value_array[isat] = value;
              Q_sat_time_array[isat] = time;
              isat++;
          break;
          case 3:
              sigma_ghost_value_array[ighostErr] = value;
              sigma_ghost_time_array[ighostErr] = time;
              ighostErr++;
          break;
          case 4:
              sigma_sat_value_array[isatErr] = value;
              sigma_sat_time_array[isatErr] = time;
              isatErr++;
          break;
          case 5:
              sigma_ghost_stat_value_array[ighostErrStat] = value;
              ighostErrStat++;
          break;
          case 6:
              sigma_ghost_syst_value_array[ighostErrSyst] = value;
              ighostErrSyst++;
          break;
          case 7:
              sigma_sat_stat_value_array[isatErrStat] = value;
              isatErrStat++;
          break;
          case 8:
              sigma_sat_syst_value_array[isatErrSyst] = value;
              isatErrSyst++;
          break;
          default: cout << "unknown variable" << endl; continue;
          break;
        } //end switch
      } //end if isInsideScan
    } // end if fillInfo
    else continue;
  } //end while loop

  //Create tree file
  TFile *ghost_file = new TFile(ghost_file_name, "RECREATE");
  TTree *ghost_tree = new TTree("GhostCharge","GhostCharge");
  TTree *satellite_tree = new TTree("SatelliteCharge","SatelliteCharge");

  UInt_t Q_ghost_time         = 0;
  UInt_t Q_sat_time           = 0;
  UInt_t sigma_ghost_time     = 0;
  UInt_t sigma_sat_time       = 0;
  Double_t Q_ghost_value      = 0;
  Double_t Q_sat_value        = 0;
  Double_t sigma_ghost_value  = 0;
  Double_t sigma_sat_value    = 0;
  Double_t sigma_ghost_value_stat = 0;
  Double_t sigma_ghost_value_syst = 0;
  Double_t sigma_sat_value_stat   = 0;
  Double_t sigma_sat_value_syst   = 0;

  ghost_tree->Branch("Q_ghost_time",&Q_ghost_time,"Q_ghost_time/i");
  ghost_tree->Branch("Q_ghost_value",&Q_ghost_value,"Q_ghost_value/D");
  ghost_tree->Branch("sigma_ghost_time",&sigma_ghost_time,"sigma_ghost_time/i");
  ghost_tree->Branch("sigma_ghost_value",&sigma_ghost_value,"sigma_ghost_value/D");
  ghost_tree->Branch("sigma_ghost_value_stat",&sigma_ghost_value_stat,"sigma_ghost_value_stat/D");
  ghost_tree->Branch("sigma_ghost_value_syst",&sigma_ghost_value_syst,"sigma_ghost_value_syst/D");

  satellite_tree->Branch("Q_sat_time",&Q_sat_time,"Q_sat_time/i");
  satellite_tree->Branch("Q_sat_value",&Q_sat_value,"Q_sat_value/D");
  satellite_tree->Branch("sigma_sat_time",&sigma_sat_time,"sigma_sat_time/i");
  satellite_tree->Branch("sigma_sat_value",&sigma_sat_value,"sigma_sat_value/D");
  satellite_tree->Branch("sigma_sat_value_stat",&sigma_sat_value_stat,"sigma_sat_value_stat/D");
  satellite_tree->Branch("sigma_sat_value_syst",&sigma_sat_value_syst,"sigma_sat_value_syst/D");


  for (Int_t i=0; i<ighost; i++){
    Q_ghost_time = Q_ghost_time_array[i];
    Q_ghost_value = Q_ghost_value_array[i];
    sigma_ghost_time = sigma_ghost_time_array[i];
    sigma_ghost_value = sigma_ghost_value_array[i];
    Q_sat_time = Q_sat_time_array[i];
    Q_sat_value = Q_sat_value_array[i];
    sigma_sat_time = sigma_sat_time_array[i];
    sigma_sat_value = sigma_sat_value_array[i];
    sigma_ghost_value_stat = sigma_ghost_stat_value_array[i];
    sigma_ghost_value_syst = sigma_ghost_syst_value_array[i];
    sigma_sat_value_stat = sigma_sat_stat_value_array[i];
    sigma_sat_value_syst = sigma_sat_syst_value_array[i];
    ghost_tree->Fill();
    satellite_tree->Fill();
  }


  data_file.close();
  //save
  ghost_tree->Write();
  satellite_tree->Write();
  ghost_file->Close();

}

void Create_ghostSatellite_tree(Int_t Fill) {

  // initialize
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // find indices for start and end of scans
  cout << "Determining the start and end times of each scan, may take some time" << endl;
  Find_start_and_end_of_scans();
  SetScanTimes();

  //beam B1_A
  Get_ghostSatellite_tree(Fill,0);
  //beam B1_B
  Get_ghostSatellite_tree(Fill,1);
  //beam B2_A
  Get_ghostSatellite_tree(Fill,2);

}
