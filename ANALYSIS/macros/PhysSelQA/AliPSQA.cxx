/////////////////////////////////////
// Created by: Kevin McDermott     //
// email: kmcderm3@nd.edu          //
// CERN Summer Student 2012        //
// University of Notre Dame du Lac //
//                                 // 
// Revision: 1.0                   //
// Created on: August 6, 2012      //
/////////////////////////////////////

#include "AliPSQA.h"
ClassImp(AliPSQA)

AliPSQA::AliPSQA():
TObject(),
  fRunNumbers(0)
{
  // Default constructor settings for private data members
  fQAcycle = 0; 
  fLocalMode = kTRUE;
  fGridPath = ""; 
  fPSInput = ""; 
  fLocalPath = "";
  fROOTInput = "";
  fOutRootName = "";
  fCacheDirectory = "";
  fRootOutDirectory = "";
  fTxtOutDirectory = "";
  fListBadFiles = new TList;
  fNGoodRuns = 0; //Files that are not empty
  fNTrigClasses = 0;
  fNTrigChannels = 0;
  fTxtFileName = "";
  fSaveCache = kFALSE;
  fDrawAllTGraphs = kFALSE;
  fSaveQAValToTxtFile = kFALSE;
  fSaveQAToROOTFile = kFALSE;

  // Initialize the macro
  InitializeRuns("");
  InitializeTriggers("","");
  InitializePlots("");

  InitializeRowTrigSet();
  InitializeFilled(); // Set once per object
  InitializeValidRun();
}

//________________________________________________________________________________________________
AliPSQA::~AliPSQA(){
  delete[] fRunNumbers;
  delete[] fTrigClasses;
  delete[] fTrigChannels;
	for(Int_t i=0;i<15;i++){
		for(Int_t j=0;j<30;j++){
			for(Int_t k=0;k<2;k++){
				for(Int_t l=0;l<fgkNTrigLogic;l++){
					delete[]  fGraphs[i][j][k][l];
				}
			}
		}
	}
 // delete[] fGraphs;
  delete[] fPlots;
  delete[] fQAColumnsNumeratorP;
  delete[] fQAColumnsDenominatorP;
  delete[] fQAColumnsNumeratorM;
  delete[] fQAColumnsDenominatorM;
  delete fListBadFiles;
}

//________________________________________________________________________________________________
void AliPSQA::InitializeRuns(const Char_t * listOfRuns)
{
  // without vectors, this process is quite contrived... Why doesn't ROOT have STL vectors??????

  Int_t run_num; //run number from list
  ifstream input_runs; //text file object
  input_runs.open(listOfRuns, ios::in );  //open the text file

  Int_t NRuns = 0; //number of runs to be processed 
  while(input_runs >> run_num) { // count number of runs
    NRuns++;
  }

  fNRuns = NRuns; // Set private data member to total runs to be processed
  input_runs.close(); // reset the file

  fRunNumbers = new Int_t [fNRuns]; // initialize private data member with number of runs
  input_runs.open(listOfRuns, ios::in ); // open the ifstream again
  for (Int_t irun = 0; irun < fNRuns; irun++){ 
    input_runs >> fRunNumbers[irun]; // fill this array with run numbers
  }

  input_runs.close();

  if (fNRuns == 0) {
    Printf("Could not initialize run List !");
    return;
  }
  
}

//________________________________________________________________________________________________
void AliPSQA::InitializeTriggers(const Char_t * listOfTriggerNames, const Char_t * listOfTrigChannels){

  TString trignames; //trigger names from list
  ifstream input_tn; //text file object
  input_tn.open(listOfTriggerNames, ios::in );  //open the text file
  Int_t Ntn = 0; //number of trigger names to be used
  
  while ( input_tn >> trignames ) { // count number of trigger names
    Ntn++;
  }
  
  fNTrigClasses = Ntn; // Set private data member to total trigger names to be used
  input_tn.close(); // reset the file

  fTrigClasses = new TString [fNTrigClasses]; // initialize private data member with number of trig names
  input_tn.open(listOfTriggerNames, ios::in );
  for (Int_t itrigclass = 0; itrigclass < fNTrigClasses; itrigclass++) { // read in trig names, fill private data members
    input_tn >> fTrigClasses[itrigclass];
  }

  input_tn.close();  

  TString trigchan; //trigger names from list
  ifstream input_tc; //text file object
  input_tc.open(listOfTrigChannels, ios::in );  //open the text file
  Int_t Ntc = 0; //number of trigger names to be used
  
  while ( input_tc >> trigchan ) { // count number of trigger names
    Ntc++;
  }
  
  fNTrigChannels = Ntc; // Set private data member to total trigger names to be used
  input_tc.close(); // reset the file

  fTrigChannels = new TString [fNTrigChannels]; // initialize private data member with number of trig names
  input_tc.open(listOfTrigChannels, ios::in );
  for (Int_t itrigchannel = 0; itrigchannel < fNTrigChannels; itrigchannel++) { // read in trig names, fill private data members
    input_tc >> fTrigChannels[itrigchannel]; 
  }

  input_tc.close();  
}

//________________________________________________________________________________________________
void AliPSQA::InitializePlots(const Char_t * listOfPlots){   // Allocate memory for all fGraphs

  // First Dimension is the name of the TGraphErrors from list file

  TString plotnames; //trigger names from list
  ifstream input_pn; //text file object
  input_pn.open(listOfPlots, ios::in );  //open the text file
  Int_t Npn = 0; //number of trigger names to be used
  
  while ( input_pn >> plotnames ) { // count number of trigger names
    Npn++;
  }
  
  fNPlots = Npn; // Set private data member to total trigger names to be used
  input_pn.close(); // reset the file

  fPlots = new TString [fNPlots]; // initialize private data member with number of trig names
  input_pn.open(listOfPlots, ios::in );
	
	cout << "Number of Plots : " << fNPlots << endl;
	
  for (Int_t iplot = 0; iplot < fNPlots; iplot++) { // read in trig names, fill private data members
    input_pn >> fPlots[iplot];
  }

  input_pn.close();

  // Set the Names of the Plots
  fQAColumnsNumeratorP = new TString [fNPlots];
  fQAColumnsDenominatorP = new TString [fNPlots];
  for (Int_t iplot = 0; iplot < fNPlots; iplot++){
    SetQAColumnLabelForPlots(iplot); // Set fQAColumn for each pair of variables to be analyzed, keeping the "_" for plot names
  }
  
  // Set the Names of the columns to be used, see GetQAVals and PlotQA
  fQAColumnsNumeratorM = new TString [fNPlots];
  fQAColumnsDenominatorM = new TString [fNPlots];
  for (Int_t iplot = 0; iplot < fNPlots; iplot++){
    SetQAColumnLabelForMatching(iplot); // Set fQAColumn for each pair of variables to be analyzed, removing all "_", needed for string matching of fHistStats
  }

  // Allocate Memory for all trending Plots

  for (Int_t i = 0; i < fNPlots; i++){
    for (Int_t j = 0; j < fNTrigClasses; j++) {
      for (Int_t k = 0; k < fNTrigChannels; k++) {
	for (Int_t l = 0; l < fgkNTrigLogic; l++){ // Settings for Graphs 
	  fGraphs[i][j][k][l] = new TGraphErrors(fNRuns);
	  fGraphs[i][j][k][l]->SetName(Form("%s_over_%s_for_%s_%s_Trigger_%02i",fQAColumnsNumeratorP[i].Data(),fQAColumnsDenominatorP[i].Data(),fTrigClasses[j].Data(),fTrigChannels[k].Data(), l));
	  fGraphs[i][j][k][l]->SetTitle(Form("%s over %s for %s [%s] Trigger %02i",fQAColumnsNumeratorM[i].Data(),fQAColumnsDenominatorM[i].Data(),fTrigClasses[j].Data(),fTrigChannels[k].Data(), l));
	  fGraphs[i][j][k][l]->GetXaxis()->SetTitle("Run Number");
	  fGraphs[i][j][k][l]->GetYaxis()->SetTitle(Form("%s/%s",fQAColumnsNumeratorM[i].Data(),fQAColumnsDenominatorM[i].Data()));
	  fGraphs[i][j][k][l]->SetMarkerStyle(6);
		
		cout << "Names : " << fGraphs[i][j][k][l]->GetName() << endl;
	}	
      }
    }
  }
  
}

//________________________________________________________________________________________________
void AliPSQA::SetQAColumnLabelForPlots(Int_t iplot){

  Ssiz_t start  = fPlots[iplot].First('h'); // h in fGraph
  Ssiz_t finish = fPlots[iplot].Last('O'); // O in Over, first O may be from FO
  Ssiz_t length = fPlots[iplot].Length();

  TString label[2];
  label[0] = TString(fPlots[iplot]( start + 1, finish - start - 1)); // Shift one over from the 'h'
  label[1] = TString(fPlots[iplot]( finish + 4, length - finish - 4 ));  //Shift five over from the 'O' 

  fQAColumnsNumeratorP[iplot]   = label[0];
  fQAColumnsDenominatorP[iplot] = label[1];
}

//________________________________________________________________________________________________
void AliPSQA::SetQAColumnLabelForMatching(Int_t iplot){
  fQAColumnsNumeratorM[iplot] = fQAColumnsNumeratorP[iplot].Copy();
  fQAColumnsDenominatorM[iplot] =  fQAColumnsDenominatorP[iplot].Copy();
  
  if(fQAColumnsNumeratorM[iplot].Contains("_")){ // Plot name contains _ in name
    fQAColumnsNumeratorM[iplot].ReplaceAll("_",1," ",1);
  }
  if(fQAColumnsDenominatorM[iplot].Contains("_")){ // Plot name contains _ in name
    fQAColumnsDenominatorM[iplot].ReplaceAll("_",1," ",1);
  }
}

//________________________________________________________________________________________________
void AliPSQA::InitializeRowTrigSet(){ // Reset the check for each unique trig combo to false for each run
  for (Int_t i = 0; i < fNTrigClasses; i++) {
    for (Int_t j = 0; j < fNTrigChannels; j++) {
      for (Int_t k = 0; k < fgkNTrigLogic; k++){
	fRowTrigSet[i][j][k] = kFALSE;
      }
    }
  }
}

//________________________________________________________________________________________________
void AliPSQA::InitializeFilled(){ // Initialize check to see if a trending plot should be written to a ROOT file
  for (Int_t i = 0; i < fNPlots; i++){
    for (Int_t j = 0; j < fNTrigClasses; j++) {
      for (Int_t k = 0; k < fNTrigChannels; k++) {
	for (Int_t l = 0; l < fgkNTrigLogic; l++){
	  fFilled[i][j][k][l] = kFALSE;
	}
      }
    }
  }
}

//________________________________________________________________________________________________
void AliPSQA::InitializeValidRun(){ // Reset the check for a valid run for each run
  for (Int_t i = 0; i < fNTrigClasses; i++) {
    for (Int_t j = 0; j < fNTrigChannels; j++) {
      for (Int_t k = 0; k < fgkNTrigLogic; k++){
	fValidRun[i][j][k] = kFALSE;
      }
    }
  }
}
 
//________________________________________________________________________________________________
void AliPSQA::CacheFiles(){
  
  for(Int_t irun = 0; irun < fNRuns; irun++){

    fCacheDirectory        = Form("%s/%09i/%s",fCachePath.Data(),fRunNumbers[irun],fPSInput.Data());
    TString cacheFile      = Form("%s/%s",fCacheDirectory.Data(),fROOTInput.Data());
  
    // Check first to see if the root file exists, have to declare a dummy FileStat_t to use gSystem->GetPathInfo
    // GetPathInfo returns 0 if path/file exists, 1 if it does not

    FileStat_t dummy_filestat;

    if(gSystem->GetPathInfo(cacheFile.Data(), dummy_filestat) == 1){ // file does not exist
      // make this directory if cache directory does not exist, to prevent overwriting cached files already stored
      if (gSystem->GetPathInfo(fCacheDirectory.Data(), dummy_filestat) == 1) { // cache directory does not exist
	MakeDir(fCacheDirectory);
      }
    
      TString gridFile = ""; // full path for grid file
      // Get the grid path of the corresponding file to be cached
      if(fPSInput.Contains("AOD")){ //Include QA Cycle for AODS
	gridFile.Form("alien://%s/%09d/%s%3.3d/%s",fGridPath.Data(),fRunNumbers[irun],fPSInput.Data(),fQAcycle,fROOTInput.Data());
      }
      else{ // Do not Include QA Cycle for AODS
	gridFile.Form("alien://%s/%09d/%s/%s",fGridPath.Data(),fRunNumbers[irun],fPSInput.Data(),fROOTInput.Data());
      }

      // copy to the directory 
      TString cpCache = "alien_cp ";
      cpCache+=gridFile;
      cpCache+=" "; // copy to the current directory to directory specified
      cpCache+=fCacheDirectory; 
      gSystem->Exec(cpCache); // copy to the directory 
    }
  }
}

//________________________________________________________________________________________________

void AliPSQA::MakeDir(TString dir){ // have to pass a private data member, kinda ugly, but otherwise use sentinel values or copy this function three times...
  TString mkDir = "mkdir -p ";
  mkDir += dir.Data();
  gSystem->Exec(mkDir); // mkdir for the cache/output
}

//________________________________________________________________________________________________
Bool_t AliPSQA::ComputeQAPerRun(){
  for(Int_t irun = 0; irun < fNRuns; irun++){
    
    // Reset RowTrigSet to kFALSE for every run
    
    InitializeRowTrigSet();
    InitializeValidRun();
    
    TFile *rootfile = 0; // root file to be used
    TString filename = ""; // full path directory of root file to be used

    filename = GetRootFileName(irun, filename);
    Printf("\nBegin reading: %s", filename.Data());
    rootfile=TFile::Open(filename.Data());

	  cout << "filename processed = " << filename << endl;
	  
    // If the file is not available, continue
    if(!rootfile){
      Printf("File %d is not available.\n",fRunNumbers[irun]);
      fListBadFiles->Add(new TObjString(Form("No root file found in path: %s",filename.Data())));
      continue;
    }

    // Get the relevant histogram, fHistStatistics
    TH2D * hStats = (TH2D*) rootfile->Get("fHistStatistics");
    if(!hStats){
      fListBadFiles->Add(new TObjString(Form("fHistStats not found in path: %s",filename.Data())));
      continue;
    }

    // Create output to be stored in a txt file, to be used for webpage
    std::ofstream outfile; // to be used to make a txtfile
    OpenTxtFile(irun, outfile);

    Bool_t oneFilledGraphPerRun = kTRUE; // Check to see if run has any entries

    // Loop over all possible values of & and * to generate data for each run
    
    for(Int_t itrigclass = 0; itrigclass < fNTrigClasses; itrigclass++){
      for(Int_t itrigchannel = 0; itrigchannel < fNTrigChannels; itrigchannel++){
        for(Int_t itriglogic = 0; itriglogic < fgkNTrigLogic; itriglogic++){
            fValidRun[itrigclass][itrigchannel][itriglogic] = ProcessQA(itrigclass, itrigchannel, itriglogic, irun, hStats, outfile);
            // Set boolean once to check if run contains any data
            if( (oneFilledGraphPerRun == kFALSE) && (fValidRun[itrigclass][itrigchannel][itriglogic] == kTRUE) ) {
                oneFilledGraphPerRun = kTRUE;
            }
        }
      }
    }
    
    if (oneFilledGraphPerRun == kTRUE){
      cout << "Wrote " << Form("%s/%s_%09i.txt",fTxtOutDirectory.Data(),fTxtFileName.Data(),fRunNumbers[irun]) << endl;
      fNGoodRuns++;
    }
    else {
      fListBadFiles->Add(new TObjString(Form("No events found in path %s",filename.Data())));
    }
    // Close the txt file, need a new one every run
    outfile.close();
  }// End loop over runs

  if (fSaveQAToROOTFile){
    SaveQAPlots(); // Will write all good graphs (after being cleaned up) to a root file.  Will make output directory if it does not exist
  }
      
  cout << "Files Statistics" << endl;
  cout << " Total Runs [" << fNRuns << "]" << endl;
  cout << " Good  Runs [" << fNGoodRuns << "]" <<  endl;
  
  if (fNGoodRuns != 0){
    return kTRUE;
  }
  else{
    return kFALSE;
  }
}// End ComputeQAPerRun()

//________________________________________________________________________________________________
TString AliPSQA::GetRootFileName(Int_t irun, TString filename){
  // check if QA file is available and open it
  // try different QA train outputs
  if(!fLocalMode){ // Run on Alien
    if(fPSInput.Contains("AOD")){ //Include QA Cycle for AODS
      filename.Form("alien://%s/%09d/%s%3.3d/%s",fGridPath.Data(),fRunNumbers[irun],fPSInput.Data(),fQAcycle,fROOTInput.Data());
    }
    else{ // Do not Include QA Cycle for AODS
      filename.Form("alien://%s/%09d/%s/%s",fGridPath.Data(),fRunNumbers[irun],fPSInput.Data(),fROOTInput.Data());
    }
  }
  else{ // Run locally, either cached files, or some other directory
    filename.Form("%s/%09d/%s/%s",fLocalPath.Data(),fRunNumbers[irun],fPSInput.Data(),fROOTInput.Data());
  }
  return filename;
}

//________________________________________________________________________________________________
void AliPSQA::OpenTxtFile(Int_t irun, std::ofstream & outfile){ // Create a text file per run of qavals
  TString txtfile = "";
  txtfile.Form("%s/%s_%09i.txt",fTxtOutDirectory.Data(),fTxtFileName.Data(),fRunNumbers[irun]);
  outfile.open(txtfile.Data(), ios::out);
}

//________________________________________________________________________________________________

Bool_t AliPSQA::ProcessQA(Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, Int_t irun, TH2D* hStats, std::ofstream & outfile){
  
  Int_t ibiny = SearchYLabels(itrigclass, itrigchannel, itriglogic, hStats);
  if (ibiny < 0){ // ifibiny = -1, (no labels match) returns false for ProcessQA, skip this combo of * and & 
    return kFALSE;
  }
  
  Bool_t isNotTotallyEmpty = kFALSE; // Check to see if qa is totally empty, assume totally empty
  Double_t qavals[2]; // qavals[0] = All event hits, used to check if run is empty, qavals[1] = actual qa value being processed for a given run, and trig combo

  for (Int_t iplot = 0; iplot < fNPlots; iplot++){ // Loop over all plots to get qa vals and plot them
   
    Int_t ibinx[2] = {-1,-1}; // if ibinx = -1 (if no labels match), then move to next column
    SearchXLabels(iplot, ibinx, hStats); 
    if ( (ibinx[0] < 0) || (ibinx[1] < 0) ){ 
      continue;
    }

    // Get QA Vals if all labels are valid

    GetQAVals(hStats, ibinx, ibiny, qavals);
    if ( (qavals[0] < 1e-6) || (qavals[1] < 1e-6) ){ // if no events in PS, then continue
      continue;
    }
    else{ // otherwise, fill plots and save to txt file
      if (isNotTotallyEmpty == kFALSE){
	isNotTotallyEmpty = kTRUE;
      }
      
      // Part of loop to save the output data

      if (fSaveQAToROOTFile){ 
	PlotQA(iplot, itrigclass, itrigchannel, itriglogic, irun, qavals);
      }
      
      if (fSaveQAValToTxtFile){ // if true, save a text file of all qavals per run
	// First Check to make sure path exists to save data
	FileStat_t dummy_filestat;
	if (gSystem->GetPathInfo(fTxtOutDirectory.Data(), dummy_filestat) == 1) { // output directory does not exist
	  MakeDir(fTxtOutDirectory); // therefore make it
	}
	SaveQAValToTxtFile(iplot, itrigclass, itrigchannel, itriglogic, qavals, outfile);
      }
    }
  }
  return isNotTotallyEmpty;
}

//________________________________________________________________________________________________
Int_t AliPSQA::SearchYLabels(Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, TH2D* hStats){ // returns the ybin to get data
  // Search for each unique combination of * and & along y-axis of fHistStats
  Int_t trigMask[3]; // trigMask[0] = trig class, reduced &val; trigMask[1] = trig channel/partition, from &val if 30th bit is on; trigMask[2] = *val, trig logic

  Int_t nbiny = hStats->GetNbinsY();
  
  for(Int_t ibiny = 1; ibiny <= nbiny; ibiny++){ // Y-Axis for fHistStats starts at bin 1
    TString label = hStats->GetYaxis()->GetBinLabel(ibiny);
    // Check the trig logic and trig mask

    GetTrigMask(label, trigMask);
    //parse the label to return whether fast or slow --> return trigMaskChannel

    if( (trigMask[0]==itrigclass) && (trigMask[1]==itrigchannel) && (trigMask[2]==itriglogic) ) { // set for unique trig class, channel, and trig number
      
      if (fRowTrigSet[itrigclass][itrigchannel][itriglogic] == kFALSE){
	fRowTrigSet[itrigclass][itrigchannel][itriglogic] = kTRUE;
	return ibiny;  // set ybin to get data; 
      }
    }
  }	
  return -1; // returns -1 if there is no matching y label
}

//________________________________________________________________________________________________
void AliPSQA::GetTrigMask(TString label, Int_t * trigMask){
  // Ssiz_t -> int returned for length of string, at this index

  Ssiz_t start  = label.First('&');
  Ssiz_t finish = label.First('*');
  Ssiz_t length = label.Length();
  
  TString maskLabel( label( start+1, finish - start - 1 ) ); // returns the substring between the & and *
  Int_t   maskInt = maskLabel.Atoi(); // Convert string to int (the int is in decimal form)
  
  TString logicLabel( label( finish+1, length - finish) );
  Int_t   triglogic = logicLabel.Atoi(); // trig logic number, from *val

  // convert int to binary check for what channel it is and see what class/mask

  ////////////////////////////////

  // Converts int to string with the binary representation --> code snippet taken from online

  TString binary = DecInttToBinTString(maskInt);

  ///////////////////////////////

  Ssiz_t secondbit = binary.Last('1');
  Ssiz_t binlength = binary.Length();
  
  Int_t trigclass; // trig class name
  Int_t trigchannel; // fast or regular?

  if (binlength != 31){ //Regular channel, 30th bit for kFast is off
    trigclass   = binlength - 1; // equivalent statements = secondbit - 1 = firstbit - 1 = length - firstbit
    trigchannel = 0;
  }
  else{ //Fast channel, if more than just regular vs fast is specified, must modify parsing and provide more robust else-ifs
    trigclass   = binlength - secondbit -1;
    trigchannel = 1;
  }
  trigMask[0] = trigclass;
  trigMask[1] = trigchannel;
  trigMask[2] = triglogic;
}

//________________________________________________________________________________________________
TString AliPSQA::DecInttToBinTString(Int_t maskInt){
  if ( maskInt == 0 ) return "0";
  if ( maskInt == 1 ) return "1";
  
  if ( maskInt % 2 == 0 ){
    return (DecInttToBinTString(maskInt / 2) + "0");
  }
  else{
    return (DecInttToBinTString(maskInt / 2) + "1");
  }
}

//________________________________________________________________________________________________
void AliPSQA::SearchXLabels(Int_t iplot, Int_t * ibinx, TH2D * hStats){

  // Search x-axis for label matching qa column val to be used
  Int_t nbinx = hStats->GetNbinsX(); // Total nbins in x-axis
  
  for (Int_t jbinx = 1; jbinx <= nbinx; jbinx++){ // look for the column value for numerator
    TString label = hStats->GetXaxis()->GetBinLabel(jbinx);
    if (label.CompareTo(fQAColumnsNumeratorM[iplot].Data(), TString::kExact) == 0){ // CompareTo returns 0 if the two TStrings are the same
      ibinx[0] = jbinx;
      break;
    }
  }

  for (Int_t jbinx = 1; jbinx <= nbinx; jbinx++){ // look for the column value for denominator
    TString label = hStats->GetXaxis()->GetBinLabel(jbinx);
    if (label.CompareTo(fQAColumnsDenominatorM[iplot].Data(), TString::kExact) == 0){ // CompareTo returns 0 if the two TStrings are the same
      ibinx[1] = jbinx;
      break;
    }
  }
}

//________________________________________________________________________________________________
void AliPSQA::GetQAVals(TH2D * hStats, Int_t * ibinx, Int_t ibiny, Double_t * qavals){ 
  qavals[0] = hStats->GetBinContent(ibinx[0],ibiny); // QAval of the first column called, numerator 
  qavals[1] = hStats->GetBinContent(ibinx[1],ibiny); // QAval demoninator for ratio
}

//________________________________________________________________________________________________
void AliPSQA::PlotQA(Int_t iplot, Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic, Int_t irun, Double_t * qavals){
  fGraphs[iplot][itrigclass][itrigchannel][itriglogic]->SetPoint(irun, fRunNumbers[irun], qavals[0]/qavals[1]); // QA Val for the plot
  fGraphs[iplot][itrigclass][itrigchannel][itriglogic]->SetPointError(irun, 0., GetError(qavals, iplot)); // QA val error for correlated subsets
  
  if (fFilled[iplot][itrigclass][itrigchannel][itriglogic] == kFALSE){
    fFilled[iplot][itrigclass][itrigchannel][itriglogic] = kTRUE; // need this boolean for later to only save filled graphs
  }
}

//________________________________________________________________________________________________
Double_t AliPSQA::GetError(Double_t * qavals, Int_t iplot){
  Double_t error;
  if ( fQAColumnsDenominatorM[iplot].CompareTo("Trigger class", TString::kExact) == 0) { // Set error if one QA val is subset of other
    error = UncNRatioCorrelated(qavals);
  }
  else{ // Set error for correlated sets, but not subsets of each other
    error = ErrorRatioNotSubsets(qavals);
  }
  return error;
}

//________________________________________________________________________________________________
Double_t AliPSQA::UncNRatioCorrelated(Double_t * qavals)
{
  Double_t eff = qavals[0] / qavals[1];
  Double_t err = TMath::Sqrt( eff*(1.0-eff) ) / TMath::Sqrt( qavals[1] );
  return err;
}

//________________________________________________________________________________________________
Double_t AliPSQA::ErrorRatioNotSubsets(Double_t * qavals)
{
  Double_t qa0err = TMath::Power((TMath::Sqrt(qavals[0]) / qavals[0]), 2);
  Double_t qa1err = TMath::Power((TMath::Sqrt(qavals[1]) / qavals[1]), 2);
  Double_t err = TMath::Sqrt(qa0err + qa1err) * (qavals[0] / qavals[1]);
  return err;
}

//________________________________________________________________________________________________
void AliPSQA::SaveQAValToTxtFile(Int_t iplot, Int_t itrigclass, Int_t itrigchannel, Int_t itriglogic,  Double_t * qavals, std::ofstream & outfile){
  outfile << "PS: \"" << fTrigClasses[itrigclass].Data() << " [" << fTrigChannels[itrigchannel].Data() << "] Trigger " << itriglogic << "\":" << fQAColumnsNumeratorP[iplot].Data() << "/" << fQAColumnsDenominatorP[iplot].Data() << ", " << qavals[0]/qavals[1] << ", " << GetError(qavals, iplot) << endl; // store qavals in txt file
}

//________________________________________________________________________________________________
void AliPSQA::SaveQAPlots(){ // Return a boolean to see if QA check produced any plots

  // First check to see if output directory exists, otherwise create it
  FileStat_t dummy_filestat;
  if (gSystem->GetPathInfo(fRootOutDirectory.Data(), dummy_filestat) == 1) { // output directory does not exist
    MakeDir(fRootOutDirectory);
  }

  // Save the output in a root file only if fFilled for that plot == kTRUE
  
  TFile * rootfile = new TFile(Form("%s/%s",fRootOutDirectory.Data(),fOutRootName.Data()),"RECREATE");
  rootfile->cd();
 
  for (Int_t iplot = 0; iplot < fNPlots; iplot++){
    for (Int_t itrigclass = 0; itrigclass < fNTrigClasses; itrigclass++){
      for (Int_t itrigchannel = 0; itrigchannel < fNTrigChannels; itrigchannel++){
	for (Int_t itriglogic = 0; itriglogic < fgkNTrigLogic; itriglogic++){
		cout << "File name : " << fGraphs[iplot][itrigclass][itrigchannel][itriglogic]->GetName() << " filled? yes/no " << fFilled[iplot][itrigclass][itrigchannel][itriglogic] << endl; 
	  if (fFilled[iplot][itrigclass][itrigchannel][itriglogic] == kTRUE){
		  
	    CleanUpGraphs(fGraphs[iplot][itrigclass][itrigchannel][itriglogic]);
	    fGraphs[iplot][itrigclass][itrigchannel][itriglogic]->Write();
		  
		  //cout << "Written to file : " << fGraphs[iplot][itrigclass][itrigchannel][itriglogic]->GetName() << endl; 
	  }
	} 
      }
    }
  }
  cout << "Wrote " << Form("%s/%s",fRootOutDirectory.Data(),fOutRootName.Data()) << endl;
  rootfile->Close();
  delete rootfile;
}

//________________________________________________________________________________________________
void AliPSQA::CleanUpGraphs(TGraphErrors * plot){ // Clean up default values of graphs to be saved in ROOT file

  // Have to unfortunately change the looping variables, scary, I know, but otherwise the file will not be cleaned up correctly.  See documentation on ROOT page for RemovePoint and GetPoint.

  Double_t x;
  Double_t y;
  for (Int_t irun = 0; irun < plot->GetN(); irun++){
    plot->GetPoint(irun, x, y);
    if ( (TMath::Abs(x) < 1e-6) && (TMath::Abs(y) < 1e-6) ){
      plot->RemovePoint(irun);
      irun--;
    }
  }
}

//________________________________________________________________________________________________
void AliPSQA::DrawAllTGraphs(){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  // do not use scientific notation for run number
  TGaxis::SetMaxDigits(6);
  
  TCanvas * c1;
  
  for (Int_t i = 0; i < fNPlots; i++){
    for (Int_t j = 0; j < fNTrigClasses; j++){
      for (Int_t k = 0; k < fNTrigChannels; k++){
	for (Int_t l = 0; l < fgkNTrigLogic; l++){
	  if (fFilled[i][j][k][l] == kTRUE){
	    c1 = new TCanvas;
	    c1->cd();
	    fGraphs[i][j][k][l]->SetName(Form("%s_over_%s_for_%s_%s_Trigger_%02i",fQAColumnsNumeratorP[i].Data(),fQAColumnsDenominatorP[i].Data(),fTrigClasses[j].Data(),fTrigChannels[k].Data(), l));
	    fGraphs[i][j][k][l]->SetTitle(Form("%s over %s for %s [%s] Trigger %02i",fQAColumnsNumeratorM[i].Data(),fQAColumnsDenominatorM[i].Data(),fTrigClasses[j].Data(),fTrigChannels[k].Data(), l));
	    fGraphs[i][j][k][l]->GetXaxis()->SetTitle("Run Number");
	    fGraphs[i][j][k][l]->GetYaxis()->SetTitle(Form("%s/%s",fQAColumnsNumeratorM[i].Data(),fQAColumnsDenominatorM[i].Data()));
	     
	    fGraphs[i][j][k][l]->Draw("AP");
	  }
	} 
      }
    }
  }
}

