/// \file runAnalysisBC.C
/// \ingroup EMCALOfflineMacros
/// \brief Main function to run the bad channel analysis
///
///
/// See https://twiki.cern.ch/twiki/bin/view/ALICE/BadChannelAnalysis for general documentation
///
/// ---------------------
///  Running the macro
/// ---------------------
/// use root -b to speed up (no canvas drawn)                                                  <br>
/// root [1] .L $ALICE_WORK_DIR/../ali-master/AliPhysics/PWGPP/EMCAL/BCMacros/runAnalysisBC.C++                                                    <br>
/// root [2] runAnalysisBC(0,"LHC15o","Train_771","INT7",244918,"","AllRuns.txt")  	           //for merging to one runblock   <br>
/// root [2] runAnalysisBC(-1,"LHC15o","Train_771","INT7",244918,"244918_INT7Filtered.root","")  //for single files              <br>
///
/// \author Eliane Epple <eliane.epple@yale.edu>, Yale University
///
/// \date June 29, 2017


// --- ROOT system ---
#include <TStopwatch.h>
#include <TROOT.h>

// --- ANALYSIS system ---
#include "BadChannelAna.h" //include when compile

///definition of methods
///________________________________________________________________________
void runAnalysisBC(Int_t nversion = -1, TString period = "LHC15n", TString train = "Train_603", TString trigger= "INT7", Int_t runNum= 245683, TString externalFile= "",TString listName="runList.txt",TString workDir=".")
{
  gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output
  gROOT->SetBatch(1); //..Prevent ROOT from stealing focus when plotting

  std::vector<Int_t> badcellsManual;

  //..15o block 1
  //badcellsManual.insert(badcellsManual.end(),{14655,14622,14640,14728,14726,14593,14599,14600,14645,14646,14759,14776});

  //..If nversion=-1
  //..this is detected as a run-by-run analysis
  if(nversion == -1)nversion=runNum; //..If you do the analysis run by run - this might be helpful

  TStopwatch* watch = new TStopwatch();
  watch->Start();

  BadChannelAna* Analysis;
  Analysis=new BadChannelAna(period,train,trigger,runNum,nversion,workDir,listName);

  //. . . . . . . . . . . . . . . . . . . . . . . .
  //. . Settings
  //. . . . . . . . . . . . . . . . . . . . . . . .
  Analysis->SetExternalMergedFile(externalFile);
  //Analysis->AddManualMasking(badcellsManual);
  //Analysis->AddMaskSM(16);                    //..switch off entire SMs
  //Analysis->AddMaskSM(9);                   //..switch off entire SMs
  //Analysis->AddMaskSM(11);                  //..switch off entire SMs
  //Analysis->SetLowerBound(1);               //..If the Emin of the energy range (Emin-Emax) is higher than X GeV then dont apply a lower cut on the distribution
  Analysis->SetLowerBound(0.4);               //..If the Emin of the energy range (Emin-Emax) is higher than X GeV then dont apply a lower cut on the distribution
  //Analysis->SetStartEndCell(0,12288);       //..only EMCal
  //Analysis->SetStartEndCell(12288,17664);   //..only DCal
  //Analysis->SetQAChecks(1);                 //..1= Perform QA checks - takes a long time! Saves all good cells for cross check to pdf
  //Analysis->SetPrintOutput(1);              //..1= prints more information about excluded cells
  //Analysis->SetTrackCellRecord(1);          //..1= prints non-zero flag elements thoughout the routine
  //Analysis->SetExternalBadMap("Version286350_mSM11/LHC18d_INT7_Histograms_V286350.root");

  //. . . . . . . . . . . . . . . . . . . . . . . .
  //. . Add different period analyses
  //. . . . . . . . . . . . . . . . . . . . . . . .
  Double_t sigmaNHits=5.5;
  Double_t sigmaE_hit=4.5;
  //..a little stricter
  sigmaNHits=5.0;
  sigmaE_hit=4.5;
  //..the range of sigmas should be selected such
  //..that one does not cut into the natural fluctuation over the modules
  Analysis->AddPeriodAnalysis(5, sigmaNHits,0.1,0.3);  // hits in cell in range Emin Emax
  Analysis->AddPeriodAnalysis(1, sigmaE_hit,0.1,0.3);  // energy/hit in range Emin Emax
  Analysis->AddPeriodAnalysis(5, sigmaNHits,0.2,0.5);  // hits in cell range Emin Emax
  Analysis->AddPeriodAnalysis(1, sigmaE_hit,0.2,0.5);  // energy/hit in range Emin Emax
  Analysis->AddPeriodAnalysis(5, sigmaNHits,0.5,1.0);  // hits in cell range Emin Emax
  Analysis->AddPeriodAnalysis(1, sigmaE_hit,0.5,1.0);  // energy/hit in range Emin Emax
  Analysis->AddPeriodAnalysis(5, sigmaNHits,1.0,4.0);  // hits in cell range Emin Emax
  Analysis->AddPeriodAnalysis(1, sigmaE_hit,1.0,4.0);  // mean energy in range Emin Emax
  Analysis->AddPeriodAnalysis(5, sigmaNHits,1.0,10.0); // hits in cell in range Emin Emax
  Analysis->AddPeriodAnalysis(1, sigmaE_hit,1.0,10.0); // energy/hit in range Emin Emax
  Analysis->AddPeriodAnalysis(5, sigmaNHits-1,0.11,0.29);  // hits in cell in range Emin Emax


  //..special test for extra high energy fluctuations
  Analysis->AddPeriodAnalysis(1, 4.0,3.0,40.0); // energy/hit in cell in range Emin Emax
  //	Analysis->AddPeriodAnalysis(1, 80.0,3.0,40.0); // IN CASE something fails - energy/hit in cell in range Emin Emax
  //	Analysis->AddPeriodAnalysis(5, 4.0,3.0,40.0); // hits in cell range Emin Emax
  //Analysis->AddPeriodAnalysis(1, 4.5,3.0,5.0);  // mean energy in range Emin Emax - cliff
  //Analysis->AddPeriodAnalysis(5, 5.5,3.0,5.0);  // hits in cell range Emin Emax   - cliff

  //..Add a cut on the time distribution
  //	Analysis->AddPeriodAnalysis(3, 6,-20,+20);  // hits with time -20-20ns/all times
  Analysis->AddPeriodAnalysis(3, 4,550,750); // hits with time 550-750ns/all times
  //	Analysis->AddPeriodAnalysis(3, 16,550,750); // IN CASE there are some problems with the fit hits with time 550-750ns/all times

  //*test time stuff*/	Analysis->AddPeriodAnalysis(3, 6,-20,+20);// energy/hit in range Emin Emax

  //. . . . . . . . . . . . . . . . . . . . . . . .
  //. . Start the bad channel analysis
  //. . . . . . . . . . . . . . . . . . . . . . . .
  Bool_t mergeOnly=0;//.. =1 do only merge and filter
  Analysis->Run(mergeOnly);

  watch->Stop();
  cout<<"Finished BC analysis "<<watch->RealTime()/60<<" min"<<endl;
  cout<<"Canvases will pop up now - please wait"<<endl;

  //if(watch)    delete watch;
  //if(Analysis) delete Analysis; //..Delete Causes problems
}
