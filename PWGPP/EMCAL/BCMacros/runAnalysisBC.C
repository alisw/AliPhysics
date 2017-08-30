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
/// root [1] .L runAnalysisBC.C++                                                    <br>
/// root [2] Run_BadChannel(1,"LHC15o","Train_771","INT7",244918,"","GloballyGood.txt")           //for merging to one runblock   <br>
/// root [2] Run_BadChannel(-1,"LHC15o","Train_771","INT7",244918,"244918_INT7Filtered.root","")  //for single files              <br>
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
void Run_BadChannel(Int_t nversion = -1, TString period = "LHC15n", TString train = "Train_603", TString trigger= "AnyINT", Int_t runNum= 245683, TString externalFile= "",TString listName="runList.txt",TString workDir=".")
{
	gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;"); //..to supress a lot of standard output
    std::vector<Int_t> badcellsManual;

    //..LHC16l
/*    badcellsManual.insert(badcellsManual.end(),{7120,7121,7122,7123,7124,7125,7126,7127,7128,7129,7130,7131,7132,7133,7134,7135});
    badcellsManual.insert(badcellsManual.end(),{7168,7169,7170,7171,7172,7173,7174,7175,7176,7177,7178,7180,7181,7182,7183});
    badcellsManual.insert(badcellsManual.end(),{7537,7872,7873,7874,7875,7876,7877,7878,7879,7880,7881,7882,7883,7884,7885,7886,7887});
    badcellsManual.insert(badcellsManual.end(),{7920,7921,7922,7923,7924,7925,7926,7927,7928,7929,7930,7931,7932,7933,7934,7935});
    badcellsManual.insert(badcellsManual.end(),{7968,7969,7970,7971,7972,7973,7974,7975,7976,7977,7978,7979,7980,7981,7982,7983});
    badcellsManual.insert(badcellsManual.end(),{8016,8017,8018,8019,8020,8021,8022,8023,8024,8025,8026,8027,8028,8029,8030,8031});
    badcellsManual.insert(badcellsManual.end(),{8656,8657,8658,8659,8660,8661,8662,8663,8664,8665,8666,8667,8668,8669,8670,8671});
    badcellsManual.insert(badcellsManual.end(),{8704,8705,8706,8707,8708,8709,8710,8711,8712,8713,8714,8715,8716,8717,8718,8719});
    badcellsManual.insert(badcellsManual.end(),{8753,8754,8755,8756,8757,8758,8759,8760,8761,8762,8763,8764,8765,8766,8767});
    badcellsManual.insert(badcellsManual.end(),{8800,8801,8802,8803,8804,8805,8806,8807,8808,8809,8810,8811,8812,8813,8814,8815});
    badcellsManual.insert(badcellsManual.end(),{10134,10811,11630,11904,11905,11906,11907,11908,11909,11910,11911,11912,11913,11914,11915,11916,11917,11918,11919});
    badcellsManual.insert(badcellsManual.end(),{12032,12033,12034,12035,12036,12037,12038,12039,12040,12041,12042,12043,12044,12045,12046,12047});
    badcellsManual.insert(badcellsManual.end(),{12170,12172,12622,12870,12876,13973,14264,14265,14266,14267,14320,14400,14553,14621,14722,14980});
	badcellsManual.insert(badcellsManual.end(),{15298,15476,16477,16503,16505});
*/
    //..15o block 1
    //badcellsManual.insert(badcellsManual.end(),{14655,14622,14640,14728,14726,14593,14599,14600,14645,14646,14759,14776});
/*
 * //..15o block 2
    badcellsManual.insert(badcellsManual.end(),{6644,6655,10140,12036,12037,12038,12039,12040,12041,12926,13067,13066,13125});
    badcellsManual.insert(badcellsManual.end(),{13133,13483,13971,13978,14116,14118,14122,14411,14593,14599,14600,14606,14699});
	badcellsManual.insert(badcellsManual.end(),{14728,15158,15462,16309});
*/  //..15o block 3
	/*badcellsManual.insert(badcellsManual.end(),{292,294,297,301,13483, 13975, 14116, 14320, 14326,14593,14597,14621,14657,14671,14705});
	badcellsManual.insert(badcellsManual.end(),{14707,14716,14717,14728,14740,14748,14752,14759});
	 */
    //..15o block 4
	/*badcellsManual.insert(badcellsManual.end(),{3472,3473,3474,3475,3476,3477,3478,3479,3480,3481,3482,3483,3484,3485,3486,3487});
	badcellsManual.insert(badcellsManual.end(),{3520,3521,3522,3523,3524,3525,3526,3527,3528,3529,3530,3531,3532,3533,3534,3535});
	badcellsManual.insert(badcellsManual.end(),{3665,3666,3667,3668,3669,3670,3671,3672,3673,3674,3675,3676,3677,3678,3679});
	badcellsManual.insert(badcellsManual.end(),{3712,3713,3714,3715,3716,3717,3718,3719,3720,3721,3722,3723,3724,3725,3726,3727});
	badcellsManual.insert(badcellsManual.end(),{8848,8849,8850,8851,8852,8853,8854,8855,8856,8857,8858,8859,8860,8861,8862,8863});
	badcellsManual.insert(badcellsManual.end(),{8896,8897,8898,8899,8900,8901,8902,8903,8904,8905,89106,8907,8908,8909,8910});
	badcellsManual.insert(badcellsManual.end(),{11906,11907,11908,11909,11910,11911,11912,11913,11914,11915,11916,11917,11918,11919});
	badcellsManual.insert(badcellsManual.end(),{12034,12035,12036,12037,12038,12039,12040,12041,12042,12043,12044,12045,13469,13483,16427,16430});
*/

	//..If nversion=-1
	//..this is detected as a run-by-run analysis
	if(nversion == -1)nversion=runNum; //..If you do the analysis run by run - this might be helpful

	TStopwatch watch;
	watch.Start();

	BadChannelAna* Analysis;
	Analysis=new BadChannelAna(period,train,trigger,runNum,nversion,workDir,listName);

	//..Settings
	Analysis->SetExternalMergedFile(externalFile);
    Analysis->AddManualMasking(badcellsManual);
    //Analysis->SetStartEndCell(0,12288);     //..only EMCal
    //Analysis->SetStartEndCell(12288,17664);//..only DCal
    //Analysis->SetQAChecks(1);     //1= Perform QA checks - takes a long time! Prints all good cells for cross check
	//Analysis->SetPrintOutput(1);  //1= prints more information about excluded cells

	//. . . . . . . . . . . . . . . . . . . . . . . .
	//. . Add different period analyses
	//. . . . . . . . . . . . . . . . . . . . . . . .
	//..the range of sigmas should be selected such
	//..that one does not cut into the natural fluctuation over the modules
	Analysis->AddPeriodAnalysis(2, 5.5,0.1,0.3);  // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.0,0.1,0.3);  // energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.2,0.5);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.2,0.5);  // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.5,0.3,0.6);  //neu* hits in cell range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 4.5,0.3,0.6);  //neu* energy/hit in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,0.5,1.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,0.5,1.0);  // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.5,1.0,2.0);  //neu* hits in cell range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 4.5,1.0,2.0);  //neu* mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,1.0,4.0);  // hits in cell range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,4.0);  // mean energy in range Emin Emax
	Analysis->AddPeriodAnalysis(2, 5.5,1.0,10.0); // hits in cell in range Emin Emax
	Analysis->AddPeriodAnalysis(1, 4.5,1.0,10.0); // energy/hit in range Emin Emax


	//..special test for extra high energy fluctuations
	Analysis->AddPeriodAnalysis(1, 4.0,3.0,40.0); // energy/hit in cell in range Emin Emax
//	Analysis->AddPeriodAnalysis(2, 4.0,3.0,40.0); // hits in cell range Emin Emax
	//Analysis->AddPeriodAnalysis(1, 4.5,3.0,5.0);  // mean energy in range Emin Emax - cliff
	//Analysis->AddPeriodAnalysis(2, 5.5,3.0,5.0);  // hits in cell range Emin Emax   - cliff
	//Analysis->AddPeriodAnalysis(1, 4.5,5,20);     // energy/hit in range Emin Emax
	//Analysis->AddPeriodAnalysis(2, 5.5,5,20);     // hits in range Emin Emax

	//*test time stuff*/	Analysis->AddPeriodAnalysis(3, 6,-20,+20);// energy/hit in range Emin Emax

	//..Start the bad channel analysis
	Bool_t mergeOnly=0;//.. =1 do only merge and filter
	Analysis->Run(mergeOnly);

	watch.Stop();
	watch.Print();
}
