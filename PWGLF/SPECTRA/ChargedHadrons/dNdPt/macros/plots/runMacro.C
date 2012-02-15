//------------------------------------------------------------------------------
// runMacro.C
//
// steering macro to create figures for paper 3
// settings are in settings.C
//------------------------------------------------------------------------------


{
//
// Init
//
using namespace std;
gROOT->Reset();
gROOT->SetStyle("Plain");

//
// load Settings & define Variables
//
gROOT->LoadMacro("settings.C");
gROOT->LoadMacro("defineVariables.C");
cout << "---------------------------------------------------------" << endl;
cout << "using fit function (nsd) " << endl;
cout << fitNsd->GetExpFormula() << endl;
cout << "---------------------------------------------------------" << endl;


//
// graphics and plot options
//
gStyle->SetTextFont(textFont);
gStyle->SetTitleFont(titleFont);
gStyle->SetTitleFont(titleFont,"xy");
gStyle->SetLabelFont(labelFont,"xyz");
gStyle->SetLabelSize(labelSize);
gStyle->SetTitleSize(titleSize);
gStyle->SetTitleFontSize(titleFontSize);
gStyle->SetMarkerSize(markerSize);
gStyle->SetHatchesSpacing(0.8);
gStyle->SetHatchesLineWidth(2.0);

//
// load macros
//
gROOT->LoadMacro("divide.C");
gROOT->LoadMacro("setAttrib.C");
gROOT->LoadMacro("logoPrelim.C");

gROOT->LoadMacro("readAliceNsd.C");
gROOT->LoadMacro("readAliceInel.C");
gROOT->LoadMacro("readAliceYield.C");
gROOT->LoadMacro("readAtlas.C");
gROOT->LoadMacro("readCms.C");
gROOT->LoadMacro("readUa1.C");
gROOT->LoadMacro("readPhojet.C");
gROOT->LoadMacro("readPythia109.C");
gROOT->LoadMacro("readPythia306.C");
gROOT->LoadMacro("readPythia320.C");

gROOT->LoadMacro("makePlotsAlice3.C");
gROOT->LoadMacro("makeCompNSD.C");
gROOT->LoadMacro("makeCompYield.C");
gROOT->LoadMacro("makeCompInel.C");

gROOT->LoadMacro("storeOutput.C");


//
// read data
//
readAliceNsd();
readAliceInel();
readAliceYield();

readAtlas();
readCms();
readUa1();
readPhojet();
readPythia109();
readPythia306();
readPythia320();

//
// pt range to plot
//
Double_t minPt = 0.1;
Double_t maxPt = 10;

//
// generate plots & store output
//
makePlotsAlice3(); // figure 2 in paper
makeCompYield(); // figure 3 (b)
makeCompInel(); // figure 5

// different pt range for atlas comparison
maxPt = 20; 
makeCompNSD(); // figure 3 (a)

storeOutput();


}
