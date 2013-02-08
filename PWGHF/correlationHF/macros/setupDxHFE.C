//-*- Mode: C++ -*-
// $Id$

/// @file   setupDxHFE.C
/// @author Matthias.Richter@cern.ch
/// @date   2012-09-12
/// @brief  Setup the environment specifically for running DxHFE analysis
///
/// Helper macro to initialize the environment for DxHFE. The macro can just
/// prepend other macros like run-single-task in the command line.
/// Usage:
/// aliroot -b -q -l setupDxHFE.C'("localAodDirectory", nofDirectories)' run-single-task.C'(...)'
/// aliroot -b -q -l setupDxHFE.C'("lhcPeriod", "mcPeriod")' run-single-task.C'(...)'
///
/// Example:
/// aliroot -b -q -l setupDxHFE.C run-single-task.C'(...)'
///
/// The macro has the following tasks:
/// - load the necessary libraries, in order to have those library names also
///   available for the alien handler initialization, a specific configuration
///   object is created
/// - setting a default analysis name via a configuration object
/// - the optional parameter 'localAodDirectory' allows to create an input chain from
///   local AODs; either a single AliAOD.root, or a folder containing directories
///   named "1, 2, ...", the number of directories is specified as parameter
///   nofDirectories
/// - loading the runs defined by AddGoodRuns of vertexingHF using lhcPeriod and
///   optional mcPeriod

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment specific for DxHFE
//
const char* analysisName="DxHFECorrelation";
const char* includePath="-I$ALICE_ROOT/PWGHF/correlationHF -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/hfe";
const char* libraryDependencies=
  "libSTEERBase "
  "libESD "
  "libAOD "
  "libANALYSIS "
  "libANALYSISalice "
  "libPWGflowBase.so "
  "libPWGflowTasks.so "
  "libCORRFW.so "
  "libPWGHFvertexingHF.so "
  "libPWGHFhfe.so "
  ;
// Note: "libPWGHFcorrelationHF.so " shouldn't be added here. If the library is
// loaded already, compilation of source files by Cint does not have any effect.
// The already existing class implementations are not overidden. To just load
// all necessary libraries to make the classes available in aliroot use macro
// LoadLibraries.C

void setupDxHFE(const char* localAodDirectory, int nofDirectories, const char* lhcPeriod, const char* mcProd="")
{
  //
  // adding include path and libraries
  //
  gSystem->AddIncludePath(includePath);
  TString libraries=libraryDependencies;
  TObjArray* pTokens=libraries.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      if (gSystem->Load(pTokens->At(i)->GetName())==0) {
	cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }
  libraries="";

  //
  // allow run-single-task to fetch the analysis name and library names
  //
  if (gDirectory) gDirectory->Add(new TNamed("analysis_name", analysisName));
  if (gDirectory) gDirectory->Add(new TNamed("analysis_libraries", libraryDependencies));

  if (lhcPeriod) {
    //
    // setting up the runs for the dpecified period
    //
    TString alienHandlerName(analysisName); alienHandlerName+="Handler";
    AliAnalysisAlien* alienHandler=new AliAnalysisAlien(alienHandlerName);
    gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/AddGoodRuns.C");
    int nruns=AddGoodRuns(alienHandler, lhcPeriod, mcProd);
    if (nruns<=0) {
      ::Error("setupDxHFE.C", Form("can not find any good runs for period %s", lhcPeriod));
      return;
    }
    gDirectory->Add(alienHandler);
    ::Info("setupDxHFE.C", Form("setting up alien plugin '%s' for period %s\n>>>>> please use '%s' as input parameter for run-single-task.C <<<<<<", alienHandlerName.Data(), lhcPeriod, alienHandlerName.Data()));

  } else if (localAodDirectory) {
    //
    // create AOD tree from local files
    // the created object is added automatically to gDirectory and can be fetched
    // from there later
    //
    gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/MakeAODInputChain.C");
    TString aodPathName(localAodDirectory);
    if (!aodPathName.EndsWith("/")) aodPathName+="/";
    aodPathName+="AliAOD.root";
    if (gSystem->AccessPathName(aodPathName)==0) {
      // Create a chain with one set of AliAOD.root and AliAOD.VertexingHF.root. The set needs 
      // to be located in the same folder as you run from (physically or linked)
      ::Info("setupDxHFE.C", Form("make chain from single chunk %s", aodPathName.Data()));
      TChain* chain = MakeAODInputChain(localAodDirectory ,1, -1);
    } else {
      // Assume several folders containing different AODs. 
      // The AODs need to be in folders named 1, 2,...
      aodPathName=localAodDirectory;
      if (!aodPathName.EndsWith("/")) aodPathName+="/";
      ::Info("setupDxHFE.C", Form("make chain from directory %s", aodPathName.Data()));
      chain=MakeAODInputChain(aodPathName, 1, nofDirectories);
    }
    ::Info("setupDxHFE.C", Form("local AOD chain: %d entries", chain->GetEntries()));
  }
}

void setupDxHFE(const char* lhcPeriod=NULL, const char* mcProd="")
{
  // Grid mode with optional calling of AddGoodRuns for specified
  // period
  setupDxHFE(NULL, 0, lhcPeriod, mcProd);
}

void setupDxHFE(const char* localAodDirectory, int nofDirectories)
{
  // local mode for AOD data
  setupDxHFE(localAodDirectory, nofDirectories, NULL);
}
