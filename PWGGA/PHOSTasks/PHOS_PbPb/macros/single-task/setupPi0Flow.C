//-*- Mode: C++ -*-
// $Id$

/// @file   setupPi0Flow.C
/// @author Matthias.Richter@cern.ch, modified for PHOSPi0Flow by Henrik.Qvigstad@cern.ch
/// @date   2012-09-12
/// @brief  Setup the environment specifically for running Pi0Flow analysis
///
/// Helper macro to initialize the environment for Pi0Flow. The macro can just
/// prepend other macros like run-single-task in the command line.
/// Usage:
/// aliroot -b -q -l setupPi0Flow.C'("localAodDirectory", nFiles)' run-single-task.C'(...)'
/// aliroot -b -q -l setupPi0Flow.C'("lhcPeriod", "mcPeriod")' run-single-task.C'(...)'
///
/// Example:
/// aliroot -b -q -l setupPi0Flow.C run-single-task.C'(...)'
///
/// The macro has the following tasks:
/// - load the necessary libraries, in order to have those library names also
///   available for the alien handler initialization, a specific configuration
///   object is created
/// - setting a default analysis name via a configuration object
/// - the optional parameter 'localAodDirectory' allows to create an input chain from
///   local AODs; either a single AliAOD.root, or a folder containing directories
///   named "1, 2, ...", the number of directories is specified as parameter
///   nFiles
/// - loading the runs defined by AddGoodRuns of vertexingHF using lhcPeriod and
///   optional mcPeriod

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment specific for Pi0Flow
//
const char* analysisName="PHOSPi0Flow";
const char* includePath="-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS -I$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb";
const char* libraryDependencies="libCore.so libTree.so libGeom.so libVMC.so libPhysics.so libMinuit.so libGui.so libXMLParser.so libMinuit2.so libProof.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libANALYSISalice.so libCDB.so libRAWDatabase.so libSTEER.so libCORRFW.so libPHOSUtils.so libPHOSbase.so libPHOSpi0Calib.so libPHOSrec.so libPHOSshuttle.so libPHOSsim.so libPWGGAPHOSTasks.so libTENDER.so libTRDbase.so libVZERObase.so libVZEROrec.so libTENDERSupplies.so";

void setupPi0Flow(const char* localAodDirectory, int nFiles, const char* lhcPeriod, const char* mcProd="")
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
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/macros/single-task/AddGoodRuns.C");
    int nruns=AddGoodRuns(alienHandler, lhcPeriod, mcProd);
    if (nruns<=0) {
      ::Error("setupPi0Flow.C", Form("can not find any good runs for period %s", lhcPeriod));
      return;
    }
    gDirectory->Add(alienHandler);
    ::Info("setupPi0Flow.C", Form("setting up alien plugin '%s' for period %s\n>>>>> please use '%s' as input parameter for run-single-task.C <<<<<<", alienHandlerName.Data(), lhcPeriod, alienHandlerName.Data()));

  } else if (localAodDirectory) {
    //
    // create AOD tree from local files
    // the created object is added automatically to gDirectory and can be fetched
    // from there later
    //
    //gROOT->LoadMacro("$ALICE_ROOT/PWGHF/vertexingHF/MakeAODInputChain.C");
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPb/macros/single-task/MakeAODInputChain.C");
    TChain* chain = 0;
    chain = MakeAODInputChain(localAodDirectory, nFiles);
    // TString aodPathName(localAodDirectory);
    // if (!aodPathName.EndsWith("/")) aodPathName+="/";
    // aodPathName+="AliAOD.root";
    TChain* chain;
    // if (gSystem->AccessPathName(aodPathName)==0) {
    //   // Create a chain with one set of AliAOD.root and AliAOD.VertexingHF.root. The set needs 
    //   // to be located in the same folder as you run from (physically or linked)
    //   ::Info("setupPi0Flow.C", Form("make chain from single chunk %s", aodPathName));
    //   chain = MakeAODInputChain(localAodDirectory ,1, -1);
    // } else {
    //   // Assume several folders containing different AODs. 
    //   // The AODs need to be in folders named 1, 2,...
    //   ::Info("setupPi0Flow.C", Form("make chain from directory %s", localAodDirectory));
    //   chain=MakeAODInputChain(localAodDirectory, nFiles);
    // }
    ::Info("setupPi0Flow.C", Form("local AOD chain: %d entries", chain->GetEntries()));
  }
}

void setupPi0Flow(const char* lhcPeriod=NULL, const char* mcProd="")
{
  // Grid mode with optional calling of AddGoodRuns for specified
  // period
  setupPi0Flow(NULL, 0, lhcPeriod, mcProd);
}

void setupPi0Flow(const char* localAodDirectory, int nFiles)
{
  // local mode for AOD data
  setupPi0Flow(localAodDirectory, nFiles, NULL);
}
