// $Id$
/**
 * @file sim-hlt-rawddl.C
 * @brief Publishing of RAW DDL data into the HLTOUT.
 *
 * Example macro to produce ddl raw data blocks in the HLTOUT. The HLT
 * chain is run embedded into AliRoot simulation. 
 *
 * The example publishes the ITSSDD data into the HLTOUT according to
 * the configuration in the conf-hlt-rawddl.C macro. The normal digit
 * to raw conversion is used to have the raw data available.
 * \b Note: if you want to change this example for another detector,
 * you need to change the detector in both this macro and the
 * configuration.
 *
 * Usage: aliroot -b -q sim-hlt-rawddl.C | tee sim-hlt-rawddl.log
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro assumes the data to be already simulated. If it should run
 * within the initial simulation, comment the corresponding functions
 * below (SetRunGeneration etc.)
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void sim_hlt_rawddl() {  

  AliSimulation sim;

  // switch of simulation and data generation
  // comment all that stuff to also simulate the events and data
  sim.SetRunGeneration(kFALSE);
  sim.SetMakeDigits("");
  sim.SetMakeSDigits("");
  sim.SetMakeDigitsFromHits("");
  sim.SetMakeTrigger("");

  // write the raw data for the ITS since we want to publish those
  // write HLT raw data since we want to replace the original	
  // detector data from the HLTOUT
  sim.SetWriteRawData("HLT");
  
  // the normal simulation sets the specific storage for the GRP entry
  if (gSystem->AccessPathName("GRP/GRP/Data")) {
    cerr << "*********************************************************" << endl;
    cerr << "error: no GRP entry found in the currect directory, simulation might be incomplete. Skip setting specific storage for GRP entry" << endl;
    cerr << "*********************************************************" << endl << endl;
  } else {
    sim.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  }
  
  // set the options for the HLT simulation:
  // libAliHLTUtil.so libAliHLTSample.so
  //     loads the specified libraries since the HLT chain will use components
  //     from those two
  // loglevel=
  //     the internal logging level in the HLT, use 0x7c for higher verbosity
  // config=<file>
  //     the configuration to be run
  // chains=<chain>
  //     run the specified chains, defined in the configuration macro
  // rawfile=<path>
  //     in this case we want to just forward the DDL data blocks
  //     to the HLTOUT. We need to set the location of the formerly generated
  //     rawfiles with the rawfile
 
  sim.SetRunHLT("libAliHLTUtil.so libAliHLTSample.so loglevel=0x7c rawfile=./ config=$ALICE_ROOT/HLT/exa/conf-hlt-rawddl.C chains=publisher");
  sim.Run();
}
