// $Id$
/*
 * @file sim-hlt-tpc.C
 * @brief HLT Conformal mapping tracker embedded into AliRoot simulation.
 *
 * Example macro to run the HLT Conformal mapping tracker embedded into
 * AliRoot simulation. The reconstruction is done from the TPC digits.
 *
 * Usage: aliroot -b -q sim-hlt-tpc.C | tee sim-hlt-tpc.log
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * The macro asumes the data to be already simulated. If it should run
 * within the initial simulation, comment the corresponding functions
 * below (SetRunGeneration etc.)
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tpc
 */
{
  // this is just a tool to switch the logging systems
  AliHLTLogging log;
  //log.SwitchAliLog(0);

  AliSimulation sim;

  // switch of simulation and data generation
  // comment all that stuff to also simulate the events and data
  sim.SetRunGeneration(kFALSE);
  sim.SetMakeDigits("");
  sim.SetMakeSDigits("");
  sim.SetMakeDigitsFromHits("");
  sim.SetMakeTrigger("");

  // set the options for the HLT simulation
  sim.SetRunHLT("libAliHLTUtil.so libAliHLTTPC.so loglevel=0x7c "
		"config=$ALICE_ROOT/HLT/TPCLib/macros/conf-tpc-esd.C chains=sink-esd,sink-clusters");
  sim.Run();
}
