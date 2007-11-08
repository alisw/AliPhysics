// $Id$
/*
 * Example macro to run the HLT Conformal mapping tracker embedded into
 * AliRoot simulation. The reconstruction is done from the TPC digits.
 *
 * aliroot -b -q sim-hlt-tpc.C | tee sim-hlt-tpc.log
 *
 * The chain to be run is defined by the macro given to the parameter
 * 'config='
 *
 * Matthias.Richter@ift.uib.no
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
  //sim.SetRunHLT("libAliHLTTPC.so loglevel=0x3c config=conf-tpc-writer.C chains=sink1");
  sim.SetRunHLT("libAliHLTTPC.so loglevel=0x3c config=conf-tpc-esd.C chains=esd-writer");
  sim.Run();
}
