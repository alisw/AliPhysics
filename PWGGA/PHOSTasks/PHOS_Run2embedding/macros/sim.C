/*
 * AliDPG - ALICE Experiment Data Preparation Group
 * Simulation steering script
 *
 */

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

void sim() 
{

  // run number
  Int_t runNumber = -1;
  if (gSystem->Getenv("CONFIG_RUN"))
    runNumber = atoi(gSystem->Getenv("CONFIG_RUN"));
  if (runNumber <= 0) {
    printf("Invalid run number: %d \n", runNumber);
    abort();
  }

  // number of events configuration
  Int_t nev = 200;
  if (gSystem->Getenv("CONFIG_NEVENTS"))
    nev = atoi(gSystem->Getenv("CONFIG_NEVENTS"));

  // simulation configuration
  gROOT->LoadMacro("Sim/SimulationConfig.C");
  Int_t simulationConfig = kSimulationDefault;
  if (gSystem->Getenv("CONFIG_SIMULATION")) {
    Bool_t valid = kFALSE;
    for (Int_t isim = 0; isim < kNSimulations; isim++)
      if (strcmp(gSystem->Getenv("CONFIG_SIMULATION"), SimulationName[isim]) == 0) {
        simulationConfig = isim;
        valid = kTRUE;
        break;
      }
    if (!valid) {
      printf(">>>>> Unknown simulation configuration: %s \n", gSystem->Getenv("CONFIG_SIMULATION"));
      abort();
    }
  }

  /* initialisation */
  AliSimulation sim("Sim/Config.C");

  /* configuration */
  SimulationConfig(sim, simulationConfig, runNumber);

  /* run */
  sim.Run(nev);

}

