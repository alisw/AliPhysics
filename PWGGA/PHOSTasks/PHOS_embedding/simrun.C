// simrun.C
{
// set job and simulation variables as :
// root.exe -b -q simrun.C  --run <x> --event <y> --process <kPythia6/kPhojet/kPythia6ATLAS_Flat/kPythia6D6T> --field <kNoField/k5kG> --energy <900/2360/10000>

  cout<<">>>>> SIMULATION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q sim.C > sim.log 2>&1");
  gSystem->Exec("mv syswatch.log simwatch.log");
  gSystem->Exec("rm *ITS*.root  *TPC*.root *TOF*.root *TRD*.root *PMD*.root ") ;
  cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q rec.C > rec.log 2>&1");
  gSystem->Exec("mv syswatch.log recwatch.log");
  cout<<">>>>> AOD <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q CreateAOD.C > aod.log 2>&1");

}
