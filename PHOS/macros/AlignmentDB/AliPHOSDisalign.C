/* $Id$ */

// Run simulations with PHOS alignment parameters read from CDB.
// CDB storage should be defined by environment variable CDB_PATH before 
// the aliroot session, otherwise the default hard-coded geometry will be used
// The alignment object is located in the local file
// $CDB_PATH/PHOS/Alignment/Geometry/RunXXX.root
//
// Author: Yuri Kharlov
// 7 March 2006

void AliPHOSDisalign(Int_t nevents=1)
{
  AliSimulation sim ; 
  sim.SetMakeSDigits("PHOS") ;
  sim.SetMakeDigits("PHOS") ;
  sim.Run(nevents) ;  
}
