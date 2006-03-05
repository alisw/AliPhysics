/* $Id$ */

// Run simulations with PHOS alignment parameters read from CDB
// located in the local file
// InitAlignDB/PHOS/Align/IdealGeometry2007/RunXXX.root
//
// Author: Yuri Kharlov
// 4 March 2006

void AliPHOSDisalign(char *path="InitAlignDB", Int_t nevents=1)
{

  //Load alignment database into aliroot session
  
  TString localPath="local://";
  localPath += path;
  AliPHOSAlignData* alignda  = (AliPHOSAlignData*)(AliCDBManager::Instance()
    ->GetStorage(localPath)->Get("PHOS/Alignment/Geometry",0)
    ->GetObject());
  alignda->Print();

  // Create PHOS geometry instance with alignment object read from CDB

  AliPHOSGeometry* geom = AliPHOSGeometry::GetInstance("IHEP","",alignda);

  // Run simulations

  AliSimulation sim ; 
//   sim.SetMakeSDigits("PHOS") ;
//   sim.SetMakeDigits("PHOS") ;
  sim.Run(nevents) ;  
}
