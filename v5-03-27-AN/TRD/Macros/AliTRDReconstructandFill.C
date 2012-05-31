#if !defined( __CINT__) || defined(__MAKECINT__)


#include <Riostream.h>
#include <TSystem.h>
#include "AliReconstruction.h"
#include "../TRD/AliTRDCalibra.h"
#include <TStopwatch.h>

#endif



void AliTRDReconstructandFill() 
{
  //
  // This macro fills 2d histo or vectors during the reconstruction
  // If it is vectors, it fits them directly after the reconstruction
  // and writes the result in the file coeftest.root
  //
 
  TStopwatch timer;
  timer.Start();

  ////Set the parameters of AliTRDCalibra***************
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();

  ////What do you want to use?
  calibra->SetMITracking(); //Offline tracking
  //calibra->Setmcmtracking();
  
  
  ////Do you want to try the correction due to the angles of the tracks for mcm tracklets?
  calibra->SetMcmCorrectAngle();
  
  ////What do you want to fill?
  calibra->SetCH2dOn();//relative gain calibration
  calibra->SetPH2dOn();//drift velocity and time0 calibration
  calibra->SetPRF2dOn();//Pad Response Function calibration


  ////How do you want to store the infos?
  calibra->SetVector2d();//vector method
  calibra->SetHisto2d();//2Dhistograms

  ////Which mode do you want?
  calibra->SetNz(2,2);//For the PRF z direction
  calibra->SetNrphi(2,2);//For The PRF rphi direction
  calibra->SetNz(0,0);//For the gain z direction
  calibra->SetNrphi(0,0);//For the gain rphi direction
  calibra->SetNz(1,3);//For the drift velocity and time0 z direction
  calibra->SetNrphi(1,3);//For the drift velocity and time 0 rphi direction

  ////How many bins?
  calibra->SetNumberBinCharge(100);
  calibra->SetNumberBinPRF(20);
  
  
  ////Do you want to accept more tracks?
  calibra->SetProcent(1.2);//For the gain if one group has a signal above 1.2 the other group then fill
  calibra->SetDifference(10);//For the drift velocity if one group has at least 10 time bins then fill
  calibra->SetNumberClusters(18);//For mcm tracklets only fill only with tracklet with at least 18 clusters
  
  ////Do you want to take only the middle pad for gain or Vdrift?
  //calibra->SetTraMaxPad();
  
  //Do you want to apply more strict cut on the clusters for the PRF calibration?
  calibra->SetThresholdClusterPRF1(2);//The neighbors pads must have a signal smaller than 2 ADC counts
  calibra->SetThresholdClusterPRF2(10);//The 3 pads in the cluster must have a signal above 10 ADC counts
  
  
  ////What do you want to write?
  calibra->SetWrite(0);//For the gain
  calibra->SetWrite(1);//For the average pulse height
  calibra->SetWrite(2);//For the PRF
  
  
  ////If you want to change the name of the file where it is stored (not very good)
  //calibra->SetWriteName("test.root");
  
  
  //Begin the reconstruction
  AliReconstruction rec;
  rec.SetGAliceFile("galice.root"); 
  rec.SetLoadAlignFromCDB(kFALSE);
  rec.SetRunHLTTracking(kFALSE);
  rec.SetFillESD("");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.Run();
  timer.Stop();
  timer.Print();
  calibra->Write2d();    


  TStopwatch timerfit;
  timerfit.Start();
  ////Fit directly after having filling****

  ////Do you want to try with less statistics?
  calibra->SetMinEntries(10);//If there is at least 10 entries in the histo, it will fit

  ////Do you want to write the result?
  calibra->SetWriteCoef(0);//gain
  calibra->SetWriteCoef(1);//time 0 and drift velocity
  calibra->SetWriteCoef(2);//PRF

  ////Do you want to change the name of the file (TRD.coefficient.root)?
  calibra->SetWriteNameCoef("coeftest.root");

  ////Do you want to see something?
  calibra->SetDebug(1);

  ////Do you want to fit?
  calibra->FitPHOnline(); 
  calibra->FitCHOnline(); 
  calibra->FitPRFOnline();  
  
  
  timerfit.Stop();
  timerfit.Print();

}
