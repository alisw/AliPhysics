//   It is recommended to train not all momenta at on eturn. The time consumption
//   for a training loop is very large! 

//   Please be careful if you want to train the networks. The networks will be saved
//   in a file ("./Networks_%d", date(yyyymmdd)). You should create seperate directories 
//   where the training for the several momentum bins should be done. Otherwise the 
//   networks will be deleted by the following momentum bin training.

#ifndef __CINT__
#include "TSystem.h"
#include "$ALICE_ROOT/TRD/qaRec/AliTRDpidRefMaker.h"
#endif

Int_t runTraining(){

  gSystem -> Load("libANALYSIS.so");
  gSystem -> Load("libTRDqaRec.so");

  AliTRDpidRefMaker *ref = new AliTRDpidRefMaker();
  ref->SetDebugLevel(2);

//   Sets the momentum bin that should be trained:
//   Bin  0 =  0.6 GeV/c
//   Bin  1 =  0.8 GeV/c
//   Bin  2 =  1.0 GeV/c
//   Bin  3 =  1.5 GeV/c
//   Bin  4 =  2.0 GeV/c
//   Bin  5 =  3.0 GeV/c
//   Bin  6 =  4.0 GeV/c
//   Bin  7 =  5.0 GeV/c
//   Bin  8 =  6.0 GeV/c
//   Bin  9 =  8.0 GeV/c
//   Bin 10 = 10.0 GeV/c
//   Bin 11 = All
  ref->SetTrainMomBin(1);

//   Sets if the network should be trained
  ref->SetDoTraining(1);

//   Sets if the training should be continued or an untrained network will be taken.
//   If the training should be continued with a pretrained network, 
//   please uncomment this line
//   ref->SetContinueTraining(1);

//   Sets the path for old networks. If the training should be continued 
//   with a pretrained network, please uncomment this line and set the date (yyyymmdd).
//   ref->SetTrainPath(20081022);

//   Do the training
  ref->PostProcess();
  
  return 1;

}
