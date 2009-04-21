//   This macro helps to decide if the training was successful and in addition
//   if further training of the network is necessary. The shape of the training
//   progress should have a first deep stall and reaching a plateau. After several
//   training epochs the plateau will become a sloping line. If this slope results 
//   again in a flat plateau the training is probably done.
//   The networks were saved using the date they were trained (yyyymmdd).
// 
//   Please enter the date of the network you want to monitor and the momentum bin..
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

#ifndef __CINT__
#include "TSystem.h"
#include "qaRec/AliTRDpidRefMaker.h"
#endif

Int_t MonitorTraining(Int_t bin, Int_t Date){

  gSystem -> Load("libANALYSIS.so");
  gSystem -> Load("libTRDqaRec.so");
  
  AliTRDpidRefMaker *ref = new AliTRDpidRefMaker();
  ref->SetDebugLevel(2);

  ref->SetDate(Date);
  ref->MakeTrainingLists();
  ref->MonitorTraining(bin);
  
  return 1;

}
