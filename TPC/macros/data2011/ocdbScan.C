/// \file ocdbScan.C

TChain * chain0=0;
TChain * chain1=0;

void Init(){  
  //  ls /hera/alice/local/benchmark/TestFor2011/r64765/000*/cpass1/dcsTime.root > /u/miranov/cpass1Calib.list
  //  ls /hera/alice/local/benchmark/TestFor2011/r64765/000*/cpass0/dcsTime.root > /u/miranov/cpass0Calib.list
  
  chain0 = AliXRDPROOFtoolkit::MakeChain("/u/miranov/cpass0Calib.list","dcs",0,10000);
  chain1 = AliXRDPROOFtoolkit::MakeChain("/u/miranov/cpass1Calib.list","dcs",0,10000);
}


void CheckHVCorrection(){
  //
  //
  // 0.)
  chain0->Draw("gainMIP:ptrel0","vdriftITS<0.01");
  // 1.)
  TStatToolkit::MakeGraphSparse(chain0,"gainMIP:ptrel0","vdriftITS<0.01",25,2)->Draw("alp");

}
