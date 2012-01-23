

using namespace std;
void printCuts(TString str){



  TString s = str(0,1);
  Int_t goodId= s.Atoi();
  s = str(1,1);
  Int_t v0FinderType= s.Atoi();
  s = str(2,1);
  Int_t eProbCut= s.Atoi();
  s = str(3,1);
  Int_t ededxSigmaCut= s.Atoi();
  s = str(4,1);
  Int_t pidedxSigmaCut= s.Atoi();
  s = str(5,1);
  Int_t piMomdedxSigmaCut= s.Atoi();
  s = str(6,1);
  Int_t chi2GammaCut= s.Atoi();
  s = str(7,1);
  Int_t singlePtCut= s.Atoi();
  s = str(8,1);
  Int_t clsTPCCut= s.Atoi();
  s = str(9,1);
  Int_t etaCut= s.Atoi();

  /*  cout<<"etaCut: "<<etaCut<<endl;
  cout<<"clsTPCCut: "<<clsTPCCut<<endl;
  cout<<"singlePtCut: "<<singlePtCut<<endl;
  cout<<"chi2GammaCut: "<<chi2GammaCut<<endl;
  cout<<"piMomdedxSigmaCut: "<<piMomdedxSigmaCut<<endl;
  cout<<"pidedxSigmaCut: "<<pidedxSigmaCut <<endl;
  cout<<"ededxSigmaCut: "<<ededxSigmaCut <<endl;
  cout<<"eProbCut: "<< eProbCut<<endl;
  cout<<"v0FinderType: "<<v0FinderType <<endl;
  cout<<"goodId: "<<goodId <<endl;
  */

  if(goodId !=9){
    cout<<"Analysis Cut Selection too short or does not start with 9"<<endl;
    return;
  }

  switch (v0FinderType){
  case 0:  // on fly V0 finder
    cout<<"Using the on fly v0 finder"<<endl;
    break;
  case 1:  // offline V0 finder
    cout<<"Using the offline v0 finder"<<endl;
    break;
  default:
    return;
  }
  switch(eProbCut){
  case 0:  // 0.
    cout<<"Prob electron is 0.000"<<endl;
    //    kGCprobElectron = 0.000;
    break;
  case 1:  // 0.001
    cout<<"Prob electron is 0.001"<<endl;
    //    kGCprobElectron = 0.001;
    break;
  case 2:  // 0.01
    cout<<"Prob electron is 0.01"<<endl;
    //    kGCprobElectron = 0.01;
    break;
  default:
    return;
  }

  switch(ededxSigmaCut){
  case 0: // -10,10
    cout<<"Using sigma dedx [-10,10]"<<endl;
    //    kGCPIDnSigmaBelowElectronLine=-10;
    //    kGCPIDnSigmaAboveElectronLine=10;
    break;
  case 1: // -5,5 
    cout<<"Using sigma dedx [-5,5]"<<endl;
    //    kGCPIDnSigmaBelowElectronLine=-5;
    //   kGCPIDnSigmaAboveElectronLine=5;
    break;
  case 2: // -3,5
    cout<<"Using sigma dedx [-3,5]"<<endl;
    //    kGCPIDnSigmaBelowElectronLine=-3;
    //    kGCPIDnSigmaAboveElectronLine=5;
    break;
  default:
    return;
  }
  
  switch(pidedxSigmaCut){
  case 0:  // -10
    cout<<"using pidedxsigmacut: -10"<<endl;
    //   kGCPIDnSigmaAbovePionLine=-10;
    break;
  case 1:   // 0
    cout<<"using pidedxsigmacut: 0"<<endl;
    //    kGCPIDnSigmaAbovePionLine=0;
    break;
  case 2:  // 1
    cout<<"using pidedxsigmacut: 1"<<endl;
    //   kGCPIDnSigmaAbovePionLine=1;
    break;
  default:
    return;
  }
  
  switch(piMomdedxSigmaCut){
  case 0:  // 0.5 GeV
    cout<<"piMomdedxSigmaCut: 0.5"<<endl;
    //  kGCPIDMinPnSigmaAbovePionLine=0.5;
    break;
  case 1:  // 1. GeV
    cout<<"piMomdedxSigmaCut: 1"<<endl;
    //   kGCPIDMinPnSigmaAbovePionLine=1.;
    break;
  case 2:  // 1.5 GeV
    cout<<"piMomdedxSigmaCut: 1.5"<<endl;
    //   kGCPIDMinPnSigmaAbovePionLine=1.5;
    break;
  default:
    return;
  }
  
  switch(chi2GammaCut){
  case 0: // 100
    cout<<"chi2CutConversion = 100."<<endl;
    //   kGCchi2CutConversion = 100.;
    break;
  case 1:  // 50
    cout<<"chi2CutConversion = 50."<<endl;
    //   kGCchi2CutConversion = 50.;
    break;
  case 2:  // 30
    cout<<"chi2CutConversion = 30."<<endl;
    //   kGCchi2CutConversion = 30.;
    break;
  default:
    return;
  }

  switch(singlePtCut){
  case 0: // 0.050 GeV
    cout<<"kGCsingleptCut = 0,050"<<endl;
    //   kGCsingleptCut = 0.050;
    break;
  case 1:  // 0.100 GeV
    cout<<"kGCsingleptCut = 0,100"<<endl;
    //   kGCsingleptCut = 0.100;
    break;
  case 2:  // 0.150 GeV
    cout<<"kGCsingleptCut = 0,150"<<endl;
    //   kGCsingleptCut = 0.150;
    break;
  case 3:  // 0.200 GeV
    cout<<"kGCsingleptCut = 0,200"<<endl;
    //   kGCsingleptCut = 0.200;
    break;
  default:
    return;
 }

  switch(clsTPCCut){
  case 0: // 0 
    cout<<"kGCminClsTPCCut = 0"<<endl;
    //   kGCminClsTPCCut= 0.;
    break;
  case 1:  // 70 
    //   kGCminClsTPCCut= 70.;
    cout<<"kGCminClsTPCCut = 70"<<endl;
    break;
  case 2:  // 80 
    //   kGCminClsTPCCut= 80.;
    cout<<"kGCminClsTPCCut = 80"<<endl;
    break;
  case 3:  // 100 
    //   kGCminClsTPCCut= 100.;
    cout<<"kGCminClsTPCCut = 100"<<endl;
    break;
  default:
    return;
  }

  switch(etaCut){
  case 0: // 0.9 
    cout<<"eta 0.9"<<endl;
    //   kGCetaCut    = 0.9;
    //    kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  case 1:  // 1.2
    cout<<"eta 1.2"<<endl;
    //   kGCetaCut    = 1.2;
    //   kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  case 2:  // 1.4
    cout<<"eta 1.4"<<endl;
    //    kGCetaCut    = 1.4;
    //   kGCLineCutZRSlope = tan(2*atan(exp(-kGCetaCut)));
    break;
  default:
    return;
  }


}
