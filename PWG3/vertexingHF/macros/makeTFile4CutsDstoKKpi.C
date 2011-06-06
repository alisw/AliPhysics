#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsDstoKKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//macro to make a .root file which contains an AliRDHFCutsDstoKKpi with loose set of cuts (for significance maximization) and TParameter with the tighest value of these cuts
//Needed for AliAnalysisTaskSEDs, AliCFTaskVertexingHF3Prong, AliAnalysisTaskSESignificance

//Use:
//Set hard coded commented with //set this!!

//.L makeTFile4CutsDstoKKpi.C
// makeInputAliAnalysisTaskSEDs()
// makeInputAliAnalysisTaskSESignificanceMaximization()

//similar macros for the other D mesons

//Author: Chiara Bianchin, cbianchi@pd.infn.it
//        Giacomo Ortona, ortona@to.infn.it
//        Renu Bala [Dplus Analysis and CF]

//Modified for Ds meson: G.M. Innocenti, innocent@to.infn.it

void makeInputAliAnalysisTaskSEDs(){

 //  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -I$ALICE_ROOT/PWG3/vertexingH/macros -g"); 

    
    AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    //default
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    //esdTrackCuts->SetMinNClustersITS(4); // default is 5
    esdTrackCuts->SetMinNClustersTPC(70);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny); 
    // default is kBoth, otherwise kAny
    esdTrackCuts->SetMinDCAToVertexXY(0.);
    esdTrackCuts->SetPtRange(0.3,1.e10);
    
    
    
    const Int_t nptbins=7;
    Float_t* ptbins;
    ptbins=new Float_t[nptbins+1];
    ptbins[0]=2.;
    ptbins[1]=3.;
    ptbins[2]=4.;
    ptbins[3]=5.;
    ptbins[4]=6.;
    ptbins[5]=8.;
    ptbins[6]=12.;
    ptbins[7]=999999999999.;
    
    
    const Int_t nvars=16;
   
    
    Float_t** prodcutsval;
    prodcutsval=new Float_t*[nvars];
    for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}  
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      prodcutsval[0][ipt]=0.2;
      prodcutsval[1][ipt]=0.3;
      prodcutsval[2][ipt]=0.3;
      prodcutsval[3][ipt]=0.;
      prodcutsval[4][ipt]=0.;
      prodcutsval[5][ipt]=0.005;
      prodcutsval[6][ipt]=0.06;
      prodcutsval[7][ipt]=0.0;
      prodcutsval[8][ipt]=0.;
      prodcutsval[9][ipt]=0.7;
      prodcutsval[10][ipt]=0.;
      prodcutsval[11][ipt]=1000.0;
      prodcutsval[12][ipt]=0.1;
      prodcutsval[13][ipt]=0.1;
      prodcutsval[14][ipt]=0.;
      prodcutsval[15][ipt]=1.;
      
    }
    
    
    Float_t** anacutsval;
    anacutsval=new Float_t*[nvars];
  
    for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}  
    for(Int_t ipt=0;ipt<nptbins;ipt++){
      
      anacutsval[0][ipt]=0.2;
      anacutsval[1][ipt]=0.3;
      anacutsval[2][ipt]=0.3;
      anacutsval[3][ipt]=0.;
      anacutsval[4][ipt]=0.;
      anacutsval[5][ipt]=0.005;
      anacutsval[6][ipt]=0.06;
      anacutsval[7][ipt]=0.0;
      anacutsval[8][ipt]=0.;
      anacutsval[9][ipt]=0.7;
      anacutsval[10][ipt]=0.;
      anacutsval[11][ipt]=1000.0;
      anacutsval[12][ipt]=0.1;
      anacutsval[13][ipt]=0.1;
      anacutsval[14][ipt]=0.;
      anacutsval[15][ipt]=1.;
      
   
    }
    
        
    /*    
    Spiegazione della selezione di ALiRDHFCutsDstoKKpi      condizione di rejection

    0           "inv. mass [GeV]",                          invmassDS-massDspdg>fCutsRD
    1			"pTK [GeV/c]",                              pTK<fCutsRd
    2			"pTPi [GeV/c]",                             pTPi<fCutsRd
    3			"d0K [cm]",                                 d0K<fCutsRd
    4			"d0Pi [cm]",                                d0Pi<fCutsRd
    5			"dist12 [cm]",                              dist12<fCutsRd
    6			"sigmavert [cm]",                           sigmavert>fCutsRd
    7			"decLen [cm]",                              decLen<fCutsRD
    8			"ptMax [GeV/c]",                            ptMax<fCutsRD
    9			"cosThetaPoint",                            CosThetaPoint<fCutsRD
    10			"Sum d0^2 (cm^2)",                          sumd0<fCutsRD
    11			"dca [cm]",                                 dca(i)>fCutsRD
    12			"inv. mass (Mphi-MKK) [GeV]",               invmass-pdg>fCutsRD
    13			"inv. mass (MKo*-MKpi) [GeV]"};             invmass-pdg>fCutsRD
    14    		"Abs(CosineKpiPhiRFrame)^3",
	15  		"CosPiDsLabFrame"};
    */
 
    AliRDHFCutsDstoKKpi *prodcuts = new AliRDHFCutsDstoKKpi();
    prodcuts->SetName("ProdCuts");
    prodcuts->SetPtBins(nptbins+1,ptbins);
    prodcuts->SetCuts(nvars,nptbins,prodcutsval);
    
    
    AliRDHFCutsDstoKKpi* analysiscuts=new AliRDHFCutsDstoKKpi();
    analysiscuts->SetName("AnalysisCuts");
    analysiscuts->SetTitle("Cuts for Ds Analysis and CF");
    analysiscuts->SetPtBins(nptbins+1,ptbins);
    analysiscuts->SetCuts(nvars,nptbins,anacutsval);
    analysiscuts->AddTrackCuts(esdTrackCuts);
    analysiscuts->SetUsePID(kTRUE);
    cout<<"This is the odject I'm going to save:"<<nptbins<<endl;
    
    analysiscuts->PrintAll();
    TFile* fout=new TFile("DstoKKpiCuts.root","recreate");   
    fout->cd();
    prodcuts->Write();
    analysiscuts->Write();
    fout->Close();
    
    
}


void makeInputAliAnalysisTaskSESignificanceMaximization(){
  
  AliRDHFCutsDstoKKpi* RDHFDstoKKpi=new AliRDHFCutsDstoKKpi();
  RDHFDstoKKpi->SetUsePID(kTRUE);
  RDHFDstoKKpi->SetName("loosercuts");
  RDHFDstoKKpi->SetTitle("Cuts for significance maximization");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  RDHFDstoKKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=16;

  const Int_t nptbins=7;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=8.;
  ptbins[6]=12.;
  ptbins[7]=999999999999.;
  
  RDHFDstoKKpi->SetPtBins(nptbins+1,ptbins);
    
  //setting my cut values
    
  
  
  
  Float_t** prodcutsval;
  prodcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}  
  for(Int_t ipt=0;ipt<nptbins;ipt++){
      prodcutsval[0][ipt]=0.2;
      prodcutsval[1][ipt]=0.3;
      prodcutsval[2][ipt]=0.3;
      prodcutsval[3][ipt]=0.;
      prodcutsval[4][ipt]=0.;
      prodcutsval[5][ipt]=0.005;
      prodcutsval[6][ipt]=0.06;
      prodcutsval[7][ipt]=0.0;
      prodcutsval[8][ipt]=0.;
      prodcutsval[9][ipt]=0.9;
      prodcutsval[10][ipt]=0.;
      prodcutsval[11][ipt]=1000.0;
      prodcutsval[12][ipt]=0.03;
      prodcutsval[13][ipt]=0.1;
      prodcutsval[14][ipt]=0.;
      prodcutsval[15][ipt]=1.;
  }

  RDHFDstoKKpi->SetCuts(nvars,nptbins,prodcutsval);

  Int_t nvarsforopt=RDHFDstoKKpi->GetNVarsForOpt();
  //Int_t nvarsforopt=2;
  const Int_t dim=4; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFDstoKKpi->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFDstoKKpi->GetVarNames();
      for(Int_t i=0;i<nvars;i++){
	cout<<names[i]<<" for opt? (y/n)"<<endl;
	cin>>answer;
	if(answer=="y") {
	  boolforopt[i]=kTRUE;
	  checktrue++;
	}
	else boolforopt[i]=kFALSE;
      }
      if (checktrue!=dim) {
	cout<<"Error! You set "<<checktrue<<" kTRUE instead of "<<dim<<endl;
	return;
      }
      RDHFDstoKKpi->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];
 
  //number of steps for each variable is 4 now
  //set this!!
  
  /*    
   Spiegazione della selezione di ALiRDHFCutsDstoKKpi      condizione di rejection

   0           "inv. mass [GeV]",                          invmassDS-massDspdg>fCutsRD
   1			"pTK [GeV/c]",                              pTK<fCutsRd
   2			"pTPi [GeV/c]",                             pTPi<fCutsRd
   3			"d0K [cm]",                                 d0K<fCutsRd
   4			"d0Pi [cm]",                                d0Pi<fCutsRd
   5			"dist12 [cm]",                              dist12<fCutsRd
   6			"sigmavert [cm]",                           sigmavert>fCutsRd
   7			"decLen [cm]",                              decLen<fCutsRD
   8			"ptMax [GeV/c]",                            ptMax<fCutsRD
   9			"cosThetaPoint",                            CosThetaPoint<fCutsRD
   10			"Sum d0^2 (cm^2)",                          sumd0<fCutsRD
   11			"dca [cm]",                                 dca(i)>fCutsRD
   12			"inv. mass (Mphi-MKK) [GeV]",               invmass-pdg>fCutsRD
   13			"inv. mass (MKo*-MKpi) [GeV]"};             invmass-pdg>fCutsRD
   14    		"Abs(CosineKpiPhiRFrame)^3",
   15  	    	"CosPiDsLabFrame"}
   */  
  
  // tighterval[var][ptbin]
  
  //sigmavert [cm]
  
   
  tighterval[0][0]=0.0;
  tighterval[0][1]=0.0;
  tighterval[0][2]=0.0;
  tighterval[0][3]=0.0;
  tighterval[0][4]=0.0;
  
  
   
  tighterval[1][0]=0.05;
  tighterval[1][1]=0.05;
  tighterval[1][2]=0.05;
  tighterval[1][3]=0.05;
  tighterval[1][4]=0.05;
  
 
  //costhetapoint
  
  tighterval[2][0]=1.;
  tighterval[2][1]=1.;
  tighterval[2][2]=1.;
  tighterval[2][3]=1.;
  tighterval[2][4]=1.;
  
  //inv mass phi meson
  
  tighterval[3][0]=0.002;
  tighterval[3][1]=0.002;
  tighterval[3][2]=0.002;
  tighterval[3][3]=0.002;
  tighterval[3][4]=0.002;
  

  TString name=""; 
  Int_t arrdim=dim*nptbins;
  TClonesArray max("TParameter<float>",arrdim);
  for (Int_t ipt=0;ipt<nptbins;ipt++){
    for(Int_t ival=0;ival<dim;ival++){
      name=Form("par%dptbin%d",ival,ipt);
      new(max[ipt*dim+ival])TParameter<float>(name.Data(),tighterval[ival][ipt]);
    }
  }

  TFile* fout=new TFile("cuts4SignifMaximDs.root","recreate");   //set this!! 
  fout->cd();
  RDHFDstoKKpi->Write();
  max.Write();
  fout->Close();

}
