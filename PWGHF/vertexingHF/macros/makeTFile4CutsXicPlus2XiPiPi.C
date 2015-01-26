#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsXicPlustoXiPiPifromAODtracks.h"
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>

void makeInputAliAnalysisTaskSEXicPlus2XiPiPi()
{

  AliRDHFCutsXicPlustoXiPiPifromAODtracks* RDHFXicPlus=new AliRDHFCutsXicPlustoXiPiPifromAODtracks();

  RDHFXicPlus->SetName("XicPlusAnalysisCuts");
  RDHFXicPlus->SetTitle("Analysis cuts for Xic analysis");
  RDHFXicPlus->SetUsePhysicsSelection(kTRUE);
	RDHFXicPlus->SetTriggerClass("");
	//RDHFXicPlus->SetUseInt7TriggerPP2012();
	RDHFXicPlus->ResetMaskAndEnableMBTrigger();
	RDHFXicPlus->EnableCentralTrigger();
	RDHFXicPlus->EnableSemiCentralTrigger();

  const Int_t nptbins=5;
  RDHFXicPlus->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]= 5.;
  ptbins[2]= 6.;
  ptbins[3]= 8.;
  ptbins[4]=12.;
  ptbins[5]=100.;
  RDHFXicPlus->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=13;

  Float_t** anacutsarray;
	anacutsarray=new Float_t*[nvars];
	for(Int_t ic=0;ic<nvars;ic++)
	{
		anacutsarray[ic]=new Float_t[nptbins];
	}

  for(Int_t ipt2=0;ipt2<nptbins;ipt2++)
	{
	 anacutsarray[ 0][ipt2]=0.25;  // Tolerance Invariant mass of Xic
	 anacutsarray[ 1][ipt2]=3.00;  // Min pT of Xic
	 anacutsarray[ 2][ipt2]=0.008;  // Tolerance Mass of Xi
	 anacutsarray[ 3][ipt2]=0.010;  // Tolerance Mass of Lambda
	 anacutsarray[ 4][ipt2]=0.2;  // Max Dca between pi-pi
	 anacutsarray[ 5][ipt2]=2.0;  //Max Dca between pi-casc
	 anacutsarray[ 6][ipt2]=0.15;  // Max d0 of pion
	 anacutsarray[ 7][ipt2]=99999.;  // Max d0 of Xi
	 anacutsarray[ 8][ipt2]=0.00;  // Min d0 of Xi bachelor
	 anacutsarray[ 9][ipt2]=0.00;  // Min d0 of Xi Lambda
	 anacutsarray[10][ipt2]=0.8;  // Min Xic Pointing angle
	 anacutsarray[11][ipt2]=0.00;  // Min Xic Decay Length XY
	 anacutsarray[12][ipt2]=0.30;  // Min pT Bachelor
  }

  RDHFXicPlus->SetCuts(nvars,nptbins,anacutsarray);

  RDHFXicPlus->SetUsePID(kFALSE);
	RDHFXicPlus->SetPIDStrategy(AliRDHFCutsXicPlustoXiPiPifromAODtracks::kNSigmaCuts);
  AliAODPidHF* pidObj=new AliAODPidHF();
  Double_t sigmaspi[5]={3.,0.,0.,0.,0.};
  pidObj->SetSigma(sigmaspi);
  pidObj->SetTPC(kTRUE);
  RDHFXicPlus->SetPidHF(pidObj);
	RDHFXicPlus->GetPidHF()->SetUseCombined(kTRUE);
	RDHFXicPlus->GetPidHF()->SetUseDefaultPriors(kTRUE);
	RDHFXicPlus->GetPidHF()->SetCombDetectors(AliAODPidHF::kTPCTOF);
	RDHFXicPlus->SetCombinedPIDThreshold(0.02);

	RDHFXicPlus->SetMaxVtxZ(10.);

	//The following cuts are used before object creation
	RDHFXicPlus->SetProdTrackPtMin(0.2);
	RDHFXicPlus->SetProdTrackEtaRange(0.8);
	RDHFXicPlus->SetProdUseAODFilterBit(kTRUE);
	RDHFXicPlus->SetProdMassTolLambda(0.010);
	RDHFXicPlus->SetProdMassTolXi(0.010);
	RDHFXicPlus->SetProdRfidMinV0(0.6);
	RDHFXicPlus->SetProdRfidMaxV0(100.);
	RDHFXicPlus->SetProdRfidMinXi(0.6);
	RDHFXicPlus->SetProdRfidMaxXi(100.);
	RDHFXicPlus->SetProdLikeSignDcaMax(0.5);
	RDHFXicPlus->SetProdRoughMassTol(0.25);
	RDHFXicPlus->SetProdRoughPtMin(3.0);
  RDHFXicPlus->SetUsePID(kTRUE);

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFXicPlus->PrintAll();
  TFile* fout=new TFile("XicPlusCuts.root","RECREATE"); 
  fout->cd();
  RDHFXicPlus->Write();
  fout->Close();
  delete fout;

}
