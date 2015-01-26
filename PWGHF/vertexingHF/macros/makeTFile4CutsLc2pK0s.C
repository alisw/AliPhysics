#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsLctopK0sfromAODtracks.h"
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>

void makeInputAliAnalysisTaskSELc2pK0s()
{

  AliRDHFCutsLctopK0sfromAODtracks* RDHFLc=new AliRDHFCutsLctopK0sfromAODtracks();

  RDHFLc->SetName("LcAnalysisCuts");
  RDHFLc->SetTitle("Analysis cuts for Lc analysis");
  RDHFLc->SetUsePhysicsSelection(kTRUE);
	RDHFLc->SetTriggerClass("");
	//RDHFLc->SetUseInt7TriggerPP2012();
	RDHFLc->ResetMaskAndEnableMBTrigger();
	RDHFLc->EnableCentralTrigger();
	RDHFLc->EnableSemiCentralTrigger();

  const Int_t nptbins=8;
  RDHFLc->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]= 2.;
  ptbins[2]= 3.;
  ptbins[3]= 4.;
  ptbins[4]= 5.;
  ptbins[5]= 6.;
  ptbins[6]= 8.;
  ptbins[7]=12.;
  ptbins[8]=100.;
  RDHFLc->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=7;

  Float_t** anacutsarray;
	anacutsarray=new Float_t*[nvars];
	for(Int_t ic=0;ic<nvars;ic++)
	{
		anacutsarray[ic]=new Float_t[nptbins];
	}

  for(Int_t ipt2=0;ipt2<nptbins;ipt2++)
	{
	 anacutsarray[ 0][ipt2]=0.25;  // Tolerance Invariant mass of Lc
	 anacutsarray[ 1][ipt2]=2.00;  // Min Lc pT
	 anacutsarray[ 2][ipt2]=0.3;  // Min Bach pT
	 anacutsarray[ 3][ipt2]=0.5;  // Max Bach d0
	 anacutsarray[ 4][ipt2]=1.5;  // Max V0 d0
	 anacutsarray[ 5][ipt2]=0.010;  // Tolrelance mass K0s
	 anacutsarray[ 6][ipt2]=0.000;  // Decay Length XY
  }

	RDHFLc->SetUseOnTheFlyV0(kFALSE);

  RDHFLc->SetCuts(nvars,nptbins,anacutsarray);

  RDHFLc->SetUsePID(kTRUE);
	RDHFLc->SetPIDStrategy(AliRDHFCutsLctopK0sfromAODtracks::kNSigmaCuts);
  AliAODPidHF* pidObj=new AliAODPidHF();
  Double_t sigmaspi[5]={3.,0.,0.,3.,0.};
  pidObj->SetSigma(sigmaspi);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetMatch(1);
  RDHFLc->SetPidHF(pidObj);
	RDHFLc->GetPidHF()->SetUseCombined(kTRUE);
	RDHFLc->GetPidHF()->SetUseDefaultPriors(kTRUE);
	RDHFLc->GetPidHF()->SetCombDetectors(AliAODPidHF::kTPCTOF);
	RDHFLc->SetCombinedPIDThreshold(0.02);

	RDHFLc->SetMaxVtxZ(10.);

	RDHFLc->SetUseOnTheFlyV0(kFALSE);
	RDHFLc->SetProdRoughMassTol(0.25);
	RDHFLc->SetProdRoughPtMin(2.0);
	RDHFLc->SetProdTrackPtMin(0.3);
	RDHFLc->SetProdTrackEtaRange(0.8);
	RDHFLc->SetProdUseAODFilterBit(kTRUE);
	RDHFLc->SetProdV0MassTolK0s(0.010);
	RDHFLc->SetProdV0PtMin(0.3);
	RDHFLc->SetProdV0CosPointingAngleToPrimVtxMin(0.99);
	RDHFLc->SetProdV0DcaDaughtersMax(1.5);
	RDHFLc->SetProdV0DaughterEtaRange(0.8);
	RDHFLc->SetProdV0DaughterPtMin(0.1);
	RDHFLc->SetProdV0DaughterTPCClusterMin(70);

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFLc->PrintAll();
  TFile* fout=new TFile("LcCuts.root","RECREATE"); 
  fout->cd();
  RDHFLc->Write();
  fout->Close();
  delete fout;

}
