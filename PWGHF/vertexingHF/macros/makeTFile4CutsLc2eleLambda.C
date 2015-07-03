#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsLctoeleLambdafromAODtracks.h"
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>

void makeInputAliAnalysisTaskSELc2eleLambda()
{

  AliRDHFCutsLctoeleLambdafromAODtracks* RDHFLc=new AliRDHFCutsLctoeleLambdafromAODtracks();

  RDHFLc->SetName("eleLambdaAnalysisCuts");
  RDHFLc->SetTitle("Analysis cuts for eleLambda analysis");
  RDHFLc->SetUsePhysicsSelection(kTRUE);
	RDHFLc->SetTriggerClass("");
	//RDHFLc->SetUseInt7TriggerPP2012();
	RDHFLc->ResetMaskAndEnableMBTrigger();
	//RDHFLc->EnableCentralTrigger();
	//RDHFLc->EnableSemiCentralTrigger();

  const Int_t nptbins=1;
  RDHFLc->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]=100.;
  RDHFLc->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=2;

  Float_t** anacutsarray;
	anacutsarray=new Float_t*[nvars];
	for(Int_t ic=0;ic<nvars;ic++)
	{
		anacutsarray[ic]=new Float_t[nptbins];
	}

  for(Int_t ipt2=0;ipt2<nptbins;ipt2++)
	{
	 anacutsarray[ 0][ipt2]=10.;  // Upper limit on Inv Mass
	 anacutsarray[ 1][ipt2]=0.;  // Lower limit on cos(Opening angle)
  }

	RDHFLc->SetUseOnTheFlyV0(kFALSE);

  RDHFLc->SetCuts(nvars,nptbins,anacutsarray);

  RDHFLc->SetUsePID(kTRUE);
	RDHFLc->SetPIDStrategy(AliRDHFCutsLctoeleLambdafromAODtracks::kNSigmaCustomizedCuts);
  AliAODPidHF* pidObj=new AliAODPidHF();
	//used in kNSigma
  Double_t sigmaspi[5]={3.,0.,0.,3.,0.};
  pidObj->SetSigma(sigmaspi);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetMatch(1);
	//used in kNSigmaCombinedCuts
	RDHFLc->SetSigmaElectronTPCRange(-0.5,3.0);
	RDHFLc->SetSigmaElectronTOFRange(-3.0,3.0);
	//RDHFLc->SetSigmaElectronTPCRange(-99999.,99999.);
	//RDHFLc->SetSigmaElectronTOFRange(-99999.,99999.);
	//used both in kNsigma and kNSigmacustomizedcuts
  RDHFLc->SetPidHF(pidObj);
	RDHFLc->SetExcludePionTPC(kFALSE);
	RDHFLc->SetExcludeProtonTPC(kFALSE);
	RDHFLc->SetExcludeKaonTPC(kFALSE);
	RDHFLc->SetExcludenSigmaPionTPC(3.0);
	RDHFLc->SetExcludenSigmaProtonTPC(3.0);
	RDHFLc->SetExcludenSigmaKaonTPC(3.0);

	RDHFLc->SetMaxVtxZ(10.);

	RDHFLc->SetProdUseAODFilterBit(kTRUE);
  RDHFLc->SetProdTrackTPCNclsPIDMin(80);//could not find in esdtrackcutsw
  RDHFLc->SetProdTrackTPCNclsRatioMin(0.6);//could not find in esdtrackcutsw

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(100);
	//esdTrackCuts->SetMinNCrossedRowsTPC(70);
	//esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);//kBoth
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.6,0.6);
  esdTrackCuts->SetPtRange(0.5,1.e10);
	esdTrackCuts->SetMaxDCAToVertexXY(1.0);
	esdTrackCuts->SetMaxDCAToVertexZ(2.0);
	esdTrackCuts->SetDCAToVertex2D(kTRUE);
	esdTrackCuts->SetMaxChi2PerClusterITS(36);
  RDHFLc->AddTrackCuts(esdTrackCuts);

	RDHFLc->SetUseOnTheFlyV0(kFALSE);
	RDHFLc->SetProdV0MassTolLambda(0.006);
	RDHFLc->SetProdV0PtMin(0.0);
	RDHFLc->SetProdV0CosPointingAngleToPrimVtxMin(0.995);//AN 0.995(seems too strong for this analysis)
	RDHFLc->SetProdV0DcaDaughtersMax(1.0);
	RDHFLc->SetProdV0DaughterEtaRange(0.8);
	RDHFLc->SetProdV0DaughterPtMin(0.1);
	RDHFLc->SetProdV0DaughterTPCClusterMin(70);
  RDHFLc->SetProdRfidMinV0(0.5);
  RDHFLc->SetProdRfidMaxV0(9000.);
  RDHFLc->SetProdDcaV0PrToPrimVertexMin(0.06);//default 0.06
  RDHFLc->SetProdDcaV0PiToPrimVertexMin(0.06);//default 0.06
  RDHFLc->SetProdV0ProperDecayLengthMax(30.);
  RDHFLc->SetProdMassRejK0s(0.01);//default 0.01


  RDHFLc->SetUseLambdaPID(kTRUE);
  AliAODPidHF* pidObjpr=new AliAODPidHF();
  Double_t sigmasproton[5]={5.,0.,0.,0.,0.};
  pidObjpr->SetSigma(sigmasproton);
  pidObjpr->SetTPC(kTRUE);
  RDHFLc->SetPidProton(pidObjpr);
  AliAODPidHF* pidObjpi=new AliAODPidHF();
  Double_t sigmaspion[5]={5.,0.,0.,0.,0.};
  pidObjpi->SetSigma(sigmaspion);
  pidObjpi->SetTPC(kTRUE);
  RDHFLc->SetPidPion(pidObjpi);

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFLc->PrintAll();
  TFile* fout=new TFile("eleLambdaCuts.root","RECREATE"); 
  fout->cd();
  RDHFLc->Write();
  fout->Close();
  delete fout;

}
