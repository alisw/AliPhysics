#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsXictoeleXifromAODtracks.h"
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>

void makeInputAliAnalysisTaskSEXic2eleXi()
{

  AliRDHFCutsXictoeleXifromAODtracks* RDHFXic=new AliRDHFCutsXictoeleXifromAODtracks();

  RDHFXic->SetName("eleXiAnalysisCuts");
  RDHFXic->SetTitle("Analysis cuts for eleXi analysis");
  RDHFXic->SetUsePhysicsSelection(kTRUE);
	RDHFXic->SetTriggerClass("");
	//RDHFXic->SetUseInt7TriggerPP2012();
	RDHFXic->ResetMaskAndEnableMBTrigger();
	//RDHFXic->EnableCentralTrigger();
	//RDHFXic->EnableSemiCentralTrigger();

  const Int_t nptbins=2;
  RDHFXic->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]=100.;
  RDHFXic->SetPtBins(nptbins+1,ptbins);

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

  RDHFXic->SetCuts(nvars,nptbins,anacutsarray);

  RDHFXic->SetUsePID(kTRUE);
	RDHFXic->SetPIDStrategy(AliRDHFCutsXictoeleXifromAODtracks::kNSigmaCustomizedCuts);
  AliAODPidHF* pidObj=new AliAODPidHF();
	//used in kNSigma
  Double_t sigmaspi[5]={3.,0.,0.,3.,0.};
  pidObj->SetSigma(sigmaspi);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetMatch(1);
	//used in kNSigmaCombinedCuts
	RDHFXic->SetSigmaElectronTPCRange(-0.5,3.0);
	RDHFXic->SetSigmaElectronTOFRange(-3.0,3.0);
	//RDHFXic->SetSigmaElectronTPCRange(-99999.,99999.);
	////RDHFXic->SetSigmaElectronTOFRange(-99999.,99999.);
	//used both in kNsigma and kNSigmacustomizedcuts
  RDHFXic->SetPidHF(pidObj);
	RDHFXic->SetExcludePionTPC(kFALSE);
	RDHFXic->SetExcludeProtonTPC(kFALSE);
	RDHFXic->SetExcludeKaonTPC(kFALSE);
	RDHFXic->SetExcludenSigmaPionTPC(3.0);
	RDHFXic->SetExcludenSigmaProtonTPC(3.0);
	RDHFXic->SetExcludenSigmaKaonTPC(3.0);
	//RDHFXic->SetBachelorType(0);//0: electron, 2: pion, 3: kaon, 4: proton (used in PID)

	RDHFXic->SetMaxVtxZ(10.);

	RDHFXic->SetProdUseAODFilterBit(kTRUE);
  RDHFXic->SetProdTrackTPCNclsPIDMin(80);//could not find in esdtrackcutsw
  RDHFXic->SetProdTrackTPCNclsRatioMin(0.6);//could not find in esdtrackcutsw

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
  RDHFXic->AddTrackCuts(esdTrackCuts);


  RDHFXic->SetProdMassTolLambda(0.0075);
  RDHFXic->SetProdMassTolXi(0.015);
  RDHFXic->SetProdMassRejOmega(0.000);
  RDHFXic->SetProdRfidMinV0(2.67);
  RDHFXic->SetProdRfidMaxV0(9000.);
  RDHFXic->SetProdRfidMinXi(0.38);
  RDHFXic->SetProdRfidMaxXi(9000.);
  RDHFXic->SetProdCascProperDecayLengthMax(999999.);
  RDHFXic->SetProdDcaXiDaughtersMax(1.68);
  RDHFXic->SetProdDcaV0DaughtersMax(1.18);
  RDHFXic->SetProdDcaBachToPrimVertexMin(0.0204);
  RDHFXic->SetProdDcaV0ToPrimVertexMin(0.03);
  RDHFXic->SetProdDcaV0PrToPrimVertexMin(0.073);
  RDHFXic->SetProdDcaV0PiToPrimVertexMin(0.073);
  RDHFXic->SetProdXiCosineOfPoiningAngleMin(0.983);
  RDHFXic->SetProdV0CosineOfPoiningAngleXiMin(0.9826);
  RDHFXic->SetProdCascNTPCClustersMin(0.);

  RDHFXic->SetUseCascadePID(kTRUE);
  AliAODPidHF* pidObjcascpi=new AliAODPidHF();
  Double_t sigmascascpi[5]={4.,0.,0.,0.,0.};
  pidObjcascpi->SetSigma(sigmascascpi);
  pidObjcascpi->SetTPC(kTRUE);
  RDHFXic->SetPidCascPi(pidObjcascpi);

  AliAODPidHF* pidObjcascpr=new AliAODPidHF();
  Double_t sigmascascpr[5]={4.,0.,0.,0.,0.};
  pidObjcascpr->SetSigma(sigmascascpr);
  pidObjcascpr->SetTPC(kTRUE);
  RDHFXic->SetPidCascPr(pidObjcascpr);

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFXic->PrintAll();
  TFile* fout=new TFile("eleXiCuts.root","RECREATE"); 
  fout->cd();
  RDHFXic->Write();
  fout->Close();
  delete fout;

}
