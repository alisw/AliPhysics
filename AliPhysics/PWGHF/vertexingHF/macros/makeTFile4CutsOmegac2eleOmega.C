#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsOmegactoeleOmegafromAODtracks.h"
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>

void makeInputAliAnalysisTaskSEOmegac2eleOmega()
{

  AliRDHFCutsOmegactoeleOmegafromAODtracks* RDHFOmegac=new AliRDHFCutsOmegactoeleOmegafromAODtracks();

  RDHFOmegac->SetName("eleOmegaAnalysisCuts");
  RDHFOmegac->SetTitle("Analysis cuts for eleOmega analysis");
  RDHFOmegac->SetUsePhysicsSelection(kTRUE);
	RDHFOmegac->SetTriggerClass("");
	//RDHFOmegac->SetUseInt7TriggerPP2012();
	RDHFOmegac->ResetMaskAndEnableMBTrigger();
	//RDHFOmegac->EnableCentralTrigger();
	//RDHFOmegac->EnableSemiCentralTrigger();

  const Int_t nptbins=2;
  RDHFOmegac->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]=100.;
  RDHFOmegac->SetPtBins(nptbins+1,ptbins);

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

  RDHFOmegac->SetCuts(nvars,nptbins,anacutsarray);

  RDHFOmegac->SetUsePID(kTRUE);
	RDHFOmegac->SetPIDStrategy(AliRDHFCutsOmegactoeleOmegafromAODtracks::kNSigmaCustomizedCuts);
  AliAODPidHF* pidObj=new AliAODPidHF();
	//used in kNSigma
  Double_t sigmaspi[5]={3.,0.,0.,3.,0.};
  pidObj->SetSigma(sigmaspi);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetMatch(1);
	//used in kNSigmaCombinedCuts
	RDHFOmegac->SetSigmaElectronTPCRange(-0.5,3.0);
	RDHFOmegac->SetSigmaElectronTOFRange(-3.0,3.0);
	//RDHFOmegac->SetSigmaElectronTPCRange(-99999.,99999.);
	////RDHFOmegac->SetSigmaElectronTOFRange(-99999.,99999.);
	//used both in kNsigma and kNSigmacustomizedcuts
  RDHFOmegac->SetPidHF(pidObj);
	RDHFOmegac->SetExcludePionTPC(kFALSE);
	RDHFOmegac->SetExcludeProtonTPC(kFALSE);
	RDHFOmegac->SetExcludeKaonTPC(kFALSE);
	RDHFOmegac->SetExcludenSigmaPionTPC(3.0);
	RDHFOmegac->SetExcludenSigmaProtonTPC(3.0);
	RDHFOmegac->SetExcludenSigmaKaonTPC(3.0);
	//RDHFOmegac->SetBachelorType(0);//0: electron, 2: pion, 3: kaon, 4: proton (used in PID)

	RDHFOmegac->SetMaxVtxZ(10.);

	RDHFOmegac->SetProdUseAODFilterBit(kTRUE);
  RDHFOmegac->SetProdTrackTPCNclsPIDMin(80);//could not find in esdtrackcutsw
  RDHFOmegac->SetProdTrackTPCNclsRatioMin(0.6);//could not find in esdtrackcutsw

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
  RDHFOmegac->AddTrackCuts(esdTrackCuts);


  RDHFOmegac->SetProdMassTolLambda(0.0075);
  RDHFOmegac->SetProdMassTolOmega(0.015);
  RDHFOmegac->SetProdMassRejXi(0.000);
  RDHFOmegac->SetProdRfidMinV0(2.768);
  RDHFOmegac->SetProdRfidMaxV0(9000.);
  RDHFOmegac->SetProdRfidMinOmega(0.79);
  RDHFOmegac->SetProdRfidMaxOmega(9000.);
  RDHFOmegac->SetProdCascProperDecayLengthMax(999999.);
  RDHFOmegac->SetProdDcaOmegaDaughtersMax(1.033);
  RDHFOmegac->SetProdDcaV0DaughtersMax(1.116);
  RDHFOmegac->SetProdDcaBachToPrimVertexMin(0.049);
  RDHFOmegac->SetProdDcaV0ToPrimVertexMin(0.088);
  RDHFOmegac->SetProdDcaV0PrToPrimVertexMin(0.073);
  RDHFOmegac->SetProdDcaV0PiToPrimVertexMin(0.073);
  RDHFOmegac->SetProdXiCosineOfPoiningAngleMin(0.9821);
  RDHFOmegac->SetProdV0CosineOfPoiningAngleXiMin(0.951);
  RDHFOmegac->SetProdCascNTPCClustersMin(0.);

  RDHFOmegac->SetUseCascadePID(kTRUE);
  AliAODPidHF* pidObjcascpi=new AliAODPidHF();
  Double_t sigmascascpi[5]={4.,0.,0.,0.,0.};
  pidObjcascpi->SetSigma(sigmascascpi);
  pidObjcascpi->SetTPC(kTRUE);
  RDHFOmegac->SetPidCascPi(pidObjcascpi);

  AliAODPidHF* pidObjcascpr=new AliAODPidHF();
  Double_t sigmascascpr[5]={4.,0.,0.,0.,0.};
  pidObjcascpr->SetSigma(sigmascascpr);
  pidObjcascpr->SetTPC(kTRUE);
  RDHFOmegac->SetPidCascPr(pidObjcascpr);

  AliAODPidHF* pidObjcascka=new AliAODPidHF();
  Double_t sigmascascka[5]={4.,0.,0.,0.,0.};
  pidObjcascka->SetSigma(sigmascascka);
  pidObjcascka->SetTPC(kTRUE);
  RDHFOmegac->SetPidCascKa(pidObjcascka);

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFOmegac->PrintAll();
  TFile* fout=new TFile("eleOmegaCuts.root","RECREATE"); 
  fout->cd();
  RDHFOmegac->Write();
  fout->Close();
  delete fout;

}
