//Class to extract data to do ITS+TPC global Spectra
//Autor Marek Chojnacki
//emali Marek.Chojnacki@cern.ch
//Used on 2009 data
//last line of comments


#include "AliAnalysisChargedHadronSpectraITSTruncatedMeanTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"


//#include "AliESDtrack.h"

#include "Riostream.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
//#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TH3F.h"
//#include "TMCProcess.h"
#include "AliVEvent.h"

#include "AliESDtrackCuts.h"
//#include "AliESDpidCuts.h"
//#include "AliESDpid.h"
#include "AliCentrality.h"
#include "AliESDUtils.h"
#include "AliMultiplicity.h"

   class     AliMCEventHandler;
     class   Riostream;

using namespace std;
ClassImp(AliAnalysisChargedHadronSpectraITSTruncatedMeanTask)

//________________________________________________________________________
AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::AliAnalysisChargedHadronSpectraITSTruncatedMeanTask(const char *name) 
:AliAnalysisTaskSE(name),fESD(0),fCuts(0),fCutsMul(0),fMC(0),
fLowMultiplicity(-1),fUpMultiplicity(-1),fLowCentrality(-10.0),fUpCentrality(-10.0),fSPD(0),fUsePilerejection(0),
fYCut(100.0),fsigmacut(3.0),fnsigmaxy(7.0),fnsigmaz(5.0),fchargeCut(0.0),
fCorrectSDD(0),fCorrectSSD(0),fHIsettings(0),fdovertexrescuts(0),
 fK0weight(0),flambdaweight(0),fAntilambdaweight(0),
fHistStats(0),fHistZVertexBeforeCut(0),fHistZVertexAfterCut(0),fHistXYVertexBeforeCut(0),fHistXYVertexAfterCut(0),
fHistPhiPtBeforeCuts(0),fHistPhiPtAfterCuts(0),fHistEtaPtBeforeCuts(0),fHistEtaPtAfterCuts(0),fHistDCABeforeCuts(0),fHistDCAAfterCuts(0),
fHistPminusTPCinPAfterCuts(0),fHistPminusTPCinPglobalAfterCuts(0),
fHistMydEPpositive(0),fHistMydETPCinPpositive(0),fHistMydETPCinPglobalpositive(0),
fHistMydEPnegative(0),fHistMydETPCinPnegative(0),fHistMydETPCinPglobalnegative(0),
fHistL3dEP(0),fHistL4dEP(0),fHistL5dEP(0),fHistL6dEP(0),fHistL3dETPCinP(0),
 fHistL4dETPCinP(0),fHistL5dETPCinP(0),fHistL6dETPCinP(0),fHistwhichhasmin(0),fHistMysignalminusESD(0),
fHistminsignalifPionP(0),fHistminsignalifKaonP(0),fHistminsignalifProtonP(0),fHistminsignalifAntiPionP(0),fHistminsignalifAntiKaonP(0),fHistminsignalifAntiProtonP(0),
fDCAXYZforcleanPions(0),fDCAXYZforcleanAntiPions(0),fDCAXYZforcleanProtons(0),fDCAXYZforcleanAntiProtons(0),
fDCAXYZOpenforcleanPions(0),fDCAXYZOpenforcleanAntiPions(0),fDCAXYZOpenforcleanProtons(0),fDCAXYZOpenforcleanAntiProtons(0),
fHistNtrackwithstandardcuts(0),fHistNtrackwithITSPIDcuts(0),
fHistSignalinTPCKaonforstandardcuts(0),fHistSignalinTPCKaonforITSPIDcuts(0),fHistSignalinTPCAntiKaonforstandardcuts(0),fHistSignalinTPCAntiKaonforITSPIDcuts(0),
fHistSignalinTPCProtonforstandardcuts(0),fHistSignalinTPCProtonforITSPIDcuts(0),fHistSignalinTPCAntiProtonforstandardcuts(0),fHistSignalinTPCAntiProtonforITSPIDcuts(0),
fHistStandartMul(0),fHistMytrackMul(0),fHistStandartMulvSPD2(0),
fHistminsignalifPionPPrimary(0),fHistminsignalifKaonPPrimary(0),fHistminsignalifProtonPPrimary(0),fHistminsignalifProtonPPrimaryfake(0),
fHistminsignalifAntiPionPPrimary(0),fHistminsignalifAntiKaonPPrimary(0),fHistminsignalifAntiProtonPPrimary(0),fHistminsignalifAntiProtonPPrimaryfake(0),
fHistminsignalifPionPSecondary(0),fHistminsignalifKaonPSecondary(0),
fHistminsignalifProtonPSecondaryWD(0),fHistminsignalifProtonPSecondaryHI(0),fHistminsignalifProtonPSecondaryRest(0),
fHistminsignalifProtonPSecondaryWDfake(0),fHistminsignalifProtonPSecondaryHIfake(0),
fHistminsignalifAntiPionPSecondary(0),fHistminsignalifAntiKaonPSecondary(0),
fHistminsignalifAntiProtonPSecondaryWD(0),fHistminsignalifAntiProtonPSecondaryHI(0), fHistminsignalifAntiProtonPSecondaryRest(0),
fHistminsignalifAntiProtonPSecondaryWDfake(0),fHistminsignalifAntiProtonPSecondaryHIfake(0),
fHistminsignalifMuEPositiveP(0),fHistminsignalifMuENegativeP(0),
fHistminsignalifPionPrimaryfake(0),fHistminsignalifKaonPrimaryfake(0),fHistminsignalifAntiPionPrimaryfake(0),fHistminsignalifAntiKaonPrimaryfake(0),
fHistminsignalifPionSecondaryfake(0),fHistminsignalifKaonSecondaryfake(0),fHistminsignalifAntiPionSecondaryfake(0),fHistminsignalifAntiKaonSecondaryfake(0),
fHistminsignalifPionPMCPrimary(0),fHistminsignalifKaonPMCPrimary(0),fHistminsignalifProtonPMCPrimary(0),
fHistminsignalifAntiPionPMCPrimary(0),fHistminsignalifAntiKaonPMCPrimary(0),fHistminsignalifAntiProtonPMCPrimary(0),
fHistminsignalifPionPMCPrimaryBeforeEventCuts(0),fHistminsignalifKaonPMCPrimaryBeforeEventCuts(0),fHistminsignalifProtonPMCPrimaryBeforeEventCuts(0),
fHistminsignalifAntiPionPMCPrimaryBeforeEventCuts(0),fHistminsignalifAntiKaonPMCPrimaryBeforeEventCuts(0),fHistminsignalifAntiProtonPMCPrimaryBeforeEventCuts(0),
fHistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex(0),fHistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex(0),fHistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex(0),
fHistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex(0),fHistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex(0),fHistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex(0),
fHistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ(0),fHistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ(0),fHistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ(0),
fHistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ(0),fHistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ(0),fHistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ(0),
fDCAXYZforcleanPionsMCPrimary(0),fDCAXYZforcleanAntiPionsMCPrimary(0),fDCAXYZforcleanProtonsMCPrimary(0),fDCAXYZforcleanAntiProtonsMCPrimary(0),
 fDCAXYZforcleanPionsWD(0),fDCAXYZforcleanAntiPionsWD(0),fDCAXYZforcleanProtonsWD(0), fDCAXYZforcleanAntiProtonsWD(0),fDCAXYZforcleanPionsHI(0),fDCAXYZforcleanAntiPionsHI(0),
fDCAXYZforcleanProtonsHI(0),fDCAXYZforcleanAntiProtonsHI(0),fDCAXYZforcleanPionsMEPrimary(0),fDCAXYZforcleanAntiPionsMEPrimary(0),fDCAXYZforcleanPionsMESecondary(0),fDCAXYZforcleanAntiPionsMESecondary(0),fDCAXYZforcleanPionsR(0),fDCAXYZforcleanAntiPionsR(0),fDCAXYZforcleanProtonsR(0),fDCAXYZforcleanAntiProtonsR(0),
fDCAXYZOpenforcleanPionsMCPrimary(0),fDCAXYZOpenforcleanAntiPionsMCPrimary(0),fDCAXYZOpenforcleanProtonsMCPrimary(0),fDCAXYZOpenforcleanAntiProtonsMCPrimary(0),
 fDCAXYZOpenforcleanPionsWD(0),fDCAXYZOpenforcleanAntiPionsWD(0),fDCAXYZOpenforcleanProtonsWD(0), fDCAXYZOpenforcleanAntiProtonsWD(0),fDCAXYZOpenforcleanPionsHI(0),fDCAXYZOpenforcleanAntiPionsHI(0),
fDCAXYZOpenforcleanProtonsHI(0),fDCAXYZOpenforcleanAntiProtonsHI(0),fDCAXYZOpenforcleanPionsMEPrimary(0),fDCAXYZOpenforcleanAntiPionsMEPrimary(0),fDCAXYZOpenforcleanPionsMESecondary(0),fDCAXYZOpenforcleanAntiPionsMESecondary(0),fDCAXYZOpenforcleanPionsR(0),fDCAXYZOpenforcleanAntiPionsR(0),fDCAXYZOpenforcleanProtonsR(0),fDCAXYZOpenforcleanAntiProtonsR(0),
fElectronsource(0),fAntiElectronsource(0), 
fMuonsource(0),fAntiMuonsource(0),
fPionNTPCClusters(0),fAntiPionNTPCClusters(0),fKaonNTPCClusters(0),fAntiKaonNTPCClusters(0),fProtonNTPCClusters(0),fAntiProtonNTPCClusters(0),
fPionchi2(0),fAntiPionchi2(0),fKaonchi2(0),fAntiKaonchi2(0),fProtonchi2(0),fAntiProtonchi2(0),
fTracksCutmonitoring(0),fParticlesCutmonitoring(0),fVertexshift(0),fPtESDminusPtMCvPtESDafterallcuts(0),fPtESDminusPtMCvPtESDafterTPCcuts(0),fMulESDMulMCVz(0),
fTPCPIDCUT(0), fESDpid(0),fPrimaryElectronsMother(0),
flist(0)
{
	//Constructor
	 fESDpid=new AliESDpid();
	// fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
	 fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28949/50., 2.74095e+01, TMath::Exp(-3.21763e+01), 2.44026, 6.58800); 
	 
	
	 
	fCutsMul=new AliESDtrackCuts("Mul","Mul");
	fCutsMul->SetMinNClustersTPC(70);
	fCutsMul->SetMaxChi2PerClusterTPC(4);
	fCutsMul->SetAcceptKinkDaughters(kFALSE);
	fCutsMul->SetRequireTPCRefit(kTRUE);
	// ITS
	fCutsMul->SetRequireITSRefit(kTRUE);
	fCutsMul->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
						AliESDtrackCuts::kAny);
	fCutsMul->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
	
	fCutsMul->SetMaxDCAToVertexZ(2);
	fCutsMul->SetDCAToVertex2D(kFALSE);
	fCutsMul->SetRequireSigmaToVertex(kFALSE);
	
	fCutsMul->SetEtaRange(-0.8,+0.8);
	fCutsMul->SetPtRange(0.15, 1e10);
	 
	fdcaxypar[0]=0.0050;
	fdcaxypar[1]=0.0060;
	fdcaxypar[2]=0.9;
	
	fdcazpar[0]=0.0146;
	fdcazpar[1]=0.0070;
	fdcazpar[2]=1.114758;
	fdcazpar[3]=0.0216;

	
	
	
 	Printf("end of AliAnalysisChargedHadronSpectraITSTruncatedMeanTask");
	 DefineOutput(1, TList::Class());
}


//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::UserCreateOutputObjects()
{
	//UserCreateOutputObject
	Printf("AliAnalysisChargedHadronSpectraITSTruncatedMeanTask UserCreateOutputObjects");
	  flist=new TList();
	 flist->SetOwner();
	Float_t ptmax=2.0;
	Float_t etamax=1.0;
	
	Int_t netabins=100;
	Int_t nptbins=40;
	Double_t dcamax=3.7;
	const  Int_t ndec=2;
	Int_t startvalue=-1;
	const  Int_t npredec=50;
	Double_t tabx[ndec*npredec+1];
	for (Int_t i=0;i<ndec;i++)
	{
		for (Int_t j=0;j<npredec;j++)
		{
			tabx[npredec*i+j]=TMath::Power(10,((Double_t)i)+((Double_t)startvalue)+((Double_t)j)/((Double_t)npredec));
		}	
	}
  	tabx[ndec*npredec]=TMath::Power(10,ndec+startvalue);
	
	const Int_t  ny=600;
	const Double_t jump=1.5;
	const Double_t starty=0.0; 
	
	
	Int_t kPtBins=30;
	Double_t binsPtDummy[kPtBins+1];
	binsPtDummy[0]=0.0;
	for(int i=1;i<=kPtBins+1;i++)
	{
		if(binsPtDummy[i-1]+0.05<1.01)
			binsPtDummy[i]=binsPtDummy[i-1]+0.05;
		else
			binsPtDummy[i]=binsPtDummy[i-1]+0.1;	
	} 
	//{0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
	
	
     	fHistStats=new TH1F("StatsHist","StatsHist",10,-0.5,9.5);
	fHistStats->GetXaxis()->SetBinLabel(1,"Nevents");
	fHistStats->GetXaxis()->SetBinLabel(2,"Nevents_physics");
	fHistStats->GetXaxis()->SetBinLabel(3,"Nevents_physics_with_good_SPD_vertex");
	fHistStats->GetXaxis()->SetBinLabel(4,"Nevents_physics_with_good_SPD_vertex_and_|z|<10.0");
	fHistStats->GetXaxis()->SetBinLabel(5,"N_tracks_with_3_or_4_in_SSD_SDD");
	fHistStats->GetXaxis()->SetBinLabel(6,"N_tracks_with_3_or_4_with_goodq_in_SSD_SDD");
	fHistStats->GetXaxis()->SetBinLabel(7,"e_in_pi");
	fHistStats->GetXaxis()->SetBinLabel(8,"mu_in_pi");
	fHistStats->GetXaxis()->SetBinLabel(9,"MC_event");
	fHistStats->GetXaxis()->SetBinLabel(10,"MC_event_with_z<10.0");
	flist->Add(fHistStats);
	
	fHistZVertexBeforeCut=new TH1F("HistZVertexBeforeCut","ZVertex;z[cm];N_{counts}",400,-20,20);
	flist->Add(fHistZVertexBeforeCut);
	fHistZVertexAfterCut=new TH1F("HistZVertexAfterCut","ZVertex;z[cm];N_{counts}",400,-20,20);
	flist->Add(fHistZVertexAfterCut);
	fHistXYVertexBeforeCut=new TH2F("HistXYVertexBeforeCut","XYVertex;x[cm];y[cm];N_{conuts}",100,-0.4,0.4,100,-0.4,0.4); 
	flist->Add(fHistXYVertexBeforeCut);
	fHistXYVertexAfterCut=new TH2F("HistXYVertexAfterCut","XYVertex;x[cm];y[cm];N_{conuts}",100,-0.4,0.4,100,-0.4,0.4); 
	flist->Add(fHistXYVertexAfterCut);
	
	fHistPhiPtBeforeCuts=new  TH2F("HistPhiPtBeforeCuts",";#phi;pt[GeV/c]",70,0,2.0*TMath::Pi(),nptbins,-1.0*ptmax,ptmax);
	flist->Add(fHistPhiPtBeforeCuts);
	
	fHistPhiPtAfterCuts=new  TH2F("HistPhiPtAfterCuts",";#phi;pt[GeV/c]",70,0,2.0*TMath::Pi(),nptbins,-1.0*ptmax,ptmax);
	flist->Add(fHistPhiPtAfterCuts);
	
	fHistEtaPtBeforeCuts=new  TH2F("HistEtaPtBeforeCuts",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtBeforeCuts);
	
	fHistEtaPtAfterCuts=new  TH2F("HistEtaPtAfterCuts",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtAfterCuts);
	
	fHistDCABeforeCuts=new  TH2F("HistDCABeforeCut",";dcaxy[cm];dcaz[cm]",100,-1.0*dcamax,dcamax,100,-1.0*dcamax,1.0*dcamax);
	flist->Add(fHistDCABeforeCuts);
	
	fHistDCAAfterCuts=new  TH2F("HistDCAAfterCut",";dcaxy[cm];dcaz[cm]",100,-1.0*dcamax,dcamax,100,-1.0*dcamax,1.0*dcamax);
	flist->Add(fHistDCAAfterCuts);
	
	
	fHistPminusTPCinPAfterCuts= new TH2F("HistPminusTPCinPVPTPCinAfterCuts",";P-PTPCin [GeV/c];PTPCin",100,-0.5,0.5,ndec*npredec,tabx); 
	flist->Add(fHistPminusTPCinPAfterCuts);
	
	fHistPminusTPCinPglobalAfterCuts= new TH2F("HistPminusTPCinPVPTPCinglobalAfterCuts",";P-PTPCinglobal [GeV/c];PTPCin",100,-0.5,0.5,ndec*npredec,tabx); 
	flist->Add(fHistPminusTPCinPglobalAfterCuts);
	
	fHistMydEPpositive=new TH2F("HistMydEPpositive",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydEPpositive);
		
	fHistMydETPCinPpositive=new TH2F("HistMydETPCinPpositive",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydETPCinPpositive);

	fHistMydETPCinPglobalpositive=new TH2F("HistMydETPCinPglobalpositive",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydETPCinPglobalpositive);
	
	fHistMydEPnegative=new TH2F("HistMydEPnegative",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydEPnegative);
		
	fHistMydETPCinPnegative=new TH2F("HistMydETPCinPnegative",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydETPCinPnegative);

	fHistMydETPCinPglobalnegative=new TH2F("HistMydETPCinPglobalnegative",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistMydETPCinPglobalnegative);
	
	
	fHistL3dEP=new TH2F("HistL3dEP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL3dEP);
	
	fHistL4dEP=new TH2F("HistL4dEP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL4dEP);
	
	fHistL5dEP=new TH2F("HistL5dEP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL5dEP);
	
	fHistL6dEP=new TH2F("HistL6dEP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL6dEP);

	fHistL3dETPCinP=new TH2F("HistL3dETPCinP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL3dETPCinP);
	
	fHistL4dETPCinP=new TH2F("HistL4dETPCinP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL4dETPCinP);
	
	fHistL5dETPCinP=new TH2F("HistL5dETPCinP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL5dETPCinP);
	
	fHistL6dETPCinP=new TH2F("HistL6dETPCinP",";P[GeV/c];dE[in 300#mum]",ndec*npredec,tabx,ny,starty,ny*jump+starty);
	flist->Add(fHistL6dETPCinP);
	


	fHistwhichhasmin=new TH2F("Histwhichhasmin","Histwhichhasmin;L;Q",4,-0.5,3.5,100,0,1000);
	fHistwhichhasmin->GetXaxis()->SetBinLabel(1,"SDD1");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(2,"SDD2");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(3,"SSD1");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(4,"SSD2");
	flist->Add(fHistwhichhasmin);
	fHistMysignalminusESD=new TH1F("HistMysignalminus","HistMysignalminus;my-ESD;N",100,-0.2,0.2);
	flist->Add(fHistMysignalminusESD);
	
	 
	fHistminsignalifPionP=new TH2F("HistminsignalifPionP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifPionP);
	fHistminsignalifKaonP=new TH2F("HistminsignalifKaonP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifKaonP);
	fHistminsignalifProtonP=new TH2F("HistminsignalifProtonP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonP);
	
	
	fHistminsignalifAntiPionP=new TH2F("HistminsignalifAntiPionP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiPionP);
	fHistminsignalifAntiKaonP=new TH2F("HistminsignalifAntiKaonP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiKaonP);
	fHistminsignalifAntiProtonP=new TH2F("HistminsignalifAntiProtonP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonP);
	
	Int_t kDCABins=20+2+40;
	Double_t binsDCADummy[62+1]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,-0.25,-0.2,-0.19,-0.18,-0.17,-0.16,-0.15,-0.14,-0.13,-0.12,-0.11,-0.10,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};
	

	
	
	fDCAXYZforcleanPions=new TH3F("fDCAXYZforcleanPions",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPions);
	fDCAXYZforcleanAntiPions=new TH3F("fDCAXYZforcleanAntiPions",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPions);
	fDCAXYZforcleanProtons=new TH3F("fDCAXYZforcleanProtons",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanProtons);
	fDCAXYZforcleanAntiProtons=new TH3F("fDCAXYZforcleanAntiProtons",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiProtons);
	
	
	fDCAXYZOpenforcleanPions=new TH3F("fDCAXYZOpenforcleanPions",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPions);
	fDCAXYZOpenforcleanAntiPions=new TH3F("fDCAXYZOpenforcleanAntiPions",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPions);
	fDCAXYZOpenforcleanProtons=new TH3F("fDCAXYZOpenforcleanProtons",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanProtons);
	fDCAXYZOpenforcleanAntiProtons=new TH3F("fDCAXYZOpenforcleanAntiProtons",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiProtons);
	
	fHistNtrackwithstandardcuts=new TH2F ("fHistNtrackwithstandardcuts",";Pt[GeV/c];type;counts",kPtBins,binsPtDummy,3,-0.5,2.5);
	flist->Add(fHistNtrackwithstandardcuts);
	fHistNtrackwithITSPIDcuts=new TH2F ("fHistNtrackwithITSPIDcuts",";Pt[GeV/c];type;counts",kPtBins,binsPtDummy,3,-0.5,2.5);
	flist->Add(fHistNtrackwithITSPIDcuts);
	
	
	fHistSignalinTPCKaonforstandardcuts= new TH2F("fHistSignalinTPCKaonforstandardcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCKaonforstandardcuts);
	fHistSignalinTPCKaonforITSPIDcuts= new TH2F("fHistSignalinTPCKaonforITSPIDcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCKaonforITSPIDcuts);
	
	fHistSignalinTPCAntiKaonforstandardcuts= new TH2F("fHistSignalinTPCAntiKaonforstandardcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCAntiKaonforstandardcuts);
	fHistSignalinTPCAntiKaonforITSPIDcuts= new TH2F("fHistSignalinTPCAntiKaonforITSPIDcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCAntiKaonforITSPIDcuts);
	
	fHistSignalinTPCProtonforstandardcuts= new TH2F("fHistSignalinTPCProtonforstandardcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCProtonforstandardcuts);
	fHistSignalinTPCProtonforITSPIDcuts= new TH2F("fHistSignalinTPCProtonforITSPIDcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCProtonforITSPIDcuts);
	
	
	fHistSignalinTPCAntiProtonforstandardcuts= new TH2F("fHistSignalinTPCAntiProtonforstandardcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCAntiProtonforstandardcuts);
	fHistSignalinTPCAntiProtonforITSPIDcuts= new TH2F("fHistSignalinTPCAntiProtonforITSPIDcuts",";Pt[GeV/c];signal",kPtBins,binsPtDummy,100,-1.0,1.0) ;
	flist->Add(fHistSignalinTPCAntiProtonforITSPIDcuts);
	
	fPionNTPCClusters=new TH2F("fPionNTPCClusters","fPionNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160);
	flist->Add(fPionNTPCClusters);
	fAntiPionNTPCClusters=new TH2F("fAntiPionNTPCClusters","fAntiPionNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160);
	flist->Add(fAntiPionNTPCClusters);
	fKaonNTPCClusters=new TH2F("fKaonNTPCClusters","fKaonNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160); 
	flist->Add(fKaonNTPCClusters);
	fAntiKaonNTPCClusters=new TH2F("fAntiKaonNTPCClusters","fAntiKaonNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160);
	flist->Add(fAntiKaonNTPCClusters);
	fProtonNTPCClusters=new TH2F("fProtonNTPCClusters","fProtonNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160); 
	flist->Add(fProtonNTPCClusters);
	fAntiProtonNTPCClusters=new TH2F("fAntiProtonNTPCClusters","fAntiProtonNTPCClusters;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,60,160);
	flist->Add(fAntiProtonNTPCClusters);
	
	
	fPionchi2=new TH2F("fPionchi2","fPionchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6);
	flist->Add(fPionchi2);
	fAntiPionchi2=new TH2F("fAntiPionchi2","fAntiPionchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6);
	flist->Add(fAntiPionchi2);
	fKaonchi2=new TH2F("fKaonchi2","fKaonchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6); 
	flist->Add(fKaonchi2);
	fAntiKaonchi2=new TH2F("fAntiKaonchi2","fAntiKaonchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6);
	flist->Add(fAntiKaonchi2);
	fProtonchi2=new TH2F("fProtonchi2","fProtonchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6); 
	flist->Add(fProtonchi2);
	fAntiProtonchi2=new TH2F("fAntiProtonchi2","fAntiProtonchi2;Pt [GeV/c];NClusters",kPtBins,binsPtDummy,20,0,6);
	flist->Add(fAntiProtonchi2);
	
	if(fHIsettings)
	{
		fHistStandartMul=new TH1F("fHistStandartMul",";Ntracks;counts",300,0,3000);
		flist->Add(fHistStandartMul);
		fHistMytrackMul=new TH1F("fHistMytrackMul",";Ntracks;counts",300,0,3000);
		flist->Add(fHistMytrackMul);
		fHistStandartMulvSPD2=new TH2F("fHistStandartMulvSPD2",";Ntracks;nSPD2;counts",300,0,3000,300,0,3000);
		flist->Add(fHistStandartMulvSPD2);
	}
	else
	{
		fHistStandartMul=new TH1F("fHistStandartMul",";Ntracks;counts",300,0,300);
		flist->Add(fHistStandartMul);
		fHistMytrackMul=new TH1F("fHistMytrackMul",";Ntracks;counts",300,0,300);
		flist->Add(fHistMytrackMul);
		fHistStandartMulvSPD2=new TH2F("fHistStandartMulvSPD2",";Ntracks;nSPD2;counts",300,0,300,300,0,300);
		flist->Add(fHistStandartMulvSPD2);
	}
	fTracksCutmonitoring=new TH2F("fTracksCutmonitoring",";cut;pt[GeV/c];N_{entries}",4,0.5,4.5,kPtBins,binsPtDummy);	
	fTracksCutmonitoring->GetXaxis()->SetBinLabel(1,"TPCin");
	fTracksCutmonitoring->GetXaxis()->SetBinLabel(2,"standard");
	fTracksCutmonitoring->GetXaxis()->SetBinLabel(3,"ITSpid");
	fTracksCutmonitoring->GetXaxis()->SetBinLabel(4,"DCA");
	flist->Add(fTracksCutmonitoring);
	
	
	
	
	if(!fMC)
	{
		
		Printf("end of CreateOutputObjects no MC");
		PostData(1,flist);
		return;
	}
	
	

	fHistminsignalifPionPPrimary=new TH2F("HistminsignalifPionPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifPionPPrimary);
	fHistminsignalifKaonPPrimary=new TH2F("HistminsignalifKaonPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifKaonPPrimary);
	fHistminsignalifProtonPPrimary=new TH2F("HistminsignalifProtonPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPPrimary);
	fHistminsignalifProtonPPrimaryfake=new TH2F("HistminsignalifProtonPPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPPrimaryfake);
		
	fHistminsignalifAntiPionPPrimary=new TH2F("HistminsignalifAntiPionPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiPionPPrimary);
	fHistminsignalifAntiKaonPPrimary=new TH2F("HistminsignalifAntiKaonPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiKaonPPrimary);
	fHistminsignalifAntiProtonPPrimary=new TH2F("HistminsignalifAntiProtonPPrimary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPPrimary);
	fHistminsignalifAntiProtonPPrimaryfake=new TH2F("HistminsignalifAntiProtonPPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPPrimaryfake);
	
	fHistminsignalifPionPSecondary=new TH2F("HistminsignalifPionPSecondary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifPionPSecondary);
	fHistminsignalifKaonPSecondary=new TH2F("HistminsignalifKaonPSecondary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifKaonPSecondary);
	fHistminsignalifProtonPSecondaryWD=new TH2F("HistminsignalifProtonPSecondaryWD",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPSecondaryWD);
	fHistminsignalifProtonPSecondaryHI=new TH2F("HistminsignalifProtonPSecondaryHI",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPSecondaryHI);
	fHistminsignalifProtonPSecondaryRest=new TH2F("HistminsignalifProtonPSecondaryRest",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPSecondaryRest);
	fHistminsignalifProtonPSecondaryWDfake=new TH2F("HistminsignalifProtonPSecondaryWDfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPSecondaryWDfake);
	fHistminsignalifProtonPSecondaryHIfake=new TH2F("HistminsignalifProtonPSecondaryHIfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifProtonPSecondaryHIfake);

		
	fHistminsignalifAntiPionPSecondary=new TH2F("HistminsignalifAntiPionPSecondary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiPionPSecondary);
	fHistminsignalifAntiKaonPSecondary=new TH2F("HistminsignalifAntiKaonPSecondary",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiKaonPSecondary);
	fHistminsignalifAntiProtonPSecondaryWD=new TH2F("HistminsignalifAntiProtonPSecondaryWD",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPSecondaryWD);
	fHistminsignalifAntiProtonPSecondaryHI=new TH2F("HistminsignalifAntiProtonPSecondaryHI",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPSecondaryHI);
	fHistminsignalifAntiProtonPSecondaryRest=new TH2F("HistminsignalifAntiProtonPSecondaryRest",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPSecondaryRest);	
	fHistminsignalifAntiProtonPSecondaryWDfake=new TH2F("HistminsignalifAntiProtonPSecondaryWDfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPSecondaryWDfake);
	fHistminsignalifAntiProtonPSecondaryHIfake=new TH2F("HistminsignalifAntiProtonPSecondaryHIfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiProtonPSecondaryHIfake);
	
	fHistminsignalifMuEPositiveP=new TH2F("HistminsignalifMuEPositiveP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifMuEPositiveP);
	fHistminsignalifMuENegativeP=new TH2F("HistminsignalifMuENegativeP",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifMuENegativeP);
	
	
	fHistminsignalifPionPrimaryfake=new TH2F("HistminsignalifPionPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifPionPrimaryfake);
  	fHistminsignalifKaonPrimaryfake=new TH2F("HistminsignalifKaonPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
 	flist->Add(fHistminsignalifKaonPrimaryfake);
	
	fHistminsignalifAntiPionPrimaryfake=new TH2F("HistminsignalifAntiPionPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiPionPrimaryfake);
  	fHistminsignalifAntiKaonPrimaryfake=new TH2F("HistminsignalifAntiKaonPrimaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
 	flist->Add(fHistminsignalifAntiKaonPrimaryfake);
	
	
	fHistminsignalifPionSecondaryfake=new TH2F("HistminsignalifPionSecondaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifPionSecondaryfake);
  	fHistminsignalifKaonSecondaryfake=new TH2F("HistminsignalifKaonSecondaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
 	flist->Add(fHistminsignalifKaonSecondaryfake);
	
	fHistminsignalifAntiPionSecondaryfake=new TH2F("HistminsignalifAntiPionSecondaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
	flist->Add(fHistminsignalifAntiPionSecondaryfake);
  	fHistminsignalifAntiKaonSecondaryfake=new TH2F("HistminsignalifAntiKaonSecondaryfake",";Pt[GeV/c];log(dE_real)-log(dE_fit)",kPtBins,binsPtDummy,ny,-4,4);
 	flist->Add(fHistminsignalifAntiKaonSecondaryfake);
	
	
	fHistminsignalifPionPMCPrimary=new TH1F("HistminsignalifPionPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifPionPMCPrimary);
	fHistminsignalifKaonPMCPrimary=new TH1F("HistminsignalifKaonPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifKaonPMCPrimary);
	fHistminsignalifProtonPMCPrimary=new TH1F("HistminsignalifProtonPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifProtonPMCPrimary);
		
	fHistminsignalifAntiPionPMCPrimary=new TH1F("HistminsignalifAntiPionPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiPionPMCPrimary);
	fHistminsignalifAntiKaonPMCPrimary=new TH1F("HistminsignalifAntiKaonPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiKaonPMCPrimary);
	fHistminsignalifAntiProtonPMCPrimary=new TH1F("HistminsignalifAntiProtonPMCPrimary",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiProtonPMCPrimary);
	
	
	fHistminsignalifPionPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifPionPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifPionPMCPrimaryBeforeEventCuts);
	fHistminsignalifKaonPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifKaonPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifKaonPMCPrimaryBeforeEventCuts);
	fHistminsignalifProtonPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifProtonPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifProtonPMCPrimaryBeforeEventCuts);
		
	fHistminsignalifAntiPionPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifAntiPionPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiPionPMCPrimaryBeforeEventCuts);
	fHistminsignalifAntiKaonPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifAntiKaonPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiKaonPMCPrimaryBeforeEventCuts);
	fHistminsignalifAntiProtonPMCPrimaryBeforeEventCuts=new TH1F("HistminsignalifAntiProtonPMCPrimaryBeforeEventCuts",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiProtonPMCPrimaryBeforeEventCuts);
	
	fHistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex);
	fHistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex);
	fHistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex);
		
	fHistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex);
	fHistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex);
	fHistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex=new TH1F("HistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex);
	
	fHistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ);
	fHistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ);
	fHistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ);
		
	fHistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ);
	fHistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ);
	fHistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ=new TH1F("HistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ",";Pt[GeV/c]",kPtBins,binsPtDummy);
	flist->Add(fHistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ);
	
	
	
	fDCAXYZforcleanPionsMCPrimary=new TH3F("fDCAXYZforcleanPionsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsMCPrimary);
	fDCAXYZforcleanAntiPionsMCPrimary=new TH3F("fDCAXYZforcleanAntiPionsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsMCPrimary);
	fDCAXYZforcleanProtonsMCPrimary=new TH3F("fDCAXYZforcleanProtonsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanProtonsMCPrimary);
	fDCAXYZforcleanAntiProtonsMCPrimary=new TH3F("fDCAXYZforcleanAntiProtonsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiProtonsMCPrimary);
	
	fDCAXYZforcleanPionsWD=new TH3F("fDCAXYZforcleanPionsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsWD);
	fDCAXYZforcleanAntiPionsWD=new TH3F("fDCAXYZforcleanAntiPionsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsWD);
	//Secondrary Protons weak deacy
	
	fDCAXYZforcleanProtonsWD=new TH3F("fDCAXYZforcleanProtonsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanProtonsWD);
	fDCAXYZforcleanAntiProtonsWD=new TH3F("fDCAXYZforcleanAntiProtonsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiProtonsWD);

	fDCAXYZforcleanPionsHI=new TH3F("fDCAXYZforcleanPionsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsHI);
	fDCAXYZforcleanAntiPionsHI=new TH3F("fDCAXYZforcleanAntiPionsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsHI);
	//Secondrary Protons Hadronic
	fDCAXYZforcleanProtonsHI=new TH3F("fDCAXYZforcleanProtonsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanProtonsHI);
	fDCAXYZforcleanAntiProtonsHI=new TH3F("fDCAXYZforcleanAntiProtonsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiProtonsHI);
	//Secondrary Pions mu el
	fDCAXYZforcleanPionsMEPrimary=new TH3F("fDCAXYZforcleanPionsMEPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsMEPrimary);
	fDCAXYZforcleanAntiPionsMEPrimary=new TH3F("fDCAXYZforcleanAntiPionsMEPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsMEPrimary);
	fDCAXYZforcleanPionsMESecondary=new TH3F("fDCAXYZforcleanPionsMESecondary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsMESecondary);
	fDCAXYZforcleanAntiPionsMESecondary=new TH3F("fDCAXYZforcleanAntiPionsMESecondary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsMESecondary);
	
	fDCAXYZforcleanPionsR=new TH3F("fDCAXYZforcleanPionsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanPionsR);
	fDCAXYZforcleanAntiPionsR=new TH3F("fDCAXYZforcleanAntiPionsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiPionsR);
	//Secondrary Protons Hadronic
	fDCAXYZforcleanProtonsR=new TH3F("fDCAXYZforcleanProtonsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanProtonsR);
	fDCAXYZforcleanAntiProtonsR=new TH3F("fDCAXYZforcleanAntiProtonsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZforcleanAntiProtonsR);
	
	
	
	fDCAXYZOpenforcleanPionsMCPrimary=new TH3F("fDCAXYZOpenforcleanPionsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsMCPrimary);
	fDCAXYZOpenforcleanAntiPionsMCPrimary=new TH3F("fDCAXYZOpenforcleanAntiPionsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsMCPrimary);
	fDCAXYZOpenforcleanProtonsMCPrimary=new TH3F("fDCAXYZOpenforcleanProtonsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanProtonsMCPrimary);
	fDCAXYZOpenforcleanAntiProtonsMCPrimary=new TH3F("fDCAXYZOpenforcleanAntiProtonsMCPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiProtonsMCPrimary);
	
	fDCAXYZOpenforcleanPionsWD=new TH3F("fDCAXYZOpenforcleanPionsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsWD);
	fDCAXYZOpenforcleanAntiPionsWD=new TH3F("fDCAXYZOpenforcleanAntiPionsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsWD);
	//Secondrary Protons weak deacy
	
	fDCAXYZOpenforcleanProtonsWD=new TH3F("fDCAXYZOpenforcleanProtonsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanProtonsWD);
	fDCAXYZOpenforcleanAntiProtonsWD=new TH3F("fDCAXYZOpenforcleanAntiProtonsWD",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiProtonsWD);

	fDCAXYZOpenforcleanPionsHI=new TH3F("fDCAXYZOpenforcleanPionsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsHI);
	fDCAXYZOpenforcleanAntiPionsHI=new TH3F("fDCAXYZOpenforcleanAntiPionsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsHI);
	//Secondrary Protons Hadronic
	fDCAXYZOpenforcleanProtonsHI=new TH3F("fDCAXYZOpenforcleanProtonsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanProtonsHI);
	fDCAXYZOpenforcleanAntiProtonsHI=new TH3F("fDCAXYZOpenforcleanAntiProtonsHI",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiProtonsHI);
	//Secondrary Pions mu el
	
	fDCAXYZOpenforcleanPionsMEPrimary=new TH3F("fDCAXYZOpenforcleanPionsMEPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsMEPrimary);
	fDCAXYZOpenforcleanAntiPionsMEPrimary=new TH3F("fDCAXYZOpenforcleanAntiPionsMEPrimary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsMEPrimary);
	fDCAXYZOpenforcleanPionsMESecondary=new TH3F("fDCAXYZOpenforcleanPionsMESecondary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsMESecondary);
	fDCAXYZOpenforcleanAntiPionsMESecondary=new TH3F("fDCAXYZOpenforcleanAntiPionsMESecondary",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsMESecondary);
	
	fDCAXYZOpenforcleanPionsR=new TH3F("fDCAXYZOpenforcleanPionsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanPionsR);
	fDCAXYZOpenforcleanAntiPionsR=new TH3F("fDCAXYZOpenforcleanAntiPionsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiPionsR);
	//Secondrary Protons Hadronic
	fDCAXYZOpenforcleanProtonsR=new TH3F("fDCAXYZOpenforcleanProtonsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanProtonsR);
	fDCAXYZOpenforcleanAntiProtonsR=new TH3F("fDCAXYZOpenforcleanAntiProtonsR",";Pt[GeV/c];dcaxy[cm];dcaz[cm]",kPtBins,binsPtDummy,kDCABins,binsDCADummy,kDCABins,binsDCADummy);
	flist->Add(fDCAXYZOpenforcleanAntiProtonsR);
	
	
	
	fElectronsource=new TH2F("fElectronsource","electrons;Pt [GeV/c];Process",kPtBins,binsPtDummy, kMaxMCProcess,0,kMaxMCProcess);
	flist->Add(fElectronsource);
	fAntiElectronsource=new TH2F("fAntiElectronsource","positrons;Pt [GeV/c];Process",kPtBins,binsPtDummy, kMaxMCProcess,0,kMaxMCProcess);
	flist->Add(fAntiElectronsource);
	fMuonsource=new TH2F("fMuonsource","electrons;Pt [GeV/c];Process",kPtBins,binsPtDummy, kMaxMCProcess,0,kMaxMCProcess);
	flist->Add(fMuonsource);
	 fAntiMuonsource=new TH2F("fAntiMuonsource","positrons;Pt [GeV/c];Process",kPtBins,binsPtDummy, kMaxMCProcess,0,kMaxMCProcess);
	flist->Add(fAntiMuonsource);
	
	fPrimaryElectronsMother=new TH1F("fPrimaryElectronsMother",";pdg code",4990,10.5,5000.5);
	flist->Add(fPrimaryElectronsMother);
	
	Double_t type[13]={-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5};
	Double_t cutlevel[10]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
	
	fParticlesCutmonitoring=new TH3F("fParticlesCutmonitoring",";particle;cut;Pt [GeV/c]",12,type,9,cutlevel,kPtBins,binsPtDummy);

	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(1,"pion");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(2,"kaon");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(3,"proton");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(4,"antipion");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(5,"antikaon");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(6,"antiproton");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(7,"pionfake");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(8,"kaonfake");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(9,"protonfake");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(10,"antipionfake");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(11,"antikaonfake");
	fParticlesCutmonitoring->GetXaxis()->SetBinLabel(12,"antiprotonfake");
	
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(1,"TPCin");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(2,"TPCrefit");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(3,"nTPCclu");	
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(4,"chi2");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(5,"ITSrefit");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(6,"SPDany");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(7,"standard");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(8,"ITSpid");
	fParticlesCutmonitoring->GetYaxis()->SetBinLabel(9,"DCA");
	flist->Add(fParticlesCutmonitoring);
	
	fVertexshift=new TH3F("fVertexshift",";#delta_{x};#delta_{y};#delta_{z}",50,-0.06,0.06,50,-0.06,0.06,50,-2,2);
	flist->Add(fVertexshift);
	
	Double_t deltapttpc[41];
	Double_t deltaptall[41];
	for(int i=0;i<41;i++)
	{
		deltapttpc[i]=-0.8+i*(1.6/40);
		deltaptall[i]=-0.2+i*(0.4/40);
	}
	fPtESDminusPtMCvPtESDafterallcuts= new TH3F("fPtESDminusPtMCvPtESDafterallcuts",";#delta_{PtESD-PtMC};PtESD;type",40,deltaptall,kPtBins,binsPtDummy,12,type);
	flist->Add(fPtESDminusPtMCvPtESDafterallcuts);
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(1,"pion");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(2,"kaon");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(3,"proton");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(4,"antipion");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(5,"antikaon");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(6,"antiproton");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(7,"pionfake");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(8,"kaonfake");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(9,"protonfake");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(10,"antipionfake");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(11,"antikaonfake");
	fPtESDminusPtMCvPtESDafterallcuts->GetZaxis()->SetBinLabel(12,"antiprotonfake");
	
	fPtESDminusPtMCvPtESDafterTPCcuts= new TH3F("fPtESDminusPtMCvPtESDafterTPCcuts",";#delta_{PtESD-PtMC};PtESD;type",40,deltapttpc,kPtBins,binsPtDummy,12,type);
	flist->Add(fPtESDminusPtMCvPtESDafterTPCcuts);
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(1,"pion");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(2,"kaon");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(3,"proton");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(4,"antipion");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(5,"antikaon");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(6,"antiproton");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(7,"pionfake");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(8,"kaonfake");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(9,"protonfake");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(10,"antipionfake");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(11,"antikaonfake");
	fPtESDminusPtMCvPtESDafterTPCcuts->GetZaxis()->SetBinLabel(12,"antiprotonfake");
	
	
	fMulESDMulMCVz=new TH3F("fMulESDMulMCVz",";NtracksESD;NparticlesMC;Vrt_z ",50,0,50,100,0,100,20,-10,10);
	flist->Add(fMulESDMulMCVz);
	PostData(1,  flist);
	Printf("end of CreateOutputObjects with MC");
}
//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::LocalInit() 
{
	//LocalInit
	Printf("end of LocalInit");
}

//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::UserExec(Option_t *) 
{
	//UserExec
	Bool_t isphysevent=0;
	Bool_t isgoodvertex=0;
	Bool_t isvxerteinZ=0;
	 fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	if (!fESD) 
	{
		Printf("ERROR: fESD not available");
		return;
	}
	
	Float_t refmultiplicity=fCutsMul->CountAcceptedTracks(fESD);
	if(!fSPD)
	{
		if(fLowMultiplicity>-1)
		{
			if(refmultiplicity<fLowMultiplicity)
				return;
		}
		if(fUpMultiplicity>-1)
		{
			if(refmultiplicity>fUpMultiplicity)
				return;
		}
	}	
	AliStack* stack=0x0;
	Double_t mcXvertex=0.0;
	Double_t mcYvertex=0.0;
	Double_t mcZvertex=0.0;
	
	if(fMC)
	{
		AliMCEvent* mcEvent  = (AliMCEvent*) MCEvent();
		//Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
		stack = mcEvent->Stack();
		mcXvertex=mcEvent->GetPrimaryVertex()->GetX();
		mcYvertex=mcEvent->GetPrimaryVertex()->GetY();
		mcZvertex=mcEvent->GetPrimaryVertex()->GetZ();
	}	
	
	
	fHistStats->Fill(0);
	//Event selection 
	//if( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()==0)	
	  UInt_t isSelected = 0;
	 if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))
	  	isSelected=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	//Printf("Mask_selection %u %u", isSelected,AliVEvent::kMB);  
        if(!(isSelected&AliVEvent::kMB))
        {
		 isphysevent=0;
		// Printf("No phys event.........\n");
	}	
	else
	{
		 isphysevent=1;
		 if(!fHIsettings)
			fHistStats->Fill(1);
	}
	if(isphysevent&&fHIsettings)
	{	
		AliCentrality *centrality = fESD->GetCentrality();
		if(!(fLowCentrality<0.0)&&fUpCentrality>0.0)
		{
			if(!centrality->IsEventInCentralityClass(fLowCentrality,fUpCentrality,"V0M"))
				return;
			else
				fHistStats->Fill(1);
		}	
	}
	//Printf("MC 1");
	//Good vertex	
	const AliESDVertex *vertex = 0x0;
	if(isphysevent)
	{
		vertex = fESD->GetPrimaryVertexTracks();
		if(vertex->GetNContributors()<1) 
		{
			// SPD vertex
			vertex = fESD->GetPrimaryVertexSPD();
			if(vertex->GetNContributors()<1) 
			{
				//Printf("No good  Vertex.........\n");
				isgoodvertex=0;
			}
			else
			{
				isgoodvertex=1;
				//fHistStats->Fill(2);	
				fHistZVertexBeforeCut->Fill(vertex ->GetZ());
				fHistXYVertexBeforeCut->Fill(vertex ->GetX(),vertex ->GetY());
			}	
		}
		else
		{
			isgoodvertex=1;	
			//fHistStats->Fill(2);	
			fHistZVertexBeforeCut->Fill(vertex ->GetZ());
			fHistXYVertexBeforeCut->Fill(vertex ->GetX(),vertex ->GetY()); 
		}
		if(isgoodvertex&&fUsePilerejection)
		{
			if(fESD->IsPileupFromSPDInMultBins())
				isgoodvertex=0;
		}
		if(isgoodvertex)
		{	
			if(TMath::Abs(vertex ->GetZ())>10.0)
			{
				//Printf("No good  Z of Vertex.........\n");
				isvxerteinZ=0;
			}
			else
				isvxerteinZ=1;
		}	
	}
	
	
	if(fdovertexrescuts&&fMC)
	{
		cout<<TMath::Abs(vertex->GetX()-mcXvertex)<<" "<<TMath::Abs(vertex->GetY()-mcYvertex)<<" "<<TMath::Abs(vertex->GetZ()-mcZvertex)<<endl;
		if(TMath::Abs(vertex->GetX()-mcXvertex)>0.015||TMath::Abs(vertex->GetY()-mcYvertex)>0.015||TMath::Abs(vertex->GetZ()-mcZvertex)>0.15)
			isvxerteinZ=0;	
	}  
	Float_t spdCorr=-1.0;
	if(isgoodvertex)
	{
		const AliMultiplicity *mult = fESD->GetMultiplicity();
		Float_t nClusters[6]={0.0,0.0,0.0,0.0,0.0,0.0};
		for(Int_t ilay=0; ilay<6; ilay++)
		{
			nClusters[ilay] = (Float_t)mult->GetNumberOfITSClusters(ilay);
		} 
		spdCorr = AliESDUtils::GetCorrSPD2(nClusters[1],vertex->GetZ());
		if(fSPD)
		{
			if(fLowMultiplicity>-1)
			{
				if(((Int_t)spdCorr)<fLowMultiplicity)
				{
					PostData(1,  flist);				
					return;
				}	
			}
			if(fUpMultiplicity>-1)
			{
				if(((Int_t)spdCorr)>fUpMultiplicity)
				{
					PostData(1,  flist);
					return;
				}	
			}		
		}
		fHistStats->Fill(2);
	}
	
	//Printf("MC 2");
	Int_t fMCmult=0;
	if(stack&&fMC)//Looping over MC information of all events
	{
		Float_t minpt=0.0;
		Float_t maxpt=0.0;
		Float_t mineta=0.0;
		Float_t maxeta=0.0;
		//Printf("MC 12");
		fCutsMul->GetPtRange(minpt,maxpt);
		fCutsMul->GetEtaRange(mineta,maxeta);
		fHistStats->Fill(8);
		if(TMath::Abs(mcZvertex)<10.0)
			fHistStats->Fill(9);
		for (int imc=0;imc<stack->GetNtrack();imc++)
		{
			if(!(stack->IsPhysicalPrimary(imc)))
				continue;
			TParticle *particleMC = stack->Particle(imc);
			if(!particleMC)
				continue;
			Int_t pdgcodeMC = particleMC->GetPdgCode();
			if(!(pdgcodeMC==211||pdgcodeMC==-211||pdgcodeMC==321||pdgcodeMC==-321||pdgcodeMC==2212||pdgcodeMC==-2212))	
				continue;
			if(particleMC->Pt()>minpt&&particleMC->Pt()<maxpt&&particleMC->Eta()>mineta&&particleMC->Eta()<maxeta)
				fMCmult++;	
			if (TMath::Abs(particleMC->Y())>fYCut)
				continue;
			if (particleMC->Pt()>2.0)
				continue;
			//Printf("%d aa",imc);			
		//	Printf("MC 22");
			if(pdgcodeMC==211)	
				fHistminsignalifPionPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			if(pdgcodeMC==-211)	
				fHistminsignalifAntiPionPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			if(pdgcodeMC==321)
				fHistminsignalifKaonPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			if(pdgcodeMC==-321)
				fHistminsignalifAntiKaonPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			if(pdgcodeMC==2212)
				fHistminsignalifProtonPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			if(pdgcodeMC==-2212)
				fHistminsignalifAntiProtonPMCPrimaryBeforeEventCuts->Fill(particleMC->Pt());
			
			if(TMath::Abs(mcZvertex)<10.0)
			{		
				if(pdgcodeMC==211)	
					fHistminsignalifPionPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
				if(pdgcodeMC==-211)	
					fHistminsignalifAntiPionPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
				if(pdgcodeMC==321)
					fHistminsignalifKaonPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
				if(pdgcodeMC==-321)
					fHistminsignalifAntiKaonPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
				if(pdgcodeMC==2212)
					fHistminsignalifProtonPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
				if(pdgcodeMC==-2212)
					fHistminsignalifAntiProtonPMCPrimaryBeforeEventCutswithgoodZvertex->Fill(particleMC->Pt());
			}
			if(!isphysevent)
				continue;
			if(!isgoodvertex)
				continue;	
			if(pdgcodeMC==211)	
				fHistminsignalifPionPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());
			if(pdgcodeMC==-211)	
				fHistminsignalifAntiPionPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());
			if(pdgcodeMC==321)
				fHistminsignalifKaonPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());
			if(pdgcodeMC==-321)
				fHistminsignalifAntiKaonPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());
			if(pdgcodeMC==2212)
				fHistminsignalifProtonPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());
			if(pdgcodeMC==-2212)
				fHistminsignalifAntiProtonPMCPrimaryAfterEventCutsBeforeVertexZ->Fill(particleMC->Pt());	
			if(!isvxerteinZ)
				continue;		
			if(pdgcodeMC==211)	
			{
				fHistminsignalifPionPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-211)
			{	
				fHistminsignalifAntiPionPMCPrimary->Fill(particleMC->Pt());
			}		
			if(pdgcodeMC==321)
			{	
				fHistminsignalifKaonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-321)
			{	
				fHistminsignalifAntiKaonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==2212)
			{	
				fHistminsignalifProtonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-2212)
			{	
				fHistminsignalifAntiProtonPMCPrimary->Fill(particleMC->Pt());
			}		
							
		}
	}
	else if(fMC)
	    return;
	else 
		Printf("Data mode \n"); 
	     		 
	if(!(isphysevent&&isgoodvertex&&isvxerteinZ))
	{
		//Printf("No Good event.........\n");
		PostData(1,  flist);
		//Printf("end of Exec");
		return;
	}
	
	

	fHistStandartMulvSPD2->Fill(refmultiplicity,spdCorr);
	fHistStats->Fill(3);
		
	fHistZVertexAfterCut->Fill(vertex ->GetZ());
	fHistXYVertexAfterCut->Fill(vertex ->GetX(),vertex ->GetY()); 
	
	if(fMC)
	{
		fVertexshift->Fill(vertex->GetX()-mcXvertex,vertex->GetY()-mcYvertex,vertex->GetZ()-mcZvertex);
		fMulESDMulMCVz->Fill(refmultiplicity,fMCmult,vertex->GetZ());
	}		
	if(fCuts==0)
	{
		//Printf("No CUTS Defined.........\n");
		PostData(1,  flist);
		//Printf("end of Exec");
		return;
	}
	
	//Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
	 Int_t nTracks=fESD->GetNumberOfTracks();
	
	Int_t mynumberoftracks=0;
	 AliESDtrack *trackESD=0;
	 
	const Float_t pionmass=AliPID::ParticleMass(2);
	const Float_t kaonmass=AliPID::ParticleMass(3);
	const Float_t protonmass=AliPID::ParticleMass(4);

 	for(int tr1=0;tr1<nTracks;tr1++)
  	{	
		
		trackESD=fESD->GetTrack(tr1);	
		//fHistStats->Fill(2);
		
		Double_t pt=trackESD->Pt()*trackESD->GetSign();
		Double_t p=trackESD->P();
		Double_t eta=trackESD->Eta();
		Double_t phi=trackESD->Phi();
		Float_t dcaxy=0.0; 
		Float_t dcaz=0.0;
		trackESD->GetImpactParameters(dcaxy,dcaz);
		Double_t pz=trackESD->Pz();
		 UShort_t nTPCclusters=trackESD->GetTPCNcls();
		 Float_t chi2=trackESD->GetTPCchi2();
		 if(nTPCclusters>0)
		 	chi2=chi2/((Float_t)nTPCclusters);
		else
			chi2=-1.0;
		if(!trackESD->IsOn(AliESDtrack::kTPCin))
			continue;			
		//Y assumtion
		Float_t yforpion=0.5*TMath::Log((TMath::Sqrt(pionmass*pionmass+p*p)+pz)/(TMath::Sqrt(pionmass*pionmass+p*p)-pz));
		Float_t yforkaon=0.5*TMath::Log((TMath::Sqrt(kaonmass*kaonmass+p*p)+pz)/(TMath::Sqrt(kaonmass*kaonmass+p*p)-pz));
		Float_t yforproton=0.5*TMath::Log((TMath::Sqrt(protonmass*protonmass+p*p)+pz)/(TMath::Sqrt(protonmass*protonmass+p*p)-pz));
		
		if(TMath::Abs(yforpion)>fYCut&&TMath::Abs(yforkaon)>fYCut&&TMath::Abs(yforproton)>fYCut) //go trought one y cut
			continue;
		Int_t label=-1;
		if(fMC)
			label=trackESD->GetLabel();
		//if(label<0)	
	//	Printf("label %d %f %f %f %f %d %f %f\n",label,p,pt,eta,chi2,nTPCclusters,dcaxy,dcaz);	
		Int_t pdgcode=0;
		Int_t primary=0;
		Double_t chargeMC=1.0;
		Float_t etaMC=10.0;
		Float_t ptMC=10.0;
		Int_t   uniqueID=-1;
		Int_t pdgcodefake=0;
		Int_t primaryfake=0;
		TParticle *particle2=0x0;
		if(label>=0&&stack&&fMC)
		{
			primary=stack->IsPhysicalPrimary(TMath::Abs(label));
			particle2 = stack->Particle(TMath::Abs(label));
			pdgcode=particle2->GetPdgCode();
			chargeMC=particle2->GetPDG(0)->Charge()/3.0;
			etaMC=particle2->Eta();
			ptMC=particle2->Pt();
			uniqueID=particle2->GetUniqueID();
		}
		if(label<0&&stack&&fMC)
		{
			primaryfake=stack->IsPhysicalPrimary(TMath::Abs(label));
			particle2 = stack->Particle(TMath::Abs(label));
			pdgcodefake=particle2->GetPdgCode();
			uniqueID=particle2->GetUniqueID();
			
		}	
		
		Int_t typeParticle=-10;
		if((primaryfake||primary))
		{
			
			if((pdgcodefake==211||pdgcode==211)&&TMath::Abs(yforpion)<fYCut)
				typeParticle=0;
			if((pdgcodefake==321||pdgcode==321)&&TMath::Abs(yforkaon)<fYCut)
				typeParticle=1;
			if((pdgcodefake==2212||pdgcode==2212)&&TMath::Abs(yforproton)<fYCut)
				typeParticle=2;
			if((pdgcodefake==-211||pdgcode==-211)&&TMath::Abs(yforpion)<fYCut)
				typeParticle=3;
			if((pdgcodefake==-321||pdgcode==-321)&&TMath::Abs(yforkaon)<fYCut)
				typeParticle=4;
			if((pdgcodefake==-2212||pdgcode==-2212)&&TMath::Abs(yforproton)<fYCut)
				typeParticle=5;
			
			if(primaryfake)	
				typeParticle+=6;
		}	
		
		fTracksCutmonitoring->Fill(1,TMath::Abs(pt));
		if(fMC)
		{
			fParticlesCutmonitoring->Fill(typeParticle,1,TMath::Abs(pt));	
			if(trackESD->IsOn(AliESDtrack::kTPCrefit))
			{
				fParticlesCutmonitoring->Fill(typeParticle,2,TMath::Abs(pt));
				if(nTPCclusters>70)
				{
					fParticlesCutmonitoring->Fill(typeParticle,3,TMath::Abs(pt));
					if(chi2<4.0)
					{
						fParticlesCutmonitoring->Fill(typeParticle,4,TMath::Abs(pt));
						fPtESDminusPtMCvPtESDafterTPCcuts->Fill(TMath::Abs(pt)-particle2->Pt(),TMath::Abs(pt),typeParticle);
						if(trackESD->IsOn(AliESDtrack::kITSrefit))
						{
							fParticlesCutmonitoring->Fill(typeParticle,5,TMath::Abs(pt));
							if(trackESD->HasPointOnITSLayer(0)||trackESD->HasPointOnITSLayer(1))
							{
								fParticlesCutmonitoring->Fill(typeParticle,6,TMath::Abs(pt));
							}
						}
					}		
					
				}	
			}
		}	
			
			
		fHistPhiPtBeforeCuts->Fill(phi,pt);//phi pt
		fHistEtaPtBeforeCuts->Fill(eta,pt);
		fHistDCABeforeCuts->Fill(dcaxy,dcaz);	
		
		

			
		//standart cuts 	
		if(fCuts->AcceptTrack(trackESD)==kFALSE)
			continue;
		fTracksCutmonitoring->Fill(2,TMath::Abs(pt));	
		if(fMC)
			fParticlesCutmonitoring->Fill(typeParticle,7,TMath::Abs(pt));
		//Tpc pid cut for debug	
		Double_t pinTPC=trackESD->GetTPCInnerParam()->GetP();//momentum in primary vertex taken from TPC tracking
		Double_t pinTPCglobal=trackESD->GetInnerParam()->GetP();//momentum at the inner  wall of the TPC taken from global tracking
		Float_t sigKaon= fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
		Float_t sigProton= fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
		Double_t tpcSignal =trackESD ->GetTPCsignal();
		
		if(fTPCPIDCUT)
		{
			if(fTPCPIDCUT->AcceptTrack(trackESD,fESD)==kFALSE)
				continue;
		}
				
		Bool_t cutDCA=SelectOnImpPar(trackESD);	
		
		//fHistStats->Fill(4);
		Double_t tmpQESD[4]={-1.0,-1.0,-1.0,-1.0};
		trackESD->GetITSdEdxSamples(tmpQESD);
		if(fCorrectSDD&&fMC)
			CorrectSDD(tmpQESD);	
		if(fCorrectSSD&&fMC)
			CorrectSSD(tmpQESD);	
		Int_t nSSDSDD=0;
		Int_t nSSDSDDall=0;
		
		
		
		if(TMath::Abs(yforpion)<fYCut&&cutDCA) 
			fHistNtrackwithstandardcuts->Fill(TMath::Abs(pt),0);
		if(TMath::Abs(yforkaon)<fYCut&&cutDCA) 
		{
			fHistNtrackwithstandardcuts->Fill(TMath::Abs(pt),1);
			if(pt>0.0)
				fHistSignalinTPCKaonforstandardcuts->Fill(pt,(tpcSignal-sigKaon)/sigKaon);
			else
				fHistSignalinTPCAntiKaonforstandardcuts->Fill(TMath::Abs(pt),(tpcSignal-sigKaon)/sigKaon);	
		}
		if(TMath::Abs(yforproton)<fYCut&&cutDCA) 
		{
			fHistNtrackwithstandardcuts->Fill(TMath::Abs(pt),2);
			if(pt>0.0)
				fHistSignalinTPCProtonforstandardcuts->Fill(pt,(tpcSignal-sigProton)/sigProton);
			else
				fHistSignalinTPCAntiProtonforstandardcuts->Fill(TMath::Abs(pt),(tpcSignal-sigProton)/sigProton);	
		}		
		
			
		for (int iL=0;iL<4;iL++)
		{
			if(tmpQESD[iL]>fchargeCut)
				nSSDSDD++;
			if(tmpQESD[iL]>0.0)
				nSSDSDDall++;	
		}	
		if(nSSDSDDall>=3)
			fHistStats->Fill(4);
		//ITS PId cut
		if(nSSDSDD<3)
		{
		//	cout<<"BAD "<<stack->Particle(TMath::Abs(trackESD->GetLabel()))->GetPdgCode()<<endl;
			continue;
		}	
		fTracksCutmonitoring->Fill(3,TMath::Abs(pt));
		if(fMC)
			fParticlesCutmonitoring->Fill(typeParticle,8,TMath::Abs(pt));
		if(cutDCA)
		{
			fTracksCutmonitoring->Fill(4,TMath::Abs(pt));
			if(fMC)
			{
				fParticlesCutmonitoring->Fill(typeParticle,9,TMath::Abs(pt));	
				fPtESDminusPtMCvPtESDafterallcuts->Fill(TMath::Abs(pt)-particle2->Pt(),TMath::Abs(pt),typeParticle);
			}		
		}
		if(TMath::Abs(yforpion)<fYCut&&cutDCA) 
		{
			fHistNtrackwithITSPIDcuts->Fill(TMath::Abs(pt),0);
			if(pt>0.0)
			{
				 fPionNTPCClusters->Fill(pt,nTPCclusters);
				 fPionchi2->Fill(pt,chi2);
			}
			else
			{
				fAntiPionNTPCClusters->Fill(TMath::Abs(pt),nTPCclusters);
				fAntiPionchi2->Fill(TMath::Abs(pt),chi2);
			}		
		}	
		if(TMath::Abs(yforkaon)<fYCut&&cutDCA) 
		{
			fHistNtrackwithITSPIDcuts->Fill(TMath::Abs(pt),1);
			if(pt>0.0)
				fHistSignalinTPCKaonforITSPIDcuts->Fill(pt,(tpcSignal-sigKaon)/sigKaon);
			else
				fHistSignalinTPCAntiKaonforITSPIDcuts->Fill(TMath::Abs(pt),(tpcSignal-sigKaon)/sigKaon);	
		}	
		if(TMath::Abs(yforproton)<fYCut&&cutDCA) 
		{
			fHistNtrackwithITSPIDcuts->Fill(TMath::Abs(pt),2);
			if(pt>0.0)
				fHistSignalinTPCProtonforITSPIDcuts->Fill(pt,(tpcSignal-sigProton)/sigProton);
			else
				fHistSignalinTPCAntiProtonforITSPIDcuts->Fill(TMath::Abs(pt),(tpcSignal-sigProton)/sigProton);	
		}
		fHistStats->Fill(5);				
		fHistPhiPtAfterCuts->Fill(phi,pt);
		fHistEtaPtAfterCuts->Fill(eta,pt);
		fHistDCAAfterCuts->Fill(dcaxy,dcaz);
		fHistPminusTPCinPAfterCuts->Fill(pinTPC-p,pinTPC);	
		fHistPminusTPCinPglobalAfterCuts->Fill(pinTPC-p,pinTPCglobal);
		if(tmpQESD[0]>fchargeCut)
		{
			fHistL3dEP->Fill(p,tmpQESD[0]);
			fHistL3dETPCinP->Fill(pinTPC,tmpQESD[0]);
		}	
		if(tmpQESD[1]>fchargeCut)
		{
			fHistL4dEP->Fill(p,tmpQESD[1]);
			fHistL4dETPCinP->Fill(pinTPC,tmpQESD[1]);
		}	
		if(tmpQESD[2]>fchargeCut)
		{
			fHistL5dEP->Fill(p,tmpQESD[2]);
			fHistL5dETPCinP->Fill(pinTPC,tmpQESD[2]);
		}
		if(tmpQESD[3]>fchargeCut)
		{
			fHistL6dEP->Fill(p,tmpQESD[3]);
			fHistL6dETPCinP->Fill(pinTPC,tmpQESD[3]);
		}
		Float_t myITSsignal=0.0;
		Float_t minITSsignal=0.0;
		Int_t whichLmin=-1;
		Int_t nosignaL=-1;
		if(nSSDSDD==3)
		{
			Double_t tmp2QESD[3]={-1.0,-1.0,-1.0};
			Int_t iLnotZero=0;
			for (int iL=0;iL<4;iL++)
			{
				if(tmpQESD[iL]>fchargeCut)
				{
					tmp2QESD[iLnotZero]=tmpQESD[iL];
					iLnotZero++;
				}
				else
					nosignaL=iL;
			}
			whichLmin=TMath::LocMin(3,tmp2QESD);
			if(nosignaL>-1&&nosignaL<=whichLmin)
				whichLmin++;
			minITSsignal=TMath::MinElement(3,tmp2QESD);
			myITSsignal=MyITSsignalusing3points(tmp2QESD);
		}	
		if(nSSDSDD==4)
		{
			myITSsignal=MyITSsignalusing4points(tmpQESD);
			whichLmin=TMath::LocMin(4,tmpQESD);	
			minITSsignal=TMath::MinElement(4,tmpQESD);
		}		
	
		if(whichLmin==0)
			fHistwhichhasmin->Fill(0.0,tmpQESD[0]);
		if(whichLmin==1)
			fHistwhichhasmin->Fill(1.0,tmpQESD[1]);
		if(whichLmin==2)
			fHistwhichhasmin->Fill(2.0,tmpQESD[2]);
		if(whichLmin==3)
			fHistwhichhasmin->Fill(3.0,tmpQESD[3]);
		if(pt>0.0)
		{			
			fHistMydEPpositive->Fill(p,myITSsignal);
			fHistMydETPCinPglobalpositive->Fill(pinTPCglobal,myITSsignal);
			fHistMydETPCinPpositive->Fill(pinTPC,myITSsignal);
		}
		else
		{
			fHistMydEPnegative->Fill(p,myITSsignal);
			fHistMydETPCinPglobalnegative->Fill(pinTPCglobal,myITSsignal);
			fHistMydETPCinPnegative->Fill(pinTPC,myITSsignal);
		}	
		Float_t signaltouse=myITSsignal;	
			
		
		Float_t itspidsignalforpions=TMath::Log(signaltouse)-TMath::Log(fESDpid->GetITSResponse().Bethe(p,pionmass,kFALSE));
		Float_t itspidsignalforkaons=TMath::Log(signaltouse)-TMath::Log(fESDpid->GetITSResponse().Bethe(p,kaonmass,kFALSE));
		Float_t itspidsignalforprotons=TMath::Log(signaltouse)-TMath::Log(fESDpid->GetITSResponse().Bethe(p,protonmass,kFALSE));
		if(cutDCA)
		{
			mynumberoftracks++;
			if(nSSDSDD==4)
				fHistMysignalminusESD->Fill((signaltouse-trackESD->GetITSsignal())/signaltouse);	
		}
		//Printf("Select on clean \n");	
		if(TMath::Abs(yforpion)<=fYCut)
		{
			Float_t weight=1.0;
			if(fMC)
				weight=GetWeight(label,stack);	
			if(pt>0.0)
			{
				if(cutDCA)
				{
					if(fMC)
						fHistminsignalifPionP->Fill(pt,itspidsignalforpions,weight);
					else
						fHistminsignalifPionP->Fill(pt,itspidsignalforpions);
					if(itspidsignalforpions>-0.5&&itspidsignalforpions<0.2) //select on clean
					{
						fDCAXYZforcleanPions->Fill(pt,dcaxy,dcaz);
						if(fMC)
						{
							if(primary&&pdgcode==211)
								fDCAXYZforcleanPionsMCPrimary->Fill(pt,dcaxy,dcaz);
							else if(!primary&&pdgcode==211&&uniqueID==kPDecay)
								fDCAXYZforcleanPionsWD->Fill(pt,dcaxy,dcaz);
							else if(!primary&&pdgcode==211&&uniqueID==kPHadronic)
								fDCAXYZforcleanPionsHI->Fill(pt,dcaxy,dcaz);
							else if(primary&&(pdgcode==-11||pdgcode==-13))
								fDCAXYZforcleanPionsMEPrimary->Fill(pt,dcaxy,dcaz);
							else if(!primary&&(pdgcode==-11||pdgcode==-13))
								fDCAXYZforcleanPionsMESecondary->Fill(pt,dcaxy,dcaz);
							else
								fDCAXYZforcleanPionsR->Fill(pt,dcaxy,dcaz);
					
						}	  	  	   
					}	 //select on clean
					if(primary)
					{	
						if(pdgcode==211)
							fHistminsignalifPionPPrimary->Fill(pt,itspidsignalforpions);
					}
					else 
					{
						if(pdgcode==211)							
							fHistminsignalifPionPSecondary->Fill(pt,itspidsignalforpions,weight);
					//cout<<pdgcode<<" "<<	uniqueID<<"  "<<kPDecay<<" "<<kPHadronic<<endl;
					}						
					if(pdgcode==-11||pdgcode==-13)
					{
						fHistminsignalifMuEPositiveP->Fill(pt,itspidsignalforpions);
						if(!primary)
						{
							if(pdgcode==-11)
							{
								fHistStats->Fill(6);
								fAntiElectronsource->Fill(pt,uniqueID);
							}	
							else if(pdgcode==-13)	
							{
								fHistStats->Fill(7);
								fAntiMuonsource->Fill(pt,uniqueID);
							}
						}
						else if(primary&&pdgcode==-11)
						{
							Printf("%d Mom",particle2->GetFirstMother());
							if(particle2->GetFirstMother()>-1)
								fPrimaryElectronsMother->Fill(TMath::Abs(stack->Particle(particle2->GetFirstMother())->GetPdgCode()));
							else
								fPrimaryElectronsMother->Fill(-1);
							fAntiElectronsource->Fill(pt,0);
						}			
						else if(primary&&pdgcode==-13)
							fAntiMuonsource->Fill(pt,0);
					}
					if(pdgcodefake==211)
					{
						if(primaryfake)
							fHistminsignalifPionPrimaryfake->Fill(pt,itspidsignalforpions);
						else
							fHistminsignalifPionSecondaryfake->Fill(pt,itspidsignalforpions);	
					}	
				}				
				if(itspidsignalforpions>-0.5&&itspidsignalforpions<0.2) //select on clean
				{
					fDCAXYZOpenforcleanPions->Fill(pt,dcaxy,dcaz);
					if(fMC)
					{
						if(primary&&pdgcode==211)
							fDCAXYZOpenforcleanPionsMCPrimary->Fill(pt,dcaxy,dcaz);
						else if(!primary&&pdgcode==211&&uniqueID==kPDecay)	 
							fDCAXYZOpenforcleanPionsWD->Fill(pt,dcaxy,dcaz);
						else if(!primary&&pdgcode==211&&uniqueID==kPHadronic)
							fDCAXYZOpenforcleanPionsHI->Fill(pt,dcaxy,dcaz);
						else if(primary&&(pdgcode==-11||pdgcode==-13))	
							fDCAXYZOpenforcleanPionsMEPrimary->Fill(pt,dcaxy,dcaz);
						else if(!primary&&(pdgcode==-11||pdgcode==-13))	
							fDCAXYZOpenforcleanPionsMESecondary->Fill(pt,dcaxy,dcaz);	
						else
							fDCAXYZOpenforcleanPionsR->Fill(pt,dcaxy,dcaz);					
					}	  	  	   
				}			
			}
			else 	
			{
				if(cutDCA)
				{
					if(fMC)
						fHistminsignalifAntiPionP->Fill(TMath::Abs(pt),itspidsignalforpions,weight);
					else	
						fHistminsignalifAntiPionP->Fill(TMath::Abs(pt),itspidsignalforpions);
					if(itspidsignalforpions>-0.5&&itspidsignalforpions<0.2)//select on clean
					{
						fDCAXYZforcleanAntiPions->Fill(TMath::Abs(pt),dcaxy,dcaz);
						if(fMC)
						{
							if(primary&&pdgcode==-211)
								fDCAXYZforcleanAntiPionsMCPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(!primary&&pdgcode==-211&&uniqueID==kPDecay)	 
								fDCAXYZforcleanAntiPionsWD->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(!primary&&pdgcode==-211&&uniqueID==kPHadronic)
								fDCAXYZforcleanAntiPionsHI->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(primary&&(pdgcode==11||pdgcode==13))	
								fDCAXYZforcleanAntiPionsMEPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(!primary&&(pdgcode==11||pdgcode==13))	
								fDCAXYZforcleanAntiPionsMESecondary->Fill(TMath::Abs(pt),dcaxy,dcaz);	
							else
								fDCAXYZforcleanAntiPionsR->Fill(TMath::Abs(pt),dcaxy,dcaz);	 
						}				
					}	//select on clean
					if(primary)
					{
						if(pdgcode==-211)
							fHistminsignalifAntiPionPPrimary->Fill(TMath::Abs(pt),itspidsignalforpions);
					}
					else 
					{
						if(pdgcode==-211)
							fHistminsignalifAntiPionPSecondary->Fill(TMath::Abs(pt),itspidsignalforpions,weight);
						//cout<<pdgcode<<" "<<	uniqueID<<"  "<<kPDecay<<" "<<kPHadronic<<endl;
					}
					if(pdgcode==11||pdgcode==13)
					{
						fHistminsignalifMuENegativeP->Fill(TMath::Abs(pt),itspidsignalforpions);	
						if(!primary)
						{
							if(pdgcode==11)
							{
								fHistStats->Fill(6);
								fElectronsource->Fill(TMath::Abs(pt),uniqueID);
							}	
							else if(pdgcode==13)	
							{
								fHistStats->Fill(7);
								fMuonsource->Fill(TMath::Abs(pt),uniqueID);
							}
						}
						else if(primary&&pdgcode==11)
						{
							Printf("%d Mom",particle2->GetFirstMother());
							if(particle2->GetFirstMother()>-1)
								fPrimaryElectronsMother->Fill(TMath::Abs(stack->Particle(particle2->GetFirstMother())->GetPdgCode()));
							else
								fPrimaryElectronsMother->Fill(-1);
							fElectronsource->Fill(TMath::Abs(pt),0);
						}	
						else if(primary&&pdgcode==13)
							fMuonsource->Fill(TMath::Abs(pt),0);
					}
					if(pdgcodefake==-211)
					{
						if(primaryfake)
							fHistminsignalifAntiPionPrimaryfake->Fill(TMath::Abs(pt),itspidsignalforpions);
						else
							fHistminsignalifAntiPionSecondaryfake->Fill(TMath::Abs(pt),itspidsignalforpions);	
					}	
				}			
				if(itspidsignalforpions>-0.5&&itspidsignalforpions<0.2)//select on clean
				{
					fDCAXYZOpenforcleanAntiPions->Fill(TMath::Abs(pt),dcaxy,dcaz);
					if(fMC)
					{
						if(primary&&pdgcode==-211)
							fDCAXYZOpenforcleanAntiPionsMCPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(!primary&&pdgcode==-211&&uniqueID==kPDecay)	 
							fDCAXYZOpenforcleanAntiPionsWD->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(!primary&&pdgcode==-211&&uniqueID==kPHadronic)
							fDCAXYZOpenforcleanAntiPionsHI->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(primary&&(pdgcode==11||pdgcode==13))	
							fDCAXYZOpenforcleanAntiPionsMEPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(!primary&&(pdgcode==11||pdgcode==13))	
							fDCAXYZOpenforcleanAntiPionsMESecondary->Fill(TMath::Abs(pt),dcaxy,dcaz);	
						else
							fDCAXYZOpenforcleanAntiPionsR->Fill(TMath::Abs(pt),dcaxy,dcaz);	 
					}				
				}	//sel								
			}	
		}	
		if(TMath::Abs(yforkaon)<=fYCut&&cutDCA)	
		{
			if(pt>0.0)
			{
				fHistminsignalifKaonP->Fill(pt,itspidsignalforkaons);
				if((itspidsignalforkaons>-0.2)&&(itspidsignalforkaons<0.2))
				{
					fKaonNTPCClusters->Fill(pt,nTPCclusters);
					fKaonchi2->Fill(pt,chi2);
				}		
				if(primary)
				{
					if(pdgcode==321)
						fHistminsignalifKaonPPrimary->Fill(pt,itspidsignalforkaons);
				}
				else 
				{
					if(pdgcode==321)
						fHistminsignalifKaonPSecondary->Fill(pt,itspidsignalforkaons);
				}
				if(pdgcodefake==321)
				{
					if(primaryfake)
						fHistminsignalifKaonPrimaryfake->Fill(pt,itspidsignalforkaons);
					else
						fHistminsignalifKaonSecondaryfake->Fill(pt,itspidsignalforkaons);
				}
				
				
			}	
			else
			{
				fHistminsignalifAntiKaonP->Fill(TMath::Abs(pt),itspidsignalforkaons);
				if((itspidsignalforkaons>-0.2)&&(itspidsignalforkaons<0.2))
				{
					fAntiKaonNTPCClusters->Fill(TMath::Abs(pt),nTPCclusters);
					fAntiKaonchi2->Fill(TMath::Abs(pt),chi2);
				}
				if(primary)
				{
					if(pdgcode==-321)
						fHistminsignalifAntiKaonPPrimary->Fill(TMath::Abs(pt),itspidsignalforkaons);
				}
				else 
				{
					if(pdgcode==-321)
						fHistminsignalifAntiKaonPSecondary->Fill(TMath::Abs(pt),itspidsignalforkaons);
				}	
				if(pdgcodefake==-321)
				{
					if(primaryfake)
						fHistminsignalifAntiKaonPrimaryfake->Fill(TMath::Abs(pt),itspidsignalforkaons);
					else
						fHistminsignalifAntiKaonSecondaryfake->Fill(TMath::Abs(pt),itspidsignalforkaons);
					
				}				
			}		
		}	
		if(TMath::Abs(yforproton)<=fYCut)
		{	
			Float_t weight=1.0;
			if(fMC)
				weight=GetWeight(label,stack);	
			if(pt>0.0)
			{
				if(cutDCA)
				{
					if(fMC)
						fHistminsignalifProtonP->Fill(pt,itspidsignalforprotons,weight);
					else	
						fHistminsignalifProtonP->Fill(pt,itspidsignalforprotons);
					
					//if(((itspidsignalforprotons))>(TMath::Abs(pt)<0.45?-0.2:0.0))&&(itspidsignalforprotons))<0.2)//select on
					// clean
					if((itspidsignalforprotons>-0.2)&&(itspidsignalforprotons<0.5))
					{
						fDCAXYZforcleanProtons->Fill(pt,dcaxy,dcaz);
						fProtonNTPCClusters->Fill(pt,nTPCclusters);
						fProtonchi2->Fill(pt,chi2);
						if(fMC)
						{	
							if(primary&&pdgcode==2212)
								fDCAXYZforcleanProtonsMCPrimary->Fill(pt,dcaxy,dcaz);
							else if(!primary&&!primaryfake&&(pdgcode==2212||pdgcodefake==2212)&&uniqueID==kPDecay)	 
								fDCAXYZforcleanProtonsWD->Fill(pt,dcaxy,dcaz);
							else if(!primary&&!primaryfake&&(pdgcode==2212||pdgcodefake==2212)&&uniqueID==kPHadronic)
								fDCAXYZforcleanProtonsHI->Fill(pt,dcaxy,dcaz);
							else 
								fDCAXYZforcleanProtonsR->Fill(pt,dcaxy,dcaz);
						}	 
					}//select on clean
					if(primary)
					{
						if(pdgcode==2212)
							fHistminsignalifProtonPPrimary->Fill(pt,itspidsignalforprotons);
					}
					else  if(primaryfake)
					{
						if(pdgcodefake==2212)
							fHistminsignalifProtonPPrimaryfake->Fill(pt,itspidsignalforprotons);
					}
					else
					{
						if(pdgcode==2212&&uniqueID==kPDecay)		
							fHistminsignalifProtonPSecondaryWD->Fill(pt,itspidsignalforprotons,weight);
						else if(pdgcode==2212&&uniqueID==kPHadronic)	
							fHistminsignalifProtonPSecondaryHI->Fill(pt,itspidsignalforprotons);
						else	if(pdgcodefake==2212&&uniqueID==kPDecay)
							fHistminsignalifProtonPSecondaryWDfake->Fill(pt,itspidsignalforprotons,weight);
						else if(pdgcodefake==2212&&uniqueID==kPHadronic)
							fHistminsignalifProtonPSecondaryHIfake->Fill(pt,itspidsignalforprotons);
						else	 if(fMC)
							fHistminsignalifProtonPSecondaryRest->Fill(pt,itspidsignalforprotons);
					}
				}
				
				if((itspidsignalforprotons>-0.2)&&(itspidsignalforprotons<0.5))
				{
					fDCAXYZOpenforcleanProtons->Fill(pt,dcaxy,dcaz);
					if(fMC)
					{	
						if(primary&&pdgcode==2212)
							fDCAXYZOpenforcleanProtonsMCPrimary->Fill(pt,dcaxy,dcaz);
						else if(!primary&&!primaryfake&&(pdgcode==2212||pdgcodefake==2212)&&uniqueID==kPDecay)	 
							fDCAXYZOpenforcleanProtonsWD->Fill(pt,dcaxy,dcaz);
						else if(!primary&&!primaryfake&&(pdgcode==2212||pdgcodefake==2212)&&uniqueID==kPHadronic)
							fDCAXYZOpenforcleanProtonsHI->Fill(pt,dcaxy,dcaz);
						else 
							fDCAXYZOpenforcleanProtonsR->Fill(pt,dcaxy,dcaz);
					}	 
				}//select on clean
					
			}		
			else
			{
				if(cutDCA)
				{
					if(fMC)
						fHistminsignalifAntiProtonP->Fill(TMath::Abs(pt),itspidsignalforprotons,weight);
					else	
						fHistminsignalifAntiProtonP->Fill(TMath::Abs(pt),itspidsignalforprotons);
					if((itspidsignalforprotons>-0.2)&&(itspidsignalforprotons<0.5))
					{//select on clean
						fDCAXYZforcleanAntiProtons->Fill(TMath::Abs(pt),dcaxy,dcaz);
						
						if(fMC)
						{
							if(primary&&pdgcode==-2212)
								fDCAXYZforcleanAntiProtonsMCPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(!primary&&!primaryfake&&(pdgcode==-2212||pdgcodefake==-2212)&&uniqueID==kPDecay)	 
								fDCAXYZforcleanAntiProtonsWD->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(!primary&&!primaryfake&&(pdgcode==-2212||pdgcodefake==-2212)&&uniqueID==kPHadronic)
								fDCAXYZforcleanAntiProtonsHI->Fill(TMath::Abs(pt),dcaxy,dcaz);
							else if(fMC)
								fDCAXYZforcleanAntiProtonsR->Fill(TMath::Abs(pt),dcaxy,dcaz);		 
						}  
					}//select on clean
					if(primary)
					{
						if(pdgcode==-2212)
							fHistminsignalifAntiProtonPPrimary->Fill(TMath::Abs(pt),itspidsignalforprotons);
					}
					else if(primaryfake)
					{
						if(pdgcodefake==-2212)
							fHistminsignalifAntiProtonPPrimaryfake->Fill(TMath::Abs(pt),itspidsignalforprotons);
					}
					else 
					{
						if(pdgcode==-2212&&uniqueID==kPDecay)
							fHistminsignalifAntiProtonPSecondaryWD->Fill(TMath::Abs(pt),itspidsignalforprotons,weight);
						else if(pdgcode==-2212&&uniqueID==kPHadronic)	
							fHistminsignalifAntiProtonPSecondaryHI->Fill(TMath::Abs(pt),itspidsignalforprotons);	
						else if(pdgcodefake==-2212&&uniqueID==kPDecay)
							fHistminsignalifAntiProtonPSecondaryWDfake->Fill(TMath::Abs(pt),itspidsignalforprotons,weight);
						else if(pdgcodefake==-2212&&uniqueID==kPHadronic)	
							fHistminsignalifAntiProtonPSecondaryHIfake->Fill(TMath::Abs(pt),itspidsignalforprotons);
						else if(fMC)
							fHistminsignalifAntiProtonPSecondaryRest->Fill(TMath::Abs(pt),itspidsignalforprotons);	
					}
				}
				if((itspidsignalforprotons>-0.2)&&(itspidsignalforprotons<0.5))
				{//select on clean
					fDCAXYZOpenforcleanAntiProtons->Fill(TMath::Abs(pt),dcaxy,dcaz);
					fAntiProtonNTPCClusters->Fill(TMath::Abs(pt),nTPCclusters);
					fAntiProtonchi2->Fill(TMath::Abs(pt),chi2);
					if(fMC)
					{
						if(primary&&pdgcode==-2212)
							fDCAXYZOpenforcleanAntiProtonsMCPrimary->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(!primary&&!primaryfake&&(pdgcode==-2212||pdgcodefake==-2212)&&uniqueID==kPDecay)	 
							fDCAXYZOpenforcleanAntiProtonsWD->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(!primary&&!primaryfake&&(pdgcode==-2212||pdgcodefake==-2212)&&uniqueID==kPHadronic)
							fDCAXYZOpenforcleanAntiProtonsHI->Fill(TMath::Abs(pt),dcaxy,dcaz);
						else if(fMC)
							fDCAXYZOpenforcleanAntiProtonsR->Fill(TMath::Abs(pt),dcaxy,dcaz);		 
					}
				}  
			}		
		}			
	}	
	fHistStandartMul->Fill(refmultiplicity);
	fHistMytrackMul->Fill(mynumberoftracks);
	
		
	// Post output data.
	Printf("Done..........\n");
	PostData(1,  flist);
	//Printf("....................Done!\n");
	//Printf("end of Exec");
}      

//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::Terminate(Option_t *) 
{
	//Terminate
	if(fCuts)
		fCuts->Dump();
	Printf("YCut=%f",fYCut);
	Printf("nsigma=%f",fsigmacut);
	Printf("DCA cut xy sigma =%f  zsigma=%f", fnsigmaxy, fnsigmaz);	
	Printf("ChargeCut=%f ", fchargeCut);
	Printf("DCAxy parameters %f  %f %f",fdcaxypar[0],fdcaxypar[1],fdcaxypar[2]);
	Printf("DCAz parameters %f  %f %f %f",fdcazpar[0],fdcazpar[1],fdcazpar[2],fdcazpar[3]);
	if(fTPCPIDCUT)
		fTPCPIDCUT->Dump();
	if(fMC)
		Printf("MC On\n");
	if(fCorrectSDD)
		Printf("correct SDD On\n");
  	if(fCorrectSSD)
		Printf("correct SSD On\n");
	if(fK0weight) 
	{
		Printf("weigth for pions");
		fK0weight->Print("All");
	}
	if(flambdaweight)
	{
		Printf("weigth for protons");
		flambdaweight->Print("All");
	}
	if(fAntilambdaweight)
	{
		Printf("weigth for antiprotons");
		fAntilambdaweight->Print("All");
	}
	Printf("Mul low %d Mul up %d",fLowMultiplicity, fUpMultiplicity);
	Printf("cent low %f cent up %f",fLowCentrality,fUpCentrality);
	if(fdovertexrescuts)
		Printf("Veretx resolution cut");
	Printf("end of Terminate");
}
//___________________________________________________
Float_t  AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::MyITSsignalusing4points(Double_t* tmpQESD) const
{
	//dE signal in case of 4 points in ITS
	Int_t indexes[4]={-1,-1,-1,-1};
	TMath::Sort(4,tmpQESD,indexes,0);
	return	0.5*(tmpQESD[indexes[0]]+tmpQESD[indexes[1]]);	
}
//________________________________________________________
Float_t  AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::MyITSsignalusing3points( Double_t* tmpQESD) const
{
	//dE signal in case of 3 points in ITS
	 Int_t indexes[3]={-1,-1,-1};
	TMath::Sort(3,tmpQESD,indexes,0);
	//cout<<tmpQESD[indexes[0]]<<" "<<tmpQESD[indexes[1]]<<" "<<tmpQESD[indexes[2]]<<endl;
	return	(tmpQESD[indexes[0]]+tmpQESD[indexes[1]]*0.5)/1.5;	
}			 
//____________________________________________________________________________________________________
  void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SetFunctionParam( Double_t * const par)
  {
  	fESDpid->GetITSResponse().SetBetheBlochParamsITSTPC(par);
  }
  //_____________________________________________________________________________________________________
  void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::CorrectSDD(Double_t *tmpQESD) const
  {
  //correction of SDD signal 
  	if(tmpQESD[0]>0.0)
	{
		tmpQESD[0]=tmpQESD[0]*3.34/5.43;
  		if(tmpQESD[0]<30.0)
			tmpQESD[0]=-1.0;
	}	
	if(tmpQESD[1]>0.0)
	{
		tmpQESD[1]=tmpQESD[1]*3.34/5.43;
  		if(tmpQESD[1]<30.0)
			tmpQESD[1]=-1.0;
	}	
  }
    //_____________________________________________________________________________________________________
  void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::CorrectSSD(Double_t *tmpQESD) const
  {
  	 //Correction of SSD signal
  	tmpQESD[2]=(85.0/77.0)*tmpQESD[2];
	tmpQESD[3]=(85.0/77.0)*tmpQESD[3];	
  }
  //_______________________________________________________________________________________________________
  Bool_t AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SelectOnImpPar(AliESDtrack* t)  const
  {
 	//
 	// cut on transverse impact parameter
 	//
 	Float_t d0z0[2],covd0z0[3];
 	t->GetImpactParameters(d0z0,covd0z0);
 	Float_t sigma= fdcaxypar[0]+fdcaxypar[1]/TMath::Power(t->Pt(),fdcaxypar[2]);
 	Float_t d0max = fnsigmaxy*sigma;
 	//
 	Float_t sigmaZ = fdcazpar[0]+fdcazpar[1]/TMath::Power(t->Pt(),fdcazpar[2]);
 	if (t->Pt() > 1) 
		sigmaZ = fdcazpar[3];
	 Float_t d0maxZ = fnsigmaz*sigmaZ;
 	//
 	if(TMath::Abs(d0z0[0]) < d0max && TMath::Abs(d0z0[1]) < d0maxZ) //error 
		return kTRUE;
 	return kFALSE;
}
//__________________________________________________________________________________________________
    Float_t AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::GetWeight(Int_t t,AliStack* const stack) const
    {
    
    	//Get weigth for pions protons and antiprotons
	if(stack->IsPhysicalPrimary(TMath::Abs(t)))
		return 1.0;	
	TParticle *particleMC = stack->Particle(TMath::Abs(t));
	Int_t pdgcodeMC = particleMC->GetPdgCode();
	if(TMath::Abs(pdgcodeMC)!=211&&TMath::Abs(pdgcodeMC)!=2212)
		return 1.0;
	if(!stack->IsPhysicalPrimary(TMath::Abs(particleMC->GetFirstMother())))
		return 1.0;	
	TParticle *particleMother=stack->Particle(TMath::Abs(particleMC->GetFirstMother()));
	Int_t pdgcodeMother = particleMother->GetPdgCode();
	Float_t motherpt=particleMother-> Pt();
	if(TMath::Abs(pdgcodeMC)==211&&pdgcodeMother==310&&fK0weight)
		return fK0weight->Eval(motherpt);
	else if (pdgcodeMother==3122&&flambdaweight)
		return flambdaweight->Eval(motherpt);
	else if(pdgcodeMother==-3122&&fAntilambdaweight)
		return fAntilambdaweight->Eval(motherpt);	
	return 1.0;	
    }
  //________________________________________________________________________________________________  
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SetDCA2010()
{
	//setting the DCA for 2010
	fdcaxypar[0]=0.0026;
	fdcaxypar[1]=0.005;
	fdcaxypar[2]=1.01;
	
	fdcazpar[0]=1000000.0;
	fdcazpar[1]=0.0;
	fdcazpar[2]=1.0;
	fdcazpar[3]=1000000.0;
}
//______________________________________________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SetCentralityCut(Float_t low, Float_t up)
{
	//centrality cut setter
	if((up>low)&&(!(low<0.0))&&(!(up>100.0)))
	{
		SetHImode();
		fLowCentrality=low;
		fUpCentrality=up;
	}
}
//_____________________________________________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SetMultiplicityCut(Int_t low, Int_t up)
{
	//mulyiplicty cut setter
	if((!(up>low))&&low>=0&&up>=0)
	{ 
		fLowMultiplicity=-1;
		fUpMultiplicity=-1;
	}
	else
	{
		fLowMultiplicity=low;
		fUpMultiplicity=up;
	}	
}	
