//Class to extract data to do ITS+TPC global Spectra
//Autor Marek Chojnacki
//emali Marek.Chojnacki@cern.ch
//Used on 2009 data
//last line of comments


#include "AliAnalysisChargedHadronSpectraITSTruncatedMeanTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"


//#include "AliESDtrack.h"

//#include "Riostream.h"
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

   class     AliMCEventHandler;
     class   Riostream;

using namespace std;
ClassImp(AliAnalysisChargedHadronSpectraITSTruncatedMeanTask)

//________________________________________________________________________
AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::AliAnalysisChargedHadronSpectraITSTruncatedMeanTask(const char *name) 
:AliAnalysisTaskSE(name),fESD(0),fCuts(0),fMC(0),
fLowMultiplicity(-1),fUpMultiplicity(-1),
fYCut(100.0),fsigmacut(3.0),fnsigmaxy(7.0),fnsigmaz(5.0),fchargeCut(0.0),
fCorrectSDD(0),fCorrectSSD(0),
 fK0weight(0),flambdaweight(0),fAntilambdaweight(0),
fHistStats(0),fHistPhiPtBeforeCuts(0),fHistPhiPtAfterCuts(0),fHistEtaPtBeforeCuts(0),fHistEtaPtAfterCuts(0),fHistDCABeforeCuts(0),fHistDCAAfterCuts(0),
fHistPminusTPCinPAfterCuts(0),fHistPminusTPCinPglobalAfterCuts(0),
fHistMydEPpositive(0),fHistMydETPCinPpositive(0),fHistMydETPCinPglobalpositive(0),
fHistMydEPnegative(0),fHistMydETPCinPnegative(0),fHistMydETPCinPglobalnegative(0),
fHistL3dEP(0),fHistL4dEP(0),fHistL5dEP(0),fHistL6dEP(0),fHistL3dETPCinP(0),
 fHistL4dETPCinP(0),fHistL5dETPCinP(0),fHistL6dETPCinP(0),fHistEtaPtPions(0), fHistEtaPtKaons(0),fHistEtaPtProtons(0),fHistwhichhasmin(0),
fHistminsignalforPionP(0),fHistminsignalforKaonP(0),fHistminsignalforProtonP(0),
fHistminsignalifPionP(0),fHistminsignalifKaonP(0),fHistminsignalifProtonP(0),fHistminsignalifAntiPionP(0),fHistminsignalifAntiKaonP(0),fHistminsignalifAntiProtonP(0),
fDCAXYZforcleanPions(0),fDCAXYZforcleanAntiPions(0),fDCAXYZforcleanProtons(0),fDCAXYZforcleanAntiProtons(0),
fDCAXYZOpenforcleanPions(0),fDCAXYZOpenforcleanAntiPions(0),fDCAXYZOpenforcleanProtons(0),fDCAXYZOpenforcleanAntiProtons(0),
fHistNtrackwithstandardcuts(0),fHistNtrackwithITSPIDcuts(0),
fHistSignalinTPCKaonforstandardcuts(0),fHistSignalinTPCKaonforITSPIDcuts(0),fHistSignalinTPCAntiKaonforstandardcuts(0),fHistSignalinTPCAntiKaonforITSPIDcuts(0),
fHistSignalinTPCProtonforstandardcuts(0),fHistSignalinTPCProtonforITSPIDcuts(0),fHistSignalinTPCAntiProtonforstandardcuts(0),fHistSignalinTPCAntiProtonforITSPIDcuts(0),
fHistStandartMul(0),fHistMytrackMul(0),
fHistEtaPtPionsMC(0),fHistEtaPtKaonsMC(0),fHistEtaPtProtonsMC(0),
 fHistEtaPtPionsMCDET(0),fHistEtaPtKaonsMCDET(0),fHistEtaPtProtonsMCDET(0),
 fHistEtaPtPionsCon(0),fHistEtaPtKaonsCon(0),fHistEtaPtProtonsCon(0),
fHistEtaPtPionsConPID(0),fHistEtaPtKaonsConPID(0),fHistEtaPtProtonsConPID(0),
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
fPionNTPCClusters(0),fAntiPionNTPCClusters(0),
fTPCPIDCUT(0), fESDpid(0),fPrimaryElectronsMother(0),
flist(0)
{
	//Constructor
	 fESDpid=new AliESDpid();
	 fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,2.63394e+01,5.04114e-11,2.12543e+00,4.88663e+00);
	 
	 
	
	 
	 for(int i=0;i<2;i++){ flinearpar[i]=0.0;}
	for(int i=0;i<5;i++){ fpar[i]=0.0;}
 	Printf("end of AliAnalysisChargedHadronSpectraITSTruncatedMeanTask");
	 DefineOutput(1, TList::Class());
}


//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::UserCreateOutputObjects()
{
	//UserCreateOutputObject
	Printf("AliAnalysisChargedHadronSpectraITSTruncatedMeanTask UserCreateOutputObjects");
	  flist=new TList();
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
	
	fHistEtaPtPions=new  TH2F("HistEtaPtPions",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtPions);
	fHistEtaPtKaons=new  TH2F("HistEtaPtKaons",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtKaons);
	fHistEtaPtProtons=new  TH2F("HistEtaPtProtons",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtProtons);

	fHistwhichhasmin=new TH1F("Histwhichhasmin","Histwhichhasmin",4,-0.5,3.5);
	fHistwhichhasmin->GetXaxis()->SetBinLabel(1,"SDD1");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(2,"SDD2");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(3,"SSD1");
	fHistwhichhasmin->GetXaxis()->SetBinLabel(4,"SSD2");
	flist->Add(fHistwhichhasmin);
	

	
	fHistminsignalforPionP=new TH2F("HistminsignalforPionP",";P[GeV/c];dE[in 300#mum]",kPtBins,binsPtDummy,ny,starty,ny*jump+starty);
	flist->Add(fHistminsignalforPionP);
	fHistminsignalforKaonP=new TH2F("HistminsignalforKaonP",";P[GeV/c];dE[in 300#mum]",kPtBins,binsPtDummy,ny,starty,ny*jump+starty);
	flist->Add(fHistminsignalforKaonP);
	fHistminsignalforProtonP=new TH2F("HistminsignalforProtonP",";P[GeV/c];dE[in 300#mum]",kPtBins,binsPtDummy,ny,starty,ny*jump+starty);
	flist->Add(fHistminsignalforProtonP);
	 
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
	
	fHistStandartMul=new TH1F("fHistStandartMul",";Ntracks;counts",300,0,300);
	flist->Add(fHistStandartMul);
	fHistMytrackMul=new TH1F("fHistMytrackMul",";Ntracks;counts",300,0,300);
	flist->Add(fHistMytrackMul);
	
	if(!fMC)
	{
		
		Printf("end of CreateOutputObjects no MC");
		PostData(1,flist);
		return;
	}
	
	
	fHistEtaPtPionsMC=new  TH2F("HistEtaPtPionsMC",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtPionsMC);
	fHistEtaPtKaonsMC=new  TH2F("HistEtaPtKaonsMC",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtKaonsMC);
	fHistEtaPtProtonsMC=new  TH2F("HistEtaPtProtonsMC",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtProtonsMC);
	
	fHistEtaPtPionsMCDET=new  TH2F("HistEtaPtPionsMCDET",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtPionsMCDET);
	fHistEtaPtKaonsMCDET=new  TH2F("HistEtaPtKaonsMCDET",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtKaonsMCDET);
	fHistEtaPtProtonsMCDET=new  TH2F("HistEtaPtProtonsMCDET",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtProtonsMCDET);
	
	fHistEtaPtPionsCon=new  TH2F("HistEtaPtPionsCon",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtPionsCon);
	fHistEtaPtKaonsCon=new  TH2F("HistEtaPtKaonsCon",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtKaonsCon);
	fHistEtaPtProtonsCon=new  TH2F("HistEtaPtProtonsCon",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtProtonsCon);
	
	fHistEtaPtPionsConPID=new  TH2F("HistEtaPtPionsConPID",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtPionsConPID);
	fHistEtaPtKaonsConPID=new  TH2F("HistEtaPtKaonsConPID",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtKaonsConPID);
	fHistEtaPtProtonsConPID=new  TH2F("HistEtaPtProtonsConPID",";#eta;pt[GeV/c]",netabins,-1.0*etamax,etamax,nptbins,-1.0*ptmax,ptmax);//eta pt
	flist->Add(fHistEtaPtProtonsConPID);

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
	 fESD = dynamic_cast<AliESDEvent*> (InputEvent());
	if (!fESD) 
	{
		Printf("ERROR: fESD not available");
		return;
	}
	
	Float_t refmultiplicity=AliESDtrackCuts::GetReferenceMultiplicity(fESD,1);
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
	AliStack* stack=0x0;
	Double_t mcZvertex=0.0;
	if(fMC)
	{
		AliMCEvent* mcEvent  = (AliMCEvent*) MCEvent();
		Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
		stack = mcEvent->Stack();
		mcZvertex=mcEvent->GetPrimaryVertex()->GetZ();
	}	
	
	if(stack)//Looping over MC information of all events
	{
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
			if (TMath::Abs(particleMC->Y())>fYCut)
				continue;
			if (particleMC->Pt()>2.0)
				continue;		
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
			
			if(TMath::Abs(mcZvertex)>10.0)
				continue;
			
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
	}	
	fHistStats->Fill(0);
	//Event selection 
	//if( ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()==0)	
	  UInt_t isSelected = 0;
	 if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))
	  	isSelected=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
	Printf("Mask_selection %d", isSelected);  
        if(!(isSelected&AliVEvent::kMB))
        {
		Printf("Not Physics event.........\n");
		PostData(1,  flist);
		Printf("end of Exec");
		return;
	}	
	fHistStats->Fill(1);
	//Good vertex	
	const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
	if(vertex->GetNContributors()<1) 
	{
  		// SPD vertex
  		vertex = fESD->GetPrimaryVertexSPD();
  		if(vertex->GetNContributors()<1) 
		{
   			Printf("No good  Vertex.........\n");
			PostData(1,  flist);
			Printf("end of Exec");
			return;
  		}
	}
	if(stack)
	{
		for (int imc=0;imc<stack->GetNtrack();imc++)
		{
			if(!(stack->IsPhysicalPrimary(imc)))
				continue;
			TParticle *particleMC = stack->Particle(imc);
			if(!particleMC)
				continue;
			Int_t pdgcodeMC = particleMC->GetPdgCode();
			if (TMath::Abs(particleMC->Y())>fYCut)
				continue;
			if (particleMC->Pt()>2.0)
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
		}
	}
	fHistStats->Fill(2);
	if(TMath::Abs(vertex ->GetZ())>10.0)
	{
		Printf("No good  Z of Vertex.........\n");
		PostData(1,  flist);
		Printf("end of Exec");
		return;
	}
	fHistStats->Fill(3);
	
	if(fCuts==0)
	{
		Printf("No CUTS Defined.........\n");
		PostData(1,  flist);
		Printf("end of Exec");
		return;
	}
	
	Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
	 Int_t nTracks=fESD->GetNumberOfTracks();
	
	Int_t mynumberoftracks=0;
	 AliESDtrack *trackESD=0;
	 
	const Float_t pionmass=0.13957;
	const Float_t kaonmass=0.493677;
	const Float_t protonmass=0.938272;

 	for(int tr1=0;tr1<nTracks;tr1++)
  	{	
		
		trackESD=fESD->GetTrack(tr1);	
		//fHistStats->Fill(2);
		
		Double_t pt=trackESD->Pt()*trackESD->GetSign();
		Double_t p=trackESD->P();
		Double_t eta=trackESD->Eta();
		Double_t phi=trackESD->Phi();
		Float_t dcaxy, dcaz;
		Double_t pz=trackESD->Pz();
		 UShort_t nTPCclusters=trackESD->GetTPCNcls();
		fHistPhiPtBeforeCuts->Fill(phi,pt);//phi pt
		fHistEtaPtBeforeCuts->Fill(eta,pt);
		fHistDCABeforeCuts->Fill(dcaxy,dcaz);		
		//standart cuts 	
		if(fCuts->AcceptTrack(trackESD)==kFALSE)
			continue;
		//Tpc pid cut for debug	
		Double_t pinTPC=trackESD->GetTPCInnerParam()->GetP();//momentum in primary vertex taken from TPC tracking
		Double_t pinTPCglobal=trackESD->GetInnerParam()->GetP();//momentum at the inner  wall of the TPC taken from global tracking
		Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
		Float_t sigProton     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
		Double_t tpcSignal =trackESD ->GetTPCsignal();
		trackESD->GetImpactParameters(dcaxy,dcaz);
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
			CorrectSDD(tmpQESD	);	
		if(fCorrectSSD&&fMC)
			CorrectSSD(tmpQESD	);	
		Int_t nSSDSDD=0;
		Int_t nSSDSDDall=0;
		
		//Y assumtion
		Float_t yforpion=0.5*TMath::Log((TMath::Sqrt(pionmass*pionmass+p*p)+pz)/(TMath::Sqrt(pionmass*pionmass+p*p)-pz));
		Float_t yforkaon=0.5*TMath::Log((TMath::Sqrt(kaonmass*kaonmass+p*p)+pz)/(TMath::Sqrt(kaonmass*kaonmass+p*p)-pz));
		Float_t yforproton=0.5*TMath::Log((TMath::Sqrt(protonmass*protonmass+p*p)+pz)/(TMath::Sqrt(protonmass*protonmass+p*p)-pz));
		
		if(TMath::Abs(yforpion)>fYCut&&TMath::Abs(yforkaon)>fYCut&&TMath::Abs(yforproton)>fYCut) //go trought one y cut
			continue;
		
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
		if(TMath::Abs(yforpion)<fYCut&&cutDCA) 
		{
			fHistNtrackwithITSPIDcuts->Fill(TMath::Abs(pt),0);
			if(pt>0.0)
				 fPionNTPCClusters->Fill(pt,nTPCclusters);
			else
				fAntiPionNTPCClusters->Fill(TMath::Abs(pt),nTPCclusters);	
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
			Double_t tmp2QESD[3];
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
			fHistwhichhasmin->Fill(0);
		if(whichLmin==1)
			fHistwhichhasmin->Fill(1);
		if(whichLmin==2)
			fHistwhichhasmin->Fill(2);
		if(whichLmin==3)
			fHistwhichhasmin->Fill(3);
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
			
		
		Float_t itspidsignalforpions=TMath::Log(signaltouse)-TMath::Log(BBparametrization(p/pionmass,fpar));
		Float_t itspidsignalforkaons=TMath::Log(signaltouse)-TMath::Log(BBparametrization(p/kaonmass,fpar));
		Float_t itspidsignalforprotons=TMath::Log(signaltouse)-TMath::Log(BBparametrization(p/protonmass,fpar));
		if(cutDCA)
			mynumberoftracks++;
		
		Int_t pa= TypeofParticle(pinTPC,signaltouse);			
		Int_t label=-1;
		if(fMC)
			label=trackESD->GetLabel();
		Int_t pdgcode=0;
		Int_t primary=0;
		Double_t chargeMC=1.0;
		Float_t etaMC=10.0;
		Float_t ptMC=10.0;
		Int_t   uniqueID=-1;
		Int_t pdgcodefake=0;
		Int_t primaryfake=0;
		
		
		TParticle *particle2=0x0;
			
		if(label>0&&stack&&fMC)
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
		if(pa==0&&TMath::Abs(yforpion)<=fYCut&&cutDCA)
		{
			fHistEtaPtPions->Fill(eta,pt);
			fHistminsignalforPionP->Fill(pinTPC,signaltouse);
			if(stack&&fMC)
			{
				if(!primary)
					fHistEtaPtPionsCon->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)!=211)	
					fHistEtaPtPionsConPID->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)==211)	
					fHistEtaPtPionsMCDET->Fill(eta,pt);
			}
		}	
		if(pa==1&&TMath::Abs(yforkaon)<=fYCut&&cutDCA)
		{
			fHistEtaPtKaons->Fill(eta,pt);
			fHistminsignalforKaonP->Fill(pinTPC,signaltouse);
			if(stack&&fMC)
			{	
				if(!primary)
					fHistEtaPtKaonsCon->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)!=321)		
					fHistEtaPtKaonsConPID->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)==321)	
					fHistEtaPtKaonsMCDET->Fill(eta,pt);	
			}
		}	
		if(pa==2&&TMath::Abs(yforproton)<=fYCut&&cutDCA)
		{
			fHistEtaPtProtons->Fill(eta,pt);
			fHistminsignalforProtonP->Fill(pinTPC,signaltouse);
			if(stack&&fMC)
			{
				if(!primary)
					fHistEtaPtProtonsCon->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)!=2212)	
					fHistEtaPtProtonsConPID->Fill(etaMC,chargeMC*ptMC);
				if(primary&&TMath::Abs(pdgcode)==2212)	
					fHistEtaPtProtonsMCDET->Fill(eta,pt);		
			}
		}
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
						if(primary&&pdgcode==-11)
							fPrimaryElectronsMother->Fill(stack->Particle(particle2->GetFirstMother())->GetPdgCode());
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
						if(primary&&pdgcode==11)
							fPrimaryElectronsMother->Fill(stack->Particle(particle2->GetFirstMother())->GetPdgCode());
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
	
		
	if(stack)//Looping over MC information of all events which passed ESD cuts
	{
		for (int imc=0;imc<stack->GetNtrack();imc++)
		{
			if(!(stack->IsPhysicalPrimary(imc)))
				continue;
			TParticle *particleMC = stack->Particle(imc);
			if(!particleMC)
				continue;
			Int_t pdgcodeMC = particleMC->GetPdgCode();
			if (TMath::Abs(particleMC->Y())>fYCut)
				continue;
			if (particleMC->Pt()>2.0)
				continue;		
			if(pdgcodeMC==211)	
			{
				fHistEtaPtPionsMC->Fill(particleMC->Eta(),particleMC->Pt());
				fHistminsignalifPionPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-211)
			{	
				fHistEtaPtPionsMC->Fill(particleMC->Eta(),-1.0*particleMC->Pt());
				fHistminsignalifAntiPionPMCPrimary->Fill(particleMC->Pt());
			}		
			if(pdgcodeMC==321)
			{	
				fHistEtaPtKaonsMC->Fill(particleMC->Eta(),particleMC->Pt());
				fHistminsignalifKaonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-321)
			{	
				fHistEtaPtKaonsMC->Fill(particleMC->Eta(),-1.0*particleMC->Pt());
				fHistminsignalifAntiKaonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==2212)
			{	
				fHistEtaPtProtonsMC->Fill(particleMC->Eta(),particleMC->Pt());
				fHistminsignalifProtonPMCPrimary->Fill(particleMC->Pt());
			}	
			if(pdgcodeMC==-2212)
			{	
				fHistEtaPtProtonsMC->Fill(particleMC->Eta(),-1.0*particleMC->Pt());
				fHistminsignalifAntiProtonPMCPrimary->Fill(particleMC->Pt());
			}	
			
		}
	}	
	// Post output data.
	Printf("Done..........\n");
	PostData(1,  flist);
	Printf("....................Done!\n");
	Printf("end of Exec");
}      

//________________________________________________________________________
void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::Terminate(Option_t *) 
{
	//Terminate
	if(fCuts)
		fCuts->Dump();
	Printf("BB parameters %f  %f %f %f %f",fpar[0],fpar[1],fpar[2],fpar[3],fpar[4]);
	Printf("linear parameters a=%f  b=%f ",flinearpar[0],flinearpar[1]);
	Printf("YCut=%f",fYCut);
	Printf("nsigma=%f",fsigmacut);
	Printf("DCA cut xy sigma =%f  zsigma=%f", fnsigmaxy, fnsigmaz);	
	Printf("ChargeCut=%f ", fchargeCut);
	
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
//______________________________________________________________
  Int_t AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::TypeofParticle(Float_t mom,Float_t signal) const
  {
  	//nsigma cut show which type
  	const Float_t pionmass=0.13957;
	 const Float_t kaonmass=0.493677;
	 const Float_t protonmass=0.938272;
  	if(mom<0.15||mom>1.1)
		return -1;
	Float_t bg[3]={mom/pionmass,mom/kaonmass,mom/protonmass};
	Float_t nsigma[3]={4.0,4.0,4.0};
	for(int i=0;i<3;i++)
	{
		Float_t peak=BBparametrization(bg[i],fpar);
		Float_t rms= flinearpar[0]*peak+flinearpar[1];
		nsigma[i]=TMath::Abs((peak-signal)/rms);		
	}
	Int_t pa=TMath::LocMin(3,nsigma);	
	if(nsigma[pa]<fsigmacut)
	{		
		return pa;
	}		
	else
		return -1;
}			 
  //______________________________________________________________________________________
   Float_t AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::BBparametrization(Float_t x,const Float_t* par) const
{
	//BBparametrization
	
	Float_t bg=x;
//	cout<<"bg  "<<x[0]<<endl;
//	cout<<par[0]<<"  "<<par[1]<<"  "<<par[2]<<"  "<<par[3]<<"  "<<par[4]<<"  "<<endl;
	Float_t beta = bg/TMath::Sqrt(1.+ bg*bg);
	Float_t gamma=bg/beta;
	
	Float_t eff=1.0;
	if(bg<par[2])
		 eff=(bg-par[3])*(bg-par[3])+par[4];
	else
		eff=(par[2]-par[3])*(par[2]-par[3])+par[4];
	return (par[1]+2.0*TMath::Log(gamma)-beta*beta)*(par[0]/(beta*beta))*eff;	

}
//____________________________________________________________________________________________________
  void AliAnalysisChargedHadronSpectraITSTruncatedMeanTask::SetFunctionParam( Float_t * const par)
  {
  	//setter for BB parameters
  	for(int i=0;i<5;i++)
		fpar[i]=par[i];
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
 	Float_t sigma= 0.0050+0.0060/TMath::Power(t->Pt(),0.9);
 	Float_t d0max = fnsigmaxy*sigma;
 	//
 	Float_t sigmaZ = 0.0146+0.0070/TMath::Power(t->Pt(),1.114758);
 	if (t->Pt() > 1) 
		sigmaZ = 0.0216;
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
    
