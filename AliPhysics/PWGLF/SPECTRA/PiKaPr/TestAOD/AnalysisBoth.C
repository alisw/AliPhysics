#if !defined (__CINT__) || (defined(__MAKECINT__))

#include "AliSpectraBothHistoManager.h"
#include "AliSpectraBothEventCuts.h"
#include "AliSpectraBothTrackCuts.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFormula.h"	
#include "TMath.h"
#include "TList.h"
#include "TCanvas.h"
#include "TFractionFitter.h"
#include  <Riostream.h>
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <TROOT.h>
#include <QAPlotsBoth.C>
#include <TDatabasePDG.h>
#include <TDirectoryFile.h>
#include <TLatex.h>
#include <TGraph.h>
#endif

#if defined(__CINT__)
TString Particle[]={"Pion","Kaon","Proton"};
#endif

class AliSpectraBothHistoManager;
class AliSpectraBothEventCuts; 
class AliSpectraBothTrackCuts;
TString Charge[]={"Pos","Neg"};
TString Sign[]={"Plus","Minus"};
//TString Particle[]={"Pion","Kaon","Proton"};
TString symboles[]={"#pi^{+}","K^{+}","p","pi^{-}","K^{-}","#bar{p}"}; 
AliSpectraBothHistoManager* managerdata=0x0;
AliSpectraBothEventCuts* ecutsdata=0x0; 
AliSpectraBothTrackCuts* tcutsdata=0x0;
	
AliSpectraBothHistoManager* managermc=0x0;
AliSpectraBothEventCuts* ecutsmc=0x0; 
AliSpectraBothTrackCuts* tcutsmc=0x0;

Float_t TOFMatchingScalling[2]={-1,-1};
Int_t Color[3]={1,2,4};
Int_t Marker[6]={20,21,22,24,25,26};
Double_t Range[3]={0.3,0.3,0.5}; // LowPt range for pi k p


Double_t FitRange[2]={-2.0,2.0};
Double_t CutRange[2]={-3.0,3.0};
Double_t minptformaterial[6]={0.0,0.2,0.0,0.0,0.2,0.0};
Double_t maxptformaterial[6]={0.0,0.6,1.3,0.0,0.6,0.0};
Double_t minptforWD[6]={0.2,100.0,0.3,0.2,100.0,0.3};
Double_t maxptforWD[6]={1.5,-100.0,2.0,1.5,-100.0,2.0};
Double_t minRanges[3]={0.3,0.3,0.45};
Double_t maxRanges[3]={1.5,1.2,2.2};
Double_t TOFPIDsignalmatching[]={-1.0,-1.0,-1.0};
Double_t fMaxContaminationPIDMC=0.2;
TString filenames[]={"eff.root","pid.root","sec.root"};

enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
enum {
 kdodca=0x1, //dca fits are made 
 kgeantflukaKaon=0x2,// geant fluka correction is used for kaons 
 kgeantflukaProton=0x4, // geant fluka correction is used for protons and antiprotons
 knormalizationtoeventspassingPhySel=0x8,// spectra are divided by number of events passing physic selection   
 kveretxcorrectionandbadchunkscorr=0x10, // correction for difference in z vertex distribution in data and MC and correction for bad chunks is applied
 kmcisusedasdata=0x20, // the result of the looping over MC is used as data input 
 kdonotusedcacuts=0x40, // allows to use the constant dca cut for all pt bins not the pt dependet defined in stardrad track cuts 2011
 kuseprimaryPIDcont=0x80, //pid contamination is calculated using only primiary particle in this case K should use dca fits 
 knormalizationwithbin0integralsdata=0x100, // the normalization factor is calcualte using integral over z vertex distributions (in this case reconstructed vertex disitrbution uses z vertex for data) 
 knormalizationwithbin0integralsMC=0x200, //in this case reconstructed vertex disitrbution uses z vertex for data, those to options will be use only if knormalizationtoeventspassingPhySel is not set
 kuserangeonfigfile=0x400, // use of config file for dca fit settings
 kskipconcutonspectra=0x800, //do not use conPID<02 cut  useful for syst. studies
 kuseTOFmatchingcorrection=0x1000, // if set tof matching correction is applied.
 kuseTOFcorrforPIDsignalmatching=0x2000, // rescale the for spectra by the factor given in config files	
 kuseeffcorrectionfromfile=0x4000, //use the efficiency from the file specfied in config file			
 kusePIDcontaminatiofromfile=0x8000, //use the PID contamination from the file specfied in config file
 kuseseccontaminatiofromfile=0x10000, //use the secondary contamination from the file specfied in config file
 kusespecialbinninginDCAfits=0x20000, //use constum binning in the dca fit 
 kusePIDfits=0x40000 //raw yields from fits
 							
};	

Bool_t OpenFile(TString dirname, TString outputname, Bool_t mcflag,Bool_t mcasdata=false);
 void GetMCTruth(TH1F** MCTruth);
void GetPtHistFromPtDCAhisto(TString hnamein, TString hnameout, AliSpectraBothHistoManager* hman,TH1F** histo,TFormula* dcacutxy);
void CleanHisto(TH1F* h, Float_t minV, Float_t maxV,TH1* contpid=0x0);
void DCACorrectionMarek(AliSpectraBothHistoManager* hman_data, AliSpectraBothHistoManager* hman_mc,TFormula* fun,TFile *fout,TH1F** hcon,TH1F** hconWD,TH1F** hconMat,TH1F** hprimary,Bool_t binning=false);
void RecomputeErrors(TH1* h);
void SetBintoOne(TH1* h);
void GetCorrectedSpectra(TH1F* corr,TH1F* raw,TH1F* eff, TH1F* con);
void GetCorrectedSpectraLeonardo(TH1F* spectra,TH1F* correction, TH1F* hprimaryData,TH1F* hprimaryMC);
void GFCorrection(TH1F **Spectra,Float_t tofpt,UInt_t options);
void MatchingTOFEff(TH1F** Spectra, TList* list=0x0);
void Divideby2pipt(TH1* hist);
void SubHistWithFullCorr(TH1F* h1, TH1F* h2, Float_t factor1=1.0, Float_t factor2=1.0);
void TOFMatchingForNch(TH1* h);
void TOFPIDsignalmatchingApply(TH1* h, Float_t factor);
void CalculateDoubleCounts(TH1* doubleconunts,TH1F** rawspectra,Int_t ipar, Bool_t dataflag);
void CopyCorrectionFromFile(TString filename,TString correctionname,TH1F** corrtab);
void RawYieldFromFits(AliSpectraBothHistoManager* hman,TH1F** histo,AliSpectraBothTrackCuts* cuts,TList* lfits,TH1F** primaryfractionfromfit);
void RawYieldFromFitsForParticleType(TH2F* nsigTOF,TH2F* nsigTPC,TH1F* histo,Int_t icharge, Int_t ipart,Float_t pttof,TList* lfits,TH1F* primaryfractionfromfit,TH2F* nsigtrueTOF,TH2F* nsigtrueTPC);
void Rescallesecondarycontaimation(TH1F* contfit, TH1F* contWDfit , TH1F* contMatfit ,TH1F* primaryfit, TH1F* primaryfractionfromfit);

TH1F* GetOneHistFromPtDCAhisto(TString name,TString hnameout,AliSpectraBothHistoManager* hman,TFormula* dcacutxy);
TF1* TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge);
Double_t TrackingPtGeantFlukaCorrectionNull(Double_t pTmc);
Double_t TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc);
Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc);
TF1* TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge);
Double_t MatchingPtGeantFlukaCorrectionNull(Double_t pTmc);
Double_t MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc);
Double_t MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc);
Double_t eta2y(Double_t pt, Double_t mass, Double_t eta);
TH1* GetSumAllCh(TH1F** spectra, Double_t* mass,Double_t etacut);
Short_t DCAfitsettings (Float_t pt, Int_t type);
Float_t Normaliztionwithbin0integrals(UInt_t options);
Bool_t ReadConfigFile(TString configfile);
TH1F* ReBinDCAHisto(TH1* h);



void AnalysisBoth (UInt_t options,TString outdate, TString outnamedata, TString outnamemc="",TString configfile="",TString customoutfilename="")
{
	cout<<"A"<<endl;
	gStyle->SetOptStat(0);	
	TH1::AddDirectory(kFALSE);
	/*#if defined(__CINT__)
	gSystem->Load("libCore");
	gSystem->Load("libPhysics");
	gSystem->Load("libTree");
	gSystem->Load("libMatrix");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libOADB");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libTender");
	gSystem->Load("libCORRFW");
	gSystem->Load("libPWGTools");
	gSystem->Load("libPWGLFspectra");
  	
  	gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/PiKaPr/TestAOD/QAPlotsBoth.C");
	#endif*/

	
	cout<<"A"<<endl;
	Double_t mass[3];
	mass[0]   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
	mass[1]   = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
	mass[2] = TDatabasePDG::Instance()->GetParticle("proton")->Mass();

	TFormula* dcacutxy=0x0;
	if(!(options&kdonotusedcacuts))
	{
	
		AliESDtrackCuts* esdtrackcuts= AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE);
		TString formulastring(esdtrackcuts->GetMaxDCAToVertexXYPtDep());	
		formulastring.ReplaceAll("pt","x");
		dcacutxy=new TFormula("dcacutxy",formulastring.Data());
	}
	if(options&kuserangeonfigfile)
		if(!ReadConfigFile(configfile))
			return;		
	TList* lout=new TList();
	

	TString indirname=Form("/output/train%s",outdate.Data());
	//TString indirname("/output/train24072012");
	if(outnamemc.Length()==0)
	outnamemc=outnamedata;
	cout<<indirname.Data()<<" "<<outnamemc.Data()<<endl;
	// Do the job 


	OpenFile(indirname,outnamemc,true);
	OpenFile(indirname,outnamedata,false,((Bool_t)(options&kmcisusedasdata)));
	if(!managermc||!managerdata)
	{
		cout<<managermc<<" "<<managerdata<<endl;
		return;	
	}
	TH1F* rawspectradata[6];
	TH1F* rawspectramc[6];
	TH1F* MCTruth[6];
	TH1F* eff[6];
	TH1F* contallMC[6];
	TH1F* contPID[6];
	TH1F* contWD[6];
	TH1F* contMat[6];
	TH1F* confinal[6];
	TH1F* contPIDpri[6];
	TH1F* contSecMC[6];
	
	TH1F* contfit[12];
	TH1F* contWDfit[12];
	TH1F* contMatfit[12];
	TH1F* primaryfit[12];

	
	
	TH1F* spectra[6];
	TH1F* spectraLeonardo[6];
	
	TH1F* corrLeonardo[6];
	TH1F* primaryfractionfromfitdata[6];
	TH1F* primaryfractionfromfitmc[6];



        TH1F* doubleconuntsdata[3]; 
	 TH1F* doubleconuntsMC[3]; 
	//GetSpectra(managerdata,rawspectradata,true);
	//GetSpectra(managermc,rawspectramc,true,true);
	TList* lfits=new TList();
	GetPtHistFromPtDCAhisto("hHistPtRecTrue","conPID",managermc,contPID,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigmaSecondaryWeakDecay","conWD",managermc,contWD,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigmaSecondaryMaterial","conMat",managermc,contMat,dcacutxy);
	GetPtHistFromPtDCAhisto("hHistPtRecSigmaPrimary","conPIDprimary",managermc,contPIDpri,dcacutxy);
	if((options&kusePIDfits)==kusePIDfits)
	{
		GetPtHistFromPtDCAhisto("hHistPtRecPrimary","eff",managermc,eff,dcacutxy);
		TString tmpnames[]={"PionPlus","KaonPlus","ProtonPlus","PionMinus","KaonMinus","ProtonMinus"};
		for (int i=0;i<6;i++)
		{
			rawspectramc[i]=(TH1F*)eff[i]->Clone(Form("SpectraMC%s",tmpnames[i].Data()));
			rawspectradata[i]=(TH1F*)eff[i]->Clone(Form("SpectraDATA%s",tmpnames[i].Data()));
			rawspectramc[i]->SetTitle(Form("SpectraMC%s",tmpnames[i].Data()));
			rawspectradata[i]->SetTitle(Form("SpectraDATA%s",tmpnames[i].Data()));
			rawspectramc[i]->Reset();
			rawspectradata[i]->Reset();
			
			primaryfractionfromfitdata[i]=(TH1F*)eff[i]->Clone(Form("primaryfractionfromfitdata%s",tmpnames[i].Data()));
			primaryfractionfromfitmc[i]=(TH1F*)eff[i]->Clone(Form("primaryfractionfromfitmc%s",tmpnames[i].Data()));
			primaryfractionfromfitdata[i]->SetTitle(Form("primaryfractionfromfitdata%s",tmpnames[i].Data()));
			primaryfractionfromfitmc[i]->SetTitle(Form("primaryfractionfromfitmc%s",tmpnames[i].Data()));
			primaryfractionfromfitdata[i]->Reset();
			primaryfractionfromfitmc[i]->Reset();

			
		
		}	
		RawYieldFromFits(managermc,rawspectramc,tcutsmc,lfits,primaryfractionfromfitmc);
		RawYieldFromFits(managerdata,rawspectradata,tcutsdata,lfits,primaryfractionfromfitdata);

		

	}
	else
	{
		GetPtHistFromPtDCAhisto("hHistPtRecTruePrimary","eff",managermc,eff,dcacutxy);
		GetPtHistFromPtDCAhisto("hHistPtRecSigma","SpectraMC",managermc,rawspectramc,dcacutxy);
		GetPtHistFromPtDCAhisto("hHistPtRecSigma","SpectraDATA",managerdata,rawspectradata,dcacutxy);

	}

	
	Double_t neventsmcall = 1 ;  //if loop over MC is done after or befor events cuts this will be changed 
	Double_t neventsdata =  1;
	Double_t neventsmc =  1;

	//Normaliztion of MCtruth depends if the loop was done after of before ESD event cuts.
	//In currect code this cannot be check on the level of macro.
	//If the loop was done before MC should be done to all processed events (NumberOfProcessedEvents())
	//If loop was done after MC should be normalized to all accepted events (NumberOfEvents()) 
	// The option one will be alaways use.
	
	neventsmcall= ecutsmc->NumberOfProcessedEvents();
	if(managermc->GetEventStatHist())
	{
		if(managermc->GetEventStatHist()->GetBinContent(1)!=neventsmcall)
			cout<<"merging problem MC"<<endl;
		neventsmcall=managermc->GetEventStatHist()->GetBinContent(1);	
	}

		


	if(options&knormalizationtoeventspassingPhySel)
	{
		//neventsmcall= ecutsmc->NumberOfProcessedEvents();
		 neventsdata=ecutsdata->NumberOfPhysSelEvents();
		 neventsmc=ecutsmc->NumberOfPhysSelEvents();
		 if(managerdata->GetEventStatHist())
		 {
			if(managerdata->GetEventStatHist()->GetBinContent(2)!=neventsdata)
				cout<<"merging problem data"<<endl;
			neventsdata=managerdata->GetEventStatHist()->GetBinContent(2);
		 }	
		 if(managermc->GetEventStatHist())
		 {
			if(managermc->GetEventStatHist()->GetBinContent(2)!=neventsmc)
				cout<<"merging problem MC"<<endl;
			neventsmc=managermc->GetEventStatHist()->GetBinContent(2);	
		 }
		

	}
	else if ((options&knormalizationwithbin0integralsdata)||(options&knormalizationwithbin0integralsMC))
	{
		neventsdata=Normaliztionwithbin0integrals(options);
		neventsmc=ecutsmc->NumberOfPhysSelEvents();
		if(managermc->GetEventStatHist())
		{
			if(managermc->GetEventStatHist()->GetBinContent(2)!=neventsmc)
				cout<<"merging problem MC"<<endl;
			neventsmc=managermc->GetEventStatHist()->GetBinContent(2);	
		 }
		



	}
	else
	{
		neventsdata=ecutsdata->NumberOfEvents(); //number of accepted events
		 neventsmc=ecutsmc->NumberOfEvents();
		neventsmcall= ecutsmc->NumberOfEvents();
		if(managerdata->GetEventStatHist())
		{
			if(managerdata->GetEventStatHist()->GetBinContent(3)!=neventsdata)
				cout<<"merging problem data"<<endl;
			neventsdata=managerdata->GetEventStatHist()->GetBinContent(3);
		}	
		if(managermc->GetEventStatHist())
		{
			if(managermc->GetEventStatHist()->GetBinContent(3)!=neventsmc)
				cout<<"merging problem MC"<<endl;
			neventsmc=managermc->GetEventStatHist()->GetBinContent(3);
			neventsmcall=managermc->GetEventStatHist()->GetBinContent(3);
		}





	}
	GetMCTruth(MCTruth);
	cout<<neventsdata<<" Events"<<endl;
        cout<< neventsmc<<" Events "<<endl;
       if(neventsdata<1)
       {
               cout<<"No events DATA"<<endl;
               return;
       }
	if(neventsmc<1&&(((options&kusePIDcontaminatiofromfile)!=kusePIDcontaminatiofromfile)||((options&kuseseccontaminatiofromfile)!=kuseseccontaminatiofromfile)||((options&kuseseccontaminatiofromfile)!=kuseseccontaminatiofromfile)))
       {
               cout<<"No events MC"<<endl;
               return;
       }
	else if(neventsmc<1)
	{
		neventsmc=1; //dumpy normalization
		neventsmcall=1; //dumpy normalization

	}

	cout<<"A "<<neventsmc<<endl;
	TH1F* allgen=(TH1F*)((TH1F*)managermc->GetPtHistogram1D("hHistPtGen",1,1))->Clone();
	allgen->SetName("AllGen");
	TH1F* allrecMC=GetOneHistFromPtDCAhisto("hHistPtRec","rawallMC",managermc,dcacutxy);
	TH1F* alleff=GetOneHistFromPtDCAhisto("hHistPtRecPrimary","effall",managermc,dcacutxy);
	TH1F* allrecdata=GetOneHistFromPtDCAhisto("hHistPtRec","rawalldata",managerdata,dcacutxy);
	
        TH1F* muons[4]; 
	muons[0]=GetOneHistFromPtDCAhisto("hHistPtRecTrueMuonPlus","MuonPlusAll",managermc,dcacutxy);
	muons[1]=GetOneHistFromPtDCAhisto("hHistPtRecTruePrimaryMuonPlus","MuonPlusPrmiary",managermc,dcacutxy);
	muons[2]=GetOneHistFromPtDCAhisto("hHistPtRecTrueMuonMinus","MuonMinusAll",managermc,dcacutxy);
	muons[3]=GetOneHistFromPtDCAhisto("hHistPtRecTruePrimaryMuonMinus","MuonMinusPrmiary",managermc,dcacutxy);
	
	
	
	TH1F* spectraall=(TH1F*)allrecdata->Clone("recNch");
	spectraall->SetTitle("recNch");
	TH1F* contall=(TH1F*)allrecMC->Clone("contall");
	contall->SetTitle("contall");
	//contall->Add(alleff,-1);
	SubHistWithFullCorr(contall,alleff);
	alleff->Divide(alleff,allgen,1,1,"B");
	contall->Divide(contall,allrecMC,1,1,"B");
	
	GetCorrectedSpectra(spectraall,allrecdata,alleff,contall);
	Divideby2pipt(spectraall);

	allrecdata->Scale(1./neventsdata,"width");
	allgen->Scale(1./neventsmcall,"width");
	allrecMC->Scale(1./neventsmc,"width");
	spectraall->Scale(1./neventsdata,"width");


	lout->Add(allgen);
	lout->Add(allrecMC);
	lout->Add(alleff);
	lout->Add(allrecdata);
	lout->Add(spectraall);
	lout->Add(contall);
	
	for (int i=0;i<6;i++)
	{
	
		
		TString tmpname(rawspectramc[i]->GetTitle());
		tmpname.ReplaceAll("SpectraMC","%s");
		contallMC[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contallMC"));
		contallMC[i]->SetTitle(Form(tmpname.Data(),"contallMC"));
		contfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contfit"));
		contfit[i]->SetTitle(Form(tmpname.Data(),"contfit"));
		contWDfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contWDfit"));
		contWDfit[i]->SetTitle(Form(tmpname.Data(),"contWDfit"));
		contMatfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contMatfit"));
		contMatfit[i]->SetTitle(Form(tmpname.Data(),"contMatfit"));
		primaryfit[i]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"primaryfit"));
		primaryfit[i]->SetTitle(Form(tmpname.Data(),"primaryfit"));
		contfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contfitonMC"));
		contfit[i+6]->SetTitle(Form(tmpname.Data(),"contfitonMC"));
		contWDfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contWDfitonMC"));
		contWDfit[i+6]->SetTitle(Form(tmpname.Data(),"contWDfitonMC"));
		contMatfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"contMatfitonMC"));
		contMatfit[i+6]->SetTitle(Form(tmpname.Data(),"contMatfitonMC"));
		primaryfit[i+6]=(TH1F*)rawspectramc[i]->Clone(Form(tmpname.Data(),"primaryfitMC"));
		primaryfit[i+6]->SetTitle(Form(tmpname.Data(),"primaryfitMC"));
		if(i==0||i==3)
		{
			muons[i/3]->Divide(muons[i/3],rawspectramc[i],1,1,"B");
			muons[i/3+1]->Divide(muons[i/3+1],rawspectramc[i],1,1,"B");

		}	
	
		
			
		contfit[i]->Reset();
		contWDfit[i]->Reset();
		contMatfit[i]->Reset();
		primaryfit[i]->Reset();
		

		contfit[i+6]->Reset();
		contWDfit[i+6]->Reset();
		contMatfit[i+6]->Reset();
		primaryfit[i+6]->Reset();
		
		SetBintoOne(primaryfit[i]);
		SetBintoOne(primaryfit[i+6]);		
		spectra[i]=(TH1F*)rawspectradata[i]->Clone(Form(tmpname.Data(),"SpectraFinal"));
		spectra[i]->SetTitle(Form(tmpname.Data(),"SpectraFinal"));
		spectraLeonardo[i]=(TH1F*)rawspectradata[i]->Clone(Form(tmpname.Data(),"SpectraFinalLeonardo"));
		spectraLeonardo[i]->SetTitle(Form(tmpname.Data(),"SpectraFinalLeonardo"));

		corrLeonardo[i]=(TH1F*)MCTruth[i]->Clone(Form(tmpname.Data(),"CorrFactLeonardo"));
		corrLeonardo[i]->SetTitle(Form(tmpname.Data(),"CorrFactLeonardo"));			
		corrLeonardo[i]->Divide(corrLeonardo[i],rawspectramc[i],1,1,"B");
		
		
		if((options&kusePIDcontaminatiofromfile)==kusePIDcontaminatiofromfile)
		{
			CopyCorrectionFromFile(filenames[1],"conPID",contPID);
			CopyCorrectionFromFile(filenames[1],"conPIDprimary",contPIDpri);

		}
		else
		{
			//contallMC[i]->Add(eff[i],-1.0);
			SubHistWithFullCorr(contallMC[i],eff[i]);
			//RecomputeErrors(contallMC[i]);
			contallMC[i]->Sumw2(); 
			contallMC[i]->Divide(contallMC[i],rawspectramc[i],1,1,"B");
			// contamintaion from PID but only primaries
			//contPIDpri[i]->Add(eff[i],-1.0);
			SubHistWithFullCorr(contPIDpri[i],eff[i]);
			//RecomputeErrors(contPIDpri[i]);
			contPIDpri[i]->Divide(contPIDpri[i],rawspectramc[i],1,1,"B");
		}

		if((options&kuseeffcorrectionfromfile)==kuseeffcorrectionfromfile)
			CopyCorrectionFromFile(filenames[0],"eff",eff);
		else
			eff[i]->Divide(eff[i],MCTruth[i],1,1,"B");
		cout<<"EFF0"<<endl;	
		eff[i]->Print("all");
		contPID[i]->Sumw2();
		rawspectramc[i]->Sumw2();
		//contPID[i]->Add(contPID[i],rawspectramc[i],-1,1);
		SubHistWithFullCorr(contPID[i],rawspectramc[i]);
		contPID[i]->Scale(-1.0);

		//RecomputeErrors(contPID[i]);
		contPID[i]->ResetStats();
		contPID[i]->Sumw2();
		contPID[i]->Divide(contPID[i],rawspectramc[i],1,1,"B");
	
		if((options&kuseprimaryPIDcont)==kuseprimaryPIDcont)
			confinal[i]=(TH1F*)contPIDpri[i]->Clone(Form(tmpname.Data(),"confinal"));
		else	
			confinal[i]=(TH1F*)contPID[i]->Clone(Form(tmpname.Data(),"confinal"));
		confinal[i]->SetTitle(Form(tmpname.Data(),"confinal"));

		contSecMC[i]=(TH1F*)contWD[i]->Clone(Form(tmpname.Data(),"conSecMC"));
		contSecMC[i]->Add(contMat[i]);
		contWD[i]->Divide(contWD[i],rawspectramc[i],1,1,"B");
		contMat[i]->Divide(contMat[i],rawspectramc[i],1,1,"B");
		contSecMC[i]->Divide(contSecMC[i],rawspectramc[i],1,1,"B");
		if(i>2)
		{
			doubleconuntsdata[i-3]=(TH1F*)rawspectradata[i]->Clone(Form("DoublecountsDATA%s",Particle[i-3].Data()));
			doubleconuntsdata[i-3]->Reset();
			doubleconuntsMC[i-3]=(TH1F*)rawspectramc[i]->Clone(Form("DoublecountsMC%s",Particle[i-3].Data()));
			doubleconuntsMC[i-3]->Reset();

			CalculateDoubleCounts(doubleconuntsdata[i-3],rawspectradata,i-3,kTRUE);
			CalculateDoubleCounts(doubleconuntsMC[i-3],rawspectramc,i-3,kFALSE);

			
	
		}
	

	
		rawspectradata[i]->Scale(1./neventsdata,"width");
		rawspectramc[i]->Scale(1./neventsmc,"width");
		MCTruth[i]->Scale(1./neventsmcall,"width");
		spectraLeonardo[i]->Scale(1./neventsdata,"width");
		
	
		lout->Add(rawspectradata[i]);
		lout->Add(rawspectramc[i]);
		lout->Add(MCTruth[i]);
		lout->Add(eff[i]);
		lout->Add(contallMC[i]);
		lout->Add(contPID[i]);
		lout->Add(contWD[i]);
		lout->Add(contMat[i]);
		lout->Add(contfit[i]);
		lout->Add(contWDfit[i]);
		lout->Add(contMatfit[i]);
 		lout->Add(primaryfit[i]);	
		lout->Add(contfit[i+6]);
		lout->Add(contWDfit[i+6]);
		lout->Add(contMatfit[i+6]);
		lout->Add(primaryfit[i+6]);	
		lout->Add(spectra[i]);
		lout->Add(spectraLeonardo[i]);
		lout->Add(confinal[i]);
		lout->Add(contPIDpri[i]);
		lout->Add(contSecMC[i]);
		if(i>2)
		{
			lout->Add(doubleconuntsdata[i-3]);
			lout->Add(doubleconuntsMC[i-3]);
		}
	}
	lout->Add(muons[0]);
	lout->Add(muons[1]);
	lout->Add(muons[2]);
	lout->Add(muons[3]);

	outdate.ReplaceAll("/","_");
	configfile.ReplaceAll(".","_");
	TFile* fout=0x0;
	if(customoutfilename.Length()>0)
		fout=new TFile(customoutfilename.Data(),"RECREATE");
	else if(configfile.Length()>0&&(options&kuserangeonfigfile))
		fout=new TFile(Form("./results/ResMY_%s_%s_%#X_%s.root",outnamemc.Data(),outdate.Data(),options,configfile.Data()),"RECREATE");
	else
		fout=new TFile(Form("./results/ResMY_%s_%s_%#X.root",outnamemc.Data(),outdate.Data(),options),"RECREATE");
	if (((options&kdodca)==kdodca)&&((options&kuseseccontaminatiofromfile)!=kuseseccontaminatiofromfile))
		DCACorrectionMarek(managerdata,managermc,dcacutxy,fout,contfit,contWDfit,contMatfit,primaryfit,((options&kusespecialbinninginDCAfits)==kusespecialbinninginDCAfits));
	else if ((options&kuseseccontaminatiofromfile)==kuseseccontaminatiofromfile)
	{
		CopyCorrectionFromFile(filenames[2],"contfit",contfit);
		CopyCorrectionFromFile(filenames[2],"contWDfit",contWDfit);
		CopyCorrectionFromFile(filenames[2],"contMatfit",contMatfit);
		CopyCorrectionFromFile(filenames[2],"primaryfit",primaryfit);
		cout<<"CONT FIT0"<<endl;
		contfit[0]->Print("all");

	}
	else 
		cout<<"Secondary from MC"<<endl;
		

	for (int i=0;i<6;i++)
	{
			if(((options&kdodca)==kdodca)||((options&kuseseccontaminatiofromfile)==kuseseccontaminatiofromfile))
			{
				cout<<"CONT FIT"<<endl;
				contfit[i]->Print("all");
				if(((options&kusePIDfits)==kusePIDfits)&&!((options&kuseseccontaminatiofromfile)==kuseseccontaminatiofromfile))
				{
					confinal[i]->Reset();	
					if(!((options&kuseseccontaminatiofromfile)==kuseseccontaminatiofromfile))					
						Rescallesecondarycontaimation(contfit[i],contWDfit[i],contMatfit[i],primaryfit[i],primaryfractionfromfitdata[i]);
				}	
				if((options&kuseprimaryPIDcont)||(i!=1&&i!=4)) //if we do not use cont PId only for primary  and this is a kaon that do not use fit
					confinal[i]->Add(contfit[i]);
				GetCorrectedSpectra(spectra[i],rawspectradata[i],eff[i],confinal[i]);
			}
			else
			{
				GetCorrectedSpectra(spectra[i],rawspectradata[i],eff[i],contallMC[i]);	
			}
			GetCorrectedSpectraLeonardo(spectraLeonardo[i],corrLeonardo[i],primaryfit[i],primaryfit[i+6]);
			if(options&kskipconcutonspectra)
				continue;
			if(options&kuseprimaryPIDcont)
			{
				CleanHisto(spectra[i],-1,100,contPIDpri[i]);
				CleanHisto(spectraLeonardo[i],-1,100,contPIDpri[i]);		
			}
			else
			{
				CleanHisto(spectra[i],-1,100,contPID[i]);
				CleanHisto(spectraLeonardo[i],-1,100,contPID[i]);
			}
			// Apply correction for wrongly simulated TOF signal in MC 
			if(options&kuseTOFcorrforPIDsignalmatching)
				TOFPIDsignalmatchingApply(spectra[i],TOFPIDsignalmatching[i%3]);
							
	}
	
	GFCorrection(spectra,tcutsdata->GetPtTOFMatching(),options);
	GFCorrection(spectraLeonardo,tcutsdata->GetPtTOFMatching(),options);
	
	Double_t ycut=tcutsdata->GetY();
	Double_t etacut=tcutsdata->GetEtaMax()-tcutsdata->GetEtaMin();
	if(etacut<0.00001)
	{
		cout<<"using eta window 1.6"<<endl; 
		etacut=1.6;

	}

	TH1F* allch=(TH1F*)GetSumAllCh(spectra,mass,etacut);
	lout->Add(allch);	
       	if(options&kuseTOFmatchingcorrection)
	{	
		MatchingTOFEff(spectra,lout);
	 	MatchingTOFEff(spectraLeonardo);
		TOFMatchingForNch(spectraall);
	}
//	lout->ls();
	fout->cd();	
	TList* listqa=new TList();
	TList* canvaslist=new TList();
	Float_t vertexcorrection=1.0;
	Float_t corrbadchunksvtx=1.0;
	if ((options&knormalizationwithbin0integralsdata)||(options&knormalizationwithbin0integralsMC))
		corrbadchunksvtx=QAPlotsBoth(managerdata,managermc,ecutsdata,ecutsmc,tcutsdata,tcutsmc,listqa,canvaslist,0);
	else
		corrbadchunksvtx=QAPlotsBoth(managerdata,managermc,ecutsdata,ecutsmc,tcutsdata,tcutsmc,listqa,canvaslist,1);
	if (options&kveretxcorrectionandbadchunkscorr)
		vertexcorrection=corrbadchunksvtx;
	cout<<" VTX corr="<<vertexcorrection<<endl;
	if(TMath::Abs(ycut)>0.0)
		vertexcorrection=vertexcorrection/(2.0*ycut);
	for (int i=0;i<6;i++)
	{
		spectra[i]->Scale(vertexcorrection);
		spectraLeonardo[i]->Scale(vertexcorrection);
		if(TMath::Abs(ycut)>0.0)
		{
			rawspectradata[i]->Scale(1.0/(2.0*ycut));
                	rawspectramc[i]->Scale(1.0/(2.0*ycut));
                	MCTruth[i]->Scale(1.0/(2.0*ycut));
		}
	}	
	allch->Scale(vertexcorrection);
	spectraall->Scale(vertexcorrection/etacut);

	//spectraall->Scale(1.0/1.6);
	lout->Write("output",TObject::kSingleKey);	
	listqa->Write("outputQA",TObject::kSingleKey);
	canvaslist->Write("outputcanvas",TObject::kSingleKey);
	if(lfits->GetEntries()>0)
		lfits->Write("PIDfits",TObject::kSingleKey);
	fout->Close();
	//Normaliztionwithbin0integrals();

}

Bool_t   OpenFile(TString dirname,TString outputname, Bool_t mcflag, Bool_t mcasdata)
{
	

	TString nameFile = Form("./%s/AnalysisResults%s.root",dirname.Data(),(mcflag?"MC":"DATA"));
	TFile *file = TFile::Open(nameFile.Data());
	if(!file)
	{
		cout<<"no file"<<endl;
		return false;
	}	
	TString sname=Form("OutputBothSpectraTask_%s_%s",(mcflag?"MC":"Data"),outputname.Data());
	if(mcasdata)
	{
		cout<<"using MC as data "<<endl;
		sname=Form("OutputBothSpectraTask_%s_%s","MC",outputname.Data());
	}
	file->ls();
	TDirectoryFile *dir=(TDirectoryFile*)file->Get(sname.Data());
	if(!dir)
	{
	//	cout<<"no dir "<<sname.Data()<<endl;	
		if(mcasdata)
		{
			cout<<"using MC as data "<<endl;
			sname=Form("OutputAODSpectraTask_%s_%s","MC",outputname.Data());
		}
		else	
			sname=Form("OutputAODSpectraTask_%s_%s",(mcflag?"MC":"Data"),outputname.Data());
	//	cout<<"trying "<<sname.Data()<<endl;
		dir=(TDirectoryFile*)file->Get(sname.Data());
		if(!dir)
		{
			cout<<"no dir "<<sname.Data()<<endl;
			return false;
		}
	}
	cout << " -- Info about " <<(mcflag?"MC":"DATA") <<" -- "<< endl;
	if(mcflag)
	{
		managermc= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");	
		ecutsmc = (AliSpectraBothEventCuts*) dir->Get("Event Cuts");
		tcutsmc = (AliSpectraBothTrackCuts*) dir->Get("Track Cuts");
		ecutsmc->PrintCuts();
		tcutsmc->PrintCuts();
		if(!managermc||!ecutsmc||!tcutsmc)
			return false;
		if(managermc->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries()!=ecutsmc->GetHistoCuts()->GetBinContent(3))
			cout<<"Please check MC file possible problem with merging"<<" "<<managermc->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries()<<" "<<ecutsmc->GetHistoCuts()->GetBinContent(3)<<endl;
	}
	else
	{
		managerdata= (AliSpectraBothHistoManager*) dir->Get("SpectraHistos");	
		ecutsdata = (AliSpectraBothEventCuts*) dir->Get("Event Cuts");
		tcutsdata = (AliSpectraBothTrackCuts*) dir->Get("Track Cuts");
		ecutsdata->PrintCuts();
		tcutsdata->PrintCuts();
		if(!managerdata||!ecutsdata||!tcutsdata)
			return false;
		if(managerdata->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries()!=ecutsdata->GetHistoCuts()->GetBinContent(3))
			cout<<"Please check DATA file possible problem with merging"<<" "<<managerdata->GetGenMulvsRawMulHistogram("hHistGenMulvsRawMul")->GetEntries()<<" "<<ecutsdata->GetHistoCuts()->GetBinContent(3)<<endl;

	}
	return true;
}

 void GetMCTruth(TH1F** MCTruth)
 {
	for(Int_t icharge=0;icharge<2;icharge++)
	{
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			TString hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
			MCTruth[index]=(TH1F*)((TH1F*)managermc->GetPtHistogram1D(hname.Data(),1,1))->Clone();
			MCTruth[index]->SetName(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
			MCTruth[index]->SetTitle(Form("MCTruth_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
			MCTruth[index]->Sumw2(); 
		}
	}
}

TH1F* GetOneHistFromPtDCAhisto(TString name,TString hnameout,AliSpectraBothHistoManager* hman,TFormula* dcacutxy)
{
			TH1F* histo =(TH1F*)((TH1F*) hman->GetPtHistogram1D(name.Data(),-1,-1))->Clone();
			histo->SetName(hnameout.Data());
			histo->SetTitle(hnameout.Data());
		  
			if(dcacutxy)
			{
				for(int ibin=1;ibin<histo->GetNbinsX();ibin++)
				{
					Double_t lowedge=histo->GetBinLowEdge(ibin);
					Float_t cut=dcacutxy->Eval(lowedge);
					TH1F* dcahist=(TH1F*)hman->GetDCAHistogram1D(name.Data(),lowedge,lowedge);
//					Float_t inyield=dcahist->Integral(dcahist->GetXaxis()->FindBin(-1.0*cut),dcahist->GetXaxis()->FindBin(cut));
					Float_t testyield=0.0;
					Float_t testerror=0.0; 	
					for (int itest=dcahist->GetXaxis()->FindBin(-1.0*cut);itest<=dcahist->GetXaxis()->FindBin(cut);itest++)
					{
						testyield+=dcahist->GetBinContent(itest);
						testerror+=dcahist->GetBinError(itest)*dcahist->GetBinError(itest);
					}
//					cout<<"corr data "<<histo->GetBinContent(ibin)<<" "<<inyield<<" "<<dcahist->Integral()<<" "<<hnameout.Data()<<endl;
//					cout<<"test dca "<<lowedge<<" "<<dcacutxy->Eval(lowedge)<<" "<<dcacutxy->Eval(histo->GetXaxis()->GetBinUpEdge(ibin))<<" "<<dcahist->GetBinLowEdge(dcahist->GetXaxis()->FindBin(-1.0*cut))<<" "<<dcahist->GetXaxis()->GetBinUpEdge(dcahist->GetXaxis()->FindBin(-1.0*cut))<<endl;

//					cout<<testyield<<" "<<TMath::Sqrt(testerror)<<" error2 "<<inyield<<" "<<TMath::Sqrt(inyield)<<endl;
						
					//histo->SetBinContent(ibin,inyield);
					//histo->SetBinError(ibin,TMath::Sqrt(inyield));
					histo->SetBinContent(ibin,testyield);
					histo->SetBinError(ibin,TMath::Sqrt(testerror));

				}
			}
			histo->Sumw2();
			return histo;
}



	
void GetPtHistFromPtDCAhisto(TString hnamein, TString hnameout, AliSpectraBothHistoManager* hman,TH1F** histo,TFormula* dcacutxy)
{
	//Float_t min[3]={0.3,0.3,0.45};
	//Float_t max[3]={1.5,1.2,2.2};
	for(Int_t icharge=0;icharge<2;icharge++)
	{
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			Printf("Getting %s",hnamein.Data());
			TString nameinfinal=Form("%s%s%s",hnamein.Data(),Particle[ipart].Data(),Sign[icharge].Data());
			TString nameoutfinal=Form("%s%s%s",hnameout.Data(),Particle[ipart].Data(),Sign[icharge].Data());
			
			
			histo[index]=GetOneHistFromPtDCAhisto(nameinfinal,nameoutfinal,hman,dcacutxy);
			CleanHisto(histo[index],minRanges[ipart],maxRanges[ipart]);
		}
	} 
}
void CleanHisto(TH1F* h, Float_t minV, Float_t maxV,TH1* contpid)
{
	for (int i=0;i<=h->GetNbinsX();i++)
	{	
		if(h->GetXaxis()->GetBinCenter(i)<minV||h->GetXaxis()->GetBinCenter(i)>maxV)
		{
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}	
		if(contpid)
		{
			if(contpid->GetBinContent(i)>fMaxContaminationPIDMC)
			{
				h->SetBinContent(i,0);
				h->SetBinError(i,0);
			}
		}
	}
}


void DCACorrectionMarek(AliSpectraBothHistoManager* hman_data, AliSpectraBothHistoManager* hman_mc,TFormula* fun,TFile *fout,TH1F** hcon,TH1F** hconWD,TH1F** hconMat,TH1F** hprimary, Bool_t binning)
{
  printf("\n\n-> DCA Correction");  
  
 
  Printf("\n DCACorr");
  TString sample[2]={"data","mc"};
  ofstream debug("debugDCA.txt");
  TList* listofdcafits=new TList();
  for(Int_t icharge=0;icharge<2;icharge++)
  {
		for(Int_t ipart=0;ipart<3;ipart++)
		{
			Int_t index=ipart+3*icharge;
			for(Int_t isample=0;isample<2;isample++)
			{

				for(Int_t ibin_data=7;ibin_data<40;ibin_data++)
				{	
						
					Double_t lowedge=hcon[index]->GetBinLowEdge(ibin_data);
					Double_t binwidth=hcon[index]->GetBinWidth(ibin_data);
					debug<<"NEW "<<Particle[ipart].Data()<<" "<<Sign[icharge].Data()<<" "<<lowedge<<endl;
					if(fun)
					{						
						CutRange[1]=fun->Eval(lowedge);
						CutRange[0]=-1.0*CutRange[1];
					}	
					debug<<"cut  "<<CutRange[1]<<" "<<CutRange[0]<<endl;		
	     				Short_t fitsettings=DCAfitsettings(lowedge+0.5*binwidth,index);
					debug<<"settings "<< fitsettings<<" "<<isample<<endl;
					cout<<"MC FIT "<<endl;
					if(fitsettings==0)
						continue;	
						
					TCanvas *cDCA=new TCanvas(Form("cDCA%d%s%s%sbin%d",index,sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),ibin_data),Form("cDCA%d%s%s%sbin%d",index,sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),ibin_data),1700,1500);
					cDCA->SetMargin(0.1,0.02,0.1,0.02);
					TLegend* Leg1 = new TLegend(0.7,0.7,0.88,0.88,"","NDC");
					Leg1->SetFillStyle(kFALSE);
					Leg1->SetLineColor(kWhite);
					Leg1->SetBorderSize(0);

					TH1F *hToFit =0x0;	
					if(isample==0)
						hToFit =(TH1F*) ((TH1F*)hman_data->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					else if(isample==1)
						hToFit =(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					else 
						return;
					debug<<Particle[ipart].Data()<<" "<<Sign[icharge].Data()<<" "<<lowedge<<endl;
					TH1F *hmc1=0x0;
					if(hman_mc->GetIncludecorrectlyidentifiedinMCtemplates())
						hmc1=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					else
						hmc1=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaPrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					TH1F *hmc2=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryWeakDecay%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					TH1F *hmc3=(TH1F*) ((TH1F*)hman_mc->GetDCAHistogram1D(Form("hHistPtRecSigmaSecondaryMaterial%s%s",Particle[ipart].Data(),Sign[icharge].Data()),lowedge,lowedge))->Clone();
					Double_t minentries=100;
					debug<<hToFit->GetEntries()<<" "<<hmc1->GetEntries()<<" "<<hmc2->GetEntries()<<" "<<hmc3->GetEntries()<<" "<<minentries<<endl;
					debug<<((fitsettings&0x1)&&hmc2->GetEntries()<=minentries)<<" "<<((fitsettings&0x2)&&hmc3->GetEntries()<=minentries)<<endl;
					cout<<hToFit->GetEntries()<<" "<<hmc1->GetEntries()<<" "<<hmc2->GetEntries()<<" "<<hmc3->GetEntries()<<endl;
                                        cout<<((fitsettings&0x1)&&hmc2->GetEntries()<=minentries)<<" "<<((fitsettings&0x2)&&hmc3->GetEntries()<=minentries)<<endl;
					if(hToFit->GetEntries()<=minentries || hmc1->GetEntries()<=minentries)
					{
						delete 	hToFit;
						delete 	hmc1;
						delete 	hmc2;
						delete 	hmc3;
						continue;
					}	
					if(((fitsettings&0x1)&&hmc2->GetEntries()<=minentries))
						fitsettings=fitsettings-1;
					if(((fitsettings&0x2)&&hmc3->GetEntries()<=minentries))
						fitsettings=fitsettings-2;
					if(fitsettings==0)
					{
						delete 	hToFit;
						delete 	hmc1;
						delete 	hmc2;
						delete 	hmc3;
						continue;

					}	

					hToFit->Sumw2();
					hmc1->Sumw2();
					hmc2->Sumw2();
					hmc3->Sumw2();
					
					Float_t corrforrebinning[4]={1.0,1.0,1.0,1.0};	
					Int_t binCutRange[]={hToFit->GetXaxis()->FindBin(CutRange[0]),hToFit->GetXaxis()->FindBin(CutRange[1]),hmc1->GetXaxis()->FindBin(CutRange[0]),hmc1->GetXaxis()->FindBin(CutRange[1])};			

					
					if(hmc3->GetNbinsX()>100||hToFit->GetNbinsX()>100)
					{
					
						corrforrebinning[0]=hToFit->Integral(binCutRange[0],binCutRange[1]);
						corrforrebinning[1]=hmc1->Integral(binCutRange[2],binCutRange[3]);
						corrforrebinning[2]=hmc2->Integral(binCutRange[2],binCutRange[3]);
						corrforrebinning[3]=hmc3->Integral(binCutRange[2],binCutRange[3]);

						//hToFit->Rebin(30);
						//hmc1->Rebin(30);
						//hmc2->Rebin(30);
						//hmc3->Rebin(30);
						if(binning)
						{
							hToFit=ReBinDCAHisto(hToFit);
							hmc1=ReBinDCAHisto(hmc1);
							hmc2=ReBinDCAHisto(hmc2);
							hmc3=ReBinDCAHisto(hmc3);
						}
						else
						{
							Int_t datarebinning=30;
							if(managerdata->GetNRebin()>0)
								datarebinning=30/managerdata->GetNRebin();
							Int_t mcrebinning=30;
							if(managermc->GetNRebin()>0)
								mcrebinning=30/managermc->GetNRebin();
							if(isample==0)
								hToFit->Rebin(datarebinning);
							else
								hToFit->Rebin(mcrebinning);
							cout<<isample<<" Rebin "<<mcrebinning<<" "<<datarebinning<<endl;
							hmc1->Rebin(mcrebinning);
							hmc2->Rebin(mcrebinning);
							hmc3->Rebin(mcrebinning);

							
						}

						if(hToFit->GetXaxis()->GetNbins()!=hmc1->GetXaxis()->GetNbins())
						{
							delete 	hToFit;
							delete 	hmc1;
							delete 	hmc2;
							delete 	hmc3;
							return;
						}
						binCutRange[0]=hmc1->GetXaxis()->FindBin(CutRange[0]);
						binCutRange[1]=hmc1->GetXaxis()->FindBin(CutRange[1]);


						//after rebbing we lose resolution of the dca this correction also us to do obtain inside used dca

						if(hToFit->Integral(binCutRange[0],binCutRange[1])>0.0)
							corrforrebinning[0]=corrforrebinning[0]/hToFit->Integral(binCutRange[0],binCutRange[1]);
						else	
							corrforrebinning[0]=1.0;
						if(hmc1->Integral(binCutRange[0],binCutRange[1])>0.0)
							corrforrebinning[1]=corrforrebinning[1]/hmc1->Integral(binCutRange[0],binCutRange[1]);
						else	
							corrforrebinning[1]=1.0;
						if(hmc2->Integral(binCutRange[0],binCutRange[1])>0.0)
							corrforrebinning[2]=corrforrebinning[2]/hmc2->Integral(binCutRange[0],binCutRange[1]);
						else	
							corrforrebinning[2]=1.0;
						if(hmc3->Integral(binCutRange[0],binCutRange[1])>0.0)
							corrforrebinning[3]=corrforrebinning[3]/hmc3->Integral(binCutRange[0],binCutRange[1]);
						else	
							corrforrebinning[3]=1.0;
							


					}

					cDCA->cd();
					gPad->SetGridy();
					gPad->SetGridx();
					gPad->SetLogy();
	 
					TObjArray *mc=0x0;
					Int_t Npar=3;	
					if(fitsettings==3)
						mc = new TObjArray(3);        // MC histograms are put in this array
					else
					{
						mc = new TObjArray(2);
						Npar=2;
					}
					mc->Add(hmc1);
					if(fitsettings&0x1)
						mc->Add(hmc2);
					if(fitsettings&0x2)
						mc->Add(hmc3);
					TFractionFitter* fit = new TFractionFitter(hToFit,mc); // initialise
					fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					if(fitsettings&0x1)
						fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					if(fitsettings&0x2)
						fit->Constrain(1+(fitsettings&0x1),0.0,1.0);               // constrain fraction 1 to be between 0 and 1
					
					Int_t binFitRange[]={hmc1->GetXaxis()->FindBin(FitRange[0]),hmc1->GetXaxis()->FindBin(FitRange[1])};			
					fit->SetRangeX(binFitRange[0],binFitRange[1]);
					hToFit->GetXaxis()->SetRange(binFitRange[0],binFitRange[1]);
				//	hToFit->SetTitle(Form("DCA distr - %s %s %s %lf",sample[isample].Data(),Particle[ipart].Data(),Sign[icharge].Data(),lowedge));
					hToFit->SetTitle("");
					Int_t nfits=0;
					Int_t status = 1;
					while (status!=0&&nfits<20)
					{
						for (int iparm = 0; iparm < Npar; ++iparm)
                                        	{
                                                	TString name("frac"); name += iparm;
                                                	if(iparm==0)
                                                	        fit->GetFitter()->SetParameter(iparm,name.Data(),0.85-0.02*nfits,0.01,0.0,1.0);
                                                	else if (Npar==2)
                                                	        fit->GetFitter()->SetParameter(iparm,name.Data(),0.15+0.02*nfits,0.01,0.0,1.0);
                                                	else
                                                        	fit->GetFitter()->SetParameter(iparm,name.Data(),0.075+0.01*nfits,0.01,0.0,1.0);

                                        	}
						status = fit->Fit(); // perform the fit
						nfits++;
					}
					cout << "fit status: " << status << endl;
					debug<<"fit status: " << status << endl;
		
					if (status == 0)
					{ 	
						Double_t v1=0.0,v2=0.0,v3=0.0;
						Double_t ev1=0.0,ev2=0.0,ev3=0.0;
						Double_t cov=0.0;
						                   
						// check on fit status
						TH1F* result = (TH1F*) fit->GetPlot();
						TH1F* PrimMCPred=(TH1F*)fit->GetMCPrediction(0);
						
						TH1F* secStMCPred=0X0;
					
						TH1F* secMCPred=0x0;
	    
						//first method, use directly the fit result
						fit->GetResult(0,v1,ev1);

						if(fitsettings&0x1)
						{
							fit->GetResult(1,v2,ev2);
							secStMCPred=(TH1F*)fit->GetMCPrediction(1);

						}
						if(fitsettings&0x2)
						{
							fit->GetResult(1+(fitsettings&0x1),v3,ev3);
							secMCPred=(TH1F*)fit->GetMCPrediction(1+(fitsettings&0x1));
							if(fitsettings&0x1)
								cov=fit->GetFitter()->GetCovarianceMatrixElement(1,2);


						}
						debug<<v1<<" "<<ev1<<" "<<v2<<" "<<ev2<<" "<<v3<<" "<<ev3<<" "<<" "<<cov<<" "<<endl;
					
	    					// becuase dca cut range is not a fit range the results from TFractionFitter should be rescale
						// This can be done in two ways or use input histograms or output histograms
						// The difference between those two methods should be on the level of statistical error	 	
	   					// I use output histograms 
						
 						// Method 1 input histo 	

						debug<<hToFit->Integral(binFitRange[0],binFitRange[1])<<" "<<hmc1->Integral(binFitRange[0],binFitRange[1])<<" "<<hmc2->Integral(binFitRange[0],binFitRange[1])<<" "<<hmc3->Integral(binFitRange[0],binFitRange[1])<<endl;
						Float_t normalizationdata=hToFit->Integral(hToFit->GetXaxis()->FindBin(CutRange[0]),hToFit->GetXaxis()->FindBin(CutRange[1]))/hToFit->Integral(binFitRange[0],binFitRange[1]);
						debug<<normalizationdata<<" "<<hToFit->Integral(hToFit->GetXaxis()->FindBin(CutRange[0]),hToFit->GetXaxis()->FindBin(CutRange[1]))<<endl;	
						normalizationdata*=corrforrebinning[0];
						debug<<normalizationdata<<endl;

						Float_t normalizationmc1=(hmc1->Integral(binCutRange[0],binCutRange[1])/hmc1->Integral(binFitRange[0],binFitRange[1]))/normalizationdata;
						Float_t normalizationmc2=0.0;
						if(fitsettings&0x1)
							normalizationmc2=(hmc2->Integral(binCutRange[0],binCutRange[1])/hmc2->Integral(binFitRange[0],binFitRange[1]))/normalizationdata;
						Float_t normalizationmc3=0.0;
						if(fitsettings&0x2)
							normalizationmc3=(hmc3->Integral(binCutRange[0],binCutRange[1])/hmc3->Integral(binFitRange[0],binFitRange[1]))/normalizationdata;

						normalizationmc1*=corrforrebinning[1];
						normalizationmc2*=corrforrebinning[2];
						normalizationmc3*=corrforrebinning[3];

						debug<<"After Nor"<<endl;
						cout<<"A"<<endl;
						debug<<v1*normalizationmc1<<" "<<ev1*normalizationmc1<<" "<<v2*normalizationmc2<<" "<<ev2*normalizationmc2<<" "<<v3*normalizationmc3<<" "<<ev3*normalizationmc3<<" "<<endl;
						cout<<"B"<<endl;
						debug<<1.0-v1*normalizationmc1<<" "<<ev1*normalizationmc1<<" "<<v2*normalizationmc2+v3*normalizationmc3<<" "<<TMath::Sqrt(ev2*ev2*normalizationmc2*normalizationmc2+ev3*ev3*normalizationmc3*normalizationmc3+cov*normalizationmc3*normalizationmc2)<<endl;
						cout<<"C"<<endl;

						debug<<"addtional info"<<endl;


					       //Method 2 output histo		

						Float_t normalizationdata1=result->Integral(binCutRange[0],binCutRange[1])/result->Integral(binFitRange[0],binFitRange[1]);
						

						// if the cut range is bigger the fit range we should calculate the normalization factor for data using the data histogram 
						// because result histogram has entries only in fits range 	 
						if(FitRange[0]>CutRange[0]||FitRange[1]<CutRange[1])	
							normalizationdata1=normalizationdata;
					
						normalizationdata1*=corrforrebinning[0];


						Float_t normalizationmc11=(PrimMCPred->Integral(binCutRange[0],binCutRange[1])/PrimMCPred->Integral(binFitRange[0],binFitRange[1]))/normalizationdata1;
						Float_t normalizationmc21=0.0;
						if(fitsettings&0x1)
							normalizationmc21=(secStMCPred->Integral(binCutRange[0],binCutRange[1])/secStMCPred->Integral(binFitRange[0],binFitRange[1]))/normalizationdata1;
						Float_t normalizationmc31=0.0;
						if(fitsettings&0x2)
							normalizationmc31=(secMCPred->Integral(binCutRange[0],binCutRange[1])/secMCPred->Integral(binFitRange[0],binFitRange[1]))/normalizationdata1;
						
						normalizationmc11*=corrforrebinning[1];
						normalizationmc21*=corrforrebinning[2];
						normalizationmc31*=corrforrebinning[3];

						debug<<"After Nor 2"<<endl;
						debug<<v1*normalizationmc11<<" "<<ev1*normalizationmc11<<" "<<v2*normalizationmc21<<" "<<ev2*normalizationmc21<<" "<<v3*normalizationmc31<<" "<<ev3*normalizationmc31<<endl;
						
						debug<<1.0-v1*normalizationmc11<<" "<<ev1*normalizationmc11<<" "<<v2*normalizationmc21+v3*normalizationmc31<<" "<<TMath::Sqrt(ev2*ev2*normalizationmc21*normalizationmc21+ev3*ev3*normalizationmc31*normalizationmc31+cov*normalizationmc31*normalizationmc21)<<endl;
					
						
						hconWD[index+6*isample]->SetBinContent(ibin_data,v2*normalizationmc21);
						hconWD[index+6*isample]->SetBinError(ibin_data,ev2*normalizationmc21);
						hconMat[index+6*isample]->SetBinContent(ibin_data,v3*normalizationmc31);
						hconMat[index+6*isample]->SetBinError(ibin_data,ev3*normalizationmc31);
						hprimary[index+6*isample]->SetBinContent(ibin_data,v1*normalizationmc11);
						hprimary[index+6*isample]->SetBinError(ibin_data,ev1*normalizationmc11);
						hcon[index+6*isample]->SetBinContent(ibin_data,v2*normalizationmc21+v3*normalizationmc31);
						hcon[index+6*isample]->SetBinError(ibin_data,TMath::Sqrt(ev2*ev2*normalizationmc21*normalizationmc21+ev3*ev3*normalizationmc31*normalizationmc31+cov*normalizationmc31*normalizationmc21));
						
						
						
						//Drawing section
						result->Scale(1.0/result->Integral(result->GetXaxis()->FindBin(FitRange[0]),result->GetXaxis()->FindBin(FitRange[1])));
						hToFit->Scale(1.0/hToFit->Integral(hToFit->GetXaxis()->FindBin(FitRange[0]),hToFit->GetXaxis()->FindBin(FitRange[1])));
						PrimMCPred->Scale(v1/PrimMCPred->Integral(PrimMCPred->GetXaxis()->FindBin(FitRange[0]),PrimMCPred->GetXaxis()->FindBin(FitRange[1])));

						hToFit->SetMinimum(0.0001);
						hToFit->GetXaxis()->SetTitle("DCA_{xy} (cm)");
						hToFit->GetYaxis()->SetTitle("1/N_{counts}(-3cm;3cm) dN/ddca_{xy} (cm)^{-1}");
						hToFit->GetYaxis()->SetTitleOffset(1.3);
						hToFit->Scale(1.0,"width");
						hToFit->SetMarkerStyle(20);
						hToFit->SetMarkerColor(kBlack);
						hToFit->SetMarkerSize(1.5);
						hToFit->DrawClone("E1x0");
						Leg1->AddEntry(hToFit,"data","p");
						result->SetTitle("Fit result");
						result->SetLineColor(kBlack);
						Leg1->AddEntry(result,"fit result","l");
						 result->Scale(1.0,"width");
						result->DrawClone("histsame");
					
						PrimMCPred->SetLineColor(kGreen+2);
						PrimMCPred->SetLineStyle(2);
						 PrimMCPred->SetLineWidth(3.0);
						Leg1->AddEntry(PrimMCPred,"primaries","l");
						PrimMCPred->Scale(1.0,"width");
						PrimMCPred->DrawClone("histsame");
						if(fitsettings&0x1)
						{

							secStMCPred->Scale(v2/secStMCPred->Integral(secStMCPred->GetXaxis()->FindBin(FitRange[0]),secStMCPred->GetXaxis()->FindBin(FitRange[1])));
							secStMCPred->SetLineColor(kRed);
							secStMCPred->SetLineWidth(3.0);

							secStMCPred->SetLineStyle(3);
							Leg1->AddEntry(secStMCPred,"weak decays","l");
							secStMCPred->Scale(1.0,"width");
							secStMCPred->DrawClone("histsame");

						}
						if(fitsettings&0x2)
						{
							
							secMCPred->Scale(v3/secMCPred->Integral(secMCPred->GetXaxis()->FindBin(FitRange[0]),secMCPred->GetXaxis()->FindBin(FitRange[1])));
							secMCPred->SetLineColor(kBlue);
							secMCPred->SetLineWidth(3.0);

							secMCPred->SetLineStyle(4);	
							Leg1->AddEntry(secMCPred,"material","l");
							secMCPred->Scale(1.0,"width");
							secMCPred->DrawClone("histsame");
	    
						}   
					}				
					else
					{
						hconWD[index+6*isample]->SetBinContent(ibin_data,0.0);
						hconWD[index+6*isample]->SetBinError(ibin_data,0.0);
						hconMat[index+6*isample]->SetBinContent(ibin_data,0.0);
						hconMat[index+6*isample]->SetBinError(ibin_data,0.0);					
						hcon[index+6*isample]->SetBinContent(ibin_data,0.0);
						hcon[index+6*isample]->SetBinError(ibin_data,0.0);
						hprimary[index+6*isample]->SetBinContent(ibin_data,1.0);
						hprimary[index+6*isample]->SetBinError(ibin_data,0.0);
					}
					Leg1->DrawClone();
					TLatex* texttitle=new TLatex();
					texttitle->SetNDC();
					texttitle->SetTextSize(0.04);
					texttitle->DrawLatex(0.12,0.92,Form("%s %.2f<#it{p}_{T} < %.2f (GeV/#it{c})",symboles[index].Data(),lowedge,binwidth+lowedge));
					listofdcafits->Add(cDCA);
					
					//cDCA->Write();
					delete 	hmc1;
					delete 	hmc2;
					delete 	hmc3;
					delete hToFit;
				}
	

			}

		}
	}
	fout->cd();
	listofdcafits->Write("DCAfits",TObject::kSingleKey);	
}

void RecomputeErrors(TH1* h)
{
	for (int i=0; i<=h->GetXaxis()->GetNbins(); i++)
	{
		cout<<h->GetBinContent(i)<<" "<<h->GetBinError(i)<<" error "<<TMath::Sqrt(h->GetBinContent(i))<<endl;
		h->SetBinError(i,TMath::Sqrt(h->GetBinContent(i)));
	}
	h->Sumw2(); 	
}
void SetBintoOne(TH1* h)
{
	for (int i=0;i<=h->GetXaxis()->GetNbins(); i++) 
	{
		h->SetBinContent(i,1);
		h->SetBinError(i,0);
	}
}


void GetCorrectedSpectra(TH1F* corr,TH1F* raw,TH1F* eff, TH1F* con)
{
	for (int i=0;i<=corr->GetXaxis()->GetNbins(); i++) 
	{
		corr->SetBinContent(i,1);
		corr->SetBinError(i,0);
	}
	
	corr->Sumw2();
	cout<<"CORR"<<endl;
	corr->Print("all");	 
	corr->Add(con,-1);
	corr->Sumw2();
	cout<<"EFF"<<endl;
	eff->Print("all"); 	
	corr->Divide(eff);
	corr->Sumw2(); 
	eff->Print("all");
	corr->Multiply(raw);
	cout<<"Final"<<endl;
	corr->Print("all");
	corr->Sumw2(); 
	corr->Print("all");
}
void GetCorrectedSpectraLeonardo(TH1F* spectra,TH1F* correction, TH1F* hprimaryData,TH1F* hprimaryMC)
{
	spectra->Sumw2(); 
	spectra->Multiply(correction);
	spectra->Sumw2(); 
	hprimaryData->Sumw2(); 
	spectra->Multiply(hprimaryData);
	hprimaryMC->Sumw2(); 
	spectra->Divide(hprimaryMC);
}

void GFCorrection(TH1F **Spectra,Float_t tofpt,UInt_t options)
{
	TFile *fGeanFlukaK=0x0;
	if (options&kgeantflukaKaon)
	{		
	 	 //Geant/Fluka Correction
	  	Printf("\nGF correction for Kaons");
	  	//Getting GF For Kaons in TPC
	  	TGraph *gGFCorrectionKaonPlus=new TGraph();
	  	gGFCorrectionKaonPlus->SetName("gGFCorrectionKaonPlus");
	  	gGFCorrectionKaonPlus->SetTitle("gGFCorrectionKaonPlus");
	  	TGraph *gGFCorrectionKaonMinus=new TGraph();
	  	gGFCorrectionKaonMinus->SetName("gGFCorrectionKaonMinus");
	  	gGFCorrectionKaonMinus->SetTitle("gGFCorrectionKaonMinus");
	 	 TString fnameGeanFlukaK="GFCorrection/correctionForCrossSection.321.root";
  	  	fGeanFlukaK= TFile::Open(fnameGeanFlukaK.Data());
	  	if (!fGeanFlukaK)
		{
			fnameGeanFlukaK="$ALICE_PHYSICS/PWGLF/SPECTRA/PiKaPr/TestAOD/correctionForCrossSection.321.root";
			fGeanFlukaK= TFile::Open(fnameGeanFlukaK.Data());
			if (!fGeanFlukaK)
				return;
		}
	  	TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
	  	TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
	  	//getting GF func for Kaons with TOF
	  	TF1 *fGFKPosTracking;
	  	fGFKPosTracking = TrackingEff_geantflukaCorrection(3,kPositive);
	  	TF1 *fGFKNegTracking;
	 	 fGFKNegTracking = TrackingEff_geantflukaCorrection(3,kNegative);
	 	 TF1 *fGFKPosMatching;
	 	 fGFKPosMatching = TOFmatchMC_geantflukaCorrection(3,kPositive);
	 	 TF1 *fGFKNegMatching;
	 	 fGFKNegMatching = TOFmatchMC_geantflukaCorrection(3,kNegative);
		 if(tcutsdata->GetUseTypeDependedTOFCut())
			tofpt=tcutsdata->GetPtTOFMatchingKaon();

	 	 for(Int_t binK=0;binK<=Spectra[1]->GetNbinsX();binK++)
	  	{
			if(Spectra[1]->GetBinCenter(binK)<tofpt)
			{//use TPC GeantFlukaCorrection
				Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(Spectra[1]->GetBinCenter(binK)));
				Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(Spectra[4]->GetBinCenter(binK)));
			//	Printf("TPC Geant/Fluka: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPos,FlukaCorrKNeg);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPos);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNeg);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPos);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNeg);
				gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),FlukaCorrKPos);
				gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),FlukaCorrKNeg);
			}
			else
			{
				gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),0);
				gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),0);
				Float_t FlukaCorrKPosTracking=fGFKPosTracking->Eval(Spectra[1]->GetBinCenter(binK));
				Float_t FlukaCorrKNegTracking=fGFKNegTracking->Eval(Spectra[1]->GetBinCenter(binK));
			//	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosTracking,FlukaCorrKNegTracking);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosTracking);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegTracking);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosTracking);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegTracking);
				Float_t FlukaCorrKPosMatching=fGFKPosMatching->Eval(Spectra[1]->GetBinCenter(binK));
				Float_t FlukaCorrKNegMatching=fGFKNegMatching->Eval(Spectra[1]->GetBinCenter(binK));
			//	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosMatching,FlukaCorrKNegMatching);
				Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosMatching);
				Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegMatching);
				Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosMatching);
				Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegMatching);
			}
	  	}
	  }
	  if(!(options&kgeantflukaProton))	
		return;
	  //Geant Fluka for P in TPC
	  Printf("\nGF correction for Protons");
	  const Int_t kNCharge=2;
	  Int_t kPos=0;
	  Int_t kNeg=1;
          TString fnameGFProtons= "GFCorrection/correctionForCrossSection.root";
	  TFile* fGFProtons = TFile::Open(fnameGFProtons.Data());
	  if (!fGFProtons)
	  { 
		fnameGFProtons="$ALICE_PHYSICS/PWGLF/SPECTRA/PiKaPr/TestAOD/correctionForCrossSection.root";
		fGFProtons = TFile::Open(fnameGFProtons.Data());
		if (!fGFProtons)
			return;
	  }
		


	  TH2D * hCorrFluka[kNCharge];
	  hCorrFluka[kPos] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionProtons");
	  hCorrFluka[kNeg] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionAntiProtons");
	  //getting GF func for Kaons with TPCTOF
	  TF1 *fGFpPosTracking;
	  fGFpPosTracking = TrackingEff_geantflukaCorrection(4,kPositive);
	  TF1 *fGFpNegTracking;
	  fGFpNegTracking = TrackingEff_geantflukaCorrection(4,kNegative);
	  TF1 *fGFpPosMatching;
	  fGFpPosMatching = TOFmatchMC_geantflukaCorrection(4,kPositive);
	  TF1 *fGFpNegMatching;
	  fGFpNegMatching = TOFmatchMC_geantflukaCorrection(4,kNegative);
	  
	 
	   Int_t nbins = Spectra[2]->GetNbinsX();
	   if(tcutsdata->GetUseTypeDependedTOFCut())
		tofpt=tcutsdata->GetPtTOFMatchingProton();

		for(Int_t ibin = 0; ibin < nbins; ibin++)
		{
			if(Spectra[2]->GetBinCenter(ibin)<tofpt)
			{//use TPC GeantFlukaCorrection
				for(Int_t icharge = 0; icharge < kNCharge; icharge++)
				{
					Int_t nbinsy=hCorrFluka[icharge]->GetNbinsY();
					Float_t pt = Spectra[2]->GetBinCenter(ibin);
					Float_t minPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(1);
					Float_t maxPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
					if (pt < minPtCorrection) 
						pt = minPtCorrection+0.0001;
					if (pt > maxPtCorrection) 
						pt = maxPtCorrection;
					Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
				
					if (correction > 0.0) 
					{// If the bin is empty this is a  0
			//			cout<<icharge<<" "<<ibin<<" "<<correction<<endl;
						Spectra[icharge*3+2]->SetBinContent(ibin,Spectra[icharge*3+2]->GetBinContent(ibin)*correction);
						Spectra[icharge*3+2]->SetBinError(ibin,Spectra[icharge*3+2]->GetBinError(ibin)*correction);
					}
					else if (Spectra[icharge*3+2]->GetBinContent(ibin) > 0.0) 
					{ 
						// If we are skipping a non-empty bin, we notify the user
						cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for " <<"protons Pt:"<< pt<< endl;
						cout << " Bin content: " << Spectra[icharge*3+2]->GetBinContent(ibin)  << endl;
					}
				}
			}
			else
			{
				Float_t FlukaCorrpPosTracking=fGFpPosTracking->Eval(Spectra[2]->GetBinCenter(ibin));
				Float_t FlukaCorrpNegTracking=fGFpNegTracking->Eval(Spectra[2]->GetBinCenter(ibin));
			//	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosTracking,FlukaCorrpNegTracking);
				Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosTracking);
				Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegTracking);
				Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosTracking);
				Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegTracking);
				Float_t FlukaCorrpPosMatching=fGFpPosMatching->Eval(Spectra[2]->GetBinCenter(ibin));
				Float_t FlukaCorrpNegMatching=fGFpNegMatching->Eval(Spectra[2]->GetBinCenter(ibin));
			//	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosMatching,FlukaCorrpNegMatching);
				Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosMatching);
				Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegMatching);
				Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosMatching);
				Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegMatching);
			}		
		}
	if(fGeanFlukaK)
	{	
		 fGeanFlukaK->Close();
		 delete fGeanFlukaK;
	}
}


///////////
TF1 *
TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

	TF1 *f = 0x0;
  if (ipart == 3 && icharge == kNegative) {
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t
TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t
TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t
TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}
///////////////////////////////////////////
TF1 *
TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge)
{
 TF1 *f = 0x0;

  if (ipart == 3 && icharge == kNegative) {
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "MatchingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "MatchingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge].Data()), "MatchingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}


Double_t
MatchingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t 
MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  Float_t ptTPCoutP =pTmc*(1-6.81059e-01*TMath::Exp(-pTmc*4.20094));
  return (TMath::Power(1 - 0.129758*TMath::Exp(-ptTPCoutP*0.679612),0.07162/0.03471));
}

Double_t
MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  Float_t ptTPCoutK=pTmc*(1- 3.37297e-03/pTmc/pTmc - 3.26544e-03/pTmc);
  return TMath::Min((TMath::Power(0.972865 + 0.0117093*ptTPCoutK,0.07162/0.03471)), 1.);
}
void MatchingTOFEff(TH1F** Spectra, TList* list)
{
	  if(TOFMatchingScalling[0]<0.0&&TOFMatchingScalling[1]<0.0)
	  {
		TH1F *hMatcEffPos_data=(TH1F*)tcutsdata->GetHistoNMatchedPos();
		  hMatcEffPos_data->Sumw2();
		  //hMatcEffPos_data->Divide((TH1F*)tcutsdata->GetHistoNSelectedPos());
		  hMatcEffPos_data->Divide(hMatcEffPos_data,(TH1F*)tcutsdata->GetHistoNSelectedPos(),1,1,"B");
		  hMatcEffPos_data->SetTitle("Matching Eff Pos - data");
		  TH1F *hMatcEffNeg_data=(TH1F*)tcutsdata->GetHistoNMatchedNeg();
		  hMatcEffNeg_data->Sumw2();
		  //hMatcEffNeg_data->Divide((TH1F*)tcutsdata->GetHistoNSelectedNeg());
		  hMatcEffNeg_data->Divide(hMatcEffNeg_data,(TH1F*)tcutsdata->GetHistoNSelectedNeg(),1,1,"B");
		  hMatcEffNeg_data->SetTitle("Matching Eff Neg - data");
		  TH1F *hMatcEffPos_mc=(TH1F*)tcutsmc->GetHistoNMatchedPos();
		  hMatcEffPos_mc->Sumw2();
		  //hMatcEffPos_mc->Divide((TH1F*)tcutsmc->GetHistoNSelectedPos());
		  hMatcEffPos_mc->Divide(hMatcEffPos_mc,(TH1F*)tcutsmc->GetHistoNSelectedPos(),1,1,"B");
		  hMatcEffPos_mc->SetTitle("Matching Eff Pos - mc");
		  TH1F *hMatcEffNeg_mc=(TH1F*)tcutsmc->GetHistoNMatchedNeg();
		  hMatcEffNeg_mc->Sumw2();
		  //hMatcEffNeg_mc->Divide((TH1F*)tcutsmc->GetHistoNSelectedNeg());
		  hMatcEffNeg_mc->Divide(hMatcEffNeg_mc,(TH1F*)tcutsmc->GetHistoNSelectedNeg(),1,1,"B");
		  hMatcEffNeg_mc->SetTitle("Matching Eff Neg - mc");


		  hMatcEffPos_data->Divide(hMatcEffPos_mc);
		  hMatcEffNeg_data->Divide(hMatcEffNeg_mc);
		  hMatcEffPos_data->SetName("MatchingTOFPos");
		  hMatcEffNeg_data->SetName("MatchingTOFNeg");
		  
		  
		  TF1 *pol0MatchPos_data=new TF1("pol0MatchPos_data","pol0",.6,5);
		  hMatcEffPos_data->Fit("pol0MatchPos_data","MNR");
		  TF1 *pol0MatchNeg_data=new TF1("pol0MatchNeg_data","pol0",.6,5);
		  hMatcEffNeg_data->Fit("pol0MatchNeg_data","MNR");
		
		list->Add(hMatcEffPos_data);
		  list->Add(hMatcEffNeg_data);
		
			
		  TOFMatchingScalling[0]=pol0MatchPos_data->GetParameter(0);
		  TOFMatchingScalling[1]=pol0MatchNeg_data->GetParameter(0);

	  }
	  //Correction spectra for matching efficiency
	  //For the moment I'm using the inclusive correction
	  for(Int_t ipart=0;ipart<3;ipart++)
	  {
		  
		for(Int_t ibin=1;ibin<Spectra[ipart]->GetNbinsX();ibin++)
		{
			Float_t ptspectra=Spectra[ipart]->GetBinCenter(ibin);
			if(ptspectra<tcutsdata->GetPtTOFMatching())
				continue;
		  //Spectra[ipart]->SetBinContent(ibin,( Spectra[ipart]->GetBinContent(ibin)/hMatcEffPos_data->GetBinContent(hMatcEffPos_data->FindBin(ptspectra))));
		  //Spectra[ipart+3]->SetBinContent(ibin,( Spectra[ipart+3]->GetBinContent(ibin)/hMatcEffNeg_data->GetBinContent(hMatcEffNeg_data->FindBin(ptspectra))));
			Spectra[ipart]->SetBinContent(ibin,( Spectra[ipart]->GetBinContent(ibin)/TOFMatchingScalling[0]));
			Spectra[ipart+3]->SetBinContent(ibin,( Spectra[ipart+3]->GetBinContent(ibin)/TOFMatchingScalling[1]));
		}
	  }
}
Double_t eta2y(Double_t pt, Double_t mass, Double_t eta)
{
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

TH1* GetSumAllCh(TH1F** spectra, Double_t* mass,Double_t etacut)
{
	TH1F* allch=(TH1F*)spectra[0]->Clone("allCh");
	allch->Reset();
	for (int i=0;i<6;i++)
	{
		Double_t masstmp=mass[i%3];
		for (int j=1;j<spectra[i]->GetXaxis()->GetNbins();j++)
		{
			Float_t value=spectra[i]->GetBinContent(j);
			Float_t error=spectra[i]->GetBinError(j);
			if(value>0.0)
			{
				Float_t pt=spectra[i]->GetXaxis()->GetBinCenter(j);
				Float_t mt2=masstmp*masstmp+pt*pt;		
				Float_t mt=TMath::Sqrt(mt2);
				Float_t maxy=eta2y(pt,masstmp,etacut/2.0);
				Float_t conver=maxy*(TMath::Sqrt(1-masstmp*masstmp/(mt2*TMath::CosH(maxy)*TMath::CosH(maxy)))+TMath::Sqrt(1-masstmp*masstmp/(mt2*TMath::CosH(0.0)*TMath::CosH(0.0))));
				conver=conver/etacut;
				//cout<<maxy<<" "<<conver<<" "<<masstmp<<""<<spectra[i]->GetName()<<endl;
				Float_t bincontent=allch->GetBinContent(j);
				Float_t binerror=allch->GetBinError(j);
				bincontent+=conver*value;
				binerror=TMath::Sqrt(binerror*binerror+conver*conver*error*error);
				allch->SetBinContent(j,bincontent);
				allch->SetBinError(j,binerror);
			}
			
		}
	}
	Divideby2pipt(allch);
	allch->SetTitle("N_{ch};p_{T} (GeV/c);1/(2 #pi p_{T})dN/p_{T} |#eta|<0.8");		
	return allch;
}

void Divideby2pipt(TH1* hist)
{

	for (int i=1;i<hist->GetXaxis()->GetNbins();i++)
	{
		Float_t value=hist->GetBinContent(i);
		Float_t error=hist->GetBinError(i);
		Float_t pt=hist->GetXaxis()->GetBinCenter(i);
		hist->SetBinContent(i,value/(2.0*TMath::Pi()*pt));
		hist->SetBinError(i,error/(2.0*TMath::Pi()*pt));

	}
}

Short_t DCAfitsettings (Float_t pt, Int_t type)
{
	Short_t value=0x0;
	if (pt<maxptformaterial[type]&&pt>minptformaterial[type])
		value=value+2;
	if (pt<maxptforWD[type]&&pt>minptforWD[type])
		value=value+1;
	return value;	

} 

Float_t Normaliztionwithbin0integrals(UInt_t options)
{
	
	TH1F* bin0mcRec=(TH1F*)ecutsmc->GetHistoVtxGenerated()->Clone("Bin0_rec");
	TH1F* bin0mcMC=(TH1F*)ecutsmc->GetHistoVtxGenerated()->Clone("Bin0_MC");

	TH1F* vertexmc=ecutsmc->GetHistoVtxAftSelwithoutZvertexCut(); 
	TH1F* vertexmcMCz=ecutsmc->GetHistoVtxAftSelwithoutZvertexCutusingMCz(); 
	TH1F* vertexdata=ecutsdata->GetHistoVtxAftSelwithoutZvertexCut();

	TH1I* histodata=ecutsdata->GetHistoCuts();
	TH1I* histomc=ecutsmc->GetHistoCuts();

	Float_t dataevents=(Float_t)histodata->GetBinContent(3);
	//cout<<histodata->GetBinContent(2)<<endl;
	Float_t databin0events=((Float_t)histodata->GetBinContent(2))-((Float_t)histodata->GetBinContent(4));	

	bin0mcRec->Sumw2();
	bin0mcMC->Sumw2();
		
	bin0mcRec->Add(vertexmc,-1);
	bin0mcMC->Add(vertexmcMCz,-1);
	
	bin0mcRec->Divide(vertexmc);
	bin0mcMC->Divide(vertexmcMCz);
	
	bin0mcRec->Multiply(vertexdata);
	bin0mcMC->Multiply(vertexdata);
	
	Float_t bin0mcRecN=0.0;
	Float_t bin0mcMCN=0.0;

	for (int i=0;i<=bin0mcRec->GetXaxis()->GetNbins();i++)
	{
		bin0mcRecN+=bin0mcRec->GetBinContent(i);
		bin0mcMCN+=bin0mcMC->GetBinContent(i);

	}
	bin0mcRec->Scale(databin0events/bin0mcRecN);
	bin0mcMC->Scale(databin0events/bin0mcMCN);		
	
	Int_t binmin=bin0mcRec->GetXaxis()->FindBin(-10);
	Int_t binmax=bin0mcRec->GetXaxis()->FindBin(10)-1;
	cout<<	bin0mcRec->GetXaxis()->GetBinLowEdge(binmin)<<" "<<bin0mcRec->GetXaxis()->GetBinUpEdge(binmax)<<endl;
	cout<<bin0mcRecN<<" "<<bin0mcMCN<<" "<<databin0events<<endl;	
	cout<<dataevents<<" normalization "<<dataevents+bin0mcRec->Integral(binmin,binmax)<<" "<<dataevents+bin0mcMC->Integral(binmin,binmax)<<endl;
	cout<<histodata->GetBinContent(2)<<" "<<histodata->GetBinContent(4)<<endl;
	if ((options&knormalizationwithbin0integralsdata)==knormalizationwithbin0integralsdata)
		return 	dataevents+bin0mcRec->Integral(binmin,binmax);
 	else if ((options&knormalizationwithbin0integralsMC)==knormalizationwithbin0integralsMC)   
		return dataevents+bin0mcMC->Integral(binmin,binmax) ;
	else
		return 1;		
}
 

Bool_t ReadConfigFile(TString configfile)
{
	ifstream infile(configfile.Data());
	if(infile.is_open()==false)
		return false;
	TString namesofSetting[]={"CutRangeMin","CutRangeMax","FitRangeMin","FitRangeMax","MinMatPionPlus","MaxMatPionPlus","MinMatKaonPlus","MaxMatKaonPlus","MinMatProtonPlus","MaxMatProtonPlus","MinMatPionMinus","MaxMatPionMinus","MinMatKaonMinus","MaxMatKaonMinus","MinMatProtonMinus","MaxMatProtonMinus","MinWDPionPlus","MaxWDPionPlus","MinWDKaonPlus","MaxWDKaonPlus","MinWDProtonPlus","MaxWDProtonPlus","MinWDPionMinus","MaxWDPionMinus","MinWDKaonMinus","MaxWDKaonMinus","MinWDProtonMinus","MaxWDProtonMinus","MaxContaminationPIDMC","MinPions","MaxPions","MinKaons","MaxKaons","MinProtons","MaxProtons","TOFPIDsignalmatchPion","TOFPIDsignalmatchKaon","TOFPIDsignalmatchProton","NamefileEff","NameFilePIDcon","NameFileSeccon"};	

	char buffer[256];
	while (infile.eof()==false)
	{
		buffer[0]='#'; 
		while (buffer[0]=='#'&&infile.eof()==false)
			infile.getline(buffer,256);
		TString tmpstring(buffer);
		cout<<buffer<<endl;
		if(tmpstring.Contains(namesofSetting[0]))
			CutRange[0]=(tmpstring.Remove(0,namesofSetting[0].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[1]))
			CutRange[1]=(tmpstring.Remove(0,namesofSetting[1].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[2]))
			FitRange[0]=(tmpstring.Remove(0,namesofSetting[2].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[3]))
			FitRange[1]=(tmpstring.Remove(0,namesofSetting[3].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[4]))
			minptformaterial[0]=(tmpstring.Remove(0,namesofSetting[4].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[5]))
			maxptformaterial[0]=(tmpstring.Remove(0,namesofSetting[5].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[6]))
			minptformaterial[1]=(tmpstring.Remove(0,namesofSetting[6].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[7]))
			maxptformaterial[1]=(tmpstring.Remove(0,namesofSetting[7].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[8]))
			minptformaterial[2]=(tmpstring.Remove(0,namesofSetting[8].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[9]))
			maxptformaterial[2]=(tmpstring.Remove(0,namesofSetting[9].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[10]))
			minptformaterial[3]=(tmpstring.Remove(0,namesofSetting[10].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[11]))
			maxptformaterial[3]=(tmpstring.Remove(0,namesofSetting[11].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[12]))
			minptformaterial[4]=(tmpstring.Remove(0,namesofSetting[12].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[13]))
			maxptformaterial[4]=(tmpstring.Remove(0,namesofSetting[13].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[14]))
			minptformaterial[5]=(tmpstring.Remove(0,namesofSetting[14].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[15]))
			maxptformaterial[5]=(tmpstring.Remove(0,namesofSetting[15].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[16]))
			minptforWD[0]=(tmpstring.Remove(0,namesofSetting[16].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[17]))
			maxptforWD[0]=(tmpstring.Remove(0,namesofSetting[17].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[18]))
			minptforWD[1]=(tmpstring.Remove(0,namesofSetting[18].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[19]))
			maxptforWD[1]=(tmpstring.Remove(0,namesofSetting[19].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[20]))
			minptforWD[2]=(tmpstring.Remove(0,namesofSetting[20].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[21]))
			maxptforWD[2]=(tmpstring.Remove(0,namesofSetting[21].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[22]))
			minptforWD[3]=(tmpstring.Remove(0,namesofSetting[22].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[23]))
			maxptforWD[3]=(tmpstring.Remove(0,namesofSetting[23].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[24]))
			minptforWD[4]=(tmpstring.Remove(0,namesofSetting[24].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[25]))
			maxptforWD[4]=(tmpstring.Remove(0,namesofSetting[25].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[26]))
			minptforWD[5]=(tmpstring.Remove(0,namesofSetting[26].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[27]))
			maxptforWD[5]=(tmpstring.Remove(0,namesofSetting[27].Length()+1)).Atof();				
		else if (tmpstring.Contains(namesofSetting[28]))
			fMaxContaminationPIDMC=(tmpstring.Remove(0,namesofSetting[28].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[29]))
			minRanges[0]=(tmpstring.Remove(0,namesofSetting[29].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[30]))
			maxRanges[0]=(tmpstring.Remove(0,namesofSetting[30].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[31]))
			minRanges[1]=(tmpstring.Remove(0,namesofSetting[31].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[32]))
			maxRanges[1]=(tmpstring.Remove(0,namesofSetting[32].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[33]))
			minRanges[2]=(tmpstring.Remove(0,namesofSetting[33].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[34]))
			maxRanges[2]=(tmpstring.Remove(0,namesofSetting[34].Length()+1)).Atof();		
	        else if (tmpstring.Contains(namesofSetting[35]))
			TOFPIDsignalmatching[0]=(tmpstring.Remove(0,namesofSetting[35].Length()+1)).Atof();		
		else if (tmpstring.Contains(namesofSetting[36]))
			TOFPIDsignalmatching[1]=(tmpstring.Remove(0,namesofSetting[36].Length()+1)).Atof();
		else if (tmpstring.Contains(namesofSetting[37]))
			TOFPIDsignalmatching[2]=(tmpstring.Remove(0,namesofSetting[37].Length()+1)).Atof();	
		else if (tmpstring.Contains(namesofSetting[38]))
			filenames[0]=tmpstring.Remove(0,namesofSetting[38].Length()+1);
		else if (tmpstring.Contains(namesofSetting[39]))
			filenames[1]=tmpstring.Remove(0,namesofSetting[39].Length()+1);
		else if (tmpstring.Contains(namesofSetting[40]))
			filenames[2]=tmpstring.Remove(0,namesofSetting[40].Length()+1);
		else 
			continue;


//	Double_t minRanges[3]={0.3,0.3,0.45};
//Double_t maxRanges[3]={1.5,1.2,2.2};
//Double_t fMaxContaminationPIDMC=0.2;



	}
	for(int i=0;i<6;i++)
		cout<<minptformaterial[i]<<" "<<maxptformaterial[i]<<" "<<minptforWD[i]<<" "<<maxptforWD[i]<<endl;
	cout<<FitRange[0]<<" "<<FitRange[1]<<" "<<CutRange[0]<<CutRange[1]<<endl;
	if(FitRange[0]>=FitRange[1])	
	{
		cout<<"A"<<endl;				
		return false;
	}
	if(CutRange[0]>=CutRange[1])
	{	
		cout<<"B"<<endl;				
		return false;
	}
	for(int i=0;i<6;i++)
	{
		if((minptformaterial[i]>maxptformaterial[i]&&minptformaterial[i]>0.0)||minptformaterial[i]<0.0||maxptformaterial[i]<0.0)
		{
			cout<<"C"<<endl;
			return false;
		}
		if((minptforWD[i]>maxptforWD[i]&&minptforWD[i]>0.0)||minptforWD[i]<0.0||maxptforWD[i]<0.0)
		{
			cout<<"D"<<endl;
			return false;
		}
	}
	for(int i=0;i<3;i++)
		if(minRanges[i]>maxRanges[i])
			return false;
		

	return true;
}

void SubHistWithFullCorr(TH1F* h1, TH1F* h2, Float_t factor1, Float_t factor2)
{
	if(h1->GetNbinsX()!=h2->GetNbinsX())
		return;
	for (int i=0;i<=h1->GetNbinsX();i++)
	{
		Float_t tmpvalue=factor1*h1->GetBinContent(i)-factor2*h2->GetBinContent(i);
		Float_t tmperror=TMath::Abs(factor1*factor1*h1->GetBinError(i)*h1->GetBinError(i)-factor2*factor2*h2->GetBinError(i)*h2->GetBinError(i));
		h1->SetBinContent(i,tmpvalue);
		h1->SetBinError(i,TMath::Sqrt(tmperror));
	}		
	
}

void TOFMatchingForNch(TH1* h)
{
	 if(TOFMatchingScalling[0]>0.0&&TOFMatchingScalling[1]>0.0)
	 {
		Float_t factor=0.5*TOFMatchingScalling[0]+0.5*TOFMatchingScalling[1];
		for(Int_t ibin=1;ibin<h->GetNbinsX();ibin++)
		{
			Float_t ptspectra=h->GetBinCenter(ibin);
			if(ptspectra<tcutsdata->GetPtTOFMatching())
				continue;
			h->SetBinContent(ibin,(h->GetBinContent(ibin)/factor));
		}

	 }
	  else
		return;			


}
void TOFPIDsignalmatchingApply(TH1* h, Float_t factor)
{
	if(factor<0.0)
		return;
	for(Int_t ibin=1;ibin<h->GetNbinsX();ibin++)
	{
		Float_t ptspectra=h->GetBinCenter(ibin);
		if(ptspectra<tcutsdata->GetPtTOFMatching())
			continue;
		h->SetBinContent(ibin,(h->GetBinContent(ibin)*factor));

	}

}
void CalculateDoubleCounts(TH1* doubleconunts,TH1F** rawspectra,Int_t ipar, Bool_t dataflag)
{
	TH2F* tmphist=0x0;	
	if (dataflag)
	 	tmphist=(TH2F*)managerdata->GetGenMulvsRawMulHistogram("hHistDoubleCounts");	
	else
		tmphist=(TH2F*)managermc->GetGenMulvsRawMulHistogram("hHistDoubleCounts");

	if(!tmphist)	
		return;
	TH1F* tmphist1=(TH1F*)rawspectra[ipar]->Clone("test");
	tmphist1->Add(rawspectra[ipar+3]);
	doubleconunts->Add(tmphist->ProjectionX(doubleconunts->GetName(),1,1));
	if(ipar!=2)
		doubleconunts->Add(tmphist->ProjectionX("pi+k",2,2));			
	if(ipar!=1)
		doubleconunts->Add(tmphist->ProjectionX("pi+p",3,3));
	if(ipar!=0)
		doubleconunts->Add(tmphist->ProjectionX("k+p",4,4));
	doubleconunts->Divide(doubleconunts,tmphist1,1,1,"B");

	delete tmphist1;
	

}

void CopyCorrectionFromFile(TString filename,TString correctionname,TH1F** corrtab)
{
	TFile* ftmp=TFile::Open(filename.Data());
	TString tmp("tmp");
	TList* ltmp=(TList*)ftmp->Get("output");

	for(Int_t icharge=0;icharge<2;icharge++)
        {
                for(Int_t ipart=0;ipart<3;ipart++)
                {
                        Int_t index=ipart+3*icharge;
                        tmp.Form("%s%s%s",correctionname.Data(),Particle[ipart].Data(),Sign[icharge].Data());
			TH1* histtmp=(TH1*)ltmp->FindObject(tmp.Data());
			if(!histtmp)
				continue;
			for(Int_t ibin=1;ibin<=histtmp->GetNbinsX();ibin++)
			{
				Float_t content=histtmp->GetBinContent(ibin);
				Float_t error=histtmp->GetBinError(ibin);
				Int_t bin=corrtab[index]->FindBin(histtmp->GetBinCenter(ibin));
				if(content>0.0)
					cout<<"TEST1 "<<content<<" "<<error<<endl;
				corrtab[index]->SetBinContent(bin,content);
				corrtab[index]->SetBinError(bin,error);
				
			}  
		}
	}
	delete ltmp;
	ftmp->Close();
}

TH1F* ReBinDCAHisto(TH1* h)
{
	TString name=h->GetName();
	name+="Rebin";
	Int_t kDCABins=88;
	Double_t binsDCADummy[]={-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0};
	
	TH1F* hout=new TH1F(name.Data(),";;dcaxy[cm];",kDCABins,binsDCADummy);
	Int_t j=1;
	Float_t sum=0.0;
	for (int i=1;i<=h->GetXaxis()->GetNbins();i++)
	{
		if(TMath::Abs(h->GetXaxis()->GetBinUpEdge(i)-hout->GetXaxis()->GetBinUpEdge(j))>0.0002)
		{
			sum+=h->GetBinContent(i);	
		}
		else
		{
			sum+=h->GetBinContent(i);
			hout->SetBinContent(j,sum);
			hout->SetBinError(j,TMath::Sqrt(sum));
			j++;
			sum=0.0;	
		}
	}
	hout->SetMarkerStyle(h->GetMarkerStyle());
	hout->SetMarkerSize(h->GetMarkerSize());
	hout->SetMarkerColor(h->GetMarkerColor());
	hout->SetLineStyle(h->GetLineStyle());
        //hout->SetLineSize(h->GetLineSize());
        hout->SetLineColor(h->GetLineColor());

	delete h;
	hout->Sumw2();
	return hout;

}
void RawYieldFromFits(AliSpectraBothHistoManager* hman,TH1F** histo,AliSpectraBothTrackCuts* cuts, TList* lfits,TH1F** primaryfractionfromfit)
{
	
	for(Int_t ipart=0;ipart<3;ipart++)
	{
		Float_t pttof=cuts->GetPtTOFMatchingPion();
		if(ipart==1)
			pttof=cuts->GetPtTOFMatchingKaon();
		if(ipart==2)
			pttof=cuts->GetPtTOFMatchingProton();	

		if(!hman->GetNSigHistogram(Form("hHistNSig%sPtTOF",Particle[ipart].Data())))
					continue;
		TH2F *nsigTOF= (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSig%sPtTOF",Particle[ipart].Data())))->Clone();
		TH2F *nsigTPC= (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSig%sPtTPC",Particle[ipart].Data())))->Clone();
		TH2F *nsigtrueTOF=0x0;
		TH2F *nsigtrueTPC=0x0;
		TString nametmp(histo[0]->GetName());
		if(nametmp.Contains("MC"))
		{
			nsigtrueTOF= (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSigTrue%sPtTOF",Particle[ipart].Data())))->Clone();
			nsigtrueTPC= (TH2F*)((TH2F*)hman->GetNSigHistogram(Form("hHistNSigTrue%sPtTPC",Particle[ipart].Data())))->Clone();

		}
		
		for(Int_t icharge=0;icharge<2;icharge++)
		{
			Int_t index=ipart+3*icharge;
			RawYieldFromFitsForParticleType(nsigTOF,nsigTPC,histo[index],icharge,ipart,pttof,lfits,primaryfractionfromfit[index],nsigtrueTOF,nsigtrueTPC); 
			
					
		}
	} 

}			
void RawYieldFromFitsForParticleType(TH2F* nsigTOF,TH2F* nsigTPC,TH1F* histo,Int_t icharge, Int_t ipart,Float_t pttof,TList* lfits,TH1F* primaryfractionfromfit,TH2F* nsigtrueTOF,TH2F* nsigtrueTPC)
{
	Int_t  charge =(icharge==0?1:-1);
	Float_t ptmin=minRanges[ipart]*charge;
	Float_t ptmax=maxRanges[ipart]*charge;
	Int_t startbin=nsigTOF->GetXaxis()->FindBin(ptmin+charge*0.0001);
	Int_t stopbin=nsigTOF->GetXaxis()->FindBin(ptmax-charge*0.0001);
	
	for (int i=startbin;i!=(stopbin+charge);i=i+charge)	
	{
		cout<<i<<" "<<icharge<<" "<<ipart<<" "<<histo->GetTitle()<<endl;
		TH1F *nsig_data_Proj1=0x0;
		TH1F *nsigtrue_data_Proj1=0x0;
		if(pttof>TMath::Abs(nsigTOF->GetXaxis()->GetBinCenter(i)))
		{	
			nsig_data_Proj1=(TH1F*)nsigTPC->ProjectionY(Form("%s[%.2f,%.2f] %s",nsigTPC->GetName(),nsigTPC->GetXaxis()->GetBinLowEdge(i),nsigTPC->GetXaxis()->GetBinUpEdge(i),histo->GetName()),i,i);
			if(nsigtrueTPC)
			{
				nsigtrue_data_Proj1=(TH1F*)nsigtrueTPC->ProjectionY(Form("%s[%.2f,%.2f] %s",nsigtrueTPC->GetName(),nsigtrueTPC->GetXaxis()->GetBinLowEdge(i),nsigtrueTPC->GetXaxis()->GetBinUpEdge(i),histo->GetName()),i,i);
				lfits->Add(nsigtrue_data_Proj1);
			}	
		}
		else
		{
			nsig_data_Proj1=(TH1F*)nsigTOF->ProjectionY(Form("%s[%.2f,%.2f] %s",nsigTOF->GetName(),nsigTOF->GetXaxis()->GetBinLowEdge(i),nsigTOF->GetXaxis()->GetBinUpEdge(i),histo->GetName()),i,i);
			if(nsigtrueTOF)
			{
				nsigtrue_data_Proj1=(TH1F*)nsigtrueTOF->ProjectionY(Form("%s[%.2f,%.2f] %s",nsigtrueTOF->GetName(),nsigtrueTOF->GetXaxis()->GetBinLowEdge(i),nsigtrueTOF->GetXaxis()->GetBinUpEdge(i),histo->GetName()),i,i);
				lfits->Add(nsigtrue_data_Proj1);
			}	
		}

		TF1* f=new TF1(Form("FitFun_%d_%d_%d",icharge,ipart,i),"gausn",-3,3);
		nsig_data_Proj1->Fit(Form("FitFun_%d_%d_%d",icharge,ipart,i),"R");
		
		lfits->Add(nsig_data_Proj1);
		lfits->Add(primaryfractionfromfit);
		Float_t yield=f->GetParameter(0)/ nsig_data_Proj1->GetXaxis()->GetBinWidth(1);
		Float_t yielderror=f->GetParError(0)/nsig_data_Proj1->GetXaxis()->GetBinWidth(1);
		Int_t bintofill=histo->GetXaxis()->FindBin(TMath::Abs(nsigTOF->GetXaxis()->GetBinCenter(i)));
		cout<<yield<<" "<<yielderror<<" "<<nsig_data_Proj1->GetXaxis()->GetBinWidth(1)<<" "<<bintofill<<endl;
		histo->SetBinContent(bintofill,yield);
		histo->SetBinError(bintofill,yielderror);
		
		Float_t fract=yield/nsig_data_Proj1->Integral(nsig_data_Proj1->GetXaxis()->FindBin(-2.999),nsig_data_Proj1->GetXaxis()->FindBin(2.999));
		Float_t fracterror=fract*yielderror/yield;
		primaryfractionfromfit->SetBinContent(bintofill,fract);
		primaryfractionfromfit->SetBinError(bintofill,fracterror);


		TVirtualFitter*  fitter=TVirtualFitter::GetFitter();
		delete fitter;
	}
}
void Rescallesecondarycontaimation(TH1F* contfit, TH1F* contWDfit, TH1F* contMatfit, TH1F* primaryfit, TH1F* primaryfractionfromfit)
{
		for (int i=1;i<=contfit->GetXaxis()->GetNbins();i++)
		{
			Float_t vcontfit=contfit->GetBinContent(i);
			Float_t vcontWDfit=contWDfit->GetBinContent(i);
			Float_t vcontMatfit=contMatfit->GetBinContent(i);
			Float_t vprimaryfit=primaryfit->GetBinContent(i);
			Float_t vprimaryfractionfromfit=primaryfractionfromfit->GetBinContent(i);
			if(vcontfit>0.0)
			{
				vprimaryfit=vprimaryfractionfromfit*vprimaryfit;
				vcontfit=vcontfit/(vcontfit+vprimaryfit);
				vcontWDfit=vcontWDfit/(vcontfit+vprimaryfit);
				vcontMatfit=vcontWDfit/(vcontfit+vprimaryfit);
				
				contfit->SetBinError(i,contfit->GetBinError(i)*vcontfit/contfit->GetBinContent(i));
				contfit->SetBinContent(i,vcontfit);
				if(contWDfit->GetBinContent(i)>0.0)
					contWDfit->SetBinError(i,contWDfit->GetBinError(i)*vcontWDfit/contWDfit->GetBinContent(i));
				contWDfit->SetBinContent(i,vcontWDfit);
				if(contMatfit->GetBinContent(i)>0.0)
					contMatfit->SetBinError(i,contMatfit->GetBinError(i)*vcontMatfit/contMatfit->GetBinContent(i));
				contMatfit->SetBinContent(i,vcontMatfit);
				if(primaryfit->GetBinContent(i)>0.0)
					primaryfit->SetBinError(i,primaryfit->GetBinError(i)*vprimaryfit/primaryfit->GetBinContent(i));
				primaryfit->SetBinContent(i,vprimaryfit);

			}

		}
}
