/*************************************************************************
 *
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

//==================================================================
// Simple class for di-charged jet analyses.
// by Beomkyu KIM
//==================================================================
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TList.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TVector.h>
#include <TSystem.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TRandom.h>
#include "THistManager.h"
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskDijet.h"
#include <TRandom3.h>
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliEmcalPythiaInfo.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"


const Double_t pionmass = AliPID::ParticleMass(AliPID::kPion);
const Double_t pi = TMath::Pi();

AliAnalysisTaskDijet::AliBSDiJetTask() 
  : AliAnalysisTaskEmcalJet("AliAnalysisTaskDijet",kTRUE)
    , fOutput(0)
{
  DefineOutput (1, TList::Class());
}

  AliAnalysisTaskDijet::AliBSDiJetTask(const char *name) 
  : AliAnalysisTaskEmcalJet(name,kTRUE)
  , fOutput(0)
    , fhJetPt(0)
{
  DefineOutput (1, TList::Class());
}
  AliAnalysisTaskDijet::AliBSDiJetTask(const char *name, const char *option) 
  : AliAnalysisTaskEmcalJet(name,kTRUE)
  , fOutput(0)
  , fOption(option)
    , fhJetPt(0)
{
  DefineOutput (1, TList::Class());
}

  AliAnalysisTaskDijet::AliBSDiJetTask(const AliBSDiJetTask& ap) 
  : AliAnalysisTaskEmcalJet(ap.fName,kTRUE)
  , fOutput(ap.fOutput)
    , fhJetPt(ap.fhJetPt)
{
}
AliAnalysisTaskDijet& AliBSDiJetTask::operator = (const AliBSDiJetTask& ap)
{

  this->~AliAnalysisTaskDijet();
  new(this) AliAnalysisTaskDijet(ap);
  return *this;
}

AliAnalysisTaskDijet::~AliBSDiJetTask()
{
  delete fOutput;
  delete fhJetPt; 
  delete fBSRandom;
}

void AliAnalysisTaskDijet::UserCreateOutputObjects(){
  fBSRandom = new TRandom3;
  fBSRandom->SetSeed();
  //===============================
  // BINS
  //===============================
  auto binJetMass    = AxisFix( "JetMass", 10, 0, 500 );
  auto binInvM       = AxisFix( "InvM", 30,0,300);
  auto binCent       = AxisFix( "Cent", 1, 0, 100 );
  auto binJetPtCut   = AxisVar( "JetPtCut"
	,{0, 5, 10, 15, 20, 30, 50, 70, 100, 150, 200,  1000}  );
  auto binDiJetSel   = AxisStr( "DiJetSel", { "LS","LS-PI2","L-PI2-S" } );
  fNDiJetSelection   = binDiJetSel.GetNbins();
  auto binType       = AxisStr( "Type", { "DATA","Mixing-Simple" } );
  fNType             = binType.GetNbins();
  if( fIsAA ){
    binCent = AxisFix( "Cent", 10, 0, 100 );
  }
  auto binLog1k     = AxisLog("Log1k",500,0.1,1000,0);
  auto binLog3c     = AxisVar("Log3c",{0,5,7,9,12,16,21,28,36,45,57,70,85,99,115,132,150,169,190,212,235,500});
  auto binAsim      = AxisFix("Asim", 100 ,0,1);


  //===============================
  // HISTOGRAMS
  //===============================
  fHistos = new THistManager("jethists");

  const int nbins=100;
  Double_t logbins[nbins+1];
  Double_t low= 0.1;
  Double_t high=300;
  Double_t logbw= (log(high)-log(low))/nbins;
  for(int ij=0;ij<=nbins;ij++) logbins[ij]=low*exp(ij*logbw);

  fHistos -> CreateTH1("hJetMulti","",101,-1,100);
  fHistos -> CreateTH1("jetpt","jetpt",nbins,logbins);
  fHistos -> CreateTH1("jetptresum","jetptresum",nbins,logbins);
  fHistos -> CreateTProfile("jetptresumRatio","jetptresumRatio",nbins,logbins);
  fHistos -> CreateTH1("jetinv","jetinv",nbins,logbins);
  fHistos -> CreateTH1("zvtx","zvtx",60,-30,30);
  fHistos -> CreateTH1("nCentrality","nCentrality", 10,0,100);
  fHistos -> CreateTH2("trketaphi","trketaphi",20,-1,1,90,0,TMath::TwoPi());

  CreateTHnSparse ("hPt","",2,{binCent, binLog3c},"s");
  CreateTHnSparse( "hJetPtLPt","",3,{binCent
      ,AxisLog("JetPt",500,0.1,100,0)
      ,AxisLog("LetPt",500,0.1,100,0)
      }, "s" );
  CreateTHnSparse( "hJetPtLeading","",5,{
      binType, binDiJetSel,binCent,binJetPtCut,binLog3c}, "s" );
  CreateTHnSparse( "hJetPtSubLeading", "", "hJetPtLeading" , "s" );
  CreateTHnSparse( "hJetMassLeading", "", "hJetPtLeading" , "s" );
  CreateTHnSparse( "hJetMassSubLeading", "", "hJetPtLeading" , "s" );
  CreateTHnSparse( "hDiJetInvM", "DiJetMass", 5, {
      binType, binDiJetSel, binCent, binJetPtCut, binLog3c }, "s" );
  CreateTHnSparse( "hDiJetInvMRes", "DiJetMass Res", 5, {
      binDiJetSel, binCent, binJetPtCut, binLog3c, binLog3c }, "s" );
  CreateTHnSparse( "hDiJetDPhi", "DiJet #Delta#Phi", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, AxisFix( "",100, -pi, pi ) },"s");
  CreateTHnSparse( "hDiJetDPhi_0_2pi", "DiJet #Delta#Phi", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, AxisFix( "",100, 0, 2*pi ) },"s");
  CreateTHnSparse( "hDiJetPtPair", "DiJet PtPair", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, binLog3c  },"s");
  CreateTHnSparse( "hDiJetPtPairKine", "DiJet PtPair", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, binLog3c  },"s");
  CreateTHnSparse( "hDiJetPtAsim", "DiJet PtPair", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, binAsim  },"s");
  CreateTHnSparse( "hDiJetEAsim", "DiJet PtPair", 6, {
      binType, binDiJetSel, binCent, binJetPtCut, binInvM, binAsim  },"s");
  CreateTHnSparse( "hDiJetPtPairRes", "DiJet PtPair Res Matrix", 6, {
      binDiJetSel, binCent, binJetPtCut, binInvM, binLog3c, binLog3c  },"s");
  CreateTHnSparse( "hDiJetInvMPtPairRes", "DiJet InvM PtPair Res Matrix", 7, {
      binDiJetSel, binCent, binJetPtCut, binInvM,binInvM, binLog3c, binLog3c  },"s");

  PostData(1, fHistos->GetListOfHistograms());

  //===============================
  // For Sure
  //===============================
  std::cout<< "DEBUG4 IsAA?"<< (fIsAA?"AA":"pp")<<std::endl;
  std::cout<<"NBins of Cent : "<<binCent.GetNbins()<<"\t"<<binCent.GetXmin()<<"\t"<<binCent.GetXmax()<<endl;
  if( fNDiJetSelection != kBDiJetSelEnd-1 ){
    cout<<"fNDiJetSelection("<<fNDiJetSelection
      <<") is not match with kBDiJetSelEnd("<<kBDiJetSelEnd<<")"<<endl;
    gSystem->Exit(1);
  }
  if( fNType != kBTypeEnd-1 ){
    cout<<"fNType("<<fNType
      <<") is not match with kBTypeEnd("<<kBTypeEnd<<")"<<endl;
    gSystem->Exit(1);
  }
  fDijetPtPair ={0,0,0}; 
  fDijetInvM ={0,0,0}; 
  fDijetapt ={0,0,0};
  invmscale = new TF1("invmscale","7.44183e-01+5.55996e-03*x",0,250) ;
  
}

//________________________________________________________________________

Bool_t AliAnalysisTaskDijet::Run(){
	using TMath::Abs;

	AliCentrality *Cent = InputEvent()->GetCentrality();
	if( ! Cent ) return false;
	fCent = Cent->GetCentralityPercentile("V0M");

	//TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
	//if(!tree) return kFALSE;

	//TFile *curfile = tree->GetCurrentFile();

	//if(!curfile) return kFALSE;

	//TString fCurrFileName = TString(curfile->GetName());
	//Float_t fXsec,fTrials;
	//Int_t pthard;
	//PythiaInfoFromFile(	fCurrFileName.Data() , fXsec, fTrials, pthard);
  //cout<<"fCurrFileName : "<<fCurrFileName.Data()<<endl;
 	//cout<<"cross section : "<<fXsec<<endl;
  //cout<<"ntrials : "<<fTrials<<endl;
	//Double_t sf = 1./fTrials;
  //cout<<fXsec<<" "<<fTrials<<" "<<pthard<<endl;

	//https://twiki.cern.ch/twiki/bin/viewauth/ALICE/JetMCProductionsCrossSections
  Double_t sf = 1.;
	Int_t NTrials = -1;
	Double_t XSection =-1;
  TLorentzVector p6,p7,p67;
	if(fIsMC){
		if (fOption.Contains("AOD")){
			AliVEvent *event = InputEvent();
			if (!event) {
				Printf("ERROR: Could not retrieve event");
				return kFALSE;
			}
			AliAODMCHeader *cHeaderAOD  = dynamic_cast<AliAODMCHeader*>
				(event->FindListObject(AliAODMCHeader::StdBranchName()));
			TList *genHeaders         = cHeaderAOD->GetCocktailHeaders();
			NTrials = -1;
			XSection =-1;
			AliGenEventHeader* gh       = 0;
			for(Int_t i = 0; i<genHeaders->GetEntries();i++){
				gh   = (AliGenEventHeader*)genHeaders->At(i);
				TString GeneratorName   = gh->GetName();
				if (GeneratorName.CompareTo("Pythia") == 0){
					AliGenPythiaEventHeader* gPythia = dynamic_cast<AliGenPythiaEventHeader*>(gh);
					NTrials = gPythia->Trials();
					XSection = gPythia->GetXsection();
				}
			}
			sf = XSection/NTrials;
		}
    else {
			// ----------------------------------------------------------------------
			// Pointer to a MC event-------------------------------------------------
			AliMCEvent *mcEvent = MCEvent();
			AliStack *stack = mcEvent->Stack();
			// ----------------------------------------------------------------------
      AliGenPythiaEventHeader*  gPythia =
				dynamic_cast<AliGenPythiaEventHeader*>(mcEvent->GenEventHeader());
			NTrials = gPythia->Trials();
			XSection = gPythia->GetXsection();
			sf = XSection/NTrials;
      Int_t nPrim  = stack->GetNprimary();
			for (Int_t i = 0; i < nPrim; i++){
				TParticle* p = stack->Particle(i);
        if (!p) continue;
			}
      //FillTHnSparse( "hDiJetPtPairKine",{ 1, 1, 0, spt, p67.M(),p67.Pt()  }, 1/p67.Pt()*sf); 
    }
	}
 


  // Fill z_vertex and cut z_vertex range
  //const AliVVertex* vtx = InputEvent()->GetPrimaryVertexTPC();
  const AliVVertex* vtx = InputEvent()->GetPrimaryVertex();
  //const AliVVertex* vtx = InputEvent()->GetPrimaryVertexSPD();
  if (vtx->GetNContributors()<1) return kFALSE;
  Double_t zvtx = vtx->GetZ();
  fHistos -> FillTH1 ("zvtx",zvtx);
  if (Abs(zvtx)>10) return kFALSE;
  fHistos -> FillTH1 ("nCentrality",fCent);

	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	AliAnalysisTaskDijet *task = (AliBSDiJetTask*)man -> GetTask("DiJetV5x20_LPt5_kANY_kinech");

  //Double1D mctrueptpair;
	//if (task) mctrueptpair = task ->GetDijetPtPair();

  //===============================
  // RESUM JETS
  //===============================

  auto jetContainer = GetJetContainer(0);
  auto trkContainer = GetParticleContainer(0);
  

  // =============================
  // QA plots
  // =============================
  for (auto trk : trkContainer->accepted()){
    fHistos->FillTH2("trketaphi",trk->Eta(),trk->Phi());
    FillTHnSparse( "hPt", {fCent, trk->Pt()}, 1./trk->Pt());
  }

  //=== NOTE : We will not use sj[0] because bin in THnSparse is always begin with 1.
  TLorentzVector2D sj( fNDiJetSelection+1, TLorentzVector1D(2) );
  TLorentzVector1D sumJets;
  std::vector<Double_t> lpts;
  int ij=0;
  for( auto j : jetContainer->accepted() ){
    //fhJetPt->Fill( fCent, j->Pt() );

    TLorentzVector sum (0,0,0,0);
    double sumpt=0;
    double lpt = 0;
    for( int it=0; it<j->GetNumberOfTracks(); it++ ) {
      auto trk =  j->TrackAt(it,trkContainer->GetArray());
      //cout<<"track selected : "<<aodsel->IsSelected(trk)<<endl;
      TLorentzVector temp;
      temp.SetXYZM(trk->Px(),trk->Py(),trk->Pz(),pionmass);
      sum+=temp;
      sumpt+=trk->Pt();
      if( lpt < temp.Pt() )  lpt = temp.Pt();
    }
    //=== Jet CUT
    if( Abs(sum.Eta())>0.4 ) continue;
    if( lpt < fLeadingParticlePtMin ) continue;

    //=== Store Selected Jet
    sumJets.push_back(sum);
    lpts.push_back(lpt);

    //=== DiJetSelection 1 : LS = Leading-SubLeading
    if (sj[kBLS][0].Pt() < sum.Pt() ) sj[kBLS][0]=sum;
    else if( sj[kBLS][1].Pt() < sum.Pt() ) sj[kBLS][1]=sum;

    //==== FILL
    fHistos->FillTH1("jetptresum",sumpt,sf);
    fHistos->FillProfile("jetptresumRatio",sumpt, j->Pt()/sumpt);
    fHistos->FillTH1("jetinv",sum.M());
    FillTHnSparse("hJetPtLPt",{fCent,sum.Pt(),lpt},1./lpt);

    ij++;
  }

  //=== DiJetSelection 2 : LS-PI2 = Leading-Subreading & dPhi>pi/2
  if( Abs(sj[kBLS][0].DeltaPhi(sj[kBLS][1])) > pi/2. ) 
    sj[kBLS_PI2] = sj[kBLS];
  else {
    sj[kBLS_PI2][0] = TLorentzVector(0,0,0,0);
    sj[kBLS_PI2][1] = TLorentzVector(0,0,0,0);
  }

  //=== DiJetSelection 3 : L-PI2-S = Leading - Subleading in opposite hemi-spear
  auto tj = sj[kBL_PI2_S][0] = sj[kBLS][0]; // Leading is same as 1:LS
  sj[kBL_PI2_S][1] = TLorentzVector(0,0,0,0);
  for( auto jet: sumJets ){
    if( Abs(tj.DeltaPhi( jet )) < pi/2. ) continue; 
    if( sj[kBL_PI2_S][1].Pt() < jet.Pt() ) sj[kBL_PI2_S][1] = jet;
  }

  //==== FILL
  fHistos->FillTH1("hJetMulti", sumJets.size() );

  //===============================
  // ptPair : leading - subleading
  //===============================

  fDijetPtPair[0] = 1;
  fDijetPtPair[1] = 1;
  fDijetPtPair[2] = 1;
  fDijetInvM[0] = 0;
  fDijetInvM[1] = 0;
  fDijetInvM[2] = 0;
  fDijetapt[0] = -1;
  fDijetapt[1] = -1;
  fDijetapt[2] = -1;
  for( int itype=kBTypeBegin;itype<kBTypeEnd;itype++ ){
    for( int ids=kBDiJetSelBegin;ids<kBDiJetSelEnd;ids++ ){ // NOTE : begin with 1
      auto j = sj[ids];

      //=== SKIP Empty DiJet
      if( j[0].E() < 1e-4 || j[1].E() < 1e-4 ) continue;
      Double_t nphi=0;
      if( itype == kBMixingSimple ){   // MixingSimple
        //== rotate j[1] with random phi and eta(-0.4~0.4)
        if( ids == kBLS )    // 1:LS
          nphi = fBSRandom->Uniform(-pi, pi );
        else              // others
          nphi = fBSRandom->Uniform( pi/2., pi*3./2. ) + j[0].Phi();

        auto ntheta = fBSRandom->Uniform(2*TMath::ATan( TMath::Exp( -0.4 )), 2*TMath::ATan( TMath::Exp( 0.4 )));
        TVector3 nv;nv.SetPtThetaPhi( 1,ntheta, nphi ); 
        auto pratio = j[1].P()/nv.Mag();
        j[1].SetXYZM( nv.Px()*pratio, nv.Py()*pratio, nv.Pz()*pratio, j[1].M() );
      }
      auto type     = Double_t( itype );
      auto diJetSel = Double_t( ids );
      auto dijet    = j[0] + j[1];
      auto invM     = dijet.M();
      auto dPhi     = j[0].DeltaPhi(j[1]);
      auto dPhiA    = Abs(dPhi);
      auto dPhi_0_2pi  = TVector2::Phi_0_2pi(dPhi);
      auto tpt      = j[0].Pt();
      auto apt      = j[1].Pt();
      auto ptpair   = dijet.Pt(); 
      auto testdphi = nphi-j[0].Phi();
      auto ptAsim   = (tpt-apt)/(tpt+apt);
      auto eAsim    = (j[0].E()-j[1].E())/(j[0].E()+j[1].E());
			if (itype == kBData) { 
 				fDijetPtPair[ids-1] = ptpair;
        fDijetInvM[ids-1]=invM;
        fDijetapt[ids-1]=apt;
 			}
			Double_t massscale = invmscale->Eval(invM);
			if (itype == kBData && task){	
        Double_t truthptpair =  (task->GetDijetPtPair()).at(ids-1);
        Double_t truthinvM =  (task->GetDijetInvM()).at(ids-1);
				FillTHnSparse( "hDiJetPtPairRes",  {diJetSel,fCent, apt, invM,ptpair, truthptpair}); 
				FillTHnSparse( "hDiJetInvMRes", { diJetSel, fCent, apt, invM, truthinvM });
				FillTHnSparse( "hDiJetInvMPtPairRes", { diJetSel,fCent,apt,invM,truthinvM,ptpair,truthptpair});
			}
			FillTHnSparse( "hJetPtLeading",      { type, diJetSel, fCent, apt, tpt},sf/tpt );
			FillTHnSparse( "hJetPtSubLeading",   { type, diJetSel, fCent, apt, apt},sf/apt);
			FillTHnSparse( "hJetMassLeading",    { type, diJetSel, fCent, apt, tpt},sf/tpt );
			FillTHnSparse( "hJetMassSubLeading", { type, diJetSel, fCent, apt, apt},sf/apt );
			FillTHnSparse( "hDiJetInvM",         { type, diJetSel, fCent, apt, invM },sf);
			FillTHnSparse( "hDiJetDPhi",         { type, diJetSel, fCent, apt, invM, dPhi },sf);
			FillTHnSparse( "hDiJetDPhi_0_2pi",   { type, diJetSel, fCent, apt, invM, dPhi_0_2pi },sf);
			FillTHnSparse( "hDiJetPtPair",       { type, diJetSel, fCent, apt, invM, ptpair },sf/ptpair); 
			FillTHnSparse( "hDiJetPtAsim",       { type, diJetSel, fCent, apt, invM, ptAsim },sf/ptAsim);
			FillTHnSparse( "hDiJetEAsim",        { type, diJetSel, fCent, apt, invM, eAsim  },sf);

    }
  }

  PostData(1, fHistos->GetListOfHistograms());;
  //PostData(1, fOutput);
  return kTRUE;
}

THnSparse * AliAnalysisTaskDijet::CreateTHnSparse(TString name, TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}

THnSparse * AliAnalysisTaskDijet::CreateTHnSparse(TString name, TString title, TString templ, Option_t * opt){
  auto o = fHistos->FindObject(templ);
  if( !o ) {
    cout<<"ERROR: no "<<templ<<endl;
    gSystem->Exit(1);
  }
  auto ht = dynamic_cast<THnSparse*>( o ); 
  const TAxis * axises[ht->GetNdimensions()];
  for( int i=0;i<ht->GetNdimensions();i++ ) axises[i]= ht->GetAxis(i);
  auto h= fHistos->CreateTHnSparse(name, title, ht->GetNdimensions(), axises,opt );
  return h;
}

Long64_t AliAnalysisTaskDijet::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}

/*
   Long64_t AliAnalysisTaskDijet::Fill( TString name, std::vector<Double_t> x, Double_t w ){
   return FillTHnSparse( name, x, w );
   }
   */

Long64_t AliAnalysisTaskDijet::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}

TAxis AliAnalysisTaskDijet::AxisFix( TString name, int nbin, Double_t xmin, Double_t xmax ){ 
  TAxis axis(nbin, xmin, xmax);axis.SetName(name);
  return axis; 
}

TAxis AliAnalysisTaskDijet::AxisStr( TString name, std::vector<TString> bin ){ 
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax; 
}

TAxis AliAnalysisTaskDijet::AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}

TAxis AliAnalysisTaskDijet::AxisLog( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}
