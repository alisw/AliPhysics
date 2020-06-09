#include "TFile.h"
#include "TSystem.h"
#include "TGrid.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliCentrality.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliVMultiplicity.h"
#include "AliVVZERO.h"
#include "AliJetContainer.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <TVector3.h>
#include <TVectorT.h>
#include "AliJJet.h"
#include "AliAnalysisTaskGenMCRidge.h"
#include <AliDirList.h>
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalContainer.h"
#include "AliStack.h"
#include "AliEmcalTrackSelection.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalJet.h"
#include <TClonesArray.h>
#include <TList.h>
#include <TProfile.h>
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskRhoSparse.h"
#include "AliRhoParameter.h"
#include "AliGenEventHeader.h"
#include "AliPWG0Helper.h"
#include "AliAnalysisDataContainer.h"
#include <TDatabasePDG.h>
#include "AliJCDijetAna.h"

using namespace std;

const Double_t pi = TMath::Pi();



//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge()
:AliAnalysisTaskEmcalJet("AliAnalysisTaskGenMCRidge"),
 fOption(), fEMpool(),
 	fPtHardMin(0),
	fPtHardMax(0),
    TagThisEvent(),
 	fJFJTask(NULL),
	fJFJTaskName("JFJTask")
{
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge
(     
      const char *name
    , const char *option
)
:AliAnalysisTaskEmcalJet(name),
 fOption(option), fEMpool(),
 	fPtHardMin(0),
	fPtHardMax(0),
	TagThisEvent(),
 	fJFJTask(NULL),
	fJFJTaskName("JFJTask")
{
    DefineOutput (1, AliDirList::Class());
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge
(
      const AliAnalysisTaskGenMCRidge& ap
)
:fOption(ap.fOption), fEMpool(ap.fEMpool),
   fPtHardMin(ap.fPtHardMin),
   fPtHardMax(ap.fPtHardMax),
   TagThisEvent(),
   fJFJTask(ap.fJFJTask),
   fJFJTaskName(ap.fJFJTaskName)
{
    DefineOutput (1, AliDirList::Class());
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge& AliAnalysisTaskGenMCRidge::operator =
(
      const AliAnalysisTaskGenMCRidge& ap
)
{
    // assignment operator

    DefineOutput (1, AliDirList::Class());
    this->~AliAnalysisTaskGenMCRidge();
    new(this) AliAnalysisTaskGenMCRidge(ap);
    return *this;
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge::~AliAnalysisTaskGenMCRidge()
{
    delete fOutput;
}

void AliAnalysisTaskGenMCRidge::UserCreateOutputObjects()
{
	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	// Add JFJ fastjet standalone task.
	fJFJTask = (AliJFJTask*)(man->GetTask( fJFJTaskName));

	fOutput = new AliDirList();
	fOutput->SetOwner();

	fHistos = new THistManager("Ridgehists");
	
	Double1D varmultbin = {0, 35, 80, 105, 1000};
	Double1D ptbin = {0.1,0.2,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,14.0,100};
	Double1D ltpttrackbin = {0.2, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 13.0, 20.0};
	Double1D jetptbin = {0, 10, 20, 30, 40, 50, 60, 80, 100, 1e5 };

	Double1D LongRangeEta = {
	-1.8,-1.7,-1.6,-1.5,
	 1.5, 1.6, 1.7, 1.8};

	Double1D verzbin = {-10,-8,-6,-4,-2,0,2,4,6,8,10};

	binPhi = AxisFix("phi",32,-0.5*pi-pi/32.0, 1.5*pi-pi/32.0);
	binEta = AxisFix("eta",80,-4.0,4.0);
	if( fOption.Contains("OnlyOneDim") ){
		binEta = AxisVar("eta",LongRangeEta);
	}

	binCent = AxisVar("Cent",varmultbin);
	binTPt = AxisVar("TPt",ptbin);
        binAPt = AxisVar("APt",ptbin);
	binLHPt = AxisVar("LHPT",ltpttrackbin);
	binJetPt = AxisVar("JetPT",jetptbin);
	binNtrig = AxisFix("Ntrig",1,0.5,1.5);
	binZ = AxisVar("Z",verzbin);

	CreateTHnSparse("hNTrigLH","hNTrigLH",4,{binCent,binTPt,binNtrig,binLHPt},"s");
	CreateTHnSparse("hRidgeLH","hRidgeLH",6,{binCent,binPhi,binEta,binTPt,binAPt,binLHPt},"s");
	CreateTHnSparse("hRidgeMixingLH","hRidgeMixingLH",6,{binCent,binPhi,binEta,binTPt,binAPt,binLHPt},"s");

	CreateTHnSparse("hNTrigJet","hNTrigJet",4,{binCent,binTPt,binNtrig,binJetPt},"s");
	CreateTHnSparse("hRidgeJet","hRidgeJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetPt},"s");
	CreateTHnSparse("hRidgeMixingJet","hRidgeMixingJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetPt},"s");

	fHistos->CreateTH1("hNChargedMultiplixity","hNChargedMultiplixity",1000,0,1000);

	fEMpool.resize(binCent.GetNbins(),vector<eventpool> (binZ.GetNbins()));

	fOutput -> Add( fHistos->GetListOfHistograms() );
	PostData(1, fOutput);

}
//___________________________________________________________________
void AliAnalysisTaskGenMCRidge::Exec(Option_t *)
{
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	if( fMcHandler ) fMcEvent = (AliMCEvent*)fMcHandler->MCEvent();
	else{ return; }

	if( !fMcEvent ) return;
	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if( !fStack ) return;

	fMult = 0.0;
	fMult = this->GetMultiplicity( fStack );
	fHistos -> FillTH1("hNChargedMultiplixity",fMult,1);
	
	fZ = 0.0;
	fLHPt = 0.0;
	fJetPt = 0.0;
	for(int iE=0;iE<5;iE++) TagThisEvent[iE]=kFALSE; // init

	if( this->GetProperTracks( fStack ) ) this->GetCorrelations();

	PostData(1, fOutput);
}

double AliAnalysisTaskGenMCRidge::GetMultiplicity( AliStack* stack ){
 double ntrk = 0;
 for(int i=0;i<stack->GetNprimary();i++){
	TParticle* track = stack->Particle(i);
	if( !track ) continue;
	if( !AliPWG0Helper::IsPrimaryCharged(track, stack->GetNprimary()) ) continue;

	Int_t pdg = track->GetPdgCode();
	Double_t chp = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
	if( chp == 0 ) continue;
	if( !( ( track->Eta() >  2.8 && track->Eta() <  5.1 ) ||
	       ( track->Eta() > -3.7 && track->Eta() < -1.7 ) ) ) continue;

	ntrk++;
 }
 return ntrk; 
}

bool AliAnalysisTaskGenMCRidge::GetProperTracks( AliStack* stack ){

 tracklist *etl;
 eventpool *ep;

 if( binCent.FindBin( fMult ) >= 1 && binZ.FindBin( fZ ) >= 1 &&
	binCent.FindBin( fMult ) <= binCent.GetNbins() &&
	binZ.FindBin( fZ ) <= binZ.GetNbins() ){

	ep = &fEMpool[binCent.FindBin( fMult )-1][binZ.FindBin( fZ )-1];
	ep -> push_back( tracklist() ); //
	etl = &(ep->back());
 }

 goodTracks.clear();
 for(int i=0;i<stack->GetNprimary();i++){
	TParticle* track = stack->Particle(i);
	if( !track ) continue;
	if( !AliPWG0Helper::IsPrimaryCharged(track, stack->GetNprimary()) ) continue;

	Int_t pdg = track->GetPdgCode();
	Double_t chp = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
	if( chp == 0 ) continue;
	if( fabs( track->Eta() ) > 0.9 ) continue;

	goodTracks.push_back( (TParticle*)track -> Clone() );

	etl->push_back( (TParticle*)track -> Clone() );
 }

 if( !goodTracks.size() ) ep->pop_back();
 if( ep->size() > fbookingsize ){
	for (auto it: ep->front()) delete it;
	ep->pop_front();
 }

 return (bool)goodTracks.size();
}


void AliAnalysisTaskGenMCRidge::GetCorrelations(){
 NTracksPerPtBin.clear();
 NTracksPerPtBin.resize( binTPt.GetNbins() );

 for(int i=0;i<binTPt.GetNbins();i++){
	NTracksPerPtBin[i] = 0 ;
 }

 TParticle* track;
 for(int i=0;i<goodTracks.size();i++){
	track = (TParticle*)goodTracks.at(i);
	if( !track ) continue;
	if( fLHPt < track->Pt() ) fLHPt = track->Pt();
	if( binTPt.FindBin( track->Pt() )-1 >= 0 ){
		NTracksPerPtBin[ binTPt.FindBin( track->Pt() )-1 ] += 1.0;
	}
 }

 for(int i=0;i<binTPt.GetNbins();i++){
	FillTHnSparse("hNTrigLH",  {fMult,binTPt.GetBinCenter(i+1),1.0,fLHPt }, NTracksPerPtBin[i] );
	FillTHnSparse("hNTrigJet", {fMult,binTPt.GetBinCenter(i+1),1.0,fJetPt}, NTracksPerPtBin[i] );
 }

 TParticle* trackTrig;
 TParticle* trackAssoc;
 for(int i=0;i<goodTracks.size()-1;i++){
	trackTrig = (TParticle*)goodTracks.at(i);
	if( !trackTrig ) continue;
	for(int j=i+1;j<goodTracks.size();j++){
		trackAssoc = (TParticle*)goodTracks.at(j);
		if( !trackAssoc ) continue;

		double DeltaEta = trackTrig->Eta() - trackAssoc->Eta();
		double DeltaPhi = trackTrig->Phi() - trackAssoc->Phi();

		if( trackTrig->Pt() < trackAssoc->Pt() ){
			DeltaEta *= -1.0;
			DeltaPhi *= -1.0;
		}
		if( DeltaPhi > 1.5*pi ) DeltaPhi -= 2.0*pi;

		FillTHnSparse("hRidgeLH",  {fMult, DeltaPhi, DeltaEta, max(trackTrig->Pt(),trackAssoc->Pt()), min(trackTrig->Pt(),trackAssoc->Pt()), fLHPt},  1.0 );
		FillTHnSparse("hRidgeJet", {fMult, DeltaPhi, DeltaEta, max(trackTrig->Pt(),trackAssoc->Pt()), min(trackTrig->Pt(),trackAssoc->Pt()), fJetPt}, 1.0 );
	}
 }

 tracklist trackpool;

 int epsize=1;
 if( binCent.FindBin( fMult ) >= 1 && binZ.FindBin( fZ ) >= 1 &&
	binCent.FindBin( fMult ) <= binCent.GetNbins() &&
	binZ.FindBin( fZ ) <= binZ.GetNbins() ){

	eventpool &ep = fEMpool[binCent.FindBin( fMult )-1][binZ.FindBin( fZ )-1];
	epsize = ep.size();

	if (ep.size() < fbookingsize  ) return;
	int n = 0;
	for (auto pool: ep){
		if (n == (ep.size() -1 )) continue;
		for (auto track: pool) trackpool.push_back((TParticle*)track);
		n++;
	}
 }

 TParticle* trackMixing;
 for(int i=0;i<goodTracks.size();i++){
	trackTrig = (TParticle*)goodTracks.at(i);
	if( !trackTrig ) continue;
	for(int j=0;j<trackpool.size();j++){
		trackMixing = (TParticle*)trackpool.at(j);
		if( !trackMixing ) continue;

		double DeltaEta = trackTrig->Eta() - trackMixing->Eta();
		double DeltaPhi = trackTrig->Phi() - trackMixing->Phi();

		if( trackTrig->Pt() < trackMixing->Pt() ){
			DeltaEta *= -1.0;
			DeltaPhi *= -1.0;
		}
		if( DeltaPhi > 1.5*pi ) DeltaPhi -= 2.0*pi;


		FillTHnSparse("hRidgeMixingLH", {fMult, DeltaPhi, DeltaEta, max(trackTrig->Pt(),trackMixing->Pt()), min(trackTrig->Pt(),trackMixing->Pt()), fLHPt},  1.0 );
		FillTHnSparse("hRidgeMixingJet",{fMult, DeltaPhi, DeltaEta, max(trackTrig->Pt(),trackMixing->Pt()), min(trackTrig->Pt(),trackMixing->Pt()), fJetPt}, 1.0 );
	}
 }

}

void AliAnalysisTaskGenMCRidge::FinishTaskOutput()
{
// fOutput->Write();
}

void AliAnalysisTaskGenMCRidge::Terminate(Option_t*)
{

}


void AliAnalysisTaskGenMCRidge::RhoSparse(AliJetContainer *ktContainer, AliJetContainer * aktContainer , Int_t numberofexcludingjets) {
        AliEmcalJet *leading = nullptr;
        AliEmcalJet *subleading = nullptr;
        Int_t n = 0;

        n = 0;
        Int_t njetacc = 0;
        static Double_t rhovec[999];
        Double_t TotaljetAreaPhys=0;
        Double_t TotalTPCArea=2*TMath::Pi()*0.9;
        for (auto iBg : ktContainer->accepted_momentum())
        {
                if (n < numberofexcludingjets) {
                        n++;
                        continue;
                }
                auto bgjet = iBg.second;

                Bool_t matched = false;
                for (auto iSg : aktContainer->accepted_momentum())
                {
                        auto sgjet = iSg.second;
                        matched = (isOverlapping(bgjet, sgjet)) ? true : false;
                }
                if (bgjet -> GetNumberOfTracks()>0 && bgjet->Pt()>0.1 ){
                        rhovec[njetacc] = bgjet->Pt() / bgjet->Area();
                        TotaljetAreaPhys += bgjet->Area();
                        njetacc++;
                }
                n++;
        }

        Double_t OccCorr=1;
        OccCorr = TotaljetAreaPhys/TotalTPCArea;
        if (njetacc > 0 ) RHO = TMath::Median(njetacc, rhovec);
        else{ RHO = 0; }

        RHO *= OccCorr;
}
Bool_t AliAnalysisTaskGenMCRidge::isOverlapping(AliEmcalJet *jet1, AliEmcalJet *jet2)
{
        for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
        {
                Int_t jet1Track = jet1->TrackAt(i);
                for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
                {
                        Int_t jet2Track = jet2->TrackAt(j);
                        if (jet1Track == jet2Track)
                                return kTRUE;
                }
        }
        return kFALSE;
}




//___________________________________________________________________
THnSparse * AliAnalysisTaskGenMCRidge::CreateTHnSparse(TString name
        , TString title, Int_t ndim, std::vector<TAxis> bins, Option_t * opt){
  const TAxis * axises[bins.size()];
  for( UInt_t i=0;i<bins.size();i++ ) axises[i]= &bins[i];
  THnSparse * h= fHistos->CreateTHnSparse(name, title, ndim, axises,opt );
  return h;
}
//___________________________________________________________________
THnSparse * AliAnalysisTaskGenMCRidge::CreateTHnSparse(TString name
        , TString title, TString templ, Option_t * opt){
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
//___________________________________________________________________
Long64_t AliAnalysisTaskGenMCRidge::FillTHnSparse( TString name, std::vector<Double_t> x, Double_t w ){
  auto hsparse = dynamic_cast<THnSparse*>( fHistos->FindObject(name) );
  if(! hsparse ){
    cout<<"ERROR : no "<<name<<endl;
    exit(1);
  }
  return FillTHnSparse( hsparse, x, w );
}
//___________________________________________________________________
Long64_t AliAnalysisTaskGenMCRidge::FillTHnSparse( THnSparse *h, std::vector<Double_t> x, Double_t w ){
  if( int(x.size()) != h->GetNdimensions() ){
    cout<<"ERROR : wrong sized of array while Fill "<<h->GetName()<<endl;
    exit(1);
  }
  return h->Fill( &x.front(), w );
}
//___________________________________________________________________
TAxis AliAnalysisTaskGenMCRidge::AxisFix
        ( TString name, int nbin, Double_t xmin, Double_t xmax ){
                TAxis axis(nbin, xmin, xmax);axis.SetName(name);
                return axis;
}
//___________________________________________________________________
TAxis AliAnalysisTaskGenMCRidge::AxisStr( TString name, std::vector<TString> bin ){
  TAxis ax = AxisFix( name, bin.size(), 0.5, bin.size()+0.5);
  UInt_t i=1;
  for( auto blabel : bin )
    ax.SetBinLabel( i++, blabel );
  return ax;
}

TAxis AliAnalysisTaskGenMCRidge::AxisVar( TString name, std::vector<Double_t> bin ){
  TAxis axis( bin.size()-1, &bin.front() ) ;axis.SetName(name);
  return axis;
}
//___________________________________________________________________
TAxis AliAnalysisTaskGenMCRidge::AxisLog
        ( TString name, int nbin, Double_t xmin, Double_t xmax, Double_t xmin0){
  int binoffset = ( xmin0<0 || (xmin-xmin0)<1e-9) ? 0 : 1;
  std::vector<Double_t> bin(nbin+1+binoffset,0);
  double logBW3 = (log(xmax)-log(xmin))/nbin;
  for(int ij=0;ij<=nbin;ij++) bin[ij+binoffset]=xmin*exp(ij*logBW3);
  TAxis axis( nbin, &bin.front() ) ;
  axis.SetName(name);
  return axis;
}

//___________________________________________________________________
// pT_Max leading particle originally but just pt
//___________________________________________________________________
void AliAnalysisTaskGenMCRidge::ESETagging(int iESE, double pT_max) {
	double Ljetpt = -999;
	double subLjetpt = -999;
	double minSubLeadingJetPt = 5.0;
	double asym = -999;
	double InvM = -999;
#if !defined(__CINT__) && !defined(__MAKECINT__)
	AliJCDijetAna *fJFJAna = (AliJCDijetAna*)fJFJTask->GetJCDijetAna();
	vector<vector<fastjet::PseudoJet>> jets = fJFJAna->GetJets();
	int iBGSubtr = 1; // private AliJCDijetAna::iBGSubtr
	jets.at(iBGSubtr) = fastjet::sorted_by_pt(jets.at(iBGSubtr));
	vector<fastjet::PseudoJet> psjets = jets.at(iBGSubtr); // only selected jet
	fastjet::PseudoJet psLjet = jets.at(iBGSubtr).at(0);
	fastjet::PseudoJet pssubLjet = jets.at(iBGSubtr).at(1);
	Ljetpt = psLjet.pt();
	subLjetpt = pssubLjet.pt();
	asym = (Ljetpt - subLjetpt)/(Ljetpt + subLjetpt);
	InvM = ( psLjet + pssubLjet).m();
	//if(fDebugMode) cout << Form("LPpt=%.1f:LPjet=%.1f:subJet=%.1f:DiJetAsym=%.1f",pT_max,Ljetpt,subLjetpt,asym) << endl;
	switch(iESE) {
		  case 0: // Leading particle
				TagThisEvent[iESE] = kTRUE;
				break;
		  case 1: // Leading particle
				if( pT_max > fPtHardMin && pT_max < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 2: // jet
				for (unsigned ijet = 0; ijet<psjets.size(); ijet++){
					double ptt = psjets.at(ijet).pt();
					if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				}
				break;
		  case 3:  // Leading jet	
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 4: // di-jet
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax && subLjetpt > minSubLeadingJetPt ) TagThisEvent[iESE] = kTRUE;
				break;					
	} // end of switch

	fJetPt = 0.0;
	for(int i=0;i<jets.at(iBGSubtr).size();i++){
		if( jets.at(iBGSubtr).at(i).eta() > 0.4 ) continue;
		if( jets.at(iBGSubtr).at(i).pt() > fJetPt ) fJetPt = jets.at(iBGSubtr).at(i).pt();
	}
#endif 
}
