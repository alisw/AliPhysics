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

using namespace std;

const Double_t pi = TMath::Pi();



//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge()
:AliAnalysisTaskEmcalJet("AliAnalysisTaskGenMCRidge")

{
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge
(     
      const char *name
    , const char *option
)
:AliAnalysisTaskEmcalJet(name)
{
    DefineOutput (1, AliDirList::Class());
}

//___________________________________________________________________
AliAnalysisTaskGenMCRidge::AliAnalysisTaskGenMCRidge
(
      const AliAnalysisTaskGenMCRidge& ap
)
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
    delete fTrigger;
    delete fTrackCuts;
    delete fRunTable;
    delete fOutput;
}

void AliAnalysisTaskGenMCRidge::UserCreateOutputObjects()
{

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

	binPhi = AxisFix("phi",32,-0.5*pi-pi/32.0, 1.5*pi-pi/32.0);
	binEta = AxisFix("eta",80,-4.0,4.0);
	if( fOption.Contains("OnlyOneDim") ){
		binEta = AxisVar("eta",LongRangeEta);
	}

	binCent = AxisVar("Cent",varmultbin);
	binTPt = AxisVar("TPt",ptbin);
        binAPt = AxisVar("TPt",ptbin);
	binLHPt = AxisVar("LHPT",ltpttrackbin);
	binJetPt = AxisVar("JetPT",jetptbin);


	CreateTHnSparse("hNTrigLH","hNTrigLH",6,{binCent,binPhi,binEta,binTPt,binAPt,binLHPt},"s");
	CreateTHnSparse("hRidgeLH","hRidgeLH",6,{binCent,binPhi,binEta,binTPt,binAPt,binLHPt},"s");
	CreateTHnSparse("hRidgeMixingLH","hRidgeMixingLH",6,{binCent,binPhi,binEta,binTPt,binAPt,binLHPt},"s");

	CreateTHnSparse("hNTrigJet","hNTrigJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetPt},"s");
	CreateTHnSparse("hRidgeJet","hRidgeJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetPt},"s");
	CreateTHnSparse("hRidgeMixingJet","hRidgeMixingJet",6,{binCent,binPhi,binEta,binTPt,binAPt,binJetPt},"s");



	fHistos->CreateTH1("hNChargedMultiplixity","hNChargedMultiplixity",1000,0,1000);

	fOutput -> Add( fHistos->GetListOfHistograms() );
	PostData(1, fOutput);

}
//___________________________________________________________________
void AliAnalysisTaskGenMCRidge::Exec(Option_t *)
{
	fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

	if( fMcHandler ) fMcEvent = fMcHandler->MCEvent();
	else{ return; }

	if( !fMcEvent ) return;
	fStack = ((AliMCEvent*)fMcEvent)->Stack();
	if( !fStack ) return;

	PostData(1, fListOfObjects);
	return;

}


void AliAnalysisTaskGenMCRidge::FinishTaskOutput()
{
    //fOutput->Write();
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
