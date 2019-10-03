//---------------------------------------------------------------------------------
// The example of usage of AliKFParticle & AliKFVertex classes for V0 analysis
// .
// @author  S.Gorbunov, I.Kisel
// @version 1.0
// @since   13.05.07
// 
// The AliKFParticleTest macro contains a toy V0 finder for ESD tracks.
// At the first step, event primary vertex is reconstructed. 
// At the second step, ideal PID hypothesis are assigned to all the particles.
// At the third step, V0 candidates are constructed for each pair 
// of positive and negative particles.
// V0 candidate considered as good when it passes Chi^2 cut and 
// it is placed >= 3 Sigma away from the primary vertex.
// Invariant mass distribution for all good V0 candidates is plotted.
//
//  -= Copyright &copy ALICE HLT Group =-
//_________________________________________________________________________________


#ifndef ALIKFPARTICLETEST_H
#define ALIKFPARTICLETEST_H

#include "TH1.h"
#include "TCanvas.h"
#include "AliESD.h"
#include "AliAnalysisTaskRL.h"

class AliKFParticleTest: public AliAnalysisTaskRL {
 public:
  AliKFParticleTest(const char *name);
  virtual ~AliKFParticleTest() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESD *fESD; //ESD object
  TH1D   *fHistoMassAll;
  TH1D   *fHistoMassSignal;
  TCanvas *fCanvas;
  void DrawV0(); 

  TObjArray * fOutputContainer; // ! output data container

  ClassDef(AliKFParticleTest, 0); // example of analysis
};


#include "Riostream.h"
#include "TChain.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TParticle.h"
#include "TRandom.h"
#include "AliESD.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliTracker.h"
#include "TLatex.h"
#include "TDatabasePDG.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"


ClassImp(AliKFParticleTest);


//________________________________________________________________________
AliKFParticleTest::AliKFParticleTest(const char *name) :AliAnalysisTaskRL(name,""), fESD(0) {
  // Constructor.
  // Input slot #0 works with an Ntuple
  DefineInput(0, TChain::Class());
 // Output slot #0 writes into a TH1 container
  DefineOutput(0, TObjArray::Class());
}

//___________________________________________________________________________
void AliKFParticleTest::ConnectInputData(Option_t *) {
  // Initialize branches.
  printf("   ConnectInputData of task %s\n", GetName());
  if (!fESD) {
    char ** address = (char **)GetBranchAddress(0, "ESD");
    if (address) fESD = (AliESD*)(*address);
    if (!fESD) {
      fESD = new AliESD();
      SetBranchAddress(0, "ESD", &fESD);
    }
  }

}

//___________________________________________________________________________
void AliKFParticleTest::CreateOutputObjects() {
  printf("   CreateOutputObjects of task %s\n", GetName());

  if ( !gROOT->IsBatch() ) fCanvas = new TCanvas();
  else fCanvas = 0;

  fHistoMassAll = new TH1D("massAll","AliKFParticle test", 500,0,3);
  fHistoMassAll->SetXTitle("V^{0} invariant mass [GeV]");
  fHistoMassAll->SetLineColor(kGreen);
  fHistoMassAll->SetFillColor(kGreen);

  fHistoMassSignal = new TH1D("massSignal","V^{0} signal", 500,0,3);
  fHistoMassSignal->SetXTitle("V^{0} invariant mass [GeV]");
  fHistoMassSignal->SetLineColor(kBlue);
  fHistoMassSignal->SetFillColor(kBlue);

  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->Add(fHistoMassAll);
  fOutputContainer->Add(fHistoMassSignal);
}


void AliKFParticleTest::DrawV0() 
{
  //* Draw the invariant mass histogram

  if ( gROOT->IsBatch() ) return;

  const Int_t histoPDG [4]= {22,310,3122,421};
  const Char_t* histoName[4]= {"#gamma","K^{0}_{s}","#Lambda","D^{0}"};
  TLatex histoText[4];

  if( !fCanvas ) return;
  fCanvas->Clear();
  fCanvas->cd();

  for( Int_t i=0; i<4; i++ ){
    histoText[i].SetTextColor(kBlue);
  }

  fHistoMassAll->Draw();
  fHistoMassSignal->Draw("same");

  for( Int_t i=0; i<4; i++ ){
    Double_t mass = TDatabasePDG::Instance()->GetParticle(histoPDG[i])->Mass();    
    Int_t bin = fHistoMassSignal->FindBin(mass) -5;
    if( bin<0 ) bin =0;
    Double_t max = 0;
    for( Int_t j=bin; j<bin+10; j++ ) 
    if( max<fHistoMassSignal->GetBinContent(j) ) max = fHistoMassSignal->GetBinContent(j);
    if(max>0) histoText[i].DrawLatex( mass+.05, max, histoName[i] );
  }
  fCanvas->Update();
}

//________________________________________________________________________
void AliKFParticleTest::Exec(Option_t *) {

  static Int_t iEvent = 0;
  if( ++iEvent%100 ==0 ) cout<<"Event "<<iEvent<<endl;

  // Get input data
  TTree *tinput = (TTree*)GetInputData(0);
  Long64_t ientry = tinput->GetReadEntry();
  if (AliAnalysisTaskRL::GetEntry(ientry) == kFALSE) {
    printf("Couldn't get event from the runLoader\n");
    return;
  }  
  if (!fESD) return;

  Double_t Bz = fESD->GetMagneticField();
  AliKFParticle::SetField( Bz );

  if (ientry==0) Notify ();

 //cout <<"Event "<<ientry<<endl;

  AliStack* stack = GetStack();
  if (!stack) {
    AliDebug(AliLog::kError, "Stack not available");
    return;
  }


  class TESDTrackInfo
  {
  public:
    TESDTrackInfo(){}
    AliKFParticle fParticle; //* assigned KFParticle
    Bool_t fPrimUsedFlag;    //* flag shows that the particle was used for primary vertex fit
    Bool_t fOK;              //* is the track good enough
    Int_t mcPDG;             //* Monte Carlo PDG code
    Int_t mcMotherID;        //* Monte Carlo ID of its mother
  };

  Int_t nESDTracks=fESD->GetNumberOfTracks();
  if( nESDTracks>1000 ) nESDTracks=1000;

  TESDTrackInfo ESDTrackInfo[1000]; //* parallel to ESD tracks

  //* Fill ESDTrackInfo array

  for (Int_t iTr=0; iTr<nESDTracks; iTr++)
    {   
      TESDTrackInfo &info = ESDTrackInfo[iTr];
      info.fOK = 0;
      info.fPrimUsedFlag = 0;
      info.mcPDG = -1;
      info.mcMotherID = -1;

      //* track quality check

      AliESDtrack *pTrack = fESD->GetTrack(iTr);    
      if( !pTrack  ) continue;
      if (pTrack->GetKinkIndex(0)>0) continue;
      if ( !( pTrack->GetStatus()&AliESDtrack::kITSrefit ) ) continue;
      Int_t indi[12];
      if( pTrack->GetITSclusters(indi) <5 ) continue;
      Int_t PDG = ( pTrack->Get1Pt() <0 ) ?321 :211;

      //* take MC PDG  
      { 
	Int_t mcID = TMath::Abs(pTrack->GetLabel());
	TParticle * part = stack->Particle(TMath::Abs(mcID));
	if( part ){
	  info.mcPDG = part->GetPdgCode();
	  PDG = info.mcPDG;
	  if( mcID>=0 ) info.mcMotherID = part->GetFirstMother();
	}    
      }

      //* Construct KFParticle for the track

      info.fParticle = AliKFParticle( *pTrack, PDG );
      info.fOK = 1;   
    }

  //* Find event primary vertex
  
  AliKFVertex primVtx;  
  {
    const AliKFParticle * vSelected[1000]; //* Selected particles for vertex fit
    Int_t vIndex[1000];                    //* Indices of selected particles
    Bool_t vFlag[1000];                    //* Flags returned by the vertex finder

    Int_t nSelected = 0;
    for( Int_t i = 0; i<nESDTracks; i++){ 
      if(ESDTrackInfo[i].fOK ){
	vSelected[nSelected] = &(ESDTrackInfo[i].fParticle);
	vIndex[nSelected] = i;
	nSelected++;
      }
    }
    primVtx.ConstructPrimaryVertex( vSelected, nSelected, vFlag, 3. );
    for( Int_t i = 0; i<nSelected; i++){ 
      if( vFlag[i] ) ESDTrackInfo[vIndex[i]].fPrimUsedFlag = 1;
    }
    if( primVtx.GetNDF() <1 ) return; //* Less then two tracks in primary vertex 
  }

  //* V0 finder

  for( Int_t iTr = 0; iTr<nESDTracks; iTr++ ){ //* first daughter

    TESDTrackInfo &info = ESDTrackInfo[iTr];
    if( !info.fOK ) continue;    

    for( Int_t jTr = iTr+1; jTr<nESDTracks; jTr++ ){  //* second daughter
      TESDTrackInfo &jnfo = ESDTrackInfo[jTr];
      if( !jnfo.fOK ) continue;
      
      //* check for different charge

      if( info.fParticle.GetQ() == jnfo.fParticle.GetQ() ) continue;      

      //* construct V0 mother

      AliKFParticle V0( info.fParticle, jnfo.fParticle );     

      //* check V0 Chi^2

      if( V0.GetNDF()<1 ) continue;
      if( TMath::Sqrt(TMath::Abs(V0.GetChi2()/V0.GetNDF())) >3. ) continue;

      //* subtruct daughters from primary vertex 

      AliKFVertex primVtxCopy = primVtx;
       
      if( info.fPrimUsedFlag ) primVtxCopy -= info.fParticle;
      if( jnfo.fPrimUsedFlag ) primVtxCopy -= jnfo.fParticle;

      //* Check V0 Chi^2 deviation from primary vertex 

      if( V0.GetDeviationFromVertex( primVtxCopy ) >3. ) continue;

      //* Add V0 to primary vertex to improve the primary vertex resolution

      primVtxCopy += V0;      

      //* Set production vertex for V0

      V0.SetProductionVertex( primVtxCopy );

      //* Check chi^2 for a case

      if( TMath::Sqrt( TMath::Abs(V0.GetChi2()/V0.GetNDF()) >3. )) continue;

      //* Get V0 decay length with estimated error

      Double_t length, sigmaLength;
      if( V0.GetDecayLength( length, sigmaLength ) ) continue;

      //* Reject V0 if it decays too close to the primary vertex

      if( length  <3.*sigmaLength ) continue;

      //* Get V0 invariant mass and plot it

      Double_t mass, sigmaMass;
      if( V0.GetMass( mass, sigmaMass ) ) continue;
      fHistoMassAll->Fill(mass);            

      //* Fill signal histograms using MC information

      if( info.mcMotherID==jnfo.mcMotherID && info.mcMotherID>=0 ){
	TParticle *mother = stack->Particle(info.mcMotherID);
	if( mother && TMath::Abs(mother->GetPdgCode())!=21 ){	 
	  if( mother->GetNDaughters()==2 ){
	    fHistoMassSignal->Fill(mass);	
	  }
	  cout<<"PDG V0,pi,pj, ndaughters, mc mass, reco mass = "<<mother->GetPdgCode()<<","<<info.mcPDG<<","<<jnfo.mcPDG<<", "
	      << mother->GetNDaughters()<<", "<<mother->GetMass()<<", "<<mass<<endl;
	}
      }
    }
  }
  if( iEvent %1000 == 0 || (iEvent ==200)) DrawV0();     

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(0, fOutputContainer);
}      


//________________________________________________________________________
void AliKFParticleTest::Terminate(Option_t *) 
{
  DrawV0();
  fOutputContainer = (TObjArray*)GetOutputData(0);
  //fHistoMassAll=(TH1D*)fOutputContainer->At(0) ;
}

#endif
