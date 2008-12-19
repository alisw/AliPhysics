
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



#if !defined( __CINT__) || defined(__MAKECINT__)

  //gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");

  #include <Riostream.h>
  #include <TTree.h>
  #include <TFile.h>
  #include <TCanvas.h>
  #include <TH1D.h>
  #include <TLatex.h>
  #include <TDatabasePDG.h>
  #include <TParticle.h>

  #include "AliRun.h"
  #include "AliStack.h"
  #include "AliESDEvent.h"
  #include "AliTracker.h"
  #include "AliKFParticle.h"
  #include "AliKFVertex.h"

#endif


Bool_t old_file = 1;
TCanvas *canvas=0;
TH1D  *histoMass[3]={0,0,0};
Int_t histoPDG [4]= {22,310,3122,421};
char* histoName[4]= {"#gamma","K^{0}_{s}","#Lambda","D^{0}"};
TLatex *histoText[4]={0,0,0,0};

void DrawV0() 
{
  //* Draw the invariant mass histogram

  if( !canvas ) return;
  histoMass[0]->Draw();
  for( Int_t i=1; i<3; i++ ) histoMass[i]->Draw("same");
  for( Int_t i=0; i<4; i++ ){
    Double_t mass = TDatabasePDG::Instance()->GetParticle(histoPDG[i])->Mass();    
    Int_t bin = histoMass[2]->FindBin(mass) -5;
    if( bin<0 ) bin =0;
    Double_t max = 0;
    for( Int_t j=bin; j<bin+10; j++ ) 
    if( max<histoMass[2]->GetBinContent(j) ) max = histoMass[2]->GetBinContent(j);
    histoText[i]->SetY(max);
    histoText[i]->SetX(mass+.05);
    if(max>0) histoText[i]->Draw();
  }
  if( canvas ) canvas->Update();
}

void StartV0()
{
  //* (Re)create histograms and fill them from temporary file
  
  TDirectory *curr = gDirectory;
  canvas = new TCanvas();
  histoMass[0] = new TH1D("massAll","AliKFParticle test", 500,0,3);
  histoMass[0]->SetXTitle("V^{0} invariant mass [GeV]");
  histoMass[1] = new TH1D("massMulti","V^{+-} contributions", 500,0,3);
  histoMass[2] = new TH1D("massV0","V^{0} signal", 500,0,3);
  histoMass[0]->SetLineColor(8);
  histoMass[0]->SetFillColor(8);
  histoMass[1]->SetLineColor(kMagenta);
  histoMass[1]->SetFillColor(kMagenta);
  histoMass[2]->SetLineColor(kBlue);
  histoMass[2]->SetFillColor(kBlue);
  for( Int_t i=0; i<4; i++ ){
    histoText[i] = new TLatex(0,0,histoName[i]);
    histoText[i]->SetTextColor(kBlue);
  }
  TFile* file = new TFile("AliKFParticleTest.root", "READ" );
  if( !old_file && file ){
    file->cd();
    for( Int_t i=0; i<3; i++ ){
      TObject *old = file->FindObjectAny(histoMass[i]->GetName());
      if( old && histoMass[i] ) histoMass[i]->Add( (TH1D*) old);
    }
    file->Close();
  }
  curr->cd();
  DrawV0();
}

void EndV0()
{
  //* End of gAlice -> store all histos in temporary file, clean memory

  TDirectory *curr = gDirectory;
  TFile* file = new TFile("AliKFParticleTest.root", "RECREATE" );
  old_file = 0;
  file->cd();
  for( Int_t i=0; i<3; i++ ) histoMass[i]->Write();
  file->Close();
  for( Int_t i=0; i<3; i++ ) delete histoMass[i];
  for( Int_t i=0; i<4; i++ ) delete histoText[i];
  delete canvas;
  curr->cd();
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

void RunV0(  AliESDEvent *event )
{
  //* V0 finder

  static Int_t iEvent = 0;
  cout<<"Event "<<iEvent++<<endl;

  if( !gAlice ) return;
  AliRunLoader *rl = AliRunLoader::GetRunLoader(); 
 
  AliStack *stack = rl->Stack();
  if( !stack ) return;  

  Int_t nESDTracks=event->GetNumberOfTracks();
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

      AliESDtrack *pTrack = event->GetTrack(iTr);    
      if( !pTrack  ) continue;
      if (pTrack->GetKinkIndex(0)>0) continue;
      if ( !( pTrack->GetStatus()&AliESDtrack::kITSrefit ) ) continue;
      Int_t indi[12];
      if( pTrack->GetITSclusters(indi) <5 ) continue;
      Int_t PDG = ( pTrack->GetSigned1Pt() <0 ) ?321 :211;

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
      histoMass[0]->Fill(mass);            

      //* Fill signal histograms using MC information

      if( info.mcMotherID==jnfo.mcMotherID && info.mcMotherID>=0 ){
	TParticle *mother = stack->Particle(info.mcMotherID);
	if( mother && TMath::Abs(mother->GetPdgCode())!=21 ){
	  histoMass[1]->Fill(mass);
	  if( mother->GetNDaughters()==2 ){
	    histoMass[2]->Fill(mass);	
	  }
	  cout<<"PDG V0,pi,pj, ndaughters, mc mass, reco mass = "<<mother->GetPdgCode()<<","<<info.mcPDG<<","<<jnfo.mcPDG<<", "
	      << mother->GetNDaughters()<<", "<<mother->GetMass()<<", "<<mass<<endl;
	}
      }
    }
  }
  if( iEvent %100 == 0 || (iEvent<100 && iEvent %10==0)) DrawV0();     
}



Int_t AliKFParticleTest(Int_t n1=0,Int_t n2=1000,char *dire="/d/alice10/sma/my_v4-05-Release/pp/"){
  //* Main macro
    
  //  LOOP  OVER  SERIES  OF  DIRECTORIES
  
  for (Int_t ifi=n1; ifi<=n2; ifi++) {
   
    char nstring[5], filename[100], esdfile[100];
    sprintf(nstring,"%3.3d",ifi);
    sprintf(filename,"%s%s/galice.root",dire,nstring);
    sprintf(esdfile,"%s%s/AliESDs.root",dire,nstring);
    
    cout <<"  Opening "<<filename<<"\nand ESD "<<esdfile<<endl;
    
    if (gAlice) {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
    
    AliRunLoader *rl = AliRunLoader::Open(filename);
    if ( !rl ) { 
      ::Error("AliKFParticleTest.C","Can not open session !");
      continue;
    }
    if (rl->LoadgAlice()) {
      ::Error("AliKFParticleTest.C","LoadgAlice returned error !");
      delete rl;
      continue;
    }
    if (rl->LoadHeader()) {
      ::Error("AliKFParticleTest.C","LoadHeader returned error !");
      delete rl;
      continue;
    }
    rl->LoadKinematics();
    AliTracker::SetFieldMap(gAlice->Field(),1);
    AliKFParticle::SetField( AliTracker::GetBz() );

    //---------------------------------------//
    //                                       //
    //               ESD  file               //
    //                                       //
    //---------------------------------------//
    
    // Open file with the ESD
    TFile *ef=TFile::Open(esdfile);

    //Check if the file could be opened
    if (!ef || !ef->IsOpen()) {cerr<<"Error open AliESDs.root !\n"; continue ;}

    //create event object
    AliESD *event = new AliESDEvent;

    //Set pointer to the esd tree in the file
    TTree* tree = (TTree*) ef->Get("esdTree");
    
    //check if the tree exists
    if (!tree) {cerr<<"no ESD tree found\n"; continue;};
    
    //Set pointer to the esd object in the tree
    event->ReadFromTree(event);
    
    //Number of events
    Int_t nevents=tree->GetEntriesFast();
    cout << "Number of events: " << nevents << endl;

    StartV0();
    for (Int_t iev=0; iev<nevents; iev++){
      tree->GetEvent(iev);
      rl->GetEvent(iev);
      RunV0(event);
    }
    EndV0();
    delete event;
    ef->Close();
  }
  StartV0();
  return 0;
}
