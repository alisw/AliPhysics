#include "Analyze.h"
#include <TMath.h>
#include <iostream>

// This macro reads in a starlight output file and creates histograms of 
// the p_T and rapidity of the daughters, as well as the p_T, rapidity and 
// mass of the parent.  The macro assumes there are only two daughter 
// tracks.  Currently, the macro only accomodates electrons, muons or 
// pions as daughter particles.  The histograms for the daughter particles 
// are called fPt2, fPt2, fRap1, and fRap2.  Parent histograms are created 
// for each possible daughter species (e.g., parent p_T histograms are 
// created with the names fPtEl, fPtMu, and fPtPi), but only the ones 
// corresponding to  the actual daughter particle are filled. The 
// histograms are saved in a  file called histograms.root.
//
// To use this macro, modify the file ana.C to call your input file 
// (as downloaded, it calls slight.out) and the number of events you wish to 
// process (as downloaded, it processes 20 events).  Then open root and type 
// ".x Ana.C" .

using namespace std;
Analyze::Analyze() :
  fInfile("slight.out"),
  fNEvents(1)
{
  //Constructor
  
  //Creating histograms
  fPtEl = new TH1F("PtEl", "Transverse momentum e+/e-", 100, 0, 2.);
  fPtMu = new TH1F("PtMu", "Transverse momentum mu+/mu-", 100, 0, 2.);
  fPtPi = new TH1F("PtPi", "Transverse momentum pi+/pi-", 100, 0, 2.);
  
  fRapEl = new TH1F("RapEl", "Rapidity e+/e-", 200, -10, 10);
  fRapMu = new TH1F("RapMu", "Rapidity mu+/mu-", 200, -10, 10);
  fRapPi = new TH1F("RapPi", "Rapidity pi+/pi-", 200, -10, 10);

  fInvMassEl = new TH1F("InvMassEl", "Invariant mass", 100, 0, 5);
  fInvMassMu = new TH1F("InvMassMu", "Invariant mass", 100, 0, 5);
  fInvMassPi = new TH1F("InvMassPi", "Invariant mass", 100, 0, 5);
  
  fPt1 = new TH1F("fPt1", "Transverse momentum track 1", 100, 0, 2.);
  fPt2 = new TH1F("fPt2", "Transverse momentum track 2", 100, 0, 2.);

  fRap1 = new TH1F("fRap1", "Rapidity track 1", 200, -10, 10);
  fRap2 = new TH1F("fRap2", "Rapidity track 2", 200, -10, 10);
}

Analyze::Analyze(TString infile, Int_t nEvents) :
  fInfile(infile),
  fNEvents(nEvents)
{
  //Special constructor

  //Creating histograms
  fPtEl = new TH1F("PtEl", "Transverse momentum e+/e-", 100, 0, 2.);
  fPtMu = new TH1F("PtMu", "Transverse momentum mu+/mu-", 100, 0, 2.);
  fPtPi = new TH1F("PtPi", "Transverse momentum pi+/pi-", 100, 0, 2.);
  
  fRapEl = new TH1F("RapEl", "Rapidity e+/e-", 200, -10, 10);
  fRapMu = new TH1F("RapMu", "Rapidity mu+/mu-", 200, -10, 10);
  fRapPi = new TH1F("RapPi", "Rapidity pi+/pi-", 200, -10., 10.);

  fInvMassEl = new TH1F("InvMassEl", "Invariant mass", 100, 0, 5);
  fInvMassMu = new TH1F("InvMassMu", "Invariant mass", 100, 0, 5);
  fInvMassPi = new TH1F("InvMassPi", "Invariant mass", 100, 0., 5.);
 

  fPt1 = new TH1F("fPt1", "Transverse momentum track 1", 100, 0, 2.);
  fPt2 = new TH1F("fPt2", "Transverse momentum track 2", 100, 0, 2.);

  fRap1 = new TH1F("fRap1", "Rapidity track 1", 200, -10, 10);
  fRap2 = new TH1F("fRap2", "Rapidity track 2", 200, -10, 10);
}

Analyze::~Analyze()
{
  //Destructor
  delete fPtEl;
  delete fPtMu;
  delete fPtPi;
 
  delete fRapEl;
  delete fRapMu;
  delete fRapPi;
  
  delete fInvMassEl;
  delete fInvMassMu;
  delete fInvMassPi;
  
  delete fPt1;
  delete fPt2;

  delete fRap1;
  delete fRap2;  
 
}

Int_t Analyze::Init()
{
  
  
  cout << "Opening textfile " << fInfile << endl;
  if( !(filelist=fopen(fInfile.Data(),"r")) ){
    cout<<"Couldn't open input file: "<<fInfile<<endl;
    return -1;
  }
  cout << "Done opening textfile" << endl;
  return 0;
}

Int_t Analyze::NextEvent()
{
  char linelabel[20];
  int i1=0;
  int i2=0;
  int i3=0;
  double x1=0.0;
  double x2=0.0;
  double x3=0.0;
  double x4=0.0;
  int ntrk=0;
  int nvtx=0;
  // Event line 
  fscanf(filelist,"%s %d %d %d ",linelabel,&i1,&ntrk,&i2);
  fNParticles = ntrk;
  cout<<linelabel<<"  "<<i1<<"  "<<ntrk<<"  "<<i2<<"   "<<fNParticles<<endl;
  // Vertex line 
  fscanf(filelist,"%s %lf %lf %lf %lf %d %d %d %d",
         linelabel,&x1,&x2,&x3,&x4,&i1,&i2,&i3,&nvtx);
  cout<<linelabel<<"  "<<x1<<"  "<<x2<<"  "<<x3<<"  "<<x4<<"  "<<i1<<"  "<<i2<<"  "<<i3<<"  "<<nvtx<<endl;
  if(ntrk != nvtx)cout<<"ERROR: ntrk = "<<ntrk<<"  nvtx = "<<nvtx<<endl;
  //
  return fNParticles;
}

TParticle* Analyze::NextParticle()
{
  char tracklabel[20];
  int i1=0;
  int i2=0;
  int i3=0;
  int i4=0;
  Int_t idpart = 0;
  Double_t px = 0.0;
  Double_t py = 0.0;
  Double_t pz = 0.0;
  Double_t ep = 0.0;
  
  cout<<"In NextParticle: fNparticles = "<<fNParticles<<endl;

  // for ( int itk=0; itk < fNParticles; itk++){


    fscanf(filelist,"%s %d %le %le %le %d %d %d %d",
         tracklabel,&i1,&px,&py,&pz,&i2,&i3,&i4,&idpart);
    cout<<"   "<<tracklabel<<"  "<<i1<<"  "<<px<<"  "<<py<<"  "<<pz<<"  "<<i2<<"  "<<i3<<"  "<<i4<<"  "<<idpart<<endl;
  
    TParticle *particle = 
      new TParticle(idpart, 0, -1, -1, -1, -1, px, py, pz, ep, 0., 0., 0., 0.);
      if(idpart == 11 || idpart == -11){particle->SetCalcMass(0.00051099907);}
      else if(idpart == 13 || idpart == -13){particle->SetCalcMass(0.105658389);}
      else if(idpart == 211 || idpart == -211){particle->SetCalcMass(0.13956995);}
      else {cout << "unknown daughter!  please edit the code to accomodate"<< endl;
	      exit(0);
	   }

  return particle;
}

void Analyze::doAnalysis()
{

  Int_t check = Init();
  if(check < 0) return;
  Double_t mass;
  //Doing the analysis:  loop over events
  for(Int_t ev = 0; ev < fNEvents; ev++){
    	   
    const Int_t ntracks = NextEvent();
    //Array of TLorentzVectors. One vector for each track
    TLorentzVector* vecArr[ntracks];
      Int_t idpart = 0;
    //Looping over the tracks of the event
    for(Int_t tr = 0; tr < ntracks; tr++){
      //Getting the TParticle for the track
      TParticle *part = NextParticle();
      mass = part->GetCalcMass();
      idpart = part->GetPdgCode();
      Double_t energy = TMath::Sqrt(mass*mass+part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
      //Creating a new TLorentzVector and setting px, py, pz and E.
      vecArr[tr] = new TLorentzVector;
      vecArr[tr]->SetPxPyPzE(part->Px(), part->Py(), part->Pz(), energy); 
      cout << "particle " << tr << ": px: " << part->Px() << " py: " << part->Py() << " pz: " << part->Pz() << " Energy: " << energy << endl;
    }
  
    // Fill the individual track histograms

    fPt1->Fill(vecArr[0]->Pt());
    fPt2->Fill(vecArr[1]->Pt());

    fRap1->Fill(vecArr[0]->Rapidity());
    fRap2->Fill(vecArr[1]->Rapidity());
    
    //Creating a new TLorentzVector, which is the sum of the elements in vecArr
    TLorentzVector sum;
    for(Int_t i = 0; i < ntracks; i++){
      sum += *vecArr[i];
    }
    //Filling the parent histograms, depending on daughter particle type
    
    if(idpart== 11 || idpart== -11){
      fPtEl->Fill(sum.Pt());
      cout << "sum.Rapidity: " << sum.Rapidity() << endl;
      cout << "sum.M(): " << sum.M() << endl;
      fRapEl->Fill(sum.Rapidity());
      fInvMassEl->Fill(sum.M());
    }
    else if(idpart == 13 || idpart == -13){
      fPtMu->Fill(sum.Pt());
      cout << "sum.Rapidity: " << sum.Rapidity() << endl;
      cout << "sum.M(): " << sum.M() << endl;
      fRapMu->Fill(sum.Rapidity());
      fInvMassMu->Fill(sum.M());
    }
    else if(idpart == 211 || idpart == -211){
      fPtPi->Fill(sum.Pt());
      cout << "sum.Rapidity: " << sum.Rapidity() << endl;
      cout << "sum.M(): " << sum.M() << endl;
      fRapPi->Fill(sum.Rapidity());
      fInvMassPi->Fill(sum.M());
    }
  
  }

  //Writing the histograms to file
  TFile file("histograms.root", "RECREATE");
  fPtEl->Write();
  fRapEl->Write();
  fInvMassEl->Write();
  fPtMu->Write();
  fRapMu->Write();
  fInvMassMu->Write();
  fPtPi->Write();
  fRapPi->Write();
  fInvMassPi->Write();
  fPt1->Write();
  fPt2->Write();
  fRap1->Write();
  fRap2->Write();
  
}
