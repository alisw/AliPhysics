/*******************************************************************************
 *                                                                             *
 *    HYDJET++ , event generator under the ROOT FRAMEWORK for simulation of    *
 *    relativistic heavy ion AA collisions as the superposition of soft,       *
 *    hydro-type state and hard, multi-parton state.                           *
 *                                                                             *
 *     The main routine is written in the object-oriented C++ language         *        
 *     under the ROOT environment. The hard, multi-partonic part of            *  
 *     HYDJET++ event is identical to the hard part of Fortran-written         *
 *     HYDJET (PYTHIA6.4xx + PYQUEN1.5) and is included in the generator       *
 *     structure as the separate directory. The soft part of HYDJET++          * 
 *     event represents the "thermal" hadronic state obtained with the         *
 *     parameterization Bjorken-like of freeze-out hypersurface and            *
 *     includes longitudinal, radial and elliptic flow effects and             *
 *     decays of hadronic resonances. The corresponding fast                   * 
 *     Monte-Carlo simulation procedure (C++ code) FAST MC is adapted.         *
 * --------------------------------------------------------------              *
 *     Web-page:                                                               *
 *    http://cern.ch/lokhtin/hydjet++                                          *   
 *     --------------------------------------------------------------          *  
 *                                                                             *                                                                             *
 *                                                                             *
 * This program is a free software; you can use and redistribute it freely.    *  
 * Any publication of results obtained using this code must reference          * 
 *                                                                             *
 *                                                                             * 
 *                                                                             *
 *      Main reference for HYDJET++:                                           *
 *     I.P. Lokhtin, L.V. Malinina, S.V. Petrushanko, A.M. Snigirev,           *
 *     I. Arsene, K. Tywoniuk, submitted to Comp. Phys. Comm.                  *
 *                                                                             * 
 *     Reference for HYDJET and PYQUEN:                                        *
 *     I.P. Lokhtin, A.M. Snigirev, Eur. Phys. J. C 46 (2006) 211;             *
 *     http://cern.ch/lokhtin/hydro/hydjet.html                                * 
 *     http://cern.ch/lokhtin/pyquen.                                          *  
 *                                                                             *    
 *     Reference for PYTHIA6.4:                                                *
 *     T.Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026;                *
 *     http://home.thep.lu.se/~torbjorn/Pythia.html.                           * 
 *                                                                             * 
 *     References for FAST MC:                                                 *  
 *     N.S. Amelin, R. Lednicky, T.A. Pocheptsov, I.P. Lokhtin,                * 
 *     L.V. Malinina, A.M. Snigirev, Iu.A. Karpenko and Yu.M. Sinyukov,        * 
 *     Phys. Rev. C 74 (2006) 064901;                                          *
 *     N.S. Amelin, I. Arsene, L. Bravina, Iu.A. Karpenko, R. Lednicky,        *  
 *     I.P. Lokhtin, L.V. Malinina, A.M. Snigirev and Yu.M. Sinyukov,          *  
 *     Phys. Rev. C 77 (2008) 014903;                                          *
 *     http://uhkm.jinr.ru.                                                    *   
 *                                                                             *
 *     Reference for nuclear shadowing model:                                  *
 *     K. Tywoniuk, I.C. Arsene, L. Bravina, A. Kaidalov and                   *
 *     E. Zabrodin, Phys. Lett. B 657 (2007) 170.                              *
 *                                                                             * 
 *       version 2.0:                                                          *
 *                                                                             *
 *     Igor Lokhtin, SINP MSU, Moscow, RU                                      *
 *     e-mail: Igor.Lokhtin@cern.ch                                            *
 *                                                                             *
 *     Ludmila Malinina, SINP MSU, Moscow, RU                                  *   
 *     e-mail: malinina@lav01.sinp.msu.ru                                      * 
 *                                                                             *
 *******************************************************************************/ 

#include <iostream> 
#include <fstream>
#include <vector>
#include <time.h>

#include <TNtuple.h>
#include <TError.h>
#include <TTree.h>
#include <TFile.h>

#include "InitialState.h"
#include "InitialStateHydjet.h"

#include <TRandom.h>
#include "Particle.h"


#include "HYJET_COMMONS.h"
extern HYIPARCommon HYIPAR;
extern HYFPARCommon HYFPAR;
extern HYJPARCommon HYJPAR;
extern HYPARTCommon HYPART;
extern SERVICECommon SERVICE;


//Main program:
//reads input parameters from file "RunInputHydjet" ;
//calculates particle densities and average initial multiplicities and writes them
//in output file "multiplicities.txt";
//creates trees (tree with direct hadrons and hadrons after resonance decays)
//with space-time and momentum-energy information of produced hadrons;
//writes trees in file "RunOutput.root".


int main()

{

  clock_t start;
  start = clock();
 
//new
  time_t  now;
  struct tm  *ts;
  char       buf[80];
         
 // Get the current time
   time(&now);
              
 // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
    ts = localtime(&now);
    strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts);
    printf("%s\n", buf);
                                    
 
 
 
  TFile *outputFile=new TFile("RunOutput.root", "RECREATE"); 

  //SET MAXIMAl VALUE OF PARTICLE MULTIPLICITY!!!
  const Int_t kMax = 500000; 
 
  Int_t npart;
  //define event number
  Int_t nev;

//define event characteristics: 
//total event multiplicity, number of produced hadrons in hard part/soft part
  Int_t Ntot, Npyt, Nhyd;

// number of jets, number of binary collisions, number of participants :
  Int_t Njet, Nbcol, Npart;
  
 //impact parameter 
  Float_t Bgen, Sigin, Sigjet;
 
  //define hadron characteristic vectors
  std::vector<Int_t> pdg(kMax); //pdg encodings
  std::vector<Int_t> Mpdg(kMax);//pdg encodings for mother hadrons
  std::vector<Int_t> type(kMax);//type of particle: 0-from hydro or decays, 1111 -from jets
  std::vector<Float_t> Px(kMax);//x-hadron momentum component,[GeV/c]
  std::vector<Float_t> Py(kMax);//y-hadron momentum component,[GeV/c]
  std::vector<Float_t> Pz(kMax);//z-hadron momentum component,[GeV/c]
  std::vector<Float_t> E(kMax); //hadron total energy,[GeV]  
  std::vector<Float_t> X(kMax);//x-hadron coordinate component,[fm]
  std::vector<Float_t> Y(kMax);//y-hadron coordinate component,[fm]
  std::vector<Float_t> Z(kMax);//z-hadron coordinate component,[fm]
  std::vector<Float_t> T(kMax);//hadron time,[fm/c] 
   
   
   
  TTree *ti=new TTree("ti","Initial");

  ti->Branch("nev",&nev,"nev/I");
  ti->Branch("Bgen",&Bgen,"Bgen/F");
  ti->Branch("Sigin",&Sigin,"Sigin/F");
  ti->Branch("Sigjet",&Sigjet,"Sigjet/F");
  ti->Branch("Ntot",&Ntot,"Ntot/I");
  ti->Branch("Nhyd",&Nhyd,"Nhyd/I");
  ti->Branch("Npyt",&Npyt,"Npyt/I");
  ti->Branch("Njet",&Njet,"Njet/I");  
  ti->Branch("Nbcol",&Nbcol,"Nbcol/I");
  ti->Branch("Npart",&Npart,"Npart/I");
  ti->Branch("Px",&Px[0],"Px[Ntot]/F");
  ti->Branch("Py",&Py[0],"Py[Ntot]/F");
  ti->Branch("Pz",&Pz[0],"Pz[Ntot]/F");
  ti->Branch("E",&E[0],"E[Ntot]/F");  
  ti->Branch("X",&X[0],"X[Ntot]/F");
  ti->Branch("Y",&Y[0],"Y[Ntot]/F");
  ti->Branch("Z",&Z[0],"Z[Ntot]/F");
  ti->Branch("T",&T[0],"T[Ntot]/F");
  ti->Branch("pdg",&pdg[0],"pdg[Ntot]/I");
  ti->Branch("Mpdg",&Mpdg[0],"Mpdg[Ntot]/I");

  TTree *td=new TTree("td","After decays");
  td->Branch("nev",&nev,"nev/I");
  td->Branch("Bgen",&Bgen,"Bgen/F");
  td->Branch("Sigin",&Sigin,"Sigin/F");
  td->Branch("Sigjet",&Sigjet,"Sigjet/F");
  td->Branch("Ntot",&Ntot,"Ntot/I");
  td->Branch("Nhyd",&Nhyd,"Nhyd/I");
  td->Branch("Npyt",&Npyt,"Npyt/I");
  td->Branch("Njet",&Njet,"Njet/I");  
  td->Branch("Nbcol",&Nbcol,"Nbcol/I");
  td->Branch("Npart",&Npart,"Npart/I");
  td->Branch("Px",&Px[0],"Px[Ntot]/F");
  td->Branch("Py",&Py[0],"Py[Ntot]/F");
  td->Branch("Pz",&Pz[0],"Pz[Ntot]/F");
  td->Branch("E",&E[0],"E[Ntot]/F");  
  td->Branch("X",&X[0],"X[Ntot]/F");
  td->Branch("Y",&Y[0],"Y[Ntot]/F");
  td->Branch("Z",&Z[0],"Z[Ntot]/F");
  td->Branch("T",&T[0],"T[Ntot]/F");
  td->Branch("pdg",&pdg[0],"pdg[Ntot]/I");
  td->Branch("Mpdg",&Mpdg[0],"Mpdg[Ntot]/I");
  td->Branch("type",&type[0],"type[Ntot]/I");
  
  
  InitialState *FASTMC;
    FASTMC = new InitialStateHydjet();
 
   
  if(!FASTMC->ReadParams()) {
    Error("RunHadronSource::main", "No initial model parameters found!!\n");
    return 0;
  }
 
  if(!FASTMC->MultIni()) {
    Error("RunHadronSource::main", "Initial multiplicities are zero!!\n");
    return 0;
  }

  ParticleAllocator allocator;
  List_t source;
  List_t secondaries;
  std::cout << "Generating " << FASTMC->GetNev() << " events" << std::endl;
  std::cout << "Starting the event loop" << std::endl;
    
  // Set Random Number seed 
  Int_t sseed =0;               //Set 0 to use the current time 
  gRandom->SetSeed(sseed); 
  std::cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<std::endl;  
  
  // Loop over events  
  for(Int_t ev = 0; ev < FASTMC->GetNev(); ++ev) {
    nev = ev;
    // Initialize the source
    FASTMC->Initialize(source, allocator);

      Npart = HYFPAR.npart;      
      Bgen = HYFPAR.bgen;
      Njet = HYJPAR.njet;
      Nbcol = HYFPAR.nbcol;
      if(ev==0) { 
            Sigin=HYJPAR.sigin;
            Sigjet=HYJPAR.sigjet;
         }
     std::cout<<"in RunHadronSource: ev"<<ev<<" Njet "<<Njet<<" Nbcol "<<Nbcol<<" Npart "<<Npart<<std::endl;

    if(source.empty()) {
      Error("RunHadronSource::main", "Source is not initialized!!");
      return 0;
    }
    
    // Run the decays
    if(FASTMC->GetTime() >= 0.)
      FASTMC->Evolve(source, secondaries, allocator, FASTMC->GetWeakDecayLimit());
   
    std::cout << "event #" << ev << "\r" << std::flush;
    Ntot = 0;
    Npyt = 0;
    Nhyd = 0;
    Ntot = 0;

    LPIT_t it;
    LPIT_t e;

    // Fill the source tree

    Ntot = 0; Nhyd=0; Npyt=0;      
    for(it = source.begin(), e = source.end(); it != e; ++it) {
      TVector3 pos(it->Pos().Vect());
      TVector3 mom(it->Mom().Vect());
      Float_t m1 = it->TableMass();
      pdg[Ntot] = it->Encoding();
      Mpdg[Ntot] = -1;
      Px[Ntot] = mom[0];
      Py[Ntot] = mom[1];
      Pz[Ntot] = mom[2];
      E[Ntot] = TMath::Sqrt(mom.Mag2() + m1*m1);
      X[Ntot] = pos[0];
      Y[Ntot] = pos[1];
      Z[Ntot] = pos[2];
      T[Ntot] = it->T();
      type[Ntot] = it->GetType();
      if(type[Ntot]==0)Nhyd++;
      if(type[Ntot]==1)Npyt++;            
      Ntot++;
      if(Ntot > kMax)
        Error("in main:", "Ntot is too large %d", Ntot);
    }
    ti->Fill();
    
    std::cout<<"ti Ntot= " <<Ntot<<" Npyt= "<<Npyt<<" Nhyd=  "<<Nhyd<<std::endl;
    
    // Fill the decayed tree
    Ntot = 0; Nhyd=0; Npyt=0;      
    for(it = secondaries.begin(), e = secondaries.end(); it != e; ++it) {
      TVector3 pos(it->Pos().Vect());
      TVector3 mom(it->Mom().Vect());
      Float_t m1 = it->TableMass();
      pdg[Ntot] = it->Encoding();
      Mpdg[Ntot] = it->GetLastMotherPdg();
      Px[Ntot] = mom[0];
      Py[Ntot] = mom[1];
      Pz[Ntot] = mom[2];
      E[Ntot] =  TMath::Sqrt(mom.Mag2() + m1*m1);
      X[Ntot] = pos[0];
      Y[Ntot] = pos[1];
      Z[Ntot] = pos[2];
      T[Ntot] = it->T();
      type[Ntot] = it->GetType();
      if(type[Ntot]==0)Nhyd++;
      if(type[Ntot]==1)Npyt++;            
      Ntot++;
      if(Ntot > kMax)
        Error("in main:", "Ntot is too large %d", Ntot);
    }
    td->Fill();      
    
    std::cout<<"td Ntot= " <<Ntot<<" Npyt= "<<Npyt<<" Nhyd=  "<<Nhyd<<std::endl;

    allocator.FreeList(source);
    allocator.FreeList(secondaries);
  }



  
  // Close the output file by getting it from the tree object to avoid crashes when
  // the output file is automatically switched by root.
  TFile *saveFile = ti->GetCurrentFile();
  saveFile->cd();
  ti->Write();
  td->Write();
  saveFile->Close();
  //  ti->GetCurrentFile()->Close();
  
  clock_t stop;
  stop = clock();
  std::cout << "*********************************************" << std::endl;
  std::cout << "Execution time: " << (stop - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "*********************************************" << std::endl;


//new
  time_t  now1;
  struct tm  *ts1;
  char       buf1[80];
         
 // Get the current time
   time(&now1);
              
 // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
    ts1 = localtime(&now1);
    strftime(buf1, sizeof(buf1), "%a %Y-%m-%d %H:%M:%S %Z", ts1);
    printf("%s\n", buf1);
     
                                    
  


  //return 0;
}
