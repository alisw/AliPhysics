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
#include <TH1D.h>
#include <TFile.h>

#include "InitialState.h"
#include "InitialStateHydjet.h"


#include <TRandom.h>

#include "Particle.h"
//#include "HYJET_COMMONS.h"
//extern SERVICECommon SERVICE;


//Main program:
//reads input parameters from file "RunInputBjorken" or "RunInputHubble";
//calculates particle densities and average initial multiplicities and writes them
//in output file "multiplicities.txt";
//creates trees (tree with direct hadrons and hadrons after resonance decays)
//with space-time and momentum-energy information of produced hadrons;
//writes trees in file "RunOutput.root".

Int_t main() {

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
  //define hadron number
  Int_t ntot;
  //define event number
  Int_t nev;
  //define hadron characteristic vectors
  std::vector<Int_t> pdg(kMax); //pdg encodings
  std::vector<Int_t> Mpdg(kMax);//pdg encodings for mother hadrons
  std::vector<Int_t> type(kMax);//type: 0-from hydro or decay, 1111 from jets
  std::vector<Float_t> Px(kMax);//x-hadron momentum component,[GeV/c]
  std::vector<Float_t> Py(kMax);//y-hadron momentum component,[GeV/c]
  std::vector<Float_t> Pz(kMax);//z-hadron momentum component,[GeV/c]
  std::vector<Float_t> E(kMax); //hadron total energy,[GeV]  
  std::vector<Float_t> X(kMax);//x-hadron coordinate component,[fm]
  std::vector<Float_t> Y(kMax);//y-hadron coordinate component,[fm]
  std::vector<Float_t> Z(kMax);//z-hadron coordinate component,[fm]
  std::vector<Float_t> T(kMax);//hadron time,[fm/c] 

   TH1D *hpt1 = new TH1D("hpt1", "hpt1", 100, 0., 20.);
   TH1D *hpt1j = new TH1D("hpt1j", "hpt1j", 100, 0., 20.);
   TH1D *hpt1h = new TH1D("hpt1h", "hpt1h", 100, 0., 20.);

   TH1D *hv2 = new TH1D("hv2", "hv2", 100, 0.0, 10.);
   TH1D *hv0 = new TH1D("hv0", "hv0", 100, 0.0, 10.);

   TH1D *hy = new TH1D("hy", "hy", 51, -5.1, 5.1);
   TH1D *hyjets = new TH1D("hyjets", "hyjets", 51, -5.1, 5.1);
   TH1D *hyhydro = new TH1D("hyhydro", "hyhydro", 51, -5.1, 5.1);


   double pdg1, Mpdg1, Px1, Py1, E1, Z1, Pz1, pt, phi, v2, eta;
   int type1;

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
    
  
  // Loop over events  
  for(Int_t ev = 0; ev < FASTMC->GetNev(); ++ev) {
    nev = ev;
    // Initialize the source
    FASTMC->Initialize(source, allocator);
    if(source.empty()) {
      Error("RunHadronSource::main", "Source is not initialized!!");
      //return 0;
      continue;  
    }
    
    // Run the decays //fDecay
    if(FASTMC->GetTime() >= 0.) 
      FASTMC->Evolve(source, secondaries, allocator, FASTMC->GetWeakDecayLimit());
   
    std::cout << "event #" << ev << "\r" << std::flush;
//    npart = 0;
    LPIT_t it;
    LPIT_t e;
    
    // Fill the decayed tree
//    npart = 0;      
    
    for(it = secondaries.begin(), e = secondaries.end(); it != e; ++it) {
      TVector3 pos(it->Pos().Vect());
      TVector3 mom(it->Mom().Vect());
      Float_t m1 = it->TableMass();
      pdg1 = it->Encoding();
      Mpdg1 = it->GetLastMotherPdg();
      Px1 = mom[0];
      Py1 = mom[1];
      Pz1 = mom[2];
      E1 =  TMath::Sqrt(mom.Mag2() + m1*m1);
      type1 = it->GetType();
      if(pdg1==211 && abs(0.5*log((E1+Pz1)/(E1-Pz1)))<1.) {
      hpt1->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));
         }
      
      if(pdg1==211 && abs(0.5*log((E1+Pz1)/(E1-Pz1)))<1. && type1==0) hpt1h->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));
      if(pdg1==211 && abs(0.5*log((E1+Pz1)/(E1-Pz1)))<1. && type1==1)hpt1j->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));

      if(((abs(pdg1)==211)||(abs(pdg1)==321)||(abs(pdg1)==2212)) 
       && (abs(0.5*log((E1+Pz1)/(E1-Pz1)))<1.0)){
       pt = TMath::Sqrt(Px1*Px1+Py1*Py1);      
       phi = TMath::ATan2(Py1,Px1);
       v2 = TMath::Cos(2*phi);       
       hv2->Fill(pt,v2);
       hv0->Fill(pt,1.);
       }
       
       if((abs(pdg1)==211)||(abs(pdg1)==321)||(abs(pdg1)==2212)){    
       eta=0.5*TMath::Log((sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)+Pz1)/(sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)-Pz1));
       if(type1==1)hyjets->Fill(eta);
       if(type1==0)hyhydro->Fill(eta);
       hy->Fill(eta);
         }

     // npar++;
     // if(npart > kMax)
    //    Error("in main:", "npart is too large %d", npart);


    }
     
    allocator.FreeList(source);
    allocator.FreeList(secondaries);
 
 
  }
  
  hpt1->Write();
  hpt1h->Write();
  hpt1j->Write();
  hv2->Write();
  hv0->Write();
  hyhydro->Write();
  hyjets->Write();
  hy->Write();
  
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
    
    

  return 0;
}
