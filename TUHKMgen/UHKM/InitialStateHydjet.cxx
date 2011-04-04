//expanding localy equilibated fireball with volume hadron radiation
//thermal part: Blast wave model, Bjorken-like parametrization
//hyght-pt: PYTHIA + jet quenching model PYQUEN
//                                                                           
//         HYDJET++ 
//         version 1.0:  
//         InitialStateHydjet is the modified InitialStateBjorken 
//         The high-pt part related with PYTHIA-PYQUEN is included       
//         InitialStateBjorken (FASTMC) was used.
//
//
//         
//         InitialStateBjorken           
//         version 2.0: 
//         Ludmila Malinina  malinina@lav01.sinp.msu.ru,   SINP MSU/Moscow and JINR/Dubna
//         Ionut Arsene  i.c.arsene@fys.uio.no,            Oslo University                                                
//                     June 2007
//        
//         version 1.0:                                                               
//         Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//         amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
//                           November. 2, 2005 

                     


#include <iostream> 
#include <fstream>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TMath.h>

#include "InitialStateHydjet.h"
#include "RandArrayFunction.h"
#include "GrandCanonical.h"
#include "StrangePotential.h"
#include "Particle.h"
#include "ParticlePDG.h"
#include "HYJET_COMMONS.h"

extern "C" void  hyevnt_();
extern "C" void  myini_();
extern HYIPARCommon HYIPAR;
extern HYFPARCommon HYFPAR;
extern HYJPARCommon HYJPAR;
extern HYPARTCommon HYPART;
extern SERVICECommon SERVICE;

using std::cout;
using std::endl;

class ParticleAllocator;
class TRandom3;

// declaration of the static member fLastIndex
Int_t Particle::fgLastIndex;

//_________________________________________________________________________________
void InitialStateHydjet::Initialize(List_t &source, ParticleAllocator & allocator) {
  // Generate initial particles from the soft and hard components

  // Initialize the static "last index variable"
  Particle::InitIndexing(); 

  //----- high-pt part------------------------------
  TLorentzVector partJMom, partJPos, zeroVec;

  // run a HYDJET event
  hyevnt_(); 

  fBgen  =  HYFPAR.bgen * HYIPAR.RA;
  fNpart =  HYFPAR.npart;
  fNcoll =  HYFPAR.nbcol;

  
      
  if(fParams.fNhsel != 0) {   
    //get number of particles in jets
    Int_t numbJetPart = HYPART.njp;

    for(Int_t i = 0; i <numbJetPart; i++) {
      Int_t pdg = Int_t(HYPART.ppart[i][1]);
      Double_t px = HYPART.ppart[i][2];
      Double_t py = HYPART.ppart[i][3];
      Double_t pz = HYPART.ppart[i][4];
      Double_t e =  HYPART.ppart[i][5];
      Double_t vx = HYPART.ppart[i][6];
      Double_t vy = HYPART.ppart[i][7];
      Double_t vz = HYPART.ppart[i][8];
      Double_t vt = HYPART.ppart[i][9];    
      ParticlePDG *partDef = fDatabase->GetPDGParticle(pdg);
      Int_t type =1;                //from jet
      if(partDef) {
	partJMom.SetXYZT(px, py, pz, e);
	partJPos.SetXYZT(vx, vy, vz, vt);
	Particle *particle=new Particle(partDef, partJPos, partJMom, 0, 0, type, -1, zeroVec, zeroVec);
	particle->SetIndex();
	allocator.AddParticle(*particle, source);
	delete particle;
      }
    }
  }       //nhsel !=0 not only hydro!             

         
  //----------HYDRO part------------------------------------------------
  if(fParams.fNhsel < 3) {
    const Double_t  weightMax = 2*TMath::CosH(fParams.fUmax);
    const Int_t nBins = 100;
    Double_t probList[nBins];
    RandArrayFunction arrayFunctDistE(nBins);
    RandArrayFunction arrayFunctDistR(nBins);
    
    TLorentzVector partPos, partMom, n1, p0;
    TVector3 vec3;
    //set maximal hadron energy
    const Double_t eMax = 5.;  
    //-------------------------------------
    // get impact parameter    
    
    //effective volume for central     
    double dYl= 2 * fParams.fYlmax; //uniform distr. [-Ylmax; Ylmax]  
    if(fParams.fEtaType >0) dYl = TMath::Sqrt(2 * TMath::Pi()) * fParams.fYlmax ;  //Gaussian distr. 
    Double_t volEffcent = 2 * TMath::Pi() * fParams.fTau * dYl * 
    (fParams.fR * fParams.fR)/TMath::Power((fParams.fUmax),2)*
    ((fParams.fUmax)*TMath::SinH((fParams.fUmax))-TMath::CosH((fParams.fUmax))+ 1);
 
    //effective volume for non-central Simpson2 
    Double_t volEffnoncent = fParams.fTau * dYl * SimpsonIntegrator2(0., 2.*TMath::Pi());
    fVolEff = volEffcent * HYFPAR.npart/HYFPAR.npart0;

    Double_t coeffRB = TMath::Sqrt(volEffcent * HYFPAR.npart/HYFPAR.npart0/volEffnoncent);
    Double_t coeffR1 = HYFPAR.npart/HYFPAR.npart0;
    coeffR1 = TMath::Power(coeffR1, 0.333333);

    double veff=fVolEff;
    
    //------------------------------------
    //cycle on particles types
    for(Int_t i = 0; i < fParams.fNPartTypes; ++i) {
      Double_t mparam = fParams.fPartMult[2 * i] * veff;
      Int_t multiplicity = gRandom->Poisson(mparam);
      const Int_t encoding = fParams.fPartEnc[i];

      if(multiplicity > 0) {
	ParticlePDG *partDef = fDatabase->GetPDGParticle(encoding);
	if(!partDef) {
	  Error("InitialStateHydjet::Initialize", "No particle with encoding %d", encoding);
	  continue;
	}
	//no charm now !
	if(TMath::Abs(partDef->GetCharmQNumber())>0 || TMath::Abs(partDef->GetCharmAQNumber())>0){
	  continue;
	}

	//compute chemical potential for single f.o. mu==mu_ch
	//compute chemical potential for thermal f.o.                
	Double_t mu = fParams.fPartMu[2 * i];
	
	//choose Bose-Einstein or Fermi-Dirac statistics
	const Double_t d    = !(Int_t(2*partDef->GetSpin()) & 1) ? -1 : 1;
	const Double_t mass = partDef->GetMass();                	 
	
	//prepare histogram to sample hadron energy: 
	Double_t h = (eMax - mass) / nBins;
	Double_t x = mass + 0.5 * h;
        Int_t ii; 
	for(ii = 0; ii < nBins; ++ii) {
	  if(x>=mu && fParams.fThFO>0)probList[ii] = x * TMath::Sqrt(x * x - mass * mass) / 
					(TMath::Exp((x - mu) / (fParams.fThFO)) + d);
	  if(x>=mu && fParams.fThFO<=0)probList[ii] = x * TMath::Sqrt(x * x - mass * mass) / 
					 (TMath::Exp((x - mu) / (fParams.fT)) + d);								 
	  if(x<mu)probList[ii] = 0.; 
	  x += h;
	}
	arrayFunctDistE.PrepareTable(probList);

	//prepare histogram to sample hadron transverse radius: 
	h = (fParams.fR) / nBins;
	x =  0.5 * h;
	Double_t param = (fParams.fUmax) / (fParams.fR);
	for(ii = 0; ii < nBins; ++ii) {
	  probList[ii] = x * TMath::CosH(param*x);
	  x += h;
	}
	arrayFunctDistR.PrepareTable(probList);

	//loop over hadrons, assign hadron coordinates and momenta
	Double_t weight = 0., yy = 0., px0 = 0., py0 = 0., pz0 = 0.;
	Double_t e = 0., x0 = 0., y0 = 0., z0 = 0., t0 = 0., etaF = 0.; 
	Double_t r, rB, phiF;
      
	for(Int_t j = 0; j < multiplicity; ++j) {               
	  do {
	    fParams.fEtaType <=0 ? etaF = fParams.fYlmax * (2. * gRandom->Rndm() - 1.) 
	      : etaF = (fParams.fYlmax) * (gRandom->Gaus());                                             
	    n1.SetXYZT(0.,0.,TMath::SinH(etaF),TMath::CosH(etaF));  
	    if(TMath::Abs(etaF)>5.)continue;
	    
	    rB = fParams.fR * coeffRB * coeffR1;
	        
	    Double_t rho = TMath::Sqrt(gRandom->Rndm());
	    Double_t phi = TMath::TwoPi() * gRandom->Rndm();
	    Double_t rx =  TMath::Sqrt(1-fParams.fEpsilon)*rB; 
	    Double_t ry =  TMath::Sqrt(1+fParams.fEpsilon)*rB;
	    
	    x0 = rx * rho * TMath::Cos(phi);
	    y0 = ry * rho * TMath::Sin(phi);
	    r = TMath::Sqrt(x0*x0+y0*y0);
	    phiF = TMath::Abs(TMath::ATan(y0/x0));
	    
	    if(x0<0&&y0>0)phiF = TMath::Pi()-phiF;
	    if(x0<0&&y0<0)phiF = TMath::Pi()+phiF;
	    if(x0>0&&y0<0)phiF = 2.*TMath::Pi()-phiF;
	    
	    //proper time with emission duration                                                               
	    Double_t tau = coeffR1 * fParams.fTau +  sqrt(2.) * fParams.fSigmaTau * coeffR1 * (gRandom->Gaus()); 	  
	    z0 = tau  * TMath::SinH(etaF);                                                                             
	    t0 = tau  * TMath::CosH(etaF);
	  
	    Double_t rhou = fParams.fUmax * r / rB;

        
	    Double_t uxf = TMath::SinH(rhou)*TMath::Sqrt(1+fParams.fDelta)*TMath::Cos(phiF); 
	    Double_t uyf = TMath::SinH(rhou)*TMath::Sqrt(1-fParams.fDelta)*TMath::Sin(phiF);
	    Double_t utf = TMath::CosH(etaF) * TMath::CosH(rhou) * 
	      TMath::Sqrt(1+fParams.fDelta*TMath::Cos(2*phiF)*TMath::TanH(rhou)*TMath::TanH(rhou));
	    Double_t uzf = TMath::SinH(etaF) * TMath::CosH(rhou) * 
	      TMath::Sqrt(1+fParams.fDelta*TMath::Cos(2*phiF)*TMath::TanH(rhou)*TMath::TanH(rhou));
        
	    vec3.SetXYZ(uxf / utf, uyf / utf, uzf / utf);
	    n1.Boost(-vec3); 
    
	    yy = weightMax * gRandom->Rndm();        
               
	    Double_t php0 = TMath::TwoPi() * gRandom->Rndm();
	    Double_t ctp0 = 2. * gRandom->Rndm() - 1.;
	    Double_t stp0 = TMath::Sqrt((1.-ctp0)*(1.+ctp0)); 
	    e = mass + (eMax - mass) * arrayFunctDistE(); 
	    Double_t pp0 = TMath::Sqrt((e-mass)*(e+mass));
	    px0 = pp0 * stp0 * TMath::Sin(php0); 
	    py0 = pp0 * stp0 * TMath::Cos(php0);
	    pz0 = pp0 * ctp0;
	    p0.SetXYZT(px0, py0, pz0, e);
	    
	    //weight for rdr          
	    weight = (n1 * p0) /e;  // weight for rdr gammar: weight = (n1 * p0) / n1[3] / e; 
	  } while(yy >= weight); 
	
          partMom.SetXYZT(px0, py0, pz0, e);
	  partPos.SetXYZT(x0, y0, z0, t0);
	  partMom.Boost(vec3);
	  Int_t type =0; //hydro

	  Particle *particle=new Particle(partDef, partPos, partMom, 0., 0, type, -1, zeroVec, zeroVec);
	  particle->SetIndex();
	  allocator.AddParticle(*particle, source);
	  delete particle;
        } //nhsel==4 , no hydro part
      }
    } 
  }
  
}

//_________________________________________________________________________________
Bool_t InitialStateHydjet::ReadParams() {     
  // Read parameters from an input file in ascii 
 
  Float_t par[200] = {0.};
  Int_t i = 0; 
  std::string s(40,' '); 
  std::ifstream input("RunInputHydjet");
  if (!input) {
    Error("Ukm::ReadParams", "Cannot open RunInputHydjet");
    return kFALSE;
  }
  
  while (std::getline(input, s)) {
    input>>par[i];
    if (i < 140) 
      std::cout<<s<<"     =  "<<par[i]<<std::endl;
    ++i;
    std::getline(input,s);
  }
  
  std::cout<<"\nFor output use the files RunOutput.root  \n\n"<< std::endl; 
  
  fParams.fNevnt  = Int_t(par[0]); //number of events
  fParams.fSqrtS  = par[1];        //cms energy per nucleon
  fParams.fAw     = par[2];        // atomic number of colliding nuclei
  fParams.fIfb    = Int_t(par[3]);      // flag of type of centrality generation (=0 is fixed by fBfix, not 0 
                                         //impact parameter is generated in each event between fBfmin 
                                          //and fBmax according with Glauber model (f-la 30)
  fParams.fBmin = par[4];         //minimum impact parameter in units of nuclear radius RA 
  fParams.fBmax = par[5];         //maximum impact parameter in units of nuclear radius RA
  fParams.fBfix = par[6];         //fix impact parameter in units of nuclear radius RA
  
  fParams.fSeed = Int_t(par[7]);         //parameter to set the random nuber seed (=0 the current time is used
                                   //to set the random generator seed, !=0 the value fSeed is 
                                   //used to set the random generator seed and then the state of random
                                   //number generator in PYTHIA MRPY(1)=fSeed
       
  fParams.fT         = par[8];     //chemical freeze-out temperature in GeV    
  fParams.fMuB       = par[9];     //baryon potential 
  fParams.fMuS       = par[10];    //strangeness potential 
  fParams.fMuI3      = par[11];    //isospin potential   
  fParams.fThFO      = par[12];    //thermal freeze-out temperature T^th in GeV
  fParams.fMu_th_pip = par[13];    // effective chemical potential of positivly charged pions at thermal in GeV 

       
  fParams.fTau       = par[14];     //proper time value
  fParams.fSigmaTau  = par[15];     //its standart deviation (emission duration)
  fParams.fR         = par[16];     //maximal transverse radius 
  fParams.fYlmax     = par[17];     //maximal longitudinal rapidity 
  fParams.fUmax      = par[18];     //maximal transverse velocity multiplaed on \gamma_r 
  fParams.fDelta     = par[19];     //momentum asymmetry parameter
  fParams.fEpsilon   = par[20];     //coordinate asymmetry parameter
  
  fParams.fDecay      = Int_t(par[21]);    // flag to switch on/off hadron decays<0: decays off,>=0: decays on, (default: 0)
  fParams.fWeakDecay  = Int_t(par[22]);    //flag to switch on/off weak hadron decays <0: decays off, >0: decays on, (default: 0)
  
  fParams.fEtaType   = Int_t(par[23]);     // flag to choose rapidity distribution, if fEtaType<=0, 
                                     //then uniform rapidity distribution in [-fYlmax,fYlmax] if fEtaType>0,
                                     //then Gaussian with dispertion = fYlmax 
  
  fParams.fTMuType   = Int_t(par[24]);     // flag to use calculated chemical freeze-out temperature,
                                     //baryon potential and strangeness potential as a function of fSqrtS 

  fParams.fCorrS     = par[25];     // flag and value to include strangeness supression factor    
  fParams.fNhsel = Int_t(par[26]);         //flag to switch on/off jet and hydro-state production (0: jet
                                     // production off and hydro on, 1: jet production on and jet quenching
                                     // off and hydro on, 2: jet production on and jet quenching on and
                                     // hydro on, 3: jet production on and jet quenching off and hydro
                                     // off, 4: jet production on and jet quenching on and hydro off
 
  fParams.fIshad= Int_t(par[27]);         //flag to switch on/off impact parameter dependent nuclear
                                    // shadowing for gluons and light sea quarks (u,d,s) (0: shadowing off,
                                    // 1: shadowing on for fAw=207, 197, 110, 40, default: 1
  
  fParams.fPtmin = par[28];       //minimal transverse momentum transfer p_T of hard
                                   // parton-parton scatterings in GeV (the PYTHIA parameter ckin(3)=fPtmin)
   
  //  PYQUEN energy loss model parameters:
 
  fParams.fT0 = par[29];          // initial temperature (in GeV) of QGP for
                                   //central Pb+Pb collisions at mid-rapidity (initial temperature for other
                                  //centralities and atomic numbers will be calculated automatically) (allowed range is 0.2<fT0<2) 
  
  fParams.fTau0= par[30];        //proper QGP formation time in fm/c (0.01<fTau0<10)
  fParams.fNf= Int_t(par[31]);          //number of active quark flavours N_f in QGP fNf=0, 1,2 or 3 
  fParams.fIenglu= Int_t(par[32]);      // flag to fix type of in-medium partonic energy loss 
                                           //(0: radiative and collisional loss, 1: radiative loss only, 2:
                                  //collisional loss only) (default: 0);
  fParams.fIanglu= Int_t(par[33]);      //flag to fix type of angular distribution of in-medium emitted
                                  // gluons (0: small-angular, 1: wide-angular, 2:collinear) (default: 0).


  //PYTHIA parameters:
  Int_t jj;
  for (Int_t j = 0; j <25; ++j) {
    jj= 35+j;
    SERVICE.parPYTH[j]=par[jj];
  } 

  // Set Random Number seed 
        
  gRandom->SetSeed(fParams.fSeed); //Set 0 to use the current time
  //to send seed in PYTHIA
  SERVICE.iseed_fromC=gRandom->GetSeed(); 
  std::cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<std::endl;  
  
  fParams.fNPartTypes = 0;         //counter of hadron species
  return kTRUE; 
}

//_________________________________________________________________________________
Bool_t InitialStateHydjet::MultIni() {
  // Calculate average multiplicities, chemical potentials (if necessary),
  // initialize pyquen 

  //check and redefine input parameters
  if(fParams.fTMuType>0 &&  fParams.fSqrtS > 2.24) {
    if(fParams.fSqrtS < 2.24){
      Error("InitialStateHydjet::MultIni", "SqrtS<2.24 not allowed with fParams.fTMuType>0");
      return 0;
    }
    
    //sqrt(s) = 2.24 ==> T_kin = 0.8 GeV
    //see J. Cleymans, H. Oeschler, K. Redlich,S. Wheaton, Phys Rev. C73 034905 (2006)
    fParams.fMuB = 1.308/(1. + fParams.fSqrtS*0.273);
    fParams.fT = 0.166 - 0.139*fParams.fMuB*fParams.fMuB - 0.053*fParams.fMuB*fParams.fMuB*
      fParams.fMuB*fParams.fMuB;
    fParams.fMuI3 = 0.;
    fParams.fMuS = 0.;
    //create strange potential object and set strangeness density 0
    StrangePotential psp(0., fDatabase);
    psp.SetBaryonPotential(fParams.fMuB);
    psp.SetTemperature(fParams.fT);
    //compute strangeness potential
    if(fParams.fMuB > 0.01)
      fParams.fMuS = psp.CalculateStrangePotential();
    
    //if user choose fYlmax larger then allowed by kinematics at the specified beam energy sqrt(s)     
    if(fParams.fYlmax > TMath::Log(fParams.fSqrtS/0.94)){
      Error("InitialStateHydjet::MultIni", "fParams.fYlmax > TMath::Log(fParams.fSqrtS/0.94)!!! ");
      return 0;
    }
    
    
    if(fParams.fCorrS <= 0.) {
      //see F. Becattini, J. Mannien, M. Gazdzicki, Phys Rev. C73 044905 (2006)
      fParams.fCorrS = 1. - 0.386* TMath::Exp(-1.23*fParams.fT/fParams.fMuB);
      
    }
    std::cout<<"The phenomenological f-la J. Cleymans et al. PRC73 034905 (2006) for Tch mu_B was used." << std::endl;
    std::cout<<"The simulation will be done with the calculated parameters:" << std::endl;
    std::cout<<"Baryon chemical potential = "<<fParams.fMuB<< " [GeV]" << std::endl;
    std::cout<<"Strangeness chemical potential = "<<fParams.fMuS<< " [GeV]" << std::endl;
    std::cout<<"Isospin chemical potential = "<<fParams.fMuI3<< " [GeV]" << std::endl;
    std::cout<<"Strangeness suppression parameter = "<<fParams.fCorrS << std::endl;
    std::cout<<"Eta_max = "<<fParams.fYlmax<<  std::endl;
    std::cout << std::endl;
    
  }
  
  
  //initialisation of high-pt part 
  
  HYJPAR.nhsel = fParams.fNhsel;
  HYJPAR.ptmin = fParams.fPtmin;
  HYJPAR.ishad = fParams.fIshad;
  HYIPAR.bminh = fParams.fBmin;
  HYIPAR.bmaxh = fParams.fBmax;
  HYIPAR.AW = fParams.fAw;
  
  HYPYIN.ifb = fParams.fIfb;
  HYPYIN.bfix = fParams.fBfix;
  HYPYIN.ene = fParams.fSqrtS;
  
  PYQPAR.T0 = fParams.fT0;
  PYQPAR.tau0 = fParams.fTau0;
  PYQPAR.nf = fParams.fNf;
  PYQPAR.ienglu = fParams.fIenglu;
  PYQPAR.ianglu = fParams.fIanglu;
  
  
  myini_();  //
  
  
  // calculation of  multiplicities of different particle species
  // according to the grand canonical approach
  GrandCanonical gc(15, fParams.fT, fParams.fMuB, fParams.fMuS, fParams.fMuI3);
  GrandCanonical gcCh(15, fParams.fT, fParams.fMuB, fParams.fMuS, fParams.fMuI3);
  GrandCanonical gcPiTh(15, fParams.fThFO, 0., 0., fParams.fMu_th_pip);
  GrandCanonical gcTh0(15, fParams.fThFO, 0., 0., 0.);
  
  //effective volume for central     
  double dYl= 2 * fParams.fYlmax; //uniform distr. [-Ylmax; Ylmax]  
  if (fParams.fEtaType >0) dYl = TMath::Sqrt(2 * TMath::Pi()) * fParams.fYlmax ;  //Gaussian distr.                                                                            
  fVolEff = 2 * TMath::Pi() * fParams.fTau * dYl * (fParams.fR * fParams.fR)/TMath::Power((fParams.fUmax),2) * 
    ((fParams.fUmax)*TMath::SinH((fParams.fUmax))-TMath::CosH((fParams.fUmax))+ 1);
  
  Double_t particleDensityPiCh=0;
  Double_t particleDensityPiTh=0;
  
  if(fParams.fThFO != fParams.fT && fParams.fThFO > 0){
    particleDensityPiCh = gcCh.ParticleNumberDensity(fDatabase->GetPDGParticle(211));
    particleDensityPiTh = gcPiTh.ParticleNumberDensity(fDatabase->GetPDGParticle(211));
  }

  if(fDatabase->GetNParticles()>=kNPartTypes) {
    cout << "InitialStateHydjet::MultIni(): ERROR Particle definitions in the PDG database exceeds the hardcoded limit of " << kNPartTypes << endl;
    cout << "      There is either an error with reading the particles file or you might need to increase the maximum allowed definitions" << endl;
  }
  for(Int_t particleIndex = 0; particleIndex<fDatabase->GetNParticles() && particleIndex<kNPartTypes; particleIndex++) {
    ParticlePDG *currParticle = fDatabase->GetPDGParticleByIndex(particleIndex);
    Int_t encoding = currParticle->GetPDG();
    //strangeness supression
    Double_t gammaS = 1;
    Int_t s = Int_t(currParticle->GetStrangeness());
    if(encoding == 333) 
      s = 2;
    if(fParams.fCorrS < 1. && s != 0)
      gammaS = TMath::Power(fParams.fCorrS,-TMath::Abs(s));
    //average densities      
    Double_t particleDensity = gc.ParticleNumberDensity(currParticle)/gammaS;
    
    //compute chemical potential for single f.o. mu==mu_ch
    Double_t mu = fParams.fMuB  * Int_t(currParticle->GetBaryonNumber()) + 
      fParams.fMuS  * Int_t(currParticle->GetStrangeness()) +
      fParams.fMuI3 * Int_t(currParticle->GetElectricCharge());
    
    //thermal f.o.
    if(fParams.fThFO != fParams.fT && fParams.fThFO > 0){
      Double_t particleDensityCh = gcCh.ParticleNumberDensity(currParticle);
      Double_t particleDensityTh0 = gcTh0.ParticleNumberDensity(currParticle);
      Double_t numbDensBolt = particleDensityPiTh*particleDensityCh/particleDensityPiCh;               
      mu = fParams.fThFO*TMath::Log(numbDensBolt/particleDensityTh0);
      if(abs(encoding)==211 || encoding==111)mu= fParams.fMu_th_pip; 
      particleDensity = numbDensBolt;         
    }
    
    // set particle number density to zero for some species
    // photons
    if(abs(encoding)==22)
      particleDensity=0;
    // K0L and K0S
    if(abs(encoding)==130 || abs(encoding)==310) {
      particleDensity=0;
    }    
    
    if(particleDensity > 0.) {
      //      outMult<<encoding<< "         " <<particleDensity<< "      "<<mu<<std::endl;
      fParams.fPartEnc[fParams.fNPartTypes] = encoding;
      fParams.fPartMult[2 * fParams.fNPartTypes] = particleDensity;
      fParams.fPartMu[2 * fParams.fNPartTypes] = mu;
      ++fParams.fNPartTypes;
      if(fParams.fNPartTypes > 1000)
	Error("in Bool_t MultIni:", "fNPartTypes is too large %d", fParams.fNPartTypes);
    }
  }
  return kTRUE;
}

//_________________________________________________________________________________
Double_t InitialStateHydjet::SimpsonIntegrator2(Double_t a, Double_t b) {
  // Simpson integration
  Int_t nsubIntervals=10000;
  Double_t h = (b - a)/nsubIntervals; //0-pi, phi
  Double_t s=0;
  Double_t x = 0; //phi
  for(Int_t j = 1; j < nsubIntervals; j++) {
    x += h; // phi
    Double_t e = fParams.fEpsilon;
    Double_t rSB = fParams.fR; //test: podstavit' *coefff_RB
    Double_t rB = rSB *(TMath::Sqrt(1-e*e)/TMath::Sqrt(1+e*TMath::Cos(2*x))); //f-la7 rB    
    Double_t sr = SimpsonIntegrator(0,rB,x);
    s += sr;
  }
  return s*h;
  
}

//_________________________________________________________________________________
Double_t InitialStateHydjet::SimpsonIntegrator(Double_t a, Double_t b, Double_t phi) {
  // Simpson integration
  Int_t nsubIntervals=100;
  Double_t h = (b - a)/nsubIntervals;
  Double_t s = F2(phi,a + 0.5*h);
  Double_t t = 0.5*(F2(phi,a) + F2(phi,b));
  Double_t x = a;
  Double_t y = a + 0.5*h;
  for(Int_t i = 1; i < nsubIntervals; i++) {
    x += h;
    y += h;
    s += F2(phi,y);
    t += F2(phi,x);
  }	
  t += 2.0*s;
  return t*h/3.0;
}


//f2=f(phi,r)
//_________________________________________________________________________________
Double_t InitialStateHydjet::F2(Double_t x, Double_t y) {
  // formula
  Double_t rSB = fParams.fR; //test: podstavit' *coefff_RB
  Double_t rhou =  fParams.fUmax * y / rSB;
  Double_t ff = y*TMath::CosH(rhou)*
    TMath::Sqrt(1+fParams.fDelta*TMath::Cos(2*x)*TMath::TanH(rhou)*TMath::TanH(rhou));
//n_mu u^mu f-la 20
  return ff;
}

//_________________________________________________________________________________
Double_t InitialStateHydjet::MidpointIntegrator2(Double_t a, Double_t b) {
  // Perform integration through the mid-point method
  Int_t nsubIntervals=2000; 
  Int_t nsubIntervals2=1; 
  Double_t h = (b - a)/nsubIntervals; //0-pi , phi
  Double_t h2 = (fParams.fR)/nsubIntervals; //0-R maximal rB ?
  Double_t x = a + 0.5*h;
  Double_t y = 0;
  Double_t t = F2(x,y);                    
  Double_t e = fParams.fEpsilon;
  for(Int_t j = 1; j < nsubIntervals; j++) {
    x += h; // integr  phi
    Double_t rSB = fParams.fR; //test: podstavit' *coefff_RB
    Double_t  rB = rSB *(TMath::Sqrt(1-e*e)/TMath::Sqrt(1+e*TMath::Cos(2*x))); //f-la7 rB
    nsubIntervals2 = Int_t(rB / h2)+1;
    // integr R 
    y=0;
    for(Int_t i = 1; i < nsubIntervals2; i++) 
      t += F2(x,(y += h2));
  }
  return t*h*h2;
}

