//#include <TClonesArray.h> 

#include <TDatabasePDG.h>
#include <TFile.h>
#include "AliConst.h"
#include "AliLog.h"
#include "AliGenMUONLMR.h" 
#include "AliMC.h" 
#include "AliRun.h" 
#include "AliGenEventHeader.h"

ClassImp(AliGenMUONLMR)


AliGenMUONLMR::AliGenMUONLMR () : AliGenMC(), fNMuMin(2), fGenSingleProc(-1), fCosTheta(0), fRhoLineShape(0), fHMultMu(0), fHNProc(0) { 
  //
  // default constructor 
  //
  // initialize pt and y distributions according to a fit to 
  // Pythia simulation at sqrt(s) = 7 TeV
  printf ("Using AliGenMUONLMR as generator\n");
  for (Int_t ipart=0; ipart < fgkNpart; ipart++) fScaleMult[ipart] = 1; 
  fScaleMult[kPionLMR] = 0; // set pion multiplicity to zero 
  fScaleMult[kKaonLMR] = 0; // set kaon multiplicity to zero
  Int_t pdg[7] = {211, 321, 221, 113, 223, 333, 331}; 
  const char* fptname[7] = {"fPtPion","fPtKaon","fPtEta","fPtRho","fPtOmega","fPtPhi","fPtEtaPrime"};
  const char* fyname[7] = {"fYPion","fYKaon","fYEta","fYRho","fYOmega","fYPhi","fYEtaPrime"}; 
  const char* fnname[7] = {"fMultPion","fMultKaon","fMultEta","fMultRho","fMultOmega","fMultPhi","fMultEtaPrime"};
  const char* fdname[2] = {"fDecPion","fDecKaon"};
  Double_t ptparam[7][3] = { {1,0.427,2.52}, // pions from Pythia
			     {1,0.58,2.57},  // kaons from Pythia
 			     {1,0.641,2.62},   // eta from Pythia
			     {1,1.2,2.5},    // rho+omega from ALICE muon  
			     {1,1.2,2.5},    // rho+omega from ALICE muon  
			     {1,1.03,2.5},   // phi from ALICE muon  
			     {1,0.72,2.5}};  // etaPrime from Pythia
  
  Double_t yparam[7][3] = { {1, 0.8251, 3.657}, // pions from pythia
			    {1, 1.83, 2.698},   // kaons from pythia
			    {1, 1.169, 3.282},  // eta from pythia
			    {1, 1.234, 3.264},  // rho from pythia
			    {1, 1.311, 3.223},  // omega from pythia
			    {1, 2.388, 2.129},  // phi from pythia
			    {1, 1.13,3.3}};     // eta prime from pythia

  // multiplicity parameters from pythia
  Double_t nparam[7][9] = { {353.582, 6.76263, 1.66979, 998.445, 9.73281, 12.6704, 175.187, 29.08, 40.2531},
			    {1.e4, 0.2841,  0,0,0,0,0,0,0},
			    {1.e4, 0.2647,  0,0,0,0,0,0,0},
			    {7055, 0.1786,  0,0,0,0,0,0,0},
			    {7500, 0.1896,  0,0,0,0,0,0,0},
			    {5.e04, 1.167,  0,0,0,0,0,0,0}, 
			    {2.9e04,0.714,  0,0,0,0,0,0,0}}; 

  Double_t ctau[2] = {7.8045, 3.712};

  for (Int_t i=0; i<fgkNpart; i++) { 
    fPDG[i] = pdg[i]; 
    if (i!=0) { 
      fMult[i] = new TF1(fnname[i],"[0]*exp(-[1]*x)",0,30);
      fMult[i]->SetParameters(nparam[i][0],nparam[i][1]);  
    }
    else { 
      fMult[i] = new TF1(fnname[i],"gaus(0)+gaus(3)+gaus(6)",0,150);
      for (Int_t j=0; j<9; j++) fMult[i]->SetParameter(j,nparam[i][j]);
    }
  
    fPt[i] = new TF1(fptname[i],AliGenMUONLMR::PtDistr,0,20,3);
    fPt[i]->SetParameters(ptparam[i][0], ptparam[i][1], ptparam[i][2]);  
    fY[i] = new TF1(fyname[i],AliGenMUONLMR::YDistr,-10,10,3);
    fY[i]->SetParameters(yparam[i][0], yparam[i][1], yparam[i][2]); 
  }

  for(Int_t i = 0; i<2; i++){
    fDecay[i] = new TF1(fdname[i],"exp(-x/[0])",0,150);
    fDecay[i]->SetParameter(0,ctau[i]);
  }
  
  for (Int_t ipart = 0; ipart < fgkNpart; ipart++) { 
    fParticle[ipart] = new TParticle(); 
    fParticle[ipart]->SetPdgCode(fPDG[ipart]); 
  }

  TDatabasePDG *pdgdb = TDatabasePDG::Instance(); 
  Double_t mumass = pdgdb->GetParticle(13)->Mass();
  fMu[0] = new TParticle(); 
  fMu[0]->SetPdgCode(-13); 
  fMu[0]->SetCalcMass(mumass); 
  fMu[1] = new TParticle(); 
  fMu[1]->SetPdgCode(13); 
  fMu[1]->SetCalcMass(mumass); 
  
  // function for polarized theta distributions
  fCosTheta = new TF1 ("fCosTheta","1+[0]*x*x",-1,1);
  fCosTheta->SetParameter(0,1);

  // Dalitz decays 
  Int_t nbins = 1000;
  Double_t xmin = 0, xmax = 2; 
  fDalitz[0] = new TH1F("hDalitzEta","",nbins,xmin,xmax);
  fDalitz[1] = new TH1F("hDalitzOmega","",nbins,xmin,xmax);
  fDalitz[2] = new TH1F("hDalitzEtaPrime","",nbins,xmin,xmax);
  
  Double_t meta   = pdgdb->GetParticle("eta")->Mass(); 
  Double_t momega = pdgdb->GetParticle("omega")->Mass(); 
  Double_t metaPrime = pdgdb->GetParticle("eta'")->Mass(); 
  Double_t mpi0   = pdgdb->GetParticle("pi0")->Mass(); 
  Double_t md3 = 0, mres = 0; 
  
  for (Int_t index = 0; index < 3; index++) { 
    if (index == 0) { 
      mres = meta; 
      md3 = 0; 
    }
    else if (index == 1) { 
      mres = momega; 
      md3 = mpi0; 
    }
    else if (index == 2) { 
      mres = metaPrime; 
      md3 = 0; 
    }
    Double_t delta   = md3 * md3 / (mres * mres);
    Double_t epsilon = mumass * mumass / (mres * mres);
    Int_t nbins0 = fDalitz[index]->GetNbinsX();
    Double_t xmin0 = fDalitz[index]->GetXaxis()->GetXmin(); 
    Double_t deltax =  fDalitz[index]->GetBinWidth(1);
    Double_t xd = xmin0 - deltax/2.; 
    for (Int_t ibin = 0; ibin< nbins0; ibin++) { 
      Double_t dalval = 0; 
      xd += deltax; 
      if (xd > 4. *epsilon) { 
	Double_t bracket = TMath::Power(1. + xd/(1. - delta),2)      
	  - 4. * xd / ((1. - delta) * (1. - delta));
	if (bracket > 0) { 
	  dalval = TMath::Power(bracket,1.5) /xd *
	    TMath::Sqrt(1 - 4 * epsilon / xd) * (1 + 2 * epsilon / xd) * 
	    FormFactor(xd * mres * mres, index);
	  fDalitz[index]->Fill(xd,dalval); 
	}
      }
    }
  }

  fRhoLineShape = new TF1("fRhoLineShape",RhoLineShapeNew,0,2,2); 
  fHMultMu = new TH1D("fHMultMu","Muon multiplicity",20,-0.5,19.5); 
  fHNProc = new TH1D("fHNProc","Number of gen. evts. per process in 4 pi",9,-0.5,8.5); 
}

//-----------------------------------------------------------

AliGenMUONLMR::~AliGenMUONLMR()
{
  // Default destructor
  for (Int_t i=0; i<7; i++) { 
    delete fPt[i]; 
    delete fY[i]; 
    delete fMult[i]; 
    delete fParticle[i]; 
  }    
  
  for (Int_t i=0; i<2; i++) { 
    delete fDecay[i]; 
    delete fMu[i]; 
  }

  for (Int_t i=0; i<3; i++) delete fDalitz[i]; 

  delete fCosTheta; fCosTheta = 0;  
  delete fRhoLineShape; fRhoLineShape = 0;  
  delete fHMultMu; fHMultMu = 0;   
  delete fHNProc;  fHNProc = 0;   
}

//-----------------------------------------------------------

void AliGenMUONLMR::FinishRun(){ 
  // save some histograms to an output file 
  Int_t nbins = fHNProc->GetNbinsX(); 
  for (Int_t ibin=1; ibin <= nbins; ibin++) printf ("ibin = %d nEvProc = %g\n",
						      ibin,fHNProc->GetBinContent(ibin));
    TFile *fout = new TFile("AliGenMUONLMR_histos.root","recreate"); 
  fHMultMu->Write(); 
  fHNProc->Write(); 
  fout->Close(); 
}

//-----------------------------------------------------------

Double_t AliGenMUONLMR::YDistr(const Double_t *px, const Double_t *par){ 
  // function for rapidity distribution: plateau at par[0] +
  // gaussian tails centered at par[1] and with par[2]=sigma  
  Double_t x = TMath::Abs(px[0]);
  Double_t func = 0;
  if (x<par[1]) func = par[0]; 
  else { 
    Double_t z = (x-par[1])/(par[2]); 
    func = par[0] * TMath::Exp(-0.5 * z * z); 
  }
  return func; 
}

//-----------------------------------------------------------

Double_t AliGenMUONLMR::PtDistr(const Double_t *px, const Double_t *par){
  // pt distribution: power law 
  Double_t x = px[0];
  Double_t func = par[0] * x / TMath::Power((1+(x/par[1])*(x/par[1])),par[2]); 
  return func; 
}

//-----------------------------------------------------------

void AliGenMUONLMR::Generate() {
  //
  // generate the low mass resonances and their decays according to  
  // the multiplicity parameterized by pythia and BR from PDG  
  // rapidity distributions parametrized from pythia 
  // pt distributions from data (or pythia for etaprime) 
  //
  Double_t pxPushed[100], pyPushed[100], pzPushed[100], ePushed[100]; 
  Int_t nmuons = -1, npartPushed = 0, pdgPushed[100]; 
  Double_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
  Double_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
  Double_t time0;              // Time0  of the generated parent particle
  // Calculating vertex position per event
  for (Int_t j=0;j<3;j++) origin0[j]=fOrigin[j];
  time0 = fTimeOrigin;
  if(fVertexSmear==kPerEvent) {
    Vertex();
    for (Int_t j=0;j<3;j++) origin0[j]=fVertex[j];
    time0 = fTime;
  }
  
  printf ("In Generate()\n");
  TParticle *mother; 
  TDatabasePDG* pdg = TDatabasePDG::Instance();

  Double_t pt, y, phi, mass, px, py, pz, ene, mt; 

  const Int_t nproc = 9; 
  Int_t idRes[nproc] = {kEtaLMR, kEtaLMR, kRhoLMR, kOmegaLMR, kOmegaLMR, kPhiLMR, kEtaPrimeLMR, kPionLMR, kKaonLMR}; 
  Double_t BR[nproc] = {5.8e-6, 3.1e-4, 4.55e-5, 7.28e-4, 1.3e-4, 2.86e-4, 1.04e-4, 1, 0.6344};
  //  Double_t BR[nproc] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  Int_t idDec[nproc] = {0, 1, 0, 0, 1, 0, 1, 2, 2}; // 0:2body, 1:Dalitz, 2:pi/K 
  Int_t mult[nproc] = {0,0,0,0,0,0,0,0,0}; 

  while (nmuons < fNMuMin) { 

    nmuons = 0; 
    npartPushed = 0; 
    for (Int_t iproc=0; iproc<nproc; iproc++) { 
      if (fGenSingleProc == -1) { 
	mult[iproc] = Int_t(fMult[idRes[iproc]]->GetRandom()*fScaleMult[idRes[iproc]]); 
      }
      else { 
	if (iproc==fGenSingleProc) { 
	  mult[iproc] = 1; 
	  BR[iproc] = 1;
	} 
	else { 
	  mult[iproc] = 0; 
	  BR[iproc] = 0;
	}
      }
    }
    
    if (fGenSingleProc == -1) { 
      mult[1] = mult[0]; 
      mult[4] = mult[3]; 
    }
    
    for (Int_t iproc = 0; iproc < nproc; iproc++) { 
//       printf ("Multiplicity for process %d is %d\n",iproc,mult[iproc]); 
      for (Int_t imult=0; imult<mult[iproc]; imult++) { 
	if (gRandom->Rndm() < BR[iproc]) { 
	  fHNProc->Fill(iproc); 
	  Int_t ipart = idRes[iproc]; 
	  pt  = fPt[ipart]->GetRandom(); 
	  y   = fY[ipart]->GetRandom(); 
	  phi = gRandom->Rndm() * 2 * TMath::Pi(); 
	  mass = pdg->GetParticle(fPDG[ipart])->Mass(); 
	  px  = pt * TMath::Cos(phi); 
	  py  = pt * TMath::Sin(phi); 
	  mt  = TMath::Sqrt(pt * pt + mass * mass);
	  pz  = mt * TMath::SinH(y); 
	  ene = mt * TMath::CosH(y); 
	
	  mother = fParticle[ipart]; 
	  mother->SetMomentum(px,py,pz,ene); 
	  mother->SetCalcMass(mass);
	  if (!KinematicSelection(mother,0)) continue; 

	  Bool_t hasDecayed = kTRUE;
	  if (idDec[iproc] == 0) Decay2Body(mother);
	  else if (idDec[iproc] == 1) DalitzDecay(mother); 
	  else DecayPiK(mother,hasDecayed); 
	  if (!hasDecayed) continue; 
	  Bool_t isMu0Acc = KinematicSelection(fMu[0],1); 
	  Bool_t isMu1Acc = KinematicSelection(fMu[1],1); 
	  Bool_t isMuFromPiKAcc = kTRUE; 

	  if (idDec[iproc] == 2) isMuFromPiKAcc = (mother->GetPdgCode()>0) ? isMu0Acc : isMu1Acc;
	  // mother 
	  if ((idDec[iproc]  < 2 && (isMu0Acc || isMu1Acc)) || 
	      (idDec[iproc] == 2 && isMuFromPiKAcc)) { 
	    pdgPushed[npartPushed] = mother->GetPdgCode(); 
	    pxPushed[npartPushed] = mother->Px(); 
	    pyPushed[npartPushed] = mother->Py(); 
	    pzPushed[npartPushed] = mother->Pz();
	    ePushed[npartPushed] = mother->Energy(); 
	    npartPushed++; 
	    if (isMu0Acc && (idDec[iproc] < 2 || mother->GetPdgCode() > 0)) { 
	      pdgPushed[npartPushed] = fMu[0]->GetPdgCode(); 
	      pxPushed[npartPushed] = fMu[0]->Px(); 
	      pyPushed[npartPushed] = fMu[0]->Py(); 
	      pzPushed[npartPushed] = fMu[0]->Pz();
	      ePushed[npartPushed] = fMu[0]->Energy(); 
	      npartPushed++; 
	      nmuons++; 
	    }
	    
	    if (isMu1Acc && (idDec[iproc] < 2 || mother->GetPdgCode() < 0)) { 
	      pdgPushed[npartPushed] = fMu[1]->GetPdgCode(); 
	      pxPushed[npartPushed] = fMu[1]->Px(); 
	      pyPushed[npartPushed] = fMu[1]->Py(); 
	      pzPushed[npartPushed] = fMu[1]->Pz();
	      ePushed[npartPushed] = fMu[1]->Energy(); 
	      npartPushed++; 
	      nmuons++; 
	    }
	  } 
	} // end if BR
      } // end loop on multiplicity 
    }  // end loop on process 
    fHMultMu->Fill(nmuons); 
  } // keep on generating until at least a muon is created in the event
  
  Int_t ntmother = 0, ntchild =0; 
  for (Int_t ipart = 0; ipart < npartPushed; ipart++) { 
    if (TMath::Abs(pdgPushed[ipart]) != 13) { // particle is not a muon, hence it's a mother
      PushTrack(0,-1,pdgPushed[ipart],
		pxPushed[ipart],pyPushed[ipart],pzPushed[ipart],ePushed[ipart],
		origin0[0],origin0[1],origin0[2],time0,
		polar[0],polar[1],polar[2],
		kPPrimary,ntmother,1,11);
      KeepTrack(ntmother); 
    }
    else { 
      PushTrack(1,ntmother,pdgPushed[ipart],
		pxPushed[ipart],pyPushed[ipart],pzPushed[ipart],ePushed[ipart],
		origin0[0],origin0[1],origin0[2],time0,
		polar[0],polar[1],polar[2],
		kPDecay,ntchild,1,1);
      KeepTrack(ntchild); 
    }
  }
  SetHighWaterMark(ntchild); 
  AliGenEventHeader* header = new AliGenEventHeader("LMR");
  header->SetPrimaryVertex(fVertex);
  header->SetInteractionTime(fTime);
  header->SetNProduced(fNprimaries);
  AddHeader(header); 
}

//------------------------------------------------------------------

void AliGenMUONLMR::Decay2Body(const TParticle *mother){ 
  // performs decay in two muons of the low mass resonances
  Double_t md1 = fMu[0]->GetMass(); 
  Int_t pdg = mother->GetPdgCode(); 
  Double_t mres =0; 
  // if mother is a rho, extract the mass from its line shape
  // otherwise consider the resonance mass 
  if (pdg == 113) mres = fRhoLineShape->GetRandom(); 
  else mres = mother->GetCalcMass(); 
  //  while (mres < md1 + md2) mres =  fDsigmaDm[res]->GetRandom();
  // energies and momenta in rest frame 
  Double_t e1 = mres / 2.;
  Double_t p1 = TMath::Sqrt((e1 + md1)*(e1 - md1)); 
  // orientation in decaying particle rest frame
  Double_t costheta = gRandom->Rndm() * 2 - 1;
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2. * TMath::Pi() * gRandom->Rndm(); 
  Double_t px1      = p1 * sintheta * TMath::Cos(phi); 
  Double_t py1      = p1 * sintheta * TMath::Sin(phi); 
  Double_t pz1      = p1 * costheta; 

  // boost muons into lab frame 

  TLorentzVector vmother, v1, v2;
  //  TLorentzVector boosted1, boosted2;   
  vmother.SetPxPyPzE(mother->Px(),mother->Py(),mother->Pz(),mother->Energy());
  v1.SetPxPyPzE(px1,py1,pz1,e1); 
  v2.SetPxPyPzE(-px1,-py1,-pz1,e1); 

  TVector3 betaParent = (1./vmother.E())*vmother.Vect(); // beta = p/E
  v1.Boost(betaParent);
  v2.Boost(betaParent);

//   TLorentzVector vtot = v1 + v2; 
//   printf ("mother: %g   %g    %g     %g\n",vmother.Px(), vmother.Py(), vmother.Pz(), vmother.E());
//   printf ("vtot  : %g   %g    %g     %g\n",vtot.Px(), vtot.Py(), vtot.Pz(), vtot.E());

  fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  fMu[1]->SetMomentum(v2.Px(),v2.Py(),v2.Pz(),v2.E());
} 

//------------------------------------------------------------------

void AliGenMUONLMR::DecayPiK(TParticle *mother, Bool_t &hasDecayed){ 
  // performs decays of pions and kaons
  Double_t md1 = fMu[0]->GetMass(); 
  // extract the mass from the resonance's line shape
  Double_t mres = mother->GetMass(); 
  // choose the pi/k sign, assuming 50% probabilities for both signs
  Int_t sign = (gRandom->Rndm() > 0.5) ? 1 : -1;
  mother->SetPdgCode(sign * TMath::Abs(mother->GetPdgCode())); 

  // energies and momenta in rest frame 
  Double_t e1 = (mres*mres + md1*md1)/(2*mres);
  Double_t p1 = TMath::Sqrt((e1 + md1)*(e1 - md1)); 
  // orientation in decaying particle rest frame
  Double_t costheta = gRandom->Rndm() * 2 - 1;
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2. * TMath::Pi() * gRandom->Rndm(); 
  Double_t px1      = p1 * sintheta * TMath::Cos(phi); 
  Double_t py1      = p1 * sintheta * TMath::Sin(phi); 
  Double_t pz1      = p1 * costheta;  

  // boost muons into lab frame 
  TLorentzVector vmother, v1;
  vmother.SetPxPyPzE(mother->Px(),mother->Py(),mother->Pz(),mother->Energy());
  v1.SetPxPyPzE(px1,py1,pz1,e1); 

  TVector3 betaParent = (1./vmother.E())*vmother.Vect(); // beta = p/E
  v1.Boost(betaParent);  
  if (mother->GetPdgCode()>0) fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  else fMu[1]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());

  Int_t idmother = 0; 
  if (TMath::Abs(mother->GetPdgCode())== 211) {
      idmother = 0; 
  } else if (TMath::Abs(mother->GetPdgCode())== 321) {
      idmother = 1; 
  } else {
      AliWarning("Mother particle is not a pion or kaon \n");
  }
  
  Double_t gammaRes = mother->Energy()/mres;
  Double_t zResCM = fDecay[idmother]->GetRandom();
  Double_t zResLab = gammaRes*zResCM;  
  if(zResLab > 0.938) hasDecayed = 0; // 0.938: distance from IP to absorber + lambda_i
  else hasDecayed = 1;

} 

//-------------------------------------------------------------------

void AliGenMUONLMR::DalitzDecay(const TParticle *mother){
  //
  // perform dalitz decays of eta, omega and etaprime 
  //
  //in the rest frame of the virtual photon:
  Double_t mres = mother->GetCalcMass(); 
  Double_t mumass  = fMu[0]->GetMass(); 
  Double_t md3  = 0;  // unless differently specified, third particle is a photon 
  if (mother->GetPdgCode() == 223) md3 = 0.134977; // if mother is an omega, third particle is a pi0

  Int_t index = 0; 
  if (TMath::Abs(mother->GetPdgCode())== 221) {
    // eta
      index = 0; 
  } else if (TMath::Abs(mother->GetPdgCode())== 223) {
    // omega
      index = 1; 
  } else if (mother->GetPdgCode() == 331) {
    // eta'
      index = 2; 
  } else {
      AliWarning("Mother particle is not a eta, eta' or omega \n");
  }

  Int_t flag = 0; 
  Double_t xd=0, mvirt2=0; 
  Double_t countIt = 0;
  while (flag==0) {  
    xd       = fDalitz[index]->GetRandom(); 
    mvirt2   = xd * mres * mres;   // mass of virtual photon 
    // check kinematics 
    if (mres - md3 > TMath::Sqrt(mvirt2) && TMath::Sqrt(mvirt2)/2. > mumass) flag=1;
    if (++countIt>1E11) {
      mvirt2 =  mres * mres * 0.998; 
      break;
    }
  }  

  //
  //        Generate muons in virtual photon rest frame. 
  //        z axis is the virt. photon direction (before boost)  
  //

  Double_t e1 = TMath::Sqrt(mvirt2)/2.; // energy of mu1 in the virtual photon frame
  Double_t psquare = (e1 + mumass)*(e1 - mumass); 
  if (psquare<0) {
    printf("Error in AliGenMUONLMR::DalitzDecay: sqrt of psquare = %f put to 0\n",psquare); 
    psquare = 0;
  }
  Double_t p1 = TMath::Sqrt(psquare);
  //theta angle between the pos. muon and the virtual photon 
  Double_t costheta = fCosTheta->GetRandom();
  if (costheta>1)  costheta = 1; 
  if (costheta<-1) costheta = -1; 
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2 * TMath::Pi() * gRandom->Rndm();
  Double_t sinphi   = TMath::Sin(phi);
  Double_t cosphi   = TMath::Cos(phi);

  // fill 4-vectors of leptons in the virtual photon frame

  Double_t px1 = p1*sintheta*cosphi; 
  Double_t py1 = p1*sintheta*sinphi; 
  Double_t pz1 = p1*costheta; 
  Double_t px2 = -p1*sintheta*cosphi; 
  Double_t py2 = -p1*sintheta*sinphi; 
  Double_t pz2 = -p1*costheta; 
  Double_t e2  = e1; 

  fMu[0]->SetMomentum(px1,py1,pz1,e1); 
  fMu[1]->SetMomentum(px2,py2,pz2,e2); 

  // calculate components of non-dilepton in CMS of parent resonance 

  Double_t e3 = (mres * mres + md3 * md3 - mvirt2) / (2.*mres);
  Double_t psquare3 = (e3 + md3)*(e3 - md3); 
  if (psquare3<0) {
    printf("Error in AliGenMUONLMR::DalitzDecay: sqrt of psquare3 = %f put to 0\n",psquare3); 
    psquare3 = 0;
  }
  Double_t p3 = TMath::Sqrt(psquare3);
  Double_t costheta2 = 2.* gRandom->Rndm() - 1.;   // angle between virtual photon and resonance
  if (costheta2>1)  costheta2 = 1; 
  if (costheta2<-1) costheta2 = -1; 
  Double_t sintheta2 = TMath::Sqrt((1. + costheta2)*(1. - costheta2));
  Double_t phi2      = 2 * TMath::Pi() * gRandom->Rndm();
  Double_t sinphi2   = TMath::Sin(phi2);
  Double_t cosphi2   = TMath::Cos(phi2);
  Double_t px3 = p3*sintheta2*cosphi2; 
  Double_t py3 = p3*sintheta2*sinphi2; 
  Double_t pz3 = p3*costheta2; 
  TLorentzVector v3(px3,py3,pz3,e3); 

  sintheta2 = -sintheta2;
  cosphi2   = -cosphi2;
  sinphi2   = -sinphi2;

  Double_t px1new = px1*costheta2*cosphi2 - py1*sinphi2 + pz1*sintheta2*cosphi2; 
  Double_t py1new = px1*costheta2*sinphi2 + py1*cosphi2 + pz1*sintheta2*sinphi2; 
  Double_t pz1new =-px1*sintheta2                       + pz1*costheta2; 
  Double_t px2new = px2*costheta2*cosphi2 - py2*sinphi2 + pz2*sintheta2*cosphi2; 
  Double_t py2new = px2*costheta2*sinphi2 + py2*cosphi2 + pz2*sintheta2*sinphi2; 
  Double_t pz2new =-px2*sintheta2                       + pz2*costheta2; 

  fMu[0]->SetMomentum(px1new,py1new,pz1new,e1); 
  fMu[1]->SetMomentum(px2new,py2new,pz2new,e2); 

  Double_t evirt = mres - e3; 
  Double_t pxvirt = -px3;
  Double_t pyvirt = -py3;
  Double_t pzvirt = -pz3;
  TLorentzVector vvirt(pxvirt,pyvirt,pzvirt,evirt); 
  TVector3 betaVirt = (1./evirt) * vvirt.Vect(); // virtual photon beta in res frame

  TLorentzVector v1(px1,py1,pz1,e1); 
  TLorentzVector v2(px2,py2,pz2,e2);

  // boost the muons in the frame where the resonance is at rest 

  v1.Boost(betaVirt); 
  v2.Boost(betaVirt); 

  // boost muons and third particle in lab frame

  TLorentzVector vmother(mother->Px(), mother->Py(), mother->Pz(), mother->Energy());  
  TVector3 resBetaLab = (1./vmother.E())*vmother.Vect(); // eta beta in lab frame
  v1.Boost(resBetaLab); 
  v2.Boost(resBetaLab); 
  v3.Boost(resBetaLab); 
  vvirt.Boost(resBetaLab); 

  fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  fMu[1]->SetMomentum(v2.Px(),v2.Py(),v2.Pz(),v2.E());
//   part3->SetMomentum(v3.Px(),v3.Py(),v3.Pz(),v3.E());

//   TLorentzVector vtot = v1 + v2 + v3; 
//   TLorentzVector vdimu = v1 + v2; 
//   printf ("mother: %g   %g    %g     %g\n",vmother.Px(), vmother.Py(), vmother.Pz(), vmother.E());
//   printf ("vtot  : %g   %g    %g     %g\n",vtot.Px(), vtot.Py(), vtot.Pz(), vtot.E());
//   printf ("vvirt : %g   %g    %g     %g\n",vvirt.Px(), vvirt.Py(), vvirt.Pz(), vvirt.E());
//   printf ("vdimu : %g   %g    %g     %g\n",vdimu.Px(), vdimu.Py(), vdimu.Pz(), vdimu.E());

}

//------------------------------------------------------------------

Double_t AliGenMUONLMR::FormFactor(Double_t q2, Int_t decay){ 
  //  Calculates the form factor for Dalitz decays A->B+l+l
  //  Returns: |F(q^2)|^2
  //
  //  References: L.G. Landsberg, Physics Reports 128 No.6 (1985) 301-376. 
 
  Double_t ff2, mass2;
  Double_t n2, n4, m2; 
  // Lepton-G
  
  Double_t lambda2inv = 0;
  switch (decay) { 
  case 0:   // eta -> mu mu gamma  
  // eta   -> l+ l- gamma: pole approximation
    lambda2inv = 1.95; 
    mass2 = fParticle[kEtaLMR]->GetMass() * fParticle[kEtaLMR]->GetMass(); 
    if (q2 < mass2) ff2 = TMath::Power(1./(1.-lambda2inv*q2),2);
    else ff2 = 0; 
    break;
  case 1:   // omega -> mu mu pi0 
    // omega -> l+ l- pi0: pole approximation
    mass2 = fParticle[kOmegaLMR]->GetMass() * fParticle[kOmegaLMR]->GetMass(); 
    lambda2inv = 2.26; 
    if (q2 < mass2) ff2 = TMath::Power(1./(1.-lambda2inv*q2),2);
    else ff2 = 0; 
    break;
  case 2:   // etaPrime -> mu mu gamma 
    mass2 = fParticle[kEtaPrimeLMR]->GetMass() * fParticle[kEtaPrimeLMR]->GetMass(); 
    // eta'  -> l+ l- gamma: Breit-Wigner fitted to data
    n2 = 0.764 * 0.764; 
    n4 = n2 * n2; 
    m2 = 0.1020 * 0.1020;
    if (q2 < mass2) ff2 = n4 / (TMath::Power(n2-q2,2) + m2 * n2); 
    else ff2 = 0; 
    break;
  default:
    printf ("FormFactor: Decay not found\n"); 
    return 0; 
    break; 
  }
  return ff2; 
}

//____________________________________________________________

Double_t AliGenMUONLMR::RhoLineShapeNew(const Double_t *x, const Double_t* /*para*/){
  //new parameterization implemented by Hiroyuki Sako (GSI)
  Double_t mass = *x;
  double r, GammaTot;
  Double_t mRho    = TDatabasePDG::Instance()->GetParticle("rho0")->Mass();
  Double_t mPi     = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();
  Double_t mMu     = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t Gamma0  = TDatabasePDG::Instance()->GetParticle("rho0")->Width();

  const double Norm = 0.0744416*1.01;  
  // 0.0744416 at m = 0.72297
  // is the max number with Norm=1 (for rho)
  
  double mThreshold = 2.*mPi;

  const double T = 0.170; // Assumption of pi+ temperature [GeV/c^2]
  //const double T = 0.11; // Taken from fit to pi+ temperature [GeV/c^2]
  // with Reference: LEBC-EHS collab., Z. Phys. C 50 (1991) 405

  if (mass < mThreshold) {
    r = 0.;
    return r;
  }

  double k = sqrt(0.25*mass*mass-(mThreshold/2)*(mThreshold/2));
  double k0 = sqrt(0.25*mRho*mRho-(mThreshold/2)*(mThreshold/2));

  GammaTot = (k/k0)*(k/k0)*(k/k0)*(mRho/mass)*(mRho/mass)*Gamma0;

  double FormFactor2 = 1/((mass*mass-mRho*mRho)*(mass*mass-mRho*mRho)+
			  mass*mass*GammaTot*GammaTot);

  r = pow(mass,1.5)*pow((1-mThreshold*mThreshold/(mass*mass)),1.5)*
    ((mass*mass+2*mMu*mMu)/(mass*mass))*(pow((mass*mass-4*mMu*mMu),0.5)/mass)*FormFactor2
    *exp(-mass/T)/Norm;

  return r;
}


