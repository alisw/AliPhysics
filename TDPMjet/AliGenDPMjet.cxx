/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


// Generator using DPMJET as an external generator
// The main DPMJET options are accessable for the user through this interface.
// Uses the TDPMjet implementation of TGenerator.

#include <TDPMjet.h>
#include <TRandom.h>
#include <TArrayI.h>
#include <TParticle.h>
#include <TGraph.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticleClassPDG.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include "AliRunLoader.h"
#include "AliGenDPMjet.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliRun.h"
#include "AliDpmJetRndm.h"
#include "AliHeader.h"
#include "AliStack.h"
#include "AliMC.h"

ClassImp(AliGenDPMjet)

//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet()
    :AliGenMC(), 
     fBeamEn(2750.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(1),
     fDecayAll(0),
     fGenImpPar(0.),
     fProcess(kDpmMb),
     fTriggerMultiplicity(0),
     fTriggerMultiplicityEta(0),
     fTriggerMultiplicityPtMin(0),
     fkTuneForDiff(0),
     fProcDiff(0)
{
// Constructor
    fEnergyCMS = 5500.;
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}


//______________________________________________________________________________
AliGenDPMjet::AliGenDPMjet(Int_t npart)
    :AliGenMC(npart),
     fBeamEn(2750.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(1),
     fDecayAll(0),
     fGenImpPar(0.),
     fProcess(kDpmMb),
     fTriggerMultiplicity(0),
     fTriggerMultiplicityEta(0),
     fTriggerMultiplicityPtMin(0),
     fkTuneForDiff(0),
     fProcDiff(0)
{
// Default PbPb collisions at 5. 5 TeV
//
    fEnergyCMS = 5500.;
    fName = "DPMJET";
    fTitle= "Particle Generator using DPMJET";
    SetTarget();
    SetProjectile();
    fVertex.Set(3);
    AliDpmJetRndm::SetDpmJetRandom(GetRandom());
}

AliGenDPMjet::AliGenDPMjet(const AliGenDPMjet &/*Dpmjet*/)
    :AliGenMC(),
     fBeamEn(2750.),
     fMinImpactParam(0.),
     fMaxImpactParam(5.),
     fICentr(0),
     fSelectAll(0),
     fFlavor(0),
     fTrials(0),
     fSpectators(1),
     fSpecn(0),
     fSpecp(0),
     fDPMjet(0),
     fNoGammas(0),
     fLHC(0),
     fPi0Decay(1),
     fDecayAll(0),
     fGenImpPar(0.),
     fProcess(kDpmMb),
     fTriggerMultiplicity(0),
     fTriggerMultiplicityEta(0),
     fTriggerMultiplicityPtMin(0),
     fkTuneForDiff(0),
     fProcDiff(0)
{
    // Dummy copy constructor
    fEnergyCMS = 5500.;
}

//______________________________________________________________________________
AliGenDPMjet::~AliGenDPMjet()
{
// Destructor
}
//______________________________________________________________________________
void AliGenDPMjet::Init()
{
// Initialization
    
    SetMC(new TDPMjet(fProcess, fAProjectile, fZProjectile, fATarget, fZTarget, 
		      fBeamEn,fEnergyCMS));

    fDPMjet=(TDPMjet*) fMCEvGen;
    //
    // **** Flag to force central production
    // fICentr=1. central production forced 
    // fICentr<0 && fICentr>-100 -> bmin = fMinImpactParam, bmax = fMaxImpactParam	  
    // fICentr<-99 -> fraction of x-sec. = XSFRAC		  
    // fICentr=-1. -> evaporation/fzc suppressed		  
    // fICentr<-1. -> evaporation/fzc suppressed		  
    if (fAProjectile == 1 && TMath::Abs(fZProjectile == 1)) fDPMjet->SetfIdp(1);
    
    fDPMjet->SetfFCentr(fICentr);  
    fDPMjet->SetbRange(fMinImpactParam, fMaxImpactParam); 
    fDPMjet->SetPi0Decay(fPi0Decay);
    fDPMjet->SetDecayAll(fDecayAll);
//
//  Initialize DPMjet  
//    
    fDPMjet->Initialize();
}


//______________________________________________________________________________
void AliGenDPMjet::Generate()
{
// Generate one event

  Double_t polar[3]    =   {0,0,0};
  Double_t origin[3]   =   {0,0,0};
  Double_t p[4]        =   {0};
  Float_t tof;

//  converts from mm/c to s
  const Float_t kconv = 0.001/2.999792458e8;
  Int_t nt  = 0;
  Int_t jev = 0;
  Int_t kf, ks, imo;
  kf = 0;
  fTrials = 0;
  //  Set collision vertex position 
  if (fVertexSmear == kPerEvent) Vertex();
  
  while(1)
  {
//    Generate one event
// --------------------------------------------------------------------------
      fSpecn = 0;  
      fSpecp = 0;
// --------------------------------------------------------------------------
      fDPMjet->GenerateEvent();
      
      fTrials++;

      fDPMjet->ImportParticles(&fParticles,"All");      
      if (fLHC) Boost();

      // Temporaneo
      fGenImpPar = fDPMjet->GetBImpac();

      if(TMath::Abs(fXingAngleY) > 1.e-10) BeamCrossAngle();

      Int_t np = fParticles.GetEntriesFast();
      //
      // Multiplicity Trigger
      if (fTriggerMultiplicity > 0) {
	Int_t multiplicity = 0;
	for (Int_t i = 0; i < np; i++) {
	  TParticle *  iparticle = (TParticle *) fParticles.At(i);
	
	  Int_t statusCode = iparticle->GetStatusCode();
	
	  // Initial state particle
	  if (statusCode != 1)
	    continue;
	  // eta cut
	  if (fTriggerMultiplicityEta > 0 && TMath::Abs(iparticle->Eta()) > fTriggerMultiplicityEta)
	    continue;
	  // pt cut
	  if (iparticle->Pt() < fTriggerMultiplicityPtMin) 
	    continue;
	  
	  TParticlePDG* pdgPart = iparticle->GetPDG();
	  if (pdgPart && pdgPart->Charge() == 0)
	    continue;
	  ++multiplicity;
	}
	//
	//
	if (multiplicity < fTriggerMultiplicity) continue;
	Printf("Triggered on event with multiplicity of %d >= %d", multiplicity, fTriggerMultiplicity);
      }    


      if(fkTuneForDiff && (TMath::Abs(fEnergyCMS - 900) < 1)) {
	if(!CheckDiffraction() ) continue;
      }


      Int_t nc = 0;
      if (np == 0) continue;

      Int_t i;
      Int_t* newPos     = new Int_t[np];
      Int_t* pSelected  = new Int_t[np];

      for (i = 0; i<np; i++) {
	  newPos[i]    = i;
	  pSelected[i] = 0;
      }
      
//      First select parent particles

      for (i = 0; i<np; i++) {
	  TParticle *iparticle = (TParticle *) fParticles.At(i);

// Is this a parent particle ?

	  if (Stable(iparticle)) continue;

	  Bool_t  selected             =  kTRUE;
	  Bool_t  hasSelectedDaughters =  kFALSE;
	  
	  kf = iparticle->GetPdgCode();
	  if (kf == 92 || kf == 99999) continue;
	  ks = iparticle->GetStatusCode();
// No initial state partons
          if (ks==21) continue;
	    
	  if (!fSelectAll) selected = KinematicSelection(iparticle, 0) && 
			       SelectFlavor(kf);

	  
	  hasSelectedDaughters = DaughtersSelection(iparticle);


// Put particle on the stack if it is either selected or 
// it is the mother of at least one seleted particle

	  if (selected || hasSelectedDaughters) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop parents

// Now select the final state particles


      for (i=0; i<np; i++) {
	  TParticle *iparticle = (TParticle *) fParticles.At(i);

// Is this a final state particle ?

	  if (!Stable(iparticle)) continue;
      
	  Bool_t  selected =  kTRUE;
	  kf = iparticle->GetPdgCode();
	  ks = iparticle->GetStatusCode();

// --------------------------------------------------------------------------
// Count spectator neutrons and protons (ks == 13, 14)
	  if(ks == 13 || ks == 14){
	      if(kf == kNeutron) fSpecn += 1;
	      if(kf == kProton)  fSpecp += 1;
	  }
// --------------------------------------------------------------------------

	  if (!fSelectAll) {
	      selected = KinematicSelection(iparticle,0)&&SelectFlavor(kf);
	      if (!fSpectators && selected) selected = (ks == 13 || ks == 14);
	  }

// Put particle on the stack if selected

	  if (selected) {
	      nc++;
	      pSelected[i] = 1;
	  } // selected
      } // particle loop final state

// Write particles to stack

      for (i = 0; i<np; i++) {
	  TParticle *  iparticle = (TParticle *) fParticles.At(i);
	  Bool_t  hasMother   = (iparticle->GetFirstMother()>=0);
	  if (pSelected[i]) {
	      
	      kf   = iparticle->GetPdgCode();	      
	      ks   = iparticle->GetStatusCode();	      
	      
	      p[0] = iparticle->Px();
	      p[1] = iparticle->Py();
	      p[2] = iparticle->Pz();
	      p[3] = iparticle->Energy();
	      origin[0] = fVertex[0]+iparticle->Vx()/10; // [cm]
	      origin[1] = fVertex[1]+iparticle->Vy()/10; // [cm]
	      origin[2] = fVertex[2]+iparticle->Vz()/10; // [cm]
		    
	      tof = fTime + kconv*iparticle->T();
	      
	      imo = -1;
	      TParticle* mother = 0;
	      if (hasMother) {
		  imo = iparticle->GetFirstMother();
		  mother = (TParticle *) fParticles.At(imo);
		  imo = (mother->GetPdgCode() != 92 && mother->GetPdgCode() != 99999) ? newPos[imo] : -1;
	      } // if has mother   


	      
	      Bool_t tFlag = (fTrackIt && (ks == 1));
	      PushTrack(tFlag, imo, kf, 
			p[0], p[1], p[2], p[3], 
			origin[0], origin[1], origin[2], tof,
			polar[0], polar[1], polar[2],
			kPNoProcess, nt, 1., ks);
	      KeepTrack(nt);
	      newPos[i] = nt;
	  } // if selected
      } // particle loop
      delete[] newPos;
      delete[] pSelected;
      if (nc>0) {
	  jev += nc;
	  if (jev >= fNpart || fNpart == -1) {
	      break;
	  }
      }
  } // event loop
  MakeHeader();
  SetHighWaterMark(nt);
}

//______________________________________________________________________________
Bool_t AliGenDPMjet::DaughtersSelection(TParticle* iparticle)
{
//
// Looks recursively if one of the daughters has been selected
//
//    printf("\n Consider daughters %d:",iparticle->GetPdgCode());
    Int_t imin = -1;
    Int_t imax = -1;
    Int_t i;
    Bool_t hasDaughters = (iparticle->GetFirstDaughter() >=0);
    Bool_t selected = kFALSE;
    if (hasDaughters) {
	imin = iparticle->GetFirstDaughter();
	imax = iparticle->GetLastDaughter();       
	for (i = imin; i <= imax; i++){
	    TParticle *  jparticle = (TParticle *) fParticles.At(i);	
	    Int_t ip = jparticle->GetPdgCode();
	    if (KinematicSelection(jparticle,0)&&SelectFlavor(ip)) {
		selected=kTRUE; break;
	    }
	    if (DaughtersSelection(jparticle)) {selected=kTRUE; break; }
	}
    } else {
	return kFALSE;
    }
    return selected;
}



//______________________________________________________________________________
Bool_t AliGenDPMjet::SelectFlavor(Int_t pid)
{
// Select flavor of particle
// 0: all
// 4: charm and beauty
// 5: beauty
    Bool_t res = 0;
    
    if (fFlavor == 0) {
	res = kTRUE;
    } else {
	Int_t ifl = TMath::Abs(pid/100);
	if (ifl > 10) ifl/=10;
	res = (fFlavor == ifl);
    }
//
//  This part if gamma writing is inhibited
    if (fNoGammas) 
	res = res && (pid != kGamma && pid != kPi0);
//
    return res;
}

//______________________________________________________________________________
Bool_t AliGenDPMjet::Stable(TParticle*  particle)
{
// Return true for a stable particle
//
    
//    if (particle->GetFirstDaughter() < 0 ) return kTRUE;
    if (particle->GetStatusCode() == 1) return kTRUE;
    else return kFALSE;

}

//______________________________________________________________________________
void AliGenDPMjet::MakeHeader()
{
//  printf("MakeHeader %13.3f \n", fDPMjet->GetBImpac());
// Builds the event header, to be called after each event
    AliGenEventHeader* header = new AliGenDPMjetEventHeader("DPMJET");
    ((AliGenDPMjetEventHeader*) header)->SetNProduced(fDPMjet->GetNumStablePc());
    ((AliGenDPMjetEventHeader*) header)->SetImpactParameter(fDPMjet->GetBImpac());
    ((AliGenDPMjetEventHeader*) header)->SetTotalEnergy(fDPMjet->GetTotEnergy());
    ((AliGenDPMjetEventHeader*) header)->SetParticipants(fDPMjet->GetProjParticipants(), 
    							 fDPMjet->GetTargParticipants());

    if(fProcDiff>0){
    ((AliGenDPMjetEventHeader*) header)->SetProcessType(fProcDiff);
    }
    else 
      ((AliGenDPMjetEventHeader*) header)->SetProcessType(fDPMjet->GetProcessCode());

    // Bookkeeping for kinematic bias
    ((AliGenDPMjetEventHeader*) header)->SetTrials(fTrials);
    // Event Vertex
    header->SetPrimaryVertex(fVertex);
    header->SetInteractionTime(fTime);
    gAlice->SetGenEventHeader(header);    
    AddHeader(header);
}

void AliGenDPMjet::AddHeader(AliGenEventHeader* header)
{
    // Add header to container or runloader
    if (fContainer) {
        fContainer->AddHeader(header);
    } else {
        AliRunLoader::Instance()->GetHeader()->SetGenEventHeader(header);
    }
}


//______________________________________________________________________________
AliGenDPMjet& AliGenDPMjet::operator=(const  AliGenDPMjet& /*rhs*/)
{
// Assignment operator
    return *this;
}


void AliGenDPMjet::FinishRun()
{
    // Print run statistics
    fDPMjet->Dt_Dtuout();
}



Bool_t AliGenDPMjet::CheckDiffraction()
{

  //  printf("AAA\n");

   Int_t np = fParticles.GetEntriesFast();

   Int_t iPart1=-1;
   Int_t iPart2=-1;

   Double_t y1 = 1e10;
   Double_t y2 = -1e10;

  const Int_t kNstable=20;
  const Int_t pdgStable[20] = {
    22,             // Photon
    11,             // Electron
    12,             // Electron Neutrino 
    13,             // Muon 
    14,             // Muon Neutrino
    15,             // Tau 
    16,             // Tau Neutrino
    211,            // Pion
    321,            // Kaon
    311,            // K0
    130,            // K0s
    310,            // K0l
    2212,           // Proton 
    2112,           // Neutron
    3122,           // Lambda_0
    3112,           // Sigma Minus
    3222,           // Sigma Plus
    3312,           // Xsi Minus 
    3322,           // Xsi0
    3334            // Omega
  };
    
     for (Int_t i = 0; i < np; i++) {
	TParticle *  part = (TParticle *) fParticles.At(i);
	
	Int_t statusCode = part->GetStatusCode();
	
	// Initial state particle
	if (statusCode != 1)
	  continue;

	Int_t pdg = TMath::Abs(part->GetPdgCode());
	Bool_t isStable = kFALSE;
	for (Int_t i1 = 0; i1 < kNstable; i1++) {
	  if (pdg == pdgStable[i1]) {
	    isStable = kTRUE;
	    break;
	  }
	}
	if(!isStable) 
	  continue;

	Double_t y = part->Y();

	if (y < y1)
	  {
	    y1 = y;
	    iPart1 = i;
	  }
	if (y > y2)
	{
	  y2 = y;
	  iPart2 = i;
	}
     }

     if(iPart1<0 || iPart2<0) return kFALSE;

     y1=TMath::Abs(y1);
     y2=TMath::Abs(y2);

     TParticle *  part1 = (TParticle *) fParticles.At(iPart1);
     TParticle *  part2 = (TParticle *) fParticles.At(iPart2);

     Int_t pdg1 = part1->GetPdgCode();
     Int_t pdg2 = part2->GetPdgCode();


     Int_t iPart = -1;
     if (pdg1 == 2212 && pdg2 == 2212)
       {
	 if(y1 > y2) 
	   iPart = iPart1;
	 else if(y1 < y2) 
	   iPart = iPart2;
	 else {
	   iPart = iPart1;
	   if((AliDpmJetRndm::GetDpmJetRandom())->Uniform(0.,1.)>0.5) iPart = iPart2;
	 }
       }
     else if (pdg1 == 2212)
       iPart = iPart1;
     else if (pdg2 == 2212)
       iPart = iPart2;





     Double_t M=-1.;
     if(iPart>0) {
       TParticle *  part = (TParticle *) fParticles.At(iPart);
       Double_t E= part->Energy();
       Double_t P= part->P();
       M= TMath::Sqrt((fEnergyCMS-E-P)*(fEnergyCMS-E+P));
     }

const Int_t nbin=120;
Double_t bin[]={
1.080000, 1.274258, 1.468516, 1.662773, 1.857031, 2.051289, 
2.245547, 2.439805, 2.634062, 2.828320, 3.022578, 3.216836, 
3.411094, 3.605352, 3.799609, 3.993867, 4.188125, 4.382383, 
4.576641, 4.770898, 4.965156, 5.547930, 6.130703, 6.713477, 
7.296250, 7.879023, 8.461797, 9.044570, 9.627344, 10.210117, 
10.792891, 11.375664, 11.958437, 12.541211, 13.123984, 13.706758, 
14.289531, 14.872305, 15.455078, 16.037852, 16.620625, 17.203398, 
17.786172, 18.368945, 18.951719, 19.534492, 20.117266, 20.700039, 
21.282812, 21.865586, 22.448359, 23.031133, 23.613906, 24.196680, 
24.779453, 25.362227, 25.945000, 26.527773, 27.110547, 27.693320, 
28.276094, 28.858867, 29.441641, 30.024414, 30.607187, 31.189961, 
31.772734, 32.355508, 32.938281, 33.521055, 34.103828, 34.686602, 
35.269375, 35.852148, 36.434922, 37.017695, 37.600469, 38.183242, 
38.766016, 39.348789, 39.931562, 40.514336, 41.097109, 41.679883, 
42.262656, 42.845430, 43.428203, 44.010977, 44.593750, 45.176523, 
45.759297, 46.342070, 46.924844, 47.507617, 48.090391, 48.673164, 
49.255937, 49.838711, 50.421484, 57.220508, 64.019531, 70.818555, 
77.617578, 84.416602, 91.215625, 98.014648, 104.813672, 111.612695, 
118.411719, 125.210742, 132.009766, 138.808789, 145.607812, 152.406836, 
159.205859, 166.004883, 172.803906, 179.602930, 186.401953, 193.200977, 
200.000000};
Double_t w[]={
1.000000, 0.367136, 0.239268, 0.181139, 0.167470, 0.160072, 
0.147832, 0.162765, 0.176103, 0.156382, 0.146040, 0.143375, 
0.134038, 0.126747, 0.123152, 0.119424, 0.113839, 0.109433, 
0.107180, 0.104690, 0.096427, 0.090603, 0.083706, 0.077206, 
0.074603, 0.069698, 0.067315, 0.064980, 0.063560, 0.059573, 
0.058712, 0.057581, 0.055944, 0.055442, 0.053272, 0.051769, 
0.051672, 0.049284, 0.048980, 0.048797, 0.047434, 0.047039, 
0.046395, 0.046227, 0.044288, 0.044743, 0.043772, 0.043902, 
0.042771, 0.043232, 0.042222, 0.041668, 0.041988, 0.040858, 
0.039672, 0.040069, 0.040274, 0.039438, 0.039903, 0.039083, 
0.038741, 0.038182, 0.037664, 0.038610, 0.038759, 0.038688, 
0.038039, 0.038220, 0.038145, 0.037445, 0.036765, 0.037333, 
0.036753, 0.036405, 0.036339, 0.037659, 0.036139, 0.036706, 
0.035393, 0.037136, 0.036570, 0.035234, 0.036832, 0.035560, 
0.035509, 0.035579, 0.035100, 0.035471, 0.035421, 0.034494, 
0.035596, 0.034935, 0.035810, 0.034324, 0.035355, 0.034323, 
0.033486, 0.034622, 0.034805, 0.034419, 0.033946, 0.033927, 
0.034224, 0.033942, 0.034088, 0.034190, 0.034620, 0.035294, 
0.035650, 0.035378, 0.036028, 0.035933, 0.036753, 0.037171, 
0.037528, 0.037985, 0.039589, 0.039359, 0.040269, 0.040755};

 Double_t wSD=1.;
 Double_t wDD=0.100418;
 Double_t wND=0.050277;

 if(M>-1 && M<bin[0]) return kFALSE;
 if(M>bin[nbin]) M=-1;

 Int_t procType=fDPMjet->GetProcessCode();//fPythia->GetMSTI(1);
 Int_t proc0=2;
 if(procType== 7) proc0=1;
 if(procType== 5 || procType== 6) proc0=0;


 // printf("M = %f   bin[nbin] = %f\n",M, bin[nbin]);

 Int_t proc=2;
 if(M>0) proc=0;
 else if(proc0==1) proc=1;

 if(proc==0 && (AliDpmJetRndm::GetDpmJetRandom())->Uniform(0.,1.) > wSD) return kFALSE;
 if(proc==1 && (AliDpmJetRndm::GetDpmJetRandom())->Uniform(0.,1.) > wDD) return kFALSE;
 if(proc==2 && (AliDpmJetRndm::GetDpmJetRandom())->Uniform(0.,1.) > wND) return kFALSE;


    //     if(proc==1 || proc==2) return kFALSE;

    if(proc!=0) {
      if(proc0!=0) fProcDiff = procType;
      else       fProcDiff = 1;
      return kTRUE;
    }

    Int_t ibin=nbin-1;
    for(Int_t i=1; i<=nbin; i++) 
      if(M<=bin[i]) {
	ibin=i-1;
	//	printf("Mi> %f && Mi< %f\n", bin[i-1], bin[i]);
	break;
      }

    //    printf("w[ibin] = %f\n", w[ibin]);

    if((AliDpmJetRndm::GetDpmJetRandom())->Uniform(0.,1.)> w[ibin]) return kFALSE;

    //    printf("iPart = %d\n", iPart);

    if(iPart==iPart1) fProcDiff=5;
    else if(iPart==iPart2) fProcDiff=6;
    else {
      printf("EROOR:  iPart!=iPart1 && iPart!=iPart2\n");

    }

    return kTRUE;
}


//______________________________________________________________________________
