/******************************************************************************
 *                      T H E R M I N A T O R                                 *
 *                   THERMal heavy-IoN generATOR                              *
 *                           version 1.0                                      *
 *                                                                            *
 * Authors of the model: Wojciech Broniowski, Wojciech.Broniowski@ifj.edu.pl, *
 *                       Wojciech Florkowski, Wojciech.Florkowski@ifj.edu.pl  *
 * Authors of the code:  Adam Kisiel, kisiel@if.pw.edu.pl                     *
 *                       Tomasz Taluc, ttaluc@if.pw.edu.pl                    *
 * Code designers: Adam Kisiel, Tomasz Taluc, Wojciech Broniowski,            *
 *                 Wojciech Florkowski                                        *
 *                                                                            *
 * For the detailed description of the program and furhter references         * 
 * to the description of the model plesase refer to: nucl-th/0504047,         *
 * accessibile at: http://www.arxiv.org/nucl-th/0504047                       *
 *                                                                            *
 * Homepage: http://hirg.if.pw.edu.pl/en/therminator/                         *
 *                                                                            *
 * This code can be freely used and redistributed. However if you decide to   *
 * make modifications to the code, please contact the authors, especially     *
 * if you plan to publish the results obtained with such modified code.       *
 * Any publication of results obtained using this code must include the       *
 * reference to nucl-th/0504047 and the published version of it, when         *
 * available.                                                                 *
 *                                                                            *
 *****************************************************************************/
#include "THGlobal.h"
#include "ParticleDecayer.h"
#include "DecayChannel.h"
#include "DecayTable.h"
#include <TMath.h>
#include <TDatime.h>

inline Double_t ParticleDecayer::BreitWigner(Double_t Mass, Double_t Gamma) const{
  Double_t x,y;

  y=mRandom->Rndm();
  x=Mass+Gamma/2*TMath::Tan(TMath::Pi()*(y-0.5));

  return x;
}

ParticleDecayer::ParticleDecayer(ParticleDB *aDB)
{
  TDatime dat;

  mDB = aDB;
  mRandom = new TRandom2();
  mRandom->SetSeed(dat.Get() / 2 * 3);
}

ParticleDecayer::~ParticleDecayer()
{
  delete mRandom;
}

void
ParticleDecayer::DecayParticle(Particle *aFather, Particle** aDaughter1, Particle** aDaughter2, Particle** aDaughter3)
{
  ParticleType *tType = aFather->GetParticleType();
#ifdef _RESCALE_CHANNELS_
  int tChannelIndex = (aFather->GetParticleType())->GetTable()->ChooseDecayChannel(mRandom->Rndm());
#else
  Double_t tProb = mRandom->Rndm();
  int tChannelIndex = (aFather->GetParticleType())->GetTable()->ChooseDecayChannelOrNot(tProb);
  if (tChannelIndex == -1) {
    (*aDaughter1) =  NULL;
    (*aDaughter2) =  NULL;
    (*aDaughter3) =  NULL; 
    
    DecayTable *tab = (aFather->GetParticleType())->GetTable();
    PRINT_DEBUG_3("Not decaying " << (aFather->GetParticleType())->GetName() << " for prob " << tProb);
    for (int tIter=0; tIter<=tab->GetChannelCount(); tIter++) {
      PRINT_DEBUG_3(mDB->GetParticleType(tab->GetDecayChannel(tIter)->GetParticle1())->GetName() << " ");
      PRINT_DEBUG_3(mDB->GetParticleType(tab->GetDecayChannel(tIter)->GetParticle2())->GetName() << " ");
      if (tab->GetDecayChannel(tIter)->GetParticle3()>-1)
	PRINT_DEBUG_3(mDB->GetParticleType(tab->GetDecayChannel(tIter)->GetParticle3())->GetName() << " ");
      PRINT_DEBUG_3(tab->GetDecayChannel(tIter)->GetBranchingRatio() << endl);
    }
    
    return;
  }

#endif
  const DecayChannel *tChan = (aFather->GetParticleType())->GetTable()->GetDecayChannel(tChannelIndex);

  if (tChan->Is3Particle())
    {
      // It is a 3-body decay channel
      ParticleType *tType1 = mDB->GetParticleType(tChan->GetParticle1());
      ParticleType *tType2 = mDB->GetParticleType(tChan->GetParticle2());
      ParticleType *tType3 = mDB->GetParticleType(tChan->GetParticle3());
      
      double tE = aFather->GetEnergy();
      double tM = aFather->GetMass();		/*MCH uncommented*/
      //      double tFatherMass=aFather->GetMass();
      //      double tFatherGamma=tType->GetGamma();
      double tM1 = tType1->GetMass();
      double tM2 = tType2->GetMass();
      double tM3 = tType3->GetMass();
/* MCH commented begin
      double tM;
      do {
	tM=BreitWigner(tFatherMass,tFatherGamma);
      }
      while (tM1+tM2+tM3>tM);
MCH commented end*/
      double tES1, tES2, tP1, tP2, tCos12, tZ;
      
      do 
	{
	  // Generate E1 and E2 with the Momnte-Carlo method
	  do 
	    {
	      tES1 = mRandom->Rndm() * (tM - tM2 - tM3 - tM1) + tM1;
	      tES2 = mRandom->Rndm() * (tM - tM1 - tM3 - tM2) + tM2;
	    }
	  while (tES1+tES2 > tM); // The sum of both energies must be smaller than the resonance mass
	  
	  tP1  = TMath::Sqrt(tES1*tES1 - tM1*tM1);
	  tP2  = TMath::Sqrt(tES2*tES2 - tM2*tM2);
	  
	  tZ = tM - tES1 - tES2;
	  tZ *= tZ;
	  tCos12 = (tZ - tP1*tP1 - tP2*tP2 - tM3*tM3)/(2*tP1*tP2);
	}
      while ((tCos12 < -1.0) || (tCos12 > 1.0)); // Cos Theta must exist (be within -1.0 to 1.0 )
      
      double tTime;
      if (tType->GetGamma() == 0.0)
	tTime = 1.0e10;
      else {
	double tTau0 = tE/(tType->GetMass()*tType->GetGamma());
	// When it decays
	tTime = -tTau0*TMath::Log(mRandom->Rndm());
      }
      
      // Decay coordinates
      double rxr = aFather->rx + (aFather->px/tE)*tTime;
      double ryr = aFather->ry + (aFather->py/tE)*tTime;
      double rzr = aFather->rz + (aFather->pz/tE)*tTime;
      double rtr = aFather->rt + tTime;
      
      double tPxr2 = tP2 * TMath::Sqrt(1-tCos12*tCos12);
      double tPzr2 = tP2*tCos12;
      double tPxr3 = - tPxr2;
      double tPzr3 = - (tP1 + tPzr2);
      double tP3 = TMath::Hypot(tPxr3, tPzr3);
      double tES3 = TMath::Hypot(tM3, tP3);

      // Generating Euler angles
      double tPhi = mRandom->Rndm() * 2 * TMath::Pi();
      double tKsi = mRandom->Rndm() * 2 * TMath::Pi();
      double tCosTh = mRandom->Rndm() * 2.0 - 1.0;

      double sp = TMath::Sin(tPhi);
      double cp = TMath::Cos(tPhi);
      double sk = TMath::Sin(tKsi);
      double ck = TMath::Cos(tKsi);
      double st = TMath::Sqrt(1.0-tCosTh*tCosTh);
      double ct = tCosTh;

      // Rotating the whole system
      double tPxp1 = - st*ck * tP1;
      double tPyp1 = st*sk * tP1;
      double tPzp1 = ct * tP1;
      
      double tPxp2 = (cp*ct*ck - sp*sk)  * tPxr2 + (-st*ck) * tPzr2;
      double tPyp2 = (-cp*ct*sk - sp*ck) * tPxr2 + (st*sk)  * tPzr2;
      double tPzp2 = cp*st               * tPxr2 + ct       * tPzr2;

      double tPxp3 = (cp*ct*ck - sp*sk)  * tPxr3 + (-st*ck) * tPzr3;
      double tPyp3 = (-cp*ct*sk - sp*ck) * tPxr3 + (st*sk)  * tPzr3;
      double tPzp3 = cp*st               * tPxr3 + ct       * tPzr3;
      
      double tVx = aFather->px/aFather->GetEnergy();
      double tVy = aFather->py/aFather->GetEnergy();
      double tVz = aFather->pz/aFather->GetEnergy();
      
       tES1 = TMath::Sqrt(tM1*tM1+tPxp1*tPxp1+tPyp1*tPyp1+tPzp1*tPzp1);
       tES2 = TMath::Sqrt(tM2*tM2+tPxp2*tPxp2+tPyp2*tPyp2+tPzp2*tPzp2);
       tES3 = TMath::Sqrt(tM3*tM3+tPxp3*tPxp3+tPyp3*tPyp3+tPzp3*tPzp3);

      double tV2 = tVx*tVx + tVy*tVy + tVz*tVz;
      double tGamma = TMath::Power(1-tV2,-0.5);
      
      // Boosting by the parent velocity
      double tVP = tVx*tPxp1 + tVy*tPyp1 + tVz*tPzp1;
      double tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

      double tPx1 =   tPxp1 + (tgvp + tGamma * tES1) * tVx;
      double tPy1 =   tPyp1 + (tgvp + tGamma * tES1) * tVy;
      double tPz1 =   tPzp1 + (tgvp + tGamma * tES1) * tVz;
  
      tVP = tVx*tPxp2 + tVy*tPyp2 + tVz*tPzp2;
      tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

      double tPx2 =   tPxp2 + (tgvp + tGamma * tES2) * tVx;
      double tPy2 =   tPyp2 + (tgvp + tGamma * tES2) * tVy;
      double tPz2 =   tPzp2 + (tgvp + tGamma * tES2) * tVz;
  
      tVP = tVx*tPxp3 + tVy*tPyp3 + tVz*tPzp3;
      tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;

      double tPx3 =   tPxp3 + (tgvp + tGamma * tES3) * tVx;
      double tPy3 =   tPyp3 + (tgvp + tGamma * tES3) * tVy;
      double tPz3 =   tPzp3 + (tgvp + tGamma * tES3) * tVz;
      
      tES1 = TMath::Sqrt(tM1*tM1+tPx1*tPx1+tPy1*tPy1+tPz1*tPz1);
      tES2 = TMath::Sqrt(tM2*tM2+tPx2*tPx2+tPy2*tPy2+tPz2*tPz2);
      tES3 = TMath::Sqrt(tM3*tM3+tPx3*tPx3+tPy3*tPy3+tPz3*tPz3);
      
      (*aDaughter1) = new Particle(tType1, 
				   tPx1, tPy1, tPz1,
				   rxr, ryr, rzr, 
				   rtr);
      (*aDaughter2) = new Particle(tType2, 
				   tPx2, tPy2, tPz2,
				   rxr, ryr, rzr, 
				   rtr);
      (*aDaughter3) = new Particle(tType3, 
				   tPx3, tPy3, tPz3,
				   rxr, ryr, rzr, 
				   rtr);

      aFather->SetDecayed();
    }
  else
    {
      // It is a regular two-body decay channel
      ParticleType *tType1 = mDB->GetParticleType(tChan->GetParticle1());
      ParticleType *tType2 = mDB->GetParticleType(tChan->GetParticle2());
      
      double tE = aFather->GetEnergy();
      double tM = aFather->GetMass();		/*MCH uncommented*/
      //      double tFatherMass=aFather->GetMass();
      //      double tFatherGamma=tType->GetGamma();
      double tM1 = tType1->GetMass();
      double tM2 = tType2->GetMass();
/* MCH commented begin
      double tM;
      do {
	tM=BreitWigner(tFatherMass,tFatherGamma);
      }
      while (tM1+tM2>tM);
MCH commented end*/

      double tTime;
      if (tType->GetGamma() == 0.0)
	tTime = 1.0e10;
      else {
	double tTau0 = tE/(tType->GetMass()*tType->GetGamma());
	// When it decays
	tTime = -tTau0*TMath::Log(mRandom->Rndm());
      }
      
      // Decay coordinates
      double rxr = aFather->rx + (aFather->px/tE)*tTime;
      double ryr = aFather->ry + (aFather->py/tE)*tTime;
      double rzr = aFather->rz + (aFather->pz/tE)*tTime;
      double rtr = aFather->rt + tTime;
      
      // Decay energy
      double tMC1 = (tM*tM - (tM1+tM2)*(tM1+tM2));
      double tMC2 = (tM*tM - (tM1-tM2)*(tM1-tM2));
      double tMom = TMath::Sqrt(tMC1*tMC2)/(2*tM);
      double tPhi = mRandom->Rndm() * 2 * TMath::Pi();
      double tCosTh = mRandom->Rndm()*2.0-1.0;
      
      double tPtr = tMom*TMath::Sqrt(1-tCosTh*tCosTh);
      double tPxr1 = tPtr*TMath::Cos(tPhi);
      double tPyr1 = tPtr*TMath::Sin(tPhi);
      double tPzr1 = tMom*tCosTh;
      
      double tVx = aFather->px/aFather->GetEnergy();
      double tVy = aFather->py/aFather->GetEnergy();
      double tVz = aFather->pz/aFather->GetEnergy();
      
      double tES1 = TMath::Sqrt(tM1*tM1+tPxr1*tPxr1+tPyr1*tPyr1+tPzr1*tPzr1);
      double tES2 = TMath::Sqrt(tM2*tM2+tPxr1*tPxr1+tPyr1*tPyr1+tPzr1*tPzr1);
      
      double tV2 = tVx*tVx + tVy*tVy + tVz*tVz;
      double tGamma = TMath::Power(1-tV2,-0.5);
      double tVP = tVx*tPxr1 + tVy*tPyr1 + tVz*tPzr1;
      double tgvp = (tGamma - 1.0) * (1.0/tV2) * tVP;
      
      double tPx1 =   tPxr1 + (tgvp + tGamma * tES1) * tVx;
      double tPy1 =   tPyr1 + (tgvp + tGamma * tES1) * tVy;
      double tPz1 =   tPzr1 + (tgvp + tGamma * tES1) * tVz;
  
      double tPx2 = - tPxr1 + (-tgvp + tGamma * tES2) * tVx;
      double tPy2 = - tPyr1 + (-tgvp + tGamma * tES2) * tVy;
      double tPz2 = - tPzr1 + (-tgvp + tGamma * tES2) * tVz;

      (*aDaughter1) = new Particle(tType1, 
				   tPx1, tPy1, tPz1,
				   rxr, ryr, rzr, 
				   rtr);
      (*aDaughter2) = new Particle(tType2, 
				   tPx2, tPy2, tPz2,
				   rxr, ryr, rzr, 
				   rtr);
      (*aDaughter3) = 0;

      aFather->SetDecayed();
  
    }
}

void ParticleDecayer::SeedSet(int aSeed)
{
  //  mRandom->SetSeed2(aSeed, aSeed * 11 % 9);
  mRandom->SetSeed(aSeed);
}

