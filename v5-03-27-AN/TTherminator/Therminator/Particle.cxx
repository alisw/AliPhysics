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
#include "Particle.h"
#include "ParticleType.h"
#include <TMath.h>
#define FIELDWIDTH 15

Particle::Particle()
{
  mPartType  = 0;
  px = 0; py = 0; pz = 0;
  rx = 0; ry = 0; rz = 0; rt = 0;
  mDecayed = 0;
  mHasFather = -1;
}

Particle::Particle(double aRapidity, double aPt, double aPhip, 
		   double aAlfam, double aRho, double aPhis, double aTau,
		   ParticleType *aType)
{
  mPartType = aType;
  mDecayed = 0;
  mHasFather = -1;
  
  px = aPt*TMath::Cos(aPhip);
  py = aPt*TMath::Sin(aPhip);
  double tMt = TMath::Hypot(GetMass(),aPt);
  pz = tMt*TMath::SinH(aRapidity);
  
  rx = aRho*TMath::Cos(aPhis);
  ry = aRho*TMath::Sin(aPhis);

// New method of calculating rz, rt
  rz = aTau*TMath::SinH(aAlfam);
  rt = aTau*TMath::CosH(aAlfam);
}

Particle::Particle(ParticleType *aType, 
		   double aPx, double aPy, double aPz, 
		   double aRx, double aRy, double aRz,
		   double aTime)
{
  mPartType = aType;
  mDecayed = 0;
  mHasFather = -1;
  
  px = aPx;
  py = aPy;
  pz = aPz;

  rx = aRx;
  ry = aRy;
  rz = aRz;
  rt = aTime;
}


Particle::Particle(const Particle& aParticle)
{
  mPartType = aParticle.GetParticleType();
  mDecayed = 0;
  mHasFather = aParticle.GetFather();
  
  px = aParticle.px;
  py = aParticle.py;
  pz = aParticle.pz;
  
  rx = aParticle.rx;
  ry = aParticle.ry;
  rz = aParticle.rz;
  rt = aParticle.rt;
}

Particle::~Particle()
{
  /* no-op */
}

void   
Particle::WriteParticle(ostream *aOuts)
{
  aOuts->flags(ios::left | ios::scientific);
  aOuts->width(10);
  (*aOuts) << mPartType->GetPDGCode();
  aOuts->width(FIELDWIDTH);
  (*aOuts) << px;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << py;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << pz;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << GetEnergy();
  aOuts->width(FIELDWIDTH);
  (*aOuts) << GetMass();
  aOuts->width(FIELDWIDTH);
  (*aOuts) << rx;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << ry;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << rz;
  aOuts->width(FIELDWIDTH);
  (*aOuts) << rt;
//  (*aOuts) << mPartType->GetCharge() << "   ";
//  if (mHasFather>-1) (*aOuts) << mHasFather << "   ";
  aOuts->width(6);
  (*aOuts) << (mHasFather);
  aOuts->width(2);
  (*aOuts) << (mDecayed);
}


double 
Particle::Pt()
{
  return TMath::Hypot(px, py);
}

double 
Particle::Rapidity()
{
  double tE = TMath::Sqrt(px*px+py*py+pz*pz+GetMass()*GetMass());
  return 0.5*TMath::Log((tE-pz)/(tE+pz));
}

ParticleType* 
Particle::GetParticleType() const
{
  return mPartType;
}

int 
Particle::HadDecayed()
{
  return mDecayed;
}

double 
Particle::GetMass()
{
  return mPartType->GetMass();
}

double        
Particle::GetEnergy()
{
  return TMath::Sqrt(GetMass()*GetMass()+px*px+py*py+pz*pz);
}


double 
Particle::GetI3()
{
  return mPartType->GetI3();
}

double        
Particle::GetBarionN()
{
  return mPartType->GetBarionN();
}

double        
Particle::GetStrangeness()
{
  return mPartType->GetStrangeness();
}

void
Particle::SetDecayed()
{
  mDecayed = 1;
}

void          
Particle::SetFather(int aFather)
{
  mHasFather = aFather;
}

int           
Particle::GetFather() const
{
  return mHasFather;
}

