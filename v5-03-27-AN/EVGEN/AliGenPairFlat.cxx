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


// Generator for particle pairs in a preset
// kinematic range 
// ranges can be set for invariant mass, pair pT, pair rapidity
// and pair azimuth

//
// Comments and suggestions: markus.konrad.kohler@cern.ch
//

#include "TPDGCode.h"
#include <TDatabasePDG.h>
#include "AliGenEventHeader.h"
#include "TF1.h"

#include "AliGenPairFlat.h"

ClassImp(AliGenPairFlat)

//_____________________________________________________________________________
AliGenPairFlat::AliGenPairFlat()
    :AliGenerator(), 
     fPairNpart(10),
     fPairYMin(-10.),
     fPairYMax(10.),
     fPairPhiMin(0.),
     fPairPhiMax(TMath::TwoPi()),
     fPairPtMin(0.001),
     fPairPtMax(50.),
     fPairMassMin(0.),
     fPairMassMax(10.),
     fLegPdg1(11),
     fLegPdg2(-11),
     fAlpha(0.),
     fDebug(0),
     fPol(0)
{
  //
  // Default constructor
  //
}


//_____________________________________________________________________________
AliGenPairFlat::~AliGenPairFlat()
{
  //
  // Destructor
  //
    delete fPol; 

}

//_____________________________________________________________________________
void AliGenPairFlat::Generate()
{
  //
  // Generate a random pair of particles
  //
  
    Float_t polar[3]= {0,0,0};
    Float_t origin[3];
    Float_t time;

    Float_t p[3];
    Int_t i, j, nt;
    Double_t phi, y, mass, mt, e, pt;
    Float_t weight = 1.;
    Double_t  pt1, pt2, y1, y2;

    TLorentzVector mother, dau1, dau2;
    Float_t random[6];

    fPol = new TF1("fPol","1.+[0]*x*x",-1.,1.);
    fPol->SetParameter(0, fAlpha);


    for (j=0;j<3;j++) origin[j]=fOrigin[j];
    time = fTimeOrigin;
    if(fVertexSmear==kPerEvent) {
	Vertex();
	for (j=0;j<3;j++) origin[j]=fVertex[j];
	time = fTime;
     }

	if(fDebug == 2){
	printf("\n\n------------------GENERATOR SETTINGS------------------\n\n");
	printf("You choosed for the mother the Mass range %f - %f     \n",fPairMassMin,fPairMassMax);
	printf("You choosed for the mother the transverse Momentum range %f - %f \n",fPairPtMin,fPairPtMax);
	printf("You choosed for the mother the Phi range %f - %f \n",fPairPhiMin,fPairPhiMax);
	printf("The Particle will decay in (pdg) %i and %i\n \n",fLegPdg1, fLegPdg2);
	printf("The rest Mass of first daughter (%s) is %f \n",
	TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg1))->GetTitle(),
	TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg1))->Mass());
	printf("The rest Mass of second daughter (%s) is %f \n",
	TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg2))->GetTitle(),
	TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg2))->Mass());
	printf("polarization factor is alpha == %lf \n",fAlpha);
	printf("vertex is at x == %f || y == %f || z == %f   \n",origin[0],origin[1],origin[2]);
	printf("\n----------------------------------------------------------\n");
	}

    for(i=0;i<fPairNpart;i++) {

	// mother properties
	Rndm(random,4);
	mass  	= fPairMassMin+random[0]*(fPairMassMax-fPairMassMin);
	pt	= fPairPtMin+random[1]*(fPairPtMax-fPairPtMin);
	y	= fPairYMin+random[2]*(fPairYMax-fPairYMin);
	phi	= fPairPhiMin+random[3]*(fPairPhiMax-fPairPhiMin);

        mt 	= TMath::Sqrt(pt*pt + mass*mass);
	p[0] 	= pt*TMath::Cos(phi);
	p[1] 	= pt*TMath::Sin(phi);
	p[2] 	= mt*TMath::SinH(y);
	e 	= mt*TMath::CosH(y);

	mother.SetPxPyPzE(p[0],p[1],p[2],e); 
	

         if (fDebug == 2)    printf("p = (%+11.4e,%+11.4e,%+11.4e) GeV || pt = %+11.4e || y =  %+11.4e || weight=%+11.4e || E= %+11.4e\n",p[0],p[1],p[2],pt, y,weight,e);

	//decay procedure
	if(!Decay(mother,dau1,dau2,fPol))continue;

	pt1 = dau1.Pt();
	pt2 = dau2.Pt();
	y1 = dau1.Rapidity();
	y2 = dau2.Rapidity();

	//first daughter
	p[0] = dau1.Px();
	p[1] = dau1.Py();
	p[2] = dau1.Pz();
	e = dau1.E();

	if(fVertexSmear==kPerTrack) {
	    Rndm(random,6);
	    for (j=0;j<3;j++) {
		origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }

	    Rndm(random,2);
	    time = fTimeOrigin + fOsigma[2]/TMath::Ccgs()*
	      TMath::Cos(2*random[0]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[1]));
	}

	PushTrack(fTrackIt,-1,fLegPdg1,p,origin,polar,time,kPPrimary,nt, weight);
	if (fDebug == 2) printf("dau1-->id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV || pt = %+11.4e || weight=%+11.4e || E= %+11.4e\n",fLegPdg1,p[0],p[1],p[2],pt1,weight,e);

	//second daughter
	p[0] = dau2.Px();
	p[1] = dau2.Py();
	p[2] = dau2.Pz();
	e = dau2.E();

	if(fVertexSmear==kPerTrack) {
	    Rndm(random,6);
	    for (j=0;j<3;j++) {
		origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }

	    Rndm(random,2);
	    time = fTimeOrigin + fOsigma[2]/TMath::Ccgs()*
	      TMath::Cos(2*random[0]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[1]));
	}

	PushTrack(fTrackIt,-1,fLegPdg2,p,origin,polar,time,kPPrimary,nt,weight);
      if (fDebug == 2) printf("dau2-->id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV || pt = %+11.4e || weight=%+11.4e || E= %+11.4e\n",fLegPdg2,p[0],p[1],p[2],pt2, weight,e);

    }//particle loop

}

//_____________________________________________________________________________
void AliGenPairFlat::Init()
{
	//
	// Init
	//

	printf("AliGenPairFlat::Init() not implemented!!!\n");

}


//_____________________________________________________________________________
Bool_t AliGenPairFlat::Decay(const TLorentzVector& mother, TLorentzVector &dau1, TLorentzVector &dau2, TF1* polfactor)
{
	//
	// decay procedure
	//

	Double_t mp, md1, md2, ed1, ed2, pd1, pd2;
	Double_t costheta, sintheta, phi;

	mp = mother.M();
	md1 = TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg1))->Mass();
	md2 = TDatabasePDG::Instance()->GetParticle(TMath::Abs(fLegPdg2))->Mass();

	if(mp < md1+md2){
	printf("Daughters are heavier than Mother!! Check Kinematics!! \n");
	return 0;
	}

  ed1 = (mp*mp + md1*md1 - md2*md2)/(2.*mp);
  ed2 = mp-ed1;
  pd1 = TMath::Sqrt((ed1+md1)*(ed1-md1));
  pd2 = TMath::Sqrt((mp*mp-(md1+md2)*(md1+md2))*(mp*mp -(md1-md2)*(md1-md2)))/(2.*mp);

  costheta = polfactor->GetRandom(); //polarization
  sintheta = TMath::Sqrt((1.+costheta)*(1.-costheta));
  phi      = 2.0*TMath::Pi()*gRandom->Rndm();

	dau1.SetPxPyPzE(pd1*sintheta*TMath::Cos(phi), 
		        pd1*sintheta*TMath::Sin(phi), 
		        pd1*costheta,
			ed1); 

	dau2.SetPxPyPzE((-1.)*pd1*sintheta*TMath::Cos(phi),
			(-1.)*pd1*sintheta*TMath::Sin(phi), 
			(-1.)*pd1*costheta,
			ed2);

  	TVector3 boost = mother.BoostVector();

  	dau1.Boost(boost);
  	dau2.Boost(boost);

  return 1;

}



