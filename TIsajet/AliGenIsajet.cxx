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

/* $Id$ */

//
// Wrapper for ISAJET generator.
// It is using ISAJET to comunicate with fortarn code
// 
//	
// 
//

#include <Riostream.h>

//#include "TSystem.h"
//#include "TUnixSystem.h"
#include <TParticle.h>
#include "TIsajet.h"
#include "AliMC.h"
#include "AliIsajetRndm.h"
#include "AliGenIsajet.h"
#include "AliRun.h"



ClassImp(AliGenIsajet)

 AliGenIsajet::AliGenIsajet() :
   AliGenMC(),
   fKineBias(1),
   fTrials(0)
 {  
  
  }
//____________________________________________________________________________
AliGenIsajet::AliGenIsajet(Int_t npart) 
  :AliGenMC(npart),
   fKineBias(1),
   fTrials(0)
{
  
    fIsajet = new TIsajet();
  AliIsajetRndm::SetIsajetRandom(GetRandom());
  
}

AliGenIsajet::AliGenIsajet(const AliGenIsajet &  Isajet) 
  :AliGenMC(Isajet),
   fKineBias(1),
   fTrials(0)
{
  Isajet.Copy(*this);
  
}
//____________________________________________________________________________
AliGenIsajet::~AliGenIsajet() 
{ 
  //
  // Standard destructor
  //
  if (fIsajet) delete fIsajet;
}

//____________________________________________________________________________
void AliGenIsajet::Init() 
{
  //
  // Generator initialisation method
  //
 cout<<"isajet4"<<endl;
   fIsajet->SetECM(2000);
 cout<<"isajet5"<<endl;
  fIsajet->SetJobtype("MINBIAS");
 cout<<"isajet6"<<endl;
 fIsajet->SetIDIN(0,1);
 cout<<"isajet7"<<endl;
  fIsajet->SetIDIN(1,1);
 cout<<"isajet8"<<endl;
   fIsajet->Initialise();
 cout<<"isajet9"<<endl;
}

//____________________________________________________________________________
void AliGenIsajet::Generate() 
{

  Float_t polar[3] =   {0,0,0};
  Float_t origin[3]=   {0,0,0};
  Float_t origin0[3]=  {0,0,0};
  Float_t p[4], random[6];
  Int_t j, kf, ks, imo;
  Int_t nt=0; 
  Int_t jev=0;

  fTrials=0;
  
 static TClonesArray *particles;
 //cout<<"generate"<<endl;
 if(!particles) particles=new TClonesArray("TParticle",10000);
 const Float_t kconv=0.001/2.999792458e8;

  VertexInternal();

  origin[0] = fVertex[0];
  origin[1] = fVertex[1];
  origin[2] = fVertex[2];
  //cout<<"generate1"<<endl;

   while(1)
          {
      //cout<<"generate2"<<endl;
  fIsajet->GenerateEvent();
  //cout<<"generate3"<<endl;
  fTrials++;
  fIsajet->ImportParticles(particles,"All");
  //cout<<"generate4"<<endl;
  Int_t np = particles->GetEntriesFast()-1;
  //cout<<"generate51"<<endl;
   if (np == 0 ) continue;

       Int_t nc=0;
       Int_t * newPos = new Int_t[np];
       //cout<<"generate52"<<endl;

        for (Int_t i = 0; i<np; i++) *(newPos+i)=-1;
	//cout<<"generate5"<<endl;
  for (Int_t i = 0; i<np; i++) 
 {
     TParticle *  iparticle       = (TParticle *) particles->At(i);
     //cout<<"generate6"<<endl;
     imo = iparticle->GetFirstMother();
     kf        = iparticle->GetPdgCode();
     ks        = iparticle->GetStatusCode();
             if (ks != 3 &&
              KinematicSelection(iparticle,0))
          {
              nc++;
     p[0]=iparticle->Px();
     p[1]=iparticle->Py();
     p[2]=iparticle->Pz();
     p[3]=iparticle->Energy();
     Float_t tof = kconv*iparticle->T();
     Int_t   iparent = (imo > -1) ? newPos[imo] : -1;
     Int_t   trackIt = (ks == 1) && fTrackIt;
     PushTrack(trackIt, iparent, kf,
                          p[0], p[1], p[2], p[3],
                          origin[0], origin[1], origin[2],
                          tof,
                          polar[0], polar[1], polar[2],
                          kPPrimary, nt,1.,ks);
      KeepTrack(nt);
        newPos[i]=nt;
          }
  }// end of for: particle loop
      if (newPos) delete[] newPos;
      // MakeHeader();
        if (nc > 0) {
            jev+=nc;
            if (jev >= fNpart || fNpart == -1) {
                fKineBias=Float_t(fNpart)/Float_t(fTrials);
                break;
	    }}
 
  }
 SetHighWaterMark(nt);
}
AliGenIsajet& AliGenIsajet::operator=(const  AliGenIsajet& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return (*this);
}

//____________________________________________________________________________
//____________________________________________________________________________

