/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                          *
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

/*
$Log$
Revision 1.1.2.1  2000/06/11 20:16:05  barbera
New: Fast simulation class for the ITS, class as part of new ITS code
structure.

*/

#include <TParticle.h>
#include "AliITS.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITSstatistics.h"

ClassImp(AliITSsimulationFastPoints)

AliITSsimulationFastPoints::AliITSsimulationFastPoints()
{
  //constructor
  fSx = new AliITSstatistics(2);
  fSz = new AliITSstatistics(2);
}

//----------------------------------------------------------
AliITSsimulationFastPoints::~AliITSsimulationFastPoints()
{
  //destructor
  delete fSx;
  delete fSz;

}

//-------------------------------------------------------------
void AliITSsimulationFastPoints::CreateFastRecPoints(AliITSmodule *mod){
  // Fast points simulator for all of the ITS.
  Int_t   nhit,h,trk,ifirst;
  Float_t x,y,z,t,e;// local coordinate (cm) and time of flight, and dedx.
  Float_t x1,y1,z1;
  AliITShit *hit;

  fSx->Reset(); // Start out with things clearly zeroed
  fSz->Reset(); // Start out with things clearly zeroed
  e = 0.; // Start out with things clearly zeroed
  Double_t weight=1.;
  nhit = mod->GetNhits();
  ifirst = 1;
  for(h=0;h<nhit;h++){
    hit = mod->GetHit(h);
    hit->GetPositionL(x,y,z,t);
    if(ifirst) {x1=x;y1=y;z1=z;}
    e += hit->GetIonization();
    trk = hit->GetTrack();
    fSx->AddValue((Double_t)x,weight);
    fSz->AddValue((Double_t)z,weight);
    ifirst = 0;
    if(hit->StatusExiting()||  // leaving volume
       hit->StatusDisappeared()|| // interacted/decayed...
       hit->StatusStop() // dropped below E cuts.
       ){ // exiting track, write out RecPoint.
      //      if(fSz->GetRMS()>1.E-1) {
      //	TParticle *part = hit->GetParticle();
      //	printf("idpart %d energy %f \n",part->GetPdgCode(),part->Energy());
      //	printf("diffx=%e diffy=%e diffz=%e\n",x-x1,y-y1,z-z1);
      //      }
      switch (mod->GetLayer()){
      case 1: case 2:  // SPDs
	AddSPD(e,mod,trk);
	break;
      case 3: case 4:  // SDDs
	AddSDD(e,mod,trk);
	break;
      case 5: case 6:  // SSDs
	AddSSD(e,mod,trk);
	break;
      } // end switch
      fSx->Reset();
      fSz->Reset();
      e = 0.;
      ifirst = 1;
      continue;
    }// end if
  } // end for h
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::AddSPD(Float_t &e,
					 AliITSmodule *mod,Int_t trackNumber){
  const Float_t kcmTomicron = 1.0e4;
  //  const Float_t kdEdXtoQ = ;
  const Float_t kRMSx = 12.0; // microns ITS TDR Table 1.3
  const Float_t kRMSz = 70.0; // microns ITS TDR Table 1.3
  Float_t a1,a2; // general float.
  AliITSRecPoint rpSPD;
  Int_t *trk = rpSPD.GetTracks();

  trk[0] = trackNumber;
  trk[1] = 0; trk[2] = 0;
  rpSPD.SetX(kcmTomicron*fSx->GetMean());
  rpSPD.SetZ(kcmTomicron*fSz->GetMean());
  rpSPD.SetdEdX(0.0);
  rpSPD.SetQ(1.0);
  a1 = kcmTomicron*fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSPD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSPD.SetSigmaX2(a1);
  a2 = kcmTomicron*fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e\n",a2);
  rpSPD.SetSigmaZ2(a2);
  rpSPD.SetProbability(1.0);

  (mod->GetITS())->AddRecPoint(rpSPD);
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::AddSDD(Float_t &e,
					 AliITSmodule *mod,Int_t trackNumber){

  const Float_t kcmTomicron = 1.0e4;
  const Float_t kdEdXtoQ = 2.778e+8; // Boris Batyuna June 10 2000.
  const Float_t kRMSx = 38.0; // microns ITS TDR Table 1.3
  const Float_t kRMSz = 28.0; // microns ITS TDR Table 1.3
  Float_t a1,a2; // general float.
  AliITSRecPoint rpSDD;
  Int_t *trk = rpSDD.GetTracks();

  trk[0] = trackNumber;
  trk[1] = 0; trk[2] = 0;
  rpSDD.SetX(kcmTomicron*fSx->GetMean());
  rpSDD.SetZ(kcmTomicron*fSz->GetMean());
  rpSDD.SetdEdX(e);
  rpSDD.SetQ(kdEdXtoQ*e);
  a1 = kcmTomicron*fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSDD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSDD.SetSigmaX2(a1);
  a2 = kcmTomicron*fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e\n",a2);
  rpSDD.SetSigmaZ2(a2);
  rpSDD.SetProbability(1.0);

  (mod->GetITS())->AddRecPoint(rpSDD);
}
//_______________________________________________________________________
void AliITSsimulationFastPoints::AddSSD(Float_t &e,
					 AliITSmodule *mod,Int_t trackNumber){

  const Float_t kcmTomicron = 1.0e4;
  const Float_t kdEdXtoQ = 2.778e+8; // Boris Batyuna June 10 2000.
  const Float_t kRMSx = 20.0; // microns ITS TDR Table 1.3
  const Float_t kRMSz = 830.0; // microns ITS TDR Table 1.3
  Float_t a1,a2; // general float.
  AliITSRecPoint rpSSD;
  Int_t *trk = rpSSD.GetTracks();

  trk[0] = trackNumber;
  trk[1] = 0; trk[2] = 0;
  rpSSD.SetX(kcmTomicron*fSx->GetMean());
  rpSSD.SetZ(kcmTomicron*fSz->GetMean());
  rpSSD.SetdEdX(e);
  rpSSD.SetQ(kdEdXtoQ*e);
  a1 = kcmTomicron*fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSSD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSSD.SetSigmaX2(a1);
  a2 = kcmTomicron*fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e RMSx=%e RMSz=%e\n",a2,fSx->GetRMS(),fSz->GetRMS());
  rpSSD.SetSigmaZ2(a2);
  rpSSD.SetProbability(1.0);

  (mod->GetITS())->AddRecPoint(rpSSD);
}
//_______________________________________________________________________
