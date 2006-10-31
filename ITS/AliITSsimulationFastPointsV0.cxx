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
 ***********************************************************************/

/* $Id$ */
/////////////////////////////////////////////////////////
//  fast simulation V0                                 //
//                                                     //
//                                                     //
/////////////////////////////////////////////////////////

#include "AliITS.h"
#include "AliITShit.h"
#include "AliITSRecPoint.h"
#include "AliITSmodule.h"
#include "AliITSsimulationFastPointsV0.h"
#include "AliITSstatistics.h"

ClassImp(AliITSsimulationFastPointsV0)

AliITSsimulationFastPointsV0::AliITSsimulationFastPointsV0():
fNrecp(0),
fSx(0),
fSz(0){
  // default constructor
}
AliITSsimulationFastPointsV0::AliITSsimulationFastPointsV0(const char *dataType):
fNrecp(0),
fSx(0),
fSz(0){
  //constructor
  Info("AliITSsimulationFastPointsV0","Standard constructor %s",dataType);
  fSx = new AliITSstatistics(2);
  fSz = new AliITSstatistics(2);
}
//______________________________________________________________________
AliITSsimulationFastPointsV0::AliITSsimulationFastPointsV0(const AliITSsimulationFastPointsV0 & rec):AliITSsimulation(rec),
fNrecp(rec.fNrecp),
fSx(rec.fSx),
fSz(rec.fSz){
    // Copy constructor. 

}
//______________________________________________________________________
AliITSsimulationFastPointsV0& AliITSsimulationFastPointsV0::operator=(const AliITSsimulationFastPointsV0&  /*source*/){
    // Assignment operator
  Error("operator=","Assignment operator not allowed");
  return *this;
}

//----------------------------------------------------------
AliITSsimulationFastPointsV0::~AliITSsimulationFastPointsV0()
{
  //destructor
  if(fSx) delete fSx;
  if(fSz) delete fSz;

}

//-------------------------------------------------------------
void AliITSsimulationFastPointsV0::CreateFastRecPoints(AliITSmodule *mod,Int_t module,TRandom *rndm, TClonesArray* recp) {
  // Fast points simulator for all of the ITS.

  Int_t   nhit,h,trk,ifirst;
  Float_t x,y,z,t,e;// local coordinate (cm) and time of flight, and dedx.
  Float_t x1,y1,z1;
  AliITShit *hit;

  if(rndm!=0) module=0; // fix unsued parameter warnings.
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
	AddSPD(e,mod,trk,recp);
	break;
      case 3: case 4:  // SDDs
	AddSDD(e,mod,trk,recp);
	break;
      case 5: case 6:  // SSDs
	AddSSD(e,mod,trk,recp);
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
void AliITSsimulationFastPointsV0::AddSPD(Float_t &e,
					  AliITSmodule* mod,Int_t trackNumber,TClonesArray* recp){
  //
  TClonesArray &pt=*recp;
  const Float_t kmicronTocm = 1.0e-4;
  //  const Float_t kdEdXtoQ = ;
  const Float_t kRMSx = 12.0*kmicronTocm; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 70.0*kmicronTocm; // microns->cm ITS TDR Table 1.3
  Float_t a1,a2; // general float.
  
  AliITSgeom* gm = mod->GetITS()->GetITSgeom();

  AliITSRecPoint rpSPD(gm);

  rpSPD.SetLabel(trackNumber,0);
  rpSPD.SetLabel(-3,1);
  rpSPD.SetLabel(-3,2);
  rpSPD.SetXZ(mod->GetIndex(),fSx->GetMean(),fSz->GetMean());
  rpSPD.SetdEdX(e);
  rpSPD.SetQ(1.0);
  a1 = fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSPD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSPD.SetSigmaDetLocX2(a1);
  a2 = fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e\n",a2);
  rpSPD.SetSigmaZ2(a2);
  rpSPD.SetDetectorIndex(mod->GetIndex());
  rpSPD.SetLayer(mod->GetLayer());
  //(mod->GetITS())->AddRecPoint(rpSPD);
  new (pt[fNrecp]) AliITSRecPoint(rpSPD);
  fNrecp++;
}
//_______________________________________________________________________
void AliITSsimulationFastPointsV0::AddSDD(Float_t &e,
					  AliITSmodule* mod,Int_t trackNumber,TClonesArray* recp){
  //
  TClonesArray &pt=*recp;
  const Float_t kmicronTocm = 1.0e-4;
  const Float_t kdEdXtoQ = 2.778e+8; 
  const Float_t kRMSx = 38.0*kmicronTocm; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 28.0*kmicronTocm; // microns->cm ITS TDR Table 1.3
  Float_t a1,a2; // general float.
  AliITSgeom* gm = mod->GetITS()->GetITSgeom();

  AliITSRecPoint rpSDD(gm);

  rpSDD.SetLabel(trackNumber,0);
  rpSDD.SetLabel(-3,1);
  rpSDD.SetLabel(-3,2);
  rpSDD.SetXZ(mod->GetIndex(),fSx->GetMean(),fSz->GetMean());
  rpSDD.SetdEdX(e);
  rpSDD.SetQ(kdEdXtoQ*e);
  a1 = fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSDD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSDD.SetSigmaDetLocX2(a1);
  a2 = fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e\n",a2);
  rpSDD.SetSigmaZ2(a2);
  rpSDD.SetDetectorIndex(mod->GetIndex());
  rpSDD.SetLayer(mod->GetLayer());
  new (pt[fNrecp]) AliITSRecPoint(rpSDD);
  fNrecp++;

}
//_______________________________________________________________________
void AliITSsimulationFastPointsV0::AddSSD(Float_t &e,
					  AliITSmodule* mod,Int_t trackNumber,TClonesArray* recp){
  // 
  TClonesArray &pt=*recp;
  const Float_t kmicronTocm = 1.0e-4;
  const Float_t kdEdXtoQ = 2.778e+8;
  const Float_t kRMSx = 20.0*kmicronTocm;  // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 830.0*kmicronTocm; // microns->cm ITS TDR Table 1.3
  Float_t a1,a2; // general float.

  AliITSgeom* gm = mod->GetITS()->GetITSgeom();
 
  AliITSRecPoint rpSSD(gm);

  rpSSD.SetLabel(trackNumber,0);
  rpSSD.SetLabel(-3,1);
  rpSSD.SetLabel(-3,2);
  rpSSD.SetXZ(mod->GetIndex(),fSx->GetMean(),fSz->GetMean());
  rpSSD.SetdEdX(e);
  rpSSD.SetQ(kdEdXtoQ*e);
  a1 = fSx->GetRMS(); a1 *= a1; a1 += kRMSx*kRMSx;
  //  if(a1>1.E5) printf("addSSD: layer=%d track #%d dedx=%e sigmaX2= %e ",
  //		    mod->GetLayer(),trackNumber,e,a1);
  rpSSD.SetSigmaDetLocX2(a1);
  a2 = fSz->GetRMS(); a2 *= a2; a2 += kRMSz*kRMSz;
  //  if(a1>1.E5) printf(" sigmaZ2= %e RMSx=%e RMSz=%e\n",a2,fSx->GetRMS(),fSz->GetRMS());
  rpSSD.SetSigmaZ2(a2);
  rpSSD.SetDetectorIndex(mod->GetIndex());
  rpSSD.SetLayer(mod->GetLayer());
  new (pt[fNrecp]) AliITSRecPoint(rpSSD);
  fNrecp++;


}
//_______________________________________________________________________

