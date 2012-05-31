// $Id$

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
#include "AliHLTITSVertexerZ.h"
#include <TTree.h>
#include<TString.h>
#include<TH1.h>
#include<TMath.h>
#include <AliRun.h>
#include <AliITS.h>
#include "AliITSLoader.h"
#include <AliITSgeom.h>
#include <AliITSgeomTGeo.h>
#include <AliITSRecPoint.h>
#include <AliITSclusterV2.h>

//-------------------------------------------------------------------------
//                Implementation of the HLT ITS vertexer class
//
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

ClassImp(AliHLTITSVertexerZ)

AliHLTITSVertexerZ::AliHLTITSVertexerZ():
  AliITSVertexerZ(),
  fZCombf(0),
  fStepFine(0)
{
  // Constructor in case that there is no runloader
  SetBinWidthFine();
}

AliHLTITSVertexerZ::AliHLTITSVertexerZ(Float_t x0, Float_t y0):
  AliITSVertexerZ(x0,y0),
  fZCombf(0),
  fStepFine(0)
{
  // Standard Constructor
  SetBinWidthFine();
}

AliHLTITSVertexerZ::~AliHLTITSVertexerZ()
{
  // Destructor
  if (fZCombf) delete fZCombf;
}

//______________________________________________________________________
AliESDVertex* AliHLTITSVertexerZ::FindVertexForCurrentEvent(AliITSgeom* /* geom */,TTree *tR){
  // Defines the AliESDVertex for the current event

  fCurrentVertex = 0;

  Double_t lc[3]; for(Int_t ii=0; ii<3; ii++) lc[ii]=0.;
  Double_t gc[3]; for(Int_t ii=0; ii<3; ii++) gc[ii]=0.;
  //Float_t lc2[3]; for(Int_t ii=0; ii<3; ii++) lc2[ii]=0.;
  //Float_t gc2[3]; for(Int_t ii=0; ii<3; ii++) gc2[ii]=0.;

  TClonesArray dummy("AliITSclusterV2",10000), *clusters=&dummy;
  TBranch *branch;
  branch = tR->GetBranch("Clusters");
  branch->SetAddress(&clusters);

  Int_t nrpL1 = 0;
  Int_t nrpL2 = 0;
  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    if(module%4==0 || module%4==3)continue;
    //   cout<<"Procesing module "<<module<<" ";
    tR->GetEvent(module);
    //    cout<<"Number of clusters "<<clusters->GetEntries()<<endl;
    nrpL1+= clusters->GetEntriesFast();
  }
  for(Int_t module= fFirstL2; module<=fLastL2;module++){
    tR->GetEvent(module);
    nrpL2+= clusters->GetEntriesFast();
  }
  if(nrpL1 == 0 || nrpL2 == 0){
    ResetHistograms();
    return fCurrentVertex;
  }

  Int_t nPhiBins = (Int_t)(TMath::Pi()/fDiffPhiMax);
  Float_t phiBinSize = 2*TMath::Pi()/(Float_t)nPhiBins;

  Int_t maxind1 = 2*nrpL1/nPhiBins;
  Float_t **zc1 = new Float_t *[nPhiBins];
  Float_t **phi1 = new Float_t *[nPhiBins];
  Float_t **r1 = new Float_t *[nPhiBins];
  Int_t *ind1 = new Int_t [nPhiBins];
  Int_t maxind2 = 2*nrpL2/nPhiBins;
  Float_t **zc2 = new Float_t *[nPhiBins];
  Float_t **phi2 = new Float_t *[nPhiBins];
  Float_t **r2 = new Float_t *[nPhiBins];
  Int_t *ind2 = new Int_t [nPhiBins];
  for(Int_t i=0;i<nPhiBins;i++) {
    zc1[i] = new Float_t [maxind1];
    phi1[i] = new Float_t [maxind1];
    r1[i] = new Float_t [maxind1];
    zc2[i] = new Float_t [maxind2];
    phi2[i] = new Float_t [maxind2];
    r2[i] = new Float_t [maxind2];
  }
  
  Float_t yshift = 0;
  Float_t zshift[4] = {-10.708000, -3.536000, 3.536000, 10.708000};

  yshift = 0.248499;
  memset(ind1,0,nPhiBins*sizeof(Int_t));
  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    if(module%4==0 || module%4==3)continue;
    tR->GetEvent(module);
    Int_t nrecp1 = clusters->GetEntriesFast();
    for(Int_t j=0;j<nrecp1;j++){
      AliITSclusterV2 *recp = (AliITSclusterV2*)clusters->UncheckedAt(j);
      lc[0]=-recp->GetY()+yshift;
      lc[2]=-recp->GetZ()+zshift[module%4];
      AliITSgeomTGeo::LocalToGlobal(module,lc,gc);
      //      geom->LtoG(module,lc,gc);
      gc[0]-=GetNominalPos()[0];
      gc[1]-=GetNominalPos()[1];
      Float_t xc1,yc1;
      xc1=gc[0];
      yc1=gc[1];
      Float_t phi = TMath::ATan2(gc[1],gc[0]);
      if(phi<0)phi=2*TMath::Pi()+phi;
      Int_t bin = (Int_t)(phi/phiBinSize);
      if(bin>=nPhiBins || bin<0) bin = 0;
      Int_t ind = ind1[bin];
      if(ind<maxind1) {
	phi1[bin][ind] = phi;
	zc1[bin][ind]=gc[2]/fStepFine;
	r1[bin][ind]=sqrt(xc1*xc1+yc1*yc1);
	ind1[bin]++;
      }
    }
    clusters->Delete();
  }
  yshift = 3.096207;
  memset(ind2,0,nPhiBins*sizeof(Int_t));
  for(Int_t module= fFirstL2; module<=fLastL2;module++){
    tR->GetEvent(module);
    Int_t nrecp2 = clusters->GetEntriesFast();
    for(Int_t j=0;j<nrecp2;j++){
      AliITSclusterV2 *recp = (AliITSclusterV2*)clusters->UncheckedAt(j);
      lc[0]=recp->GetY()+yshift;
      lc[2]=-recp->GetZ()+zshift[module%4];
      AliITSgeomTGeo::LocalToGlobal(module,lc,gc);
      //      geom->LtoG(module,lc,gc);
      gc[0]-=GetNominalPos()[0];
      gc[1]-=GetNominalPos()[1];
      Float_t xc2,yc2;
      xc2=gc[0];
      yc2=gc[1];
      Float_t phi = TMath::ATan2(gc[1],gc[0]);
      if(phi<0)phi=2*TMath::Pi()+phi;
      Int_t bin = (Int_t)(phi/phiBinSize+0.5);
      if(bin>=nPhiBins || bin<0) bin = 0;
      Int_t ind = ind2[bin];
      if(ind<maxind2) {
	phi2[bin][ind] = phi;
	zc2[bin][ind]=gc[2]/fStepFine;
	r2[bin][ind]=sqrt(xc2*xc2+yc2*yc2);
	ind2[bin]++;
      }
    }
    clusters->Delete();
  }
  Int_t nbinfine = static_cast<Int_t>((fHighLim-fLowLim)/fStepFine);
  Float_t lowz = fLowLim/fStepFine;
  Int_t *harray = new Int_t[nbinfine];
  memset(harray,0,nbinfine*sizeof(Int_t));
  for(Int_t ibin=0;ibin<nPhiBins;ibin++) {
    Float_t *pphi1 = phi1[ibin];
    Float_t *pr1 = r1[ibin];
    Float_t *pzc1 = zc1[ibin];
    for(Int_t i=0;i<ind1[ibin];i++){
      for(Int_t jbin=ibin;jbin<=(ibin+1);jbin++) {
	Int_t jbin2 = jbin;
	if(jbin==nPhiBins) jbin2 = 0;
	Float_t *pphi2 = phi2[jbin2];
	Float_t *pr2 = r2[jbin2];
	Float_t *pzc2 = zc2[jbin2];
	for(Int_t j=0;j<ind2[jbin2];j++){
	  Float_t diff = TMath::Abs(pphi2[j]-pphi1[i]);
	  if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff;
	  if(diff<fDiffPhiMax){
	    Float_t zr0=(pr2[j]*pzc1[i]-pr1[i]*pzc2[j])/(pr2[j]-pr1[i]);
	    Int_t bin = (Int_t)(zr0-lowz);
	    if(bin>=0 && bin<nbinfine){
	      harray[bin]++;
	    }
	  }
	}
      }
    }
  }
  for(Int_t i=0;i<nPhiBins;i++) {
    delete [] zc1[i];
    delete [] phi1[i];
    delete [] r1[i];
    delete [] zc2[i];
    delete [] phi2[i];
    delete [] r2[i];
  }
  delete [] zc1;
  delete [] phi1;
  delete [] r1;
  delete [] ind1;
  delete [] zc2;
  delete [] phi2;
  delete [] r2;
  delete [] ind2;

  Int_t nbinratio = (Int_t)(fStepCoarse/fStepFine+0.5);
  Int_t nbincoarse = nbinfine/nbinratio;

  if(fZCombc)delete fZCombc;
  fZCombc = new TH1F("fZCombc","Z",nbincoarse,fLowLim,fLowLim+nbincoarse*fStepCoarse);
  if(fZCombf)delete fZCombf;
  fZCombf = new TH1F("fZCombf","Z",nbinfine,fLowLim,fLowLim+nbinfine*fStepFine);
  Stat_t contents=0;
  Int_t counter=0;
  for(Int_t bin=0;bin<nbinfine;bin++) {
    fZCombf->SetBinContent(bin+1,(Stat_t)harray[bin]);
    fZCombf->SetBinError(bin+1,TMath::Sqrt((Stat_t)harray[bin]));
    contents+=(Stat_t)harray[bin];
    counter++;
    if(counter==nbinratio) {
      Int_t binc=bin/nbinratio; 
      fZCombc->SetBinContent(binc+1,contents);
      fZCombc->SetBinError(binc+1,TMath::Sqrt(contents));
      contents=0;
      counter=0;
    }
  }

  delete [] harray;

  if(fZCombc->GetEntries()==0){
    Warning("FindVertexForCurrentEvent","Insufficient number of rec. points\n");
    ResetHistograms();
    return fCurrentVertex;
  }
  //  else {
  //    cout<<"Number of entries in hist. "<<fZCombc->GetEntries()<<endl;
  //  }
  Int_t bi = fZCombc->GetMaximumBin();
  Float_t centre = fZCombc->GetBinCenter(bi);
  Int_t n1 = static_cast<Int_t>((centre-fZCombc->GetBinWidth(bi)-fZCombf->GetBinLowEdge(0))/fZCombf->GetBinWidth(0));
  Int_t n2 = static_cast<Int_t>((centre+fZCombc->GetBinWidth(bi)-fZCombf->GetBinLowEdge(0))/fZCombf->GetBinWidth(0));
  Int_t niter = 0;
  Bool_t goon = kTRUE;
  Int_t num;
  while(goon){
    fZFound = 0.;
    fZsig = 0.;
    num=0;
    for(Int_t n=n1;n<=n2;n++){
      fZFound+=fZCombf->GetBinCenter(n)*fZCombf->GetBinContent(n);
      num+=static_cast<Int_t>(fZCombf->GetBinContent(n));
      fZsig+=fZCombf->GetBinCenter(n)*fZCombf->GetBinCenter(n)*fZCombf->GetBinContent(n);
    }
    if(num<2){
      fZsig = 0.;
    }
    else {
      Float_t radi =  fZsig/(num-1)-fZFound*fZFound/num/(num-1);
      if(radi>0.)fZsig=TMath::Sqrt(radi);
      else fZsig=0.;
      fZFound/=num;
    }
    goon = TMath::Abs(TMath::Abs(fZFound-fZCombf->GetBinCenter(n1))-TMath::Abs(fZFound-fZCombf->GetBinCenter(n2)))>fTolerance;
    n1 = static_cast<Int_t>((fZFound-fZCombc->GetBinWidth(bi)-fZCombf->GetBinLowEdge(0))/fZCombf->GetBinWidth(0));
    n2 = static_cast<Int_t>((fZFound+fZCombc->GetBinWidth(bi)-fZCombf->GetBinLowEdge(0))/fZCombf->GetBinWidth(0));
    niter++;
    if(niter>=10){
      goon = kFALSE;
      Warning("FindVertexForCurrentEvent","The procedure dows not converge\n");
    }
  }
  //  cout<<"Numer of Iterations "<<niter<<endl<<endl;
  fCurrentVertex = new AliESDVertex(fZFound,fZsig,num);
  fCurrentVertex->SetTitle("vertexer: HLT");
  ResetHistograms();
  PrintStatus();
  return fCurrentVertex;
}
