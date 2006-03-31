/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
#include <AliITSVertexerZ.h>
#include <TString.h>
#include "AliITSLoader.h"
#include<TBranch.h>
#include<TClonesArray.h>
#include<TH1.h>
#include<TTree.h>
#include <AliITSgeom.h>
#include "AliITSDetTypeRec.h"
#include <AliITSRecPoint.h>

/////////////////////////////////////////////////////////////////
// this class implements a fast method to determine
// the Z coordinate of the primary vertex
// for p-p collisions it seems to give comparable or better results
// with respect to what obtained with AliITSVertexerPPZ
// It can be used successfully with Pb-Pb collisions
////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexerZ)



//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ():AliITSVertexer() {
  // Default constructor
  SetDiffPhiMax(0);
  fX0 = 0.;
  fY0 = 0.;
  SetFirstLayerModules(0);
  SetSecondLayerModules(0);
  fZFound = 0;
  fZsig = 0.;
  fZCombc = 0;
  fZCombf = 0;
  SetLowLimit(0.);
  SetHighLimit(0.);
  SetBinWidthCoarse(0.);
  SetBinWidthFine(0.);
  SetTolerance(0.);
}

//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ(TString fn, Float_t x0, Float_t y0):AliITSVertexer(fn) {
  // Standard Constructor
  SetDiffPhiMax();
  fX0 = x0;
  fY0 = y0;
  SetFirstLayerModules();
  SetSecondLayerModules();
  fZFound = 0;
  fZsig = 0.;
  fZCombc = 0;
  fZCombf = 0;
  SetLowLimit();
  SetHighLimit();
  SetBinWidthCoarse();
  SetBinWidthFine();
  SetTolerance();

}

//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ(const AliITSVertexerZ &vtxr) : AliITSVertexer(vtxr) {
  // Copy constructor
  // Copies are not allowed. The method is protected to avoid misuse.
  Error("AliITSVertexerZ","Copy constructor not allowed\n");
}

//______________________________________________________________________
AliITSVertexerZ& AliITSVertexerZ::operator=(const AliITSVertexerZ& /* vtxr */){
  // Assignment operator
  // Assignment is not allowed. The method is protected to avoid misuse.
  Error("= operator","Assignment operator not allowed\n");
  return *this;
}


//______________________________________________________________________
AliITSVertexerZ::~AliITSVertexerZ() {
  // Default Destructor
  //fITS = 0;
  if(fZCombc)delete fZCombc;
  if(fZCombf)delete fZCombf;
}

//______________________________________________________________________
AliESDVertex* AliITSVertexerZ::FindVertexForCurrentEvent(Int_t evnumber){
  // Defines the AliESDVertex for the current event
  
  fCurrentVertex = 0;
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  TDirectory * olddir = gDirectory;
  rl->CdGAFile();
  AliITSgeom* geom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
  olddir->cd();

  itsLoader->LoadRecPoints();
  rl->GetEvent(evnumber);

  AliITSDetTypeRec detTypeRec;

  TTree *tR = itsLoader->TreeR();
  detTypeRec.SetTreeAddressR(tR);
  TClonesArray *itsRec  = 0;
  Float_t lc[3]; for(Int_t ii=0; ii<3; ii++) lc[ii]=0.;
  Float_t gc[3]; for(Int_t ii=0; ii<3; ii++) gc[ii]=0.;
  Float_t lc2[3]; for(Int_t ii=0; ii<3; ii++) lc2[ii]=0.;
  Float_t gc2[3]; for(Int_t ii=0; ii<3; ii++) gc2[ii]=0.;

  itsRec = detTypeRec.RecPoints();
  TBranch *branch;
  branch = tR->GetBranch("ITSRecPoints");

  Int_t nbinfine = static_cast<Int_t>((fHighLim-fLowLim)/fStepFine);
  Int_t nbincoarse = static_cast<Int_t>((fHighLim-fLowLim)/fStepCoarse);
  if(fZCombc)delete fZCombc;
  fZCombc = new TH1F("fZCombc","Z",nbincoarse,fLowLim,fLowLim+nbincoarse*fStepCoarse);
  if(fZCombf)delete fZCombf;
  fZCombf = new TH1F("fZCombf","Z",nbinfine,fLowLim,fLowLim+nbinfine*fStepFine);

  Int_t nrpL1 = 0;
  Int_t nrpL2 = 0;
  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    if(module%4==0 || module%4==3)continue;
    branch->GetEvent(module);
    nrpL1+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  for(Int_t module= fFirstL2; module<=fLastL2;module++){
    branch->GetEvent(module);
    nrpL2+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  if(nrpL1 == 0 || nrpL2 == 0){
    ResetHistograms();
    return fCurrentVertex;
  }
  Float_t *xc1 = new Float_t [nrpL1];
  Float_t *yc1 = new Float_t [nrpL1];
  Float_t *zc1 = new Float_t [nrpL1];
  Float_t *phi1 = new Float_t [nrpL1];
  Float_t *xc2 = new Float_t [nrpL2];
  Float_t *yc2 = new Float_t [nrpL2];
  Float_t *zc2 = new Float_t [nrpL2];
  Float_t *phi2 = new Float_t [nrpL2];
  Int_t ind = 0;
  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    if(module%4==0 || module%4==3)continue;
    branch->GetEvent(module);
    Int_t nrecp1 = itsRec->GetEntries();
    for(Int_t j=0;j<nrecp1;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)itsRec->At(j);
      lc[0]=recp->GetDetLocalX();
      lc[2]=recp->GetDetLocalZ();
      geom->LtoG(module,lc,gc);
      gc[0]-=fX0;
      gc[1]-=fY0;
      xc1[ind]=gc[0];
      yc1[ind]=gc[1];
      zc1[ind]=gc[2];
      phi1[ind] = TMath::ATan2(gc[1],gc[0]);
      if(phi1[ind]<0)phi1[ind]=2*TMath::Pi()+phi1[ind];
      ind++;
    }
    detTypeRec.ResetRecPoints();
  }
  ind = 0;
  for(Int_t module= fFirstL2; module<=fLastL2;module++){
    branch->GetEvent(module);
    Int_t nrecp2 = itsRec->GetEntries();
    for(Int_t j=0;j<nrecp2;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)itsRec->At(j);
      lc[0]=recp->GetDetLocalX();
      lc[2]=recp->GetDetLocalZ();
      geom->LtoG(module,lc,gc);
      gc[0]-=fX0;
      gc[1]-=fY0;
      xc2[ind]=gc[0];
      yc2[ind]=gc[1];
      zc2[ind]=gc[2];
      phi2[ind] = TMath::ATan2(gc[1],gc[0]);
      if(phi2[ind]<0)phi2[ind]=2*TMath::Pi()+phi2[ind];
      ind++;
    }
    detTypeRec.ResetRecPoints();
  }
  for(Int_t i=0;i<nrpL1;i++){
    Float_t r1=TMath::Sqrt(xc1[i]*xc1[i]+yc1[i]*yc1[i]);
    for(Int_t j=0;j<nrpL2;j++){
      Float_t diff = TMath::Abs(phi2[j]-phi1[i]);
      if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff;
      if(diff<fDiffPhiMax){
	Float_t r2=TMath::Sqrt(xc2[j]*xc2[j]+yc2[j]*yc2[j]);
	Float_t zr0=(r2*zc1[i]-r1*zc2[j])/(r2-r1);
	fZCombf->Fill(zr0);
	fZCombc->Fill(zr0);
      }
    }
  }
  delete [] xc1;
  delete [] yc1;
  delete [] zc1;
  delete [] phi1;
  delete [] xc2;
  delete [] yc2;
  delete [] zc2;
  delete [] phi2;
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
  fCurrentVertex->SetTitle("vertexer: B");
  ResetHistograms();
  return fCurrentVertex;
}

//_____________________________________________________________________
void AliITSVertexerZ::ResetHistograms(){
  // delete TH1 data members
  if(fZCombc)delete fZCombc;
  if(fZCombf)delete fZCombf;
  fZCombc = 0;
  fZCombf = 0;
}

//______________________________________________________________________
void AliITSVertexerZ::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsLoader->ReloadRecPoints();
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    cout<<"Processing event "<<i<<endl;
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);
    if(fCurrentVertex){
      WriteCurrentVertex();
    }
    else {
      if(fDebug>0){
	cout<<"Vertex not found for event "<<i<<endl;
	cout<<"fZFound = "<<fZFound<<", fZsig= "<<fZsig<<endl;
      }
    }
  }
}

//________________________________________________________
void AliITSVertexerZ::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";
  cout <<" First layer first and last modules: "<<fFirstL1<<", ";
  cout <<fLastL1<<endl;
  cout <<" Second layer first and last modules: "<<fFirstL2<<", ";
  cout <<fLastL2<<endl;
  cout <<" Max Phi difference: "<<fDiffPhiMax<<endl;
  cout <<"Limits for Z histograms: "<<fLowLim<<"; "<<fHighLim<<endl;
  cout <<"Bin sizes for coarse and fine z histos "<<fStepCoarse<<"; "<<fStepFine<<endl;
  cout <<" Current Z "<<fZFound<<"; Z sig "<<fZsig<<endl;
  cout <<" Debug flag: "<<fDebug<<endl;
  cout <<"First event to be processed "<<fFirstEvent;
  cout <<"\n Last event to be processed "<<fLastEvent<<endl;
 
 
  cout <<"=======================================================\n";
}

