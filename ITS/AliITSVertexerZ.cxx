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
#include "AliITSVertexerZ.h"
#include<TBranch.h>
#include<TClonesArray.h>
#include<TFile.h>
#include<TH1.h>
#include <TString.h>
#include<TTree.h>
#include "AliITSLoader.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSZPoint.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"

/////////////////////////////////////////////////////////////////
// this class implements a fast method to determine
// the Z coordinate of the primary vertex
// for p-p collisions it seems to give comparable or better results
// with respect to what obtained with AliITSVertexerPPZ
// It can be used successfully with Pb-Pb collisions
////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexerZ)



//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ():AliITSVertexer(),
fFirstL1(0),
fLastL1(0),
fFirstL2(0),
fLastL2(0),
fDiffPhiMax(0),
fX0(0.),
fY0(0.),
fZFound(0),
fZsig(0.),
fZCombc(0),
fLowLim(0.),
fHighLim(0.),
fStepCoarse(0),
fTolerance(0.),
fMaxIter(0),
fWindowWidth(0) {
  // Default constructor
  SetDiffPhiMax();
  SetFirstLayerModules();
  SetSecondLayerModules();
  SetLowLimit();
  SetHighLimit();
  SetBinWidthCoarse();
  SetTolerance();
  SetPPsetting();
  ConfigIterations();
  SetWindowWidth();
}

//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ(TString fn, Float_t x0, Float_t y0):AliITSVertexer(fn),
fFirstL1(0),
fLastL1(0),
fFirstL2(0),
fLastL2(0),
fDiffPhiMax(0),
fX0(x0),
fY0(y0),
fZFound(0),
fZsig(0.),
fZCombc(0),
fLowLim(0.),
fHighLim(0.),
fStepCoarse(0),
fTolerance(0.),
fMaxIter(0),
fWindowWidth(0) {
  // Standard Constructor
  SetDiffPhiMax();
  SetFirstLayerModules();
  SetSecondLayerModules();
  SetLowLimit();
  SetHighLimit();
  SetBinWidthCoarse();
  SetTolerance();
  SetPPsetting();
  ConfigIterations();
  SetWindowWidth();

}

//______________________________________________________________________
AliITSVertexerZ::AliITSVertexerZ(const AliITSVertexerZ &vtxr) : AliITSVertexer(vtxr),
fFirstL1(vtxr.fFirstL1),
fLastL1(vtxr.fLastL1),
fFirstL2(vtxr.fFirstL2),
fLastL2(vtxr.fLastL2),
fDiffPhiMax(vtxr.fDiffPhiMax),
fX0(vtxr.fX0),
fY0(vtxr.fY0),
fZFound(vtxr.fZFound),
fZsig(vtxr.fZsig),
fZCombc(vtxr.fZCombc),
fLowLim(vtxr.fLowLim),
fHighLim(vtxr.fHighLim),
fStepCoarse(vtxr.fStepCoarse),
fTolerance(vtxr.fTolerance),
fMaxIter(vtxr.fMaxIter),
fWindowWidth(vtxr.fWindowWidth){
  // Copy constructor

}

//______________________________________________________________________
AliITSVertexerZ& AliITSVertexerZ::operator=(const AliITSVertexerZ&  vtxr ){
  // Assignment operator

  this->~AliITSVertexerZ();
  new(this) AliITSVertexerZ(vtxr);
  return *this;
}

//______________________________________________________________________
AliITSVertexerZ::~AliITSVertexerZ() {
  // Destructor
  delete fZCombc;
}

//______________________________________________________________________
void AliITSVertexerZ::ConfigIterations(Int_t noiter,Float_t *ptr){
  // configure the iterative procedure to gain efficiency for
  // pp events with very low multiplicity
  Float_t defaults[5]={0.05,0.1,0.2,0.3,0.5};
  fMaxIter=noiter;
  if(noiter>5){
    Error("ConfigIterations","Maximum number of iterations is 5\n");
    fMaxIter=5;
  }
  for(Int_t j=0;j<5;j++)fPhiDiffIter[j]=defaults[j];
  if(ptr)for(Int_t j=0;j<fMaxIter;j++)fPhiDiffIter[j]=ptr[j];
}

//______________________________________________________________________
Int_t AliITSVertexerZ::GetPeakRegion(TH1F*h, Int_t &binmin, Int_t &binmax) const {
  // Finds a region around a peak in the Z histogram
  // Case of 2 peaks is treated 
  Int_t imax=h->GetNbinsX();
  Float_t maxval=0;
  Int_t bi1=h->GetMaximumBin();
  Int_t bi2=0;
  for(Int_t i=imax;i>=1;i--){
    if(h->GetBinContent(i)>maxval){
      maxval=h->GetBinContent(i);
      bi2=i;
    }
  }
  Int_t npeaks=0;

  if(bi1==bi2){
    binmin=bi1-3;
    binmax=bi1+3;
    npeaks=1;
  }else{
    TH1F *copy = new TH1F(*h);
    copy->SetBinContent(bi1,0.);
    copy->SetBinContent(bi2,0.);
    Int_t l1=TMath::Max(bi1-3,1);
    Int_t l2=TMath::Min(bi1+3,h->GetNbinsX());
    Float_t cont1=copy->Integral(l1,l2);
    Int_t ll1=TMath::Max(bi2-3,1);
    Int_t ll2=TMath::Min(bi2+3,h->GetNbinsX());
    Float_t cont2=copy->Integral(ll1,ll2);
    if(cont1>cont2){
      binmin=l1;
      binmax=l2;
      npeaks=1;
    }
    if(cont2>cont1){
      binmin=ll1;
      binmax=ll2;
      npeaks=1;
    }
    if(cont1==cont2){
      binmin=l1;
      binmax=ll2;
      if(bi2-bi1==1) npeaks=1;
      else npeaks=2;
    }  
    delete copy;
  }    
  return npeaks;
}
//______________________________________________________________________
AliESDVertex* AliITSVertexerZ::FindVertexForCurrentEvent(Int_t evnumber){
  // Defines the AliESDVertex for the current event
  VertexZFinder(evnumber);
  Int_t ntrackl=0;
  for(Int_t iteraz=0;iteraz<fMaxIter;iteraz++){
    if(fCurrentVertex) ntrackl=fCurrentVertex->GetNContributors();
    if(!fCurrentVertex || ntrackl==0 || ntrackl==-1){
      Float_t diffPhiMaxOrig=fDiffPhiMax;
      fDiffPhiMax=GetPhiMaxIter(iteraz);
      VertexZFinder(evnumber);
      fDiffPhiMax=diffPhiMaxOrig;
    }
  }
  FindMultiplicity(evnumber);
  return fCurrentVertex;
}  




//______________________________________________________________________
void AliITSVertexerZ::VertexZFinder(Int_t evnumber){
  // Defines the AliESDVertex for the current event
  fCurrentVertex = 0;
  AliRunLoader *rl =AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader = (AliITSLoader*)rl->GetLoader("ITSLoader");
  AliITSgeom* geom = itsLoader->GetITSgeom();
  itsLoader->LoadRecPoints();
  rl->GetEvent(evnumber);

  AliITSDetTypeRec detTypeRec;

  TTree *tR = itsLoader->TreeR();
  detTypeRec.SetTreeAddressR(tR);
  TClonesArray *itsRec  = 0;
  // lc and gc are local and global coordinates for layer 1
  Float_t lc[3]; for(Int_t ii=0; ii<3; ii++) lc[ii]=0.;
  Float_t gc[3]; for(Int_t ii=0; ii<3; ii++) gc[ii]=0.;
  // lc2 and gc2 are local and global coordinates for layer 2
  Float_t lc2[3]; for(Int_t ii=0; ii<3; ii++) lc2[ii]=0.;
  Float_t gc2[3]; for(Int_t ii=0; ii<3; ii++) gc2[ii]=0.;

  itsRec = detTypeRec.RecPoints();
  TBranch *branch;
  branch = tR->GetBranch("ITSRecPoints");

  Int_t nrpL1 = 0;
  Int_t nrpL2 = 0;
  // By default fFirstL1=0 and fLastL1=79
  // This loop counts the number of recpoints on layer1 (central modules)
  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    // Keep only central modules
    //    if(module%4==0 || module%4==3)continue;
    //   cout<<"Procesing module "<<module<<" ";
    branch->GetEvent(module);
    //    cout<<"Number of clusters "<<clusters->GetEntries()<<endl;
    nrpL1+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  //By default fFirstL2=80 and fLastL2=239
  //This loop counts the number of RP on layer 2
  for(Int_t module= fFirstL2; module<=fLastL2;module++){
    branch->GetEvent(module);
    nrpL2+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  if(nrpL1 == 0 || nrpL2 == 0){
    ResetHistograms();
    itsLoader->UnloadRecPoints();
    fCurrentVertex = new AliESDVertex(0.,5.3,-2);
    return;
  }
  // The vertex finding is attempted only if the number of RP is !=0 on
  // both layers
  Float_t *xc1 = new Float_t [nrpL1]; // coordinates of the L1 Recpoints
  Float_t *yc1 = new Float_t [nrpL1];
  Float_t *zc1 = new Float_t [nrpL1];
  Float_t *phi1 = new Float_t [nrpL1];
  Float_t *err1 = new Float_t [nrpL1];
  Float_t *xc2 = new Float_t [nrpL2]; // coordinates of the L1 Recpoints
  Float_t *yc2 = new Float_t [nrpL2];
  Float_t *zc2 = new Float_t [nrpL2];
  Float_t *phi2 = new Float_t [nrpL2];
  Float_t *err2 = new Float_t [nrpL2];
  Int_t ind = 0;// running index for RP
  // Force a coarse bin size of 200 microns if the number of clusters on layer 2
  // is low
  if(nrpL2<fPPsetting[0])SetBinWidthCoarse(fPPsetting[1]);
  // By default nbincoarse=(10+10)/0.01=2000
  Int_t nbincoarse = static_cast<Int_t>((fHighLim-fLowLim)/fStepCoarse);
  if(fZCombc)delete fZCombc;
  fZCombc = new TH1F("fZCombc","Z",nbincoarse,fLowLim,fLowLim+nbincoarse*fStepCoarse);

  // Loop on modules of layer 1 

  for(Int_t module= fFirstL1; module<=fLastL1;module++){
    //    if(module%4==0 || module%4==3)continue;
    branch->GetEvent(module);
    Int_t nrecp1 = itsRec->GetEntries();
    for(Int_t j=0;j<nrecp1;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)itsRec->At(j);
      // Local coordinates of this recpoint
      lc[0]=recp->GetDetLocalX();
      lc[2]=recp->GetDetLocalZ();
      geom->LtoG(module,lc,gc);
      // Global coordinates of this recpoints
      gc[0]-=fX0; // Possible beam offset in the bending plane
      gc[1]-=fY0; //   "               "
      xc1[ind]=gc[0];
      yc1[ind]=gc[1];
      zc1[ind]=gc[2];
      // azimuthal angle is computed in the interval 0 --> 2*pi
      phi1[ind] = TMath::ATan2(gc[1],gc[0]);
      if(phi1[ind]<0)phi1[ind]=2*TMath::Pi()+phi1[ind];
      err1[ind]=recp->GetSigmaZ2();
      ind++;
    }
    detTypeRec.ResetRecPoints();
  }
  ind = 0; // the running index is reset for Layer 2
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
      err2[ind]=recp->GetSigmaZ2();
      ind++;
    }
    detTypeRec.ResetRecPoints();
  }
 
/* Test the ffect of mutiple scatternig on error. Negligible
  // Multiple scattering
  Float_t beta=1.,pmed=0.875; //pmed=875 MeV (for tracks with dphi<0.01 rad)
  Float_t beta2=beta*beta;
  Float_t p2=pmed*pmed;
  Float_t rBP=3; //Beam Pipe radius = 3cm
  Float_t dBP=0.08/35.3; // 800 um of Be
  Float_t dL1=0.01; //approx. 1% of radiation length  
  Float_t theta2BP=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(dBP);
  Float_t theta2L1=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(dL1);
*/
  TClonesArray *points = new TClonesArray("AliITSZPoint",nrpL1*nrpL2);
  TClonesArray &pts = *points;
  Int_t nopoints =0;
  for(Int_t i=0;i<nrpL1;i++){ // loop on L1 RP
    Float_t r1=TMath::Sqrt(xc1[i]*xc1[i]+yc1[i]*yc1[i]); // radius L1 RP
    for(Int_t j=0;j<nrpL2;j++){ // loop on L2 RP
      Float_t diff = TMath::Abs(phi2[j]-phi1[i]); // diff in azimuth
      if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff; //diff<pi
      if(diff<fDiffPhiMax){ // cut on 10 milliradians by def.
	Float_t r2=TMath::Sqrt(xc2[j]*xc2[j]+yc2[j]*yc2[j]); // radius L2 RP
//	Float_t tgth=(zc2[j]-zc1[i])/(r2-r1); // slope
	Float_t zr0=(r2*zc1[i]-r1*zc2[j])/(r2-r1); //Z @ null radius
	Float_t ezr0q=(r2*r2*err1[i]+r1*r1*err2[j])/(r2-r1)/(r2-r1); //error on Z @ null radius
	/*
        ezr0q+=r1*r1*(1+tgth*tgth)*theta2L1/2; // multiple scattering in layer 1
        ezr0q+=rBP*rBP*(1+tgth*tgth)*theta2BP/2; // multiple scattering in beam pipe
	*/
	new(pts[nopoints++])AliITSZPoint(zr0,ezr0q);

	fZCombc->Fill(zr0);
      }
    }
  }
  delete [] xc1;
  delete [] yc1;
  delete [] zc1;
  delete [] phi1;
  delete [] err1;
  delete [] xc2;
  delete [] yc2;
  delete [] zc2;
  delete [] phi2;
  delete [] err2;

  points->Sort();

  Double_t contents = fZCombc->GetEntries()- fZCombc->GetBinContent(0)-fZCombc->GetBinContent(nbincoarse+1);
  if(contents<1.){
    //    Warning("FindVertexForCurrentEvent","Insufficient number of rec. points\n");
    ResetHistograms();
    itsLoader->UnloadRecPoints();
    fCurrentVertex = new AliESDVertex(0.,5.3,-1);
    return;
  }

  TH1F *hc = fZCombc;

  
  if(hc->GetBinContent(hc->GetMaximumBin())<3)hc->Rebin(3);
  Int_t binmin,binmax;
  Int_t nPeaks=GetPeakRegion(hc,binmin,binmax);   
  if(nPeaks==2)AliWarning("2 peaks found");
  Float_t zm =0.;
  Float_t ezm =0.;
  Float_t lim1 = hc->GetBinLowEdge(binmin);
  Float_t lim2 = hc->GetBinLowEdge(binmax)+hc->GetBinWidth(binmax);

  if(nPeaks ==1 && (lim2-lim1)<fWindowWidth){
    Float_t c=(lim1+lim2)/2.;
    lim1=c-fWindowWidth/2.;
    lim2=c+fWindowWidth/2.;
  }
  Int_t niter = 0, ncontr=0;
  do {
    // symmetrization
    if(zm  !=0.){
      Float_t semilarg=TMath::Min((lim2-zm),(zm-lim1));
      lim1=zm - semilarg;
      lim2=zm + semilarg;
    }

    zm=0.;
    ezm=0.;
    ncontr=0;
    for(Int_t i =0; i<points->GetEntriesFast(); i++){
      AliITSZPoint* p=(AliITSZPoint*)points->UncheckedAt(i);
      if(p->GetZ()>lim1 && p->GetZ()<lim2){
        Float_t deno = p->GetErrZ();
        zm+=p->GetZ()/deno;
        ezm+=1./deno;
        ncontr++;
      }
    }
    zm/=ezm;
    ezm=TMath::Sqrt(1./ezm);
    niter++;
  } while(niter<10 && TMath::Abs((zm-lim1)-(lim2-zm))>fTolerance);
  fCurrentVertex = new AliESDVertex(zm,ezm,ncontr);
  fCurrentVertex->SetTitle("vertexer: B");
  ResetHistograms();
  itsLoader->UnloadRecPoints();
  return;
}

//_____________________________________________________________________
void AliITSVertexerZ::ResetHistograms(){
  // delete TH1 data members
  if(fZCombc)delete fZCombc;
  fZCombc = 0;
}

//______________________________________________________________________
void AliITSVertexerZ::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsLoader->ReloadRecPoints();
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    //  cout<<"Processing event "<<i<<endl;
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);
    if(fCurrentVertex){
      WriteCurrentVertex();
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
  cout <<"Bin sizes for coarse z histos "<<fStepCoarse<<endl;
  cout <<" Current Z "<<fZFound<<"; Z sig "<<fZsig<<endl;
  cout <<" Debug flag: "<<fDebug<<endl;
  cout <<"First event to be processed "<<fFirstEvent;
  cout <<"\n Last event to be processed "<<fLastEvent<<endl;
  if(fZCombc){
    cout<<"fZCombc exists - entries="<<fZCombc->GetEntries()<<endl;
  }
  else{
    cout<<"fZCombc does not exist\n";
  }
 
  cout <<"=======================================================\n";
}

