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
#include "AliITSVertexerZD.h"
#include<TBranch.h>
#include<TClonesArray.h>
#include<TH1.h>
#include <TString.h>
#include<TTree.h>
#include "AliESDVertex.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSZPoint.h"
#include <iostream>
#include "AliLog.h"

/////////////////////////////////////////////////////////////////
// this class implements a fast method to determine
// the Z coordinate of the primary vertex
// it is based on the 2 layers of SDD detectors
// wider Z acceptance w.r.t. AliITSVertexerZ
// lower resolution
// origin masera@to.infn.it
////////////////////////////////////////////////////////////////

using std::endl;
using std::cout;
ClassImp(AliITSVertexerZD)



//______________________________________________________________________
AliITSVertexerZD::AliITSVertexerZD():AliITSVertexer(),
fFirstL1(0),
fLastL1(0),
fFirstL2(0),
fLastL2(0),
fDiffPhiMax(0),
fZFound(0),
fZsig(0.),
fZCombc(0),
fLowLim(0.),
fHighLim(0.),
fStepCoarse(0),
fTolerance(0.),
fMaxIter(0),
fWindowWidth(0),
fSearchForPileup(kTRUE)
{
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
AliITSVertexerZD::AliITSVertexerZD(Float_t x0, Float_t y0):AliITSVertexer(),
fFirstL1(0),
fLastL1(0),
fFirstL2(0),
fLastL2(0),
fDiffPhiMax(0),
fZFound(0),
fZsig(0.),
fZCombc(0),
fLowLim(0.),
fHighLim(0.),
fStepCoarse(0),
fTolerance(0.),
fMaxIter(0),
fWindowWidth(0),
fSearchForPileup(kTRUE)
{
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
  SetVtxStart((Double_t)x0,(Double_t)y0,0.);

}

//______________________________________________________________________
AliITSVertexerZD::~AliITSVertexerZD() {
  // Destructor
  delete fZCombc;
}

//______________________________________________________________________
void AliITSVertexerZD::ConfigIterations(Int_t noiter,Float_t *ptr){
  // configure the iterative procedure to gain efficiency for
  // pp events with very low multiplicity
  Float_t defaults[5]={0.02,0.05,0.1,0.2,0.3};
  fMaxIter=noiter;
  if(noiter>5){
    Error("ConfigIterations","Maximum number of iterations is 5\n");
    fMaxIter=5;
  }
  for(Int_t j=0;j<5;j++)fPhiDiffIter[j]=defaults[j];
  if(ptr)for(Int_t j=0;j<fMaxIter;j++)fPhiDiffIter[j]=ptr[j];
}

//______________________________________________________________________
Int_t AliITSVertexerZD::GetPeakRegion(TH1F*h, Int_t &binmin, Int_t &binmax){
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
AliESDVertex* AliITSVertexerZD::FindVertexForCurrentEvent(TTree *itsClusterTree){
  // Defines the AliESDVertex for the current event
  VertexZFinder(itsClusterTree);
  Int_t ntrackl=0;
  for(Int_t iteraz=0;iteraz<fMaxIter;iteraz++){
    if(fCurrentVertex) ntrackl=fCurrentVertex->GetNContributors();
    if(!fCurrentVertex || ntrackl==0 || ntrackl==-1){
      Float_t diffPhiMaxOrig=fDiffPhiMax;
      fDiffPhiMax=GetPhiMaxIter(iteraz);
      VertexZFinder(itsClusterTree);
      fDiffPhiMax=diffPhiMaxOrig;
    }
  }
  if(fComputeMultiplicity) FindMultiplicity(itsClusterTree);
  return fCurrentVertex;
}  

//______________________________________________________________________
void AliITSVertexerZD::VertexZFinder(TTree *itsClusterTree){
  // Defines the AliESDVertex for the current event
  //debug  printf("===========  VertexZFinder has been called \n");
  fCurrentVertex = 0;
  Double_t startPos[3]={GetNominalPos()[0],GetNominalPos()[1],GetNominalPos()[2]};
  Double_t startCov[6]={GetNominalCov()[0],GetNominalCov()[1],GetNominalCov()[2],
			GetNominalCov()[3],GetNominalCov()[4],GetNominalCov()[5]};
  ResetVertex();
  TClonesArray *itsRec  = 0;
  // lc1 and gc1 are local and global coordinates for layer 1
  Float_t gc1[3]={0.,0.,0.}; // ; for(Int_t ii=0; ii<3; ii++) gc1[ii]=0.;
  // lc2 and gc2 are local and global coordinates for layer 2
  Float_t gc2[3]={0.,0.,0.}; //; for(Int_t ii=0; ii<3; ii++) gc2[ii]=0.;
  AliITSRecPointContainer* rpcont=AliITSRecPointContainer::Instance();
  rpcont->FetchClusters(0,itsClusterTree);
  if(!rpcont->IsSDDActive()){
    AliWarning("Null pointer for RecPoints branch, vertex not calculated");
    ResetHistograms();
    fCurrentVertex = new AliESDVertex(startPos,startCov,99999.,-2);
    return;
  }

  Int_t nrpL1 = 0;
  Int_t nrpL2 = 0;
  nrpL1=rpcont->GetNClustersInLayerFast(3);
  nrpL2=rpcont->GetNClustersInLayerFast(4);

  if(nrpL1 == 0 || nrpL2 == 0){
    AliDebug(1,Form("No RecPoints in at least one SDD layer (%d %d)",nrpL1,nrpL2));
    ResetHistograms();
    fCurrentVertex = new AliESDVertex(startPos,startCov,99999.,-2);
    return;
  }
  // Force a coarse bin size of 200 microns if the number of clusters on layer 2
  // is low
  if(nrpL2<fPPsetting[0])SetBinWidthCoarse(fPPsetting[1]);
  // By default nbincoarse=(10+10)/0.01=2000
  Int_t nbincoarse = static_cast<Int_t>((fHighLim-fLowLim)/fStepCoarse);
  if(fZCombc)delete fZCombc;
  fZCombc = new TH1F("fZCombc","Z",nbincoarse,fLowLim,fLowLim+nbincoarse*fStepCoarse);


  Int_t nEntriesMod[260];
  TClonesArray* recpArr[260];
  for(Int_t modul=0; modul<260; ++modul) {
      recpArr[modul]=rpcont->UncheckedGetClusters(modul+240);
      nEntriesMod[modul]=recpArr[modul]->GetEntriesFast();
  }
      
  Int_t maxdim=TMath::Min(nrpL1*nrpL2,50000);  // temporary; to limit the size in PbPb
  static TClonesArray points("AliITSZPoint",maxdim);
  Int_t nopoints =0;
  for(Int_t modul1= fFirstL1; modul1<=fLastL1;modul1++){   // Loop on modules of layer 1
    //    UShort_t ladder=int(modul1/6)+1;  // ladders are numbered starting from 1
    TClonesArray *prpl1=recpArr[modul1-240]; //rpcont->UncheckedGetClusters(modul1);
    Int_t nrecp1 = nEntriesMod[modul1-240]; //prpl1->GetEntries();
    for(Int_t j1=0;j1<nrecp1;j1++){
      AliITSRecPoint *recp1 = (AliITSRecPoint*)prpl1->At(j1);
      recp1->GetGlobalXYZ(gc1);
      gc1[0]-=GetNominalPos()[0]; // Possible beam offset in the bending plane
      gc1[1]-=GetNominalPos()[1]; //   "               "
      Float_t phi1 = TMath::ATan2(gc1[1],gc1[0]);
      if(phi1<0)phi1+=TMath::TwoPi();
      for(Int_t modul2=fFirstL2;modul2<=fLastL2;modul2++){
	  itsRec=recpArr[modul2-240]; // rpcont->UncheckedGetClusters(modul2);
	  Int_t nrecp2 = nEntriesMod[modul2-240]; // itsRec->GetEntries();
	  for(Int_t j2=0;j2<nrecp2;j2++){
	    AliITSRecPoint *recp2 = (AliITSRecPoint*)itsRec->At(j2);
	    recp2->GetGlobalXYZ(gc2);
	    gc2[0]-=GetNominalPos()[0];
	    gc2[1]-=GetNominalPos()[1];
	    Float_t phi2 = TMath::ATan2(gc2[1],gc2[0]);
	    if(phi2<0)phi2+=TMath::TwoPi();

	    Float_t diff = TMath::Abs(phi2-phi1); 
	    if(diff>TMath::Pi())diff=TMath::TwoPi()-diff;
	    if(diff<fDiffPhiMax){
	      //=================   DEBUG
	      /*	      Bool_t goodMatch = kFALSE;
	      for(Int_t kl1=0;kl1<3;kl1++){
		Int_t label1=recp1->GetLabel(kl1);
		if(label1>0){
		  for(Int_t kl2=0;kl2<3;kl2++){
		    if(label1 == recp2->GetLabel(kl2))goodMatch = kTRUE;
		  }
		}
	      }
	      */
	      //=================   END DEBUG
	      Float_t r1=TMath::Sqrt(gc1[0]*gc1[0]+gc1[1]*gc1[1]);
	      Float_t zc1=gc1[2];
	      Float_t erz1=recp1->GetSigmaZ2();
	      Float_t r2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
	      //debug	      printf("r1 =%f ; r2=%f \n",r1,r2);
	      Float_t zc2=gc2[2];
	      Float_t erz2=recp2->GetSigmaZ2();
	      //	Float_t tgth=(zc2[j]-zc1[i])/(r2-r1); // slope (used for multiple scattering)
	      Float_t zr0=(r2*zc1-r1*zc2)/(r2-r1); //Z @ null radius
	      Float_t ezr0q=(r2*r2*erz1+r1*r1*erz2)/((r2-r1)*(r2-r1)); //error on Z @ null radius
	      /*
	      if(goodMatch)printf("Good tracklet: z0=%f, ezr0q=%f\n",zr0,ezr0q);
	      else printf("Fake tracklet: z0=%f, ezr0q=%f\n",zr0,ezr0q);
	      if(goodMatch){
		printf("Labels: ");
		for(Int_t kkk=0;kkk<3;kkk++)printf("%d ",recp1->GetLabel(kkk));
		for(Int_t kkk=0;kkk<3;kkk++)printf("%d ",recp1->GetLabel(kkk));
	      }
	      printf("\n");
	      if(nopoints<maxdim) new(points[nopoints++])AliITSZPoint(zr0,ezr0q,goodMatch);
	      */
	      if(nopoints<maxdim) new(points[nopoints++])AliITSZPoint(zr0,ezr0q);  
	      fZCombc->Fill(zr0);
	    }
	  }
      }
    }
  }
  points.Sort();

  Double_t contents = fZCombc->GetEntries()- fZCombc->GetBinContent(0)-fZCombc->GetBinContent(nbincoarse+1);
  if(contents<1.){
    //    Warning("FindVertexForCurrentEvent","Insufficient number of rec. points\n");
    ResetHistograms();
    fCurrentVertex = new AliESDVertex(startPos,startCov,99999.,-1);
    points.Clear();
    return;
  }

  TH1F *hc = fZCombc;

  
  if(hc->GetBinContent(hc->GetMaximumBin())<3)hc->Rebin(4);
  Int_t binmin,binmax;
  Int_t nPeaks=GetPeakRegion(hc,binmin,binmax);   
  if(nPeaks==2)AliDebug(2,"2 peaks found");
  Float_t zm =0.;
  Float_t ezm =0.;
  Float_t lim1 = hc->GetBinLowEdge(binmin);
  Float_t lim2 = hc->GetBinLowEdge(binmax)+hc->GetBinWidth(binmax);
  Float_t widthSR=lim2-lim1;

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
    for(Int_t i =0; i<points.GetEntriesFast(); i++){
      AliITSZPoint* p=(AliITSZPoint*)points.UncheckedAt(i);
      if(p->GetZ()>lim1 && p->GetZ()<lim2){
        Float_t deno = p->GetErrZ();
        zm+=p->GetZ()/deno;
        ezm+=1./deno;
        ncontr++;
      }
    }
    if(ezm>0) {
      zm/=ezm;
      ezm=TMath::Sqrt(1./ezm);
    }
    niter++;
  } while(niter<10 && TMath::Abs((zm-lim1)-(lim2-zm))>fTolerance);
  if(nPeaks==2) ezm=widthSR;
  Double_t position[3]={GetNominalPos()[0],GetNominalPos()[1],zm};
  Double_t covmatrix[6]={GetNominalCov()[0],0.,GetNominalCov()[2],0.,0.,ezm};
  fCurrentVertex = new AliESDVertex(position,covmatrix,99999.,ncontr);
  fCurrentVertex->SetTitle("vertexer: Z");
  fCurrentVertex->SetDispersion(fDiffPhiMax);
  fNoVertices=1;
  points.Clear();
  if(fSearchForPileup && ncontr>fMinTrackletsForPilup){ 
    Float_t secPeakPos;
    Int_t ncontr2=FindSecondPeak(fZCombc,binmin,binmax,secPeakPos);
    if(ncontr2>=fMinTrackletsForPilup){ 
      fIsPileup=kTRUE;
      fNoVertices=2;
      fZpuv=secPeakPos;
      fNTrpuv=ncontr2;
      AliESDVertex secondVert(secPeakPos,0.1,ncontr2);
      fVertArray = new AliESDVertex[2];
      fVertArray[0]=(*fCurrentVertex);
      fVertArray[1]=secondVert;
    }
  }
  if(fNoVertices==1){
    fVertArray = new AliESDVertex[1];
    fVertArray[0]=(*fCurrentVertex);	  
  }
  
  ResetHistograms();
  return;
}

//_____________________________________________________________________
Int_t AliITSVertexerZD::FindSecondPeak(TH1F* h, Int_t binmin,Int_t binmax, Float_t& secPeakPos){ 
  // Resets bin contents between binmin and binmax and then search 
  // for a second peak position 
  for(Int_t i=binmin-1;i<=binmax+1;i++){
    h->SetBinContent(i,0.);
  }
  Int_t secPeakBin=h->GetMaximumBin();
  secPeakPos=h->GetBinCenter(secPeakBin);
  Int_t secPeakCont=(Int_t)h->GetBinContent(secPeakBin);
  secPeakCont+=(Int_t)h->GetBinContent(secPeakBin-1);
  secPeakCont+=(Int_t)h->GetBinContent(secPeakBin+1);  
  secPeakCont+=(Int_t)h->GetBinContent(secPeakBin-2);
  secPeakCont+=(Int_t)h->GetBinContent(secPeakBin+2);  
  return secPeakCont;
}

//_____________________________________________________________________
void AliITSVertexerZD::ResetHistograms(){
  // delete TH1 data members
  if(fZCombc)delete fZCombc;
  fZCombc = 0;
}

//________________________________________________________
void AliITSVertexerZD::PrintStatus() const {
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
  if(fZCombc){
    cout<<"fZCombc exists - entries="<<fZCombc->GetEntries()<<endl;
  }
  else{
    cout<<"fZCombc does not exist\n";
  }
 
  cout <<"=======================================================\n";
}

