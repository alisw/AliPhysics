/**************************************************************************
 * Copyright(c) 2006-2008, ALICE Experiment at CERN, All rights reserved. *
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
#include <TTree.h>
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliStrLine.h"
#include "AliTracker.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
#include "AliITSgeomTGeo.h"
#include "AliVertexerTracks.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
/////////////////////////////////////////////////////////////////
// this class implements a method to determine
// the 3 coordinates of the primary vertex
// for p-p collisions 
// It can be used successfully with Pb-Pb collisions
////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexer3D)

/* $Id$ */

//______________________________________________________________________
AliITSVertexer3D::AliITSVertexer3D():AliITSVertexer(),
fLines("AliStrLine",1000),
fVert3D(),
fCoarseDiffPhiCut(0.),
fCoarseMaxRCut(0.),
fMaxRCut(0.),
fZCutDiamond(0.),
fMaxZCut(0.),
fDCAcut(0.),
fDiffPhiMax(0.),
fMeanPSelTrk(0.),
fMeanPtSelTrk(0.),
fUsedCluster(kMaxCluPerMod*kNSPDMod),
fZHisto(0),
fDCAforPileup(0.),
fPileupAlgo(0)
{
  // Default constructor
  SetCoarseDiffPhiCut();
  SetCoarseMaxRCut();
  SetMaxRCut();
  SetZCutDiamond();
  SetMaxZCut();
  SetDCACut();
  SetDiffPhiMax();
  SetMeanPSelTracks();
  SetMeanPtSelTracks();
  SetMinDCAforPileup();
  SetPileupAlgo();
  Float_t binsize=0.02; // default 200 micron
  Int_t nbins=static_cast<Int_t>(1+2*fZCutDiamond/binsize);
  fZHisto=new TH1F("hz","",nbins,-fZCutDiamond,-fZCutDiamond+binsize*nbins);
}

//______________________________________________________________________
AliITSVertexer3D::~AliITSVertexer3D() {
  // Destructor
  fLines.Clear("C");
  if(fZHisto) delete fZHisto;
}

//______________________________________________________________________
void AliITSVertexer3D::ResetVert3D(){
  //
  fVert3D.SetXv(0.);
  fVert3D.SetYv(0.);
  fVert3D.SetZv(0.);
  fVert3D.SetDispersion(0.);
  fVert3D.SetNContributors(0);
  fIsPileup=kFALSE;
  fNTrpuv=-2;
  fZpuv=-99999.;
  fUsedCluster.ResetAllBits(0);
}
//______________________________________________________________________
AliESDVertex* AliITSVertexer3D::FindVertexForCurrentEvent(TTree *itsClusterTree){
  // Defines the AliESDVertex for the current event
  ResetVert3D();
  AliDebug(1,"FindVertexForCurrentEvent - 3D - PROCESSING NEXT EVENT");
  fLines.Clear();

  Int_t nolines = FindTracklets(itsClusterTree,0);
  fCurrentVertex = 0;
  if(nolines>=2){
    Int_t rc=Prepare3DVertex(0);
    if(rc==0) fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0);
  /*  uncomment to debug
      printf("Vertex found in first iteration:\n");
      fVert3D.Print();
      printf("Start second iteration\n");
  end of debug lines  */
    if(fVert3D.GetNContributors()>0){
      fLines.Clear("C");
      nolines = FindTracklets(itsClusterTree,1);
      if(nolines>=2){
	rc=Prepare3DVertex(1);
	if(rc==0) fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0);
      }
    }
  /*  uncomment to debug 
      printf("Vertex found in second iteration:\n");
      fVert3D.Print();
      end of debug lines  */ 
 
    Float_t vRadius=TMath::Sqrt(fVert3D.GetXv()*fVert3D.GetXv()+fVert3D.GetYv()*fVert3D.GetYv());
    if(vRadius<GetPipeRadius() && fVert3D.GetNContributors()>0){
      Double_t position[3]={fVert3D.GetXv(),fVert3D.GetYv(),fVert3D.GetZv()};
      Double_t covmatrix[6];
      fVert3D.GetCovMatrix(covmatrix);
      Double_t chi2=99999.;
      Int_t    nContr=fVert3D.GetNContributors();
      fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
      fCurrentVertex->SetTitle("vertexer: 3D");
      fCurrentVertex->SetName("SPDVertex3D");
      fCurrentVertex->SetDispersion(fVert3D.GetDispersion());

      switch(fPileupAlgo){
      case 0: PileupFromZ(); break;
      case 1: FindOther3DVertices(itsClusterTree); break;
      default: AliError("Wrong pileup algorithm"); break;
      }
    }
  }
  if(!fCurrentVertex){
    AliITSVertexerZ vertz(GetNominalPos()[0],GetNominalPos()[1]);
    vertz.SetDetTypeRec(GetDetTypeRec());
    AliDebug(1,"Call Vertexer Z\n");
    vertz.SetLowLimit(-fZCutDiamond);
    vertz.SetHighLimit(fZCutDiamond);
    AliESDVertex* vtxz = vertz.FindVertexForCurrentEvent(itsClusterTree);
    if(vtxz){
      Double_t position[3]={GetNominalPos()[0],GetNominalPos()[1],vtxz->GetZv()};
      Double_t covmatrix[6];
      vtxz->GetCovMatrix(covmatrix);
      Double_t chi2=99999.;
      Int_t    nContr=vtxz->GetNContributors();
      fCurrentVertex = new AliESDVertex(position,covmatrix,chi2,nContr);    
      fCurrentVertex->SetTitle("vertexer: Z");
      fCurrentVertex->SetName("SPDVertexZ");
      delete vtxz;
      //      printf("Vertexer Z success\n");
    }

  }
  FindMultiplicity(itsClusterTree);
  return fCurrentVertex;
}  

//______________________________________________________________________
Int_t AliITSVertexer3D::FindTracklets(TTree *itsClusterTree, Int_t optCuts){
  // All the possible combinations between recpoints on layer 1and 2 are
  // considered. Straight lines (=tracklets)are formed. 
  // The tracklets are processed in Prepare3DVertex

  TTree *tR = itsClusterTree;
  fDetTypeRec->SetTreeAddressR(tR);
  TClonesArray *itsRec  = 0;
  if(optCuts==0) fZHisto->Reset();
  // lc1 and gc1 are local and global coordinates for layer 1
  //  Float_t lc1[3]={0.,0.,0.};
  Float_t gc1[3]={0.,0.,0.};
  // lc2 and gc2 are local and global coordinates for layer 2
  //  Float_t lc2[3]={0.,0.,0.};
  Float_t gc2[3]={0.,0.,0.};

  itsRec = fDetTypeRec->RecPoints();
  TBranch *branch;
  branch = tR->GetBranch("ITSRecPoints");

  // Set values for cuts
  Float_t xbeam=GetNominalPos()[0]; 
  Float_t ybeam=GetNominalPos()[1];
  Float_t zvert=0.;
  Float_t deltaPhi=fCoarseDiffPhiCut;
  Float_t deltaR=fCoarseMaxRCut;
  Float_t dZmax=fZCutDiamond;
  if(optCuts==1){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    zvert=fVert3D.GetZv();
    deltaPhi = fDiffPhiMax; 
    deltaR=fMaxRCut;
    dZmax=fMaxZCut;
  }else if(optCuts==2){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    deltaPhi = fDiffPhiMax; 
    deltaR=fMaxRCut;
  }

  Int_t nrpL1 = 0;    // number of rec points on layer 1
  Int_t nrpL2 = 0;    // number of rec points on layer 2

  // By default irstL1=0 and lastL1=79
  Int_t firstL1 = AliITSgeomTGeo::GetModuleIndex(1,1,1);
  Int_t lastL1 = AliITSgeomTGeo::GetModuleIndex(2,1,1)-1;
  for(Int_t module= firstL1; module<=lastL1;module++){  // count number of recopints on layer 1
    branch->GetEvent(module);
    nrpL1+= itsRec->GetEntries();
    fDetTypeRec->ResetRecPoints();
  }
  //By default firstL2=80 and lastL2=239
  Int_t firstL2 = AliITSgeomTGeo::GetModuleIndex(2,1,1);
  Int_t lastL2 = AliITSgeomTGeo::GetModuleIndex(3,1,1)-1;
  for(Int_t module= firstL2; module<=lastL2;module++){  // count number of recopints on layer 2
    branch->GetEvent(module);
    nrpL2+= itsRec->GetEntries();
    fDetTypeRec->ResetRecPoints();
  }
  if(nrpL1 == 0 || nrpL2 == 0){
    return -1;
  }
  AliDebug(1,Form("RecPoints on Layer 1,2 = %d, %d\n",nrpL1,nrpL2));

  Double_t a[3]={xbeam,ybeam,0.}; 
  Double_t b[3]={xbeam,ybeam,10.};
  AliStrLine zeta(a,b,kTRUE);
  Float_t bField=AliTracker::GetBz()/10.; //T
  SetMeanPPtSelTracks(bField);

  Int_t nolines = 0;
  // Loop on modules of layer 1
  for(Int_t modul1= firstL1; modul1<=lastL1;modul1++){   // Loop on modules of layer 1
    if(!fUseModule[modul1]) continue;
    UShort_t ladder=int(modul1/4)+1; // ladders are numbered starting from 1
    branch->GetEvent(modul1);
    Int_t nrecp1 = itsRec->GetEntries();
    static TClonesArray prpl1("AliITSRecPoint",nrecp1);
    prpl1.SetOwner();
    for(Int_t j=0;j<nrecp1;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)itsRec->At(j);
      new(prpl1[j])AliITSRecPoint(*recp);
    }
    fDetTypeRec->ResetRecPoints();
    for(Int_t j=0;j<nrecp1;j++){
      if(j>kMaxCluPerMod) continue;
      UShort_t idClu1=modul1*kMaxCluPerMod+j;
      if(fUsedCluster.TestBitNumber(idClu1)) continue;
      AliITSRecPoint *recp1 = (AliITSRecPoint*)prpl1.At(j);
      // Local coordinates of this recpoint
      /*
      lc[0]=recp1->GetDetLocalX();
      lc[2]=recp1->GetDetLocalZ();
      */
      recp1->GetGlobalXYZ(gc1);
      Double_t phi1 = TMath::ATan2(gc1[1]-ybeam,gc1[0]-xbeam);
      if(phi1<0)phi1=2*TMath::Pi()+phi1;
      for(Int_t ladl2=0 ; ladl2<fLadOnLay2*2+1;ladl2++){
	for(Int_t k=0;k<4;k++){
	  Int_t ladmod=fLadders[ladder-1]+ladl2;
 	  if(ladmod>AliITSgeomTGeo::GetNLadders(2)) ladmod=ladmod-AliITSgeomTGeo::GetNLadders(2);
	  Int_t modul2=AliITSgeomTGeo::GetModuleIndex(2,ladmod,k+1);
	  if(!fUseModule[modul2]) continue;
	  branch->GetEvent(modul2);
	  Int_t nrecp2 = itsRec->GetEntries();
	  for(Int_t j2=0;j2<nrecp2;j2++){
	    if(j2>kMaxCluPerMod) continue;
	    UShort_t idClu2=modul2*kMaxCluPerMod+j2;
	    if(fUsedCluster.TestBitNumber(idClu2)) continue;
	    AliITSRecPoint *recp2 = (AliITSRecPoint*)itsRec->At(j2);
	    /*
	    lc2[0]=recp2->GetDetLocalX();
	    lc2[2]=recp2->GetDetLocalZ();
	    */
	    recp2->GetGlobalXYZ(gc2);
	    Double_t phi2 = TMath::ATan2(gc2[1]-ybeam,gc2[0]-xbeam);
	    if(phi2<0)phi2=2*TMath::Pi()+phi2;
	    Double_t diff = TMath::Abs(phi2-phi1); 
	    if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff; 
	    if(optCuts==0 && diff<fDiffPhiMax){
	      Float_t r1=TMath::Sqrt(gc1[0]*gc1[0]+gc1[1]*gc1[1]);
	      Float_t zc1=gc1[2];
	      Float_t r2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
	      Float_t zc2=gc2[2];
	      Float_t zr0=(r2*zc1-r1*zc2)/(r2-r1); //Z @ null radius
	      fZHisto->Fill(zr0);
	    }
	    if(diff>deltaPhi)continue;
	    AliStrLine line(gc1,gc2,kTRUE);
	    Double_t cp[3];
	    Int_t retcode = line.Cross(&zeta,cp);
	    if(retcode<0)continue;
	    Double_t dca = line.GetDCA(&zeta);
	    if(dca<0.) continue;
	    if(dca>deltaR)continue;
	    Double_t deltaZ=cp[2]-zvert;
	    if(TMath::Abs(deltaZ)>dZmax)continue;

	    if(nolines == 0){
	      if(fLines.GetEntriesFast()>0)fLines.Clear("C");
	    }
	    Float_t cov[6];
	    recp2->GetGlobalCov(cov);


	    Float_t rad1=TMath::Sqrt(gc1[0]*gc1[0]+gc1[1]*gc1[1]);
	    Float_t rad2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
	    Float_t factor=(rad1+rad2)/(rad2-rad1); //factor to account for error on tracklet direction 

	    Float_t curvErr=0;
	    if(bField>0.00001){
	      Float_t curvRadius=fMeanPtSelTrk/(0.3*bField)*100; //cm 
	      Float_t dRad=TMath::Sqrt(TMath::Power((gc1[0]-gc2[0]),2)+TMath::Power((gc1[1]-gc2[1]),2));
	      Float_t aux=dRad/2.+rad1;
	      curvErr=TMath::Sqrt(curvRadius*curvRadius-dRad*dRad/4.)-TMath::Sqrt(curvRadius*curvRadius-aux*aux); //cm
	    }

	    Float_t sigmasq[3];
	    sigmasq[0]=(cov[0]+curvErr*curvErr/2.)*factor*factor;
	    sigmasq[1]=(cov[3]+curvErr*curvErr/2.)*factor*factor;
	    sigmasq[2]=cov[5]*factor*factor;

	    // Multiple scattering
 	    Float_t beta=1.;
 	    Float_t beta2=beta*beta;
 	    Float_t p2=fMeanPSelTrk*fMeanPSelTrk;
 	    Float_t rBP=GetPipeRadius();
 	    Float_t dBP=0.08/35.3; // 800 um of Be
 	    Float_t dL1=0.01; //approx. 1% of radiation length  
 	    Float_t theta2BP=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(dBP);
 	    Float_t theta2L1=14.1*14.1/(beta2*p2*1e6)*TMath::Abs(dL1);
 	    Float_t thetaBP=TMath::Sqrt(theta2BP);
 	    Float_t thetaL1=TMath::Sqrt(theta2L1);
//  	    Float_t geomfac[3];
//  	    geomfac[0]=sin(phi1)*sin(phi1);
//  	    geomfac[1]=cos(phi1)*cos(phi1);
//  	    Float_t tgth=(gc2[2]-gc1[2])/(rad2-rad1);	    
//  	    geomfac[2]=1+tgth*tgth;
 	    for(Int_t ico=0; ico<3;ico++){    
// 	      printf("Error on coord. %d due to cov matrix+curvErr=%f\n",ico,sigmasq[ico]);
// 	      //	      sigmasq[ico]+=rad1*rad1*geomfac[ico]*theta2L1/2; // multiple scattering in layer 1
// 	      //  sigmasq[ico]+=rBP*rBP*geomfac[ico]*theta2BP/2; // multiple scattering in beam pipe
 	      sigmasq[ico]+=TMath::Power(rad1*TMath::Tan(thetaL1),2)/3.;
 	      sigmasq[ico]+=TMath::Power(rBP*TMath::Tan(thetaBP),2)/3.;

// 	      printf("Multipl. scatt. contr %d = %f (LAY1), %f (BP)\n",ico,rad1*rad1*geomfac[ico]*theta2L1/2,rBP*rBP*geomfac[ico]*theta2BP/2);
// 	      printf("Total error on coord %d = %f\n",ico,sigmasq[ico]);
 	    }
	    Float_t wmat[9]={1.,0.,0.,0.,1.,0.,0.,0.,1.};
	    if(sigmasq[0]!=0.) wmat[0]=1./sigmasq[0];
	    if(sigmasq[1]!=0.) wmat[4]=1./sigmasq[1];
	    if(sigmasq[2]!=0.) wmat[8]=1./sigmasq[2];
	    new(fLines[nolines++])AliStrLine(gc1,sigmasq,wmat,gc2,kTRUE,idClu1,idClu2);

	  }
	  fDetTypeRec->ResetRecPoints();
	}
      }
    }
    prpl1.Clear();
  }
  if(nolines == 0)return -2;
  return nolines;
}

//______________________________________________________________________
Int_t  AliITSVertexer3D::Prepare3DVertex(Int_t optCuts){
  // Finds the 3D vertex information using tracklets
  Int_t retcode = -1;

  Float_t xbeam=GetNominalPos()[0];
  Float_t ybeam=GetNominalPos()[1];
  Float_t zvert=0.;
  Float_t deltaR=fCoarseMaxRCut;
  Float_t dZmax=fZCutDiamond;
  if(optCuts==1){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    zvert=fVert3D.GetZv();
    deltaR=fMaxRCut;
    dZmax=fMaxZCut;
  }else if(optCuts==2){
    xbeam=fVert3D.GetXv();
    ybeam=fVert3D.GetYv();
    deltaR=fMaxRCut;
  }

  Int_t nbr=50;
  Float_t rl=-fCoarseMaxRCut;
  Float_t rh=fCoarseMaxRCut;
  Int_t nbz=100;
  Float_t zl=-fZCutDiamond;
  Float_t zh=fZCutDiamond;
  Float_t binsizer=(rh-rl)/nbr;
  Float_t binsizez=(zh-zl)/nbz;
  TH3F *h3d = new TH3F("h3d","xyz distribution",nbr,rl,rh,nbr,rl,rh,nbz,zl,zh);
  Int_t nbrcs=25;
  Int_t nbzcs=50;
  TH3F *h3dcs = new TH3F("h3dcs","xyz distribution",nbrcs,rl,rh,nbrcs,rl,rh,nbzcs,zl,zh);

  // cleanup of the TCLonesArray of tracklets (i.e. fakes are removed)
  Int_t *validate = new Int_t [fLines.GetEntriesFast()];
  for(Int_t i=0; i<fLines.GetEntriesFast();i++)validate[i]=0;
  for(Int_t i=0; i<fLines.GetEntriesFast()-1;i++){
    AliStrLine *l1 = (AliStrLine*)fLines.At(i);
    for(Int_t j=i+1;j<fLines.GetEntriesFast();j++){
      AliStrLine *l2 = (AliStrLine*)fLines.At(j);
      Float_t dca=l1->GetDCA(l2);
      if(dca > fDCAcut || dca<0.00001) continue;
      Double_t point[3];
      Int_t retc = l1->Cross(l2,point);
      if(retc<0)continue;
      Double_t deltaZ=point[2]-zvert;
     if(TMath::Abs(deltaZ)>dZmax)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(rad>fCoarseMaxRCut)continue;
      Double_t deltaX=point[0]-xbeam;
      Double_t deltaY=point[1]-ybeam;
      Double_t raddist=TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);
      if(raddist>deltaR)continue;
      validate[i]=1;
      validate[j]=1;
      h3d->Fill(point[0],point[1],point[2]);
      h3dcs->Fill(point[0],point[1],point[2]);
    }
  }



  Int_t numbtracklets=0;
  for(Int_t i=0; i<fLines.GetEntriesFast();i++)if(validate[i]>=1)numbtracklets++;
  if(numbtracklets<2){delete [] validate; delete h3d; delete h3dcs; return retcode; }

  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    if(validate[i]<1)fLines.RemoveAt(i);
  }
  fLines.Compress();
  AliDebug(1,Form("Number of tracklets (after compress)%d ",fLines.GetEntriesFast()));
  delete [] validate;

  //        Find peaks in histos

  Double_t peak[3]={0.,0.,0.};
  Int_t ntrkl,ntimes;
  FindPeaks(h3d,peak,ntrkl,ntimes);  
  delete h3d;

  if(optCuts==0 && ntrkl<=2){
    ntrkl=0;
    ntimes=0;
    FindPeaks(h3dcs,peak,ntrkl,ntimes);  
    binsizer=(rh-rl)/nbrcs;
    binsizez=(zh-zl)/nbzcs;
    if(ntrkl==1 || ntimes>1){delete h3dcs; return retcode;}
  }
  delete h3dcs;

  //         Second selection loop

  Float_t bs=(binsizer+binsizez)/2.;
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *l1 = (AliStrLine*)fLines.At(i);
    if(l1->GetDistFromPoint(peak)>2.5*bs)fLines.RemoveAt(i);
  }
  fLines.Compress();
  AliDebug(1,Form("Number of tracklets (after 2nd compression) %d",fLines.GetEntriesFast()));

  if(fLines.GetEntriesFast()>1){
    //  find a first candidate for the primary vertex
    fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0); 
    // make a further selection on tracklets based on this first candidate
    fVert3D.GetXYZ(peak);
    AliDebug(1,Form("FIRST V candidate: %f ; %f ; %f",peak[0],peak[1],peak[2]));
    for(Int_t i=0; i<fLines.GetEntriesFast();i++){
      AliStrLine *l1 = (AliStrLine*)fLines.At(i);
      if(l1->GetDistFromPoint(peak)> fDCAcut)fLines.RemoveAt(i);
    }
    fLines.Compress();
    AliDebug(1,Form("Number of tracklets (after 3rd compression) %d",fLines.GetEntriesFast()));
    if(fLines.GetEntriesFast()>1) retcode=0; // this new tracklet selection is used
    else retcode =1; // the previous tracklet selection will be used
  }
  else {
    retcode = 0;
  }
  return retcode;  
}

//________________________________________________________
void AliITSVertexer3D::SetMeanPPtSelTracks(Float_t fieldTesla){
  // Sets mean values of P and Pt based on the field
  if(TMath::Abs(fieldTesla-0.5)<0.01){
    SetMeanPSelTracks(0.885);
    SetMeanPtSelTracks(0.630);
  }else if(TMath::Abs(fieldTesla-0.4)<0.01){
    SetMeanPSelTracks(0.805);
    SetMeanPtSelTracks(0.580);
  }else if(TMath::Abs(fieldTesla-0.2)<0.01){
    SetMeanPSelTracks(0.740);
    SetMeanPtSelTracks(0.530);
  }else if(fieldTesla<0.00001){
    SetMeanPSelTracks(0.730);
    SetMeanPtSelTracks(0.510);
  }else{
    SetMeanPSelTracks();
    SetMeanPtSelTracks();
  }
}

//________________________________________________________
void AliITSVertexer3D::FindPeaks(TH3F* histo, Double_t *peak, Int_t &nOfTracklets, Int_t &nOfTimes){
  // Finds bin with max contents in 3D histo of tracket intersections
  TAxis *xax = histo->GetXaxis();  
  TAxis *yax = histo->GetYaxis();
  TAxis *zax = histo->GetZaxis();
  peak[0]=0.;
  peak[1]=0.;
  peak[2]=0.;
  nOfTracklets = 0;
  nOfTimes=0;
  for(Int_t i=xax->GetFirst();i<=xax->GetLast();i++){
    Float_t xval = xax->GetBinCenter(i);
    for(Int_t j=yax->GetFirst();j<=yax->GetLast();j++){
      Float_t yval = yax->GetBinCenter(j);
      for(Int_t k=zax->GetFirst();k<=zax->GetLast();k++){
	Float_t zval = zax->GetBinCenter(k);
	Int_t bc =(Int_t)histo->GetBinContent(i,j,k);
	if(bc>nOfTracklets){
	  nOfTracklets = bc;
	  peak[2] = zval;
	  peak[1] = yval;
	  peak[0] = xval;
	  nOfTimes = 1;
	}else if(bc==nOfTracklets){
	  nOfTimes++;
	}
      }
    }
  }
}

//________________________________________________________
void AliITSVertexer3D::MarkUsedClusters(){
  // Mark clusters of tracklets used in vertex claulation
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *lin = (AliStrLine*)fLines.At(i);
    Int_t idClu1=lin->GetIdPoint(0);
    Int_t idClu2=lin->GetIdPoint(1);
    fUsedCluster.SetBitNumber(idClu1);
    fUsedCluster.SetBitNumber(idClu2);
  }
}
//________________________________________________________
Int_t AliITSVertexer3D::RemoveTracklets(){
  // Remove trackelts close to first found vertex
  Double_t vert[3]={fVert3D.GetXv(),fVert3D.GetYv(),fVert3D.GetZv()};
  Int_t nRemoved=0;
  for(Int_t i=0; i<fLines.GetEntriesFast();i++){
    AliStrLine *lin = (AliStrLine*)fLines.At(i);
    if(lin->GetDistFromPoint(vert)<fDCAforPileup){
      Int_t idClu1=lin->GetIdPoint(0);
      Int_t idClu2=lin->GetIdPoint(1);
      fUsedCluster.SetBitNumber(idClu1);
      fUsedCluster.SetBitNumber(idClu2);
      fLines.RemoveAt(i);
      ++nRemoved;
    }
  }
  fLines.Compress();
  return nRemoved;
}
//________________________________________________________
void AliITSVertexer3D::FindOther3DVertices(TTree *itsClusterTree){
  // pileup identification based on 3D vertexing with not used clusters
  MarkUsedClusters();
  fLines.Clear("C");
  Int_t nolines = FindTracklets(itsClusterTree,2);
  if(nolines>=2){
    Int_t nr=RemoveTracklets();
    nolines-=nr;
    if(nolines>=2){
      Int_t rc=Prepare3DVertex(2);
      if(rc==0){ 
	fVert3D=AliVertexerTracks::TrackletVertexFinder(&fLines,0);
	if(fVert3D.GetNContributors()>=fMinTrackletsForPilup){
	  fIsPileup=kTRUE;
	  fZpuv=fVert3D.GetZv();
	  fNTrpuv=fVert3D.GetNContributors();
	}
      }
    }
  }
}
//______________________________________________________________________
void AliITSVertexer3D::PileupFromZ(){
  // Calls the pileup algorithm of ALiITSVertexerZ
  Int_t binmin, binmax;
  Int_t nPeaks=AliITSVertexerZ::GetPeakRegion(fZHisto,binmin,binmax);   
  if(nPeaks==2)AliWarning("2 peaks found");
  Int_t firstPeakCont=0;
  Float_t firstPeakPos=0.;
  for(Int_t i=binmin-1;i<=binmax+1;i++){
    firstPeakCont+=fZHisto->GetBinContent(i);
    firstPeakPos+=fZHisto->GetBinContent(i)*fZHisto->GetBinCenter(i);
  }
  firstPeakPos/=firstPeakCont;
  Int_t ncontr2=0;
  if(firstPeakCont>fMinTrackletsForPilup){     
    Float_t secPeakPos;
    ncontr2=AliITSVertexerZ::FindSecondPeak(fZHisto,binmin,binmax,secPeakPos);
    if(ncontr2>=fMinTrackletsForPilup){ 
      fIsPileup=kTRUE;
      fZpuv=secPeakPos;
      fNTrpuv=ncontr2;
    }
  }
}
//________________________________________________________
void AliITSVertexer3D::PrintStatus() const {
  // Print current status
  printf("=======================================================\n");
  printf("Loose cut on Delta Phi %f\n",fCoarseDiffPhiCut);
  printf("Cut on tracklet DCA to Z axis %f\n",fCoarseMaxRCut);
  printf("Cut on tracklet DCA to beam axis %f\n",fMaxRCut);
  printf("Cut on diamond (Z) %f\n",fZCutDiamond);
  printf("Cut on DCA - tracklet to tracklet and to vertex %f\n",fDCAcut);
  printf("Max Phi difference: %f\n",fDiffPhiMax);
  printf("Pileup algo: %d\n",fPileupAlgo);
  printf("Min DCA to 1st vetrtex for pileup: %f\n",fDCAforPileup);
  printf("=======================================================\n");
}
