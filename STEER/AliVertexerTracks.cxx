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

//-----------------------------------------------------------------
//    Implementation of the vertexer from ESD tracks
//
// Origin: AliITSVertexerTracks
//         A.Dainese, Padova, 
//         andrea.dainese@pd.infn.it
//         M.Masera,  Torino, 
//         massimo.masera@to.infn.it
// Moved to STEER and adapted to ESD tracks: 
//          F.Prino,  Torino, prino@to.infn.it
//-----------------------------------------------------------------

//---- Root headers --------
#include <TTree.h>
//---- AliRoot headers -----
#include "AliStrLine.h"
#include "AliVertexerTracks.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp(AliVertexerTracks)


//----------------------------------------------------------------------------
AliVertexerTracks::AliVertexerTracks():
 TObject(),fVert() 
{
//
// Default constructor
//
  SetVtxStart();
  SetMinTracks();
  fDCAcut=0;
  fAlgo=1;
}
//-----------------------------------------------------------------------------
AliVertexerTracks::AliVertexerTracks(Double_t xStart, Double_t yStart):
 TObject(),fVert() 
{
//
// Standard constructor
//
  SetVtxStart(xStart,yStart);
  SetMinTracks();
  fDCAcut=0;
  fAlgo=1;
}
//-----------------------------------------------------------------------------
AliVertexerTracks::~AliVertexerTracks() {
  // Default Destructor
  // The objects poited by the following pointers are not owned
  // by this class and are not deleted
}
//----------------------------------------------------------------------------
Int_t AliVertexerTracks::PrepareTracks(TTree &trkTree) {
//
// Propagate tracks to initial vertex position and store them in a TObjArray
//
  Double_t alpha,xlStart;
  Int_t    nTrks    = 0;

  Int_t    nEntries = (Int_t)trkTree.GetEntries();
  if(!fTrkArray.IsEmpty()) fTrkArray.Clear();
  fTrkArray.Expand(nEntries);

  for(Int_t i=0; i<nEntries; i++) {
    AliESDtrack *track = new AliESDtrack; 
    trkTree.SetBranchAddress("tracks",&track);
    trkTree.GetEvent(i);

    // propagate track to vtxSeed
    alpha  = track->GetAlpha();
    xlStart = fNominalPos[0]*TMath::Cos(alpha)+fNominalPos[1]*TMath::Sin(alpha);
    track->PropagateTo(xlStart,GetField());   // to vtxSeed
    fTrkArray.AddLast(track);
    nTrks++; 
  }
  return nTrks;
} 
//----------------------------------------------------------------------------
AliVertex* AliVertexerTracks::VertexForSelectedTracks(TTree *trkTree) {
//
// Return vertex from tracks in trkTree
//

  // get tracks and propagate them to initial vertex position
  Int_t nTrks = PrepareTracks(*trkTree);
  if(nTrks < fMinTracks) {
    printf("TooFewTracks\n");
    Double_t vtx[3]={0,0,0};
    fVert.SetXYZ(vtx);
    fVert.SetDispersion(999);
    fVert.SetNContributors(-5);
    return &fVert;
  }
 
  // Set initial vertex position from ESD
  if(fAlgo==1)  StrLinVertexFinderMinDist(1);
  if(fAlgo==2)  StrLinVertexFinderMinDist(0);
  if(fAlgo==3)  HelixVertexFinder();
  if(fAlgo==4)  VertexFinder(1);
  if(fAlgo==5)  VertexFinder(0);
  return &fVert;
}

//----------------------------------------------------------------------------
AliESDVertex* AliVertexerTracks::FindVertex(const AliESD *event) {
//
// This is a simple wrapping (by Jouri.Belikov@cern.ch) over the original
// code by the authors of this class.
//
  Int_t nt=event->GetNumberOfTracks(), nacc=0;
  while (nt--) {
    AliESDtrack *t=event->GetTrack(nt);
    if ((t->GetStatus()&AliESDtrack::kITSrefit)==0) continue;
    fTrkArray.AddLast(t);
    nacc++;
  }

  // get tracks and propagate them to initial vertex position
  if(nacc < fMinTracks) {
    printf("TooFewTracks\n");
    Double_t vtx[3]={0,0,0};
    fVert.SetXYZ(vtx);
    fVert.SetDispersion(999);
    fVert.SetNContributors(-5);
  } else 
    switch (fAlgo) {
       case 1: StrLinVertexFinderMinDist(1); break;
       case 2: StrLinVertexFinderMinDist(0); break;
       case 3: HelixVertexFinder();          break;
       case 4: VertexFinder(1);              break;
       case 5: VertexFinder(0);              break;
       default: printf("Wrong algorithm\n"); break;  
    }

  fTrkArray.Clear();
  return &fVert;
}


//----------------------------------------------------------------------------
AliVertex* AliVertexerTracks::VertexForSelectedTracks(TObjArray *trkArray) {
//
// Return vertex from array of tracks
//

  // get tracks and propagate them to initial vertex position
  Int_t nTrks = trkArray->GetEntriesFast();
  if(nTrks < fMinTracks) {
    printf("TooFewTracks\n");
    Double_t vtx[3]={0,0,0};
    fVert.SetXYZ(vtx);
    fVert.SetDispersion(999);
    fVert.SetNContributors(-5);
    return &fVert;
  }
  TTree *trkTree = new TTree("TreeT","tracks");
  AliESDtrack *esdTrack = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdTrack);
  for(Int_t i=0; i<nTrks; i++){
    esdTrack = (AliESDtrack*)trkArray->At(i);
    trkTree->Fill();
  }
    
  AliVertex *vtx =  VertexForSelectedTracks(trkTree);
  delete trkTree;
  return vtx;
}

//---------------------------------------------------------------------------
void AliVertexerTracks::VertexFinder(Int_t optUseWeights) {

  // Get estimate of vertex position in (x,y) from tracks DCA
 
  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();
  Double_t aver[3]={0.,0.,0.};
  Double_t aversq[3]={0.,0.,0.};
  Double_t sigmasq[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliESDtrack *track1;
  AliESDtrack *track2;
  Double_t alpha,mindist;
  Double_t field=GetField();

  for(Int_t i=0; i<nacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    alpha=track1->GetAlpha();
    mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];

    Double_t pos[3]; track1->GetXYZAt(mindist,field,pos);
    Double_t dir[3]; track1->GetPxPyPzAt(mindist,field,dir);
    AliStrLine *line1 = new AliStrLine(pos,dir);
   
    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliESDtrack*)fTrkArray.At(j);
      alpha=track2->GetAlpha();
      mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];

      Double_t pos[3]; track2->GetXYZAt(mindist,field,pos);
      Double_t dir[3]; track2->GetPxPyPzAt(mindist,field,dir);
      AliStrLine *line2 = new AliStrLine(pos,dir);

      Double_t distCA=line2->GetDCA(line1);
      if(fDCAcut<=0 || (fDCAcut>0&&distCA<fDCAcut)){
	Double_t pnt1[3],pnt2[3],crosspoint[3];

	if(optUseWeights<=0){
	  Int_t retcode = line2->Cross(line1,crosspoint);
	  if(retcode>=0){
	    ncombi++;
	    for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	    for(Int_t jj=0;jj<3;jj++)aversq[jj]+=(crosspoint[jj]*crosspoint[jj]);
	  }
	}
	if(optUseWeights>0){
	  Int_t retcode = line1->CrossPoints(line2,pnt1,pnt2);
	  if(retcode>=0){
	    Double_t alpha, cs, sn;
	    alpha=track1->GetAlpha();
	    cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);	  
	    Double_t sx1=sn*sn*track1->GetSigmaY2(), sy1=cs*cs*track1->GetSigmaY2();
	    alpha=track2->GetAlpha();
	    cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	    Double_t sx2=sn*sn*track2->GetSigmaY2(), sy2=cs*cs*track2->GetSigmaY2();
	    Double_t sz1=track1->GetSigmaZ2(), sz2=track2->GetSigmaZ2();
	    Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
	    Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
	    Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
	    crosspoint[0]=wx1*pnt1[0] + wx2*pnt2[0]; 
	    crosspoint[1]=wy1*pnt1[1] + wy2*pnt2[1]; 
	    crosspoint[2]=wz1*pnt1[2] + wz2*pnt2[2];
	  
	    ncombi++;
	    for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	    for(Int_t jj=0;jj<3;jj++)aversq[jj]+=(crosspoint[jj]*crosspoint[jj]);
	  }
	}
      }
      delete line2;
    }
    delete line1;
  }
  if(ncombi>0){
    for(Int_t jj=0;jj<3;jj++){
      initPos[jj] = aver[jj]/ncombi;
      aversq[jj]/=ncombi;
      sigmasq[jj]=aversq[jj]-initPos[jj]*initPos[jj];
      sigma+=sigmasq[jj];
    }
    sigma=TMath::Sqrt(TMath::Abs(sigma));
  }
  else {
    Warning("VertexFinder","Finder did not succed");
    sigma=999;
  }
  fVert.SetXYZ(initPos);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(ncombi);
}
//---------------------------------------------------------------------------
void AliVertexerTracks::HelixVertexFinder() {

  // Get estimate of vertex position in (x,y) from tracks DCA


  Double_t initPos[3];
  initPos[2] = 0.;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  Double_t field=GetField();

  Int_t nacc = (Int_t)fTrkArray.GetEntriesFast();

  Double_t aver[3]={0.,0.,0.};
  Double_t averquad[3]={0.,0.,0.};
  Double_t sigmaquad[3]={0.,0.,0.};
  Double_t sigma=0;
  Int_t ncombi = 0;
  AliESDtrack *track1;
  AliESDtrack *track2;
  Double_t distCA;
  Double_t x, par[5];
  Double_t alpha, cs, sn;
  Double_t crosspoint[3];
  for(Int_t i=0; i<nacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    

    for(Int_t j=i+1; j<nacc; j++){
      track2 = (AliESDtrack*)fTrkArray.At(j);

      distCA=track2->PropagateToDCA(track1,field);

      if(fDCAcut<=0 ||(fDCAcut>0&&distCA<fDCAcut)){
	track1->GetExternalParameters(x,par);
	alpha=track1->GetAlpha();
	cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	Double_t x1=x*cs - par[0]*sn;
	Double_t y1=x*sn + par[0]*cs;
	Double_t z1=par[1];
	Double_t sx1=sn*sn*track1->GetSigmaY2(), sy1=cs*cs*track1->GetSigmaY2(); 
	track2->GetExternalParameters(x,par);
	alpha=track2->GetAlpha();
	cs=TMath::Cos(alpha); sn=TMath::Sin(alpha);
	Double_t x2=x*cs - par[0]*sn;
	Double_t y2=x*sn + par[0]*cs;
	Double_t z2=par[1];
	Double_t sx2=sn*sn*track2->GetSigmaY2(), sy2=cs*cs*track2->GetSigmaY2();
	Double_t sz1=track1->GetSigmaZ2(), sz2=track2->GetSigmaZ2();
	Double_t wx1=sx2/(sx1+sx2), wx2=1.- wx1;
	Double_t wy1=sy2/(sy1+sy2), wy2=1.- wy1;
	Double_t wz1=sz2/(sz1+sz2), wz2=1.- wz1;
	crosspoint[0]=wx1*x1 + wx2*x2; 
	crosspoint[1]=wy1*y1 + wy2*y2; 
	crosspoint[2]=wz1*z1 + wz2*z2;

	ncombi++;
	for(Int_t jj=0;jj<3;jj++)aver[jj]+=crosspoint[jj];
	for(Int_t jj=0;jj<3;jj++)averquad[jj]+=(crosspoint[jj]*crosspoint[jj]);
      }
    }
      
  }
  if(ncombi>0){
    for(Int_t jj=0;jj<3;jj++){
      initPos[jj] = aver[jj]/ncombi;
      averquad[jj]/=ncombi;
      sigmaquad[jj]=averquad[jj]-initPos[jj]*initPos[jj];
      sigma+=sigmaquad[jj];
    }
    sigma=TMath::Sqrt(TMath::Abs(sigma));
  }
  else {
    Warning("HelixVertexFinder","Finder did not succed");
    sigma=999;
  }
  fVert.SetXYZ(initPos);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(ncombi);
}
//---------------------------------------------------------------------------
void AliVertexerTracks::StrLinVertexFinderMinDist(Int_t optUseWeights){

  // Calculate the point at minimum distance to prepared tracks 
  
  Double_t initPos[3];
  initPos[2] = 0.;
  Double_t sigma=0;
  for(Int_t i=0;i<2;i++)initPos[i]=fNominalPos[i];
  const Int_t knacc = (Int_t)fTrkArray.GetEntriesFast();
  Double_t field=GetField();

  AliESDtrack *track1;
  Double_t (*vectP0)[3]=new Double_t [knacc][3];
  Double_t (*vectP1)[3]=new Double_t [knacc][3];
  
  Double_t sum[3][3];
  Double_t dsum[3]={0,0,0};
  for(Int_t i=0;i<3;i++)
    for(Int_t j=0;j<3;j++)sum[i][j]=0;
  for(Int_t i=0; i<knacc; i++){
    track1 = (AliESDtrack*)fTrkArray.At(i);
    Double_t alpha=track1->GetAlpha();
    Double_t mindist = TMath::Cos(alpha)*fNominalPos[0]+TMath::Sin(alpha)*fNominalPos[1];

    Double_t pos[3]; track1->GetXYZAt(mindist,field,pos);
    Double_t dir[3]; track1->GetPxPyPzAt(mindist,field,dir);
    AliStrLine *line1 = new AliStrLine(pos,dir);

    Double_t p0[3],cd[3];
    line1->GetP0(p0);
    line1->GetCd(cd);
    Double_t p1[3]={p0[0]+cd[0],p0[1]+cd[1],p0[2]+cd[2]};
    vectP0[i][0]=p0[0];
    vectP0[i][1]=p0[1];
    vectP0[i][2]=p0[2];
    vectP1[i][0]=p1[0];
    vectP1[i][1]=p1[1];
    vectP1[i][2]=p1[2];
    
    Double_t matr[3][3];
    Double_t dknow[3];
    if(optUseWeights==0)GetStrLinDerivMatrix(p0,p1,matr,dknow);
    if(optUseWeights==1){
      Double_t sigmasq[3];
      sigmasq[0]=track1->GetSigmaY2();
      sigmasq[1]=track1->GetSigmaY2();
      sigmasq[2]=track1->GetSigmaZ2();
      GetStrLinDerivMatrix(p0,p1,sigmasq,matr,dknow);
    }

    for(Int_t iii=0;iii<3;iii++){
      dsum[iii]+=dknow[iii]; 
      for(Int_t lj=0;lj<3;lj++) sum[iii][lj]+=matr[iii][lj];
    }
    delete line1;
  }
  
  Double_t vett[3][3];
  Double_t det=GetDeterminant3X3(sum);
  
   if(det!=0){
     for(Int_t zz=0;zz<3;zz++){
       for(Int_t ww=0;ww<3;ww++){
	 for(Int_t kk=0;kk<3;kk++) vett[ww][kk]=sum[ww][kk];
       }
       for(Int_t kk=0;kk<3;kk++) vett[kk][zz]=dsum[kk];
       initPos[zz]=GetDeterminant3X3(vett)/det;
     }


     for(Int_t i=0; i<knacc; i++){
       Double_t p0[3]={0,0,0},p1[3]={0,0,0};
       for(Int_t ii=0;ii<3;ii++){
	 p0[ii]=vectP0[i][ii];
	 p1[ii]=vectP1[i][ii];
       }
       sigma+=GetStrLinMinDist(p0,p1,initPos);
     }

     sigma=TMath::Sqrt(sigma);
   }else{
    Warning("StrLinVertexFinderMinDist","Finder did not succed");
    sigma=999;
  }
  delete vectP0;
  delete vectP1;
  fVert.SetXYZ(initPos);
  fVert.SetDispersion(sigma);
  fVert.SetNContributors(knacc);
}
//_______________________________________________________________________
Double_t AliVertexerTracks::GetDeterminant3X3(Double_t matr[][3]){
  //
  Double_t det=matr[0][0]*matr[1][1]*matr[2][2]-matr[0][0]*matr[1][2]*matr[2][1]-matr[0][1]*matr[1][0]*matr[2][2]+matr[0][1]*matr[1][2]*matr[2][0]+matr[0][2]*matr[1][0]*matr[2][1]-matr[0][2]*matr[1][1]*matr[2][0];
 return det;
}
//____________________________________________________________________________
void AliVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t m[][3],Double_t *d){

  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t kk=x12*x12+y12*y12+z12*z12;
  m[0][0]=2-2/kk*x12*x12;
  m[0][1]=-2/kk*x12*y12;
  m[0][2]=-2/kk*x12*z12;
  m[1][0]=-2/kk*x12*y12;
  m[1][1]=2-2/kk*y12*y12;
  m[1][2]=-2/kk*y12*z12;
  m[2][0]=-2/kk*x12*z12;
  m[2][1]=-2*y12*z12;
  m[2][2]=2-2/kk*z12*z12;
  d[0]=2*p0[0]-2/kk*p0[0]*x12*x12-2/kk*p0[2]*x12*z12-2/kk*p0[1]*x12*y12;
  d[1]=2*p0[1]-2/kk*p0[1]*y12*y12-2/kk*p0[0]*x12*y12-2/kk*p0[2]*z12*y12;
  d[2]=2*p0[2]-2/kk*p0[2]*z12*z12-2/kk*p0[0]*x12*z12-2/kk*p0[1]*z12*y12;

}
//____________________________________________________________________________
void AliVertexerTracks::GetStrLinDerivMatrix(Double_t *p0,Double_t *p1,Double_t *sigmasq,Double_t m[][3],Double_t *d){
  //
  Double_t x12=p1[0]-p0[0];
  Double_t y12=p1[1]-p0[1];
  Double_t z12=p1[2]-p0[2];

  Double_t den= x12*x12*sigmasq[1]*sigmasq[2]+y12*y12*sigmasq[0]*sigmasq[2]+z12*z12*sigmasq[0]*sigmasq[1];

  Double_t kk= 2*(x12*x12/sigmasq[0]+y12*y12/sigmasq[1]+z12*z12/sigmasq[2]);

  Double_t cc[3];
  cc[0]=-x12/sigmasq[0];
  cc[1]=-y12/sigmasq[1];
  cc[2]=-z12/sigmasq[2];

  Double_t ww=(-p0[0]*x12*sigmasq[1]*sigmasq[2]-p0[1]*y12*sigmasq[0]*sigmasq[2]-p0[2]*z12*sigmasq[0]*sigmasq[1])/den;

  Double_t ss= -p0[0]*cc[0]-p0[1]*cc[1]-p0[2]*cc[2];

  Double_t aa[3];
  aa[0]=x12*sigmasq[1]*sigmasq[2]/den;
  aa[1]=y12*sigmasq[0]*sigmasq[2]/den;
  aa[2]=z12*sigmasq[0]*sigmasq[1]/den;

  m[0][0]=aa[0]*(aa[0]*kk+2*cc[0])+2*cc[0]*aa[0]+2/sigmasq[0];
  m[0][1]=aa[1]*(aa[0]*kk+2*cc[0])+2*cc[1]*aa[0];
  m[0][2]=aa[2]*(aa[0]*kk+2*cc[0])+2*cc[2]*aa[0];

  m[1][0]=aa[0]*(aa[1]*kk+2*cc[1])+2*cc[0]*aa[1];
  m[1][1]=aa[1]*(aa[1]*kk+2*cc[1])+2*cc[1]*aa[1]+2/sigmasq[1];
  m[1][2]=aa[2]*(aa[1]*kk+2*cc[1])+2*cc[2]*aa[1];

  m[2][0]=aa[0]*(aa[2]*kk+2*cc[2])+2*cc[0]*aa[2];
  m[2][1]=aa[1]*(aa[2]*kk+2*cc[2])+2*cc[1]*aa[2];
  m[2][2]=aa[2]*(aa[2]*kk+2*cc[2])+2*cc[2]*aa[2]+2/sigmasq[2];

  d[0]=-ww*(aa[0]*kk+2*cc[0])-2*ss*aa[0]+2*p0[0]/sigmasq[0];
  d[1]=-ww*(aa[1]*kk+2*cc[1])-2*ss*aa[1]+2*p0[1]/sigmasq[1];
  d[2]=-ww*(aa[2]*kk+2*cc[2])-2*ss*aa[2]+2*p0[2]/sigmasq[2];

}
//_____________________________________________________________________________
Double_t AliVertexerTracks::GetStrLinMinDist(Double_t *p0,Double_t *p1,Double_t *x0){
  //
  Double_t x12=p0[0]-p1[0];
  Double_t y12=p0[1]-p1[1];
  Double_t z12=p0[2]-p1[2];
  Double_t x10=p0[0]-x0[0];
  Double_t y10=p0[1]-x0[1];
  Double_t z10=p0[2]-x0[2];
  return ((x10*x10+y10*y10+z10*z10)*(x12*x12+y12*y12+z12*z12)-(x10*x12+y10*y12+z10*z12)*(x10*x12+y10*y12+z10*z12))/(x12*x12+y12*y12+z12*z12);
}
