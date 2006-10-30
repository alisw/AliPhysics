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
#include <AliESDVertex.h>
#include <AliITSVertexer3D.h>
#include <AliStrLine.h>
#include <AliVertexerTracks.h>
#include <Riostream.h>
#include <TH3F.h>
#include <TTree.h>
#include<TClonesArray.h>
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliITSgeom.h"
#include "AliITSDetTypeRec.h"
#include "AliITSRecPoint.h"
/////////////////////////////////////////////////////////////////
// this class implements a method to determine
// the 3 coordinates of the primary vertex
// for p-p collisions 
// It can be used successfully with Pb-Pb collisions
////////////////////////////////////////////////////////////////

ClassImp(AliITSVertexer3D)

//______________________________________________________________________
AliITSVertexer3D::AliITSVertexer3D():AliITSVertexer(),
fLines(),
fVert3D(),
fCoarseDiffPhiCut(0.),
fMaxRCut(0.),
fZCutDiamond(0.),
fDCAcut(0.),
fDiffPhiMax(0.) {
  // Default constructor
  SetCoarseDiffPhiCut();
  SetMaxRCut();
  SetZCutDiamond();
  SetDCAcut();
  SetDiffPhiMax();

}

//______________________________________________________________________
AliITSVertexer3D::AliITSVertexer3D(TString fn): AliITSVertexer(fn),
fLines(),
fVert3D(),
fCoarseDiffPhiCut(0.),
fMaxRCut(0.),
fZCutDiamond(0.),
fDCAcut(0.),
fDiffPhiMax(0.)  {
  // Standard constructor
  fLines = new TClonesArray("AliStrLine",1000);
  SetCoarseDiffPhiCut();
  SetMaxRCut();
  SetZCutDiamond();
  SetDCAcut();
  SetDiffPhiMax();
}

//______________________________________________________________________
AliITSVertexer3D::~AliITSVertexer3D() {
  // Destructor
  fVert3D = 0;  // this object is not owned by this class
  if(fLines){
    fLines->Delete();
    delete fLines;
  }
}

//______________________________________________________________________
AliESDVertex* AliITSVertexer3D::FindVertexForCurrentEvent(Int_t evnumber){
  // Defines the AliESDVertex for the current event
  if (fDebug) cout<<"\n \n FindVertexForCurrentEvent - 3D - PROCESSING EVENT "<<evnumber<<endl;
  fVert3D = 0;
  if(fLines)fLines->Clear();

  Int_t nolines = FindTracklets(evnumber,0);
  fCurrentVertex = 0;
  if(nolines<2)return fCurrentVertex;
  Int_t rc=Prepare3DVertex();
  if(rc==0)Find3DVertex();
  if(fVert3D){
    if(fLines) fLines->Delete();
    nolines = FindTracklets(evnumber,1);
    if(nolines<2)return fCurrentVertex;
    rc=Prepare3DVertex();
    if(rc==0)Find3DVertex();
  }
  
  if(fVert3D){
    fCurrentVertex = new AliESDVertex();
    fCurrentVertex->SetTitle("vertexer: 3D");
    fCurrentVertex->SetName("Vertex");
    fCurrentVertex->SetXv(fVert3D->GetXv());
    fCurrentVertex->SetYv(fVert3D->GetYv());
    fCurrentVertex->SetZv(fVert3D->GetZv());
    fCurrentVertex->SetDispersion(fVert3D->GetDispersion());
    fCurrentVertex->SetNContributors(fVert3D->GetNContributors());
  }
  FindMultiplicity(evnumber);
  return fCurrentVertex;
}  

//______________________________________________________________________
Int_t AliITSVertexer3D::FindTracklets(Int_t evnumber, Int_t optCuts){
  // All the possible combinations between recpoints on layer 1and 2 are
  // considered. Straight lines (=tracklets)are formed. 
  // The tracklets are processed in Prepare3DVertex
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

  // Set values for cuts
  Float_t xbeam=0., ybeam=0.;
  Float_t deltaPhi=fCoarseDiffPhiCut;
  if(optCuts){
    xbeam=fVert3D->GetXv();
    ybeam=fVert3D->GetYv();
    deltaPhi = fDiffPhiMax; 
  }

  Int_t nrpL1 = 0;    // number of rec points on layer 1
  Int_t nrpL2 = 0;    // number of rec points on layer 2

  // By default irstL1=0 and lastL1=79
  Int_t irstL1 = geom->GetStartDet(0);
  Int_t lastL1 = geom->GetModuleIndex(2,1,1)-1;
  for(Int_t module= irstL1; module<=lastL1;module++){  // count number of recopints on layer 1
    branch->GetEvent(module);
    nrpL1+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  //By default irstL2=80 and lastL2=239
  Int_t irstL2 = geom->GetModuleIndex(2,1,1);
  Int_t lastL2 = geom->GetLastDet(0);
  for(Int_t module= irstL2; module<=lastL2;module++){  // count number of recopints on layer 2
    branch->GetEvent(module);
    nrpL2+= itsRec->GetEntries();
    detTypeRec.ResetRecPoints();
  }
  if(nrpL1 == 0 || nrpL2 == 0){
    itsLoader->UnloadRecPoints();
    return -1;
  }
  if(fDebug) printf("RecPoints on Layer 1,2 = %d, %d\n",nrpL1,nrpL2);

  Double_t a[3]={xbeam,ybeam,0.}; 
  Double_t b[3]={xbeam,ybeam,10.};
  AliStrLine zeta(a,b,kTRUE);

  Int_t nolines = 0;
  // Loop on modules of layer 1
  for(Int_t modul1= irstL1; modul1<=lastL1;modul1++){   // Loop on modules of layer 1
    branch->GetEvent(modul1);
    Int_t nrecp1 = itsRec->GetEntries();
    TClonesArray *prpl1 = new TClonesArray("AliITSRecPoint",nrecp1);
    prpl1->SetOwner();
    TClonesArray &rpl1 = *prpl1;
    for(Int_t j=0;j<nrecp1;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)itsRec->At(j);
      new(rpl1[j])AliITSRecPoint(*recp);
    }
    detTypeRec.ResetRecPoints();
    for(Int_t j=0;j<nrecp1;j++){
      AliITSRecPoint *recp = (AliITSRecPoint*)prpl1->At(j);
      // Local coordinates of this recpoint
      lc[0]=recp->GetDetLocalX();
      lc[2]=recp->GetDetLocalZ();
      geom->LtoG(modul1,lc,gc); // global coordinates
      Double_t phi1 = TMath::ATan2(gc[1]-ybeam,gc[0]-xbeam);
      if(phi1<0)phi1=2*TMath::Pi()+phi1;
      for(Int_t modul2= irstL2; modul2<=lastL2;modul2++){
	branch->GetEvent(modul2);
	Int_t nrecp2 = itsRec->GetEntries();
	for(Int_t j2=0;j2<nrecp2;j2++){
	  recp = (AliITSRecPoint*)itsRec->At(j2);
	  lc2[0]=recp->GetDetLocalX();
	  lc2[2]=recp->GetDetLocalZ();
	  geom->LtoG(modul2,lc2,gc2);
	  Double_t phi2 = TMath::ATan2(gc2[1]-ybeam,gc2[0]-xbeam);
	  if(phi2<0)phi2=2*TMath::Pi()+phi2;
          Double_t diff = TMath::Abs(phi2-phi1); 
          if(diff>TMath::Pi())diff=2.*TMath::Pi()-diff; 
          if(diff>deltaPhi)continue;
	  AliStrLine line(gc,gc2,kTRUE);
	  Double_t cp[3];
	  Int_t retcode = line.Cross(&zeta,cp);
	  if(retcode<0)continue;
	  Double_t dca = line.GetDCA(&zeta);
	  if(dca<0.) continue;
	  if(dca>fMaxRCut)continue;
	  if(cp[2]>fZCutDiamond || cp[2]<-fZCutDiamond)continue;
	  MakeTracklet(gc,gc2,nolines);
	}
	detTypeRec.ResetRecPoints();
      }
    }
    delete prpl1;
  }
  if(nolines == 0)return -2;
  return nolines;
}

//______________________________________________________________________
void AliITSVertexer3D::FindVertices(){
  // computes the vertices of the events in the range FirstEvent - LastEvent
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliITSLoader* itsLoader =  (AliITSLoader*) rl->GetLoader("ITSLoader");
  itsLoader->ReloadRecPoints();
  for(Int_t i=fFirstEvent;i<=fLastEvent;i++){
    rl->GetEvent(i);
    FindVertexForCurrentEvent(i);
    if(fCurrentVertex){
      WriteCurrentVertex();
    }
  }
}


//______________________________________________________________________
void AliITSVertexer3D::Find3DVertex(){
  // Determines the vertex position 
  // same algorithm as in AliVertexerTracks::StrLinVertexFinderMinDist(0)
  // adapted to pure AliStrLine objects

  Int_t knacc = fLines->GetEntriesFast();  
  Double_t (*vectP0)[3]=new Double_t [knacc][3];
  Double_t (*vectP1)[3]=new Double_t [knacc][3];
  
  Double_t sum[3][3];
  Double_t dsum[3]={0,0,0};
  for(Int_t i=0;i<3;i++)
    for(Int_t j=0;j<3;j++)sum[i][j]=0;
  for(Int_t i=0; i<knacc; i++){
    AliStrLine *line1 = (AliStrLine*)fLines->At(i); 

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
    AliVertexerTracks::GetStrLinDerivMatrix(p0,p1,matr,dknow);

    for(Int_t iii=0;iii<3;iii++){
      dsum[iii]+=dknow[iii]; 
      for(Int_t lj=0;lj<3;lj++) sum[iii][lj]+=matr[iii][lj];
    }
  }
  
  Double_t vett[3][3];
  Double_t det=AliVertexerTracks::GetDeterminant3X3(sum);
  Double_t initPos[3];
  Double_t sigma = 0.;
  for(Int_t i=0;i<3;i++)initPos[i]=0.;

  Int_t goodvert=0;

  if(det!=0){
    for(Int_t zz=0;zz<3;zz++){
      for(Int_t ww=0;ww<3;ww++){
	for(Int_t kk=0;kk<3;kk++) vett[ww][kk]=sum[ww][kk];
      }
      for(Int_t kk=0;kk<3;kk++) vett[kk][zz]=dsum[kk];
      initPos[zz]=AliVertexerTracks::GetDeterminant3X3(vett)/det;
    }

    Double_t rvert=TMath::Sqrt(initPos[0]*initPos[0]+initPos[1]*initPos[1]);
    Double_t zvert=initPos[2];
    if(rvert<fMaxRCut && TMath::Abs(zvert)<fZCutDiamond){
      for(Int_t i=0; i<knacc; i++){
	Double_t p0[3]={0,0,0},p1[3]={0,0,0};
	for(Int_t ii=0;ii<3;ii++){
	  p0[ii]=vectP0[i][ii];
	  p1[ii]=vectP1[i][ii];
	}
	sigma+=AliVertexerTracks::GetStrLinMinDist(p0,p1,initPos);
      }
      sigma=TMath::Sqrt(sigma);
      goodvert=1;
    }
  }
  delete vectP0;
  delete vectP1;
  if(!goodvert){
    Warning("Find3DVertex","Finder did not succed");
    for(Int_t i=0;i<3;i++)initPos[i]=0.;
    sigma=999;
    knacc=-1;
  }
  if(fVert3D){
    fVert3D-> SetXYZ(initPos);
    fVert3D->SetDispersion(sigma);
    fVert3D->SetNContributors(knacc);
  }else{
    fVert3D = new AliVertex(initPos,sigma,knacc);
  }
  if(fDebug) fVert3D->Print();
}



//______________________________________________________________________
Int_t  AliITSVertexer3D::Prepare3DVertex(){
  // Finds the 3D vertex information using tracklets
  Int_t retcode = -1;

  Int_t nbr=50;
  Float_t rl=-fMaxRCut;
  Float_t rh=fMaxRCut;
  Int_t nbz=100;
  Float_t zl=-fZCutDiamond;
  Float_t zh=fZCutDiamond;
  Float_t binsizer=(rh-rl)/nbr;
  Float_t binsizez=(zh-zl)/nbz;
  TH3F *h3d = new TH3F("h3d","xyz distribution",nbr,rl,rh,nbr,rl,rh,nbz,zl,zh);

 
  // cleanup of the TCLonesArray of tracklets (i.e. fakes are removed)
  Int_t *validate = new Int_t [fLines->GetEntriesFast()];
  for(Int_t i=0; i<fLines->GetEntriesFast();i++)validate[i]=0;
  //  Int_t itot=fLines->GetEntriesFast()-1-1;
  for(Int_t i=0; i<fLines->GetEntriesFast()-1;i++){
    //if(validate[i]==1)continue;
    AliStrLine *l1 = (AliStrLine*)fLines->At(i);
    for(Int_t j=i+1;j<fLines->GetEntriesFast();j++){
      AliStrLine *l2 = (AliStrLine*)fLines->At(j);
      if(l1->GetDCA(l2) > fDCAcut) continue;
      Double_t point[3];
      Int_t retc = l1->Cross(l2,point);
      if(retc<0)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(TMath::Abs(point[2])>fZCutDiamond)continue;
      if(rad>fMaxRCut)continue;
      validate[i]++;
      validate[j]++;
      h3d->Fill(point[0],point[1],point[2]);
    }
  }



  for(Int_t i=0; i<fLines->GetEntriesFast()-1;i++){ // remove tracklets sharing one recpoint
    if(validate[i]==0)continue;
    AliStrLine *l1 = (AliStrLine*)fLines->At(i);
    for(Int_t j=i+1;j<fLines->GetEntriesFast();j++){
    if(validate[j]==0)continue;
      AliStrLine *l2 = (AliStrLine*)fLines->At(j);
      if(l1->GetDCA(l2) > 0.00001) continue;
      Double_t point[3];
      Int_t retc = l1->Cross(l2,point);
      if(retc<0)continue;
      Double_t rad=TMath::Sqrt(point[0]*point[0]+point[1]*point[1]);
      if(rad<fMaxRCut)continue;
      validate[i]--;
      validate[j]--;
    }
  }

  Int_t numbtracklets=0;
  for(Int_t i=0; i<fLines->GetEntriesFast();i++)if(validate[i]>=1)numbtracklets++;
  if(numbtracklets<2){delete [] validate; delete h3d; return retcode; }

  for(Int_t i=0; i<fLines->GetEntriesFast();i++){
    if(validate[i]<1)fLines->RemoveAt(i);
  }
  fLines->Compress();
  if (fDebug) cout<<"Number of tracklets (after compress) "<<fLines->GetEntriesFast()<<endl;
  delete [] validate;


  // finds region of origin"
  Int_t nobi = 0;
  Float_t cont = 0.;
  Float_t cont2 = 0.;
  TAxis *xax = h3d->GetXaxis();  
  TAxis *yax = h3d->GetYaxis();
  TAxis *zax = h3d->GetZaxis();
  for(Int_t i=xax->GetFirst();i<=xax->GetLast();i++){
    for(Int_t j=yax->GetFirst();j<=yax->GetLast();j++){
      for(Int_t k=zax->GetFirst();k<=zax->GetLast();k++){
	nobi++;
	cont+=h3d->GetBinContent(i,j,k);
	cont2+=h3d->GetBinContent(i,j,k)*h3d->GetBinContent(i,j,k);
      }
    }
  }
  cont/=nobi;
  cont2/=nobi;
  cont2 = TMath::Sqrt((cont2 - cont*cont)/(nobi-1));
  Double_t peak[3]={0.,0.,0.};
  Float_t contref = 0.;
  for(Int_t i=xax->GetFirst();i<=xax->GetLast();i++){
    Float_t xval = xax->GetBinCenter(i);
    for(Int_t j=yax->GetFirst();j<=yax->GetLast();j++){
      Float_t yval = yax->GetBinCenter(j);
      for(Int_t k=zax->GetFirst();k<=zax->GetLast();k++){
	Float_t bc = h3d->GetBinContent(i,j,k);
	Float_t zval = zax->GetBinCenter(k);
	if(bc>contref){
	  contref = bc;
	  peak[2] = zval;
	  peak[1] = yval;
	  peak[0] = xval;
	}
      }
    }
  }
  delete h3d;

  //         Second selection loop
  Float_t bs=(binsizer+binsizez)/2.;
  for(Int_t i=0; i<fLines->GetEntriesFast()-1;i++){
    AliStrLine *l1 = (AliStrLine*)fLines->At(i);
    if(l1->GetDistFromPoint(peak)>2.5*bs)fLines->RemoveAt(i);
  }
  fLines->Compress();
  if (fDebug) cout<<"Number of tracklets (after 2nd compression) "<<fLines->GetEntriesFast()<<endl;


  if(fLines->GetEntriesFast()>1){
    Find3DVertex();   //  find a first candidate for the primary vertex
    // make a further selection on tracklets based on this first candidate
    fVert3D->GetXYZ(peak);
    if (fDebug) cout<<"FIRST V candidate: "<<peak[0]<<"; "<<peak[1]<<"; "<<peak[2]<<endl;
    for(Int_t i=0; i<fLines->GetEntriesFast()-1;i++){
      AliStrLine *l1 = (AliStrLine*)fLines->At(i);
      if(l1->GetDistFromPoint(peak)> fDCAcut)fLines->RemoveAt(i);
    }
    fLines->Compress();
    if (fDebug) cout<<"Number of tracklets (after 3rd compression) "<<fLines->GetEntriesFast()<<endl;
    if(fLines->GetEntriesFast()>1){
      retcode = 0;   // this last tracklet selection will be used
      delete fVert3D;
      fVert3D = 0;
    }
    else {
	retcode =1; // the previous tracklet selection will be used
    }
  }
  else {
    retcode = 0;
  }
  return retcode;  
}

 //______________________________________________________________________
void  AliITSVertexer3D::MakeTracklet(Double_t *pA, Double_t *pB, Int_t &nolines) {
  // makes a tracklet
  TClonesArray &lines = *fLines;
  if(nolines == 0){
    if(fLines->GetEntriesFast()>0)fLines->Clear();
  }
  if(fLines->GetEntriesFast()==fLines->GetSize()){
    Int_t newsize=(Int_t) 1.5*fLines->GetEntriesFast();
    fLines->Expand(newsize);
  }

  new(lines[nolines++])AliStrLine(pA,pB,kTRUE);
}

 //______________________________________________________________________
void  AliITSVertexer3D::MakeTracklet(Float_t *pA, Float_t *pB, Int_t &nolines) {// Makes a tracklet
  //
  Double_t a[3],b[3];
  for(Int_t i=0;i<3;i++){
    a[i] = pA[i];
    b[i] = pB[i];
  }
  MakeTracklet(a,b,nolines);
}

//________________________________________________________
void AliITSVertexer3D::PrintStatus() const {
  // Print current status
  cout <<"=======================================================\n";
  cout << "Loose cut on Delta Phi "<<fCoarseDiffPhiCut<<endl;
  cout << "Cut on tracklet DCA to Z axis "<<fMaxRCut<<endl;
  cout << "Cut on diamond (Z) "<<fZCutDiamond<<endl;
  cout << "Cut on DCA - tracklet to tracklet and to vertex "<<fDCAcut<<endl;
  cout <<" Max Phi difference: "<<fDiffPhiMax<<endl;

 
  cout <<"=======================================================\n";
}
