
/*
  This test macro to test AliTracker functionality.
  Additional focus on the tracking of  cosmic tracks .

  Test:
  1. Simulate random tracks  - AliTracker used to propagate track
  2. Missalign and smear space points
  3. Fit the tracks

  4. Visualize results/residual/pulls using tree draw

  The input MC paramters are stored together with reconstructed output in the tree.


  Usage:
  .L AliTrackerTest.C+g
  AliTrackerTest(10000)
  TFile f("fit.root");

  Make a simple plots:
  //Pulls in P3 parameter - tangent lambda
  fit->Draw("(seedF0.fP[3]-mcU.fP[3])/sqrt(seedF0.fC[9])>>hisSP3(100,-5,5)","")
  //Pulls of seeding for seed
  fit->Draw("(seedF0.fP[4]-mcU.fP[4])/sqrt(seedF0.fC[14])>>hisSP4(100,-5,5)","")
  //
  //Pulls in P3 parameter - tangent lambda
  fit->Draw("(trackF0.fP[3]-mcU.fP[3])/sqrt(trackF02.fC[9])>>hisP3(100,-5,5)","")
  //Pulls in P4 parameter
  fit->Draw("(trackF0.fP[4]-mcU.fP[4])/sqrt(trackF02.fC[14])>>hisP4(100,-5,5)","")

  //
  //
  // seeding procedure test (errors scaled 0.01 cm)
  fit->Draw("(seedF0.Pt()-mcD.Pt())/mcD.Pt()^2:mcD.Pt()>>hisSPt(10,0,400,100,-0.005,0.005)","abs(seedF0.fP[4])<3&&abs(pointsU.fNPoints+pointsD.fNPoints)>300","colz");
  //"standart track
  fit->Draw("(trackF0.Pt()-mcU.Pt())/mcU.Pt()^2:mcU.Pt()>>his0Pt(10,0,400,100,-0.005,0.005)","abs(trackF0.fP[4])<3&&abs(pointsU.fNPoints+pointsD.fNPoints)>300","colz");
  // misaligned track - up down
  fit->Draw("(trackF02.Pt()-mcU.Pt())/mcU.Pt()^2:mcU.Pt()>>his2Pt(10,0,400,100,-0.005,0.005)","abs(trackF02.fP[4])<3&&abs(pointsU.fNPoints+pointsD.fNPoints)>300","colz");
  hisSPt->FitSlicesY();
  his0Pt->FitSlicesY();
  his2Pt->FitSlicesY();
  
*/

//ROOT includes
#include <iostream>
#include "TRandom.h"
#include "TRandom3.h"
#include "TNtuple.h" 
#include "TStopwatch.h" 
#include "TDatabasePDG.h" 
#include "TMath.h" 
#include "TGeoManager.h" 
#include "TClonesArray.h" 
#include "TTree.h" 
#include "TFile.h" 
//AliRoot includes   
#include "AliGeomManager.h" 
#include "AliMagF.h" 
#include "AliESDVertex.h" 
#include "AliExternalTrackParam.h" 
#include "TTreeStream.h" 
#include "AliTrackPointArray.h" 
#include "AliTrackerBase.h"
#include "AliTracker.h"

//
// 
#include "refitTrack.C" 

void TestRotateMI(AliExternalTrackParam *paramIn, TTreeSRedirector *pcstream);


Int_t simulateCosmicTrack(AliExternalTrackParam &paramU, AliExternalTrackParam &paramD,
			  AliTrackPointArray &pointsU, AliTrackPointArray &pointsD,
			  AliTrackPointArray &pointsF){  
  //
  // Toy - Simulate cosmic muon track 
  // space points genearted with the step 1 cm - at given radius
  //
  // Return value    - number of points
  //        paramU   - parameters at the beginning of track
  //        paramD   - parameters at the end of track
  //        pointsU  - points in upper half of 
  //        pointsD  - points in lower part
  //        pointsF  - all points along trajectory
  // 
  Double_t rTPC1=250;
  Double_t rTPC0=80;
  Double_t kMuon = TDatabasePDG::Instance()->GetParticle("mu+")->Mass();
  Double_t kMaxSnp = 0.85;  
  Double_t xyz[3]={0,0,0};
  Double_t xyzdw[3]={0,0,0};
  Double_t pxyzup[3]={0,0,0};
  Double_t pxyzdw[3]={0,0,0};
  Double_t pxyz[3]={0,0,0};
  Double_t cv[21];
  for (Int_t i=0; i<21; i++) cv[i]=0;
  Float_t covPoint[6]={0,0,0,0,0,0};
  UShort_t volid=0;
  //
  static TClonesArray arrup("AliTrackPoint",500);
  static TClonesArray arrdw("AliTrackPoint",500);
  static TClonesArray arr("AliTrackPoint",500);
  if (arrup[0]==0){
    // init point arrays
    for (Int_t i=0;i<500;i++){
      new (arrup[i]) AliTrackPoint;
      new (arrdw[i]) AliTrackPoint;
      new (arr[i])   AliTrackPoint;
    }
  }

  //
  xyz[0]=(gRandom->Rndm()-0.5)*250;
  xyz[1]=250;
  xyz[2]=(gRandom->Rndm()-0.5)*100;
  //
  Double_t pt  = gRandom->Exp(8.);
  if (gRandom->Rndm()>0.3)  pt= gRandom->Exp(50.);
  if (gRandom->Rndm()>0.6)  pt= gRandom->Rndm()*400;
  //  Double_t pt  = gRandom->Exp(60.);
  Double_t pz  = gRandom->Gaus(0,0.3)*pt;
  Double_t phi = gRandom->Gaus(0,0.3);
  pxyz[0] = pt*TMath::Sin(phi);
  pxyz[1] = pt*TMath::Cos(phi);
  pxyz[2] = pz;
  Short_t sign= (gRandom->Rndm()>0.5)? -1:1;
  //
  pxyzup[0]=pxyz[0];
  pxyzup[1]=pxyz[1];
  pxyzup[2]=pxyz[2];
  //  
  AliExternalTrackParam lparamU(xyz, pxyzup,cv,sign);
  Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
  lparamU.Rotate(alpha);
  paramU=lparamU;
  //
  //
  //
  Double_t txyzup[3]={0,0,0};
  Double_t txyzdw[3]={0,0,0};
  //
  Int_t npointsUp=0;
  AliTrackerBase::PropagateTrackTo(&lparamU,rTPC1,kMuon,3,kTRUE,kMaxSnp);
  paramU=lparamU;
  for (Double_t r=rTPC1; r>0; r-=1){
    Bool_t status =	 AliTrackerBase::PropagateTrackTo(&lparamU,r,kMuon,3,kTRUE,kMaxSnp);
    if (!status) break;
    if (TMath::Abs(lparamU.GetSnp())>kMaxSnp) break;
    if (r<rTPC0) continue;
    lparamU.GetXYZ(txyzup);
    new (arrup[npointsUp]) AliTrackPoint(txyzup[0],txyzup[1],txyzup[2],covPoint,volid,0.,0.,0.,0);
    npointsUp++;
  }
  //
  lparamU.GetXYZ(txyzup);
  Double_t bz=AliTrackerBase::GetBz(txyzup);
  Double_t maxd=100000;
  Double_t pvtx[3]={0,0,0};
  Double_t sigmavtx[3]={0.01,0.01,3000};
  AliESDVertex *vtx= new AliESDVertex(pvtx,sigmavtx,"vertex");
  Double_t dz[2]={0,0};
  Double_t cvtx[3]={0.0,0.0,0};
  lparamU.PropagateToDCA(vtx,bz,maxd,dz,cvtx); 
  //
  // make step to other side
  lparamU.PropagateTo(-30,bz);  
  lparamU.GetXYZ(xyzdw); 
  lparamU.GetPxPyPz(pxyzdw);
  // invert the sign of the momentum
  pxyzdw[0]=-pxyzdw[0];
  pxyzdw[1]=-pxyzdw[1];
  pxyzdw[2]=-pxyzdw[2];
  Short_t sign2=-sign;
  AliExternalTrackParam lparamD(xyzdw,pxyzdw,cv,sign2);
  lparamU.GetXYZ(xyzdw); 
  lparamU.GetPxPyPz(pxyzdw);
  Double_t alphadw = TMath::ATan2(xyzdw[1],xyzdw[0]);
  lparamD.Rotate(alphadw);//I have to rotate gobal to local coordenate externalparam
  Double_t radius0=TMath::Sqrt(xyzdw[1]*xyzdw[1]+xyzdw[0]*xyzdw[0]);
  Int_t npointsDown=0;
  for (Double_t r=radius0; r<rTPC1; r+=1){ 
    Bool_t status =  AliTrackerBase::PropagateTrackTo(&lparamD,r,kMuon,3,kTRUE,0.99); 
    if (!status) continue;
    if (TMath::Abs(lparamD.GetSnp())>kMaxSnp) continue;
    if(r>rTPC0){
      lparamD.GetXYZ(txyzdw);
      new (arrdw[npointsDown]) AliTrackPoint(txyzdw[0],txyzdw[1],txyzdw[2],covPoint,volid,0.,0.,0.,0);
      npointsDown++;	
    }
  }
  //
  // Fill MC point arrays
  //  
  Int_t npoints=npointsUp+npointsDown;
  AliTrackPointArray lpointsF(npoints);
  AliTrackPointArray lpointsU(npointsUp);
  AliTrackPointArray lpointsD(npointsDown);
  for (Int_t i=0; i<npointsUp;i++){
    AliTrackPoint *point = (AliTrackPoint*)arrup[i];
    lpointsF.AddPoint(i,point);
    lpointsU.AddPoint(i,point);
  }
  //
  for (Int_t i=0; i<npointsDown;i++){
    AliTrackPoint *point = (AliTrackPoint*)arrdw[i];
    lpointsF.AddPoint(i+npointsUp,point);
    lpointsD.AddPoint(i,point);
  }
  // export points
  AliTrackerBase::PropagateTrackTo(&lparamD,rTPC1,kMuon,3,kTRUE,kMaxSnp);
  paramD=lparamD;
  //
  pointsU=lpointsU;
  pointsD=lpointsD;
  pointsF=lpointsF;
  return npoints;
}


AliTrackPointArray * SmearPoints(AliTrackPointArray &pointArray, Double_t sigmaY, Double_t sigmaZ, Double_t shiftY,Double_t shiftZ){
  //
  // Smear ideal points form the simulation
  // 1. Smear the input points (in the "local frame")
  // 2. Assing corresponding covariance matrix
  // 3. Add systematic shift - shift y and shift z
  
  Int_t  npoints=pointArray.GetNPoints();
  AliTrackPointArray *outputArray = new  AliTrackPointArray(npoints);
  Double_t xyz[3]={0,0,0}; 
  Float_t covPoint[6]={0,0,0, sigmaY*sigmaY,0,sigmaZ*sigmaZ};  //covaraince at the local frame
  //
  //
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    AliTrackPoint pointIn;
    pointArray.GetPoint(pointIn,ipoint);
    Double_t alpha = TMath::ATan2(pointIn.GetY(),pointIn.GetX());
    AliTrackPoint pr = pointIn.Rotate(alpha);
    xyz[0]=pr.GetX();                           //local x
    xyz[1]=pr.GetY()+gRandom->Gaus(0,sigmaY);   //local y
    xyz[2]=pr.GetZ()+gRandom->Gaus(0,sigmaZ);   //local z
    if (pointIn.GetY()>0) xyz[1]+=shiftY;
    if (pointIn.GetY()>0) xyz[2]+=shiftZ;
    //
    pr.SetXYZ(xyz[0],xyz[1],xyz[2],covPoint);  // set covariance matrix
    AliTrackPoint      pg= pr.Rotate(-alpha);
    AliTrackPoint prCheck= pg.Rotate(alpha);
    outputArray->AddPoint(ipoint,&pg);    
  }
  return outputArray;
}



AliExternalTrackParam *  MakeSeed(AliTrackPointArray &pointArray, Int_t seedDelta){
  //
  // Example: creation of seed
  // Make seed for array of track points
  // seedDelta - gap between seeding points
  //
  Int_t  npoints=pointArray.GetNPoints();
  if(npoints<=3) return 0;   //not enough points to make a trac
  if (npoints-2*seedDelta-1<0) return 0;
  AliTrackPoint   point1;
  AliTrackPoint   point2;
  AliTrackPoint   point3;
  pointArray.GetPoint(point1,npoints-1);
  pointArray.GetPoint(point2,npoints-seedDelta-1);
  pointArray.GetPoint(point3,npoints-2*seedDelta-1);
  //
  AliExternalTrackParam * trackParam = AliTrackerBase::MakeSeed(point1, point2, point3);
  return trackParam;
}


void AliTrackerTest(Int_t ntracks) {
  //   
  // 
  //
  TGeoManager::Import("./geometry.root");
  AliGeomManager::LoadGeometry("./geometry.root");
  AliMagF::BMap_t smag = AliMagF::k5kG;
  Double_t        kMuon = TDatabasePDG::Instance()->GetParticle("mu+")->Mass();
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., smag));
  //
  
  TTreeSRedirector *pcstream = new TTreeSRedirector("fit.root"); 
  AliExternalTrackParam paramU;
  AliExternalTrackParam paramD;
  AliTrackPointArray    pointsU; 
  AliTrackPointArray    pointsD;
  AliTrackPointArray    pointsF;  

  for (Int_t i=0;i<ntracks;i++) {
    if (i%10==0) printf("Track\t%d\n",i);
    Int_t npoints = simulateCosmicTrack(paramU, paramD, pointsU, pointsD, pointsF);  
    if (npoints<10) continue;
    //
    // point for test of seeding - with scaled errors
    AliTrackPointArray* pointsFS=SmearPoints(pointsF, 0.01,0.01,0,0);
    AliTrackPointArray* pointsDS=SmearPoints(pointsD, 0.01,0.01,0,0);
    AliTrackPointArray* pointsUS=SmearPoints(pointsU, 0.01,0.01,0,0);
    //test of seeding routines
    AliExternalTrackParam *seedF0  =  MakeSeed(*pointsFS,50);
    AliExternalTrackParam *seedD0  =  MakeSeed(*pointsDS,50);
    AliExternalTrackParam *seedU0  =  MakeSeed(*pointsUS,50);
    // points smeeared according TPC resolution
    AliTrackPointArray* pointsF0=SmearPoints(pointsF, 0.1,0.1,0,0);
    AliTrackPointArray* pointsD0=SmearPoints(pointsD, 0.1,0.1,0,0);
    AliTrackPointArray* pointsU0=SmearPoints(pointsU, 0.1,0.1,0,0);
    // points smeared according TPCresolution - + systematic shift
    AliTrackPointArray* pointsF02=SmearPoints(pointsF, 0.1,0.1,0.2,0.2);
    AliTrackPointArray* pointsD02=SmearPoints(pointsD, 0.1,0.1,0.2,0.2);
    AliTrackPointArray* pointsU02=SmearPoints(pointsU, 0.1,0.1,0.2,0.2);
    //test of tracking    
    AliExternalTrackParam *trackF0  =  MakeSeed(*pointsF0,10);
    AliExternalTrackParam *trackD0  =  MakeSeed(*pointsD0,10);
    AliExternalTrackParam *trackU0  =  MakeSeed(*pointsU0,10);
    AliExternalTrackParam *trackF02 =  MakeSeed(*pointsF02,10);
    AliExternalTrackParam *trackD02 =  MakeSeed(*pointsD02,10);
    AliExternalTrackParam *trackU02 =  MakeSeed(*pointsU02,10);
    //
    //
    if (!trackF0) continue;
    if (!trackD0) continue;
    if (!trackU0) continue;
    if (!seedF0) continue;
    if (!seedD0) continue;
    if (!seedU0) continue;    
    AliTrackerBase::FitTrack(trackF0, pointsF0,kMuon,3);
    AliTrackerBase::FitTrack(trackU0, pointsU0,kMuon,3);
    AliTrackerBase::FitTrack(trackD0, pointsD0,kMuon,3);
    AliTrackerBase::FitTrack(trackF02, pointsF02,kMuon,3);
    AliTrackerBase::FitTrack(trackD02, pointsD02,kMuon,3);
    AliTrackerBase::FitTrack(trackU02, pointsU02,kMuon,3);    
    //
    TestRotateMI(trackF0, pcstream);
    if (trackF0&&trackD0&&trackU0){
      (*pcstream)<<"fit"<<	
	"pointsD.="<<&pointsD<<
	"pointsD0.="<<pointsD0<<
	"pointsD02.="<<pointsD02<<
	//
	"pointsU.="<<&pointsU<<
	"pointsU0.="<<pointsU0<<
	"pointsU02.="<<pointsU02<<
	//
	"pointsF.="<<&pointsF<<
	"pointsF0.="<<pointsF0<<
	"pointsF02.="<<pointsF02<<
	//
	"mcU.="<<&paramU<<             //MC track down
	"mcD.="<<&paramD<<             //MC track top
	// track fit  - 0 misalignemnt
	"seedF0.="<<seedF0<<       //full track
	"seedD0.="<<seedD0<<       //down track
	"seedU0.="<<seedU0<<       //up track
	//
	"trackF0.="<<trackF0<<       //full track
	"trackD0.="<<trackD0<<       //down track
	"trackU0.="<<trackU0<<       //up track
	// track fit  - 0.2 cm  misalignemnt
	"trackF02.="<<trackF02<<       //full track
	"trackD02.="<<trackD02<<       //down track
	"trackU02.="<<trackU02<<       //up track
	"\n";
    }
  }
  delete pcstream; 
  
}

void TestRotateMI(AliExternalTrackParam *paramIn, TTreeSRedirector *pcstream){
  //
  // test of rotation function
  // rotate by 360 degrees
  // dump the state vector to the tree after each rotation 
  AliExternalTrackParam param(*paramIn);
  AliExternalTrackParam paramMI(*paramIn);
  for (Int_t idiv=0; idiv<=18; idiv++){
    Double_t alphaRot = paramIn->GetAlpha()+2*TMath::Pi()*idiv/18.;
    param.Rotate(alphaRot);
    paramMI.RotateMI(alphaRot);
    (*pcstream)<<"rotateTest"<<
      "pIn.="<<paramIn<<
      "pRot.="<<&param<<
      "pRotMI.="<<&paramMI<<
      "idiv="<<idiv<<
      "\n";
  }
  /*
    to check:
    
   */
}
