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
/*
$Log$
Revision 1.1  2002/02/18 09:06:54  kowal2
TPC parametrization by Andrea Dainese

*/
/*************************************************************************
 *                                                                        *
 * This class builds AliTPCtrack objects from generated tracks to feed    *
 * ITS tracking (V2). The AliTPCtrack is built from its first hit in      *
 * the TPC. The track is assigned a Kalman-like covariance matrix         *
 * depending on its pT and pseudorapidity and track parameters are        *
 * smeared according this covariance matrix.                              *
 * Output file contains sorted tracks, ready for matching with ITS        *
 *                                                                        *
 * Test macro is: AliBarrelRec_TPCparam.C                                 * 
 *                                                                        *
 *  Origin: Andrea Dainese, Padova - e-mail: andrea.dainese@pd.infn.it    * 
 *                                                                        *
 **************************************************************************/
#include "AliTPCtrackerParam.h"
#include "alles.h"
#include "AliMagF.h"
#include "AliTPCtrack.h"
#include "TMatrixD.h"
#include "AliKalmanTrack.h"

ClassImp(AliTPCtrackerParam)

//-----------------------------------------------------------------
AliTPCtrackerParam::AliTPCtrackerParam(const Int_t coll,const Double_t Bz) 
{
  fColl = coll; // collision code (0: PbPb6000)
  fBz = Bz;     // value of the z component of L3 field (Tesla)
  
}
//-----------------------------------------------------------------
AliTPCtrackerParam::~AliTPCtrackerParam() 
{}
//-----------------------------------------------------------------

Int_t AliTPCtrackerParam::BuildTPCtracks(const TFile *inp, TFile *out, Int_t n)
{
 
  if(fColl!=0) { 
    cerr<<"AliTPCtrackerParam::BuildTPCtracks:  Invalid collision!"<<endl<<
                "Available:  0   ->   PbPb6000"<<endl; return 0; 
  }
  if(fBz!=0.4) {
    cerr<<"AliTPCtrackerParam::BuildTPCtracks:  Invalid field!"<<endl<<
                "Available:  0.4"<<endl; return 0;
  }

  TFile *infile=(TFile*)inp;

  AliKalmanTrack::SetConvConst(100/0.299792458/fBz);

  // Get gAlice object from file
  if(!(gAlice=(AliRun*)infile->Get("gAlice"))) {
    cerr<<"gAlice has not been found on galice.root !\n";
    return 1;
  }

  // Set random number generator seed
  TDatime t;
  gRandom->SetSeed(t.Get());


  // loop over first n events in file
  for(Int_t evt=0; evt<n; evt++){
    cerr<<"+++\n+++ Processing event "<<evt<<"\n+++\n";

    AliTPCtrack *tpctrack=0;

    // tree for TPC tracks
    Char_t tname[100];
    sprintf(tname,"TreeT_TPC_%d",evt);
    TTree *tracktree = new TTree(tname,"Tree with TPC tracks");
    tracktree->Branch("tracks","AliTPCtrack",&tpctrack,20000,0);

    // array for TPC tracks
    TObjArray tarray(20000);


    Int_t nparticles=gAlice->GetEvent(evt);   

    Bool_t done[500000];
    for(Int_t l=0; l<500000; l++) { done[l]=kFALSE; }


    // Get TPC detector 
    AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
    Int_t ver = TPC->IsVersion(); 
    cerr<<"+++ TPC version "<<ver<<" has been found !\n";
    AliTPCParam *digp=(AliTPCParam*)infile->Get("75x40_100x60");
    if(!digp) { cerr<<"TPC parameters have not been found !\n"; return 1; }
    TPC->SetParam(digp);

    // Get TreeH with hits
    TTree *TH=gAlice->TreeH(); 
    Int_t ntracks=(Int_t)TH->GetEntries();
    cerr<<"+++\n+++ Number of particles in event "<<evt<<":  "<<nparticles<<"\n+++\n+++ Number of \"primary tracks\" (entries in TreeH): "<<ntracks<<"\n+++\n\n";

    TParticle* Part;
    AliTPChit* tpcHit;
    Double_t hPx,hPy,hPz,hPt,xg,yg,zg,xl,yl,zl;
    Double_t Alpha;
    Float_t cosAlpha,sinAlpha;
    Int_t label,Pdg,charge,bin;
    Int_t tracks=0;
    //Int_t N1=0,N2=0;

    // loop over entries in TreeH
    for(Int_t i=0; i<ntracks; i++) {
      if(i%1000==0) cerr<<"  --- Processing primary track "<<i<<" of "<<ntracks<<" ---\r";
      TPC->ResetHits();
      TH->GetEvent(i);
      // Get FirstHit
      tpcHit=(AliTPChit*)TPC->FirstHit(-1);
      for( ; tpcHit; tpcHit=(AliTPChit*)TPC->NextHit() ) {
	if(tpcHit->fQ !=0.) continue;
	// Get particle momentum at hit
	hPx=tpcHit->X(); hPy=tpcHit->Y(); hPz=tpcHit->Z();
	hPt=TMath::Sqrt(hPx*hPx+hPy*hPy);
	// reject hits with Pt<mag*0.45 GeV/c
	if(hPt<(fBz*0.45)) continue;

	// Get track label
	label=tpcHit->Track();
	// check if this track has already been processed
	if(done[label]) continue;
	// electric charge
	Part = gAlice->Particle(label);
	Pdg = Part->GetPdgCode();
	if(Pdg>200 || Pdg==-11 || Pdg==-13) { charge=1; }
	else if(Pdg<-200 || Pdg==11 || Pdg==13) { charge=-1; }
	else continue;


	if((tpcHit=(AliTPChit*)TPC->NextHit())==0) break;
	if(tpcHit->fQ != 0.) continue;
	// Get global coordinates of hit
	xg=tpcHit->X(); yg=tpcHit->Y(); zg=tpcHit->Z();
	if(TMath::Sqrt(xg*xg+yg*yg)>90.) continue;

	// Get TPC sector, Alpha angle and local coordinates
	// printf("Sector %d\n",tpcHit->fSector);
	digp->AdjustCosSin(tpcHit->fSector,cosAlpha,sinAlpha);
	Alpha = TMath::ATan2(sinAlpha,cosAlpha);
	xl = xg*cosAlpha + yg*sinAlpha;
	yl =-xg*sinAlpha + yg*cosAlpha;
	zl = zg;
	//printf("Alpha %f   xl %f  yl %f  zl %f\n",Alpha,xl,yl,zl);

	// reject tracks which are not in the TPC acceptance
	if(TMath::Abs(zl+(244.-xl)*hPz/hPt)>252.) continue;

	// Get bin in pT,eta
	bin = GetBin(hPt,Part->Eta());

	// Apply selection according to TPC efficiency
	//if(abs(Pdg)==211) N1++;
	if(!SelectedTrack(Pdg,hPt,Part->Eta())) continue; 
	//if(abs(Pdg)==211) N2++;

	// Mark track as "done"
	done[label]=kTRUE; tracks++;
 
	// create AliTPCtrack object
	tpctrack = BuildTrack(Alpha,xl,yl,zl,hPx,hPy,hPz,hPt,charge,label);

	// put track in array
	tarray.AddLast(tpctrack);

      }
 
    } // loop over entries in TreeH

    TObjArray newtarray(20000);

    // assing covariance matrixes and smear track parameters
    CookTracks(tarray,newtarray);

    // sort array with TPC tracks (decreasing pT)
    newtarray.Sort();


    Int_t arrentr = newtarray.GetEntriesFast();
    //printf("\n  %d  \n\n",arrentr);
    for(Int_t l=0; l<arrentr; l++) {
      tpctrack=(AliTPCtrack*)newtarray.UncheckedAt(l);
      tracktree->Fill();
    }

    // write the tree with tracks in the output file
    out->cd();
    tracktree->Write();

    delete tracktree;

    printf("\n\n+++\n+++ Number of TPC tracks: %d\n+++\n",tracks);
    //printf("Average Eff: %f\n",(Float_t)N2/N1);

  } // loop on events

  return 0;
}

//-----------------------------------------------------------------
AliTPCtrack* AliTPCtrackerParam::BuildTrack(Double_t alpha,Double_t x,
					    Double_t y,Double_t z,Double_t Px,
					    Double_t Py,Double_t Pz,Double_t Pt,
					    Int_t ch,Int_t lab) 
{  
  Double_t xref = x;
  Double_t xx[5],cc[15];
  cc[0]=cc[2]=cc[5]=cc[9]=cc[14]=10.;
  cc[1]=cc[3]=cc[4]=cc[6]=cc[7]=cc[8]=cc[10]=cc[11]=cc[12]=cc[13]=0.;
  
  // Magnetic field
  TVector3 Bfield(0.,0.,fBz);
  
  
  // radius [cm] of track projection in (x,y) 
  Double_t rho = Pt*100/0.299792458/Bfield.Z();
  // center of track projection in local reference frame
  TVector3 MomH,PosH;


  // position (local) and momentum (local) at the hit
  // in the bending plane (z=0)
  PosH.SetXYZ(x,y,0.);
  MomH.SetXYZ(Px*TMath::Cos(alpha)+Py*TMath::Sin(alpha),-Px*TMath::Sin(alpha)+Py*TMath::Cos(alpha),0.);
  TVector3 Vrho = MomH.Cross(Bfield);
  Vrho *= ch;
  Vrho.SetMag(rho);

  TVector3 Vcenter = PosH+Vrho;

  Double_t x0 = Vcenter.X();

  // fX     = xref         X-coordinate of this track (reference plane)
  // fAlpha = Alpha        Rotation angle the local (TPC sector) 
  // fP0    = YL           Y-coordinate of a track
  // fP1    = ZG           Z-coordinate of a track
  // fP2    = C*x0         x0 is center x in rotated frame
  // fP3    = Tgl          tangent of the track momentum dip angle
  // fP4    = C            track curvature
  xx[0] = y;
  xx[1] = z;
  xx[3] = Pz/Pt;
  xx[4] = -ch/rho;
  xx[2] = xx[4]*x0;

  // create the object AliTPCtrack    
  AliTPCtrack* track = new AliTPCtrack(0,xx,cc,xref,alpha);
  // set the label
  track->SetLabel(lab);

  return track;
}

//-----------------------------------------------------------------
Bool_t AliTPCtrackerParam::SelectedTrack(Int_t pdg,Double_t pt,Double_t eta) 
{
  Double_t eff=0.;
  
  //eff computed with | zl+(244-xl)*pz/pt | < 252
  Double_t effPi[27] = {0.724587,0.743389,0.619273,0.798477,0.812036,0.823195,0.771437,0.775826,0.784136,0.809071,0.762001,0.774576,0.848834,0.787201,0.792548,0.942089,0.951631,0.951085,0.960885,0.971451,0.969103,0.983245,0.978939,0.988706,0.990852,0.985679,0.993606};
  Double_t effK[18]  = {0.377934,0.363962,0.321721,0.518784,0.547459,0.517878,0.612704,0.619101,0.620894,0.733411,0.732128,0.750373,0.790630,0.806565,0.791353,0.967486,0.970483,0.974527};
  Double_t effP[15]  = {0.131173,0.165114,0.229658,0.365357,0.412989,0.483297,0.454614,0.505173,0.658615,0.694753,0.730661,0.815680,0.873461,0.887227,0.899324};
  Double_t effEl[15] = {0.835549,0.853746,0.718207,0.835230,0.831489,0.862222,0.757783,0.747301,0.824096,0.867949,0.871891,0.808480,0.890625,0.911765,0.973684};
  Double_t effMu[15] = {0.553486,0.641392,0.609932,0.591126,0.706729,0.750755,0.747952,0.729051,0.760849,0.898810,0.737500,0.830357,0.735294,0.800000,0.882353};
 

  if(TMath::Abs(pdg)==211)  eff = LinearInterpolation(9,effPi,pt,eta);
  if(TMath::Abs(pdg)==321)  eff = LinearInterpolation(6,effK,pt,eta);
  if(TMath::Abs(pdg)==2212) eff = LinearInterpolation(5,effP,pt,eta);
  if(TMath::Abs(pdg)==11)   eff = LinearInterpolation(5,effEl,pt,eta);
  if(TMath::Abs(pdg)==13)   eff = LinearInterpolation(5,effMu,pt,eta);

  if(gRandom->Rndm() < eff) return kTRUE;

  return kFALSE;
}

//-----------------------------------------------------------------
Double_t AliTPCtrackerParam::LinearInterpolation(Int_t ptBins,Double_t *value,
						 Double_t pt,Double_t eta) 
{
  Double_t kValue=0,kValue1=0,kValue2=0;
  Int_t kEtaSide = (TMath::Abs(eta)<.45 ? 0 : 1);
  Double_t kEta[3]={0.15,0.45,0.75};
  Double_t kPt[9]={0.244,0.390,0.676,1.190,2.36,4.,6.,10.,20.};
  if(ptBins==6) kPt[5]=10.;

  for(Int_t i=0; i<ptBins; i++) {
    if(pt<kPt[i]) {
      if(i==0) i=1;
      kValue1 = value[3*i-3+kEtaSide]+(value[3*i-3+kEtaSide]-value[3*i+kEtaSide])/(kPt[i-1]-kPt[i])*(pt-kPt[i-1]);
      kValue2 = value[3*i-3+kEtaSide+1]+(value[3*i-3+kEtaSide+1]-value[3*i+kEtaSide+1])/(kPt[i-1]-kPt[i])*(pt-kPt[i-1]);

      kValue = kValue1+(kValue1-kValue2)/(kEta[kEtaSide]-kEta[kEtaSide+1])*(eta-kEta[kEtaSide]);
      break;
    }
    if(i==ptBins-1) { 
      i=ptBins-2;
      kValue1 = value[3*i-3+kEtaSide]+(value[3*i-3+kEtaSide]-value[3*i+kEtaSide])/(kPt[i-1]-kPt[i])*(pt-kPt[i-1]);
      kValue2 = value[3*i-3+kEtaSide+1]+(value[3*i-3+kEtaSide+1]-value[3*i+kEtaSide+1])/(kPt[i-1]-kPt[i])*(pt-kPt[i-1]);

      kValue = kValue1+(kValue1-kValue2)/(kEta[kEtaSide]-kEta[kEtaSide+1])*(eta-kEta[kEtaSide]);
      break;
    }
  }

  return kValue;
}

//-----------------------------------------------------------------
Int_t AliTPCtrackerParam::GetBin(Double_t pt,Double_t eta) {

  if(TMath::Abs(eta)<0.3) {
    if(pt<0.3)            return 0;
    if(pt>=0.3 && pt<0.5) return 3;
    if(pt>=0.5 && pt<1.)  return 6;
    if(pt>=1. && pt<1.5)  return 9;
    if(pt>=1.5 && pt<3.)  return 12;
    if(pt>=3. && pt<5.)   return 15;
    if(pt>=5. && pt<7.)   return 18;
    if(pt>=7. && pt<15.)  return 21;
    if(pt>=15.)           return 24;
  }
  if(TMath::Abs(eta)>=0.3 && TMath::Abs(eta)<0.6) {
    if(pt<0.3)            return 1;
    if(pt>=0.3 && pt<0.5) return 4;
    if(pt>=0.5 && pt<1.)  return 7;
    if(pt>=1. && pt<1.5)  return 10;
    if(pt>=1.5 && pt<3.)  return 13;
    if(pt>=3. && pt<5.)   return 16;
    if(pt>=5. && pt<7.)   return 19;
    if(pt>=7. && pt<15.)  return 22;
    if(pt>=15.)           return 25;
  }
  if(TMath::Abs(eta)>=0.6) {
    if(pt<0.3)            return 2;
    if(pt>=0.3 && pt<0.5) return 5;
    if(pt>=0.5 && pt<1.)  return 8;
    if(pt>=1. && pt<1.5)  return 11;
    if(pt>=1.5 && pt<3.)  return 14;
    if(pt>=3. && pt<5.)   return 17;
    if(pt>=5. && pt<7.)   return 20;
    if(pt>=7. && pt<15.)  return 23;
    if(pt>=15.)           return 26;
  }

  return -1;

}

//-----------------------------------------------------------------

TMatrixD AliTPCtrackerParam::GetSmearingMatrix(Double_t* cc,Double_t pt,Double_t eta) 
{
  TMatrixD CovMat(5,5);

  CovMat(0,0)=cc[0];
  CovMat(1,0)=cc[1];  CovMat(0,1)=CovMat(1,0);
  CovMat(1,1)=cc[2];
  CovMat(2,0)=cc[3];  CovMat(0,2)=CovMat(2,0);
  CovMat(2,1)=cc[4];  CovMat(1,2)=CovMat(2,1);
  CovMat(2,2)=cc[5];
  CovMat(3,0)=cc[6];  CovMat(0,3)=CovMat(3,0);
  CovMat(3,1)=cc[7];  CovMat(1,3)=CovMat(3,1);
  CovMat(3,2)=cc[8];  CovMat(2,3)=CovMat(3,2);
  CovMat(3,3)=cc[9];
  CovMat(4,0)=cc[10]; CovMat(0,4)=CovMat(4,0);
  CovMat(4,1)=cc[11]; CovMat(1,4)=CovMat(4,1);
  CovMat(4,2)=cc[12]; CovMat(2,4)=CovMat(4,2);
  CovMat(4,3)=cc[13]; CovMat(3,4)=CovMat(4,3);
  CovMat(4,4)=cc[14];

  TMatrixD ScaleMat(5,5);
  for(Int_t k=0;k<5;k++) {
    for(Int_t l=0;l<5;l++) {
      ScaleMat(k,l)=0.;
    }
  }

  
  Double_t s00[27]={1.3553,1.2973,1.2446,1.3428,1.3031,1.2345,1.3244,1.2681,1.2046,1.3046,1.2430,1.1606,1.2462,1.2104,1.1082,1.2207,1.1189,1.0789,1.2079,1.1300,1.0426,1.1502,1.1122,1.0559,1.1419,1.1072,1.0208};
  Double_t s11[27]={1.0890,1.1463,1.2313,1.1149,1.1597,1.2280,1.1225,1.1584,1.2329,1.1149,1.1550,1.2369,1.1095,1.1353,1.2050,1.0649,1.0858,1.1546,1.0663,1.0672,1.1340,1.0416,1.0738,1.0945,1.0663,1.0654,1.0909};
  Double_t s22[27]={1.1709,1.1367,1.0932,1.2247,1.1832,1.1247,1.2470,1.2017,1.1441,1.2202,1.1653,1.1050,1.1819,1.1269,1.0583,1.1546,1.0621,0.9928,1.1305,1.0512,0.9576,1.0995,1.0445,0.9884,1.0968,1.0368,0.9574};
  Double_t s33[27]={0.9720,0.9869,1.0271,1.0030,1.0223,1.0479,1.0164,1.0305,1.0559,1.0339,1.0450,1.0686,1.0450,1.0507,1.0784,1.0334,1.0208,1.0863,1.0252,1.0114,1.0835,0.9854,1.0144,1.0507,1.0124,1.0159,1.0464};
  Double_t s44[27]={1.1104,1.0789,1.0479,1.1597,1.1234,1.0728,1.2087,1.1602,1.1041,1.1942,1.1428,1.0831,1.1572,1.1036,1.0431,1.1296,1.0498,0.9844,1.1145,1.0266,0.9489,1.0959,1.0450,0.9875,1.0775,1.0266,0.9406};


  ScaleMat(0,0) = LinearInterpolation(9,s00,pt,eta); 
  ScaleMat(1,1) = LinearInterpolation(9,s11,pt,eta); 
  ScaleMat(2,2) = LinearInterpolation(9,s22,pt,eta); 
  ScaleMat(3,3) = LinearInterpolation(9,s33,pt,eta); 
  ScaleMat(4,4) = LinearInterpolation(9,s44,pt,eta); 

  TMatrixD Mat(ScaleMat,TMatrixD::kMult,CovMat);
  TMatrixD CovMatSmear(Mat,TMatrixD::kMult,ScaleMat);


  return CovMatSmear;
}

//-----------------------------------------------------------------
void AliTPCtrackerParam::SmearTrack(Double_t* xx,Double_t* xxsm,TMatrixD cov) {

  // build triangular matrix Bmat
  TMatrixD* Bmat = new TMatrixD(5,5);

  for(Int_t j=0;j<5;j++){
    Double_t accum = 0;
    for(Int_t k=0;k<j;k++){
      accum += (*Bmat)(j,k)* (*Bmat)(j,k);
    }
    (*Bmat)(j,j)=TMath::Sqrt(TMath::Abs(cov(j,j)-accum));
    for(Int_t i=j+1;i<5;i++){
      accum = 0;
      for(Int_t k=0;k<j;k++){
	accum+=(*Bmat)(i,k)* (*Bmat)(j,k);
      }
      (*Bmat)(i,j) = (cov(i,j)-accum) / (*Bmat)(j,j);
    }
  }

  // get array of numbers with normal distribution
  TArrayD norm;
  norm.Set(5);
  for(Int_t l=0;l<5;l++) {
    norm[l]=gRandom->Gaus();
  }

  // use Bmat matrix to generate correlated numbers
  TArrayD corr;
  corr.Set(5);
  for(Int_t i=0;i<5;i++){
    corr[i]=0;
    for(Int_t j=0;j<=i;j++) corr[i] += (*Bmat)(i,j)*norm[j];
  }

  delete Bmat;

  for(Int_t l=0;l<5;l++) {
    xxsm[l] = xx[l]+corr[l];
  }

  return;
}

//-----------------------------------------------------------------
void AliTPCtrackerParam::CookTracks(TObjArray& tarray,TObjArray& newtarray) 
{
  TString* s = new TString("$ALICE_ROOT/TPC/CovMatrixDB_");
  if(fColl==0) s->Append("PbPb6000");
  if(fBz==0.4) s->Append("_B0.4T.root");  
  // open file with matrixes DB
  TFile* DBfile = new TFile(s->Data());  

  AliTPCtrack* track=new AliTPCtrack;
  Int_t entr = (Int_t)tarray.GetEntriesFast();

  for(Int_t k=0; k<entr; k++) {
    track=(AliTPCtrack*)tarray.UncheckedAt(k);  

    Double_t pt = 1/TMath::Abs(track->Get1Pt());
    Double_t eta = -TMath::Log(TMath::Tan(0.25*TMath::Pi()-0.5*TMath::ATan(track->GetTgl())));

    Int_t bin = GetBin(pt,eta);

    // string with tree name
    TString str("Tree_bin");
    str+=bin;
 
    // get the right tree
    TTree* CovTree = (TTree*)DBfile->Get(str.Data());
    TArrayF* matrix = new TArrayF; 
    CovTree->SetBranchAddress("matrixes",&matrix);

    // get random entry from the tree
    Int_t Entries = (Int_t)CovTree->GetEntries();
    CovTree->GetEvent(gRandom->Integer(Entries));

    Double_t xref,alpha,xx[5],xxsm[5],cc[15];
    Int_t lab;
    // get P and Cosl from track
    Double_t Cosl = TMath::Cos(TMath::ATan(track->GetTgl()));
    Double_t P = 1/TMath::Abs(track->Get1Pt())/Cosl;
    
    // get covariance matrix from regularized matrix
    cc[0]=matrix->At(0)*(1.588e-3+1.911e-4/TMath::Power(P,1.5));
    cc[1]=matrix->At(1);
    cc[2]=matrix->At(2)*(1.195e-3+0.8102e-3/P);
    cc[3]=matrix->At(3)*(7.275e-5+1.181e-5/TMath::Power(P,1.5));
    cc[4]=matrix->At(4);
    cc[5]=matrix->At(5)*(5.163e-6+1.138e-6/TMath::Power(P,2.)/Cosl);
    cc[6]=matrix->At(6);
    cc[7]=matrix->At(7)*(1.176e-5+1.175e-5/TMath::Power(P,1.5)/Cosl/Cosl/Cosl);
    cc[8]=matrix->At(8);
    cc[9]=matrix->At(9)*(1.289e-7+4.636e-7/TMath::Power(P,1.7)/Cosl/Cosl/Cosl/Cosl);
    cc[10]=matrix->At(10)*(4.056e-7+4.379e-8/TMath::Power(P,1.5));
    cc[11]=matrix->At(11);
    cc[12]=matrix->At(12)*(3.049e-8+8.244e-9/TMath::Power(P,2.)/Cosl);
    cc[13]=matrix->At(13);
    cc[14]=matrix->At(14)*(1.847e-10+5.822e-11/TMath::Power(P,2.)/Cosl/Cosl/Cosl);
  
    TMatrixD CovMatSmear(5,5);
    
    CovMatSmear = GetSmearingMatrix(cc,pt,eta);

    // get track original parameters
    xref=track->GetX();
    alpha=track->GetAlpha();
    xx[0]=track->GetY();
    xx[1]=track->GetZ();
    xx[2]=track->GetX()*track->GetC()-track->GetSnp();
    xx[3]=track->GetTgl();
    xx[4]=track->GetC();
    lab=track->GetLabel();

    // use smearing matrix to smear the original parameters
    SmearTrack(xx,xxsm,CovMatSmear);

    AliTPCtrack* tpctrack = new AliTPCtrack(0,xxsm,cc,xref,alpha);
    tpctrack->SetLabel(lab);

    // fill the array
    newtarray.AddLast(tpctrack);
  
  }

  DBfile->Close();

  return; 
}
//-----------------------------------------------------------------












