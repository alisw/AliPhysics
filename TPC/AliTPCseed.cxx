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




//-----------------------------------------------------------------
//           Implementation of the TPC seed class
//        This class is used by the AliTPCtrackerMI class
//      Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch
//-----------------------------------------------------------------
#include "TClonesArray.h"
#include "AliTPCseed.h"

ClassImp(AliTPCseed)



AliTPCseed::AliTPCseed():AliTPCtrack(){
  //
  fRow=0; 
  fRemoval =0; 
  for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
  for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;
  for (Int_t i=0;i<3;i++)   fKinkIndexes[i]=0;
  for (Int_t i=0;i<5;i++)   fTPCr[i]=0.2;
  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  =0;
  fRemoval = 0;
  fSort =0;
  fFirstPoint =0;
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
  fEsd =0;
  fCircular = 0;  // not curling track
}
AliTPCseed::AliTPCseed(const AliTPCseed &s):AliTPCtrack(s){
  //---------------------
  // dummy copy constructor
  //-------------------------
  for (Int_t i=0;i<160;i++) fClusterPointer[i] = s.fClusterPointer[i];
  for (Int_t i=0;i<160;i++) fIndex[i] = s.fIndex[i];

  fPoints  = 0;
  fEPoints = 0;
  fCircular =0;
  fEsd =0;
}
AliTPCseed::AliTPCseed(const AliTPCtrack &t):AliTPCtrack(t){
  //
  //copy constructor
  fPoints = 0;
  fEPoints = 0;
  fNShared  =0; 
  //  fTrackPoints =0;
  fRemoval =0;
  fSort =0;
  for (Int_t i=0;i<3;i++)   fKinkIndexes[i]=t.GetKinkIndex(i);
  for (Int_t i=0;i<5;i++)   fTPCr[i]=0.2;
  for (Int_t i=0;i<160;i++) {
    fClusterPointer[i] = 0;
    Int_t index = t.GetClusterIndex(i);
    if (index>=-1){ 
      SetClusterIndex2(i,index);
    }
    else{
      SetClusterIndex2(i,-3); 
    }    
  }
  fFirstPoint =0;
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
  fCircular =0;
  fEsd =0;
}

AliTPCseed::AliTPCseed(UInt_t index,  const Double_t xx[5], const Double_t cc[15], 
                                        Double_t xr, Double_t alpha):      
  AliTPCtrack(index, xx, cc, xr, alpha) {
   //
  //
  //constructor
  fRow =0;
  for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
  for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;
  for (Int_t i=0;i<3;i++)   fKinkIndexes[i]=0;
  for (Int_t i=0;i<5;i++)   fTPCr[i]=0.2;

  fPoints = 0;
  fEPoints = 0;
  fNFoundable =0;
  fNShared  = 0;
  //  fTrackPoints =0;
  fRemoval =0;
  fSort =0;
  fFirstPoint =0;
  //  fHelixIn = new TClonesArray("AliHelix",0);
  //fHelixOut = new TClonesArray("AliHelix",0);
  fNoCluster =0;
  fBSigned = kFALSE;
  fSeed1 =-1;
  fSeed2 =-1;
  fCurrentCluster =0;
  fCurrentSigmaY2=0;
  fCurrentSigmaZ2=0;
  fEsd =0;
}

AliTPCseed::~AliTPCseed(){
  //
  // destructor
  if (fPoints) delete fPoints;
  fPoints =0;
  if (fEPoints) delete fEPoints;
  fEPoints = 0;
  fNoCluster =0;
}

AliTPCTrackerPoint * AliTPCseed::GetTrackPoint(Int_t i)
{
  //
  // 
  return &fTrackPoints[i];
}

void AliTPCseed::RebuildSeed()
{
  //
  // rebuild seed to be ready for storing
  AliTPCclusterMI cldummy;
  cldummy.SetQ(0);
  AliTPCTrackPoint pdummy;
  pdummy.GetTPoint().fIsShared = 10;
  for (Int_t i=0;i<160;i++){
    AliTPCclusterMI * cl0 = fClusterPointer[i];
    AliTPCTrackPoint *trpoint = (AliTPCTrackPoint*)fPoints->UncheckedAt(i);     
    if (cl0){
      trpoint->GetTPoint() = *(GetTrackPoint(i));
      trpoint->GetCPoint() = *cl0;
      trpoint->GetCPoint().SetQ(TMath::Abs(cl0->GetQ()));
    }
    else{
      *trpoint = pdummy;
      trpoint->GetCPoint()= cldummy;
    }
    
  }

}


Double_t AliTPCseed::GetDensityFirst(Int_t n)
{
  //
  //
  // return cluster for n rows bellow first point
  Int_t nfoundable = 1;
  Int_t nfound      = 1;
  for (Int_t i=fLastPoint-1;i>0&&nfoundable<n; i--){
    Int_t index = GetClusterIndex2(i);
    if (index!=-1) nfoundable++;
    if (index>0) nfound++;
  }
  if (nfoundable<n) return 0;
  return Double_t(nfound)/Double_t(nfoundable);

}


void AliTPCseed::GetClusterStatistic(Int_t first, Int_t last, Int_t &found, Int_t &foundable, Int_t &shared, Bool_t plus2)
{
  // get cluster stat.  on given region
  //
  found       = 0;
  foundable   = 0;
  shared      =0;
  for (Int_t i=first;i<last; i++){
    Int_t index = GetClusterIndex2(i);
    if (index!=-1) foundable++;
    if (fClusterPointer[i]) {
      found++;
    }
    else 
      continue;

    if (fClusterPointer[i]->IsUsed(10)) {
      shared++;
      continue;
    }
    if (!plus2) continue; //take also neighborhoud
    //
    if ( (i>0) && fClusterPointer[i-1]){
      if (fClusterPointer[i-1]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    if ( fClusterPointer[i+1]){
      if (fClusterPointer[i+1]->IsUsed(10)) {
	shared++;
	continue;
      }
    }
    
  }
  //if (shared>found){
    //Error("AliTPCseed::GetClusterStatistic","problem\n");
  //}
}





void AliTPCseed::Reset(Bool_t all)
{
  //
  //
  SetNumberOfClusters(0);
  fNFoundable = 0;
  SetChi2(0);
  ResetCovariance();
  /*
  if (fTrackPoints){
    for (Int_t i=0;i<8;i++){
      delete [] fTrackPoints[i];
    }
    delete fTrackPoints;
    fTrackPoints =0;
  }
  */

  if (all){   
    for (Int_t i=0;i<200;i++) SetClusterIndex2(i,-3);
    for (Int_t i=0;i<160;i++) fClusterPointer[i]=0;
  }

}


void AliTPCseed::Modify(Double_t factor)
{

  //------------------------------------------------------------------
  //This function makes a track forget its history :)  
  //------------------------------------------------------------------
  if (factor<=0) {
    ResetCovariance();
    return;
  }
  fC00*=factor;
  fC10*=0;  fC11*=factor;
  fC20*=0;  fC21*=0;  fC22*=factor;
  fC30*=0;  fC31*=0;  fC32*=0;  fC33*=factor;
  fC40*=0;  fC41*=0;  fC42*=0;  fC43*=0;  fC44*=factor;
  SetNumberOfClusters(0);
  fNFoundable =0;
  SetChi2(0);
  fRemoval = 0;
  fCurrentSigmaY2 = 0.000005;
  fCurrentSigmaZ2 = 0.000005;
  fNoCluster     = 0;
  //fFirstPoint = 160;
  //fLastPoint  = 0;
}




Int_t  AliTPCseed::GetProlongation(Double_t xk, Double_t &y, Double_t & z) const
{
  //-----------------------------------------------------------------
  // This function find proloncation of a track to a reference plane x=xk.
  // doesn't change internal state of the track
  //-----------------------------------------------------------------
  
  Double_t x1=fX, x2=x1+(xk-x1), dx=x2-x1;

  if (TMath::Abs(fP4*xk - fP2) >= 0.999) {   
    return 0;
  }

  //  Double_t y1=fP0, z1=fP1;
  Double_t c1=fP4*x1 - fP2, r1=sqrt(1.- c1*c1);
  Double_t c2=fP4*x2 - fP2, r2=sqrt(1.- c2*c2);
  
  y = fP0;
  z = fP1;
  //y += dx*(c1+c2)/(r1+r2);
  //z += dx*(c1+c2)/(c1*r2 + c2*r1)*fP3;
  
  Double_t dy = dx*(c1+c2)/(r1+r2);
  Double_t dz = 0;
  //
  Double_t delta = fP4*dx*(c1+c2)/(c1*r2 + c2*r1);
  /*
  if (TMath::Abs(delta)>0.0001){
    dz = fP3*TMath::ASin(delta)/fP4;
  }else{
    dz = dx*fP3*(c1+c2)/(c1*r2 + c2*r1);
  }
  */
  //  dz =  fP3*AliTPCFastMath::FastAsin(delta)/fP4;
  dz =  fP3*TMath::ASin(delta)/fP4;
  //
  y+=dy;
  z+=dz;
  

  return 1;  
}


//_____________________________________________________________________________
Double_t AliTPCseed::GetPredictedChi2(const AliCluster *c) const 
{
  //-----------------------------------------------------------------
  // This function calculates a predicted chi2 increment.
  //-----------------------------------------------------------------
  //Double_t r00=c->GetSigmaY2(), r01=0., r11=c->GetSigmaZ2();
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;
  r00+=fC00; r01+=fC10; r11+=fC11;

  Double_t det=r00*r11 - r01*r01;
  if (TMath::Abs(det) < 1.e-10) {
    //Int_t n=GetNumberOfClusters();
    //if (n>4) cerr<<n<<" AliKalmanTrack warning: Singular matrix !\n";
    return 1e10;
  }
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  
  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  
  return (dy*r00*dy + 2*r01*dy*dz + dz*r11*dz)/det;
}


//_________________________________________________________________________________________


Int_t AliTPCseed::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the sector - for given sector according z
  //-----------------------------------------------------------------
  AliTPCseed *t=(AliTPCseed*)o;

  if (fSort == 0){
    if (t->fRelativeSector>fRelativeSector) return -1;
    if (t->fRelativeSector<fRelativeSector) return 1;
    Double_t z2 = t->GetZ();
    Double_t z1 = GetZ();
    if (z2>z1) return 1;
    if (z2<z1) return -1;
    return 0;
  }
  else {
    Float_t f2 =1;
    f2 = 1-20*TMath::Sqrt(t->fC44)/(TMath::Abs(t->GetC())+0.0066);
    if (t->fBConstrain) f2=1.2;

    Float_t f1 =1;
    f1 = 1-20*TMath::Sqrt(fC44)/(TMath::Abs(GetC())+0.0066);

    if (fBConstrain)   f1=1.2;
 
    if (t->GetNumberOfClusters()*f2 <GetNumberOfClusters()*f1) return -1;
    else return +1;
  }
}




//_____________________________________________________________________________
Int_t AliTPCseed::Update(const AliCluster *c, Double_t chisq, UInt_t /*index*/) {
  //-----------------------------------------------------------------
  // This function associates a cluster with this track.
  //-----------------------------------------------------------------
  Double_t r00=fErrorY2, r01=0., r11=fErrorZ2;

  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11 - r01*r01;
  Double_t tmp=r00; r00=r11/det; r11=tmp/det; r01=-r01/det;

  Double_t k00=fC00*r00+fC10*r01, k01=fC00*r01+fC10*r11;
  Double_t k10=fC10*r00+fC11*r01, k11=fC10*r01+fC11*r11;
  Double_t k20=fC20*r00+fC21*r01, k21=fC20*r01+fC21*r11;
  Double_t k30=fC30*r00+fC31*r01, k31=fC30*r01+fC31*r11;
  Double_t k40=fC40*r00+fC41*r01, k41=fC40*r01+fC41*r11;

  Double_t dy=c->GetY() - fP0, dz=c->GetZ() - fP1;
  Double_t cur=fP4 + k40*dy + k41*dz, eta=fP2 + k20*dy + k21*dz;
  if (TMath::Abs(cur*fX-eta) >= 0.9) {
    return 0;
  }

  fP0 += k00*dy + k01*dz;
  fP1 += k10*dy + k11*dz;
  fP2  = eta;
  fP3 += k30*dy + k31*dz;
  fP4  = cur;

  Double_t c01=fC10, c02=fC20, c03=fC30, c04=fC40;
  Double_t c12=fC21, c13=fC31, c14=fC41;

  fC00-=k00*fC00+k01*fC10; fC10-=k00*c01+k01*fC11;
  fC20-=k00*c02+k01*c12;   fC30-=k00*c03+k01*c13;
  fC40-=k00*c04+k01*c14; 

  fC11-=k10*c01+k11*fC11;
  fC21-=k10*c02+k11*c12;   fC31-=k10*c03+k11*c13;
  fC41-=k10*c04+k11*c14; 

  fC22-=k20*c02+k21*c12;   fC32-=k20*c03+k21*c13;
  fC42-=k20*c04+k21*c14; 

  fC33-=k30*c03+k31*c13;
  fC43-=k40*c03+k41*c13; 

  fC44-=k40*c04+k41*c14; 

  Int_t n=GetNumberOfClusters();
  //  fIndex[n]=index;
  SetNumberOfClusters(n+1);
  SetChi2(GetChi2()+chisq);

  return 1;
}



//_____________________________________________________________________________
Float_t AliTPCseed::CookdEdx(Double_t low, Double_t up,Int_t i1, Int_t i2, Bool_t onlyused) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[200];
  Float_t angular[200];
  Float_t weight[200];
  Int_t index[200];
  //Int_t nc = 0;
  //  TClonesArray & arr = *fPoints; 
  Float_t meanlog = 100.;
  
  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  for (Int_t of =0; of<4; of++){    
    for (Int_t i=of+i1;i<i2;i+=4)
      {
	Int_t index = fIndex[i];
	if (index<0||index&0x8000) continue;

	//AliTPCTrackPoint * point = (AliTPCTrackPoint *) arr.At(i);
	AliTPCTrackerPoint * point = GetTrackPoint(i);
	//AliTPCTrackerPoint * pointm = GetTrackPoint(i-1);
	//AliTPCTrackerPoint * pointp = 0;
	//if (i<159) pointp = GetTrackPoint(i+1);

	if (point==0) continue;
	AliTPCclusterMI * cl = fClusterPointer[i];
	if (cl==0) continue;	
	if (onlyused && (!cl->IsUsed(10))) continue;
	if (cl->IsUsed(11)) {
	  fNShared++;
	  continue;
	}
	Int_t   type   = cl->GetType();
	//if (point->fIsShared){
	//  fNShared++;
	//  continue;
	//}
	//if (pointm) 
	//  if (pointm->fIsShared) continue;
	//if (pointp) 
	//  if (pointp->fIsShared) continue;

	if (type<0) continue;
	//if (type>10) continue;       
	//if (point->GetErrY()==0) continue;
	//if (point->GetErrZ()==0) continue;

	//Float_t ddy = (point->GetY()-cl->GetY())/point->GetErrY();
	//Float_t ddz = (point->GetZ()-cl->GetZ())/point->GetErrZ();
	//if ((ddy*ddy+ddz*ddz)>10) continue; 


	//	if (point->GetCPoint().GetMax()<5) continue;
	if (cl->GetMax()<5) continue;
	Float_t angley = point->GetAngleY();
	Float_t anglez = point->GetAngleZ();

	Float_t rsigmay2 =  point->GetSigmaY();
	Float_t rsigmaz2 =  point->GetSigmaZ();
	/*
	Float_t ns = 1.;
	if (pointm){
	  rsigmay +=  pointm->GetTPoint().GetSigmaY();
	  rsigmaz +=  pointm->GetTPoint().GetSigmaZ();
	  ns+=1.;
	}
	if (pointp){
	  rsigmay +=  pointp->GetTPoint().GetSigmaY();
	  rsigmaz +=  pointp->GetTPoint().GetSigmaZ();
	  ns+=1.;
	}
	rsigmay/=ns;
	rsigmaz/=ns;
	*/

	Float_t rsigma = TMath::Sqrt(rsigmay2*rsigmaz2);

	Float_t ampc   = 0;     // normalization to the number of electrons
	if (i>64){
	  //	  ampc = 1.*point->GetCPoint().GetMax();
	  ampc = 1.*cl->GetMax();
	  //ampc = 1.*point->GetCPoint().GetQ();	  
	  //	  AliTPCClusterPoint & p = point->GetCPoint();
	  //	  Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.6)) - TMath::Abs(p.GetY()/0.6)+0.5);
	  // Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);
	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	  Float_t zfactor = (AliTPCReconstructor::GetCtgRange()-0.0004*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 
	}
	else{ 
	  //ampc = 1.0*point->GetCPoint().GetMax(); 
	  ampc = 1.0*cl->GetMax(); 
	  //ampc = 1.0*point->GetCPoint().GetQ(); 
	  //AliTPCClusterPoint & p = point->GetCPoint();
	  // Float_t dy = TMath::Abs(Int_t( TMath::Abs(p.GetY()/0.4)) - TMath::Abs(p.GetY()/0.4)+0.5);
	  //Float_t iz =  (250.0-TMath::Abs(p.GetZ())+0.11)/0.566;
	  //Float_t dz = 
	  //  TMath::Abs( Int_t(iz) - iz + 0.5);

	  //ampc *= 1.15*(1-0.3*dy);
	  //ampc *= 1.15*(1-0.3*dz);
	  //	Float_t zfactor = (1.02-0.000*TMath::Abs(point->GetCPoint().GetZ()));
	  //ampc               *=zfactor; 

	}
	ampc *= 2.0;     // put mean value to channel 50
	//ampc *= 0.58;     // put mean value to channel 50
	Float_t w      =  1.;
	//	if (type>0)  w =  1./(type/2.-0.5); 
	//	Float_t z = TMath::Abs(cl->GetZ());
	if (i<64) {
	  ampc /= 0.6;
	  //ampc /= (1+0.0008*z);
	} else
	  if (i>128){
	    ampc /=1.5;
	    //ampc /= (1+0.0008*z);
	  }else{
	    //ampc /= (1+0.0008*z);
	  }
	
	if (type<0) {  //amp at the border - lower weight
	  // w*= 2.;
	  
	  continue;
	}
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	amp[nc[of]]        = ampc;
	angular[nc[of]]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	weight[nc[of]]     = w;
	nc[of]++;
      }
    
    TMath::Sort(nc[of],amp,index,kFALSE);
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;
    //meanlog = amp[index[Int_t(nc[of]*0.33)]];
    meanlog = 50;
    for (Int_t i=int(nc[of]*low+0.5);i<int(nc[of]*up+0.5);i++){
      Float_t ampl      = amp[index[i]]/angular[index[i]];
      ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
      //
      sumw    += weight[index[i]]; 
      sumamp  += weight[index[i]]*ampl;
      sumamp2 += weight[index[i]]*ampl*ampl;
      norm[of]    += angular[index[i]]*weight[index[i]];
    }
    if (sumw<1){ 
      SetdEdx(0);  
    }
    else {
      norm[of] /= sumw;
      mean[of]  = sumamp/sumw;
      sigma[of] = sumamp2/sumw-mean[of]*mean[of];
      if (sigma[of]>0.1) 
	sigma[of] = TMath::Sqrt(sigma[of]);
      else
	sigma[of] = 1000;
      
    mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
    //mean  *=(1-0.02*(sigma/(mean*0.17)-1.));
    //mean *=(1-0.1*(norm-1.));
    }
  }

  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //  mean[0]*= (1-0.05*(sigma[0]/(0.01+mean[1]*0.18)-1));
  //  mean[1]*= (1-0.05*(sigma[1]/(0.01+mean[0]*0.18)-1));

  
  //  dedx = (mean[0]* TMath::Sqrt((1.+nc[0]))+ mean[1]* TMath::Sqrt((1.+nc[1])) )/ 
  //  (  TMath::Sqrt((1.+nc[0]))+TMath::Sqrt((1.+nc[1])));

  Int_t norm2 = 0;
  Int_t norm3 = 0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i];
      fSdEdx    += sigma[i]*(nc[i]-2);
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i];
      norm3     += nc[i]-2;
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return 0;
  }
  //  Float_t dedx1 =dedx;
  /*
  dedx =0;
  for (Int_t i =0;i<4;i++){
    if (nc[i]>2&&nc[i]<1000){
      mean[i]   = mean[i]*(1-0.12*(sigma[i]/(fSdEdx)-1.));
      dedx      += mean[i] *nc[i];
    }
    fDEDX[i]  = mean[i];                
  }
  dedx /= norm2;
  */

  
  SetdEdx(dedx);
  return dedx;
}
Double_t AliTPCseed::Bethe(Double_t bg){
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  //
  Double_t bg2=bg*bg;
  Double_t bethe;
  if (bg<3.5e1) 
    bethe=(1.+ bg2)/bg2*(log(5940*bg2) - bg2/(1.+ bg2));
  else // Density effect ( approximately :) 
    bethe=1.15*(1.+ bg2)/bg2*(log(3.5*5940*bg) - bg2/(1.+ bg2));
  return bethe/11.091;
}

void AliTPCseed::CookPID()
{
  //
  // cook PID information according dEdx
  //
  Double_t fRange = 10.;
  Double_t fRes   = 0.1;
  Double_t fMIP   = 47.;
  //
  Int_t ns=AliPID::kSPECIES;
  Double_t sumr =0;
  for (Int_t j=0; j<ns; j++) {
    Double_t mass=AliPID::ParticleMass(j);
    Double_t mom=P();
    Double_t dedx=fdEdx/fMIP;
    Double_t bethe=Bethe(mom/mass); 
    Double_t sigma=fRes*bethe;
    if (sigma>0.001){
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	fTPCr[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
	sumr+=fTPCr[j];
	continue;
      }
      fTPCr[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      sumr+=fTPCr[j];
    }
    else{
      fTPCr[j]=1.;
      sumr+=fTPCr[j];
    }
  }
  for (Int_t j=0; j<ns; j++) {
    fTPCr[j]/=sumr;           //normalize
  }
}

/*
void AliTPCseed::CookdEdx2(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------

  Float_t amp[200];
  Float_t angular[200];
  Float_t weight[200];
  Int_t index[200];
  Bool_t inlimit[200];
  for (Int_t i=0;i<200;i++) inlimit[i]=kFALSE;
  for (Int_t i=0;i<200;i++) amp[i]=10000;
  for (Int_t i=0;i<200;i++) angular[i]= 1;;
  

  //
  Float_t meanlog = 100.;
  Int_t indexde[4]={0,64,128,160};

  Float_t amean     =0;
  Float_t asigma    =0;
  Float_t anc       =0;
  Float_t anorm     =0;

  Float_t mean[4]  = {0,0,0,0};
  Float_t sigma[4] = {1000,1000,1000,1000};
  Int_t nc[4]      = {0,0,0,0};
  Float_t norm[4]    = {1000,1000,1000,1000};
  //
  //
  fNShared =0;

  //  for (Int_t of =0; of<3; of++){    
  //  for (Int_t i=indexde[of];i<indexde[of+1];i++)
  for (Int_t i =0; i<160;i++)
    {
	AliTPCTrackPoint * point = GetTrackPoint(i);
	if (point==0) continue;
	if (point->fIsShared){
	  fNShared++;	  
	  continue;
	}
	Int_t   type   = point->GetCPoint().GetType();
	if (type<0) continue;
	if (point->GetCPoint().GetMax()<5) continue;
	Float_t angley = point->GetTPoint().GetAngleY();
	Float_t anglez = point->GetTPoint().GetAngleZ();
	Float_t rsigmay =  point->GetCPoint().GetSigmaY();
	Float_t rsigmaz =  point->GetCPoint().GetSigmaZ();
	Float_t rsigma = TMath::Sqrt(rsigmay*rsigmaz);

	Float_t ampc   = 0;     // normalization to the number of electrons
	if (i>64){
	  ampc =  point->GetCPoint().GetMax();
	}
	else{ 
	  ampc = point->GetCPoint().GetMax(); 
	}
	ampc *= 2.0;     // put mean value to channel 50
	//	ampc *= 0.565;     // put mean value to channel 50

	Float_t w      =  1.;
	Float_t z = TMath::Abs(point->GetCPoint().GetZ());
	if (i<64) {
	  ampc /= 0.63;
	} else
	  if (i>128){
	    ampc /=1.51;
	  }	  	
	if (type<0) {  //amp at the border - lower weight	  	  
	  continue;
	}
	if (rsigma>1.5) ampc/=1.3;  // if big backround
	angular[i]    = TMath::Sqrt(1.+angley*angley+anglez*anglez);
	amp[i]        = ampc/angular[i];
	weight[i]     = w;
	anc++;
    }

  TMath::Sort(159,amp,index,kFALSE);
  for (Int_t i=int(anc*low+0.5);i<int(anc*up+0.5);i++){      
    inlimit[index[i]] = kTRUE;  // take all clusters
  }
  
  //  meanlog = amp[index[Int_t(anc*0.3)]];
  meanlog =10000.;
  for (Int_t of =0; of<3; of++){    
    Float_t sumamp=0;
    Float_t sumamp2=0;
    Float_t sumw=0;    
   for (Int_t i=indexde[of];i<indexde[of+1];i++)
      {
	if (inlimit[i]==kFALSE) continue;
	Float_t ampl      = amp[i];
	///angular[i];
	ampl              = meanlog*TMath::Log(1.+ampl/meanlog);
	//
	sumw    += weight[i]; 
	sumamp  += weight[i]*ampl;
	sumamp2 += weight[i]*ampl*ampl;
	norm[of]    += angular[i]*weight[i];
	nc[of]++;
      }
   if (sumw<1){ 
     SetdEdx(0);  
   }
   else {
     norm[of] /= sumw;
     mean[of]  = sumamp/sumw;
     sigma[of] = sumamp2/sumw-mean[of]*mean[of];
     if (sigma[of]>0.1) 
       sigma[of] = TMath::Sqrt(sigma[of]);
     else
       sigma[of] = 1000;      
     mean[of] = (TMath::Exp(mean[of]/meanlog)-1)*meanlog;
   }
  }
    
  Float_t dedx =0;
  fSdEdx =0;
  fMAngular =0;
  //
  Int_t norm2 = 0;
  Int_t norm3 = 0;
  Float_t www[3] = {12.,14.,17.};
  //Float_t www[3] = {1.,1.,1.};

  for (Int_t i =0;i<3;i++){
    if (nc[i]>2&&nc[i]<1000){
      dedx      += mean[i] *nc[i]*www[i]/sigma[i];
      fSdEdx    += sigma[i]*(nc[i]-2)*www[i]/sigma[i];
      fMAngular += norm[i] *nc[i];    
      norm2     += nc[i]*www[i]/sigma[i];
      norm3     += (nc[i]-2)*www[i]/sigma[i];
    }
    fDEDX[i]  = mean[i];             
    fSDEDX[i] = sigma[i];            
    fNCDEDX[i]= nc[i]; 
  }

  if (norm3>0){
    dedx   /=norm2;
    fSdEdx /=norm3;
    fMAngular/=norm2;
  }
  else{
    SetdEdx(0);
    return;
  }
  //  Float_t dedx1 =dedx;
  
  dedx =0;
  Float_t norm4 = 0;
  for (Int_t i =0;i<3;i++){
    if (nc[i]>2&&nc[i]<1000&&sigma[i]>3){
      //mean[i]   = mean[i]*(1+0.08*(sigma[i]/(fSdEdx)-1.));
      dedx      += mean[i] *(nc[i])/(sigma[i]);
      norm4     += (nc[i])/(sigma[i]);
    }
    fDEDX[i]  = mean[i];                
  }
  if (norm4>0) dedx /= norm4;
  

  
  SetdEdx(dedx);
    
  //mi deDX

}
*/
