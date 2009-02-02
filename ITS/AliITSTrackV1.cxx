//   ITS Track Class
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// It contain all the usefull information for the track and the method to calculate, modify or extract them
// The track is mainly caracterized by the state vector of elements (fX0, fX1, fX2, fX3, fX4) and the
// corresponding covariance matrix of elements (C00, C10, ..... C44) that is triangular
//
#include <Riostream.h>
#include <TMath.h>
#include <TVector.h>
#include <TObjArray.h>
#include "AliRun.h"
#include "AliITSRad.h"
#include "AliITSTrackV1.h"
#include "AliGenerator.h"
//#include "AliMagF.h"


ClassImp(AliITSTrackV1)

AliITSTrackV1::AliITSTrackV1() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// default constructor   
 
  fTPCtrack = 0;
  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=fC44=0.;
  flistCluster = 0;   
  fNumClustInTrack =0;
  fChi2=-1;
  flabel =0;
  fLayer = -1; 
  fClusterInTrack = 0; 
  frtrack=0.;
  fnoclust=0;
  fMass=0.13956995; //a pion by default
  fFieldFactor = 0.0;
  fdEdx = 0.;                          // oggi
  Int_t ia=0;                          // oggi
  for( ia=0; ia<4; ia++) fcor[ia]=0.;  // oggi
  
}
AliITSTrackV1::AliITSTrackV1(Double_t fieldfactor) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// default constructor   
 
  fTPCtrack = 0;
  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=fC44=0.;
  flistCluster = 0;   
  fNumClustInTrack =0;
  fChi2=-1;
  flabel =0; 
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);
  fLayer = -1; 
  fClusterInTrack = 0; 
  frtrack=0.;
  fnoclust=0;     
  fd2.ResizeTo(6);
  ftgl2.ResizeTo(6); 
  fdtgl.ResizeTo(6);
  fMass=0.13956995; //a pion by default
  fdEdx = 0.;
  Int_t ia=0;
  for( ia=0; ia<4; ia++) fcor[ia]=0.;   

  
//////////////////////////////////////// gets magnetic field factor ////////////////////////////////

 // AliMagF * fieldPointer = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField());
  // fFieldFactor =(Double_t)fieldPointer-> SolenoidField()/10/.2;
    fFieldFactor = fieldfactor;
  //cout<< " field factor = "<<fFieldFactor<<"\n"; getchar();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
   
}


 
AliITSTrackV1::AliITSTrackV1(const AliITSTrackV1 &cobj) : TObject(cobj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// copy constructor    

  fClusterInTrack = new TMatrix(6,9);
  Int_t i,j;
  for(i=0; i<6; i++){
  for(j=0; j<9; j++) (*fClusterInTrack)(i,j)=-1.;   
  }  
  flistCluster = new TObjArray; 
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);
  fVertex = cobj.fVertex;
  fErrorVertex = cobj.fErrorVertex;
  flabel = cobj.flabel;
  fLayer=cobj.fLayer;
  fTPCtrack = cobj.fTPCtrack;
  fNumClustInTrack = cobj.fNumClustInTrack; 
  fChi2= cobj.fChi2;
  fX0=cobj.fX0; fX1=cobj.fX1; fX2=cobj.fX2; fX3=cobj.fX3; fX4=cobj.fX4; 
  frtrack=cobj.frtrack;
  fDv=cobj.fDv;
  fZv=cobj.fZv;
  fsigmaDv=cobj.fsigmaDv;
  fsigmaZv=cobj.fsigmaZv;
  fd2.ResizeTo(6);
  ftgl2.ResizeTo(6); 
  fdtgl.ResizeTo(6); 
  fd2=cobj.fd2;
  ftgl2=cobj.ftgl2;
  fdtgl=cobj.fdtgl;
  fnoclust=cobj.fnoclust; 
  fdEdx = cobj.fdEdx;
  Int_t ia=0;
  for( ia=0; ia<4; ia++) fcor[ia]=cobj.fcor[ia];    

    
  fC00=cobj.fC00; fC10=cobj.fC10; fC11=cobj.fC11; fC20=cobj.fC20; fC21=cobj.fC21;
  fC22=cobj.fC22; fC30=cobj.fC30; fC31=cobj.fC31; fC32=cobj.fC32; fC33=cobj.fC33; 
  fC40=cobj.fC40; fC41=cobj.fC41; fC42=cobj.fC42; fC43=cobj.fC43; fC44=cobj.fC44; 
 
  *fClusterInTrack = *cobj.fClusterInTrack;
  
  fFieldFactor=cobj.fFieldFactor;
  fMass=cobj.fMass; 
   
  for(i=0; i<cobj.flistCluster->GetSize(); i++) 
    flistCluster->AddLast(cobj.flistCluster->At(i));
 
}

AliITSTrackV1::AliITSTrackV1(AliTPCtrack &obj, Double_t fieldfactor)
{ 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// special constructor to convert a TPC track into an ITS track

//////////////////////////////////////// gets magnetic field factor ////////////////////////////////

   // AliMagF * fieldPointer = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField());
  // fFieldFactor =(Double_t)fieldPointer-> SolenoidField()/10/.2;
    fFieldFactor = fieldfactor;
 // cout<< " field factor dentro alitrack = "<<fFieldFactor<<"\n";/* getchar();*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////


  fTPCtrack = &obj;
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);

  fd2.ResizeTo(6);
  ftgl2.ResizeTo(6); 
  fdtgl.ResizeTo(6); 
  AliGenerator *gener = gAlice->Generator();
  Float_t vxg,vyg,vzg;
  gener->GetOrigin(vxg,vyg,vzg);
  
 
  fVertex(0)=(Double_t)vxg;
  fVertex(1)=(Double_t)vyg; 
  fVertex(2)=(Double_t)vzg;    
  
  fLayer = 7;
  fClusterInTrack = new TMatrix(6,9);
  
  Int_t i,j;
  for(i=0; i<6; i++){
  for(j=0; j<9; j++) (*fClusterInTrack)(i,j)=-1.;   
  }  
  flistCluster = new TObjArray; 
  fNumClustInTrack = 0;
  fnoclust=0;
  fdEdx = 0.;
  Int_t ia=0;
  for( ia=0; ia<4; ia++) fcor[ia]=0.;          
  LmTPC(); 

}


AliITSTrackV1::~AliITSTrackV1() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
//destructor
 
  if(flistCluster) {
    flistCluster->Delete();
    delete flistCluster;
  }
  if(fClusterInTrack) delete fClusterInTrack;

}     

void AliITSTrackV1::PutCElements(Double_t C00, Double_t C10, Double_t C11, Double_t C20, Double_t C21, 
Double_t C22, Double_t C30, Double_t C31, Double_t C32, Double_t C33, Double_t C40, 
Double_t C41, Double_t C42, Double_t C43, Double_t C44){
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// puts elements of covariance matrix

  fC00=C00; fC10=C10; fC11=C11; fC20=C20; fC21=C21; fC22=C22; fC30=C30; fC31=C31; fC32=C32; fC33=C33;
  fC40=C40; fC41=C41; fC42=C42; fC43=C43; fC44=C44; 
}
  
void AliITSTrackV1::GetCElements(Double_t &C00, Double_t &C10, Double_t &C11, Double_t &C20, Double_t &C21, 
Double_t &C22, Double_t &C30, Double_t &C31, Double_t &C32, Double_t &C33, Double_t &C40, 
Double_t &C41, Double_t &C42, Double_t &C43, Double_t &C44) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// gets elements of covariance matrix


  C00=fC00; C10=fC10; C11=fC11; C20=fC20; C21=fC21; C22=fC22; C30=fC30; C31=fC31; C32=fC32; C33=fC33;
  C40=fC40; C41=fC41; C42=fC42; C43=fC43; C44=fC44; 

}

void AliITSTrackV1::GetXElements(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3, Double_t &X4) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// gets the elements of the state vector
  X0=fX0; X1=fX1; X2=fX2; X3=fX3; X4=fX4; 
}

void AliITSTrackV1::PutXElements(Double_t X0, Double_t X1, Double_t X2, Double_t X3, Double_t X4){
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// puts the elements of the state vector
  fX0=X0; fX1=X1; fX2=X2; fX3=X3; fX4=X4; 
}   

void AliITSTrackV1::LmTPC() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// Transform the TPC state vector from TPC-local to master and build a new state vector ITS-type
// The covariance matrix is also modified accordingly      

   
  Double_t alpha = fTPCtrack->GetAlpha(); 
 
  //printf("LmTPC: alpha %f\n",alpha); 

  Double_t yTPC   = fTPCtrack->GetY();
  Double_t zTPC   = fTPCtrack->GetZ();
  Double_t cTPC   = fTPCtrack->GetC();
  Double_t etaTPC = fTPCtrack->GetEta();
  Double_t tglTPC = fTPCtrack->GetTgl();
  

  Double_t xm, ym, zm;
  Double_t sina = TMath::Sin(alpha);
  Double_t cosa = TMath::Cos(alpha);
  Double_t xl= fTPCtrack->GetX();
  xm = xl * cosa - yTPC*sina;
  ym = xl * sina + yTPC*cosa;
  zm = zTPC;  
  //cout<<" xl e alpha = "<<xl<<" "<<alpha<<"\n"; getchar(); 

  Double_t x0m,y0m;

  ///////////////////////////////////// determine yo //////////////////////////////////////////////////
  
  Double_t vxl=fVertex(0)*cosa+fVertex(1)*sina;
  Double_t vyl= -fVertex(0)*sina+fVertex(1)*cosa;
  Double_t xo,yo, signy;
  Double_t r = 1./cTPC;
  xo =  etaTPC / cTPC;
  Double_t yo1, yo2, diffsq1, diffsq2;  
  yo1 = yTPC +  TMath::Sqrt((r-(xl-xo))*(r+(xl-xo))); 
  yo2 = yTPC -  TMath::Sqrt((r-(xl-xo))*(r+(xl-xo)));   
  diffsq1=TMath::Abs((yo1-vyl)*(yo1-vyl)+((xo-vxl)-r)*((xo-vxl)+r));
  diffsq2=TMath::Abs((yo2-vyl)*(yo2-vyl)+((xo-vxl)-r)*((xo-vxl)+r));
  if(diffsq1<diffsq2) {yo=yo1; signy=1.;} else {yo=yo2; signy=-1.;};
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  x0m = xo * cosa - yo * sina;
  y0m = xo * sina + yo * cosa;  

  frtrack=TMath::Sqrt(xm*xm+ym*ym);
	
  Double_t phi=TMath::ATan2(ym,xm);  if(phi<0) phi=2.*TMath::Pi()+phi;     

  fX0=phi;
  fX1=zm;
  fX3=tglTPC;
  fX4=cTPC;

  
  Double_t dd=TMath::Sqrt((x0m-fVertex(0))*(x0m-fVertex(0))+(y0m-fVertex(1))*(y0m-fVertex(1)));
  Double_t signdd;
  if (r>0) signdd=1.; else signdd=-1.;
  fX2=signdd*dd-r;
  //cout<<" fvertex = "<<fVertex(0)<<" "<<fVertex(1)<<" "<<fVertex(2)<<"\n";
      
  Double_t cov[15];
  fTPCtrack->GetCovariance(cov);

  Double_t  dfidy, dDdy, dDdC, dDdeta;

  dfidy=(xm*cosa+ym*sina)/(frtrack*frtrack);
  dDdy=signdd*((y0m-fVertex(1))*cosa-(x0m-fVertex(0))*sina)/dd;
  Double_t dyodr=signy*(r+(xl-xo)*etaTPC)/TMath::Sqrt((r-(xl-xo))*(r+(xl-xo)));
  Double_t dyomdr=sina*etaTPC+cosa*dyodr;
  Double_t dxomdr=cosa*etaTPC-sina*dyodr;
  Double_t ddddR=((x0m-fVertex(0))*dxomdr+(y0m-fVertex(1))*dyomdr)/dd;
  dDdC=-r*r*(signdd*ddddR-1.);
  Double_t dyoldxol=signy*(xl-xo)/TMath::Sqrt((r-(xl-xo))*(r+(xl-xo)));
  Double_t dxomdeta=r*(cosa-sina*dyoldxol);
  Double_t dyomdeta=r*(sina+cosa*dyoldxol);
  dDdeta=signdd*((x0m-fVertex(0))*dxomdeta+(y0m-fVertex(1))*dyomdeta)/dd;
  
  Double_t f00=dfidy;
  Double_t f20=dDdy;    
  Double_t f22=dDdC;    
  Double_t f23=dDdeta;
  
  Double_t t00=cov[0]*f00;
  Double_t t02=cov[0]*f20+cov[6]*f22+cov[3]*f23;
  Double_t t20=cov[6]*f00;
  Double_t t22=cov[6]*f20+cov[9]*f22+cov[8]*f23;
  
  fC00=f00*t00;
  fC10=cov[1]*f00;
  fC11=cov[2];
  fC20=f20*t00+f22*t20+f23*cov[3]*f00;
  fC21=f20*cov[1]+f22*cov[7]+f23*cov[4];
  fC22=f20*t02+f22*t22+f23*(cov[3]*f20+cov[8]*f22+cov[5]*f23);
  fC30=cov[10]*f00;
  fC31=cov[11];
  fC32=cov[10]*f20+cov[13]*f22+cov[12]*f23;
  fC33=cov[14];
  fC40=t20;
  fC41=cov[7];
  fC42=t22;
  fC43=cov[13];
  fC44=cov[9];
  
   
}


AliITSTrackV1 &AliITSTrackV1::operator=(AliITSTrackV1 obj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// assignement operator

  if (flistCluster) {
    flistCluster->Delete();
    delete flistCluster;
  }
  delete fClusterInTrack;
  fClusterInTrack = new TMatrix(6,9);
  flistCluster = new TObjArray; 
  flabel = obj.flabel;
  fTPCtrack = obj.fTPCtrack; 
  fNumClustInTrack = obj.fNumClustInTrack; 
  fChi2= obj.fChi2;
  fVertex=obj.fVertex;
  fErrorVertex=obj.fErrorVertex;
  fX0=obj.fX0; fX1=obj.fX1; fX2=obj.fX2; fX3=obj.fX3; fX4=obj.fX4;  
  fLayer=obj.fLayer;
  frtrack=obj.frtrack;
  fDv=obj.fDv;
  fZv=obj.fZv;
  fsigmaDv=obj.fsigmaDv;
  fsigmaZv=obj.fsigmaZv;
  fd2=obj.fd2;
  ftgl2=obj.ftgl2; 
  fdtgl=obj.fdtgl;
  fnoclust=obj.fnoclust; 
  
  fC00=obj.fC00; fC10=obj.fC10; fC11=obj.fC11; fC20=obj.fC20; fC21=obj.fC21;
  fC22=obj.fC22; fC30=obj.fC30; fC31=obj.fC31; fC32=obj.fC32; fC33=obj.fC33; 
  fC40=obj.fC40; fC41=obj.fC41; fC42=obj.fC42; fC43=obj.fC43; fC44=obj.fC44;
   
  fMass=obj.fMass; 
  fdEdx = obj.fdEdx; 
  Int_t ia=0;
  for( ia=0; ia<4; ia++) fcor[ia]=obj.fcor[ia];   

  
  *fClusterInTrack = *obj.fClusterInTrack;
  Int_t i;
  for(i=0; i<obj.flistCluster->GetSize(); i++) flistCluster->AddLast(obj.flistCluster->At(i));

  return *this;
  
}

void AliITSTrackV1::PutCluster(Int_t layerc, TVector vecclust) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// put information for clusters
  
  (*fClusterInTrack)(layerc,0) = vecclust(0);
  (*fClusterInTrack)(layerc,1) = vecclust(1);
  (*fClusterInTrack)(layerc,2) = vecclust(2);
  (*fClusterInTrack)(layerc,3) = vecclust(3);
  (*fClusterInTrack)(layerc,4) = vecclust(4);  
  (*fClusterInTrack)(layerc,5) = vecclust(5);
  (*fClusterInTrack)(layerc,6) = vecclust(6);
  (*fClusterInTrack)(layerc,7) = vecclust(7);
  (*fClusterInTrack)(layerc,8) = vecclust(8);

}


void AliITSTrackV1::GetClusters() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// prints the clusters belonging to the current track

  TMatrix a(*fClusterInTrack);
  TMatrix b(6,3);
  Int_t i;
  for(i=0;i<6; i++){
   b(i,0)=a(i,6); b(i,1)=a(i,7); b(i,2)=a(i,8); 
  }
  a.Print();
 // b.Print(); 
 
}


TVector AliITSTrackV1::GetLabTrack(Int_t lay) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// gets the label of the track

  TVector vecLabel(3);
  vecLabel(0)=( (Float_t) (*fClusterInTrack)(lay,6) );
  vecLabel(1)=( (Float_t) (*fClusterInTrack)(lay,7) );
  vecLabel(2)=( (Float_t) (*fClusterInTrack)(lay,8) );
  return vecLabel;
}

void AliITSTrackV1::Search(TVector VecTotLabref, Long_t &labref, Int_t &freq){
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// define label

  Int_t vecfreq[18];
  
  Int_t i,j;
  for(i=0; i<18; i++) vecfreq[i]=0;
  
  for(i=0; i<18; i++) {
    for(j=0; j<18; j++) {
     // if(VecTotLabref(i) == 0.) VecTotLabref(i)=-3.;  //commentato il 5-3-2001
      if( (VecTotLabref(i)>=0.) && (VecTotLabref(i)==VecTotLabref(j)) ) vecfreq[i]++;    
    }  
  }
  Int_t imax=-1000;
  Long_t  labdefault= (Long_t)1000000.;
  freq=0;
  for(i=0; i<18; i++) {
    if(vecfreq[i]>freq) {freq=vecfreq[i]; imax=i;}  
  }
  if(imax<0) labref=labdefault; else labref=(Long_t) VecTotLabref(imax);
} 


void AliITSTrackV1::Propagation(Double_t rk) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
//Propagation of track 
  Double_t duepi=2.*TMath::Pi();
  Double_t rkm1=frtrack;
//cout<<" rk e rkm1 dentro Propagation   "<<rk<<" "<<rkm1<<"\n";
  
  //
  Double_t aAk=ArgA(rk), aAkm1=ArgA(rkm1);
  Double_t ak=Arga(rk), akm1=Arga(rkm1);
  fX0+=TMath::ASin(aAk)-TMath::ASin(aAkm1);

  if(fX0>duepi) fX0-=duepi;
  if(fX0<0.) fX0+=duepi;
	
  Double_t tgl=fX3;
  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2;
  fX1+=tgl/cy*(TMath::ASin(ak)-TMath::ASin(akm1));  
  frtrack=rk;

	
  Double_t bk=ArgB(rk), bkm1=ArgB(rkm1);	
  Double_t ck=ArgC(rk), ckm1=ArgC(rkm1);  

  Double_t f02=ck/TMath::Sqrt((1.-aAk)*(1.+aAk)) - ckm1/TMath::Sqrt((1.-aAkm1)*(1.+aAkm1));
  Double_t f04=bk/TMath::Sqrt((1.-aAk)*(1.+aAk)) - bkm1/TMath::Sqrt((1.-aAkm1)*(1.+aAkm1));  	
  Double_t f12=tgl*d*(1./rk - 1./rkm1);
  Double_t f13=rk - rkm1;
  

  Double_t c00=fC00;
  Double_t c10=fC10;
  Double_t c11=fC11;
  Double_t c20=fC20;
  Double_t c21=fC21;
  Double_t c22=fC22;
  Double_t c30=fC30;
  Double_t c31=fC31;
  Double_t c32=fC32;
  Double_t c33=fC33;
  Double_t c40=fC40;
  Double_t c41=fC41;
  Double_t c42=fC42;
  Double_t c43=fC43;
  Double_t c44=fC44;
  
  Double_t r10=c10+c21*f02+c41*f04;
  Double_t r20=c20+c22*f02+c42*f04;
  Double_t r30=c30+c32*f02+c43*f04;
  Double_t r40=c40+c42*f02+c44*f04;
  Double_t r21=c21+c22*f12+c32*f13;
  Double_t r31=c31+c32*f12+c33*f13;									
  Double_t r41=c41+c42*f12+c43*f13;

  fC00=c00+c20*f02+c40*f04+f02*r20+f04*r40;
  fC10=r10+f12*r20+f13*r30;
  fC11=c11+c21*f12+c31*f13+f12*r21+f13*r31;
  fC20=r20;
  fC21=r21;
  fC30=r30;
  fC31=r31;
  fC40=r40;
  fC41=r41;
  
}

void AliITSTrackV1::AddEL(Double_t signdE, Bool_t flagtot, Double_t mass) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
//  add energy loss
// AliITSRad *rl was passed as argument. Now rl has been commented out

  mass=fMass;  
  
  TVector s(6);  
  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
  //s(0)=0.0026+0.00283*2.; s(1)=0.018*2.; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
  //0.00277 is added in the first layer to take into account the energy loss in the beam pipe

  //for(int k=0; k<6; k++) cout<<s(k)<<" "; cout<<"\n";
  Int_t k;
  //for(k=0; k<6; k++) s(k)=s(k)*1.6;    //forint
  for(k=0; k<6; k++) s(k)=s(k)*1.7;    //forint
  
  Double_t phi=fX0;
  
  if(phi<0.174 ) s(5)=s(5)+0.012; 
  if(phi>6.1 ) s(5)=s(5)+0.012; // to take into account rail 
  if(phi>2.96 && phi<3.31 ) s(5)=s(5)+0.012;   
  
  /* 
  Double_t tgl=fX3;
  Double_t theta=((TMath::Pi())/2.)-TMath::ATan(tgl);
  //phi*=180./TMath::Pi();
  //theta*=180./TMath::Pi();
  //Double_t rad90=(TMath::Pi())/2.;
  Double_t rad40=(TMath::Pi())*40./180.;
  Double_t rad100=(TMath::Pi())*100/180;
  Double_t rad360=(TMath::Pi())*2.;
  Int_t imax=rl->Getimax();
  Int_t jmax=rl->Getjmax();
  Int_t i=(Int_t) ( (theta-rad40)/rad100*imax);
  Int_t j=(Int_t) ( phi/rad360*jmax );
  //Int_t i=(Int_t)( ((theta-((TMath::Pi())/4.))/((TMath::Pi())/2.))*imax );
  //Int_t j=(Int_t)( (phi/((TMath::Pi())*2.))*jmax );
  if(i<0) i=0;
  if(i>=imax) i=imax-1;
  if(j<0) j=0;
  if(j>=jmax) j=jmax-1;
  
  s(0) = 0.0028/TMath::Sin(theta)+( rl->GetRadMatrix1() )(i,j);   // 0.0028 takes into account the beam pipe
  s(1) = ( rl->GetRadMatrix2() )(i,j);
  s(2) = ( rl->GetRadMatrix3() )(i,j);
  s(3) = ( rl->GetRadMatrix4() )(i,j);
  s(4) = ( rl->GetRadMatrix5() )(i,j);
  s(5) = ( rl->GetRadMatrix6() )(i,j);
  

  */  
  
  //for(int k=0; k<6; k++) cout<<s(k)<<" "; getchar();
    
  //if(phi>60) {cout<<" phi = "<<phi<<"\n"; getchar();}
  //if(theta<45 || theta>135) {cout<<" theta = "<<theta<<"\n"; getchar();}
  //cout<<" dentro AddEl: phi, theta = "<<phi<<" "<<theta<<"\n"; getchar(); 
    
  Double_t cl=1.+fX3*fX3;  // cl=1/(cosl)**2 = 1 + (tgl)**2
  Double_t sqcl=TMath::Sqrt(cl);
  Double_t pt=GetPt();
     
  Double_t p2=pt*pt*cl;
  Double_t e=TMath::Sqrt(p2+mass*mass);
  Double_t beta2=p2/(p2+mass*mass);
  
  Double_t dE;
  if(flagtot) {
    Double_t stot=s(0)+s(1)+s(2)+s(3)+s(4)+s(5);
    dE=0.153/beta2*(log(5940*beta2/(1-beta2)) - beta2)*stot*21.82*sqcl;
  } else {
    dE=0.153/beta2*(log(5940*beta2/(1-beta2)) - beta2)*s(fLayer-1)*21.82*sqcl;
  }   
  dE=signdE*dE/1000.; 
	  
  e+=dE;
  Double_t p=TMath::Sqrt((e-mass)*(e+mass));   
  Double_t sign=1.;
  if(fX4 < 0.) sign=-1.; 
  pt=sign*p/sqcl; 
  Double_t cc=(0.299792458*0.2*fFieldFactor)/(pt*100.);
  fX4=cc;
  
}

void  AliITSTrackV1::Correct(Double_t rk) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// correct track to take into account real geometry detector

  Double_t duepi=2.*TMath::Pi();
  Double_t rkm1=frtrack;
  Double_t aAk=ArgA(rk), aAkm1=ArgA(rkm1);
  Double_t ak=Arga(rk), akm1=Arga(rkm1);

  fX0+=TMath::ASin(aAk)-TMath::ASin(aAkm1);
  if(fX0>duepi) fX0-=duepi;
  if(fX0<0.) fX0+=duepi;
	
  Double_t tgl=fX3;
  Double_t c=fX4;
  Double_t cy=c/2;
  fX1+=tgl/cy*(TMath::ASin(ak)-TMath::ASin(akm1));
  frtrack=rk;
  		
}

void AliITSTrackV1::AddMS(Double_t mass) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it       
//////////   Modification of the covariance matrix to take into account multiple scattering  ///////////     
   
  mass=fMass;
 
  TVector s(6);

  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
 //  s(0)=0.0026+0.00283*2.; s(1)=0.018*2.; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
//0.00277 is added in the first layer to take into account the energy loss in the beam pipe 

  Int_t k;
  //for(k=0; k<6; k++) s(k)=s(k)*1.6;   // forint
  for(k=0; k<6; k++) s(k)=s(k)*1.7;   // forint
  
  Double_t phi=fX0;
 if(phi<0.174 ) s(5)=s(5)+0.012; //aggiunta provvisoria
  if(phi>6.1 ) s(5)=s(5)+0.012; //aggiunta provvisoria  
   if(phi>2.96 && phi< 3.31) s(5)=s(5)+0.012; //aggiunta provvisoria 
  
      
  Double_t tgl=fX3;
  /*
  Double_t theta=((TMath::Pi())/2.)-TMath::ATan(tgl);
  Double_t rad40=(TMath::Pi())*40./180.;      // rivedere
  Double_t rad100=(TMath::Pi())*100/180;
  Double_t rad360=(TMath::Pi())*2.;
  Int_t imax=rl->Getimax();
  Int_t jmax=rl->Getjmax();
  Int_t i=(Int_t) ( (theta-rad40)/rad100*imax);
  Int_t j=(Int_t) ( phi/rad360*jmax);

  if(i<0) i=0;
  if(i>=imax) i=imax-1;
  if(j<0) j=0;
  if(j>=jmax) j=jmax-1;
  
  s(0) = 0.0028/TMath::Sin(theta)+( rl->GetRadMatrix1() )(i,j);   // 0.0028 takes into account the beam pipe
  s(1) = ( rl->GetRadMatrix2() )(i,j);
  s(2) = ( rl->GetRadMatrix3() )(i,j);
  s(3) = ( rl->GetRadMatrix4() )(i,j);
  s(4) = ( rl->GetRadMatrix5() )(i,j);
  s(5) = ( rl->GetRadMatrix6() )(i,j);
   */   
  //Double_t mass=0.1396;
  Int_t layer=(Int_t)GetLayer();
  
  Double_t cosl=TMath::Cos(TMath::ATan(tgl));	
  Double_t d=fX2;
  Double_t c=fX4;
  Double_t cy=c/2.;
  Double_t q20=1./(cosl*cosl);	 
  Double_t q30=c*tgl;
   
  Double_t q40=cy*(frtrack*frtrack-d*d)/(1.+ 2.*cy*d);
  Double_t dd=d+cy*d*d-cy*frtrack*frtrack;
  Double_t dprova=frtrack*frtrack - dd*dd;
  Double_t q41=0.;
  if(dprova>0.) q41=-1./cosl*TMath::Sqrt(dprova)/(1.+ 2.*cy*d);
  	 	
  Double_t p2=(GetPt()*GetPt())/(cosl*cosl);
  Double_t beta2=p2/(p2+mass*mass);
 // Double_t theta2=14.1*14.1/(beta2*p2*1.e6)*(s(layer-1)/cosl);
   Double_t theta2=14.1*14.1/(beta2*p2*1.e6)*(s(layer-1)/TMath::Abs(cosl));

  fC22+=theta2*(q40*q40+q41*q41);
  fC32+=theta2*q20*q40;
  fC33+=theta2*q20*q20;
  fC42+=theta2*q30*q40;
  fC43+=theta2*q30*q20;
  fC44+=theta2*q30*q30;
    
}
void AliITSTrackV1::PrimaryTrack() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// calculation of part of covariance matrix for vertex constraint

  Double_t rlayer[6];
	
  rlayer[0]=4.; rlayer[1]=7.;  rlayer[2]=14.9;  rlayer[3]=23.8;  
  rlayer[4]=39.1;  rlayer[5]=43.6;
  	
  Double_t cy=fX4/2.;	
  Double_t tgl=(fX1-fZv)*cy/TMath::ASin(cy*frtrack);
  Double_t frtrack=1.;	
  fX0=0.;
  fX1=frtrack*tgl; 
  fX2=fDv;
  fX3=tgl; 

  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=0.;

  AddEL(1.,1);
  fLayer=0;
  Int_t i;
  for (i=0; i<6; i++) {
    Propagation(rlayer[i]);
    fLayer++;
    fd2(i)=fC22;
    ftgl2(i)=fC33;
    fdtgl(i)=fC32; 
    AddMS();    
    AddEL(-1,0); 	   
  } 		
}	
 
  
Int_t AliITSTrackV1::DoNotCross(Double_t rk) const{
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// determine if the track cross a layer

  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2.;
  return (TMath::Abs((cy*rk+(1.+cy*d)*d/rk)/(1.+2.*cy*d))>=1.)?1:0;
}
  

Double_t AliITSTrackV1::ArgA(Double_t rk) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// quantity usefull in Propagation  
  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2.;
  Double_t arg=(cy*rk + (1 + cy*d)*d/rk)/(1.+ 2.*cy*d);		
  if (TMath::Abs(arg) < 1.) return arg;
  //cout<<"class AliITSTrack: ArgA out of range !\n";/* getchar();*/
  return (arg>0) ? 0.99999999999 : -0.9999999999;
}
   
Double_t AliITSTrackV1::Arga(Double_t rk) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// quantity usefull in Propagation 
  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2.;
  Double_t arg=(rk*rk - d*d)/(1.+ 2.*cy*d);		
  if (arg<0.) {/*cout<<"class AliITSTrack: Arga out of range !\n";*/ arg=0.;} 
  return cy*TMath::Sqrt(arg);
}
   	
Double_t AliITSTrackV1::ArgB(Double_t rk) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// quantity usefull in Propagation 
  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2.;	   
  return (rk*rk - d*d)/(rk*(1.+ 2.*cy*d)*(1.+ 2.*cy*d));
}
   
Double_t AliITSTrackV1::ArgC(Double_t rk) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// quantity usefull in Propagation 
  Double_t c=fX4;
  Double_t d=fX2;
  Double_t cy=c/2.;		
  return  (1./rk - 2.*cy*ArgA(rk)/(1.+ 2.*cy*d));
}


Double_t AliITSTrackV1::GetPredChi2(Double_t m[2], Double_t sigma[2] ) const {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// This function calculates a predicted chi2 increment.
 
  Double_t r00=sigma[0], r01=0., r11=sigma[1];
  r00+=fC00; r01+=fC10; r11+=fC11;
  Double_t det=r00*r11-r01*r01;
  if(TMath::Abs(det) < 1.e-15) {cout<<" Problems on matrix in GetPredChi2 "<<det<<"\n";
  return 1e10;}
  Double_t tmp=r00; r00=r11; r11=tmp; r01=-r01;
  Double_t dphi=m[0]-fX0; 
  Double_t dz=m[1]-fX1;
  Double_t chi2 = (dphi*r00*dphi +2.*r01*dphi*dz + dz*r11*dz)/det;
  return chi2;
  
}
