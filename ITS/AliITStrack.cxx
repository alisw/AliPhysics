#include <iostream.h>
#include <TMath.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TObjArray.h>

#include "AliITS.h"
#include "AliRun.h"
#include "AliITStrack.h"
#include "AliGenerator.h"
#include "AliITSRad.h"


ClassImp(AliITStrack)

AliITStrack::AliITStrack() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// default constructor   
 
  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=fC44=0.;
  flistCluster = new TObjArray; 
  fNumClustInTrack =0;
  fChi2=-1;
  flabel =0; 
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);
  fLayer = -1; 
  ClusterInTrack = new TMatrix(6,9);
  Int_t i;
  for(i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
                           (*ClusterInTrack)(i,8)=-1.;
  rtrack=0.; 
  d2.ResizeTo(6);
  tgl2.ResizeTo(6); 
  dtgl.ResizeTo(6);   
}


 
AliITStrack::AliITStrack(const AliITStrack &cobj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it    

  ClusterInTrack = new TMatrix(6,9);
  Int_t i;
  for(i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
                           (*ClusterInTrack)(i,8)=-1.;
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
  rtrack=cobj.rtrack;
  Dv=cobj.Dv;
  Zv=cobj.Zv;
  sigmaDv=cobj.sigmaDv;
  sigmaZv=cobj.sigmaZv;
  d2.ResizeTo(6);
  tgl2.ResizeTo(6); 
  dtgl.ResizeTo(6); 
  d2=cobj.d2;
  tgl2=cobj.tgl2;
  dtgl=cobj.dtgl;
    
  fC00=cobj.fC00; fC10=cobj.fC10; fC11=cobj.fC11; fC20=cobj.fC20; fC21=cobj.fC21;
  fC22=cobj.fC22; fC30=cobj.fC30; fC31=cobj.fC31; fC32=cobj.fC32; fC33=cobj.fC33; 
  fC40=cobj.fC40; fC41=cobj.fC41; fC42=cobj.fC42; fC43=cobj.fC43; fC44=cobj.fC44; 
 
  *ClusterInTrack = *cobj.ClusterInTrack;
 
  for(i=0; i<cobj.flistCluster->GetSize(); i++) 
    flistCluster->AddLast(cobj.flistCluster->At(i));
 
}

AliITStrack::AliITStrack(AliTPCtrack &obj)
{ 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 

  fTPCtrack = &obj;
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);

  d2.ResizeTo(6);
  tgl2.ResizeTo(6); 
  dtgl.ResizeTo(6); 
  AliGenerator *gener = gAlice->Generator();
  Float_t Vxg,Vyg,Vzg;
  gener->GetOrigin(Vxg,Vyg,Vzg);
  
 
  fVertex(0)=(Double_t)Vxg;
  fVertex(1)=(Double_t)Vyg; 
  fVertex(2)=(Double_t)Vzg;    
  
  fLayer = 7;
  //fmCovariance = new TMatrix(5,5);
  ClusterInTrack = new TMatrix(6,9);
  
  Int_t i;
  for(i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
                           (*ClusterInTrack)(i,8)=-1.;
  flistCluster = new TObjArray; 
  fNumClustInTrack = 0; 
  LmTPC(); 

}


AliITStrack::~AliITStrack() { 

  //destructor
 
  if(flistCluster) delete flistCluster; 
  if(ClusterInTrack) delete ClusterInTrack;

}     

void AliITStrack::PutCElements(Double_t C00, Double_t C10, Double_t C11, Double_t C20, Double_t C21, 
Double_t C22, Double_t C30, Double_t C31, Double_t C32, Double_t C33, Double_t C40, 
Double_t C41, Double_t C42, Double_t C43, Double_t C44){

  fC00=C00; fC10=C10; fC11=C11; fC20=C20; fC21=C21; fC22=C22; fC30=C30; fC31=C31; fC32=C32; fC33=C33;
  fC40=C40; fC41=C41; fC42=C42; fC43=C43; fC44=C44; 
}
  
void AliITStrack::GetCElements(Double_t &C00, Double_t &C10, Double_t &C11, Double_t &C20, Double_t &C21, 
Double_t &C22, Double_t &C30, Double_t &C31, Double_t &C32, Double_t &C33, Double_t &C40, 
Double_t &C41, Double_t &C42, Double_t &C43, Double_t &C44){

  C00=fC00; C10=fC10; C11=fC11; C20=fC20; C21=fC21; C22=fC22; C30=fC30; C31=fC31; C32=fC32; C33=fC33;
  C40=fC40; C41=fC41; C42=fC42; C43=fC43; C44=fC44; 

}

void AliITStrack::GetXElements(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3, Double_t &X4) {
  X0=fX0; X1=fX1; X2=fX2; X3=fX3; X4=fX4; 
}

void AliITStrack::PutXElements(Double_t X0, Double_t X1, Double_t X2, Double_t X3, Double_t X4){
  fX0=X0; fX1=X1; fX2=X2; fX3=X3; fX4=X4; 
}   

void AliITStrack::LmTPC() {

// Transform the TPC state vector from TPC-local to master and build a new state vector ITS-type
// The covariance matrix is also modified accordingly
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it      

   
  Double_t Alpha = fTPCtrack->GetAlpha(); 
 
  //printf("LmTPC: Alpha %f\n",Alpha); 

  Double_t YTPC   = fTPCtrack->GetY();
  Double_t ZTPC   = fTPCtrack->GetZ();
  Double_t CTPC   = fTPCtrack->GetC();
  Double_t EtaTPC = fTPCtrack->GetEta();
  Double_t TglTPC = fTPCtrack->GetTgl();
  

  Double_t xm, ym, zm;
  Double_t sina = TMath::Sin(Alpha);
  Double_t cosa = TMath::Cos(Alpha);
  Double_t xl= fTPCtrack->GetX();
  xm = xl * cosa - YTPC*sina;
  ym = xl * sina + YTPC*cosa;
  zm = ZTPC;  
  //cout<<" xl e alpha = "<<xl<<" "<<Alpha<<"\n"; getchar(); 

  Double_t x0m,y0m;

  ///////////////////////////////////// determine yo //////////////////////////////////////////////////
  
  Double_t Vxl=fVertex(0)*cosa+fVertex(1)*sina;
  Double_t Vyl= -fVertex(0)*sina+fVertex(1)*cosa;
  Double_t Xo,Yo, signy;
  Double_t R = 1./CTPC;
  Xo =  EtaTPC / CTPC;
  xoTPC=Xo;
  Double_t Yo1, Yo2, diffsq1, diffsq2;  
  Yo1 = YTPC +  TMath::Sqrt(R*R - (xl-Xo)*(xl-Xo)); 
  Yo2 = YTPC -  TMath::Sqrt(R*R - (xl-Xo)*(xl-Xo));   
  diffsq1=TMath::Abs((Yo1-Vyl)*(Yo1-Vyl)+(Xo-Vxl)*(Xo-Vxl)-R*R);
  diffsq2=TMath::Abs((Yo2-Vyl)*(Yo2- Vyl)+(Xo-Vxl)*(Xo-Vxl)-R*R);
  if(diffsq1<diffsq2) {Yo=Yo1; signy=1.;} else {Yo=Yo2; signy=-1.;};
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  x0m = Xo * cosa - Yo * sina;
  y0m = Xo * sina + Yo * cosa;  

  rtrack=TMath::Sqrt(xm*xm+ym*ym);
  Double_t pigre=TMath::Pi();
  Double_t phi=0.0;
  if(ym == 0. || xm == 0.) {
    if(ym == 0. && xm == 0.) {cout << "  Error in AliITStrack::LmTPC x=0 and y=0 \n"; getchar();}
    if(ym ==0. && xm>0.) phi=0.;
    if(ym==0. && xm<0.) phi=pigre;
    if(xm==0 && ym>0.) phi=pigre/2.;
    if(xm==0 && ym<0.) phi=1.5*pigre;   
  }
  else {
    if (xm>0. && ym>0.) phi=TMath::ATan(ym/xm);
    if (xm<0. && ym>0.) phi=pigre+TMath::ATan(ym/xm);
    if (xm<0. && ym<0.) phi=pigre+TMath::ATan(ym/xm); 
    if (xm>0. && ym<0.) phi=(2*pigre)+TMath::ATan(ym/xm);     
  };
  if(phi<0. || phi>(2*pigre)) {cout<<"attention error on phi in  AliITStrack:LmTPC \n"; getchar();}  

  fX0=phi;
  fX1=zm;
  fX3=TglTPC;
  fX4=CTPC;

  
  Double_t dd=TMath::Sqrt((x0m-fVertex(0))*(x0m-fVertex(0))+(y0m-fVertex(1))*(y0m-fVertex(1)));
  Double_t signdd;
  if (R>0) signdd=1.; else signdd=-1.;
  fX2=signdd*dd-R;
  //cout<<" fvertex = "<<fVertex(0)<<" "<<fVertex(1)<<" "<<fVertex(2)<<"\n";
      
  Double_t cov[15];
  fTPCtrack->GetCovariance(cov);

  Double_t  dfidy, dDdy, dDdC, dDdeta;

  dfidy=(xm*cosa+ym*sina)/(rtrack*rtrack);
  dDdy=signdd*((y0m-fVertex(1))*cosa-(x0m-fVertex(0))*sina)/dd;
  Double_t dyodr=signy*(R+(xl-Xo)*EtaTPC)/TMath::Sqrt(R*R-(xl-Xo)*(xl-Xo));
  Double_t dyomdr=sina*EtaTPC+cosa*dyodr;
  Double_t dxomdr=cosa*EtaTPC-sina*dyodr;
  Double_t ddddR=((x0m-fVertex(0))*dxomdr+(y0m-fVertex(1))*dyomdr)/dd;
  dDdC=-R*R*(signdd*ddddR-1.);
  Double_t dyoldxol=signy*(xl-Xo)/TMath::Sqrt(R*R-(xl-Xo)*(xl-Xo));
  Double_t dxomdeta=R*(cosa-sina*dyoldxol);
  Double_t dyomdeta=R*(sina+cosa*dyoldxol);
  dDdeta=signdd*((x0m-fVertex(0))*dxomdeta+(y0m-fVertex(1))*dyomdeta)/dd;
  
  Double_t F00=dfidy;
  Double_t F20=dDdy;    
  Double_t F22=dDdC;    
  Double_t F23=dDdeta;
  
  Double_t T00=cov[0]*F00;
  Double_t T02=cov[0]*F20+cov[6]*F22+cov[3]*F23;
  Double_t T20=cov[6]*F00;
  Double_t T22=cov[6]*F20+cov[9]*F22+cov[8]*F23;
  
  fC00=F00*T00;
  fC10=cov[1]*F00;
  fC11=cov[2];
  fC20=F20*T00+F22*T20+F23*cov[3]*F00;
  fC21=F20*cov[1]+F22*cov[7]+F23*cov[4];
  fC22=F20*T02+F22*T22+F23*(cov[3]*F20+cov[8]*F22+cov[5]*F23);
  fC30=cov[10]*F00;
  fC31=cov[11];
  fC32=cov[10]*F20+cov[13]*F22+cov[12]*F23;
  fC33=cov[14];
  fC40=T20;
  fC41=cov[7];
  fC42=T22;
  fC43=cov[13];
  fC44=cov[9];
  
  //cout<<" C32 e C44 = "<<fC32<<" "<<fC44<<"\n"; getchar();
   
}


AliITStrack &AliITStrack::operator=(AliITStrack obj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  

  delete flistCluster;
  delete ClusterInTrack;
  ClusterInTrack = new TMatrix(6,9);
  flistCluster = new TObjArray; 
  flabel = obj.flabel;
  fTPCtrack = obj.fTPCtrack; 
  fNumClustInTrack = obj.fNumClustInTrack; 
  fChi2= obj.fChi2;
  fVertex=obj.fVertex;
  fErrorVertex=obj.fErrorVertex;
  fX0=obj.fX0; fX1=obj.fX1; fX2=obj.fX2; fX3=obj.fX3; fX4=obj.fX4;  
  fLayer=obj.fLayer;
  rtrack=obj.rtrack;
  Dv=obj.Dv;
  Zv=obj.Zv;
  sigmaDv=obj.sigmaDv;
  sigmaZv=obj.sigmaZv;
  d2=obj.d2;
  tgl2=obj.tgl2; 
  dtgl=obj.dtgl; 
  
  fC00=obj.fC00; fC10=obj.fC10; fC11=obj.fC11; fC20=obj.fC20; fC21=obj.fC21;
  fC22=obj.fC22; fC30=obj.fC30; fC31=obj.fC31; fC32=obj.fC32; fC33=obj.fC33; 
  fC40=obj.fC40; fC41=obj.fC41; fC42=obj.fC42; fC43=obj.fC43; fC44=obj.fC44; 
  
  
  *ClusterInTrack = *obj.ClusterInTrack;
  Int_t i;
  for(i=0; i<obj.flistCluster->GetSize(); i++) flistCluster->AddLast(obj.flistCluster->At(i));

  return *this;
  
}

void AliITStrack::PutCluster(Int_t layerc, TVector vecclust) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
  
  (*ClusterInTrack)(layerc,0) = vecclust(0);
  (*ClusterInTrack)(layerc,1) = vecclust(1);
  (*ClusterInTrack)(layerc,2) = vecclust(2);
  (*ClusterInTrack)(layerc,3) = vecclust(3);
  (*ClusterInTrack)(layerc,4) = vecclust(4);  
  (*ClusterInTrack)(layerc,5) = vecclust(5);
  (*ClusterInTrack)(layerc,6) = vecclust(6);
  (*ClusterInTrack)(layerc,7) = vecclust(7);
  (*ClusterInTrack)(layerc,8) = vecclust(8);

}


void AliITStrack::GetClusters() {

  TMatrix A(*ClusterInTrack);
  TMatrix B(6,3);
  Int_t i;
  for(i=0;i<6; i++){
   B(i,0)=A(i,6); B(i,1)=A(i,7); B(i,2)=A(i,8); 
  }
  A.Print();
 // B.Print(); 
 
}


TVector AliITStrack::GetLabTrack(Int_t lay) {
  TVector VecLabel(3);
  VecLabel(0)=( (Float_t) (*ClusterInTrack)(lay,6) );
  VecLabel(1)=( (Float_t) (*ClusterInTrack)(lay,7) );
  VecLabel(2)=( (Float_t) (*ClusterInTrack)(lay,8) );
  return VecLabel;
}

void AliITStrack::Search(TVector VecTotLabref, Long_t &labref, Int_t &freq){
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


void AliITStrack::Propagation(Double_t rk) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
//Propagation of track 
  Double_t duepi=2.*TMath::Pi();
  Double_t rkm1=rtrack;
//cout<<" rk e rkm1 dentro Propagation   "<<rk<<" "<<rkm1<<"\n";
  
  //
  Double_t Ak=argA(rk), Akm1=argA(rkm1);
  Double_t ak=arga(rk), akm1=arga(rkm1);
  fX0+=TMath::ASin(Ak)-TMath::ASin(Akm1);

  if(fX0>duepi) fX0-=duepi;
  if(fX0<0.) fX0+=duepi;
	
  Double_t tgl=fX3;
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2;
  fX1+=tgl/Cy*(TMath::ASin(ak)-TMath::ASin(akm1));  
  rtrack=rk;

	
  Double_t Bk=argB(rk), Bkm1=argB(rkm1);	
  Double_t Ck=argC(rk), Ckm1=argC(rkm1);  

  Double_t F02=Ck/TMath::Sqrt(1.-Ak*Ak) - Ckm1/TMath::Sqrt(1.-Akm1*Akm1);
  Double_t F04=Bk/TMath::Sqrt(1.-Ak*Ak) - Bkm1/TMath::Sqrt(1.-Akm1*Akm1);  	
  Double_t F12=tgl*D*(1./rk - 1./rkm1);
  Double_t F13=rk - rkm1;
  

  Double_t C00=fC00;
  Double_t C10=fC10;
  Double_t C11=fC11;
  Double_t C20=fC20;
  Double_t C21=fC21;
  Double_t C22=fC22;
  Double_t C30=fC30;
  Double_t C31=fC31;      //    provare se si puo' fare a meno
  Double_t C32=fC32;
  Double_t C33=fC33;
  Double_t C40=fC40;
  Double_t C41=fC41;
  Double_t C42=fC42;
  Double_t C43=fC43;
  Double_t C44=fC44;
  
  Double_t R10=C10+C21*F02+C41*F04;
  Double_t R20=C20+C22*F02+C42*F04;
  Double_t R30=C30+C32*F02+C43*F04;
  Double_t R40=C40+C42*F02+C44*F04;
  Double_t R21=C21+C22*F12+C32*F13;
  Double_t R31=C31+C32*F12+C33*F13;
  Double_t R41=C41+C42*F12+fC43*F13;

  fC00=C00+C20*F02+C40*F04+F02*R20+F04*R40;
  fC10=R10+F12*R20+F13*R30;
  fC11=C11+C21*F12+C31*F13+F12*R21+F13*R31;
  fC20=R20;
  fC21=R21;
  fC30=R30;
  fC31=R31;
  fC40=R40;
  fC41=R41;
  
}

void AliITStrack::AddEL(AliITSRad *rl, Double_t signdE, Bool_t flagtot, Double_t mass) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
//  add energy loss

  TVector s(6);  
  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
  //0.00277 is added in the first layer to take into account the energy loss in the beam pipe

  //for(int k=0; k<6; k++) cout<<s(k)<<" "; cout<<"\n";
  for(int k=0; k<6; k++) s(k)=s(k)*1.6;
 
  
  Double_t phi=fX0;
  
  if(phi<0.174 ) s(5)=s(5)+0.012; 
   if(phi>6.1 ) s(5)=s(5)+0.012; // to take into account rail 
   if(phi>2.96 && phi<3.31 ) s(5)=s(5)+0.012;   
  
   
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
  /*
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
  Double_t E=TMath::Sqrt(p2+mass*mass);
  Double_t beta2=p2/(p2+mass*mass);
  
  Double_t dE;
  if(flagtot) {
    Double_t stot=s(0)+s(1)+s(2)+s(3)+s(4)+s(5);
    dE=0.153/beta2*(log(5940*beta2/(1-beta2)) - beta2)*stot*21.82*sqcl;
  } else {
    dE=0.153/beta2*(log(5940*beta2/(1-beta2)) - beta2)*s(fLayer-1)*21.82*sqcl;
  }   
  dE=signdE*dE/1000.; 
	  
  E+=dE;
  Double_t p=TMath::Sqrt(E*E-mass*mass);   
  Double_t sign=1.;
  if(fX4 < 0.) sign=-1.; 
  pt=sign*p/sqcl; 
  Double_t CC=(0.299792458*0.2)/(pt*100.);
  //Double_t CC=(0.299792458*0.5)/(pt*100.);
  fX4=CC;
  
}

void  AliITStrack::Correct(Double_t rk) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// correct track to take into account real geometry detector

  Double_t duepi=2.*TMath::Pi();
  Double_t rkm1=rtrack;
  Double_t Ak=argA(rk), Akm1=argA(rkm1);
  Double_t ak=arga(rk), akm1=arga(rkm1);

  fX0+=TMath::ASin(Ak)-TMath::ASin(Akm1);
  if(fX0>duepi) fX0-=duepi;
  if(fX0<0.) fX0+=duepi;
	
  Double_t tgl=fX3;
  Double_t C=fX4;
  Double_t Cy=C/2;
  fX1+=tgl/Cy*(TMath::ASin(ak)-TMath::ASin(akm1));
  rtrack=rk;
  		
}

void AliITStrack::AddMS(AliITSRad *rl) {
       
//////////   Modification of the covariance matrix to take into account multiple scattering  ///////////
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it     
   
  TVector s(6);

  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
//0.00277 is added in the first layer to take into account the energy loss in the beam pipe 

  for(int k=0; k<6; k++) s(k)=s(k)*1.6;
  
  Double_t phi=fX0;
 if(phi<0.174 ) s(5)=s(5)+0.012; //aggiunta provvisoria
  if(phi>6.1 ) s(5)=s(5)+0.012; //aggiunta provvisoria  
   if(phi>2.96 && phi< 3.31) s(5)=s(5)+0.012; //aggiunta provvisoria 
      
  Double_t tgl=fX3;
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
  /*
  s(0) = 0.0028/TMath::Sin(theta)+( rl->GetRadMatrix1() )(i,j);   // 0.0028 takes into account the beam pipe
  s(1) = ( rl->GetRadMatrix2() )(i,j);
  s(2) = ( rl->GetRadMatrix3() )(i,j);
  s(3) = ( rl->GetRadMatrix4() )(i,j);
  s(4) = ( rl->GetRadMatrix5() )(i,j);
  s(5) = ( rl->GetRadMatrix6() )(i,j);
   */   
  Double_t mass=0.1396;
  Int_t layer=(Int_t)GetLayer();
  
  Double_t cosl=TMath::Cos(TMath::ATan(tgl));	
  Double_t D=fX2;
  Double_t C=fX4;
  Double_t Cy=C/2.;
  Double_t Q20=1./(cosl*cosl);	 
  Double_t Q30=C*tgl;
   
  Double_t Q40=Cy*(rtrack*rtrack-D*D)/(1.+ 2.*Cy*D);
  Double_t dd=D+Cy*D*D-Cy*rtrack*rtrack;
  Double_t dprova=rtrack*rtrack - dd*dd;
  Double_t Q41=0.;
  if(dprova>0.) Q41=-1./cosl*TMath::Sqrt(dprova)/(1.+ 2.*Cy*D);
  	 	
  Double_t p2=(GetPt()*GetPt())/(cosl*cosl);
  Double_t beta2=p2/(p2+mass*mass);
  Double_t theta2=14.1*14.1/(beta2*p2*1.e6)*(s(layer-1)/cosl);

  fC22+=theta2*(Q40*Q40+Q41*Q41);
  fC32+=theta2*Q20*Q40;
  fC33+=theta2*Q20*Q20;
  fC42+=theta2*Q30*Q40;
  fC43+=theta2*Q30*Q20;
  fC44+=theta2*Q30*Q30;
    
}
void AliITStrack::PrimaryTrack(AliITSRad *rl) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// calculation of part of covariance matrix for vertex  constraint

  Double_t Rlayer[6];
	
  Rlayer[0]=4.; Rlayer[1]=7.;  Rlayer[2]=14.9;  Rlayer[3]=23.8;  
  Rlayer[4]=39.1;  Rlayer[5]=43.6;
  	
  Double_t Cy=fX4/2.;	
  Double_t tgl=(fX1-Zv)*Cy/TMath::ASin(Cy*rtrack);
  Double_t rtrack=1.;	
  fX0=0.;
  fX1=rtrack*tgl; 
  fX2=Dv;
  fX3=tgl; 

  fC00=fC10=fC11=fC20=fC21=fC22=fC30=fC31=fC32=fC33=fC40=fC41=fC42=fC43=0.;

  AddEL(rl,1.,1);
  fLayer=0;
  Int_t i;
  for (i=0; i<6; i++) {
    Propagation(Rlayer[i]);
    fLayer++;
    d2(i)=fC22;
    tgl2(i)=fC33;
    dtgl(i)=fC32; 
    AddMS(rl);    
    AddEL(rl,-1,0); 	   
  } 		
}	
 
  
Int_t AliITStrack::DoNotCross(Double_t rk) const{
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2.;
  return (TMath::Abs((Cy*rk+(1.+Cy*D)*D/rk)/(1.+2.*Cy*D))>=1.)?1:0;
}
  

Double_t AliITStrack::argA(Double_t rk) const {
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2.;
  Double_t arg=(Cy*rk + (1 + Cy*D)*D/rk)/(1.+ 2.*Cy*D);		
  if (TMath::Abs(arg) < 1.) return arg;
  //cout<<"class AliITSTrack: argA out of range !\n";/* getchar();*/
  return (arg>0) ? 0.99999999999 : -0.9999999999;
}
   
Double_t AliITStrack::arga(Double_t rk) const {
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2.;
  Double_t arg=(rk*rk - D*D)/(1.+ 2.*Cy*D);		
  if (arg<0.) {/*cout<<"class AliITSTrack: arga out of range !\n";*/ arg=0.;} 
  return Cy*TMath::Sqrt(arg);
}
   	
Double_t AliITStrack::argB(Double_t rk) const {
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2.;	   
  return (rk*rk - D*D)/(rk*(1.+ 2.*Cy*D)*(1.+ 2.*Cy*D));
}
   
Double_t AliITStrack::argC(Double_t rk) const {
  Double_t C=fX4;
  Double_t D=fX2;
  Double_t Cy=C/2.;		
  return  (1./rk - 2.*Cy*argA(rk)/(1.+ 2.*Cy*D));
}

Double_t AliITStrack::PhiDef(Double_t x, Double_t y){
  Double_t pigre= TMath::Pi();
  Double_t phi=10000.;
  if(y == 0. || x == 0.) {
    if(y == 0. && x == 0.) {
      cout << "  Error in AliITStracking::PhiDef x=0 and y=0 \n"; getchar();
    }
    if(y==0. && x>0.) phi=0.;
    if(y==0. && x<0.) phi=pigre;
    if(x==0 && y>0.) phi=pigre/2.;
    if(x==0 && y<0.) phi=1.5*pigre;   
  }
    else {
      if (x>0. && y>0.) phi=TMath::ATan(y/x);
      if (x<0. && y>0.) phi=pigre+TMath::ATan(y/x);
      if (x<0. && y<0.) phi=pigre+TMath::ATan(y/x); 
      if (x>0. && y<0.) phi=(2.*pigre)+TMath::ATan(y/x);     
    }
  if(phi<0. || phi>(2*pigre)) {
    cout<<" Error on phi in  AliITStracking::PhiDef \n"; getchar();
  }  
  return phi;
}

