#include <iostream.h>
#include <TMath.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TObjArray.h>

#include "AliITS.h"
#include "AliRun.h"
#include "AliITStrack.h"
#include "AliGenerator.h"


ClassImp(AliITStrack)

AliITStrack::AliITStrack() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// default constructor   
 
  fvTrack.ResizeTo(5);
  fmCovariance = new TMatrix(5,5);
  flistCluster = new TObjArray; 
  fNumClustInTrack =0;
  fChi2=-1;
  flabel =0; 
  fVertex.ResizeTo(3); 
  fErrorVertex.ResizeTo(3);
  fLayer = -1; 
  ClusterInTrack = new TMatrix(6,9);
  for(Int_t i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
                           (*ClusterInTrack)(i,8)=-1.;
  rtrack=0.; 
  //alphaprov=-50.; //provvisorio 
  d2.ResizeTo(6);
  tgl2.ResizeTo(6); 
  dtgl.ResizeTo(6);   
}


 
AliITStrack::AliITStrack(const AliITStrack &cobj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it    

  fvTrack.ResizeTo(5); 
  fmCovariance = new TMatrix(5,5);
  ClusterInTrack = new TMatrix(6,9);
  for(Int_t i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
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
  fvTrack = cobj.fvTrack;
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
  //alphaprov=cobj.alphaprov;  //provvisorio    
 
  *fmCovariance = *cobj.fmCovariance;
 
  *ClusterInTrack = *cobj.ClusterInTrack;
 
  Int_t i;
  for(i=0; i<cobj.flistCluster->GetSize(); i++) 
    flistCluster->AddLast(cobj.flistCluster->At(i));
 
}

AliITStrack::AliITStrack(AliTPCtrack &obj)
{ 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 

  fTPCtrack = &obj;
  fvTrack.ResizeTo(5); 
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
  fmCovariance = new TMatrix(5,5);
  ClusterInTrack = new TMatrix(6,9);
  for(Int_t i=0; i<6; i++) (*ClusterInTrack)(i,6)=(*ClusterInTrack)(i,7)=
                           (*ClusterInTrack)(i,8)=-1.;
  flistCluster = new TObjArray; 
  fNumClustInTrack = 0; 
  LmTPC(); 

}

AliITStrack::~AliITStrack() { 

  //destructor

  if(fmCovariance) delete fmCovariance; 
  if(flistCluster) delete flistCluster; 
  if(ClusterInTrack) delete ClusterInTrack;

}     



void AliITStrack::LmTPC() {

// Transform the TPC state vector from TPC-local to master and build a new state vector ITS-type
// The covariance matrix is also modified accordingly
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it      

  
  TVector tpctrack(5); 
  Double_t Alpha = fTPCtrack->GetAlpha(); 
 
  //printf("LmTPC: Alpha %f\n",Alpha); 

  tpctrack(0) = fTPCtrack->GetY();
  tpctrack(1) = fTPCtrack->GetZ();
  tpctrack(2) = fTPCtrack->GetC();
  tpctrack(3) = fTPCtrack->GetEta();
  tpctrack(4) = fTPCtrack->GetTgl();
  //cout<< " tpctracks = "<<tpctrack(0)<<" "<<tpctrack(1)<<" "<<tpctrack(2)<<" "<<tpctrack(3)<<" "<<tpctrack(4)<<
  //"\n ";

  Double_t xm, ym, zm;
  Double_t sina = TMath::Sin(Alpha);
  Double_t cosa = TMath::Cos(Alpha);
  Double_t xl= fTPCtrack->GetX();
  xm = xl * cosa - tpctrack(0)*sina;
  ym = xl * sina + tpctrack(0)*cosa;
  zm = tpctrack(1);  
  //cout<<" xl e alpha = "<<xl<<" "<<Alpha<<"\n"; getchar(); 

  Double_t x0m,y0m;

  ///////////////////////////////////// determine yo //////////////////////////////////////////////////
  Double_t Vxl=fVertex(0)*cosa+fVertex(1)*sina;
  Double_t Vyl= -fVertex(0)*sina+fVertex(1)*cosa;
  Double_t Xo,Yo, signy;
  Double_t R = 1./tpctrack(2);
  Xo =  tpctrack(3) / tpctrack(2);
  xoTPC=Xo;
  Double_t Yo1, Yo2, diffsq1, diffsq2;  
  Yo1 = tpctrack(0) +  TMath::Sqrt(R*R - (xl-Xo)*(xl-Xo)); 
  Yo2 = tpctrack(0) -  TMath::Sqrt(R*R - (xl-Xo)*(xl-Xo));   
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

  fvTrack(0)= phi;
  fvTrack(1)=zm;
  fvTrack(3)=tpctrack(4);
  fvTrack(4)=tpctrack(2);
  /*
  //provvisorio
   alphaprov=pigre-PhiDef(x0m,y0m)-PhiDef(fVertex(0),fVertex(1));
  if(alphaprov<0.) alphaprov=alphaprov+2.*pigre;
  //cout<<" PhiDef(x0m,y0m) ="<<PhiDef(x0m,y0m)<<"\n";
  //cout<<" PhiDef(fVertex(0),fVertex(1)) = "<<PhiDef(fVertex(0),fVertex(1))<<"\n";
  //cout<<"alphaprov = "<<alphaprov<<"\n"; getchar();
  //cout<<" fvTrack(0) prima = "<<fvTrack(0)<<"\n";
  fvTrack(0)-=alphaprov;
  if(fvTrack(0)<0.) fvTrack(0)+= 2.*pigre;
  //cout<<" fvTrack(0) dopo = "<<fvTrack(0)<<"  alphaprov = "<<alphaprov<<"\n"; 
  cout<<" fvertex = "<<fVertex(0)<<" "<<fVertex(1)<<" "<<fVertex(2)<<"\n";   
  fVertex(0)=0.;
  fVertex(1)=0.;
/////////////////////////////////////////////////////////////////////////////////////////////////// 
*/ 
  
  Double_t dd=TMath::Sqrt((x0m-fVertex(0))*(x0m-fVertex(0))+(y0m-fVertex(1))*(y0m-fVertex(1)));
  Double_t signdd;
  if (R>0) signdd=1.; else signdd=-1.;
  fvTrack(2)=signdd*dd-R;
  //cout<<" fvertex = "<<fVertex(0)<<" "<<fVertex(1)<<" "<<fVertex(2)<<"\n";
  //cout<<" fvTrack(2) ="<<fvTrack(2)<<"\n"; getchar();
  
  TMatrix localM(5,5);    
  Double_t cov[15];
  fTPCtrack->GetCovariance(cov);
  localM(0,0)=cov[0];
  localM(1,1)=cov[2];
  localM(2,2)=cov[5];
  localM(3,3)=cov[9];
  localM(4,4)=cov[14];
  localM(1,0)=localM(0,1)=cov[1];
  localM(2,0)=localM(0,2)=cov[3];
  localM(2,1)=localM(1,2)=cov[4];
  localM(3,0)=localM(0,3)=cov[6];
  localM(3,1)=localM(1,3)=cov[7];
  localM(3,2)=localM(2,3)=cov[8];
  localM(4,0)=localM(0,4)=cov[10];
  localM(4,1)=localM(1,4)=cov[11];
  localM(4,2)=localM(2,4)=cov[12];
  localM(4,3)=localM(3,4)=cov[13];

  TMatrix F(5,5);
  
  Double_t  dfidy, dDdy, dDdC, dDdeta;

  dfidy=(xm*cosa+ym*sina)/(rtrack*rtrack);
  dDdy=signdd*((y0m-fVertex(1))*cosa-(x0m-fVertex(0))*sina)/dd;
  Double_t dyodr=signy*(R+(xl-Xo)*tpctrack(3))/TMath::Sqrt(R*R-(xl-Xo)*(xl-Xo));
  Double_t dyomdr=sina*tpctrack(3)+cosa*dyodr;
  Double_t dxomdr=cosa*tpctrack(3)-sina*dyodr;
  Double_t ddddR=((x0m-fVertex(0))*dxomdr+(y0m-fVertex(1))*dyomdr)/dd;
  dDdC=-R*R*(signdd*ddddR-1.);
  Double_t dyoldxol=signy*(xl-Xo)/TMath::Sqrt(R*R-(xl-Xo)*(xl-Xo));
  Double_t dxomdeta=R*(cosa-sina*dyoldxol);
  Double_t dyomdeta=R*(sina+cosa*dyoldxol);
  dDdeta=signdd*((x0m-fVertex(0))*dxomdeta+(y0m-fVertex(1))*dyomdeta)/dd;
  F(0,0)=dfidy;
  F(0,1)=F(0,2)=F(0,3)=F(0,4)=0.;
  F(1,1)=1.;
  F(1,0)=F(1,2)=F(1,3)=F(1,4)=0.;
  F(2,0)=dDdy;    
  F(2,2)=dDdC;    
  F(2,3)=dDdeta;  
  F(2,1)=F(2,4)=0.;
  F(3,0)=F(3,1)=F(3,2)=F(3,3)=0.;
  F(3,4)=1.;
  F(4,0)=F(4,1)=F(4,3)=F(4,4)=0.;
  F(4,2)=1.;
// cout<<" Matrice F \n";  
// F.Print(); getchar();
  TMatrix Ft(TMatrix::kTransposed,F);
  TMatrix temp(localM);
  temp*=Ft;
  TMatrix B(F);
  B*=temp;                            // B=F*localM*Ft 
// cout<<" Matrice C\n";  
// B.Print(); getchar();
  *fmCovariance = B;   
}


AliITStrack &AliITStrack::operator=(AliITStrack obj) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  

  delete fmCovariance;
  delete flistCluster;
  delete ClusterInTrack;
  fmCovariance = new TMatrix(5,5);
  ClusterInTrack = new TMatrix(6,9);
  flistCluster = new TObjArray; 
  flabel = obj.flabel;
  fTPCtrack = obj.fTPCtrack;
  fvTrack.ResizeTo(5); 
  fNumClustInTrack = obj.fNumClustInTrack; 
  fChi2= obj.fChi2;
  fVertex=obj.fVertex;
  fErrorVertex=obj.fErrorVertex;
  fvTrack = obj.fvTrack; 
  fLayer=obj.fLayer;
  rtrack=obj.rtrack;
  Dv=obj.Dv;
  Zv=obj.Zv;
  sigmaDv=obj.sigmaDv;
  sigmaZv=obj.sigmaZv;
  d2=obj.d2;
  tgl2=obj.tgl2; 
  dtgl=obj.dtgl; 
  //alphaprov=obj.alphaprov; //proviisorio   
    
  *fmCovariance = *obj.fmCovariance;
  *ClusterInTrack = *obj.ClusterInTrack;
  for(Int_t i=0; i<obj.flistCluster->GetSize(); i++) flistCluster->AddLast(obj.flistCluster->At(i));

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
  for(Int_t i=0;i<6; i++){
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
      if(VecTotLabref(i) == 0.) VecTotLabref(i)=-3.;  
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
  fvTrack(0)+=TMath::ASin(Ak)-TMath::ASin(Akm1);
  if(fvTrack(0)>duepi) fvTrack(0)-=duepi;
  if(fvTrack(0)<0.) fvTrack(0)+=duepi;
	
  Double_t tgl=fvTrack(3);
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2;
  fvTrack(1)+=tgl/Cy*(TMath::ASin(ak)-TMath::ASin(akm1));
  rtrack=rk;
//cout<<"fvTrack(0) fvTrack(1) e rtrack dentro Propagation = "<<fvTrack(0)<<" "<<fvTrack(1)<<" "<<rtrack<<"\n";
//getchar();
	
  Double_t Bk=argB(rk), Bkm1=argB(rkm1);	
  Double_t Ck=argC(rk), Ckm1=argC(rkm1);
  TMatrix F(5,5);
  F(0,2)=Ck/TMath::Sqrt(1.-Ak*Ak) - Ckm1/TMath::Sqrt(1.-Akm1*Akm1);
  F(0,4)=Bk/TMath::Sqrt(1.-Ak*Ak) - Bkm1/TMath::Sqrt(1.-Akm1*Akm1);	
  F(1,2)=tgl*D*(1./rk - 1./rkm1);
  F(1,3) = rk - rkm1;
  F(0,0)=F(1,1)=F(2,2)=F(3,3)=F(4,4)=1.;	
	
  TMatrix Ft(TMatrix::kTransposed,F);
  TMatrix temp(*fmCovariance);
//cout<<" C Matrix prima propagation =\n";
//temp.Print(); getchar();         
  temp*=Ft;
  TMatrix B(F);
  B*=temp;
//cout<<" C Matrix dopo propagation =\n";                            // B=F*C*Ft 
//B.Print(); getchar();
  *fmCovariance = B;    
  	
}
void AliITStrack::AddEL(Double_t signdE, Bool_t flagtot, Double_t mass) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
//  add energy loss

  TVector s(6);  
  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
//0.00277 is added in the first layer to take into account the energy loss in the beam pipe 
    
  Double_t cl=1.+fvTrack(3)*fvTrack(3);  // cl=1/(cosl)**2 = 1 + (tgl)**2 
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
  if((fvTrack)(4) < 0.) sign=-1.; 
  pt=sign*p/sqcl;  
  Double_t CC=(0.3*0.2)/(pt*100.);
  fvTrack(4)=CC; 
    		
}

void  AliITStrack::Correct(Double_t rk) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// correct track to take into account real geometry detector

  Double_t duepi=2.*TMath::Pi();
  Double_t rkm1=rtrack;
  Double_t Ak=argA(rk), Akm1=argA(rkm1);
  Double_t ak=arga(rk), akm1=arga(rkm1);

  fvTrack(0)+=TMath::ASin(Ak)-TMath::ASin(Akm1);
  if(fvTrack(0)>duepi) fvTrack(0)-=duepi;
  if(fvTrack(0)<0.) fvTrack(0)+=duepi;
	
  Double_t tgl=fvTrack(3);
  Double_t C=fvTrack(4);
  Double_t Cy=C/2;
  fvTrack(1)+=tgl/Cy*(TMath::ASin(ak)-TMath::ASin(akm1));
  rtrack=rk;
  		
}

void AliITStrack::AddMS() {
       
//////////   Modification of the covariance matrix to take into account multiple scattering  ///////////
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it     
   
  TVector s(6);

  s(0)=0.0026+0.00283; s(1)=0.018; s(2)=0.0094; s(3)=0.0095; s(4)=0.0091; s(5)=0.0087;
//0.00277 is added in the first layer to take into account the energy loss in the beam pipe     
  Double_t mass=0.1396;
  Int_t layer=(Int_t)GetLayer();
  
  Double_t tgl=fvTrack(3);
  Double_t cosl=TMath::Cos(TMath::ATan(tgl));	
  Double_t D=fvTrack(2);
  Double_t C=fvTrack(4);
  Double_t Cy=C/2.;
  Double_t Q20=1./(cosl*cosl);	 
  Double_t Q30=C*tgl;
   
  Double_t Q40=Cy*(rtrack*rtrack-D*D)/(1.+ 2.*Cy*D);
  Double_t dd=D+Cy*D*D-Cy*rtrack*rtrack;
  Double_t Q41=-1./cosl*TMath::Sqrt(rtrack*rtrack - dd*dd)/(1.+ 2.*Cy*D);
  
  /*
  Double_t xk=rtrack*TMath::Cos(fvTrack(0));
  Double_t yk=rtrack*TMath::Sin(fvTrack(0));  
  Double_t rvertex=TMath::Sqrt((xk-fVertex(0))*(xk-fVertex(0)) + (yk-fVertex(1))*(yk-fVertex(1)));
  Double_t Q40=Cy*(rvertex*rvertex-D*D)/(1.+ 2.*Cy*D);
  Double_t dd=D+Cy*D*D-Cy*rvertex*rvertex;
  Double_t Q41=-1./cosl*TMath::Sqrt(rvertex*rvertex - dd*dd)/(1.+ 2.*Cy*D); 
  */
   
  TMatrix J(5,2);	
  J(0,0)=0.; J(0,1)=0.;
  J(1,0)=0.; J(1,1)=0.;
  J(2,0)=Q40;    
  J(2,1)=Q41;   	
  J(3,0)=Q20;  J(3,1)=0.;
  J(4,0)=Q30;   J(4,1)=0.;
		 	
  Double_t p2=(GetPt()*GetPt())/(cosl*cosl);
  Double_t beta2=p2/(p2+mass*mass);
  Double_t theta2=14.1*14.1/(beta2*p2*1.e6)*(s(layer-1)/cosl);	
	
  TMatrix Jt(TMatrix::kTransposed,J);	
  TMatrix Q(J,TMatrix::kMult,Jt);
  Q*=theta2;	
		    
  (*fmCovariance)+=Q;
  
}

void AliITStrack::PrimaryTrack() {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it  
// calculation of part of covariance matrix for vertex  constraint

  Double_t Rlayer[6];
	
  Rlayer[0]=4.; Rlayer[1]=7.;  Rlayer[2]=14.9;  Rlayer[3]=23.8;  
  Rlayer[4]=39.1;  Rlayer[5]=43.6;
  	
  Double_t Cy=fvTrack(4)/2.;	
  Double_t tgl=(fvTrack(1)-Zv)*Cy/TMath::ASin(Cy*rtrack);
  Double_t rtrack=1.;	
  fvTrack(0)=0.;
  fvTrack(1)=rtrack*tgl; 
  fvTrack(2)=Dv;
  fvTrack(3)=tgl; 
	
  TMatrix  newC(5,5);
  newC(4,4)=(*fmCovariance)(4,4);
  (*fmCovariance)=newC;
  AddEL(1.,1);
  fLayer=0;
  for (Int_t i=0; i<6; i++) {
    Propagation(Rlayer[i]);
    fLayer++;
    d2(i)=(*fmCovariance)(2,2);
    tgl2(i)=(*fmCovariance)(3,3);
    dtgl(i)=(*fmCovariance)(2,3); 
    AddMS();    
    AddEL(-1,0); 	   
  } 		
}	
 
  
Int_t AliITStrack::DoNotCross(Double_t rk) const{
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2.;
  return (TMath::Abs((Cy*rk+(1.+Cy*D)*D/rk)/(1.+2.*Cy*D))>=1.)?1:0;
}
  

Double_t AliITStrack::argA(Double_t rk) const {
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2.;
  Double_t arg=(Cy*rk + (1 + Cy*D)*D/rk)/(1.+ 2.*Cy*D);		
  if (TMath::Abs(arg) < 1.) return arg;
  //cout<<"class AliITSTrack: argA out of range !\n";/* getchar();*/
  return (arg>0) ? 1. : -1.;
}
   
Double_t AliITStrack::arga(Double_t rk) const {
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2.;
  Double_t arg=(rk*rk - D*D)/(1.+ 2.*Cy*D);		
  if (arg<0.) {/*cout<<"class AliITSTrack: arga out of range !\n";*/ arg=0.;} 
  return Cy*TMath::Sqrt(arg);
}
   	
Double_t AliITStrack::argB(Double_t rk) const {
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2.;	   
  return (rk*rk - D*D)/(rk*(1.+ 2.*Cy*D)*(1.+ 2.*Cy*D));
}
   
Double_t AliITStrack::argC(Double_t rk) const {
  Double_t C=fvTrack(4);
  Double_t D=fvTrack(2);
  Double_t Cy=C/2.;		
  return  (1./rk - 2.*Cy*argA(rk)/(1.+ 2.*Cy*D));
}
/*
Double_t AliITStrack::PhiDef(Double_t x, Double_t y){
  Double_t pigre= TMath::Pi();
  Double_t phi;
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
*/
