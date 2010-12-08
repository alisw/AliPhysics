#define AnaEveMyTree_cxx
#include "AnaEveMyTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnaEveMyTree::Loop()
{
//      Root > .L AnaEveMyTree.C
//      Root > AnaEveMyTree t
//      Root > t.Loop();       // Loop on all entries
//
// T0C -3.3 ~ -2.9
// T0A  4.5 ~  5.0
//
// V0C -1.7 ~ -2.2 ~ -2.7 ~ -3.2 ~ -3.7
// V0A  2.8 ~  3.4 ~  3.9 ~  4.5 ~  5.1
//
 float pi=acos(-1.0);
 float t0phi[24];
 for(int it=0; it<12; it++) {
  t0phi[it]=(it+1.0)*2.0*pi/12.0;
  t0phi[it+12]=it*2.0*pi/12.0;
 }
 for(int it=0; it<24; it++) {
  float phi=t0phi[it];
  t0phi[it]=atan2(sin(phi),cos(phi));
 }
 float v0phi[64];
 int is=0;
 for(int iRing = 0; iRing<4 ; iRing++){
  for(int iSec = 0; iSec<8 ; iSec++){
   v0phi[is] = 22.5 + 45. * iSec;
   is++;
  }
 }
 for(int iRing = 0; iRing<4 ; iRing++){
  for(int iSec = 0; iSec<8 ; iSec++){
   v0phi[is] = 22.5 + 45. * iSec;
   is++;
  }
 }
 for (int i=0; i<64; i++) {
  float phi=v0phi[i]*pi/180.0;
  v0phi[i]=atan2(sin(phi),cos(phi));
  cout << v0phi[i] << " ";
  if (i%8==7) cout << endl;
 }
 float z0gai[8],t0gai[24],v0gai[64];
 float mean[10][10][4][27][2];
 float widt[10][10][4][27][2];
 float four[10][10][4][27][2][8];
 float f0,f1,f2,f3,f4,f5,f6,f7;
 ifstream ifs;
 ifs.open("AnaEveMyTree.cal");
 for (int i=0; i<1; i++) {
  ifs>>f0>>f1>>f2>>f3>>f4>>f5>>f6>>f7;
  z0gai[i*8+0]=f0; z0gai[i*8+1]=f1;
  z0gai[i*8+2]=f2; z0gai[i*8+3]=f3;
  z0gai[i*8+4]=f4; z0gai[i*8+5]=f5;
  z0gai[i*8+6]=f6; z0gai[i*8+7]=f7;
  cout << f0 << " " << f1 << " " << f2 << " " << f3 << " "
       << f4 << " " << f5 << " " << f6 << " " << f7 << endl;
 }
 for (int i=0; i<3; i++) {
  ifs>>f0>>f1>>f2>>f3>>f4>>f5>>f6>>f7;
  t0gai[i*8+0]=f0; t0gai[i*8+1]=f1;
  t0gai[i*8+2]=f2; t0gai[i*8+3]=f3;
  t0gai[i*8+4]=f4; t0gai[i*8+5]=f5;
  t0gai[i*8+6]=f6; t0gai[i*8+7]=f7;
  cout << f0 << " " << f1 << " " << f2 << " " << f3 << " "
       << f4 << " " << f5 << " " << f6 << " " << f7 << endl;
 }
 for (int i=0; i<8; i++) {
  ifs>>f0>>f1>>f2>>f3>>f4>>f5>>f6>>f7;
  v0gai[i*8+0]=f0; v0gai[i*8+1]=f1;
  v0gai[i*8+2]=f2; v0gai[i*8+3]=f3;
  v0gai[i*8+4]=f4; v0gai[i*8+5]=f5;
  v0gai[i*8+6]=f6; v0gai[i*8+7]=f7;
  cout << f0 << " " << f1 << " " << f2 << " " << f3 << " "
       << f4 << " " << f5 << " " << f6 << " " << f7 << endl;
 }
 for (int ic=0; ic<10; ic++) {
  for (int iz=0; iz<10; iz++) {
   for (int ih=0; ih<4; ih++) {
    for (int id=0; id<27; id++) {
     for (int ib=0; ib<2; ib++) {
      mean[ic][iz][ih][id][ib]=0.0;
      widt[ic][iz][ih][id][ib]=1.0;
     }
     ifs >> f0 >> f1 >> f2 >> f3;
     if (f1==0.0 || f1<0.0) f1=1.0;
     if (f3==0.0 || f3<0.0) f3=1.0;
     mean[ic][iz][ih][id][0]=f0;
     widt[ic][iz][ih][id][0]=f1;
     mean[ic][iz][ih][id][1]=f2;
     widt[ic][iz][ih][id][1]=f3;
     if (ic==4 && iz==4 && ih==0) {
      cout << f0 << " " << f1 << " " << f2 << " " << f3 << endl;
     }
     for (int ib=0; ib<2; ib++) {
      for (int io=0; io<8; io++) {
       four[ic][iz][ih][id][ib][io]=0.0;
      }
      ifs >> f0 >> f1 >> f2 >> f3 >> f4 >> f5 >> f6 >> f7;
      four[ic][iz][ih][id][ib][0]=f0;
      four[ic][iz][ih][id][ib][1]=f1;
      four[ic][iz][ih][id][ib][2]=f2;
      four[ic][iz][ih][id][ib][3]=f3;
      four[ic][iz][ih][id][ib][4]=f4;
      four[ic][iz][ih][id][ib][5]=f5;
      four[ic][iz][ih][id][ib][6]=f6;
      four[ic][iz][ih][id][ib][7]=f7;
      if (ic==4 && iz==4 && ih==0) {
       cout << f0 << " " << f1 << " " << f2 << " " << f3 << " "
            << f4 << " " << f5 << " " << f6 << " " << f7 << endl;
      }
     }
    }
   }
  }
 }
 ifs.close();
 float av1=0;
 float no1=0;
 for (int i=0; i<8;  i++) {
  if (z0gai[i]>3) {
   av1+=z0gai[i];
   no1+=1.0;
  } else {
   z0gai[i]=-1.0;
  }
 }
 av1/=no1;
 float av2=0;
 float no2=0;
 for (int i=0; i<24;  i++) {
  if (t0gai[i]>3) {
   av2+=t0gai[i];
   no2+=1.0;
  } else {
   t0gai[i]=-1.0;
  }
 }
 av2/=no2;
 float av3=0;
 float no3=0;
 for (int i=0; i<64;  i++) {
  if (v0gai[i]>3) {
   av3+=v0gai[i];
   no3+=1.0;
  } else {
   v0gai[i]=-1.0;
  }
 }
 av3/=no3;
 cout << av1 << " " << av2 << " " << av3 << endl;
 TFile *tf = new TFile("AnaEveMyTree.root","recreate");
 TH2D     *cnzv = new TH2D("cnzv","",50,0.0,11000.0,50,-15.0,15.0);
 TProfile *det0g = new TProfile("det0g","", 16,0.0, 16.0   ,-600.0,6000.0);
 TH2D     *det0m = new TH2D    ("det0m","", 16,0.0, 16.0,50,-300.0,3000.0);
 TProfile *det1g = new TProfile("det1g","", 48,0.0, 48.0   , -20.0, 200.0);
 TH2D     *det1m = new TH2D    ("det1m","", 48,0.0, 48.0,50, -10.0, 100.0);
 TProfile *det2g = new TProfile("det2g","",128,0.0,128.0   ,-100.0,1000.0);
 TH2D     *det2m = new TH2D    ("det2m","",128,0.0,128.0,50, -50.0, 500.0);
 char name[80];
 TH2D     *zdcxy[2][5];
 for (int i=0; i<2; i++) {
  for (int j=0; j<5; j++) {
   sprintf(name,"zdcxy%d%d",i,j);
   zdcxy[i][j] = new TH2D(name,"",50,-2,2,50,-2,2);
  }
 }
 float rng0=11000.0;
 float rng1=av3*64.0*4.5;
 float rng2=av2*12.0*7.0; 
 float rng3=av1* 8.0*2.2;
 TH2D *cor00 = new TH2D("corr00","",50,0.0,rng0,50,0.0,rng1);
 TH2D *cor01 = new TH2D("corr01","",50,0.0,rng0,50,0.0,rng2);
 TH2D *cor02 = new TH2D("corr02","",50,0.0,rng0,50,0.0,rng3);
 TH2D *cor03 = new TH2D("corr03","",50,0.0,rng1,50,0.0,rng2);
 TH2D *cor04 = new TH2D("corr04","",50,0.0,rng1,50,0.0,rng3);
 TH2D *cor05 = new TH2D("corr05","",50,0.0,rng2,50,0.0,rng3);
 TH2D *cor10 = new TH2D("corr10","",50,0.0,rng0,50,0.0,rng1);
 TH2D *cor11 = new TH2D("corr11","",50,0.0,rng0,50,0.0,rng2);
 TH2D *cor12 = new TH2D("corr12","",50,0.0,rng0,50,0.0,rng3);
 TH2D *cor13 = new TH2D("corr13","",50,0.0,rng1,50,0.0,rng2);
 TH2D *cor14 = new TH2D("corr14","",50,0.0,rng1,50,0.0,rng3);
 TH2D *cor15 = new TH2D("corr15","",50,0.0,rng2,50,0.0,rng3);
 TProfile *ave[10][10][4][27];
 TProfile *flt[10][10][4][27];
 for (int ic=0; ic<10; ic++) {
  for (int iz=0; iz<10; iz++) {
   for (int ih=0; ih<4; ih++) {
    for (int id=0; id<27; id++) {
     sprintf(name,"ave%d%d%d%d",ic,iz,ih,id);
     ave[ic][iz][ih][id] = new TProfile(name,name,4,-0.5,3.5,-1.0E+06,1.0E+06,"S");
     sprintf(name,"flt%d%d%d%d",ic,iz,ih,id);
     flt[ic][iz][ih][id] = new TProfile(name,name,32,-0.5,31.5,-1.0E+06,1.0E+06);
    }
   }
  }
 }
 TProfile *res[4][22][2];
 TH2D     *cor[4][22][5];
 for (int ih=0; ih<4; ih++) {
  for (int id=0; id<22; id++) {
   for (int ia=0; ia<2; ia++) {
    sprintf(name,"res%d%d%d",ih,id,ia);
    res[ih][id][ia] = new TProfile(name,name,10,-0.5,9.5,-1.1,1.1);
   }
   for (int ia=0; ia<5; ia++) {
    sprintf(name,"cor%d%d%d",ih,id,ia);
    cor[ih][id][ia] = new TH2D(name,name,50,-pi,pi,50,-pi,pi);
   }
  }
 }
 TH1D *t0d[24][5],*t0e[24][5];
 for (int i=0; i<24; i++) {
  for (int j=0; j<5; j++) {
   if (j<2) {
    sprintf(name,"t0d%d_%d",i,j);
    t0d[i][j] = new TH1D(name,"",24,-pi,pi);
    sprintf(name,"t0e%d_%d",i,j);
    t0e[i][j] = new TH1D(name,"",24,-pi,pi);
   } else {
    sprintf(name,"t0d%d_%d",i,j);
    t0d[i][j] = new TH1D(name,"",24,-pi/2.,pi/2.);
    sprintf(name,"t0e%d_%d",i,j);
    t0e[i][j] = new TH1D(name,"",24,-pi/2.,pi/2.);
   }
  }
 }
 TH2D *tvcr0[10][4],*tvcr1[10][4];
 for (int i=0; i<10; i++) {
  for (int j=0; j<4; j++) {
   sprintf(name,"tvcr0_%d_%d",i,j);
   tvcr0[i][j] = new TH2D(name,"",48,0,48,24,-pi,pi);
   sprintf(name,"tvcr1_%d_%d",i,j);
   tvcr1[i][j] = new TH2D(name,"",128,0,128,24,-pi,pi);
  }
 }
 int neve[10][10];
 float buft[10][10][3][24];
 float bufv[10][10][3][64];
 for (int i=0; i<10; i++) {
  for (int j=0; j<10; j++) {
   neve[i][i]=0;
  }
 }
 if (fChain == 0) return;
 Long64_t nentries = fChain->GetEntriesFast();
 cout << nentries << endl;
 Long64_t nbytes = 0, nb = 0;
 int ievent==0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) break;
  nb = fChain->GetEntry(jentry);   nbytes += nb;
  if (t1s >   1 && t2s >  1
   && v1s >  10 && v2s > 10
   && z1s > 100 && z2s >100
   && zvt > -10 && zvt < 10
   && v2s < 0.8*ntr+1000
   && v2s > 0.8*ntr-1000
   && t1s > 0.02*ntr) {
  cnzv->Fill(ntr,zvt);
  int icen=(int)(ntr/900.0);
  if (icen<0) icen=0;
  if (icen>9) icen=9;
  int izvt=(int)((zvt+10.0)/2.0);
  if (izvt<0) izvt=0;
  if (izvt>9) izvt=9;

  float sumxy[4][27][4];
  for (int ih=0; ih<4; ih++) {
   for (int id=0; id<27; id++) {
    for (int iv=0; iv<4; iv++) {
     sumxy[ih][id][iv]=0;
    }
   }
  }
  float z0co[8];
  float z0sm=0;
  float z0x[8] = {-1,  1, -1,  1,  1, -1,  1, -1};
  float z0y[8] = {-1, -1,  1,  1, -1, -1,  1,  1};
  for (int i=0; i<8; i++) {
   float gain=z0gai[i];
   if (gain<0) gain=0;
   else gain=av1/gain;
   if (i<4) { // for C-side
    det0g->Fill(i+0.5,z1nt[i+1]);
    det0m->Fill(i+0.5,z1nt[i+1]);
    det0g->Fill(i+8.5,z1nt[i+1]*gain);
    det0m->Fill(i+8.5,z1nt[i+1]*gain);
    z0co[i]=z1nt[i+1]*gain;
    z0sm+=z1nt[i+1]*gain;
   } else { // for A-side
    det0g->Fill(i+0.5,z2nt[i-3]);
    det0m->Fill(i+0.5,z2nt[i-3]);
    det0g->Fill(i+8.5,z2nt[i-3]*gain);
    det0m->Fill(i+8.5,z2nt[i-3]*gain);
    z0co[i]=z2nt[i-3]*gain;
    z0sm+=z2nt[i-3]*gain;
   }
   sumxy[0][i/4][0]+=z0co[i]*z0x[i];
   sumxy[0][i/4][1]+=z0co[i]*z0y[i];
   sumxy[0][i/4][2]+=z0co[i];
   float sgn=-1.0*(1-(int)(i/4)); // sign flip for C-side
   sumxy[0][2][0]+=z0co[i]*z0x[i]*sgn;
   sumxy[0][2][1]+=z0co[i]*z0y[i]*sgn;
   sumxy[0][2][2]+=z0co[i];
  }
  for (int i=0; i<3; i++) {
   sumxy[0][i][0]/=sumxy[0][i][2];
   sumxy[0][i][1]/=sumxy[0][i][2];
  }
  zdcxy[0][icen/2]->Fill(sumxy[0][0][0],sumxy[0][1][0]);
  zdcxy[1][icen/2]->Fill(sumxy[0][0][1],sumxy[0][1][1]);
  float t0co[24];
  float t0sm=0;
  for (int i=0; i<24; i++) {
   float gain=t0gai[i];
   if (gain<0) gain=0;
   else gain=av2/gain;
   if (t0tim[i]>100 && t0amp[i]>0) {
    det1g->Fill(i+0.5,t0amp[i]);
    det1m->Fill(i+0.5,t0amp[i]);
    det1g->Fill(i+24.5,t0amp[i]*gain);
    det1m->Fill(i+24.5,t0amp[i]*gain);
    t0co[i]=t0amp[i]*gain;
    t0sm+=t0amp[i]*gain;
   } else {
    t0co[i]=0.0;
   }
  }
  float v0co[64];
  float v0sm=0;
  for (int i=0; i<64; i++) {
   float gain=v0gai[i];
   if (gain<0) gain=0;
   else gain=av3/gain;
   if (v0mul[i]>0) {
    det2g->Fill(i+0.5,v0mul[i]);
    det2m->Fill(i+0.5,v0mul[i]);
    det2g->Fill(i+64.5,v0mul[i]*gain);
    det2m->Fill(i+64.5,v0mul[i]*gain);
    v0co[i]=v0mul[i]*gain;
    v0sm+=v0mul[i]*gain;
   } else {
    v0co[i]=0.0;
   }
  }
  for (int i=0; i<14; i++) {
   int id=(int)(i/12);
   for (int ih=0; ih<4; ih++) {
    float sgn=1;
    if (ih%2==0) sgn=-1.0*(1-id); // sign flip for C-side
    sumxy[ih][3+id][0]+=t0co[i]*cos((ih+1.0)*t0phi[i]);
    sumxy[ih][3+id][1]+=t0co[i]*sin((ih+1.0)*t0phi[i]);
    sumxy[ih][3+id][2]+=t0co[i];
    sumxy[ih][5][0]+=t0co[i]*cos((ih+1.0)*t0phi[i]);
    sumxy[ih][5][1]+=t0co[i]*sin((ih+1.0)*t0phi[i]);
    sumxy[ih][5][2]+=t0co[i];
    sumxy[ih][6][0]+=t0co[i]*cos((ih+1.0)*t0phi[i])*sgn;
    sumxy[ih][6][1]+=t0co[i]*sin((ih+1.0)*t0phi[i])*sgn;
    sumxy[ih][6][2]+=t0co[i];
   }
  }
  for (int i=0; i<64; i++) {
   int id=(int)(i/32);
   int jd=(int)(i/8);
   int kd=jd%4;
   for (int ih=0; ih<4; ih++) {
    float sgn=1;
    if (ih%2==0) sgn=-1.0*(1-id); // sign flip for C-side
    sumxy[ih][7+id][0]+=v0co[i]*cos((ih+1.0)*v0phi[i]);
    sumxy[ih][7+id][1]+=v0co[i]*sin((ih+1.0)*v0phi[i]);
    sumxy[ih][7+id][2]+=v0co[i];
    sumxy[ih][9][0]+=v0co[i]*cos((ih+1.0)*v0phi[i]);
    sumxy[ih][9][1]+=v0co[i]*sin((ih+1.0)*v0phi[i]);
    sumxy[ih][9][2]+=v0co[i];
    sumxy[ih][10][0]+=v0co[i]*cos((ih+1.0)*v0phi[i])*sgn;
    sumxy[ih][10][1]+=v0co[i]*sin((ih+1.0)*v0phi[i])*sgn;
    sumxy[ih][10][2]+=v0co[i];
    sumxy[ih][11+4*kd+id][0]+=v0co[i]*cos((ih+1.0)*v0phi[i]);
    sumxy[ih][11+4*kd+id][1]+=v0co[i]*sin((ih+1.0)*v0phi[i]);
    sumxy[ih][11+4*kd+id][2]+=v0co[i];
    sumxy[ih][13+4*kd][0]+=v0co[i]*cos((ih+1.0)*v0phi[i]);
    sumxy[ih][13+4*kd][1]+=v0co[i]*sin((ih+1.0)*v0phi[i]);
    sumxy[ih][13+4*kd][2]+=v0co[i];
    sumxy[ih][14+4*kd][0]+=v0co[i]*cos((ih+1.0)*v0phi[i])*sgn;
    sumxy[ih][14+4*kd][1]+=v0co[i]*sin((ih+1.0)*v0phi[i])*sgn;
    sumxy[ih][14+4*kd][2]+=v0co[i];
   }
  }
  if (icen>0 && icen<7) {
   for (int i=0; i<24; i++) {
    float amp1=t0co[i]/av2;
    float phi1=t0phi[i];
    for (int j=0; j<64; j++) {
     int id=(int)(j/32);
     float amp2=v0co[j]/av3;
     float phi2=v0phi[j];
     float dphi=phi2-phi1;
     dphi=atan2(sin(dphi),cos(dphi));
     t0d[i][id]->Fill(phi2,amp1*amp2);
     t0e[i][id]->Fill(dphi,amp1*amp2);
    }
   }
  }
  int ieve=neve[icen][izvt];
  for (int i=0; i<24; i++) buft[icen][izvt][ieve][i]=t0co[i]/av2;
  for (int i=0; i<64; i++) bufv[icen][izvt][ieve][i]=v0co[i]/av3;
  neve[icen][izvt]++;
  if (neve[icen][izvt]==3) {
   neve[icen][izvt]=0;
   for (int ie=0; ie<3; ie++) {
    for (int je=0; je<3; je++) {
     int mix=1;
     if (ie==je) mix=0;
     for (int ip=0; ip<24; ip++) {
      int iq=(int)(ip/12);
      float phi1=t0phi[ip];
      float amp1=buft[icen][izvt][ie][ip];
      for (int jp=0; jp<64; jp++) {
       int jq=(int)(jp/32);
       float phi2=v0phi[jp];
       float amp2=bufv[icen][izvt][je][jp];
       float dphi=phi2-phi1;
       dphi=atan2(sin(dphi),cos(dphi));
       tvcr0[icen][iq*2+jq]->Fill(mix*24+ip+0.5,dphi,amp1*amp2);
       tvcr1[icen][iq*2+jq]->Fill(mix*64+jp+0.5,dphi,amp1*amp2);
      }
     }
    }
   }
  }

// r.p. flattening // -----------------------------------
  int calFlag=1;
  for (int ih=0; ih<4; ih++) {
   for (int id=0; id<27; id++) {
    for (int ib=0; ib<2; ib++) {
     if (calFlag>0) ave[icen][izvt][ih][id]->Fill(ib+0.0,sumxy[ih][id][ib]);
     float sxy=sumxy[ih][id][ib];
     float mxy=mean[icen][izvt][ih][id][ib];
     float wxy=widt[icen][izvt][ih][id][ib];
     sumxy[ih][id][ib]=(sxy-mxy)/wxy;
     if (calFlag>0) ave[icen][izvt][ih][id]->Fill(ib+2.0,sumxy[ih][id][ib]);
    }
    if (sumxy[ih][id][2]>0.0) {
     sumxy[ih][id][3]=atan2(sumxy[ih][id][1],sumxy[ih][id][0])/(ih+1.0);
     float psi=sumxy[ih][id][3]*(ih+1.0);
     float dp=0.0;
     for (int io=0; io<8; io++) {
      float cc=cos((io+1.0)*psi);
      float ss=sin((io+1.0)*psi);
      if (calFlag>0) flt[icen][izvt][ih][id]->Fill(io+0.0,cc);
      if (calFlag>0) flt[icen][izvt][ih][id]->Fill(io+8.0,ss);
      float aa=four[icen][izvt][ih][id][0][io]; // mean cos
      float bb=four[icen][izvt][ih][id][1][io]; // mean sin
      dp+=(aa*ss-bb*cc)*2.0/(io+1.0);
     }
     psi+=dp;
     for (int io=0; io<8; io++) {
      float cc=cos((io+1.0)*psi);
      float ss=sin((io+1.0)*psi);
      if (calFlag>0) flt[icen][izvt][ih][id]->Fill(io+16.0,cc);
      if (calFlag>0) flt[icen][izvt][ih][id]->Fill(io+24.0,ss);
     }
     sumxy[ih][id][3]=atan2(sin(psi),cos(psi))/(ih+1.0);
    } else {
     sumxy[ih][id][3]=-9999.9;
    }
   }
  }

  if (icen>0 && icen<7) {
   for (int i=0; i<24; i++) {
    float amp1=t0co[i]/av2;
    float phi1=t0phi[i];
    for (int id=7; id<10; id++) {
     float phi2=sumxy[1][id][3];
     if (phi2>-9000) {
      float dphi=phi2-phi1;
      dphi=atan2(sin(2.0*dphi),cos(2.0*dphi))/2.0;
      t0d[i][id-5]->Fill(phi2,amp1);
      t0e[i][id-5]->Fill(dphi,amp1);
     }
    }
   }
  }

// r.p. correlation // -----------------------------------
  for (int ih=0; ih<4; ih++) {
   for (int id=0; id<22; id++) {
    int id1=0;
    int id2=0;
    int ih1=ih;
    int ih2=ih;
    if (id==0) {id1=0; id2=1; ih1=0; ih2=0;} // z0c-z0c (1st)
    if (id==1) {id1=3; id2=4;}               // t0c-t0a
    if (id==2) {id1=7; id2=8;}               // v0c-v0a
    if (id==3) {id1=2; id2=5; ih1=0;}        // zd1-t0ac  6 for t0ac flip
    if (id==4) {id1=2; id2=9; ih1=0;}        // zd1-v0ac 10 for v0ac flip
    if (id==5) {id1=3; id2=7;}               // t0c-v0c
    if (id==6) {id1=3; id2=8;}               // t0c-v0a
    if (id==7) {id1=4; id2=7;}               // t0a-v0c
    if (id==8) {id1=4; id2=8;}               // t0a-v0a
    if (id==9) {id1=5; id2=9;}               // t0ac-v0ac
    if (id==10) {id1=11; id2=12;}            // v0c0-v0a0
    if (id==11) {id1=15; id2=16;}            // v0c1-v0a1
    if (id==12) {id1=19; id2=20;}            // v0c2-v0a2
    if (id==13) {id1=23; id2=24;}            // v0c3-v0a3
    if (id==14) {id1=3; id2=12;}             // t0c-v0a0
    if (id==15) {id1=3; id2=16;}             // t0c-v0a1
    if (id==16) {id1=3; id2=20;}             // t0c-v0a2
    if (id==17) {id1=3; id2=24;}             // t0c-v0a3
    if (id==18) {id1=4; id2=11;}             // t0a-v0c0
    if (id==19) {id1=4; id2=15;}             // t0a-v0c1
    if (id==20) {id1=4; id2=19;}             // t0a-v0c2
    if (id==21) {id1=4; id2=23;}             // t0a-v0c3
    if (sumxy[ih1][id1][3]>-9000 && sumxy[ih2][id2][3]>-9000) {
     float dph=sumxy[ih1][id1][3]-sumxy[ih2][id2][3];
     res[ih][id][0]->Fill((float)icen,cos((ih+1.0)*dph));
     res[ih][id][1]->Fill((float)icen,sin((ih+1.0)*dph));
     float phi1=sumxy[ih1][id1][3];
     float phi2=sumxy[ih2][id2][3];
     if (ih>ih1) phi1*=(ih+1.0);
     else        phi1*=(ih1+1.0);
     if (ih>ih2) phi2*=(ih+1.0);
     else        phi2*=(ih2+1.0);
     phi1=atan2(sin(phi1),cos(phi1));
     phi2=atan2(sin(phi2),cos(phi2));
     cor[ih][id][icen/2]->Fill(phi1,phi2);
    }
   }
  }

  corr00->Fill(1.0*ntr,v1s+v2s);
  corr01->Fill(1.0*ntr,t1s+t2s);
  corr02->Fill(1.0*ntr,z1s+z2s);
  corr03->Fill(v1s+v2s,t1s+t2s);
  corr04->Fill(v1s+v2s,z1s+z2s);
  corr05->Fill(t1s+t2s,z1s+z2s);
  corr10->Fill(1.0*ntr,v0sm);
  corr11->Fill(1.0*ntr,t0sm);
  corr12->Fill(1.0*ntr,z0sm);
  corr13->Fill(v0sm,t0sm);
  corr14->Fill(v0sm,z0sm);
  corr15->Fill(t0sm,z0sm);
  if (ievent%1000==0) cout << "eve " << jentry << " " << ievent << " : "
                      << jev << " " << iev << " : " << ntr << " : "
                      << z1s+z2s << " " << t1s+t2s << " " << v1s+v2s << " : "
                      << z0sm << " " << t0sm << " " << v0sm << endl;
  ievent++;
 }
 }
 ofstream ofs;
 ofs.open("AnaEveMyTree.cal");
 for (int i=0; i<8; i++) {
  ofs << det0g->GetBinContent(i+1) << " ";
  if (i%8==7) ofs << endl;
  cout << det0g->GetBinContent(i+1) << " ";
  if (i%8==7) cout << endl;
 }
 for (int i=0; i<24; i++) {
  ofs << det1g->GetBinContent(i+1) << " ";
  if (i%8==7) ofs << endl;
  cout << det1g->GetBinContent(i+1) << " ";
  if (i%8==7) cout << endl;
 }
 for (int i=0; i<64; i++) {
  ofs << det2g->GetBinContent(i+1) << " ";
  if (i%8==7) ofs << endl;
  cout << det2g->GetBinContent(i+1) << " ";
  if (i%8==7) cout << endl;
 }
 for (int ic=0; ic<10; ic++) {
  for (int iz=0; iz<10; iz++) {
   for (int ih=0; ih<4; ih++) {
    for (int id=0; id<27; id++) {
     for (int ib=0; ib<2; ib++) {
      ofs << ave[ic][iz][ih][id]->GetBinContent(ib+1) << " ";
      ofs << ave[ic][iz][ih][id]->GetBinError  (ib+1) << " ";
      if (ic==4 && iz==4 && ih==0) {
       cout << ave[ic][iz][ih][id]->GetBinContent(ib+1) << " ";
       cout << ave[ic][iz][ih][id]->GetBinError  (ib+1) << " ";
      }
     }
     ofs << endl;
     if (ic==4 && iz==4 && ih==0) cout << endl;
     for (int ib=0; ib<2; ib++) {
      for (int io=0; io<8; io++) {
       ofs << flt[ic][iz][ih][id]->GetBinContent(ib*8+io+1) << " ";
       if (ic==4 && iz==4 && ih==0) {
        cout << flt[ic][iz][ih][id]->GetBinContent(ib*8+io+1) << " ";
       }
      }
      ofs << endl;
      if (ic==4 && iz==4 && ih==0) cout << endl;
     }
    }
   }
  }
 }
 ofs.close();
 tf->Write();
 tf->Close();
}
