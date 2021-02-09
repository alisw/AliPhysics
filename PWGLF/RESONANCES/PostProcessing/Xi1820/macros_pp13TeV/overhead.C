char s0[500],s1[500];//charachter strings used for various purposes
char gDIR[500];//stores the name of the working directory
char color[100];//plot colors stored here
int style;//plot styles are stored here

//int cb,ncb=1,cbmin[20]={0},cbmax[20]={400};//defines multiplicity binning, currently only 1 bin//original code
int cb,ncb=4;//4
//int cb,ncb=1;//2
//float cbmin[20]={0,0},cbmax[20]={100,0.1};//defines multiplicity binning, currently only 1 bin//previous version
//float cbmin[20]={0,0,15,40},cbmax[20]={100,15,40,100};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,0,15,50},cbmax[20]={100,15,50,100};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,0,5,15},cbmax[20]={100,5,15,100};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,0,5,10},cbmax[20]={15,5,10,15};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,0,10,20},cbmax[20]={100,10,20,100};//defines multiplicity binning, currently only 1 bin
float cbmin[20]={0,0,10,30},cbmax[20]={100,10,30,100};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,20,40,60},cbmax[20]={20,40,60,100};//defines multiplicity binning, currently only 1 bin
//float cbmin[20]={0,10,20,30,40,50,60,70},cbmax[20]={10,20,30,40,50,60,70,80};//defines multiplicity binning, currently only 1 bin
///float cbmin[20]={0,0,30,50},cbmax[20]={100,30,50,100};//defines multiplicity binning, currently only 1 bin pPb502


int pb,npb;
float ptmin[100],ptmax[100],ptb[100];//defines pT binning

float massmin[100]={1.6};
float massmax[100]={2.4};

int jDC,nDC=9;//"DC" stands for "decay channel"; keeps track of the different decay channels
int Kpikx=0,Kpik0=1,Kkxk0=2,Kpkx=3,Kpk0=4,Klpi=5,Klkx=6,Klk0=7,Klp=8;//the different channels
char DC0[100],DC1[500];//these strings are filled differently depending on the decay channel

int jBS,nBS=6;//"BS" stands for "background subtraction"; keeps track of the type of combinatorial background being used
int Ksum=0,Kgm=1,Kmix=2;//different types of combinatorial background
char BS0[100],BS1[100],BS2[500];//these strings are filled differently depending on the combinatorial background

int jNR,nNR=1;//"NR" stands for "normalization region"; keeps track of the normalization region used for mixed-event background
double nr1,nr2,nr3,nr4;//stores the values of the normalization region boundaries; it is possible to have a normalization region with a hole in the middle (e.g., surrounding a peak), which is why there are four values for the boundaries

//list of things to declare
bool ptbins(int);
void StyleColorCentrality(int);
bool TitlesDC();
bool TitlesDC2(int);
int GetDC(char*);
bool TitlesBS();
bool TitlesBS2(int);
bool NRbins();
bool NRbins2(int);
double normalize_mix(TH1D*,TH1D*,TH1D*,TH1D*);
double standard_err(double,double,double,double);
double standard_err_sqrt(double,double);
void AKColor(char*,char*);
//list of things to declare

bool ptbins(int m=0){//fills information on pT binning, which can in priciple be different for each mutliplicity bin or decay channel (but is currently all the same)
  npb=3;//5
  ///double temp[100]={0.,0.5,1.,2.,4.,10.};
    //double temp[100]={0.,0.5,1.5,3.,5.,10.};
    ////double temp[100]={0.,1.5,2.5,4.,7.,10.};
    double temp[100]={1.5,2.5,4.,20.};//Current
    /////double temp[100]={1.5,2.5,4.,20.};//Current
    //double temp[100]={1.5,2.5,4.,7.};

    
  for(pb=0;pb<npb;pb++){ptmin[pb]=temp[pb]; ptmax[pb]=temp[pb+1]; ptb[pb]=temp[pb];}
  //ptmin[pb]=0.; ptmax[pb]=30.;
    ptmin[pb]=1.5; ptmax[pb]=20.;//30
  ptb[pb]=ptmax[pb-1];
  if(!m) npb++;

  return true;
}


//void StyleColorCentrality(){StyleColorCentrality(cb);}
void StyleColorCentrality(){//sets styles and colors for different multiplicity bins (not much use right now)
  if(!cb){style=8; AKColor((char*)"black",color); return;}

  cerr<<"unknown value for cb "<<cb<<endl;
  style=1; AKColor((char*)"brown",color);
  return;
}


bool TitlesDC(){return TitlesDC2(jDC);}
bool TitlesDC2(int j=0){//fills information for decay channels
  if(j==Kpikx){
    sprintf(DC0,"pikx"); sprintf(DC1,"X#rightarrow#pi^{#pm}K^{#mp}");
    nNR=3;
    return true;
  }else if(j==Kpik0){
    sprintf(DC0,"pik0"); sprintf(DC1,"X#rightarrow#pi^{#pm}K^{0}_{S}");
    nNR=3;
    return true;
  }else if(j==Kkxk0){
    sprintf(DC0,"kxk0"); sprintf(DC1,"X#rightarrowK^{#pm}K^{0}_{S}");
    nNR=3;
    return true;
  }else if(j==Kpkx){
    sprintf(DC0,"pkx"); sprintf(DC1,"X#rightarrowp(#bar{p})K^{#mp}");
    nNR=3;
    return true;
  }else if(j==Kpk0){
    sprintf(DC0,"pk0"); sprintf(DC1,"X#rightarrowp(#bar{p})K^{0}_{S}");
    nNR=3;
    return true;
  }else if(j==Klpi){
    sprintf(DC0,"Lambdapi"); sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})#pi^{#pm}");
    nNR=2;
    return true;
  }else if(j==Klkx){
    sprintf(DC0,"Lambdakx"); sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})K^{#pm}");
    nNR=4;
      ////nNR=1;
    return true;
  }else if(j==Klk0){
    sprintf(DC0,"Lambdak0"); sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})K^{0}_{S}");
    nNR=4;
      ////nNR=1;
    return true;
  }else if(j==Klp){
    sprintf(DC0,"Lambdap"); sprintf(DC1,"X#rightarrow#Lambdap(#bar{#Lambda}#bar{p})");
    nNR=3;
    return true;
  }

  style=8; AKColor((char*)"black",color);
  return false;
}


int GetDC(char* name){//identifies the decay channel based on the string passed
  int j=-1;
  for(jDC=0;jDC<nDC && TitlesDC();jDC++) if(!strcmp(name,DC0)) break;
  if(jDC<nDC) j=jDC;
  else{cerr<<"unknown decay channel "<<name<<endl; return -1;}
  return j;
}

bool TitlesBS(){return TitlesBS2(jBS);}
bool TitlesBS2(int j=0){//fills information for the combinatorial background
  if(j==Ksum){
    sprintf(BS0,"sum"); sprintf(BS1,"Like Charge (Sum) "); sprintf(BS2,"n(--)+n(++)");
    style=4; AKColor((char*)"dgreen",color);
    return true;
  }else if(j==Kgm){
    sprintf(BS0,"gm"); sprintf(BS1,"Like Charge (2GM) "); sprintf(BS2,"2#sqrt{n(--)n(++)}");
    style=8; AKColor((char*)"blue",color);
    return true;
  }else if(j==Kmix){
    sprintf(BS0,"mix"); sprintf(BS1,"Mixed Event "); sprintf(BS2,"Mixed Event");
    style=29; AKColor((char*)"red",color);
    return true;
  }

  cerr<<"unknown value for jBS "<<j<<endl;
  sprintf(BS0,""); sprintf(BS1,""); sprintf(BS2,""); style=1; AKColor((char*)"brown",color);
  return false;
}


bool NRbins(){return NRbins2(jNR);}
bool NRbins2(int j=0){//fills information for the different normalization regions
  nr1=nr2=nr3=nr4=-1.;
  if(jDC==Kpikx){
    if(!jNR){nr1=1.2; nr2=1.3;}
    else if(jNR==1){nr1=1.55; nr2=1.65;}
    else if(jNR==2){nr1=1.85; nr2=1.9;}
  }else if(jDC==Kpik0){
    if(!jNR){nr1=1.2; nr2=1.3;}
    else if(jNR==1){nr1=1.55; nr2=1.65;}
    else if(jNR==2){nr1=1.85; nr2=1.9;}
  }else if(jDC==Kkxk0){
    if(!jNR){nr1=1.5; nr2=1.6;}
    else if(jNR==1){nr1=2.05; nr2=1.2;}
    else if(jNR==2){nr1=2.6; nr2=2.7;}
  }else if(jDC==Kpkx){
    if(!jNR){nr1=1.55; nr2=1.6;}
    else if(jNR==1){nr1=1.9; nr2=2.;}
    else if(jNR==2){nr1=2.7; nr2=2.8;}
  }else if(jDC==Kpk0){
    if(!jNR){nr1=1.55; nr2=1.6;}
    else if(jNR==1){nr1=1.9; nr2=2.;}
    else if(jNR==2){nr1=2.7; nr2=2.8;}
  }else if(jDC==Klpi){
    if(!jNR){nr1=1.8; nr2=1.9;}
    else if(jNR==1){nr1=2.4; nr2=2.5;}
  }else if(jDC==Klkx){
      /*
    if(!jNR){nr1=1.63; nr2=1.66;}
    else if(jNR==1){nr1=1.85; nr2=1.95;}
    else if(jNR==2){nr1=2.3; nr2=2.4;}
    else if(jNR==3){nr1=2.8; nr2=2.9;}
       */

      if(!jNR){nr1=1.70; nr2=1.74;}
      else if(jNR==1){nr1=1.85; nr2=1.95;}
      else if(jNR==2){nr1=2.0; nr2=2.1;}
      else if(jNR==3){nr1=2.3; nr2=2.4;}

      ////if(!jNR){nr1=2.0; nr2=2.1;}
  }else if(jDC==Klk0){
      /*
    if(!jNR){nr1=1.63; nr2=1.66;}
    else if(jNR==1){nr1=1.85; nr2=1.95;}
    else if(jNR==2){nr1=2.3; nr2=2.4;}
    else if(jNR==3){nr1=2.8; nr2=2.9;}
       */

      if(!jNR){nr1=1.70; nr2=1.74;}
      else if(jNR==1){nr1=1.85; nr2=1.95;}
      else if(jNR==2){nr1=2.0; nr2=2.1;}
      else if(jNR==3){nr1=2.3; nr2=2.4;}
      ////if(!jNR){nr1=2.0; nr2=2.1;}
  }else if(jDC==Klp){
    if(!jNR){nr1=2.3; nr2=2.4;}
    else if(jNR==1){nr1=2.8; nr2=2.9;}
    else if(jNR==2){nr1=3.5; nr2=3.6;}
  }else return false;

  if(nr1<0.){cerr<<"unknown normalization range "<<jNR<<endl; return false;}
  return true;
}


double normalize_mix(TH1D* u,TH1D* b,TH1D* n,TH1D* s){
  //normalizes the mixed-event background (b) to the signal distribution (u) in the chosen normalization region; fills n with the normalized background and s with the background-subtracted distribution 
  n->Reset();
  s->Reset();

  int b1=n->GetXaxis()->FindBin(nr1+1.e-4);
  int b2=n->GetXaxis()->FindBin(nr2-1.e-4);
  int b3=n->GetXaxis()->FindBin(nr3+1.e-4);
  int b4=n->GetXaxis()->FindBin(nr4-1.e-4);

  double iu=0.,ib=0.,nmix,nmix_unc;
  if(b1>=1 && b2>=1){
    iu+=u->Integral(b1,b2);
    ib+=b->Integral(b1,b2);
  }
  if(b3>=1 && b4>=1){
    iu+=u->Integral(b3,b4);
    ib+=b->Integral(b3,b4);
  }

  if(iu<1.e-10){cerr<<"error in normalize_mix(): unlike-charge histogram has 0 integral"<<endl; return 0.;}
  if(ib<1.e-10){cerr<<"error in normalize_mix(): mixed-event histogram has 0 integral"<<endl; return 0.;}

  nmix=ib/iu;
  nmix_unc=standard_err_sqrt(ib,iu);
  n->SetBinContent(0,nmix);
  n->SetBinError(0,nmix_unc);

  int j;
  for(j=1;j<=s->GetNbinsX();j++){
    n->SetBinContent(j,b->GetBinContent(j)/nmix);//generate background
    n->SetBinError(j,b->GetBinError(j)/nmix);
    s->SetBinContent(j,u->GetBinContent(j)-n->GetBinContent(j));
    s->SetBinError(j,sqrt(pow(u->GetBinError(j),2)+pow(n->GetBinError(j),2)));//subtract background from unlike-charge
  }

  return nmix;
}


double standard_err(double A,double B,double EA,double EB){
  double ER=0.;
  if(B>0.) ER=sqrt(EA*EA/(B*B)+A*A*EB*EB/(B*B*B*B));
  return ER;
}


double standard_err_sqrt(double A,double B){
  return standard_err(A,B,sqrt(abs(A)),sqrt(abs(B)));
}


void AKColor(char* name,char* out){
  if(!strcmp(name,"red")) sprintf(out,"#ff0000");
  if(!strcmp(name,"dred")) sprintf(out,"#800000");
  if(!strcmp(name,"orange")) sprintf(out,"#ffaa00");
  if(!strcmp(name,"burntorange")) sprintf(out,"#cc5500");
  if(!strcmp(name,"yellow")) sprintf(out,"#dddd00");
  if(!strcmp(name,"trueyellow")) sprintf(out,"#ffff00");
  if(!strcmp(name,"green") || !strcmp(name,"lgreen")) sprintf(out,"#00dd00");
  if(!strcmp(name,"truegreen")) sprintf(out,"#00ffff");
  if(!strcmp(name,"dgreen")) sprintf(out,"#238e23");
  if(!strcmp(name,"teal") || !strcmp(name,"bluegreen")) sprintf(out,"#008080");
  if(!strcmp(name,"cyan")) sprintf(out,"#00dddd");
  if(!strcmp(name,"truecyan")) sprintf(out,"#00ffff");
  if(!strcmp(name,"blue")) sprintf(out,"#0000ff");
  if(!strcmp(name,"lblue")) sprintf(out,"#aaaaff");
  if(!strcmp(name,"purple")) sprintf(out,"#aa00ff");
  if(!strcmp(name,"magenta")) sprintf(out,"#ff00ff");
  if(!strcmp(name,"pink")) sprintf(out,"#ffc0cb");
  if(!strcmp(name,"rose")) sprintf(out,"#ddaaaa");
  if(!strcmp(name,"brown")) sprintf(out,"#964b00");
  if(!strcmp(name,"black")) sprintf(out,"#000000");
  if(!strcmp(name,"gray")) sprintf(out,"#7f7f7f");
  return;
}
