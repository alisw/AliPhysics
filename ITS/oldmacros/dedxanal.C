//=========================================
int HMax=80,PMax=6.;
int NStat=15000;
Float_t pmin,pmax;

char tit[]="--------------------------";                                                
TH1F pions=TH1F("Qpi",tit,100,0.5,PMax);   
TH1F *qpi =&pions;     

TH2F qplotP, qplotKa, qplotPi, qplotE;
TH1F *q1plot;
TH2F *qplot;
void dedxhis(Int_t pcod=0,Float_t pmin=0,Float_t pmax=1300)
{
     if(!q1plot)q1plot =  new TH1F("Qhis","Particle charge",100,0.5,PMax);
     if(!qplot){  qplot =  new TH2F("Qtrm","Qtrm vs Pmom",100,0,1300,100,0,13);
     //     TH2F qplotP, qplotKa, qplotPi, qplotE;
     qplot->Copy(qplotP); 
     qplotP.SetMarkerStyle(8); qplotP.SetMarkerColor(kBlack); qplotP.SetMarkerSize(.2);
     qplot->Copy(qplotKa); 
     qplotKa.SetMarkerStyle(8); qplotKa.SetMarkerColor(kRed); qplotKa.SetMarkerSize(.2);
     qplot->Copy(qplotPi); 
     qplotPi.SetMarkerStyle(8); qplotPi.SetMarkerColor(kBlue); qplotPi.SetMarkerSize(.2);
     qplot->Copy(qplotE); 
     qplotE.SetMarkerStyle(8); qplotE.SetMarkerColor(kGreen); qplotE.SetMarkerSize(.2);
     }else{
         qplotP.Reset();
	  qplotPi.Reset(); qplotKa.Reset(); qplotE.Reset();
     }
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
    if(pcod==0) 
      { if(tabv(0)>=1 ){ 
	//qplot->Fill( Pmom,Qtr );
	if(tabv(11)==2212)qplotP.Fill( Pmom,Qtr );
	if(tabv(11)== 321)qplotKa.Fill( Pmom,Qtr );
	if(tabv(11)== 211)qplotPi.Fill( Pmom,Qtr );
	if(tabv(11)==  11)qplotE.Fill( Pmom,Qtr );
      } 
      }
	else { if(tabv(11)==pcod && tabv(0)>=1 )qplot->Fill( Pmom,Qtr );}
	
	//if( tabv(0)>=2 && tabv(11)==pcod )qplot->Fill( Pmom,Qtr );
	//if(tabv(0)>=1)qplot->Fill( Pmom,Qtr );
	if(Pmom>pmin&&Pmom<pmax)q1plot->Fill(Qtr);
    }
// ------------------------------------------------------------
//    c1=new TCanvas("c1","  ",10,10,700,500);
//    pad1 =new TPad("p1","  " ,0  ,  0,0,1,21);
//    pad1 =new TPad("p1","  " ,0  ,  0,.3,1,21);
//    pad2 =new TPad("p2","  " ,.35,  0, 1,1,0);  
//    pad2->Draw();  pad2->cd();
//.....................
    //if(pcod==0){ qplotP.Draw(); qplotKa.Draw("same"); qplotPi.Draw("same"); qplotE.Draw("same"); }

    if(pcod==0){

      //qplotP->GetXaxis()->SetTitleSize(0.05);
      //qplotP->GetYaxis()->SetTitleSize(0.05);
      //gPad->SetFillColor(kWhite);    // b.b.
      gStyle->SetOptStat(0);
      qplotP->SetXTitle("(Mev/c)");
      qplotP->SetYTitle("(mips)");
      qplotP.Draw(); qplotKa.Draw("same"); qplotPi.Draw("same");
      qplotE.Draw("same");
      TText *text = new TText(800.,11.,"Protons");
      text->SetTextSize(0.05);
      text->SetTextColor(kBlack);
      text->Draw();
      TText *text = new TText(800.,10.,"Kaons");
      text->SetTextSize(0.05);
      text->SetTextColor(kRed);
      text->Draw();
      TText *text = new TText(800.,9.,"Pions");
      text->SetTextSize(0.05);
      text->SetTextColor(kBlue);
      text->Draw();
      TText *text = new TText(800.,8.,"Electrons");
      text->SetTextSize(0.05);
      text->SetTextColor(kGreen);
      text->Draw();

    }
    else{
      qplot->Draw();}
    c1->Range(0,0,1300,10);
    gStyle->SetLineColor(kRed);
    gStyle->SetLineWidth(2);
	TLine *lj[3],*lk[3]; 
	for(Int_t j=0;j<3;j++){
		Float_t x1,x2,y1,y2,xx1,xx2,yy1,yy2;
		x1=pid->cut[j+1][0]; x2=pid->cut[j+2][0];
		y1=y2=pid->cut[j+2][2];
	    lj[j]=new TLine(1000*x1,y1,1000*x2,y2);
            //lj[j]->Draw();
	    if(j==0){yy1=10.;}else{yy1=lj[j-1]->GetY1();}
	    yy2=lj[j]->GetY1();
	    xx1=xx2=x1;
	    lk[j]=new TLine(1000*xx1,yy1,1000*xx2,yy2);
            //lk[j]->Draw();
	}
	//Draw pions-kaons cuts.
//..................................	
	TLine *mj[7],*mk[7]; 
	for(Int_t j=0;j<7;j++){
		Float_t x1,x2,y1,y2,xx1,xx2,yy1,yy2;
		x1=pid->cut[j+2][0]; x2=pid->cut[j+3][0];
		y1=y2=pid->cut[j+3][5];
	    mj[j]=new TLine(1000*x1,y1,1000*x2,y2);
            //mj[j]->Draw();
	    if(j==0){yy1=10.;}else{yy1=mj[j-1]->GetY1();}
	    yy2=mj[j]->GetY1();
	    xx1=xx2=x1;
	    mk[j]=new TLine(1000*xx1,yy1,1000*xx2,yy2);
            //mk[j]->Draw();
	}
	//Draw kaons-protons cuts.	
    gStyle->SetLineWidth(1.);
}
 
void qhispi(Float_t pmin=0,Float_t pmax=1300)
{
    char tit[]="--------------------------";
    sprintf(tit,"%s%.1f %.1f","Pmom ",pmin,pmax);
//     TH1F *qpi =  new TH1F("Qpi",tit,100,0.5,PMax);
//     TH1F *nhit =  new TH1F("Nhit",tit,100,0,6);
     TH2F *nhit =  new TH2F("Nhit",tit,100,0,6,100,0,1300);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
	if(Pmom>pmin&&Pmom<pmax&&abs(tabv(11))==211)qpi->Fill(Qtr);
//	if(tabv(11)==2212)nhit->Fill(tabv(0));
	if(tabv(11)==211)nhit->Fill(tabv(0),Pmom);
    }
// ------------------------------------------------------------
    int pistat=qpi->GetEntries();
    //qpi->SetMaximum(pistat+1);
    qpi->SetMaximum(HMax);
    qpi->SetFillColor(42);
    qpi->Fit("gaus"); 
    TF1 *fun=qpi->GetFunction("gaus");
    fun->SetLineWidth(.1);
    qpi->SetLineWidth(.1);
    qpi->Draw();
//---------------------------
    if(0){
	gaus1=new TF1("g1","gaus",.5,2);
	qpi->SetMaximum(-1111);
	gaus1->SetParameter(0, qpi->GetMaximum() );
	qpi->SetMaximum(HMax);
	gaus1->SetParameter(1,1);
	gaus1->SetParameter(2,.12);
	gaus1->SetLineWidth(.1);
	gaus1->Draw("same");
    }
}

void qhiska(Float_t pmin=0,Float_t pmax=1300)
{
     TH1F *qka =  new TH1F("Qka","Particle charge",100,0.5,PMax);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
//	qplot->Fill( Pmom,Qtr );
	if(Pmom>pmin&&Pmom<pmax&&tabv(11)==321)qka->Fill(Qtr);
    }
// ------------------------------------------------------------
    qka->SetFillColor(0);
    qka->SetLineColor(kRed);
//...................
    qka->Fit("gaus","0"); 
    TF1 *fun=qka->GetFunction("gaus");
    fun->SetLineWidth(.1);
    qka->SetLineWidth(.1);
    qka->Draw("same");
//...................
//    qka->Draw("same");
//    gaus_k(pmin,pmax,qka);

}

void qhisp(Float_t pmin=0,Float_t pmax=1300)
{
     TH1F *qpr =  new TH1F("Qpr","Particle charge",100,0.5,PMax);
// -------------------------------------------------------------
    qpr->Clear();
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
//	qplot->Fill( Pmom,Qtr );
	if(Pmom>pmin&&Pmom<pmax&&tabv(11)==2212)qpr->Fill(Qtr);
    }
// ------------------------------------------------------------
    qpr->SetFillColor(16);
//...................
    qpr->Fit("gaus","0"); 
    TF1 *fun=qpr->GetFunction("gaus");
    fun->SetLineWidth(.1);
    qpr->SetLineWidth(.1);
    qpr->Draw("same");
//...................
//    qpr->Draw("same");
//    gaus_p(pmin,pmax,qpr);
}
//---------------------------------------

void gaus_k(Float_t pmin=450,Float_t pmax=470,TH1F *qka){
    Float_t xmean,xsig;
    if(pmin==410&&pmax==470){xmean=pid->cut[5][3]; xsig=pid->cut[5][4];}
    else 
    if(pmin==470&&pmax==530){xmean=pid->cut[6][3]; xsig=pid->cut[6][4];}
    else
    if(pmin==730&&pmax==830){xmean=pid->cut[10][3]; xsig=pid->cut[10][4];}
    else
    if(pmin==830&&pmax==930){xmean=pid->cut[11][3]; xsig=pid->cut[11][4];}
    else
    if(pmin==930&&pmax==1030){xmean=pid->cut[12][3]; xsig=pid->cut[12][4];}
    else{return;}
	gaus1=new TF1("g1","gaus",xmean-3*xsig,xmean+3*xsig);
	qka->SetMaximum(-1111);
	gaus1->SetParameter(0, qka->GetMaximum() );
	qka->SetMaximum(HMax);
	gaus1->SetParameter(1,xmean);
	gaus1->SetParameter(2,xsig);
	gaus1->SetLineWidth(.1);
	gaus1->Draw("same");
}
//---------------------------------------
void gaus_p(Float_t pmin=730,Float_t pmax=830,TH1F *qp){
    Float_t xmean,xsig;
    if(pmin==730&&pmax==830){xmean=pid->cut[10][5]; xsig=pid->cut[10][6];}
    else
    if(pmin==830&&pmax==930){xmean=pid->cut[11][5]; xsig=pid->cut[11][6];}
    else
    if(pmin==930&&pmax==1030){xmean=pid->cut[12][5]; xsig=pid->cut[12][6];}
    else{return;}
	gaus1=new TF1("g1","gaus",xmean-3*xsig,xmean+3*xsig);
	qp->SetMaximum(-1111);
	gaus1->SetParameter(0, qp->GetMaximum() );
	qp->SetMaximum(HMax);
	gaus1->SetParameter(1,xmean);
	gaus1->SetParameter(2,xsig);
	gaus1->SetLineWidth(.1);
	gaus1->Draw("same");
}

//----b.b.---------------------------------
void fitpi(Float_t pmin=0,Float_t pmax=1300,char *tit)
{
  cout<<"fitpi: NStat, PMax, pmin, pmax ="<<NStat<<","<<PMax<<","<<pmin<<","<<pmax<<endl;
     TH1F *q1fit =  new TH1F("Qfit",tit,100,0.5,PMax);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
	if(Pmom>pmin&&Pmom<pmax&&tabv(11)==211)q1fit->Fill(Qtr);
    }
    q1fit->SetFillColor(0);
    q1fit->SetLineColor(kRed);
    q1fit->Draw();
    q1fit->Fit("gaus");
}

//---------------------------------------
void fitka(Float_t pmin=0,Float_t pmax=1300,char *tit)
{
  //TH1F *q1fit =  new TH1F("Qfit",tit,100,0.5,PMax);
     TH1F *q1fit =  new TH1F("Qfit",tit,100,0.5,10.);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
	if(Pmom>pmin&&Pmom<pmax&&tabv(11)==321)q1fit->Fill(Qtr);
    }
    q1fit->SetFillColor(0);
    q1fit->SetLineColor(kRed);
    q1fit->Draw();
    q1fit->Fit("gaus");
}
//---------------------------------------
void fitp(Float_t pmin=0,Float_t pmax=1300,char *tit)
{
     TH1F *q1fit =  new TH1F("Qfit",tit,100,0.5,PMax);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
	if(Pmom>pmin&&Pmom<pmax&&tabv(11)==2212)q1fit->Fill(Qtr);
    }
    q1fit->SetFillColor(0);
    q1fit->Draw();
    q1fit->Fit("gaus");
}
//---------------------------------------
void effall(){
 
      eff(211,211,"",3);
      eff(321,321,"same",2);
     eff(2212,2212,"same",1);
}
//---------------------------------------
void eff(int pc=211,int pc2=211,char *opt="",int color=2)
{
    const Int_t HSize=100;
     TH1F *efpi =  new TH1F("Effpi","Eff of PID",HSize,0.025,1400);
     TH1F *pmom =  new TH1F("Pmom" ,"Eff of PID",HSize,0.025,1400);
     TH1F *heff =  new TH1F("Eff%" ,"Eff of PID",HSize,0.025,1400);
// -------------------------------------------------------------
    Float_t Qtr,Pmom;
    Int_t   Pcode;
    for(Int_t i=0;i<NStat;i++){
	TVector tabv(*( pid->GetVec(i) ));
	    Qtr=tabv(6);Pmom=tabv(10);
	Pcode=  (  pid->GetPcode(Qtr,Pmom/1000.) );
    if(tabv(0)>=4){
	if(Pcode==pc&&tabv(11)==pc2)   efpi->Fill(Pmom);
	if(tabv(11)==pc2)    pmom->Fill(Pmom);
	}
    }
// ------------------------------------------------------------
    for(int i=1;i<=HSize;i++){
	if(  pmom->fArray[i] > 0 )
	heff->fArray[i] = color * ( efpi->fArray[i] )/( pmom->fArray[i] );
    }
    heff->SetLineColor(color);
    if(color==1)heff->SetLineWidth(2)else heff->SetLineWidth(.5);
    heff->Draw(opt);
return;
    pmom->SetFillColor(16);    pmom->Draw(); 
    efpi->SetFillColor(0); //42 
    efpi->Draw("same");
}
//---------------------------------------
Int_t NRange=0;
    Float_t xpmin[]={410,470,730,830,930};
    Float_t xpmax[]={470,530,830,930,1030};
void all(){
    PMax=3; HMax=500; NStat=99000;
  Float_t pmin,pmax;
  pmin=xpmin[NRange];pmax=xpmax[NRange];NRange++;if(NRange==5)NRange=0;
   qhispi(pmin,pmax); qhisp(pmin,pmax); qhiska(pmin,pmax);
}
//---------------------------------------
void fitkall(){
 PMax=3;
 Float_t pmin,pmax;
 //pmin=xpmin[NRange];pmax=xpmax[NRange];NRange++;if(NRange==5)NRange=0;
 pmin=200;pmax=300; // b.b.
 char str[]="                  ";
 sprintf(str,"Kaons %d - %d MeV/c",pmin,pmax);
 gStyle->SetOptFit();
 fitka(pmin,pmax,str);
}

//--b.b.-------------------------------------
void fitpiall(){
 PMax=3;
 NRange=3; // b.b.
 Float_t pmin,pmax;
 cout<<"fitpiall: NRange ="<<NRange<<endl;
 //pmin=xpmin[NRange];pmax=xpmax[NRange];NRange++;if(NRange==5)NRange=0;
 pmin=100;pmax=200; // b.b.
 char str[]="                  ";
 sprintf(str,"Pions %d - %d MeV/c",pmin,pmax);
 gStyle->SetOptFit();
 cout<<"fitpiall: pmin, pmax ="<<pmin<<","<<pmax<<endl;
 fitpi(pmin,pmax,str);
}

//---------------------------------------
void fitpall(){
 PMax=3;
 Float_t pmin,pmax;
 if(NRange==0)NRange=2;
 pmin=xpmin[NRange];pmax=xpmax[NRange];NRange++;if(NRange==5)NRange=2;
 char str[]="                        ";
 sprintf(str,"Protons %d - %d MeV/c",pmin,pmax);
 gStyle->SetOptFit();
 fitp(pmin,pmax,str);
}
//-------------
void newcuts(){
//             

    pid->SetCut(3,.3, 0   , 2.5,  2.5  ,  9,  9,  10   ); //200-300

    pid->SetCut(5,.47, 1  , 0.12,   1.98 , 1.17 ,  2.5,  10   );//410-470
    pid->SetCut(6,.53, 1  , 0.12,   1.73 , 0.15 ,  2.5,  10   );//470-530
    pid->SetCut(7,.59, 0  , 0,      1.18 , 1.125 ,  2.5,  10   );//530-590

    pid->SetCut(8,.65, 0  , 0,   1.18,  1.125 ,  2.3,  10   );//590-650
}
//---------------------------------------
void qhisall(Float_t pmin=0.25,Float_t pmax=.700)
//void qhisall()
{
    qhispi(pmin,pmax);
    qhisp(pmin,pmax);
    qhiska(pmin,pmax);
}
//--------------------------------------
void pcode(){
    TH1F *his =  new TH1F("pcode","tit",100,10,2300);                                      
Float_t Pcod;                                                                       
for(Int_t i=0;i<NStat;i++){                                                             
    TVector tabv(*( pid->GetVec(i) ));                                                  
    Pcod=tabv(11); if(Pcod>0)  his->Fill(Pcod);                            
}                                                                                       
his->SetFillColor(0);his->SetLineColor(kRed);his->Draw();                                                                          
}

void signal(){                                                                               
    TH1F *hisR =  new TH1F("Rec signal","tit",100,0,13);    
    TH1F *hisH =  new TH1F("Hit signal","tit",100,0,13);    
    Float_t xx;                                                                               
    for(Int_t i=0;i<NStat;i++){                                                                 
        TVector tabv(*( pid->GetVec(i) ));                                                      
	    xx=tabv(1); if(tabv(0)>0)  hisR->Fill(xx);                                             
    }                                                                                           
    hisR->SetFillColor(0);hisR->SetLineColor(kRed);hisR->Draw();         
    hisH->SetFillColor(0);hisH->SetLineColor(kBlue);hisH->Draw("same");                                   
}     
//-------------------------------------

void pmom(){                                                                               
    TH1F *his =  new TH1F("pmom","tit",100,0,1000);    
    Float_t xx;                                                                               
    for(Int_t i=0;i<NStat;i++){                                                                 
        TVector tabv(*( pid->GetVec(i) ));                                                      
	    xx=tabv(10); if(tabv(0)>0)  his->Fill(xx);                                             
    }                                                                                           
    his->SetFillColor(0);his->SetLineColor(kRed);his->Draw();                                   
}     

// Print next track
int jtrack=0;
void track(){
    if(jtrack==NStat)jtrack=0;
    for(Int_t i=jtrack;i<NStat;i++){                                                                 
        TVector tabv(*( pid->GetVec(i) ));                                                      
	     if(tabv(0)>0) 
	         {  jtrack=i;break; }                                             
    }                                                                                           
cout<<"Track No "<<jtrack<<endl;
pid->Print(jtrack++);
}

// Fill histogram for track number
//--------------------------------
void tracks(){
   TH1F *his =  new TH1F("tracks","tit",100,0,15000);    
    Float_t xx;                                                                               
    for(Int_t i=0;i<NStat;i++){                                                                 
        TVector tabv(*( pid->GetVec(i) ));                                                      
	     if(tabv(0)>0) 
	         {  his->Fill(i); }                                             
    }                                                                                           
    his->SetFillColor(0);his->SetLineColor(kRed);his->Draw();      
}
// Fill pid table with reconstructed tracks


#ifndef __CINT__
  #include <iostream.h>
  #include <fstream.h>
 
  #include "AliRun.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliITStrackV2.h"
  #include "AliITSclusterV2.h"
 
  #include "TFile.h"
  #include "TTree.h"
  #include "TH1.h"
  #include "TObjArray.h"
  #include "TStyle.h"
  #include "TCanvas.h"
  #include "TLine.h"
  #include "TText.h"
#endif
 
struct GoodTrackITS {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};



void filltab_tracks(){

  cerr<<"Filling of track table...\n";

  Int_t event=0;
Int_t good_tracks_its(GoodTrackITS *gt, const Int_t max, const Int_t event);

  const Int_t MAX=15000;
  Int_t nentr=0; TObjArray tarray(2000);
  {/* Load tracks */
    TFile *tf=TFile::Open("AliITStracksV2.root");
    if (!tf->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 3;}
    char tname[100]; sprintf(tname,"TreeT_ITS_%d",event);
    TTree *tracktree=(TTree*)tf->Get(tname);
    if (!tracktree) {cerr<<"Can't get a tree with ITS tracks !\n"; return 4;}
    TBranch *tbranch=tracktree->GetBranch("tracks");
    nentr=(Int_t)tracktree->GetEntries();
    for (Int_t i=0; i<nentr; i++) {
      AliITStrackV2 *iotrack=new AliITStrackV2;
      tbranch->SetAddress(&iotrack);
      tracktree->GetEvent(i);
      tarray.AddLast(iotrack);
    }
    delete tracktree; //Thanks to Mariana Bondila
    tf->Close();
  }

  /* Generate a list of "good" tracks */
  GoodTrackITS gt[MAX];
  Int_t ngood=0;
  Float_t ConvRadAng=180./TMath::Pi();
  ifstream in("good_tracks_its");
  if (in) {
    cerr<<"Reading good tracks...\n";
    while (in>>gt[ngood].lab>>gt[ngood].code>>
	   gt[ngood].px>>gt[ngood].py>>gt[ngood].pz>>
	   gt[ngood].x >>gt[ngood].y >>gt[ngood].z) {
      ngood++;
      cerr<<ngood<<'\r';
      if (ngood==MAX) {
	cerr<<"Too many good tracks !\n";
	break;
      }
    }
    if (!in.eof()) cerr<<"Read error (good_tracks_its) !\n";
  } else {
    cerr<<"Marking good tracks (this will take a while)...\n";
    ngood=good_tracks_its(gt,MAX,event);
    ofstream out("good_tracks_its");
    if (out) {
      for (Int_t ngd=0; ngd<ngood; ngd++)
	out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
	   <<gt[ngd].px<<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
	   <<gt[ngd].x <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<endl;
    } else cerr<<"Can not open file (good_tracks_its) !\n";
    out.close();
  }
  cerr<<"Number of good tracks : "<<ngood<<endl;
  TH1F *hp=new TH1F("hp","PHI resolution",50,-20.,20.); hp->SetFillColor(4);
  TH1F *hl=new TH1F("hl","LAMBDA resolution",50,-20,20);hl->SetFillColor(4);
  TH1F *hpt=new TH1F("hpt","Relative Pt resolution",30,-10.,10.);
  hpt->SetFillColor(2);
  TH1F *hmpt=new TH1F("hmpt","Transverse impact parameter",30,-300,300);
  hmpt->SetFillColor(6);
  TH1F *hz=new TH1F("hz","Longitudinal impact parameter",30,-300,300);
  //hmpt->SetFillColor(6);

  AliITStrackV2 *trk=(AliITStrackV2*)tarray.UncheckedAt(0);
  Double_t pmin=0.1*(100/0.299792458/0.2/trk->GetConvConst());
  Double_t pmax=6.0+pmin;

  TH1F *hgood=new TH1F("hgood","Good tracks",30,pmin,pmax);
  TH1F *hfound=new TH1F("hfound","Found tracks",30,pmin,pmax);
  TH1F *hfake=new TH1F("hfake","Fake tracks",30,pmin,pmax);
  TH1F *hg=new TH1F("hg","",30,pmin,pmax); //efficiency for good tracks
  hg->SetLineColor(4); hg->SetLineWidth(2);
  TH1F *hf=new TH1F("hf","Efficiency for fake tracks",30,pmin,pmax);
  hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

  TH1F *hptw=new TH1F("hptw","Weghted pt",30,pmax,pmin);
  
  while (ngood--) {
    Int_t lab=gt[ngood].lab, tlab=-1;
    Double_t pxg=gt[ngood].px, pyg=gt[ngood].py, pzg=gt[ngood].pz;
    Double_t ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
    Double_t pg=1000.*TMath::Sqrt(pxg*pxg+pyg*pyg+pzg*pzg); // b.b.


    if (ptg<pmin) continue;

    hgood->Fill(ptg);

    AliITStrackV2 *track=0;
    Int_t j;
    for (j=0; j<nentr; j++) {
      track=(AliITStrackV2*)tarray.UncheckedAt(j);
      tlab=track->GetLabel();
      if (lab==TMath::Abs(tlab)) break;
    }
    if (j==nentr) {
    cerr<<"Track "<<lab<<" was not found !\n";
    continue;
    }
    //track->Propagate(track->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
    track->PropagateTo(3.,0.0028,65.19);

  track->PropagateToVertex();

  if (lab==tlab) hfound->Fill(ptg);
  else { hfake->Fill(ptg); cerr<<lab<<" fake\n";}

  Double_t xv,par[5]; track->GetExternalParameters(xv,par);
  Float_t phi=TMath::ASin(par[2]) + track->GetAlpha();
  if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
  if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
  Float_t lam=TMath::ATan(par[3]);
  Float_t pt_1=TMath::Abs(par[4]);
  Double_t phig=TMath::ATan2(pyg,pxg);
  //hp->Fill((phi - phig)*1000.);

  Double_t lamg=TMath::ATan2(pzg,ptg);
  //hl->Fill((lam - lamg)*1000.);

  Double_t d=10000*track->GetD();
  //hmpt->Fill(d);

  //hptw->Fill(ptg,TMath::Abs(d));

  Double_t z=10000*track->GetZ();
  //hz->Fill(z);



// ------------------------------------------------ b.b. -----
    Int_t nev=0;
    Float_t mom=1000.*(1./(pt_1*TMath::Cos(lam)));
    Float_t dedx=track->GetdEdx();
    //hep->Fill(mom,dedx,1.);
    Int_t pcode=gt[ngood].code;
    pid->SetEdep(10000*nev+ngood,dedx);
    pid->SetPmom(10000*nev+ngood,mom);
    pid->SetPcod(10000*nev+ngood,abs(pcode));
    //cout<<"!!!!! pcode, pmod, dedx ="<<pcode<<","<<mom<<","<<dedx<<endl;
 
  }
 
  pid->Tab();

  Stat_t ng=hgood->GetEntries();
  //cerr<<"Good tracks "<<ng<<endl;
  Stat_t nf=hfound->GetEntries();
  Stat_t nfak=hfake->GetEntries();
  if (ng!=0)
    cerr<<"Integral efficiency is about "<<nf/ng*100.<<" %\n";
  //delete gAlice; //b.b.
  cout<<"Done with nfound(good): "<<nf<<" + nfake: "<<nfak<<" = "<<nf+nfak<<" tracks."<<endl;
// -------------------- b.b. -------------------
 
  //return 0;
}

Int_t good_tracks_its(GoodTrackITS *gt, const Int_t max, const Int_t event) {
  if (gAlice) {delete gAlice; gAlice=0;}

  TFile *file=TFile::Open("galice.root");
  if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}
  if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
    cerr<<"gAlice have not been found on galice.root !\n";
    exit(5);
  }

  Int_t np=gAlice->GetEvent(event);

  Int_t *good=new Int_t[np];
  Int_t k;
  for (k=0; k<np; k++) good[k]=0;

  AliITS *ITS=(AliITS*)gAlice->GetDetector("ITS");
  if (!ITS) {
    cerr<<"can't get ITS !\n"; exit(8);
  }
  AliITSgeom *geom=ITS->GetITSgeom();
  if (!geom) {
    cerr<<"can't get ITS geometry !\n"; exit(9);
  }

  TFile *cf=TFile::Open("AliITSclustersV2.root");
  if (!cf->IsOpen()){
    cerr<<"Can't open AliITSclustersV2.root !\n"; exit(6);
  }
  char cname[100]; sprintf(cname,"TreeC_ITS_%d",event);
  TTree *cTree=(TTree*)cf->Get(cname);
  if (!cTree) {
    cerr<<"Can't get cTree !\n"; exit(7);
  }
  TBranch *branch=cTree->GetBranch("Clusters");
  if (!branch) {
    cerr<<"Can't get clusters branch !\n"; exit(8);
  }
  TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
  branch->SetAddress(&clusters);

  Int_t entr=(Int_t)cTree->GetEntries();
  for (k=0; k<entr; k++) {
    cTree->GetEvent(k);
    Int_t ncl=clusters->GetEntriesFast(); if (ncl==0) continue;
    Int_t lay,lad,det;  geom->GetModuleId(k,lay,lad,det);
    if (lay<1 || lay>6) {
      cerr<<"wrong layer !\n"; exit(10);
    }
    while (ncl--) {
      AliITSclusterV2 *pnt=(AliITSclusterV2*)clusters->UncheckedAt(ncl);
      Int_t l0=pnt->GetLabel(0);
      Int_t l1=pnt->GetLabel(1);
      Int_t l2=pnt->GetLabel(2);
      Int_t mask=1<<(lay-1);
      if (l0>=0) good[l0]|=mask;
      if (l1>=0) good[l1]|=mask;
      if (l2>=0) good[l2]|=mask;
    }
  }
  clusters->Delete(); delete clusters;
  delete cTree; //Thanks to Mariana Bondila
  cf->Close();

  ifstream in("good_tracks_tpc");
  if (!in) {
    cerr<<"can't get good_tracks_tpc !\n"; exit(11);
  }
  Int_t nt=0;
  Double_t px,py,pz,x,y,z;
  Int_t code,lab;
  while (in>>lab>>code>>px>>py>>pz>>x>>y>>z) {
    if (good[lab] != 0x3F) continue;
    TParticle *p = (TParticle*)gAlice->Particle(lab);
    gt[nt].lab=lab;
    gt[nt].code=p->GetPdgCode();
    //**** px py pz - in global coordinate system
    gt[nt].px=p->Px(); gt[nt].py=p->Py(); gt[nt].pz=p->Pz();
    gt[nt].x=gt[nt].y=gt[nt].z=0.;
    nt++;
    if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
  }

  delete[] good;

  delete gAlice; gAlice=0;
  file->Close();

  return nt;
}


//----------------------------------------
void filltab(){
  cout<<"Fill tab..";
  int nev=0;
  int track,pcode;
  float signal,pmom;
  for(int t=0;t<10000;t++){
    track=t; signal=5.; pmom=1200.; pcode=321;
       pid->SetEdep(10000*nev+track,signal);   
       pid->SetPmom(10000*nev+track,pmom);                                     
       pid->SetPcod(10000*nev+track,abs(pcode));
  }
  pid->Tab();
  cout<<"Done."<<endl;
}
//
void loadpid(){
  if(pid==0){                                                                  
            TFile *f=new TFile("pidhit.root");                                       
	    AliITSPid *pid=(AliITSPid*)f->Get("AliITSPid");
  if(pid)                                                                      
	{cout<<"Load PID object from PIDHIT.ROOT file"<<endl;}                   
    else{cout<<"ERROR load PID object from PIDHIT.ROOT file"<<endl;}         
 }                                                                            
 pid->Print(0);       
}

void quit(){
    gROOT->Reset();
    gROOT->ProcessLine(".q");
}
//---------------------------------------
















