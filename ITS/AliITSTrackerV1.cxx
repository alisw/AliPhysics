#include <iostream.h>
#include <TROOT.h>
#include <TMath.h>
#include <TRandom.h>
#include <TBranch.h>
#include <TVector.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>


#include "TParticle.h"
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSgeom.h"
#include "AliITSmodule.h"
#include "AliITSRecPoint.h"
#include "AliMC.h"
#include "AliKalmanTrack.h" 
#include "AliMagF.h"


#include "AliITSRad.h" 
#include "AliITStrack.h"
#include "AliITSiotrack.h"
#include "AliITStracking.h"   
#include "../TPC/AliTPC.h"
#include "../TPC/AliTPCParam.h"
#include "../TPC/AliTPCtracker.h"

#include "AliITSgeoinfo.h"
#include "AliITSTrackerV1.h"



ClassImp(AliITSTrackerV1)


//________________________________________________________________

AliITSTrackerV1::AliITSTrackerV1(AliITS* IITTSS) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// default constructor   
  ITS = IITTSS;

}

//________________________________________________________________

AliITStrack  AliITSTrackerV1::Tracking(AliITStrack &track, AliITStrack *reference,TObjArray *fastpoints, 
Int_t **vettid, Bool_t flagvert,  AliITSRad *rl, AliITSgeoinfo *geoinfo) { 
										
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
  										    
  TList *list= new TList();   

  AliITStrack tr(track);
  
  list->AddLast(&tr);
  
  Double_t Pt=(tr).GetPt();
  //cout << "\n Pt = " << Pt <<"\n";  //stampa

  AliITStracking obj(list, reference, ITS, fastpoints,TMath::Abs(Pt),vettid, flagvert, rl, geoinfo);  // nuova ITS 
  list->Delete();
  delete list;

  Int_t itot=-1;
  TVector VecTotLabref(18);
  Int_t lay, k;
  for(lay=5; lay>=0; lay--) {
    TVector VecLabref(3); 
    VecLabref=(*reference).GetLabTrack(lay);
    Float_t ClustZ=(*reference).GetZclusterTrack( lay);   
    for(k=0; k<3; k++){ 
		Int_t lpp=(Int_t)VecLabref(k);
		if(lpp>=0) {
		  TParticle *p=(TParticle*) gAlice->Particle(lpp);
		  Int_t pcode=p->GetPdgCode();
		  if(pcode==11) VecLabref(k)=p->GetFirstMother();
		}    
    itot++; VecTotLabref(itot)=VecLabref(k);
    if(VecLabref(k)==0. && ClustZ == 0.) VecTotLabref(itot) =-3.; }  
  }
  Long_t labref;
  Int_t freq;  
  (*reference).Search(VecTotLabref, labref, freq);
    
  //if(freq < 6) labref=-labref;        // cinque - sei
  if(freq < 5) labref=-labref;        // cinque - sei	
  (*reference).SetLabel(labref);

  return *reference; 

}



//________________________________________________________________



void AliITSTrackerV1::DoTracking(Int_t evNumber, Int_t min_t, Int_t max_t, TFile *file, Bool_t flagvert) {
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
//   ex macro for tracking ITS

  //Loop variables

  Int_t i;

  printf("begin DoTracking - file %p\n",file);
  
///////////////////////////////////////  gets information on geometry ///////////////////////////////////  
  AliITSgeoinfo *geoinfo = new AliITSgeoinfo;

  AliITSgeom *g1 = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom(); 
  
  Int_t ll=1, dd=1;
  TVector det(9);
  
  //cout<<" nlad ed ndet \n";
  for(i=0; i<6; i++) {
    geoinfo->Nlad[i]=g1->GetNladders(i+1);
    geoinfo->Ndet[i]=g1->GetNdetectors(i+1);
	 //cout<<geoinfo->Nlad[i]<<" "<<geoinfo->Ndet[i]<<"\n"; 
  }
  //getchar();

  //cout<<" raggio medio = ";
  for(i=0; i<6; i++) {  
    g1->GetCenterThetaPhi(i+1,ll,dd,det);
    geoinfo->Avrad[i]=TMath::Sqrt(det(0)*det(0)+det(1)*det(1));
	 //cout<<geoinfo->Avrad[i]<<" ";
  }
  //cout<<"\n"; getchar();
  
  geoinfo->Detx[0] = ((AliITSgeomSPD*)(g1->GetShape(1, ll, dd)))->GetDx();
  geoinfo->Detz[0] = ((AliITSgeomSPD*)(g1->GetShape(1, ll, dd)))->GetDz();
  
  geoinfo->Detx[1] = ((AliITSgeomSPD*)(g1->GetShape(2, ll, dd)))->GetDx();
  geoinfo->Detz[1] = ((AliITSgeomSPD*)(g1->GetShape(2, ll, dd)))->GetDz();
  
  geoinfo->Detx[2] = ((AliITSgeomSDD*)(g1->GetShape(3, ll, dd)))->GetDx();
  geoinfo->Detz[2] = ((AliITSgeomSDD*)(g1->GetShape(3, ll, dd)))->GetDz();
  
  geoinfo->Detx[3] = ((AliITSgeomSDD*)(g1->GetShape(4, ll, dd)))->GetDx();
  geoinfo->Detz[3] = ((AliITSgeomSDD*)(g1->GetShape(4, ll, dd)))->GetDz();
  
  geoinfo->Detx[4] = ((AliITSgeomSSD*)(g1->GetShape(5, ll, dd)))->GetDx();
  geoinfo->Detz[4] = ((AliITSgeomSSD*)(g1->GetShape(5, ll, dd)))->GetDz();
  
  geoinfo->Detx[5] = ((AliITSgeomSSD*)(g1->GetShape(6, ll, dd)))->GetDx();
  geoinfo->Detz[5] = ((AliITSgeomSSD*)(g1->GetShape(6, ll, dd)))->GetDz();
  
  //cout<<"    Detx     Detz\n";
  //for(Int_t la=0; la<6; la++) cout<<"    "<<geoinfo->Detx[la]<<"     "<<geoinfo->Detz[la]<<"\n";
  //getchar();
//////////////////////////////////////////////////////////////////////////////////////////////////////////  
  
  //const char *pname="75x40_100x60";
  
 Int_t imax=200,jmax=450;
  AliITSRad *rl = new AliITSRad(imax,jmax);
  //cout<<" dopo costruttore AliITSRad\n"; getchar();
    
  struct GoodTrack {
    Int_t lab,code;
    Float_t px,py,pz,x,y,z,pxg,pyg,pzg,ptg;
    Bool_t flag;
  };
  

  gAlice->GetEvent(0);
 
    AliKalmanTrack *kkprov;
    kkprov->SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());  

  TFile *cf=TFile::Open("AliTPCclusters.root");  
  AliTPCParam *digp= (AliTPCParam*)cf->Get("75x40_100x60");
  if (!digp) { cerr<<"TPC parameters have not been found !\n"; getchar();}
  
   AliTPCtracker *tracker = new AliTPCtracker(digp);  
 
 // Load clusters
   tracker->LoadInnerSectors();
   tracker->LoadOuterSectors();
       
  
  GoodTrack gt[15000];
  Int_t ngood=0;
  ifstream in("itsgood_tracks");

  cerr<<"Reading itsgood tracks...\n";
  while (in>>gt[ngood].lab>>gt[ngood].code
	  >>gt[ngood].px >>gt[ngood].py>>gt[ngood].pz
	  >>gt[ngood].x  >>gt[ngood].y >>gt[ngood].z
	  >>gt[ngood].pxg  >>gt[ngood].pyg >>gt[ngood].pzg
	  >>gt[ngood].ptg >>gt[ngood].flag) {
    ngood++;
    cerr<<ngood<<'\r';
    if (ngood==15000) {
      cerr<<"Too many good tracks !\n";
      break;
    }
  }
  if (!in.eof()) cerr<<"Read error (itsgood_tracks) !\n";
  
  
// Load tracks
  TFile *tf=TFile::Open("AliTPCtracks.root"); 
  if (!tf->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return ;}
  TObjArray tracks(200000);
   TTree *tracktree=(TTree*)tf->Get("TPCf"); 
   if (!tracktree) {cerr<<"Can't get a tree with TPC tracks !\n";}   
  TBranch *tbranch=tracktree->GetBranch("tracks");
  Int_t nentr=(Int_t)tracktree->GetEntries();
  Int_t kk;
   AliTPCtrack *iotracktpc=0;    
  for (kk=0; kk<nentr; kk++) {
    iotracktpc=new AliTPCtrack; 
    tbranch->SetAddress(&iotracktpc);
    tracktree->GetEvent(kk);    
    tracker->CookLabel(iotracktpc,0.1);       
    tracks.AddLast(iotracktpc);         
  }  
   delete tracker;      
  tf->Close();


  Int_t nt = tracks.GetEntriesFast();
  cerr<<"Number of found tracks "<<nt<<endl;
  
  TVector DataOut(9);
  Int_t kkk=0;
  
  Double_t ptg=0.,pxg=0.,pyg=0.,pzg=0.;

  //////////////////////////////  good tracks definition in TPC  ////////////////////////////////
      
  ofstream out1 ("AliITSTrag.out");
  for (i=0; i<ngood; i++) out1 << gt[i].ptg << "\n";
  out1.close();


  TVector vec(5);
  TTree *TR=gAlice->TreeR();
  Int_t nent=(Int_t)TR->GetEntries();  
  TClonesArray  *recPoints = ITS->RecPoints();  
  
  Int_t numbpoints;
  Int_t totalpoints=0;
  Int_t *np = new Int_t[nent];
  Int_t **vettid = new Int_t* [nent];
  Int_t mod;
  
  for (mod=0; mod<nent; mod++) {
    vettid[mod]=0;
    ITS->ResetRecPoints();  
    //gAlice->TreeR()->GetEvent(mod+1); //first entry in TreeR is empty
    gAlice->TreeR()->GetEvent(mod); //first entry in TreeR is empty
    numbpoints = recPoints->GetEntries();
    totalpoints+=numbpoints;
    np[mod] = numbpoints;
  //cout<<" mod = "<<mod<<"   numbpoints = "<<numbpoints<<"\n"; getchar();
    vettid[mod] = new Int_t[numbpoints];
    Int_t ii;
    for (ii=0;ii<numbpoints; ii++) *(vettid[mod]+ii)=0;
  }

  AliTPCtrack *track=0;

     
  if(min_t < 0) {min_t = 0; max_t = nt-1;}   

/*
  ///////////////////////////////// Definition of vertex end its error ////////////////////////////
  ////////////////////////// In the future it will be given by a method ///////////////////////////
  Double_t Vx=0.;
  Double_t Vy=0.;
  Double_t Vz=0.;
  
  Float_t sigmavx=0.0050;      // 50  microns
  Float_t sigmavy=0.0050;      // 50  microns
  Float_t sigmavz=0.010;       // 100 microns

  //Vx+=gRandom->Gaus(0,sigmavx);  Vy+=gRandom->Gaus(0,sigmavy);  Vz+=gRandom->Gaus(0,sigmavz);
  TVector vertex(3), ervertex(3)
  vertex(0)=Vx; vertex(1)=Vy; vertex(2)=Vz;
  ervertex(0)=sigmavx;  ervertex(1)=sigmavy;  ervertex(2)=sigmavz;
  /////////////////////////////////////////////////////////////////////////////////////////////////
*/      
 

  TTree tracktree1("TreeT","Tree with ITS tracks");
  AliITSiotrack *iotrack=0;
  tracktree1.Branch("ITStracks","AliITSiotrack",&iotrack,32000,0);

  ofstream out ("AliITSTra.out");
   
  Int_t j;       
  for (j=min_t; j<=max_t; j++) {     
    track=(AliTPCtrack*)tracks.UncheckedAt(j);
    Int_t flaglab=0;
    if (!track) continue;
    ////// elimination of not good tracks ////////////	 
    Int_t ilab=TMath::Abs(track->GetLabel());
    Int_t iii;
    for (iii=0;iii<ngood;iii++) {
	 //cout<<" ilab, gt[iii].lab = "<<ilab<<" "<<gt[iii].lab<<"\n"; getchar();
      if (ilab==gt[iii].lab) { 
	flaglab=1;
	ptg=gt[iii].ptg; 
	pxg=gt[iii].pxg;
	pyg=gt[iii].pyg;
	pzg=gt[iii].pzg;	
	break;
      }
    }
	 //cout<<" j flaglab =  " <<j<<" "<<flaglab<<"\n";  getchar();
    if (!flaglab) continue;  
	 //cout<<" j =  " <<j<<"\n";  getchar();
 
 	 	 
	 ////// new propagation to the end of TPC //////////////
    Double_t xk=77.415;
    track->PropagateTo(xk, 28.94, 1.204e-3);	 //Ne    
	 xk -=0.01;
    track->PropagateTo(xk, 44.77, 1.71);	 //Tedlar
	 xk -=0.04;
    track->PropagateTo(xk, 44.86, 1.45);	 //Kevlar
	 xk -=2.0;
    track->PropagateTo(xk, 41.28, 0.029);	 //Nomex	 
    xk-=16;
    track->PropagateTo(xk,36.2,1.98e-3); //C02
	 xk -=0.01;
    track->PropagateTo(xk, 24.01, 2.7);	 //Al	 
	 xk -=0.01;
    track->PropagateTo(xk, 44.77, 1.71);	 //Tedlar
	 xk -=0.04;
    track->PropagateTo(xk, 44.86, 1.45);	 //Kevlar
	 xk -=0.5;
    track->PropagateTo(xk, 41.28, 0.029);	 //Nomex	 	 	 	 	 	 
	    	     
       /////////////////////////////////////////////////////////////// 	 		 
 
   ///////////////////////////////////////////////////////////////
    AliITStrack trackITS(*track);
    AliITStrack result(*track);
    AliITStrack primarytrack(*track); 
    
///////////////////////////////////////////////////////////////////////////////////////////////
	 TVector Vgeant(3);
	 Vgeant=result.GetVertex(); 
			  
  // Definition of Dv and Zv for vertex constraint	
     Double_t sigmaDv=0.0050;  Double_t sigmaZv=0.010;	
    //Double_t sigmaDv=0.0015;  Double_t sigmaZv=0.0015;				  
	Double_t uniform= gRandom->Uniform();
	Double_t signdv;
	if(uniform<=0.5) signdv=-1.;
	   else
		 signdv=1.;
	 
	Double_t Vr=TMath::Sqrt(Vgeant(0)*Vgeant(0)+ Vgeant(1)*Vgeant(1));
	  Double_t Dv=gRandom->Gaus(signdv*Vr,(Float_t)sigmaDv); 
    Double_t Zv=gRandom->Gaus(Vgeant(2),(Float_t)sigmaZv);
				
  //cout<<" Dv e Zv = "<<Dv<<" "<<Zv<<"\n";				
    trackITS.SetDv(Dv);  trackITS.SetZv(Zv);
    trackITS.SetsigmaDv(sigmaDv); trackITS.SetsigmaZv(sigmaZv); 
    result.SetDv(Dv);  result.SetZv(Zv);
    result.SetsigmaDv(sigmaDv); result.SetsigmaZv(sigmaZv);
    primarytrack.SetDv(Dv);  primarytrack.SetZv(Zv);
    primarytrack.SetsigmaDv(sigmaDv); primarytrack.SetsigmaZv(sigmaZv); 				 				

/////////////////////////////////////////////////////////////////////////////////////////////////  	 
	 	
    primarytrack.PrimaryTrack(rl);
    TVector  d2=primarytrack.Getd2();
    TVector  tgl2=primarytrack.Gettgl2();
    TVector  dtgl=primarytrack.Getdtgl();
    trackITS.Setd2(d2); trackITS.Settgl2(tgl2);  trackITS.Setdtgl(dtgl); 
    result.Setd2(d2); result.Settgl2(tgl2);  result.Setdtgl(dtgl); 	   
	 /*	          	 
    trackITS.SetVertex(vertex); trackITS.SetErrorVertex(ervertex);
    result.SetVertex(vertex);   result.SetErrorVertex(ervertex);   
    */
                         
    Tracking(trackITS,&result,recPoints,vettid, flagvert,rl,geoinfo);  
	     
    // cout<<" progressive track number = "<<j<<"\r";
   // cout<<j<<"\r";
    Int_t NumofCluster=result.GetNumClust();  
   // cout<<" progressive track number = "<<j<<"\n";    // stampa
    Long_t labITS=result.GetLabel();
   // cout << " ITS track label = " << labITS << "\n"; 	// stampa	    
    int lab=track->GetLabel();		    
    //cout << " TPC track label = " << lab <<"\n";      // stampa
	 
	     
//propagation to vertex
	
    Double_t rbeam=3.;
     
    result.Propagation(rbeam);
       	 
	 Double_t C00,C10,C11,C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44;
	 result.GetCElements(C00,C10,C11,C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44);
	 	 
    Double_t pt=TMath::Abs(result.GetPt());
    Double_t Dr=result.GetD();
    Double_t Z=result.GetZ();
    Double_t tgl=result.GetTgl();
    Double_t C=result.GetC();
    Double_t Cy=C/2.;
    Double_t Dz=Z-(tgl/Cy)*TMath::ASin(result.arga(rbeam));
	 Dz-=Vgeant(2);
	  
	 // cout<<" Dr e dz alla fine = "<<Dr<<" "<<Dz<<"\n"; getchar();
    Double_t phi=result.Getphi();
    Double_t phivertex = phi - TMath::ASin(result.argA(rbeam));
    Double_t duepi=2.*TMath::Pi();	 
    if(phivertex>duepi) phivertex-=duepi;
    if(phivertex<0.) phivertex+=duepi;
    Double_t Dtot=TMath::Sqrt(Dr*Dr+Dz*Dz);
	 
//////////////////////////////////////////////////////////////////////////////////////////    	
  
    Int_t idmodule,idpoint;
	 if(NumofCluster >=5)  {            // cinque - sei
	 //if(NumofCluster ==6)  {            // cinque - sei 


      AliITSiotrack outtrack;

      iotrack=&outtrack;

      iotrack->SetStatePhi(phi);
      iotrack->SetStateZ(Z);
      iotrack->SetStateD(Dr);
      iotrack->SetStateTgl(tgl);
      iotrack->SetStateC(C);
		Double_t radius=result.Getrtrack();
		iotrack->SetRadius(radius);
		Int_t charge;
		if(C>0.) charge=-1;  else charge=1;
		iotrack->SetCharge(charge);
		


      iotrack->SetCovMatrix(C00,C10,C11,C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44);  

      Double_t px=pt*TMath::Cos(phivertex);
      Double_t py=pt*TMath::Sin(phivertex);
      Double_t pz=pt*tgl;
		
      Double_t xtrack=Dr*TMath::Sin(phivertex);
      Double_t ytrack=Dr*TMath::Cos(phivertex);
      Double_t ztrack=Dz+Vgeant(2);


      iotrack->SetPx(px);
      iotrack->SetPy(py);
      iotrack->SetPz(pz);
      iotrack->SetX(xtrack);
      iotrack->SetY(ytrack);
      iotrack->SetZ(ztrack);
      iotrack->SetLabel(labITS);
		
      Int_t il;		
		for(il=0;il<6; il++){
		  iotrack->SetIdPoint(il,result.GetIdPoint(il));
		  iotrack->SetIdModule(il,result.GetIdModule(il));		
		}
      tracktree1.Fill();

   //cout<<" labITS = "<<labITS<<"\n";
	//cout<<" phi z Dr tgl C = "<<phi<<" "<<Z<<" "<<Dr<<" "<<tgl<<" "<<C<<"\n";  getchar();	   

     DataOut(kkk) = ptg; kkk++; DataOut(kkk)=labITS; kkk++; DataOut(kkk)=lab; kkk++;		

      for (il=0;il<6;il++) {
        idpoint=result.GetIdPoint(il);
        idmodule=result.GetIdModule(il);
	*(vettid[idmodule]+idpoint)=1; 
	iotrack->SetIdPoint(il,idpoint);
        iotrack->SetIdModule(il,idmodule);
      }
      
      Double_t difpt= (pt-ptg)/ptg*100.;  	            	                
      DataOut(kkk)=difpt; kkk++;                                             
      Double_t lambdag=TMath::ATan(pzg/ptg);
      Double_t   lam=TMath::ATan(tgl);      
      Double_t diflam = (lam - lambdag)*1000.;
      DataOut(kkk) = diflam; kkk++;	    	  			    
      Double_t phig=TMath::ATan2(pyg,pxg);  if(phig<0) phig=2.*TMath::Pi()+phig;       
      Double_t phi=phivertex;
        
      Double_t difphi = (phi - phig)*1000.;
      DataOut(kkk)=difphi; kkk++;
      DataOut(kkk)=Dtot*1.e4; kkk++;
      DataOut(kkk)=Dr*1.e4; kkk++;
      DataOut(kkk)=Dz*1.e4; kkk++; 
      Int_t r;
      for (r=0; r<9; r++) { out<<DataOut(r)<<" ";}
      out<<"\n";
      kkk=0;  
		
	    
    } // end if on NumofCluster
  //gObjectTable->Print();    // stampa memoria     
  }  //  end for (int j=min_t; j<=max_t; j++)
  
  out.close();  
 
 
  static Bool_t first=kTRUE;
  static TFile *tfile;

	if(first) {
	    tfile=new TFile("itstracks.root","RECREATE");
	    //cout<<"I have opened itstracks.root file "<<endl;
	}	    
	first=kFALSE;
	tfile->cd();
	tfile->ls();

   char hname[30];
   sprintf(hname,"TreeT%d",evNumber);

  tracktree1.Write(hname);


  
	    TTree *fAli=gAlice->TreeK();
            TFile *fileAli=0;
	    
	    if (fAli) fileAli =fAli->GetCurrentFile();
	    fileAli->cd();
     
  ////////////////////////////////////////////////////////////////////////////////////////////////

  printf("delete vectors\n");
  if(np) delete [] np;
  if(vettid) delete [] vettid;
  
  if(geoinfo) delete geoinfo;
  
}
