#include "iostream.h"


  void AliITSComparisonV1(Int_t evNumber1=0,Int_t evNumber2=0) {

  const char *filename="itstracks.root";
  
  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);
   
   ifstream in("itsgood_tracks"); 
//
//   Loop over events 
//

   char tname[30];

     
  struct GoodTrack {
   Int_t fEventN;  
    Int_t lab,code;
    Float_t px,py,pz,x,y,z,pxg,pyg,pzg,ptg;
    Bool_t flag;
  };   
    Int_t i;
  GoodTrack gt[15000]; 
  Int_t ngood=0;  
  cerr<<"Reading itsgood tracks...\n";
  while (in>>gt[ngood].fEventN>>gt[ngood].lab>>gt[ngood].code
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
   ofstream out1 ("AliITSTrag.out");
  Int_t countpos=0,countneg=0;

   
  //for (i=0; i<ngood; i++) out1 << gt[i].ptg << "\n";
  for (i=0; i<ngood; i++) {
    out1 << gt[i].ptg << "\n";
    Int_t codpar=gt[i].code;
    if(codpar==2212 || codpar==-11 || codpar==-13 || codpar==211 || codpar==321 || codpar==3222
	 || codpar==213 || codpar==323 || codpar==10323 || codpar==3224 || codpar==2224 || codpar==2214
	 || codpar==-1114 || codpar==-3112 || codpar==-3312 || codpar==3224 || codpar==-3114 || codpar==-3314
	 || codpar==411 || codpar==431 || codpar==413 || codpar==433 || codpar==-15 || codpar==4232
	 || codpar==4222 || codpar==4322 || codpar==4422 || codpar==4412 || codpar==4432 || codpar==4224 
	 ||codpar==4214 || codpar==4324 || codpar==4424 || codpar==4414 || codpar==4434 || codpar==4444)
    countpos++;
	 if(codpar==-2212 || codpar==11 || codpar==13 || codpar==-211 || codpar==-321 || codpar==3112
	 || codpar==-213 || codpar==-323 || codpar==-10323 || codpar==3114 || codpar==1114 || codpar==-2224
	 || codpar==-2214 || codpar==33112 || codpar==-3222 || codpar==3114 || codpar==3314 || codpar==3334 
	 || codpar==-3224 || codpar==-411 || codpar==-431 || codpar==-413 || codpar==-433 || codpar==15 
	 || codpar==-4422 || codpar==-4432 || codpar==-4214 || codpar==-4324 || codpar==-4424 || codpar==-4434 
	 || codpar==-444)
	 countneg++; 		   
  }
  out1.close();
   ofstream out ("AliITSTra.out"); 
///   definition of nm of good track for each event
   TVector ntrackevtpc(evNumber2-evNumber1 +1);
   Int_t nchange=-1;
   Int_t nmev=-100;
   
   
    Int_t i;
    for( i =0; i<ngood ; i++){ 
    if(gt[i].fEventN != nmev ){
    nmev=gt[i].fEventN;
    nchange++; 
    }
    if(nmev == gt[i].fEventN) ntrackevtpc(nchange)++;
    }  
//////////////////////////////////////////////////////////////     
    
   for (Int_t nev=evNumber1; nev<= evNumber2; nev++) {


   sprintf(tname,"TreeT%d",nev);
   TTree *tracktree=(TTree*)file->Get(tname);
   TBranch *tbranch=tracktree->GetBranch("ITStracks");
   cout<<" nev = "<<nev<<"\n";   
	//cout<<" open the file \n"; 
	
   Int_t nentr=tracktree->GetEntries();

   TObjArray tarray(nentr);
  // AliITSIOTrack *iotrack=0;
   printf("nentr %d\n",nentr);
	
   for (Int_t i=0; i<nentr; i++) {
      AliITSIOTrack *iotrack=new AliITSIOTrack;
      // tarray.AddAt(new AliITSIOTrack,i);
      // iotrack=(AliITSiotrack*)tarray.UncheckedAt(i);
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
	tarray.AddLast(iotrack);
   }
	//file->Close();		 
	
	  AliITSIOTrack *iotrack;
   for (Int_t i=0; i<nentr; i++) {
     AliITSIOTrack *iotrack=new AliITSIOTrack;      	
	 iotrack=(AliITSIOTrack*)tarray.UncheckedAt(i);
	 if(!iotrack) continue;
     Int_t labITS=iotrack->GetLabel();
     Int_t labTPC=iotrack->GetTPCLabel();  
     //Int_t labTPC=labITS;  //provvisoria  
	  Double_t phistate=iotrack->GetStatePhi();
	  Double_t tgl=iotrack->GetStateTgl();	
	  Double_t Zstate=iotrack->GetStateZ();
	  Double_t Dr=iotrack->GetStateD();		  
	  Double_t C=iotrack->GetStateC();
	  Double_t px=iotrack->GetPx();
	  Double_t py=iotrack->GetPy();	  
	  Double_t pz=iotrack->GetPz();	
	  Double_t pt=TMath::Sqrt(px*px+py*py);
	  Int_t c = iotrack->GetCharge(); 
	  Double_t x=iotrack->GetX();
	  Double_t y=iotrack->GetY();
	  Double_t z= iotrack->GetZ(); 
	//  Double_t Dz=z;   //non e' vero bisogna levare vertice
	 // Double_t Dtot= TMath::Sqrt(Dr*Dr+Dz*Dz);
	  
     // cout<<" track label = "<<label<<"\n";
     // cout<<" phi z D tanl C = "<<phistate<<" "<<Zstate<<" "<<Dr<<" "<<tgl<<" "<<C<<"\n"; 	  
	  
  //  Int_t ilab=TMath::Abs(iotrack->GetLabel());
    Int_t flaglab=0;   
    Int_t iii=0;
   Double_t ptg=0.,pxg=0.,pyg=0.,pzg=0.;
	Double_t xo=0., yo=0., zo=0.;  	 

    Int_t mingood=0, maxgood=0, jj=0;
    if(nev==evNumber1) mingood=0;
    else
    {for(jj=evNumber1; jj<nev; jj++) mingood += ntrackevtpc(jj);}
     for(jj=evNumber1; jj<=nev; jj++) maxgood+= ntrackevtpc(jj);    

   Int_t ilab=TMath::Abs(labTPC);
   // for (iii=0;iii<ngood;iii++) {
   for(iii=mingood; iii<maxgood; iii++){
	 //cout<<" ilab, gt[iii].lab = "<<ilab<<" "<<gt[iii].lab<<"\n"; getchar();
      if (ilab==gt[iii].lab) { 
	flaglab=1;
	ptg=gt[iii].ptg; 
	pxg=gt[iii].pxg;
	pyg=gt[iii].pyg;
	pzg=gt[iii].pzg;
	xo=gt[iii].x;
	yo=gt[iii].y;
	zo=gt[iii].z;			
	break;
      }
    }   

    if (!flaglab) continue;  
    //cout<<" j =  " <<j<<"\n";  getchar();  	  		  		    
   TVector dataOut(10);
    Int_t kkk=0;
    
    dataOut(kkk) = ptg; kkk++; dataOut(kkk)=labITS; kkk++; dataOut(kkk)=labTPC; kkk++;	
  
      ///////////////////////////////
      Double_t difpt= (pt-ptg)/ptg*100.;
      //cout<<" difpt = "<<difpt<<"\n"; getchar();  	            	                
      dataOut(kkk)=difpt; kkk++;                                             
      Double_t lambdag=TMath::ATan(pzg/ptg);
      Double_t   lam=TMath::ATan(tgl);      
      Double_t diflam = (lam - lambdag)*1000.;
      dataOut(kkk) = diflam; kkk++;	    	  			    
      Double_t phig=TMath::ATan2(pyg,pxg);  if(phig<0) phig=2.*TMath::Pi()+phig;       
     // Double_t phi=phivertex;
       Double_t phi=TMath::ACos(px/pt);
     Double_t duepi=2.*TMath::Pi();	 
     if(phi>duepi) phi-=duepi;
     if(phi<0.) phi+=duepi;      
      Double_t signC=0.; 
      if(c>0) signC=1.; else signC=-1.;
 	  Double_t Dz=z-zo;   // vertex subtraction
	  Double_t Dtot= TMath::Sqrt(Dr*Dr+Dz*Dz);       
      Double_t difphi = (phi - phig)*1000.;
      dataOut(kkk)=difphi; kkk++;
      dataOut(kkk)=Dtot*1.e4; kkk++;
      dataOut(kkk)=Dr*1.e4; kkk++;
      dataOut(kkk)=Dz*1.e4; kkk++;
      dataOut(kkk)=signC; kkk++; 
      Int_t r;
      for (r=0; r<10; r++) { out<<dataOut(r)<<" ";}
      out<<"\n";
      kkk=0;        	    

    delete iotrack;		 
   }  
  //out.close(); 
   }   // event loop 
  out.close();
  in.close();   
  file->Close();      
}

