#include <iostream.h>
#include <TMath.h> 
#include <TList.h> 
#include <TTree.h> 
#include <TVector.h>
#include <TMatrix.h>
#include <TObjectTable.h>

#include "AliITStracking.h"
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITStrack.h"

ClassImp(AliITStracking)
 

AliITStracking::AliITStracking(TList *trackITSlist, AliITStrack *reference, 
                AliITS *aliITS, TObjArray *rpoints, Double_t Ptref, Int_t **vettid, Bool_t flagvert,  AliITSRad *rl) {										 
///////////////////////   This function perform the tracking in ITS detectors /////////////////////
///////////////////////     reference is a pointer to the final best track    ///////////////////// 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// The authors  thank   Mariana Bondila to have help them to resolve some problems.  July-2000                                                      

  Rlayer[0]=4.; Rlayer[1]=7.;  Rlayer[2]=14.9;  Rlayer[3]=23.8;  Rlayer[4]=39.1;  Rlayer[5]=43.6;
  
  Int_t index;   
  for(index =0; index<trackITSlist->GetSize(); index++) {
    AliITStrack *trackITS = (AliITStrack *) trackITSlist->At(index);

    if((*trackITS).GetLayer()==7) reference->SetChi2(10.223e140);
   // cout <<" Layer inizio = "<<(*trackITS).GetLayer()<<"\n";
   //  cout<<"fvtrack =" <<"\n";
   //  cout << (*trackITS)(0) << " "<<(*trackITS)(1)<<" "<<(*trackITS)(2)<<" "<<(*trackITS)(3)<<" "<<(*trackITS)(4)<<"\n";
   //  cout<< " rtrack = "<<(*trackITS).Getrtrack()<<"\n";
   //  cout<< " Pt = "<<(*trackITS).GetPt()<<"\n";
   //  getchar();    
    Double_t Chi2Now, Chi2Ref;
    if((*trackITS).GetLayer()==1 ) {
      Chi2Now = trackITS->GetChi2();
      Float_t NumClustNow = trackITS->GetNumClust();
      if(trackITS->GetNumClust()) Chi2Now /= (Double_t )trackITS->GetNumClust();
      Chi2Ref = reference->GetChi2();
      Float_t NumClustRef = reference->GetNumClust();	
      if(reference->GetNumClust()) Chi2Ref /= (Double_t )reference->GetNumClust();
      //cout<<" Chi2Now and Chi2Ref = "<<Chi2Now<<" "<<Chi2Ref<<"\n";
		if( NumClustNow > NumClustRef ) {*reference = *trackITS;} 
      if((NumClustNow == NumClustRef )&& (Chi2Now < Chi2Ref))  {*reference = *trackITS;}
      continue;	
    }
    Float_t NumClustNow = trackITS->GetNumClust();
    if(NumClustNow) { 
      Chi2Now = trackITS->GetChi2();
      Chi2Now/=NumClustNow;
      //cout<<" Chi2Now =  "<<Chi2Now<<"\n"; 
    /*
    // if(Ptref > 0.6 && Chi2Now > 20.) continue; 
    if(Ptref > 0.6 && Chi2Now > 30.) continue; 	  
    if((Ptref <= 0.6 && Ptref>0.2)&& Chi2Now > 15.) continue;        
     // if(Chi2Now>5.) continue; 
      //if(Chi2Now>15.) continue;     
     // if(Ptref <= 0.2 && Chi2Now > 10.) continue;  
     if(Ptref <= 0.2 && Chi2Now > 8.) continue; 
     */
     if(Ptref > 1.0 && Chi2Now > 30.) continue; 
     if((Ptref >= 0.6 && Ptref<=1.0) && Chi2Now > 40.) continue;     	  
     if((Ptref <= 0.6 && Ptref>0.2)&& Chi2Now > 40.) continue;            
     if(Ptref <= 0.2 && Chi2Now > 8.) continue;        		 		    	       	 	 	 
    }
     	         
    Int_t layerInit = (*trackITS).GetLayer();
    Int_t layernew = layerInit - 2;  // -1 for new layer, -1 for matrix index 
					  
    Int_t NLadder[]= {20, 40, 14, 22, 34, 38}; 
    Int_t NDetector[]= {4,  4,   6,  8, 23, 26}; 
				 		
    TList listoftrack;    	 
    Int_t ladp, ladm, detp,detm,ladinters,detinters; 	
    Int_t layerfin=layerInit-1;
    Double_t Rfin=Rlayer[layerfin-1];
    // cout<<"Prima di intersection \n";

    Int_t  outinters=NewIntersection(*trackITS, Rfin, layerfin, ladinters, detinters);
	 	 
   // cout<<" outinters = "<<outinters<<"\n";
   //  cout<<" Layer ladder detector intersection ="<<layerfin<<" "<<ladinters<<" "<<detinters<<"\n";
   //  cout << " phiinters zinters = "<<(*trackITS)(0) << " "<<(*trackITS)(1)<<"\n"; getchar();
			 
    if(outinters==-1) continue;
	 
    Int_t flaghit=0;   	 	        
    if(outinters==0){   
      TVector Touclad(9), Toucdet(9);	 
      Int_t lycur=layerfin;                                            
      ladp=ladinters+1;
      ladm=ladinters-1;
      if(ladm <= 0) ladm=NLadder[layerfin-1];  
      if(ladp > NLadder[layerfin-1]) ladp=1;		
      detp=detinters+1;
      detm=detinters-1;
      Int_t idetot=1;
      Touclad(0)=ladinters; Touclad(1)=ladm; Touclad(2)=ladp;
      Touclad(3)=ladinters; Touclad(4)=ladm; Touclad(5)=ladp;
      Touclad(6)=ladinters; Touclad(7)=ladm; Touclad(8)=ladp;
      Toucdet(0)=detinters; Toucdet(1)=detinters; Toucdet(2)=detinters;
      if(detm > 0 && detp <= NDetector[layerfin-1]) {
        idetot=9;
        Toucdet(3)=detm; Toucdet(4)=detm; Toucdet(5)=detm;	   
        Toucdet(6)=detp; Toucdet(7)=detp; Toucdet(8)=detp;
      }
	 
      if(detm > 0 && detp > NDetector[layerfin-1]) {
        idetot=6;
        Toucdet(3)=detm; Toucdet(4)=detm; Toucdet(5)=detm;
      }
	 
      if(detm <= 0 && detp <= NDetector[layerfin-1]) {
        idetot=6;
        Toucdet(3)=detp; Toucdet(4)=detp; Toucdet(5)=detp;
      }
      Int_t iriv; 	
      for (iriv=0; iriv<idetot; iriv++) {  //for on detectors
        AliITSgeom *g1 = aliITS->GetITSgeom();  
        TVector CTF(9);
        g1->GetCenterThetaPhi(layerInit-1,(Int_t)Touclad(iriv),(Int_t)Toucdet(iriv),CTF);

        // cout<<" layer, ladder, det, xo, yo, zo = "<<layerInit-1<<" "<<(Int_t)Touclad(iriv)<<
        // " "<<(Int_t)Toucdet(iriv)<<" "<<CTF(0)<<" "<<CTF(1)<<" "<<CTF(2)<< " "<<CTF(6)<<"\n"; getchar(); 

        ////////////////////////////////////////////////////////////////////////////////////////////////

        /*** Rec points sorted by module *****/
        /**************************************/

        Int_t index;
        AliITSRecPoint *recp;
        AliITSgeom *geom = aliITS->GetITSgeom();
        index = geom->GetModuleIndex(lycur,Touclad(iriv),Toucdet(iriv));
        Int_t lay,lad,det;
        geom->GetModuleId(index,lay,lad,det);
        aliITS->ResetRecPoints();
        //gAlice->TreeR()->GetEvent(index+1); //first entry in TreeR is empty
        gAlice->TreeR()->GetEvent(index); //first entry in TreeR is empty

        Int_t npoints=rpoints->GetEntries();
        Int_t *indlist=new Int_t[npoints+1];
        Int_t counter=0;
        Int_t ind;
        for (ind=0; ind<=npoints; ind++) {
          indlist[ind]=-1;
	       if (*(vettid[index]+ind)==0) {
              indlist[counter]=ind;
	           counter++;
	      }
        }

        ind=-1;
	
        for(;;) { 
          ind++;
          if(indlist[ind] < 0) recp=0;
	       else recp = (AliITSRecPoint*)rpoints->UncheckedAt(indlist[ind]);

	  if((!recp)  )  break;	
	  TVector cluster(3),vecclust(9);
	  vecclust(6)=vecclust(7)=vecclust(8)=-1.;
	  Double_t sigma[2];
	    //modificata 8-3-2001                
	  // set veclust in global
	  Float_t global[3], local[3];
	  local[0]=recp->GetX();
	  local[1]=0.;
	  local[2]= recp->GetZ();
          AliITSgeom *g1 = aliITS->GetITSgeom();
          Int_t play = lycur;
          Int_t plad = TMath::Nint(Touclad(iriv));   
          Int_t pdet = TMath::Nint(Toucdet(iriv)); 		
          g1->LtoG(play,plad,pdet,local,global); 
	  
          vecclust(0)=global[0];
          vecclust(1)=global[1];
          vecclust(2)=global[2];
	  
	  /*
	  ////  modificato 8-3-2001
          vecclust(0)=recp->GetRhit();;
          vecclust(1)=recp->Getphi();
          vecclust(2)=recp->GetZglobal();
	  */
	  	  	  	  	  	     
          vecclust(3) = (float)recp->fTracks[0]; 
          vecclust(4) = (float)indlist[ind];
          vecclust(5) = (float)index;
          vecclust(6) = (float)recp->fTracks[0];
          vecclust(7) = (float)recp->fTracks[1];
          vecclust(8) = (float)recp->fTracks[2];
     
          sigma[0] = (Double_t)  recp->GetSigmaX2();	   
	  sigma[1] = (Double_t) recp->GetSigmaZ2();
	   //commentato 8-3-2001		 
         //now we are in r,phi,z in global
          cluster(0) = TMath::Sqrt(vecclust(0)*vecclust(0)+vecclust(1)*vecclust(1));//r hit
          cluster(1) = PhiDef(vecclust(0),vecclust(1));    // phi hit
          cluster(2) = vecclust(2);                   // z hit
	  
	  /*
	 //modificato 8-3-2001
          cluster(0) = vecclust(0);//r hit
          cluster(1) = vecclust(1);    // phi hit
          cluster(2) = vecclust(2);                   // z hit	
	  */	  	
	 // cout<<" layer = "<<play<<"\n";
	 // cout<<" cluster prima = "<<vecclust(0)<<" "<<vecclust(1)<<" "
	 // <<vecclust(2)<<"\n"; getchar();    
          //cluster(1)= cluster(1)-trackITS->Getalphaprov();  //provvisorio;
			 //if(cluster(1)<0.) cluster(1)+=2.*TMath::Pi(); //provvisorio
			 //cout<<" cluster(1) dopo = "<<cluster(1)<< " alphaprov = "<<trackITS->Getalphaprov()<<"\n"; 			 
          Float_t sigmatotphi, sigmatotz;
   		  		  
          //Float_t epsphi=3.2, epsz=3.; 
	    Float_t epsphi=5.0, epsz=5.0;              
          if(Ptref<0.2) {epsphi=3.; epsz=3.;}
		  		  
          Double_t Rtrack=(*trackITS).Getrtrack();
          Double_t sigmaphi=sigma[0]/(Rtrack*Rtrack);
          sigmatotphi=epsphi*TMath::Sqrt(sigmaphi + (*trackITS).GetSigmaphi());
	  	 
          sigmatotz=epsz*TMath::Sqrt(sigma[1] + (*trackITS).GetSigmaZ());
  //cout<<"cluster e sigmatotphi e track = "<<cluster(0)<<" "<<cluster(1)<<" "<<sigmatotphi<<" "<<vecclust(3)<<"\n";
  //if(vecclust(3)==481) getchar();
	       if(cluster(1)<6. && (*trackITS).Getphi()>6.) cluster(1)=cluster(1)+(2.*TMath::Pi());
	       if(cluster(1)>6. && (*trackITS).Getphi()<6.) cluster(1)=cluster(1)-(2.*TMath::Pi());			 	  
          if(TMath::Abs(cluster(1)-(*trackITS).Getphi()) > sigmatotphi) continue;
			// cout<<" supero sigmaphi \n";      
          AliITStrack *newTrack = new AliITStrack((*trackITS));
          (*newTrack).SetLayer((*trackITS).GetLayer()-1); 
 
	       if (TMath::Abs(Rtrack-cluster(0))/Rtrack>1e-6) 
                                               (*newTrack).Correct(Double_t(cluster(0)));	
			//cout<<" cluster(2) e (*newTrack).GetZ() = "<<cluster(2)<<" "<<	(*newTrack).GetZ()<<"\n";									  	
          if(TMath::Abs(cluster(2)-(*newTrack).GetZ()) > sigmatotz){ 
             delete newTrack;
             continue;}
       
          if(iriv == 0) flaghit=1;
 
          (*newTrack).AddMS(rl);  // add the multiple scattering matrix to the covariance matrix 
	  (*newTrack).AddEL(rl,1.,0);
	   	 
	  Double_t sigmanew[2];
	  sigmanew[0]= sigmaphi;
	  sigmanew[1]=sigma[1];

	  if(flagvert) 	 
	    KalmanFilterVert(newTrack,cluster,sigmanew);  
	  else		    	  		 		  
	    KalmanFilter(newTrack,cluster,sigmanew);
	       
      		  
          (*newTrack).PutCluster(layernew, vecclust);
           newTrack->AddClustInTrack();            
 		 		          		  		
	   listoftrack.AddLast(newTrack);

        }   // end of for(;;) on rec points 

        delete [] indlist;
  
      }  // end of for on detectors
     
    }//end if(outinters==0) 
  
    if(flaghit==0 || outinters==-2) {
      AliITStrack *newTrack = new AliITStrack(*trackITS);	 
      (*newTrack).SetLayer((*trackITS).GetLayer()-1); 
      (*newTrack).AddMS(rl);  // add the multiple scattering matrix to the covariance matrix  
      (*newTrack).AddEL(rl,1.,0);  	      
				    		  
      listoftrack.AddLast(newTrack);	  
    }	
 	        

    //gObjectTable->Print();   // stampa memoria
	 
    AliITStracking(&listoftrack, reference, aliITS, rpoints,Ptref,vettid,flagvert,rl);          
    listoftrack.Delete();
  } // end of for on tracks

  //gObjectTable->Print();   // stampa memoria

}   


Int_t AliITStracking::NewIntersection(AliITStrack &track, Double_t rk,Int_t layer, Int_t &ladder, Int_t &detector) { 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// Found the intersection and the detector 

  if(track.DoNotCross(rk)){ /*cout<< " Do not cross \n";*/ return -1;} 
  track.Propagation(rk);
  Double_t zinters=track.GetZ();
  Double_t phinters=track.Getphi();
  //phinters = phinters+track.Getalphaprov(); //provvisorio
  //if(phinters>2.*3.14) phinters=phinters-2.*3.14; //provvisorio
  //cout<<"zinters = "<<zinters<<"  phinters = "<<phinters<<"\n";

  //////////////////////////////////      limits for Geometry 5      /////////////////////////////
  
  Int_t NLadder[]= {20, 40, 14, 22, 34, 38};
  Int_t NDetector[]= {4,  4,   6,  8, 23, 26}; 

  Float_t Detx[]= {0.64, 0.64, 3.509, 3.509, 3.65, 3.65 };
  Float_t Detz[]= {4.19, 4.19, 3.75 , 3.75 , 2   , 2    };
  ////////////////////////////////////////////////////////////////////////////////////////////////  
  
  TVector det(9);
  TVector ListDet(2);
  TVector DistZCenter(2);  
  AliITSgeom *g1 = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
  Int_t iz=0; 
  Double_t epsz=1.2;
  Double_t epszpixel=0.05;

  Int_t iD;
  for(iD = 1; iD<= NDetector[layer-1]; iD++) {
    g1->GetCenterThetaPhi(layer,1,iD,det);
    Double_t zmin=det(2)-Detz[layer-1];
    if(iD==1) zmin=det(2)-(Detz[layer-1])*epsz;		
    Double_t zmax=det(2)+Detz[layer-1];
    if(iD==NDetector[layer-1]) zmax=det(2)+(Detz[layer-1])*epsz;
    //added to take into account problem on drift
    if(layer == 4 || layer==3) zmin=zmin-epszpixel; zmax=zmax+epszpixel;
    //cout<<"zmin zinters zmax det(2)= "<<zmin<<" "<<zinters<<" "<<zmax<<" "<<det(2)<<"\n";	
    if(zinters > zmin && zinters <= zmax) { 
      if(iz>1) {cout<< " Errore su iz in NewIntersection \n"; getchar();}
      else {
        ListDet(iz)= iD; DistZCenter(iz)=TMath::Abs(zinters-det(2)); iz++;
      }
    }			      
  }
  
  if(iz==0) {/* cout<< " No detector along Z \n";*/ return -2;}
  detector=Int_t (ListDet(0));
  if(iz>1 && (DistZCenter(0)>DistZCenter(1)))   detector=Int_t (ListDet(1));
  
  AliITSgeom *g2 = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom(); 
  Float_t global[3];
  Float_t local[3];
  TVector ListLad(2);
  TVector DistphiCenter(2);
  Int_t ip=0;
  Double_t pigre=TMath::Pi();
  
  Int_t iLd;   
  for(iLd = 1; iLd<= NLadder[layer-1]; iLd++) {
          g1->GetCenterThetaPhi(layer,iLd,detector,det);
  Double_t phidet=PhiDef(Double_t(det(0)),Double_t(det(1)));
  // cout<<" layer phidet e det(6) = "<< layer<<" "<<phidet<<" "<<det(6)<<"\n"; getchar();
  Double_t xmin,ymin,xmax,ymax;	
 // Double_t phiconfr=0.0;
  //cout<<" phiconfr inizio =  "<<phiconfr <<"\n"; getchar();  
  local[1]=local[2]=0.;  
  local[0]= -(Detx[layer-1]);
  if(layer==1)    local[0]= (Detx[layer-1]);  //take into account different reference system
  g2->LtoG(layer,iLd,detector,local,global);
  xmax=global[0]; ymax=global[1];
  local[0]= (Detx[layer-1]);
  if(layer==1)    local[0]= -(Detx[layer-1]);  //take into account different reference system  
  g2->LtoG(layer,iLd,detector,local,global);
  xmin=global[0]; ymin=global[1];
  Double_t phimin=PhiDef(xmin,ymin);
  Double_t phimax=PhiDef(xmax,ymax);
  //cout<<" xmin ymin = "<<xmin<<" "<<ymin<<"\n";
  // cout<<" xmax ymax = "<<xmax<<" "<<ymax<<"\n";  
  // cout<<" iLd phimin phimax ="<<iLd<<" "<<phimin<<" "<<phimax<<"\n";

  Double_t phiconfr=phinters;
 if(phimin>phimax ){    
     if(phimin <5.5) {cout<<" Error in NewIntersection for phi \n"; getchar();}
    phimin=phimin-(2.*pigre);
    if(phinters>(1.5*pigre)) phiconfr=phinters-(2.*pigre); 
    if(phidet>(1.5*pigre)) phidet=phidet-(2.*pigre);
  }              
  //  cout<<" phiconfr finale = "<<phiconfr<<"\n"; getchar(); 
  if(phiconfr>phimin && phiconfr<= phimax) {
    if(ip>1) {
      cout<< " Errore su ip in NewIntersection \n"; getchar();
    }
      else  {
        ListLad(ip)= iLd; DistphiCenter(ip)=TMath::Abs(phiconfr-phidet); ip++;
      }  
    }
  }
  if(ip==0) { cout<< " No detector along phi \n"; getchar();}
  ladder=Int_t (ListLad(0));
  if(ip>1 && (DistphiCenter(0)>DistphiCenter(1)))   ladder=Int_t (ListLad(1));       

  return 0;
}


Double_t AliITStracking::PhiDef(Double_t x, Double_t y){
  Double_t pigre= TMath::Pi();
  Double_t phi=0.0;
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


void AliITStracking::KalmanFilter(AliITStrack *newTrack,TVector &cluster,Double_t sigma[2]){ 
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// Kalman filter without vertex constraint


  ////////////////////////////// Evaluation of the measurement vector /////////////////////////////////////  

  Double_t m[2];
  Double_t rk,phik,zk;
  rk=cluster(0);   phik=cluster(1);  zk=cluster(2);
  m[0]=phik;    m[1]=zk; 
       
  ///////////////////////////////////// Evaluation of the error matrix V  ///////////////////////////////          

  Double_t V00=sigma[0];
  Double_t V11=sigma[1];
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  
  Double_t Cin00,Cin10,Cin20,Cin30,Cin40,Cin11,Cin21,Cin31,Cin41,Cin22,Cin32,Cin42,Cin33,Cin43,Cin44;
			    
  newTrack->GetCElements(Cin00,Cin10,Cin11,Cin20,Cin21,Cin22,Cin30,Cin31,Cin32,Cin33,Cin40,
                         Cin41,Cin42,Cin43,Cin44); //get C matrix
  			  
  Double_t Rold00=Cin00+V00;
  Double_t Rold10=Cin10;
  Double_t Rold11=Cin11+V11;
  
//////////////////////////////////// R matrix inversion  ///////////////////////////////////////////////
  
  Double_t det=Rold00*Rold11-Rold10*Rold10;
  Double_t R00=Rold11/det;
  Double_t R10=-Rold10/det;
  Double_t R11=Rold00/det;

////////////////////////////////////////////////////////////////////////////////////////////////////////			    

  Double_t K00=Cin00*R00+Cin10*R10;
  Double_t K01=Cin00*R10+Cin10*R11;
  Double_t K10=Cin10*R00+Cin11*R10;  
  Double_t K11=Cin10*R10+Cin11*R11;
  Double_t K20=Cin20*R00+Cin21*R10;  
  Double_t K21=Cin20*R10+Cin21*R11;  
  Double_t K30=Cin30*R00+Cin31*R10;  
  Double_t K31=Cin30*R10+Cin31*R11;  
  Double_t K40=Cin40*R00+Cin41*R10;
  Double_t K41=Cin40*R10+Cin41*R11;
  
  Double_t X0,X1,X2,X3,X4;
  newTrack->GetXElements(X0,X1,X2,X3,X4);     // get the state vector
  
  Double_t savex0=X0, savex1=X1;
  
  X0+=K00*(m[0]-savex0)+K01*(m[1]-savex1);
  X1+=K10*(m[0]-savex0)+K11*(m[1]-savex1);
  X2+=K20*(m[0]-savex0)+K21*(m[1]-savex1);
  X3+=K30*(m[0]-savex0)+K31*(m[1]-savex1);
  X4+=K40*(m[0]-savex0)+K41*(m[1]-savex1);
  
  Double_t C00,C10,C20,C30,C40,C11,C21,C31,C41,C22,C32,C42,C33,C43,C44;
  
  C00=Cin00-K00*Cin00-K01*Cin10;
  C10=Cin10-K00*Cin10-K01*Cin11;
  C20=Cin20-K00*Cin20-K01*Cin21;
  C30=Cin30-K00*Cin30-K01*Cin31;
  C40=Cin40-K00*Cin40-K01*Cin41;
  
  C11=Cin11-K10*Cin10-K11*Cin11;
  C21=Cin21-K10*Cin20-K11*Cin21;
  C31=Cin31-K10*Cin30-K11*Cin31;
  C41=Cin41-K10*Cin40-K11*Cin41;
  
  C22=Cin22-K20*Cin20-K21*Cin21;
  C32=Cin32-K20*Cin30-K21*Cin31;
  C42=Cin42-K20*Cin40-K21*Cin41;

  C33=Cin33-K30*Cin30-K31*Cin31;
  C43=Cin43-K30*Cin40-K31*Cin41;
  
  C44=Cin44-K40*Cin40-K41*Cin41;
  
  newTrack->PutXElements(X0,X1,X2,X3,X4);               // put the new state vector
   
  newTrack->PutCElements(C00,C10,C11,C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44); // put in track the
                                                                                       // new cov matrix  
  Double_t VMCold00=V00-C00;
  Double_t VMCold10=-C10;
  Double_t VMCold11=V11-C11;
  
///////////////////////////////////// Matrix VMC inversion  ////////////////////////////////////////////////
  
  det=VMCold00*VMCold11-VMCold10*VMCold10;
  Double_t VMC00=VMCold11/det;
  Double_t VMC10=-VMCold10/det;
  Double_t VMC11=VMCold00/det;

////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Double_t chi2=(m[0]-X0)*( VMC00*(m[0]-X0) + 2.*VMC10*(m[1]-X1) ) +                    
                (m[1]-X1)*VMC11*(m[1]-X1);	 
   	 
  newTrack->SetChi2(newTrack->GetChi2()+chi2);
   
} 


void AliITStracking::KalmanFilterVert(AliITStrack *newTrack,TVector &cluster,Double_t sigma[2]){
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 
// Kalman filter with vertex constraint

  ////////////////////////////// Evaluation of the measurement vector m ///////////////  

  Double_t m[4];
  Double_t rk,phik,zk;
  rk=cluster(0);   phik=cluster(1);  zk=cluster(2);
  m[0]=phik;    m[1]=zk;
 
  Double_t CC=(*newTrack).GetC();
  Double_t Zv=(*newTrack).GetZv(); 
  Double_t Dv=(*newTrack).GetDv();
  Double_t Cy=CC/2.;
  Double_t tgl= (zk-Zv)*Cy/TMath::ASin(Cy*rk);
  m[2]=Dv;    m[3]=tgl;

  ///////////////////////////////////// Evaluation of the error matrix V  //////////////
  Int_t Layer=newTrack->GetLayer();
  Double_t V00=sigma[0];
  Double_t V11=sigma[1];
  Double_t V31=sigma[1]/rk;
  Double_t sigmaDv=newTrack->GetsigmaDv();
  Double_t V22=sigmaDv*sigmaDv  + newTrack->Getd2(Layer-1);
  Double_t V32=newTrack->Getdtgl(Layer-1);
  Double_t sigmaZv=newTrack->GetsigmaZv();  
  Double_t V33=(sigma[1]+sigmaZv*sigmaZv)/(rk*rk) + newTrack->Gettgl2(Layer-1);
  ///////////////////////////////////////////////////////////////////////////////////////
  
  Double_t Cin00,Cin10,Cin11,Cin20,Cin21,Cin22,Cin30,Cin31,Cin32,Cin33,Cin40,Cin41,Cin42,Cin43,Cin44;
			    
  newTrack->GetCElements(Cin00,Cin10,Cin11,Cin20,Cin21,Cin22,Cin30,Cin31,Cin32,Cin33,Cin40,
                         Cin41,Cin42,Cin43,Cin44); //get C matrix
			  
  Double_t R[4][4];
  R[0][0]=Cin00+V00;
  R[1][0]=Cin10;
  R[2][0]=Cin20;
  R[3][0]=Cin30;
  R[1][1]=Cin11+V11;
  R[2][1]=Cin21;
  R[3][1]=Cin31+sigma[1]/rk;
  R[2][2]=Cin22+sigmaDv*sigmaDv+newTrack->Getd2(Layer-1);
  R[3][2]=Cin32+newTrack->Getdtgl(Layer-1);
  R[3][3]=Cin33+(sigma[1]+sigmaZv*sigmaZv)/(rk*rk) + newTrack->Gettgl2(Layer-1);
  
  R[0][1]=R[1][0]; R[0][2]=R[2][0]; R[0][3]=R[3][0]; R[1][2]=R[2][1]; R[1][3]=R[3][1]; 
  R[2][3]=R[3][2];

/////////////////////  Matrix R inversion ////////////////////////////////////////////
 
  const Int_t n=4;
  Double_t big, hold;
  Double_t d=1.;
  Int_t ll[n],mm[n];

  Int_t i,j,k;
  
  for(k=0; k<n; k++) {
    ll[k]=k;
    mm[k]=k;
    big=R[k][k];
    for(j=k; j<n ; j++) {
      for (i=j; i<n; i++) {
        if(TMath::Abs(big) < TMath::Abs(R[i][j]) ) { big=R[i][j]; ll[k]=i; mm[k]=j; }
      }
    }    
//
    j= ll[k];
    if(j > k) {
      for(i=0; i<n; i++) { hold=-R[k][i]; R[k][i]=R[j][i]; R[j][i]=hold; }
      
    }
//
    i=mm[k];
    if(i > k ) { 
      for(j=0; j<n; j++) { hold=-R[j][k]; R[j][k]=R[j][i]; R[j][i]=hold; }
    }
//
    if(!big) {
      d=0.;
      cout << "Singular matrix\n"; 
    }
    for(i=0; i<n; i++) {
      if(i == k) { continue; }    
      R[i][k]=R[i][k]/(-big);
    }   
//
    for(i=0; i<n; i++) {
      hold=R[i][k];
      for(j=0; j<n; j++) {
        if(i == k || j == k) { continue; }
	R[i][j]=hold*R[k][j]+R[i][j];
      }
    }
//  
    for(j=0; j<n; j++) {
      if(j == k) { continue; }
      R[k][j]=R[k][j]/big;
    }
//
    d=d*big;
//
    R[k][k]=1./big;        
  } 
//  
  for(k=n-1; k>=0; k--) {
    i=ll[k];
    if(i > k) {
      for (j=0; j<n; j++) {hold=R[j][k]; R[j][k]=-R[j][i]; R[j][i]=hold;}
    }  
    j=mm[k];
    if(j > k) {
      for (i=0; i<n; i++) {hold=R[k][i]; R[k][i]=-R[j][i]; R[j][i]=hold;}
      }
  }
//////////////////////////////////////////////////////////////////////////////////


  Double_t K00=Cin00*R[0][0]+Cin10*R[1][0]+Cin20*R[2][0]+Cin30*R[3][0];
  Double_t K01=Cin00*R[1][0]+Cin10*R[1][1]+Cin20*R[2][1]+Cin30*R[3][1];
  Double_t K02=Cin00*R[2][0]+Cin10*R[2][1]+Cin20*R[2][2]+Cin30*R[3][2];
  Double_t K03=Cin00*R[3][0]+Cin10*R[3][1]+Cin20*R[3][2]+Cin30*R[3][3];
  Double_t K10=Cin10*R[0][0]+Cin11*R[1][0]+Cin21*R[2][0]+Cin31*R[3][0];  
  Double_t K11=Cin10*R[1][0]+Cin11*R[1][1]+Cin21*R[2][1]+Cin31*R[3][1];
  Double_t K12=Cin10*R[2][0]+Cin11*R[2][1]+Cin21*R[2][2]+Cin31*R[3][2];
  Double_t K13=Cin10*R[3][0]+Cin11*R[3][1]+Cin21*R[3][2]+Cin31*R[3][3];
  Double_t K20=Cin20*R[0][0]+Cin21*R[1][0]+Cin22*R[2][0]+Cin32*R[3][0];  
  Double_t K21=Cin20*R[1][0]+Cin21*R[1][1]+Cin22*R[2][1]+Cin32*R[3][1];  
  Double_t K22=Cin20*R[2][0]+Cin21*R[2][1]+Cin22*R[2][2]+Cin32*R[3][2];
  Double_t K23=Cin20*R[3][0]+Cin21*R[3][1]+Cin22*R[3][2]+Cin32*R[3][3];
  Double_t K30=Cin30*R[0][0]+Cin31*R[1][0]+Cin32*R[2][0]+Cin33*R[3][0];  
  Double_t K31=Cin30*R[1][0]+Cin31*R[1][1]+Cin32*R[2][1]+Cin33*R[3][1];  
  Double_t K32=Cin30*R[2][0]+Cin31*R[2][1]+Cin32*R[2][2]+Cin33*R[3][2];  
  Double_t K33=Cin30*R[3][0]+Cin31*R[3][1]+Cin32*R[3][2]+Cin33*R[3][3];
  Double_t K40=Cin40*R[0][0]+Cin41*R[1][0]+Cin42*R[2][0]+Cin43*R[3][0];
  Double_t K41=Cin40*R[1][0]+Cin41*R[1][1]+Cin42*R[2][1]+Cin43*R[3][1];
  Double_t K42=Cin40*R[2][0]+Cin41*R[2][1]+Cin42*R[2][2]+Cin43*R[3][2];  
  Double_t K43=Cin40*R[3][0]+Cin41*R[3][1]+Cin42*R[3][2]+Cin43*R[3][3];
  
  Double_t X0,X1,X2,X3,X4;
  newTrack->GetXElements(X0,X1,X2,X3,X4);     // get the state vector
  
  Double_t savex0=X0, savex1=X1, savex2=X2, savex3=X3;
  
  X0+=K00*(m[0]-savex0)+K01*(m[1]-savex1)+K02*(m[2]-savex2)+
      K03*(m[3]-savex3);
  X1+=K10*(m[0]-savex0)+K11*(m[1]-savex1)+K12*(m[2]-savex2)+
      K13*(m[3]-savex3);
  X2+=K20*(m[0]-savex0)+K21*(m[1]-savex1)+K22*(m[2]-savex2)+
      K23*(m[3]-savex3);
  X3+=K30*(m[0]-savex0)+K31*(m[1]-savex1)+K32*(m[2]-savex2)+
      K33*(m[3]-savex3);
  X4+=K40*(m[0]-savex0)+K41*(m[1]-savex1)+K42*(m[2]-savex2)+
      K43*(m[3]-savex3);       

  Double_t C00,C10,C20,C30,C40,C11,C21,C31,C41,C22,C32,C42,C33,C43,C44;
  
  C00=Cin00-K00*Cin00-K01*Cin10-K02*Cin20-K03*Cin30;
  C10=Cin10-K00*Cin10-K01*Cin11-K02*Cin21-K03*Cin31;
  C20=Cin20-K00*Cin20-K01*Cin21-K02*Cin22-K03*Cin32;
  C30=Cin30-K00*Cin30-K01*Cin31-K02*Cin32-K03*Cin33;
  C40=Cin40-K00*Cin40-K01*Cin41-K02*Cin42-K03*Cin43;
  
  C11=Cin11-K10*Cin10-K11*Cin11-K12*Cin21-K13*Cin31;
  C21=Cin21-K10*Cin20-K11*Cin21-K12*Cin22-K13*Cin32;
  C31=Cin31-K10*Cin30-K11*Cin31-K12*Cin32-K13*Cin33;
  C41=Cin41-K10*Cin40-K11*Cin41-K12*Cin42-K13*Cin43;
  
  C22=Cin22-K20*Cin20-K21*Cin21-K22*Cin22-K23*Cin32;
  C32=Cin32-K20*Cin30-K21*Cin31-K22*Cin32-K23*Cin33;
  C42=Cin42-K20*Cin40-K21*Cin41-K22*Cin42-K23*Cin43;

  C33=Cin33-K30*Cin30-K31*Cin31-K32*Cin32-K33*Cin33;
  C43=Cin43-K30*Cin40-K31*Cin41-K32*Cin42-K33*Cin43;
  
  C44=Cin44-K40*Cin40-K41*Cin41-K42*Cin42-K43*Cin43;
  
  newTrack->PutXElements(X0,X1,X2,X3,X4);               // put the new state vector
  
  newTrack->PutCElements(C00,C10,C11,C20,C21,C22,C30,C31,C32,C33,C40,C41,C42,C43,C44); // put in track the
                                                                                       // new cov matrix
  
  Double_t VMC[4][4];
  
  VMC[0][0]=V00-C00; VMC[1][0]=-C10; VMC[2][0]=-C20; VMC[3][0]=-C30;
  VMC[1][1]=V11-C11; VMC[2][1]=-C21; VMC[3][1]=V31-C31;
  VMC[2][2]=V22-C22; VMC[3][2]=V32-C32;
  VMC[3][3]=V33-C33;
  
  VMC[0][1]=VMC[1][0]; VMC[0][2]=VMC[2][0]; VMC[0][3]=VMC[3][0];
  VMC[1][2]=VMC[2][1]; VMC[1][3]=VMC[3][1];
  VMC[2][3]=VMC[3][2];
  

/////////////////////// VMC matrix inversion ///////////////////////////////////  
 
  d=1.;
  
  for(k=0; k<n; k++) {
    ll[k]=k;
    mm[k]=k;
    big=VMC[k][k];
    for(j=k; j<n ; j++) {
      for (i=j; i<n; i++) {
        if(TMath::Abs(big) < TMath::Abs(VMC[i][j]) ) { big=VMC[i][j]; ll[k]=i; mm[k]=j; }
      }
    }    
//
    j= ll[k];
    if(j > k) {
      for(i=0; i<n; i++) { hold=-VMC[k][i]; VMC[k][i]=VMC[j][i]; VMC[j][i]=hold; }
      
    }
//
    i=mm[k];
    if(i > k ) { 
      for(j=0; j<n; j++) { hold=-VMC[j][k]; VMC[j][k]=VMC[j][i]; VMC[j][i]=hold; }
    }
//
    if(!big) {
      d=0.;
      cout << "Singular matrix\n"; 
    }
    for(i=0; i<n; i++) {
      if(i == k) { continue; }    
      VMC[i][k]=VMC[i][k]/(-big);
    }   
//
    for(i=0; i<n; i++) {
      hold=VMC[i][k];
      for(j=0; j<n; j++) {
        if(i == k || j == k) { continue; }
	VMC[i][j]=hold*VMC[k][j]+VMC[i][j];
      }
    }
//  
    for(j=0; j<n; j++) {
      if(j == k) { continue; }
      VMC[k][j]=VMC[k][j]/big;
    }
//
    d=d*big;
//
    VMC[k][k]=1./big;        
  } 
//  
  for(k=n-1; k>=0; k--) {
    i=ll[k];
    if(i > k) {
      for (j=0; j<n; j++) {hold=VMC[j][k]; VMC[j][k]=-VMC[j][i]; VMC[j][i]=hold;}
    }  
    j=mm[k];
    if(j > k) {
      for (i=0; i<n; i++) {hold=VMC[k][i]; VMC[k][i]=-VMC[j][i]; VMC[j][i]=hold;}
      }
  }


////////////////////////////////////////////////////////////////////////////////

  Double_t chi2=(m[0]-X0)*( VMC[0][0]*(m[0]-X0) + 2.*VMC[1][0]*(m[1]-X1) + 
                   2.*VMC[2][0]*(m[2]-X2)+ 2.*VMC[3][0]*(m[3]-X3) ) +
                (m[1]-X1)* ( VMC[1][1]*(m[1]-X1) + 2.*VMC[2][1]*(m[2]-X2)+ 
		   2.*VMC[3][1]*(m[3]-X3) ) +
                (m[2]-X2)* ( VMC[2][2]*(m[2]-X2)+ 2.*VMC[3][2]*(m[3]-X3) ) +
                (m[3]-X3)*VMC[3][3]*(m[3]-X3);	
 	 
  newTrack->SetChi2(newTrack->GetChi2()+chi2);
   
} 

