//------------------------------------------------------------//
// Class for identification of pions,kaons and protons in ITS //
// Prior particles population (probabilities) are taken from  //
// Hijing event generator.                                   //
//------------------------------------------------------------//
// #include <stdlib.h>
#include "AliITSPid.h"
//#include "TMath.h"
#include <Riostream.h>
#include <TClonesArray.h>
#include <TVector.h>
#include "AliKalmanTrack.h"
#include "AliITSIOTrack.h"
#include "AliITStrackV2.h"
#include <TF1.h>

ClassImp(AliITSPid)


AliITSPid::AliITSPid(const AliITSPid &source) : TObject(source),
fMxtrs(source.fMxtrs),
fTrs(source.fTrs),
fWpi(source.fWpi),
fWk(source.fWk),
fWp(source.fWp),
fRpik(source.fRpik),
fRppi(source.fRppi),
fRpka(source.fRpka),
fRp(source.fRp),
fPcode(source.fPcode),
fSigmin(source.fSigmin),
fSilent(source.fSilent),
fCutKa(source.fCutKa),
fCutPr(source.fCutPr),
fggpi(source.fggpi),
fggka(source.fggka),
fggpr(source.fggpr){
    // Copy constructor. This is a function which is not allowed to be

  
}
  
//______________________________________________________________________
AliITSPid& AliITSPid::operator=(const AliITSPid& source){
  // Assignment operator. This is a function which is not allowed to be
  this->~AliITSPid();
  new(this) AliITSPid(source);
  return *this;
}

//
  Float_t AliITSPid::Qtrm(Int_t track) {
//
// This calculates truncated mean signal for given track.
//
    TVector q(*( this->GetVec(track)  ));
    Int_t ml=(Int_t)q(0);
    if(ml<1)return 0.;
    if(ml>6)ml=6;
    float vf[6];
    Int_t nl=0; for(Int_t i=0;i<ml;i++){if(q(i)>fSigmin){vf[nl]=q(i+1);nl++;}}
    if(nl==0)return 0.;
    switch(nl){
      case  1:q(6)=q(1); break;
      case  2:q(6)=(q(1)+q(2))/2.; break;
      default:
	for(int fi=0;fi<2;fi++){
	  Int_t swap;
	  do{ swap=0; float qmin=vf[fi];
	  for(int j=fi+1;j<nl;j++)
	    if(qmin>vf[j]){qmin=vf[j];vf[j]=vf[fi];vf[fi]=qmin;swap=1;};
	  }while(swap==1);
	}
        q(6)= (vf[0]+vf[1])/2.;
        break;
    }
    for(Int_t i=0;i<nl;i++){q(i+1)=vf[i];} this->SetVec(track,q);
    return (q(6));
}

Float_t AliITSPid::Qtrm(Float_t qarr[6],Int_t narr) const{
//
//This calculates truncated mean signal for given signal set.
//
  Float_t q[6],qm,qmin;
  Int_t nl,ml;
  if(narr>0&&narr<7){ml=narr;}else{return 0;};
  nl=0; for(Int_t i=0;i<ml;i++){if(qarr[i]>fSigmin){q[nl]=qarr[i];nl++;}}
  if(nl==0)return 0.;
    switch(nl){
      case  1:qm=q[0]; break;
      case  2:qm=(q[0]+q[1])/2.; break;
      default:
	Int_t swap;
	for(int fi=0;fi<2;fi++){
	  do{ swap=0; qmin=q[fi];
	  for(int j=fi+1;j<nl;j++)
	    if(qmin>q[j]){qmin=q[j];q[j]=q[fi];q[fi]=qmin;swap=1;};
	  }while(swap==1);
	}
        qm= (q[0]+q[1])/2.;
        break;
    }
    return qm;
}

Int_t	AliITSPid::Wpik(Float_t pm,Float_t q){
  //Calcutates probabilityes of pions and kaons
  //Returns particle code for dominant probability.
  Double_t par[6];
  for(int i=0;i<6;i++){par[i]=fGGpi[i]->Eval(pm);}
  fggpi->SetParameters(par);

  for(int i=0;i<3;i++){par[i]=fGGka[i]->Eval(pm);}
  fggka->SetParameters(par);

  Float_t ppi=fggpi->Eval(q);
  Float_t pka=fggka->Eval(q);
  Float_t p=ppi+pka;
  /*
  if(!fSilent){
    fggka->Print();
    fggpi->Print();
    if(p>0)cout<<" ppi,pka="<<ppi/p<<"  "<<pka/p<<endl;
  }
  */

  if(p>0){
    ppi=ppi/p; 
    pka=pka/p;
    fWp=0.; fWpi=ppi; fWk=pka;
    if( pka>ppi){return fPcode=321;}else{return fPcode=211;}
  }else{return 0;}
}
//-----------------------------------------------------------
Int_t	AliITSPid::Wpikp(Float_t pm,Float_t q){
  //
  //Calcutates probabilityes of pions,kaons and protons.
  //Returns particle code for dominant probability.
  Double_t par[6];
  for(int i=0;i<6;i++){par[i]=fGGpi[i]->Eval(pm);}
  fggpi->SetParameters(par);

  for(int i=0;i<3;i++){par[i]=fGGka[i]->Eval(pm);}
  fggka->SetParameters(par);

  for(int i=0;i<3;i++){par[i]=fGGpr[i]->Eval(pm);}
  fggpr->SetParameters(par);

  Float_t p,ppi,pka,ppr;
  if( q>(fggpr->GetParameter(1)+fggpr->GetParameter(2)) )
      { p=1.0; ppr=1.0; ppi=pka=0.0;
    }else{ 
    ppi=fggpi->Eval(q);
    pka=fggka->Eval(q);
    ppr=fggpr->Eval(q);
    p=ppi+pka+ppr;
  }
  if(p>0){
    ppi=ppi/p; 
    pka=pka/p;
    ppr=ppr/p;
    fWp=ppr; fWpi=ppi; fWk=pka;
    //if(!fSilent)cout<<" ppi,pka,ppr="<<ppi<<"  "<<pka<<" "<<ppr<<endl;

   if( ppi>pka&&ppi>ppr )
           {return fPcode=211;}
   else{ if(pka>ppr){return fPcode=321;}else{return fPcode=2212;}
   }

  }else{return 0;}
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(TClonesArray* rps,Float_t pm)
{
  //Returns particle code
  Info("GetPcode","method not implemented - Inputs TClonesArray *%x , Float_t %f",rps,pm); 
    return 0;    
}
//-----------------------------------------------------------
Int_t   AliITSPid::GetPcode(AliKalmanTrack *track)
{
  //Returns particle code for given track.
      Double_t xk,par[5]; track->GetExternalParameters(xk,par);
      Float_t phi=TMath::ASin(par[2]) + track->GetAlpha();
      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
      Float_t lam=TMath::ATan(par[3]); 
      Float_t pt1=TMath::Abs(par[4]);
      Float_t mom=1./(pt1*TMath::Cos(lam));
      Float_t dedx=track->GetPIDsignal();
    Int_t pcode=GetPcode(dedx/40.,mom);
//    cout<<"TPCtrack dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
    return pcode?pcode:211;
    }
//------------------------------------------------------------
Int_t   AliITSPid::GetPcode(AliITSIOTrack *track)
{
  //Returns particle code for given track(V1).
    Double_t px,py,pz;
    px=track->GetPx();
    py=track->GetPy();
    pz=track->GetPz();
    Float_t mom=TMath::Sqrt(px*px+py*py+pz*pz);
//???????????????????
    // Float_t dedx=1.0;
    Float_t dedx=track->GetdEdx();
//???????????????????    
    Int_t pcode=GetPcode(dedx,mom);
//    cout<<"ITSV1 dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
return pcode?pcode:211;
}
//-----------------------------------------------------------
Int_t   AliITSPid::GetPcode(AliITStrackV2 *track)
{
//Returns particle code for given track(V2).
  if(track==0)return 0;
  //      track->Propagate(track->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
      track->PropagateTo(3.,0.0028,65.19);
      //track->PropagateToVertex();          Not needed. (I.B.)
    Double_t xk,par[5]; track->GetExternalParameters(xk,par);
    Float_t lam=TMath::ATan(par[3]);
    Float_t pt1=TMath::Abs(par[4]);
    Float_t mom=0.;
    if( (pt1*TMath::Cos(lam))!=0. ){ mom=1./(pt1*TMath::Cos(lam)); }else{mom=0.;};
    Float_t dedx=track->GetdEdx();
//    cout<<"lam,pt1,mom,dedx="<<lam<<","<<pt1<<","<<mom<<","<<dedx<<endl;
    Int_t pcode=GetPcode(dedx,mom);
//    cout<<"ITS V2 dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
return pcode?pcode:211;
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(Float_t q,Float_t pm)
{
  //Returns particle code for given signal and momentum.
    fWpi=fWk=fWp=0.;     fPcode=0;

    if ( pm<=0.400 )
	{ if( q<fCutKa->Eval(pm) )
	    {return Pion();}
	else{ if( q<fCutPr->Eval(pm) )
		{return Kaon();}
	    else{return Proton();}
	    } 
	}
    if ( pm<=0.750 )
	if ( q>fCutPr->Eval(pm)  )
	    {return Proton();} else {return Wpik(pm,q);};
    if( pm<=1.10 ){ return Wpikp(pm,q); }
    return fPcode;    
}
//-----------------------------------------------------------
void	AliITSPid::SetCut(Int_t n,Float_t pm,Float_t pilo,Float_t pihi,
			Float_t klo,Float_t khi,Float_t plo,Float_t phi)
{
  // The cut-table initializer method.
    fCut[n][0]=pm;
    fCut[n][1]=pilo;
    fCut[n][2]=pihi;
    fCut[n][3]=klo;
    fCut[n][4]=khi;
    fCut[n][5]=plo;
    fCut[n][6]=phi;
    return ;    
}
//------------------------------------------------------------
void AliITSPid::SetVec(Int_t ntrack,const TVector& info) const
{
  //Store track info in tracls table
TClonesArray& arr=*fTrs;
    new( arr[ntrack] ) TVector(info);
}
//-----------------------------------------------------------
TVector* AliITSPid::GetVec(Int_t ntrack) const
{
  //Get given track from track table 
TClonesArray& arr=*fTrs;
    return (TVector*)arr[ntrack];
}
//-----------------------------------------------------------
void AliITSPid::SetEdep(Int_t track,Float_t Edep)
{
  //Set dEdx for given track
    TVector xx(0,11);
    if( ((TVector*)fTrs->At(track))->IsValid() )
	{TVector yy( *((TVector*)fTrs->At(track)) );xx=yy; }
    Int_t j=(Int_t)xx(0); if(j>4)return;
    xx(++j)=Edep;xx(0)=j;
    TClonesArray &arr=*fTrs;
    new(arr[track])TVector(xx);
}
//-----------------------------------------------------------
void AliITSPid::SetPmom(Int_t track,Float_t Pmom)
{
  //Set momentum for given track
    TVector xx(0,11);
    if( ((TVector*)fTrs->At(track))->IsValid() )
	{TVector yy( *((TVector*)fTrs->At(track)) );xx=yy; }
    xx(10)=Pmom;
    TClonesArray &arr=*fTrs;
    new(arr[track])TVector(xx);
}
//-----------------------------------------------------------
void AliITSPid::SetPcod(Int_t track,Int_t partcode)
{
  //Set particle code for given track
    TVector xx(0,11);
    if( ((TVector*)fTrs->At(track))->IsValid() )
	{TVector yy( *((TVector*)fTrs->At(track)) );xx=yy; }
    if(xx(11)==0)
	{xx(11)=partcode; fMxtrs++;
	TClonesArray &arr=*fTrs;
	new(arr[track])TVector(xx);
	}
}
//-----------------------------------------------------------
void AliITSPid::Print(Int_t track)
{
  //Prints information for given track
cout<<fMxtrs<<" tracks in AliITSPid obj."<<endl;
    if( ((TVector*)fTrs->At(track))->IsValid() )
	{TVector xx( *((TVector*)fTrs->At(track)) );
	 xx.Print();
	 }
    else 
	{cout<<"No data for track "<<track<<endl;return;}
}
//-----------------------------------------------------------
void AliITSPid::Tab(void)
{
  //Make PID for tracks stored in tracks table
if(fTrs->GetEntries()==0){cout<<"No entries in TAB"<<endl;return;}
cout<<"------------------------------------------------------------------------"<<endl;
cout<<"Nq"<<"   q1  "<<"   q2  "<<"   q3  "<<"   q4  "<<"   q5   "<<
      " Qtrm    "    <<"  Wpi  "<<"  Wk   "<<"  Wp  "<<"Pmom  "<<endl;
cout<<"------------------------------------------------------------------------"<<endl;
for(Int_t i=0;i<fTrs->GetEntries();i++)
{
  TVector xx( *((TVector*)fTrs->At(i)) );     
    if( xx.IsValid() && xx(0)>0 )
	{
	    TVector xx( *((TVector*)fTrs->At(i)) );
	    if(xx(0)>=2)
		{
//       1)Calculate Qtrm	
		    xx(6)=(this->Qtrm(i));

	    	 }else{
		     xx(6)=xx(1);
		 }
//	 2)Calculate Wpi,Wk,Wp
  	    this->GetPcode(xx(6),xx(10)/1000.);
	    xx(7)=GetWpi();
	    xx(8)=GetWk();
	    xx(9)=GetWp();
//       3)Print table
	    if(xx(0)>0){
//		    cout<<xx(0)<<" ";
		    for(Int_t j=1;j<11;j++){
		      if(i<7){ cout.width(7);cout.precision(4);cout<<xx(j);}
                      if(i>7){ cout.width(7);cout.precision(5);cout<<xx(j);}
		    }
		    cout<<endl;
		}
//	  4)Update data in TVector
	    TClonesArray &arr=*fTrs;
	    new(arr[i])TVector(xx);	 
	}
    else 
      {/*cout<<"No data for track "<<i<<endl;*/}
}// End loop for tracks
}
void AliITSPid::Reset(void)
{
  //Reset tracks table
  for(Int_t i=0;i<fTrs->GetEntries();i++){
    TVector xx(0,11);
    TClonesArray &arr=*fTrs;
    new(arr[i])TVector(xx);
  }
}
//-----------------------------------------------------------
AliITSPid::AliITSPid(Int_t ntrack):
fMxtrs(0),
fTrs(0),
fWpi(0),
fWk(0),
fWp(0),
fRpik(0),
fRppi(0),
fRpka(0),
fRp(0),
fPcode(0),
fSigmin(0.01),
fSilent(0),
fCutKa(0),
fCutPr(0),
fggpi(0),
fggka(0),
fggpr(0){
  //Constructor for AliITSPid class
 
    fTrs = new TClonesArray("TVector",ntrack);
    TClonesArray &arr=*fTrs;
    for(Int_t i=0;i<ntrack;i++)new(arr[i])TVector(0,11);
     //   
    fCutKa=new TF1("fcutka","pol4",0.05,0.4);
    Double_t ka[5]={25.616, -161.59, 408.97, -462.17, 192.86};
    fCutKa->SetParameters(ka);
    //
    fCutPr=new TF1("fcutpr","[0]/x/x+[1]",0.05,1.1);
    Double_t pr[2]={0.70675,0.4455};
    fCutPr->SetParameters(pr);
    //
    //---------- signal fit ----------
{//Pions
fGGpi[0]=new TF1("fp1pi","pol4",0.34,1.2);
  Double_t parpi0[10]={ -1.9096471071e+03, 4.5354331545e+04, -1.1860738840e+05,
   1.1405329025e+05, -3.8289694496e+04  };
  fGGpi[0]->SetParameters(parpi0);
fGGpi[1]=new TF1("fp2pi","[0]/x/x+[1]",0.34,1.2);
  Double_t parpi1[10]={ 1.0791668283e-02, 9.7347716496e-01  };
  fGGpi[1]->SetParameters(parpi1);
fGGpi[2]=new TF1("fp3pi","[0]/x/x+[1]",0.34,1.2);
  Double_t parpi2[10]={ 5.8191602279e-04, 9.7285601334e-02  };
  fGGpi[2]->SetParameters(parpi2);
fGGpi[3]=new TF1("fp4pi","pol4",0.34,1.2);
  Double_t parpi3[10]={ 6.6267353195e+02, 7.1595101104e+02, -5.3095111914e+03,
   6.2900977606e+03, -2.2935862292e+03  };
  fGGpi[3]->SetParameters(parpi3);
fGGpi[4]=new TF1("fp5pi","[0]/x/x+[1]",0.34,1.2);
  Double_t parpi4[10]={ 9.0419011783e-03, 1.1628922525e+00  };
  fGGpi[4]->SetParameters(parpi4);
fGGpi[5]=new TF1("fp6pi","[0]/x/x+[1]",0.34,1.2);
  Double_t parpi5[10]={ 1.8324872519e-03, 2.1503968838e-01  };
  fGGpi[5]->SetParameters(parpi5);
}//End Pions
{//Kaons
fGGka[0]=new TF1("fp1ka","pol4",0.24,1.2);
  Double_t parka0[20]={
  -1.1204243395e+02,4.6716191428e+01,2.2584059281e+03,
  -3.7123338009e+03,1.6003647641e+03  };
  fGGka[0]->SetParameters(parka0);
fGGka[1]=new TF1("fp2ka","[0]/x/x+[1]",0.24,1.2);
  Double_t parka1[20]={
  2.5181172905e-01,8.7566001814e-01  };
  fGGka[1]->SetParameters(parka1);
fGGka[2]=new TF1("fp3ka","pol6",0.24,1.2);
  Double_t parka2[20]={
  8.6236021573e+00,-7.0970427531e+01,2.4846827669e+02,
  -4.6094401290e+02,4.7546751408e+02,-2.5807112462e+02,
  5.7545491696e+01  };
  fGGka[2]->SetParameters(parka2);
}//End Kaons
{//Protons
fGGpr[0]=new TF1("fp1pr","pol4",0.4,1.2);
  Double_t parpr0[10]={
  6.0150106543e+01,-8.8176206410e+02,3.1222644604e+03,
  -3.5269200901e+03,1.2859128345e+03  };
  fGGpr[0]->SetParameters(parpr0);
fGGpr[1]=new TF1("fp2pr","[0]/x/x+[1]",0.4,1.2);
  Double_t parpr1[10]={
  9.4970837607e-01,7.3573504201e-01  };
  fGGpr[1]->SetParameters(parpr1);
fGGpr[2]=new TF1("fp3pr","[0]/x/x+[1]",0.4,1.2);
  Double_t parpr2[10]={
  1.2498403757e-01,2.7845072306e-02  };
  fGGpr[2]->SetParameters(parpr2);
}//End Protons
    //----------- end fit -----------

    fggpr=new TF1("ggpr","gaus",0.4,1.2);
    fggpi=new TF1("ggpi","gaus+gaus(3)",0.4,1.2);
    fggka=new TF1("ggka","gaus",0.4,1.2);

    //-------------------------------------------------
const int kInf=10;
//         Ncut Pmom   pilo  pihi    klo    khi     plo    phi
//       cut[j] [0]    [1]    [2]    [3]    [4]     [5]    [6]
//----------------------------------------------------------------
    SetCut(  1, 0.12 ,  0.  ,  0.  , kInf  , kInf   , kInf  , kInf  );
    SetCut(  2, 0.20 ,  0.  ,  6.0 , 6.0  , kInf   , kInf  , kInf  );
    SetCut(  3, 0.30 ,  0.  ,  3.5 , 3.5  , 9.0   , 9.0  , kInf  );
    SetCut(  4, 0.41 ,  0.  ,  1.9 , 1.9  , 4.0   , 4.0  , kInf  );
//----------------------------------------------------------------
    SetCut(  5, 0.47 , 0.935, 0.139, 1.738 , 0.498  , 3.5  , kInf  );  //410-470
    SetCut(  6, 0.53 , 0.914, 0.136, 1.493 , 0.436  , 3.0  , kInf  );  //470-530
//----------------------------------------------------------------    
    SetCut(  7, 0.59 , 0.895, 0.131, 1.384 , 0.290 , 2.7  , kInf  );    //530-590
    SetCut(  8, 0.65 , 0.887, 0.121, 1.167 , 0.287 , 2.5  , kInf  );     //590-650
    SetCut(  9, 0.73 , 0.879, 0.120, 1.153 , 0.257 , 2.0  , kInf  );     //650-730
//----------------------------------------------------------------    
    SetCut( 10, 0.83 , 0.880, 0.126, 1.164 , 0.204 , 2.308 , 0.297 );       //730-830
    SetCut( 11, 0.93 , 0.918, 0.145, 1.164 , 0.204 , 2.00 , 0.168 );        //830-930
    SetCut( 12, 1.03 , 0.899, 0.128, 1.164 , 0.204  ,1.80 , 0.168);
    //------------------------ pi,K ---------------------
    fAprob[0][0]=1212;     fAprob[1][0]=33.;   // fAprob[0][i] - const for pions,cut[i+5] 
    fAprob[0][1]=1022;     fAprob[1][1]=46.2 ; // fAprob[1][i] -           kaons
    //---------------------------------------------------
    fAprob[0][2]= 889.7;   fAprob[1][2]=66.58; fAprob[2][2]=14.53;
    fAprob[0][3]= 686.;    fAprob[1][3]=88.8;  fAprob[2][3]=19.27;   
    fAprob[0][4]= 697.;    fAprob[1][4]=125.6; fAprob[2][4]=28.67;
    //------------------------ pi,K,p -------------------
    fAprob[0][5]= 633.7;   fAprob[1][5]=100.1;   fAprob[2][5]=37.99;   // fAprob[2][i] -  protons
    fAprob[0][6]= 469.5;   fAprob[1][6]=20.74;   fAprob[2][6]=25.43;
    fAprob[0][7]= 355.;    fAprob[1][7]=
                          355.*(20.74/469.5);  fAprob[2][7]=34.08;
}
//End AliITSPid.cxx
