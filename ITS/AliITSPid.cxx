#include "AliITSPid.h"
#include "TMath.h"
#include "AliITSIOTrack.h"
#include <Riostream.h>
ClassImp(AliITSPid)

Float_t AliITSPid::qcorr(Float_t xc)
{
    assert(0);
  Float_t fcorr;
  fcorr=( 0.766 +0.9692*xc -1.267*xc*xc )*( 1.-TMath::Exp(-xc*64.75) );
  if(fcorr<=0.1)fcorr=0.1;
return qtot/fcorr;
}

Float_t AliITSPid::qtrm(Int_t track)
{
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

Float_t AliITSPid::qtrm(Float_t qarr[6],Int_t narr)
{
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
//-----------------------------------------------------------
//inline
Int_t	AliITSPid::wpik(Int_t nc,Float_t q)
{
  //         return pion();   

    Float_t qmpi,qmk,sigpi,sigk,dpi,dk,ppi,pk;
    Float_t appi,apk;
    qmpi =cut[nc][1];
    sigpi=cut[nc][2];
    qmk  =cut[nc][3];
    sigk =cut[nc][4];
    appi = aprob[0][nc-5];
    apk  = aprob[1][nc-5];
//    cout<<"qmpi,sigpi,qmk,sigk="<<qmpi<<"  "<<sigpi<<"  "<<qmk<<"  "<<sigk<<endl;
//    cout<<"appi,apk="<<appi<<","<<apk<<endl;
    Float_t dqpi=(q-qmpi)/sigpi;
    Float_t dqk =(q-qmk )/sigk;
    if( dqk<-1. )return pion();
    dpi =TMath::Abs(dqpi);
    dk  =TMath::Abs(dqk);
    //ppi =1.- TMath::Erf(dpi);  // +0.5;
    //pk  =1.- TMath::Erf(dk);   // +0.5;
    ppi=appi*TMath::Gaus(q,qmpi,sigpi)
        /(appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk));
    pk = apk*TMath::Gaus(q,qmk, sigk )
        /(appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk));

//    Float_t rpik=ppi/(pk+0.0000001); 
//	cout<<"q,dqpi,dqk, wpik: ppi,pk,rpik="
//	<<q<<"  "<<dqpi<<"  "<<dqk<<"  "<<ppi<<"  "<<pk<<"  "<<rpik<<endl;

    fWp=0.; fWpi=ppi; fWk=pk;
    if( pk>ppi){return kaon();}else{return pion();}
}
//-----------------------------------------------------------
//inline
Int_t	AliITSPid::wpikp(Int_t nc,Float_t q)
{
   Float_t qmpi,qmk,qmp,sigpi,sigk,sigp,ppi,pk,pp;
   Float_t appi,apk,app;

    qmpi =cut[nc][1];
    sigpi=cut[nc][2];
    qmk  =cut[nc][3];
    sigk =cut[nc][4];
    qmp  =cut[nc][5];
    sigp =cut[nc][6];

    //appi = apk = app = 1.;
    appi = aprob[0][nc-5];
    apk  = aprob[1][nc-5];
    app  = aprob[2][nc-5];

    //ppi =TMath::Erf( TMath::Abs((q-qmpi)/sigpi) )+0.5;
    //pka =TMath::Erf( TMath::Abs((q-qmk )/sigk)  )+0.5;
    //ppr =TMath::Erf( TMath::Abs((q-qmp )/sigp)  )+0.5;

    ppi=appi*TMath::Gaus(q,qmpi,sigpi)
      /( appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk)+app*TMath::Gaus(q,qmp,sigp) );
    pk = apk*TMath::Gaus(q,qmk, sigk )
      /( appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk)+app*TMath::Gaus(q,qmp,sigp) );
    pp = app*TMath::Gaus(q,qmp, sigp )
      /( appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk)+app*TMath::Gaus(q,qmp,sigp) );

    fWp=pp; fWpi=ppi; fWk=pk;

//cout<<" wpikp: mid,sig pi,k,p="<<qmpi<<" "<<sigpi<<";   "<<qmk<<" "<<sigk<<";   "
//    <<qmp<<" "<<sigp<<"; "<<endl;
// cout<<" aprob: "<<appi<<"  "<<apk<<"  "<<app<<endl;
//cout<<" ppi,pk,pp="<<ppi<<"  "<<pk<<"  "<<pp<<endl;

    if( ppi>pk&&ppi>pp )  { return pion(); }
    if(pk>pp){return kaon();}else{return proton();}
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(TClonesArray* rps,Float_t pm)
{
    return 0;    
}
//-----------------------------------------------------------
Int_t   AliITSPid::GetPcode(AliTPCtrack *track)
{
      Double_t xk,par[5]; track->GetExternalParameters(xk,par);
      Float_t phi=TMath::ASin(par[2]) + track->GetAlpha();
      if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
      if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
      Float_t lam=TMath::ATan(par[3]); 
      Float_t pt_1=TMath::Abs(par[4]);
      Float_t mom=1./(pt_1*TMath::Cos(lam));
      Float_t dedx=track->GetdEdx();
    Int_t pcode=GetPcode(dedx/40.,mom);
//    cout<<"TPCtrack dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
    return pcode?pcode:211;
    }
//------------------------------------------------------------
Int_t   AliITSPid::GetPcode(AliITSIOTrack *track)
{
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
  if(track==0)return 0;
  //      track->Propagate(track->GetAlpha(),3.,0.1/65.19*1.848,0.1*1.848);
      track->PropagateTo(3.,0.0028,65.19);
      track->PropagateToVertex();
    Double_t xk,par[5]; track->GetExternalParameters(xk,par);
    Float_t lam=TMath::ATan(par[3]);
    Float_t pt_1=TMath::Abs(par[4]);
    Float_t mom=0.;
    if( (pt_1*TMath::Cos(lam))!=0. ){ mom=1./(pt_1*TMath::Cos(lam)); }else{mom=0.;};
    Float_t dedx=track->GetdEdx();
//    cout<<"lam,pt_1,mom,dedx="<<lam<<","<<pt_1<<","<<mom<<","<<dedx<<endl;
    Int_t pcode=GetPcode(dedx,mom);
//    cout<<"ITS V2 dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
return pcode?pcode:211;
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(Float_t q,Float_t pm)
{
    fWpi=fWk=fWp=0.;     fPcode=0;
//1)---------------------- 0-120 MeV/c --------------
    if ( pm<=cut[1][0] )
	{ return pion(); }
//2)----------------------120-200 Mev/c ------------- 
    if ( pm<=cut[2][0] )
	{ if( q<cut[2][2] ){ return pion(); } else { return  kaon();} }
//3)----------------------200-300 Mev/c -------------
    if ( pm<=cut[3][0] )
	if( q<cut[3][2])
		{ return pion(); }
	     else
		{ if ( q<=cut[3][5] ) {return kaon();} else {return proton();}}
//4)----------------------300-410 Mev/c -------------
    if ( pm<=cut[4][0] )
	if( q<cut[4][2] )
		{ return pion(); }
	    else
		{ if( q<=cut[4][4] ) {return kaon();} else {return proton(); }}
//5)----------------------410-470 Mev/c -------------
    if ( pm<=cut[5][0] )
	if ( q>cut[5][5] ) {return proton();} else {return wpik(5,q);};
//6)----------------------470-530 Mev/c -------------
    if ( pm<=cut[6][0] )
	if ( q>cut[6][5] ) {return proton();} else {return wpik(6,q);};
//7)----------------------530-590 Mev/c -------------
    if ( pm<=cut[7][0] )
	if ( q<=cut[7][5] ) {return wpik(7,q);} else {return proton();};
//8)----------------------590-650 Mev/c -------------
    if ( pm<=cut[8][0] )
	if ( q<=cut[8][5] ) {return wpik(8,q);} else {return proton();}; 
//9)----------------------650-730 Mev/c -------------
    if ( pm<=cut[9][0] )
	if ( q<=cut[9][5] ) {return wpik(9,q);} else {return proton();};
//10)----------------------730-830 Mev/c -------------
    if( pm<=cut[10][0] ){ return wpikp(10,q); }
//11)----------------------830-930 Mev/c -------------
    if( pm<=cut[11][0] ){ return wpikp(11,q); }
//12)----------------------930-1030 Mev/c -------------
    if( pm<=cut[12][0] ){ return wpikp(12,q); }

    return fPcode;    
}
//-----------------------------------------------------------
void	AliITSPid::SetCut(Int_t n,Float_t pm,Float_t pilo,Float_t pihi,
			Float_t klo,Float_t khi,Float_t plo,Float_t phi)
{
    cut[n][0]=pm;
    cut[n][1]=pilo;
    cut[n][2]=pihi;
    cut[n][3]=klo;
    cut[n][4]=khi;
    cut[n][5]=plo;
    cut[n][6]=phi;
    return ;    
}
//------------------------------------------------------------
void AliITSPid::SetVec(Int_t ntrack,TVector info)
{
TClonesArray& arr=*trs;
    new( arr[ntrack] ) TVector(info);
}
//-----------------------------------------------------------
TVector* AliITSPid::GetVec(Int_t ntrack)
{
TClonesArray& arr=*trs;
    return (TVector*)arr[ntrack];
}
//-----------------------------------------------------------
void AliITSPid::SetEdep(Int_t track,Float_t Edep)
{
    TVector xx(0,11);
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector yy( *((TVector*)trs->At(track)) );xx=yy; }
    Int_t j=(Int_t)xx(0); if(j>4)return;
    xx(++j)=Edep;xx(0)=j;
    TClonesArray &arr=*trs;
    new(arr[track])TVector(xx);
}
//-----------------------------------------------------------
void AliITSPid::SetPmom(Int_t track,Float_t Pmom)
{
    TVector xx(0,11);
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector yy( *((TVector*)trs->At(track)) );xx=yy; }
    xx(10)=Pmom;
    TClonesArray &arr=*trs;
    new(arr[track])TVector(xx);
}
//-----------------------------------------------------------
void AliITSPid::SetPcod(Int_t track,Int_t partcode)
{
    TVector xx(0,11);
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector yy( *((TVector*)trs->At(track)) );xx=yy; }
    if(xx(11)==0)
	{xx(11)=partcode; mxtrs++;
	TClonesArray &arr=*trs;
	new(arr[track])TVector(xx);
	}
}
//-----------------------------------------------------------
void AliITSPid::Print(Int_t track)
{cout<<mxtrs<<" tracks in AliITSPid obj."<<endl;
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector xx( *((TVector*)trs->At(track)) );
	 xx.Print();
	 }
    else 
	{cout<<"No data for track "<<track<<endl;return;}
}
//-----------------------------------------------------------
void AliITSPid::Tab(void)
{
if(trs->GetEntries()==0){cout<<"No entries in TAB"<<endl;return;}
cout<<"------------------------------------------------------------------------"<<endl;
cout<<"Nq"<<"   q1  "<<"   q2  "<<"   q3  "<<"   q4  "<<"   q5   "<<
      " Qtrm    "    <<"  Wpi  "<<"  Wk   "<<"  Wp  "<<"Pmom  "<<endl;
cout<<"------------------------------------------------------------------------"<<endl;
for(Int_t i=0;i<trs->GetEntries();i++)
{
  TVector xx( *((TVector*)trs->At(i)) );     
    if( xx.IsValid() && xx(0)>0 )
	{
	    TVector xx( *((TVector*)trs->At(i)) );
	    if(xx(0)>=2)
		{
//       1)Calculate Qtrm	
		    xx(6)=(this->qtrm(i));

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
	    TClonesArray &arr=*trs;
	    new(arr[i])TVector(xx);	 
	}
    else 
      {/*cout<<"No data for track "<<i<<endl;*/}
}// End loop for tracks
}
void AliITSPid::Reset(void)
{
  for(Int_t i=0;i<trs->GetEntries();i++){
    TVector xx(0,11);
    TClonesArray &arr=*trs;
    new(arr[i])TVector(xx);
  }
}
//-----------------------------------------------------------
AliITSPid::AliITSPid(Int_t ntrack)
{
  fSigmin=0.01;
    trs = new TClonesArray("TVector",ntrack);
    TClonesArray &arr=*trs;
    for(Int_t i=0;i<ntrack;i++)new(arr[i])TVector(0,11);
    mxtrs=0;
const int inf=10;
//         Ncut Pmom   pilo  pihi    klo    khi     plo    phi
//       cut[j] [0]    [1]    [2]    [3]    [4]     [5]    [6]
//----------------------------------------------------------------
    SetCut(  1, 0.12 ,  0.  ,  0.  , inf  , inf   , inf  , inf  );
    SetCut(  2, 0.20 ,  0.  ,  6.0 , 6.0  , inf   , inf  , inf  );
    SetCut(  3, 0.30 ,  0.  ,  3.5 , 3.5  , 9.0   , 9.0  , inf  );
    SetCut(  4, 0.41 ,  0.  ,  1.9 , 1.9  , 4.0   , 4.0  , inf  );
//----------------------------------------------------------------
    SetCut(  5, 0.47 , 0.935, 0.139, 1.738 , 0.498  , 3.5  , inf  );  //410-470
    SetCut(  6, 0.53 , 0.914, 0.136, 1.493 , 0.436  , 3.0  , inf  );  //470-530
//----------------------------------------------------------------    
    SetCut(  7, 0.59 , 0.895, 0.131, 1.384 , 0.290 , 2.7  , inf  );    //530-590
    SetCut(  8, 0.65 , 0.887, 0.121, 1.167 , 0.287 , 2.5  , inf  );     //590-650
    SetCut(  9, 0.73 , 0.879, 0.120, 1.153 , 0.257 , 2.0  , inf  );     //650-730
//----------------------------------------------------------------    
    SetCut( 10, 0.83 , 0.880, 0.126, 1.164 , 0.204 , 2.308 , 0.297 );       //730-830
    SetCut( 11, 0.93 , 0.918, 0.145, 1.164 , 0.204 , 2.00 , 0.168 );        //830-930
    SetCut( 12, 1.03 , 0.899, 0.128, 1.164 , 0.204  ,1.80 , 0.168);
    //------------------------ pi,K ---------------------
    aprob[0][0]=1212;     aprob[1][0]=33.;   // aprob[0][i] - const for pions,cut[i+5] 
    aprob[0][1]=1022;     aprob[1][1]=46.2 ; // aprob[1][i] -           kaons
    //---------------------------------------------------
    aprob[0][2]= 889.7;   aprob[1][2]=66.58; aprob[2][2]=14.53;
    aprob[0][3]= 686.;    aprob[1][3]=88.8;  aprob[2][3]=19.27;   
    aprob[0][4]= 697.;    aprob[1][4]=125.6; aprob[2][4]=28.67;
    //------------------------ pi,K,p -------------------
    aprob[0][5]= 633.7;   aprob[1][5]=100.1;   aprob[2][5]=37.99;   // aprob[2][i] -  protons
    aprob[0][6]= 469.5;   aprob[1][6]=20.74;   aprob[2][6]=25.43;
    aprob[0][7]= 355.;    aprob[1][7]=
                          355.*(20.74/469.5);  aprob[2][7]=34.08;
}
//-----------------------------------------------------------

