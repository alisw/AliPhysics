#include "AliITSPid.h"
#include "TMath.h"
#include <iostream.h>
ClassImp(AliITSPid)

//-----------------------------------------------------------
Float_t AliITSPid::qcorr(Float_t xc)
{
    assert(0);
  Float_t fcorr;
  fcorr=( 0.766 +0.9692*xc -1.267*xc*xc )*( 1.-TMath::Exp(-xc*64.75) );
  if(fcorr<=0.1)fcorr=0.1;
return qtot/fcorr;
}

//-----------------------------------------------------------
#include "vector.h"
#include <algorithm>
Float_t AliITSPid::qtrm(Int_t track)
{
    TVector q(*( this->GetVec(track)  ));
    Int_t nl=(Int_t)q(0);
    if(nl<1)return 0.;
    vector<float> vf(nl);
    switch(nl){
      case  1:q(6)=q(1); break;
      case  2:q(6)=(q(1)+q(2))/2.; break;
      default:
	for(Int_t i=0;i<nl;i++){vf[i]=q(i+1);}
        sort(vf.begin(),vf.end());
        q(6)= (vf[0]+vf[1])/2.;
        break;
    }
    for(Int_t i=0;i<nl;i++){q(i+1)=vf[i];} this->SetVec(track,q);
    return (q(6));
}
//-----------------------------------------------------------
inline
Int_t	AliITSPid::wpik(Int_t nc,Float_t q)
{
       return pion();   

    Float_t qmpi,qmk,sigpi,sigk,dpi,dk,ppi,pk;
    qmpi =cut[nc][1];
    sigpi=cut[nc][2];
    qmk  =cut[nc][3];
    sigk =cut[nc][4];
    Float_t dqpi=(q-qmpi)/sigpi;
    Float_t dqk =(q-qmk )/sigk;
    if( dqk<-1. )return pion();
    dpi =fabs(dqpi);
    dk  =fabs(dqk);
    ppi =1.- TMath::Erf(dpi);  // +0.5;
    pk  =1.- TMath::Erf(dk);   // +0.5;
    Float_t rpik=ppi/(pk+0.0000001); 
	cout<<"q,dqpi,dqk, wpik: ppi,pk,rpik="
	<<q<<" "<<dqpi<<" "<<dqk<<" "<<ppi<<" "<<pk<<" "<<rpik<<endl;
    if( rpik>=1000. ){return pion();}
    if( rpik>0.001 ){fWk=1./(rpik+1.); fWpi=1.-fWk;return 0;}
    return pion();    
}
//-----------------------------------------------------------
inline
Int_t	AliITSPid::wpikp(Int_t nc,Float_t q)
{
       return pion();   

   Float_t qmpi,qmk,qmp,sigpi,sigk,sigp,ppi,pka,ppr,rppi,rpka;

    qmpi =cut[nc][1];
    sigpi=cut[nc][2];
    qmk  =cut[nc][3];
    sigk =cut[nc][4];
    qmp  =cut[nc][5];
    sigp =cut[nc][6];

    ppi =TMath::Erf( TMath::Abs((q-qmpi)/sigpi) )+0.5;
    pka =TMath::Erf( TMath::Abs((q-qmk )/sigk)  )+0.5;
    ppr =TMath::Erf( TMath::Abs((q-qmp )/sigp)  )+0.5;

    rppi=ppr/ppi;
    rpka=ppr/pka;
    Float_t rp;
    if( rppi>rpka )  { rp=rpka; }
	         else{ rp=rppi; }
    if( rp>=1000. ){ return proton(); }
    if( rp>0.001  ){ fWp=rp/(rp+1.);return proton(); }
    
    return 0;    
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(TClonesArray* rps,Float_t pm)
{
    return 0;    
}
//-----------------------------------------------------------
Int_t	AliITSPid::GetPcode(Float_t q,Float_t pm)
{
    fWpi=fWk=fWp=-1.;     fPcode=0;
//1)
    if ( pm<=cut[1][0] )
	{ return pion(); }
//2)
    if ( pm<=cut[2][0] )
	{ if( q<cut[2][2] ){ return pion(); } else { return  kaon();} }
//3)
    if ( pm<=cut[3][0] )
	if( q<cut[3][2])
		{ return pion(); }
	     else
		{ if ( q<=cut[3][5] ) {return kaon();} else {return proton();}}
//4)
    if ( pm<=cut[4][0] )
	if( q<cut[4][2] )
		{ return pion(); }
	    else
		{ if( q<=cut[4][4] ) {return kaon();} else {return proton(); }}
//5)
    if ( pm<=cut[5][0] )
	if ( q>cut[5][5] ) {return proton();} else {return wpik(5,q);};
//6)
    if ( pm<=cut[6][0] )
	if ( q>cut[6][5] ) {return proton();} else {return wpik(6,q);};
//7)
    if ( pm<=cut[7][0] )
	if ( q<=cut[7][5] ) {return fPcode;} else {return proton();};
//8)
    if ( pm<=cut[8][0] )
	if ( q<=cut[8][5] ) {return fPcode;} else {return proton();}; 
//9)
    if ( pm<=cut[9][0] )
	if ( q<=cut[9][5] ) {return fPcode;} else {return proton();};
//10)
    if( pm<=cut[10][0] ){ return wpikp(10,q); }
//11)
    if( pm<=cut[11][0] ){ return wpikp(11,q); }
//12)
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
//-----------------------------------------------------------
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
		    cout<<xx(0)<<" ";
		    for(Int_t j=1;j<11;j++){
		        cout.width(7);cout.precision(5);cout<<xx(j);
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
    SetCut(  5, 0.47 ,  1.  ,  0.12, 1.98 , 0.17  , 3.5  , inf  );
    SetCut(  6, 0.53 ,  1.  ,  0.12, 1.75 , 0.16  , 3.0  , inf  );
//----------------------------------------------------------------    
    SetCut(  7, 0.59 ,  0.  ,  0.  , 1.18 , 0.125 , 2.7  , inf  );
    SetCut(  8, 0.65 ,  0.  ,  0.  , 1.18 , 0.125 , 2.5  , inf  );
    SetCut(  9, 0.73 ,  0.  ,  0.  , 0.   , 0.125 , 2.0  , inf  );
//----------------------------------------------------------------    
    SetCut( 10, 0.83 ,  0.  ,  0.  , 1.25 , 0.13  , 2.14 , 0.20 );
    SetCut( 11, 0.93 ,  0.  ,  0.  , 1.18 , 0.125 , 1.88 , 0.18 );
    SetCut( 12, 1.03 ,  0.  ,  0.  , 1.13 , 0.12  , 1.68 , 0.155);
}
//-----------------------------------------------------------

