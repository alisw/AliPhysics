/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include "AliTPCPid.h"
#include "TMath.h"
//#include "AliITSIOTrack.h"
#include "Riostream.h"
ClassImp(AliTPCPid)
// Correction 13.01.2003 Z.S.,Dubna
//            22.01.2003
//------------------------------------------------------------
Float_t AliTPCPid::qcorr(Float_t xc)
{
    assert(0);
  Float_t fcorr;
  fcorr=( 0.766 +0.9692*xc -1.267*xc*xc )*( 1.-TMath::Exp(-xc*64.75) );
  if(fcorr<=0.1)fcorr=0.1;
return qtot/fcorr;
}
//-----------------------------------------------------------
Float_t AliTPCPid::qtrm(Int_t track)
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
}// --- End AliTPCPid::qtrm ---

Float_t AliTPCPid::qtrm(Float_t qarr[6],Int_t narr)
{
  //..................
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
}// --- End AliTPCPid::qtrm  ---

Int_t	AliTPCPid::wpik(Int_t nc,Float_t q)
{
    Float_t qmpi,qmk,sigpi,sigk,dpi,dk,ppi,pk;
    Float_t appi,apk;
    qmpi =cut[nc][1];
    sigpi=cut[nc][2];
    qmk  =cut[nc][3];
    sigk =cut[nc][4];
    appi = aprob[0][nc-5];
    apk  = aprob[1][nc-5];
    if( !fSilent ){
    cout<<"qmpi,sigpi,qmk,sigk="<<qmpi<<"  "<<sigpi<<"  "<<qmk<<"  "<<sigk<<endl;
    cout<<"appi,apk="<<appi<<","<<apk<<endl;
    }
    Float_t dqpi=(q-qmpi)/sigpi;
    Float_t dqk =(q-qmk )/sigk;
    dpi =TMath::Abs(dqpi);
    dk  =TMath::Abs(dqk);
    Double_t dn=appi*TMath::Gaus(q,qmpi,sigpi)+apk*TMath::Gaus(q,qmk,sigk);
    if(dn>0.){
      ppi=appi*TMath::Gaus(q,qmpi,sigpi)/dn;
      pk = apk*TMath::Gaus(q,qmk, sigk )/dn;
    }else{fWpi=1;return pion();}
    Float_t rpik=ppi/(pk+0.0000001); 
    if( !fSilent )
	cout<<"q,dqpi,dqk, wpik: ppi,pk,rpik="
	<<q<<"  "<<dqpi<<"  "<<dqk<<"  "<<ppi<<"  "<<pk<<"  "<<rpik<<endl;

    fWp=0.; fWpi=ppi; fWk=pk;
    if( pk>ppi){return kaon();}else{return pion();}
}
//-----------------------------------------------------------
//inline
Int_t	AliTPCPid::wpikp(Int_t nc,Float_t q)
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

    Double_t dn= appi*TMath::Gaus(q,qmpi,sigpi)
      +apk*TMath::Gaus(q,qmk,sigk)+app*TMath::Gaus(q,qmp,sigp);
    if(dn>0.){
    ppi=appi*TMath::Gaus(q,qmpi,sigpi)/dn;
    pk = apk*TMath::Gaus(q,qmk, sigk )/dn;
    pp = app*TMath::Gaus(q,qmp, sigp )/dn;
     }
    else{fWpi=1;return pion();}
    fWp=pp; fWpi=ppi; fWk=pk;
    if( !fSilent ){
cout<<" wpikp: mid,sig pi,k,p="<<qmpi<<" "<<sigpi<<";   "<<qmk<<" "<<sigk<<";   "
    <<qmp<<" "<<sigp<<"; "<<endl;
cout<<" aprob: "<<appi<<"  "<<apk<<"  "<<app<<endl;
cout<<" ppi,pk,pp="<<ppi<<"  "<<pk<<"  "<<pp<<endl;
     }
    if( ppi>pk&&ppi>pp )  { return pion(); }
    if(pk>pp){return kaon();}else{return proton();}
}
//-----------------------------------------------------------
Int_t	AliTPCPid::GetPcode(TClonesArray* /*rps*/,Float_t /*pm*/)
{
    return 0;    
}
//-----------------------------------------------------------
Int_t   AliTPCPid::GetPcode(AliTPCtrack *track)
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
    cout<<"TPCtrack dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
    return pcode?pcode:211;
    }
//-----------------------------------------------------------
Int_t   AliTPCPid::GetPcode(AliITStrackV2 *track)
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
    cout<<"lam,pt_1,mom,dedx="<<lam<<","<<pt_1<<","<<mom<<","<<dedx<<endl;
    Int_t pcode=GetPcode(dedx,mom);
    cout<<"ITS V2 dedx,mom,pcode="<<dedx<<","<<mom<<","<<pcode<<endl;
return pcode?pcode:211;
}
//-----------------------------------------------------------
Int_t	AliTPCPid::GetPcode(Float_t q,Float_t pm)
{
    fWpi=fWk=fWp=0.;     fPcode=0;
//1)---------------------- 0-120 MeV/c --------------
    if ( pm<=cut[1][0] )
	{ fWpi=1.; return pion(); }
//2)----------------------120-200 Mev/c ( 29.7.2002 ) ------------- 
    if ( pm<=cut[2][0] )
	{ if( q<fCutKa->Eval(pm) ){fWpi=1.; return pion(); } 
	                     else {fWk =1.; return  kaon();} }
//3)----------------------200-400 Mev/c -------------
    if ( pm<=cut[3][0] )
	if( q<fCutKa->Eval(pm) )
		{ fWpi=1.;return pion(); }
	     else
		{ if ( q<=fCutPr->Eval(pm) ) 
			     {fWk=1.;return kaon();} 
			else {fWp=1.;return proton();}
		}
//4)----------------------400-450 Mev/c -------------
    if ( pm<=cut[4][0] )
	if( q<fCutKaTune*fCutKa->Eval(pm) )
		{ fWpi=1.;return pion(); }
	    else
		{ if( q<fCutPr->Eval(pm) ) 
                  {fWk=1.;return kaon();} else {fWp=1.;return proton(); }
                }
//5)----------------------450-500 Mev/c -------------
    if ( pm<=cut[5][0] )
	if ( q>fCutPr->Eval(pm) )
           {fWp=1.;return proton();} else {return wpik(5,q);};
//6)----------------------500-550 Mev/c -------------
    if ( pm<=cut[6][0] )
	if ( q>fCutPr->Eval(pm) )
           {fWp=1.;return proton();} else {return wpik(6,q);};
//7)----------------------550-600 Mev/c -------------
    if ( pm<=cut[7][0] )
	if ( q>fCutPr->Eval(pm) )
           {fWp=1.;return proton();} else {return wpik(7,q);};
//8)----------------------600-650 Mev/c -------------
    if ( pm<=cut[8][0] )
      if ( q>fCutPrTune*fCutPr->Eval(pm) ){fWp=1.;return proton();} 
                                     else {return wpik(8,q);};
//9)----------------------650-730 Mev/c -------------
    if ( pm<=cut[9][0] )
      if ( q>fCutPrTune*fCutPr->Eval(pm) ){fWp=1.;return proton();}
                                     else {return wpik(9,q);};
//10)----------------------730-830 Mev/c -------------
    if( pm<=cut[10][0] )
      if ( q>fCutPrTune*fCutPr->Eval(pm) ){fWp=1.;return proton();}
                                     else {return wpik(10,q);};
//11)----------------------830-930 Mev/c -------------
    if( pm<=cut[11][0] ){ return wpikp(11,q); }
//12)----------------------930-1030 Mev/c -------------
    if( pm<=cut[12][0] )
     { return wpikp(12,q); };

    return fPcode;    
}
//-----------------------------------------------------------
void	AliTPCPid::SetCut(Int_t n,Float_t pm,Float_t pilo,Float_t pihi,
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
void AliTPCPid::SetVec(Int_t ntrack,TVector info)
{
TClonesArray& arr=*trs;
    new( arr[ntrack] ) TVector(info);
}
//-----------------------------------------------------------
TVector* AliTPCPid::GetVec(Int_t ntrack)
{
TClonesArray& arr=*trs;
    return (TVector*)arr[ntrack];
}
//-----------------------------------------------------------
void AliTPCPid::SetEdep(Int_t track,Float_t Edep)
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
void AliTPCPid::SetPmom(Int_t track,Float_t Pmom)
{
    TVector xx(0,11);
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector yy( *((TVector*)trs->At(track)) );xx=yy; }
    xx(10)=Pmom;
    TClonesArray &arr=*trs;
    new(arr[track])TVector(xx);
}
//-----------------------------------------------------------
void AliTPCPid::SetPcod(Int_t track,Int_t partcode)
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
void AliTPCPid::Print(Int_t track)
{cout<<mxtrs<<" tracks in AliITSPid obj."<<endl;
    if( ((TVector*)trs->At(track))->IsValid() )
	{TVector xx( *((TVector*)trs->At(track)) );
	 xx.Print();
	 }
    else 
	{cout<<"No data for track "<<track<<endl;return;}
}
//-----------------------------------------------------------
void AliTPCPid::Tab(void)
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
void AliTPCPid::Reset(void)
{
  for(Int_t i=0;i<trs->GetEntries();i++){
    TVector xx(0,11);
    TClonesArray &arr=*trs;
    new(arr[i])TVector(xx);
  }
}
//-----------------------------------------------------------
AliTPCPid::AliTPCPid(Int_t ntrack)
{
    trs = new TClonesArray("TVector",ntrack);
    TClonesArray &arr=*trs;
    for(Int_t i=0;i<ntrack;i++)new(arr[i])TVector(0,11);
    mxtrs=0;

    //fCutKa = new TF1("fkaons","[0]/x/x+[1]",0.1,1.2);
    //fCutPr = new TF1("fprotons","[0]/x/x +[1]",0.2,1.2);
    TF1 *f_rmska=0;
    
    f_rmska = new TF1("x_frmska","1.46-7.82*x+16.78*x^2-15.53*x^3+5.24*x^4 ",
	        0.1,1.2);
    fCutKa = new TF1("fkaons",
	   "1.25+0.044/x/x+1.25+0.044*x-13.87*x^2+22.37*x^3-10.05*x^4-2.5*x_frmska",
	   0.1,1.2);
    fCutPr = new TF1("fprotons",
                     "0.83*(15.32-56.094*x+89.962*x^2-66.1856*x^3+18.4052*x^4)",
		     0.2,1.2); 
    fCutKaTune=1.1; // 0.92; 
    fCutPrTune=1.0; //0.80;
    
const int inf=10;
//         Ncut Pmom   pilo  pihi    klo    khi     plo    phi
//       cut[j] [0]    [1]    [2]    [3]    [4]     [5]    [6]
//----------------------------------------------------------------
    SetCut(  1, 0.120,  0.  ,  0.  , inf  , inf   , inf  , inf  );
    SetCut(  2, 0.200,  0.  ,  6.0 , 6.0  , inf   , inf  , inf  ); //120-200
    SetCut(  3, 0.400,  0.  ,  3.5 , 3.5  , 9.0   , 9.0  , inf  ); //200-400
    SetCut(  4, 0.450,  0.  ,  1.9 , 1.9  , 4.0   , 4.0  , inf  ); //400-450
//----------------------------------------------------------------
    SetCut(  5, 0.500, 0.976, 0.108, 1.484 , 0.159  , 3.5  , inf  );  //450-500
    SetCut(  6, 0.550, 0.979, 0.108, 1.376 , 0.145  , 3.0  , inf  );  //500-550
//----------------------------------------------------------------    
    SetCut(  7, 0.600, 0.984, 0.111, 1.295 , 0.146 , 2.7  , inf  );   //550-600
    SetCut(  8, 0.650, 0.989, 0.113, 1.239 , 0.141 , 2.5  , inf  );   //600-650
    SetCut(  9, 0.730, 0.995, 0.109, 1.172 , 0.132 , 2.0  , inf  );   //650-730
//----------------------------------------------------------------    
    SetCut( 10, 0.830, 1.008, 0.116, 1.117 , 0.134 , 1.703, 0.209 ); //730-830
    SetCut( 11, 0.930, 1.019, 0.115, 1.072 , 0.121 , 1.535, 0.215 ); //830-930
    SetCut( 12, 1.230, 1.035, 0.117, 1.053 , 0.140  ,1.426, 0.270); //930-1030
    //------------------------ pi,K ---------------------
    aprob[0][0]=33219.;   aprob[1][0]=1971.;   // aprob[0][i] -    pions 
    aprob[0][1]=28828.;   aprob[1][1]=1973.; // aprob[1][i] -    kaons
    //---------------------------------------------------
    aprob[0][2]=24532.;   aprob[1][2]=1932.; aprob[2][2]=1948.;
    aprob[0][3]=20797.;   aprob[1][3]=1823.; aprob[2][3]=1970.;   
    aprob[0][4]=27017.;   aprob[1][4]=2681.; aprob[2][4]=2905.;
    //------------------------ pi,K,p -------------------
    aprob[0][5]= 24563.;    aprob[1][5]=2816.;  aprob[2][5]=3219.;  
    aprob[0][6]= 16877.;    aprob[1][6]=2231.;  aprob[2][6]=2746.;
    aprob[0][7]= 11557.;    aprob[1][7]=1681;   aprob[2][7]=2190.;

    fSilent=kTRUE;
}
//-----------------------------------------------------------



