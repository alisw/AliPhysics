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


// Implementation of the class to calculate the parton energy loss
// Based on the "BDMPS" quenching weights by C.A.Salgado and U.A.Wiedemann
// References:
// C.A.Salgado and U.A.Wiedemann, Phys.Rev.D68 (2003) 014008 [hep-ph/0302184]
// A.Dainese, Eur.Phys.J.C, in press, [nucl-ex/0312005]             
//
//
//            Origin:  C. Loizides   constantinos.loizides@cern.ch
//                     A. Dainese    andrea.dainese@pd.infn.it            
//
//=================== Added by C. Loizides 27/03/04 ===========================
//
// Added support for k-Quenching, where wc=I1*k and R=2I1^2/I0*k
// (see the AliFastGlauber class for definition of I0/I1)
//-----------------------------------------------------------------------------

#include <Riostream.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TLegend.h>
#include "AliQuenchingWeights.h"

ClassImp(AliQuenchingWeights)

// conversion from fm to GeV^-1: 1 fm = fmGeV GeV^-1
const Double_t AliQuenchingWeights::fgkConvFmToInvGeV = 1./0.197; 

// maximum value of R
const Double_t AliQuenchingWeights::fgkRMax = 1.e6; 

// hist binning
const Int_t AliQuenchingWeights::fgkBins = 1300; 
const Double_t AliQuenchingWeights::fgkMaxBin = 1.3; 

// counter for histogram labels
Int_t AliQuenchingWeights::fgCounter = 0; 


AliQuenchingWeights::AliQuenchingWeights() 
    : TObject(),
      fInstanceNumber(fgCounter++),
      fMultSoft(kTRUE), 
      fECMethod(kDefault),
      fQTransport(1.),
      fMu(1.),
      fK(4.e5),
      fLengthMax(20),
      fLengthMaxOld(0),
      fHistos(0),
      fHisto(0),
      fTablesLoaded(kFALSE)
{
  //default constructor 

}

AliQuenchingWeights::AliQuenchingWeights(const AliQuenchingWeights& a) 
    : TObject(),
      fInstanceNumber(fgCounter++),
      fMultSoft(kTRUE), 
      fECMethod(kDefault),
      fQTransport(1.),
      fMu(1.),
      fK(4.e5),
      fLengthMax(20),
      fLengthMaxOld(0),
      fHistos(0),
      fHisto(0),
      fTablesLoaded(kFALSE)
{
  // copy constructor 

  fTablesLoaded=kFALSE;
  fHistos=0;
  fLengthMaxOld=0;
  fMultSoft=a.GetMultSoft();; 
  fMu=a.GetMu();
  fK=a.GetK();
  fQTransport=a.GetQTransport();
  fECMethod=(kECMethod)a.GetECMethod();
  fLengthMax=a.GetLengthMax();
  fInstanceNumber=fgCounter++;
  Char_t name[100];
  snprintf(name,100, "hhistoqw_%d",fInstanceNumber);
  fHisto = new TH1F(name,"",fgkBins,0.,fgkMaxBin);
  for(Int_t bin=1;bin<=fgkBins;bin++) 
      fHisto->SetBinContent(bin,0.);

  //Missing in the class is the pathname
  //to the tables, can be added if needed
}

AliQuenchingWeights::~AliQuenchingWeights()
{
  Reset();
  delete fHisto;
}

void AliQuenchingWeights::Init()
{
//    Initialization
    if (fHisto) return;
    fHisto = new TH1F(Form("hhistoqw_%d",fInstanceNumber),"",fgkBins,0.,fgkMaxBin);
    for(Int_t bin=1;bin<=fgkBins;bin++) 
	fHisto->SetBinContent(bin,0.);
}

void AliQuenchingWeights::Reset()
{
  //reset tables if there were used

  if(!fHistos) return;
  for(Int_t l=0;l<4*fLengthMaxOld;l++){
    delete fHistos[0][l];
    delete fHistos[1][l];
  }
  delete[] fHistos;
  fHistos=0;
  fLengthMaxOld=0;
}

void AliQuenchingWeights::SetECMethod(kECMethod type)
{
  //set energy constraint method

  fECMethod=type;
  if(fECMethod==kDefault)
    Info("SetECMethod","Energy Constraint Method set to DEFAULT:\nIf (sampled energy loss > parton energy) then sampled energy loss = parton energy.");
  else if(fECMethod==kReweight)
    Info("SetECMethod","Energy Constraint Method set to REWEIGHT:\nRequire sampled energy loss <= parton energy.");
  else Info("SetECMethod","Energy Constraint Method set to REWEIGHTCONT:\nRequire sampled energy loss <= parton energy (only implemented for FAST method.");
}

Int_t AliQuenchingWeights::InitMult(const Char_t *contall,const Char_t *discall) 
{
  // read in tables for multiple scattering approximation
  // path to continuum and to discrete part

  fTablesLoaded = kFALSE;
  fMultSoft=kTRUE;
  
  Char_t fname[1024];
  snprintf(fname,1024, "%s",gSystem->ExpandPathName(contall));
  //PH  ifstream fincont(fname);
  fstream fincont(fname,ios::in);
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fincont.rdbuf()->is_open()) return -1;
#else
  if(!fincont.is_open()) return -1;
#endif

  Int_t nn=0; //quarks
  while(fincont>>fxx[nn]>>fcaq[0][nn]>>fcaq[1][nn]>>fcaq[2][nn]>>fcaq[3][nn]>>
	fcaq[4][nn]>>fcaq[5][nn]>>fcaq[6][nn]>>fcaq[7][nn]>>fcaq[8][nn]>>
	fcaq[9][nn]>>fcaq[10][nn]>>fcaq[11][nn]>>fcaq[12][nn]>>fcaq[13][nn]>>
	fcaq[14][nn]>>fcaq[15][nn]>>fcaq[16][nn]>>fcaq[17][nn]>>fcaq[18][nn]>>
	fcaq[19][nn]>>fcaq[20][nn]>>fcaq[21][nn]>>fcaq[22][nn]>>fcaq[23][nn]>>
	fcaq[24][nn]>>fcaq[25][nn]>>fcaq[26][nn]>>fcaq[27][nn]>>fcaq[28][nn]>>
	fcaq[29][nn]>>fcaq[30][nn]>>fcaq[31][nn]>>fcaq[32][nn]>>fcaq[33][nn]) 
    {
      nn++;
      if(nn==261) break;
    }

  nn=0;       //gluons
  while(fincont>>fxxg[nn]>>fcag[0][nn]>>fcag[1][nn]>>fcag[2][nn]>>fcag[3][nn]>>
	fcag[4][nn]>>fcag[5][nn]>>fcag[6][nn]>>fcag[7][nn]>>fcag[8][nn]>>
	fcag[9][nn]>>fcag[10][nn]>>fcag[11][nn]>>fcag[12][nn]>>fcag[13][nn]>>
        fcag[14][nn]>>fcag[15][nn]>>fcag[16][nn]>>fcag[17][nn]>>fcag[18][nn]>>
	fcag[19][nn]>>fcag[20][nn]>>fcag[21][nn]>>fcag[22][nn]>>fcag[23][nn]>>
	fcag[24][nn]>>fcag[25][nn]>>fcag[26][nn]>>fcag[27][nn]>>fcag[28][nn]>>
	fcag[29][nn]>>fcag[30][nn]>>fcag[31][nn]>>fcag[32][nn]>>fcag[33][nn]) 
  {
    nn++;
    if(nn==261) break;
  }
  fincont.close();

  snprintf(fname,1024, "%s",gSystem->ExpandPathName(discall));
  //PH  ifstream findisc(fname); 
  fstream findisc(fname,ios::in); 
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!findisc.rdbuf()->is_open()) return -1;
#else
  if(!findisc.is_open()) return -1;
#endif

  nn=0; //quarks
  while(findisc>>frrr[nn]>>fdaq[nn]) {
    nn++;
    if(nn==34) break;
  }
  nn=0; //gluons
  while(findisc>>frrrg[nn]>>fdag[nn]) {
    nn++;
    if(nn==34) break;
  }
  findisc.close();
  fTablesLoaded = kTRUE;
  return 0;
}

/*
C***************************************************************************
C	Quenching Weights for Multiple Soft Scattering
C		February 10, 2003
C
C	Refs:
C
C  Carlos A. Salgado and Urs A. Wiedemann, hep-ph/0302184.                 
C
C  Carlos A. Salgado and Urs A. Wiedemann Phys.Rev.Lett.89:092303,2002.
C 
C
C   This package contains quenching weights for gluon radiation in the
C   multiple soft scattering approximation.
C
C   swqmult returns the quenching weight for a quark (ipart=1) or 
C   a gluon (ipart=2) traversing a medium with transport coeficient q and
C   length L. The input values are rrrr=0.5*q*L^3 and xxxx=w/wc, where
C   wc=0.5*q*L^2 and w is the energy radiated. The output values are
C   the continuous and discrete (prefactor of the delta function) parts
C   of the quenching weights.
C	
C   In order to use this routine, the files cont.all and disc.all need to be
C   in the working directory. 
C
C   An initialization of the tables is needed by doing call initmult before
C   using swqmult.
C
C   Please, send us any comment:
C
C	urs.wiedemann@cern.ch
C	carlos.salgado@cern.ch
C
C
C-------------------------------------------------------------------

      SUBROUTINE swqmult(ipart,rrrr,xxxx,continuous,discrete)
*
      REAL*8           xx(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqua/    xx, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglu/    xxg, dag, cag, rrrg

      REAL*8           rrrr,xxxx, continuous, discrete
      REAL*8           rrin, xxin
      INTEGER          nrlow, nrhigh, nxlow, nxhigh
      REAL*8           rrhigh, rrlow, rfraclow, rfrachigh
      REAL*8           xfraclow, xfrachigh
      REAL*8           clow, chigh
*

      continuous=0.d0
      discrete=0.d0

      rrin = rrrr
      xxin = xxxx
*
      do 666, nr=1,34
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow = rrr(nr)
            nrlow = nr
            nrhigh = nr-1
            goto 665
         endif
 666     enddo
 665     continue
*
      rfraclow = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
      if (rrin.gt.10000d0) then
         rfraclow = dlog(rrhigh/rrin)/dlog(rrhigh/rrlow)
         rfrachigh = dlog(rrin/rrlow)/dlog(rrhigh/rrlow)
      endif
*
      if (ipart.eq.1.and.rrin.ge.rrr(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (ipart.ne.1.and.rrin.ge.rrrg(1)) then
         nrlow=1
         nrhigh=1
         rfraclow=1
         rfrachigh=0
      endif

      if (xxxx.ge.xx(261)) go to 245

      nxlow = int(xxin/0.01) + 1
      nxhigh = nxlow + 1
      xfraclow = (xx(nxhigh)-xxin)/0.01
      xfrachigh = (xxin - xx(nxlow))/0.01
*
      if(ipart.eq.1) then
      clow = xfraclow*caq(nrlow,nxlow)+xfrachigh*caq(nrlow,nxhigh)
      chigh = xfraclow*caq(nrhigh,nxlow)+xfrachigh*caq(nrhigh,nxhigh)
      else
      clow = xfraclow*cag(nrlow,nxlow)+xfrachigh*cag(nrlow,nxhigh)
      chigh = xfraclow*cag(nrhigh,nxlow)+xfrachigh*cag(nrhigh,nxhigh)
      endif

      continuous = rfraclow*clow + rfrachigh*chigh

245   continue

      if(ipart.eq.1) then
      discrete = rfraclow*daq(nrlow) + rfrachigh*daq(nrhigh)
      else
      discrete = rfraclow*dag(nrlow) + rfrachigh*dag(nrhigh)
      endif
*
      END

      subroutine initmult
      REAL*8           xxq(400), daq(34), caq(34,261), rrr(34)
      COMMON /dataqua/    xxq, daq, caq, rrr
*
      REAL*8           xxg(400), dag(34), cag(34,261), rrrg(34)
      COMMON /dataglu/    xxg, dag, cag, rrrg
*
      OPEN(UNIT=20,FILE='contnew.all',STATUS='OLD',ERR=90)
      do 110 nn=1,261
      read (20,*) xxq(nn), caq(1,nn), caq(2,nn), caq(3,nn),
     +     caq(4,nn), caq(5,nn), caq(6,nn), caq(7,nn), caq(8,nn),
     +     caq(9,nn), caq(10,nn), caq(11,nn), caq(12,nn), 
     +     caq(13,nn),
     +     caq(14,nn), caq(15,nn), caq(16,nn), caq(17,nn), 
     +     caq(18,nn),
     +     caq(19,nn), caq(20,nn), caq(21,nn), caq(22,nn), 
     +     caq(23,nn),
     +     caq(24,nn), caq(25,nn), caq(26,nn), caq(27,nn), 
     +     caq(28,nn),
     +     caq(29,nn), caq(30,nn), caq(31,nn), caq(32,nn), 
     +     caq(33,nn), caq(34,nn)
 110     continue
      do 111 nn=1,261
      read (20,*) xxg(nn), cag(1,nn), cag(2,nn), cag(3,nn),
     +     cag(4,nn), cag(5,nn), cag(6,nn), cag(7,nn), cag(8,nn),
     +     cag(9,nn), cag(10,nn), cag(11,nn), cag(12,nn), 
     +     cag(13,nn),
     +     cag(14,nn), cag(15,nn), cag(16,nn), cag(17,nn), 
     +     cag(18,nn),
     +     cag(19,nn), cag(20,nn), cag(21,nn), cag(22,nn), 
     +     cag(23,nn),
     +     cag(24,nn), cag(25,nn), cag(26,nn), cag(27,nn), 
     +     cag(28,nn),
     +     cag(29,nn), cag(30,nn), cag(31,nn), cag(32,nn), 
     +     cag(33,nn), cag(34,nn)
 111     continue
      close(20)
*
      OPEN(UNIT=21,FILE='discnew.all',STATUS='OLD',ERR=91)
      do 112 nn=1,34
      read (21,*) rrr(nn), daq(nn)
 112     continue
      do 113 nn=1,34
      read (21,*) rrrg(nn), dag(nn)
 113     continue
      close(21)
*
      goto 888
 90   PRINT*, 'input - output error' 
 91   PRINT*, 'input - output error #2' 
 888  continue

      end


=======================================================================

   Adapted to ROOT macro by A. Dainese - 13/07/2003
   Ported to class by C. Loizides - 12/02/2004
   New version for extended R values added - 06/03/2004 
*/

Int_t AliQuenchingWeights::CalcMult(Int_t ipart, Double_t rrrr,Double_t xxxx,
                              Double_t &continuous,Double_t &discrete) const
{
  // Calculate Multiple Scattering approx. 
  // weights for given parton type, 
  // rrrr=0.5*q*L^3 and xxxx=w/wc, wc=0.5*q*L^2

  //set result to zero
  continuous=0.;
  discrete=0.;

  //read-in data before first call
  if(!fTablesLoaded){
    Error("CalcMult","Tables are not loaded.");
    return -1;
  }
  if(!fMultSoft){
    Error("CalcMult","Tables are not loaded for Multiple Scattering.");
    return -1;
  }

  Double_t rrin = rrrr;
  Double_t xxin = xxxx;

  if(xxin>fxx[260]) return -1;
  Int_t nxlow     = (Int_t)(xxin/0.01) + 1;
  Int_t nxhigh    = nxlow + 1;
  Double_t xfraclow  = (fxx[nxhigh-1]-xxin)/0.01;
  Double_t xfrachigh = (xxin - fxx[nxlow-1])/0.01;

  //why this?
  if(rrin<=frrr[33]) rrin = 1.05*frrr[33]; // AD
  if(rrin>=frrr[0])  rrin = 0.95*frrr[0];  // AD

  Int_t nrlow=0,nrhigh=0;
  Double_t rrhigh=0,rrlow=0;
  for(Int_t nr=1; nr<=34; nr++) {
    if(rrin<frrr[nr-1]) {
      rrhigh = frrr[nr-1];
    } else {
      rrhigh = frrr[nr-1-1];
      rrlow  = frrr[nr-1];
      nrlow  = nr;
      nrhigh = nr-1;
      break;
    }
  }

  rrin = rrrr; // AD

  Double_t rfraclow  = (rrhigh-rrin)/(rrhigh-rrlow);
  Double_t rfrachigh = (rrin-rrlow)/(rrhigh-rrlow);

  if(rrin>1.e4){
    rfraclow = TMath::Log2(rrhigh/rrin)/TMath::Log2(rrhigh/rrlow);    
    rfrachigh = TMath::Log2(rrin/rrlow)/TMath::Log2(rrhigh/rrlow);
  }
  if((ipart==1) && (rrin>=frrr[0]))
  {
    nrlow=1;
    nrhigh=1;
    rfraclow=1.;
    rfrachigh=0.;
  }
  if((ipart==2) && (rrin>=frrrg[0]))
  {
    nrlow=1;
    nrhigh=1;
    rfraclow=1.;
    rfrachigh=0.;
  }

  //printf("R = %f,\nRlow = %f, Rhigh = %f,\nRfraclow = %f, Rfrachigh = %f\n",rrin,rrlow,rrhigh,rfraclow,rfrachigh); // AD

  Double_t clow=0,chigh=0;
  if(ipart==1) {
    clow  = xfraclow*fcaq[nrlow-1][nxlow-1]+xfrachigh*fcaq[nrlow-1][nxhigh-1];
    chigh = xfraclow*fcaq[nrhigh-1][nxlow-1]+xfrachigh*fcaq[nrhigh-1][nxhigh-1];
  } else {
    clow  = xfraclow*fcag[nrlow-1][nxlow-1]+xfrachigh*fcag[nrlow-1][nxhigh-1];
    chigh = xfraclow*fcag[nrhigh-1][nxlow-1]+xfrachigh*fcag[nrhigh-1][nxhigh-1];
  }

  continuous = rfraclow*clow + rfrachigh*chigh;
  //printf("rfraclow %f, clow %f, rfrachigh %f, chigh %f,\n continuous %f\n",
  //rfraclow,clow,rfrachigh,chigh,continuous);

  if(ipart==1) {
    discrete = rfraclow*fdaq[nrlow-1] + rfrachigh*fdaq[nrhigh-1];
  } else {
    discrete = rfraclow*fdag[nrlow-1] + rfrachigh*fdag[nrhigh-1];
  }

  return 0;
}

Int_t AliQuenchingWeights::InitSingleHard(const Char_t *contall,const Char_t *discall) 
{
  // read in tables for Single Hard Approx.
  // path to continuum and to discrete part

  fTablesLoaded = kFALSE;
  fMultSoft=kFALSE;
  
  Char_t fname[1024];
  snprintf(fname, 1024, "%s",gSystem->ExpandPathName(contall));
  //PH  ifstream fincont(fname);
  fstream fincont(fname,ios::in);
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fincont.rdbuf()->is_open()) return -1;
#else
  if(!fincont.is_open()) return -1;
#endif

  Int_t nn=0; //quarks
  while(fincont>>fxx[nn]>>fcaq[0][nn]>>fcaq[1][nn]>>fcaq[2][nn]>>fcaq[3][nn]>>
	fcaq[4][nn]>>fcaq[5][nn]>>fcaq[6][nn]>>fcaq[7][nn]>>fcaq[8][nn]>>
	fcaq[9][nn]>>fcaq[10][nn]>>fcaq[11][nn]>>fcaq[12][nn]>>
	fcaq[13][nn]>>
	fcaq[14][nn]>>fcaq[15][nn]>>fcaq[16][nn]>>fcaq[17][nn]>>
	fcaq[18][nn]>>
	fcaq[19][nn]>>fcaq[20][nn]>>fcaq[21][nn]>>fcaq[22][nn]>>
	fcaq[23][nn]>>
	fcaq[24][nn]>>fcaq[25][nn]>>fcaq[26][nn]>>fcaq[27][nn]>>
	fcaq[28][nn]>>
	fcaq[29][nn]) 
    {
      nn++;
      if(nn==261) break;
    }

  nn=0;       //gluons
  while(fincont>>fxxg[nn]>>fcag[0][nn]>>fcag[1][nn]>>fcag[2][nn]>>fcag[3][nn]>>
	fcag[4][nn]>>fcag[5][nn]>>fcag[6][nn]>>fcag[7][nn]>>fcag[8][nn]>>
	fcag[9][nn]>>fcag[10][nn]>>fcag[11][nn]>>fcag[12][nn]>>
	fcag[13][nn]>>
	fcag[14][nn]>>fcag[15][nn]>>fcag[16][nn]>>fcag[17][nn]>>
	fcag[18][nn]>>
	fcag[19][nn]>>fcag[20][nn]>>fcag[21][nn]>>fcag[22][nn]>>
	fcag[23][nn]>>
	fcag[24][nn]>>fcag[25][nn]>>fcag[26][nn]>>fcag[27][nn]>>
	fcag[28][nn]>>
	fcag[29][nn]) {
    nn++;
    if(nn==261) break;
  }
  fincont.close();

  snprintf(fname, 1024, "%s",gSystem->ExpandPathName(discall));
  //PH  ifstream findisc(fname); 
  fstream findisc(fname,ios::in); 
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!findisc.rdbuf()->is_open()) return -1;
#else
  if(!findisc.is_open()) return -1;
#endif

  nn=0; //quarks
  while(findisc>>frrr[nn]>>fdaq[nn]) {
    nn++;
    if(nn==30) break;
  }
  nn=0; //gluons
  while(findisc>>frrrg[nn]>>fdag[nn]) {
    nn++;
    if(nn==30) break;
  }
  findisc.close();

  fTablesLoaded = kTRUE;
  return 0;
}

/*
C***************************************************************************
C       Quenching Weights for Single Hard Scattering
C               February 20, 2003
C
C       Refs:
C
C  Carlos A. Salgado and Urs A. Wiedemann, hep-ph/0302184. 
C
C  Carlos A. Salgado and Urs A. Wiedemann Phys.Rev.Lett.89:092303,2002.
C  
C
C   This package contains quenching weights for gluon radiation in the
C   single hard scattering approximation.
C
C   swqlin returns the quenching weight for a quark (ipart=1) or
C   a gluon (ipart=2) traversing a medium with Debye screening mass mu and
C   length L. The input values are rrrr=0.5*mu^2*L^2 and xxxx=w/wc, where
C   wc=0.5*mu^2*L and w is the energy radiated. The output values are
C   the continuous and discrete (prefactor of the delta function) parts
C   of the quenching weights.
C
C   In order to use this routine, the files contlin.all and disclin.all 
C   need to be in the working directory.
C
C   An initialization of the tables is needed by doing call initlin before
C   using swqlin.
C
C   Please, send us any comment:
C
C       urs.wiedemann@cern.ch
C       carlos.salgado@cern.ch
C
C
C-------------------------------------------------------------------


      SUBROUTINE swqlin(ipart,rrrr,xxxx,continuous,discrete)
*
      REAL*8           xx(400), dalq(30), calq(30,261), rrr(30)
      COMMON /datalinqua/    xx, dalq, calq, rrr
*
      REAL*8           xxlg(400), dalg(30), calg(30,261), rrrlg(30)
      COMMON /datalinglu/    xxlg, dalg, calg, rrrlg

      REAL*8           rrrr,xxxx, continuous, discrete
      REAL*8           rrin, xxin
      INTEGER          nrlow, nrhigh, nxlow, nxhigh
      REAL*8           rrhigh, rrlow, rfraclow, rfrachigh
      REAL*8           xfraclow, xfrachigh
      REAL*8           clow, chigh
*
      rrin = rrrr
      xxin = xxxx
*
      nxlow = int(xxin/0.038) + 1
      nxhigh = nxlow + 1
      xfraclow = (xx(nxhigh)-xxin)/0.038
      xfrachigh = (xxin - xx(nxlow))/0.038
*
      do 666, nr=1,30
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow = rrr(nr)
            nrlow = nr
            nrhigh = nr-1
            goto 665
         endif
 666     enddo
 665     continue
*
      rfraclow = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
*
      if(ipart.eq.1) then
      clow = xfraclow*calq(nrlow,nxlow)+xfrachigh*calq(nrlow,nxhigh)
      chigh = xfraclow*calq(nrhigh,nxlow)+xfrachigh*calq(nrhigh,nxhigh)
      else
      clow = xfraclow*calg(nrlow,nxlow)+xfrachigh*calg(nrlow,nxhigh)
      chigh = xfraclow*calg(nrhigh,nxlow)+xfrachigh*calg(nrhigh,nxhigh)
      endif

      continuous = rfraclow*clow + rfrachigh*chigh

      if(ipart.eq.1) then
      discrete = rfraclow*dalq(nrlow) + rfrachigh*dalq(nrhigh)
      else
      discrete = rfraclow*dalg(nrlow) + rfrachigh*dalg(nrhigh)
      endif
*
      END

      subroutine initlin
      REAL*8           xxlq(400), dalq(30), calq(30,261), rrr(30)
      COMMON /datalinqua/    xxlq, dalq, calq, rrr
*
      REAL*8           xxlg(400), dalg(30), calg(30,261), rrrlg(30)
      COMMON /datalinglu/    xxlg, dalg, calg, rrrlg
*
      OPEN(UNIT=20,FILE='contlin.all',STATUS='OLD',ERR=90)
      do 110 nn=1,261
      read (20,*) xxlq(nn), calq(1,nn), calq(2,nn), calq(3,nn),
     +     calq(4,nn), calq(5,nn), calq(6,nn), calq(7,nn), calq(8,nn),
     +     calq(9,nn), calq(10,nn), calq(11,nn), calq(12,nn), 
     +     calq(13,nn),
     +     calq(14,nn), calq(15,nn), calq(16,nn), calq(17,nn), 
     +     calq(18,nn),
     +     calq(19,nn), calq(20,nn), calq(21,nn), calq(22,nn), 
     +     calq(23,nn),
     +     calq(24,nn), calq(25,nn), calq(26,nn), calq(27,nn), 
     +     calq(28,nn),
     +     calq(29,nn), calq(30,nn)
 110     continue
      do 111 nn=1,261
      read (20,*) xxlg(nn), calg(1,nn), calg(2,nn), calg(3,nn),
     +     calg(4,nn), calg(5,nn), calg(6,nn), calg(7,nn), calg(8,nn),
     +     calg(9,nn), calg(10,nn), calg(11,nn), calg(12,nn), 
     +     calg(13,nn),
     +     calg(14,nn), calg(15,nn), calg(16,nn), calg(17,nn), 
     +     calg(18,nn),
     +     calg(19,nn), calg(20,nn), calg(21,nn), calg(22,nn), 
     +     calg(23,nn),
     +     calg(24,nn), calg(25,nn), calg(26,nn), calg(27,nn), 
     +     calg(28,nn),
     +     calg(29,nn), calg(30,nn)
 111     continue
      close(20)
*
      OPEN(UNIT=21,FILE='disclin.all',STATUS='OLD',ERR=91)
      do 112 nn=1,30
      read (21,*) rrr(nn), dalq(nn)
 112     continue
      do 113 nn=1,30
      read (21,*) rrrlg(nn), dalg(nn)
 113     continue
      close(21)
*
      goto 888
 90   PRINT*, 'input - output error' 
 91   PRINT*, 'input - output error #2' 
 888  continue

      end

=======================================================================

   Ported to class by C. Loizides - 17/02/2004

*/

Int_t AliQuenchingWeights::CalcSingleHard(Int_t ipart, Double_t rrrr,Double_t xxxx,
                                    Double_t &continuous,Double_t &discrete) const
{
  // calculate Single Hard approx. 
  // weights for given parton type, 
  // rrrr=0.5*mu^2*L^2 and xxxx=w/wc, wc=0.5*mu^2*L

  // read-in data before first call
  if(!fTablesLoaded){
    Error("CalcSingleHard","Tables are not loaded.");
    return -1;
  }
  if(fMultSoft){
    Error("CalcSingleHard","Tables are not loaded for Single Hard Scattering.");
    return -1;
  }

  Double_t rrin = rrrr;
  Double_t xxin = xxxx;

  Int_t nxlow     = (Int_t)(xxin/0.038) + 1;
  Int_t nxhigh    = nxlow + 1;
  Double_t xfraclow  = (fxx[nxhigh-1]-xxin)/0.038;
  Double_t xfrachigh = (xxin - fxx[nxlow-1])/0.038;

  //why this?
  if(rrin<=frrr[29]) rrin = 1.05*frrr[29]; // AD
  if(rrin>=frrr[0])  rrin = 0.95*frrr[0];  // AD

  Int_t nrlow=0,nrhigh=0;
  Double_t rrhigh=0,rrlow=0;
  for(Int_t nr=1; nr<=30; nr++) {
    if(rrin<frrr[nr-1]) {
      rrhigh = frrr[nr-1];
    } else {
      rrhigh = frrr[nr-1-1];
      rrlow  = frrr[nr-1];
      nrlow  = nr;
      nrhigh = nr-1;
      break;
    }
  }

  rrin = rrrr; // AD

  Double_t rfraclow  = (rrhigh-rrin)/(rrhigh-rrlow);
  Double_t rfrachigh = (rrin-rrlow)/(rrhigh-rrlow);

  //printf("R = %f,\nRlow = %f, Rhigh = %f,\nRfraclow = %f, Rfrachigh = %f\n",rrin,rrlow,rrhigh,rfraclow,rfrachigh); // AD

  Double_t clow=0,chigh=0;
  if(ipart==1) {
    clow  = xfraclow*fcaq[nrlow-1][nxlow-1]+xfrachigh*fcaq[nrlow-1][nxhigh-1];
    chigh = xfraclow*fcaq[nrhigh-1][nxlow-1]+xfrachigh*fcaq[nrhigh-1][nxhigh-1];
  } else {
    clow  = xfraclow*fcag[nrlow-1][nxlow-1]+xfrachigh*fcag[nrlow-1][nxhigh-1];
    chigh = xfraclow*fcag[nrhigh-1][nxlow-1]+xfrachigh*fcag[nrhigh-1][nxhigh-1];
  }

  continuous = rfraclow*clow + rfrachigh*chigh;
  //printf("rfraclow %f, clow %f, rfrachigh %f, chigh %f,\n continuous %f\n",
  //	 rfraclow,clow,rfrachigh,chigh,continuous);

  if(ipart==1) {
    discrete = rfraclow*fdaq[nrlow-1] + rfrachigh*fdaq[nrhigh-1];
  } else {
    discrete = rfraclow*fdag[nrlow-1] + rfrachigh*fdag[nrhigh-1];
  }

  return 0;
}

Int_t AliQuenchingWeights::CalcMult(Int_t ipart, 
                              Double_t w,Double_t qtransp,Double_t length,
                              Double_t &continuous,Double_t &discrete) const
{
  // Calculate Multiple Scattering approx. 
  // weights for given parton type, 
  // rrrr=0.5*q*L^3 and xxxx=w/wc, wc=0.5*q*L^2

  Double_t wc=CalcWC(qtransp,length);
  Double_t rrrr=CalcR(wc,length);
  Double_t xxxx=w/wc;
  return CalcMult(ipart,rrrr,xxxx,continuous,discrete);
}

Int_t AliQuenchingWeights::CalcSingleHard(Int_t ipart, 
                                    Double_t w,Double_t mu,Double_t length,
                                    Double_t &continuous,Double_t &discrete) const
{
  // calculate Single Hard approx. 
  // weights for given parton type, 
  // rrrr=0.5*mu^2*L^2 and xxxx=w/wc, wc=0.5*mu^2*L

  Double_t wcbar=CalcWCbar(mu,length);
  Double_t rrrr=CalcR(wcbar,length);
  Double_t xxxx=w/wcbar;
  return CalcSingleHard(ipart,rrrr,xxxx,continuous,discrete);
}

Double_t AliQuenchingWeights::CalcR(Double_t wc, Double_t l) const 
{ 
  //calculate r value and 
  //check if it is less then maximum

  Double_t r = wc*l*fgkConvFmToInvGeV;
  if(r >= fgkRMax) {
    Warning("CalcR","Value of r = %.2f; should be less than %.2f", r, fgkRMax);
    return fgkRMax-1;
  }  
  return r;
}

Double_t AliQuenchingWeights::CalcRk(Double_t k, Double_t I0, Double_t I1) const
{ 
  //calculate R value and 
  //check if it is less then maximum

  Double_t r = fgkRMax-1;
  if(I0>0)
    r = 2*I1*I1/I0*k;
  if(r>=fgkRMax) {
    Warning("CalcRk","Value of r = %.2f; should be less than %.2f",r,fgkRMax);
    return fgkRMax-1;
  }  
  return r;
}

Double_t AliQuenchingWeights::GetELossRandom(Int_t ipart, Double_t length, Double_t e) const
{
  // return DeltaE for MS or SH scattering
  // for given parton type, length and energy
  // Dependant on ECM (energy constraint method)
  // e is used to determine where to set bins to zero.

  if(!fHistos){
    Fatal("GetELossRandom","Call SampleEnergyLoss method before!");
    return -1000.;
  }
  if((ipart<1) || (ipart>2)) {
    Fatal("GetELossRandom","ipart =%d; but has to be 1 (quark) or 2 (gluon)",ipart);
    return -1000.;
  }

  Int_t l=GetIndex(length);
  if(l<=0) return 0.;

  if(fECMethod==kReweight){
    Double_t ret = 2.*e;
    Int_t ws=0;
    while(ret>e){
      ret=fHistos[ipart-1][l-1]->GetRandom(); 
      if(++ws==1e6){
	Warning("GetELossRandom",
                "Stopping reweighting; maximum loss assigned after 1e6 trials.");
	return e;
      }
    }
    return ret;
  }
  //kDefault
  Double_t ret=fHistos[ipart-1][l-1]->GetRandom();
  if(ret>e) return e;
  return ret;
}

Double_t AliQuenchingWeights::CalcQuenchedEnergy(Int_t ipart, Double_t length, Double_t e) const
{
  //return quenched parton energy
  //for given parton type, length and energy

  Double_t loss=GetELossRandom(ipart,length,e);
  return e-loss;
}

Double_t AliQuenchingWeights::GetELossRandom(Int_t ipart, TH1F *hell, Double_t e) const
{
  // return DeltaE for MS or SH scattering
  // for given parton type, length distribution and energy
  // Dependant on ECM (energy constraint method)
  // e is used to determine where to set bins to zero.

  if(!hell){
    Warning("GetELossRandom","Pointer to length distribution is NULL.");
    return 0.;
  }
  Double_t ell=hell->GetRandom();
  return GetELossRandom(ipart,ell,e);
}

Double_t AliQuenchingWeights::CalcQuenchedEnergy(Int_t ipart, TH1F *hell, Double_t e)  const
{
  //return quenched parton energy
  //for given parton type, length distribution and energy

  Double_t loss=GetELossRandom(ipart,hell,e);
  return e-loss;
}

Double_t AliQuenchingWeights::GetELossRandomK(Int_t ipart, Double_t I0, Double_t I1, Double_t e)
{
  // return DeltaE for new dynamic version
  // for given parton type, I0 and I1 value and energy
  // Dependant on ECM (energy constraint method)
  // e is used to determine where to set bins to zero.

  // read-in data before first call
  if(!fTablesLoaded){
    Fatal("GetELossRandomK","Tables are not loaded.");
    return -1000.;
  }
  if((ipart<1) || (ipart>2)) {
    Fatal("GetELossRandomK","ipart =%d; but has to be 1 (quark) or 2 (gluon)",ipart);
    return -1000.;
  }

  Double_t r=CalcRk(I0,I1);
  if(r<0.){
    Fatal("GetELossRandomK","R should not be negative");
    return -1000.;
  }
  Double_t wc=CalcWCk(I1);
  if(wc<=0.){
    Fatal("GetELossRandomK","wc should be greater than zero");
    return -1000.;
  }
  if(SampleEnergyLoss(ipart,r)!=0){
    Fatal("GetELossRandomK","Could not sample energy loss");
    return -1000.;
  }

  if(fECMethod==kReweight){
    Double_t ret = 2.*e;
    Int_t ws=0;
    while(ret>e){
      ret=fHisto->GetRandom(); 
      if(++ws==1e6){
	Warning("GetELossRandomK",
                "Stopping reweighting; maximum loss assigned after 1e6 trials.");
	return e;
      }
    }
    return ret;
  }

  //kDefault
  Double_t ret=fHisto->GetRandom()*wc;
  if(ret>e) return e;
  return ret;
}

Double_t AliQuenchingWeights::CalcQuenchedEnergyK(Int_t ipart, Double_t I0, Double_t I1, Double_t e)
{
  //return quenched parton energy
  //for given parton type, I0 and I1 value and energy

  Double_t loss=GetELossRandomK(ipart,I0,I1,e);
  return e-loss;
}

Double_t AliQuenchingWeights::GetELossRandomKFast(Int_t ipart, Double_t I0, Double_t I1, Double_t e)
{
  // return DeltaE for new dynamic version
  // for given parton type, I0 and I1 value and energy
  // Dependant on ECM (energy constraint method)
  // e is used to determine where to set bins to zero.
  // method is optimized and should only be used if 
  // all parameters are well within the bounds.
  // read-in data tables before first call 

  Double_t r=CalcRk(I0,I1);
  if(r<=0.){
    return 0.;
  }

  Double_t wc=CalcWCk(I1);
  if(wc<=0.){
    return 0.;
  }

  return GetELossRandomKFastR(ipart,r,wc,e);
}

Double_t AliQuenchingWeights::GetELossRandomKFastR(Int_t ipart, Double_t r, Double_t wc, Double_t e)
{
  // return DeltaE for new dynamic version
  // for given parton type, R and wc value and energy
  // Dependant on ECM (energy constraint method)
  // e is used to determine where to set bins to zero.
  // method is optimized and should only be used if 
  // all parameters are well within the bounds.
  // read-in data tables before first call 

  if(r>=fgkRMax) {
    r=fgkRMax-1;
  }  
  
  Double_t discrete=0.;
  Double_t continuous=0.;
  Int_t bin=1;
  Double_t xxxx = fHisto->GetBinCenter(bin);
  if(fMultSoft)
    CalcMult(ipart,r,xxxx,continuous,discrete);
  else
    CalcSingleHard(ipart,r,xxxx,continuous,discrete);

  if(discrete>=1.0) {
    return 0.; //no energy loss
  }
  if (!fHisto) Init();
  
  fHisto->SetBinContent(bin,continuous);
  Int_t kbinmax=fHisto->FindBin(e/wc);
  if(kbinmax>=fgkBins) kbinmax=fgkBins-1;
  if(kbinmax==1) return e; //maximum energy loss

  if(fMultSoft) {
    for(bin=2; bin<=kbinmax; bin++) {
      xxxx = fHisto->GetBinCenter(bin);
      CalcMult(ipart,r,xxxx,continuous,discrete);
      fHisto->SetBinContent(bin,continuous);
    }
  } else {
    for(bin=2; bin<=kbinmax; bin++) {
      xxxx = fHisto->GetBinCenter(bin);
      CalcSingleHard(ipart,r,xxxx,continuous,discrete);
      fHisto->SetBinContent(bin,continuous);
    }
  }

  if(fECMethod==kReweight){
    fHisto->SetBinContent(kbinmax+1,0);
    fHisto->Fill(0.,discrete*fgkBins/fgkMaxBin);
  } else if (fECMethod==kReweightCont) {
    fHisto->SetBinContent(kbinmax+1,0);
    const Double_t kdelta=fHisto->Integral(1,kbinmax);
    fHisto->Scale(1./kdelta*(1-discrete));
    fHisto->Fill(0.,discrete);
  } else {
    const Double_t kdelta=fHisto->Integral(1,kbinmax);
    Double_t val=discrete*fgkBins/fgkMaxBin;
    fHisto->Fill(0.,val);
    fHisto->SetBinContent(kbinmax+1,(1-discrete)*fgkBins/fgkMaxBin-kdelta);
  }
  for(bin=kbinmax+2; bin<=fgkBins; bin++) {
    fHisto->SetBinContent(bin,0);
  }
  //cout << kbinmax << " " << discrete << " " << fHisto->Integral() << endl;
  Double_t ret=fHisto->GetRandom()*wc;
  if(ret>e) return e;
  return ret;
}

Double_t AliQuenchingWeights::CalcQuenchedEnergyKFast(Int_t ipart, Double_t I0, Double_t I1, Double_t e)
{
  //return quenched parton energy (fast method)
  //for given parton type, I0 and I1 value and energy

  Double_t loss=GetELossRandomKFast(ipart,I0,I1,e);
  return e-loss;
}

Double_t AliQuenchingWeights::GetDiscreteWeight(Int_t ipart, Double_t I0, Double_t I1)
{
  // return discrete weight

  Double_t r=CalcRk(I0,I1);
  if(r<=0.){
    return 1.;
  }
  return GetDiscreteWeightR(ipart,r);
}

Double_t AliQuenchingWeights::GetDiscreteWeightR(Int_t ipart, Double_t r)
{
  // return discrete weight

  if(r>=fgkRMax) {
    r=fgkRMax-1;
  }  

  Double_t discrete=0.;
  Double_t continuous=0.;
  Int_t bin=1;
  Double_t xxxx = fHisto->GetBinCenter(bin);
  if(fMultSoft)
    CalcMult(ipart,r,xxxx,continuous,discrete);
  else
    CalcSingleHard(ipart,r,xxxx,continuous,discrete);
  return discrete;
}

void AliQuenchingWeights::GetZeroLossProb(Double_t &p,Double_t &prw,Double_t &prwcont,
					  Int_t ipart,Double_t I0,Double_t I1,Double_t e)
{
  //calculate the probabilty that there is no energy
  //loss for different ways of energy constraint
  p=1.;prw=1.;prwcont=1.;
  Double_t r=CalcRk(I0,I1);
  if(r<=0.){
    return;
  }
  Double_t wc=CalcWCk(I1);
  if(wc<=0.){
    return;
  }
  GetZeroLossProbR(p,prw,prwcont,ipart,r,wc,e);
}

void AliQuenchingWeights::GetZeroLossProbR(Double_t &p,Double_t &prw,Double_t &prwcont,
					   Int_t ipart, Double_t r,Double_t wc,Double_t e)
{
  //calculate the probabilty that there is no energy
  //loss for different ways of energy constraint
  if(r>=fgkRMax) {
    r=fgkRMax-1;
  }  

  Double_t discrete=0.;
  Double_t continuous=0.;
  if (!fHisto) Init();
  Int_t kbinmax=fHisto->FindBin(e/wc);
  if(kbinmax>=fgkBins) kbinmax=fgkBins-1;
  if(fMultSoft) {
    for(Int_t bin=1; bin<=kbinmax; bin++) {
      Double_t xxxx = fHisto->GetBinCenter(bin);
      CalcMult(ipart,r,xxxx,continuous,discrete);
      fHisto->SetBinContent(bin,continuous);
    }
  } else {
    for(Int_t bin=1; bin<=kbinmax; bin++) {
      Double_t xxxx = fHisto->GetBinCenter(bin);
      CalcSingleHard(ipart,r,xxxx,continuous,discrete);
      fHisto->SetBinContent(bin,continuous);
    }
  }

  //non-reweighted P(Delta E = 0)
  const Double_t kdelta=fHisto->Integral(1,kbinmax);
  Double_t val=discrete*fgkBins/fgkMaxBin;
  fHisto->Fill(0.,val);
  fHisto->SetBinContent(kbinmax+1,(1-discrete)*fgkBins/fgkMaxBin-kdelta);
  Double_t hint=fHisto->Integral(1,kbinmax+1);
  p=fHisto->GetBinContent(1)/hint;

  // reweighted
  hint=fHisto->Integral(1,kbinmax);
  prw=fHisto->GetBinContent(1)/hint;

  Double_t xxxx = fHisto->GetBinCenter(1);
  CalcMult(ipart,r,xxxx,continuous,discrete);
  fHisto->SetBinContent(1,continuous);
  hint=fHisto->Integral(1,kbinmax);
  fHisto->Scale(1./hint*(1-discrete));
  fHisto->Fill(0.,discrete);
  prwcont=fHisto->GetBinContent(1);
}

Int_t AliQuenchingWeights::SampleEnergyLoss() 
{
  // Has to be called to fill the histograms
  //
  // For stored values fQTransport loop over 
  // particle type and length = 1 to fMaxLength (fm)
  // to fill energy loss histos
  //
  //    Take histogram of continuous weights 
  //    Take discrete_weight
  //    If discrete_weight > 1, put all channels to 0, except channel 1 
  //    Fill channel 1 with discrete_weight/(1-discrete_weight)*integral 

  // read-in data before first call
  if(!fTablesLoaded){
    Error("SampleEnergyLoss","Tables are not loaded.");
    return -1;
  }

  if(fMultSoft) {
    Int_t lmax=CalcLengthMax(fQTransport);
    if(fLengthMax>lmax){
      Info("SampleEnergyLoss","Maximum length changed from %d to %d;\nin order to have R < %.f",fLengthMax,lmax,fgkRMax);
      fLengthMax=lmax;
    }
  } else {
      Warning("SampleEnergyLoss","Maximum length not checked,\nbecause SingeHard is not yet tested.");
  }

  Reset();
  fHistos=new TH1F**[2];
  fHistos[0]=new TH1F*[4*fLengthMax];
  fHistos[1]=new TH1F*[4*fLengthMax];
  fLengthMaxOld=fLengthMax; //remember old value in case 
                            //user wants to reset

  Int_t medvalue=0;
  Char_t meddesc[100];
  if(fMultSoft) {
    medvalue=(Int_t)(fQTransport*1000.);
    snprintf(meddesc, 100, "MS");
  } else {
    medvalue=(Int_t)(fMu*1000.);
    snprintf(meddesc, 100, "SH");
  }

  for(Int_t ipart=1;ipart<=2;ipart++){
    for(Int_t l=1;l<=4*fLengthMax;l++){
      Char_t hname[100];
      snprintf(hname, 100, "hDisc-ContQW_%s_%d_%d_%d_%d",meddesc,fInstanceNumber,ipart,medvalue,l);
      Double_t len=l/4.;
      Double_t wc = CalcWC(len);
      fHistos[ipart-1][l-1] = new TH1F(hname,hname,fgkBins,0.,fgkMaxBin*wc);
      fHistos[ipart-1][l-1]->SetXTitle("#Delta E [GeV]");
      fHistos[ipart-1][l-1]->SetYTitle("p(#Delta E)");
      fHistos[ipart-1][l-1]->SetLineColor(4);

      Double_t rrrr = CalcR(wc,len);
      Double_t discrete=0.;
      // loop on histogram channels
      for(Int_t bin=1; bin<=fgkBins; bin++) {
	Double_t xxxx = fHistos[ipart-1][l-1]->GetBinCenter(bin)/wc;
	Double_t continuous;
	if(fMultSoft)
	  CalcMult(ipart,rrrr,xxxx,continuous,discrete);
	else
	  CalcSingleHard(ipart,rrrr,xxxx,continuous,discrete);
	fHistos[ipart-1][l-1]->SetBinContent(bin,continuous);
      }
      // add discrete part to distribution
      if(discrete>=1.)
	for(Int_t bin=2;bin<=fgkBins;bin++) 
	  fHistos[ipart-1][l-1]->SetBinContent(bin,0.);
      else {
	Double_t val=discrete/(1.-discrete)*fHistos[ipart-1][l-1]->Integral(1,fgkBins);
	fHistos[ipart-1][l-1]->Fill(0.,val);
      }
      Double_t hint=fHistos[ipart-1][l-1]->Integral(1,fgkBins);
      fHistos[ipart-1][l-1]->Scale(1./hint);
    }
  }
  return 0;
}

Int_t AliQuenchingWeights::SampleEnergyLoss(Int_t ipart, Double_t r)
{
  // Sample energy loss directly for one particle type
  // choses R (safe it and keep it until next call of function)

  // read-in data before first call
  if(!fTablesLoaded){
    Error("SampleEnergyLoss","Tables are not loaded.");
    return -1;
  }

  Double_t discrete=0.;
  Double_t continuous=0;;
  Int_t bin=1;
  if (!fHisto) Init();
  Double_t xxxx = fHisto->GetBinCenter(bin);
  if(fMultSoft)
    CalcMult(ipart,r,xxxx,continuous,discrete);
  else
    CalcSingleHard(ipart,r,xxxx,continuous,discrete);

  if(discrete>=1.) {
    fHisto->SetBinContent(1,1.);
    for(bin=2;bin<=fgkBins;bin++) 
      fHisto->SetBinContent(bin,0.);
    return 0;
  }

  fHisto->SetBinContent(bin,continuous);
  for(bin=2; bin<=fgkBins; bin++) {
    xxxx = fHisto->GetBinCenter(bin);
    if(fMultSoft)
      CalcMult(ipart,r,xxxx,continuous,discrete);
    else
      CalcSingleHard(ipart,r,xxxx,continuous,discrete);
    fHisto->SetBinContent(bin,continuous);
  }

  Double_t val=discrete/(1.-discrete)*fHisto->Integral(1,fgkBins);
  fHisto->Fill(0.,val);
  Double_t hint=fHisto->Integral(1,fgkBins);
  if(hint!=0)
    fHisto->Scale(1./hint);
  else {
    //cout << discrete << " " << hint << " " << continuous << endl;
    return -1;
  }
  return 0;
}

const TH1F* AliQuenchingWeights::GetHisto(Int_t ipart,Double_t length) const
{
  //return quenching histograms 
  //for ipart and length

  if(!fHistos){
    Fatal("GetELossRandom","Call SampleEnergyLoss method before!");
    return 0;
  }
  if((ipart<1) || (ipart>2)) {
    Fatal("GetELossRandom","ipart =%d; but has to be 1 (quark) or 2 (gluon)",ipart);
    return 0;
  }

  Int_t l=GetIndex(length);
  if(l<=0) return 0;
  return fHistos[ipart-1][l-1];
}

TH1F* AliQuenchingWeights::ComputeQWHisto(Int_t ipart,Double_t medval,Double_t length) const 
{
  // ipart = 1 for quark, 2 for gluon
  // medval a) qtransp = transport coefficient (GeV^2/fm)
  //        b) mu      = Debye mass (GeV)
  // length = path length in medium (fm)
  // Get from SW tables:
  // - continuous weight, as a function of dE/wc

  Double_t wc = 0;
  Char_t meddesc[100];
  if(fMultSoft) {
    wc=CalcWC(medval,length);
    snprintf(meddesc, 100, "MS");
  } else {
    wc=CalcWCbar(medval,length);
    snprintf(meddesc, 100, "SH");
  }

  Char_t hname[100];
  snprintf(hname, 100, "hContQWHisto_%s_%d_%d_%d",meddesc,ipart,
                (Int_t)(medval*1000.),(Int_t)length);

  TH1F *hist = new TH1F("hist",hname,fgkBins,0.,fgkMaxBin*wc);
  hist->SetXTitle("#Delta E [GeV]");
  hist->SetYTitle("p(#Delta E)");
  hist->SetLineColor(4);

  Double_t rrrr = CalcR(wc,length);
  //loop on histogram channels
  for(Int_t bin=1; bin<=fgkBins; bin++) {
    Double_t xxxx = hist->GetBinCenter(bin)/wc;
    Double_t continuous,discrete;
    Int_t ret=0;
    if(fMultSoft) ret=CalcMult(ipart,rrrr,xxxx,continuous,discrete);
    else ret=CalcSingleHard(ipart,rrrr,xxxx,continuous,discrete);
    if(ret!=0){
      delete hist;
      return 0;
    };
    hist->SetBinContent(bin,continuous);
  }
  return hist;
}

TH1F* AliQuenchingWeights::ComputeQWHistoX(Int_t ipart,Double_t medval,Double_t length) const 
{
  // ipart = 1 for quark, 2 for gluon
  // medval a) qtransp = transport coefficient (GeV^2/fm)
  //        b) mu      = Debye mass (GeV)
  // length = path length in medium (fm)
  // Get from SW tables:
  // - continuous weight, as a function of dE/wc

  Double_t wc = 0;
  Char_t meddesc[100];
  if(fMultSoft) {
    wc=CalcWC(medval,length);
    snprintf(meddesc, 100, "MS");
  } else {
    wc=CalcWCbar(medval,length);
    snprintf(meddesc, 100, "SH");
  }

  Char_t hname[100];
  snprintf(hname, 100, "hContQWHistox_%s_%d_%d_%d",meddesc,ipart,
                (Int_t)(medval*1000.),(Int_t)length);

  TH1F *histx = new TH1F("histx",hname,fgkBins,0.,fgkMaxBin);
  histx->SetXTitle("x = #Delta E/#omega_{c}");
  if(fMultSoft)
    histx->SetYTitle("p(#Delta E/#omega_{c})");
  else
    histx->SetYTitle("p(#Delta E/#bar#omega_{c})");
  histx->SetLineColor(4);

  Double_t rrrr = CalcR(wc,length);
  //loop on histogram channels
  for(Int_t bin=1; bin<=fgkBins; bin++) {
    Double_t xxxx = histx->GetBinCenter(bin);
    Double_t continuous,discrete;
    Int_t ret=0;
    if(fMultSoft) ret=CalcMult(ipart,rrrr,xxxx,continuous,discrete);
    else ret=CalcSingleHard(ipart,rrrr,xxxx,continuous,discrete);
    if(ret!=0){
      delete histx;
      return 0;
    };
    histx->SetBinContent(bin,continuous);
  }
  return histx;
}

TH1F* AliQuenchingWeights::ComputeQWHistoX(Int_t ipart,Double_t r) const 
{
  // compute P(E) distribution for
  // given ipart = 1 for quark, 2 for gluon 
  // and R

  Char_t meddesc[100];
  if(fMultSoft) {
    snprintf(meddesc, 100, "MS");
  } else {
    snprintf(meddesc, 100, "SH");
  }

  Char_t hname[100];
  snprintf(hname, 100, "hQWHistox_%s_%d_%.2f",meddesc,ipart,r);
  TH1F *histx = new TH1F("histx",hname,fgkBins,0.,fgkMaxBin);
  histx->SetXTitle("x = #Delta E/#omega_{c}");
  if(fMultSoft)
    histx->SetYTitle("p(#Delta E/#omega_{c})");
  else
    histx->SetYTitle("p(#Delta E/#bar#omega_{c})");
  histx->SetLineColor(4);

  Double_t rrrr = r;
  Double_t continuous=0.,discrete=0.;
  //loop on histogram channels
  for(Int_t bin=1; bin<=fgkBins; bin++) {
    Double_t xxxx = histx->GetBinCenter(bin);
    Int_t ret=0;
    if(fMultSoft) ret=CalcMult(ipart,rrrr,xxxx,continuous,discrete);
    else ret=CalcSingleHard(ipart,rrrr,xxxx,continuous,discrete);
    if(ret!=0){
      delete histx;
      return 0;
    };
    histx->SetBinContent(bin,continuous);
  }

  //add discrete part to distribution
  if(discrete>=1.)
    for(Int_t bin=2;bin<=fgkBins;bin++) 
      histx->SetBinContent(bin,0.);
  else {
    Double_t val=discrete/(1.-discrete)*histx->Integral(1,fgkBins);
    histx->Fill(0.,val);
  }
  Double_t hint=histx->Integral(1,fgkBins);
  if(hint!=0) histx->Scale(1./hint);

  return histx;
}

TH1F* AliQuenchingWeights::ComputeELossHisto(Int_t ipart,Double_t medval,Double_t l,Double_t e) const
{
  // compute energy loss histogram for 
  // parton type, medium value, length and energy

  AliQuenchingWeights *dummy=new AliQuenchingWeights(*this);
  if(fMultSoft){
    dummy->SetQTransport(medval);
    dummy->InitMult();
  } else {
    dummy->SetMu(medval);
    dummy->InitSingleHard();
  }
  dummy->SampleEnergyLoss();

  Char_t name[100];
  Char_t hname[100];
  if(ipart==1){
    snprintf(name, 100, "Energy Loss Distribution - Quarks;E_{loss} (GeV);#"); 
    snprintf(hname,100, "hLossQuarks"); 
  } else {
    snprintf(name, 100, "Energy Loss Distribution - Gluons;E_{loss} (GeV);#"); 
    snprintf(hname, 100, "hLossGluons"); 
  }

  TH1F *h = new TH1F(hname,name,250,0,250);
  for(Int_t i=0;i<100000;i++){
    //if(i % 1000 == 0) cout << "." << flush;
    Double_t loss=dummy->GetELossRandom(ipart,l,e);
    h->Fill(loss);
  }
  h->SetStats(kTRUE);
  delete dummy;
  return h;
}

TH1F* AliQuenchingWeights::ComputeELossHisto(Int_t ipart,Double_t medval,TH1F *hEll,Double_t e) const
{
  // compute energy loss histogram for 
  // parton type, medium value, 
  // length distribution and energy

  AliQuenchingWeights *dummy=new AliQuenchingWeights(*this);
  if(fMultSoft){
    dummy->SetQTransport(medval);
    dummy->InitMult();
  } else {
    dummy->SetMu(medval);
    dummy->InitSingleHard();
  }
  dummy->SampleEnergyLoss();

  Char_t name[100];
  Char_t hname[100];
  if(ipart==1){
    snprintf(name, 100, "Energy Loss Distribution - Quarks;E_{loss} (GeV);#"); 
    snprintf(hname,100, "hLossQuarks"); 
  } else {
    snprintf(name,100, "Energy Loss Distribution - Gluons;E_{loss} (GeV);#"); 
    snprintf(hname, 100, "hLossGluons"); 
  }

  TH1F *h = new TH1F(hname,name,250,0,250);
  for(Int_t i=0;i<100000;i++){
    //if(i % 1000 == 0) cout << "." << flush;
    Double_t loss=dummy->GetELossRandom(ipart,hEll,e);
    h->Fill(loss);
  }
  h->SetStats(kTRUE);
  delete dummy;
  return h;
}

TH1F* AliQuenchingWeights::ComputeELossHisto(Int_t ipart,Double_t r) const 
{
  // compute energy loss histogram for 
  // parton type and given R

  TH1F *dummy = ComputeQWHistoX(ipart,r);
  if(!dummy) return 0;

  Char_t hname[100];
  snprintf(hname, 100, "hELossHistox_%d_%.2f",ipart,r);
  TH1F *histx = new TH1F("histxr",hname,fgkBins,0.,fgkMaxBin);
  for(Int_t i=0;i<100000;i++){
    //if(i % 1000 == 0) cout << "." << flush;
    Double_t loss=dummy->GetRandom();
    histx->Fill(loss);
  }
  delete dummy;
  return histx;
}

Double_t AliQuenchingWeights::GetMeanELoss(Int_t ipart,Double_t medval,Double_t l) const
{
  // compute average energy loss for 
  // parton type, medium value, length and energy

  TH1F *dummy = ComputeELossHisto(ipart,medval,l);
  if(!dummy) return 0;
  Double_t ret=dummy->GetMean();
  delete dummy;
  return ret;
}

Double_t AliQuenchingWeights::GetMeanELoss(Int_t ipart,Double_t medval,TH1F *hEll) const
{
  // compute average energy loss for 
  // parton type, medium value, 
  // length distribution and energy

  TH1F *dummy = ComputeELossHisto(ipart,medval,hEll);
  if(!dummy) return 0;
  Double_t ret=dummy->GetMean();
  delete dummy;
  return ret;
}

Double_t  AliQuenchingWeights::GetMeanELoss(Int_t ipart,Double_t r) const 
{
  // compute average energy loss over wc 
  // for parton type and given R

  TH1F *dummy = ComputeELossHisto(ipart,r);
  if(!dummy) return 0;
  Double_t ret=dummy->GetMean();
  delete dummy;
  return ret;
}

void AliQuenchingWeights::PlotDiscreteWeights(Double_t len,Double_t qm) const
{
  // plot discrete weights for given length

  TCanvas *c; 
  if(fMultSoft) 
    c = new TCanvas("cdiscms","Discrete Weight for Multiple Scattering",0,0,500,400);
  else 
    c = new TCanvas("cdiscsh","Discrete Weight for Single Hard Scattering",0,0,500,400);
  c->cd();

  TH2F *hframe = new TH2F("hdisc","",2,0,qm+.1,2,0,1.25);
  hframe->SetStats(0);
  if(fMultSoft) 
    hframe->SetXTitle("#hat{q} [GeV^{2}/fm]");
  else
    hframe->SetXTitle("#mu [GeV]");
  //hframe->SetYTitle("Probability #Delta E = 0 , p_{0}");
  hframe->SetYTitle("p_{0} (discrete weight)");
  hframe->Draw();

  Int_t points=(Int_t)qm*4;
  TGraph *gq=new TGraph(points);
  Int_t i=0;
  if(fMultSoft) {
    for(Double_t q=0.05;q<=qm+.05;q+=0.25){
      Double_t disc,cont;
      CalcMult(1,1.0,q,len,cont,disc);
      gq->SetPoint(i,q,disc);i++;
    }
  } else {
    for(Double_t m=0.05;m<=qm+.05;m+=0.25){
      Double_t disc,cont;
      CalcSingleHard(1,1.0,m,len,cont, disc);
      gq->SetPoint(i,m,disc);i++;
    }
  }
  gq->SetMarkerStyle(20);
  gq->SetMarkerColor(1);
  gq->SetLineStyle(1);
  gq->SetLineColor(1);
  gq->Draw("l");

  TGraph *gg=new TGraph(points);
  i=0;
  if(fMultSoft){
    for(Double_t q=0.05;q<=qm+.05;q+=0.25){
      Double_t disc,cont;
      CalcMult(2,1.0,q,len,cont,disc);
      gg->SetPoint(i,q,disc);i++;
    }
  } else {
    for(Double_t m=0.05;m<=qm+.05;m+=0.25){
      Double_t disc,cont;
      CalcSingleHard(2,1.0,m,len,cont,disc);
      gg->SetPoint(i,m,disc);i++;
    }
  }
  gg->SetMarkerStyle(24);
  gg->SetMarkerColor(2);
  gg->SetLineStyle(2);
  gg->SetLineColor(2);
  gg->Draw("l");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "L = %.1f fm",len);
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark projectile","l");
  l1a->AddEntry(gg,"gluon projectile","l");
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotContWeights(Int_t itype,Double_t ell) const
{
  // plot continous weights for 
  // given parton type and length

  Float_t medvals[3];
  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft) {
    if(itype==1)
      snprintf(title, 1024, "Cont. Weight for Multiple Scattering - Quarks");
    else if(itype==2)
      snprintf(title, 1024, "Cont. Weight for Multiple Scattering - Gluons");
    else return;
    medvals[0]=4;medvals[1]=1;medvals[2]=0.5;
    snprintf(name, 1024, "ccont-ms-%d",itype);
  } else {
    if(itype==1)
      snprintf(title, 1024, "Cont. Weight for Single Hard Scattering - Quarks");
    else if(itype==2)
      snprintf(title, 1024, "Cont. Weight for Single Hard Scattering - Gluons");
    else return;
    medvals[0]=2;medvals[1]=1;medvals[2]=0.5;
    snprintf(name, 1024, "ccont-ms-%d",itype);
  }

  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH1F *h1=ComputeQWHisto(itype,medvals[0],ell); 
  h1->SetName("h1");
  h1->SetTitle(title); 
  h1->SetStats(0);
  h1->SetLineColor(1);
  h1->DrawCopy();
  TH1F *h2=ComputeQWHisto(itype,medvals[1],ell); 
  h2->SetName("h2");
  h2->SetLineColor(2);
  h2->DrawCopy("SAME");
  TH1F *h3=ComputeQWHisto(itype,medvals[2],ell); 
  h3->SetName("h3");
  h3->SetLineColor(3);
  h3->DrawCopy("SAME");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "L = %.1f fm",ell);
  l1a->AddEntry(h1,label,"");
  if(fMultSoft) {
    snprintf(label, 100, "#hat{q} = %.1f GeV^{2}/fm",medvals[0]);
    l1a->AddEntry(h1,label,"pl");
    snprintf(label, 100, "#hat{q} = %.1f GeV^{2}/fm",medvals[1]);
    l1a->AddEntry(h2,label,"pl");
    snprintf(label, 100, "#hat{q} = %.1f GeV^{2}/fm",medvals[2]);
    l1a->AddEntry(h3, label,"pl");
  } else {
    snprintf(label, 100, "#mu = %.1f GeV",medvals[0]);
    l1a->AddEntry(h1,label,"pl");
    snprintf(label, 100, "#mu = %.1f GeV",medvals[1]);
    l1a->AddEntry(h2,label,"pl");
    snprintf(label, 100, "#mu = %.1f GeV",medvals[2]);
    l1a->AddEntry(h3,label,"pl");
  }
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotContWeightsVsL(Int_t itype,Double_t medval) const
{
  // plot continous weights for 
  // given parton type and medium value

  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft) {
    if(itype==1)
      snprintf(title,1024, "Cont. Weight for Multiple Scattering - Quarks");
    else if(itype==2)
      snprintf(title,1024, "Cont. Weight for Multiple Scattering - Gluons");
    else return;
    snprintf(name,1024, "ccont2-ms-%d",itype);
  } else {
    if(itype==1)
      snprintf(title, 1024, "Cont. Weight for Single Hard Scattering - Quarks");
    else if(itype==2)
      snprintf(title, 1024, "Cont. Weight for Single Hard Scattering - Gluons");
    else return;
    snprintf(name, 1024, "ccont2-sh-%d",itype);
  }
  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH1F *h1=ComputeQWHisto(itype,medval,8); 
  h1->SetName("h1");
  h1->SetTitle(title); 
  h1->SetStats(0);
  h1->SetLineColor(1);
  h1->DrawCopy();
  TH1F *h2=ComputeQWHisto(itype,medval,5); 
  h2->SetName("h2");
  h2->SetLineColor(2);
  h2->DrawCopy("SAME");
  TH1F *h3=ComputeQWHisto(itype,medval,2); 
  h3->SetName("h3");
  h3->SetLineColor(3);
  h3->DrawCopy("SAME");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  if(fMultSoft)
    snprintf(label, 100, "#hat{q} = %.1f GeV^{2}/fm",medval);
  else
    snprintf(label, 100, "#mu = %.1f GeV",medval);

  l1a->AddEntry(h1,label,"");
  l1a->AddEntry(h1,"L = 8 fm","pl");
  l1a->AddEntry(h2,"L = 5 fm","pl");
  l1a->AddEntry(h3,"L = 2 fm","pl");
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotAvgELoss(Double_t len,Double_t qm,Double_t e) const
{
  // plot average energy loss for given length
  // and parton energy 

  if(!fTablesLoaded){
    Error("PlotAvgELoss","Tables are not loaded.");
    return;
  }

  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft){ 
    snprintf(title, 1024, "Average Energy Loss for Multiple Scattering");
    snprintf(name, 1024, "cavgelossms");
  } else {
    snprintf(title, 1024, "Average Energy Loss for Single Hard Scattering");
    snprintf(name, 1024, "cavgelosssh");
  }

  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH2F *hframe = new TH2F("avgloss","",2,0,qm+.1,2,0,100);
  hframe->SetStats(0);
  if(fMultSoft) 
    hframe->SetXTitle("#hat{q} [GeV^{2}/fm]");
  else
    hframe->SetXTitle("#mu [GeV]");
  hframe->SetYTitle("<E_{loss}> [GeV]");
  hframe->Draw();

  TGraph *gq=new TGraph(20);
  Int_t i=0;
  for(Double_t v=0.05;v<=qm+.05;v+=0.25){
    TH1F *dummy=ComputeELossHisto(1,v,len,e);
    Double_t avgloss=dummy->GetMean();
    gq->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gq->SetMarkerStyle(21);
  gq->Draw("pl");

  Int_t points=(Int_t)qm*4;
  TGraph *gg=new TGraph(points);
  i=0;
  for(Double_t v=0.05;v<=qm+.05;v+=0.25){
    TH1F *dummy=ComputeELossHisto(2,v,len,e);
    Double_t avgloss=dummy->GetMean();
    gg->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gg->SetMarkerStyle(20);
  gg->SetMarkerColor(2);
  gg->Draw("pl");

  TGraph *gratio=new TGraph(points);
  for(i=0;i<points;i++){
    Double_t x,y,x2,y2;
    gg->GetPoint(i,x,y);
    gq->GetPoint(i,x2,y2);
    if(y2>0)
      gratio->SetPoint(i,x,y/y2*10/2.25);
    else gratio->SetPoint(i,x,0);
  }
  gratio->SetLineStyle(4);
  gratio->Draw();
  TLegend *l1a = new TLegend(0.15,0.60,0.50,0.90);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "L = %.1f fm",len);
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark projectile","pl");
  l1a->AddEntry(gg,"gluon projectile","pl");
  l1a->AddEntry(gratio,"gluon/quark/2.25*10","pl");
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotAvgELoss(TH1F *hEll,Double_t e) const
{
  // plot average energy loss for given
  // length distribution and parton energy

  if(!fTablesLoaded){
    Error("PlotAvgELossVs","Tables are not loaded.");
    return;
  }

  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft){ 
    snprintf(title, 1024, "Average Energy Loss for Multiple Scattering");
    snprintf(name, 1024, "cavgelossms2");
  } else {
    snprintf(title, 1024, "Average Energy Loss for Single Hard Scattering");
    snprintf(name, 1024, "cavgelosssh2");
  }

  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH2F *hframe = new TH2F("havgloss",title,2,0,5.1,2,0,100);
  hframe->SetStats(0);
  if(fMultSoft) 
    hframe->SetXTitle("#hat{q} [GeV^{2}/fm]");
  else
    hframe->SetXTitle("#mu [GeV]");
  hframe->SetYTitle("<E_{loss}> [GeV]");
  hframe->Draw();

  TGraph *gq=new TGraph(20);
  Int_t i=0;
  for(Double_t v=0.05;v<=5.05;v+=0.25){
    TH1F *dummy=ComputeELossHisto(1,v,hEll,e);
    Double_t avgloss=dummy->GetMean();
    gq->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gq->SetMarkerStyle(20);
  gq->Draw("pl");

  TGraph *gg=new TGraph(20);
  i=0;
  for(Double_t v=0.05;v<=5.05;v+=0.25){
    TH1F *dummy=ComputeELossHisto(2,v,hEll,e);
    Double_t avgloss=dummy->GetMean();
    gg->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gg->SetMarkerStyle(24);
  gg->Draw("pl");

  TGraph *gratio=new TGraph(20);
  for(i=0;i<20;i++){
    Double_t x,y,x2,y2;
    gg->GetPoint(i,x,y);
    gq->GetPoint(i,x2,y2);
    if(y2>0)
      gratio->SetPoint(i,x,y/y2*10/2.25);
    else gratio->SetPoint(i,x,0);
  }
  gratio->SetLineStyle(4);
  //gratio->Draw();

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "<L> = %.2f fm",hEll->GetMean());
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark","pl");
  l1a->AddEntry(gg,"gluon","pl");
  //l1a->AddEntry(gratio,"gluon/quark/2.25*10","pl");
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotAvgELossVsL(Double_t e)  const
{
  // plot average energy loss versus ell

  if(!fTablesLoaded){
    Error("PlotAvgELossVsEll","Tables are not loaded.");
    return;
  }

  Char_t title[1024];
  Char_t name[1024];
  Float_t medval;
  if(fMultSoft){ 
    snprintf(title, 1024, "Average Energy Loss for Multiple Scattering");
    snprintf(name, 1024,  "cavgelosslms");
    medval=fQTransport;
  } else {
    snprintf(title, 1024, "Average Energy Loss for Single Hard Scattering");
    snprintf(name, 1024,  "cavgelosslsh");
    medval=fMu;
  }

  TCanvas *c = new TCanvas(name,title,0,0,600,400);
  c->cd();
  TH2F *hframe = new TH2F("avglossell",title,2,0,fLengthMax,2,0,250);
  hframe->SetStats(0);
  hframe->SetXTitle("length [fm]");
  hframe->SetYTitle("<E_{loss}> [GeV]");
  hframe->Draw();

  TGraph *gq=new TGraph((Int_t)fLengthMax*4);
  Int_t i=0;
  for(Double_t v=0.25;v<=fLengthMax;v+=0.25){
    TH1F *dummy=ComputeELossHisto(1,medval,v,e);
    Double_t avgloss=dummy->GetMean();
    gq->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gq->SetMarkerStyle(20);
  gq->Draw("pl");

  TGraph *gg=new TGraph((Int_t)fLengthMax*4);
  i=0;
  for(Double_t v=0.25;v<=fLengthMax;v+=0.25){
    TH1F *dummy=ComputeELossHisto(2,medval,v,e);
    Double_t avgloss=dummy->GetMean();
    gg->SetPoint(i,v,avgloss);i++;
    delete dummy;
  }
  gg->SetMarkerStyle(24);
  gg->Draw("pl");

  TGraph *gratio=new TGraph((Int_t)fLengthMax*4);
  for(i=0;i<=(Int_t)fLengthMax*4;i++){
    Double_t x,y,x2,y2;
    gg->GetPoint(i,x,y);
    gq->GetPoint(i,x2,y2);
    if(y2>0)
      gratio->SetPoint(i,x,y/y2*100/2.25);
    else gratio->SetPoint(i,x,0);
  }
  gratio->SetLineStyle(1);
  gratio->SetLineWidth(2);
  gratio->Draw();
  TLegend *l1a = new TLegend(0.15,0.65,.60,0.85);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  if(fMultSoft) 
    snprintf(label, 100, "#hat{q} = %.2f GeV^{2}/fm",medval);
  else
    snprintf(label, 100, "#mu = %.2f GeV",medval);
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark","pl");
  l1a->AddEntry(gg,"gluon","pl");
  l1a->AddEntry(gratio,"gluon/quark/2.25*100","pl");
  l1a->Draw();

  TF1 *f=new TF1("f","100",0,fLengthMax);
  f->SetLineStyle(4);
  f->SetLineWidth(1);
  f->Draw("same");
  c->Update();
}

void AliQuenchingWeights::PlotAvgELossVsPt(Double_t medval,Double_t len) const
{
  // plot relative energy loss for given
  // length and parton energy versus pt

  if(!fTablesLoaded){
    Error("PlotAvgELossVsPt","Tables are not loaded.");
    return;
  }

  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft){
    snprintf(title, 1024, "Relative Energy Loss for Multiple Scattering");
    snprintf(name, 1024, "cavgelossvsptms");
  } else {
    snprintf(title, 1024, "Relative Energy Loss for Single Hard Scattering");
    snprintf(name, 1024, "cavgelossvsptsh");
  }

  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH2F *hframe = new TH2F("havglossvspt",title,2,0,100,2,0,1);
  hframe->SetStats(0);
  hframe->SetXTitle("p_{T} [GeV]");
  hframe->SetYTitle("<E_{loss}>/p_{T} [GeV]");
  hframe->Draw();

  TGraph *gq=new TGraph(40);
  Int_t i=0;
  for(Double_t pt=2.5;pt<=100.05;pt+=2.5){
    TH1F *dummy=ComputeELossHisto(1,medval,len,pt);
    Double_t avgloss=dummy->GetMean();
    gq->SetPoint(i,pt,avgloss/pt);i++;
    delete dummy;
  }
  gq->SetMarkerStyle(20);
  gq->Draw("pl");

  TGraph *gg=new TGraph(40);
  i=0;
  for(Double_t pt=2.5;pt<=100.05;pt+=2.5){
    TH1F *dummy=ComputeELossHisto(2,medval,len,pt);
    Double_t avgloss=dummy->GetMean();
    gg->SetPoint(i,pt,avgloss/pt);i++;
    delete dummy;
  }
  gg->SetMarkerStyle(24);
  gg->Draw("pl");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "L = %.1f fm",len);
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark","pl");
  l1a->AddEntry(gg,"gluon","pl");
  l1a->Draw();

  c->Update();
}

void AliQuenchingWeights::PlotAvgELossVsPt(Double_t medval,TH1F *hEll) const
{
  // plot relative energy loss for given
  // length distribution and parton energy versus pt

  if(!fTablesLoaded){
    Error("PlotAvgELossVsPt","Tables are not loaded.");
    return;
  }

  Char_t title[1024];
  Char_t name[1024];
  if(fMultSoft){
    snprintf(title, 1024, "Relative Energy Loss for Multiple Scattering");
    snprintf(name,  1024, "cavgelossvsptms2");
  } else {
    snprintf(title, 1024, "Relative Energy Loss for Single Hard Scattering");
    snprintf(name,  1024, "cavgelossvsptsh2");
  }
  TCanvas *c = new TCanvas(name,title,0,0,500,400);
  c->cd();
  TH2F *hframe = new TH2F("havglossvspt",title,2,0,100,2,0,1);
  hframe->SetStats(0);
  hframe->SetXTitle("p_{T} [GeV]");
  hframe->SetYTitle("<E_{loss}>/p_{T} [GeV]");
  hframe->Draw();

  TGraph *gq=new TGraph(40);
  Int_t i=0;
  for(Double_t pt=2.5;pt<=100.05;pt+=2.5){
    TH1F *dummy=ComputeELossHisto(1,medval,hEll,pt);
    Double_t avgloss=dummy->GetMean();
    gq->SetPoint(i,pt,avgloss/pt);i++;
    delete dummy;
  }
  gq->SetMarkerStyle(20);
  gq->Draw("pl");

  TGraph *gg=new TGraph(40);
  i=0;
  for(Double_t pt=2.5;pt<=100.05;pt+=2.5){
    TH1F *dummy=ComputeELossHisto(2,medval,hEll,pt);
    Double_t avgloss=dummy->GetMean();
    gg->SetPoint(i,pt,avgloss/pt);i++;
    delete dummy;
  }
  gg->SetMarkerStyle(24);
  gg->Draw("pl");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  snprintf(label, 100, "<L> = %.2f fm",hEll->GetMean());
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark","pl");
  l1a->AddEntry(gg,"gluon","pl");
  l1a->Draw();

  c->Update();
}

Int_t AliQuenchingWeights::GetIndex(Double_t len) const
{
  //get the index according to length
  if(len>fLengthMax) len=fLengthMax;

  Int_t l=Int_t(len/0.25);
  if((len-l*0.25)>0.125) l++;
  return l;
}

