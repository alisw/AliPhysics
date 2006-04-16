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

/*
$Log$
Revision 1.6  2006/04/11 15:28:32  hristov
Checks on cache status before deleting calibration objects (A.Colla)

Revision 1.5  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.4  2006/03/31 11:26:46  arcelli
 changing CDB Ids according to standard convention

Revision 1.3  2006/03/28 14:57:02  arcelli
updates to handle new V5 geometry & some re-arrangements

Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFcalib.h"
#include "AliRun.h"
#include <TTask.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include "AliTOF.h"
#include "AliTOFcalibESD.h"
#include "AliESD.h"
#include <TObject.h>
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliESDtrack.h"
#include "AliTOFChannel.h"
#include "AliTOFChSim.h"
#include "AliTOFGeometryV5.h"
#include "TClonesArray.h"
#include "AliTOFCal.h"
#include "TRandom.h"
#include "AliTOFcluster.h"
#include "TList.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"

extern TROOT *gROOT;
extern TStyle *gStyle;

ClassImp(AliTOFcalib)

const Int_t AliTOFcalib::fgkchannel = 5000;
//_______________________________________________________________________
AliTOFcalib::AliTOFcalib():TTask("AliTOFcalib",""){ 
  //TOF Calibration Class ctor
  fArrayToT = 0x0;
  fArrayTime = 0x0;
  fESDsel = 0x0;
  AliTOFGeometry *geom=new AliTOFGeometryV5();
  AliInfo("V5 TOF Geometry is taken as the default");
  fNSector = geom->NSectors();
  fNPlate  = geom->NPlates();
  fNStripA = geom->NStripA();
  fNStripB = geom->NStripB();
  fNStripC = geom->NStripC();
  fNpadZ = geom->NpadZ();
  fNpadX = geom->NpadX();
  fNChannels = fNSector*(2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX; //generalized version
  fTOFCal = 0x0;
  fTOFSimCal = 0x0;
  fTOFSimToT = 0x0;
  delete geom;
}
//_______________________________________________________________________
AliTOFcalib::AliTOFcalib(AliTOFGeometry *geom):TTask("AliTOFcalib",""){ 
  //TOF Calibration Class ctor, taking the TOF geometry as input
  fArrayToT = 0x0;
  fArrayTime = 0x0;
  fESDsel = 0x0;
  fNSector = geom->NSectors();
  fNPlate  = geom->NPlates();
  fNStripA = geom->NStripA();
  fNStripB = geom->NStripB();
  fNStripC = geom->NStripC();
  fNpadZ = geom->NpadZ();
  fNpadX = geom->NpadX();
  fNChannels = fNSector*(2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX; //generalized version
  fTOFCal = 0x0;
  fTOFSimCal = 0x0;
  fTOFSimToT = 0x0;
}
//____________________________________________________________________________ 

AliTOFcalib::AliTOFcalib(const AliTOFcalib & calib):TTask("AliTOFcalib","")
{
  //TOF Calibration Class copy ctor
  fNSector = calib.fNSector;
  fNPlate = calib.fNPlate;
  fNStripA = calib.fNStripA;
  fNStripB = calib.fNStripB;
  fNStripC = calib.fNStripC;
  fNpadZ = calib.fNpadZ;
  fNpadX = calib.fNpadX;
  fNChannels = calib.fNChannels;
  fArrayToT = calib.fArrayToT;
  fArrayTime = calib.fArrayTime;
  fTOFCal=calib.fTOFCal;
  fTOFSimCal = calib.fTOFSimCal;
  fTOFSimToT=calib.fTOFSimToT;
}

//____________________________________________________________________________ 

AliTOFcalib::~AliTOFcalib()
{
  //TOF Calibration Class dtor
  delete fArrayToT;
  delete fArrayTime;

  if(!(AliCDBManager::Instance()->GetCacheFlag())){ // CDB objects must NOT be deleted if cache is active!
  	delete fTOFCal;
  	delete fTOFSimCal;
  }
  
  delete fESDsel;
}
//__________________________________________________________________________

TF1* AliTOFcalib::SetFitFunctions(TH1F *histo)
{
  //Define Fit Functions for Slewing Correction
  TF1 * fpol[3];
  const Int_t knbins = histo->GetNbinsX();
  Float_t delta = histo->GetBinWidth(1);  //all the bins have the same width
  Double_t max = histo->GetBinLowEdge(knbins)+delta;
  max = 15;
  fpol[0]=new TF1("poly3","pol3",5,max);
  fpol[1]=new TF1("poly4","pol4",5,max);
  fpol[2]=new TF1("poly5","pol5",5,max);
  char npoly[10];
  Double_t chi[3]={1E6,1E6,1E6};
  Int_t ndf[3]={-1,-1,-1};
  Double_t nchi[3]={1E6,1E6,1E6};
  Double_t bestchi=1E6;
  TF1 * fGold=0x0;
  Int_t nonzero =0;
  Int_t numberOfpar =0;
  for (Int_t j=0; j<knbins; j++){
    if (histo->GetBinContent(j)!=0) {
      nonzero++;
    }
  }
  Int_t norderfit = 0;
  if (nonzero<=4) {
    AliError(" Too few points in the histo. No fit performed.");
    return 0x0;
  }
  else if (nonzero<=5) {
    norderfit = 3;
    AliInfo(" Only 3rd order polynomial fit possible.");
  }
  else if (nonzero<=6) {
    norderfit = 4;
    AliInfo(" Only 3rd and 4th order polynomial fit possible.");
  }
  else {
    norderfit = 5;
    AliInfo(" All 3rd, 4th and 5th order polynomial fit possible.");
  }
  for (Int_t ifun=norderfit-3;ifun<norderfit-2;ifun++){
    sprintf(npoly,"poly%i",ifun+3);
    histo->Fit(npoly, "ERN", " ", 5.,14.);
    chi[ifun] = fpol[ifun]->GetChisquare();
    ndf[ifun] = fpol[ifun]->GetNDF();
    nchi[ifun] = (Double_t)chi[ifun]/ndf[ifun];
    if (nchi[ifun]<bestchi) {
      bestchi=nchi[ifun];
      fGold = fpol[ifun];
      numberOfpar = fGold->GetNpar();
    }
  }
  fGold=fpol[2];  //Gold fit function set to pol5 in any case
  histo->Fit(fGold,"ER " ," ",5.,15.);
  return fGold;
}
//____________________________________________________________________________

void AliTOFcalib::SelectESD(AliESD *event) 
{
  //track selection for Calibration
  Float_t lowerMomBound=0.8; // [GeV/c] default value Pb-Pb
  Float_t upperMomBound=1.8 ; // [GeV/c] default value Pb-Pb
  Int_t ntrk =0;
  Int_t ngoodtrkfinalToT = 0;
  ntrk=event->GetNumberOfTracks();
  fESDsel = new TObjArray(ntrk);
  fESDsel->SetOwner();
  TObjArray  uCdatatemp(ntrk);
  Int_t ngoodtrk = 0;
  Int_t ngoodtrkfinal = 0;
  Float_t mintime =1E6;
  for (Int_t itrk=0; itrk<ntrk; itrk++) {
    AliESDtrack *t=event->GetTrack(itrk);
    //track selection: reconstrution to TOF:
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) {
      continue;
    }
    //IsStartedTimeIntegral
    if ((t->GetStatus()&AliESDtrack::kTIME)==0) {
      continue;
    }
    Double_t time=t->GetTOFsignal();	
    time*=1.E-3; // tof given in nanoseconds
    if(time <= mintime)mintime=time;
    Double_t mom=t->GetP();
    if (!(mom<=upperMomBound && mom>=lowerMomBound))continue;
    UInt_t assignedTOFcluster=t->GetTOFcluster();//index of the assigned TOF cluster, >0 ?
    if(assignedTOFcluster==0){ // not matched
      continue;
    }
    AliTOFcalibESD *unc = new AliTOFcalibESD;
    unc->CopyFromAliESD(t);
    Double_t c1[15]; 
    unc->GetExternalCovariance(c1);
    uCdatatemp.Add(unc);
    ngoodtrk++;
  }
  for (Int_t i = 0; i < ngoodtrk ; i ++){
    AliTOFcalibESD *unc = (AliTOFcalibESD*)uCdatatemp.At(i);
    if((unc->GetTOFsignal()-mintime*1.E3)<5.E3){
      fESDsel->Add(unc);
      ngoodtrkfinal++;
      ngoodtrkfinalToT++;
    }
  }
  fESDsel->Sort();
}
//_____________________________________________________________________________

void AliTOFcalib::CombESDId()
{
  //track PID for calibration
  Float_t t0offset=0;
  Float_t loffset=0;
  Int_t   ntracksinset=6;
  Float_t exptof[6][3];
  Float_t momentum[6]={0.,0.,0.,0.,0.,0.};
  Int_t   assparticle[6]={3,3,3,3,3,3};
  Float_t massarray[3]={0.13957,0.493677,0.9382723};
  Float_t timeofflight[6]={0.,0.,0.,0.,0.,0.};
  Float_t beta[6]={0.,0.,0.,0.,0.,0.};
  Float_t texp[6]={0.,0.,0.,0.,0.,0.};
  Float_t sqMomError[6]={0.,0.,0.,0.,0.,0.};
  Float_t sqTrackError[6]={0.,0.,0.,0.,0.,0.};
  Float_t tracktoflen[6]={0.,0.,0.,0.,0.,0.};
  Float_t timeResolution   = 0.90e-10; // 90 ps by default	
  Float_t timeresolutioninns=timeResolution*(1.e+9); // convert in [ns]
  Float_t timezero[6]={0.,0.,0.,0.,0.,0.};
  Float_t weightedtimezero[6]={0.,0.,0.,0.,0.,0.};
  Float_t besttimezero[6]={0.,0.,0.,0.,0.,0.};
  Float_t bestchisquare[6]={0.,0.,0.,0.,0.,0.};
  Float_t bestweightedtimezero[6]={0.,0.,0.,0.,0.,0.};
  Float_t bestsqTrackError[6]={0.,0.,0.,0.,0.,0.};

  Int_t nelements = fESDsel->GetEntries();
  Int_t nset= (Int_t)(nelements/ntracksinset);
  for (Int_t i=0; i< nset; i++) {   

    AliTOFcalibESD **unc=new AliTOFcalibESD*[ntracksinset];
    for (Int_t itrk=0; itrk<ntracksinset; itrk++) {
      Int_t index = itrk+i*ntracksinset;
      AliTOFcalibESD *element=(AliTOFcalibESD*)fESDsel->At(index);
      unc[itrk]=element;
    }
    
    for (Int_t j=0; j<ntracksinset; j++) {
      AliTOFcalibESD *element=unc[j];
      Double_t mom=element->GetP();
      Double_t time=element->GetTOFsignal()*1.E-3; // in ns	
      Double_t exptime[10]; 
      element->GetIntegratedTimes(exptime);
      Double_t toflen=element->GetIntegratedLength()/100.;  // in m
      timeofflight[j]=time+t0offset;
      tracktoflen[j]=toflen+loffset;
      exptof[j][0]=exptime[2]*1.E-3+0.005;
      exptof[j][1]=exptime[3]*1.E-3+0.005;
      exptof[j][2]=exptime[4]*1.E-3+0.005;
      momentum[j]=mom;
    }
    Float_t t0best=999.;
    Float_t et0best=999.;
    Float_t chisquarebest=999.;
    for (Int_t i1=0; i1<3;i1++) {
      beta[0]=momentum[0]/sqrt(massarray[i1]*massarray[i1]+momentum[0]*momentum[0]);
      texp[0]=exptof[0][i1];
      for (Int_t i2=0; i2<3;i2++) { 
	beta[1]=momentum[1]/sqrt(massarray[i2]*massarray[i2]+momentum[1]*momentum[1]);
	texp[1]=exptof[1][i2];
	for (Int_t i3=0; i3<3;i3++) {
	  beta[2]=momentum[2]/sqrt(massarray[i3]*massarray[i3]+momentum[2]*momentum[2]);
	  texp[2]=exptof[2][i3];
	  for (Int_t i4=0; i4<3;i4++) {
	    beta[3]=momentum[3]/sqrt(massarray[i4]*massarray[i4]+momentum[3]*momentum[3]);
	    texp[3]=exptof[3][i4];
	    
	    for (Int_t i5=0; i5<3;i5++) {
	      beta[4]=momentum[4]/sqrt(massarray[i5]*massarray[i5]+momentum[4]*momentum[4]);
	      texp[4]=exptof[4][i5];
	      for (Int_t i6=0; i6<3;i6++) {
		beta[5]=momentum[5]/sqrt(massarray[i6]*massarray[i6]+momentum[5]*momentum[5]);
		texp[5]=exptof[5][i6];
	
		Float_t sumAllweights=0.;
		Float_t meantzero=0.;
		Float_t emeantzero=0.;
		
		for (Int_t itz=0; itz<ntracksinset;itz++) {
		  sqMomError[itz]=
		    ((1.-beta[itz]*beta[itz])*0.025)*
		    ((1.-beta[itz]*beta[itz])*0.025)*
		    (tracktoflen[itz]/
		     (0.299792*beta[itz]))*
		    (tracktoflen[itz]/
		     (0.299792*beta[itz])); 
		  sqTrackError[itz]=
		    (timeresolutioninns*
		     timeresolutioninns
		     +sqMomError[itz]); 
		  
		  timezero[itz]=texp[itz]-timeofflight[itz];		    
		  weightedtimezero[itz]=timezero[itz]/sqTrackError[itz];
		  sumAllweights+=1./sqTrackError[itz];
		  meantzero+=weightedtimezero[itz];
		  
		} // end loop for (Int_t itz=0; itz<15;itz++)
		
		meantzero=meantzero/sumAllweights; // it is given in [ns]
		emeantzero=sqrt(1./sumAllweights); // it is given in [ns]
		
		// calculate chisquare
		
		Float_t chisquare=0.;		
		for (Int_t icsq=0; icsq<ntracksinset;icsq++) {
		  chisquare+=(timezero[icsq]-meantzero)*(timezero[icsq]-meantzero)/sqTrackError[icsq];
		} // end loop for (Int_t icsq=0; icsq<15;icsq++) 
		//		cout << " chisquare " << chisquare << endl;
		
		Int_t npion=0;
		if(i1==0)npion++;
		if(i2==0)npion++;
		if(i3==0)npion++;
		if(i4==0)npion++;
		if(i5==0)npion++;
		if(i6==0)npion++;
		
	     	if(chisquare<=chisquarebest  && ((Float_t) npion/ ((Float_t) ntracksinset)>0.3)){
		  //  if(chisquare<=chisquarebest){
		  
		  for(Int_t iqsq = 0; iqsq<ntracksinset; iqsq++) {
		    bestsqTrackError[iqsq]=sqTrackError[iqsq]; 
		    besttimezero[iqsq]=timezero[iqsq]; 
		    bestweightedtimezero[iqsq]=weightedtimezero[iqsq]; 
		    bestchisquare[iqsq]=(timezero[iqsq]-meantzero)*(timezero[iqsq]-meantzero)/sqTrackError[iqsq]; 
		  }
		  
		  assparticle[0]=i1;
		  assparticle[1]=i2;
		  assparticle[2]=i3;
		  assparticle[3]=i4;
		  assparticle[4]=i5;
		  assparticle[5]=i6;
		  
		  chisquarebest=chisquare;
     		  t0best=meantzero;
		  et0best=emeantzero;
		} // close if(dummychisquare<=chisquare)
	      } // end loop on i6
	    } // end loop on i5
	  } // end loop on i4
	} // end loop on i3
      } // end loop on i2
    } // end loop on i1


    Float_t confLevel=999;
    if(chisquarebest<999.){
      Double_t dblechisquare=(Double_t)chisquarebest;
      confLevel=(Float_t)TMath::Prob(dblechisquare,ntracksinset-1); 
    }
    // assume they are all pions for fake sets
    if(confLevel<0.01 || confLevel==999. ){
      for (Int_t itrk=0; itrk<ntracksinset; itrk++)assparticle[itrk]=0;
    }
    for (Int_t itrk=0; itrk<ntracksinset; itrk++) {
      Int_t index = itrk+i*ntracksinset;
      AliTOFcalibESD *element=(AliTOFcalibESD*)fESDsel->At(index);
      element->SetCombID(assparticle[itrk]);
    }
  }
}

//_____________________________________________________________________________

void AliTOFcalib::CalibrateESD(){
  //Calibrate selected ESD times
  Int_t nelements = fESDsel->GetEntries();
  Int_t *number=new Int_t[fNChannels];
  fArrayToT = new AliTOFArray(fNChannels);
  fArrayTime = new AliTOFArray(fNChannels);
  for (Int_t i=0; i<fNChannels; i++){
    number[i]=0;
    fArrayToT->AddArray(i, new TArrayF(fgkchannel));
    TArrayF * parrToT = fArrayToT->GetArray(i);
    TArrayF & refaToT = * parrToT;
    fArrayTime->AddArray(i, new TArrayF(fgkchannel));
    TArrayF * parrTime = fArrayToT->GetArray(i);
    TArrayF & refaTime = * parrTime;
    for (Int_t j = 0;j<AliTOFcalib::fgkchannel;j++){
      refaToT[j]=0.;      //ToT[i][j]=j;
      refaTime[j]=0.;      //Time[i][j]=j;
    }
  }
  
  for (Int_t i=0; i< nelements; i++) {
    AliTOFcalibESD *element=(AliTOFcalibESD*)fESDsel->At(i);
    Int_t ipid = element->GetCombID();
    Double_t etime = 0;   //expected time
    Double_t expTime[10]; 
    element->GetIntegratedTimes(expTime);
    if (ipid == 0) etime = expTime[2]*1E-3; //ns
    else if (ipid == 1) etime = expTime[3]*1E-3; //ns
    else if (ipid == 2) etime = expTime[4]*1E-3; //ns
    else AliError("No pid from combinatorial algo for this track");
    Double_t mtime = (Double_t)element->GetTOFsignal()*1E-3;  //measured time
    Double_t mToT = (Double_t) element->GetToT();  //measured ToT, ns
    //select the correspondent channel with its simulated ToT spectrum
    //summing up everything, index = 0 for all channels:
    Int_t index = element->GetTOFCalChannel();
    Int_t index2 = number[index];
    TArrayF * parrToT = fArrayToT->GetArray(index);
    TArrayF & refaToT = * parrToT;
    refaToT[index2] = (Float_t)mToT;
    TArrayF * parrTime = fArrayTime->GetArray(index);
    TArrayF & refaTime = * parrTime;
    refaTime[index2] = (Float_t)(mtime-etime);
    number[index]++;
  }

  for (Int_t i=0;i<1;i++){
    TH1F * hProf = Profile(i);
    TF1* fGold = SetFitFunctions(hProf);
    Int_t nfpar = fGold->GetNpar();
    Float_t par[6];    
    for(Int_t kk=0;kk<6;kk++){
      par[kk]=0;
    }
    for (Int_t kk = 0; kk< nfpar; kk++){
      par[kk]=fGold->GetParameter(kk);
    }
    if (!fTOFCal) {
      AliTOFGeometry *geom=new AliTOFGeometryV5();
      fTOFCal = new AliTOFCal(geom);
      fTOFCal->CreateArray();
      delete geom;
    }
    AliTOFChannel * calChannel = fTOFCal->GetChannel(i);
    calChannel->SetSlewPar(par);
  }
  delete[] number;
}

//___________________________________________________________________________

TH1F* AliTOFcalib::Profile(Int_t ich)
{
  //Prepare histograms for Slewing Correction
  const Int_t knbinToT = 650;
  Int_t nbinTime = 400;
  Float_t minTime = -10.5; //ns
  Float_t maxTime = 10.5; //ns
  Float_t minToT = 7.5; //ns
  Float_t maxToT = 40.; //ns
  Float_t deltaToT = (maxToT-minToT)/knbinToT;
  Double_t mTime[knbinToT+1],mToT[knbinToT+1],meanTime[knbinToT+1], meanTime2[knbinToT+1],vToT[knbinToT+1], vToT2[knbinToT+1],meanToT[knbinToT+1],meanToT2[knbinToT+1],vTime[knbinToT+1],vTime2[knbinToT+1],xlow[knbinToT+1],sigmaTime[knbinToT+1];
  Int_t n[knbinToT+1], nentrx[knbinToT+1];
  Double_t sigmaToT[knbinToT+1];
  for (Int_t i = 0; i < knbinToT+1 ; i++){
    mTime[i]=0;
    mToT[i]=0;
    n[i]=0;
    meanTime[i]=0;
    meanTime2[i]=0;
    vToT[i]=0;
    vToT2[i]=0;
    meanToT[i]=0;
    meanToT2[i]=0;
    vTime[i]=0;
    vTime2[i]=0;
    xlow[i]=0;
    sigmaTime[i]=0;
    sigmaToT[i]=0;
    n[i]=0;
    nentrx[i]=0;
  }
  TH2F* hSlewing = new TH2F("hSlewing", "hSlewing", knbinToT, minToT, maxToT, nbinTime, minTime, maxTime);
  TArrayF * parrToT = fArrayToT->GetArray(ich);
  TArrayF & refaToT = * parrToT;
  TArrayF * parrTime = fArrayTime->GetArray(ich);
  TArrayF & refaTime = * parrTime;
  for (Int_t j = 0; j < AliTOFcalib::fgkchannel; j++){
    if (refaToT[j] == 0) continue; 
    Int_t nx = (Int_t)((refaToT[j]-minToT)/deltaToT)+1;
    if ((refaToT[j] != 0) && (refaTime[j] != 0)){
      vTime[nx]+=refaTime[j];
      vTime2[nx]+=(refaTime[j])*(refaTime[j]);
      vToT[nx]+=refaToT[j];
      vToT2[nx]+=refaToT[j]*refaToT[j];
      nentrx[nx]++;
      hSlewing->Fill(refaToT[j],refaTime[j]);
    }
  }
  Int_t nbinsToT=hSlewing->GetNbinsX();
  if (nbinsToT != knbinToT) {
    AliError("Profile :: incompatible numbers of bins");
    return 0x0;
  }

  Int_t usefulBins=0;
  TH1F *histo = new TH1F("histo", "1D Time vs ToT", nbinsToT, minToT, maxToT);
  for (Int_t i=1;i<=nbinsToT;i++){
    if (nentrx[i]!=0){
    n[usefulBins]+=nentrx[i];
    if (n[usefulBins]==0 && i == nbinsToT) {
      break;
    }
    meanTime[usefulBins]+=vTime[i];
    meanTime2[usefulBins]+=vTime2[i];
    meanToT[usefulBins]+=vToT[i];
    meanToT2[usefulBins]+=vToT2[i];
    if (n[usefulBins]<20 && i!=nbinsToT) continue; 
    mTime[usefulBins]=meanTime[usefulBins]/n[usefulBins];
    mToT[usefulBins]=meanToT[usefulBins]/n[usefulBins];
    sigmaTime[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
				   *(meanTime2[usefulBins]-meanTime[usefulBins]
				     *meanTime[usefulBins]/n[usefulBins]));
    if ((1./n[usefulBins]/n[usefulBins]
	 *(meanToT2[usefulBins]-meanToT[usefulBins]
	   *meanToT[usefulBins]/n[usefulBins]))< 0) {
      AliError(" too small radical" );
      sigmaToT[usefulBins]=0;
    }
    else{       
      sigmaToT[usefulBins]=TMath::Sqrt(1./n[usefulBins]/n[usefulBins]
				     *(meanToT2[usefulBins]-meanToT[usefulBins]
				       *meanToT[usefulBins]/n[usefulBins]));
    }
    usefulBins++;
    }
  }
  for (Int_t i=0;i<usefulBins;i++){
    Int_t binN = (Int_t)((mToT[i]-minToT)/deltaToT)+1;
    histo->Fill(mToT[i],mTime[i]);
    histo->SetBinError(binN,sigmaTime[i]);
  } 
  return histo;
}
//_____________________________________________________________________________

void AliTOFcalib::CorrectESDTime()
{
  //Calculate the corrected TOF time
  Int_t nelements = fESDsel->GetEntries();
  for (Int_t i=0; i< nelements; i++) {
    AliTOFcalibESD *element=(AliTOFcalibESD*)fESDsel->At(i);
    Int_t index = element->GetTOFCalChannel();
    Float_t tToT = element->GetToT();
    //select the correspondent channel with its simulated ToT spectrum
    //summing up everything, index = 0 for all channels:
    Int_t ipid = element->GetCombID();
    Double_t etime = 0;   //expected time
    Double_t expTime[10]; 
    element->GetIntegratedTimes(expTime);
    if (ipid == 0) etime = expTime[2]*1E-3; //ns
    else if (ipid == 1) etime = expTime[3]*1E-3; //ns
    else if (ipid == 2) etime = expTime[4]*1E-3; //ns
    Float_t par[6];
    if (!fTOFCal) {
      AliTOFGeometry *geom=new AliTOFGeometryV5();
      fTOFCal = new AliTOFCal(geom);
      fTOFCal->CreateArray();
      delete geom;
    }
    AliTOFChannel * calChannel = fTOFCal->GetChannel(index);
    for (Int_t j = 0; j<6; j++){
      par[j]=calChannel->GetSlewPar(j);
    }
    Float_t timeCorr=0;
    timeCorr= par[0]+par[1]*tToT+par[2]*tToT*tToT+par[3]*tToT*tToT*tToT+par[4]*tToT*tToT*tToT*tToT+par[5]*tToT*tToT*tToT*tToT*tToT;
  }
}
//_____________________________________________________________________________

void AliTOFcalib::CorrectESDTime(AliESD *event){
  //Calculate the corrected TOF time

  Int_t ntrk =0;
  ntrk=event->GetNumberOfTracks();
  for (Int_t itrk=0; itrk<ntrk; itrk++) {
    AliESDtrack *t=event->GetTrack(itrk);
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) {
      continue;
    }
    //IsStartedTimeIntegral
    if ((t->GetStatus()&AliESDtrack::kTIME)==0) {
      continue;
    }
    UInt_t assignedTOFcluster=t->GetTOFcluster();//index of the assigned TOF cluster, >0 ?
    if(assignedTOFcluster==0){ // not matched
      continue;
    }
    Int_t index = t->GetTOFCalChannel();
    if (!fTOFCal) {
      AliTOFGeometry *geom=new AliTOFGeometryV5();
      fTOFCal = new AliTOFCal(geom);
      fTOFCal->CreateArray();
      delete geom;
    }
    AliTOFChannel * calChannel = fTOFCal->GetChannel(index);
    Float_t par[6];
    for (Int_t j = 0; j<6; j++){
      par[j]=calChannel->GetSlewPar(j);
    }
    Float_t tToT = t->GetTOFsignalToT();
    Float_t timeCorr =0; 
    timeCorr=par[0]+par[1]*tToT+par[2]*tToT*tToT+par[3]*tToT*tToT*tToT+par[4]*tToT*tToT*tToT*tToT+par[5]*tToT*tToT*tToT*tToT*tToT;
  }
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write calibration parameters to the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "Par" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId id(out,minrun,maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  if (!fTOFCal) {
    AliTOFGeometry *geom=new AliTOFGeometryV5();
    fTOFCal = new AliTOFCal(geom);
    fTOFCal->CreateArray();
    delete geom;
  }
  man->Put(fTOFCal,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFCal *cal){
  //Write calibration parameters to the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "Par" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId id(out,minrun,maxrun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  man->Put(cal,id,md);
  delete md;
}
//_____________________________________________________________________________

void AliTOFcalib::ReadParFromCDB(Char_t *sel, Int_t nrun)
{
  //Read calibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "Par" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry = man->Get(out,nrun);
  AliTOFCal *cal =(AliTOFCal*)entry->GetObject();
  fTOFCal = cal;
}
//_____________________________________________________________________________
void AliTOFcalib::WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Sim miscalibration parameters to the CDB


  //for the time being, only one spectrum is used
  TFile *spFile = new TFile("$ALICE_ROOT/TOF/data/spectrum.root","read");
  TH1F * hToT;
  // Retrieve ToT Spectrum
  spFile->GetObject("ToT",hToT);

  fTOFSimToT=hToT;
  
  // Retrieve Time over TOT dependence 
  
  TH1F * h = (TH1F*)spFile->Get("TimeToTFit");
  TList * list = (TList*)h->GetListOfFunctions();
  TF1* f = (TF1*)list->At(0);
  Float_t par[6] = {0,0,0,0,0,0};
  Int_t npar=f->GetNpar();
  for (Int_t ipar=0;ipar<npar;ipar++){
    par[ipar]=f->GetParameter(ipar);
  }
  if (!fTOFSimCal) {
    AliTOFGeometry *geom=new AliTOFGeometryV5();
    fTOFSimCal = new AliTOFCal(geom);
    fTOFSimCal->CreateArray();
    delete geom;
  }
  for(Int_t iTOFch=0; iTOFch<fTOFSimCal->NPads();iTOFch++){
    AliTOFChannel * calChannel = fTOFSimCal->GetChannel(iTOFch);
    calChannel->SetSlewPar(par);
  }

  // Store them in the CDB

  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  Char_t *sel1 = "SimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId id1(out,minrun,maxrun);
  man->Put(fTOFSimCal,id1,md);
  Char_t *sel2 = "SimHisto" ;
  sprintf(out,"%s/%s",sel,sel2); 
  AliCDBId id2(out,minrun,maxrun);
  man->Put(fTOFSimToT,id2,md);
  delete md;
}

//_____________________________________________________________________________
void AliTOFcalib::WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun, AliTOFCal *cal, TH1F * histo){
  //Write Sim miscalibration parameters to the CDB

  fTOFSimToT=histo;
  fTOFSimCal=cal;  
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Chiara Zampolli");
  Char_t *sel1 = "SimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId id1(out,minrun,maxrun);
  man->Put(fTOFSimCal,id1,md);
  Char_t *sel2 = "SimHisto" ;
  sprintf(out,"%s/%s",sel,sel2); 
  AliCDBId id2(out,minrun,maxrun);
  man->Put(fTOFSimToT,id2,md);
  delete md;
}
//_____________________________________________________________________________
void AliTOFcalib::ReadSimParFromCDB(Char_t *sel, Int_t nrun)
{
  //Read miscalibration parameters from the CDB
  AliCDBManager *man = AliCDBManager::Instance();
  if(!man->IsDefaultStorageSet())man->SetDefaultStorage("local://$ALICE_ROOT");
  Char_t *sel1 = "SimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry1 = man->Get(out,nrun);
  AliTOFCal *cal =(AliTOFCal*)entry1->GetObject();
  fTOFSimCal=cal;
  Char_t *sel2 = "SimHisto" ;
  sprintf(out,"%s/%s",sel,sel2); 
  AliCDBEntry *entry2 = man->Get(out,nrun);
  TH1F *histo =(TH1F*)entry2->GetObject();
  fTOFSimToT=histo;
}
//_____________________________________________________________________________

Int_t AliTOFcalib::GetIndex(Int_t *detId)
{
  //Retrieve calibration channel index 
  Int_t isector = detId[0];
  if (isector >= fNSector)
    AliError(Form("Wrong sector number in TOF (%d) !",isector));
  Int_t iplate = detId[1];
  if (iplate >= fNPlate)
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
  Int_t istrip = detId[2];
  Int_t ipadz = detId[3];
  Int_t ipadx = detId[4];
  Int_t stripOffset = 0;
  switch (iplate) {
  case 0:
    stripOffset = 0;
    break;
  case 1:
    stripOffset = fNStripC;
    break;
  case 2:
    stripOffset = fNStripC+fNStripB;
    break;
  case 3:
    stripOffset = fNStripC+fNStripB+fNStripA;
    break;
  case 4:
    stripOffset = fNStripC+fNStripB+fNStripA+fNStripB;
    break;
  default:
    AliError(Form("Wrong plate number in TOF (%d) !",iplate));
    break;
  };

  Int_t idet = ((2*(fNStripC+fNStripB)+fNStripA)*fNpadZ*fNpadX)*isector +
               (stripOffset*fNpadZ*fNpadX)+
               (fNpadZ*fNpadX)*istrip+
	       (fNpadX)*ipadz+
	        ipadx;
  return idet;
}

