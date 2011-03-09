/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////
//            Implementation of the ITS clusterer V2 class                //
//                                                                        //
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch            //
//          Last revision: 13-05-09 Enrico Fragiacomo                     //
//                                  enrico.fragiacomo@ts.infn.it          //
//                                                                        //
///////////////////////////////////////////////////////////////////////////

#include "AliITSClusterFinderV2SSD.h"

#include <Riostream.h>
#include <TGeoGlobalMagField.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeomTGeo.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSSD.h"
#include <TClonesArray.h>
#include <TCollection.h>
#include "AliITSdigitSSD.h"
#include "AliITSReconstructor.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSSD.h"

Short_t *AliITSClusterFinderV2SSD::fgPairs = 0x0;
Int_t    AliITSClusterFinderV2SSD::fgPairsSize = 0;
const Float_t  AliITSClusterFinderV2SSD::fgkThreshold = 5.;

const Float_t AliITSClusterFinderV2SSD::fgkCosmic2008StripShifts[16][9] = 
  {{-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35},  // DDL 512
   {-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35,-0.35},  // DDL 513
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15},  // DDL 514
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15},  // DDL 515
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},  // DDL 516
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},  // DDL 517
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15},  // DDL 518
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15},  // DDL 519
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.25,-0.15},  // DDL 520
   {-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15,-0.15},  // DDL 521
   {-0.10,-0.10,-0.10,-0.40,-0.40,-0.40,-0.10,-0.10,-0.45},  // DDL 522
   {-0.10,-0.10,-0.10,-0.35,-0.35,-0.35,-0.10,-0.35,-0.50},  // DDL 523
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},  // DDL 524
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},  // DDL 525
   { 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35},  // DDL 526
   { 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45}}; // DDL 527

ClassImp(AliITSClusterFinderV2SSD)


  AliITSClusterFinderV2SSD::AliITSClusterFinderV2SSD(AliITSDetTypeRec* dettyp):AliITSClusterFinder(dettyp),fLastSSD1(AliITSgeomTGeo::GetModuleIndex(6,1,1)-1), fLorentzShiftP(0), fLorentzShiftN(0)
{
//Default constructor
  static AliITSRecoParam *repa = NULL;  
  if(!repa){
    repa = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
    if(!repa){
      repa = AliITSRecoParam::GetHighFluxParam();
      AliWarning("Using default AliITSRecoParam class");
    }
  }

  if (repa->GetCorrectLorentzAngleSSD()) {
    AliMagF* field = dynamic_cast<AliMagF*>(TGeoGlobalMagField::Instance()->GetField()); 
    if (field == 0) {
      AliError("Cannot get magnetic field from TGeoGlobalMagField");
    }
    else {
      Float_t Bfield = field->SolenoidField();
      // NB: spatial shift has opposite sign for lay 5 and 6, but strip numbering also changes direction, so no sign-change 
      // Shift due to ExB on drift N-side, units: strip width 
      fLorentzShiftP = -repa->GetTanLorentzAngleElectronsSSD() * 150.e-4/95.e-4 * Bfield / 5.0;
      // Shift due to ExB on drift P-side, units: strip width 
      fLorentzShiftN = -repa->GetTanLorentzAngleHolesSSD() * 150.e-4/95.e-4 * Bfield / 5.0;
      AliDebug(1,Form("Bfield %f Lorentz Shift P-side %f N-side %f",Bfield,fLorentzShiftN,fLorentzShiftP));
    }
  }
}
 
//______________________________________________________________________
AliITSClusterFinderV2SSD::AliITSClusterFinderV2SSD(const AliITSClusterFinderV2SSD &cf) : AliITSClusterFinder(cf), fLastSSD1(cf.fLastSSD1), fLorentzShiftP(cf.fLorentzShiftP), fLorentzShiftN(cf.fLorentzShiftN)
{
  // Copy constructor
}

//______________________________________________________________________
AliITSClusterFinderV2SSD& AliITSClusterFinderV2SSD::operator=(const AliITSClusterFinderV2SSD&  cf ){
  // Assignment operator

  this->~AliITSClusterFinderV2SSD();
  new(this) AliITSClusterFinderV2SSD(cf);
  return *this;
}


void AliITSClusterFinderV2SSD::FindRawClusters(Int_t mod){

  //Find clusters V2
  SetModule(mod);
  FindClustersSSD(fDigits);

}

void AliITSClusterFinderV2SSD::FindClustersSSD(TClonesArray *alldigits) {
  //------------------------------------------------------------
  // Actual SSD cluster finder
  //------------------------------------------------------------
  Int_t smaxall=alldigits->GetEntriesFast();
  if (smaxall==0) return;


  //---------------------------------------
  // load recoparam and calibration
  // 
  static AliITSRecoParam *repa = NULL;  
  if(!repa){
    repa = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
    if(!repa){
      repa = AliITSRecoParam::GetHighFluxParam();
      AliWarning("Using default AliITSRecoParam class");
    }
  }

  AliITSCalibrationSSD* cal = (AliITSCalibrationSSD*)GetResp(fModule);
  Float_t gain=0;
  Float_t noise=0;
  //---------------------------------------


  //------------------------------------
  // fill the digits array with zero-suppression condition
  // Signal is converted in KeV
  //
  TObjArray digits;
  for (Int_t i=0;i<smaxall; i++){
    AliITSdigitSSD *d=(AliITSdigitSSD*)alldigits->UncheckedAt(i);

    if(d->IsSideP()) noise = cal->GetNoiseP(d->GetStripNumber());  
    else noise = cal->GetNoiseN(d->GetStripNumber());
    if (d->GetSignal()<3.*noise) continue;

    if(d->IsSideP()) gain = cal->GetGainP(d->GetStripNumber());  
    else gain = cal->GetGainN(d->GetStripNumber());

    Float_t q=gain*d->GetSignal(); //
    q=cal->ADCToKeV(q); // converts the charge in KeV from ADC units
    d->SetSignal(Int_t(q));

    digits.AddLast(d);
  }
  Int_t smax = digits.GetEntriesFast();
  if (smax==0) return;
  //------------------------------------

  
  const Int_t kMax=1000;
  Int_t np=0, nn=0; 
  Ali1Dcluster pos[kMax], neg[kMax];
  Float_t y=0., q=0., qmax=0.; 
  Int_t lab[4]={-2,-2,-2,-2};  
  Bool_t flag5 = 0;
  
  /*
  cout<<"-----------------------------"<<endl;
  cout<<"this is module "<<fModule;
  cout<<endl;
  cout<<endl;
  */
  Int_t layer = 4;
  if (fModule>fLastSSD1) 
    layer = 5;

  //--------------------------------------------------------
  // start 1D-clustering from the first digit in the digits array
  //
  AliITSdigitSSD *d=(AliITSdigitSSD*)digits.UncheckedAt(0);
  q += d->GetSignal();
  y += d->GetCoord2()*d->GetSignal();
  qmax=d->GetSignal();
  lab[0]=d->GetTrack(0); lab[1]=d->GetTrack(1); lab[2]=d->GetTrack(2);
  
  if(d->IsSideP()) {
    noise = cal->GetNoiseP(d->GetStripNumber());  
    gain = cal->GetGainP(d->GetStripNumber());
  }
  else {
    noise = cal->GetNoiseN(d->GetStripNumber());
    gain = cal->GetGainN(d->GetStripNumber());
  }  
  noise*=gain;
  noise=cal->ADCToKeV(noise); // converts noise in KeV from ADC units

  if(qmax>fgkThreshold*noise) flag5=1; // seed for the cluster

  /*
  cout<<d->GetSignal()<<" "<<noise<<" "<<flag5<<" "<<
    d->GetCoord1()<<" "<<d->GetCoord2()<<endl;
  */

  Int_t curr=d->GetCoord2();
  Int_t flag=d->GetCoord1();

  // Note: the first side which will be processed is supposed to be the 
  // P-side which is neg
  Int_t *n=&nn;
  Ali1Dcluster *c=neg;
  if(flag) {n=&np; c=pos;} // in case we have only Nstrips (P was bad!)

  Int_t nd=1;
  Int_t milab[10];
  for (Int_t ilab=0;ilab<10;ilab++){
    milab[ilab]=-2;
  }
  milab[0]=d->GetTrack(0); milab[1]=d->GetTrack(1); milab[2]=d->GetTrack(2);


  //----------------------------------------------------------
  // search for neighboring digits
  //
  for (Int_t s=1; s<smax; s++) {
      d=(AliITSdigitSSD*)digits.UncheckedAt(s);      
      Int_t strip=d->GetCoord2();

      // if digits is not a neighbour or side did not change 
      // and at least one of the previous digits met the seed condition
      // then creates a new 1D cluster
      if ( ( ((strip-curr) > 1) || (flag!=d->GetCoord1()) ) ) {

	if(flag5) {
	  //cout<<"here1"<<endl;
	  Float_t dLorentz = 0;
	  if (!flag) { // P-side is neg clust
	    dLorentz = fLorentzShiftN;
	  }
	  else { // N-side is p clust
	    dLorentz = fLorentzShiftP;
	  }
         c[*n].SetY(y/q+dLorentz);
         c[*n].SetQ(q);
         c[*n].SetNd(nd);
	 CheckLabels2(milab);
         c[*n].SetLabels(milab);

	 if(repa->GetUseUnfoldingInClusterFinderSSD()==kTRUE) {
	   // Note: fUseUnfoldingInClusterFinderSSD=kFALSE by default in RecoParam

	   //Split suspiciously big cluster
	   if (nd>4&&nd<25) {
	     c[*n].SetY(y/q-0.25*nd+dLorentz);
	     c[*n].SetQ(0.5*q);
	     (*n)++;
	     if (*n==kMax) {
	       Error("FindClustersSSD","Too many 1D clusters !");
	       return;
	     }
	     c[*n].SetY(y/q+0.25*nd+dLorentz);
	     c[*n].SetQ(0.5*q);
	     c[*n].SetNd(nd);
	     c[*n].SetLabels(milab);
	   }	 
	   
	 } // unfolding is on

         (*n)++;
         if (*n==kMax) {
          Error("FindClustersSSD","Too many 1D clusters !");
          return;
         }

	} // flag5 set

	 // reset everything
         y=q=qmax=0.;
         nd=0;
	 flag5=0;
         lab[0]=lab[1]=lab[2]=-2;
	 for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;

	 // if side changed from P to N, switch to pos 1D clusters
	 // (if for some reason the side changed from N to P then do the opposite)
         if (flag!=d->GetCoord1()) 
	   { if(!flag) {n=&np; c=pos;} else {n=&nn; c=neg;} }

      } // end create new 1D cluster from previous neighboring digits

      // continues adding digits to the previous cluster 
      // or start a new one 
      flag=d->GetCoord1();
      q += d->GetSignal();
      y += d->GetCoord2()*d->GetSignal();
      nd++;

      if(d->IsSideP()) {
	noise = cal->GetNoiseP(d->GetStripNumber());  
	gain = cal->GetGainP(d->GetStripNumber());
      }
      else {
	noise = cal->GetNoiseN(d->GetStripNumber());
	gain = cal->GetGainN(d->GetStripNumber());
      }
      noise*=gain;
      noise=cal->ADCToKeV(noise); // converts the charge in KeV from ADC units

      if(d->GetSignal()>fgkThreshold*noise) flag5=1;

      /*
  cout<<d->GetSignal()<<" "<<noise<<" "<<flag5<<" "<<
    d->GetCoord1()<<" "<<d->GetCoord2()<<endl;
      */

      if (d->GetSignal()>qmax) {
         qmax=d->GetSignal();
         lab[0]=d->GetTrack(0); lab[1]=d->GetTrack(1); lab[2]=d->GetTrack(2);
      }
      for (Int_t ilab=0;ilab<10;ilab++) {
	if (d->GetTrack(ilab)>=0) AddLabel(milab, (d->GetTrack(ilab))); 
      }
      curr=strip;


  } // loop over digits, no more digits in the digits array


  // add the last 1D cluster 
  if(flag5) {

    //	cout<<"here2"<<endl;
    Float_t dLorentz = 0;
    if (!flag) { // P-side is neg clust
      dLorentz = fLorentzShiftN;
    }
    else { // N-side is p clust
      dLorentz = fLorentzShiftP;
    }
    
    c[*n].SetY(y/q + dLorentz);
    c[*n].SetQ(q);
    c[*n].SetNd(nd);
    c[*n].SetLabels(lab);
    
    if(repa->GetUseUnfoldingInClusterFinderSSD()==kTRUE) {
      
      //Split suspiciously big cluster
      if (nd>4 && nd<25) {
	c[*n].SetY(y/q-0.25*nd + dLorentz);
	c[*n].SetQ(0.5*q);
	(*n)++;
	if (*n==kMax) {
	  Error("FindClustersSSD","Too many 1D clusters !");
	  return;
	}
	c[*n].SetY(y/q+0.25*nd + dLorentz);
	c[*n].SetQ(0.5*q);
	c[*n].SetNd(nd);
	c[*n].SetLabels(lab);
      }
    } // unfolding is on
    
    (*n)++;
    if (*n==kMax) {
      Error("FindClustersSSD","Too many 1D clusters !");
      return;
    }

  } // if flag5 last 1D cluster added


  //------------------------------------------------------
  // call FindClustersSSD to pair neg and pos 1D clusters 
  // and create recpoints from the crosses
  // Note1: neg are Pside and pos are Nside!!
  // Note2: if there are no Pside digits nn=0 (bad strips!!) (same for Nside)
  //
  //  cout<<nn<<" Pside and "<<np<<" Nside clusters"<<endl;
  
  AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();   
  if (nn*np > 0) { 
    TClonesArray* clusters = rpc->UncheckedGetClusters(fModule);
    clusters->Clear();
    FindClustersSSD(neg, nn, pos, np, clusters);
    TIter itr(clusters);
    AliITSRecPoint *irp;
    while ((irp = (AliITSRecPoint*)itr.Next()))  fDetTypeRec->AddRecPoint(*irp);
  }
  //-----------------------------------------------------
}


void AliITSClusterFinderV2SSD::RawdataToClusters(AliRawReader* rawReader){

  //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStreamSSD inputSSD(rawReader);
  FindClustersSSD(&inputSSD);
  
}


void AliITSClusterFinderV2SSD::FindClustersSSD(AliITSRawStreamSSD* input) 
{
  //------------------------------------------------------------
  // Actual SSD cluster finder for raw data
  //------------------------------------------------------------

  AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();
  static AliITSRecoParam *repa = NULL;
  if(!repa){
    repa = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
    if(!repa){
      repa = AliITSRecoParam::GetHighFluxParam();
      AliWarning("Using default AliITSRecoParam class");
    }
  }
  Int_t nClustersSSD = 0;
  const Int_t kNADC = 12;
  const Int_t kMaxADCClusters = 1000;

  Int_t strips[kNADC][2][kMaxADCClusters][2]; // [ADC],[side],[istrip], [0]=istrip [1]=signal
  Int_t nStrips[kNADC][2];

  for( int i=0; i<kNADC; i++ ){
    nStrips[i][0] = 0;
    nStrips[i][1] = 0;
  }

  Int_t ddl = -1;
  Int_t ad = -1;
  
  //*
  //* Loop over modules DDL+AD
  //*
  
  while (kTRUE) {

    bool next = input->Next();
    
    //* 
    //* Continue if corrupted input
    //*

    if( (!next)&&(input->flag) ){
     AliWarning(Form("ITSClustersFinderSSD: Corrupted data: warning from RawReader"));
      continue; 
    }

    Int_t newDDL = input->GetDDL(); 
    Int_t newAD = input->GetAD();

    if( next ){
      if( newDDL<0 || newDDL>15 ){
	AliWarning(Form("ITSClustersFinderSSD: Corrupted data: wrong DDL number (%d)",newDDL));
	continue;
      }
      
      if( newAD<1 || newAD>9 ){
	AliWarning(Form("ITSClustersFinderSSD: Corrupted data: wrong AD number (%d)",newAD));
	continue;
      }
    }

    bool newModule = ( !next || ddl!= newDDL || ad!=newAD );

    if( newModule && ddl>=0 && ad>=0 ){ 

      //*
      //* Reconstruct the previous block of 12 modules --- actual clusterfinder
      //* 
      //cout<<endl;
      for( int adc = 0; adc<kNADC; adc++ ){
	
	//* 1D clusterfinder

	Ali1Dcluster clusters1D[2][kMaxADCClusters]; // per ADC, per side
	Int_t nClusters1D[2] = {0,0};
	//int nstat[2] = {0,0};
	fModule = AliITSRawStreamSSD::GetModuleNumber(ddl, (ad - 1)  * 12 + adc );
	
	if( fModule<0 ){
//	  AliWarning(Form("ITSClustersFinderSSD: Corrupted data: module (ddl %d ad %d adc %d) not found in the map",ddl,ad,adc));
//CM channels are always present even everything is suppressed 
	  continue;
	}
	
	Int_t layer = 4;
	if (fModule>fLastSSD1) 
	  layer = 5;

	AliITSCalibrationSSD* cal = (AliITSCalibrationSSD*)fDetTypeRec->GetCalibrationModel(fModule);
	if( !cal ){
	  AliWarning(Form("ITSClustersFinderSSD: No calibration found for module (ddl %d ad %d adc %d)",ddl,ad,adc));	    
	  continue;
	}

	Float_t dStrip = 0;

	if( repa->GetUseCosmicRunShiftsSSD()) {  // Special condition for 2007/2008 cosmic data
	  dStrip = fgkCosmic2008StripShifts[ddl][ad-1];
	  if (TMath::Abs(dStrip) > 1.5){
	    AliWarning(Form("Indexing error in Cosmic calibration: ddl = %d, dStrip %f\n",ddl,dStrip));
	    dStrip = 0;
	  }	
	}
	
	for( int side=0; side<=1; side++ ){

	  Int_t lab[3]={-2,-2,-2};
	  Float_t q = 0.;
	  Float_t y = 0.;
	  Int_t nDigits = 0;
	  Int_t ostrip = -2;
	  Bool_t snFlag = 0;

	  Float_t dLorentz = 0;
	  if (side==0) { // P-side is neg clust
	    dLorentz = fLorentzShiftN;
	  }
	  else { // N-side is pos clust
	    dLorentz = fLorentzShiftP;
	  }
	  
	  Int_t n = nStrips[adc][side];
	  for( int istr = 0; istr<n+1; istr++ ){
	    
	    bool stripOK = 1;
	    Int_t strip=0;
	    Float_t signal=0.0, noise=0.0, gain=0.0;
	    
	    if( istr<n ){
	      strip = strips[adc][side][istr][0];
	      signal = strips[adc][side][istr][1];
	      
	      //cout<<"strip "<<adc<<" / "<<side<<": "<<strip<<endl;

	      if( cal ){
		noise = side ?cal->GetNoiseN(strip) :cal->GetNoiseP(strip); 
		gain = side ?cal->GetGainN(strip) :cal->GetGainP(strip);	 
		stripOK = ( noise>=1. && signal>=3.0*noise 
			    //&& !cal->IsPChannelBad(strip) 
			    );
	      }
	    } else stripOK = 0; // end of data

	    bool newCluster = ( (abs(strip-ostrip)!=1) || !stripOK );	  
	        
	    if( newCluster ){

	      //* Store the previous cluster

	      if( nDigits>0 && q>0 && snFlag ){

		if (nClusters1D[side] >= kMaxADCClusters-1 ) {
		  AliWarning("HLT ClustersFinderSSD: Too many 1D clusters !");
		}else {
		  
		  Ali1Dcluster &cluster = clusters1D[side][nClusters1D[side]++];
		  cluster.SetY( y / q + dStrip + dLorentz);
		  cluster.SetQ(q);
		  cluster.SetNd(nDigits);
		  cluster.SetLabels(lab);
		  //cout<<"cluster 1D side "<<side<<": y= "<<y<<" q= "<<q<<" d="<<dStrip<<" Y="<<cluster.GetY()<<endl;
		  //Split suspiciously big cluster

		  if( repa->GetUseUnfoldingInClusterFinderSSD()
		      && nDigits > 4 && nDigits < 25 
		      ){
		    cluster.SetY(y/q + dStrip - 0.25*nDigits + dLorentz);	    
		    cluster.SetQ(0.5*q);	  
		    Ali1Dcluster& cluster2 = clusters1D[side][nClusters1D[side]++];
		    cluster2.SetY(y/q + dStrip + 0.25*nDigits + dLorentz);	    
		    cluster2.SetQ(0.5*q);
		    cluster2.SetNd(nDigits);
		    cluster2.SetLabels(lab);	  
		  } // unfolding is on	  	
		}
	      }
	      y = q = 0.;
	      nDigits = 0;
	      snFlag = 0;

	    } //* End store the previous cluster

	    if( stripOK ){ // add new signal to the cluster
//	      signal = (Int_t) (signal * gain); // signal is corrected for gain
	      if( signal>fgkThreshold*noise) snFlag = 1;
	      signal = signal * gain; // signal is corrected for gain
//	      if( cal ) signal = (Int_t) cal->ADCToKeV( signal ); // signal is  converted in KeV  	  
	      if( cal ) signal = cal->ADCToKeV( signal ); // signal is  converted in KeV  	  
	      q += signal;	  // add digit to current cluster
	      y += strip * signal;	  
	      nDigits++;
	      //nstat[side]++;
	      ostrip = strip;

	    }
	  } //* end loop over strips

	} //* end loop over ADC sides
	

  	//* 2D clusterfinder
      	if( nClusters1D[0] && nClusters1D[1] && fModule>=0 ){
           TClonesArray* clusters = rpc->UncheckedGetClusters(fModule);
	       FindClustersSSD( clusters1D[0], nClusters1D[0], clusters1D[1], nClusters1D[1], clusters); 
           Int_t nClustersn = clusters->GetEntriesFast();
	       nClustersSSD += nClustersn;
	}

	//cout<<"SG: "<<ddl<<" "<<ad<<" "<<adc<<": strips "<<nstat[0]<<"+"<<nstat[1]<<", clusters 1D= "<<nClusters1D[0]<<" + "<<nClusters1D[1]<<", 2D= "<<clusters.size()<<endl;

      }//* end loop over adc
      
    }//* end of reconstruction of previous block of 12 modules
    
    if( newModule ){
      
      //*
      //* Clean up arrays and set new module
      //* 
      
      for( int i=0; i<kNADC; i++ ){
	nStrips[i][0] = 0;
	nStrips[i][1] = 0;
      }     
      ddl = newDDL;
      ad = newAD;
    } 
    

    //*
    //* Exit main loop when there is no more input
    //* 

    if( !next ) break;
    
    //* 
    //* Fill the current strip information
    //* 

    Int_t adc = input->GetADC(); 
    if( adc<0 || adc>=kNADC+2 || (adc>5&&adc<8) ){
      AliWarning(Form("HLT ClustersFinderSSD: Corrupted data: wrong adc number (%d)", adc));
      continue;
    }

    if( adc>7 ) adc-= 2; // shift ADC numbers 8-13 to 6-11
    
    Bool_t side = input->GetSideFlag();
    Int_t strip = input->GetStrip();
    Int_t signal = input->GetSignal();
    

    //cout<<"SSD: "<<ddl<<" "<<ad<<" "<<adc<<" "<<side<<" "<<strip<<" : "<<signal<<endl;

    if( strip>767 ){    
      AliWarning(Form("HLT ClustersFinderSSD: Corrupted data: wrong strip number (ddl %d ad %d adc %d side %d, strip %d", 
		      ddl, ad, adc, side,strip) );
      continue;
    }
    if (strip < 0) continue;
    
    int &n = nStrips[adc][side];
    if( n >0 ){
      Int_t oldStrip = strips[adc][side][n-1][0];

      if( strip==oldStrip ){
	AliWarning(Form("HLT ClustersFinderSSD: Corrupted data: duplicated signal: ddl %d ad %d adc %d, side %d, strip %d", 
			ddl, ad, adc, side, strip ));
	continue;
      }
    }
    strips[adc][side][n][0] = strip;
    strips[adc][side][n][1] = signal;    
    n++;

    //cout<<"SSD: "<<input->GetDDL()<<" "<<input->GetAD()<<" "
    //<<input->GetADC()<<" "<<input->GetSideFlag()<<" "<<((int)input->GetStrip())<<" "<<strip<<" : "<<input->GetSignal()<<endl;

  } //* End main loop over the input
  
  AliDebug(1,Form("found clusters in ITS SSD: %d", nClustersSSD));
}


void AliITSClusterFinderV2SSD::
FindClustersSSD(const Ali1Dcluster* neg, Int_t nn, 
		const Ali1Dcluster* pos, Int_t np,
		const TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SSD cluster finder
  //------------------------------------------------------------

  const TGeoHMatrix *mT2L=AliITSgeomTGeo::GetTracking2LocalMatrix(fModule);

  //---------------------------------------
  // load recoparam
  // 
  static AliITSRecoParam *repa = NULL;  
  if(!repa){
    repa = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
    if(!repa){
      repa = AliITSRecoParam::GetHighFluxParam();
      AliWarning("Using default AliITSRecoParam class");
    }
  }

//  TClonesArray &cl=*clusters;
  
  AliITSsegmentationSSD *seg = static_cast<AliITSsegmentationSSD*>(fDetTypeRec->GetSegmentationModel(2));
  if (fModule>fLastSSD1) 
    seg->SetLayer(6);
  else 
    seg->SetLayer(5);

  Float_t hwSSD = seg->Dx()*1e-4/2;
  Float_t hlSSD = seg->Dz()*1e-4/2;

  Int_t idet=fNdet[fModule];
  Int_t ncl=0;

  //
  Int_t *cnegative = new Int_t[np];  
  Int_t *cused1 = new Int_t[np];
  Int_t *negativepair = new Int_t[10*np];
  Int_t *cpositive = new Int_t[nn];  
  Int_t *cused2 = new Int_t[nn];  
  Int_t *positivepair = new Int_t[10*nn];  
  for (Int_t i=0;i<np;i++) {cnegative[i]=0; cused1[i]=0;}
  for (Int_t i=0;i<nn;i++) {cpositive[i]=0; cused2[i]=0;}
  for (Int_t i=0;i<10*np;i++) {negativepair[i]=0;}
  for (Int_t i=0;i<10*nn;i++) {positivepair[i]=0;}

  if ((np*nn) > fgPairsSize) {

    if (fgPairs) delete [] fgPairs;
    fgPairsSize = 4*np*nn;
    fgPairs = new Short_t[fgPairsSize];
  }
  memset(fgPairs,0,sizeof(Short_t)*np*nn);

  //
  // find available pairs
  //
  Int_t ncross = 0;
  for (Int_t i=0; i<np; i++) {
    Float_t yp=pos[i].GetY(); 
    if ( (pos[i].GetQ()>0) && (pos[i].GetQ()<3) ) continue;
    for (Int_t j=0; j<nn; j++) {
      if ( (neg[j].GetQ()>0) && (neg[j].GetQ()<3) ) continue;
      Float_t yn=neg[j].GetY();

      Float_t xt, zt;
      seg->GetPadCxz(yn, yp, xt, zt);
      //cout<<yn<<" "<<yp<<" "<<xt<<" "<<zt<<endl;
      
      if (TMath::Abs(xt)<hwSSD)
      if (TMath::Abs(zt)<hlSSD) {
	Int_t in = i*10+cnegative[i];
	Int_t ip = j*10+cpositive[j];
	if ((in < 10*np) && (ip < 10*nn)) {
	  negativepair[in] =j;  //index
	  positivepair[ip] =i;
	  cnegative[i]++;  //counters
	  cpositive[j]++;
	  ncross++;	
	  fgPairs[i*nn+j]=100;
	}
	else
	  AliError(Form("Index out of range: ip=%d, in=%d",ip,in));
      }
    }
  }

  if (!ncross) {
    delete [] cnegative;
    delete [] cused1;
    delete [] negativepair;
    delete [] cpositive;
    delete [] cused2;
    delete [] positivepair;
    return;
  }
//why not to allocate memorey here?  if(!clusters) clusters = new TClonesArray("AliITSRecPoint", ncross);
  
  /* //
  // try to recover points out of but close to the module boundaries 
  //
  for (Int_t i=0; i<np; i++) {
    Float_t yp=pos[i].GetY(); 
    if ( (pos[i].GetQ()>0) && (pos[i].GetQ()<3) ) continue;
    for (Int_t j=0; j<nn; j++) {
      if ( (neg[j].GetQ()>0) && (neg[j].GetQ()<3) ) continue;
      // if both 1Dclusters have an other cross continue
      if (cpositive[j]&&cnegative[i]) continue;
      Float_t yn=neg[j].GetY();
      
      Float_t xt, zt;
      seg->GetPadCxz(yn, yp, xt, zt);
      
      if (TMath::Abs(xt)<hwSSD+0.1)
      if (TMath::Abs(zt)<hlSSD+0.15) {
	// tag 1Dcluster (eventually will produce low quality recpoint)
	if (cnegative[i]==0) pos[i].SetNd(100);  // not available pair
	if (cpositive[j]==0) neg[j].SetNd(100);  // not available pair
	Int_t in = i*10+cnegative[i];
	Int_t ip = j*10+cpositive[j];
	if ((in < 10*np) && (ip < 10*nn)) {
	  negativepair[in] =j;  //index
	  positivepair[ip] =i;
	  cnegative[i]++;  //counters
	  cpositive[j]++;	
	  fgPairs[i*nn+j]=100;
	}
	else
	  AliError(Form("Index out of range: ip=%d, in=%d",ip,in));
      }
    }
  }
  */

  //
  Float_t lp[6];
  Int_t milab[10];
  Double_t ratio;
  

  if(repa->GetUseChargeMatchingInClusterFinderSSD()==kTRUE) {


    //
    // sign gold tracks
    //
    for (Int_t ip=0;ip<np;ip++){
      Float_t xbest=1000,zbest=1000,qbest=0;
      //
      // select gold clusters
      if ( (cnegative[ip]==1) && cpositive[negativepair[10*ip]]==1){ 
	Float_t yp=pos[ip].GetY(); 
	Int_t j = negativepair[10*ip];      

	if( (pos[ip].GetQ()==0) && (neg[j].GetQ() ==0) ) { 
	  // both bad, hence continue;	  
	  // mark both as used (to avoid recover at the end)
	  cused1[ip]++; 
	  cused2[j]++;
	  continue;
	}

	ratio = (pos[ip].GetQ()-neg[j].GetQ())/(pos[ip].GetQ()+neg[j].GetQ());
	//cout<<"ratio="<<ratio<<endl;

	// charge matching (note that if posQ or negQ is 0 -> ratio=1 and the following condition is met
	if (TMath::Abs(ratio)>0.2) continue; // note: 0.2=3xsigma_ratio calculated in cosmics tests

	//
	Float_t yn=neg[j].GetY();
	
	Float_t xt, zt;
	seg->GetPadCxz(yn, yp, xt, zt);
	
	xbest=xt; zbest=zt; 

	
	qbest=0.5*(pos[ip].GetQ()+neg[j].GetQ());
	if( (pos[ip].GetQ()==0)||(neg[j].GetQ()==0)) qbest*=2; // in case of bad strips on one side keep all charge from the other one
	
	{
	  Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  lp[0]=trk[1];
	  lp[1]=trk[2];
	}
	lp[4]=qbest;        //Q
	for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	for (Int_t ilab=0;ilab<3;ilab++){
	  milab[ilab] = pos[ip].GetLabel(ilab);
	  milab[ilab+3] = neg[j].GetLabel(ilab);
	}
	//
	CheckLabels2(milab);
	milab[3]=(((ip<<10) + j)<<10) + idet; // pos|neg|det
	Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	AliITSRecPoint * cl2;
    cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);
	  
    cl2->SetChargeRatio(ratio);    	
    cl2->SetType(1);
    fgPairs[ip*nn+j]=1;

    if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
      cl2->SetType(2);
      fgPairs[ip*nn+j]=2;
    }

    if(pos[ip].GetQ()==0) cl2->SetType(3);
    if(neg[j].GetQ()==0) cl2->SetType(4);

    cused1[ip]++;
    cused2[j]++;

  	ncl++;
      }
    }
    
    for (Int_t ip=0;ip<np;ip++){
      Float_t xbest=1000,zbest=1000,qbest=0;
      //
      //
      // select "silber" cluster
      if ( cnegative[ip]==1 && cpositive[negativepair[10*ip]]==2){
	Int_t in  = negativepair[10*ip];
	Int_t ip2 = positivepair[10*in];
	if (ip2==ip) ip2 =  positivepair[10*in+1];
	Float_t pcharge = pos[ip].GetQ()+pos[ip2].GetQ();
	


	ratio = (pcharge-neg[in].GetQ())/(pcharge+neg[in].GetQ());
	if ( (TMath::Abs(ratio)<0.2) && (pcharge!=0) ) {
	  //if ( (TMath::Abs(pcharge-neg[in].GetQ())<30) && (pcharge!=0) ) { // 
	  
	  //
	  // add first pair
	  if ( (fgPairs[ip*nn+in]==100)&&(pos[ip].GetQ() ) ) {  //
	    
	    Float_t yp=pos[ip].GetY(); 
	    Float_t yn=neg[in].GetY();
	    
	    Float_t xt, zt;
	    seg->GetPadCxz(yn, yp, xt, zt);
	    
	    xbest=xt; zbest=zt; 

	    qbest =pos[ip].GetQ();
	    Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	    mT2L->MasterToLocal(loc,trk);
	    lp[0]=trk[1];
	    lp[1]=trk[2];
	    
	    lp[4]=qbest;        //Q
	    for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	    for (Int_t ilab=0;ilab<3;ilab++){
	      milab[ilab] = pos[ip].GetLabel(ilab);
	      milab[ilab+3] = neg[in].GetLabel(ilab);
	    }
	    //
	    CheckLabels2(milab);
	    ratio = (pos[ip].GetQ()-neg[in].GetQ())/(pos[ip].GetQ()+neg[in].GetQ());
	    milab[3]=(((ip<<10) + in)<<10) + idet; // pos|neg|det
	    Int_t info[3] = {pos[ip].GetNd(),neg[in].GetNd(),fNlayer[fModule]};
	    
	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	AliITSRecPoint * cl2;
	      cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);
	      cl2->SetChargeRatio(ratio);    	
	      cl2->SetType(5);
	      fgPairs[ip*nn+in] = 5;
	      if ((pos[ip].GetNd()+neg[in].GetNd())>6){ //multi cluster
		cl2->SetType(6);
		fgPairs[ip*nn+in] = 6;
	      }	    
	    ncl++;
	  }
	  
	  
	  //
	  // add second pair
	  
	  //	if (!(cused1[ip2] || cused2[in])){  //
	  if ( (fgPairs[ip2*nn+in]==100) && (pos[ip2].GetQ()) ) {
	    
	    Float_t yp=pos[ip2].GetY();
	    Float_t yn=neg[in].GetY();
	    
	    Float_t xt, zt;
	    seg->GetPadCxz(yn, yp, xt, zt);
	    
	    xbest=xt; zbest=zt; 

	    qbest =pos[ip2].GetQ();
	    
	    Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	    mT2L->MasterToLocal(loc,trk);
	    lp[0]=trk[1];
	    lp[1]=trk[2];
	    
	    lp[4]=qbest;        //Q
	    for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	    for (Int_t ilab=0;ilab<3;ilab++){
	      milab[ilab] = pos[ip2].GetLabel(ilab);
	      milab[ilab+3] = neg[in].GetLabel(ilab);
	    }
	    //
	    CheckLabels2(milab);
	    ratio = (pos[ip2].GetQ()-neg[in].GetQ())/(pos[ip2].GetQ()+neg[in].GetQ());
	    milab[3]=(((ip2<<10) + in)<<10) + idet; // pos|neg|det
	    Int_t info[3] = {pos[ip2].GetNd(),neg[in].GetNd(),fNlayer[fModule]};
	    
	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	    AliITSRecPoint * cl2;
	      cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);
	      
	      cl2->SetChargeRatio(ratio);    	
	      cl2->SetType(5);
	      fgPairs[ip2*nn+in] =5;
	      if ((pos[ip2].GetNd()+neg[in].GetNd())>6){ //multi cluster
		cl2->SetType(6);
		fgPairs[ip2*nn+in] =6;
	      }
	    ncl++;
	  }
	  
	  cused1[ip]++;
	  cused1[ip2]++;
	  cused2[in]++;

	} // charge matching condition

      } // 2 Pside cross 1 Nside
    } // loop over Pside clusters
    
    
      
      //  
    for (Int_t jn=0;jn<nn;jn++){
      if (cused2[jn]) continue;
      Float_t xbest=1000,zbest=1000,qbest=0;
      // select "silber" cluster
      if ( cpositive[jn]==1 && cnegative[positivepair[10*jn]]==2){
	Int_t ip  = positivepair[10*jn];
	Int_t jn2 = negativepair[10*ip];
	if (jn2==jn) jn2 =  negativepair[10*ip+1];
	Float_t pcharge = neg[jn].GetQ()+neg[jn2].GetQ();
	//
	

	ratio = (pcharge-pos[ip].GetQ())/(pcharge+pos[ip].GetQ());
	if ( (TMath::Abs(ratio)<0.2) && (pcharge!=0) ) {

	  /*
	if ( (TMath::Abs(pcharge-pos[ip].GetQ())<30) &&  // charge matching 
	     (pcharge!=0) ) { // reject combinations of bad strips
	  */


	  //
	  // add first pair
	  //	if (!(cused1[ip]||cused2[jn])){
	  if ( (fgPairs[ip*nn+jn]==100) && (neg[jn].GetQ()) ) {  //
	    
	    Float_t yn=neg[jn].GetY(); 
	    Float_t yp=pos[ip].GetY();

	    Float_t xt, zt;
	    seg->GetPadCxz(yn, yp, xt, zt);
	    
	    xbest=xt; zbest=zt; 

	    qbest =neg[jn].GetQ();

	    {
	      Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	      mT2L->MasterToLocal(loc,trk);
	      lp[0]=trk[1];
	      lp[1]=trk[2];
          }
	  
	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip].GetLabel(ilab);
	    milab[ilab+3] = neg[jn].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip].GetQ()-neg[jn].GetQ())/(pos[ip].GetQ()+neg[jn].GetQ());
	  milab[3]=(((ip<<10) + jn)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	  AliITSRecPoint * cl2;
	    cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);

	    cl2->SetChargeRatio(ratio);    	
	    cl2->SetType(7);
	    fgPairs[ip*nn+jn] =7;
	    if ((pos[ip].GetNd()+neg[jn].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      fgPairs[ip*nn+jn]=8;
	    }
	  ncl++;
	}
	//
	// add second pair
	//	if (!(cused1[ip]||cused2[jn2])){
	if ( (fgPairs[ip*nn+jn2]==100)&&(neg[jn2].GetQ() ) ) {  //

	  Float_t yn=neg[jn2].GetY(); 
	  Double_t yp=pos[ip].GetY(); 

	  Float_t xt, zt;
	  seg->GetPadCxz(yn, yp, xt, zt);
	  
	  xbest=xt; zbest=zt; 

	  qbest =neg[jn2].GetQ();

          {
          Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
          mT2L->MasterToLocal(loc,trk);
          lp[0]=trk[1];
          lp[1]=trk[2];
          }

	  lp[4]=qbest;        //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++){
	    milab[ilab] = pos[ip].GetLabel(ilab);
	    milab[ilab+3] = neg[jn2].GetLabel(ilab);
	  }
	  //
	  CheckLabels2(milab);
	  ratio = (pos[ip].GetQ()-neg[jn2].GetQ())/(pos[ip].GetQ()+neg[jn2].GetQ());
	  milab[3]=(((ip<<10) + jn2)<<10) + idet; // pos|neg|det
	  Int_t info[3] = {pos[ip].GetNd(),neg[jn2].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	  AliITSRecPoint * cl2;
	    cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);


	    cl2->SetChargeRatio(ratio);    	
	    fgPairs[ip*nn+jn2]=7;
	    cl2->SetType(7);
	    if ((pos[ip].GetNd()+neg[jn2].GetNd())>6){ //multi cluster
	      cl2->SetType(8);
	      fgPairs[ip*nn+jn2]=8;
	    }
	  ncl++;
	}
	cused1[ip]++;
	cused2[jn]++;
	cused2[jn2]++;

	} // charge matching condition

      } // 2 Nside cross 1 Pside
    } // loop over Pside clusters

  
    
    for (Int_t ip=0;ip<np;ip++){

      if(cused1[ip]) continue;


      Float_t xbest=1000,zbest=1000,qbest=0;
      //
      // 2x2 clusters
      //
      if ( (cnegative[ip]==2) && cpositive[negativepair[10*ip]]==2){ 
	Float_t minchargediff =4.;
	Float_t minchargeratio =0.2;

	Int_t j=-1;
	for (Int_t di=0;di<cnegative[ip];di++){
	  Int_t   jc = negativepair[ip*10+di];
	  Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	  ratio = (pos[ip].GetQ()-neg[jc].GetQ())/(pos[ip].GetQ()+neg[jc].GetQ()); 
	  //if (TMath::Abs(chargedif)<minchargediff){
	  if (TMath::Abs(ratio)<0.2){
	    j =jc;
	    minchargediff = TMath::Abs(chargedif);
	    minchargeratio = TMath::Abs(ratio);
	  }
	}
	if (j<0) continue;  // not proper cluster      
	

	Int_t count =0;
	for (Int_t di=0;di<cnegative[ip];di++){
	  Int_t   jc = negativepair[ip*10+di];
	  Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+3.) count++;
	}
	if (count>1) continue;  // more than one "proper" cluster for positive
	//
	
	count =0;
	for (Int_t dj=0;dj<cpositive[j];dj++){
	  Int_t   ic  = positivepair[j*10+dj];
	  Float_t chargedif = pos[ic].GetQ()-neg[j].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+3.) count++;
	}
	if (count>1) continue;  // more than one "proper" cluster for negative
	
	Int_t jp = 0;
	
	count =0;
	for (Int_t dj=0;dj<cnegative[jp];dj++){
	  Int_t   ic = positivepair[jp*10+dj];
	  Float_t chargedif = pos[ic].GetQ()-neg[jp].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+4.) count++;
	}
	if (count>1) continue;   
	if (fgPairs[ip*nn+j]<100) continue;
	//
	


	//almost gold clusters
	Float_t yp=pos[ip].GetY(); 
	Float_t yn=neg[j].GetY();      
	Float_t xt, zt;
	seg->GetPadCxz(yn, yp, xt, zt);	
	xbest=xt; zbest=zt; 
	qbest=0.5*(pos[ip].GetQ()+neg[j].GetQ());
	{
	  Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  lp[0]=trk[1];
	  lp[1]=trk[2];
	}
	lp[4]=qbest;        //Q
	for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	for (Int_t ilab=0;ilab<3;ilab++){
	  milab[ilab] = pos[ip].GetLabel(ilab);
	  milab[ilab+3] = neg[j].GetLabel(ilab);
	}
	//
	CheckLabels2(milab);
        if ((neg[j].GetQ()==0)&&(pos[ip].GetQ()==0)) continue; // reject crosses of bad strips!!
	ratio = (pos[ip].GetQ()-neg[j].GetQ())/(pos[ip].GetQ()+neg[j].GetQ());
	milab[3]=(((ip<<10) + j)<<10) + idet; // pos|neg|det
	Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	AliITSRecPoint * cl2;
	  cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);
	  	  
	  cl2->SetChargeRatio(ratio);    	
	  cl2->SetType(10);
	  fgPairs[ip*nn+j]=10;
	  if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	    cl2->SetType(11);
	    fgPairs[ip*nn+j]=11;
	  }
	  cused1[ip]++;
	  cused2[j]++;      
	ncl++;
	
      } // 2X2
    } // loop over Pside 1Dclusters



    for (Int_t ip=0;ip<np;ip++){

      if(cused1[ip]) continue;


      Float_t xbest=1000,zbest=1000,qbest=0;
      //
      // manyxmany clusters
      //
      if ( (cnegative[ip]<5) && cpositive[negativepair[10*ip]]<5){ 
	Float_t minchargediff =4.;
	Int_t j=-1;
	for (Int_t di=0;di<cnegative[ip];di++){
	  Int_t   jc = negativepair[ip*10+di];
	  Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff){
	    j =jc;
	    minchargediff = TMath::Abs(chargedif);
	  }
	}
	if (j<0) continue;  // not proper cluster      
	
	Int_t count =0;
	for (Int_t di=0;di<cnegative[ip];di++){
	  Int_t   jc = negativepair[ip*10+di];
	  Float_t chargedif = pos[ip].GetQ()-neg[jc].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+3.) count++;
	}
	if (count>1) continue;  // more than one "proper" cluster for positive
	//
	
	count =0;
	for (Int_t dj=0;dj<cpositive[j];dj++){
	  Int_t   ic  = positivepair[j*10+dj];
	  Float_t chargedif = pos[ic].GetQ()-neg[j].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+3.) count++;
	}
	if (count>1) continue;  // more than one "proper" cluster for negative
	
	Int_t jp = 0;
	
	count =0;
	for (Int_t dj=0;dj<cnegative[jp];dj++){
	  Int_t   ic = positivepair[jp*10+dj];
	  Float_t chargedif = pos[ic].GetQ()-neg[jp].GetQ();
	  if (TMath::Abs(chargedif)<minchargediff+4.) count++;
	}
	if (count>1) continue;   
	if (fgPairs[ip*nn+j]<100) continue;
	//
	
	//almost gold clusters
	Float_t yp=pos[ip].GetY(); 
	Float_t yn=neg[j].GetY();
      

	Float_t xt, zt;
	seg->GetPadCxz(yn, yp, xt, zt);
	
	xbest=xt; zbest=zt; 

	qbest=0.5*(pos[ip].GetQ()+neg[j].GetQ());

	{
	  Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  lp[0]=trk[1];
	  lp[1]=trk[2];
	}
	lp[4]=qbest;        //Q
	for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	for (Int_t ilab=0;ilab<3;ilab++){
	  milab[ilab] = pos[ip].GetLabel(ilab);
	  milab[ilab+3] = neg[j].GetLabel(ilab);
	}
	//
	CheckLabels2(milab);
        if ((neg[j].GetQ()==0)&&(pos[ip].GetQ()==0)) continue; // reject crosses of bad strips!!
	ratio = (pos[ip].GetQ()-neg[j].GetQ())/(pos[ip].GetQ()+neg[j].GetQ());
	milab[3]=(((ip<<10) + j)<<10) + idet; // pos|neg|det
	Int_t info[3] = {pos[ip].GetNd(),neg[j].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	AliITSRecPoint * cl2;
	  cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);
	  	  
	  cl2->SetChargeRatio(ratio);    	
	  cl2->SetType(12);
	  fgPairs[ip*nn+j]=12;
	  if ((pos[ip].GetNd()+neg[j].GetNd())>6){ //multi cluster
	    cl2->SetType(13);
	    fgPairs[ip*nn+j]=13;
	  }
	  cused1[ip]++;
	  cused2[j]++;      
	ncl++;
	
      } // manyXmany
    } // loop over Pside 1Dclusters
    
  } // use charge matching
  
  
  // recover all the other crosses
  //  
  for (Int_t i=0; i<np; i++) {
    Float_t xbest=1000,zbest=1000,qbest=0;
    Float_t yp=pos[i].GetY(); 
    if ((pos[i].GetQ()>0)&&(pos[i].GetQ()<3)) continue;
    for (Int_t j=0; j<nn; j++) {
    //    for (Int_t di = 0;di<cpositive[i];di++){
    //  Int_t j = negativepair[10*i+di];
      if ((neg[j].GetQ()>0)&&(neg[j].GetQ()<3)) continue;

      if ((neg[j].GetQ()==0)&&(pos[i].GetQ()==0)) continue; // reject crosses of bad strips!!

      if (cused2[j]||cused1[i]) continue;      
      if (fgPairs[i*nn+j]>0 &&fgPairs[i*nn+j]<100) continue;
      ratio = (pos[i].GetQ()-neg[j].GetQ())/(pos[i].GetQ()+neg[j].GetQ());      
      Float_t yn=neg[j].GetY();
      
      Float_t xt, zt;
      seg->GetPadCxz(yn, yp, xt, zt);
      
      if (TMath::Abs(xt)<hwSSD)
      if (TMath::Abs(zt)<hlSSD) {
	xbest=xt; zbest=zt; 

        qbest=0.5*(pos[i].GetQ()+neg[j].GetQ());

        {
        Double_t loc[3]={xbest,0.,zbest},trk[3]={0.,0.,0.};
        mT2L->MasterToLocal(loc,trk);
        lp[0]=trk[1];
        lp[1]=trk[2];
        }
        lp[4]=qbest;        //Q
	for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	for (Int_t ilab=0;ilab<3;ilab++){
	  milab[ilab] = pos[i].GetLabel(ilab);
	  milab[ilab+3] = neg[j].GetLabel(ilab);
	}
	//
	CheckLabels2(milab);
	milab[3]=(((i<<10) + j)<<10) + idet; // pos|neg|det
	Int_t info[3] = {pos[i].GetNd(),neg[j].GetNd(),fNlayer[fModule]};

	lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	// out-of-diagonal element of covariance matrix
 	if( (info[0]==1) && (info[1]==1) ) lp[5]=-0.00012;
	else if ( (info[0]>1) && (info[1]>1) ) { 
	  lp[2]=2.63e-06;    // 0.0016*0.0016;  //SigmaY2
	  lp[3]=0.0065;      // 0.08*0.08;   //SigmaZ2
	  lp[5]=-6.48e-05;
	}
	else {
	  lp[2]=4.80e-06;      // 0.00219*0.00219
	  lp[3]=0.0093;        // 0.0964*0.0964;
	  if (info[0]==1) {
	    lp[5]=-0.00014;
	  }
	  else { 
	    lp[2]=2.79e-06;    // 0.0017*0.0017; 
	    lp[3]=0.00935;     // 0.967*0.967;
	    lp[5]=-4.32e-05;
	  }
	}

	AliITSRecPoint * cl2;
	  cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);

	  cl2->SetChargeRatio(ratio);
	  cl2->SetType(100+cpositive[j]+cnegative[i]);	  

	  if(pos[i].GetQ()==0) cl2->SetType(200+cpositive[j]+cnegative[i]);
	  if(neg[j].GetQ()==0) cl2->SetType(300+cpositive[j]+cnegative[i]);
      	ncl++;
      }
    }
  }



  if(repa->GetUseBadChannelsInClusterFinderSSD()==kTRUE) {
    
    //---------------------------------------------------------
    // recover crosses of good 1D clusters with bad strips on the other side
    // Note1: at first iteration skip modules with a bad side (or almost), (would produce too many fake!) 
    // Note2: for modules with a bad side see below 
    
    AliITSCalibrationSSD* cal = (AliITSCalibrationSSD*)GetResp(fModule);
    Int_t countPbad=0, countNbad=0;
    for(Int_t ib=0; ib<768; ib++) {
      if(cal->IsPChannelBad(ib)) countPbad++;
      if(cal->IsNChannelBad(ib)) countNbad++;
    }
    //  AliInfo(Form("module %d has %d P- and %d N-bad strips",fModule,countPbad,countNbad));
    
    if( (countPbad<100) && (countNbad<100) ) { // no bad side!!
      
      for (Int_t i=0; i<np; i++) { // loop over Nside 1Dclusters with no crosses
	if(cnegative[i]) continue; // if intersecting Pside clusters continue;
	
	//      for(Int_t ib=0; ib<768; ib++) { // loop over all Pstrips
	for(Int_t ib=15; ib<753; ib++) { // loop over all Pstrips
	  
	  if(cal->IsPChannelBad(ib)) { // check if strips is bad
	    Float_t yN=pos[i].GetY();	
	    Float_t xt, zt;
	    seg->GetPadCxz(1.*ib, yN, xt, zt);	
	    
	    //----------
	    // bad Pstrip is crossing the Nside 1Dcluster -> create recpoint
	    // 
	    if ( (TMath::Abs(xt)<hwSSD) && (TMath::Abs(zt)<hlSSD) ) {
	      Double_t loc[3]={xt,0.,zt},trk[3]={0.,0.,0.};
	      mT2L->MasterToLocal(loc,trk);
	      lp[0]=trk[1];
	      lp[1]=trk[2];        
	      lp[4]=pos[i].GetQ(); //Q
	      for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	      for (Int_t ilab=0;ilab<3;ilab++) milab[ilab] = pos[i].GetLabel(ilab);	  
	      CheckLabels2(milab);
	      milab[3]=( (i<<10) << 10 ) + idet; // pos|neg|det
	      Int_t info[3] = {pos[i].GetNd(),0,fNlayer[fModule]};
	      
	      lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	      lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	      lp[5]=-0.00012; // out-of-diagonal element of covariance matrix
	      if (info[0]>1) {
		lp[2]=4.80e-06;
		lp[3]=0.0093;
		lp[5]=0.00014;
	      }
	      	      
	      AliITSRecPoint * cl2;
		cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);	    
		cl2->SetChargeRatio(1.);
	      cl2->SetType(50);	  
	      ncl++;
	    } // cross is within the detector
	    //
	    //--------------
	    
	  } // bad Pstrip
	  
	} // end loop over Pstrips
	
      } // end loop over Nside 1D clusters
      
      for (Int_t j=0; j<nn; j++) { // loop over Pside 1D clusters with no crosses
	if(cpositive[j]) continue;
	
	//      for(Int_t ib=0; ib<768; ib++) { // loop over all Nside strips
	for(Int_t ib=15; ib<753; ib++) { // loop over all Nside strips
	  
	  if(cal->IsNChannelBad(ib)) { // check if strip is bad
	    Float_t yP=neg[j].GetY();	
	    Float_t xt, zt;
	    seg->GetPadCxz(yP, 1.*ib, xt, zt);	
	    
	    //----------
	    // bad Nstrip is crossing the Pside 1Dcluster -> create recpoint
	    // 
	    if ( (TMath::Abs(xt)<hwSSD) && (TMath::Abs(zt)<hlSSD) ) {
	      Double_t loc[3]={xt,0.,zt},trk[3]={0.,0.,0.};
	      mT2L->MasterToLocal(loc,trk);
	      lp[0]=trk[1];
	      lp[1]=trk[2];        
	      lp[4]=neg[j].GetQ(); //Q
	      for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	      for (Int_t ilab=0;ilab<3;ilab++) milab[ilab] = neg[j].GetLabel(ilab);	  
	      CheckLabels2(milab);
	      milab[3]=( j << 10 ) + idet; // pos|neg|det
	      Int_t info[3]={0,(Int_t)neg[j].GetNd(),fNlayer[fModule]};
	      
	      lp[2]=4.968e-06;     // 0.00223*0.00223;  //SigmaY2
	      lp[3]=0.012;         // 0.110*0.110;  //SigmaZ2
	      lp[5]=-0.00012; // out-of-diagonal element of covariance matrix
	      if (info[0]>1) {
		lp[2]=2.79e-06;
		lp[3]=0.00935;
		lp[5]=-4.32e-05;
	      }
	      
	      AliITSRecPoint * cl2;
		cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);	    
		cl2->SetChargeRatio(1.);
		cl2->SetType(60);	  
	      ncl++;
	    } // cross is within the detector
	    //
	    //--------------
	    
	  } // bad Nstrip
	} // end loop over Nstrips
      } // end loop over Pside 1D clusters
      
    } // no bad sides 
    
    //---------------------------------------------------------
    
    else if( (countPbad>700) && (countNbad<100) ) { // bad Pside!!
      
      for (Int_t i=0; i<np; i++) { // loop over Nside 1Dclusters with no crosses
	if(cnegative[i]) continue; // if intersecting Pside clusters continue;
	
	Float_t xt, zt;
	Float_t yN=pos[i].GetY();	
	Float_t yP=0.;
	if (seg->GetLayer()==5) yP = yN + (7.6/1.9);
	else yP = yN - (7.6/1.9);
	seg->GetPadCxz(yP, yN, xt, zt);	
	
	if ( (TMath::Abs(xt)<hwSSD) && (TMath::Abs(zt)<hlSSD) ) {
	  Double_t loc[3]={xt,0.,zt},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  lp[0]=trk[1];
	  lp[1]=trk[2];        
	  lp[4]=pos[i].GetQ(); //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++) milab[ilab] = pos[i].GetLabel(ilab);	  
	  CheckLabels2(milab);
	  milab[3]=( (i<<10) << 10 ) + idet; // pos|neg|det
	  Int_t info[3] = {(Int_t)pos[i].GetNd(),0,fNlayer[fModule]};
	  
	  lp[2]=0.00098;    // 0.031*0.031;  //SigmaY2
	  lp[3]=1.329;      // 1.15*1.15;  //SigmaZ2
	  lp[5]=-0.0359;
	  if(info[0]>1) lp[2]=0.00097;

	  AliITSRecPoint * cl2;
	    cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);	    
	    cl2->SetChargeRatio(1.);
	    cl2->SetType(70);	  
	  ncl++;
	} // cross is within the detector
	//
	//--------------
	
      } // end loop over Nside 1D clusters
      
    } // bad Pside module
    
    else if( (countNbad>700) && (countPbad<100) ) { // bad Nside!!
      
      for (Int_t j=0; j<nn; j++) { // loop over Pside 1D clusters with no crosses
	if(cpositive[j]) continue;
	
	Float_t xt, zt;
	Float_t yP=neg[j].GetY();	
	Float_t yN=0.;
	if (seg->GetLayer()==5) yN = yP - (7.6/1.9);
	else yN = yP + (7.6/1.9);
	seg->GetPadCxz(yP, yN, xt, zt);	
	
	if ( (TMath::Abs(xt)<hwSSD) && (TMath::Abs(zt)<hlSSD) ) {
	  Double_t loc[3]={xt,0.,zt},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  lp[0]=trk[1];
	  lp[1]=trk[2];        
	  lp[4]=neg[j].GetQ(); //Q
	  for (Int_t ilab=0;ilab<10;ilab++) milab[ilab]=-2;
	  for (Int_t ilab=0;ilab<3;ilab++) milab[ilab] = neg[j].GetLabel(ilab);	  
	  CheckLabels2(milab);
	  milab[3]=( j << 10 ) + idet; // pos|neg|det
	  Int_t info[3] = {0,(Int_t)neg[j].GetNd(),fNlayer[fModule]};
	  
	  lp[2]=7.27e-05;   // 0.0085*0.0085;  //SigmaY2
	  lp[3]=1.33;       // 1.15*1.15;  //SigmaZ2
	  lp[5]=0.00931;
	  if(info[1]>1) lp[2]=6.91e-05;
	  
	  AliITSRecPoint * cl2;
	    cl2 = new ((*clusters)[ncl]) AliITSRecPoint(milab,lp,info);	    
	    cl2->SetChargeRatio(1.);
	    cl2->SetType(80);	  
	  ncl++;
	} // cross is within the detector
	//
	//--------------
	
      } // end loop over Pside 1D clusters
      
    } // bad Nside module
    
    //---------------------------------------------------------
    
  } // use bad channels
    
  //cout<<ncl<<" clusters for this module"<<endl;

  delete [] cnegative;
  delete [] cused1;
  delete [] negativepair;
  delete [] cpositive;
  delete [] cused2;
  delete [] positivepair;

}
