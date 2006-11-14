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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  fast TPC cluster simulation                                              //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include <TVector.h>
#include <TRandom.h>

#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliTPC.h"
#include "AliTPCParam.h"
#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"
#include "AliTPCcluster.h"
#include "AliComplexCluster.h"
#include "AliTPCFast.h"
#include "AliLog.h"

ClassImp(AliTPCFast)
  //____________________________________________________________________
AliTPCFast::AliTPCFast(const AliTPCFast &param)
              :TObject(param),fParam(0)

{
  //
  //  copy constructor - dummy
  //
  fParam = param.fParam;
}
AliTPCFast & AliTPCFast::operator =(const AliTPCFast & param)
{
  //
  // assignment operator - dummy
  //
  fParam=param.fParam;
  return (*this);
}

//_____________________________________________________________________________
void AliTPCFast::Hits2Clusters(AliRunLoader* runLoader) const
{
  //--------------------------------------------------------
  // TPC simple cluster generator from hits
  // obtained from the TPC Fast Simulator
  // The point errors are taken from the parametrization
  //--------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  // Adopted to Marian's cluster data structure by I.Belikov, CERN,
  // Jouri.Belikov@cern.ch
  //----------------------------------------------------------------
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //---------------------------------------------------------------------
  //   ALICE TPC Cluster Parameters
  //--------------------------------------------------------------------
       
  

  // Cluster width in rphi
  const Float_t kACrphi=0.18322;
  const Float_t kBCrphi=0.59551e-3;
  const Float_t kCCrphi=0.60952e-1;
  // Cluster width in z
  const Float_t kACz=0.19081;
  const Float_t kBCz=0.55938e-3;
  const Float_t kCCz=0.30428;


  AliLoader* loader = runLoader->GetLoader("TPCLoader");
  if (!loader) {
    AliError("No TPC loader found");
    return;
  }
  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  AliRun* aliRun = runLoader->GetAliRun();
  if (!aliRun) {
    AliError("Couldn't get AliRun object");
    return;
  }
  AliTPC* tpc = (AliTPC*) aliRun->GetDetector("TPC");
  if (!tpc) {
    AliError("Couldn't get TPC detector");
    return;
  }
  AliTPCParam* param = tpc->GetParam();
  if (!param) {
    AliError("No TPC parameters available");
    return;
  }

  //if(fDefaults == 0) SetDefaults();

  Float_t sigmaRphi,sigmaZ,clRphi,clZ;
  //
  TParticle *particle; // pointer to a given particle
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  Int_t sector;
  Int_t ipart;
  Float_t xyz[5];
  Float_t pl,pt,tanth,rpad,ratio;
  Float_t cph,sph;
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *tH = loader->TreeH();
  if (tH == 0x0)
    AliFatal("Can not find TreeH in folder");

  tpc->SetTreeAddress();
  
  Stat_t ntracks = tH->GetEntries();

  //Switch to the output file
  
  if (loader->TreeR() == 0x0) loader->MakeTree("R");
  
  AliDebug(1,Form("param->GetTitle() = %s",param->GetTitle()));
  
  runLoader->CdGAFile();
  //param->Write(param->GetTitle());

  AliTPCClustersArray carray;
  carray.Setup(param);
  carray.SetClusterType("AliTPCcluster");
  carray.MakeTree(loader->TreeR());

  Int_t nclusters=0; //cluster counter
  
  //------------------------------------------------------------
  // Loop over all sectors (72 sectors for 20 deg
  // segmentation for both lower and upper sectors)
  // Sectors 0-35 are lower sectors, 0-17 z>0, 17-35 z<0
  // Sectors 36-71 are upper sectors, 36-53 z>0, 54-71 z<0
  //
  // First cluster for sector 0 starts at "0"
  //------------------------------------------------------------
   
  for(Int_t isec=0;isec<param->GetNSector();isec++){
    //MI change
    param->AdjustCosSin(isec,cph,sph);
    
    //------------------------------------------------------------
    // Loop over tracks
    //------------------------------------------------------------
    
    for(Int_t track=0;track<ntracks;track++){
      tpc->ResetHits();
      tH->GetEvent(track);
      //
      //  Get number of the TPC hits
      //     
       tpcHit = (AliTPChit*)tpc->FirstHit(-1);

      // Loop over hits
      //
       while(tpcHit){
	 
	 if (tpcHit->fQ == 0.) {
           tpcHit = (AliTPChit*) tpc->NextHit();
           continue; //information about track (I.Belikov)
	 }
	 sector=tpcHit->fSector; // sector number
	 
	 if(sector != isec){
	   tpcHit = (AliTPChit*) tpc->NextHit();
	   continue; 
	 }
	 ipart=tpcHit->Track();
	 particle=aliRun->GetMCApp()->Particle(ipart);
	 pl=particle->Pz();
	 pt=particle->Pt();
	 if(pt < 1.e-9) pt=1.e-9;
	 tanth=pl/pt;
	 tanth = TMath::Abs(tanth);
	 rpad=TMath::Sqrt(tpcHit->X()*tpcHit->X() + tpcHit->Y()*tpcHit->Y());
	 ratio=0.001*rpad/pt; // pt must be in MeV/c - historical reason
	 
	 //   space-point resolutions
	 
	sigmaRphi=AliTPCcluster::SigmaY2(rpad,tanth,pt);
	sigmaZ   =AliTPCcluster::SigmaZ2(rpad,tanth   );
	
	//   cluster widths
	
	clRphi=kACrphi-kBCrphi*rpad*tanth+kCCrphi*ratio*ratio;
	clZ=kACz-kBCz*rpad*tanth+kCCz*tanth*tanth;
	
	// temporary protection
	
	if(sigmaRphi < 0.) sigmaRphi=0.4e-3;
	if(sigmaZ < 0.) sigmaZ=0.4e-3;
	if(clRphi < 0.) clRphi=2.5e-3;
	if(clZ < 0.) clZ=2.5e-5;
	
	//
	
	//
	// smearing --> rotate to the 1 (13) or to the 25 (49) sector,
	// then the inaccuracy in a X-Y plane is only along Y (pad row)!
	//
        Float_t xprim= tpcHit->X()*cph + tpcHit->Y()*sph;
	Float_t yprim=-tpcHit->X()*sph + tpcHit->Y()*cph;
	xyz[0]=gRandom->Gaus(yprim,TMath::Sqrt(sigmaRphi));   // y
          Float_t alpha=(isec < param->GetNInnerSector()) ?
	  param->GetInnerAngle() : param->GetOuterAngle();
          Float_t ymax=xprim*TMath::Tan(0.5*alpha);
          if (TMath::Abs(xyz[0])>ymax) xyz[0]=yprim; 
	xyz[1]=gRandom->Gaus(tpcHit->Z(),TMath::Sqrt(sigmaZ)); // z
          if (TMath::Abs(xyz[1])>param->GetZLength()) xyz[1]=tpcHit->Z(); 
	xyz[2]=sigmaRphi;                                     // fSigmaY2
	xyz[3]=sigmaZ;                                        // fSigmaZ2
	xyz[4]=tpcHit->fQ;                                    // q

        AliTPCClustersRow *clrow=carray.GetRow(sector,tpcHit->fPadRow);
        if (!clrow) clrow=carray.CreateRow(sector,tpcHit->fPadRow);	

        Int_t tracks[3]={tpcHit->Track(), -1, -1};
	AliTPCcluster cluster(tracks,xyz);

        clrow->InsertCluster(&cluster); nclusters++;

        tpcHit = (AliTPChit*)tpc->NextHit();
        

      } // end of loop over hits

    }   // end of loop over tracks

    Int_t nrows=param->GetNRow(isec);
    for (Int_t irow=0; irow<nrows; irow++) {
        AliTPCClustersRow *clrow=carray.GetRow(isec,irow);
        if (!clrow) continue;
        carray.StoreRow(isec,irow);
        carray.ClearRow(isec,irow);
    }

  } // end of loop over sectors  

  //  cerr<<"Number of made clusters : "<<nclusters<<"                        \n";
  loader->WriteRecPoints("OVERWRITE");
  
  
} // end of function



//_________________________________________________________________



void  AliTPCFast::Hits2ExactClusters(AliRunLoader* runLoader) const{



  AliLoader* loader = runLoader->GetLoader("TPCLoader");
  if (!loader) {
    AliError("No TPC loader found");
    return;
  }
  AliTPC* tpc = (AliTPC*) gAlice->GetDetector("TPC");
  if (!tpc) {
    AliError("Couldn't get TPC detector");
    return;
  }
  AliTPCParam* param = tpc->GetParam();
  if (!param) {
    AliError("No TPC parameters available");
    return;
  }
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *tH = loader->TreeH();
  if (tH == 0x0) { runLoader->LoadHits("TPC"); loader->TreeH();}

  tpc->SetTreeAddress();
  

  //Switch to the output file
  
  if (loader->TreeR() == 0x0) loader->MakeTree("R");
  
  AliDebug(1,Form("param->GetTitle() = %s",param->GetTitle()));
  
  runLoader->CdGAFile();
  //param->Write(param->GetTitle());

  AliTPCClustersArray carray;
  carray.Setup(param);
  carray.SetClusterType("AliTPCclusterMI");
  carray.MakeTree(loader->TreeR());
  
  //------------------------------------------------------------
  // Loop over all sectors (72 sectors for 20 deg)
  //------------------------------------------------------------
   
  for(Int_t isec=0;isec<param->GetNSector();isec++){
    Hits2ExactClustersSector(runLoader, &carray, isec);    
  } // end of loop over sectors  
  loader->WriteRecPoints("OVERWRITE");
}




//_________________________________________________________________
void AliTPCFast::Hits2ExactClustersSector(AliRunLoader* runLoader,
					  AliTPCClustersArray* clustersArray,
					  Int_t isec) const
{
  //--------------------------------------------------------
  //calculate exact cross point of track and given pad row
  //resulting values are expressed in "local" coordinata
  //the sigmaY2 and sigmaZ2 of cluster are the shape of cluster parameters
  //                        - thay are used later on for error parameterization
  //--------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marian Ivanov  GSI Darmstadt, m.ivanov@gsi.de
  //-----------------------------------------------------------------
  //
  if (clustersArray==0){    
    return;
  }
  AliLoader* loader = runLoader->GetLoader("TPCLoader");
  if (!loader) {
    AliError("No TPC loader found");
    return;
  }
  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  AliRun* aliRun = runLoader->GetAliRun();
  if (!aliRun) {
    AliError("Couldn't get AliRun object");
    return;
  }
  AliTPC* tpc = (AliTPC*) aliRun->GetDetector("TPC");
  if (!tpc) {
    AliError("Couldn't get TPC detector");
    return;
  }
  AliTPCParam* param = tpc->GetParam();
  if (!param) {
    AliError("No TPC parameters available");
    return;
  }
  //
  //
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  //  Int_t sector,nhits;
  Int_t ipart;
  const Int_t kcmaxhits=30000;
  TVector * xxxx = new TVector(kcmaxhits*4);
  TVector & xxx = *xxxx;
  Int_t maxhits = kcmaxhits;
  //construct array for each padrow
  for (Int_t i=0; i<param->GetNRow(isec);i++) 
    clustersArray->CreateRow(isec,i);
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *tH = loader->TreeH();
  if (tH == 0x0)
    AliFatal("Can not find TreeH in folder");

  tpc->SetTreeAddress();

  Stat_t ntracks = tH->GetEntries();
  //MI change
  TBranch * branch=0;
  if (tpc->GetHitType()>1) branch = tH->GetBranch("TPC2");
  else branch = tH->GetBranch("TPC");

  //------------------------------------------------------------
  // Loop over tracks
  //------------------------------------------------------------

  for(Int_t track=0;track<ntracks;track++){ 
    Bool_t isInSector=kTRUE;
    tpc->ResetHits();
    isInSector = tpc->TrackInVolume(isec,track);
    if (!isInSector) continue;
    //MI change
    branch->GetEntry(track); // get next track
    //
    // Loop over hits
    //
    Int_t currentIndex=0;
    Int_t lastrow=-1;  //last writen row

    //M.I. changes

    tpcHit = (AliTPChit*)tpc->FirstHit(-1);
    while(tpcHit){
      
      Int_t sector=tpcHit->fSector; // sector number
      if(sector != isec){
	tpcHit = (AliTPChit*) tpc->NextHit();
	continue; 
      }

      ipart=tpcHit->Track();
      
      //find row number

      Float_t  x[3]={tpcHit->X(),tpcHit->Y(),tpcHit->Z()};
      Int_t    index[3]={1,isec,0};
      Int_t    currentrow = param->GetPadRow(x,index);
      
      if (currentrow<0) {tpcHit = (AliTPChit*)tpc->NextHit(); continue;}
      if (lastrow<0) lastrow=currentrow;
      //
      //
      if (currentrow!=lastrow){
	if (currentIndex>2){
	  Float_t sumx=0;
	  Float_t sumx2=0;
	  Float_t sumx3=0;
	  Float_t sumx4=0;
	  Float_t sumy=0;
	  Float_t sumxy=0;
	  Float_t sumx2y=0;
	  Float_t sumz=0;
	  Float_t sumxz=0;
	  Float_t sumx2z=0;
	  Float_t sumq=0;
	  for (Int_t index=0;index<currentIndex;index++){
	    Float_t x,x2,x3,x4;
	    x=x2=x3=x4=xxx(index*4);
	    x2*=x;
	    x3*=x2;
	    x4*=x3;
	    sumx+=x;
	    sumx2+=x2;
	    sumx3+=x3;
	    sumx4+=x4;
	    sumy+=xxx(index*4+1);
	    sumxy+=xxx(index*4+1)*x;
	    sumx2y+=xxx(index*4+1)*x2;
	    sumz+=xxx(index*4+2);
	    sumxz+=xxx(index*4+2)*x;
	    sumx2z+=xxx(index*4+2)*x2;	 
	    sumq+=xxx(index*4+3);
	  }
	  Float_t det=currentIndex*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx3-sumx2*sumx2);
	  
	  Float_t detay=sumy*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxy*sumx4-sumx2y*sumx3)+
	    sumx2*(sumxy*sumx3-sumx2y*sumx2);
	  Float_t detaz=sumz*(sumx2*sumx4-sumx3*sumx3)-sumx*(sumxz*sumx4-sumx2z*sumx3)+
	    sumx2*(sumxz*sumx3-sumx2z*sumx2);
	  
	  Float_t detby=currentIndex*(sumxy*sumx4-sumx2y*sumx3)-sumy*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx2y-sumx2*sumxy);
	  Float_t detbz=currentIndex*(sumxz*sumx4-sumx2z*sumx3)-sumz*(sumx*sumx4-sumx2*sumx3)+
	    sumx2*(sumx*sumx2z-sumx2*sumxz);
	  
	  if (TMath::Abs(det)<0.00000000000000001){
	     tpcHit = (AliTPChit*)tpc->NextHit();
	    continue;
	  }
	
	  Float_t y=detay/det;
	  Float_t z=detaz/det;		  
	  Float_t by=detby/det; //y angle
	  Float_t bz=detbz/det; //z angle
	  sumy/=Float_t(currentIndex);
	  sumz/=Float_t(currentIndex);


	  //
	  //
	  Float_t sign = (tpcHit->Z()>0)? 1.:-1.;
	  y = (y+0.5)*param->GetPadPitchWidth(isec);
	  z = z*param->GetZWidth();
	  //
	  // y expected shape
	  Double_t sigmay2 = z*fParam->GetDiffL();
	  sigmay2 *= sigmay2;
	  Float_t angulary = by*param->GetPadPitchLength(isec,lastrow);
	  angulary*=angulary;
	  angulary/=12;
	  sigmay2 +=angulary+0.25*param->GetPadPitchWidth(isec)*param->GetPadPitchWidth(isec);
	  //
	  // z expected shape
	  Double_t sigmaz2 = z*fParam->GetDiffT();
	  sigmaz2 *= sigmaz2;
	  Float_t angularz = bz*param->GetPadPitchLength(isec,lastrow);
	  angularz*=angularz;
	  angularz/=12;
	  sigmaz2 +=angularz+0.25;
	  //
	  sigmaz2 = TMath::Min(sigmaz2,1.);
	  sigmay2 = TMath::Min(sigmay2,1.);
	  //
	  //
	  z = sign*(param->GetZLength() - z);
	  if (TMath::Abs(z)< param->GetZLength()-1){
	    AliTPCClustersRow * row = (clustersArray->GetRow(isec,lastrow));
	    if (row!=0) {
	      AliTPCclusterMI* cl = new((AliTPCclusterMI*)row->Append()) AliTPCclusterMI ;
	      //	  AliTPCclusterMI cl;	    
	      cl->SetZ(z);
	      cl->SetY(y);
	      cl->SetQ(sumq);
	      if (TMath::Abs(sumq)<1000){
		cl->SetMax(Short_t(sumq));
	      }
	      else{
		cl->SetMax(0);
	      }
	      cl->SetSigmaZ2(sigmaz2);
	      cl->SetSigmaY2(sigmay2);
	      cl->SetLabel(ipart,0);
	      cl->SetLabel(-1,1);
	      cl->SetLabel(-1,2);
	      cl->SetType(0);
	    }
	  }
	} //end of calculating cluster for given row		
	currentIndex=0;
	lastrow=currentrow;
      }  //  end of crossing rows
      //
      if ( currentIndex>=maxhits){
	maxhits+=kcmaxhits;
	xxx.ResizeTo(4*maxhits);
      }     
      xxx(currentIndex*4)=x[0];
      xxx(currentIndex*4+1)=x[1];
      xxx(currentIndex*4+2)=x[2];	
      xxx(currentIndex*4+3)=tpcHit->fQ;
      currentIndex++;	
      tpcHit = (AliTPChit*)tpc->NextHit();      
    } // end of loop over hits
  }   // end of loop over tracks 
  //write padrows to tree 
  for (Int_t ii=0; ii<param->GetNRow(isec);ii++) {
    clustersArray->StoreRow(isec,ii);    
    clustersArray->ClearRow(isec,ii);        
  }
  xxxx->Delete();
 
}


