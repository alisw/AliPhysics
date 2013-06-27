#include <iostream>
#include <AliMagF.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH3F.h>
#include <TClonesArray.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TPad.h>
#include <AliTrackPointArray.h>

#include "AliToyMCDrawer.h"
#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"
#include <AliCDBManager.h>
#include <AliTPCParam.h>
#include <AliGeomManager.h>
#include <TGeoGlobalMagField.h>
#include <AliTPCcalibDB.h>
#include <AliTPCROC.h>


ClassImp(AliToyMCDrawer);

AliToyMCDrawer::AliToyMCDrawer()
 : TObject()   
 ,fInputTree(0x0)
 ,inFile(0x0)
 ,fFileName(0x0)
 ,fEvent(0x0)
 ,fEventArray(0x0)
 ,fDispHist(0x0)
 ,fCenterTime(-1.)
 ,fDriftVel(-1.)
   
 {
   fEventArray = new TClonesArray("AliToyMCEvent");
   
   fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
   fTPCParam->ReadGeoMatrices();
   fDriftVel = fTPCParam->GetDriftV();///1000000;
   fMaxZ0=fTPCParam->GetZLength();
   fTimeRange = 2*fMaxZ0/fDriftVel;
   fIFCRadius=  83.5;
   fOFCRadius= 254.5;   
 
   fRoc  = AliTPCROC::Instance();
   fPoints = new TClonesArray("TPolyMarker3D");
   fDistPoints = new TClonesArray("TPolyMarker3D");
   
 }
//________________________________________________________________
AliToyMCDrawer::AliToyMCDrawer(const AliToyMCDrawer &drawer)
  : TObject(drawer)
  ,fInputTree(drawer.fInputTree)
  ,inFile(drawer.inFile)
  ,fFileName(drawer.fFileName)
  ,fEvent(drawer.fEvent)
  ,fEventArray(drawer.fEventArray)
  ,fDispHist(drawer.fDispHist)
  , fCenterTime(drawer.fCenterTime)
  ,fDriftVel(drawer.fDriftVel)
  ,fTPCParam(drawer.fTPCParam)
  ,fMaxZ0(drawer.fMaxZ0)
  ,fOFCRadius(drawer.fOFCRadius)
  ,fIFCRadius(drawer.fIFCRadius)
  ,fTimeRange(drawer.fTimeRange)
  ,fRoc(drawer.fRoc)
  ,fPoints(drawer.fPoints)
  ,fDistPoints(drawer.fDistPoints)
{
  //
}
//_____________________________________________________
AliToyMCDrawer& AliToyMCDrawer::operator = (const AliToyMCDrawer &drawer)
{
  //assignment operator
  if (&drawer == this) return *this;
  new (this) AliToyMCDrawer(drawer);

  return *this;
}

//________________________________________________________________
AliToyMCDrawer::~AliToyMCDrawer()
{
  //destructor
  delete fEvent;
  delete fEventArray;
  delete fDispHist;
  delete fPoints;
  delete fDistPoints; 
}
//________________________________________________________________
Int_t AliToyMCDrawer::FillEventArray(Double_t snapShotTime)
{
  if(!fFileName) {
    std::cout << "no input file provided, using default (toyMC.root)" << std::endl;
    fFileName = "toyMC.root";
  }
  
  inFile = new TFile(fFileName,"read");
  fInputTree = dynamic_cast<TTree*> (inFile->Get("toyMCtree"));
  fInputTree->SetBranchAddress("event",&fEvent);
  fEventArray->Clear();
  fCenterTime = snapShotTime;
  Int_t leftSearchIndex  = Int_t(snapShotTime/2e-5);
  if(leftSearchIndex>fInputTree->GetEntries()) leftSearchIndex = fInputTree->GetEntries()-1;
  if(leftSearchIndex<0) leftSearchIndex = 0;
  
  fInputTree->GetEvent(leftSearchIndex);
  Double_t leftTime = fEvent ->GetT0();
  Int_t firstEventIndex;
  //if fEvent exactly at snapShotTime, check if there are more events at the same time but with lower indices
  if( TMath::Abs(leftTime - snapShotTime)<5e-9) 
    {
      //  printf("close from start, lefttime = %2.2f\n",leftTime);
      firstEventIndex = leftSearchIndex+1;
      while( TMath::Abs(leftTime- snapShotTime) < 5e-9){
	firstEventIndex--;
	//	printf("close from start, lefttime = %2.2f\n",leftTime);
	if(fInputTree->GetEvent(firstEventIndex-1))
	  {
	    leftTime=fEvent->GetT0();
	    //     printf("close from start, lefttime = %2.2f\n",leftTime);
	  }
	else break;
      }


    }
  else
    {
      Int_t rightSearchIndex = leftSearchIndex;
      Double_t rightTime;
      
      Int_t direction = 1; //go right
      if(leftTime > snapShotTime) {
	direction = -1;
      }
      rightSearchIndex+= direction;
      fInputTree->GetEvent(rightSearchIndex);
      rightTime = fEvent->GetT0();
      //std::cout << "bef search " << leftTime << " " << rightTime<< " " << " " << snapShotTime << " "<< (leftTime-snapShotTime)*(rightTime-snapShotTime)  << std::endl;
      while ( (leftTime-snapShotTime)*(rightTime-snapShotTime)>0   ){
	//		printf("searching, leftime %f, righttim %f\n",leftTime,rightTime);
	rightSearchIndex += direction; 
	leftSearchIndex +=direction;
	fInputTree->GetEvent(leftSearchIndex);
	leftTime = fEvent->GetT0();
	fInputTree->GetEvent(rightSearchIndex);
	rightTime = fEvent->GetT0();
	//	printf("searching, leftime %f, righttim %f\n",leftTime,rightTime);
      }
      if (direction==-1) rightSearchIndex = leftSearchIndex;
      firstEventIndex = rightSearchIndex;

    }
  fInputTree->GetEvent(firstEventIndex);
  //std::cout <<"first event after sn time: " << fEvent->GetT0() << std::endl;
  return FillEventArray(firstEventIndex);
  
}
//________________________________________________________________
Int_t AliToyMCDrawer::FillEventArray(Int_t middleEventNbr)
{
  if(!fFileName) {
    std::cout << "no input file provided, using default (toyMC.root)" << std::endl;
    fFileName = "toyMC.root";
  }
  if(!inFile) inFile = new TFile(fFileName,"read");
  if(!fInputTree) 
    {
      fInputTree = dynamic_cast<TTree*> (inFile->Get("toyMCtree"));
      //   fInputTree->Add(fileName);
      fInputTree->SetBranchAddress("event",&fEvent);
      if(fInputTree->GetEvent(middleEventNbr)) fCenterTime = fEvent->GetT0();
    }
  Double_t centerEventTime;
  if(fInputTree->GetEvent(middleEventNbr)) 
    {
      centerEventTime = fEvent->GetT0();
      
    }
  else return 0;
  fEventArray->Clear();
  if(fInputTree->GetEvent(middleEventNbr)){
    new((*fEventArray)[fEventArray->GetEntriesFast()]) AliToyMCEvent(*fEvent);
    
    //add events after middle event
    Int_t eventIndex = middleEventNbr + 1;
    Double_t currentTime = centerEventTime;
    if(fInputTree->GetEvent(eventIndex)) {
      currentTime = fEvent->GetT0();
      while(currentTime - centerEventTime < fTimeRange/2) {
	new((*fEventArray)[fEventArray->GetEntriesFast()]) AliToyMCEvent(*fEvent);
	eventIndex++;
	if(!fInputTree->GetEvent(eventIndex)) break;
	currentTime = fEvent->GetT0();
      }
    }
    //add events before middle event
    eventIndex = middleEventNbr - 1;
    if(fInputTree->GetEvent(eventIndex)){
      currentTime = fEvent->GetT0();
      while(centerEventTime-currentTime < fTimeRange/2) {
      
	new((*fEventArray)[fEventArray->GetEntriesFast()]) AliToyMCEvent(*fEvent);
	eventIndex--;
	if(!fInputTree->GetEvent(eventIndex)) break;
	currentTime = fEvent->GetT0();
      }
    }
    return 1;
  }
  else {
    printf("Selected event number (%d) out of range!\n",middleEventNbr);
    return 0;
  }
  //return fEventArray;
  inFile->Close();
}
//________________________________________________________________
void AliToyMCDrawer::DrawEvents(Bool_t both, Bool_t before)
{
  
  fPoints->Clear();
  fDistPoints->Clear();
  DrawGeometry();
 
  
   Double_t phiCut = 1.;
   TLegend *leg = new TLegend(0.1,0.75,0.6,1,Form("Snapshot of TPC at %f micros after sim start, freq. = , bunchcr.rate = .",1000000*fCenterTime));
   leg->AddEntry((TObject*)0,Form("%2.2f<#phi<#pi-%2.2f && #pi+%2.2f<#phi<2#pi-%2.2f",phiCut,phiCut,phiCut,phiCut),"");

   
 

  for(Int_t iEvent = 0; iEvent<fEventArray->GetEntriesFast(); iEvent++){
 
 
     AliToyMCEvent *currentEvent = static_cast<AliToyMCEvent*> (fEventArray->At(iEvent));
     Double_t currentEventTime = currentEvent->GetT0();
     
     if((1000000*currentEventTime)>((1000000*fCenterTime)+0.05   ))  {
       printf("larger iEvent: %d, current %f, center %f\n",iEvent, currentEventTime,fCenterTime);
       printf("after\n");
       if ( !(both || !before)) continue;

     }
     if( currentEventTime<=fCenterTime)// && !(currentEventTime==fCenterTime &&both   ))
       { 
	 printf("smaller iEvent: %d, current %f, center %f\n",iEvent, currentEventTime,fCenterTime);
	 printf("before\n");
	 if ( !(both || before)) continue;
       }
     
     //std:: cout << "iev " << iEvent << " " << color << std::endl;
     Double_t currentEventTimeRelToCentral = fCenterTime - currentEventTime;
     
     Int_t color = ((currentEvent->GetEventNumber())%10);
     if(color == 0||color ==1) {
       color==0?color=2:color=3;
       // if(iEvent==0) color = 8;
       // if(iEvent==1) color = 9;
       
     }
     
      char bef[] = "before snapshot time";
      char aft[] = "after snapshot time";
      char *when; 
      Int_t sign = 1;
      if( 1000000*(currentEventTime-fCenterTime) > 0.5) {when = aft;
	sign = -1;
	  }
      else when = bef;
      TGraph *temp = new TGraph();
      temp->SetName(Form("temp%d",currentEvent->GetEventNumber()));
       temp->SetLineColor(color);
       temp->SetLineWidth(10);
       
       leg->AddEntry(temp,Form(" %2.2f microseconds %s (global time %2.2f) at z_vertex = %2.2f, centr = ., eventNr: %d", 1000000*sign*currentEventTimeRelToCentral, when, 1000000*currentEvent->GetT0(),  currentEvent->GetZ() ,currentEvent->GetEventNumber() ),"l");
     
       
       DrawEvent(currentEvent, fCenterTime, color);
       

       
       
  
   }
   leg->Draw("same");

}
//________________________________________________________________

void AliToyMCDrawer::DrawEvent(AliToyMCEvent *currentEvent, Double_t centerTime, Int_t color){
  if(!fDispHist) DrawGeometry();
  
  for(Int_t iTrack = 0; iTrack < currentEvent->GetNumberOfTracks(); iTrack++){
    
    //Double_t alpha = TMath::ATan2(currentEvent->GetTrack(iTrack)->GetSpacePoint(0)->GetY(),currentEvent->GetTrack(iTrack)->GetSpacePoint(0)->GetX());
    // if( (     (  !(abs(alpha)>phiCut && abs(alpha) < TMath::Pi() - phiCut)  ) )      )continue;
    // if(alpha<0) {
    // 	continue;
    // 	delete tempTrack;
    // }
    
    //if(currentEvent->GetTrack(iTrack)->GetSpacePoint(0)->GetY()<0) alpha +=TMath::Pi();
    //if( (       (alpha>phiCut && alpha < TMath::Pi() - phiCut) || ( alpha>TMath::Pi() + phiCut && alpha < TMath::TwoPi() - phiCut))   )continue;
    
    // if(alpha > 0) continue;

    //std::cout << TMath::ASin(tempTrack->GetParameter()[2]) << std::endl;
    // Double_t x = currentEvent->GetX();
    // Double_t y = currentEvent->GetY();
      // Double_t hlxpar[6];
      // currentEvent->GetTrack(iTrack)->GetHelixParameters(hlxpar,0);
      // Double_t trackPhi = hlxpar[2];
      // Double_t phiCut = 1.47;
      // std::cout << "bef phi cut "<< trackPhi<< std::endl;
      
      // std::cout << "aft phi cut " << std::endl;
      // Double_t trackEta = currentEvent->GetTrack(iTrack)->GetEta();
      //  Double_t z = currentEvent->GetZ();

     
     
      
      
      const AliToyMCTrack *tempTrack = currentEvent->GetTrack(iTrack);
      //AliToyMCTrack *tempTrack = new AliToyMCTrack(*currentEvent->GetTrack(iTrack));
      DrawTrack(tempTrack, centerTime, currentEvent->GetT0(),color);
      



     

      


       
    }




}

//________________________________________________________________
void AliToyMCDrawer::DrawTrack(const AliToyMCTrack *track,  Double_t centerTime, Double_t currentEventTime, Int_t color){
  if(!fDispHist) DrawGeometry();
  
  TPolyMarker3D *disttrackpoints = new((*fDistPoints)[fDistPoints->GetEntriesFast()]) TPolyMarker3D();
  TPolyMarker3D *trackpoints = new((*fPoints)[fPoints->GetEntriesFast()]) TPolyMarker3D();

  Double_t currentEventTimeRelToCentral = centerTime - currentEventTime;
  Int_t nDistPoints = track->GetNumberOfDistSpacePoints();
  Int_t nPoints = track->GetNumberOfSpacePoints();
    
  for(Int_t iPoint = 0; iPoint< (nDistPoints<nPoints?nPoints:nDistPoints);iPoint++){
    if(iPoint<nPoints) {
	    

      Double_t xp = track->GetSpacePoint(iPoint)->GetX();
      Double_t yp = track->GetSpacePoint(iPoint)->GetY();
      Double_t zp = track->GetSpacePoint(iPoint)->GetZ();
      Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;
      //Double_t zDrifted = (zp/TMath::Abs(zp))*fMaxZ0  -(zp/TMath::Abs(zp))* fDriftVel*(track->GetSpacePoint(iPoint)->GetTimeBin()      -centerTime   );
      Float_t xyzp[3] = {xp,yp,zp};
      AliTrackPoint p;
      p.SetXYZ(xyzp);
      Float_t tempcov[6] = {0};
      p.SetCov(tempcov);
      Int_t sec = track->GetSpacePoint(iPoint)->GetDetector();
      Double_t angle=((sec%18)*20.+10.)/TMath::RadToDeg();
      AliTrackPoint prot = p.Rotate(-angle);
      xp = prot.GetX();
      yp = prot.GetY();


	    
      if(track->GetSpacePoint(iPoint)->GetRow()!=255) {
	if(TMath::Abs(zDrifted)<fMaxZ0 && zDrifted*zp >=0 /*&&TMath::Sqrt(xp*xp + yp*yp)<fIFCRadius*/) trackpoints->SetNextPoint(zDrifted,xp,yp);
	      
      }
      else std::cout << "row == " << track->GetSpacePoint(iPoint)->GetRow() << std::endl;
      
    }
    if(iPoint<nDistPoints) {
      Double_t xpdist = track->GetDistortedSpacePoint(iPoint)->GetX();
      Double_t ypdist = track->GetDistortedSpacePoint(iPoint)->GetY();
      Double_t zpdist = track->GetDistortedSpacePoint(iPoint)->GetZ();
      //std::cout << zpdist << std::endl;
	    
      Float_t xyzpdist[3] = {xpdist,ypdist,zpdist};
      AliTrackPoint pdist;
      pdist.SetXYZ(xyzpdist);
      Float_t tempcovdist[6] = {0};
      pdist.SetCov(tempcovdist);
      Int_t secdist = track->GetDistortedSpacePoint(iPoint)->GetDetector();
      Double_t angledist=((secdist%18)*20.+10.)/TMath::RadToDeg();
      AliTrackPoint protdist = pdist.Rotate(-angledist);
      xpdist = protdist.GetX();
      ypdist = protdist.GetY();
      
      
      UInt_t sector = track->GetDistortedSpacePoint(iPoint)->GetDetector();
      UInt_t row = track->GetDistortedSpacePoint(iPoint)->GetRow();
      UInt_t pad  = track->GetDistortedSpacePoint(iPoint)->GetPad();
      
      Int_t nPads = fTPCParam->GetNPads(sector,row);
      Int_t intPad = TMath::Nint(Float_t(pad+nPads/2));

      Float_t xyz[3];
      //std::cout <<"sector: " << sector << " row: " << row << " pad: " <<pad<< std::endl;
      fRoc->GetPositionGlobal(sector,row,intPad,xyz);
      
      Double_t zDrifteddist = (zpdist/TMath::Abs(zpdist))*fMaxZ0  -(zpdist/TMath::Abs(zpdist))* fDriftVel*(track->GetDistortedSpacePoint(iPoint)->GetTimeBin()      -centerTime   );
      if(row!=255){
	      
	if(TMath::Abs(zDrifteddist)<fMaxZ0 && zDrifteddist*zpdist>=0 /*&&TMath::Sqrt(xpdist*xpdist + ypdist*ypdist)<fIFCRadius*/) {
		
	  disttrackpoints->SetNextPoint(zDrifteddist,xpdist,ypdist);
	  //if(TMath::Sqrt(xpdist*xpdist + ypdist*ypdist)<fIFCRadius) std::cout << "fMaxZ0 " << fMaxZ0 <<" inside " << xpdist << " " << "zpdist "  << zpdist << " " << "zDrifteddist "<< zDrifteddist << " " << zDrifteddist*zpdist << std::endl;
	}
      }
      else std::cout << "row == " << row << std::endl;
      //Double_t zDrifteddist =(zpdist/TMath::Abs(zpdist))*fMaxZ0  -(zpdist/TMath::Abs(zpdist))*(currentEvent->GetTrack(iTrack)->GetDistortedSpacePoint(iPoint)->GetTimeBin()- currentEvent->GetT0() )* fDriftVel;
	  
      
    } 
    //  if( (       (trackPhi>phiCut && trackPhi < TMath::Pi() - phiCut) || ( trackPhi>TMath::Pi() + phiCut && trackPhi < TMath::TwoPi() - phiCut))   ) {
    
    



  }
  if(1){
    if(trackpoints && trackpoints->GetN()>0) {
      //   trackpoints->SetMarkerColor(1+currentEvent->GetEventNumber()%9);
      //trackpoints->SetMarkerStyle(7);
      trackpoints->Draw("same");
    }
    if(disttrackpoints && disttrackpoints->GetN()>0) {
      //  
      
      disttrackpoints->SetMarkerColor(color);

      
      disttrackpoints->Draw("same");
      
    }
    
      }
  
  
  

}

//________________________________________________________________
void AliToyMCDrawer::DrawGeometry() {

  
  //delete fDispGraph;
  //fDispGraph = new TGraph2D();
  delete fDispHist;
 
  fDispHist = new TH3F("fDispHist","",100,-fMaxZ0, fMaxZ0, 100,-(fOFCRadius +10), fOFCRadius +10,100,-(fOFCRadius +10), fOFCRadius +10);
  //if(!fDispGraph) fDispGraph = new TGraph();

  //fDispGraph->Clear();
  fDispHist->SetStats(0);
  fDispHist->GetXaxis()->SetTitle("z [cm]");
  fDispHist->GetYaxis()->SetTitle("x [cm]");
  fDispHist->GetZaxis()->SetTitle("y [cm]");
  fDispHist->Draw();
  gPad->SetPhi(0);
  gPad->SetTheta(0);
  TPolyLine3D *endCap1 = new TPolyLine3D();
  TPolyLine3D *endCap2 = new TPolyLine3D();
  TPolyLine3D *cage[16] ={0x0};
  TPolyLine3D *innerCage[16] ={0x0};
  TPolyLine3D *innerEndCap1 = new TPolyLine3D();
  TPolyLine3D *innerEndCap2 = new TPolyLine3D();
  for(Int_t i = 0; i<16; i++){
    
    cage[i] = new TPolyLine3D();
    cage[i]->SetPoint(0,-fMaxZ0,fOFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    cage[i]->SetPoint(1,fMaxZ0,fOFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    innerCage[i] = new TPolyLine3D();
    innerCage[i]->SetPoint(0,-fMaxZ0,fIFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    innerCage[i]->SetPoint(1,fMaxZ0,fIFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;

    endCap1->SetPoint(i,fMaxZ0,fOFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    endCap2->SetPoint(i,-fMaxZ0,fOFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;

    innerEndCap1->SetPoint(i,fMaxZ0,fIFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    innerEndCap2->SetPoint(i,-fMaxZ0,fIFCRadius*TMath::Cos(i*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(i*TMath::TwoPi()/16)) ;
    innerCage[i]->Draw("same");
    if(!(i%2))  cage[i]->Draw("same");
    
  }
  endCap1->SetPoint(16,fMaxZ0,fOFCRadius*TMath::Cos(16*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(16*TMath::TwoPi()/16)) ;
  endCap2->SetPoint(16,-fMaxZ0,fOFCRadius*TMath::Cos(16*TMath::TwoPi()/16) ,fOFCRadius*TMath::Sin(16*TMath::TwoPi()/16)) ;

  innerEndCap1->SetPoint(16,fMaxZ0,fIFCRadius*TMath::Cos(16*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(16*TMath::TwoPi()/16)) ;
  innerEndCap2->SetPoint(16,-fMaxZ0,fIFCRadius*TMath::Cos(16*TMath::TwoPi()/16) ,fIFCRadius*TMath::Sin(16*TMath::TwoPi()/16)) ;
    

  //fDispGraph->SetTitle("ToyMC display");
  
  endCap1->Draw("same");
  endCap2->Draw("same");

  innerEndCap2->Draw("same");
  innerEndCap1->Draw("same");
  



}

