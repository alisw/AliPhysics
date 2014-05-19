#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TPolyLine3D.h>
#include <TPolyLine.h>
#include <TPolyMarker3D.h>
#include <TPolyMarker.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TPad.h>
#include <TGeoGlobalMagField.h>
#include <TCanvas.h>

#include <AliMagF.h>
#include <AliCDBManager.h>
#include <AliTPCParam.h>
#include <AliGeomManager.h>
#include <AliTPCcalibDB.h>
#include <AliTrackPointArray.h>
#include <AliCluster.h>
#include <AliLog.h>
#include <AliTPCLaserTrack.h>


#include "AliToyMCDrawer.h"
#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"

// Visualization class. To use

//  AliToyMCDrawer* draw = new AliToyMCDrawer()
//  draw->SetFileName("path/to/toyMC.root")

//  draw->FillEventArray(Int_t centerEventNumber)
//         or                 
//  draw->FillEventArray(Double_t time)
//    to display with a certain event in the center or at a certain time 

//  draw->DrawEvents(Bool_t both, Bool_t before)
//    where "both" will display events before and after the middle event and 
//    before will show also events before (after) the middle event if true (false) 
//    when "both" is false


ClassImp(AliToyMCDrawer);

AliToyMCDrawer::AliToyMCDrawer()
  : TObject()   
  ,fInputTree(0x0)
  ,fInFile(0x0)
  ,fFileName()
  ,fEvent(0x0)
  ,fEventArray(0x0)
  ,fDispHist(0x0)
  ,fCenterTime(-1.)
  ,fDriftVel(-1.)
  ,fTPCParam(0x0)
  ,fMaxZ0(0.)
  ,fIFCRadius(83.5)
  ,fOFCRadius(254.5)
  ,fTimeRange(0)
  ,fRoc(AliTPCROC::Instance())
  ,fPoints(0x0)
  ,fDistPoints(0x0)
  ,fProjectionType("XYT")
  ,fTimeZmin(0.)
  ,fTimeZmax(0.)
  ,fGlobalXmin(0.)
  ,fGlobalXmax(0.)
  ,fGlobalYmin(0.)
  ,fGlobalYmax(0.)
  {
   fEventArray = new TClonesArray("AliToyMCEvent");
   
   fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
   fTPCParam->ReadGeoMatrices();
   fDriftVel = fTPCParam->GetDriftV();
   fMaxZ0    =fTPCParam->GetZLength();
   fTimeRange = 2*fMaxZ0/fDriftVel;
 }
//________________________________________________________________
AliToyMCDrawer::AliToyMCDrawer(const AliToyMCDrawer &drawer)
  : TObject(drawer)
  ,fInputTree(drawer.fInputTree)
  ,fInFile(drawer.fInFile)
  ,fFileName(drawer.fFileName)
  ,fEvent(drawer.fEvent)
  ,fEventArray(drawer.fEventArray)
  ,fDispHist(drawer.fDispHist)
  , fCenterTime(drawer.fCenterTime)
  ,fDriftVel(drawer.fDriftVel)
  ,fTPCParam(drawer.fTPCParam)
  ,fMaxZ0(drawer.fMaxZ0)
  ,fIFCRadius(drawer.fIFCRadius)
  ,fOFCRadius(drawer.fOFCRadius)
  ,fTimeRange(drawer.fTimeRange)
  ,fRoc(drawer.fRoc)
  ,fPoints(drawer.fPoints)
  ,fDistPoints(drawer.fDistPoints)
  ,fProjectionType("XYT")
  ,fTimeZmin(0.)
  ,fTimeZmax(0.)
  ,fGlobalXmin(0.)
  ,fGlobalXmax(0.)
  ,fGlobalYmin(0.)
  ,fGlobalYmax(0.)
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
  delete fPoints;
  delete fDistPoints; 
  delete fDispHist;
}
//________________________________________________________________
Int_t AliToyMCDrawer::FillEventArray(Double_t snapShotTime)
{
  if(fFileName.IsNull()) {
    std::cout << "no input file provided, using default (toyMC.root)" << std::endl;
    fFileName = "toyMC.root";
  }
  
  fInFile = new TFile(fFileName.Data(),"read");
  gROOT->cd();
  fInputTree = dynamic_cast<TTree*> (fInFile->Get("toyMCtree"));
  fInputTree->SetBranchAddress("event",&fEvent);
  fEventArray->Clear();
  
  Int_t leftSearchIndex  = Int_t(snapShotTime/2e-5);
  if(leftSearchIndex>fInputTree->GetEntries()) {
    leftSearchIndex = fInputTree->GetEntries()-1;
    fInputTree->GetEvent(leftSearchIndex);
    snapShotTime = fEvent->GetT0();
    std::cout << "input time too large, setting time to time of last event" << std::endl;
  }
  if(leftSearchIndex<0) leftSearchIndex = 0;
  
  //fCenterTime = snapShotTime;

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
      if(!fInputTree->GetEvent(rightSearchIndex)) return FillEventArray(leftSearchIndex);
      rightTime = fEvent->GetT0();
      // std::cout << "bef search " << leftTime << " " << rightTime<< " " << " " << snapShotTime << " "<< (leftTime-snapShotTime)*(rightTime-snapShotTime)  << std::endl;
      // Int_t b;

      while ( (leftTime-snapShotTime)*(rightTime-snapShotTime)>0   ){
	//		printf("searching, leftime %f, righttim %f\n",leftTime,rightTime);
	rightSearchIndex += direction; 
	leftSearchIndex +=direction;
	fInputTree->GetEvent(leftSearchIndex);
	leftTime = fEvent->GetT0();
	if(!fInputTree->GetEvent(rightSearchIndex)) {
	  rightSearchIndex-=direction;
	  break;

	}
	rightTime = fEvent->GetT0();
		printf("searching, leftime %f, righttim %f\n",leftTime,rightTime);
      }
      if (direction==-1) rightSearchIndex = leftSearchIndex;
      firstEventIndex = rightSearchIndex;

    }
  // fInputTree->GetEvent(firstEventIndex);
  //std::cout <<"first event after sn time: " << fEvent->GetT0() << std::endl;
  return FillEventArray(firstEventIndex, snapShotTime);
  
}
//________________________________________________________________
Int_t AliToyMCDrawer::FillEventArray(Int_t middleEventNbr, Double_t snapShotTime)
{
  if(fFileName.IsNull()) {
    std::cout << "no input file provided, using default (toyMC.root)" << std::endl;
    fFileName = "toyMC.root";
  }
  if(!fInFile) fInFile = new TFile(fFileName.Data(),"read");
  gROOT->cd();
  if(!fInputTree) 
    {
      fInputTree = dynamic_cast<TTree*> (fInFile->Get("toyMCtree"));
      //   fInputTree->Add(fileName);
      fInputTree->SetBranchAddress("event",&fEvent);
    }
  Double_t centerEventTime;
  if(fInputTree->GetEvent(middleEventNbr)) 
    {
      centerEventTime = fEvent->GetT0();
      
    }
  else return 0;

  if(snapShotTime<0) fCenterTime = centerEventTime;
  else fCenterTime = snapShotTime;

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

    std::cout << " end fill" << std::endl;

    return 1;
  }
  else {
    printf("Selected event number (%d) out of range!\n",middleEventNbr);
     return 0;
  }
 
}
//________________________________________________________________
void AliToyMCDrawer::DrawEvents(Bool_t both, Bool_t before)
{

  if (fPoints){
    fPoints->Clear();
    fDistPoints->Clear();
  }
  DrawGeometry();
 
  
   Double_t phiCut = 1.;
   TLegend *leg = new TLegend(0.1,0.75,0.6,1,Form("Snapshot of TPC at %f micros after sim start, freq. = , bunchcr.rate = .",1000000*fCenterTime));
   leg->AddEntry((TObject*)0,Form("%2.2f<#phi<#pi-%2.2f && #pi+%2.2f<#phi<2#pi-%2.2f",phiCut,phiCut,phiCut,phiCut),"");

   
 

  for(Int_t iEvent = 0; iEvent<fEventArray->GetEntriesFast(); iEvent++){
 
 
     AliToyMCEvent *currentEvent = static_cast<AliToyMCEvent*> (fEventArray->At(iEvent));
     Double_t currentEventTime = currentEvent->GetT0();
     
     if((1000000*currentEventTime)>((1000000*fCenterTime)+0.05   ))  {
       printf("larger iEvent: %d, current %f, center %f, ntracks: %d\n",iEvent, currentEventTime,fCenterTime, currentEvent->GetNumberOfTracks());
       printf("after\n");
       if ( !(both || !before)) continue;

     }
     if( currentEventTime<=fCenterTime)// && !(currentEventTime==fCenterTime &&both   ))
       { 
	 printf("smaller iEvent: %d, current %f, center %f, ntracks: %d\n",iEvent, currentEventTime,fCenterTime, currentEvent->GetNumberOfTracks());
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
      if (fProjectionType.Length()==3){
        DrawTrack(tempTrack, centerTime, currentEvent->GetT0(),color);
      } else if (fProjectionType.Length()==2){
        DrawTrack2D(tempTrack, centerTime, currentEvent->GetT0(),color);
      }

    }
}

//________________________________________________________________
void AliToyMCDrawer::DrawLaserEvent(Int_t nLaserEvents, Int_t side, Int_t rod, Int_t bundle, Int_t beam)
{
  //
  //
  //

  if (fPoints) {
    fPoints->Clear();
    fDistPoints->Clear();
  }
  if (!ConnectInputTree()) return;

  AliTPCLaserTrack::LoadTracks();
  TObjArray *arr=AliTPCLaserTrack::GetTracks();
  
  if (fProjectionType.Length()==3){
    DrawGeometry();
  } else if (fProjectionType.Length()==2){
    DrawGeometry2D();
  }
  
  Int_t laserEvents=0;
  for (Int_t iev=0;iev<fInputTree->GetEntries();++iev) {
    fInputTree->GetEntry(iev);
    if (fEvent->GetEventType()!=AliToyMCEvent::kLaser) continue;
    if (fEvent->GetNumberOfTracks()!=arr->GetEntriesFast()) {
      AliError(Form("Wrong number of tracks in the laser event: %d!=%d",fEvent->GetNumberOfTracks(),arr->GetEntriesFast()));
      continue;
    }

    for (Int_t iTrack=0; iTrack<fEvent->GetNumberOfTracks();++iTrack){
      AliTPCLaserTrack *laserTrack = (AliTPCLaserTrack*)arr->UncheckedAt(iTrack);
      if (side  >-1 && laserTrack->GetSide()   != side   ) continue;
      if (rod   >-1 && laserTrack->GetRod()    != rod    ) continue;
      if (bundle>-1 && laserTrack->GetBundle() != bundle ) continue;
      if (beam  >-1 && laserTrack->GetBeam()   != beam   ) continue;
      const AliToyMCTrack *track = fEvent->GetTrack(iTrack);
      if (fProjectionType.Length()==3){
        DrawTrack(track,0,fEvent->GetT0(),kRed);
      } else if (fProjectionType.Length()==2){
        DrawTrack2D(track,0,fEvent->GetT0(),kRed);
      }
      
    }

    ++laserEvents;
    if (laserEvents==nLaserEvents) break;
  }
  
  
}

//________________________________________________________________
void AliToyMCDrawer::DrawTrack(const AliToyMCTrack *track,  Double_t centerTime, Double_t currentEventTime, Int_t color){
  if(!fDispHist) DrawGeometry();
  if (!fPoints) {
    fPoints = new TClonesArray("TPolyMarker3D");
    fDistPoints = new TClonesArray("TPolyMarker3D");
  }
  
  TPolyMarker3D *disttrackpoints = new((*fDistPoints)[fDistPoints->GetEntriesFast()]) TPolyMarker3D();
  TPolyMarker3D *trackpoints = new((*fPoints)[fPoints->GetEntriesFast()]) TPolyMarker3D();

  Double_t currentEventTimeRelToCentral = centerTime - currentEventTime;
  Int_t nDistPoints = track->GetNumberOfDistSpacePoints();
  Int_t nPoints = track->GetNumberOfSpacePoints();
    
  for(Int_t iITSPoint = 0; iITSPoint< track->GetNumberOfITSPoints();iITSPoint++){

    Double_t xp = track->GetITSPoint(iITSPoint)->GetX();
    Double_t yp = track->GetITSPoint(iITSPoint)->GetY();
    Double_t zp = track->GetITSPoint(iITSPoint)->GetZ();
    Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;
    trackpoints->SetNextPoint(zDrifted,xp,yp);
  }

  for(Int_t iPoint = 0; iPoint< (nDistPoints<nPoints?nPoints:nDistPoints);iPoint++){
    if(iPoint<nPoints) {
	    

      Double_t xp = track->GetSpacePoint(iPoint)->GetX();
      Double_t yp = track->GetSpacePoint(iPoint)->GetY();
      Double_t zp = track->GetSpacePoint(iPoint)->GetZ();
      Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;
      //Double_t zDrifted = (zp/TMath::Abs(zp))*fMaxZ0  -(zp/TMath::Abs(zp))* fDriftVel*(track->GetSpacePoint(iPoint)->GetTimeBin()      -centerTime   );
      Float_t xyzp[3] = {static_cast<Float_t>(xp),static_cast<Float_t>(yp),static_cast<Float_t>(zp)};
      AliTrackPoint p;
      p.SetXYZ(xyzp);
      Float_t tempcov[6] = {0};
      p.SetCov(tempcov);
      Int_t sec = track->GetSpacePoint(iPoint)->GetDetector();
      Double_t angle=((sec%18)*20.+10.)/TMath::RadToDeg();
      AliTrackPoint prot = p.Rotate(-angle);
      xp = prot.GetX();
      yp = prot.GetY();

//       Double_t ztime=zDrifted;
//       Double_t globx=xp;
//       Double_t globy=yp;

//       if (fProje)

	    
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
	    
      Float_t xyzpdist[3] = {static_cast<Float_t>(xpdist),static_cast<Float_t>(ypdist),static_cast<Float_t>(zpdist)};
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
  if(0){
    for(Int_t iTRDPoint = 0; iTRDPoint< track->GetNumberOfTRDPoints();iTRDPoint++){
      
      Double_t xp = track->GetTRDPoint(iTRDPoint)->GetX();
      Double_t yp = track->GetTRDPoint(iTRDPoint)->GetY();
      Double_t zp = track->GetTRDPoint(iTRDPoint)->GetZ();
      Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;
      trackpoints->SetNextPoint(zDrifted,xp,yp);
  }
}
  if(1){
    if(trackpoints && trackpoints->GetN()>0) {
      //   trackpoints->SetMarkerColor(1+currentEvent->GetEventNumber()%9);
      //trackpoints->SetMarkerStyle(6);
      trackpoints->Draw("same");
    }
    if(disttrackpoints && disttrackpoints->GetN()>0) {
      //  
      //disttrackpoints->SetMarkerStyle(6);
      disttrackpoints->SetMarkerColor(color);

      
      disttrackpoints->Draw("same");
      
    }
    
      }
  
  
  

}

//________________________________________________________________
void AliToyMCDrawer::DrawTrack2D(const AliToyMCTrack *track,  Double_t centerTime, Double_t currentEventTime, Int_t color){
  if(!fDispHist) DrawGeometry2D();
  if (!fPoints) {
    fPoints = new TClonesArray("TPolyMarker");
    fDistPoints = new TClonesArray("TPolyMarker");
  }

  Double_t timeZmin=-fMaxZ0;
  Double_t timeZmax= fMaxZ0;
  Double_t globXmin=-(fOFCRadius +10);
  Double_t globXmax=  fOFCRadius +10 ;
  Double_t globYmin=-(fOFCRadius +10);
  Double_t globYmax=  fOFCRadius +10 ;
  
//   const Double_t epsilon=.001;
  
  TString title;
  if (fTimeZmax>fTimeZmin) {
    timeZmin=fTimeZmin;
    timeZmax=fTimeZmax;
  }
  if (fGlobalXmax>fGlobalXmin) {
    globXmin=fGlobalXmin;
    globXmax=fGlobalXmax;
  }
  if (fGlobalYmax>fGlobalYmin) {
    globYmin=fGlobalYmin;
    globYmax=fGlobalYmax;
  }
  
  TPolyMarker *disttrackpoints = new((*fDistPoints)[fDistPoints->GetEntriesFast()]) TPolyMarker();
  TPolyMarker *trackpoints = new((*fPoints)[fPoints->GetEntriesFast()]) TPolyMarker();
  
  Double_t currentEventTimeRelToCentral = centerTime - currentEventTime;
  Int_t nDistPoints = track->GetNumberOfDistSpacePoints();
  Int_t nPoints = track->GetNumberOfSpacePoints();
  
  for(Int_t iITSPoint = 0; iITSPoint< track->GetNumberOfITSPoints();iITSPoint++){
    
    Double_t xp = track->GetITSPoint(iITSPoint)->GetX();
    Double_t yp = track->GetITSPoint(iITSPoint)->GetY();
    Double_t zp = track->GetITSPoint(iITSPoint)->GetZ();
    Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;

    if ( xp<globXmin || xp>globXmax || yp<globYmin || yp>globYmax || zDrifted<timeZmin || zDrifted>timeZmax ) continue;

    Double_t x=0.;
    Double_t y=0.;
    if (fProjectionType=="XT"){
      x=zDrifted;
      y=xp;
    } else if (fProjectionType=="YT"){
      x=zDrifted;
      y=yp;
    } else if (fProjectionType=="RT"){
      x=zDrifted;
      y=TMath::Sqrt(xp*xp+yp*yp);
    } else if (fProjectionType=="XY" || fProjectionType=="YX") {
      x=xp;
      y=yp;
    }
    trackpoints->SetNextPoint(x,y);
  }
  
  for(Int_t iPoint = 0; iPoint< (nDistPoints<nPoints?nPoints:nDistPoints);iPoint++){
    if(iPoint<nPoints) {
      
      
      Double_t xp = track->GetSpacePoint(iPoint)->GetX();
      Double_t yp = track->GetSpacePoint(iPoint)->GetY();
      Double_t zp = track->GetSpacePoint(iPoint)->GetZ();
      Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;

      //Double_t zDrifted = (zp/TMath::Abs(zp))*fMaxZ0  -(zp/TMath::Abs(zp))* fDriftVel*(track->GetSpacePoint(iPoint)->GetTimeBin()      -centerTime   );
      Float_t xyzp[3] = {static_cast<Float_t>(xp),static_cast<Float_t>(yp),static_cast<Float_t>(zp)};
      AliTrackPoint p;
      p.SetXYZ(xyzp);
      Float_t tempcov[6] = {0};
      p.SetCov(tempcov);
      Int_t sec = track->GetSpacePoint(iPoint)->GetDetector();
      Double_t angle=((sec%18)*20.+10.)/TMath::RadToDeg();
      AliTrackPoint prot = p.Rotate(-angle);
      xp = prot.GetX();
      yp = prot.GetY();

      if ( xp<globXmin || xp>globXmax || yp<globYmin || yp>globYmax || zDrifted<timeZmin || zDrifted>timeZmax ) continue;
      
//       Double_t ztime=zDrifted;
//       Double_t globx=xp;
//       Double_t globy=yp;
      
      Double_t x=0.;
      Double_t y=0.;
      if (fProjectionType=="XT"){
        x=zDrifted;
        y=xp;
      } else if (fProjectionType=="YT"){
        x=zDrifted;
        y=yp;
      } else if (fProjectionType=="RT"){
        x=zDrifted;
        y=TMath::Sqrt(xp*xp+yp*yp);
      } else if (fProjectionType=="XY" || fProjectionType=="YX") {
        x=xp;
        y=yp;
      }
      
      if(track->GetSpacePoint(iPoint)->GetRow()!=255) {
        if(TMath::Abs(zDrifted)<fMaxZ0 && zDrifted*zp >=0 /*&&TMath::Sqrt(xp*xp + yp*yp)<fIFCRadius*/) trackpoints->SetNextPoint(x,y);
        
      }
      else std::cout << "row == " << track->GetSpacePoint(iPoint)->GetRow() << std::endl;
      
    }
    if(iPoint<nDistPoints) {
      Double_t xpdist = track->GetDistortedSpacePoint(iPoint)->GetX();
      Double_t ypdist = track->GetDistortedSpacePoint(iPoint)->GetY();
      Double_t zpdist = track->GetDistortedSpacePoint(iPoint)->GetZ();
      
      //std::cout << zpdist << std::endl;
      
      Float_t xyzpdist[3] = {static_cast<Float_t>(xpdist),static_cast<Float_t>(ypdist),static_cast<Float_t>(zpdist)};
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

      if ( xpdist<globXmin || xpdist>globXmax || ypdist<globYmin || ypdist>globYmax || zDrifteddist<timeZmin || zDrifteddist>timeZmax ) continue;
      
      Double_t x=0.;
      Double_t y=0.;
      if (fProjectionType=="XT"){
        x=zDrifteddist;
        y=xpdist;
      } else if (fProjectionType=="YT"){
        x=zDrifteddist;
        y=ypdist;
      } else if (fProjectionType=="RT"){
        x=zDrifteddist;
        y=TMath::Sqrt(xpdist*xpdist+ypdist*ypdist);
      } else if (fProjectionType=="XY" || fProjectionType=="YX") {
        x=xpdist;
        y=ypdist;
      }
      
      if(row!=255){
        
        if(TMath::Abs(zDrifteddist)<fMaxZ0 && zDrifteddist*zpdist>=0 /*&&TMath::Sqrt(xpdist*xpdist + ypdist*ypdist)<fIFCRadius*/) {
          
          disttrackpoints->SetNextPoint(x,y);
          //if(TMath::Sqrt(xpdist*xpdist + ypdist*ypdist)<fIFCRadius) std::cout << "fMaxZ0 " << fMaxZ0 <<" inside " << xpdist << " " << "zpdist "  << zpdist << " " << "zDrifteddist "<< zDrifteddist << " " << zDrifteddist*zpdist << std::endl;
        }
      }
      else std::cout << "row == " << row << std::endl;
      //Double_t zDrifteddist =(zpdist/TMath::Abs(zpdist))*fMaxZ0  -(zpdist/TMath::Abs(zpdist))*(currentEvent->GetTrack(iTrack)->GetDistortedSpacePoint(iPoint)->GetTimeBin()- currentEvent->GetT0() )* fDriftVel;
      
      
    }
    //  if( (       (trackPhi>phiCut && trackPhi < TMath::Pi() - phiCut) || ( trackPhi>TMath::Pi() + phiCut && trackPhi < TMath::TwoPi() - phiCut))   ) {
      
      
      
      
      
  }
  if(0){
    for(Int_t iTRDPoint = 0; iTRDPoint< track->GetNumberOfTRDPoints();iTRDPoint++){
      
      Double_t xp = track->GetTRDPoint(iTRDPoint)->GetX();
      Double_t yp = track->GetTRDPoint(iTRDPoint)->GetY();
      Double_t zp = track->GetTRDPoint(iTRDPoint)->GetZ();
      Double_t zDrifted =  zp+(zp/TMath::Abs(zp))*currentEventTimeRelToCentral * fDriftVel;

      if ( xp<globXmin || xp>globXmax || yp<globYmin || yp>globYmax || zDrifted<timeZmin || zDrifted>timeZmax ) continue;
      
      Double_t x=0.;
      Double_t y=0.;
      if (fProjectionType=="XT"){
        x=zDrifted;
        y=xp;
      } else if (fProjectionType=="YT"){
        x=zDrifted;
        y=yp;
      } else if (fProjectionType=="RT"){
        x=zDrifted;
        y=TMath::Sqrt(xp*xp+yp*yp);
      } else if (fProjectionType=="XY" || fProjectionType=="YX") {
        x=xp;
        y=yp;
      }
      trackpoints->SetNextPoint(x,y);
    }
  }
  if(1){
    if(trackpoints && trackpoints->GetN()>0) {
      //   trackpoints->SetMarkerColor(1+currentEvent->GetEventNumber()%9);
      //trackpoints->SetMarkerStyle(6);
      trackpoints->Draw("same");
    }
    if(disttrackpoints && disttrackpoints->GetN()>0) {
      //
      //disttrackpoints->SetMarkerStyle(6);
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

  Double_t timeZmin=-fMaxZ0;
  Double_t timeZmax= fMaxZ0;
  Double_t globXmin=-(fOFCRadius +10);
  Double_t globXmax=  fOFCRadius +10 ;
  Double_t globYmin=-(fOFCRadius +10);
  Double_t globYmax=  fOFCRadius +10 ;

  const Double_t epsilon=.001;
  
  TString title;
  if (fTimeZmax>fTimeZmin) {
    timeZmin=fTimeZmin;
    timeZmax=fTimeZmax;
  }
  if (fGlobalXmax>fGlobalXmin) {
    globXmin=fGlobalXmin;
    globXmax=fGlobalXmax;
  }
  if (fGlobalYmax>fGlobalYmin) {
    globYmin=fGlobalYmin;
    globYmax=fGlobalYmax;
  }
  fDispHist = new TH3F("fDispHist",";#it{z} (cm); #it{x} (cm); #it{y} (cm)",
                       100, timeZmin-10, timeZmax+10,
                       100, globXmin, globXmax ,
                       100, globYmin, globYmax);
    
  //if(!fDispGraph) fDispGraph = new TGraph();

  //fDispGraph->Clear();
  fDispHist->SetStats(0);
  fDispHist->Draw();
  gPad->SetPhi(0);
  gPad->SetTheta(0);
  
  TPolyLine3D *endCap1 = 0x0;
  if (timeZmin-epsilon<-fMaxZ0)
    endCap1=new TPolyLine3D();
  printf("time: %p, %.2f,, %.2f, %.2f, %.2f, %d\n",endCap1,fTimeZmin,fTimeZmax,timeZmin-epsilon,fMaxZ0, timeZmin-epsilon<-fMaxZ0);
  TPolyLine3D *endCap2 = 0x0;
  if (timeZmax+epsilon> fMaxZ0)
    endCap2=new TPolyLine3D();

  TPolyLine3D *outerCE = 0x0;
  if (timeZmin<0 && timeZmax>0 )
    outerCE=new TPolyLine3D();
  
  TPolyLine3D *cage[18] ={0x0};
  
  TPolyLine3D *innerCage[18] ={0x0};
  
  TPolyLine3D *innerEndCap1 = 0x0;
  if (timeZmin-epsilon<-fMaxZ0)
    innerEndCap1=new TPolyLine3D();
  
  TPolyLine3D *innerEndCap2 = 0x0;
  if (timeZmax+epsilon> fMaxZ0)
    innerEndCap2=new TPolyLine3D();

  TPolyLine3D *innerCE = 0x0;
  if (timeZmin<0 && timeZmax>0 )
    innerCE=new TPolyLine3D();
  
  Int_t iPoint=0;
  Double_t angle    = 0.;
  Double_t globalX  = 0.;
  Double_t globalY  = 0.;
  Double_t globalXi = 0.;
  Double_t globalYi = 0.;
  
  for(Int_t i = 0; i<18; i++){
    angle    = i*TMath::TwoPi()/18;
    globalX  = fOFCRadius*TMath::Cos(angle);
    globalY  = fOFCRadius*TMath::Sin(angle);
    globalXi = fIFCRadius*TMath::Cos(angle);
    globalYi = fIFCRadius*TMath::Sin(angle);
    
    cage[iPoint] = new TPolyLine3D();
    cage[iPoint]->SetPoint(0,timeZmin,globalX ,globalY) ;
    cage[iPoint]->SetPoint(1, timeZmax,globalX ,globalY) ;
    innerCage[iPoint] = new TPolyLine3D();
    innerCage[iPoint]->SetPoint(0,timeZmin,globalXi ,globalYi) ;
    innerCage[iPoint]->SetPoint(1, timeZmax,globalXi ,globalYi) ;

    // only draw if inside range
    if (endCap1) { endCap1->SetPoint(i,timeZmax, globalX, globalY); }
    if (endCap2) { endCap2->SetPoint(i,timeZmin, globalX, globalY); }
    if (outerCE) { outerCE->SetPoint(i,      0., globalX, globalY); }
    
    if (innerEndCap1) { innerEndCap1->SetPoint(i, timeZmax, globalXi, globalYi); }
    if (innerEndCap2) { innerEndCap2->SetPoint(i, timeZmin, globalXi, globalYi); }
    if (innerCE)      {      innerCE->SetPoint(i,       0., globalXi, globalYi); }
    
    innerCage[iPoint]->Draw("same");
    
    if(!(i%2))
      cage[iPoint]->Draw("same");

    ++iPoint;
  }

  //
  // close endplate and CE polygons
  //
  Int_t i=18;
  angle    = i*TMath::TwoPi()/18;
  globalX  = fOFCRadius*TMath::Cos(angle);
  globalY  = fOFCRadius*TMath::Sin(angle);
  globalXi = fIFCRadius*TMath::Cos(angle);
  globalYi = fIFCRadius*TMath::Sin(angle);
  
  // only draw if inside range
  if (endCap1) { endCap1->SetPoint(i, timeZmax, globalX, globalY); }
  if (endCap2) { endCap2->SetPoint(i, timeZmin, globalX, globalY); }
  if (outerCE) { outerCE->SetPoint(i,       0., globalX, globalY); }
  
  if (innerEndCap1) { innerEndCap1->SetPoint(i, timeZmax, globalXi, globalYi); }
  if (innerEndCap2) { innerEndCap2->SetPoint(i, timeZmin, globalXi, globalYi); }
  if (innerCE)      { innerCE     ->SetPoint(i,       0., globalXi, globalYi); }
  
  if (endCap1) { endCap1->Draw("same"); }
  if (endCap2) { endCap2->Draw("same"); }
  if (outerCE) { outerCE->Draw("same"); }

  if (innerEndCap1) { innerEndCap1->Draw("same"); }
  if (innerEndCap2) { innerEndCap2->Draw("same"); }
  if (innerCE)      { innerCE     ->Draw("same");      }
  
}

//________________________________________________________________
void AliToyMCDrawer::DrawGeometry2D() {
  
  
  //delete fDispGraph;
  //fDispGraph = new TGraph2D();
  delete fDispHist;
  
  Double_t timeZmin=-fMaxZ0;
  Double_t timeZmax= fMaxZ0;
  Double_t globXmin=-(fOFCRadius +10);
  Double_t globXmax=  fOFCRadius +10 ;
  Double_t globYmin=-(fOFCRadius +10);
  Double_t globYmax=  fOFCRadius +10 ;
  
  const Double_t epsilon=.001;
  
  TString title;
  if (fTimeZmax>fTimeZmin) {
    timeZmin=fTimeZmin;
    timeZmax=fTimeZmax;
  }
  if (fGlobalXmax>fGlobalXmin) {
    globXmin=fGlobalXmin;
    globXmax=fGlobalXmax;
  }
  if (fGlobalYmax>fGlobalYmin) {
    globYmin=fGlobalYmin;
    globYmax=fGlobalYmax;
  }
  TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("cDrawer");

  if (fProjectionType=="XT"){
    if (!c) c=new TCanvas("cDrawer","Toy Drawer");
    fDispHist = new TH2F("fDispHist",";#it{z} (cm); #it{x} (cm)",
                         100, timeZmin-10, timeZmax+10,
                         100, globXmin, globXmax);
  } else if (fProjectionType=="YT"){
    if (!c) c=new TCanvas("cDrawer","Toy Drawer");
    fDispHist = new TH2F("fDispHist",";#it{z} (cm); #it{y} (cm)",
                         100, timeZmin-10, timeZmax+10,
                         100, globYmin, globYmax);
  } else if (fProjectionType=="RT"){
    if (!c) c=new TCanvas("cDrawer","Toy Drawer");
    fDispHist = new TH2F("fDispHist",";#it{z} (cm); #it{r} (cm)",
                         100, timeZmin-10, timeZmax+10,
                         100, globYmin, globYmax);
  } else if (fProjectionType=="YX"||fProjectionType=="XY"){
    if (!c) c=new TCanvas("cDrawer","Toy Drawer",800,800);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.05);
    fDispHist = new TH2F("fDispHist",";#it{x} (cm); #it{y} (cm)",
                         100, globXmin, globXmax,
                         100, globYmin, globYmax);
    fDispHist->GetYaxis()->SetTitleOffset(1.5);
  } else {
    AliError(Form("Display Format not known: %s",fProjectionType.Data()));
    return;
  }
  
  c->Clear();
  
  //if(!fDispGraph) fDispGraph = new TGraph();
  
  //fDispGraph->Clear();
  fDispHist->SetStats(0);
  fDispHist->Draw();
  gPad->SetPhi(0);
  gPad->SetTheta(0);
  
  TPolyLine *endCap1 = 0x0;
  if (timeZmin-epsilon<-fMaxZ0)
    endCap1=new TPolyLine();
  printf("time: %p, %.2f,, %.2f, %.2f, %.2f, %d\n",endCap1,fTimeZmin,fTimeZmax,timeZmin-epsilon,fMaxZ0, timeZmin-epsilon<-fMaxZ0);
  TPolyLine *endCap2 = 0x0;
  if (timeZmax+epsilon> fMaxZ0)
    endCap2=new TPolyLine();
  
  TPolyLine *outerCE = 0x0;
  if (timeZmin<0 && timeZmax>0 )
    outerCE=new TPolyLine();
  
  TPolyLine *cage[18] ={0x0};
  
  TPolyLine *innerCage[18] ={0x0};
  
  TPolyLine *innerEndCap1 = 0x0;
  if (timeZmin-epsilon<-fMaxZ0)
    innerEndCap1=new TPolyLine();
  
  TPolyLine *innerEndCap2 = 0x0;
  if (timeZmax+epsilon> fMaxZ0)
    innerEndCap2=new TPolyLine();
  
  TPolyLine *innerCE = 0x0;
  if (timeZmin<0 && timeZmax>0 )
    innerCE=new TPolyLine();
  
  Int_t iPoint=0;
  Double_t angle    = 0.;
  Double_t globalX  = 0.;
  Double_t globalY  = 0.;
  Double_t globalXi = 0.;
  Double_t globalYi = 0.;
  Double_t globalXi2= 0.;
  Double_t globalYi2= 0.;

  
  if (fProjectionType=="YX"||fProjectionType=="XY") {
    for(Int_t i = 0; i<18; i++){
      angle    = i*TMath::TwoPi()/18;
      globalX  = fOFCRadius*TMath::Cos(angle);
      globalY  = fOFCRadius*TMath::Sin(angle);
      globalXi = fIFCRadius*TMath::Cos(angle);
      globalYi = fIFCRadius*TMath::Sin(angle);
      globalXi = fIFCRadius*TMath::Cos(angle);
      globalYi = fIFCRadius*TMath::Sin(angle);
      globalXi2= 135.5*TMath::Cos(angle);
      globalYi2= 135.5*TMath::Sin(angle);
      
      if ( globalX<globXmin || globalX>globXmax || globalY<globYmin || globalY>globYmax ) {
        continue;
      }
      // in xy, abuse for sector boundaries
      cage[iPoint] = new TPolyLine();
      cage[iPoint]->SetPoint(0, globalX ,globalY) ;
      cage[iPoint]->SetPoint(1, globalXi ,globalYi) ;
      // only draw if inside range
      if (endCap1) { endCap1->SetPoint(iPoint, globalX, globalY); }

      if (innerEndCap1) { innerEndCap1->SetPoint(iPoint, globalXi, globalYi); }
      if (innerEndCap2) { innerEndCap2->SetPoint(iPoint, globalXi2, globalYi2); }
      
      cage[iPoint]->Draw("same");

      ++iPoint;
    }
  }
  
  angle    = 5*TMath::TwoPi()/18;
  globalX  = fOFCRadius*TMath::Cos(angle);
  globalY  = fOFCRadius*TMath::Sin(angle);
  globalXi = fIFCRadius*TMath::Cos(angle);
  globalYi = fIFCRadius*TMath::Sin(angle);
  globalXi = fIFCRadius*TMath::Cos(angle);
  globalYi = fIFCRadius*TMath::Sin(angle);
  globalXi2= 135.5*TMath::Cos(angle);
  globalYi2= 135.5*TMath::Sin(angle);
  if (endCap1) { endCap1->SetPoint(iPoint, 0, globalY); }
  
  if (innerEndCap1) { innerEndCap1->SetPoint(iPoint, 0, globalYi); }
  if (innerEndCap2) { innerEndCap2->SetPoint(iPoint, 0, globalYi2); }
  ++iPoint;
  //
  // close endplate and CE polygons
  //
  Int_t i=18;
  angle    = i*TMath::TwoPi()/18;
  globalX  = fOFCRadius*TMath::Cos(angle);
  globalY  = fOFCRadius*TMath::Sin(angle);
  globalXi = fIFCRadius*TMath::Cos(angle);
  globalYi = fIFCRadius*TMath::Sin(angle);
  globalXi2= 135.5*TMath::Cos(angle);
  globalYi2= 135.5*TMath::Sin(angle);
  
  // only draw if inside range
  if ( !(globalX<globXmin || globalX>globXmax || globalY<globYmin || globalY>globYmax) ) {
    if (endCap1) { endCap1->SetPoint(i, globalX, globalY); }
//     if (endCap2) { endCap2->SetPoint(i, globalX, globalY); }
//     if (outerCE) { outerCE->SetPoint(i, globalX, globalY); }

    if (innerEndCap1) { innerEndCap1->SetPoint(i, globalXi, globalYi); }
    if (innerEndCap2) { innerEndCap2->SetPoint(i, globalXi2, globalYi2); }
//     if (innerCE)      { innerCE     ->SetPoint(i, globalXi, globalYi); }

  }
  
  if (endCap1) { endCap1->Draw("same"); }
  //     if (endCap2) { endCap2->Draw("same"); }
  //     if (outerCE) { outerCE->Draw("same"); }
  
  if (innerEndCap1) { innerEndCap1->Draw("same"); }
  if (innerEndCap2) { innerEndCap2->Draw("same"); }
  //     if (innerCE)      { innerCE     ->Draw("same");      }
}

//________________________________________________________________
Bool_t AliToyMCDrawer::ConnectInputTree()
{
  //
  //
  //

  if(fFileName.IsNull()) {
    AliError("no input file provided, using default (toyMC.root)");
    fFileName = "toyMC.root";
  }
  if (fInFile && fInFile->GetName()==fFileName) return kTRUE;
  delete fInFile;
  fInFile=0x0;
  delete fInputTree;
  fInputTree=0x0;
  
  fInFile = new TFile(fFileName.Data(),"read");
  if (!fInFile || !fInFile->IsOpen() || fInFile->IsZombie() ) {
    delete fInFile;
    return kFALSE;
  }
  
  gROOT->cd();
  
  fInputTree = dynamic_cast<TTree*> (fInFile->Get("toyMCtree"));

  if (!fInputTree){
    AliError("Could not connect tree!\n");
    delete fInFile;
    fInFile=0x0;
    return kFALSE;
  }
  
  fInputTree->SetBranchAddress("event",&fEvent);
  return kTRUE;
}
