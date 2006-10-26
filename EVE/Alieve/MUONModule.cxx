/* HEAD11Jul06 */
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The main AliEVE drawing module for the MUON detector                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "EventAlieve.h"
#include "Reve/PointSet.h"
#include "Reve/RGTopFrame.h"

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliMUON.h"
#include "AliMpSegmentation.h"
#include "AliMpDEIterator.h"
#include "AliMpSectorSegmentation.h"
#include "AliMpSector.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONSegmentation.h"
#include "AliMpStationType.h"
#include "AliMpDEManager.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"

#include "MUONModule.h"

#include <TPolyLine3D.h>
#include <TMarker3DBox.h>
#include <TColor.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(MUONModule)

/**************************************************************************/
MUONModule::MUONModule(const Text_t* n, const Text_t* t, Color_t col) :
  Reve::RenderElement(col),
  QuadSet(n, t),
  fInfo(0),
  fID(-1), 
  fCath(0),
  fShowDigits(0), fShowClusters(0), fShowTracks(0),
  fFrameCol(col),
  fDetElemId(-1)
{
  //
  // Default constructor
  //

}

/**************************************************************************/
MUONModule::MUONModule(Int_t id, Int_t cath, MUONDigitsInfo* info, Bool_t dig, Bool_t clus, Color_t col ) :
  Reve::RenderElement(col),
  QuadSet(Form("M-DetElemId %d C%1d",id,cath)),
  fInfo(info),
  fID(-1), 
  fCath(0),
  fShowDigits(dig), fShowClusters(clus), fShowTracks(0),
  fFrameCol(col),
  fDetElemId(-1)
{
  //
  // Constructor
  //

  if (!fShowDigits && !fShowClusters) fShowTracks = 1;

  if (fShowClusters) SetName(Form("M-DetElemId %d",id));

  if (id/100 >= 11 && cath == 1) SetName(Form("M-DetElemId %d X",id));
  if (id/100 >= 11 && cath == 2) SetName(Form("M-DetElemId %d Y",id));
  if (id == -1 && cath == -1) SetName(Form("M-Chambers"));
  if (id >  -1 && cath == -1) SetName(Form("M-Track %d",id));
  if (cath == -2) SetName(Form("M-Chamber %2d",id));
  if (cath == -3) SetName(Form("M-Pads %d X",id));
  if (cath == -4) SetName(Form("M-Pads %d Y",id));

  SetID(id,cath);

}

/**************************************************************************/
MUONModule::MUONModule(const MUONModule &mmod) :
  Reve::RenderElement(),
  QuadSet(Form("M-DetElemId %d C%1d",mmod.fID,mmod.fCath)),
  fInfo(mmod.fInfo),
  fID(mmod.fID),
  fCath(mmod.fCath),
  fShowDigits(mmod.fShowDigits),
  fShowClusters(mmod.fShowClusters),
  fShowTracks(mmod.fShowTracks),
  fFrameCol(mmod.fFrameCol),
  fDetElemId(mmod.fDetElemId)
{
  //
  // Copy constructor
  //

}

/**************************************************************************/
MUONModule& MUONModule::operator=(const MUONModule &mmod)
{
  //
  // Assignment operator
  //

  if (this != &mmod) {

    fInfo         = mmod.fInfo;
    fID           = mmod.fID;
    fCath         = mmod.fCath;
    fShowDigits   = mmod.fShowDigits;
    fShowClusters = mmod.fShowClusters;
    fShowTracks   = mmod.fShowTracks;
    fFrameCol     = mmod.fFrameCol;
    fDetElemId    = mmod.fDetElemId;

  }

  return *this;

}

/**************************************************************************/
void MUONModule::SetID(Int_t id, Int_t cath)
{
  //
  // Select a single detector element id
  //

  static const Exc_t eH("MUOModule::SetID ");

  if(fInfo == 0)
    throw(eH + "MUONDigitsInfo not set.");

  fID   = id;
  fCath = cath;

  InitModule();

}

/**************************************************************************/
void MUONModule::InitModule()
{
  // 
  // Initialize and draw selected items
  //

  fDetElemId = 0;

  if (fShowDigits) LoadQuadsDigits();
  if (fShowClusters) LoadQuadsClusters();
  if (fShowTracks) {
    if (fCath == -1){
      LoadQuadsTracks(fID);
    } else if (fCath == -2 || fCath == -3 || fCath == -4) {
      if (fCath == -2) LoadQuadsChambers(fID,fID);
      if (fCath == -3) LoadQuadsChambers(fID,fID,-1,1);
      if (fCath == -4) LoadQuadsChambers(fID,fID,-1,2);
    }
  }
  ComputeBBox();

}

/**************************************************************************/
void MUONModule::LoadQuadsChambers(Int_t chamber1, Int_t chamber2, Int_t id, Int_t cat)
{
  //
  // Draw chambers
  //

  //printf("Draw chambers: %d %d %d %d %d \n",chamber1,chamber2,id,cat,fCath);

  Int_t fChamber;

  AliRunLoader *runLoader = Alieve::Event::AssertRunLoader();
  runLoader->LoadgAlice();

  AliMUON *pMUON = (AliMUON*)gAlice->GetModule("MUON");

  const AliMUONGeometryTransformer* kGeomTransformer = pMUON->GetGeometryTransformer();
  AliMUONSegmentation* segmentation = pMUON->GetSegmentation();

  Float_t xg1, xg2, yg1, yg2, zg1, zg2;

  for (fChamber = chamber1; fChamber <= chamber2; fChamber++) {

  if(fChamber < 5) {

    AliMpDEIterator it;
    for ( it.First(fChamber-1); ! it.IsDone(); it.Next() ) {
      
    Int_t detElemId = it.CurrentDE();

    if (id > 0 && id != detElemId) continue;
      
    //printf("Detector element ID for chamber %d (tracking): %d \n",fChamber,detElemId);
    
    AliMpSectorSegmentation * seg 
      = (AliMpSectorSegmentation *) AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, 0);
    const AliMpSector * sector = seg->GetSector();
    
    // get sector measurements
    TVector2 position  = sector->Position(); 
    TVector2 dimension = sector->Dimensions(); // half length

    //printf("Sector position: \n"); position.Print();
    //printf("Sector dimensions: \n"); dimension.Print();
    
    Float_t xlocal1 =  position.Px(); // FIXME: not really needed as it's 0 ?
    Float_t ylocal1 =  position.Py(); // FIXME: not really needed as it's 0 ?
    Float_t xlocal2 =  dimension.Px() * 2.;
    Float_t ylocal2 =  dimension.Px() * 2.;
    
    //printf("Local position and dimension: xp = %f , yp = %f xd = %f yd = %f \n",xlocal1,ylocal1,xlocal2,ylocal2);
    
    kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
    kGeomTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg2, yg2, zg2);
    
    //printf("Global position: xpg = %f , ypg = %f zpg = %f \n",xg1,yg1,zg1);
    //printf("Global dimension: xdg = %f , ydg = %f zdg = %f \n",xg2,yg2,zg2);
    
    Float_t *p;
    if (fCath != -3 && fCath != -4) {
    fQuads.push_back(Reve::Quad(fFrameCol));
    p = fQuads.back().vertices;
    /*
    p[0] = xg1;      p[1] =  yg1; p[2]  = zg1;
    p[3] = xg1;      p[4] =  yg1; p[5]  = zg1;
    p[6] = xg1+xg2;  p[7] =  yg1; p[8]  = zg1;
    p[9] = xg1+xg2;  p[10] = yg1; p[11] = zg1;
    */
    // switch x <> z
    p[0] = zg1;      p[1] =  yg1; p[2]  = xg1;
    p[3] = zg1;      p[4] =  yg1; p[5]  = xg1;
    p[6] = zg1;      p[7] =  yg1; p[8]  = xg1+xg2;
    p[9] = zg1;      p[10] = yg1; p[11] = xg1+xg2;
    
    Float_t xprev = xg1+xg2;
    Float_t yprev = yg1;
    Int_t nstep = 100;
    Float_t dstep = TMath::Pi()/2.0 / (Float_t)nstep;
    Float_t d;
    for (Int_t istep = 1; istep < (nstep+1); istep++) {

      d = istep * dstep;
      Float_t x = xg1 + xg2 * TMath::Cos(d);
      Float_t y = yg1 + yg2 * TMath::Sin(d);

      fQuads.push_back(Reve::Quad(fFrameCol));
      p = fQuads.back().vertices;
      /*
      p[0] = xprev; p[1] =  yprev; p[2]  = zg1;
      p[3] = xprev; p[4] =  yprev; p[5]  = zg1;
      p[6] = x;     p[7] =  y;     p[8]  = zg1;
      p[9] = x;     p[10] = y;     p[11] = zg1;
      */
      // switch x <> z
      p[0] = zg1;   p[1] =  yprev; p[2]  = xprev;
      p[3] = zg1;   p[4] =  yprev; p[5]  = xprev;
      p[6] = zg1;   p[7] =  y;     p[8]  = x;
      p[9] = zg1;   p[10] = y;     p[11] = x;

      xprev = x;
      yprev = y;

    }

    fQuads.push_back(Reve::Quad(fFrameCol));
    p = fQuads.back().vertices;
    /*
    p[0] = xprev;      p[1] =  yprev; p[2]  = zg1;
    p[3] = xprev;      p[4] =  yprev; p[5]  = zg1;
    p[6] = xg1;        p[7] =  yg1;   p[8]  = zg1;
    p[9] = xg1;        p[10] = yg1;   p[11] = zg1;
    */
    // switch x <> z
    p[0] = zg1;        p[1] =  yprev; p[2]  = xprev;
    p[3] = zg1;        p[4] =  yprev; p[5]  = xprev;
    p[6] = zg1;        p[7] =  yg1;   p[8]  = xg1;
    p[9] = zg1;        p[10] = yg1;   p[11] = xg1;
    }  //  end fCath != -3 , -4
    }  //  end detElemId

  }

  if (fChamber>4) {

    AliMpDEIterator it;
    for ( it.First(fChamber-1); ! it.IsDone(); it.Next() ) {
      
    Int_t detElemId = it.CurrentDE();
          
    if (id > 0 && id != detElemId) continue;
      
    //AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);
    
    //printf("Detector element ID for chamber %d (trigger, stationType %s): %d \n",fChamber,StationTypeName(stationType).Data(),detElemId);
    
    if (  segmentation->HasDE(detElemId) ) {
      const AliMpVSegmentation* seg1 
        = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, 0);
      if (!seg1) { 
	// Create mapping segmentation if old trigger segmentation
	// (not using mapping)
	seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, 0);
      }	  
      if (seg1) {  
	
	Float_t deltax = seg1->Dimensions().X();
	Float_t deltay = seg1->Dimensions().Y();
	Float_t xlocal1 =  -deltax;
	Float_t ylocal1 =  -deltay;
	Float_t xlocal2 =  +deltax;
	Float_t ylocal2 =  +deltay;
	
	//printf("Local corners: x1 = %f , y1 = %f x2 = %f y2 = %f \n",xlocal1,ylocal1,xlocal2,ylocal2);
	
	kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	kGeomTransformer->Local2Global(detElemId, xlocal2, ylocal2, 0, xg2, yg2, zg2);
	
	//printf("Global corner 1: xg1 = %f , yg1 = %f zg1 = %f \n",xg1,yg1,zg1);
	//printf("Global corner 2: xg2 = %f , yg2 = %f zg2 = %f \n",xg2,yg2,zg2);
	
	if (fCath != -3 && fCath != -4) {
	fQuads.push_back(Reve::Quad(fFrameCol));
	Float_t* p = fQuads.back().vertices;
	// the tilting direction now known...
	/*
	p[0] = xg1;  p[1] =  yg1; p[2]  = zg1;
	p[3] = xg1;  p[4] =  yg2; p[5]  = zg1;
	p[6] = xg2;  p[7] =  yg2; p[8]  = zg2;
	p[9] = xg2;  p[10] = yg1; p[11] = zg2;
	*/
	// switch x <> z
	p[0] = zg1;  p[1] =  yg1; p[2]  = xg1;
	p[3] = zg1;  p[4] =  yg2; p[5]  = xg1;
	p[6] = zg2;  p[7] =  yg2; p[8]  = xg2;
	p[9] = zg2;  p[10] = yg1; p[11] = xg2;
	}
	if (fCath == -3 || fCath == -4) {
	// draw all pads in the trigger chambers
	if (fChamber > 10) {

	  for (Int_t ic = 1; ic <= 2; ic++) {

	  if (ic != cat) continue;

	  const AliMpVSegmentation* seg3 
	    = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId, ic-1);
	  Int_t maxX = seg3->MaxPadIndexX();
	  Int_t maxY = seg3->MaxPadIndexY();
	  //printf("detElemId %d ic %d maxX %d maxY %d \n",detElemId,ic,maxX,maxY);
	  for (Int_t ix = 0; ix <= maxX; ix++) {
	    for (Int_t iy = 0; iy <= maxY; iy++) {

	      AliMpPad pad = seg3->PadByIndices(AliMpIntPair(ix,iy),kFALSE);

	      if (!pad.IsValid()) continue;

	      // get the pad position and dimensions
	      xlocal1 = pad.Position().X();
	      ylocal1 = pad.Position().Y();
	      xlocal2 = pad.Dimensions().X();
	      ylocal2 = pad.Dimensions().Y();

	      kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
	      // (no transformation for pad dimensions)
	      xg2 = xlocal2;
	      yg2 = ylocal2;
	      
	      Int_t pcolor = 0;
	      if (ic == 1) pcolor = 5;
	      if (ic == 2) pcolor = 5;

	      fQuads.push_back(Reve::Quad());
	      fQuads.back().ColorFromIdx(pcolor);
	      Float_t* p = fQuads.back().vertices;
	      // the tilting direction now known...
	      // switch x <> z
	      p[0] = zg1;         p[1] =  yg1-yg2; p[2]  = xg1-xg2;
	      p[3] = zg1;         p[4] =  yg1+yg2; p[5]  = xg1-xg2;
	      p[6] = zg1;         p[7] =  yg1+yg2; p[8]  = xg1+xg2;
	      p[9] = zg1;         p[10] = yg1-yg2; p[11] = xg1+xg2;
	      
	    }
	  }

	  }  // end "cathode" loop

	}  // end draw pads
	}  // end fCath == -3 , -4
     }  

    }
    
    }  // end detElemId
  
  }

  }  // end chamber loop

}

/**************************************************************************/
void MUONModule::LoadQuadsDigits()
{
  // 
  // Draw digits
  //

  if (fDetElemId > 0 && fDetElemId != fID) return;

  Int_t fChamber;

  AliRunLoader *runLoader = Alieve::Event::AssertRunLoader();
  runLoader->LoadgAlice();

  AliMUON *pMUON = (AliMUON*)gAlice->GetModule("MUON");

  const AliMUONGeometryTransformer* kGeomTransformer = pMUON->GetGeometryTransformer();

  fChamber = fID/100;

  /*           D I S P L A Y     D I G I T S                              */

  // Display X-Y strips for the trigger chambers

  Float_t adcmax = 1024; // default
  if (fChamber<11) adcmax = 4096;

  TClonesArray *digits;
  Int_t ndigits;
  digits  = fInfo->GetDigits(fChamber);
  ndigits = digits->GetEntriesFast(); 

  AliMUONDigit  *mdig;

  Int_t fCathode = fCath;

  for (Int_t id = 0; id < ndigits; id++) {
    mdig  = (AliMUONDigit*)digits->UncheckedAt(id);
    if (mdig->Cathode() != fCathode-1) continue;

    // get all needed parameters
    Int_t ix=mdig->PadX();
    Int_t iy=mdig->PadY();
    Int_t detElemId=mdig->DetElemId();      
    Int_t charge = mdig->Signal();
    Int_t index  = Int_t(TMath::Log(charge)/(TMath::Log(adcmax)/22));
    Int_t color  = 261+index;
    Int_t colorTrigger = 2;   
    if (color > 282) color = 282;

    if (detElemId != fID) continue;
      
    const AliMpVSegmentation* seg2 
      = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,fCathode-1);

    AliMpPad pad = seg2->PadByIndices(AliMpIntPair(ix,iy),kTRUE);

    //printf("Dig ix %d iy %d \n",ix,iy);

    // time delay information
    if (fChamber > 10) { // trigger chamber
      Int_t sumCharge = 0;
      Int_t n = mdig->Ntracks();
      for (Int_t icharge = 0; icharge < n; icharge++) {
	sumCharge = sumCharge+mdig->TrackCharge(icharge);
      }
      assert(sumCharge==mdig->Signal());
      Int_t testCharge = sumCharge-(Int_t(sumCharge/n))*n;
      if(sumCharge <= n || testCharge > 0) {
	colorTrigger = 4; 
      } else {
	colorTrigger = 5; 
      }
    }
	    
    // get the pad position and dimensions
    Float_t xlocal1 = pad.Position().X();
    Float_t ylocal1 = pad.Position().Y();
    Float_t xlocal2 = pad.Dimensions().X();
    Float_t ylocal2 = pad.Dimensions().Y();

    Float_t xg1, xg2, yg1, yg2, zg1;
	    
    kGeomTransformer->Local2Global(detElemId, xlocal1, ylocal1, 0, xg1, yg1, zg1);
    // (no transformation for pad dimensions)
    xg2 = xlocal2;
    yg2 = ylocal2;

    if (fChamber > 10) {
      color = colorTrigger;
    }

    fQuads.push_back(Reve::Quad());
    fQuads.back().ColorFromIdx(color);
    Float_t* p = fQuads.back().vertices;
    // the tilting direction now known...
    /*
    p[0] = xg1-xg2;  p[1] =  yg1-yg2; p[2]  = zg1;
    p[3] = xg1-xg2;  p[4] =  yg1+yg2; p[5]  = zg1;
    p[6] = xg1+xg2;  p[7] =  yg1+yg2; p[8]  = zg1;
    p[9] = xg1+xg2;  p[10] = yg1-yg2; p[11] = zg1;
    */
    // switch x <> z
    p[0] = zg1;         p[1] =  yg1-yg2; p[2]  = xg1-xg2;
    p[3] = zg1;         p[4] =  yg1+yg2; p[5]  = xg1-xg2;
    p[6] = zg1;         p[7] =  yg1+yg2; p[8]  = xg1+xg2;
    p[9] = zg1;         p[10] = yg1-yg2; p[11] = xg1+xg2;
    /*
    if (fChamber < 11) {
    printf("coordinates.............................. \n");
    for (Int_t j = 0; j < 12; j++) {
      printf("%f ",p[j]);
      if ((j+1)%3 == 0) printf("\n");
    }
    }
    */

  }  // end loop over digits

  LoadQuadsChambers(fChamber,fChamber,fID,fCath);

}

/**************************************************************************/
void MUONModule::LoadQuadsClusters()
{
  //
  // Draw clusters
  //

  Int_t fChamber;

  fChamber = fID/100;

  LoadQuadsChambers(fChamber,fChamber,fID);

  /*           D I S P L A Y     C L U S T E R S                           */

  TClonesArray *clusters;
  Int_t nclusters;
  clusters  = fInfo->GetClusters(fChamber);
  if (clusters == 0) return;
  nclusters = clusters->GetEntriesFast(); 

  //Float_t zpos = AliMUONConstants::DefaultChamberZ(fChamber-1);  

  AliMUONRawCluster  *cls;
  for (Int_t iCls = 0; iCls < nclusters; iCls++) {

    cls = (AliMUONRawCluster*)clusters->UncheckedAt(iCls);

    if (cls->GetDetElemId() != fID) continue;

    Float_t x = cls->GetX(0);
    Float_t y = cls->GetY(0);
    Float_t z = cls->GetZ(0);

    Float_t dx = 1.0/2.0;
    Float_t dy = 1.0/2.0;
    Float_t r = TMath::Sqrt(dx*dx+dy*dy);
    /*
    // just an empty  square
    fQuads.push_back(Reve::Quad(5));
    Float_t* p = fQuads.back().vertices;
    
    //p[0] = x-dx;  p[1] =  y-dy; p[2]  = z;
    //p[3] = x-dx;  p[4] =  y+dy; p[5]  = z;
    //p[6] = x+dx;  p[7] =  y+dy; p[8]  = z;
    //p[9] = x+dx;  p[10] = y-dy; p[11] = z;
    
    // switch x <> z
    p[0] = z;     p[1] =  y-dy; p[2]  = x-dx;
    p[3] = z;     p[4] =  y+dy; p[5]  = x-dx;
    p[6] = z;     p[7] =  y+dy; p[8]  = x+dx;
    p[9] = z;     p[10] = y-dy; p[11] = x+dx;
    */
    Int_t nstep = 100;
    Float_t dstep = 2.0*TMath::Pi() / (Float_t)nstep;
    Float_t d;
    for (Int_t istep = 1; istep < (nstep+1); istep++) {

      d = istep * dstep;
      Float_t xc = x + r * TMath::Cos(d);
      Float_t yc = y + r * TMath::Sin(d);

      fQuads.push_back(Reve::Quad(5));
      Float_t* p = fQuads.back().vertices;

      p[0] = z;     p[1] =  y;  p[2]  = x;
      p[3] = z;     p[4] =  y;  p[5]  = x;
      p[6] = z;     p[7] =  yc; p[8]  = xc;
      p[9] = z;     p[10] = yc; p[11] = xc;

   }

  }

}

/**************************************************************************/
void MUONModule::LoadQuadsTracks(Int_t id)
{
  //
  // Draw tracks
  //

  /*           D I S P L A Y     T R A C K S                           */

  TClonesArray * trackParamAtHit = 0;   
  TClonesArray *tracks;
  Int_t ntracks;
  tracks  = fInfo->GetTracks();
  if (tracks == 0) return;
  ntracks = tracks->GetEntriesFast(); 

  Float_t *p;

  AliMUONTrack *recTrack = (AliMUONTrack*)tracks->At(id);

  Float_t xRec, xRec0;
  Float_t yRec, yRec0;
  Float_t zRec, zRec0;
  
  AliMUONTrackParam *trackParam = recTrack->GetTrackParamAtVertex(); 
  xRec0  = trackParam->GetNonBendingCoor();
  yRec0  = trackParam->GetBendingCoor();
  zRec0  = trackParam->GetZ();

  Int_t nTrackHits = recTrack->GetNTrackHits();
  for (Int_t iHit = 0; iHit < nTrackHits; iHit++){
    trackParamAtHit = recTrack->GetTrackParamAtHit();
    trackParam = (AliMUONTrackParam*) trackParamAtHit->At(iHit); 
    xRec  = trackParam->GetNonBendingCoor();
    yRec  = trackParam->GetBendingCoor();
    zRec  = trackParam->GetZ();

    fQuads.push_back(Reve::Quad());
    fQuads.back().ColorFromIdx(4);
    p = fQuads.back().vertices;
    
    p[0] = zRec0; p[1] =  yRec0; p[2]  = xRec0;
    p[3] = zRec0; p[4] =  yRec0; p[5]  = xRec0;
    p[6] = zRec;  p[7] =  yRec;  p[8]  = xRec;
    p[9] = zRec;  p[10] = yRec;  p[11] = xRec;

    xRec0 = xRec;
    yRec0 = yRec;
    zRec0 = zRec;

  } // end loop rec. hits

}

