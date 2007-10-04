#include "TOFSector.h"

#include <Reve/RGTopFrame.h>

#include <AliTOFdigit.h>
#include <AliTOFGeometry.h>

#include <TStyle.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

Bool_t       TOFSector::fgStaticInitDone = kFALSE;
FrameBox*    TOFSector::fgTOFsectorFrameBox = 0;
RGBAPalette* TOFSector::fgTOFsectorPalette  = 0;

//_______________________________________________________
ClassImp(TOFSector)

/* ************************************************************************ */

TOFSector::TOFSector(const Text_t* n, const Text_t* t) :
  QuadSet(n, t),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(0x0),
  fTOFtree(0x0),
  fSector(-1),
  fDx(0), fDy(0), fDz(0),
  fAutoTrans (kTRUE),
  //fPlateFlag0(kTRUE), fPlateFlag1(kTRUE), fPlateFlag2(kTRUE), fPlateFlag3(kTRUE), fPlateFlag4(kTRUE),
  fThreshold (5), fMaxVal    (80),
  fSectorID  (0),
  fGeoManager(0)
{

  fPlateFlag = new Bool_t[5];
  for (Int_t ii=0; ii<5; ii++) fPlateFlag[ii]=kTRUE;


  fGeoManager = (TGeoManager*)gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  if (!fGeoManager) {
    printf("ERROR: no TGeo\n");
  }

}
/* ************************************************************************ */

TOFSector::TOFSector(TGeoManager *localGeoManager,
		     Int_t nSector)
  :
  QuadSet(Form("Sector%i",nSector)),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(0x0),
  fTOFtree(0x0),
  fSector(nSector),
  fDx(0), fDy(0), fDz(0),
  fAutoTrans (kTRUE),
  //fPlateFlag0(kTRUE), fPlateFlag1(kTRUE), fPlateFlag2(kTRUE), fPlateFlag3(kTRUE), fPlateFlag4(kTRUE),
  fThreshold (5), fMaxVal    (80),
  fSectorID  (nSector),
  fGeoManager(localGeoManager)
{

  fPlateFlag = new Bool_t[5];
  for (Int_t ii=0; ii<5; ii++) fPlateFlag[ii]=kTRUE;

  /*
  if (!fGeoManager) {
    printf("ERROR: no TGeo\n");
  }
  */

  InitModule();

}
/* ************************************************************************ */

TOFSector::TOFSector(TGeoManager *localGeoManager,
		     Int_t nSector,
		     TClonesArray *tofArray)
  :
  QuadSet(Form("Sector%i",nSector)),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(tofArray),
  fTOFtree(0x0),
  fSector(nSector),
  fDx(0), fDy(0), fDz(0),
  fAutoTrans (kTRUE),
  //fPlateFlag0(kTRUE), fPlateFlag1(kTRUE), fPlateFlag2(kTRUE), fPlateFlag3(kTRUE), fPlateFlag4(kTRUE),
  fThreshold (5), fMaxVal    (80),
  fSectorID  (nSector),
  fGeoManager(localGeoManager)
{

  fPlateFlag = new Bool_t[5];
  for (Int_t ii=0; ii<5; ii++) fPlateFlag[ii]=kTRUE;

  InitModule();

}
/* ************************************************************************ */

TOFSector::TOFSector(TGeoManager *localGeoManager,
		     Int_t nSector,
		     TTree *tofTree)
  :
  QuadSet(Form("Sector%i",nSector)),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(0x0),
  fTOFtree(tofTree),
  fSector(nSector),
  fDx(0), fDy(0), fDz(0),
  fAutoTrans (kTRUE),
  //fPlateFlag0(kTRUE), fPlateFlag1(kTRUE), fPlateFlag2(kTRUE), fPlateFlag3(kTRUE), fPlateFlag4(kTRUE),
  fThreshold (5), fMaxVal    (80),
  fSectorID  (nSector),
  fGeoManager(localGeoManager)
{

  fPlateFlag = new Bool_t[5];
  for (Int_t ii=0; ii<5; ii++) fPlateFlag[ii]=kTRUE;

  InitModule();

}
/* ************************************************************************ */

TOFSector::~TOFSector()
{
  /*
  fGeoManager = 0x0;
  delete fGeoManager;

  delete fTOFarray;
  fTOFarray = 0x0;
  */
  delete fPlateFlag;

}

/* ************************************************************************ */
/*
void TOFSector::SetDigitsInfo(TOFDigitsInfo* info)
{
  if(fInfo) fInfo->DecRefCount();
  fInfo = info;
  if(fInfo) fInfo->IncRefCount();

}
*/
/* ************************************************************************ */
void TOFSector::InitStatics()
{
  if (fgStaticInitDone) return;

  Float_t dx = 124.5;
  Float_t dz =  29.;
  Float_t dy = 370.6*2.;
  fgTOFsectorFrameBox = new FrameBox();

  fgTOFsectorFrameBox->SetAABox(-dx*0.5, -dy*0.5, -dz*0.5, dx, dy, dz);
  fgTOFsectorFrameBox->SetFrameColor((Color_t) 32);//31);

  //fgTOFsectorPalette  = new RGBAPalette(0, 2048); // TOT
  fgTOFsectorPalette  = new RGBAPalette(0, 8192/*1024*/); // TDC
  fgTOFsectorPalette->SetLimits(0, 8192); 

  fgStaticInitDone = kTRUE;
}

/* ************************************************************************ */
void TOFSector::InitModule()
{

  fDx = fTOFgeometry->XPad()*fTOFgeometry->NpadX();
  //fDy = fTOFgeometry->XPad()*fTOFgeometry->NpadX();
  fDz = fTOFgeometry->ZPad()*fTOFgeometry->NpadZ();

  if (!fgStaticInitDone) InitStatics();

  SetFrame(fgTOFsectorFrameBox);
  SetPalette(fgTOFsectorPalette);
  //fFrame   = fgTOFsectorFrameBox;
  //fPalette = fgTOFsectorPalette;

  LoadQuads();  
  ComputeBBox();
  SetTrans();

}

/* ************************************************************************ */
void TOFSector::LoadQuads()
{

  Reset(QT_FreeQuad, kFALSE, 32);

  //Int_t n_col = gStyle->GetNumberOfColors();

  Int_t vol[5] = {fSectorID, -1, -1, -1, -1};
  Int_t informations[4] = {-1, -1, -1, -1};
  Int_t dummy[3] = {-1, -1, -1};
  Int_t tdc = -1;
  Int_t tot = -1;

  Double_t **coord = new Double_t*[4];
  for (Int_t ii=0; ii<4; ii++) coord[ii] = new Double_t[3];

  //printf(" fTOFarray->GetEntries() = %4i \n",fTOFarray->GetEntries());

  if (fTOFtree && !fTOFarray) {
    //printf("Hello world\n");
    TClonesArray* digitsTOFnew = new TClonesArray("AliTOFdigit",  300);

    fTOFarray = new TClonesArray("AliTOFdigit",  300);
    TClonesArray &ldigits = *fTOFarray;
    Int_t newCounter = 0;

    AliTOFdigit *digs;

    fTOFtree->SetBranchAddress("TOF",&digitsTOFnew);
    fTOFtree->GetEntry(0);
    for (Int_t digitNumber=0; digitNumber<digitsTOFnew->GetEntries(); digitNumber++) {

      //if (digitNumber==digitsTOF->GetEntries()-1) printf(" Hello  4 -> %3i digit of %i \n", digitNumber+1, digitsTOF->GetEntries());
  
      digs = (AliTOFdigit*)digitsTOFnew->UncheckedAt(digitNumber);

      if (digs->GetSector()!=fSectorID) continue;

      vol[1] = digs->GetPlate();  // Plate Number (0-4)
      vol[2] = digs->GetStrip();  // Strip Number (0-14/18)
      vol[3] = digs->GetPadx();   // Pad Number in x direction (0-47)
      vol[4] = digs->GetPadz();   // Pad Number in z direction (0-1)

      informations[0] = digs->GetTdc();
      informations[1] = digs->GetAdc();
      informations[2] = digs->GetToT();
      informations[3] = digs->GetTdcND();
      new (ldigits[newCounter++]) AliTOFdigit(dummy, vol, informations);
    }
  }

  AliTOFdigit *tofDigit;
  //printf("  0x%lx\n",fTOFarray);
  for (Int_t ii=0; ii<fTOFarray->GetEntries(); ii++) {

    tofDigit = (AliTOFdigit*)fTOFarray->UncheckedAt(ii);
   
    if (fPlateFlag[tofDigit->GetPlate()]) {

      vol[1] = tofDigit->GetPlate();
      vol[2] = tofDigit->GetStrip();
      vol[3] = tofDigit->GetPadz();
      vol[4] = tofDigit->GetPadx();

      tot = tofDigit->GetToT();
      tdc = tofDigit->GetTdc();

      //if (fSector==4 && fPlate==2 && fStrip==0) printf(" %2i  %1i\n", iPadX, iPadZ);
      //if (iPadX==23 || iPadX==24) printf(" %2i  %1i %2i \n", fSector, fPlate, fStrip);
      //if (vol[0]==4 && vol[1]==2 && vol[2]==7) {

      for (Int_t kk=0; kk<4; kk++) for (Int_t jj=0; jj<3; jj++) coord[kk][jj]=0.;

      fTOFgeometry->DetToSectorRF(vol, coord);
      /*
	printf("\n");
	printf("  %1i   %2i,   %f  %f  %f \n", vol[3], vol[4], coord[0][0], coord[0][1], coord[0][2]);
	printf("  %1i   %2i,   %f  %f  %f \n", vol[3], vol[4], coord[1][0], coord[1][1], coord[1][2]);
	printf("  %1i   %2i,   %f  %f  %f \n", vol[3], vol[4], coord[2][0], coord[2][1], coord[2][2]);
	printf("  %1i   %2i,   %f  %f  %f \n", vol[3], vol[4], coord[3][0], coord[3][1], coord[3][2]);
      */
      Float_t vertices[12]={(Float_t)coord[0][0], (Float_t)coord[0][1], (Float_t)coord[0][2],
			    (Float_t)coord[1][0], (Float_t)coord[1][1], (Float_t)coord[1][2],
			    (Float_t)coord[2][0], (Float_t)coord[2][1], (Float_t)coord[2][2],
			    (Float_t)coord[3][0], (Float_t)coord[3][1], (Float_t)coord[3][2]};
      
      AddQuad(vertices);
      //AddQuad((Float_t*)coord);
      //AddQuad(coord[0], coord[1], coord[2], 2.5, 3.5);
      //AddQuad(-2.5*0.5, -3.5*0.5, 2.5, 3.5);

      // In principle could have color based on number of neigbours. We
      // can insert the time-of-flight value for each pad
      //QuadValue((Int_t)tot);
      QuadValue((Int_t)tdc);
      QuadId(tofDigit);
      
    //}
    } // closed if control on plates switched on
  } // closed loop on TOF sector digits

  RefitPlex();

  fTOFarray = 0x0;
 
}

/* ************************************************************ */
void TOFSector::SetTrans()
{
  fHMTrans.UnitTrans();

  //Int_t det[5] = {fSector, -1, -1, -1, -1};
  Char_t path[100];

  Int_t localSector = fSector;
  if (!fAutoTrans) localSector = 4;

  //fTOFgeometry->GetVolumePath(det,path);
  fTOFgeometry->GetVolumePath(localSector,path);
  fGeoManager->cd(path);
  TGeoHMatrix global = *fGeoManager->GetCurrentMatrix();
  Double_t *rotMat = global.GetRotationMatrix();
  Double_t *tr = global.GetTranslation();

  fHMTrans.SetBaseVec(1, rotMat[0], rotMat[3], rotMat[6]);
  fHMTrans.SetBaseVec(2, rotMat[1], rotMat[4], rotMat[7]);
  fHMTrans.SetBaseVec(3, rotMat[2], rotMat[5], rotMat[8]);

  fHMTrans.SetBaseVec(4, tr);
}

//-----------------------------------------------------

void TOFSector::SetSectorID(Int_t id)
{
  fSectorID = id;
  fSector   = id;
  if (fAutoTrans)
    SetTrans(); // Force repositioning.

  LoadQuads();
}

//-----------------------------------------------------

void TOFSector::SetPlate(Int_t nPlate, Bool_t r)
{

  fPlateFlag[nPlate] = r;

  //printf("   HELLO World ! %i %i %i \n", nPlate, r, fPlateFlag[nPlate]);
}

/**************************************************************************/

void TOFSector::SetThreshold(Short_t t)
{
  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  // ClearColorArray();
}

/**************************************************************************/

void TOFSector::SetMaxVal(Int_t mv)
{
  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  //ClearColorArray();
}

/**************************************************************************/

void TOFSector::QuadSelected(Int_t idx)
{
  // Override control-click from QuadSet

  QuadBase* qb   = GetQuad(idx);
  TObject* obj   = qb->fId.GetObject();
  AliTOFdigit* digs = dynamic_cast<AliTOFdigit*>(obj);
  // printf("TOFSector::QuadSelected "); Print();
  printf("  idx = %5i, value = %5d, obj = 0x%lx, digit = 0x%lx  ",
	 idx, qb->fValue, (ULong_t)obj, (ULong_t)digs);
  if (digs)
    printf("-> Sector = %2i  Plate = %1i  Strip = %2i  ToT = %3i  Tof = %5i\n",
	   fSector , digs->GetPlate(), digs->GetStrip(), digs->GetToT(), digs->GetTdc());
  else printf("\n");

}

/**************************************************************************/
