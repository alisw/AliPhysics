///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Ring Imaging Cherenkov                                                   //
//  This class contains the basic functions for the Ring Imaging Cherenkov   //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliRICHClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h> 
#include <TNode.h> 
#include <TRandom.h> 

#include <TClass.h> 
#include "AliRICH.h"
#include "AliRun.h"
#include "TGeant3.h" 


ClassImp(AliRICH)

//_____________________________________________________________________________
AliRICH::AliRICH()
{
  //
  // Default constructor for RICH
  //
  fIshunt   = 0;
  fHits     = 0;
  fMips     = 0;
  fCkovs    = 0;
  fPadhits  = 0;
  fNmips    = 0;
  fNckovs   = 0;
  fNpadhits = 0;
  
  fChslope  = 0;
  fAlphaFeed= 0;
  fSxcharge = 0;
  fIritri   = 0;
}
 
//_____________________________________________________________________________
AliRICH::AliRICH(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for RICH
  //
  fIshunt     = 0;
  fNmips      = 0;
  fNckovs     = 0;
  fNpadhits   = 0;
  //
  // Allocate space for different components
  fHits     = new TClonesArray("AliRICHhit", 100);
  fMips     = new TClonesArray("AliRICHmip", 100);
  fCkovs    = new TClonesArray("AliRICHckov", 100);
  fPadhits  = new TClonesArray("AliRICHpadhit", 100);
  //
  // Set parameters to default value
  fChslope  = 40;
  fAlphaFeed= 0.04;
  fSxcharge = 0.18;
  fIritri   = 0;
  //
  SetMarkerColor(6);
  SetMarkerStyle(20);
  SetMarkerSize(0.5);
}

//_____________________________________________________________________________
AliRICH::~AliRICH()
{
  //
  // Destructor for RICH
  //
  fIshunt   = 0;
  delete fHits;
  fMips->Delete();    delete fMips;
  fCkovs->Delete();   delete fCkovs;
  fPadhits->Delete(); delete fPadhits;
}

//_____________________________________________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a RICH hit
  //
  switch ((int) hits[0]) {
  case 0:
    {
      //
      // Simple hit
      TClonesArray &lhits = *fHits;
      new(lhits[fNhits++]) AliHit(fIshunt,track);
      AliHit *lhit = (AliHit*) fHits->AddrAt(fNhits-1);
      lhit->fX = hits[19];
      lhit->fY = hits[20];
      lhit->fZ = hits[21];
    }	break;
  case 1:
    //
    // MIP hit
    AddMipHit(track,vol,hits);
    break;
  case 2:
    //
    // Cherenkov hit
    AddCkovHit(track,vol,hits);
    break;
  case 3:
    //
    // Pad hit
    AddPadHit(track,vol,hits);
    break;
  case 4:
    // Update a mip hit
    UpdateMipHit(hits);
    break; 
  default:
    printf("Error: AliRICH::AddHit flag %d not defined./n",(int) hits[0]);
    return;
  }    
}

//_____________________________________________________________________________
void AliRICH::AddMipHit(Int_t track, Int_t *vol, Float_t *hits)
{ 
    // Adds a mip hit in the RICH.
    TClonesArray &lhits = *fMips;
    new(lhits[fNmips++]) AliRICHmip(fIshunt,track,vol,hits,
				    fNckovs,fNpadhits);
}

//_____________________________________________________________________________
void AliRICH::AddCkovHit(Int_t track, Int_t *vol, Float_t *hits)
{ 
  //
  // Adds a cerenkov hit in the RICH.
  //
  TClonesArray &lhits = *fCkovs;
  AliRICHmip *lmip = (AliRICHmip*) fMips->AddrAt(fNmips-1);
  //
  // If this ckov come from a mip update the mip lastckov entry.
  Int_t fmipslocal=-1;
  if (lmip->GetZ() != -999.0) {
    fmipslocal    =   fNmips-1;
    lmip->SetLastCkov(fNckovs);
  }
  new(lhits[fNckovs++]) AliRICHckov(fIshunt,track,vol,hits,
				    fmipslocal,fNpadhits);
}

//_____________________________________________________________________________
void AliRICH::AddPadHit(Int_t track, Int_t *vol, Float_t *hits)
{ 
  //
  // Adds pad hits in the RICH
  //
  TClonesArray &lhits = *fPadhits;
  // Update the last padhit of the respective particle:
  if ((int) hits[1]==50) { // a ckov
    ((AliRICHckov *) fCkovs->AddrAt(fNckovs-1))->SetLastpad(fNpadhits);
    new(lhits[fNpadhits++]) AliRICHpadhit(fIshunt,track,vol,hits,-1,fNckovs-1);
  }else { // a mip
    ((AliRICHmip *) fMips->AddrAt(fNmips-1))->SetLastpad(fNpadhits);
    new(lhits[fNpadhits++]) AliRICHpadhit(fIshunt,track,vol,hits,fNmips-1,-1);
  }
}

//_____________________________________________________________________________
void AliRICH::BuildGeometry()
{
  //
  // Builds a TNode geometry for event display
  //
  TNode *Node, *Top;
  
  const int kColorRICH = kGreen;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90);
  new TRotMatrix("rot994","rot994",90,-20,90,70,0,0);
  new TRotMatrix("rot995","rot995",90,0,90,90,0,0);
  new TRotMatrix("rot996","rot996",90,20,90,110,0,0);
  new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70);
  new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90);
  new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110);
  new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);
  Top->cd();
  Node = new TNode("RICH1","RICH1","S_RICH",0,471.8999,165.2599,"rot993");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH2","RICH2","S_RICH",171,470,0,"rot994");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH3","RICH3","S_RICH",0,500,0,"rot995");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH4","RICH4","S_RICH",-171,470,0,"rot996");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH5","RICH5","S_RICH",161.3999,443.3999,-165.3,"rot997");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH6","RICH6","S_RICH",0,471.8999,-165.3,"rot998");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("RICH7","RICH7","S_RICH",-161.399,443.3999,-165.3,"rot999");
  Node->SetLineColor(kColorRICH);
  fNodes->Add(Node); 
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//                       Section for Rich Step
#define MAXPH 1000

static Float_t sVloc[3];
static Float_t sVectIn[3];
static Float_t sDetot;
static Float_t sYanode[MAXPH];
static Float_t sXpad[MAXPH];
static Float_t sYpad[MAXPH];
static Int_t sNpx;
static Int_t sNpy;
static Float_t sDxp;
static Float_t sDyp;
static Float_t sDlx;
static Float_t sDly;
static Float_t sDpad;


static Int_t sNsecon;
//static Float_t sQint[2];

#define MAXSEC 2000

static Int_t sNpads;
static Int_t sIpx[MAXSEC];
static Int_t sIpy[MAXSEC];
static Int_t sIqpad[MAXSEC];
static Int_t sNphlink[MAXSEC];
static Int_t sNphoton;

static Int_t sNfeed, sNfeedd, sNdir;

static Float_t sPparent;
static Float_t sThincParent;
static Int_t sIloss[MAXPH];
static Int_t sIprod[MAXPH];
static Float_t sXphit[MAXPH];
static Float_t sYphit[MAXPH];

static Float_t sSycharge;
static Float_t sXsox, sYsox, sZsox;

static Int_t sMckov[MAXPH];
static Int_t idpartgx;
static Float_t phisx;
static Int_t nmodsx;
static Float_t psx;
static Float_t xsx;
static Float_t ysx;
static Float_t thetasx;


//_____________________________________________________________________________
void AliRICH::Init()
{
  //
  // Initialise RICH after that it has been built
  //
  const Float_t sconv=2*TMath::Sqrt(2*TMath::Log(2));
  const Float_t yK3=1.20;
  Float_t ansp;
  Int_t i;
  //
  sNpx=162;
  sNpy=162;
  sDxp=0.80;
  sDyp=0.80;
  ansp=sDyp/2;
  sDlx=sNpx*sDxp/2;
  sDly=sNpy*sDyp/2;
  sDpad=0.2;
  //
  for(i=0;i<sNpx;i++) {
    sXpad[i]=i*sDxp;
    sYpad[i]=i*sDyp;
  }
  for(i=0;i<2*sNpy+1;i++) sYanode[i]=ansp/2+i*ansp;
  //
  
  sSycharge=4*TMath::ATanH(1/TMath::Sqrt(2*yK3))/TMath::Pi()
    /(1-0.5*TMath::Sqrt(yK3))/sDpad/sconv;
  //
}


//___________________________________________________________________________
void AliRICH::StepManager()
{
  //
  // Called at every step in the RICH
  //

  AliMC* pMC = AliMC::GetMC();
  TGeant3 *geant3 = (TGeant3*) pMC;

  const Float_t xshift[3] = { 41.3, 0, -41.3 };
  static Float_t polar[3] = {0, 0, 0};
  const Int_t nrooth = 25;
  
  static Int_t ixold=-1, iyold=-1;
  
  // System generated locals 
  Int_t j, i1;
  Float_t r1, r2;
  
  // Local variables 
  Float_t ranf[2], rrhh[nrooth], phiangle, cost, vmod;
  //Int_t idpartsx;
  Int_t i;
  Float_t t, vxloc[3];
  Int_t ll, irivol[2];
  Int_t lcase;
  //Int_t iprimx;
  Int_t ix, iy;
  Float_t stwght;
  Int_t ncher;
  Float_t cophi;
  Float_t dir[3];
  Int_t ihitrak;
  Int_t medprod;
  
  Int_t nmult=0;
  //Float_t xtrig[200], ytrig[200];
  //Int_t itrig[200];

  
  
  //     ILOSS = 0    NOT LOST 
  //             1    REFLECTED FREON-QUARZ 
  //             2    REFLECTED QUARZ-METHANE 
  //             3    REFLECTED METHANE-CSI 
  //            11    ABSORBED IN FREON 
  //            12    ABSORBED IN QUARZ 
  //            13    ABSORBED IN METHANE 
  //            21    CSI QE 
  //     IPROD = 1    PRODUCED IN FREON 
  //     IPROD = 2    PRODUCED IN QUARZ 
  
  // new (changed NROOTH from 10 to 25!!!!!!!!!!!!!) 
  
  
  Int_t *idtmed = gAlice->Idtmed();
  
  //--------------------------------------------------------------------------
  
  //        MIP inside CsI 

  if (geant3->Gckine()->charge) {
    
    //        Charged particles treatment 
    if (fIritri && !geant3->Gctrak()->upwght) {
      if (geant3->Gctmed()->numed == idtmed[fIritri-1]) {
	if (geant3->Gcking()->ngkine > 0) {
	  strncpy((char *)&lcase,"HADR",4);
	  if (geant3->Gcking()->kcase == lcase) {
	    i1 = geant3->Gcking()->ngkine;
	    for (i = 1; i <= i1; ++i) {
	      pMC->Gmtod(geant3->Gckin3()->gpos[i-1], vxloc, 1);
	      pMC->Gmtod(geant3->Gcking()->gkin[i-1], dir, 2);
	      if (geant3->Gcking()->gkin[i-1][4] == 8. || 
		  geant3->Gcking()->gkin[i-1][4] == 9.) {
		++nmult;
		// Computing 2nd power 
		r1 = dir[0];
		// Computing 2nd power 
		r2 = dir[2];
		//theta = TMath::ATan2(TMath::Sqrt(r1*r1+r2*r2),dir[1]);
		//xtrig[nmult - 1] = theta;
		//ytrig[nmult - 1] = vxloc[1] + .25;
		//itrig[nmult - 1] = (Int_t) geant3->Gcking()->gkin[i-1][4];
	      }
	    }
	  }
	}
      }
      if ((geant3->Gctmed()->numed == idtmed[1006-1] &&
	   geant3->Gctrak()->inwvol == 2) || 
	  geant3->Gctrak()->istop) {
	if (!nmult) {
	  printf("NOT TRIGGERED\n");
	  sDetot = 0.;
	  sNsecon = 0;
	  sNpads = 0;
	  sNphoton = 0;
	  sNfeed = 0;
	  sNfeedd = 0;
	  sNdir = 0;
	  geant3->Gctrak()->istory = 0;
	  geant3->Gctrak()->upwght = 0.;
	  geant3->Gcflag()->ieotri = 1;
	  nmult = 0;
	  //sQint[0] = 0.;
	  //sQint[1] = 0.;
	} else {
	  printf("TRIGGERED %d\n",nmult);
	}
      }
    }
    //        MIP inside Methane 
    if (geant3->Gctmed()->numed == idtmed[1009-1]) {
      
      // new 
      //     If particle produced already Cerenkov Photons (istory=1) 
      //     update the impact point only 
      if (geant3->Gctrak()->istory == 1) {
	//     Direction of incidence and where did it hit ? 
	pMC->Gmtod(geant3->Gctrak()->vect, sVloc, 1);
	pMC->Gmtod(&geant3->Gctrak()->vect[3], dir, 2);
	phiangle = TMath::ATan2(dir[2], dir[0]);
	if (phiangle < 0.) phiangle += 2*TMath::Pi();
	i1 = nrooth;
	for (ll = 0; ll < i1; ++ll) rrhh[ll] = 0;
	irivol[0] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	irivol[1] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	// NMODSX 
	ihitrak = gAlice->CurrentTrack();
	rrhh[0] = 4.;
	// flag to say this is update 
	rrhh[1] = sVloc[0] + sDlx;
	// XSX 
	rrhh[2] = sVloc[2] + sDly;
	// YSX 
	rrhh[3] = (Float_t) geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	// NMODSX 
	// Computing 2nd power 
	r1 = dir[0];
	// Computing 2nd power 
	r2 = dir[2];
	rrhh[4] = TMath::ATan2(TMath::Sqrt(r1 * r1 + r2 * r2), dir[1]);
	// theta 
	rrhh[5] = phiangle;
	//          PRINT *, 'ITRA = ',ITRA,'ISTAK = ',ISTAK 
	AddHit(ihitrak,irivol,rrhh);
      }
      // enew 
      //        Record particle properties 
      //        If particle produced already Cerenkov Photons (istory=1)
      
      //        update the impact point only 
      if (geant3->Gctrak()->istory != 2) {
	if (!geant3->Gctrak()->istory) {
	  ++sNsecon;
	  
	  //        Is this a primary particle ? 
	  //iprimx = 1;
	  //if (geant3->Gctrak()->upwght) iprimx = 0;
	  
	  //        Where did it come from ? 
	  sXsox = geant3->Gckine()->vert[0];
	  sYsox = geant3->Gckine()->vert[1];
	  sZsox = geant3->Gckine()->vert[2];
	  
	  //        Momentum 
	  psx = geant3->Gctrak()->vect[6];
	  
	  //        Particle type and parent 
	  //idpartsx = geant3->Gckine()->ipart;
	  r1 = geant3->Gctrak()->upwght / 100.;
	  idpartgx = Int_t(r1+0.5);
	  if (!geant3->Gctrak()->upwght) {
	    sPparent = geant3->Gctrak()->vect[6];
	    sThincParent = thetasx;
	  }
	  
	  //        Direction of incidence and where did it hit ? 
	  pMC->Gmtod(geant3->Gctrak()->vect, sVloc, 1);
	  pMC->Gmtod(&geant3->Gctrak()->vect[3], dir, 2);
	  // Computing 2nd power 
	  r1 = dir[0];
	  // Computing 2nd power 
	  r2 = dir[2];
	  thetasx = TMath::ATan2(TMath::Sqrt(r1 * r1 + r2 * r2), dir[1]);
	  phisx = TMath::ATan2(dir[2], dir[0]);
	  if (phisx < 0.) phisx += 2*TMath::Pi();
	  ysx = sVloc[2] + sDly;
	  xsx = sVloc[0] + sDlx;
	  nmodsx = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	  //     new 
	  i1 = nrooth;
	  for (ll = 0; ll < i1; ++ll) rrhh[ll] = 0;
	  irivol[0] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	  irivol[1] = nmodsx;
	  ihitrak = gAlice->CurrentTrack();
	  rrhh[0] = 1.;
	  // Flag to say this is MIP 
	  rrhh[1] = (Float_t) geant3->Gckine()->ipart;
	  rrhh[2] = xsx;
	  rrhh[3] = ysx;
	  rrhh[4] = (Float_t) nmodsx;
	  // Module Number 
	  rrhh[5] = thetasx;
	  rrhh[6] = geant3->Gctrak()->tofg;
	  // in seconds 
	  rrhh[7] = (Float_t) idpartgx;
	  //     mips specific 
	  rrhh[8] = phisx;
	  rrhh[9] = psx;
	  //     charge of current particle in electron charge unit; 
	  rrhh[10] = geant3->Gckine()->charge;
	  rrhh[11] = -999.;
	  // Zo of ckov generation 
	  rrhh[12] = 0.;
	  // no ckov !!! 
	  AddHit(ihitrak, irivol, rrhh);
	  //     end of new 
	  
	  //        Earmark track as being recorded in methane gap 
	  
	  geant3->Gctrak()->istory = 2;
	}
      }
      
      //        Signal generation in methane gap 
      pMC->Gmtod(geant3->Gctrak()->vect, sVloc, 1);
      pMC->Gmtod(geant3->Gckine()->vert, vxloc, 1);
      ix = (Int_t) ((sVloc[0] + sDlx) /  sDxp);
      iy = (Int_t) ((sVloc[2] + sDly) /  sDyp);
      
      //        Is this the first step? 
      if (ixold == -1 && iyold == -1) {
	ixold = ix;
	iyold = iy;
	for(j=0;j<3;j++) sVectIn[j]=geant3->Gctrak()->vect[j];
      }
      
      //        Mip left gap 
      if (geant3->Gctrak()->inwvol == 2 || geant3->Gctrak()->istop) {
	sDetot += geant3->Gctrak()->destep;
	if (sDetot > 0.) RichIntegration();
	sDetot = 0.;
	ixold = -1;
	iyold = -1;
	
	//        Mip left current pad 
      } else if (ixold != ix || iyold != iy) {
	if (sDetot > 0.) RichIntegration();
	for(j=0;j<3;j++) sVectIn[j]=geant3->Gctrak()->vect[j];
	sDetot = geant3->Gctrak()->destep;
	ixold = ix;
	iyold = iy;
      } else {
	sDetot += geant3->Gctrak()->destep;
      }
    }
  }
  
  //        End charged particles treatment 
  //
  //        Treat photons produced in Freon and Quartz 
  if (geant3->Gckin2()->ngphot > 0 && 
      (geant3->Gctmed()->numed == idtmed[1004-1] || 
       geant3->Gctmed()->numed == idtmed[1003-1])) {
    if (!geant3->Gctrak()->upwght) {
      
      //        If it is a primary, save all generated photons 
      i1 = geant3->Gckin2()->ngphot;
      for (i = 1; i <= i1; ++i) {
	++sNphoton;
	if (sNphoton > MAXPH) {
	  sNphoton = MAXPH;
	  printf("ATTENTION NPHOTON %d\n",sNphoton);
	  continue;
	}
	
	//        Production medium 
	medprod = 1;
	if (geant3->Gctmed()->numed == idtmed[1003-1]) medprod = 2;
	//
	//        Production angle and energy 
	vmod=0;
	cost=0;
	for(j=0;j<3;j++) {
	  cost+=geant3->Gckin2()->xphot[i-1][3+j]*geant3->Gctrak()->vect[3+j];
	  vmod+=geant3->Gckin2()->xphot[i-1][3+j]*
	    geant3->Gckin2()->xphot[i-1][3+j];
	}
	cost/=sqrt(vmod);
	sIloss[sNphoton - 1] = 22;
	sIprod[sNphoton - 1] = medprod;
	sXphit[sNphoton - 1] = 0.;
	sYphit[sNphoton - 1] = 0.;
	stwght = geant3->Gctrak()->upwght;
	geant3->Gctrak()->upwght = (Float_t) sNphoton;
	geant3->Gskpho(i);
	gAlice->SetTrack(0, gAlice->CurrentTrack(), 50, 
			 &geant3->Gckin2()->xphot[i-1][3],geant3->Gckin2()->xphot[i-1],
			 polar,geant3->Gctrak()->tofg,"Cherenkov", ncher);
	sMckov[sNphoton - 1] = ncher;
	geant3->Gctrak()->upwght = stwght;
      }
    } else {
      stwght = geant3->Gctrak()->upwght;
      geant3->Gctrak()->upwght = 0.;
      geant3->Gskpho(0);
      geant3->Gctrak()->upwght = stwght;
    }
    
    //        Particle did not yet pass the methane gap 
    if (geant3->Gctrak()->istory == 0) {
      geant3->Gctrak()->istory = 1;
      ++sNsecon;
      //        Is this a primary particle ? 
      //iprimx = 1;
      //if (geant3->Gctrak()->upwght) iprimx = 0;
      
      //        Where did it come from ? 
      sXsox = geant3->Gckine()->vert[0];
      sYsox = geant3->Gckine()->vert[1];
      sZsox = geant3->Gckine()->vert[2];
      
      //        Where did it hit ? 
      pMC->Gmtod(geant3->Gctrak()->vect, sVloc, 1);
      pMC->Gmtod(&geant3->Gctrak()->vect[3], dir, 2);
      ysx = sVloc[2] + sDly;
      if (geant3->Gctmed()->numed == idtmed[1004-1]) {
	nmodsx = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
	xsx = sVloc[0] + sDlx + 
	  xshift[geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 2] - 1];
      } else if (geant3->Gctmed()->numed == idtmed[1003-1]) {
	nmodsx = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 3];
	xsx = sVloc[0] + sDlx;
      } else {
	nmodsx = 0;
      }
      
      //        Momentum and direction of incidence 
      psx = geant3->Gctrak()->vect[6];
      // Computing 2nd power 
      r1 = dir[0];
      // Computing 2nd power 
      r2 = dir[2];
      thetasx = TMath::ATan2(TMath::Sqrt(r1 * r1 + r2 * r2), dir[1]);
      phisx = TMath::ATan2(dir[2], dir[0]);
      if (phisx < 0.) phisx += 2*TMath::Pi();
      
      //        Particle type and parent 
      //idpartsx = geant3->Gckine()->ipart;
      r1 = geant3->Gctrak()->upwght / 100.;
      idpartgx = Int_t(r1+0.5);
      if (!geant3->Gctrak()->upwght) {
	sPparent = geant3->Gctrak()->vect[6];
	sThincParent = thetasx;
      }
      // new 
      for (ll = 0; ll <nrooth; ++ll) rrhh[ll] = 0;
      irivol[0] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
      irivol[1] = nmodsx;
      ihitrak = gAlice->CurrentTrack();
      rrhh[0] = 1.;
      // Flag to say that this is MIP 
      rrhh[1] = (Float_t) geant3->Gckine()->ipart;
      rrhh[2] = xsx;
      rrhh[3] = ysx;
      rrhh[4] = (Float_t) nmodsx;
      // Module Number 
      rrhh[5] = thetasx;
      rrhh[6] = geant3->Gctrak()->tofg;
      // in seconds 
      rrhh[7] = (Float_t) idpartgx;
      //     mips specific 
      rrhh[8] = phisx;
      rrhh[9] = psx;
      //     charge of current particle in electron charge unit; 
      rrhh[10] = geant3->Gckine()->charge;
      rrhh[11] = sVloc[1];
      // Zo of ckov generation 
      rrhh[12] = 1.;
      // ckov generation 
      AddHit(ihitrak, irivol, rrhh);
      // enew 
    }
  }
  
  //        Current particle is cherenkov photon 
  if (geant3->Gckine()->ipart == 50) {
    pMC->Gmtod(geant3->Gctrak()->vect, sVloc, 1);
    //         WRITE(6,* ) UPWGHT, VLOC(2), NUMED, DESTEP 
    //        Photon crosses ch4-csi boundary 
    //           take into account fresnel losses with complex refraction index
    if (geant3->Gctrak()->inwvol == 1 && geant3->Gctmed()->numed == idtmed[1006-1]) {
      
      //     fresnel losses commented out for the moment 
      //    make sure that qe correction for fresnel losses is also switched off !
      //            CALL FRESNELCSI 
      //            IF (ISTOP .EQ. 2) RETURN 
      //        Put transmission of electrodes in by hand 
      pMC->Gmtod(&geant3->Gctrak()->vect[3], dir, 2);
      cophi = TMath::Cos(TMath::ATan2(dir[0], dir[1]));
      t = (1. - .025 / cophi) * (1. - .05 /  cophi);
      pMC->Rndm(ranf, 1);
      if (ranf[0] > t) {
	if (geant3->Gctrak()->upwght && Int_t(geant3->Gctrak()->upwght+0.5)<MAXPH) 
	  sIloss[Int_t(geant3->Gctrak()->upwght+0.5) - 1] = 15;
	geant3->Gctrak()->istop = 2;
	return;
      }
    }
    
    //        Photon lost energy in CsI 
    if (geant3->Gctrak()->destep > 0. && geant3->Gctmed()->numed == idtmed[1006-1]) {
      geant3->Gctrak()->istop = 2;
      r1 = geant3->Gctrak()->upwght / 100.;
      if (Int_t(r1+0.5) > 50) {
	++sNfeedd;
      } else {
	++sNdir;
      }
      //            WRITE(6,*) 'PHOTON',UPWGHT, MAXPH 
      if (geant3->Gctrak()->upwght && Int_t(geant3->Gctrak()->upwght+0.5) < MAXPH) 
	sIloss[Int_t(geant3->Gctrak()->upwght+0.5) - 1] = 0;
      //            write(6,*) 'photon detected' 
      for(j=0;j<3;j++) sVectIn[j]=geant3->Gctrak()->vect[j];
      // new 
      // copied from miphit in Freon or Quartz 
      //        Where did it hit ? 
      pMC->Gmtod(&geant3->Gctrak()->vect[3], dir, 2);
      
      //        Momentum and direction of incidence 
      for (ll = 0; ll < nrooth; ++ll) rrhh[ll]=0;
      irivol[0] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
      // ??? 
      irivol[1] = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
      // ??? 
      ihitrak = gAlice->CurrentTrack();
      rrhh[0] = 2.;
      // Flag to say that this is CK 
      rrhh[1] = (Float_t) geant3->Gckine()->ipart;
      rrhh[2] = sVloc[0] + sDlx;
      rrhh[3] = sVloc[2] + sDly;
      rrhh[4] = (Float_t) geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 3];
      // ??? Module Number 
      // Computing 2nd power 
      r1 = dir[0];
      // Computing 2nd power 
      r2 = dir[2];
      rrhh[5] = TMath::ATan2(TMath::Sqrt(r1 * r1 + r2 * r2), dir[1]);
      // THETASX 
      rrhh[6] = geant3->Gctrak()->tofg;
      // in seconds 
      r1 = geant3->Gctrak()->upwght / 100.;
      rrhh[7] = (Float_t) Int_t(r1+0.5);
      //     ckov specific 
      // Feedback ??? 
      rrhh[8] = geant3->Gctrak()->getot;
      rrhh[9] = 0.;
      // Stop in CsI 
      AddHit(ihitrak, irivol, rrhh);
      //     end of new 
      RichIntegration();
      return;
    }
    if (geant3->Gctrak()->upwght && Int_t(geant3->Gctrak()->upwght+0.5) < MAXPH) 
      {
	//        Losses by reflection and absorption and leaving detector
	if (sIloss[Int_t(geant3->Gctrak()->upwght+0.5) - 1] == 22) {
	  i1 = geant3->Gctrak()->nmec;
	  for (i = 0; i < i1; ++i) {
	    //        Reflection loss 
	    if (geant3->Gctrak()->lmec[i] == 106) {
	      sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=10;
	      if (geant3->Gctmed()->numed == idtmed[1004-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=1;
	      
	      if (geant3->Gctmed()->numed == idtmed[1003-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=2;
	      
	      //        Absorption loss 
	    } else if (geant3->Gctrak()->lmec[i] == 101) {
	      geant3->Gctrak()->istop = 2;
	      sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=20;
	      if (geant3->Gctmed()->numed == idtmed[1004-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=11;
	      
	      if (geant3->Gctmed()->numed == idtmed[1003-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=12;
	      
	      if (geant3->Gctmed()->numed == idtmed[1005-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=13;
	      
	      if (geant3->Gctmed()->numed == idtmed[1009-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=13;
	      
	      if (geant3->Gctmed()->numed == idtmed[1001-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=14;
	      
	      //        CsI inefficiency 
	      if (geant3->Gctmed()->numed == idtmed[1006-1]) 
		sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=16;
	      
	      
	      //        Photon goes out of tracking scope 
	    } else if (geant3->Gctrak()->lmec[i] == 30) 
	      sIloss[Int_t(geant3->Gctrak()->upwght+0.5)-1]=21;
	    
	  }
	}
      }
  }
}


//_____________________________________________________________________________
void AliRICH::RichIntegration()
{
  //
  // Integrates over RICH pads
  //

  AliMC* pMC = AliMC::GetMC();
  TGeant3 *geant3 = (TGeant3*) pMC;
  
  Int_t i1, i2;
  Float_t r1;
  
  // Local variables 
  Float_t rrhh[25];
  Float_t qtot;
  Int_t ifeed;
  Float_t x0;
  Int_t ixmin, ixmax, iymin, iymax;
  Int_t ll, ix, iy;
  Float_t qp;
  Float_t source[3];
  Int_t irivol[2];
  Float_t y0a, qtot_check, xi1, xi2, yi1, yi2;
  Int_t nph = 0, iqp;
  Int_t ihitrak, modulen;
  //Int_t isec[MAXSEC];
  //     ILOSS = 0    NOT LOST 
  //             1    REFLECTED FREON-QUARZ 
  //             2    REFLECTED QUARZ-METHANE 
  //             3    REFLECTED METHANE-CSI 
  //            11    ABSORBED IN FREON 
  //            12    ABSORBED IN QUARZ 
  //            13    ABSORBED IN METHANE 
  //            21    CSI QE 
  //     IPROD = 1    PRODUCED IN FREON 
  //     IPROD = 2    PRODUCED IN QUARZ 
  //       get charge for MIP or cherenkov photon 
  
  if (geant3->Gckine()->ipart == 50) {
    GetCharge(qtot);
    modulen = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 3];
    //sQint[1] = qtot;
  } else {
    GetChargeMip(qtot);
    modulen = geant3->Gcvolu()->number[geant3->Gcvolu()->nlevel - 4];
    //sQint[0] = qtot;
  }
  //
  pMC->Gmtod(sVectIn, sVloc, 1);
  if (TMath::Abs(sVloc[0]) >= sDlx) return;
  
  if (TMath::Abs(sVloc[2]) >= sDly) return;
  
  sVloc[0] += sDlx;
  x0 = sVloc[0];
  sVloc[2] += sDly;
  //y0 = sVloc[2];
  AnodicWires(y0a);
  ixmin = (Int_t) ((x0 - fSxcharge * 5.) / sDxp) + 1;
  if (ixmin < 1) ixmin = 1;
  ixmax = (Int_t) ((x0 + fSxcharge * 5.) / sDxp) + 1;
  if (ixmax > sNpx) ixmax = sNpx;
  iymin = (Int_t) ((y0a - sSycharge * 5.) / sDyp) + 1;
  if (iymin < 1) iymin = 1;
  iymax = (Int_t) ((y0a + sSycharge * 5.) / sDyp) + 1;
  if (iymax > sNpy) iymax = sNpy;
  qtot_check = 0.;
  i1 = ixmax;
  for (ix = ixmin; ix <= i1; ++ix) {
    i2 = iymax;
    for (iy = iymin; iy <= i2; ++iy) {
      xi1 = (sXpad[ix - 1] - x0) / sDpad;
      xi2 = xi1 + sDxp / sDpad;
      yi1 = (sYpad[iy - 1] - y0a) / sDpad;
      yi2 = yi1 + sDyp / sDpad;
      qp = qtot * FMathieson(xi1, xi2) * FMathieson(yi1, yi2);
      iqp = (Int_t) TMath::Abs(qp);
      qtot_check += TMath::Abs(qp);
      
      //     FILL COMMON FOR PADS 
      
      if (iqp) {
	if (iqp > 4095) {
	  iqp = 4096;
	}
	++sNpads;
	if (sNpads > MAXSEC) {
	  //                  write(6,*) 'attention npads',npads 
	  sNpads = MAXSEC;
	}
	sIpx[sNpads - 1] = ix;
	sIpy[sNpads - 1] = iy;
	sIqpad[sNpads - 1] = iqp;
	//     TAG THE ORIGIN OF THE CHARGE DEPOSITION 
	r1 = geant3->Gctrak()->upwght / 100.;
	ifeed = Int_t(r1+0.5);
	sNphlink[sNpads - 1] = 0;
	if (geant3->Gckine()->ipart != 50) {
	  //     MIP 
	  //isec[sNpads - 1] = sNsecon;
	} else {
	  if (ifeed < 50) {
	    //     CERENKOV PHOTON 
	    
	    nph = Int_t(geant3->Gctrak()->upwght+0.5);
	    sNphlink[sNpads - 1] = nph;
	    sXphit[nph - 1] = sVloc[0];
	    sYphit[nph - 1] = sVloc[2];
	    //isec[sNpads - 1] = sNsecon + 100;
	  } else if (ifeed == 51) {
	    //     FEEDBACK FROM CERENKOV 
	    
	    //isec[sNpads - 1] = sNsecon + 300;
	  } else if (ifeed == 52) {
	    //     FEEDBACK FROM MIP 
	    
	    //isec[sNpads - 1] = sNsecon + 200;
	  }
	}
	// Generate the hit for Root IO 
	for (ll = 0; ll < 25; ++ll) rrhh[ll] = 0;
	irivol[0] = modulen;
	irivol[1] = nmodsx;
	rrhh[0] = 0.;
	// Flag to say this is a  hit 
	rrhh[1] = xsx;
	rrhh[2] = ysx;
	rrhh[3] = psx;
	rrhh[4] = thetasx;
	rrhh[5] = phisx;
	rrhh[6] = (Float_t) idpartgx;
	rrhh[7] = sXsox;
	rrhh[8] = sYsox;
	rrhh[9] = sZsox;
	rrhh[10] = (Float_t) sIpx[sNpads - 1];
	rrhh[11] = (Float_t) sIpy[sNpads - 1];
	rrhh[12] = (Float_t) sIqpad[sNpads - 1];
	ihitrak = gAlice->CurrentTrack();
	if (sNphlink[sNpads - 1] > 0) {
	  rrhh[13] = sPparent;
	  rrhh[14] = sThincParent;
	  rrhh[15] = (Float_t) sIloss[nph - 1];
	  rrhh[16] = (Float_t) sIprod[nph - 1];
	  rrhh[17] = sXphit[nph - 1];
	  rrhh[18] = sYphit[nph - 1];
	  ihitrak = sMckov[nph - 1];
	}
	//     This is the current position of the hit 
	rrhh[19] = geant3->Gctrak()->vect[0];
	rrhh[20] = geant3->Gctrak()->vect[1];
	rrhh[21] = geant3->Gctrak()->vect[2];
	//               PRINT *,' ' 
	//               PRINT *,'RXAHIT(oldhit)=[',RRHH(20),',',RRHH(21),',', 
	//     +              RRHH(22),']' 
	//          PRINT *, 'ITRA = ',ITRA,'ISTAK = ',ISTAK 
	AddHit(ihitrak, irivol, rrhh);
	// new Padhits 
	for (ll = 0; ll < 25; ++ll) rrhh[ll] = 0;
	rrhh[0] = 3.;
	rrhh[1] = (Float_t) geant3->Gckine()->ipart;
	rrhh[2] = (Float_t) sIpx[sNpads - 1];
	rrhh[3] = (Float_t) sIpy[sNpads - 1];
	rrhh[4] = (Float_t) modulen;
	rrhh[5] = -1.;
	// !!!Prod 
	rrhh[6] = (Float_t) sIqpad[sNpads - 1];
	AddHit(ihitrak, irivol, rrhh);
	// enew 
      }
    }
  }
  //      IF (IPART .NE. 50) WRITE(6,*) 'CHECK',QTOT,QTOT_CHECK 
  
  //     GENERATE FEEDBACK PHOTONS 
  
  source[0] = x0 - sDlx;
  source[1] = sVloc[1] - .2;
  source[2] = y0a - sDly;
  pMC->Gdtom(source, source, 1);
  FeedBack(source, qtot);
  return;
}

//_____________________________________________________________________________
void AliRICH::AnodicWires(Float_t &y0a)
{
  //
  // Position of anodic wires
  //
  Int_t i1;
  
  // Local variables 
  Int_t i;
  Float_t ass_i, y0, ass_i1;
  
  y0 = sVloc[2];
  y0a = -1.;
  i1 = (sNpy << 1) + 1;
  for (i = 1; i <= i1; ++i) {
    if (y0 > sYanode[i - 1] && y0 <= sYanode[i]) {
      ass_i1 = TMath::Abs(sYanode[i] - y0);
      ass_i = TMath::Abs(sYanode[i - 1] - y0);
      if (ass_i1 <= ass_i) y0a = sYanode[i];
      if (ass_i1 > ass_i) y0a = sYanode[i - 1];
      return;
    }
  }
} 

//_____________________________________________________________________________
void AliRICH::GetChargeMip(Float_t &qtot)
{
  //
  // Calculates the charge deposited by a MIP
  //

  AliMC* pMC = AliMC::GetMC();

  Int_t i1;
  
  // Local variables 
  Float_t ranf[2];
  Int_t i;
  Int_t nel;
  
  //     NUMBER OF ELECTRONS 
  
  nel = (Int_t) (sDetot * 1e9 / 26.);
  qtot = 0.;
  if (!nel) nel = 1;
  i1 = nel;
  for (i = 1; i <= i1; ++i) {
    pMC->Rndm(ranf, 1);
    qtot -= fChslope * TMath::Log(ranf[0]);
  }
} 

//_____________________________________________________________________________
void AliRICH::GetCharge(Float_t &qtot)
{
  //
  // Charge deposited
  //

  AliMC* pMC = AliMC::GetMC();

  Float_t ranf[1];
  
  pMC->Rndm(ranf, 1);
  qtot = -fChslope * TMath::Log(ranf[0]);
} 

//_____________________________________________________________________________
void AliRICH::FeedBack(Float_t *source, Float_t qtot)
{
  //
  // Generate FeedBack photons
  //

  AliMC* pMC = AliMC::GetMC();
  TGeant3 *geant3 = (TGeant3*) pMC;

  Int_t i1, j;
  Float_t r1, r2;
  
  // Local variables 
  Float_t cthf, ranf[2], phif, enfp = 0, sthf;
  Int_t i, ifeed;
  Float_t e1[3], e2[3], e3[3];
  Float_t vmod, uswop;
  Float_t fp, random;
  Float_t dir[3], phi;
  Int_t nfp;
  Float_t pol[3], supwght;
  
  //     DETERMINE NUMBER OF FEEDBACK PHOTONS 
  
  // Function Body 
  fp = fAlphaFeed * qtot;
  nfp = gRandom->Poisson(fp);
  
  //     GENERATE PHOTONS 
  
  geant3->Gckin2()->ngphot = 0;
  i1 = nfp;
  for (i = 1; i <= i1; ++i) {
    
    //     DIRECTION 
    pMC->Rndm(ranf, 2);
    cthf = ranf[0] * 2. - 1.;
    if (cthf < 0.)  continue;
    sthf = TMath::Sqrt((1. - cthf) * (cthf + 1.));
    phif = ranf[1] * 2 * TMath::Pi();
    ++geant3->Gckin2()->ngphot;
    if (geant3->Gckin2()->ngphot > 800) {
      printf("ATTENTION NGPHOT ! %d %f %d\n",
	     geant3->Gckin2()->ngphot,fp,nfp);
      break;
    }
    
    //     POSITION 
    for(j=0;j<3;j++) 
      geant3->Gckin2()->xphot[geant3->Gckin2()->ngphot-1][j]=source[j];
    
    //     ENERGY 
    pMC->Rndm(&random, 1);
    if (random <= .57) {
      enfp = 7.5e-9;
    } else if (random > .57 && random <= .7) {
      enfp = 6.4e-9;
    } else if (random > .7) {
      enfp = 7.9e-9;
    }
    geant3->Gckin2()->xphot[geant3->Gckin2()->ngphot-1][6] = enfp;
    dir[0] = sthf * TMath::Sin(phif);
    dir[1] = cthf;
    dir[2] = sthf * TMath::Cos(phif);
    pMC->Gdtom(dir, &geant3->Gckin2()->xphot[geant3->Gckin2()->ngphot-1][3], 2);
    
    //     POLARISATION 
    e1[0] = 0.;
    e1[1] = -dir[2];
    e1[2] = dir[1];
    
    e2[0] = -dir[1];
    e2[1] = dir[0];
    e2[2] = 0.;
    
    e3[0] = dir[1];
    e3[1] = 0.;
    e3[2] = -dir[0];
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e1[j];
      e1[j]=e3[j];
      e3[j]=uswop;
    }
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    if (!vmod) for(j=0;j<3;j++) {
      uswop=e2[j];
      e2[j]=e3[j];
      e3[j]=uswop;
    }
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e1[j]*e1[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e1[j]*=vmod;
    
    vmod=0;
    for(j=0;j<3;j++) vmod+=e2[j]*e2[j];
    vmod=TMath::Sqrt(1/vmod);
    for(j=0;j<3;j++) e2[j]*=vmod;
    
    pMC->Rndm(ranf, 1);
    phi = ranf[0] * 2 * TMath::Pi();
    r1 = TMath::Sin(phi);
    r2 = TMath::Cos(phi);
    for(j=0;j<3;j++) pol[j]=e1[j]*r1+e2[j]*r2;
    pMC->Gdtom(pol, &geant3->Gckin2()->xphot[geant3->Gckin2()->ngphot-1][7], 2);
    
    //     TIME OF FLIGHT 
    geant3->Gckin2()->xphot[geant3->Gckin2()->ngphot-1][10] = 0.;
    
    //     PUT PHOTON ON THE STACK AND LABEL IT AS FEEDBACK (51,52) 
    supwght = geant3->Gctrak()->upwght;
    r1 = geant3->Gctrak()->upwght / 100.;
    ifeed = Int_t(r1+0.5);
    ++sNfeed;
    if (geant3->Gckine()->ipart == 50 && ifeed != 52) {
      geant3->Gctrak()->upwght = 5100.;
    } else {
      geant3->Gctrak()->upwght = 5200.;
    }
    geant3->Gskpho(geant3->Gckin2()->ngphot);
    geant3->Gctrak()->upwght = supwght;
    
  }
  geant3->Gckin2()->ngphot = 0;
}

//_____________________________________________________________________________
Float_t AliRICH::FMathieson(Float_t lambda1, Float_t lambda2)
{
  //
  // Mathieson charge distribution
  //
  // CALCULATES INTEGRAL OF GATTI/MATHIESON CHARGE DISTRIBUTION 
  // FOR CPC AND CSC CHAMBERS 
  // INTEGRATION LIMITS LAMBDA1,LAMBDA2=  X/D WHERE D IS THE ANODE CATHODE
  //                                          SEPARATION 
  // K3: GEOMETRY PARAMETER 
  //
  Float_t u1, u2;
  //    const Float_t k1=0.2828;
  const Float_t k2=0.962;
  const Float_t k3=0.6;
  const Float_t k4=0.379;
  const Float_t sqrtk3=TMath::Sqrt(k3);
  
  
  u1 = tanh(lambda1 * k2) * sqrtk3;
  u2 = tanh(lambda2 * k2) * sqrtk3;
  
  return (atan(u2) - atan(u1)) * 2 * k4;
  
}


//_____________________________________________________________________________
void AliRICH::UpdateMipHit(Float_t *hits)
{
  //
  // Update some parameters when the current mip hits the detector.
  //
  TClonesArray &lhits = *fMips;
  AliRICHmip *lmip = (AliRICHmip *) lhits.AddrAt(fNmips-1);
  lmip->SetX     (      hits[1]);
  lmip->SetY     (      hits[2]);
  lmip->SetModule((int) hits[3]);
  lmip->SetTheta (      hits[4]);
  lmip->SetPhi   (      hits[5]);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//_____________________________________________________________________________
void AliRICH::MakeBranch(Option_t *option)
{
  //
  // Create a new branch in the current Root Tree.
  // The branch of fHits is automatically split
  //
  AliDetector::MakeBranch(option);
  
  Int_t buffersize = 4000;
  if (gAlice->TreeH()) {
    if(fMips){
      gAlice->TreeH()->Branch("RICHmip",&fMips, buffersize);
      printf("Making Branch %s for mips\n","RICHmip");
    }
    if(fCkovs){
      gAlice->TreeH()->Branch("RICHckov",&fCkovs, buffersize);
      printf("Making Branch %s for ckovs\n","RICHckov");
    }
    if(fPadhits){
      gAlice->TreeH()->Branch("RICHpadhit",&fPadhits, buffersize);
      printf("Making Branch %s for padhits\n","RICHpadhit");
    }
  }	
}

//_____________________________________________________________________________
void AliRICH::SetTreeAddress()
{
  //
  // Set branch address for the Hits and Digits Tree.
  //
  AliDetector::SetTreeAddress();
  TBranch *branch;
  TTree *treeH = gAlice->TreeH();
  if (treeH) {
    if (fMips) {
      branch = treeH->GetBranch("RICHmip");
      if (branch) branch->SetAddress(&fMips);
    }
    if (fCkovs) {
      branch = treeH->GetBranch("RICHckov");
      if (branch) branch->SetAddress(&fCkovs);
    }
    if (fPadhits) {
      branch = treeH->GetBranch("RICHpadhit");
      if (branch) branch->SetAddress(&fPadhits);
    }
  }
}

//_____________________________________________________________________________
void AliRICH::ResetHits()
{
  //
  // Reset number of mips and the mips array for RICH
  //
  AliDetector::ResetHits();
  fNmips    = 0;
  if (fMips)    fMips->Clear();
  // Reset number of Ckovs and the Ckovs array for RICH
  fNckovs   = 0;
  if (fCkovs)   fCkovs->Clear();
  // Reset number of Padhits and the Padhits array for RICH
  fNpadhits = 0;
  if (fPadhits) fPadhits->Clear();
}

//_____________________________________________________________________________
Int_t AliRICH::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from mouse RICH on the display
  // Dummy routine
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliRICH::PreTrack()
{
  //
  // Called before transporting a track
  //
  sDetot=0;
  sNsecon=0;
  sNpads=0;
  sNphoton=0;
  sNfeed=0;
  sNfeedd=0;
  sNdir=0;
}

//_____________________________________________________________________________
void AliRICH::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliRICH.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    AliDetector::Streamer(R__b);
    R__b >> fNmips;
    R__b >> fNckovs;
    R__b >> fNpadhits;
    R__b >> fMips; //diff
    R__b >> fCkovs; //diff
    R__b >> fPadhits; //diff
    R__b >> fChslope;
    R__b >> fAlphaFeed;
    R__b >> fSxcharge;
    R__b >> fIritri;
  } else {
    R__b.WriteVersion(AliRICH::IsA());
    AliDetector::Streamer(R__b);
    R__b << fNmips;
    R__b << fNckovs;
    R__b << fNpadhits;
    R__b << fMips; //diff
    R__b << fCkovs; //diff
    R__b << fPadhits; //diff
    R__b << fChslope;
    R__b << fAlphaFeed;
    R__b << fSxcharge;
    R__b << fIritri;
  }
}

ClassImp(AliRICHv1)

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Ring Imaging Cherenkov                                                   //
//  This class contains version 1 of the Ring Imaging Cherenkov              //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliRICHv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliRICHv1::AliRICHv1() : AliRICH()
{
  //
  // Default Constructor RICH for version 1
  //
}
 
//_____________________________________________________________________________
AliRICHv1::AliRICHv1(const char *name, const char *title)
       : AliRICH(name,title)
{
  //
  // Standard constructor for RICH version 1
  //
}

//_____________________________________________________________________________
AliRICHv1::~AliRICHv1()
{
  //
  // Standard destructor for RICH version 1
  //
}
 
//_____________________________________________________________________________
void AliRICHv1::CreateGeometry()
{
  //
  // Create the geometry for RICH version 1
  //
  // Modified by:  N. Colonna (INFN - BARI, Nicola.Colonna@ba.infn.it) 
  //               R.A. Fini  (INFN - BARI, Rosanna.Fini@ba.infn.it) 
  //               R.A. Loconsole (Bari University, loco@riscom.ba.infn.it) 
  //
  //Begin_Html
  /*
    <img src="gif/AliRICHv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliRICHv1Tree.gif">
  */
  //End_Html

  AliMC* pMC = AliMC::GetMC();

  Int_t *idtmed = gAlice->Idtmed();
  
  Int_t i;
  Float_t zs;
  Int_t idrotm[1099];
  Float_t par[3];
  
  // --- Define the RICH detector 
  //     External aluminium box 
  par[0] = 71.1;
  par[1] = 11.5;
  par[2] = 73.15;
  pMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
  
  //     Sensitive part of the whole RICH 
  par[0] = 64.8;
  par[1] = 11.5;
  par[2] = 66.55;
  pMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3);
  
  //     Honeycomb 
  par[0] = 63.1;
  par[1] = .188;
  par[2] = 66.55;
  pMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
  
  //     Aluminium sheet 
  par[0] = 63.1;
  par[1] = .025;
  par[2] = 66.55;
  pMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
  
  //     Quartz 
  par[0] = 63.1;
  par[1] = .25;
  par[2] = 65.5;
  pMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
  
  //     Spacers (cylinders) 
  par[0] = 0.;
  par[1] = .5;
  par[2] = .5;
  pMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);
  
  //     Opaque quartz 
  par[0] = 61.95;
  par[1] = .2;
  par[2] = 66.5;
  pMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
  
  //     Frame of opaque quartz 
  par[0] = 20.65;
  par[1] = .5;
  par[2] = 66.5;
  pMC->Gsvolu("OQUF", "BOX ", idtmed[1007], par, 3);
  
  //     Little bar of opaque quartz 
  par[0] = 63.1;
  par[1] = .25;
  par[2] = .275;
  pMC->Gsvolu("BARR", "BOX ", idtmed[1007], par, 3);
  
  //     Freon 
  par[0] = 20.15;
  par[1] = .5;
  par[2] = 65.5;
  pMC->Gsvolu("FREO", "BOX ", idtmed[1003], par, 3);
  
  //     Methane 
  par[0] = 64.8;
  par[1] = 5.;
  par[2] = 64.8;
  pMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
  
  //     Methane gap 
  par[0] = 64.8;
  par[1] = .2;
  par[2] = 64.8;
  pMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
  
  //     CsI photocathode 
  par[0] = 64.8;
  par[1] = .25;
  par[2] = 64.8;
  pMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
  
  //     Anode grid 
  par[0] = 0.;
  par[1] = .0025;
  par[2] = 20.;
  pMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);
  
  // --- Places the detectors defined with GSVOLU 
  //     Place material inside RICH 
  pMC->Gspos("SRIC", 1, "RICH", 0., 0., 0., 0, "ONLY");
  
  pMC->Gspos("ALUM", 1, "SRIC", 0., -6.075, 0., 0, "ONLY");
  pMC->Gspos("HONE", 1, "SRIC", 0., -5.862, 0., 0, "ONLY");
  pMC->Gspos("ALUM", 2, "SRIC", 0., -5.649, 0., 0, "ONLY");
  pMC->Gspos("OQUA", 1, "SRIC", 0., -5.424, 0., 0, "ONLY");
  
  AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
  
  for (i = 1; i <= 9; ++i) {
    zs = (5 - i) * 14.4;
    pMC->Gspos("SPAC", i, "FREO", 6.7, 0., zs, idrotm[1019], "ONLY");
  }
  for (i = 10; i <= 18; ++i) {
    zs = (14 - i) * 14.4;
    pMC->Gspos("SPAC", i, "FREO", -6.7, 0., zs, idrotm[1019], "ONLY");
  }
  
  pMC->Gspos("FREO", 1, "OQUF", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("OQUF", 1, "SRIC", 41.3, -4.724, 0., 0, "ONLY");
  pMC->Gspos("OQUF", 2, "SRIC", 0., -4.724, 0., 0, "ONLY");
  pMC->Gspos("OQUF", 3, "SRIC", -41.3, -4.724, 0., 0, "ONLY");
  pMC->Gspos("BARR", 1, "QUAR", 0., 0., -21.65, 0, "ONLY");
  pMC->Gspos("BARR", 2, "QUAR", 0., 0., 21.65, 0, "ONLY");
  pMC->Gspos("QUAR", 1, "SRIC", 0., -3.974, 0., 0, "ONLY");
  pMC->Gspos("GAP ", 1, "META", 0., 4.8, 0., 0, "ONLY");
  pMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
  pMC->Gspos("CSI ", 1, "SRIC", 0., 6.526, 0., 0, "ONLY");
  
  //     Place RICH inside ALICE apparatus 
  
  AliMatrix(idrotm[1000], 90., 0., 70.69, 90., 19.31, -90.);
  AliMatrix(idrotm[1001], 90., -20., 90., 70., 0., 0.);
  AliMatrix(idrotm[1002], 90., 0., 90., 90., 0., 0.);
  AliMatrix(idrotm[1003], 90., 20., 90., 110., 0., 0.);
  AliMatrix(idrotm[1004], 90., 340., 108.2, 70., 18.2, 70.);
  AliMatrix(idrotm[1005], 90., 0., 109.31, 90., 19.31, 90.);
  AliMatrix(idrotm[1006], 90., 20., 108.2, 110., 18.2, 110.);
  
  pMC->Gspos("RICH", 1, "ALIC", 0., 471.9, 165.26,     idrotm[1000], "ONLY");
  pMC->Gspos("RICH", 2, "ALIC", 171., 470., 0.,        idrotm[1001], "ONLY");
  pMC->Gspos("RICH", 3, "ALIC", 0., 500., 0.,          idrotm[1002], "ONLY");
  pMC->Gspos("RICH", 4, "ALIC", -171., 470., 0.,       idrotm[1003], "ONLY");
  pMC->Gspos("RICH", 5, "ALIC", 161.4, 443.4, -165.3,  idrotm[1004], "ONLY");
  pMC->Gspos("RICH", 6, "ALIC", 0., 471.9, -165.3,     idrotm[1005], "ONLY");
  pMC->Gspos("RICH", 7, "ALIC", -161.4, 443.4, -165.3, idrotm[1006], "ONLY");
  
}

//_____________________________________________________________________________
void AliRICHv1::DrawModule()
{
  //
  // Draw a shaded view of the Ring Imaging Cherenkov
  //

  AliMC* pMC = AliMC::GetMC();
  TGeant3 *geant3 = (TGeant3*) pMC;

  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible

  pMC->Gsatt("RICH","seen",0);
  pMC->Gsatt("SRIC","seen",0);
  pMC->Gsatt("HONE","seen",1);
  pMC->Gsatt("ALUM","seen",1);
  pMC->Gsatt("QUAR","seen",1);
  pMC->Gsatt("SPAC","seen",1);
  pMC->Gsatt("OQUA","seen",1);
  pMC->Gsatt("OQUF","seen",1);
  pMC->Gsatt("BARR","seen",1);
  pMC->Gsatt("FREO","seen",1);
  pMC->Gsatt("META","seen",1);
  pMC->Gsatt("GAP ","seen",1);
  pMC->Gsatt("CSI ","seen",1);
  pMC->Gsatt("GRID","seen",1);
  
  //
  geant3->Gdopt("hide", "on");
  geant3->Gdopt("shad", "on");
  geant3->Gsatt("*", "fill", 7);
  geant3->SetClipBox(".");
  geant3->Gdopt("hide", "on");
  geant3->Gdopt("shad", "on");
  geant3->Gsatt("*", "fill", 7);
  geant3->SetClipBox(".");
  //  geant3->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  geant3->DefaultRange();
  pMC->Gdraw("alic", 60, 50, 0, 10, 0, .03, .03);
  geant3->Gdhead(1111, "Ring Imaging Cherenkov version 1");
  geant3->Gdman(16, 6, "MAN");
  geant3->Gdopt("hide", "off");
}

//_____________________________________________________________________________
void AliRICHv1::CreateMaterials()
{
  //
  // *** DEFINITION OF AVAILABLE RICH MATERIALS *** 
  // ORIGIN    : NICK VAN EIJNDHOVEN 
  // Modified by:  N. Colonna (INFN - BARI, Nicola.Colonna@ba.infn.it) 
  //               R.A. Fini  (INFN - BARI, Rosanna.Fini@ba.infn.it) 
  //               R.A. Loconsole (Bari University, loco@riscom.ba.infn.it) 
  //
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t ppckov[14] = { 5.63e-9,5.77e-9,5.9e-9,6.05e-9,6.2e-9,6.36e-9,6.52e-9,
			 6.7e-9,6.88e-9,7.08e-9,7.3e-9,7.51e-9,7.74e-9,8e-9 };
  Float_t rindex_quarz[14] = { 1.528309,1.533333,
			       1.538243,1.544223,1.550568,1.55777,
			       1.565463,1.574765,1.584831,1.597027,
			       1.611858,1.6277,1.6472,1.6724 };
  Float_t rindex_quarzo[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
  Float_t rindex_methane[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
  Float_t rindex_gri[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
  Float_t absco_freon[14] = { 179.0987,179.0987,
			      179.0987,179.0987,179.0987,35.7,12.54,5.92,4.92,3.86,1.42,.336,.134,0. };
  Float_t absco_quarz[14] = { 20.126,16.27,13.49,11.728,9.224,8.38,7.44,7.17,
			      6.324,4.483,1.6,.323,.073,0. };
  Float_t absco_quarzo[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
			       1e-5,1e-5,1e-5,1e-5,1e-5 };
  Float_t absco_csi[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			    1e-4,1e-4,1e-4,1e-4 };
  Float_t absco_methane[14] = { 1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,
				1e6,1e6,1e6 };
  Float_t absco_gri[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			    1e-4,1e-4,1e-4,1e-4 };
  Float_t effic_all[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
  Float_t effic_csi[14] = { 4.74e-4,.00438,.009,.0182,.0282,.0653,.1141,.163,
			    .2101,.2554,.293,.376,.3861,.418 };
  Float_t effic_gri[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
  
  Float_t afre[2], agri, amet[2], aqua[2], ahon, zfre[2], zgri, zhon, 
    zmet[2], zqua[2];
  Int_t nlmatfre;
  Float_t densquao;
  Int_t nlmatmet, nlmatqua;
  Float_t wmatquao[2], rindex_freon[14];
  Int_t i;
  Float_t aquao[2], epsil, stmin, zquao[2];
  Int_t nlmatquao;
  Float_t radlal, densal, tmaxfd, deemax, stemax;
  Float_t aal, zal, radlgri, densfre, radlhon, densgri, denshon,densqua, densmet, wmatfre[2], wmatmet[2], wmatqua[2];
  
  Int_t *idtmed = gAlice->Idtmed();
  
  AliMC* pMC = AliMC::GetMC();
  TGeant3 *geant3 = (TGeant3*) pMC;

  // --- Photon energy (GeV) 
  // --- Refraction indexes 
  for (i = 0; i < 14; ++i) {
    rindex_freon[i] = ppckov[i] * .01095 * 1e9 + 1.2177;
  }
  // need to be changed 
  
  // --- Absorbtion lenghts (in cm) 
  //      DATA ABSCO_QUARZ / 
  //     &    5 * 1000000., 
  //     &    29.85,    7.34,     4.134,    1.273,    0.722, 
  //     &    0.365, 0.365, 0.365, 0.  / 
  // need to be changed 
  
  // --- Detection efficiencies (quantum efficiency for CsI) 
  // --- Define parameters for honeycomb. 
  //     Used carbon of equivalent rad. lenght 
  
  ahon    = 12.01;
  zhon    = 6.;
  denshon = 2.265;
  radlhon = 18.8;
  
  // --- Parameters to include in GSMIXT, relative to Quarz (SiO2) 
  
  aqua[0]    = 28.09;
  aqua[1]    = 16.;
  zqua[0]    = 14.;
  zqua[1]    = 8.;
  densqua    = 2.64;
  nlmatqua   = -2;
  wmatqua[0] = 1.;
  wmatqua[1] = 2.;
  
  // --- Parameters to include in GSMIXT, relative to opaque Quarz (SiO2) 
  
  aquao[0]    = 28.09;
  aquao[1]    = 16.;
  zquao[0]    = 14.;
  zquao[1]    = 8.;
  densquao    = 2.64;
  nlmatquao   = -2;
  wmatquao[0] = 1.;
  wmatquao[1] = 2.;
  
  // --- Parameters to include in GSMIXT, relative to Freon (C6F14) 
  
  afre[0]    = 12.;
  afre[1]    = 19.;
  zfre[0]    = 6.;
  zfre[1]    = 9.;
  densfre    = 1.7;
  nlmatfre   = -2;
  wmatfre[0] = 6.;
  wmatfre[1] = 14.;
  
  // --- Parameters to include in GSMIXT, relative to methane (CH4) 
  
  amet[0]    = 12.01;
  amet[1]    = 1.;
  zmet[0]    = 6.;
  zmet[1]    = 1.;
  densmet    = 7.17e-4;
  nlmatmet   = -2;
  wmatmet[0] = 1.;
  wmatmet[1] = 4.;
  
  // --- Parameters to include in GSMIXT, relative to anode grid (Cu) 
  
  agri    = 63.54;
  zgri    = 29.;
  densgri = 8.96;
  radlgri = 1.43;
  
  // --- Parameters to include in GSMATE related to aluminium sheet 
  
  aal    = 26.98;
  zal    = 13.;
  densal = 2.7;
  radlal = 8.9;
  
  AliMaterial(1, "Air     $", 14.61, 7.3, .001205, 30420., 67500);
  AliMaterial(6, "HON", ahon, zhon, denshon, radlhon, 0);
  AliMaterial(16, "CSI", ahon, zhon, denshon, radlhon, 0);
  AliMixture(20, "QUA", aqua, zqua, densqua, nlmatqua, wmatqua);
  AliMixture(21, "QUAO", aquao, zquao, densquao, nlmatquao, wmatquao);
  AliMixture(30, "FRE", afre, zfre, densfre, nlmatfre, wmatfre);
  AliMixture(40, "MET", amet, zmet, densmet, nlmatmet, wmatmet);
  AliMixture(41, "METG", amet, zmet, densmet, nlmatmet, wmatmet);
  AliMaterial(11, "GRI", agri, zgri, densgri, radlgri, 0);
  AliMaterial(50, "ALUM", aal, zal, densal, radlal, 0);
  
  tmaxfd = -10.;
  stemax = -.1;
  deemax = -.2;
  epsil  = .001;
  stmin  = -.001;
  
  AliMedium(1001, "DEFAULT MEDIUM AIR$", 1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1002, "HONEYCOMB$", 6, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1003, "QUARZO$", 20, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1004, "FREON$", 30, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1005, "METANO$", 40, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1006, "CSI$", 16, 1, ISXFLD, SXMGMX,tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1007, "GRIGLIA$", 11, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1008, "QUARZOO$", 21, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1009, "GAP$", 41, 1, ISXFLD, SXMGMX,tmaxfd, .1, -deemax, epsil, -stmin);
  AliMedium(1010, "ALUMINUM$", 50, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  
  //     Switch on delta-ray production in the methane and freon gaps 
  
  pMC->Gstpar(idtmed[1002], "LOSS", 1.);
  pMC->Gstpar(idtmed[1003], "LOSS", 1.);
  pMC->Gstpar(idtmed[1004], "LOSS", 1.);
  pMC->Gstpar(idtmed[1008], "LOSS", 1.);
  pMC->Gstpar(idtmed[1005], "LOSS", 1.);
  pMC->Gstpar(idtmed[1002], "HADR", 1.);
  pMC->Gstpar(idtmed[1003], "HADR", 1.);
  pMC->Gstpar(idtmed[1004], "HADR", 1.);
  pMC->Gstpar(idtmed[1008], "HADR", 1.);
  pMC->Gstpar(idtmed[1005], "HADR", 1.);
  pMC->Gstpar(idtmed[1002], "DCAY", 1.);
  pMC->Gstpar(idtmed[1003], "DCAY", 1.);
  pMC->Gstpar(idtmed[1004], "DCAY", 1.);
  pMC->Gstpar(idtmed[1008], "DCAY", 1.);
  pMC->Gstpar(idtmed[1005], "DCAY", 1.);
  geant3->Gsckov(idtmed[1000], 14, ppckov, absco_methane, effic_all, rindex_methane);
  geant3->Gsckov(idtmed[1001], 14, ppckov, absco_methane, effic_all, rindex_methane);
  geant3->Gsckov(idtmed[1002], 14, ppckov, absco_quarz, effic_all,rindex_quarz);
  geant3->Gsckov(idtmed[1003], 14, ppckov, absco_freon, effic_all,rindex_freon);
  geant3->Gsckov(idtmed[1004], 14, ppckov, absco_methane, effic_all, rindex_methane);
  geant3->Gsckov(idtmed[1005], 14, ppckov, absco_csi, effic_csi, rindex_methane);
  geant3->Gsckov(idtmed[1006], 14, ppckov, absco_gri, effic_gri, rindex_gri);
  geant3->Gsckov(idtmed[1007], 14, ppckov, absco_quarzo, effic_all, rindex_quarzo);
  geant3->Gsckov(idtmed[1008], 14, ppckov, absco_methane, effic_all, rindex_methane);
  geant3->Gsckov(idtmed[1009], 14, ppckov, absco_gri, effic_gri, rindex_gri);
}

ClassImp(AliRICHhit)
  
//_____________________________________________________________________________
AliRICHhit::AliRICHhit(Int_t shunt, Int_t track, Int_t *vol,
		       Float_t *hits, Int_t fNpadhits) :
    AliHit(shunt,track)
{
  //
  // Standard constructor for RICH hit
  //
  Int_t i;
  for (i=0;i<2;i++) fVolume[i] = vol[i];
  fTrack      = (int)  track;
  //Hit information
  fPart       = (int)  hits[ 1]; 
  //AliHit information, position of the hit
  fX          =        hits[ 2];       
  fY          =        hits[ 3];
  //Pad information
  fFirstpad   = fNpadhits;
  fLastpad    = fNpadhits-1;
  
  //Hit information
  fModule     = (int)  hits[ 4];
  fTheta      =        hits[ 5];
  fArrivaltime=        hits[ 6];
  fFeed       = (int)  hits[ 7];
}

ClassImp(AliRICHmip)

//_____________________________________________________________________________
AliRICHmip::AliRICHmip(Int_t shunt, Int_t track, Int_t *vol, 
		       Float_t *hits, Int_t fNckovs, Int_t fNpadhits) :
    AliRICHhit(shunt,track,vol,hits,fNpadhits)
{
  //
  // Standard constructor for RICH Mip hit
  //
  fPhi        =        hits[ 8];  
  fPs         =        hits[ 9];
  fQ          =        hits[10];    
  fZ          =        hits[11];
  //Ckov information
  if ((int) hits[12]){
    fFirstCkov  =    fNckovs;
    fLastCkov   =    fNckovs-1;
  } else {
    fFirstCkov  =    -1;
    fLastCkov   =    -1;
  }
}

ClassImp(AliRICHckov)
  
//_____________________________________________________________________________
AliRICHckov::AliRICHckov(Int_t shunt, Int_t track, Int_t *vol, 
			 Float_t *hits, Int_t fNmips, Int_t fNpadhits) :
  AliRICHhit(shunt,track,vol,hits,fNpadhits)
{
  //
  // Standard creator for RICH Cherenkov information
  //
  fEnergy     =        hits[8];
  fStop       = (int)  hits[9];
  
  //Parent info
  fParent     =         fNmips;
}

ClassImp(AliRICHpadhit)


//_____________________________________________________________________________
AliRICHpadhit::AliRICHpadhit(Int_t shunt, Int_t track, Int_t *vol,
			     Float_t *hits, Int_t fNmips, Int_t fNckovs):
  AliHit(shunt,track)
{
  //
  // Standard constructor for RICH pad hit
  //
  Int_t i;
  for (i=0;i<2;i++) fVolume[i] = vol[i];
  
  // Hit information
  fX          = (int)  hits[ 2];
  fY          = (int)  hits[ 3];
  fModule     = (int)  hits[ 4];
  fParentMip  =        fNmips;
  fParentCkov =        fNckovs;
  fProd       = (int)  hits[ 5];
  fCharge     =        hits[ 6];
  fZ          =          -999.0; // Not implemented
}


