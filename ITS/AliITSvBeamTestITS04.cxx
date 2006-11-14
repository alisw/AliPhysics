////////////////////////////////////////////////////////
// ITS geometry class and step manager for the        //
//   integrated ITS test beam of Nov. 04              //
//  Author: mercedes.lopez.noriega@cern.ch            //
////////////////////////////////////////////////////////
#include "AliRun.h"
#include "AliITSvBeamTestITS04.h"
#include <TClonesArray.h>
#include <TString.h>
#include "AliITS.h"
#include "AliITSDetTypeSim.h"
#include "AliITSgeom.h"
#include "AliITShit.h"
#include "AliITSresponseSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSSD.h"
#include "AliMagF.h"
#include "TVirtualMC.h"
#include "AliMC.h"

const Int_t AliITSvBeamTestITS04::fgkNumberOfSPD = 4;
const Int_t AliITSvBeamTestITS04::fgkNumberOfSDD = 2;
const Int_t AliITSvBeamTestITS04::fgkNumberOfSSD = 4;
// Dimension (thickness:Y (beam direction), width:X, length:Z)

const char*    AliITSvBeamTestITS04::fgSPDsensitiveVolName = "ITSspdSensitiv";
//dimensions (preliminary values from Petra (in cms))
const Double_t AliITSvBeamTestITS04::fgkSPDthickness    = 0.02;
const Double_t AliITSvBeamTestITS04::fgkSPDwidth        = 1.4; 
const Double_t AliITSvBeamTestITS04::fgkSPDlength       = 7.2;
const Double_t AliITSvBeamTestITS04::fgkSPDthickSens    = 0.02;
const Double_t AliITSvBeamTestITS04::fgkSPDwidthSens    = 1.2; 
const Double_t AliITSvBeamTestITS04::fgkSPDlengthSens   = 7.0;
//position
const Double_t AliITSvBeamTestITS04::fgkSPD0y = 23.7;
const Double_t AliITSvBeamTestITS04::fgkSPD1y = 33.7;

//===
const char*    AliITSvBeamTestITS04::fgSDDsensitiveVolName = "ITSsddSensitiv";
//dimensions (preliminary values from Ludovic (in cms))
const Double_t AliITSvBeamTestITS04::fgkSDDthickness     = 0.03;
const Double_t AliITSvBeamTestITS04::fgkSDDwidth         = 7.22;
const Double_t AliITSvBeamTestITS04::fgkSDDlength        = 8.76;
const Double_t AliITSvBeamTestITS04::fgkSDDthickSens     = 0.02998;
const Double_t AliITSvBeamTestITS04::fgkSDDwidthSens     = 7.017;
const Double_t AliITSvBeamTestITS04::fgkSDDlengthSens    = 7.497;
//position
const Double_t AliITSvBeamTestITS04::fgkSDD0y = 51.7;
const Double_t AliITSvBeamTestITS04::fgkSDD1y = 57.2;

//===
const char*    AliITSvBeamTestITS04::fgSSDsensitiveVolName = "ITSssdSensitiv";
//dimensions (final values from Javier (in cms))
const Double_t AliITSvBeamTestITS04::fgkSSDthickness    = 0.03;
const Double_t AliITSvBeamTestITS04::fgkSSDwidth        = 7.7;
const Double_t AliITSvBeamTestITS04::fgkSSDlength       = 4.4;
const Double_t AliITSvBeamTestITS04::fgkSSDthickSens    = 0.03;
const Double_t AliITSvBeamTestITS04::fgkSSDwidthSens    = 7.5;
const Double_t AliITSvBeamTestITS04::fgkSSDlengthSens   = 4.2;
//position
const Double_t AliITSvBeamTestITS04::fgkSSD0y = 73.6;
const Double_t AliITSvBeamTestITS04::fgkSSD1y = 80.6;

//===============================================================


#include <TLorentzVector.h>
#include "AliTrackReference.h"
#include "AliITSDetTypeSim.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSSD.h"
#include "AliITShit.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSCalibrationSPD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSSD.h"

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>

ClassImp(AliITSvBeamTestITS04)
    
//_____________________________________________________________
  AliITSvBeamTestITS04::AliITSvBeamTestITS04() : AliITS(),
fITSmotherVolume(0),
fNspd(fgkNumberOfSPD),
fNsdd(fgkNumberOfSDD),
fNssd(fgkNumberOfSSD),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE){
    //
    // Constructor
    //
    
  // SetNumberOfSPD(fgkNumberOfSPD);
  //  SetNumberOfSDD(fgkNumberOfSDD);
  //  SetNumberOfSSD(fgkNumberOfSSD);
    
    fIdN = 3;         
    fIdName    = new TString[fIdN];
    fIdName[0] = fgSPDsensitiveVolName;
    fIdName[1] = fgSDDsensitiveVolName;
    fIdName[2] = fgSSDsensitiveVolName;
    fIdSens    = new Int_t[fIdN];
    for(Int_t i=0; i<fIdN; i++) fIdSens[i] = 0;
    
    //for writing out geometry
    //fGeomDetOut   = kFALSE; 

    // for reading in geometry (JC)
    //fGeomDetIn = kFALSE;

    for(Int_t a=0;a<60;a++) fWrite[a] = '\0';
}

//_____________________________________________________________
AliITSvBeamTestITS04::AliITSvBeamTestITS04(const char* name,const char *title)
  : AliITS(name,title),
fITSmotherVolume(0),
fNspd(fgkNumberOfSPD),
fNsdd(fgkNumberOfSDD),
fNssd(fgkNumberOfSSD),
fGeomDetOut(kFALSE),
fGeomDetIn(kFALSE)
{
    //
    // Constructor
    //
    
    //SetNumberOfSPD(fgkNumberOfSPD);
    //SetNumberOfSDD(fgkNumberOfSDD);
    //SetNumberOfSSD(fgkNumberOfSSD);
    
    fIdN = 3;         
    fIdName    = new TString[fIdN];
    fIdName[0] = fgSPDsensitiveVolName;
    fIdName[1] = fgSDDsensitiveVolName;
    fIdName[2] = fgSSDsensitiveVolName;
    fIdSens    = new Int_t[fIdN];
    for(Int_t i=0; i<fIdN; i++) fIdSens[i] = 0;
    
    //for writing out geometry
    //fGeomDetOut   = kFALSE; // Don't write .det file

    // for reading in geometry (JC)
    //fGeomDetIn = kFALSE;

    for(Int_t a=0;a<60;a++) fWrite[a] = '\0';    
}

//__________________________________________________________________
AliITSvBeamTestITS04::~AliITSvBeamTestITS04()
{
    //
    // Destructor
    //
}

//______________________________________________________________________
void AliITSvBeamTestITS04::CreateMaterials()
{
    // Media defined here should correspond to the one defined in galice.cuts
    // This file is read in (AliMC*) fMCApp::Init() { ReadTransPar(); }
    
    // Create ITS materials
    Int_t   ifield = gAlice->Field()->Integ();
    Float_t fieldm = gAlice->Field()->Max();
    
    Float_t tmaxfdSi = 0.1;
    Float_t stemaxSi = 0.0075;
    Float_t deemaxSi = 0.1;
    Float_t epsilSi  = 1.0E-4;
    Float_t stminSi  = 0.0;
    
    Float_t tmaxfdAir = 0.1;
    Float_t stemaxAir = .10000E+01;
    Float_t deemaxAir = 0.1;
    Float_t epsilAir  = 1.0E-4;
    Float_t stminAir  = 0.0;
 
    // AIR
    Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
    Float_t zAir[4]={6.,7.,8.,18.};
    Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
    Float_t dAir = 1.20479E-3;
    
    AliMaterial(51,"ITSspdSi",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(51,"ITSspdSi",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    
    AliMaterial(1,"ITSsddSi",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    AliMedium(1,"ITSsddSi",1,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    
    //AliMaterial(?,"ITSssdSi",0.28086E+02,0.14000E+02,0.23300E+01,0.93600E+01,0.99900E+03);
    //AliMedium(?,"ITSssdSi",51,0,ifield,fieldm,tmaxfdSi,stemaxSi,deemaxSi,epsilSi,stminSi);
    
    AliMixture(5,"ITSair",aAir,zAir,dAir,4,wAir);
    AliMedium(5,"ITSair",5,0,ifield,fieldm,tmaxfdAir,stemaxAir,deemaxAir,epsilAir,stminAir);
    
//NEED TO ADD PLASTIC OF SCINTILLATORS!!

}

//______________________________________________________________________
void AliITSvBeamTestITS04::CreateGeometry()
{    
  //Creates geometry
    TGeoManager *geoManager = gGeoManager;
    TGeoVolume *vALIC = geoManager->GetTopVolume();
    
    //================================
    //       ITS mother volume
    //================================
    TGeoPcon *sITS = new TGeoPcon("ITS Top Volume",0.0,360.0,2);
    // DefineSection(section number, Z, Rmin, Rmax).
    sITS->DefineSection(0,-100.0,0.01,100.0); // Units in cms
    sITS->DefineSection(1,+100.0,0.01,100.0);
    
    TGeoMedium *air = gGeoManager->GetMedium("ITSair");
    fITSmotherVolume = new TGeoVolume("vITS",sITS,air);
    fITSmotherVolume->SetVisibility(kFALSE);
    vALIC->AddNode(fITSmotherVolume,1,0);
    
//     //Scintillators
//     TGeoMedium *plasticScint = new TGeoMedium("plasticScint",1,Plastic);
//     //First Scintillator
//     TGeoBBox *Scint1Shape = new TGeoBBox("Scint1Shape",0.5,0.1,0.5,0); //1x1cm
//     TGeoVolume *Scint1 = new TGeoVolume("Scint1",Scint1Shape,plasticScint);
//     TGeoTranslation *firstScint = new TGeoTranslation(0,0.7,0);
//     vALIC->AddNode(Scint1,2,firstScint);
//     //Second Scintillator
//     TGeoBBox *Scint2Shape = new TGeoBBox("Scint2Shape",1.,0.1,1.,0); //2x2cm
//     TGeoVolume *Scint2 = new TGeoVolume("Scint2",Scint2Shape,plasticScint);
//     TGeoTranslation *secondScint = new TGeoTranslation(0,90.,0);
//     vALIC->AddNode(Scint2,3,secondScint);
    
    AddSPDGeometry(fITSmotherVolume);
    AddSDDGeometry(fITSmotherVolume);
    AddSSDGeometry(fITSmotherVolume);
}

//______________________________________________________________________
void AliITSvBeamTestITS04::Init()
{
    // Initialize the ITS after it has been created.
    Int_t i;
    for(i=0;i<20;i++) printf("*");
    printf( " ITSbeamtest_Init " );
    for(i=0;i<20;i++) printf("*"); printf("\n");

//    // Create geometry
//    if(!fGeomDetIn) this->InitAliITSgeom();

    // Initialize AliITS
    AliITS::Init();
    for(i=0;i<40+16;i++) printf("*"); printf("\n");

}

//______________________________________________________________________
void AliITSvBeamTestITS04::InitAliITSgeom()
{    
  //initialisation of ITSgeom
    const Int_t knlayers = 6;
    Int_t nlad[knlayers], ndet[knlayers];
    
    nlad[0] = 1; ndet[0] = 2;
    nlad[1] = 1; ndet[1] = 2;
    nlad[2] = 1; ndet[2] = 1;
    nlad[3] = 1; ndet[3] = 1;
    nlad[4] = 1; ndet[4] = 2;
    nlad[5] = 1; ndet[5] = 2;

    Int_t nModTot = fNspd + fNsdd + fNssd;
    if (GetITSgeom()) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom(0,knlayers,nlad,ndet,nModTot);
    SetITSgeom(geom);
    //*** Set default shapes 
    const Float_t kDxyzSPD[] = {fgkSPDwidthSens/2, fgkSPDthickSens/2,fgkSPDlengthSens/2};  
    if(!(GetITSgeom()->IsShapeDefined(kSPD)))
	GetITSgeom()->ReSetShape(kSPD,new AliITSgeomSPD425Short(3,(Float_t *)kDxyzSPD));
    
    const Float_t kDxyzSDD[] = {fgkSDDwidthSens/2., fgkSDDthickSens/2.,fgkSDDlengthSens/2.};
    if(!(GetITSgeom()->IsShapeDefined(kSDD)))
	GetITSgeom()->ReSetShape(kSDD, new AliITSgeomSDD256(3,(Float_t *)kDxyzSDD));
    
    const Float_t kDxyzSSD[] = {fgkSSDlengthSens/2, fgkSSDthickSens/2,fgkSSDwidthSens/2};
    if(!(GetITSgeom()->IsShapeDefined(kSSD)))
	GetITSgeom()->ReSetShape(kSSD,new AliITSgeomSSD75and275(3,(Float_t *)kDxyzSSD));
    
    // Creating the matrices in AliITSgeom for each sensitive volume
    // (like in AliITSv11GeometrySDD) mln
    // Here, each layer is one detector
    
    char layerName[30];
    Int_t startMod = 0;
    
    // SPD
    for (Int_t i=0; i<fNspd;i++) {
	sprintf(layerName, "ITSspdWafer_%i",i+1);
	TGeoNode *layNode = fITSmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());	    
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    if ((i+1==2)||(i+1==4)) iDet = 2;
	    Int_t iLad = 1;
	    Int_t iLay = 1;
	    if (i+1>2) iLay = 2;
	    GetITSgeom()->CreateMatrix(startMod,iLay,iLad,iDet,kSPD,trans,rot);
	    startMod++;
	};
    };
    
    // SDD
    for (Int_t i=0; i<fNsdd;i++) {
	sprintf(layerName, "ITSsddWafer_%i",i+fNspd+1);
	TGeoNode *layNode = fITSmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    Int_t iLad = 1;
	    Int_t iLay = fNspd-1+i;
	    GetITSgeom()->CreateMatrix(startMod,iLay,iLad,iDet,kSDD,trans,rot);
	    startMod++;
	};
    };
    
    // SSD
    for (Int_t i=0; i<fNssd;i++) {
	sprintf(layerName, "ITSssdWafer_%i",i+fNspd+fNsdd+1);
	TGeoNode *layNode = fITSmotherVolume->GetNode(layerName);
	if (layNode) {
	    TGeoHMatrix layMatrix(*layNode->GetMatrix());	    
	    Double_t *trans  = layMatrix.GetTranslation();
	    Double_t *r      = layMatrix.GetRotationMatrix();
	    Double_t rot[10] = {r[0],r[1],r[2],
				r[3],r[4],r[5],
				r[6],r[7],r[8], 1.0};
	    Int_t iDet = 1;
	    if ((i+1==2)||(i+1==4)) iDet = 2;
	    Int_t iLad = 1;
	    Int_t iLay = 5;
	    if (i+1>2) iLay = 6;
	    GetITSgeom()->CreateMatrix(startMod,iLay,iLad,iDet,kSSD,trans,rot);
	    startMod++;
	};
    };
    
    return;
}

//______________________________________________________________________
void AliITSvBeamTestITS04::SetDefaults()
{
    // (from AliITSv11) mln
    
    const Float_t kconv = 1.0e+04; // convert cm to microns
    
    if(!fDetTypeSim) fDetTypeSim = new AliITSDetTypeSim();
    fDetTypeSim->SetITSgeom(GetITSgeom());
    fDetTypeSim->ResetCalibrationArray();
    fDetTypeSim->ResetSegmentation();
 
    AliITSgeomSPD *s0;
    AliITSgeomSDD *s1;
    AliITSgeomSSD *s2;
    Int_t i;
    Float_t bx[256],bz[280];

    // If fGeomDetIn is set true the geometry will
    // be initialised from file (JC)
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    AliITSgeom* geom = new AliITSgeom();
    SetITSgeom(geom);
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);

   
    // SPD

    s0 = (AliITSgeomSPD*) GetITSgeom()->GetShape(kSPD);// Get shape info.
    if (s0) {
	AliITSCalibration *resp0=new AliITSCalibrationSPD();
	SetCalibrationModel(kSPD,resp0);

	AliITSsegmentationSPD *seg0=new AliITSsegmentationSPD();
	seg0->SetDetSize(s0->GetDx()*2.*kconv, // base this on AliITSgeomSPD
			 s0->GetDz()*2.*kconv, // for now.
			 s0->GetDy()*2.*kconv);// x,z,y full width in microns.
	seg0->SetNPads(256,160);               // Number of Bins in x and z
	for(i=000;i<256;i++) bx[i] =  50.0;    // in x all are 50 microns.
	for(i=000;i<160;i++) bz[i] = 425.0;    // most are 425 microns except below
	for(i=160;i<280;i++) bz[i] =   0.0;    // Outside of detector.
	bz[ 31] = bz[ 32] = 625.0;             // first chip boundry
	bz[ 63] = bz[ 64] = 625.0;             // first chip boundry
	bz[ 95] = bz[ 96] = 625.0;             // first chip boundry
	bz[127] = bz[128] = 625.0;             // first chip boundry
	bz[160] = 425.0;                       // Set so that there is no zero pixel size for fNz.
	seg0->SetBinSize(bx,bz);               // Based on AliITSgeomSPD for now.
	SetSegmentationModel(kSPD,seg0);
	// set digit and raw cluster classes to be used
	const char *kData0=(fDetTypeSim->GetCalibrationModel(kSPD))->DataType();
	if (strstr(kData0,"real")) 
	  fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigit");
	else fDetTypeSim->SetDigitClassName(kSPD,"AliITSdigitSPD");
    }
  
    // SDD
   
    s1 = (AliITSgeomSDD*) GetITSgeom()->GetShape(kSDD);// Get shape info.
    if (s1) {
      AliITSCalibrationSDD *resp1=new AliITSCalibrationSDD("simulated");
      SetCalibrationModel(kSDD,resp1);
      AliITSsegmentationSDD *seg1=new AliITSsegmentationSDD();
      seg1->SetDetSize(s1->GetDx()*kconv, // base this on AliITSgeomSDD
		       s1->GetDz()*4.*kconv, // for now.
		       s1->GetDy()*4.*kconv); // x,z,y full width in microns.
      seg1->SetDriftSpeed(resp1->GetDriftSpeed());
      seg1->SetNPads(256,256);// Use AliITSgeomSDD for now
      SetSegmentationModel(kSDD,seg1);
      const char *kData1=(fDetTypeSim->GetCalibrationModel(kSDD))->DataType();
      const char *kopt=resp1->GetZeroSuppOption();
      if((!strstr(kopt,"2D")) && (!strstr(kopt,"1D")) || strstr(kData1,"real") ){
	fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigit");
	} else fDetTypeSim->SetDigitClassName(kSDD,"AliITSdigitSDD");
    }
    
    // SSD
    
    s2 = (AliITSgeomSSD*) GetITSgeom()->GetShape(kSSD);// Get shape info. Do it this way for now.
    if (s2) {
      AliITSCalibration *resp2=new AliITSCalibrationSSD("simulated");
      SetCalibrationModel(kSSD,resp2);

      AliITSsegmentationSSD *seg2=new AliITSsegmentationSSD();
      seg2->SetDetSize(s2->GetDx()*2.*kconv, // base this on AliITSgeomSSD
		       s2->GetDz()*2.*kconv, // for now.
		       s2->GetDy()*2.*kconv); // x,z,y full width in microns.
      seg2->SetPadSize(95.,0.); // strip x pitch in microns
      seg2->SetNPads(768,0); // number of strips on each side.
      seg2->SetAngles(0.0075,0.0275); // strip angels rad P and N side.
      seg2->SetAnglesLay5(0.0075,0.0275); // strip angels rad P and N side.
      seg2->SetAnglesLay6(0.0275,0.0075); // strip angels rad P and N side.
      SetSegmentationModel(kSSD,seg2); 
      const char *kData2=(fDetTypeSim->GetCalibrationModel(kSSD))->DataType();
      if(strstr(kData2,"real") ) fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigit");
      else fDetTypeSim->SetDigitClassName(kSSD,"AliITSdigitSSD");
    }
    
  if(fgkNTYPES>3){Warning("SetDefaults","Only the four basic detector types are initialised!");}
  return;
}

//______________________________________________________________________
void AliITSvBeamTestITS04::AddSPDGeometry(TGeoVolume *moth) const
{
  //Adds SPD geometry
    TGeoMedium *siliconSPD = gGeoManager->GetMedium("ITSspdSi");
    
    //outer volume
    TGeoBBox *waferSPDshape = new TGeoBBox("ITSspdWaferShape",fgkSPDwidth/2,fgkSPDthickness/2,fgkSPDlength/2,0);
    TGeoVolume *waferSPD = new TGeoVolume("ITSspdWafer",waferSPDshape,siliconSPD);
    //sensitive volume
    TGeoBBox *sensSPDbox = new TGeoBBox("ITSsddSensorSensBox",fgkSPDwidthSens/2,fgkSPDthickSens/2,fgkSPDlengthSens/2,0);
    TGeoVolume *sensVolSPD = new TGeoVolume(fgSPDsensitiveVolName,sensSPDbox,siliconSPD);
    waferSPD->AddNode(sensVolSPD, 1, 0); //added to outer volume
    
    //locate them in space (with respect top volume)
    TGeoTranslation *spd1tr = new TGeoTranslation(0,fgkSPD0y,fgkSPDlength/2);
    TGeoTranslation *spd2tr = new TGeoTranslation(0,fgkSPD0y,-fgkSPDlength/2);

    TGeoTranslation *spd3tr = new TGeoTranslation(0,fgkSPD1y,fgkSPDlength/2);
    TGeoTranslation *spd4tr = new TGeoTranslation(0,fgkSPD1y,-fgkSPDlength/2);
    
    //add to top volume
    moth->AddNode(waferSPD, 1, spd1tr);
    moth->AddNode(waferSPD, 2, spd2tr);
    moth->AddNode(waferSPD, 3, spd3tr);
    moth->AddNode(waferSPD, 4, spd4tr);
    
    //draw options
    waferSPD->SetLineColor(4);
    sensVolSPD->SetLineColor(4);
}


//______________________________________________________________________
void AliITSvBeamTestITS04::AddSDDGeometry(TGeoVolume *moth) const
{
  //Adds SDD geometry
    TGeoMedium *siliconSDD = gGeoManager->GetMedium("ITSsddSi");
    
    //outer volume
    TGeoBBox *waferSDDshape = new TGeoBBox("ITSsddWaferShape",fgkSDDwidth/2,fgkSDDthickness/2,fgkSDDlength/2,0);
    TGeoVolume *waferSDD = new TGeoVolume("ITSsddWafer",waferSDDshape,siliconSDD);
    //sensitive volume
    TGeoBBox *sensSDDbox = new TGeoBBox("ITSsddSensorSensBox",fgkSDDwidthSens/2,fgkSDDthickSens/2,fgkSDDlengthSens/2,0);
    TGeoVolume *sensVolSDD = new TGeoVolume(fgSDDsensitiveVolName,sensSDDbox,siliconSDD);
    waferSDD->AddNode(sensVolSDD, 1, 0); //added to outer volume
    
    //locate them in space
    TGeoTranslation *sdd1tr = new TGeoTranslation(0,fgkSDD0y,0);
    TGeoTranslation *sdd2tr = new TGeoTranslation(0,fgkSDD1y,0);
        
    //add to top volume
    moth->AddNode(waferSDD, fNspd+1, sdd1tr);
    moth->AddNode(waferSDD, fNspd+2, sdd2tr);
    
    //draw options
    waferSDD->SetLineColor(3);
    sensVolSDD->SetLineColor(3);
}


//______________________________________________________________________
void AliITSvBeamTestITS04::AddSSDGeometry(TGeoVolume *moth) const
{
  //Adds SSD geometry
    TGeoMedium *siliconSSD = gGeoManager->GetMedium("ITSspdSi"); // SSD medium still needed!!!
    
    //outer volume 
    TGeoBBox *waferSSDshape = new TGeoBBox("ITSssdWaferShape",fgkSSDwidth/2,fgkSSDthickness/2,fgkSSDlength/2,0);
    TGeoVolume *waferSSD = new TGeoVolume("ITSssdWafer",waferSSDshape,siliconSSD);
    //sensitive volume
    TGeoBBox *sensSSDbox = new TGeoBBox("ITSssdSensorSensBox",fgkSSDwidthSens/2,fgkSSDthickSens/2,fgkSSDlengthSens/2,0);
    TGeoVolume *sensVolSSD = new TGeoVolume(fgSSDsensitiveVolName,sensSSDbox,siliconSSD);
    waferSSD->AddNode(sensVolSSD, 1, 0);
    
    //locate them in space
    /* In the SSD, there was an overlap of sensitive volumes of 2.9mm = 0.29cm (0.29/2=0.145) 
       in the modules in the same plane, therefore the modules where not in the same plane in 
       the Y direction, there was a "thickness" (0.03cm) difference */
    TGeoTranslation *ssd1tr = new TGeoTranslation(0,fgkSSD0y,fgkSSDlength/2-0.145);
    TGeoTranslation *ssd2tr = new TGeoTranslation(0,fgkSSD0y+0.03,-fgkSSDlength/2+0.145);

    TGeoTranslation *ssd3tr = new TGeoTranslation(0,fgkSSD1y,fgkSSDlength/2-0.145);
    TGeoTranslation *ssd4tr = new TGeoTranslation(0,fgkSSD1y+0.03,-fgkSSDlength/2+0.145);

    //add to top volume
    moth->AddNode(waferSSD, fNspd+fNsdd+1, ssd1tr);
    moth->AddNode(waferSSD, fNspd+fNsdd+2, ssd2tr);
    moth->AddNode(waferSSD, fNspd+fNsdd+3, ssd3tr);
    moth->AddNode(waferSSD, fNspd+fNsdd+4, ssd4tr);
    
    //draw options
    waferSSD->SetLineColor(2);
    sensVolSSD->SetLineColor(2);
}

//______________________________________________________________________
void AliITSvBeamTestITS04::StepManager()
{
    // Called for every step in the ITS, then calles the AliITShit class
    // creator with the information to be recoreded about that hit.

    // "Standard" StepManager. (Similar to AliITSv11) mln
    Int_t copy, id;
    TLorentzVector position, momentum;
    static TLorentzVector position0;
    static Int_t stat0=0;
    
    if(!(this->IsActive())){
	return;
    } // end if !Active volume.
    
    if(!(gMC->TrackCharge())) return;
    
    id=gMC->CurrentVolID(copy);
    
    Bool_t sensvol = kFALSE;
    for(Int_t kk = 0; kk < fIdN; kk++)
	if(id == fIdSens[kk]) sensvol = kTRUE;
    
    if (sensvol && (gMC->IsTrackExiting())) {
	copy = fTrackReferences->GetEntriesFast();
	TClonesArray &lTR = *fTrackReferences;
	// Fill TrackReference structure with this new TrackReference.
	new(lTR[copy]) AliTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
    } // if Outer ITS mother Volume
    
    Int_t   vol[5];
    TClonesArray &lhits = *fHits;
    //
    // Track status
    vol[3] = 0;
    vol[4] = 0;
    // Fill hit structure.
    if(gMC->IsTrackInside())      vol[3] +=  1;
    if(gMC->IsTrackEntering())    vol[3] +=  2;
    if(gMC->IsTrackExiting())     vol[3] +=  4;
    if(gMC->IsTrackOut())         vol[3] +=  8;
    if(gMC->IsTrackDisappeared()) vol[3] += 16;
    if(gMC->IsTrackStop())        vol[3] += 32;
    if(gMC->IsTrackAlive())       vol[3] += 64;
    
    // Only entering charged tracks
    if(!(gMC->TrackCharge())) return;
    
    if( ((id = gMC->CurrentVolID(copy)) == fIdSens[0]) ||
	((id = gMC->CurrentVolID(copy)) == fIdSens[1]) ||
	((id = gMC->CurrentVolID(copy)) == fIdSens[2]) )
    {
	GetCurrentLayLaddDet(vol[0], vol[2], vol[1]);
	// vol[2], vol[1]) : in this order because the ladder
	// index and the det. index are exchanged in the constructor
	// of AliITShit...
    } else {
	return; // not an ITS volume?
    };
    
    gMC->TrackPosition(position);
    gMC->TrackMomentum(momentum);
    vol[4] = stat0;
    if(gMC->IsTrackEntering()){
	position0 = position;
	stat0 = vol[3];
	return;
    } // end if IsEntering
    // Fill hit structure with this new hit.
    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),
				   vol, gMC->Edep(),gMC->TrackTime(),position,
				   position0,momentum);
    //
    position0 = position;
    stat0 = vol[3];
    return;
}

//______________________________________________________________________
Int_t AliITSvBeamTestITS04::GetCurrentLayLaddDet(Int_t &lay,Int_t &ladd, Int_t &det) const
{ 
  // Function which gives the layer, ladder and det.
  // index of the current volume. To be used in
  // AliITS::StepManager()

    det  = 1;   ladd = 1;
    
    TGeoNode *node = gGeoManager->GetMother(1);
    if (!node) return kFALSE;
    Int_t nodeNum = node->GetNumber();
    
    // GetNumber() return the index recorded in the node
    
    if (nodeNum==5||nodeNum==6) {         // SDD: one layer, one detector
	lay = nodeNum-2;
    } else if (nodeNum==3||nodeNum==4) {  // SPD layer 2
	lay = 2;
	if (nodeNum==4) det = 2;
    } else if (nodeNum==1||nodeNum==2){   // SPD layer 1
	lay = 1;
	if (nodeNum==2) det = 2; 
    } else if (nodeNum==9||nodeNum==10) { // SSD layer 2
	lay = 6;
	if (nodeNum==10) det = 2;
    } else if (nodeNum==7||nodeNum==8){   // SSD layer 1
	lay = 5;
	if (nodeNum==8) det = 2; 
    };  
    
    return kTRUE;
}

//_____________________________________________________________

 Int_t AliITSvBeamTestITS04::GetNumberOfSubDet(const TString& det) const{
    
   //Get number of individual detectors
    if(det.Contains("SPD")) return fNspd;
    if(det.Contains("SDD")) return fNsdd;
    if(det.Contains("SSD")) return fNssd;
    return 0;
  }
