/// \file RegisterCorrection.C
///
/// \author marian.ivanov@cern.ch
///
/// Register primitive corrections: base functions for minimization:
/// Id numbers are associated to given primitive corrections.
/// See comments in the function headers. 
/// Used only for residuals minimization not in the reconstruction.
/// File with all primitives expected to be in the current directory:
/// filenames:  
/// TPCCorrectionPrimitives.root                  - Alignment, Quadrants, 2D symentrid and rod misalingment
/// TPCCorrectionPrimitivesROC.root               - ExB distortion due to the common ROC misalignment
/// TPCCorrectionPrimitivesSector.root            - ExB distortion due to the one sector misalingment (sector 0)
/// TPCCorrectionPrimitivesFieldCage.root         - ExB distortion due to the Field cage missalignment
///
///
/// RegisterCorrection();                         - Reserved id's  0 -999
///
/// RegisterAliTPCFCVoltError3D();                - Reserved id's  0 -99
/// RegisterAliTPCBoundaryVoltError();            - Reserved id's  100-199
/// RegisterAliTPCCalibGlobalMisalignment();      - Reserved id's  200-499
/// RegisterAliTPCExBBShape();                    - Reserved id's  500-600
/// RegisterAliTPCExBTwist();                     - Reserved id's  600-700
/// RegisterAliTPCROCVoltError3D()                - Reserved id's  700-800
/// RegisterAliTPCROCVoltError3DSector()          - Reserved id's  800-900
///
/// ~~~
/// .x ~/rootlogon.C
/// .L $ALICE_ROOT/TPC/CalibMacros/RegisterCorrection.C+
/// RegisterCorrection();
/// ~~~
///
/// Example use: 
///
/// ~~~
/// .x ~/rootlogon.C
/// .L $ALICE_ROOT/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C
/// ConfigCalibTrain(119037,"local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB/")
///
/// .L $ALICE_ROOT/TPC/CalibMacros/RegisterCorrection.C+
/// RegisterCorrection(0);
/// ~~~
///
/// See example usage of correction primitive/derivatives in file

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TObjArray.h"
#include "TTreeStream.h"
#include "TMath.h" 
#include "TGraph.h" 
#include "TRandom.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TAxis.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliTPCParamSR.h"
#include "TDatabasePDG.h"
#include "AliTPCCorrection.h"
#include "AliTPCFCVoltError3D.h"
#include "AliTPCROCVoltError3D.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCBoundaryVoltError.h"
#include "AliTPCCalibGlobalMisalignment.h"
#include "AliTPCExBBShape.h"
#include "AliTPCExBTwist.h"
#include "AliTPCCorrectionDrift.h"
#include "AliMagF.h"
#include "AliCDBEntry.h"
#include "AliTPCROC.h"
#include <TStatToolkit.h>
#include "TCut.h"
#include "TGraphErrors.h"
#include  "AliTrackerBase.h"
#include  "TGeoGlobalMagField.h"
#include  "TROOT.h"
#include "TLinearFitter.h"
#include "TStopwatch.h"
#include "AliTPCCorrectionFit.h"
#endif




TFile *fCorrections=0;                         //file with corrections
//
//models E field distortion   AliTPCFCVoltError3D
//
AliTPCFCVoltError3D *rotOFC  =0;     // fit models
AliTPCFCVoltError3D *rodOFC1 =0;
AliTPCFCVoltError3D *rodOFC2 =0;
AliTPCFCVoltError3D *rotIFC  =0;
AliTPCFCVoltError3D *rodIFC1 =0;
AliTPCFCVoltError3D *rodIFC2 =0;
//
AliTPCFCVoltError3D *rodIFCShift=0;     //common IFC shift
AliTPCFCVoltError3D *rodIFCSin=0;       //common IFC shift -sinus
AliTPCFCVoltError3D *rodOFCSin=0;       //common OFC shift -sinus
AliTPCFCVoltError3D *rodOFCShift=0;     //common OFC shift
AliTPCFCVoltError3D *rodIFCCos=0;       //common IFC shift -cos
AliTPCFCVoltError3D *rodOFCCos=0;       //common OFC shift -cos
//
TObjArray  *rodFCSideRadiusType=0;  // field cage shift, A/C side, OFC/IFC, Fourier(const,sin,cos) 
//



AliTPCROCVoltError3D *rocRotgXA=0;     // roc rotation A side - inclination in X
AliTPCROCVoltError3D *rocRotgYA=0;     // roc rotation A side - inclination in Y
AliTPCROCVoltError3D *rocRotgXC=0;     // roc rotation C side - inclination in X
AliTPCROCVoltError3D *rocRotgYC=0;     // roc rotation C side - inclination in Y
//
AliTPCROCVoltError3D *rocDzIROCA=0;      // roc shift A side - in Z
AliTPCROCVoltError3D *rocDzIROCC=0;      // roc shift C side - in Z
AliTPCROCVoltError3D *rocRotIROCA=0;      // roc rotation A side - in Z
AliTPCROCVoltError3D *rocRotIROCC=0;      // roc rotation C side - in Z
AliTPCROCVoltError3D *rocDzUDA=0;         // roc shift Up-Down A side
AliTPCROCVoltError3D *rocDzUDC=0;         // roc shift Up-Down C side
AliTPCROCVoltError3D *rocRotUDA=0;        // roc rotation updown A side
AliTPCROCVoltError3D *rocRotUDC=0;        // roc rotation updown C side

//
AliTPCROCVoltError3D *rocShiftIROCA0=0;      // IROC shift A0 side
AliTPCROCVoltError3D *rocRotIROCA0=0;        // IROC rot   A0 side
AliTPCROCVoltError3D *rocShiftOROCA0=0;      // OROC shift A0 side
AliTPCROCVoltError3D *rocRotOROCA0=0;        // OROC rot   A0 side
AliTPCROCVoltError3D *rocShiftIROCC0=0;      // IROC shift C0 side
AliTPCROCVoltError3D *rocRotIROCC0=0;        // IROC rot   C0 side
AliTPCROCVoltError3D *rocShiftOROCC0=0;      // OROC shift C0 side
AliTPCROCVoltError3D *rocRotOROCC0=0;        // OROC rot   C0 side
//
//
//
AliTPCBoundaryVoltError *boundaryVoltErrorA[8];  // boundary voltage error A side
AliTPCBoundaryVoltError *boundaryVoltErrorC[8];  // boundary voltage error C side
AliTPCExBBShape *exbShape     = 0;               // nominal correctin
AliTPCExBBShape *exbShapeT1X  = 0;               // nominal +deltaT1=0.1
AliTPCExBBShape *exbShapeT2X  = 0;               // nominal +deltaT2=0.1
AliTPCExBTwist  *twistX    = 0;
AliTPCExBTwist  *twistY    = 0;
//
AliTPCCalibGlobalMisalignment *alignRot0=0;
AliTPCCalibGlobalMisalignment *alignRot1=0;
AliTPCCalibGlobalMisalignment *alignRot2=0;
AliTPCCalibGlobalMisalignment *alignTrans0=0;
AliTPCCalibGlobalMisalignment *alignTrans1=0;
AliTPCCalibGlobalMisalignment *alignTrans2=0;
//
AliTPCCalibGlobalMisalignment *alignTrans0D[4]={0};       // delta alignemnt 4 quadrants
AliTPCCalibGlobalMisalignment *alignTrans1D[4]={0};       // 
AliTPCCalibGlobalMisalignment *alignTrans2D[4]={0};       // 
AliTPCCalibGlobalMisalignment *alignRot0D[4]={0};
AliTPCCalibGlobalMisalignment *alignRot1D[4]={0};
AliTPCCalibGlobalMisalignment *alignRot2D[4]={0};




AliTPCCorrectionDrift  *calibDrift[7]={0};


//
void RegisterAliTPCFCVoltError3D();
void RegisterAliTPCCalibGlobalMisalignment();
void RegisterAliTPCBoundaryVoltError();
void RegisterAliTPCExBShape();
void RegisterAliTPCExBTwist();
void RegisterAliTPCROCVoltError3D();
void RegisterAliTPCROCVoltError3DSector();
void RegisterAliTPCCorrectionDrift();
void RegisterAliTPCFCVoltError3DRodFCSideRadiusType();
//

void RegisterCorrection(Int_t type=0){
  /// check the presence of corrections in file
  ///
  /// gROOT->Macro("ConfigCalibTrain.C(119037)");

  if (type==1) return RegisterAliTPCROCVoltError3D();            // 3D distortion due misalignemnt of FC ....
  if (type==2) return RegisterAliTPCROCVoltError3DSector();      // 3D distortion due misalingment of Sectors
  if (type==4) return RegisterAliTPCFCVoltError3DRodFCSideRadiusType(); 
  //
  if (type==3) {                                                 // 2D distortions
    fCorrections = TFile::Open("TPCCorrectionPrimitives.root","recreate");
    RegisterAliTPCCalibGlobalMisalignment();
    RegisterAliTPCBoundaryVoltError();
    RegisterAliTPCFCVoltError3D();
    RegisterAliTPCExBShape();
    RegisterAliTPCExBTwist();
    RegisterAliTPCCorrectionDrift();
    fCorrections->Close();
    delete fCorrections;
    return;
  }
  fCorrections = TFile::Open("TPCCorrectionPrimitives.root");
  AliTPCComposedCorrection *corrField3D = (AliTPCComposedCorrection*) fCorrections->Get("TPCFCVoltError3D");
  // if not make new file
  if (!corrField3D) fCorrections = TFile::Open("TPCCorrectionPrimitives.root","update");
  if (type==0) {
    RegisterAliTPCROCVoltError3D();
    RegisterAliTPCROCVoltError3DSector();
    RegisterAliTPCFCVoltError3DRodFCSideRadiusType(); 
  }
  RegisterAliTPCCalibGlobalMisalignment();
  RegisterAliTPCBoundaryVoltError();
  RegisterAliTPCFCVoltError3D();
  RegisterAliTPCExBShape();
  RegisterAliTPCExBTwist();
  RegisterAliTPCCorrectionDrift();
  if (fCorrections) fCorrections->Close();
  //
}



void RegisterAliTPCFCVoltError3D(){
  /// Load the models from the file
  /// Or create it
  /// Register functions with following IDs:
  /// IMPORTANT: The nominal shift is in mm
  ///
  ///  rotOFC  - 0
  ///  rodOFC1 - 1
  ///  rodOFC2 - 2
  ///  rotIFC  - 3
  ///  rodIFC1 - 4
  ///  rodIFC2 - 5
  ///  rodIFCShift - 6
  ///  rodIFCSin   - 7
  ///  rodIFCCos   - 8
  ///  rodOFCShift - 9
  ///  rodOFCSin   - 10
  ///  rodOFCCos   - 11

  printf("RegisterAliTPCFCVoltError3D()");
  Int_t volt = 40; // 40 V ~  1mm
  AliTPCComposedCorrection *corrField3D = (AliTPCComposedCorrection*) fCorrections->Get("TPCFCVoltError3D");    
  if (corrField3D) { // load form file
    corrField3D->Print();
    TCollection *iter = corrField3D->GetCorrections();    
    rotOFC = (AliTPCFCVoltError3D*)iter->FindObject("rotOFC");
    rodOFC1 = (AliTPCFCVoltError3D*)iter->FindObject("rodOFC1");
    rodOFC2 = (AliTPCFCVoltError3D*)iter->FindObject("rodOFC2");
    rotIFC = (AliTPCFCVoltError3D*)iter->FindObject("rotIFC");
    rodIFC1 = (AliTPCFCVoltError3D*)iter->FindObject("rodIFC1");
    rodIFC2 = (AliTPCFCVoltError3D*)iter->FindObject("rodIFC2");
    //
    rodIFCShift = (AliTPCFCVoltError3D*)iter->FindObject("rodIFCShift");
    rodOFCShift = (AliTPCFCVoltError3D*)iter->FindObject("rodOFCShift");
    rodIFCSin = (AliTPCFCVoltError3D*)iter->FindObject("rodIFCSin");
    rodOFCSin = (AliTPCFCVoltError3D*)iter->FindObject("rodOFCSin");
    rodIFCCos = (AliTPCFCVoltError3D*)iter->FindObject("rodIFCCos");
    rodOFCCos = (AliTPCFCVoltError3D*)iter->FindObject("rodOFCCos");
  } else {    
    // OFC 
    rotOFC = new AliTPCFCVoltError3D();
    rotOFC->SetOmegaTauT1T2(0,1,1);
    rotOFC->SetRotatedClipVoltA(1,volt);
    rotOFC->SetRotatedClipVoltC(1,volt);
    //
    rodOFC1 = new AliTPCFCVoltError3D();
    rodOFC1->SetOmegaTauT1T2(0,1,1);
    rodOFC1->SetRodVoltShiftA(18,volt);
    rodOFC1->SetRodVoltShiftC(18,volt);
    //
    rodOFC2 = new AliTPCFCVoltError3D();
    rodOFC2->SetOmegaTauT1T2(0,1,1);
    rodOFC2->SetCopperRodShiftA(18,volt);
    rodOFC2->SetCopperRodShiftC(18,volt);    
    // IFC     
    rotIFC = new AliTPCFCVoltError3D();
    rotIFC->SetOmegaTauT1T2(0,1,1);
    rotIFC->SetRotatedClipVoltA(0,volt);
    rotIFC->SetRotatedClipVoltC(0,volt);
    //
    rodIFC1 = new AliTPCFCVoltError3D();
    rodIFC1->SetOmegaTauT1T2(0,1,1);
    rodIFC1->SetRodVoltShiftA(0,volt);
    rodIFC1->SetRodVoltShiftC(0,volt);
    //
    rodIFC2 = new AliTPCFCVoltError3D();
    rodIFC2->SetOmegaTauT1T2(0,1,1);
    rodIFC2->SetCopperRodShiftA(0,volt);
    rodIFC2->SetCopperRodShiftC(0,volt);
    //
    rodIFCShift = new AliTPCFCVoltError3D();
    rodIFCSin = new AliTPCFCVoltError3D();
    rodIFCCos = new AliTPCFCVoltError3D();
    rodOFCShift = new AliTPCFCVoltError3D();
    rodOFCSin = new AliTPCFCVoltError3D();
    rodOFCCos = new AliTPCFCVoltError3D();
    for (Int_t isec=0; isec<18; isec++){
      Double_t phi=TMath::Pi()*isec/9.;
      rodIFCShift->SetOmegaTauT1T2(0,1,1);
      rodIFCShift->SetRodVoltShiftA(isec,volt);
      rodIFCShift->SetRodVoltShiftC(isec,volt);
      rodIFCSin->SetOmegaTauT1T2(0,1,1);
      rodIFCSin->SetRodVoltShiftA(isec,volt*TMath::Sin(phi));
      rodIFCSin->SetRodVoltShiftC(isec,volt*TMath::Sin(phi));
      rodIFCCos->SetOmegaTauT1T2(0,1,1);
      rodIFCCos->SetRodVoltShiftA(isec,volt*TMath::Cos(phi));
      rodIFCCos->SetRodVoltShiftC(isec,volt*TMath::Cos(phi));
      //
      rodOFCShift->SetOmegaTauT1T2(0,1,1);
      rodOFCShift->SetRodVoltShiftA(18+isec,volt);
      rodOFCShift->SetRodVoltShiftC(18+isec,volt);
      rodOFCSin->SetOmegaTauT1T2(0,1,1);
      rodOFCSin->SetRodVoltShiftA(18+isec,volt*TMath::Sin(phi));
      rodOFCSin->SetRodVoltShiftC(18+isec,volt*TMath::Sin(phi));
      rodOFCCos->SetOmegaTauT1T2(0,1,1);
      rodOFCCos->SetRodVoltShiftA(18+isec,volt*TMath::Cos(phi));
      rodOFCCos->SetRodVoltShiftC(18+isec,volt*TMath::Cos(phi));
    }
    //
    //
    // Initialization of the lookup tables
    //
    printf(" ------- OFC rotated clip:\n"); rotOFC->InitFCVoltError3D();
    printf(" ------- OFC rod & strip:\n");  rodOFC1->InitFCVoltError3D();
    printf(" ------- OFC copper rod:\n");   rodOFC2->InitFCVoltError3D();
    printf(" ------- IFC rotated clip:\n"); rotIFC->InitFCVoltError3D();
    printf(" ------- IFC rod & strip:\n");  rodIFC1->InitFCVoltError3D();
    printf(" ------- IFC copper rod:\n");   rodIFC2->InitFCVoltError3D();

    printf(" ------- IFC rod & strip shift:\n");  rodIFCShift->InitFCVoltError3D();
    printf(" ------- IFC rod & strip sin:\n");    rodIFCSin->InitFCVoltError3D();
    printf(" ------- IFC rod & strip cos:\n");    rodIFCCos->InitFCVoltError3D();
    printf(" ------- OFC rod & strip shift:\n");  rodOFCShift->InitFCVoltError3D();
    printf(" ------- OFC rod & strip sin:\n");    rodOFCSin->InitFCVoltError3D();
    printf(" ------- OFC rod & strip cos:\n");    rodOFCCos->InitFCVoltError3D();

    // give names
    rotOFC->SetName("rotOFC");rotOFC->SetTitle("rotOFC");
    rodOFC1->SetName("rodOFC1");rodOFC1->SetTitle("rodOFC1");
    rodOFC2->SetName("rodOFC2");rodOFC2->SetTitle("rodOFC2");
    rotIFC->SetName("rotIFC");rotIFC->SetTitle("rotIFC");
    rodIFC1->SetName("rodIFC1");rodIFC1->SetTitle("rodIFC1");
    rodIFC2->SetName("rodIFC2");rodIFC2->SetTitle("rodIFC2");
    //
    rodIFCShift->SetName("rodIFCShift");rodIFCShift->SetTitle("rodIFCShift");
    rodIFCSin->SetName("rodIFCSin");rodIFCSin->SetTitle("rodIFCSin");
    rodIFCCos->SetName("rodIFCCos");rodIFCCos->SetTitle("rodIFCCos");
    //
    rodOFCShift->SetName("rodOFCShift");rodOFCShift->SetTitle("rodOFCShift");
    rodOFCSin->SetName("rodOFCSin");rodOFCSin->SetTitle("rodOFCSin");
    rodOFCCos->SetName("rodOFCCos");rodOFCCos->SetTitle("rodOFCCos");
    //
    // save in file
    corrField3D = new AliTPCComposedCorrection();
    TObjArray *cs = new TObjArray();
    cs->Add(rotIFC); cs->Add(rotOFC);
    cs->Add(rodIFC1); cs->Add(rodOFC1);
    cs->Add(rodIFC2); cs->Add(rodOFC2);
    cs->Add(rodIFCShift);    cs->Add(rodIFCSin);    cs->Add(rodIFCCos);
    cs->Add(rodOFCShift);    cs->Add(rodOFCSin);    cs->Add(rodOFCCos);
    //
    corrField3D->SetCorrections(cs);
    corrField3D->SetOmegaTauT1T2(0,1.,1.);
    corrField3D->Print();    
    fCorrections->cd();
    corrField3D->Write("TPCFCVoltError3D");
  }
  //
  AliTPCCorrection::AddVisualCorrection(rotOFC,0); 
  AliTPCCorrection::AddVisualCorrection(rodOFC1,1); 
  AliTPCCorrection::AddVisualCorrection(rodOFC2,2); 
  AliTPCCorrection::AddVisualCorrection(rotIFC,3); 
  AliTPCCorrection::AddVisualCorrection(rodIFC1,4); 
  AliTPCCorrection::AddVisualCorrection(rodIFC2,5); 
  // common corrections
  //
  AliTPCCorrection::AddVisualCorrection(rodIFCShift,6); 
  AliTPCCorrection::AddVisualCorrection(rodIFCSin,7); 
  AliTPCCorrection::AddVisualCorrection(rodIFCCos,8); 
  //
  AliTPCCorrection::AddVisualCorrection(rodOFCShift,9); 
  AliTPCCorrection::AddVisualCorrection(rodOFCSin,10); 
  AliTPCCorrection::AddVisualCorrection(rodOFCCos,11); 
}


void RegisterAliTPCFCVoltError3DRodFCSideRadiusType(){
  /// Load the models from the file
  /// Or create it
  /// Register functions with following IDs:
  /// IMPORTANT: The nominal shift is in mm
  ///
  /// naming convention:
  /// rodFCSide%dRadius%dType%d

  ::Info("RegisterAliTPCFCVoltError3DRodFCSideRadiusType()","Start");
  Int_t volt = 40; // 40 V ~  1mm
  TFile * fCorrectionsRodFCSideRadiusType = TFile::Open("TPCCorrectionPrimitivesFieldCage.root","update");
  if (!fCorrectionsRodFCSideRadiusType)  fCorrectionsRodFCSideRadiusType = TFile::Open("TPCCorrectionPrimitivesFieldCage.root","recreate");
  rodFCSideRadiusType= new TObjArray(100);
  //
  AliTPCComposedCorrection *corrFieldCage = (AliTPCComposedCorrection*) fCorrectionsRodFCSideRadiusType->Get("TPCRodFCSideRadiusType");    
  if (corrFieldCage) { // load from file
    TCollection *corrections = corrFieldCage->GetCorrections();    
    for (Int_t itype=0; itype<3; itype++){
      for (Int_t iside=0; iside<2; iside++){
	for (Int_t ifc=0; ifc<2; ifc++){
	  Int_t id=4*itype+2*iside+ifc;
	  AliTPCFCVoltError3D *corr =  (AliTPCFCVoltError3D *)corrections->FindObject(TString::Format("rodFCSide%dRadius%dType%d",iside,ifc,itype).Data());
	  rodFCSideRadiusType->AddLast(corr);
	  ::Info("Read AliTPCFCVoltError3DRodFCSideRadiusType", TString::Format("Init: rodFCSide%dRadius%dType%d\tId=%d",iside,ifc,itype,900+id).Data());
	  AliTPCCorrection::AddVisualCorrection(corr,900+id); 
	}
      }
    }    
    corrFieldCage->Print();
  } else {    
    for (Int_t itype=0; itype<3; itype++){
      for (Int_t iside=0; iside<2; iside++){
	for (Int_t ifc=0; ifc<2; ifc++){
	  Int_t id=4*itype+2*iside+ifc;
	  AliTPCFCVoltError3D *corr =  new AliTPCFCVoltError3D();
	  corr->SetName(TString::Format("rodFCSide%dRadius%dType%d",iside,ifc,itype).Data());
	  corr->SetTitle(TString::Format("rodFCSide%dRadius%dType%d",iside,ifc,itype).Data());
	  //
	  for (Int_t isec=0; isec<18; isec++){
	    Double_t phi=TMath::Pi()*isec/9.;
	    Int_t sectorOffset=(ifc==0) ? 0:18;
	    corr->SetOmegaTauT1T2(0,1,1);
	    Double_t sectorVoltage=0;
	    if (itype==0) sectorVoltage=volt;
	    if (itype==1) sectorVoltage=volt*TMath::Sin(phi);
	    if (itype==2) sectorVoltage=volt*TMath::Cos(phi);
	    if (iside==0){
	      corr->SetRodVoltShiftA(isec+sectorOffset,sectorVoltage);
	    }
	    if (iside==1){
	      corr->SetRodVoltShiftC(isec+sectorOffset,sectorVoltage);
	    }
	  }
	  rodFCSideRadiusType->AddLast(corr); 
	  ::Info("Generate AliTPCFCVoltError3DRodFCSideRadiusType", TString::Format("Init: rodFCSide%dRadius%dType%d",iside,ifc,itype).Data());
	  corr->InitFCVoltError3D();	  
	  AliTPCCorrection::AddVisualCorrection(corr,900+id); 
	}
      }
    }
    //
    corrFieldCage = new AliTPCComposedCorrection();
    corrFieldCage->SetCorrections(rodFCSideRadiusType);
    corrFieldCage->SetOmegaTauT1T2(0,1.,1.);
    
    corrFieldCage->Print();    
    fCorrectionsRodFCSideRadiusType->cd();
    corrFieldCage->Write("TPCRodFCSideRadiusType");    
  }
  fCorrectionsRodFCSideRadiusType->Close();
  ::Info("RegisterAliTPCFCVoltError3DRodFCSideRadiusType()","End");
}


void RegisterAliTPCCalibGlobalMisalignment(){
  /// Register primitive alignment components.
  /// Linear conbination of primitev forulas used for fit
  /// The nominal delta 1 mm in shift and 1 mrad in rotation
  /// Primitive formulas registeren in AliTPCCoreection::AddvisualCorrection
  /// 20 - deltaX
  /// 21 - deltaY
  /// 22 - deltaZ
  /// 23 - rot0 (phi)
  /// 24 - rot1 (theta)
  /// 25 - rot2

  printf("RegisterAliTPCCalibGlobalMisalignment()\n");
  TGeoHMatrix matrixX;
  TGeoHMatrix matrixY;
  TGeoHMatrix matrixZ;
  TGeoRotation rot0;
  TGeoRotation rot1;
  TGeoRotation rot2;  //transformation matrices
  TGeoRotation rot90;  //transformation matrices
  matrixX.SetDx(0.1); matrixY.SetDy(0.1); matrixZ.SetDz(0.1); //1 mm translation
  rot0.SetAngles(0.001*TMath::RadToDeg(),0,0);
  rot1.SetAngles(0,0.001*TMath::RadToDeg(),0);
  rot2.SetAngles(0,0,0.001*TMath::RadToDeg());
  //how to get rot02 ?
  rot90.SetAngles(0,90,0);
  rot2.MultiplyBy(&rot90,kTRUE);
  rot90.SetAngles(0,-90,0);
  rot2.MultiplyBy(&rot90,kFALSE);
  //
  alignRot0  =new  AliTPCCalibGlobalMisalignment;
  alignRot0->SetAlignGlobal(&rot0);
  alignRot0->SetName("alignRot0");
  alignRot1=new  AliTPCCalibGlobalMisalignment;
  alignRot1->SetAlignGlobal(&rot1);
  alignRot1->SetName("alignRot1");
  alignRot2=new  AliTPCCalibGlobalMisalignment;
  alignRot2->SetAlignGlobal(&rot2);
  alignRot2->SetName("alignRot2");
  //
  alignTrans0  =new  AliTPCCalibGlobalMisalignment;
  alignTrans0->SetAlignGlobal(&matrixX);
  alignTrans0->SetName("alignTrans0");
  alignTrans1=new  AliTPCCalibGlobalMisalignment;
  alignTrans1->SetAlignGlobal(&matrixY);
  alignTrans1->SetName("alignTrans1");
  alignTrans2=new  AliTPCCalibGlobalMisalignment;
  alignTrans2->SetAlignGlobal(&matrixZ);
  alignTrans2->SetName("alignTrans2");
  //
  

  alignRot0D[0]  =new  AliTPCCalibGlobalMisalignment;
  alignRot0D[0]->SetAlignGlobalDelta(&rot0);
  alignRot0D[0]->SetName("alignRot0D");
  alignRot1D[0]=new  AliTPCCalibGlobalMisalignment;
  alignRot1D[0]->SetAlignGlobalDelta(&rot1);
  alignRot1D[0]->SetName("alignRot1D");
  alignRot2D[0]=new  AliTPCCalibGlobalMisalignment;
  alignRot2D[0]->SetAlignGlobalDelta(&rot2);
  alignRot2D[0]->SetName("alignRot2D");
  //
  alignTrans0D[0]  =new  AliTPCCalibGlobalMisalignment;
  alignTrans0D[0]->SetAlignGlobalDelta(&matrixX);
  alignTrans0D[0]->SetName("alignTrans0D");
  alignTrans1D[0]=new  AliTPCCalibGlobalMisalignment;
  alignTrans1D[0]->SetAlignGlobalDelta(&matrixY);
  alignTrans1D[0]->SetName("alignTrans1D");
  alignTrans2D[0]=new  AliTPCCalibGlobalMisalignment;
  alignTrans2D[0]->SetAlignGlobalDelta(&matrixZ);
  alignTrans2D[0]->SetName("alignTrans2D");

  TObjArray * arrayDX = new TObjArray(72);
  TObjArray * arrayDY = new TObjArray(72);
  TObjArray * arrayDPhi = new TObjArray(72);
  //
  // Up down A side
  //
  for (Int_t isec=0; isec<72; isec++){  //A side
    TGeoHMatrix *matrixDX = new TGeoHMatrix;
    TGeoHMatrix *matrixDY = new TGeoHMatrix;
    TGeoRotation *matrixDPhi= new TGeoRotation;
    arrayDX->AddAt(matrixDX,isec);
    arrayDY->AddAt(matrixDY,isec);
    arrayDPhi->AddAt(matrixDPhi,isec);
    if (isec%36<18) matrixDX->SetDx(isec%18<9?0.1:-0.1);
    if (isec%36<18) matrixDY->SetDy(isec%18<9?0.1:-0.1);
    if (isec%36<18) matrixDPhi->SetAngles((isec%18<9?0.001:-0.001)*TMath::RadToDeg(),0,0);
  }
  alignTrans0D[1]  =new  AliTPCCalibGlobalMisalignment;
  alignTrans0D[1]->SetName("alignTrans0UDA");
  alignTrans0D[1]->SetAlignSectors((TObjArray *)(arrayDX->Clone()));
  alignTrans1D[1]  =new  AliTPCCalibGlobalMisalignment;
  alignTrans1D[1]->SetName("alignTrans1UDA");
  alignTrans1D[1]->SetAlignSectors((TObjArray *)(arrayDY->Clone()));
  alignRot0D[1]  =new  AliTPCCalibGlobalMisalignment;
  alignRot0D[1]->SetName("alignRot0UDA");
  alignRot0D[1]->SetAlignSectors((TObjArray *)(arrayDPhi->Clone()));
  //
  // Uu-down C side
  //
  for (Int_t isec=0; isec<72; isec++){  //A side
    TGeoHMatrix *matrixDX = new TGeoHMatrix;
    TGeoHMatrix *matrixDY = new TGeoHMatrix;
    TGeoRotation *matrixDPhi= new TGeoRotation;
    arrayDX->AddAt(matrixDX,isec);
    arrayDY->AddAt(matrixDY,isec);
    arrayDPhi->AddAt(matrixDPhi,isec);
    if (isec%36>=18) matrixDX->SetDx(isec%18<9?0.1:-0.1);
    if (isec%36>=18) matrixDY->SetDy(isec%18<9?0.1:-0.1);
    if (isec%36>=18) matrixDPhi->SetAngles((isec%18<9?0.001:-0.001)*TMath::RadToDeg(),0,0);
  }
  alignTrans0D[2]  =new  AliTPCCalibGlobalMisalignment;
  alignTrans0D[2]->SetName("alignTrans0UDC");
  alignTrans0D[2]->SetAlignSectors((TObjArray *)(arrayDX->Clone()));
  alignTrans1D[2]  =new  AliTPCCalibGlobalMisalignment;
  alignTrans1D[2]->SetName("alignTrans1UDC");
  alignTrans1D[2]->SetAlignSectors((TObjArray *)(arrayDY->Clone()));
  alignRot0D[2]  =new  AliTPCCalibGlobalMisalignment;
  alignRot0D[2]->SetName("alignRot0UDC");
  alignRot0D[2]->SetAlignSectors((TObjArray *)(arrayDPhi->Clone()));
  //


  //
  AliTPCCorrection::AddVisualCorrection(alignTrans0  ,200);
  AliTPCCorrection::AddVisualCorrection(alignTrans1  ,201);
  AliTPCCorrection::AddVisualCorrection(alignTrans2  ,202);
  AliTPCCorrection::AddVisualCorrection(alignRot0    ,203);
  AliTPCCorrection::AddVisualCorrection(alignRot1    ,204);
  AliTPCCorrection::AddVisualCorrection(alignRot2    ,205);
 
}


void RegisterAliTPCBoundaryVoltError(){
  /// Register phi symetric E filed distortions
  /// 100-108 - A side 0 Field
  /// 110-118 - C side 0 Field
  /// 120-128 - A side +0.5 Field
  /// 130-138 - C side +0.5 Field
  /// 140-148 - A side -0.5 Field
  /// 150-158 - C side -0.5 Field

  Double_t vdrift = 2.64; // [cm/us]   // to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t T1 = 1.0;
  Double_t T2 = 1.0;
  Double_t wtP = -10.0 * (0.5*10) * vdrift /  ezField ; 
  Double_t wtM = -10.0 * (0.5*10) * vdrift / -ezField ; 

  printf("RegisterAliTPCBoundaryVoltError()\n");
  AliTPCComposedCorrection *corrField2D = (AliTPCComposedCorrection*) fCorrections->Get("TPCFCVoltError2D");    
  //
  if (!corrField2D){
    TObjArray *array=new TObjArray(16);
    Double_t val = 40.; // 1mm
    Float_t bound0[8] = { 0, 0,0,0,0,0,0,0};
    Float_t boundAi[8] = { 0, 0,0,0,0,0,0,0};
    Float_t boundCi[8] = { 0, 0,0,0,0,0,0,0};    
    for (Int_t ipar=0; ipar<8; ipar++){
      //
      boundaryVoltErrorA[ipar] = new AliTPCBoundaryVoltError;
      boundaryVoltErrorC[ipar] = new AliTPCBoundaryVoltError;
      boundaryVoltErrorA[ipar]->SetName(Form("BoundaryVoltErrorAsidePar%d",ipar));
      boundaryVoltErrorA[ipar]->SetTitle(Form("BoundaryVoltErrorAsidePar%d",ipar));
      boundaryVoltErrorC[ipar]->SetName(Form("BoundaryVoltErrorCsidePar%d",ipar));
      boundaryVoltErrorC[ipar]->SetTitle(Form("BoundaryVoltErrorCsidePar%d",ipar));
      for (Int_t jpar=0; jpar<8; jpar++) if (ipar!=jpar){
	boundAi[jpar]=0;
	boundCi[jpar]=0;
      }
      boundAi[ipar]=val;
      boundCi[ipar]=val;
      //
      boundaryVoltErrorA[ipar]->SetBoundariesA(boundAi);
      boundaryVoltErrorA[ipar]->SetBoundariesC(bound0);
      boundaryVoltErrorA[ipar]->InitBoundaryVoltErrorDistortion();  
      boundaryVoltErrorA[ipar]->SetOmegaTauT1T2(0.,1,1); 
      //
      Float_t tempboundAi[8] = { 0, 0,0,0,0,0,-boundCi[6],-boundCi[7]};
      boundaryVoltErrorC[ipar]->SetBoundariesA(tempboundAi);
      boundaryVoltErrorC[ipar]->SetBoundariesC(boundCi);
    
      boundaryVoltErrorC[ipar]->InitBoundaryVoltErrorDistortion();  
      boundaryVoltErrorC[ipar]->SetOmegaTauT1T2(0.,1,1); 
      array->AddAt(boundaryVoltErrorA[ipar],ipar);
      array->AddAt(boundaryVoltErrorC[ipar],ipar+8);
      boundaryVoltErrorA[ipar]->Print();
      boundaryVoltErrorC[ipar]->Print();
      AliTPCCorrection::AddVisualCorrection(boundaryVoltErrorA[ipar], 100+ipar); 
      AliTPCCorrection::AddVisualCorrection(boundaryVoltErrorC[ipar], 150+ipar);     
    }
    corrField2D = new AliTPCComposedCorrection;
    corrField2D->SetCorrections(array);
    corrField2D->SetOmegaTauT1T2(0,1.,1.);
    corrField2D->Print();    
    fCorrections->cd();
    corrField2D->SetName("TPCFCVoltError2D");
    corrField2D->SetTitle("TPCFCVoltError2D");
    corrField2D->Write("TPCFCVoltError2D");
  }else{
    TObjArray *array = (TObjArray*)corrField2D->GetCorrections();
    for (Int_t ipar=0; ipar<8; ipar++){
      boundaryVoltErrorA[ipar] = (AliTPCBoundaryVoltError*) array->At(ipar);
      boundaryVoltErrorC[ipar] = (AliTPCBoundaryVoltError*) array->At(ipar+8);      
    }
  }
  //
  // Register correction
  for (Int_t ipar=0; ipar<8; ipar++){
    AliTPCCorrection::AddVisualCorrection(boundaryVoltErrorA[ipar], 100+ipar); 
    AliTPCCorrection::AddVisualCorrection(boundaryVoltErrorC[ipar], 110+ipar);     
    //
    // correction for +-0.5 T setting
    AliTPCCorrection *corrField =0; 
    corrField=(AliTPCCorrection *)boundaryVoltErrorA[ipar]->Clone();
    corrField->SetOmegaTauT1T2(wtP,T1,T2);
    AliTPCCorrection::AddVisualCorrection(corrField,120+ipar);

    corrField=(AliTPCCorrection *)boundaryVoltErrorC[ipar]->Clone();
    corrField->SetOmegaTauT1T2(wtP,T1,T2);
    AliTPCCorrection::AddVisualCorrection(corrField,130+ipar);
    // correction for +-0.5 T setting
    corrField=(AliTPCCorrection *)boundaryVoltErrorA[ipar]->Clone();
    corrField->SetOmegaTauT1T2(wtM,T1,T2);
    AliTPCCorrection::AddVisualCorrection(corrField,140+ipar);

    corrField=(AliTPCCorrection *)boundaryVoltErrorC[ipar]->Clone();
    corrField->SetOmegaTauT1T2(wtM,T1,T2);
    AliTPCCorrection::AddVisualCorrection(corrField,150+ipar);
  }
  
}



void RegisterAliTPCExBShape(){
  ///

  AliMagF *magF = new AliMagF("mag","mag");

  exbShape             = new AliTPCExBBShape;
  exbShape->SetBField(magF);
  exbShape->SetName("TPCExBShape");
  exbShape->SetTitle("TPCExBShape");
  exbShape->SetOmegaTauT1T2(0,1.,1.);
  exbShape->Print();   
  AliTPCCorrection::AddVisualCorrection(exbShape,500); 
  exbShapeT1X             = new AliTPCExBBShape;
  exbShapeT1X->SetBField(magF);
  exbShapeT1X->SetName("TPCExbShapeT1X");
  exbShapeT1X->SetTitle("TPCExbShapeT1X");
  exbShapeT1X->SetOmegaTauT1T2(0,1.2,1.);
  exbShapeT1X->Print();   
  AliTPCCorrection::AddVisualCorrection(exbShapeT1X,501); 
  exbShapeT2X             = new AliTPCExBBShape;
  exbShapeT2X->SetBField(magF);
  exbShapeT2X->SetName("TPCExbShapeT2X");
  exbShapeT2X->SetTitle("TPCExbShapeT2X");
  exbShapeT2X->SetOmegaTauT1T2(0,1.0,1.2);
  exbShapeT2X->Print();   
  AliTPCCorrection::AddVisualCorrection(exbShapeT2X,502); 
}


void RegisterAliTPCExBTwist(){
  ///

  twistX    = new  AliTPCExBTwist;
  twistY    = new  AliTPCExBTwist;
  twistX->SetXTwist(0.001);  // 1 mrad twist in x
  twistX->SetName("ExBTwistX");
  twistX->SetTitle("ExBTwistX");
  twistY->SetYTwist(0.001);  // 1 mrad twist in y
  twistY->SetName("ExBTwistY");
  twistY->SetTitle("ExBTwistY");
  twistX->SetOmegaTauT1T2(0,1.,1.);
  twistY->SetOmegaTauT1T2(0,1.,1.);      
  AliTPCCorrection::AddVisualCorrection(twistX,600); 
  AliTPCCorrection::AddVisualCorrection(twistY,601); 
}

void RegisterAliTPCCorrectionDrift(){
  /// Drift distortion/correction

  for (Int_t idrift=0; idrift<7; idrift++) {
    calibDrift[idrift]=new AliTPCCorrectionDrift;
  }
  calibDrift[0]->SetName("driftT0");
  calibDrift[0]->SetTitle("driftT0");
  calibDrift[0]->fZ0Aside=0.1; 
  calibDrift[0]->fZ0Cside=0.1;
  calibDrift[1]->SetName("driftScale0");
  calibDrift[1]->SetTitle("driftScale0");
  calibDrift[1]->fVScale0=0.001;
  calibDrift[2]->SetName("driftScaleR");
  calibDrift[2]->SetTitle("driftScaleR");
  calibDrift[2]->fVScaleR=0.001;
  calibDrift[3]->SetName("driftScaleX");
  calibDrift[3]->SetTitle("driftScaleX");
  calibDrift[3]->fVScaleX=0.001;
  calibDrift[4]->SetName("driftScaleY");
  calibDrift[4]->SetTitle("driftScaleY");
  calibDrift[4]->fVScaleY=0.001;
  //
  calibDrift[5]->SetName("driftIROCDZ");
  calibDrift[5]->SetTitle("driftIROCDZ");
  calibDrift[5]->fIROCZ0=0.1; //delta Z for IROCORC
  //
  calibDrift[6]->SetName("driftOROCDT");
  calibDrift[6]->SetTitle("driftOROCDT");
  calibDrift[6]->fOROCDZ=0.001; //delta Z for IROCORC
  

}


void RegisterAliTPCROCVoltError3D(){
  /// ROC rotation transformation
  ///       700 -709 - 0 field
  ///       710 -719 - +0.5 field
  ///       720 -729 - -0.5 field

  Double_t vdrift = 2.64; // [cm/us]   // to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t T1 = 1.0;
  Double_t T2 = 1.0;
  Double_t wtP = -10.0 * (0.5*10) * vdrift /  ezField ; 
  Double_t wtM = -10.0 * (0.5*10) * vdrift / -ezField ; 

  //
  rocRotgXA=0;     // roc rotation A side - inclination in X
  rocRotgYA=0;     // roc rotation A side - inclination in Y
  rocRotgXC=0;     // roc rotation C side - inclination in X
  rocRotgYC=0;     // roc rotation C side - inclination in Y
  rocDzIROCA=0;      // roc shift A side - in Z
  rocDzIROCC=0;      // roc shift C side - in Z
  //
  rocDzUDA=0;        // roc shift up-down - A side
  rocDzUDC=0;        // roc shift up-down - C side
  //
  rocRotIROCA=0;      // roc rot IROC A side - in Z
  rocRotIROCC=0;      // roc rot OROC C side - in Z
  rocRotUDA=0;        // roc rot updown A side
  rocRotUDC=0;        // roc rot updown C side
 //
  printf("RegisterAliTPCROCVoltError3D()");
  Double_t kAngle=0.001;
  // reference in lx
  AliTPCROC * rocInfo = AliTPCROC::Instance();
  Double_t lxRef  = (rocInfo->GetPadRowRadii(0,62)+rocInfo->GetPadRowRadii(36,0))/2;
  //
  TMatrixD matrix(72,3);
 
  AliTPCComposedCorrection *corrField3D = 0;
  TFile *fCorrectionsROC=0;
  fCorrectionsROC = TFile::Open("TPCCorrectionPrimitivesROC.root");
  corrField3D = ( AliTPCComposedCorrection *)fCorrectionsROC->Get("TPCROCVoltError3DRotationgXgY");
  //
  if (!corrField3D){
    fCorrectionsROC = TFile::Open("TPCCorrectionPrimitivesROC.root","recreate");
  }  
  if (corrField3D) { // load from file
    corrField3D->Print();
    TCollection *iter = corrField3D->GetCorrections();
   
    rocRotgXA=(AliTPCROCVoltError3D*)iter->FindObject("rocRotgXA");    
   
    rocRotgYA=(AliTPCROCVoltError3D*)iter->FindObject("rocRotgYA");  
   
    rocRotgXC=(AliTPCROCVoltError3D*)iter->FindObject("rocRotgXC");  
   
    rocRotgYC=(AliTPCROCVoltError3D*)iter->FindObject("rocRotgYC");  
    
    rocDzIROCA=(AliTPCROCVoltError3D*)iter->FindObject("rocDzIROCA");  
   
    rocDzIROCC=(AliTPCROCVoltError3D*)iter->FindObject("rocDzIROCC");  
    //
    rocDzUDA=(AliTPCROCVoltError3D*)iter->FindObject("rocDzUDA");  
   
    rocDzUDC=(AliTPCROCVoltError3D*)iter->FindObject("rocDzUDC");  

    rocRotIROCA=(AliTPCROCVoltError3D*)iter->FindObject("rocRotIROCA");  
   
    rocRotIROCC=(AliTPCROCVoltError3D*)iter->FindObject("rocRotIROCC");  

    rocRotUDA=(AliTPCROCVoltError3D*)iter->FindObject("rocRotUDA");  
   
    rocRotUDC=(AliTPCROCVoltError3D*)iter->FindObject("rocRotUDC");  
     
  } else {
    corrField3D = new AliTPCComposedCorrection;
    rocRotgXA=new AliTPCROCVoltError3D;    
    rocRotgYA=new AliTPCROCVoltError3D;  
    rocRotgXC=new AliTPCROCVoltError3D;  
    rocRotgYC=new AliTPCROCVoltError3D;  
    rocDzIROCA=new AliTPCROCVoltError3D;  
    rocDzIROCC=new AliTPCROCVoltError3D;  
    rocDzUDA=new AliTPCROCVoltError3D;  
    rocDzUDC=new AliTPCROCVoltError3D;  
    rocRotIROCA=new AliTPCROCVoltError3D;  
    rocRotIROCC=new AliTPCROCVoltError3D;  
    rocRotUDA=new AliTPCROCVoltError3D;  
    rocRotUDC=new AliTPCROCVoltError3D;  
    //
    for (Int_t isec=0; isec<72; isec++){
      Double_t secAlpha = TMath::DegToRad()*(10.+20.*(((Int_t)isec)%18));
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36<18){
        matrix(isec,0) = kAngle*TMath::Cos(secAlpha)*lxRef;
        matrix(isec,1) = kAngle*TMath::Cos(secAlpha);
        matrix(isec,2) = -kAngle*TMath::Sin(secAlpha);
      }
    }
    rocRotgXA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      Double_t secAlpha = TMath::DegToRad()*(10.+20.*(((Int_t)isec)%18));
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36<18){
        matrix(isec,0) = kAngle*TMath::Sin(secAlpha)*lxRef;
        matrix(isec,1) = kAngle*TMath::Sin(secAlpha);
        matrix(isec,2) = kAngle*TMath::Cos(secAlpha);
      }
    }
    rocRotgYA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
     Double_t secAlpha = TMath::DegToRad()*(10.+20.*(((Int_t)isec)%18));
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36>=18){
        matrix(isec,0) = kAngle*TMath::Cos(secAlpha)*lxRef;
        matrix(isec,1) = kAngle*TMath::Cos(secAlpha);
        matrix(isec,2) = -kAngle*TMath::Sin(secAlpha);
      }
    }
    rocRotgXC->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      Double_t secAlpha = TMath::DegToRad()*(10.+20.*(((Int_t)isec)%18));
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36>=18){
        matrix(isec,0) = kAngle*TMath::Sin(secAlpha)*lxRef;
        matrix(isec,1) = kAngle*TMath::Sin(secAlpha);
        matrix(isec,2) = kAngle*TMath::Cos(secAlpha);
      }
    }
    rocRotgYC->SetROCData(&matrix);

    //
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec<18){
        matrix(isec,0) = 0.1;  // 1 mm 
        matrix(isec,1) = 0;
        matrix(isec,2) = 0;
      }
    }
    rocDzIROCA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec>=18 && isec<36){
        matrix(isec,0) = 0.1;  // 1 mm 
        matrix(isec,1) = 0;
        matrix(isec,2) = 0;
      }
    }
    rocDzIROCC->SetROCData(&matrix);

    //
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36<18){
        matrix(isec,0) = (isec%18<9)? 0.05:-0.05;  // 1 mm 
        matrix(isec,1) = 0;
        matrix(isec,2) = 0;
      }
    }
    rocDzUDA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36>=18){
        matrix(isec,0) = (isec%18<9)?0.05:-0.05;  // 1 mm 
        matrix(isec,1) = 0;
        matrix(isec,2) = 0;
      }
    }
    rocDzUDC->SetROCData(&matrix);
    //
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec<18){
        matrix(isec,0) = 0;   
        matrix(isec,1) = kAngle;
        matrix(isec,2) = 0;
      }
    }
    //
    rocRotIROCA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec>=18 && isec<36){
        matrix(isec,0) = 0;
        matrix(isec,1) = kAngle;
        matrix(isec,2) = 0;
      }
    }
    rocRotIROCC->SetROCData(&matrix);
    //
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36<18){
        matrix(isec,0) = 0;   
        matrix(isec,1) = (isec%18<9)?kAngle:-kAngle;
        matrix(isec,2) = 0;
      }
    }
    //
    rocRotUDA->SetROCData(&matrix);
    //
    for (Int_t isec=0; isec<72; isec++){
      matrix(isec,0)=0; matrix(isec,1)=0; matrix(isec,2)=0;
      if (isec%36>=18){
        matrix(isec,0) = 0;
        matrix(isec,1) = (isec%18<9)?kAngle:-kAngle;
        matrix(isec,2) = 0;
      }
    }
    rocRotUDC->SetROCData(&matrix);


    //
    rocRotgXA->SetElectronArrivalCorrection(kFALSE);
    rocRotgYA->SetElectronArrivalCorrection(kFALSE);
    rocRotgXC->SetElectronArrivalCorrection(kFALSE);
    rocRotgYC->SetElectronArrivalCorrection(kFALSE);
    rocDzIROCA->SetElectronArrivalCorrection(kFALSE);
    rocDzIROCC->SetElectronArrivalCorrection(kFALSE);
    rocDzUDA->SetElectronArrivalCorrection(kFALSE);
    rocDzUDC->SetElectronArrivalCorrection(kFALSE);
    rocRotIROCA->SetElectronArrivalCorrection(kFALSE);
    rocRotIROCC->SetElectronArrivalCorrection(kFALSE);
    rocRotUDA->SetElectronArrivalCorrection(kFALSE);
    rocRotUDC->SetElectronArrivalCorrection(kFALSE);

    /* // verification plots
    rocRotgXA.CreateHistoOfZAlignment(0,500,500)->Draw("surf2"); 
    rocRotgYA.CreateHistoOfZAlignment(0,500,500)->Draw("surf2"); 
    rocRotgXC.CreateHistoOfZAlignment(1,500,500)->Draw("surf2"); 
    rocRotgYC.CreateHistoOfZAlignment(1,500,500)->Draw("surf2"); 
    */

    //
    rocRotgXA->SetName("rocRotgXA");rocRotgXA->SetTitle("rocRotgXA");
    rocRotgYA->SetName("rocRotgYA");rocRotgYA->SetTitle("rocRotgYA");
    rocRotgXC->SetName("rocRotgXC");rocRotgXC->SetTitle("rocRotgXC");
    rocRotgYC->SetName("rocRotgYC");rocRotgYC->SetTitle("rocRotgYC");
    rocDzIROCA->SetName("rocDzIROCA");rocDzIROCA->SetTitle("rocDzIROCA");
    rocDzIROCC->SetName("rocDzIROCC");rocDzIROCC->SetTitle("rocDzIROCC");
    rocDzUDA->SetName("rocDzUDA");rocDzUDA->SetTitle("rocDzUDA");
    rocDzUDC->SetName("rocDzUDC");rocDzUDC->SetTitle("rocDzUDC");
    rocRotIROCA->SetName("rocRotIROCA");rocRotIROCA->SetTitle("rocRotIROCA");
    rocRotIROCC->SetName("rocRotIROCC");rocRotIROCC->SetTitle("rocRotIROCC");
    rocRotUDA->SetName("rocRotUDA");rocRotUDA->SetTitle("rocRotUDA");
    rocRotUDC->SetName("rocRotUDC");rocRotUDC->SetTitle("rocRotUDC");
    //
    //
    TObjArray *cs = new TObjArray();
    cs->Add(rocRotgXA);
    cs->Add(rocRotgYA);
    cs->Add(rocRotgXC);
    cs->Add(rocRotgYC);
    cs->Add(rocDzIROCA);
    cs->Add(rocDzIROCC);
    cs->Add(rocDzUDA);
    cs->Add(rocDzUDC);
    cs->Add(rocRotIROCA);
    cs->Add(rocRotIROCC);
    cs->Add(rocRotUDA);
    cs->Add(rocRotUDC);
    corrField3D->SetCorrections(cs);
    corrField3D->SetOmegaTauT1T2(0,1.,1.);
    corrField3D->Print();
    fCorrectionsROC->cd();
    corrField3D->Init();
    corrField3D->Print("da");
    fCorrectionsROC->cd();
    corrField3D->Write("TPCROCVoltError3DRotationgXgY");
  }
  AliTPCCorrection::AddVisualCorrection(rocRotgXA,701); 
  AliTPCCorrection::AddVisualCorrection(rocRotgYA,702); 
  AliTPCCorrection::AddVisualCorrection(rocRotgXC,703); 
  AliTPCCorrection::AddVisualCorrection(rocRotgYC,704); 
  AliTPCCorrection::AddVisualCorrection(rocDzIROCA,705); 
  AliTPCCorrection::AddVisualCorrection(rocDzIROCC,706); 
  AliTPCCorrection::AddVisualCorrection(rocDzUDA,709); 
  AliTPCCorrection::AddVisualCorrection(rocDzUDC,710); 
  AliTPCCorrection::AddVisualCorrection(rocRotIROCA,707); 
  AliTPCCorrection::AddVisualCorrection(rocRotIROCC,708); 
  AliTPCCorrection::AddVisualCorrection(rocRotUDA,711); 
  AliTPCCorrection::AddVisualCorrection(rocRotUDC,712); 

  AliTPCCorrection *corrPlus =0; 
  //
  corrPlus=(AliTPCCorrection *)rocRotgXA->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,711);
  //
  corrPlus=(AliTPCCorrection *)rocRotgYA->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,712);
  //
  corrPlus=(AliTPCCorrection *)rocRotgXC->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,713);
  //
  corrPlus=(AliTPCCorrection *)rocRotgYC->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,714);
  //
  corrPlus=(AliTPCCorrection *)rocDzIROCA->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,715);
  //
  corrPlus=(AliTPCCorrection *)rocDzIROCC->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,716);
  //
  corrPlus=(AliTPCCorrection *)rocRotIROCA->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,717);
  //
  corrPlus=(AliTPCCorrection *)rocDzIROCC->Clone();
  corrPlus->SetOmegaTauT1T2(wtP,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrPlus,718);
  //
  //
  AliTPCCorrection *corrMinus =0; 
  //
  corrMinus=(AliTPCCorrection *)rocRotgXA->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,721);
  //
  corrMinus=(AliTPCCorrection *)rocRotgYA->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,722);
  //
  corrMinus=(AliTPCCorrection *)rocRotgXC->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,723);
  //
  corrMinus=(AliTPCCorrection *)rocRotgYC->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,724);
  //
  corrMinus=(AliTPCCorrection *)rocDzIROCA->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,725);
  //
  corrMinus=(AliTPCCorrection *)rocDzIROCC->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,726);
  //
  corrMinus=(AliTPCCorrection *)rocRotIROCA->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,727);
  //
  corrMinus=(AliTPCCorrection *)rocDzIROCC->Clone();
  corrMinus->SetOmegaTauT1T2(wtM,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corrMinus,728);
  //

  fCorrectionsROC->Close();
  delete fCorrectionsROC;
}


void RegisterAliTPCROCVoltError3DSector(){
  /// ROC rotation and shift transformation
  /// 800-819 -   0.0 Field
  /// 820-839 -  +0.5 Field
  /// 840-859 -  +0.5 Field

  rocShiftIROCA0=0;      // IROC shift A0 side
  rocRotIROCA0=0;        // IROC rot   A0 side
  rocShiftOROCA0=0;      // OROC shift A0 side
  rocRotOROCA0=0;        // OROC rot   A0 side
  rocShiftIROCC0=0;      // IROC shift C0 side
  rocRotIROCC0=0;        // IROC rot   C0 side
  rocShiftOROCC0=0;      // OROC shift C0 side
  rocRotOROCC0=0;        // OROC rot   C0 side
  Double_t vdrift = 2.64; // [cm/us]   // to be updated: per second (ideally)
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t T1 = 1.0;
  Double_t T2 = 1.0;

  //  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  //
  //
  printf("RegisterAliTPCROCVoltError3DSector()");
  Double_t kAngle=0.001;
  Double_t kDz=0.1;
  // reference in lx
  //
  TMatrixD matrix(72,3);
 
  AliTPCComposedCorrection *corrField3DSector = 0;
  TFile *fCorrectionsROC=0;
  fCorrectionsROC = TFile::Open("TPCCorrectionPrimitivesSector.root");
  corrField3DSector = ( AliTPCComposedCorrection *)fCorrectionsROC->Get("TPCROCVoltError3DSector");
  //
  if (!corrField3DSector){
    fCorrectionsROC = TFile::Open("TPCCorrectionPrimitivesSector.root","recreate");
  }  
  if (corrField3DSector) { // load from file
    corrField3DSector->Print();
    TCollection *iter = corrField3DSector->GetCorrections();
    //
    rocShiftIROCA0=(AliTPCROCVoltError3D*)iter->FindObject("rocShiftIROCA0");   // IROC shift A0 side
    rocRotIROCA0=(AliTPCROCVoltError3D*)iter->FindObject("rocRotIROCA0");       // IROC rot   A0 side
    rocShiftOROCA0=(AliTPCROCVoltError3D*)iter->FindObject("rocShiftOROCA0");   // OROC shift A0 side
    rocRotOROCA0=(AliTPCROCVoltError3D*)iter->FindObject("rocRotOROCA0");       // OROC rot   A0 side
    rocShiftIROCC0=(AliTPCROCVoltError3D*)iter->FindObject("rocShiftIROCC0");   // IROC shift C0 side
    rocRotIROCC0=(AliTPCROCVoltError3D*)iter->FindObject("rocRotIROCC0");       // IROC rot   C0 side
    rocShiftOROCC0=(AliTPCROCVoltError3D*)iter->FindObject("rocShiftOROCC0");   // OROC shift C0 side
    rocRotOROCC0=(AliTPCROCVoltError3D*)iter->FindObject("rocRotOROCC0");       // OROC rot   C0 side
     
  } else {
    corrField3DSector = new AliTPCComposedCorrection;
    rocShiftIROCA0=new AliTPCROCVoltError3D;      // IROC shift A0 side
    rocRotIROCA0=new AliTPCROCVoltError3D;        // IROC rot   A0 side
    rocShiftOROCA0=new AliTPCROCVoltError3D;      // OROC shift A0 side
    rocRotOROCA0=new AliTPCROCVoltError3D;        // OROC rot   A0 side
    rocShiftIROCC0=new AliTPCROCVoltError3D;      // IROC shift C0 side
    rocRotIROCC0=new AliTPCROCVoltError3D;        // IROC rot   C0 side
    rocShiftOROCC0=new AliTPCROCVoltError3D;      // OROC shift C0 side
    rocRotOROCC0=new AliTPCROCVoltError3D;        // OROC rot   C0 side
    //
    matrix.Zero(); matrix(0,0)=kDz; 
    rocShiftIROCA0->SetROCData(&matrix);
    matrix.Zero(); matrix(0,1)=kAngle; 
    rocRotIROCA0->SetROCData(&matrix);
    matrix.Zero(); matrix(36,0)=kDz; 
    rocShiftOROCA0->SetROCData(&matrix);
    matrix.Zero(); matrix(36,1)=kAngle; 
    rocRotOROCA0->SetROCData(&matrix);

    matrix.Zero(); matrix(18,0)=kDz; 
    rocShiftIROCC0->SetROCData(&matrix);
    matrix.Zero(); matrix(18,1)=kAngle; 
    rocRotIROCC0->SetROCData(&matrix);
    matrix.Zero(); matrix(36+18,0)=kDz; 
    rocShiftOROCC0->SetROCData(&matrix);
    matrix.Zero(); matrix(36+18,1)=kAngle; 
    rocRotOROCC0->SetROCData(&matrix);
    //
    rocShiftIROCA0->SetElectronArrivalCorrection(kFALSE);     // IROC shift A0 side
    rocRotIROCA0->SetElectronArrivalCorrection(kFALSE);        // IROC rot   A0 side
    rocShiftOROCA0->SetElectronArrivalCorrection(kFALSE);      // OROC shift A0 side
    rocRotOROCA0->SetElectronArrivalCorrection(kFALSE);        // OROC rot   A0 side
    rocShiftIROCC0->SetElectronArrivalCorrection(kFALSE);      // IROC shift C0 side
    rocRotIROCC0->SetElectronArrivalCorrection(kFALSE);        // IROC rot   C0 side
    rocShiftOROCC0->SetElectronArrivalCorrection(kFALSE);      // OROC shift C0 side
    rocRotOROCC0->SetElectronArrivalCorrection(kFALSE);        // OROC rot   C0 side

    /* // verification plots
    */
    //
    rocShiftIROCA0->SetName("rocShiftIROCA0");rocShiftIROCA0->SetTitle("rocShiftIROCA0");
    rocRotIROCA0->SetName("rocRotIROCA0");rocRotIROCA0->SetTitle("rocRotIROCA0");
    rocShiftOROCA0->SetName("rocShiftOROCA0"); rocShiftOROCA0->SetTitle("rocShiftOROCA0");
    rocRotOROCA0->SetName("rocRotOROCA0");rocRotOROCA0->SetTitle("rocRotOROCA0");
    //
    rocShiftIROCC0->SetName("rocShiftIROCC0");rocShiftIROCC0->SetTitle("rocShiftIROCC0");
    rocRotIROCC0->SetName("rocRotIROCC0");rocRotIROCC0->SetTitle("rocRotIROCC0");
    rocShiftOROCC0->SetName("rocShiftOROCC0");rocShiftOROCC0->SetTitle("rocShiftOROCC0");
    rocRotOROCC0->SetName("rocRotOROCC0");rocRotOROCC0->SetTitle("rocRotOROCC0");
    //
    //
    TObjArray *cs = new TObjArray();
    cs->Add(rocShiftIROCA0);      // IROC shift A0 side
    cs->Add(rocRotIROCA0);        // IROC rot   A0 side
    cs->Add(rocShiftOROCA0);      // OROC shift A0 side
    cs->Add(rocRotOROCA0);        // OROC rot   A0 side
    cs->Add(rocShiftIROCC0);      // IROC shift C0 side
    cs->Add(rocRotIROCC0);        // IROC rot   C0 side
    cs->Add(rocShiftOROCC0);      // OROC shift C0 side
    cs->Add(rocRotOROCC0);        // OROC rot   C0 side
    //
    corrField3DSector->SetCorrections(cs);
    corrField3DSector->SetOmegaTauT1T2(0,T1,T2);
    corrField3DSector->Print();
    //

    fCorrectionsROC->cd();
    corrField3DSector->Init();
    corrField3DSector->Print("da");
    fCorrectionsROC->cd();
    corrField3DSector->Write("TPCROCVoltError3DSector");
  }
  AliTPCCorrection::AddVisualCorrection(rocShiftIROCA0,800);      // IROC shift A0 side
  AliTPCCorrection::AddVisualCorrection(rocRotIROCA0,801);        // IROC rot   A0 side
  AliTPCCorrection::AddVisualCorrection(rocShiftOROCA0,802);      // OROC shift A0 side
  AliTPCCorrection::AddVisualCorrection(rocRotOROCA0,803);        // OROC rot   A0 side
  AliTPCCorrection::AddVisualCorrection(rocShiftIROCC0,804);      // IROC shift C0 side
  AliTPCCorrection::AddVisualCorrection(rocRotIROCC0,805);        // IROC rot   C0 side
  AliTPCCorrection::AddVisualCorrection(rocShiftOROCC0,806);      // OROC shift C0 side
  AliTPCCorrection::AddVisualCorrection(rocRotOROCC0,807);        // OROC rot   C0 side 
  //
  // Register correction for plus setting
  AliTPCCorrection *corr =0; 
  Double_t wtp = -10.0 * (0.5*10) * vdrift / ezField ; 
  Double_t wtm = -10.0 * (0.5*10) * vdrift / -ezField ; 

  corr=(AliTPCCorrection *)rocShiftIROCA0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,820);           // IROC shift A0 side + Plus field
  //
  corr=(AliTPCCorrection *)rocRotIROCA0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,821);           // IROC rot   A0 side
  //
  corr=(AliTPCCorrection *)rocShiftOROCA0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,822);            // OROC shift   A0 side
  //
  corr=(AliTPCCorrection *)rocRotOROCA0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,823);            // OROC rot   A0 side
  corr=(AliTPCCorrection *)rocShiftIROCC0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,824);           // IROC shift C0 side + Plus field
  //
  corr=(AliTPCCorrection *)rocRotIROCC0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,825);           // IROC rot   C0 side
  //
  corr=(AliTPCCorrection *)rocShiftOROCC0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,826);            // OROC shift   C0 side
  //
  corr=(AliTPCCorrection *)rocRotOROCC0->Clone();
  corr->SetOmegaTauT1T2(wtp,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,827);            // OROC rot   C0 side
  //
  corr=(AliTPCCorrection *)rocShiftIROCA0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,840);           // IROC shift A0 side + Plus field
  //
  corr=(AliTPCCorrection *)rocRotIROCA0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,841);           // IROC rot   A0 side
  //
  corr=(AliTPCCorrection *)rocShiftOROCA0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,842);            // OROC shift   A0 side
  //
  corr=(AliTPCCorrection *)rocRotOROCA0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,843);            // OROC rot   A0 side
  corr=(AliTPCCorrection *)rocShiftIROCC0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,844);           // IROC shift C0 side + Plus field
  //
  corr=(AliTPCCorrection *)rocRotIROCC0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,845);           // IROC rot   C0 side
  //
  corr=(AliTPCCorrection *)rocShiftOROCC0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,846);            // OROC shift   C0 side
  //
  corr=(AliTPCCorrection *)rocRotOROCC0->Clone();
  corr->SetOmegaTauT1T2(wtm,T1,T2);
  AliTPCCorrection::AddVisualCorrection(corr,847);            // OROC rot   C0 side
  //
  fCorrectionsROC->Close();
  delete fCorrectionsROC;
}





AliTPCComposedCorrection * MakeComposedCorrectionExB(){
  /// make composed corection for ExB scanning

  RegisterCorrection();
  //
   //
  TObjArray * corr = new TObjArray;       // primitive corrections - for fitting
  TObjArray * testCorr = new TObjArray;   // test corrections - to be used as benchmark for fitting
  //
  corr->AddLast(twistX);
  corr->AddLast(twistY);
  corr->AddLast(exbShape);
  corr->AddLast(exbShapeT1X);
  corr->AddLast(exbShapeT2X);
  //
  for (Int_t i=0; i<8; i++){
    corr->AddLast(boundaryVoltErrorA[i]);
    corr->AddLast(boundaryVoltErrorC[i]);
  }
  //
  // ROD alignment
  //
  corr->AddLast(rodIFCCos);
  corr->AddLast(rodIFCSin);
  corr->AddLast(rodOFCSin);
  corr->AddLast(rodOFCCos);
  //alignment
  corr->AddLast(alignTrans0);
  corr->AddLast(alignTrans1);
  corr->AddLast(alignTrans2);
  corr->AddLast(alignRot0);
  corr->AddLast(alignRot1);
  corr->AddLast(alignRot2);
  //
  corr->AddLast(alignTrans0D[0]);
  corr->AddLast(alignTrans1D[0]);
  corr->AddLast(alignTrans2D[0]);
  corr->AddLast(alignRot0D[0]);
  corr->AddLast(alignRot1D[0]);
  corr->AddLast(alignRot2D[0]);
  corr->AddLast(alignTrans0D[1]);
  corr->AddLast(alignTrans1D[1]);
  corr->AddLast(alignRot0D[1]);
  corr->AddLast(alignTrans0D[2]);
  corr->AddLast(alignTrans1D[2]);
  corr->AddLast(alignRot0D[2]);
  //
  // z alignment + E field distortion due z misalignment
  //
  corr->AddLast(rocRotgXA);  // A side C side
  corr->AddLast(rocRotgYA);
  corr->AddLast(rocRotgXC);
  corr->AddLast(rocRotgYC);
  corr->AddLast(rocDzIROCA);
  corr->AddLast(rocDzIROCC);
  corr->AddLast(rocRotIROCA);
  corr->AddLast(rocRotIROCC);
  corr->AddLast(rocDzUDA);
  corr->AddLast(rocDzUDC);
  corr->AddLast(rocRotUDA);
  corr->AddLast(rocRotUDC);
  //
  //
  corr->AddLast(calibDrift[0]);
  corr->AddLast(calibDrift[1]);
  corr->AddLast(calibDrift[2]);
  corr->AddLast(calibDrift[3]);
  corr->AddLast(calibDrift[4]);
  corr->AddLast(calibDrift[5]);
  corr->AddLast(calibDrift[6]);
  //
  // setup test correction 
  // 
  testCorr->AddLast(rodIFCCos);
  testCorr->AddLast(rodIFCSin);
  testCorr->AddLast(rodOFCSin);
  testCorr->AddLast(rodOFCCos);
  //
  testCorr->AddLast(twistX);
  testCorr->AddLast(twistY);
  testCorr->AddLast(alignTrans0);
  testCorr->AddLast(alignTrans1);
  testCorr->AddLast(alignTrans2);
  testCorr->AddLast(alignRot0);
  testCorr->AddLast(alignRot1);
  testCorr->AddLast(alignRot2);  
  testCorr->AddLast(rocRotgXA);  // A side C side
  testCorr->AddLast(rocRotgYA);
  testCorr->AddLast(rocRotgXC);
  testCorr->AddLast(rocRotgYC);
  testCorr->AddLast(rocDzIROCA);
  testCorr->AddLast(rocDzIROCC);
  testCorr->AddLast(rocRotIROCA);
  testCorr->AddLast(rocRotIROCC);
  Int_t entries=testCorr->GetEntries();
  TVectorD weights(entries);
  for (Int_t i=0; i<entries; i++) weights[i]=1;
  //
  AliTPCComposedCorrection *composedTest= new AliTPCComposedCorrection ;
  composedTest->SetName("FitSample");
  composedTest->SetTitle("FitSample");
  composedTest->SetCorrections(testCorr);
  composedTest->SetWeights(&weights);
  //
  corr->AddLast(composedTest);
  AliTPCComposedCorrection *cc= new AliTPCComposedCorrection ;
  cc->SetCorrections((TObjArray*)(corr));
   //cc->Init();
  cc->Print("DA"); // Print used correction classes
  cc->SetName("ComposedExB");
  TFile fexb("RegisterCorrectionExB.root","recreate");
  cc->Write("ComposedExB");
  fexb.Close();
  return cc;
}


AliTPCComposedCorrection * GetCorrectionFromFile(){
  /// Getthe appropariate correction form the closest file

  TFile * fexb= TFile::Open("RegisterCorrectionExB.root");
  if (!fexb)  fexb= TFile::Open("../RegisterCorrectionExB.root");
  if (!fexb)  fexb= TFile::Open("../../RegisterCorrectionExB.root");
  if (!fexb)  fexb= TFile::Open("../../../RegisterCorrectionExB.root");
  //
  if (!fexb)   return 0;
  TFile * fitter= TFile::Open("fitCorrection.root");

  if (!fitter)  fitter= TFile::Open("../fitCorrection.root");
  if (!fitter)  fitter= TFile::Open("../../fitCorrection.root");
  if (!fitter)  fitter= TFile::Open("../../../fitCorrection.root");
  //
  AliTPCComposedCorrection *cc=  (AliTPCComposedCorrection*) fexb->Get("ComposedExB");
  {if (!cc){
    printf("File or correction RegisterCorrectionExB.root doees not exist or corrupted\n\n\n");
    return 0;
    }}
  TObjArray * corr = (TObjArray*)(cc->GetCorrections());
  //  TObjArray * corrLocal =new TObjArray;
  // TObjArray * corrGlobal =new TObjArray;
  //
  if (fitter){
    if (fitter->GetKey("FitBoundary")) corr->AddLast(fitter->Get("FitBoundary"));
    if (fitter->GetKey("FitExBTwist"))  corr->AddLast(fitter->Get("FitExBTwist"));
    if (fitter->GetKey("FitAlignGlobal"))  corr->AddLast(fitter->Get("FitAlignGlobal"));
    if (fitter->GetKey("FitRodAlignGloba"))  corr->AddLast(fitter->Get("FitRodAlignGlobal"));
    if (fitter->GetKey("FitRocAlignGlobal"))  corr->AddLast(fitter->Get("FitRocAlignGlobal"));  
    if (fitter->GetKey("FitRocAlignZ"))  corr->AddLast(fitter->Get("FitRocAlignZ"));
    if (fitter->GetKey("FitAlignLocal"))  corr->AddLast(fitter->Get("FitAlignLocal"));
    if (fitter->GetKey("FitAlignTPC"))  corr->AddLast(fitter->Get("FitAlignTPC"));
    if (fitter->GetKey("FitAlignTOF"))  corr->AddLast(fitter->Get("FitAlignTOF"));
    if (fitter->GetKey("FitAlignTRD"))  corr->AddLast(fitter->Get("FitAlignTRD"));
  }
  return cc;
}

void TestParExample(){
  /// dz shift example: AliTPCCorrection::AddVisualCorrection(rocDzIROCA,705);
  /// => parabolic fit and helix fit agrees - once significant ammount of points used
  ///    160 point - agreement  ~2%;
  ///     80 points - agreement ~5%

  AliTPCCorrection* corr = AliTPCCorrection::GetVisualCorrection(705);  
  corr->SetOmegaTauT1T2(0.33,1,1);
  TF1 f705Par("f705","AliTPCCorrectionFit::EvalAtPar(0,0,85, x,0.1,705,0,80)",0,360);
  f705Par.SetLineColor(2);f705Par.SetNpx(500);
  TF1 f705Helix("f705Helix","AliTPCCorrectionFit::EvalAtHelix(0,0,85,x,0.1,705,0,80)",0,360);
  f705Helix.SetLineColor(4);f705Helix.SetNpx(500);
  f705Helix.Draw();
  f705Par.Draw("same");

}



void TestFitSpeed(Int_t nEvals){
  ///  test speed of helix fit/ resp. parabolic fir

  TStopwatch timerh; 
  ::Info("TestFitSpeed","Helix fit");
  for (Int_t i=0;i<nEvals; i++) AliTPCCorrectionFit::EvalAtPar(0,0,85,gRandom->Rndm(),0.1,705,0,80); 
  timerh.Print();
  TStopwatch timerP; 
  ::Info("TestFitSpeed","Parabolicfit");
  for (Int_t i=0;i<nEvals; i++) AliTPCCorrectionFit::EvalAtPar(0,0,85,gRandom->Rndm(),0.1,705,0,80); 
  timerP.Print();
  /*
    For the test system CPU comsumption is the same:
    I-TestFitSpeed: Helix fit
    Real time 0:00:03, CP time 3.310
    I-TestFitSpeed: Parabolicfit
    Real time 0:00:03, CP time 3.280
  */
}

