/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: federico ronchetti         for the ALICE HLT Project.*
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALConstants.h"
#include "AliHLTEMCALConstant.h"
#include "assert.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliHLTEMCALSharedMemoryInterface.h" 
#include "TVector3.h"

using namespace EmcalHLTConst;

ClassImp(AliHLTEMCALGeometry);
TGeoManager *gGeoManager = 0;

AliHLTEMCALGeometry::AliHLTEMCALGeometry() :
	AliHLTCaloGeometry ("EMCAL"),
	fShmPtr(0),
	fGeo(0),
	fEMCALGeometry(0)
{
	//fGeo = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
	//fGeo = new AliEMCALGeometry("EMCAL_COMPLETE","EMCAL");
	//fGeo =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
	//TGeoManager::Import("/home/fedro/work/AliRoot/test/QA/geometry.root");
	//fGeo = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
	fShmPtr = new AliHLTEMCALSharedMemoryInterface();
	GetGeometryFromCDB();
}

Int_t AliHLTEMCALGeometry::InitialiseGeometry()
{
   
   return GetGeometryFromCDB();
}


AliHLTEMCALGeometry::~AliHLTEMCALGeometry()
{

}
  
void 
AliHLTEMCALGeometry::GetGlobalCoordinates(AliHLTEMCALRecPointDataStruct &recPoint, AliHLTCaloGlobalCoordinate &globalCoord)
{

	 Int_t istrip = 0;
	 Float_t z0 = 0;
	 Float_t zb = 0;
	 Float_t z_is = 0;
	 Float_t d = 0;

	 Float_t x,y,z; // return variables in terry's RF

	 Float_t dz = fCaloConstants->GetMINCELLSTEPETA(); // base cell width in eta
	 Float_t dx = fCaloConstants->GetCELLSTEPPHI(); // base cell width in phi

	  // parameters for shower depth calculation
	 Float_t X0 = fCaloConstants->GetRADLENGTH();
	 Float_t Ecr = fCaloConstants->GetCRITICENERGY();
	 Float_t Cj;

	 Float_t teta0 = fCaloConstants->GetCELLANGLE(); //tapering angle (deg)
	 Float_t teta1; //working angle

	 Float_t L = fCaloConstants->GetCELLHEIGHT();

	 // converting to MEV
	 Float_t E = recPoint.fAmp * 1000;

	 //TVector3 v1;
	 	Double_t glob[3];
	 	Double_t loc[3];

	 if (recPoint.fZ >= 47.5 || recPoint.fZ<-0.5) {
		 Logging(kHLTLogError, "HLT", "EMCAL", "AliHLTEMCALGeometry::GetGlobalCoordinates: invalid Z recpoint ");
		 return;
	  }

	  if (recPoint.fX >= 23.5 || recPoint.fX<-0.5) {
		  Logging(kHLTLogError, "HLT", "EMCAL", "AliHLTEMCALGeometry::GetGlobalCoordinates: invalid X recpoint ");
		  return;
	  }


	  switch ( recPoint.fParticle )
	  {
	  case 0:
		  Cj = - fCaloConstants->GetCJ(); // photon
		  d = X0 * TMath::Log( E / Ecr + Cj);
		  break;

	  case 1:
		  Cj = + fCaloConstants->GetCJ(); // electron
		  d = X0 * TMath::Log( E / Ecr + Cj);
		  break;

	  case 2:
		  // hadron
		  d = 0.5 * L;
		  break;

	  default:
		  Cj = - fCaloConstants->GetCJ(); // defaulting to photon
		  d = X0 * TMath::Log( E / Ecr + Cj);
	  }

	   istrip = int ((recPoint.fZ + 0.5 ) / 2);

	   // module angle
	   teta1 = TMath::DegToRad() * istrip * teta0;

	   // calculation of module corner along z
	   // as a function of strip

	   for (Int_t is=0; is<= istrip; is++) {

	     teta1 = TMath::DegToRad() * is * teta0;

	     z_is = z_is + 2*dz*(TMath::Sin(teta1)*TMath::Tan(teta1) + TMath::Cos(teta1));

	   }

	   z0 = dz * (recPoint.fZ - 2*istrip + 0.5);
	   zb = (2*dz-z0-d*TMath::Tan(teta1))*TMath::Cos(teta1);

	   z = z_is - zb*TMath::Cos(teta1);

//	   cout << "----> istrip: " << istrip << endl;
//	   cout << "----> z0: "<< z0 << endl;
//	   cout << "----> zb: "<< zb << endl;
//	   cout << "----> corner z: "<< z_is << endl;

//	   cout << "----> teta1: "<< TMath::RadToDeg()*teta1 << endl;

	   y = d/TMath::Cos(teta1) + zb*TMath::Sin(teta1);

	   x = (recPoint.fX + 0.5)*dx;

	   // cout << "x: " << x << " y: "<< y << " z " << z << endl;

	// check coordinate origin
	loc[0] = x;
	loc[1] = y;
	loc[2] = z;

	 if(!fGeo)
   {
      Logging(kHLTLogError, "HLT", "EMCAL", "AliHLTEMCALGeometry::GetGlobalCoordinates: no geometry initialised");
      return;
   }

	 ConvertRecPointCoordinates(loc[1], loc[2], loc[3]);

	 fGeo->GetGlobal(loc, glob, recPoint.fModule);

	 globalCoord.fX = glob[0];
	 globalCoord.fY = glob[1];
	 globalCoord.fZ = glob[2];
}
 
void 
AliHLTEMCALGeometry::GetCellAbsId(UInt_t module, UInt_t y, UInt_t z, Int_t& AbsId)
{

  	if(!fGeo)
      	{
		 Logging(kHLTLogError, "HLT", "EMCAL", "AliHLTEMCALGeometry::GetCellAbsId: no geometry initialised");
		 return;

      	}
	
	AbsId = fGeo->GetAbsCellIdFromCellIndexes(module, y, z);
	
}

void AliHLTEMCALGeometry::ConvertRecPointCoordinates(Double_t &x, Double_t &y, Double_t &z) const
{
	Double_t DX = 13.869008;
	Double_t DY = 72.559998;
	Double_t DZ = 175.00000;

	 x = y - DX;  //fixme
	 y = -x + DY; //fixme
	 z = z - DZ;  //fixme

}


int
AliHLTEMCALGeometry::GetGeometryFromCDB()
{
  cout << "Getting geometry..." << endl;

 // AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliCDBPath path("GRP","Geometry","Data");
  if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry)
	{
	  if(!fEMCALGeometry)
	    {
	      delete fEMCALGeometry;
	      fEMCALGeometry = 0;
	    }

	  gGeoManager = (TGeoManager*) pEntry->GetObject();
	  cout<< "gGeoManager = 0x%x" << gGeoManager << endl;

	  if(gGeoManager)
	    {
	      fGeo = new AliEMCALGeoUtils("EMCAL_COMPLETE","EMCAL");
	      //fClusterAnalyserPtr->SetGeometry(fEMCALGeometry);
	    }

	}
      else
	{
	 	  //HLTError("can not fetch object \"%s\" from OCDB", path);
	}
    }
  return 0;
}
