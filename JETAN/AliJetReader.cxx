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
 
//------------------------------------------------------------------------
// Jet reader base class
// manages the reading of input for jet algorithms
// Authors: jgcn@mda.cinvestav.mx
//          magali.estienne@IReS.in2p3.fr
//
// **February 2011
// implemented  standard geometry (AliEMCALGeoUtils) (was AliJetDummyGeo implented separately in ESDReader and AODReader
// local2master matrices are now get from $ALICE_ROOT/OADB/PWG4/JetReconstruction/EMCALlocal2master.root
// you can choose the geometry (EMCAL_COMPLETE, EMCAL_FIRSTYEARv1, etc) via SetEMCALgeo2bLoad('Name_of_Geometry') in the Readerheader
// different options for survey(ed) matrice are provided too
// marco.bregant@subatech.in2p3.fr
//------------------------------------------------------------------------- 

// root
#include <TSystem.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TRefArray.h>
#include "TTask.h"
#include <TGeoManager.h>
//AliRoot
#include "AliLog.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"
#include "AliESDEvent.h"
#include "AliHeader.h"
#include "AliEMCALGeoUtils.h"
#include "AliEMCALEMCGeometry.h"
#include "AliJetESDFillUnitArrayTracks.h" 
#include "AliJetESDFillUnitArrayEMCalDigits.h"
#include "AliJetUnitArray.h"
#include "AliJetHadronCorrectionv1.h"
#include "AliOADBContainer.h"

ClassImp(AliJetReader)

////////////////////////////////////////////////////////////////////////
AliEMCALGeoUtils *AliJetReader::fGeom=0;

AliJetReader::AliJetReader():
  // Constructor
  fJetanOADBpath(""),
  fChain(0), 
  fTree(0), 
  fMomentumArray(new TClonesArray("TLorentzVector",4000)),
  fArrayMC(0),
  fFillUnitArray(new TTask("fillUnitArray","Fill unit array jet finder")),
  fESD(0),
  fReaderHeader(0),
  fAliHeader(0),
  fSignalFlag(0),
  fCutFlag(0),
  fUnitArray(new TClonesArray("AliJetUnitArray",60000)),
  fArrayInitialised(0),
  fFillUAFromTracks(new AliJetESDFillUnitArrayTracks()), 
  fFillUAFromEMCalDigits(new AliJetESDFillUnitArrayEMCalDigits()),
  fNumCandidate(0),
  fNumCandidateCut(0),
  fHadronCorrector(0),
  fHCorrection(0),
  fECorrection(0),
  fEFlag(kFALSE),
  fDebug(0)
{
  // Default constructor
  fSignalFlag = TArrayI();
  fCutFlag    = TArrayI();
}

////////////////////////////////////////////////////////////////////////

AliJetReader::~AliJetReader()
{
  // Destructor
  if (fMomentumArray) {
      fMomentumArray->Delete();
      delete fMomentumArray;
  }
  
  if (fUnitArray) {
      fUnitArray->Delete();
      delete fUnitArray;
  }
  
  if (fFillUnitArray) {
    delete fFillUnitArray;
  }
  if (fArrayMC) {
      fArrayMC->Delete();
      delete fArrayMC;
  }
  
}

////////////////////////////////////////////////////////////////////////

void AliJetReader::ClearArray()
{
  if (fMomentumArray)  fMomentumArray->Clear();
  if (fFillUnitArray)  fFillUnitArray->Clear();
}

Bool_t AliJetReader::SetEMCALGeometry()
{
  // 
  // Set the EMCal Geometry
  //
  
  fDebug = fReaderHeader->GetDebug();
  
  if(fGeom != 0){
	  Info(" SetEMCALGeometry:","was already done.. it's called just once !!");
     return kTRUE;
     }
     
  if(fDebug>9) cout<<"JetReader: Setting EMCALGeometry"<<endl;

//path to the OADB file


 TString  myPath=  fReaderHeader ->GetMyOADBfile();
 TString OADBfile;

 Bool_t customFile=kFALSE;

if(myPath.Length()) {
	Info(" SetEMCALGeometry","custom version of OADB file: %s",myPath.Data());
	 customFile=kTRUE;
	OADBfile=myPath;
 } else OADBfile.Form("%s/PWG4/JetReconstruction/EMCALlocal2master.root",(const char*) fJetanOADBpath);
 
 AliOADBContainer EMCALgeoCont;
	Info(" SetEMCALGeometry"," I'm going to read the matrices from %s",OADBfile.Data()); 
	TObjArray *EmcalMatrArray;
	if(fDebug>19) cout<<"array definito"<<endl;
	
	 EMCALgeoCont.InitFromFile((char*) OADBfile.Data(),"AliEMCALgeo");
	 EMCALgeoCont.GetDefaultList()->Print(); 
	
	const char*  geoType= fReaderHeader -> GetEMCALgeo2bLoad();
	if(fDebug>19)  cout<<"geometry: "<<geoType<<endl;
	
	const char*  loc2master = ((AliJetESDReaderHeader*) fReaderHeader)->GetEMCALmatrices2bLoad();
	if(fDebug>19) cout<<"matrices: "<<loc2master<<endl;
	
	
     if(fDebug>9)	cout<<"geometry type is: "<<geoType<<endl;
	 if(fDebug>9)   cout<<"survey matrices are: "<<loc2master<<endl;
	
	// some crosschecks to avoid not existing cases
	if(!(!strcmp(geoType, "EMCAL_COMPLETE") || !strcmp(geoType, "EMCAL_COMPLETEV1") || !strcmp(geoType, "EMCAL_FIRSTYEARV1")) ) 
	 Warning(" SetEMCALGeometry","%s is not a known good geometry!  either your are using an old one or everything will crash right now!",geoType);
		
	if(! (!strcmp(loc2master, "survey10") || !strcmp(loc2master, "survey11") || !strcmp(loc2master, "ideal") || !strcmp(loc2master, "test")) ) {
		Warning(" SetEMCALGeometry"," %s is not one of the allowed cases  (listed few lines above) !!----!!!", loc2master);
		loc2master="survey11";
		Warning(" SetEMCALGeometry"," to avoid crashes, your decision has been overrulled!, matrices '%s' will be used instead",loc2master);
		if(fDebug>9)   cout<<"survey matrices are (new, after overrulling): "<<loc2master<<endl;
		}
	
	
	// some warning for not so orthodox combination
	if(!strcmp(geoType, "EMCAL_COMPLETE"))
		 Warning(" SetEMCALGeometry:", "!!----!!  this geometry contains wrong tilt angles for stripmodules.. are you really sure?  please consider using EMCAL_COMPLETEV1 !! ---!! ");
	if(  !strcmp(loc2master, "survey11") && strcmp(geoType, "EMCAL_COMPLETEV1") )
		Warning(" SetEMCALGeometry:",  "!!----!! survey11 matrices should be used with EMCAL_COMPLETEV1 geometry !!---!!");
	if(  !strcmp(loc2master, "survey10") && strcmp(geoType, "EMCAL_FIRSTYEARV1") )
		Warning(" SetEMCALGeometry",  "!!----!! survey10 matrices should be used ONLY with EMCAL_FIRSTYEARV1 geometry!!");
    if(!strcmp(loc2master,"ideal"))
		Warning(" SetEMCALGeometry","!!----!! ideal matrices are without any survey (misalignment correction)... is it really what you want?");
    if(!strcmp(loc2master,"test") && !customFile)
		Warning(" SetEMCALGeometry","!!----!! 'test' matrices will be used. but it seems you didn't provide a custom version of OADB file, the default 'test' is as 'ideal', no survey (misalignment correction) !!----!!");
	if(!strcmp(loc2master,"test") && customFile)
		Info(" SetEMCALGeometry"," !!----!! 'test' matrices read from the custom file you provided     !!----!!");

	EmcalMatrArray=(TObjArray*)EMCALgeoCont.GetObject(100,(char*) loc2master);
	
 
  // Define EMCAL geometry
 
  if(fDebug>10) cout<<"which EMCALgeometry is going to be uploaded?"<<geoType<<endl; 
  fGeom = new AliEMCALGeoUtils(geoType,"EMCAL");
  
  
	for (Int_t mod=0;mod<(fGeom->GetEMCGeometry())->GetNumberOfSuperModules();mod+=1)
 {
      fGeom->SetMisalMatrix(((TGeoHMatrix*) EmcalMatrArray->At(mod)),mod);
     if(fDebug>9)  cout<<"and the matrix is: SM "<<mod<<" (to print the matrix, fDebug>11!) "<<endl;
    if(fDebug>11) { 
    cout<<"print the matrix, (will it work?)"<<endl;
    	((TGeoHMatrix*) EmcalMatrArray->At(mod))->Print();
    	cout<<"if you read that, it did!"<<endl;
    	} 
}
  
  Info("\n SetEMCALGeometry:"," EMCal Geometry set ! \n");
  
  return kTRUE;
}
