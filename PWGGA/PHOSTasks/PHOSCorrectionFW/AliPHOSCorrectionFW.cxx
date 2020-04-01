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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  PHOS tender, recalibrate PHOS clusters                                   //
//  and do track matching                                                    //
//  Author : Dmitri Peressounko (RRC KI)                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMinuit.h"

#include <AliLog.h>
#include <AliVEvent.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliAnalysisManager.h>
#include <AliTender.h>
#include <AliCDBManager.h>
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"

#include "AliPHOSCorrectionFW.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSAodCluster.h"
#include "AliOADBContainer.h"
#include "AliAODCaloCells.h"
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"

ClassImp(AliPHOSCorrectionFW)

AliPHOSCorrectionFW::AliPHOSCorrectionFW():AliAnalysisTaskSE(),
  fEvent(0x0),
  fPHOSGeo(0x0),
  fPHOSCalibData(0x0),
  fDigits(0x0),
  fClusters(0x0),
  fAODClusterArray(0x0),
  fVtx(0.,0.,0.),
  fNonlinearityParams{0},         // Parameters for non-linearity calculation
  fRunByRunCorr{0},               // Per module run-by-run correction
  fL1phase{0},                    // L1phases for PHOS DDLs (run2 only)
  fDigitsUsed{0},                 //Mark digits as already used in cluster (EMC:5*56*64 ; CPV: 5*56*128)
  fZScut(0.025),
  fNoiseMC(0.005),
  fEmcMinE(0.010),
  fClusteringThreshold(0.050),
  fEmcLocMaxCut(0.030),         // Miminal height of local maximum (in GeV)
  fW0(4.5),                   // Weight (def=4.5)
  fEcoreRadius(3.0),          // cire radius (in cm)
  fRunNumber(-1),      // Run number from data
  fDRN(-1),            // Fixed RunNumber if ncecessary
  fRecoPass(1),        // Reco pass for calibration
  fIsMC(false),
  fUsePrivateBadMap(false), 
  fUsePrivateCalib(false),
  fApplyZS(true),             // Apply zero suppression emulation in MC
  fAddNoiseMC(true),          // Smear cell energy by adding electroinc noise emulation in MC
  fToUnfold(true),            // Unfold clusters
  fPrivateOADBBadMap(""),
  fMCProduction(""),
  fNonlinearityVersion("Run2Tune")
{
  for(Int_t i=0; i<6; i++)fPHOSBadMap[i]=0x0 ;  //BadMaps
}
//_____________________________________________________
AliPHOSCorrectionFW::AliPHOSCorrectionFW(const char *name):AliAnalysisTaskSE(name),
  fEvent(0x0),
  fPHOSGeo(0x0),
  fPHOSCalibData(0x0),
  fDigits(0x0),
  fClusters(0x0),
  fAODClusterArray(0x0),
  fVtx(0.,0.,0.),
  fNonlinearityParams{0},         // Parameters for non-linearity calculation
  fRunByRunCorr{0},               // Per module run-by-run correction
  fL1phase{0},                    // L1phases for PHOS DDLs (run2 only)
  fDigitsUsed{0},                 //Mark digits as already used in cluster (EMC:5*56*64 ; CPV: 5*56*128)
  fZScut(0.025),
  fNoiseMC(0.005),
  fEmcMinE(0.010),
  fClusteringThreshold(0.050),
  fEmcLocMaxCut(0.030),         // Miminal height of local maximum (in GeV)
  fW0(4.5),                   // Weight (def=4.5)
  fEcoreRadius(3.0),          // cire radius (in cm)
  fRunNumber(-1),      // Run number from data
  fDRN(-1),            // Fixed RunNumber if ncecessary
  fRecoPass(1),        // Reco pass for calibration
  fIsMC(false),
  fUsePrivateBadMap(false), 
  fUsePrivateCalib(false),
  fApplyZS(true),             // Apply zero suppression emulation in MC
  fAddNoiseMC(true),          // Smear cell energy by adding electroinc noise emulation in MC
  fToUnfold(true),            // Unfold clusters
  fPrivateOADBBadMap(""),
  fMCProduction(""),
  fNonlinearityVersion("Run2Tune")
{
  for(Int_t i=0; i<6; i++)fPHOSBadMap[i]=0x0 ;  //BadMaps	
}
//_____________________________________________________
AliPHOSCorrectionFW::AliPHOSCorrectionFW(const AliPHOSCorrectionFW& ap):AliAnalysisTaskSE(ap),
  fEvent(0x0),
  fPHOSGeo(ap.fPHOSGeo),
  fPHOSCalibData(ap.fPHOSCalibData),
  fDigits(0x0),
  fClusters(0x0),
  fAODClusterArray(0x0),
  fVtx(ap.fVtx),
  fZScut(ap.fZScut),
  fNoiseMC(ap.fNoiseMC),
  fEmcMinE(ap.fEmcMinE),
  fClusteringThreshold(ap.fClusteringThreshold),
  fEmcLocMaxCut(ap.fEmcLocMaxCut),         // Miminal height of local maximum (in GeV)
  fW0(ap.fW0),                   // Weight (def=4.5)
  fEcoreRadius(ap.fEcoreRadius),          // cire radius (in cm)
  fRunNumber(ap.fRunNumber),      // Run number from data
  fDRN(ap.fDRN),            // Fixed RunNumber if ncecessary
  fRecoPass(ap.fRecoPass),        // Reco pass for calibration
  fIsMC(ap.fIsMC),
  fUsePrivateBadMap(ap.fUsePrivateBadMap), 
  fUsePrivateCalib(ap.fUsePrivateCalib),
  fApplyZS(ap.fApplyZS),             // Apply zero suppression emulation in MC
  fAddNoiseMC(ap.fAddNoiseMC),          // Smear cell energy by adding electroinc noise emulation in MC
  fToUnfold(ap.fToUnfold),            // Unfold clusters
  fPrivateOADBBadMap(ap.fPrivateOADBBadMap),
  fMCProduction(ap.fMCProduction),
  fNonlinearityVersion(ap.fNonlinearityVersion)
{
  if(ap.fDigits)fDigits = new TClonesArray("AliPHOSDigit",ap.fDigits->GetSize()) ;
  if(ap.fClusters) fClusters = new TObjArray(ap.fClusters->GetSize()); 
  if(ap.fAODClusterArray) fAODClusterArray= new TClonesArray("AliAODCaloCluster",ap.fAODClusterArray->GetSize()); 
  for(Int_t i=0; i<10; i++)fNonlinearityParams[i]= ap.fNonlinearityParams[i] ;       
  for(Int_t i=0; i<5; i++)fRunByRunCorr[i]=ap.fRunByRunCorr[i];
  for(Int_t i=0; i<15; i++)fL1phase[i]=ap.fL1phase[i] ;
  for(Int_t i=0; i<6; i++)fPHOSBadMap[i]=0x0 ;  //BadMaps	
}   
//_____________________________________________________
AliPHOSCorrectionFW::~AliPHOSCorrectionFW()
{
  if(fDigits) { fDigits->Delete(); delete fDigits; fDigits=0x0;}
  if(fClusters) { fClusters->Delete() ; delete fClusters; fClusters=0x0;}
  if(fAODClusterArray){fAODClusterArray->Delete(); delete fAODClusterArray; fAODClusterArray=0x0;}
}

//_____________________________________________________
void AliPHOSCorrectionFW::Init()
{   
  //
  // Initialise PHOS tender
  //
  if(!InputEvent()) return ;	
printf("InputEvent =%p\n",InputEvent()) ;	
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent()) ;
  if(!fEvent){
  	AliFatal("Can not find AOD event") ;
  }

  fRunNumber = fEvent->GetRunNumber() ;

  if(fIsMC && fDRN > 0){
    fRunNumber = fDRN;//only for single particle simulation
    AliInfo(Form("A dummy run number is set. fRunNumber is %d",fRunNumber));
  }  
   
  //Init fPHOSGeoetry 
  if(!fPHOSGeo){
  	fPHOSGeo =  AliPHOSGeometry::GetInstance() ; //Do not create own geometry if already created
    if(!fPHOSGeo){ //still does not exist
      if(fRunNumber<209122) //Run1
        fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
      else
        fPHOSGeo =  AliPHOSGeometry::GetInstance("Run2") ;      
      AliOADBContainer fPHOSGeoContainer("phosGeo");
      if(fIsMC){ //use excatly the same fPHOSGeoetry as in simulation
        fPHOSGeoContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)fPHOSGeoContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++) {
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
          ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
        } 
      }
      else{ //Use best approaximation to real fPHOSGeoetry
        fPHOSGeoContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
        TObjArray *matrixes = (TObjArray*)fPHOSGeoContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++) {
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
          ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
        }
      } 
    }
  }
  
  //Init Bad channels map
  if(!fUsePrivateBadMap){
    AliOADBContainer badmapContainer(Form("phosBadMap"));
    if(fPrivateOADBBadMap.Length()!=0){
      //Load standard bad maps file if no OADB file is force loaded
      AliInfo(Form("using custom bad channel map from %s\n",fPrivateOADBBadMap.Data()));
       badmapContainer.InitFromFile(fPrivateOADBBadMap.Data(),"phosBadMap");
    } else {
      //Load force loaded OADB file
      AliInfo("using standard bad channel map from $ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root\n");
      badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root","phosBadMap");
    }
    TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRunNumber,"phosBadMap");
    if(!maps){
      AliError(Form("Can not read Bad map for run %d. \n You may choose to use your map with ForceUsingBadMap()\n",fRunNumber)) ;
    }
    else{
      AliInfo(Form("Setting PHOS bad map with name %s \n",maps->GetName())) ;
      for(Int_t mod=0; mod<6;mod++){
        if(fPHOSBadMap[mod]) 
          delete fPHOSBadMap[mod] ;
        TH2I * h = (TH2I*)maps->At(mod) ;
        if(h) fPHOSBadMap[mod]=new TH2I(*h) ;
      }
    }
  }

  if(!fUsePrivateCalib){
    if(fIsMC){ //re/de-calibration for MC productions
      //Init recalibration
      AliOADBContainer calibContainer("phosRecalibration");
      calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCCalibrations.root","phosRecalibration");

      AliInfo(Form("Reading PHOS MC recalibration object for production %s, run=%d", fMCProduction.Data(),fRunNumber)) ;      
      TObjArray *recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration",fMCProduction.Data());
      if(!recalib){
        AliFatal(Form("Can not read calibrations for run %d and name >%s<\n. You may choose your specific calibration with ForceUsingCalibration()\n",fRunNumber,fMCProduction.Data())) ;
      }
      else{
	    //Now try to find object with proper name
	    for(Int_t i=0; i<recalib->GetEntriesFast(); i++){
	      AliPHOSCalibData * tmp = (AliPHOSCalibData*)recalib->At(i) ;
	      if(fMCProduction.CompareTo(tmp->GetName())==0){
            fPHOSCalibData = tmp ;
	        break ;
	      }
	    }
        if(!fPHOSCalibData) {
          AliFatal(Form("Can not find calibration for run %d, and name %s \n",fRunNumber, fMCProduction.Data())) ;
        }
      }
    }
    else{ //real data
      //Init recalibration
      //Check the pass1-pass2-pass3 reconstruction
      AliOADBContainer calibContainer("phosRecalibration");
      calibContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSCalibrations.root","phosRecalibration");
      TObjArray *recalib = (TObjArray*)calibContainer.GetObject(fRunNumber,"PHOSRecalibration");
      if(!recalib){
        AliFatal(Form("Can not read calibrations for run %d\n. You may choose your specific calibration with ForceUsingCalibration()\n",fRunNumber)) ;
      }
      else{
        fPHOSCalibData = (AliPHOSCalibData*)recalib->At(fRecoPass-1) ;
        if(!fPHOSCalibData) {
          AliFatal(Form("Can not find calibration for run %d, pass %d \n",fRunNumber, fRecoPass)) ;
        }
      }
    }
  }
   
  if(!fIsMC){
      //L1phase and run-by-run correction for Run2
      if(fRunNumber>209122){ //Run2
        //L1phase for Run2  
        AliOADBContainer L1Container("phosL1Calibration");
        L1Container.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSL1Calibrations.root","phosL1Calibration");
        TNamed* a= (TNamed*)L1Container.GetObject(fRunNumber);
	  if(!a){
	    AliError(Form("L1phase for run %d was not found, time calibration will be wrong!\n",fRunNumber)) ; 
          for(Int_t ii=0; ii<15; ii++)fL1phase[ii]=0;
	  }
	  else{
          const char*c=a->GetName();
          for(Int_t ii=0; ii<15; ii++)fL1phase[ii]=c[ii]-'0';
	  }
      

	  //Run-by-run correction
        AliOADBContainer runByRunContainer("phosRunByRunCalibration");
        runByRunContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSRunByRunCalibrations.root","phosRunByRunCalibration");
        TNamed* rbr= (TNamed*)runByRunContainer.GetObject(fRunNumber);
	  if(rbr){
          sscanf(rbr->GetName(),"%f,%f,%f,%f",&fRunByRunCorr[1],&fRunByRunCorr[2],&fRunByRunCorr[3],&fRunByRunCorr[4]) ;  
	  }
	  //In any case correction should not be zero
	  //If it is zero, set default and write warning
      for(Int_t mod=1; mod<5; mod++){
          if(fRunByRunCorr[mod]==0.){
              fRunByRunCorr[mod]=0.136 ; 
             AliWarning(Form("Run-by-Run correction for mod. %d is zero in run %d",mod,fRunNumber)) ;  
          }
      }        
    }
  }

  //Non-linearity correction
  if(fNonlinearityVersion==""){ //non-linearity not set by user yet
    if(fRunNumber>209122){ //Run2
      fNonlinearityVersion="Run2" ;
    }else{
      fNonlinearityVersion="Default" ;   
    }
  }


  if(!fDigits) fDigits = new TClonesArray("AliPHOSDigit",100) ;
  if(!fClusters) fClusters = new TObjArray(100) ;


  fAODClusterArray = new TClonesArray("AliAODCaloCluster",200);

}


//_____________________________________________________
void AliPHOSCorrectionFW::UserExec(Option_t * opt)
{
  //Choose PHOS clusters and recalibrate them
  //that it recalculate energy, position and distance 
  //to closest track extrapolation	
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent()) ;
  if(!fEvent){
  	AliFatal("Can not get AOD event") ;
      return ;
  }     
    
  if(!fPHOSCalibData || RunChanged()){ 
    Init();
  }

  fVtx.SetXYZ(fEvent->GetPrimaryVertex()->GetX(), 
  	          fEvent->GetPrimaryVertex()->GetY(),
              fEvent->GetPrimaryVertex()->GetZ()) ;

  //Process calibration
  ConvertDigits() ;

  Custerize() ;

  FillOutput() ;


  const TString brname_cluster = Form("CaloClusters%s",GetName());
  AliInfo(Form("name of new AODCluster branch = %s",brname_cluster.Data()));
  TObject* outClusters = fEvent->FindListObject(brname_cluster);
  if(!outClusters){
    AliInfo(Form("%s object is not found in input event. Add it.",brname_cluster.Data()));
    fAODClusterArray->SetName(brname_cluster);
    fEvent->AddObject(fAODClusterArray);
  }
}

//_____________________________________________________
void AliPHOSCorrectionFW::FillOutput(){
   //Create and fill output branch

  Int_t nclu = fClusters->GetEntriesFast() ;
  if(fAODClusterArray->GetSize()<nclu) fAODClusterArray->Expand(nclu) ;
  fAODClusterArray->Clear() ;
  for (Int_t irp = 0 ; irp < nclu ; irp++) {
    AliPHOSEmcRecPoint  *emcRP = static_cast<AliPHOSEmcRecPoint *>(fClusters->At(irp));
    AliAODCaloCluster   *ec    = new ((*fAODClusterArray)[irp]) AliAODCaloCluster() ; 
    
    TVector3 locpos;
    emcRP->GetLocalPosition(locpos) ;
    TVector3 globalPos ;
    fPHOSGeo->Local2Global(emcRP->GetPHOSMod(), (float)locpos.X(), (float)locpos.Z(), globalPos) ;
    Float_t xyz[3]={(float)globalPos.X(),(float)globalPos.Y(),(float)globalPos.Z()} ;
    ec->SetPosition(xyz);                       //rec.point position in MARS

    Int_t  cellMult   = emcRP->GetDigitsMultiplicity();

    //Primaries
    Int_t  primMult  = 0;
    Int_t *primList =  emcRP->GetPrimaries(primMult);

    Float_t energy = emcRP->GetEnergy();
    energy = CorrectNonlinearity(energy) ;

    ec->SetNCells(cellMult);
    ec->SetType(AliVCluster::kPHOSNeutral);
    ec->SetE(energy);                           //total or core particle energy
    ec->SetDispersion(emcRP->GetDispersion());  //cluster dispersion
    ec->SetM02(emcRP->GetM2x()) ;               //second moment M2x
    ec->SetM20(emcRP->GetM2z()) ;               //second moment M2z
    ec->SetNExMax(emcRP->GetNExMax());          //number of local maxima
    ec->SetTOF(emcRP->GetTime());               //Time of flight - already calibrated in EMCRecPoint
    //Eval FullDispersion
    ec->SetDispersion(TestFullLambda(energy,emcRP->GetM2z(),emcRP->GetM2x())) ;
      // //Eval CoreDispersion
      // Double_t m02=0.,m20=0.;
      // EvalLambdas(&cluPHOS,m02, m20);   
    ec->SetChi2(-1);                     //not yet implemented
//      ec->SetChi2(TestCoreLambda(clu->E(),m20,m02));                     //not yet implemented     

    //Distance to the nearest bad crystal
    ec->SetDistanceToBadChannel(emcRP->GetDistanceToBadCrystal()); 
  
    //Array of MC indeces
    ec->SetLabel(primList,primMult);
    
    Double_t dx,dz,pt;
    Int_t charge ;
    Int_t itr = FindTrackMatching(emcRP->GetPHOSMod(),&locpos, dx, dz, pt,charge) ;
    if(itr>0){
      ec->AddTrackMatched((AliAODTrack*)fEvent->GetTrack(itr));
    }
    ec->SetTrackDistance(dx,dz); 
    Double_t r=TestCPV(dx, dz, pt,charge) ;
    ec->SetEmcCpvDistance(r);    

    }

}


//_____________________________________________________
void AliPHOSCorrectionFW::ConvertDigits()
{
    //Convert PHOSCells to AliPHOSDigits which will be consumed by Clusterizer

    AliAODCaloCells * cells = fEvent->GetPHOSCells() ;

    fDigits->Clear() ;
    Int_t inList=0 ;
    
    Short_t ncell= cells->GetNumberOfCells() ;
    Short_t cellNumber;
    Double_t amplitude=0., time=0., efrac=0.;
    Int_t mclabel;
    if(fDigits->GetSize()<ncell){
      fDigits->Expand(ncell);
    }
    for(Short_t pos=0; pos<ncell; pos++){
	  cells->GetCell(pos, cellNumber, amplitude, time, mclabel, efrac) ;
	  if(cellNumber<0){ //this is CPV cells
	  	continue ;
	  }

      //Calibrate Energy and time
   	  Bool_t isHG=cells->GetHighGain(pos) ;
      amplitude=Calibrate(amplitude,cellNumber) ;
      time = CalibrateT(time,cellNumber,isHG) ;

      //Select digits'
      if(fApplyZS && amplitude<fZScut)
        continue ;
      //..............

      //TestBadMap
      Int_t relid[4] ;
      fPHOSGeo->AbsToRelNumbering(Int_t(cellNumber), relid) ; 
      if(!fPHOSBadMap[relid[0]]){
          AliError(Form("No Bad map for PHOS module %d",relid[0])) ;
          continue ;
      }
      if(fPHOSBadMap[relid[0]]->GetBinContent(relid[2],relid[3])>0){  //Bad channel
  	     continue ;
      }

      //If some modifications if necessary
      if(fAddNoiseMC){
         amplitude=TMath::Max(0.,amplitude+gRandom->Gaus(0,fNoiseMC)) ;
      }

      if(amplitude<fEmcMinE)
      	continue ;
      
      new((*fDigits)[inList]) AliPHOSDigit(mclabel, (int)cellNumber, (float)amplitude, (float)time,inList) ;  //Digits keep calibrated energy and time
      inList++ ;
    }
}
//_____________________________________________________
void AliPHOSCorrectionFW::Custerize(){
   //Repeat clusterization from collected digits


  Int_t nClusters = 0 ;

  //Mark all digits as unused yet
  Int_t nDigits=fDigits->GetEntriesFast() ;

  for(Int_t i=0; i<nDigits; i++){
    fDigitsUsed[i]=0 ;
  }
  Int_t iFirst = 0 ; //first index of digit which potentially can be a part of cluster
                     //e.g. first digit in this module, first CPV digit etc.
  AliPHOSDigit * digit ; 
  TArrayI clusterdigitslist(nDigits) ;   
  AliPHOSRecPoint * clu = 0 ; 
  for(Int_t i=0; i<nDigits; i++){
    if(fDigitsUsed[i])
      continue ;

    digit=static_cast<AliPHOSDigit*>(fDigits->At(i)) ;

    clu=0 ;

    Int_t index ;

    //is this digit so energetic that start cluster?
    if ( digit->GetEnergy() > fClusteringThreshold ) {
        Int_t iDigitInCluster = 0 ; 
        // start a new EMC RecPoint
        if(nClusters >= fClusters->GetSize()) 
          fClusters->Expand(2*nClusters+1) ;
          
        fClusters->AddAt(new  AliPHOSEmcRecPoint(""), nClusters) ;
        clu = static_cast<AliPHOSEmcRecPoint *>( fClusters->At(nClusters) ) ; 
	    nClusters++ ; 
	    clu->AddDigit(*digit, digit->GetEnergy(),digit->GetTime()) ;
        clusterdigitslist[iDigitInCluster] = digit->GetIndexInList() ;
        iDigitInCluster++ ;
        fDigitsUsed[i]=kTRUE ; 
             

        //Now scan remaining digits in list to find neigbours of our seed
        AliPHOSDigit * digitN ; 
        index = 0 ;
        while (index < iDigitInCluster){ // scan over digits already in cluster 
          digit =  static_cast<AliPHOSDigit*>( fDigits->At(clusterdigitslist[index]) )  ;      
          index++ ; 
          for(Int_t j=iFirst; j<nDigits; j++){
            if(fDigitsUsed[j]) 
              continue ;        //look through remaining digits
            digitN = static_cast<AliPHOSDigit*>( fDigits->At(j) ) ;
            Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
            switch (ineb ) {
            case -1:   //too early (e.g. previous module), do not look before j at subsequent passes
              iFirst=j ;
              break ;
            case 0 :   // not a neighbour
              break ;
            case 1 :   // are neighbours 
	          clu->AddDigit(*digitN, digitN->GetEnergy(),digitN->GetTime()) ;
              clusterdigitslist[iDigitInCluster] = j ; 
              iDigitInCluster++ ; 
              fDigitsUsed[j]=kTRUE ;
              break ;
            case 2 :   // too far from each other
              goto endOfLoop;   
            } // switch
          
          }
          endOfLoop: ; //scanned all possible neighbours for this digit
      } // loop over cluster  
    }   
  }

  //
  //UNFOLD clusters
  //
  if(fToUnfold && (nClusters > 0)){ 

    Int_t numberofNotUnfolded = nClusters ; 
    Int_t index ;   
    for(index = 0 ; index < numberofNotUnfolded ; index++){
      
      AliPHOSEmcRecPoint * emcRecPoint = static_cast<AliPHOSEmcRecPoint *>( fClusters->At(index) ) ;
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      AliPHOSDigit ** maxAt = new AliPHOSDigit*[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy,fEmcLocMaxCut,fDigits) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
        UnfoldCluster(emcRecPoint, nMax, maxAt, maxAtEnergy) ;

        fClusters->Remove(emcRecPoint); 
        fClusters->Compress() ;
        index-- ;
        nClusters -- ;
        numberofNotUnfolded-- ;
      }
      else{
        emcRecPoint->SetNExMax(1) ; //Only one local maximum
      }
      
      delete[] maxAt ; 
      delete[] maxAtEnergy ; 
    }
  } 
  // Unfolding of EMC clusters finished


  //Calculate cluster parameters
  //Evaluate position, dispersion and other RecPoint properties..
  for(Int_t index = 0; index < nClusters; index++){
    AliPHOSEmcRecPoint * rp = static_cast<AliPHOSEmcRecPoint *>( fClusters->At(index) );
    rp->Purify(fEmcMinE,fDigits) ;
    if(rp->GetMultiplicity()==0){
      fClusters->RemoveAt(index) ;
      delete rp ;
      continue;
    }

    rp->EvalAll(fDigits) ;
    rp->EvalCoreEnergy(fW0,fEcoreRadius,fDigits) ;
    rp->EvalAll(fW0,fVtx,fDigits) ;
    rp->SetDistanceToBadCrystal(DistanceToBadChannels(rp));
  }
  fClusters->Compress() ;
  fClusters->Sort() ; 

}
//____________________________________________________________________________
void  AliPHOSCorrectionFW::UnfoldCluster(AliPHOSEmcRecPoint * iniEmc, 
                                          Int_t nMax, 
                                          AliPHOSDigit ** maxAt, 
                                          Float_t * maxAtEnergy)
{ 
  // Performs the unfolding of a cluster with nMax overlapping showers 

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;

  Bool_t rv = FindFit(iniEmc, maxAt, maxAtEnergy, nPar, fitparameters) ;

  if( !rv ) {
    // Fit failed, return and remove cluster
    iniEmc->SetNExMax(-1) ;
    delete[] fitparameters ; 
    return ;
  }

  // create ufolded rec points and fill them with new energy lists
  // First calculate energy deposited in each sell in accordance with fit (without fluctuations): efit[]
  // and later correct this number in acordance with actual energy deposition

  Int_t nDigits = iniEmc->GetMultiplicity() ;  
  Float_t * efit = new Float_t[nDigits] ;
  Float_t xDigit=0.,zDigit=0. ;
  Float_t xpar=0.,zpar=0.,epar=0.  ;
  Int_t relid[4] ;
  AliPHOSDigit * digit = 0 ;
  Int_t * emcDigits = iniEmc->GetDigitsList() ;

  TVector3 vIncid ;

  Int_t iparam ;
  Int_t iDigit ;
  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = static_cast<AliPHOSDigit*>( fDigits->At(emcDigits[iDigit] ) ) ;   
    fPHOSGeo->AbsToRelNumbering(digit->GetId(), relid) ;
    fPHOSGeo->RelPosInModule(relid, xDigit, zDigit) ;
    efit[iDigit] = 0;

    iparam = 0 ;    
    while(iparam < nPar ){
      xpar = fitparameters[iparam] ;
      zpar = fitparameters[iparam+1] ;
      epar = fitparameters[iparam+2] ;
      iparam += 3 ;
//      fPHOSGeo->GetIncidentVector(fVtx,relid[0],xpar,zpar,vIncid) ;
//      efit[iDigit] += epar * ShowerShape(xDigit - xpar,zDigit - zpar,vIncid) ;
      efit[iDigit] += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
    }
  }
  
  // Now create new RecPoints and fill energy lists with efit corrected to fluctuations
  // so that energy deposited in each cell is distributed betwin new clusters proportionally
  // to its contribution to efit

  Float_t * emcEnergies = iniEmc->GetEnergiesList() ;
  Float_t ratio ;
  Int_t numberOfEmcClusters= fClusters->GetEntriesFast() ;
  iparam = 0 ;
  while(iparam < nPar ){
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    iparam += 3 ;    

    AliPHOSEmcRecPoint * emcRP = 0 ;  

      
    if(numberOfEmcClusters >= fClusters->GetSize())
        fClusters->Expand(2*numberOfEmcClusters) ;
      
    (*fClusters)[numberOfEmcClusters] = new AliPHOSEmcRecPoint("") ;
    emcRP = static_cast<AliPHOSEmcRecPoint *>( fClusters->At(numberOfEmcClusters) ) ;
    numberOfEmcClusters++ ;
    emcRP->SetNExMax((Int_t)nPar/3) ;
    
    Float_t eDigit ;
    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = static_cast<AliPHOSDigit*>( fDigits->At( emcDigits[iDigit] ) ) ; 
      fPHOSGeo->AbsToRelNumbering(digit->GetId(), relid) ;
      fPHOSGeo->RelPosInModule(relid, xDigit, zDigit) ;
      if(efit[iDigit]>0){
//      ratio = epar * ShowerShape(xDigit - xpar,zDigit - zpar,vIncid) / efit[iDigit] ; 
        ratio = epar * ShowerShape(xDigit - xpar,zDigit - zpar) / efit[iDigit] ; 
        eDigit = emcEnergies[iDigit] * ratio ;
        emcRP->AddDigit( *digit, eDigit,digit->GetTime() ) ;
      }
    } 
  }
 
  delete[] fitparameters ; 
  delete[] efit ; 
  
}

//_____________________________________________________________________________
void AliPHOSCorrectionFW::UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TList * toMinuit = static_cast<TList*>( gMinuit->GetObjectFit() ) ;

  AliPHOSEmcRecPoint * emcRP = static_cast<AliPHOSEmcRecPoint*>( toMinuit->At(0) )  ;
  TClonesArray * digits = static_cast<TClonesArray*>( toMinuit->At(1) )  ;
  // A bit buggy way to get an access to the fPHOSGeoetry
  // To be revised!
  AliPHOSGeometry *fPHOSGeo = static_cast<AliPHOSGeometry *>(toMinuit->At(2));

//  TVector3 * vtx = static_cast<TVector3*>(toMinuit->At(3)) ;  //Vertex position
  
  //  AliPHOSEmcRecPoint * emcRP = static_cast<AliPHOSEmcRecPoint *>( gMinuit->GetObjectFit() ) ; // EmcRecPoint to fit

  Int_t * emcDigits     = emcRP->GetDigitsList() ;

  Int_t nOdigits = emcRP->GetDigitsMultiplicity() ; 

  Float_t * emcEnergies = emcRP->GetEnergiesList() ;
  
//  TVector3 vInc ;
  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t efit ;    

  AliPHOSDigit * digit ;
  Int_t iDigit ;

  for( iDigit = 0 ; iDigit < nOdigits ; iDigit++) {

    digit = static_cast<AliPHOSDigit*>( digits->At( emcDigits[iDigit] ) ); 

    Int_t relid[4] ;
    Float_t xDigit ;
    Float_t zDigit ;

    fPHOSGeo->AbsToRelNumbering(digit->GetId(), relid) ;

    fPHOSGeo->RelPosInModule(relid, xDigit, zDigit) ;

     if(iflag == 2){  // calculate gradient
       Int_t iParam = 0 ;
       efit = 0 ;
       while(iParam < nPar ){
         Double_t dx = (xDigit - x[iParam]) ;
         iParam++ ; 
         Double_t dz = (zDigit - x[iParam]) ; 
         iParam++ ;          
//         fPHOSGeo->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),x[iParam-2],x[iParam-1],vInc) ;
//         efit += x[iParam] * ShowerShape(dx,dz,vInc) ;
         efit += x[iParam] * ShowerShape(dx,dz) ;
         iParam++ ;
       }
       Double_t sum = 2. * (efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 
       iParam = 0 ;
       while(iParam < nPar ){
         Double_t xpar = x[iParam] ;
         Double_t zpar = x[iParam+1] ;
         Double_t epar = x[iParam+2] ;
//         fPHOSGeo->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),xpar,zpar,vInc) ;
         Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
//         Double_t shape = sum * ShowerShape(xDigit - xpar,zDigit - zpar,vInc) ;
         Double_t shape = sum * ShowerShape(xDigit - xpar,zDigit - zpar) ;
//DP: No incident angle dependence in gradient yet!!!!!!
         Double_t r4 = dr*dr*dr*dr ;
         Double_t r295 = TMath::Power(dr,2.95) ;
         Double_t deriv =-4. * dr*dr * ( 2.32 / ( (2.32 + 0.26 * r4) * (2.32 + 0.26 * r4) ) +
                                         0.0316 * (1. + 0.0171 * r295) / ( ( 1. + 0.0652 * r295) * (1. + 0.0652 * r295) ) ) ;
         
         Grad[iParam] += epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x    
         iParam++ ; 
         Grad[iParam] += epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z         
         iParam++ ; 
         Grad[iParam] += shape ;                                  // Derivative over energy             
         iParam++ ; 
       }
     }
     efit = 0;
     iparam = 0 ;

     while(iparam < nPar ){
       Double_t xpar = x[iparam] ;
       Double_t zpar = x[iparam+1] ;
       Double_t epar = x[iparam+2] ;
       iparam += 3 ;
//       fPHOSGeo->GetIncidentVector(*vtx,emcRP->GetPHOSMod(),xpar,zpar,vInc) ;
//       efit += epar * ShowerShape(xDigit - xpar,zDigit - zpar,vInc) ;
       efit += epar * ShowerShape(xDigit - xpar,zDigit - zpar) ;
     }

     fret += (efit-emcEnergies[iDigit])*(efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 
     // Here we assume, that sigma = sqrt(E)
  }

}
//____________________________________________________________________________
Double_t AliPHOSCorrectionFW::DistanceToBadChannels(AliPHOSEmcRecPoint * rp)
{
  //For each EMC rec. point set the distance to the nearest bad crystal.
  
  //Calculate relId position of cluster center
  Int_t relid[4]={0,0,0,0} ;
  TVector3 lpos ;
  rp->GetLocalPosition(lpos) ;
  fPHOSGeo->RelPosToRelId(rp->GetPHOSMod(), lpos.X(), lpos.Z(),relid) ; 

  //Test distances
  if(fPHOSBadMap[relid[0]]->GetBinContent(relid[2],relid[3])) return 0.; //in bad channel

  Double_t sum=0.;
  Int_t ixMin = TMath::Max(1,relid[2]-1) ;
  Int_t ixMax = TMath::Min(64,relid[2]+1) ;
  Int_t izMin = TMath::Max(1,relid[3]-1) ;
  Int_t izMax = TMath::Min(56,relid[3]+1) ;
  for(Int_t ix=ixMin; ix<=ixMax; ix++){  
    for(Int_t iz=izMin; iz<=izMax; iz++){
       sum += fPHOSBadMap[relid[0]]->GetBinContent(ix,iz) ;
    }
  }
  if(sum>0) return 1.; 

  // 2 cells
  ixMin = TMath::Max(1,relid[2]-2) ;
  ixMax = TMath::Min(64,relid[2]+2) ;
  izMin = TMath::Max(1,relid[3]-2) ;
  izMax = TMath::Min(56,relid[3]+2) ;
  sum=0.;
  for(Int_t ix=ixMin; ix<=ixMax; ix++){ 
       sum += fPHOSBadMap[relid[0]]->GetBinContent(ix,izMin) + 
              fPHOSBadMap[relid[0]]->GetBinContent(ix,izMax) ;
  }
  for(Int_t iz=izMin; iz<=izMax; iz++){ 
       sum += fPHOSBadMap[relid[0]]->GetBinContent(ixMin,iz) + 
              fPHOSBadMap[relid[0]]->GetBinContent(ixMax,iz) ;
  }
  if(sum>0) return 2.; 
  else return 3. ; //More than 3 cells distance

}
//____________________________________________________________________________
Double_t  AliPHOSCorrectionFW::ShowerShape(Double_t x, Double_t z)
{ 
  // Shape of the shower (see PHOS TDR)
  // If you change this function, change also the gradient evaluation in ChiSquare()

  //for the moment we neglect dependence on the incident angle.  

  Double_t r2    = x*x + z*z ;
  Double_t r4    = r2*r2 ;
  Double_t r295  = TMath::Power(r2, 2.95/2.) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//_____________________________________________________
Int_t AliPHOSCorrectionFW::AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)const
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  //                                       =-1 are not neighbour, continue searching, but do not look before d2 next time
  // neighbours are defined as digits having at least a common vertex 
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

  Int_t relid1[4] ; 
  fPHOSGeo->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  fPHOSGeo->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module 
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){   //At least common vertex
      //    if (( relid1[2]==relid2[2] && coldiff <= 1 )  || ( relid1[3]==relid2[3] &&  rowdiff <= 1 )){ //common side
      if(relid1[1] != 0)  
      return 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
        return 2; //  Difference in row numbers is too large to look further 
    }
    return 0 ;

  } 
  else {
    if(relid1[0] > relid2[0] && relid1[1]==relid2[1] ) //we switched to the next module
      return -1 ;
    if(relid1[1] < relid2[1]) //we switched from EMC(0) to CPV(-1)
      return -1 ;
    
    return 2 ;

  }

  return 0 ; 
}
//________________________________________________________________________
Double_t AliPHOSCorrectionFW::Calibrate(Double_t e, Int_t absId){
  //Apply energy re-calibration 

  Int_t relId[4];
  fPHOSGeo->AbsToRelNumbering(absId,relId) ;
  Int_t   module = relId[0];
  Int_t   column = relId[3];
  Int_t   row    = relId[2];
  e*=fPHOSCalibData->GetADCchannelEmc(module, column, row);
  return e ;  
}


//________________________________________________________________________
Double_t AliPHOSCorrectionFW::CalibrateT(Double_t tof, Int_t absId, Bool_t isHG){
  //Apply time re-calibration separately for HG and LG channels
  //By default (if not filled) shifts are zero.  
    
  Int_t relId[4];
  fPHOSGeo->AbsToRelNumbering(absId,relId) ;
  Int_t   module = relId[0];
  Int_t   column = relId[3];
  Int_t   row    = relId[2];
  if(isHG)
    tof-=fPHOSCalibData->GetTimeShiftEmc(module, column, row);
  else{
    tof-=fPHOSCalibData->GetLGTimeShiftEmc(module, column, row);
  }
  if(!fIsMC){
    //Apply L1phase
    //First eval DDL
    const Int_t nmod=5; 
    Int_t ddl = (nmod-module) * 4 + (row-1)/16 - 6; //convert offline module numbering to online.
    //L1phase is 0 for Run1
    if(fRunNumber>209122){ //Run2
      UShort_t BC = fEvent->GetBunchCrossNumber();
      Int_t timeshift = BC%4 - fL1phase[ddl];
      if(timeshift<0) timeshift += 4; 
      tof -= timeshift*25e-9;
    }
  }
  return tof ;
}
//________________________________________________________________________
Bool_t AliPHOSCorrectionFW::RunChanged(){
   return fRunNumber != InputEvent()->GetRunNumber() ;
}
//____________________________________________________________
Double_t AliPHOSCorrectionFW::CorrectNonlinearity(Double_t en){

  //For backward compatibility, if no RecoParameters found
  if(fNonlinearityVersion=="Default"){
    return 0.0241+1.0504*en+0.000249*en*en ;
  }
  if(fNonlinearityVersion=="MC"){ //Default + some correction
    return (0.0241+1.0504*en+0.000249*en*en)*fNonlinearityParams[0]*(1+fNonlinearityParams[1]/(1.+en*en/fNonlinearityParams[2]/fNonlinearityParams[2])) ;
  }

  if(fNonlinearityVersion=="NoCorrection"){
    return en ;
  }
  if(fNonlinearityVersion=="Gustavo2005"){
    return fNonlinearityParams[0]+fNonlinearityParams[1]*en + fNonlinearityParams[2]*en*en ;
  }
  if(fNonlinearityVersion=="Henrik2010"){
    return en*(fNonlinearityParams[0]+fNonlinearityParams[1]*TMath::Exp(-en*fNonlinearityParams[2]))*(1.+fNonlinearityParams[3]*TMath::Exp(-en*fNonlinearityParams[4]))*(1.+fNonlinearityParams[6]/(en*en+fNonlinearityParams[5])) ;
  }
  if(fNonlinearityVersion=="Run2"){
     return (1.-0.08/(1.+en*en/0.055))*(0.03+6.65e-02*TMath::Sqrt(en)+en) ; 
  }
  if(fNonlinearityVersion=="Run2MC"){ //Default for Run2 + some correction
    return (1.-0.08/(1.+en*en/0.055))*(0.03+6.65e-02*TMath::Sqrt(en)+en)*fNonlinearityParams[0]*(1+fNonlinearityParams[1]/(1.+en*en/fNonlinearityParams[2]/fNonlinearityParams[2])) ;
  }
  
  if(fNonlinearityVersion=="Run2Tune"){ //Improved Run2 tune for Emin extended down to 100 MeV
    if(en<=0.) return 0.;
    const Double_t xMin=0.36; //low part of the param (optimized from pi0 peak)
    const Double_t xMax=5.17; //Upper part of the param (optimized from pi0 peak)
    
    //middle part param    
    const Double_t a= 1.02165   ; 
    const Double_t b=-2.548e-01 ; 
    const Double_t c= 6.483e-01 ;     
    const Double_t d=-0.4805    ;
    const Double_t e= 0.1275    ;

    Double_t ecorr=0.;
    if(en<xMin){
       const Double_t beta = 2.*a*sqrt(xMin)+b-d/(xMin)-2.*e/(xMin*sqrt(xMin));
       const Double_t alpha = a*xMin+b*sqrt(xMin)+c+d/sqrt(xMin)+e/xMin-beta*sqrt(xMin) ;
       ecorr= 1.0312526*(alpha+beta*sqrt(en)) ;  
    }
    else{
      if(en<xMax){
         ecorr= 1.0312526*(a*en+b*sqrt(en)+c+d/sqrt(en)+e/en) ;
      }
      else{
        const Double_t beta= b+2.*c/sqrt(xMax)+3.*d/xMax+4.*e/xMax/sqrt(xMax) ;
        const Double_t alpha = a+b/sqrt(xMax)+c/xMax+d/xMax/sqrt(xMax)+e/(xMax*xMax)-beta/sqrt(xMax) ;
        ecorr= 1.0312526*(alpha*en+beta*sqrt(en)) ;  
      }
    }
    if(ecorr>0){
      return ecorr;
    }
    else{
      return 0.;
    }
  }
  if(fNonlinearityVersion=="Run2TuneMC"){ //Improved Run2 tune for MC
    if(en<=0.) return 0.;

    const Double_t p0 = 1.031;
    const Double_t p1 = 0.51786058;
    const Double_t p2 = 0.13504396;
    const Double_t p3 = -0.14737537;
    const Double_t p4 = -0.455062;

    const Double_t Nonlin = p0+p1/en+p2/en/en+p3/TMath::Sqrt(en)+p4/en/TMath::Sqrt(en);

    return en * Nonlin;
  }
  if(fNonlinearityVersion=="Run2TuneMCNoNcell"){ //Improved Run2 tune for MC in the case of loose cluster cuts (no Ncell>2 cut)
    if(en<=0.) return 0.;
    
    const Double_t xMin=0.850;  //low part of the param (optimized from pi0 peak)
    const Double_t xMax=5.17;   //Upper part of the param (optimized from pi0 peak)
        
    //middle part param    
    const Double_t a= 1.02165   ; 
    const Double_t b=-2.548e-01 ; 
    const Double_t c= 0.6483 ;
    const Double_t d=-0.4980 ;
    const Double_t e= 0.1245 ;  
    
    Double_t ecorr=0.;
    if(en<xMin){
       const Double_t gamma = 0.150 ;
       const Double_t beta = 0.5*(0.5*b*sqrt(xMin)+c+1.5*d/sqrt(xMin)+2.*e/xMin)*TMath::Power((xMin*xMin+gamma*gamma),2)/
       (xMin*xMin*xMin);
       const Double_t alpha = (a*xMin+b*sqrt(xMin)+c+d/sqrt(xMin)+e/xMin-beta*xMin/(xMin*xMin+gamma*gamma))/xMin ;
       ecorr= 1.0328783*(alpha*en+beta*en/(en*en+gamma*gamma)) ;  
    }
    else{
      if(en<xMax){
         ecorr= 1.0328783*(a*en+b*sqrt(en)+c+d/sqrt(en)+e/en) ;
      }
      else{
        const Double_t beta= b+2.*c/sqrt(xMax)+3.*d/xMax+4.*e/xMax/sqrt(xMax) ;
        const Double_t alpha = a+b/sqrt(xMax)+c/xMax+d/xMax/sqrt(xMax)+e/(xMax*xMax)-beta/sqrt(xMax) ;
        ecorr= 1.0328783*(alpha*en+beta*sqrt(en)) ;  
      }
    }
    if(ecorr<0){
      return 0.;
    }

    return ecorr ;    
  }

  return en ;
}

//___________________________________________________________________________________________________
Int_t AliPHOSCorrectionFW::FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge){
  //Find track with closest extrapolation to cluster
  Double_t magF = fEvent->GetMagneticField();
 
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default") ;
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = fEvent->GetNumberOfTracks();
      
  //Calculate actual distance to PHOS module
  TVector3 globaPos ;
  fPHOSGeo->Local2Global(mod, 0.,0., globaPos) ;
  const Double_t rPHOS = globaPos.Pt() ; //Distance to center of  PHOS module
  const Double_t kYmax = 72.+10. ; //Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64.+10. ; //Size of the module (with some reserve) in z direction
  const Double_t kAlpha0=330./180.*TMath::Pi() ; //First PHOS module angular direction
  const Double_t kAlpha= 20./180.*TMath::Pi() ; //PHOS module angular size
  Double_t minDistance = 1.e6;


  Double_t gposTrack[3] ; 

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5*kAlmost0Field,bz) + bz;

  Double_t b[3]; 
  Int_t itr=-1 ;
  AliAODTrack *aodTrack=0x0 ;
  Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
  for (Int_t i=0; i<nt; i++) {
      aodTrack=(AliAODTrack*)fEvent->GetTrack(i);
      
      //Continue extrapolation from TPC outer surface
      AliExternalTrackParam outerParam;
      aodTrack->GetPxPyPz(pxpypz);
      aodTrack->GetXYZ(xyz);
      aodTrack->GetCovarianceXYZPxPyPz(cv);
      outerParam.Set(xyz,pxpypz,cv,aodTrack->Charge());
      
      Double_t z; 
      if(!outerParam.GetZAt(rPHOS,bz,z))
        continue ;

      if (TMath::Abs(z) > kZmax) 
        continue; // Some tracks miss the PHOS in Z

     
      //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*mod ;
      if(!outerParam.RotateParamOnly(phiMod)) continue ; //RS use faster rotation if errors are not needed 
    
      Double_t y;                       // Some tracks do not reach the PHOS
      if (!outerParam.GetYAt(rPHOS,bz,y)) continue; //    because of the bending
      
      if(TMath::Abs(y) < kYmax){
        outerParam.GetBxByBz(b) ;
        outerParam.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
        //outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
        outerParam.GetXYZ(gposTrack) ;
        TVector3 globalPositionTr(gposTrack) ;
        TVector3 localPositionTr ;
        fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,mod) ;
        Double_t ddx = locpos->X()-localPositionTr.X();
        Double_t ddz = locpos->Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
	      dx = ddx ;
  	      dz = ddz ;
	      minDistance=d2 ;
	      itr=i ;
          pt=aodTrack->Pt() ;
          charge=aodTrack->Charge() ;
        }
      }
    }//Scanned all tracks
 
   return itr ;
}

//_____________________________________________________________________________
Double_t AliPHOSCorrectionFW::TestCoreLambda(Double_t pt,Double_t l1,Double_t l2){
//Parameterization for core dispersion   
//For R=4.5
  Double_t   l1Mean  = 1.150200 + 0.097886/(1.+1.486645*pt+0.000038*pt*pt) ;
  Double_t   l2Mean = 1.574706 + 0.997966*exp(-0.895075*pt)-0.010666*pt ;
  Double_t   l1Sigma = 0.100255 + 0.337177*exp(-0.517684*pt)+0.001170*pt ;
  Double_t   l2Sigma = 0.232580 + 0.573401*exp(-0.735903*pt)-0.002325*pt ;
  Double_t   c = -0.110983 -0.017353/(1.-1.836995*pt+0.934517*pt*pt) ;

  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return R2 ;
  
}
//_____________________________________________________________________________
Double_t AliPHOSCorrectionFW::TestFullLambda(Double_t pt,Double_t l1,Double_t l2){
//Parameterization for full dispersion   
  //Parameterizatino for full dispersion
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;

  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return R2 ;
  
}
//____________________________________________________________________________
Double_t AliPHOSCorrectionFW::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Return distance to closest track in sigmas
   
   Double_t mf = fEvent->GetMagneticField(); //Positive for ++ and negative for --
  
   Double_t meanX=0;
   Double_t meanZ=0.;
   Double_t sx=0.; 
   Double_t sz=0.; 
  if(fRunNumber<209122){ //Run1
    sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
    sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  
    if(mf<0.){ //field --
      meanZ = -0.468318 ;
      if(charge>0)
        meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
      else
        meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
    }
    else{ //Field ++
      meanZ= -0.468318;
      if(charge>0)
        meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
      else
        meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
    }

  }
  else{//Run2
  
    sx = TMath::Min(5.2, 1.111 + 0.56 * TMath::Exp(-0.031 * pt*pt) + 4.8 /TMath::Power(pt+0.61,3));
    sz = TMath::Min(3.3, 1.12  + 0.35 * TMath::Exp(-0.032 * pt*pt) + 0.75/TMath::Power(pt+0.24,3));

    if(mf<0.){ //field --
      meanZ = 0.102;
      if(charge>0)
        meanX =  TMath::Min(5.8, 0.42 + 0.70 * TMath::Exp(-0.015 * pt*pt) + 35.8/TMath::Power(pt+1.41,3));
      else
        meanX = -TMath::Min(5.8, 0.17 + 0.64 * TMath::Exp(-0.019 * pt*pt) + 26.1/TMath::Power(pt+1.21,3));
    }
    else{ //Field ++
      meanZ= 0.102;
      if(charge>0)
        meanX = -TMath::Min(5.8, 0.58 + 0.68 * TMath::Exp(-0.027 * pt*pt) + 28.0/TMath::Power(pt+1.28,3));
      else
        meanX =  TMath::Min(5.8, 0.11 + 0.67 * TMath::Exp(-0.015 * pt*pt) + 29.9/TMath::Power(pt+1.29,3));
    }

  }
  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
  
}
