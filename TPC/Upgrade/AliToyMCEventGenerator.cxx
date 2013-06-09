#include <iostream>

#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TH2F.h>
#include <TGeoGlobalMagField.h>

#include <AliLog.h>
#include <AliTPCROC.h>
#include <AliTrackPointArray.h>
#include <AliTrackerBase.h>
#include <AliCDBManager.h>
#include <AliTPCParam.h>
#include <AliGeomManager.h>
#include <AliTPCcalibDB.h>
#include <AliTPCclusterMI.h>
#include <AliTPCSpaceCharge3D.h>

#include "AliToyMCEvent.h"
#include "AliToyMCTrack.h"

#include "AliToyMCEventGenerator.h"

ClassImp(AliToyMCEventGenerator);


AliToyMCEventGenerator::AliToyMCEventGenerator()
  :TObject()
  ,fTPCParam(0x0)
  ,fEvent(0x0)
  ,fSpaceCharge(0x0)
  ,fOutputFileName("toyMC.root")
  ,fOutFile(0x0)
  ,fOutTree(0x0)
{
  //TODO: Add method to set custom space charge file
  fSpaceCharge = new AliTPCSpaceCharge3D();
  fSpaceCharge->SetSCDataFileName("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz.root");
  fSpaceCharge->SetOmegaTauT1T2(0.325,1,1); // Ne CO2
  //fSpaceCharge->SetOmegaTauT1T2(0.41,1,1.05); // Ar CO2
  fSpaceCharge->InitSpaceCharge3DDistortion();
//   fSpaceCharge->CreateHistoSCinZR(0.,50,50)->Draw("surf1");
//   fSpaceCharge->CreateHistoDRPhiinZR(0,100,100)->Draw("colz");
  //!!! This should be handled by the CongiOCDB macro
  // const char* ocdb="local://$ALICE_ROOT/OCDB/";
  // AliCDBManager::Instance()->SetDefaultStorage(ocdb);
  // AliCDBManager::Instance()->SetRun(0);   
  // TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
  // AliGeomManager::LoadGeometry("");
  fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
  //std::cout<<"----------->"<<fTPCParam<<std::endl;
  fTPCParam->ReadGeoMatrices();
  //std::cout<<"----------->"<<fTPCParam<<std::endl;

}
//________________________________________________________________
AliToyMCEventGenerator::AliToyMCEventGenerator(const AliToyMCEventGenerator &gen)
  :TObject(gen)
  ,fTPCParam(gen.fTPCParam)
  ,fEvent(0x0)
  ,fSpaceCharge(gen.fSpaceCharge)
  ,fOutputFileName(gen.fOutputFileName)
  ,fOutFile(0x0)
  ,fOutTree(0x0)
{
  //
}
//________________________________________________________________
AliToyMCEventGenerator::~AliToyMCEventGenerator() 
{
  delete fSpaceCharge;
  //!!! Is fTPCParam not owned by the CDB manager?
  delete fTPCParam;
}
//________________________________________________________________
Bool_t AliToyMCEventGenerator::DistortTrack(AliToyMCTrack &trackIn, Double_t t0)
{
  //
  //
  //

  //!!! Should this complete function not be moved to AliToyMCTrack, and directly called from the constructor
  //    or setter of the track parameters?

  if(!fTPCParam) {
    fTPCParam = AliTPCcalibDB::Instance()->GetParameters();
    fTPCParam->ReadGeoMatrices();
    std::cout << "init tpc param" << std::endl;
  }
  const Int_t kNRows=AliTPCROC::Instance()->GetNRows(0)+AliTPCROC::Instance()->GetNRows(36);
  // const Double_t kRTPC0  =AliTPCROC::Instance()->GetPadRowRadii(0,0);
  // const Double_t kRTPC1  =AliTPCROC::Instance()->GetPadRowRadii(36,AliTPCROC::Instance()->GetNRows(36)-1);
  const Double_t kMaxSnp = 0.85;  
  const Double_t kSigmaY=0.1; 
  const Double_t kSigmaZ=0.1;
  //!!! why divide by 1000000? This will be in cm/us then. But I guess everything should be in sec;
  //   const Double_t kDriftVel = fTPCParam->GetDriftV()/1000000;
  const Double_t kDriftVel = fTPCParam->GetDriftV();
  const Double_t kMaxZ0=fTPCParam->GetZLength();
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
    
  const Double_t iFCRadius=  83.5; //radius constants found in AliTPCCorrection.cxx 
  const Double_t oFCRadius= 254.5;
  
  AliToyMCTrack track(trackIn);
  const Int_t nMaxPoints = kNRows+50;
  AliTrackPointArray pointArray0(nMaxPoints);
  AliTrackPointArray pointArray1(nMaxPoints);
  Double_t xyz[3];

  //!!! when does the propagation not work, how often does it happen?
  if (!AliTrackerBase::PropagateTrackTo(&track,iFCRadius,kMass,5,kTRUE,kMaxSnp)) {
    AliError(Form("Propagation to IFC: %.2f failed\n",iFCRadius));
    return kFALSE;
  }

  Int_t npoints=0;
  Float_t covPoint[6]={0,0,0, kSigmaY*kSigmaY,0,kSigmaZ*kSigmaZ};  //covariance at the local frame


  //Simulate track from inner field cage radius to outer and save points along trajectory
  Int_t nMinPoints = 40;
  //!!! changed to loop over pad rows
  //      an increment of 1cm only might loose pad rows in IROC: 4x7.5 mm2
  //      or we change this to a smaller increment
  //    Do we need this initial creation of space points at all? can't we
  //      simply loop over all rows, crate the clusters in the center of the row,
  //      distort them and fill the distorted clusters
  
  //for (Double_t radius=iFCRadius; radius<oFCRadius; radius++){
  for (Int_t irow=0; irow<159; ++irow ){
    const Int_t isec    = irow<63?0:36;
    const Int_t isecrow = irow<63?irow:irow-63;
    Double_t radius = fTPCParam->GetPadRowRadii(isec,isecrow);

    //!!! changed from return 0 to continue -> Please check
    if (!AliTrackerBase::PropagateTrackTo(&track,radius,kMass,5,kTRUE,kMaxSnp)) {
      AliError(Form("Propagation to %d: %.2f failed\n",irow,radius));
      continue;
    }
    track.GetXYZ(xyz);

    //!!! Why is this smeared
    xyz[0]+=gRandom->Gaus(0,0.000005);
    xyz[1]+=gRandom->Gaus(0,0.000005);
    xyz[2]+=gRandom->Gaus(0,0.000005);

    if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    //!!! for future we might want to incude space points outside the active area
    //      since due to the distotions they might get inside
    //      in this case it will of course make sense to make a dense
    //      space point creation first and afterwards calculate the cluster
    //      but then we need to treat cluster and distorted cluster symmetrically (see below)
    if (TMath::Abs(track.GetX())<iFCRadius) continue;
    if (TMath::Abs(track.GetX())>oFCRadius) continue;

    AliTrackPoint pIn0;                               // space point          
    AliTrackPoint pIn1;
    Int_t sector= (xyz[2]>0)? 0:18;
    pointArray0.GetPoint(pIn0,npoints);
    pointArray1.GetPoint(pIn1,npoints);
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    Float_t distPoint[3]={xyz[0],xyz[1],xyz[2]};
    fSpaceCharge->DistortPoint(distPoint, sector);
    pIn0.SetXYZ(xyz[0], xyz[1],xyz[2]);
    pIn1.SetXYZ(distPoint[0], distPoint[1],distPoint[2]);
    track.Rotate(alpha);
    AliTrackPoint prot0 = pIn0.Rotate(alpha);   // rotate to the local frame - non distoted  point
    AliTrackPoint prot1 = pIn1.Rotate(alpha);   // rotate to the local frame -     distorted point
    prot0.SetXYZ(prot0.GetX(),prot0.GetY(), prot0.GetZ(),covPoint);
    prot1.SetXYZ(prot1.GetX(),prot1.GetY(), prot1.GetZ(),covPoint);
    pIn0=prot0.Rotate(-alpha);                       // rotate back to global frame
    pIn1=prot1.Rotate(-alpha);                       // rotate back to global frame
    pointArray0.AddPoint(npoints, &pIn0);      
    pointArray1.AddPoint(npoints, &pIn1);
    npoints++;
    if (npoints>=nMaxPoints) break;
  }

  if (npoints<nMinPoints) return 0;

  //save space points and make clusters of distorted points

  //!!! Here both, clusters and distorted clusters should be in the center of the pad row
  //      So in case the distorted clusters are calculated as the COG of all space points
  //      in a row, the same should be done for the undistorted clusters! There should
  //      only be one cluster per track per row
  Int_t lastRow = 0;
  Int_t pntsInCurrentRow = 0;
  Double_t xt=0, yt=0, zt=0;
  for(Int_t iPoint = 0; iPoint<npoints; iPoint++){

    //undistorted cluster, currently we assume one space point per pad row, see comments above
    AliTPCclusterMI *tempCl = trackIn.AddSpacePoint(AliTPCclusterMI());

    const Float_t *xArr = pointArray0.GetX();
    const Float_t *yArr = pointArray0.GetY();
    const Float_t *zArr = pointArray0.GetZ();

    tempCl->SetX(xArr[iPoint]);
    tempCl->SetY(yArr[iPoint]);
    tempCl->SetZ(zArr[iPoint]);
    
    Float_t xyzUC[3] = {xArr[iPoint],yArr[iPoint],zArr[iPoint]};
    Int_t indexUC[3] = {0};
    
    Float_t padRowUC = fTPCParam->GetPadRow(xyzUC,indexUC);
    Int_t   sectorUC = indexUC[1];

    // use pad row = 255 if outside active area
    if(padRowUC < 0 || (sectorUC < 36 && padRowUC >62) || (sectorUC > 35 &&padRowUC > 95) ) padRowUC=255; //outside sensitive area

    tempCl->SetPad(xyzUC[1]);
    tempCl->SetRow(padRowUC);
    tempCl->SetDetector(sectorUC);
    tempCl->SetTimeBin(t0 + (kMaxZ0-TMath::Abs(tempCl->GetZ()))/kDriftVel); // set time as t0  + drift time from dist z

    // Distorted cluster
    const Float_t *xArrDist = pointArray1.GetX();
    const Float_t *yArrDist = pointArray1.GetY();
    const Float_t *zArrDist = pointArray1.GetZ();
    
    Float_t xyz2[3] = {xArrDist[iPoint],yArrDist[iPoint],zArrDist[iPoint]};
    Int_t index[3] = {0};

    Float_t padRow = fTPCParam->GetPadRow(xyz2,index);
    Int_t sector = index[1];

    if(padRow < 0 || (sector < 36 && padRow>62) || (sector > 35 &&padRow > 95) ) continue; //outside sensitive area

    //!!! could it be that the last row is missed?
    //    assume row=95 and there is only 1 space point then this part will not be reached I think
    if(lastRow!=padRow) {
      

      //make and save distorted cluster
      if(pntsInCurrentRow>0){
	xt/=pntsInCurrentRow;
	yt/=pntsInCurrentRow;
	zt/=pntsInCurrentRow;
	AliTPCclusterMI *tempDistCl = trackIn.AddDistortedSpacePoint(AliTPCclusterMI());
	tempDistCl->SetX(xt);
	tempDistCl->SetY(yt);
	tempDistCl->SetZ(zt);
	Float_t clxyz[3] = {xt,yt,zt};
	Int_t clindex[3] = {0};
	Int_t clRow = fTPCParam->GetPadRow(clxyz,clindex);
// 	Int_t nPads = fTPCParam->GetNPads(clindex[1], clRow);
// 	Int_t pad = TMath::Nint(clxyz[1] + nPads/2); //taken from AliTPC.cxx
        //!!! pad should be fractional
// 	tempDistCl->SetPad(pad);
        tempDistCl->SetPad(clxyz[1]);
        tempDistCl->SetRow(clRow);
	tempDistCl->SetDetector(clindex[1]);
	tempDistCl->SetTimeBin(t0 + (kMaxZ0-TMath::Abs(zt))/kDriftVel); // set time as t0  + drift time from dist z coordinate to pad plane
      }
      xt = 0;
      yt = 0;
      zt = 0;
      pntsInCurrentRow = 0;
    }

    xt+=xArrDist[iPoint];
    yt+=yArrDist[iPoint];
    zt+=zArrDist[iPoint];
    pntsInCurrentRow++;
    lastRow = padRow;
  }
 
  return 1;

  
  
}
//________________________________________________________________
Bool_t AliToyMCEventGenerator::ConnectOutputFile()
{
  //
  // Create the output file name and tree and connect the event
  //

  fOutFile = new TFile(fOutputFileName.Data(),"recreate");

  if (!fOutFile || !fOutFile->IsOpen()){
    delete fOutFile;
    fOutFile=0x0;
    return kFALSE;
  }
  
  fOutTree = new TTree("toyMCtree","Tree with toyMC simulation");
  fOutTree->Branch("event","AliToyMCEvent",&fEvent);

  return kTRUE;
}

//________________________________________________________________
Bool_t AliToyMCEventGenerator::CloseOutputFile()
{
  //
  // close the output file
  //
  if (!fOutFile) return kFALSE;
  fOutFile->Write();
  fOutFile->Close();
  delete fOutFile;
  fOutFile=0x0;

  return kTRUE;
}

//________________________________________________________________
void AliToyMCEventGenerator::FillTree()
{
  // fill the tree
  if (fOutTree) fOutTree->Fill();
}

