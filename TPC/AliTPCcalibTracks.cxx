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
//     Class to analyse tracks for calibration                               //
//     to be used as a component in AliTPCSelectorTracks                     //
//     In the constructor you have to specify name and title                 //
//     to get the Object out of a file.                                      //
//     The parameter 'clusterParam', a AliTPCClusterParam object             //
//      (needed for TPC cluster error and shape parameterization)            //
//     Normally you get this object out of the file 'TPCClusterParam.root'   //
//     In the parameter 'cuts' the cuts are specified, that decide           //
//     weather a track will be accepted for calibration or not.              //
//                                                                           //
//       
//     The data flow:
//     
/*
   Raw Data -> Local Reconstruction -> Tracking ->     Calibration -> RefData (component itself)
               Offline/HLT             Offline/HLT                    OCDB entries (AliTPCClusterParam) 
*/            

/*

How to retrive it from file (created using calibration task):

gSystem->Load("libANALYSIS");
gSystem->Load("libTPCcalib");
TFile fcalib("CalibObjects.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTracks * calibTracks = ( AliTPCcalibTracks *)array->FindObject("calibTracks");


//USAGE of debug stream example
 gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chainres = tool.MakeChain("tracks.txt","ResolCl",0,10200);
  chainres->Lookup();
*/


                                                               //
///////////////////////////////////////////////////////////////////////////////

//
// ROOT includes 
//
#include <iostream>
#include <fstream>
using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
#include <TProfile.h>

//
//#include <TPDGCode.h>
#include <TStyle.h>
#include "TLinearFitter.h"
//#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"
#include <TCanvas.h>
#include <TGraph2DErrors.h>
#include "TPostScript.h"
#include "TCint.h"

#include <TH2D.h>
#include <TF2.h>
#include <TSystem.h>
#include <TCollection.h>
#include <iostream>
#include <TLinearFitter.h>
#include <TString.h>

//
// AliROOT includes 
//
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliTPCROC.h"

#include "AliTPCParamSR.h"
#include "AliTrackPointArray.h"
#include "AliTPCcalibTracks.h"
#include "AliTPCClusterParam.h"
#include "AliTPCcalibTracksCuts.h"
#include "AliTPCCalPadRegion.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TText.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TStatToolkit.h"
#include "TCut.h"



ClassImp(AliTPCcalibTracks)


AliTPCcalibTracks::AliTPCcalibTracks():
  AliTPCcalibBase(),
  fClusterParam(0),
  fROC(0),
  fArrayAmpRow(0),
  fArrayAmp(0), 
  fArrayQDY(0), 
  fArrayQDZ(0), 
  fArrayQRMSY(0),
  fArrayQRMSZ(0),
  fArrayChargeVsDriftlength(0),
  fcalPadRegionChargeVsDriftlength(0),
  fDeltaY(0),
  fDeltaZ(0),
  fResolY(0),
  fResolZ(0),
  fRMSY(0),
  fRMSZ(0),
  fCuts(0),
  fHclus(0),
  fRejectedTracksHisto(0),
  fHclusterPerPadrow(0),
  fHclusterPerPadrowRaw(0),
  fClusterCutHisto(0),
  fCalPadClusterPerPad(0),
  fCalPadClusterPerPadRaw(0)
{ 
   // 
   // AliTPCcalibTracks default constructor
   //    
  SetDebugLevel(1);
  if (GetDebugLevel() > 0) cout << "AliTPCcalibTracks' default constructor called" << endl;  
}   



AliTPCcalibTracks::AliTPCcalibTracks(const AliTPCcalibTracks& calibTracks):
  AliTPCcalibBase(calibTracks),
  fClusterParam(0),
  fROC(0),
  fArrayAmpRow(0),
  fArrayAmp(0), 
  fArrayQDY(0), 
  fArrayQDZ(0), 
  fArrayQRMSY(0),
  fArrayQRMSZ(0),
  fArrayChargeVsDriftlength(0),
  fcalPadRegionChargeVsDriftlength(0),
  fDeltaY(0),
  fDeltaZ(0),
  fResolY(0),
  fResolZ(0),
  fRMSY(0),
  fRMSZ(0),
  fCuts(0),
  fHclus(0),
  fRejectedTracksHisto(0),
  fHclusterPerPadrow(0),
  fHclusterPerPadrowRaw(0),
  fClusterCutHisto(0),
  fCalPadClusterPerPad(0),
  fCalPadClusterPerPadRaw(0)
{
   // 
   // AliTPCcalibTracks copy constructor
   // 
  if (GetDebugLevel() > 0) cout << " ***** this is AliTPCcalibTracks' copy constructor ***** " << endl;
   
   Bool_t dirStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   
   Int_t length = -1;
   // backward compatibility: if the data member doesn't yet exist, it will not be merged
   (calibTracks.fArrayAmpRow) ? length = calibTracks.fArrayAmpRow->GetEntriesFast() : length = -1;
   fArrayAmpRow = new TObjArray(length);
   fArrayAmp = new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fArrayAmpRow->AddAt( (TProfile*)calibTracks.fArrayAmpRow->At(i)->Clone(), i);
      fArrayAmp->AddAt( ((TProfile*)calibTracks.fArrayAmp->At(i)->Clone()), i);
   }
   
   (calibTracks.fArrayQDY) ? length = calibTracks.fArrayQDY->GetEntriesFast() : length = -1;
   fArrayQDY= new TObjArray(length);
   fArrayQDZ= new TObjArray(length);
   fArrayQRMSY= new TObjArray(length);
   fArrayQRMSZ= new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fArrayQDY->AddAt( ((TH1F*)calibTracks.fArrayQDY->At(i)->Clone()), i);
      fArrayQDZ->AddAt( ((TH1F*)calibTracks.fArrayQDZ->At(i)->Clone()), i);
      fArrayQRMSY->AddAt( ((TH1F*)calibTracks.fArrayQRMSY->At(i)->Clone()), i);
      fArrayQRMSZ->AddAt( ((TH1F*)calibTracks.fArrayQRMSZ->At(i)->Clone()), i);
   }
   
   (calibTracks.fResolY) ? length = calibTracks.fResolY->GetEntriesFast() : length = -1;
   fResolY = new TObjArray(length);
   fResolZ = new TObjArray(length);
   fRMSY = new TObjArray(length);
   fRMSZ = new TObjArray(length);
   for (Int_t i = 0; i < length; i++) {
      fResolY->AddAt( ((TH1F*)calibTracks.fResolY->At(i)->Clone()), i);
      fResolZ->AddAt( ((TH1F*)calibTracks.fResolZ->At(i)->Clone()), i);
      fRMSY->AddAt( ((TH1F*)calibTracks.fRMSY->At(i)->Clone()), i);
      fRMSZ->AddAt( ((TH1F*)calibTracks.fRMSZ->At(i)->Clone()), i);
   } 
   
   (calibTracks.fArrayChargeVsDriftlength) ? length = calibTracks.fArrayChargeVsDriftlength->GetEntriesFast() : length = -1;
   (calibTracks.fArrayChargeVsDriftlength) ? fArrayChargeVsDriftlength = new TObjArray(length) : fArrayChargeVsDriftlength = 0;
   for (Int_t i = 0; i < length; i++) {
      fArrayChargeVsDriftlength->AddAt( ((TProfile*)calibTracks.fArrayChargeVsDriftlength->At(i)->Clone()), i);
   }
   
   fDeltaY =  (TH1F*)calibTracks.fDeltaY->Clone();
   fDeltaZ =  (TH1F*)calibTracks.fDeltaZ->Clone();
   fHclus = (TH1I*)calibTracks.fHclus->Clone();
   fClusterCutHisto = (TH2I*)calibTracks.fClusterCutHisto->Clone();
   fRejectedTracksHisto    = (TH1I*)calibTracks.fRejectedTracksHisto->Clone();
   fHclusterPerPadrow      = (TH1I*)calibTracks.fHclusterPerPadrow->Clone();
   fHclusterPerPadrowRaw   = (TH1I*)calibTracks.fHclusterPerPadrowRaw->Clone();
   fcalPadRegionChargeVsDriftlength = (AliTPCCalPadRegion*)calibTracks.fcalPadRegionChargeVsDriftlength->Clone();
   fCalPadClusterPerPad    = (AliTPCCalPad*)calibTracks.fCalPadClusterPerPad->Clone();
   fCalPadClusterPerPadRaw = (AliTPCCalPad*)calibTracks.fCalPadClusterPerPadRaw->Clone();

   fCuts = new AliTPCcalibTracksCuts(calibTracks.fCuts->GetMinClusters(), calibTracks.fCuts->GetMinRatio(), 
      calibTracks.fCuts->GetMax1pt(), calibTracks.fCuts->GetEdgeYXCutNoise(), calibTracks.fCuts->GetEdgeThetaCutNoise());
   SetNameTitle(calibTracks.GetName(), calibTracks.GetTitle());
   TH1::AddDirectory(dirStatus); // set status back to original status
//    cout << "+++++ end of copy constructor +++++" << endl;   // TO BE REMOVED
}


AliTPCcalibTracks & AliTPCcalibTracks::operator=(const AliTPCcalibTracks& calibTracks){
  //
  // assgnment operator
  //
  if (this != &calibTracks) {
    new (this) AliTPCcalibTracks(calibTracks);
  }
  return *this;

}


AliTPCcalibTracks::AliTPCcalibTracks(const Text_t *name, const Text_t *title, AliTPCClusterParam *clusterParam,  AliTPCcalibTracksCuts* cuts, Int_t logLevel) : 
  AliTPCcalibBase(),
  fClusterParam(0),
  fROC(0),
  fArrayAmpRow(0),
  fArrayAmp(0), 
  fArrayQDY(0), 
  fArrayQDZ(0), 
  fArrayQRMSY(0),
  fArrayQRMSZ(0),
  fArrayChargeVsDriftlength(0),
  fcalPadRegionChargeVsDriftlength(0),
  fDeltaY(0),
  fDeltaZ(0),
  fResolY(0),
  fResolZ(0),
  fRMSY(0),
  fRMSZ(0),
  fCuts(0),
  fHclus(0),
  fRejectedTracksHisto(0),
  fHclusterPerPadrow(0),
  fHclusterPerPadrowRaw(0),
  fClusterCutHisto(0),
  fCalPadClusterPerPad(0),
  fCalPadClusterPerPadRaw(0)
 {
   // 
   // AliTPCcalibTracks constructor
   // specify 'name' and 'title' of your object
   // specify 'clusterParam', (needed for TPC cluster error and shape parameterization)
   // In the parameter 'cuts' the cuts are specified, that decide           
   // weather a track will be accepted for calibration or not.              
   //
   // fDebugLevel - debug output: -1: silence, 0: default, 1: things like constructor called, 5: write fDebugStreamer, 6: waste your screen
   // 
   // All histograms are instatiated in this constructor.
   // 
   this->SetName(name);
   this->SetTitle(title);

   if (GetDebugLevel() > 0) cout << " ***** this is AliTPCcalibTracks' main constructor ***** " << endl;
   G__SetCatchException(0);     
   
   fClusterParam = clusterParam;
   if (fClusterParam){
     fClusterParam->SetInstance(fClusterParam);
   }
   else {
     Error("AliTPCcalibTracks","No cluster parametrization found! A valid clusterParam object is needed in the constructor. (To be found in 'TPCClusterParam.root'.)");
   } 
   fCuts = cuts;
   SetDebugLevel(logLevel);
   
   TH1::AddDirectory(kFALSE);
   
   char chname[1000];
   TProfile * prof1=0;
   TH1F     * his1 =0;
   fHclus = new TH1I("hclus","Number of clusters per track",160, 0, 160);     // valgrind 3
   fRejectedTracksHisto    = new TH1I("RejectedTracksHisto", "Rejected tracks, sorted by failed cut", 100, -1, 10);
   fHclusterPerPadrow      = new TH1I("fHclusterPerPadrow", " clusters per padRow, used for the resolution tree", 160, 0, 160);
   fHclusterPerPadrowRaw   = new TH1I("fHclusterPerPadrowRaw", " clusters per padRow, before cutting clusters", 160, 0, 160);
   fCalPadClusterPerPad    = new AliTPCCalPad("fCalPadClusterPerPad", "clusters per pad");
   fCalPadClusterPerPadRaw = new AliTPCCalPad("fCalPadClusterPerPadRaw", "clusters per pad, before cutting clusters");
   fClusterCutHisto = new TH2I("fClusterCutHisto", "Cutted cluster over padRow; Cut Criterium; PadRow", 5,1,5, 160,0,159);
   
   // Amplitude  - sector - row histograms 
   fArrayAmpRow = new TObjArray(72);
   fArrayAmp    = new TObjArray(72);
   fArrayChargeVsDriftlength = new TObjArray(72);
   
   for (Int_t i = 0; i < 36; i++){   
      sprintf(chname,"Amp_row_Sector%d",i);
      prof1 = new TProfile(chname,chname,63,0,64);          // valgrind 3   193,536 bytes in 354 blocks are still reachable 
      prof1->SetXTitle("Pad row");
      prof1->SetYTitle("Mean Max amplitude");
      fArrayAmpRow->AddAt(prof1,i);
      sprintf(chname,"Amp_row_Sector%d",i+36);
      prof1 = new TProfile(chname,chname,96,0,97);       // valgrind 3   3,912 bytes in 6 blocks are possibly lost
      prof1->SetXTitle("Pad row");
      prof1->SetYTitle("Mean Max  amplitude");
      fArrayAmpRow->AddAt(prof1,i+36);
      
      // amplitude
      sprintf(chname,"Amp_Sector%d",i);
      his1 = new TH1F(chname,chname,100,0,32);         // valgrind 
      his1->SetXTitle("Max Amplitude (ADC)");
      fArrayAmp->AddAt(his1,i);
      sprintf(chname,"Amp_Sector%d",i+36);
      his1 = new TH1F(chname,chname,100,0,32);         // valgrind 3   13,408,208 bytes in 229 blocks are still reachable
      his1->SetXTitle("Max Amplitude (ADC)");
      fArrayAmp->AddAt(his1,i+36);
      
      // driftlength
      sprintf(chname, "driftlengt vs. charge, ROC %i", i);
      prof1 = new TProfile(chname, chname, 25, 0, 250);
      prof1->SetYTitle("Charge");
      prof1->SetXTitle("Driftlength");
      fArrayChargeVsDriftlength->AddAt(prof1,i);
      sprintf(chname, "driftlengt vs. charge, ROC %i", i+36);
      prof1 = new TProfile(chname, chname, 25, 0, 250);
      prof1->SetYTitle("Charge");
      prof1->SetXTitle("Driftlength");
      fArrayChargeVsDriftlength->AddAt(prof1,i+36);
   }
   
   TH1::AddDirectory(kFALSE);
   
   fDeltaY = new TH1F("DeltaY","DeltaY",100,-1,1);
   fDeltaZ = new TH1F("DeltaZ","DeltaZ",100,-1,1);
   
   fResolY = new TObjArray(3);
   fResolZ = new TObjArray(3);
   fRMSY   = new TObjArray(3);
   fRMSZ   = new TObjArray(3);
   TH3F * his3D;
   //
   his3D = new TH3F("Resol Y0","Resol Y0", 5,20,250, 4, 0,1., 50, -1,1);
   fResolY->AddAt(his3D,0);	
   his3D = new TH3F("Resol Y1","Resol Y1", 5,20,250, 4, 0,1., 50, -1,1);
   fResolY->AddAt(his3D,1);
   his3D = new TH3F("Resol Y2","Resol Y2", 5,20,250, 4, 0,0.8, 50, -1,1);
   fResolY->AddAt(his3D,2);
   //
   his3D = new TH3F("Resol Z0","Resol Z0", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,0);
   his3D = new TH3F("Resol Z1","Resol Z1", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,1);
   his3D = new TH3F("Resol Z2","Resol Z2", 5,20,250, 4, 0,1, 50, -1,1);
   fResolZ->AddAt(his3D,2);
   //
   his3D = new TH3F("RMS Y0","RMS Y0", 5,20,250, 4, 0,1., 50, 0,0.8);
   fRMSY->AddAt(his3D,0);
   his3D = new TH3F("RMS Y1","RMS Y1", 5,20,250, 4, 0,1., 50, 0,0.8);
   fRMSY->AddAt(his3D,1);
   his3D = new TH3F("RMS Y2","RMS Y2", 5,20,250, 4, 0,0.8, 50, 0,0.8);
   fRMSY->AddAt(his3D,2);
   //
   his3D = new TH3F("RMS Z0","RMS Z0", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,0);
   his3D = new TH3F("RMS Z1","RMS Z1", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,1);
   his3D = new TH3F("RMS Z2","RMS Z2", 5,20,250, 4, 0,1, 50, 0,0.8);
   fRMSZ->AddAt(his3D,2);
   //
      
   TH1::AddDirectory(kFALSE);
   
   fArrayQDY = new TObjArray(300);
   fArrayQDZ = new TObjArray(300);
   fArrayQRMSY = new TObjArray(300);
   fArrayQRMSZ = new TObjArray(300);
   for (Int_t iq = 0; iq <= 10; iq++){
      for (Int_t ipad = 0; ipad < 3; ipad++){
         Int_t   bin   = GetBin(iq, ipad);
         Float_t qmean = GetQ(bin);
         char hname[200];
         sprintf(hname,"ResolY Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(hname, hname, 20,10,250, 20, 0,1.5, 100, -1,1);
         fArrayQDY->AddAt(his3D, bin);
         sprintf(hname,"ResolZ Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(hname, hname, 20,10,250, 20, 0,1.5, 100, -1,1);
         fArrayQDZ->AddAt(his3D, bin);
         sprintf(hname,"RMSY Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(hname, hname, 20,10,250, 20, 0,1.5, 100, 0,0.6);
         fArrayQRMSY->AddAt(his3D, bin);
         sprintf(hname,"RMSZ Pad%d Qmiddle%f",ipad, qmean);
         his3D = new TH3F(hname, hname, 20,10,250, 20, 0,1.5, 100, 0,0.6);
         fArrayQRMSZ->AddAt(his3D, bin);
      }
   }
   
   fcalPadRegionChargeVsDriftlength = new AliTPCCalPadRegion("fcalPadRegionChargeVsDriftlength", "TProfiles with charge vs driftlength for each pad region");
   TProfile *tempProf;
   for (UInt_t padSize = 0; padSize < 3; padSize++) {
      for (UInt_t isector = 0; isector < 36; isector++) {
         if (padSize == 0) sprintf(chname, "driftlengt vs. charge, sector %i, short pads", isector);
         if (padSize == 1) sprintf(chname, "driftlengt vs. charge, sector %i, medium  pads", isector);
         if (padSize == 2) sprintf(chname, "driftlengt vs. charge, sector %i, long  pads", isector);
         tempProf = new TProfile(chname, chname, 500, 0, 250);
         tempProf->SetYTitle("Charge");
         tempProf->SetXTitle("Driftlength");
         fcalPadRegionChargeVsDriftlength->SetObject(tempProf, isector, padSize);
      }
   }
   

   if (GetDebugLevel() > 1) cout << "AliTPCcalibTracks object sucessfully constructed: " << GetName() << endl; 
   cout << "end of main constructor" << endl; // TO BE REMOVED
}    


AliTPCcalibTracks::~AliTPCcalibTracks() {
   // 
   // AliTPCcalibTracks destructor
   // 
   
  if (GetDebugLevel() > 0) cout << "AliTPCcalibTracks' destuctor called." << endl;
   Int_t length = 0;
   if (fArrayAmpRow) length = fArrayAmpRow->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fArrayAmpRow->At(i);  
      delete fArrayAmp->At(i);  
   }
   delete fArrayAmpRow;
   delete fArrayAmp;
   
   delete fDeltaY;
   delete fDeltaZ;
   
   if (fResolY) length = fResolY->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fResolY->At(i);
      delete fResolZ->At(i);
      delete fRMSY->At(i);
      delete fRMSZ->At(i);
   }
   delete fResolY;
   delete fResolZ;
   delete fRMSY;
   delete fRMSZ;
   
   if (fArrayQDY) length = fArrayQDY->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fArrayQDY->At(i);
      delete fArrayQDZ->At(i);
      delete fArrayQRMSY->At(i);
      delete fArrayQRMSZ->At(i);
   }
   
   if (fArrayChargeVsDriftlength) length = fArrayChargeVsDriftlength->GetEntriesFast();
   for (Int_t i = 0; i < length; i++){
      delete fArrayChargeVsDriftlength->At(i);
   }
   
    
   delete fArrayQDY;
   delete fArrayQDZ;
   delete fArrayQRMSY;
   delete fArrayQRMSZ;
   delete fArrayChargeVsDriftlength;
   
  delete fHclus;
  delete fRejectedTracksHisto;
  delete fClusterCutHisto;
  delete fHclusterPerPadrow;
  delete fHclusterPerPadrowRaw;
  if (fCalPadClusterPerPad)    delete fCalPadClusterPerPad;
  if (fCalPadClusterPerPadRaw) delete fCalPadClusterPerPadRaw;
  if(fcalPadRegionChargeVsDriftlength) {
     fcalPadRegionChargeVsDriftlength->Delete();
     delete fcalPadRegionChargeVsDriftlength;
  }
}
   
  

void AliTPCcalibTracks::Process(AliTPCseed *track){
   // 
   // To be called in the selector
   // first AcceptTrack is evaluated, then calls all the following analyse functions: 
   // FillResolutionHistoLocal(track)
   // AlignUpDown(track, esd)
   // 
  if (GetDebugLevel() > 5) Info("Process","Starting to process the track...");
   Int_t accpetStatus = AcceptTrack(track);
   if (accpetStatus == 0) {
      FillResolutionHistoLocal(track);
      // AlignUpDown(track, esd);
   }
   else fRejectedTracksHisto->Fill(accpetStatus);
}



Int_t AliTPCcalibTracks::GetBin(Float_t q, Int_t pad){
  //
  // calculate bins for given q and pad type 
  // used in TObjArray
  //
  Int_t res = TMath::Max( TMath::Nint((TMath::Sqrt(q) - 3.)), 0 );  
  res *= 3;
  res += pad;
  return res;
}


Int_t AliTPCcalibTracks::GetBin(Int_t iq, Int_t pad){
  //
  // calculate bins for given iq and pad type 
  // used in TObjArray
  //
  return iq * 3 + pad;;
}


Float_t AliTPCcalibTracks::GetQ(Int_t bin){
   // 
   // returns to bin belonging charge
   // (bin / 3 + 3)^2
   // 
   Int_t bin0 = bin / 3;
   bin0 += 3;
   return bin0 * bin0;
}


Float_t AliTPCcalibTracks::GetPad(Int_t bin){
   // 
   // returns to bin belonging pad
   // bin % 3
   // 
   return bin % 3; 
}



Int_t AliTPCcalibTracks::AcceptTrack(AliTPCseed * track){
  //
  // Function, that decides wheather a given track is accepted for 
  // the analysis or not. 
  // The cuts are specified in the AliTPCcalibTracksCuts object 'fCuts'
  // Returns 0 if a track is accepted or an integer different from 0 
  // to indicate the failed cut
  //
  const Int_t   kMinClusters  = fCuts->GetMinClusters();
  const Float_t kMinRatio     = fCuts->GetMinRatio();
  const Float_t kMax1pt       = fCuts->GetMax1pt();
  const Float_t kEdgeYXCutNoise    = fCuts->GetEdgeYXCutNoise();
  const Float_t kEdgeThetaCutNoise = fCuts->GetEdgeThetaCutNoise();
  
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
  if ( TMath::Abs(track->GetY() / track->GetX()) > kEdgeYXCutNoise )
    if ( TMath::Abs(track->GetTgl()) < kEdgeThetaCutNoise ) return 1;
  if (track->GetNumberOfClusters() < kMinClusters) return 2;
  Float_t ratio = track->GetNumberOfClusters() / (track->GetNFoundable() + 1.);
  if (ratio < kMinRatio) return 3;
//   Float_t mpt = track->Get1Pt();       // Get1Pt() doesn't exist any more
  Float_t mpt = track->GetSigned1Pt();
  if (TMath::Abs(mpt) > kMax1pt) return 4;
  //if (TMath::Abs(track->GetZ())>240.) return kFALSE;
  //if (TMath::Abs(track->GetZ())<10.) return kFALSE;
  //if (TMath::Abs(track->GetTgl())>0.03) return kFALSE;
  
  if (GetDebugLevel() > 20) Info("AcceptTrack","Track has been accepted.");  
  return 0;
}


void  AliTPCcalibTracks::FillResolutionHistoLocal(AliTPCseed * track){
   //
   // fill resolution histograms - localy - tracklet in the neighborhood
   // write debug information to 'TPCSelectorDebug.root'
   // 
   // _ the main function, called during track analysis _
   // 
   // loop over all padrows along the track
   // fit tracklets (length: 13 clusters) calculate mean chi^2 for this track-fit in Y and Z direction
   // 
   // loop again over all padrows along the track
   // fit tracklet (clusters in given padrow +- kDelta padrows) 
   // with polynom of 2nd order and two polynoms of 1st order
   // take both polynoms of 1st order, calculate difference of their parameters
   // add covariance matrixes and calculate chi2 of this difference
   // if this chi2 is bigger than a given threshold, assume that the current cluster is
   // a kink an goto next padrow
   // if not:
   // fill fArrayAmpRow, array with amplitudes vs. row for given sector
   // fill fArrayAmp, array with amplitude histograms for give sector
   // fill fRMSY, fRMSZ, fArrayQRMSY and fArrayQRMSZ, fDeltaY, fDeltaZ, fResolY, fResolZ, fArrayQDY, fArrayQDY
   // 
   // write debug information to 'TPCSelectorDebug.root'
   // only for every kDeltaWriteDebugStream'th padrow to reduce data volume 
   // and to avoid redundant data
   // 

  static TLinearFitter fFitterLinY1(2,"pol1");   //
  static TLinearFitter fFitterLinZ1(2,"pol1");   // 
  static TLinearFitter fFitterLinY2(2,"pol1");   // 
  static TLinearFitter fFitterLinZ2(2,"pol1");   //
  static TLinearFitter fFitterParY(3,"pol2");    // 
  static TLinearFitter fFitterParZ(3,"pol2");    //

  fFitterLinY1.StoreData(kFALSE);
  fFitterLinZ1.StoreData(kFALSE);
  fFitterLinY2.StoreData(kFALSE);
  fFitterLinZ2.StoreData(kFALSE);
  fFitterParY.StoreData(kFALSE);
  fFitterParZ.StoreData(kFALSE);


  if (GetDebugLevel() > 5) Info("FillResolutionHistoLocal"," ***** Start of FillResolutionHistoLocal *****");
   const Int_t   kDelta    = 10;          // delta rows to fit
   const Float_t kMinRatio = 0.75;        // minimal ratio
   //   const Float_t kCutChi2  = 6.;          // cut chi2 - left right  - kink removal
   const Float_t kErrorFraction = 0.5;    // use only clusters with small interpolation error - for error param
   const Int_t   kFirstLargePad = 127;    // medium pads -> long pads
   const Float_t kLargePadSize  = 1.5;    // factor between medium and long pads' area
   const Int_t   kDeltaWriteDebugStream  = 5;  // only for every kDeltaWriteDebugStream'th padrow debug information is calulated and written to debugstream
   TVectorD paramY0(2);
   TVectorD paramZ0(2);
   TVectorD paramY1(2);
   TVectorD paramZ1(2);
   TVectorD paramY2(3);
   TVectorD paramZ2(3);
   TMatrixD matrixY0(2,2);
   TMatrixD matrixZ0(2,2);
   TMatrixD matrixY1(2,2);
   TMatrixD matrixZ1(2,2);
   
   // estimate mean error
   Int_t nTrackletsAll = 0;       // number of tracklets for given track
   Float_t csigmaY     = 0;       // mean sigma for tracklet refit in Y direction
   Float_t csigmaZ     = 0;       // mean sigma for tracklet refit in Z direction
   Int_t nClusters     = 0;       // working variable, number of clusters per tracklet
   Int_t sectorG       = -1;      // working variable, sector of tracklet, has to stay constant for one tracklet
   
   fHclus->Fill(track->GetNumberOfClusters());      // for statistics overview
   // ---------------------------------------------------------------------
   for (Int_t irow = 0; irow < 159; irow++){
      // loop over all rows along the track
      // fit tracklets (length: 13 rows) with pol2 in Y and Z direction
      // calculate mean chi^2 for this track-fit in Y and Z direction
      AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
      if (!cluster0) continue;  // no cluster found
      Int_t sector = cluster0->GetDetector();
      fHclusterPerPadrowRaw->Fill(irow);
      
      Int_t ipad = TMath::Nint(cluster0->GetPad());
      Float_t value = fCalPadClusterPerPadRaw->GetCalROC(sector)->GetValue((sector<36)?irow:irow-64, TMath::Nint(cluster0->GetPad()));
      fCalPadClusterPerPadRaw->GetCalROC(sector)->SetValue((sector<36)?irow:irow-64, ipad, value + 1 );
      
      if (sector != sectorG){
         // track leaves sector before it crossed enough rows to fit / initialization
         nClusters = 0;
         fFitterParY.ClearPoints();
         fFitterParZ.ClearPoints();
         sectorG = sector;
      }
      else {
         nClusters++;
         Double_t x = cluster0->GetX();
         fFitterParY.AddPoint(&x, cluster0->GetY(), 1);
         fFitterParZ.AddPoint(&x, cluster0->GetZ(), 1);
         //
         if ( nClusters >= kDelta + 3 ){  
         // if more than 13 (kDelta+3) clusters were added to the fitters
         // fit the tracklet, increase trackletCounter
         fFitterParY.Eval();
         fFitterParZ.Eval();
         nTrackletsAll++;
         csigmaY += fFitterParY.GetChisquare() / (nClusters - 3.);
         csigmaZ += fFitterParZ.GetChisquare() / (nClusters - 3.);
         nClusters = -1;
         fFitterParY.ClearPoints();
         fFitterParZ.ClearPoints();
         }
      }
   }      // for (Int_t irow = 0; irow < 159; irow++)
   // mean chi^2 for all tracklet fits in Y and in Z direction: 
   csigmaY = TMath::Sqrt(csigmaY / nTrackletsAll);
   csigmaZ = TMath::Sqrt(csigmaZ / nTrackletsAll);
   // ---------------------------------------------------------------------
 
   for (Int_t irow = 0; irow < 159; irow++){
      // loop again over all rows along the track
      // do analysis
      // 
      Int_t nclFound = 0;  // number of clusters in the neighborhood
      Int_t ncl0 = 0;      // number of clusters in rows < rowOfCenterCluster
      Int_t ncl1 = 0;      // number of clusters in rows > rowOfCenterCluster
      AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
      if (!cluster0) continue;
      Int_t sector = cluster0->GetDetector();
      Float_t xref = cluster0->GetX();
         
      // Make Fit
      fFitterParY.ClearPoints();
      fFitterParZ.ClearPoints();
      fFitterLinY1.ClearPoints();
      fFitterLinZ1.ClearPoints();
      fFitterLinY2.ClearPoints();
      fFitterLinZ2.ClearPoints();
      
      // fit tracklet (clusters in given padrow +- kDelta padrows) 
      // with polynom of 2nd order and two polynoms of 1st order
      // take both polynoms of 1st order, calculate difference of their parameters
      // add covariance matrixes and calculate chi2 of this difference
      // if this chi2 is bigger than a given threshold, assume that the current cluster is
      // a kink an goto next padrow
      
      for (Int_t idelta = -kDelta; idelta <= kDelta; idelta++){
         // loop over irow +- kDelta rows (neighboured rows)
         // 
         // 
         if (idelta == 0) continue;                                // don't use center cluster
         if (idelta + irow < 0 || idelta + irow > 159) continue;   // don't go out of ROC
         AliTPCclusterMI * currentCluster = track->GetClusterPointer(irow + idelta);
         if (!currentCluster) continue;
         if (currentCluster->GetType() < 0) continue;
         if (currentCluster->GetDetector() != sector) continue;
         Double_t x = currentCluster->GetX() - xref;  // x = differece: current cluster - cluster @ irow
         nclFound++;
         if (idelta < 0){
         ncl0++;
         fFitterLinY1.AddPoint(&x, currentCluster->GetY(), csigmaY);
         fFitterLinZ1.AddPoint(&x, currentCluster->GetZ(), csigmaZ);
         }
         if (idelta > 0){
         ncl1++;
         fFitterLinY2.AddPoint(&x, currentCluster->GetY(), csigmaY);
         fFitterLinZ2.AddPoint(&x, currentCluster->GetZ(), csigmaZ);
         }
         fFitterParY.AddPoint(&x, currentCluster->GetY(), csigmaY);  
         fFitterParZ.AddPoint(&x, currentCluster->GetZ(), csigmaZ);  
      }  // loop over neighbourhood for fitter filling 


      
      if (nclFound < kDelta * kMinRatio) fRejectedTracksHisto->Fill(10);
      if (nclFound < kDelta * kMinRatio) fClusterCutHisto->Fill(1, irow);
      if (nclFound < kDelta * kMinRatio) continue;    // if not enough clusters (7.5) found in neighbourhood goto next padrow
      fFitterParY.Eval();
      fFitterParZ.Eval();
      Double_t chi2 = (fFitterParY.GetChisquare() + fFitterParZ.GetChisquare()) / (2. * nclFound - 6.);
      //if (chi2 > kCutChi2) fRejectedTracksHisto->Fill(9);
      //if (chi2 > kCutChi2) fClusterCutHisto->Fill(2, irow);
      //if (chi2 > kCutChi2) continue;   // if chi^2 is too big goto next padrow
      TTreeSRedirector *cstream = GetDebugStreamer();
      if (cstream){
	(*cstream)<<"Cut9"<<
	  "chi2="<<chi2<<
	  "\n";
      }
      // REMOVE KINK
      // only when there are enough clusters (4) in each direction
      if (ncl0 > 4){
         fFitterLinY1.Eval();
         fFitterLinZ1.Eval();
      }
      if (ncl1 > 4){
         fFitterLinY2.Eval();
         fFitterLinZ2.Eval();
      }
      
      if (ncl0 > 4 && ncl1 > 4){
         fFitterLinY1.GetCovarianceMatrix(matrixY0);
         fFitterLinY2.GetCovarianceMatrix(matrixY1);
         fFitterLinZ1.GetCovarianceMatrix(matrixZ0);
         fFitterLinZ2.GetCovarianceMatrix(matrixZ1);
         fFitterLinY2.GetParameters(paramY1);
         fFitterLinZ2.GetParameters(paramZ1);
         fFitterLinY1.GetParameters(paramY0);
         fFitterLinZ1.GetParameters(paramZ0);
         paramY0 -= paramY1;
         paramZ0 -= paramZ1;
         matrixY0 += matrixY1;
         matrixZ0 += matrixZ1;
         Double_t cchi2 = 0;
         
         TMatrixD difY(2, 1, paramY0.GetMatrixArray());
         TMatrixD difYT(1, 2, paramY0.GetMatrixArray());
         matrixY0.Invert();
         TMatrixD mulY(matrixY0, TMatrixD::kMult, difY);
         TMatrixD chi2Y(difYT, TMatrixD::kMult, mulY);
         cchi2 += chi2Y(0, 0);
         
         TMatrixD difZ(2, 1, paramZ0.GetMatrixArray());
         TMatrixD difZT(1, 2, paramZ0.GetMatrixArray());
         matrixZ0.Invert();
         TMatrixD mulZ(matrixZ0, TMatrixD::kMult, difZ);
         TMatrixD chi2Z(difZT, TMatrixD::kMult, mulZ);
         cchi2 += chi2Z(0, 0);      
         
         // REMOVE KINK - TO be fixed - proper chi2 calculation for curved track to be implemented
         //if (chi2 * 0.25 > kCutChi2) fRejectedTracksHisto->Fill(8);
         //if (chi2 * 0.25 > kCutChi2) fClusterCutHisto->Fill(3, irow);
         //if (chi2 * 0.25 > kCutChi2) continue;   // if chi2 is too big goto next padrow
         // fit tracklet with polynom of 2nd order and two polynoms of 1st order
         // take both polynoms of 1st order, calculate difference of their parameters
         // add covariance matrixes and calculate chi2 of this difference
         // if this chi2 is bigger than a given threshold, assume that the current cluster is
         // a kink an goto next padrow

	 if (cstream){
	   (*cstream)<<"Cut8"<<
	     "chi2="<<cchi2<<
	     "\n";
	 }	 
      }
      
      // current padrow has no kink
      
      // get fit parameters from pol2 fit: 
      Double_t paramY[4], paramZ[4];
      paramY[0] = fFitterParY.GetParameter(0);
      paramY[1] = fFitterParY.GetParameter(1);
      paramY[2] = fFitterParY.GetParameter(2);
      paramZ[0] = fFitterParZ.GetParameter(0);
      paramZ[1] = fFitterParZ.GetParameter(1);
      paramZ[2] = fFitterParZ.GetParameter(2);    
      
      Double_t tracky = paramY[0];
      Double_t trackz = paramZ[0];
      Float_t  deltay = tracky - cluster0->GetY();
      Float_t  deltaz = trackz - cluster0->GetZ();
      Float_t  angley = paramY[1] - paramY[0] / xref;
      Float_t  anglez = paramZ[1];
      
      Float_t max = cluster0->GetMax();
      UInt_t isegment = cluster0->GetDetector() % 36;
      Int_t padSize = 0;                          // short pads
      if (cluster0->GetDetector() >= 36) {
         padSize = 1;                              // medium pads 
         if (cluster0->GetRow() > 63) padSize = 2; // long pads
      }

      // =========================================
      // wirte collected information to histograms
      // =========================================
      
      TProfile *profAmpRow =  (TProfile*)fArrayAmpRow->At(sector);
      if ( irow >= kFirstLargePad) max /= kLargePadSize;
      Double_t smax = TMath::Sqrt(max);
      profAmpRow->Fill( (Double_t)cluster0->GetRow(), smax );
      TH1F *hisAmp =  (TH1F*)fArrayAmp->At(sector);
      hisAmp->Fill(smax);
      
      // remove the following two lines one day:
      TProfile *profDriftLength = (TProfile*)fArrayChargeVsDriftlength->At(sector);
      profDriftLength->Fill( 250.-(Double_t)TMath::Abs(cluster0->GetZ()), smax );
      
      TProfile *profDriftLengthTmp = (TProfile*)(fcalPadRegionChargeVsDriftlength->GetObject(isegment, padSize));
      profDriftLengthTmp->Fill( 250.-(Double_t)TMath::Abs(cluster0->GetZ()), smax );
      
      fHclusterPerPadrow->Fill(irow);   // fill histogram showing clusters per padrow
      Int_t ipad = TMath::Nint(cluster0->GetPad());
      Float_t value = fCalPadClusterPerPad->GetCalROC(sector)->GetValue((sector<36)?irow:irow-64, TMath::Nint(cluster0->GetPad()));
      fCalPadClusterPerPad->GetCalROC(sector)->SetValue((sector<36)?irow:irow-64, ipad, value + 1 );
   
         
      TH3F * his3 = 0;
      his3 = (TH3F*)fRMSY->At(padSize);
      if (his3) his3->Fill(250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), TMath::Sqrt(cluster0->GetSigmaY2()) );
      his3 = (TH3F*)fRMSZ->At(padSize);
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), TMath::Sqrt(cluster0->GetSigmaZ2()) );
      
      his3 = (TH3F*)fArrayQRMSY->At(GetBin(cluster0->GetMax(), padSize));
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), TMath::Sqrt(cluster0->GetSigmaY2()) );
      his3 = (TH3F*)fArrayQRMSZ->At(GetBin(cluster0->GetMax(), padSize));
      if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), TMath::Sqrt(cluster0->GetSigmaZ2()) );
   
      
      // Fill resolution histograms
      Bool_t useForResol = kTRUE;
      if (fFitterParY.GetParError(0) > kErrorFraction * csigmaY) useForResol = kFALSE;
   
      if (cstream){
	Float_t zdrift = 250 - TMath::Abs(cluster0->GetZ());
	Float_t sy = cluster0->GetSigmaY2();
	Float_t sz = cluster0->GetSigmaZ2();
	(*cstream)<<"Resol0"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field	      
	  "padSize="<<padSize<<
	  "angley="<<angley<<
	  "anglez="<<anglez<<
	  "zdr="<<zdrift<<
	  "dy="<<deltay<<
	  "dz="<<deltaz<<
	  "sy="<<sy<<
	  "sz="<<sz<<
	  "\n";
      }

      if (useForResol){
         fDeltaY->Fill(deltay);
         fDeltaZ->Fill(deltaz);
         his3 = (TH3F*)fResolY->At(padSize);
         if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(angley), deltay );
         his3 = (TH3F*)fResolZ->At(padSize);
         if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()), TMath::Abs(anglez), deltaz );
         his3 = (TH3F*)fArrayQDY->At(GetBin(cluster0->GetMax(), padSize));
         if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()),TMath::Abs(angley), deltay );
         his3 = (TH3F*)fArrayQDZ->At(GetBin(cluster0->GetMax(), padSize));
         if (his3) his3->Fill( 250 - TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez), deltaz );
      }
      
      //=============================================================================================
      
      if (useForResol && nclFound > 2 * kMinRatio * kDelta 
	  && irow % kDeltaWriteDebugStream == 0 && GetDebugLevel() > 4){
	if (GetDebugLevel() > 20) Info("FillResolutionHistoLocal","Filling 'TPCSelectorDebug.root', irow = %i", irow);
         FillResolutionHistoLocalDebugPart(track, cluster0, irow, angley, anglez, nclFound, kDelta);
      }  // if (useForResol && nclFound > 2 * kMinRatio * kDelta)
   
   }    // loop over all padrows along the track: for (Int_t irow = 0; irow < 159; irow++)
}  // FillResolutionHistoLocal(...)



void AliTPCcalibTracks::FillResolutionHistoLocalDebugPart(AliTPCseed *track, AliTPCclusterMI *cluster0, Int_t irow, Float_t  angley, Float_t  anglez, Int_t nclFound, Int_t kDelta) {
   // 
   //  - debug part of FillResolutionHistoLocal - 
   // called only for every kDeltaWriteDebugStream'th padrow, to avoid to much redundant data
   // called only for GetStreamLevel() > 4
   // fill resolution trees
   //
      
   Int_t sector = cluster0->GetDetector();
   Float_t xref = cluster0->GetX();
   Int_t padSize = 0;                          // short pads
   if (cluster0->GetDetector() >= 36) {
      padSize = 1;                              // medium pads 
      if (cluster0->GetRow() > 63) padSize = 2; // long pads
   }
      
   static TLinearFitter fitY0(3, "pol2");
   static TLinearFitter fitZ0(3, "pol2");
   static TLinearFitter fitY2(5, "hyp4");
   static TLinearFitter fitZ2(5, "hyp4");
   static TLinearFitter fitY2Q(5, "hyp4");
   static TLinearFitter fitZ2Q(5, "hyp4");
   static TLinearFitter fitY2S(5, "hyp4");
   static TLinearFitter fitZ2S(5, "hyp4");
   fitY0.ClearPoints();
   fitZ0.ClearPoints();
   fitY2.ClearPoints();
   fitZ2.ClearPoints();
   fitY2Q.ClearPoints();
   fitZ2Q.ClearPoints();
   fitY2S.ClearPoints();
   fitZ2S.ClearPoints();
   
   for (Int_t idelta = -kDelta; idelta <= kDelta; idelta++){
     // loop over irow +- kDelta rows (neighboured rows)
     // 
     // 
     if (idelta == 0) continue;
     if (idelta + irow < 0 || idelta + irow > 159) continue;   // don't go out of ROC
     AliTPCclusterMI * cluster = track->GetClusterPointer(irow + idelta);
     if (!cluster) continue;
     if (cluster->GetType() < 0) continue;
     if (cluster->GetDetector() != sector) continue;
     Double_t x = cluster->GetX() - xref;
     Double_t sigmaY0 = fClusterParam->GetError0Par( 0, padSize, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(angley) );
     Double_t sigmaZ0 = fClusterParam->GetError0Par( 1, padSize, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(anglez) );
     //
     Double_t sigmaYQ = fClusterParam->GetErrorQPar( 0, padSize, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(angley), TMath::Abs(cluster->GetMax()) );
     Double_t sigmaZQ = fClusterParam->GetErrorQPar( 1, padSize, (250.0 - TMath::Abs(cluster->GetZ())), TMath::Abs(anglez), TMath::Abs(cluster->GetMax()) );
     Double_t sigmaYS = fClusterParam->GetErrorQParScaled( 0, padSize, (250.0 - TMath::Abs(cluster->GetZ())), 
                                                           TMath::Abs(angley), TMath::Abs(cluster->GetMax()) );
     Double_t sigmaZS = fClusterParam->GetErrorQParScaled( 1, padSize, (250.0 - TMath::Abs(cluster->GetZ())), 
                                                           TMath::Abs(anglez), TMath::Abs(cluster->GetMax()) );
     Float_t rmsYFactor = fClusterParam->GetShapeFactor( 0, padSize,(250.0 - TMath::Abs(cluster->GetZ())),
							 TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							 TMath::Sqrt(cluster0->GetSigmaY2()), 0 );
     Float_t rmsZFactor = fClusterParam->GetShapeFactor(0, padSize,(250.0 - TMath::Abs(cluster->GetZ())),
							TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							TMath::Sqrt(cluster0->GetSigmaZ2()),0 );
     sigmaYS  = TMath::Sqrt(sigmaYS * sigmaYS + rmsYFactor * rmsYFactor / 12.);
     sigmaZS  = TMath::Sqrt(sigmaZS * sigmaZS + rmsZFactor * rmsZFactor / 12. + rmsYFactor * rmsYFactor / 24.);
     //
     if (kDelta != 0){
       fitY0.AddPoint(&x, cluster->GetY(), sigmaY0);
       fitZ0.AddPoint(&x, cluster->GetZ(), sigmaZ0);
     }
     Double_t xxx[4];
     xxx[0] = ( (idelta+irow) % 2 == 0 ) ? 1 : 0;
     xxx[1] = x;
     xxx[2] = ( (idelta+irow) % 2 == 0 ) ? x : 0;
     xxx[3] = x * x;	
     fitY2.AddPoint(xxx, cluster->GetY(), sigmaY0);
     fitY2Q.AddPoint(xxx, cluster->GetY(), sigmaYQ);
     fitY2S.AddPoint(xxx, cluster->GetY(), sigmaYS);
     fitZ2.AddPoint(xxx, cluster->GetZ(), sigmaZ0);
     fitZ2Q.AddPoint(xxx, cluster->GetZ(), sigmaZQ);
     fitZ2S.AddPoint(xxx, cluster->GetZ(), sigmaZS);
     //
   }  // neigbouhood-loop
      //
   fitY0.Eval();
   fitZ0.Eval();
   fitY2.Eval();
   fitZ2.Eval();
   fitY2Q.Eval();
   fitZ2Q.Eval();
   fitY2S.Eval();
   fitZ2S.Eval();
   Float_t chi2Y0 = fitY0.GetChisquare() / (nclFound-3.);
   Float_t chi2Z0 = fitZ0.GetChisquare() / (nclFound-3.);
   Float_t chi2Y2 = fitY2.GetChisquare() / (nclFound-5.);
   Float_t chi2Z2 = fitZ2.GetChisquare() / (nclFound-5.);
   Float_t chi2Y2Q = fitY2Q.GetChisquare() / (nclFound-5.);
   Float_t chi2Z2Q = fitZ2Q.GetChisquare() / (nclFound-5.);
   Float_t chi2Y2S = fitY2S.GetChisquare() / (nclFound-5.);
   Float_t chi2Z2S = fitZ2S.GetChisquare() / (nclFound-5.);
   //
   static  TVectorD    parY0(3);
   static  TMatrixD    matY0(3, 3);
   static  TVectorD    parZ0(3);
   static  TMatrixD    matZ0(3, 3);
   fitY0.GetParameters(parY0);
   fitY0.GetCovarianceMatrix(matY0);
   fitZ0.GetParameters(parZ0);
   fitZ0.GetCovarianceMatrix(matZ0);
   //
   static  TVectorD    parY2(5);
   static  TMatrixD    matY2(5,5);
   static  TVectorD    parZ2(5);
   static  TMatrixD    matZ2(5,5);
   fitY2.GetParameters(parY2);
   fitY2.GetCovarianceMatrix(matY2);
   fitZ2.GetParameters(parZ2);
   fitZ2.GetCovarianceMatrix(matZ2);
   //
   static  TVectorD    parY2Q(5);
   static  TMatrixD    matY2Q(5,5);
   static  TVectorD    parZ2Q(5);
   static  TMatrixD    matZ2Q(5,5);
   fitY2Q.GetParameters(parY2Q);
   fitY2Q.GetCovarianceMatrix(matY2Q);
   fitZ2Q.GetParameters(parZ2Q);
   fitZ2Q.GetCovarianceMatrix(matZ2Q);
   static  TVectorD    parY2S(5);
   static  TMatrixD    matY2S(5,5);
   static  TVectorD    parZ2S(5);
   static  TMatrixD    matZ2S(5,5);
   fitY2S.GetParameters(parY2S);
   fitY2S.GetCovarianceMatrix(matY2S);
   fitZ2S.GetParameters(parZ2S);
   fitZ2S.GetCovarianceMatrix(matZ2S);
   Float_t sigmaY0   = TMath::Sqrt(matY0(0,0));
   Float_t sigmaZ0   = TMath::Sqrt(matZ0(0,0));
   Float_t sigmaDY0  = TMath::Sqrt(matY0(1,1));
   Float_t sigmaDZ0  = TMath::Sqrt(matZ0(1,1));
   Float_t sigmaY2   = TMath::Sqrt(matY2(1,1));
   Float_t sigmaZ2   = TMath::Sqrt(matZ2(1,1));
   Float_t sigmaDY2  = TMath::Sqrt(matY2(3,3));
   Float_t sigmaDZ2  = TMath::Sqrt(matZ2(3,3));
   Float_t sigmaY2Q  = TMath::Sqrt(matY2Q(1,1));
   Float_t sigmaZ2Q  = TMath::Sqrt(matZ2Q(1,1));
   Float_t sigmaDY2Q = TMath::Sqrt(matY2Q(3,3));
   Float_t sigmaDZ2Q = TMath::Sqrt(matZ2Q(3,3));
   Float_t sigmaY2S  = TMath::Sqrt(matY2S(1,1));
   Float_t sigmaZ2S  = TMath::Sqrt(matZ2S(1,1));
   Float_t sigmaDY2S = TMath::Sqrt(matY2S(3,3));
   Float_t sigmaDZ2S = TMath::Sqrt(matZ2S(3,3));
   
   // Error parameters
   Float_t csigmaY0 = fClusterParam->GetError0Par(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(angley));
   Float_t csigmaZ0 = fClusterParam->GetError0Par(1,padSize,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(anglez));
   //
   Float_t csigmaYQ = fClusterParam->GetErrorQPar(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						  TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
   Float_t csigmaZQ = fClusterParam->GetErrorQPar(1,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						  TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
   Float_t csigmaYS = fClusterParam->GetErrorQParScaled(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
							TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
   Float_t csigmaZS = fClusterParam->GetErrorQParScaled(1,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
							TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
   
   
   // RMS parameters
   Float_t meanRMSY = 0;
   Float_t meanRMSZ = 0;
   Int_t   nclRMS = 0;
   for (Int_t idelta = -2; idelta <= 2; idelta++){
     // loop over neighbourhood
     if (idelta+irow < 0 || idelta+irow > 159) continue;
     // 	if (idelta+irow>159) continue;
     AliTPCclusterMI * cluster = track->GetClusterPointer(irow+idelta);
     if (!cluster) continue;
     meanRMSY += TMath::Sqrt(cluster->GetSigmaY2());
     meanRMSZ += TMath::Sqrt(cluster->GetSigmaZ2());
     nclRMS++;
   }
   meanRMSY /= nclRMS; 
   meanRMSZ /= nclRMS; 
   
   Float_t rmsY      = TMath::Sqrt(cluster0->GetSigmaY2());  
   Float_t rmsZ      = TMath::Sqrt(cluster0->GetSigmaZ2());
   Float_t rmsYT     = fClusterParam->GetRMSQ(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
					      TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
   Float_t rmsZT     = fClusterParam->GetRMSQ(1,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
					      TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
   Float_t rmsYT0    = fClusterParam->GetRMS0(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
					      TMath::Abs(angley));
   Float_t rmsZT0    = fClusterParam->GetRMS0(1,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						 TMath::Abs(anglez));
   Float_t rmsYSigma = fClusterParam->GetRMSSigma(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						  TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
   Float_t rmsZSigma = fClusterParam->GetRMSSigma(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						  TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
   Float_t rmsYFactor = fClusterParam->GetShapeFactor(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						      TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
						      rmsY,meanRMSY);
   Float_t rmsZFactor = fClusterParam->GetShapeFactor(0,padSize,(250.0-TMath::Abs(cluster0->GetZ())),
						      TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
						      rmsZ,meanRMSZ);
   //
   // cluster debug
   TTreeSRedirector *cstream = GetDebugStreamer();
   if (cstream){
     (*cstream)<<"ResolCl"<<	// valgrind 3   40,000 bytes in 1 blocks are possibly 
       "run="<<fRun<<              //  run number
       "event="<<fEvent<<          //  event number
       "time="<<fTime<<            //  time stamp of event
       "trigger="<<fTrigger<<      //  trigger
       "mag="<<fMagF<<             //  magnetic field	      
       "Sector="<<sector<<
       "Cl.="<<cluster0<<
       "CSigmaY0="<<csigmaY0<<   // cluster errorY
       "CSigmaYQ="<<csigmaYQ<<   // cluster errorY - q scaled
       "CSigmaYS="<<csigmaYS<<   // cluster errorY - q scaled
       "CSigmaZ0="<<csigmaZ0<<   // 
       "CSigmaZQ="<<csigmaZQ<<
       "CSigmaZS="<<csigmaZS<<
       "shapeYF="<<rmsYFactor<<
       "shapeZF="<<rmsZFactor<<
       "rmsY="<<rmsY<<
       "rmsZ="<<rmsZ<<
       "rmsYM="<<meanRMSY<<
       "rmsZM="<<meanRMSZ<<
       "rmsYT="<<rmsYT<<
       "rmsZT="<<rmsZT<<
       "rmsYT0="<<rmsYT0<<
       "rmsZT0="<<rmsZT0<<
       "rmsYS="<<rmsYSigma<<  
       "rmsZS="<<rmsZSigma<<
       "padSize="<<padSize<<
       "Ncl="<<nclFound<<	
       "PY0.="<<&parY0<<
       "PZ0.="<<&parZ0<<
       "SigmaY0="<<sigmaY0<< 
       "SigmaZ0="<<sigmaZ0<< 
       "angley="<<angley<<
       "anglez="<<anglez<<
       "\n";
     
     //       tracklet dubug
     (*cstream)<<"ResolTr"<<	
       "run="<<fRun<<              //  run number
       "event="<<fEvent<<          //  event number
       "time="<<fTime<<            //  time stamp of event
       "trigger="<<fTrigger<<      //  trigger
       "mag="<<fMagF<<             //  magnetic field	      
       "padSize="<<padSize<<
       "IPad="<<padSize<<
       "Sector="<<sector<<
       "Ncl="<<nclFound<<	
       "chi2Y0="<<chi2Y0<<
       "chi2Z0="<<chi2Z0<<
       "chi2Y2="<<chi2Y2<<
       "chi2Z2="<<chi2Z2<<
       "chi2Y2Q="<<chi2Y2Q<<
       "chi2Z2Q="<<chi2Z2Q<<
       "chi2Y2S="<<chi2Y2S<<
       "chi2Z2S="<<chi2Z2S<<
       "PY0.="<<&parY0<<
       "PZ0.="<<&parZ0<<
       "PY2.="<<&parY2<<
       "PZ2.="<<&parZ2<<
       "PY2Q.="<<&parY2Q<<
       "PZ2Q.="<<&parZ2Q<<
       "PY2S.="<<&parY2S<<
       "PZ2S.="<<&parZ2S<<
       "SigmaY0="<<sigmaY0<< 
       "SigmaZ0="<<sigmaZ0<< 
       "SigmaDY0="<<sigmaDY0<< 
       "SigmaDZ0="<<sigmaDZ0<< 
       "SigmaY2="<<sigmaY2<< 
       "SigmaZ2="<<sigmaZ2<< 
       "SigmaDY2="<<sigmaDY2<< 
       "SigmaDZ2="<<sigmaDZ2<< 
       "SigmaY2Q="<<sigmaY2Q<< 
       "SigmaZ2Q="<<sigmaZ2Q<< 
       "SigmaDY2Q="<<sigmaDY2Q<< 
       "SigmaDZ2Q="<<sigmaDZ2Q<< 
       "SigmaY2S="<<sigmaY2S<< 
       "SigmaZ2S="<<sigmaZ2S<< 
       "SigmaDY2S="<<sigmaDY2S<< 
       "SigmaDZ2S="<<sigmaDZ2S<< 
       "angley="<<angley<<
       "anglez="<<anglez<<
       "\n";
   }  
}





TH2D * AliTPCcalibTracks::MakeDiff(TH2D * hfit, TF2 * func){
   // 
   // creates a new histogram which contains the difference between
   // the histogram hfit and the function func
   // 
  TH2D * result = (TH2D*)hfit->Clone();      // valgrind 3   40,139 bytes in 11 blocks are still reachable
  result->SetTitle(Form("%s fit residuals",result->GetTitle()));
  result->SetName(Form("%s fit residuals",result->GetName()));
  TAxis *xaxis = hfit->GetXaxis();
  TAxis *yaxis = hfit->GetYaxis();
  Double_t x[2];
  for (Int_t biny = 0; biny <= yaxis->GetNbins(); biny++) {
    x[1]  = yaxis->GetBinCenter(biny);
    for (Int_t binx = 0; binx <= xaxis->GetNbins(); binx++) {
      x[0]  = xaxis->GetBinCenter(binx);
      Int_t bin = hfit->GetBin(binx, biny);
      Double_t val = hfit->GetBinContent(bin);
//      result->SetBinContent( bin, (val - func->Eval(x[0], x[1])) / func->Eval(x[0], x[1]) );
      result->SetBinContent( bin, (val / func->Eval(x[0], x[1])) - 1 );
    }    
  }
  return result;
}


void  AliTPCcalibTracks::SetStyle() const {
   // 
   // set style, can be called by all draw functions
   // 
   gROOT->SetStyle("Plain");
   gStyle->SetFillColor(10);
   gStyle->SetPadColor(10);
   gStyle->SetCanvasColor(10);
   gStyle->SetStatColor(10);
   gStyle->SetPalette(1,0);
   gStyle->SetNumberContours(60);
}


void AliTPCcalibTracks::Draw(Option_t* opt){
   // 
   // draw-function of AliTPCcalibTracks
   // will draws some exemplaric pictures
   // 

  if (GetDebugLevel() > 6) Info("Draw", "Drawing an exemplaric picture.");
   SetStyle();
   Double_t min = 0;
   Double_t max = 0;
   TCanvas *c1 = new TCanvas();
   c1->Divide(0, 3);
   TVirtualPad *upperThird = c1->GetPad(1);
   TVirtualPad *middleThird = c1->GetPad(2);
   TVirtualPad *lowerThird = c1->GetPad(3);
   upperThird->Divide(2,0);
   TVirtualPad *upleft  = upperThird->GetPad(1);
   TVirtualPad *upright = upperThird->GetPad(2);
   middleThird->Divide(2,0);
   TVirtualPad *middleLeft  = middleThird->GetPad(1);
   TVirtualPad *middleRight = middleThird->GetPad(2);
   lowerThird->Divide(2,0);
   TVirtualPad *downLeft  = lowerThird->GetPad(1);
   TVirtualPad *downRight = lowerThird->GetPad(2);
   
   
   upleft->cd(0);
   min = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())-20;
   max = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())+20;
   fDeltaY->SetAxisRange(min, max);
   fDeltaY->Fit("gaus","q","",min, max);        // valgrind 3  7 block possibly lost   2,400 bytes in 1 blocks are still reachable
   c1->Update();
   
   upright->cd(0);
   max = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())+20;
   min = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())-20;
   fDeltaZ->SetAxisRange(min, max);
   fDeltaZ->Fit("gaus","q","",min, max);
   c1->Update();
   
   middleLeft->cd();
   fHclus->Draw(opt);
   
   middleRight->cd();
   fRejectedTracksHisto->Draw(opt);
   TPaveText *pt = new TPaveText(0.6,0.6, 0.8,0.8, "NDC");
   TText *t1 = pt->AddText("1: kEdgeThetaCutNoise");
   TText *t2 = pt->AddText("2: kMinClusters");
   TText *t3 = pt->AddText("3: kMinRatio");
   TText *t4 = pt->AddText("4: kMax1pt");
   t1 = t1; t2 = t2; t3 = t3; t4 = t4;    // avoid compiler warnings
   pt->SetToolTipText("Legend for failed cuts");
   pt->Draw();
   
   downLeft->cd();
   fHclusterPerPadrowRaw->Draw(opt);
   
   downRight->cd();
   fHclusterPerPadrow->Draw(opt);
}


void AliTPCcalibTracks::MakeReport(Int_t stat, const char* pathName){ 
   // 
   // all functions are called, that produce pictures
   // the histograms are written to the directory 'pathName'
   // 'stat' is a threshhold: only histograms with more than 'stat' entries are wirtten to file
   // 'stat' is also the number of minEntries for MakeResPlotsQTree
   // 

  if (GetDebugLevel() > 0) Info("MakeReport","Writing plots and trees to '%s'.", pathName);
   MakeAmpPlots(stat, pathName);
   MakeDeltaPlots(pathName);
   FitResolutionNew(pathName);
   FitRMSNew(pathName);
   MakeChargeVsDriftLengthPlots(pathName);
//    MakeResPlotsQ(1, 1); 
   MakeResPlotsQTree(stat, pathName);
}
   

void AliTPCcalibTracks::MakeAmpPlots(Int_t stat, const char* pathName){ 
   // 
   // creates several plots:
   // fArrayAmp.ps, fArrayAmpRow.ps and DeltaYZ.ps
   // fArrayAmp.ps: one histogram per sector, the histogram shows the charge per cluster
   // fArrayAmpRow.ps: one histogram per sector, mean max. amplitude vs. pad row with landau fit
   // DeltaYZ.ps: DeltaY and DeltaZ histogram with gaus fit
   // Empty histograms (sectors without data) are not written to file
   // the ps-files are written to the directory 'pathName', that is created if it does not exist
   // 'stat': only histograms with more than 'stat' entries are written to file.
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas();     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
   TPostScript *ps; 
   // histograms with accumulated amplitude for all IROCs and OROCs
   TH1F *allAmpHisIROC = ((TH1F*)(fArrayAmp->At(0))->Clone());
   allAmpHisIROC->SetName("Amp all IROCs");
   allAmpHisIROC->SetTitle("Amp all IROCs");
   TH1F *allAmpHisOROC = ((TH1F*)(fArrayAmp->At(36))->Clone());
   allAmpHisOROC->SetName("Amp all OROCs");
   allAmpHisOROC->SetTitle("Amp all OROCs");
   
   ps = new TPostScript("fArrayAmp.ps", 112);
   if (GetDebugLevel() > -1) cout << "creating fArrayAmp.ps..." << endl;
   for (Int_t i = 0; i < fArrayAmp->GetEntriesFast(); i++){
      if ( ((TH1F*)fArrayAmp->At(i))->GetEntries() < stat  ) continue;
      ps->NewPage();
      ((TH1F*)fArrayAmp->At(i))->Draw();
      c1->Update();              // valgrind 3
      if (i > 0 && i < 36) {
         allAmpHisIROC->Add(((TH1F*)fArrayAmp->At(i)));
         allAmpHisOROC->Add(((TH1F*)fArrayAmp->At(i+36)));
      }
   }
   ps->NewPage();
   allAmpHisIROC->Draw();
   c1->Update();              // valgrind
   ps->NewPage();
   allAmpHisOROC->Draw();
   c1->Update();
   ps->Close();
   delete ps;
   
   TH1F *his = 0;
   Double_t min = 0;
   Double_t max = 0;
   ps = new TPostScript("fArrayAmpRow.ps", 112);
   if (GetDebugLevel() > -1) cout << "creating fArrayAmpRow.ps..." << endl;
   for (Int_t i = 0; i < fArrayAmpRow->GetEntriesFast(); i++){
      his = (TH1F*)fArrayAmpRow->At(i);
      if (his->GetEntries() < stat) continue;
      ps->NewPage();
      min = TMath::Max( his->GetBinCenter(his->GetMaximumBin() )-100., 0.);
      max = his->GetBinCenter(5*his->GetMaximumBin()) + 100;
      his->SetAxisRange(min, max);
      his->Fit("pol3", "q", "", min, max);
      // his->Draw("error");    // don't use this line when you don't want to have empty pages in the ps-file
      c1->Update();
   }
   ps->Close();
   delete ps;
   delete c1;
   gSystem->ChangeDirectory("..");
}


void AliTPCcalibTracks::MakeDeltaPlots(const char* pathName){
   // 
   // creates several plots:
   // DeltaYZ.ps: DeltaY and DeltaZ histogram with gaus fit
   // the ps-files are written to the directory 'pathName', that is created if it does not exist
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas();     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
   TPostScript *ps; 
   Double_t min = 0;
   Double_t max = 0;
   
   ps = new TPostScript("DeltaYZ.ps", 112);
   if (GetDebugLevel() > -1) cout << "creating DeltaYZ.ps..." << endl;
   min = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())-20;
   max = fDeltaY->GetBinCenter(fDeltaY->GetMaximumBin())+20;
   fDeltaY->SetAxisRange(min, max);
   ps->NewPage();
   fDeltaY->Fit("gaus","q","",min, max);        // valgrind 3  7 block possibly lost   2,400 bytes in 1 blocks are still reachable
   c1->Update();
   ps->NewPage();
   max = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())+20;
   min = fDeltaZ->GetBinCenter(fDeltaZ->GetMaximumBin())-20;
   fDeltaZ->SetAxisRange(min, max);
   fDeltaZ->Fit("gaus","q","",min, max);
   c1->Update();
   ps->Close();
   delete ps;
   delete c1;
   gSystem->ChangeDirectory("..");
}


void AliTPCcalibTracks::MakeChargeVsDriftLengthPlotsOld(const char* pathName){
   // 
   // creates charge vs. driftlength plots, one TProfile for each ROC
   // is not correct like this, should be one TProfile for each sector and padsize
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas();     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
   TPostScript *ps; 
   ps = new TPostScript("chargeVsDriftlengthOld.ps", 112);
   if (GetDebugLevel() > -1) cout << "creating chargeVsDriftlength.ps..." << endl;
   TProfile *chargeVsDriftlengthAllIROCs = ((TProfile*)fArrayChargeVsDriftlength->At(0)->Clone());
   TProfile *chargeVsDriftlengthAllOROCs = ((TProfile*)fArrayChargeVsDriftlength->At(36)->Clone());
   chargeVsDriftlengthAllIROCs->SetName("allAmpHisIROC");
   chargeVsDriftlengthAllIROCs->SetTitle("charge vs. driftlength, all IROCs");
   chargeVsDriftlengthAllOROCs->SetName("allAmpHisOROC");
   chargeVsDriftlengthAllOROCs->SetTitle("charge vs. driftlength, all OROCs");
   
   for (Int_t i = 0; i < fArrayChargeVsDriftlength->GetEntriesFast(); i++) {
      ((TProfile*)fArrayChargeVsDriftlength->At(i))->Draw();
      c1->Update();
      if (i > 0 && i < 36) { 
         chargeVsDriftlengthAllIROCs->Add(((TProfile*)fArrayChargeVsDriftlength->At(i)));
         chargeVsDriftlengthAllOROCs->Add(((TProfile*)fArrayChargeVsDriftlength->At(i+36)));
      }
      ps->NewPage();
   }
   chargeVsDriftlengthAllIROCs->Draw();
   c1->Update();              // valgrind
   ps->NewPage();
   chargeVsDriftlengthAllOROCs->Draw();
   c1->Update();
   ps->Close();  
   delete ps;
   delete c1;
   gSystem->ChangeDirectory("..");
}   


void AliTPCcalibTracks::MakeChargeVsDriftLengthPlots(const char* pathName){
   // 
   // creates charge vs. driftlength plots, one TProfile for each ROC
   // under development....
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas* c1 = new TCanvas("c1", "c1", 700,(Int_t)(TMath::Sqrt(2)*700));     // valgrind 3 ???  634 bytes in 28 blocks are still reachable
//    TCanvas c1("c1", "c1", 500,(sqrt(2)*500))
   c1->Divide(0,3);
   TPostScript *ps; 
   ps = new TPostScript("chargeVsDriftlength.ps", 111);
   if (GetDebugLevel() > -1) cout << "creating chargeVsDriftlengthNew.ps..." << endl;
   
   TProfile *chargeVsDriftlengthAllShortPads  = ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,0)->Clone());
   TProfile *chargeVsDriftlengthAllMediumPads = ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,1)->Clone());
   TProfile *chargeVsDriftlengthAllLongPads   = ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,2)->Clone());
   chargeVsDriftlengthAllShortPads->SetName("allAmpHisShortPads");
   chargeVsDriftlengthAllShortPads->SetTitle("charge vs. driftlength, all sectors, short pads");
   chargeVsDriftlengthAllMediumPads->SetName("allAmpHisMediumPads");
   chargeVsDriftlengthAllMediumPads->SetTitle("charge vs. driftlength, all sectors, medium pads");
   chargeVsDriftlengthAllLongPads->SetName("allAmpHisLongPads");
   chargeVsDriftlengthAllLongPads->SetTitle("charge vs. driftlength, all sectors, long pads");
   
   for (Int_t i = 0; i < 36; i++) {
      c1->cd(1)->SetGridx();
      c1->cd(1)->SetGridy();
      ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(i,0))->Draw();
      c1->cd(2)->SetGridx();
      c1->cd(2)->SetGridy();
      ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(i,1))->Draw();
      c1->cd(3)->SetGridx();
      c1->cd(3)->SetGridy();
      ((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(i,2))->Draw();
      c1->Update();
      chargeVsDriftlengthAllShortPads->Add( (TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,0));
      chargeVsDriftlengthAllMediumPads->Add((TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,1));
      chargeVsDriftlengthAllLongPads->Add(  (TProfile*)fcalPadRegionChargeVsDriftlength->GetObject(0,2));
      ps->NewPage();
   }
   c1->cd(1)->SetGridx();
   c1->cd(1)->SetGridy();
   chargeVsDriftlengthAllShortPads->Draw();
   c1->cd(2)->SetGridx();
   c1->cd(2)->SetGridy();
   chargeVsDriftlengthAllMediumPads->Draw();
   c1->cd(3)->SetGridx();
   c1->cd(3)->SetGridy();
   chargeVsDriftlengthAllLongPads->Draw();
   c1->Update();              // valgrind
//    ps->NewPage();
   ps->Close();  
   delete ps;
   delete c1;
   gSystem->ChangeDirectory("..");
}   



void AliTPCcalibTracks::FitResolutionNew(const char* pathName){
   // 
   // calculates different resulution fits in Y and Z direction
   // the histograms are written to 'ResolutionYZ.ps'
   // writes calculated resolution to 'resol.txt'
   // all files are stored in the directory pathName
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas c;
   c.Divide(2,1); 
   if (GetDebugLevel() > -1) cout << "creating ResolutionYZ.ps..." << endl;
   TPostScript *ps = new TPostScript("ResolutionYZ.ps", 112); 
   TF2 *fres = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
   fres->SetParameter(0,0.02);
   fres->SetParameter(1,0.0054);
   fres->SetParameter(2,0.13);  
   
   TH1::AddDirectory(kTRUE);  // TH3F::FitSlicesZ() writes histograms into the current directory
   
   // create histogramw for Y-resolution
   TH3F * hisResY0 = (TH3F*)fResolY->At(0);
   hisResY0->FitSlicesZ();
   TH2D * hisResY02 = (TH2D*)gDirectory->Get("Resol Y0_2");
   TH3F * hisResY1 = (TH3F*)fResolY->At(1); 
   hisResY1->FitSlicesZ();
   TH2D * hisResY12 = (TH2D*)gDirectory->Get("Resol Y1_2");
   TH3F * hisResY2 = (TH3F*)fResolY->At(2);
   hisResY2->FitSlicesZ();
   TH2D * hisResY22 = (TH2D*)gDirectory->Get("Resol Y2_2");
    //
   ps->NewPage();
   c.cd(1);
   hisResY02->Fit(fres, "q");      // valgrind    132,072 bytes in 6 blocks are indirectly lost
   hisResY02->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY02,fres)->Draw("surf1");
   c.Update();
   //   c.SaveAs("ResolutionYPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResY12->Fit(fres, "q");
   hisResY12->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResY12,fres)->Draw("surf1");
   c.Update();
   //   c.SaveAs("ResolutionYPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResY22->Fit(fres, "q");
   hisResY22->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY22,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionYPad2.eps");
   
   // create histogramw for Z-resolution
   TH3F * hisResZ0 = (TH3F*)fResolZ->At(0);
   hisResZ0->FitSlicesZ();
   TH2D * hisResZ02 = (TH2D*)gDirectory->Get("Resol Z0_2");
   TH3F * hisResZ1 = (TH3F*)fResolZ->At(1);
   hisResZ1->FitSlicesZ();
   TH2D * hisResZ12 = (TH2D*)gDirectory->Get("Resol Z1_2");
   TH3F * hisResZ2 = (TH3F*)fResolZ->At(2);
   hisResZ2->FitSlicesZ();
   TH2D * hisResZ22 = (TH2D*)gDirectory->Get("Resol Z2_2");
   
   ps->NewPage();
   c.cd(1);
   hisResZ02->Fit(fres, "q");
   hisResZ02->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResZ02,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ12->Fit(fres, "q");
   hisResZ12->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResZ12,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ22->Fit(fres, "q");
   hisResZ22->Draw("surf1");  
   c.cd(2);
   MakeDiff(hisResZ22,fres)->Draw("surf1");
   c.Update();
//    c.SaveAs("ResolutionZPad2.eps");
   ps->Close();
   delete ps;
   
   // write calculated resoltuions to 'resol.txt'
   ofstream fresol("resol.txt");
   fresol<<"Pad 0.75 cm"<<"\n";
   hisResY02->Fit(fres, "q");                     // valgrind
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ02->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   //
   fresol<<"Pad 1.00 cm"<<1<<"\n";
   hisResY12->Fit(fres, "q");                     // valgrind
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ12->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   //
   fresol<<"Pad 1.50 cm"<<0<<"\n";
   hisResY22->Fit(fres, "q");
   fresol<<"Y\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   hisResZ22->Fit(fres, "q");
   fresol<<"Z\t"<<fres->GetParameter(0)<<"\t"<<fres->GetParameter(1)<<"\t"<<fres->GetParameter(2)<<"\n";
   
   TH1::AddDirectory(kFALSE);
   gSystem->ChangeDirectory("..");
   delete fres;
}


void AliTPCcalibTracks::FitRMSNew(const char* pathName){
   // 
   // calculates different resulution-rms fits in Y and Z direction
   // the histograms are written to 'RMS_YZ.ps'
   // writes calculated resolution rms to 'rms.txt'
   // all files are stored in the directory pathName
   // 
   
   SetStyle();
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   
   TCanvas c;        // valgrind 3   42,120 bytes in 405 blocks are still reachable   23,816 bytes in 229 blocks are still reachable
   c.Divide(2,1); 
   if (GetDebugLevel() > -1) cout << "creating RMS_YZ.ps..." << endl;
   TPostScript *ps = new TPostScript("RMS_YZ.ps", 112); 
   TF2 *frms = new TF2("fres","TMath::Sqrt([0]*[0]+[1]*[1]*x+[2]*[2]*y*y)",0,250,0,1);
   frms->SetParameter(0,0.02);
   frms->SetParameter(1,0.0054);
   frms->SetParameter(2,0.13);  
   
   TH1::AddDirectory(kTRUE);  // TH3F::FitSlicesZ() writes histograms into the current directory
   
   // create histogramw for Y-RMS   
   TH3F * hisResY0 = (TH3F*)fRMSY->At(0);
   hisResY0->FitSlicesZ();
   TH2D * hisResY02 = (TH2D*)gDirectory->Get("RMS Y0_1");
   TH3F * hisResY1 = (TH3F*)fRMSY->At(1);
   hisResY1->FitSlicesZ();
   TH2D * hisResY12 = (TH2D*)gDirectory->Get("RMS Y1_1");
   TH3F * hisResY2 = (TH3F*)fRMSY->At(2);
   hisResY2->FitSlicesZ();
   TH2D * hisResY22 = (TH2D*)gDirectory->Get("RMS Y2_1");
   //
   ps->NewPage();
   c.cd(1);
   hisResY02->Fit(frms, "qn0"); 
   hisResY02->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY02,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResY12->Fit(frms, "qn0");               // valgrind   several blocks possibly lost
   hisResY12->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResY12,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResY22->Fit(frms, "qn0");
   hisResY22->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResY22,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSYPad2.eps");
   
   // create histogramw for Z-RMS   
   TH3F * hisResZ0 = (TH3F*)fRMSZ->At(0);
   hisResZ0->FitSlicesZ();
   TH2D * hisResZ02 = (TH2D*)gDirectory->Get("RMS Z0_1");
   TH3F * hisResZ1 = (TH3F*)fRMSZ->At(1); 
   hisResZ1->FitSlicesZ();
   TH2D * hisResZ12 = (TH2D*)gDirectory->Get("RMS Z1_1");
   TH3F * hisResZ2 = (TH3F*)fRMSZ->At(2); 
   hisResZ2->FitSlicesZ();
   TH2D * hisResZ22 = (TH2D*)gDirectory->Get("RMS Z2_1");
   //
   ps->NewPage();
   c.cd(1);
   hisResZ02->Fit(frms, "qn0");         // valgrind
   hisResZ02->Draw("surf1");
   c.cd(2);
   MakeDiff(hisResZ02,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad0.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ12->Fit(frms, "qn0");
   hisResZ12->Draw("surf1"); 
   c.cd(2);
   MakeDiff(hisResZ12,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad1.eps");
   ps->NewPage();
   c.cd(1);
   hisResZ22->Fit(frms, "qn0");         // valgrind  1 block possibly lost
   hisResZ22->Draw("surf1");  
   c.cd(2);
   MakeDiff(hisResZ22,frms)->Draw("surf1");
   c.Update();
//    c.SaveAs("RMSZPad2.eps");
   
   // write calculated resoltuion rms to 'rms.txt'
   ofstream filerms("rms.txt");
   filerms<<"Pad 0.75 cm"<<"\n";
   hisResY02->Fit(frms, "qn0");         // valgrind   23 blocks indirectly lost
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ02->Fit(frms, "qn0");         // valgrind   23 blocks indirectly lost
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   //
   filerms<<"Pad 1.00 cm"<<1<<"\n";
   hisResY12->Fit(frms, "qn0");         // valgrind      3,256 bytes in 22 blocks are indirectly lost 
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ12->Fit(frms, "qn0");         // valgrind    66,036 bytes in 3 blocks are still reachable
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   //
   filerms<<"Pad 1.50 cm"<<0<<"\n";
   hisResY22->Fit(frms, "qn0");      // valgrind   40,139 bytes in 11 blocks are still reachable   330,180 bytes in 15 blocks are possibly lost
   filerms<<"Y\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";
   hisResZ22->Fit(frms, "qn0");
   filerms<<"Z\t"<<frms->GetParameter(0)<<"\t"<<frms->GetParameter(1)<<"\t"<<frms->GetParameter(2)<<"\n";

   TH1::AddDirectory(kFALSE);
   gSystem->ChangeDirectory("..");
   ps->Close();
   delete ps;
   delete frms;
}


void AliTPCcalibTracks::MakeResPlotsQTree(Int_t minEntries, const char* pathName){
   //
   // Make tree with resolution parameters
   // the result is written to 'resol.root' in directory 'pathname'
   // file information are available in fileInfo
   // available variables in the tree 'Resol':
   //  Entries: number of entries for this resolution point
   //  nbins:   number of bins that were accumulated
   //  Dim:     direction, Dim==0: y-direction, Dim==1: z-direction
   //  Pad:     padSize; short, medium and long
   //  Length:  pad length, 0.75, 1, 1.5
   //  QMean:   mean charge of current charge bin and its neighbours, Qmean<0: integrated spectra
   //  Zc:      center of middle bin in drift direction
   //  Zm:      mean dirftlength for accumulated Delta-Histograms
   //  Zs:      width of driftlength bin
   //  AngleC:  center of middle bin in Angle-Direction
   //  AngleM:  mean angle for accumulated Delta-Histograms
   //  AngleS:  width of Angle-bin
   //  Resol:   sigma for gaus fit through Delta-Histograms
   //  Sigma:   error of sigma for gaus fit through Delta Histograms
   //  MeanR:   mean of the Delta-Histogram
   //  SigmaR:  rms of the Delta-Histogram
   //  RMSm:    mean of the gaus fit through RMS-Histogram
   //  RMS:     sigma of the gaus fit through RMS-Histogram
   //  RMSe0:   error of mean of gaus fit in RMS-Histogram
   //  RMSe1:   error of sigma of gaus fit in RMS-Histogram
   //  
      
  if (GetDebugLevel() > -1) cout << " ***** this is MakeResPlotsQTree *****" << endl;
  if (GetDebugLevel() > -1) cout << "    relax, the calculation will take a while..." << endl;
  
   gSystem->MakeDirectory(pathName);
   gSystem->ChangeDirectory(pathName);
   TString kFileName = "resol.root";
   TTreeSRedirector fTreeResol(kFileName.Data());
   
   TH3F *resArray[2][3][11];
   TH3F *rmsArray[2][3][11];
  
   // load histograms from fArrayQDY and fArrayQDZ 
   // into resArray and rmsArray
   // that is all we need here
   for (Int_t idim = 0; idim < 2; idim++){
      for (Int_t ipad = 0; ipad < 3; ipad++){
         for (Int_t iq = 0; iq <= 10; iq++){
            rmsArray[idim][ipad][iq]=0;
            resArray[idim][ipad][iq]=0;
            Int_t bin = GetBin(iq,ipad); 
            TH3F *hresl = 0;
            if (idim == 0) hresl = (TH3F*)fArrayQDY->At(bin);
            if (idim == 1) hresl = (TH3F*)fArrayQDZ->At(bin);
            if (!hresl) continue;
            resArray[idim][ipad][iq] = (TH3F*) hresl->Clone();
            resArray[idim][ipad][iq]->SetDirectory(0);
            TH3F * hreslRMS = 0;
            if (idim == 0) hreslRMS = (TH3F*)fArrayQRMSY->At(bin);
            if (idim == 1) hreslRMS = (TH3F*)fArrayQRMSZ->At(bin);
            if (!hreslRMS) continue;
            rmsArray[idim][ipad][iq] = (TH3F*) hreslRMS->Clone();
            rmsArray[idim][ipad][iq]->SetDirectory(0);
         }
      }
   }
   if (GetDebugLevel() > -1) cout << "Histograms loaded, starting to proces..." << endl;
   
   //--------------------------------------------------------------------------------------------
   
   char name[200];
   Double_t qMean;
   Double_t zMean, angleMean, zCenter, angleCenter;
   Double_t zSigma, angleSigma;
   TH1D *projectionRes = new TH1D("projectionRes", "projectionRes", 50, -1, 1);
   TH1D *projectionRms = new TH1D("projectionRms", "projectionRms", 50, -1, 1);
   TF1 *fitFunction = new TF1("fitFunction", "gaus");
   Float_t entriesQ = 0;
   Int_t loopCounter = 1;
  
   for (Int_t idim = 0; idim < 2; idim++){
      // Loop y-z corrdinate
      for (Int_t ipad = 0; ipad < 3; ipad++){
         // loop pad type
         for (Int_t iq = -1; iq < 10; iq++){
            // LOOP Q
	   if (GetDebugLevel() > -1) 
               cout << "Loop-counter, this is loop " << loopCounter << " of 66, (" 
                  << (Int_t)((loopCounter)/66.*100) << "% done), " 
                  << "idim = " << idim << ", ipad = " << ipad << ", iq = " << iq << "  \r" << std::flush;
            loopCounter++;
            TH3F *hres = 0;
            TH3F *hrms = 0;
            qMean    = 0;
            entriesQ = 0;
            
            // calculate qMean
            if (iq == -1){
               // integrated spectra
               for (Int_t iql = 0; iql < 10; iql++){    
                  Int_t bin = GetBin(iql,ipad); 
                  TH3F *hresl = resArray[idim][ipad][iql];
                  TH3F *hrmsl = rmsArray[idim][ipad][iql];
                  if (!hresl) continue;
                  if (!hrmsl) continue;	    
                  entriesQ += hresl->GetEntries();
                  qMean += hresl->GetEntries() * GetQ(bin);      
                  if (!hres) {
                     hres = (TH3F*)hresl->Clone();
                     hrms = (TH3F*)hrmsl->Clone();
                  }
                  else{
                     hres->Add(hresl);
                     hrms->Add(hrmsl);
                  }
               }
               qMean /= entriesQ;
               qMean *= -1.;  // integral mean charge
            }
            else {
               // loop over neighboured Q-bins 
               // accumulate entries from neighboured Q-bins
               for (Int_t iql = iq - 1; iql <= iq + 1; iql++){		    
                  if (iql < 0) continue;
                  Int_t bin = GetBin(iql,ipad);
                  TH3F * hresl = resArray[idim][ipad][iql];
                  TH3F * hrmsl = rmsArray[idim][ipad][iql];
                  if (!hresl) continue;
                  if (!hrmsl) continue;
                  entriesQ += hresl->GetEntries(); 
                  qMean += hresl->GetEntries() * GetQ(bin);      
                  if (!hres) {
                     hres = (TH3F*) hresl->Clone();
                     hrms = (TH3F*) hrmsl->Clone();
                  }
                  else{
                     hres->Add(hresl);
                     hrms->Add(hrmsl);
                  }
               }
               qMean/=entriesQ;
            }
      
            TAxis *xAxisDriftLength = hres->GetXaxis();   // driftlength / z - axis
            TAxis *yAxisAngle       = hres->GetYaxis();   // angle axis
            TAxis *zAxisDelta       = hres->GetZaxis();   // delta axis
            TAxis *zAxisRms         = hrms->GetZaxis();   // rms axis
            
            // loop over all angle bins
            for (Int_t ibinyAngle = 1; ibinyAngle <= yAxisAngle->GetNbins(); ibinyAngle++) {
               angleCenter = yAxisAngle->GetBinCenter(ibinyAngle);
               // loop over all driftlength bins
               for (Int_t ibinxDL = 1; ibinxDL <= xAxisDriftLength->GetNbins(); ibinxDL++) {
                  zCenter    = xAxisDriftLength->GetBinCenter(ibinxDL);
                  zSigma     = xAxisDriftLength->GetBinWidth(ibinxDL);
                  angleSigma = yAxisAngle->GetBinWidth(ibinyAngle); 
                  zMean      = zCenter;      // changens, when more statistic is accumulated
                  angleMean  = angleCenter;  // changens, when more statistic is accumulated
                  
                  // create 2 1D-Histograms, projectionRes and projectionRms
                  // these histograms are delta histograms for given direction, padSize, chargeBin,
                  // angleBin and driftLengthBin
                  // later on they will be fitted with a gausian, its sigma is the resoltuion...
                  sprintf(name,"%s, zCenter: %f, angleCenter: %f", hres->GetName(), zCenter, angleCenter);
                  // TH1D * projectionRes = new TH1D(name, name, zAxisDelta->GetNbins(), zAxisDelta->GetXmin(), zAxisDelta->GetXmax());
                  projectionRes->SetNameTitle(name, name);
                  sprintf(name,"%s, zCenter: %f, angleCenter: %f", hrms->GetName(),zCenter,angleCenter);
                  // TH1D * projectionRms =  new TH1D(name, name, zAxisDelta->GetNbins(), zAxisRms->GetXmin(), zAxisRms->GetXmax());
                  projectionRms->SetNameTitle(name, name);
                  
                  projectionRes->Reset();
                  projectionRes->SetBins(zAxisDelta->GetNbins(), zAxisDelta->GetXmin(), zAxisDelta->GetXmax());
                  projectionRms->Reset();
                  projectionRms->SetBins(zAxisRms->GetNbins(), zAxisRms->GetXmin(), zAxisRms->GetXmax());
                  projectionRes->SetDirectory(0);
                  projectionRms->SetDirectory(0);
                  
                  Double_t entries = 0;
                  Int_t    nbins   = 0;   // counts, how many bins were accumulated
                  zMean     = 0;
                  angleMean = 0;
                  entries   = 0;
                  
                  // fill projectionRes and projectionRms for given dim, ipad and iq, 
                  // as well as for given angleBin and driftlengthBin
                  // if this gives not enough statistic, include neighbourhood 
                  // (angle and driftlength) successifely
                  for (Int_t dbin = 0; dbin <= 8; dbin++){              // delta-bins around centered angleBin and driftlengthBin
                     for (Int_t dbiny2 = -1; dbiny2 <= 1; dbiny2++) {   // delta-bins in angle direction
                        for (Int_t dbinx2 = -3; dbinx2 <= 3; dbinx2++){ // delta-bins in driftlength direction
                           if (TMath::Abs(dbinx2) + TMath::Abs(dbiny2) != dbin) continue;   // add each bin only one time !
                           Int_t binx2 = ibinxDL + dbinx2;                       // position variable in x (driftlength) direction
                           Int_t biny2 = ibinyAngle + dbiny2;                    // position variable in y (angle)  direction
                           if (binx2 < 1 || biny2 < 1) continue;                 // don't go out of the histogram!
                           if (binx2 >= xAxisDriftLength->GetNbins()) continue;  // don't go out of the histogram!
                           if (biny2 >= yAxisAngle->GetNbins()) continue;        // don't go out of the histogram!
                           nbins++;                                              // count the number of accumulated bins
                           // Fill resolution histo
                           for (Int_t ibin3 = 1; ibin3 < zAxisDelta->GetNbins(); ibin3++) {
                              // Int_t content = (Int_t)hres->GetBinContent(binx2, biny2, ibin3);     // unused variable
                              projectionRes->Fill(zAxisDelta->GetBinCenter(ibin3), hres->GetBinContent(binx2, biny2, ibin3));
                              entries   += hres->GetBinContent(binx2, biny2, ibin3);
                              zMean     += hres->GetBinContent(binx2, biny2, ibin3) * xAxisDriftLength->GetBinCenter(binx2);
                              angleMean += hres->GetBinContent(binx2, biny2, ibin3) * yAxisAngle->GetBinCenter(biny2);
                           }  // ibin3 loop
                           // fill RMS histo
                           for (Int_t ibin3 = 1; ibin3 < zAxisRms->GetNbins(); ibin3++) {
                              projectionRms->Fill(zAxisRms->GetBinCenter(ibin3), hrms->GetBinContent(binx2, biny2, ibin3));
                           }
                        }  //dbinx2 loop
                        if (entries > minEntries) break; // enough statistic accumulated
                     }  // dbiny2 loop
                     if (entries > minEntries) break;    // enough statistic accumulated
                  }  // dbin loop
                  if ( entries< minEntries) continue;  // when it was absolutly impossible to get enough statistic, don't write this point into the resolution tree  
                  zMean /= entries;
                  angleMean /= entries;
                  
                  if (entries > minEntries) {
                     //  when enough statistic is accumulated
                     //  fit Delta histograms with a gausian
                     //  of the gausian is the resolution (resol), its fit error is sigma
                     //  store also mean and RMS of the histogram
                     Float_t xmin     = projectionRes->GetMean() - 2. * projectionRes->GetRMS() - 0.2;
                     Float_t xmax     = projectionRes->GetMean() + 2. * projectionRes->GetRMS() + 0.2;
                     
//                      projectionRes->Fit("gaus", "q0", "", xmin, xmax);
//                      Float_t resol    = projectionRes->GetFunction("gaus")->GetParameter(2);
//                      Float_t sigma    = projectionRes->GetFunction("gaus")->GetParError(2);
                     fitFunction->SetMaximum(xmax);
                     fitFunction->SetMinimum(xmin);
                     projectionRes->Fit("fitFunction", "qN0", "", xmin, xmax);
                     Float_t resol    = fitFunction->GetParameter(2);
                     Float_t sigma    = fitFunction->GetParError(2);
                     
                     Float_t meanR    = projectionRes->GetMean();
                     Float_t sigmaR   = projectionRes->GetRMS();
                     // fit also RMS histograms with a gausian
                     // store mean and sigma of the gausian in rmsMean and rmsSigma
                     // store also the fit errors in errorRMS and errorSigma
                     xmin = projectionRms->GetMean() - 2. * projectionRes->GetRMS() - 0.2;
                     xmax = projectionRms->GetMean() + 2. * projectionRes->GetRMS() + 0.2;
                     
//                      projectionRms->Fit("gaus","q0","",xmin,xmax);
//                      Float_t rmsMean    = projectionRms->GetFunction("gaus")->GetParameter(1);
//                      Float_t rmsSigma   = projectionRms->GetFunction("gaus")->GetParameter(2);
//                      Float_t errorRMS   = projectionRms->GetFunction("gaus")->GetParError(1);
//                      Float_t errorSigma = projectionRms->GetFunction("gaus")->GetParError(2);
                     projectionRms->Fit("fitFunction", "qN0", "", xmin, xmax);
                     Float_t rmsMean    = fitFunction->GetParameter(1);
                     Float_t rmsSigma   = fitFunction->GetParameter(2);
                     Float_t errorRMS   = fitFunction->GetParError(1);
                     Float_t errorSigma = fitFunction->GetParError(2);
                    
                     Float_t length = 0.75;
                     if (ipad == 1) length = 1;
                     if (ipad == 2) length = 1.5;
                     
                     fTreeResol<<"Resol"<<
                        "Entries="<<entries<<      // number of entries for this resolution point
                        "nbins="<<nbins<<          // number of bins that were accumulated
                        "Dim="<<idim<<             // direction, Dim==0: y-direction, Dim==1: z-direction
                        "Pad="<<ipad<<             // padSize; short, medium and long
                        "Length="<<length<<        // pad length, 0.75, 1, 1.5
                        "QMean="<<qMean<<          // mean charge of current charge bin and its neighbours, Qmean<0: integrated spectra
                        "Zc="<<zCenter<<           // center of middle bin in drift direction
                        "Zm="<<zMean<<             // mean dirftlength for accumulated Delta-Histograms
                        "Zs="<<zSigma<<            // width of driftlength bin
                        "AngleC="<<angleCenter<<   // center of middle bin in Angle-Direction
                        "AngleM="<<angleMean<<     // mean angle for accumulated Delta-Histograms
                        "AngleS="<<angleSigma<<    // width of Angle-bin
                        "Resol="<<resol<<          // sigma for gaus fit through Delta-Histograms
                        "Sigma="<<sigma<<          // error of sigma for gaus fit through Delta Histograms
                        "MeanR="<<meanR<<          // mean of the Delta-Histogram
                        "SigmaR="<<sigmaR<<        // rms of the Delta-Histogram
                        "RMSm="<<rmsMean<<         // mean of the gaus fit through RMS-Histogram
                        "RMSs="<<rmsSigma<<        // sigma of the gaus fit through RMS-Histogram
                        "RMSe0="<<errorRMS<<       // error of mean of gaus fit in RMS-Histogram
                        "RMSe1="<<errorSigma<<     // error of sigma of gaus fit in RMS-Histogram
                        "\n";
                     if (GetDebugLevel() > 5) {
                        projectionRes->SetDirectory(fTreeResol.GetFile());
                        projectionRes->Write(projectionRes->GetName());
                        projectionRes->SetDirectory(0);
                        projectionRms->SetDirectory(fTreeResol.GetFile());
                        projectionRms->Write(projectionRms->GetName());
                        projectionRes->SetDirectory(0);
                     }
                  }  // if (projectionRes->GetSum() > minEntries)
               }  // for (Int_t ibinxDL = 1; ibinxDL <= xAxisDriftLength->GetNbins(); ibinxDL++)
            }  // for (Int_t ibinyAngle = 1; ibinyAngle <= yAxisAngle->GetNbins(); ibinyAngle++)
            
         }  // iq-loop
      }  // ipad-loop
   }  // idim-loop
   delete projectionRes;
   delete projectionRms;
   
//    TFile resolFile(fTreeResol.GetFile());
   TObjString fileInfo(Form("Resolution tree, minEntries = %i", minEntries));
   fileInfo.Write("fileInfo");
//    resolFile.Close();
//    fTreeResol.GetFile()->Close();
   if (GetDebugLevel() > -1) cout << endl;
   if (GetDebugLevel() > -1) cout << "MakeResPlotsQTree done, results are in '"<< kFileName.Data() <<"'." << endl;
   gSystem->ChangeDirectory("..");
}





Long64_t AliTPCcalibTracks::Merge(TCollection *collectionList) {
   // 
   // function to merge several AliTPCcalibTracks objects after PROOF calculation
   // The object's histograms are merged via their merge functions
   // Be carefull: histograms are linked to a file, switch this off by TH1::AddDirectory(kFALSE) !!!
   // 
   
  if (GetDebugLevel() > 0) cout << " *****  this is AliTPCcalibTracks::Merge(TCollection *collectionList)  *****"<< endl;  
   if (!collectionList) return 0;
   if (collectionList->IsEmpty()) return -1;
   
   if (GetDebugLevel() > 1) cout << "the collectionList contains " << collectionList->GetEntries() << " entries." << endl;     //    REMOVE THIS LINE!!!!!!!!!!!!!!!!!1
   if (GetDebugLevel() > 5) cout << " the list in the merge-function looks as follows: " << endl;
   collectionList->Print();
   
   // create a list for each data member
   TList* deltaYList = new TList;
   TList* deltaZList = new TList;
   TList* arrayAmpRowList = new TList;
   TList* rejectedTracksList = new TList;
   TList* hclusList = new TList;
   TList* clusterCutHistoList = new TList;
   TList* arrayAmpList = new TList;
   TList* arrayQDYList = new TList;
   TList* arrayQDZList = new TList;
   TList* arrayQRMSYList = new TList;
   TList* arrayQRMSZList = new TList;
   TList* arrayChargeVsDriftlengthList = new TList;
   TList* calPadRegionChargeVsDriftlengthList = new TList;
   TList* hclusterPerPadrowList = new TList;
   TList* hclusterPerPadrowRawList = new TList;
   TList* resolYList = new TList;
   TList* resolZList = new TList;
   TList* rMSYList = new TList;
   TList* rMSZList = new TList;
   
//    TList* nRowsList = new TList;
//    TList* nSectList = new TList;
//    TList* fileNoList = new TList;
   
   TIterator *listIterator = collectionList->MakeIterator();
   AliTPCcalibTracks *calibTracks = 0;
   if (GetDebugLevel() > 1) cout << "start to iterate, filling lists" << endl;    
   Int_t counter = 0;
   while ( (calibTracks = dynamic_cast<AliTPCcalibTracks*> (listIterator->Next())) ){
      // loop over all entries in the collectionList and get dataMembers into lists
      if (!calibTracks) continue;
      
      deltaYList->Add( calibTracks->GetfDeltaY() );
      deltaZList->Add( calibTracks->GetfDeltaZ() );
      arrayAmpRowList->Add(calibTracks->GetfArrayAmpRow());
      arrayAmpList->Add(calibTracks->GetfArrayAmp());
      arrayQDYList->Add(calibTracks->GetfArrayQDY());
      arrayQDZList->Add(calibTracks->GetfArrayQDZ());
      arrayQRMSYList->Add(calibTracks->GetfArrayQRMSY());
      arrayQRMSZList->Add(calibTracks->GetfArrayQRMSZ());
      resolYList->Add(calibTracks->GetfResolY());
      resolZList->Add(calibTracks->GetfResolZ());
      rMSYList->Add(calibTracks->GetfRMSY());
      rMSZList->Add(calibTracks->GetfRMSZ());
      arrayChargeVsDriftlengthList->Add(calibTracks->GetfArrayChargeVsDriftlength());
      calPadRegionChargeVsDriftlengthList->Add(calibTracks->GetCalPadRegionchargeVsDriftlength());
      hclusList->Add(calibTracks->GetfHclus());
      rejectedTracksList->Add(calibTracks->GetfRejectedTracksHisto());
      clusterCutHistoList->Add(calibTracks->GetfClusterCutHisto());
      hclusterPerPadrowList->Add(calibTracks->GetfHclusterPerPadrow());
      hclusterPerPadrowRawList->Add(calibTracks->GetfHclusterPerPadrowRaw());
      //
      if (fCalPadClusterPerPad && calibTracks->GetfCalPadClusterPerPad())
	fCalPadClusterPerPad->Add(calibTracks->GetfCalPadClusterPerPad());      
      //      fCalPadClusterPerPadRaw->Add(calibTracks->GetfCalPadClusterPerPadRaw());
      counter++;
      if (GetDebugLevel() > 5) cout << "filling lists, object " << counter << " added." << endl;
   }
   
   
   // merge data members
   if (GetDebugLevel() > 0) cout << "histogram's merge-functins are called... " << endl; 
   fDeltaY->Merge(deltaYList);
   fDeltaZ->Merge(deltaZList);
   fHclus->Merge(hclusList);
   fClusterCutHisto->Merge(clusterCutHistoList);
   fRejectedTracksHisto->Merge(rejectedTracksList);
   fHclusterPerPadrow->Merge(hclusterPerPadrowList);
   fHclusterPerPadrowRaw->Merge(hclusterPerPadrowRawList);
   
   TObjArray* objarray = 0;
   TH1* hist = 0;
   TList* histList = 0;
   TIterator *objListIterator = 0;
   
   if (GetDebugLevel() > 0) cout << "merging fArrayAmpRows..." << endl;
   // merge fArrayAmpRows
   for (Int_t i = 0; i < fArrayAmpRow->GetEntriesFast(); i++ ) {  // loop over data member, i<72
      objListIterator = arrayAmpRowList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayAmpRowList, get TObjArray, get object at position i, cast it into TProfile
         if (!objarray) continue;
         hist = (TProfile*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TProfile*)(fArrayAmpRow->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   if (GetDebugLevel() > 0) cout << "merging fArrayAmps..." << endl;
   // merge fArrayAmps
   for (Int_t i = 0; i < fArrayAmp->GetEntriesFast(); i++ ) {  // loop over data member, i<72
      TIterator *cobjListIterator = arrayAmpList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)cobjListIterator->Next() )) { 
         // loop over arrayAmpList, get TObjArray, get object at position i, cast it into TH1F
         if (!objarray) continue;
         hist = (TH1F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH1F*)(fArrayAmp->At(i)))->Merge(histList);
      delete histList;
      delete cobjListIterator;
   }
   
   if (GetDebugLevel() > 0) cout << "merging fArrayQDY..." << endl;
   // merge fArrayQDY
   for (Int_t i = 0; i < fArrayQDY->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQDYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDYList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQDY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   if (GetDebugLevel() > 0) cout << "merging fArrayQDZ..." << endl;
   // merge fArrayQDZ
   for (Int_t i = 0; i < fArrayQDZ->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQDZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQDZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   if (GetDebugLevel() > 0) cout << "merging fArrayQRMSY..." << endl;
   // merge fArrayQRMSY
   for (Int_t i = 0; i < fArrayQRMSY->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQRMSYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQRMSY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }   

   if (GetDebugLevel() > 0) cout << "merging fArrayQRMSZ..." << endl;
   // merge fArrayQRMSZ
   for (Int_t i = 0; i < fArrayQRMSZ->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayQRMSZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fArrayQRMSZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   } 
   
   if (GetDebugLevel() > 0) cout << "merging fArrayChargeVsDriftlength..." << endl;
   // merge fArrayChargeVsDriftlength
   for (Int_t i = 0; i < fArrayChargeVsDriftlength->GetEntriesFast(); i++) { // loop over data member, i < 300
      objListIterator = arrayChargeVsDriftlengthList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TProfile
         if (!objarray) continue;
         hist = (TProfile*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TProfile*)(fArrayChargeVsDriftlength->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }    
   
   if (GetDebugLevel() > 0) cout << "merging fcalPadRegionChargeVsDriftlength..." << endl;
   // merge fcalPadRegionChargeVsDriftlength
   AliTPCCalPadRegion *cpr = 0x0;
   
   /*
   TIterator *regionIterator = fcalPadRegionChargeVsDriftlength->MakeIterator();
   while (hist = (TProfile*)regionIterator->Next()) {
      // loop over all calPadRegion's in destination calibTracks object
         objListIterator = calPadRegionChargeVsDriftlengthList->MakeIterator();
         while (( cpr =  (AliTPCCalPadRegion*)objListIterator->Next() )) { 

      
      hist->Merge(...);
   }
   */
   
   for (UInt_t isec = 0; isec < 36; isec++) {
      for (UInt_t padSize = 0; padSize < 3; padSize++){
         objListIterator = calPadRegionChargeVsDriftlengthList->MakeIterator();
         histList = new TList;
         while (( cpr =  (AliTPCCalPadRegion*)objListIterator->Next() )) { 
            // loop over calPadRegionChargeVsDriftlengthList, get AliTPCCalPadRegion, get object 
            if (!cpr) continue;
            hist = (TProfile*)cpr->GetObject(isec, padSize);
            histList->Add(hist);
         }
         ((TProfile*)(fcalPadRegionChargeVsDriftlength->GetObject(isec, padSize)))->Merge(histList);
         delete histList;
         delete objListIterator;
      }
   }
  
   
        
   
   if (GetDebugLevel() > 0) cout << "starting to merge the rest: fResolY, fResolZ , fRMSY, fRMSZ..." << endl;
   // merge fResolY
   for (Int_t i = 0; i < fResolY->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = resolYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fResolY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   // merge fResolZ
   for (Int_t i = 0; i < fResolZ->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = resolZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fResolZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }

   // merge fRMSY
   for (Int_t i = 0; i < fRMSY->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = rMSYList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fRMSY->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
         
   // merge fRMSZ
   for (Int_t i = 0; i < fRMSZ->GetEntriesFast(); i++) { // loop over data member, i < 3
      objListIterator = rMSZList->MakeIterator();
      histList = new TList;
      while (( objarray =  (TObjArray*)objListIterator->Next() )) { 
         // loop over arrayQDZList, get TObjArray, get object at position i, cast it into TH3F
         if (!objarray) continue;
         hist = (TH3F*)(objarray->At(i));
         histList->Add(hist);
      }
      ((TH3F*)(fRMSZ->At(i)))->Merge(histList);
      delete histList;
      delete objListIterator;
   }
   
   delete deltaYList;
   delete deltaZList;
   delete arrayAmpRowList;
   delete arrayAmpList;
   delete arrayQDYList;
   delete arrayQDZList;
   delete arrayQRMSYList;
   delete arrayQRMSZList;
   delete resolYList;
   delete resolZList;
   delete rMSYList;
   delete rMSZList;
   delete listIterator;
   
   if (GetDebugLevel() > 0) cout << "merging done!" << endl;
   
   return 1;
}


AliTPCcalibTracks* AliTPCcalibTracks::TestMerge(AliTPCcalibTracks *ct, AliTPCClusterParam *clusterParam, Int_t nCalTracks){
   // 
   // function to test AliTPCcalibTrack::Merge:
   // in the file 'f' is a AliTPCcalibTrack object with name "calibTracks"
   // this object is appended 'nCalTracks' times to a TList
   // A new AliTPCcalibTrack object is created which merges the list
   // this object is returned
   // 
   /*
   // .L AliTPCcalibTracks.cxx+g
   .L libTPCcalib.so
   TFile f("Output.root");
   AliTPCcalibTracks* calTracks = (AliTPCcalibTracks*)f.Get("calibTracks");
   //f.Close();
   TFile clusterParamFile("/u/lbozyk/calibration/workdir/calibTracks/TPCClusterParam.root");
   AliTPCClusterParam *clusterParam  =  (AliTPCClusterParam *) clusterParamFile.Get("Param"); 
   clusterParamFile.Close();

   AliTPCcalibTracks::TestMerge(calTracks, clusterParam);
   */
   TList *list = new TList();
   if (ct == 0 || clusterParam == 0) return 0;
   cout << "making list with " << nCalTracks << " AliTPCcalibTrack objects" << endl;
   for (Int_t i = 0; i < nCalTracks; i++) {
      if (i%10==0) cout << "Adding element " << i << " of " << nCalTracks << endl;
      list->Add(new AliTPCcalibTracks(*ct));
   }
   
   // only for check at the end
   AliTPCcalibTracks* cal1 = new AliTPCcalibTracks(*ct);
   Double_t cal1Entries = ((TH1F*)cal1->GetfArrayAmpRow()->At(5))->GetEntries();
//    Double_t cal1Entries = 5; //((TH1F*)ct->GetfArrayAmpRow()->At(5))->GetEntries();

   cout  << "The list contains " << list->GetEntries() << " entries. " << endl;
  
   
   AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);
   AliTPCcalibTracks* cal = new AliTPCcalibTracks("calTracksMerged", "calTracksMerged", clusterParam, cuts, 5);
   cal->Merge(list);
   
   cout << "cal->GetfArrayAmpRow()->At(5)->Print():" << endl;
   cal->GetfArrayAmpRow()->At(5)->Print();
   Double_t calEntries = ((TH1F*)cal->GetfArrayAmpRow()->At(5))->GetEntries();
   
   cout << "cal1->GetfArrayAmpRow()->At(5))->GetEntries() = " << cal1Entries << endl;
   cout << " cal->GetfArrayAmpRow()->At(5))->GetEntries() = " <<  calEntries << endl;
   printf("That means there were %f / %f = %f AliTPCcalibTracks-Objects merged. \n", 
      calEntries, cal1Entries, ((Double_t)calEntries/cal1Entries));
   
   return cal;

}


void  AliTPCcalibTracks::MakeQPosNormAll(TTree * chainres, AliTPCClusterParam * param, Int_t maxPoints){
  //
  // Make position corrections
  // for the moment Only using debug streamer 
  // chainres  - debug tree
  // param     - parameters to be updated
  // maxPoints - maximal number of points using for fit
  // verbose   - print info flag
  //
  // Current implementation - only using debug streamers
  // 
  
  /*    
    //Defaults
    Int_t maxPoints=100000;
  */
  /*
    Usage: 
    //0. Load libraries
    gSystem->Load("libANALYSIS");
    gSystem->Load("libSTAT");
    gSystem->Load("libTPCcalib");
    

    //1. Load Parameters to be modified:
    //e.g:
    
    AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    AliCDBManager::Instance()->SetRun(0) 
    AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();

    //2. Load chain from debug streamers
    //
    //e.g
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
    AliXRDPROOFtoolkit tool;  
    TChain * chainres = tool.MakeChain("tracks.txt","ResolCl",0,10200);
    chainres->Lookup();
    //3. Do fits and store results
    // 
    AliTPCcalibTracks::MakeQPosNormAll(chainres,param,200000,0) ;
    TFile f("paramout.root","recreate");
    param->Write("clusterParam");
    f.Close();
    //4. Verification
    TFile f2("paramout.root");
    AliTPCClusterParam *param2 = (AliTPCClusterParam*)f2.Get("clusterParam");
    param2->SetInstance(param2);
    chainres->Draw("fitZ0:AliTPCClusterParam::SPosCorrection(1,0,Cl.fPad,Cl.fTimeBin,Cl.fZ,Cl.fSigmaY2,Cl.fSigmaZ2,Cl.fMax)","Cl.fDetector<36","",10000); // should be line 
    
   */


  TStatToolkit toolkit;
  Double_t chi2;
  TVectorD fitParamY0;
  TVectorD fitParamY1;
  TVectorD fitParamZ0;
  TVectorD fitParamZ1;
  TMatrixD covMatrix;
  Int_t npoints;
  
  chainres->SetAlias("dp","(-1+(Cl.fZ>0)*2)*((Cl.fPad-int(Cl.fPad))-0.5)");
  chainres->SetAlias("dt","(-1+(Cl.fZ>0)*2)*((Cl.fTimeBin-0.66-int(Cl.fTimeBin-0.66))-0.5)");
  chainres->SetAlias("sp","(sin(dp*pi)-dp*pi)");
  chainres->SetAlias("st","(sin(dt)-dt)");
  //
  chainres->SetAlias("di","sqrt(1.-abs(Cl.fZ/250.))");
  //
  //
  //
  TCut cutA("1");
  TString fstringY="";  
  //
  fstringY+="(dp)++";            //1
  fstringY+="(dp)*di++";         //2
  fstringY+="(sp)++";            //3
  fstringY+="(sp)*di++";         //4
  TString fstringZ="";  
  fstringZ+="(dt)++";            //1
  fstringZ+="(dt)*di++";         //2
  fstringZ+="(st)++";            //3
  fstringZ+="(st)*di++";         //4
  //
  // Z corrections
  //
  TString *strZ0 = toolkit.FitPlane(chainres,"(Cl.fZ-PZ0.fElements[0]):CSigmaZ0",fstringZ.Data(), "Cl.fDetector<36"+cutA, chi2,npoints,fitParamZ0,covMatrix,-1,0,maxPoints);
  printf("Z0 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  param->fPosZcor[0] = (TVectorD*) fitParamZ0.Clone();
  //
  TString *strZ1 = toolkit.FitPlane(chainres,"(Cl.fZ-PZ0.fElements[0]):CSigmaZ0",fstringZ.Data(), "Cl.fDetector>36"+cutA, chi2,npoints,fitParamZ1,covMatrix,-1,0,maxPoints);
  //
  printf("Z1 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  param->fPosZcor[1] = (TVectorD*) fitParamZ1.Clone();
  param->fPosZcor[2] = (TVectorD*) fitParamZ1.Clone();
  //
  // Y corrections
  //   
  TString *strY0 = toolkit.FitPlane(chainres,"(Cl.fY-PY0.fElements[0]):CSigmaY0",fstringY.Data(), "Cl.fDetector<36"+cutA, chi2,npoints,fitParamY0,covMatrix,-1,0,maxPoints);
  printf("Y0 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  param->fPosYcor[0] = (TVectorD*) fitParamY0.Clone();
  

  TString *strY1 = toolkit.FitPlane(chainres,"(Cl.fY-PY0.fElements[0]):CSigmaY0",fstringY.Data(), "Cl.fDetector>36"+cutA, chi2,npoints,fitParamY1,covMatrix,-1,0,maxPoints);
  //
  printf("Y1 - chi2/npoints = %f\n",TMath::Sqrt(chi2/npoints));
  param->fPosYcor[1] = (TVectorD*) fitParamY1.Clone();
  param->fPosYcor[2] = (TVectorD*) fitParamY1.Clone();
  //
  //
  //
  chainres->SetAlias("fitZ0",strZ0->Data());
  chainres->SetAlias("fitZ1",strZ1->Data());
  chainres->SetAlias("fitY0",strY0->Data());
  chainres->SetAlias("fitY1",strY1->Data());
  //  chainres->Draw("Cl.fZ-PZ0.fElements[0]","CSigmaY0<0.7&&CSigmaZ0<0.7"+cutA,"",10000);   
}



