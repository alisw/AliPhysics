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

//_________________________________________________________________________
// Cheking PHOSHistos procedure of PHOS
//*-- Author : Gines MARTINEZ  SUBATECH january 2000
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TTree.h"

// --- Standard library ---

#include <iostream>
#include <cstdio>


// --- AliRoot header files ---
#include "AliRun.h"
#include "TFile.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSv0.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "PHOSHistos.h"


void PHOSHistos (Text_t* infile, Int_t nevent, Int_t Module)  
{
//========== Opening galice.root file  
  TFile * file   = new TFile(infile); 
//========== Get AliRun object from file 
  gAlice = (AliRun*) file->Get("gAlice");//=========== Gets the PHOS object and associated geometry from the file 
  AliPHOSv0 * PHOS  = (AliPHOSv0 *)gAlice->GetDetector("PHOS");
  AliPHOSGeometry * Geom = AliPHOSGeometry::GetInstance(PHOS->GetGeometry()->GetName(),PHOS->GetGeometry()->GetTitle());
//========== Creates the Clusterizer
  AliPHOSClusterizerv1 clusterizer; 
  clusterizer.SetEmcEnergyThreshold(0.01) ; 
  clusterizer.SetEmcClusteringThreshold(0.1) ; 
  clusterizer.SetPpsdEnergyThreshold(0.00000005) ; 
  clusterizer.SetPpsdClusteringThreshold(0.0000005) ; 
  clusterizer.SetLocalMaxCut(0.03) ;
  clusterizer.SetCalibrationParameters(0., 0.00000001) ;
//========== Creates the track segment maker
  AliPHOSTrackSegmentMakerv1 tracksegmentmaker ;
//========== Creates the Reconstructioner  
  AliPHOSReconstructioner Reconstructioner(clusterizer,tracksegmentmaker);     
    
  TH1F * hEmcDigit       = new TH1F("hEmcDigit","hEmcDigit",1000,0.,10.);
  TH1F * hVetoDigit      = new TH1F("hVetoDigit","hVetoDigit",1000,0.,3.e-5);
  TH1F * hConvertorDigit = new TH1F("hConvertorDigit","hConvertorDigit",1000,0.,3.e-5);
  TH1F * hEmcCluster       = new TH1F("hEmcCluster","hEmcCluster",100,0.,10.);
  TH1F * hVetoCluster      = new TH1F("hVetoCluster","hVetoCluster",1000,0.,3.e-5);
  TH1F * hConvertorCluster = new TH1F("hConvertorCluster","hConvertorCluster",1000,0.,3.e-5);
  AliPHOSDigit * digit ;

//========== Loop on events
  Int_t ievent;
  for(ievent=0;ievent<nevent; ievent++)
  {
    //   cout << "Event " << ievent <<endl;

    Int_t RelId[4] ;
    //=========== Connects the various Tree's for evt
    gAlice->GetEvent(ievent);
    //=========== Gets the Digit TTree
    gAlice->TreeD()->GetEvent(0) ;     
    //=========== Gets the number of entries in the Digits array
    //    Int_t nId = PHOS->Digits()->GetEntries();      
    TIter next(PHOS->Digits()) ;
    Float_t Etot=0 ;
    Int_t nVeto=0 ;
    Int_t nConvertor=0 ;
    while( ( digit = (AliPHOSDigit *)next() ) )
    {
       Etot+=clusterizer.Calibrate(digit->GetAmp()) ;
       Geom->AbsToRelNumbering(digit->GetId(), RelId) ;        

       if (clusterizer.IsInEmc(digit))
       {   hEmcDigit->Fill(clusterizer.Calibrate(digit->GetAmp())) ; }
       else    
       {  
          if (RelId[1]==9) {nVeto++; hVetoDigit->Fill(clusterizer.Calibrate(digit->GetAmp()));} 
	  if (RelId[1]==25){nConvertor++; hConvertorDigit->Fill(clusterizer.Calibrate(digit->GetAmp()));}
       }
    }
     
    //    if (nVeto>1)       printf("AnaPHOSv0.C> Number of Veto entries  is %d \n",nVeto);
    //    if (nConvertor>1)  printf("AnaPHOSv0.C> Number of Convertor entries is %d \n",nConvertor);

//    cout <<"TestRec> Found  " << nId << "  digits in PHOS with total energy " << Etot << endl ;
//  
         PHOS->Reconstruction(Reconstructioner); 
        //  PHOS->EmcClusters()->Delete();
//          PHOS->PpsdClusters()->Delete(); 
//   
//    //=========== Cluster in Module
//     TClonesArray * EmcRP = PHOS->EmcClusters() ;
//     Etot = 0.; 
//     Int_t TotalNumberOfClusters = 0 ; 
//     Int_t NumberOfClusters = 0 ;
//     TIter nextemc(EmcRP) ;
//     AliPHOSEmcRecPoint * emc ;
//     while((emc = (AliPHOSEmcRecPoint *)nextemc())) 
//     {
//       TotalNumberOfClusters++ ;
//       if ( emc->GetPHOSMod() == Module )
//       { 
//         NumberOfClusters++ ; 
//         Energy = emc->GetTotalEnergy() ;
// 	hEmcCluster->Fill(Energy);   
//         Etot+=Energy ;  
//       }
//     }
//     cout << "TestRec> Found " << TotalNumberOfClusters << " EMC Clusters in PHOS" << endl ; 
//     cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " EMC Clusters " << endl ;
//     cout << "TestRec> Total energy  " <<Etot << endl ; 
// 
//    //=========== Cluster in Module PPSD Down
//     TClonesArray * PpsdRP = PHOS->PpsdClusters() ;
//     Etot = 0.; 
//     Int_t TotalNumberOfClusters = 0 ; 
//     Int_t NumberOfClusters = 0 ;
//     TIter nextPpsd(PpsdRP) ;
//     AliPHOSPpsdRecPoint * Ppsd ;
//     while((Ppsd = (AliPHOSPpsdRecPoint *)nextPpsd())) 
//     {
//       TotalNumberOfClusters++ ;
//       if ( Ppsd->GetPHOSMod() == Module )
//       { 
//         NumberOfClusters++ ; 
//         Energy = Ppsd->GetEnergy() ;
// 	hConvertorCluster->Fill(Energy) ;   
//         Etot+=Energy ;  
//         if (!Ppsd->GetUp()) Ppsd->Draw("P") ;
//       }
//     }
//     cout << "TestRec> Found " << TotalNumberOfClusters << " Ppsd Down Clusters in PHOS" << endl ; 
//     cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Down Clusters " << endl ;
//     cout << "TestRec> Total energy  " <<Etot << endl ; 
// 
//    //=========== Cluster in Module PPSD Up
//     PpsdRP = PHOS->PpsdClusters() ;
//     Etot = 0.; 
//     Int_t TotalNumberOfClusters = 0 ; 
//     Int_t NumberOfClusters = 0 ;
//     TIter nextPpsdUp(PpsdRP) ;
//     while((Ppsd = (AliPHOSPpsdRecPoint *)nextPpsdUp())) 
//     {
//       TotalNumberOfClusters++ ;
//       if ( Ppsd->GetPHOSMod() == Module )
//       { 
//         NumberOfClusters++ ; 
//         Energy = Ppsd->GetEnergy() ;
// 	hVetoCluster->Fill(Energy);   
//         Etot+=Energy ;  
//         if (Ppsd->GetUp()) Ppsd->Draw("P") ;
//       }
//     }
//     cout << "TestRec> Found " << TotalNumberOfClusters << " Ppsd Up Clusters in PHOS" << endl ; 
//     cout << "TestRec> Found in Module " << Module << "  " << NumberOfClusters << " Ppsd Up Clusters " << endl ;
//     cout << "TestRec> Total energy  " <<Etot << endl ; 

  } 
  TCanvas * cVetoDigit = new TCanvas("VetoDigit","VetoDigit");  
  hVetoDigit->Draw();
  TCanvas * cConvertorDigit = new TCanvas("ConvertorDigit","ConvertorDigit");  
  hConvertorDigit->Draw();
  TCanvas * cEmcDigit = new TCanvas("EmcDigit","EmcDigit");  
  hEmcDigit->Draw();
  
  
}


