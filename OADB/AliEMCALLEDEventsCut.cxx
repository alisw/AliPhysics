
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

#include "TH2F.h"

#include "AliVEvent.h"

#include "AliLog.h"

#include "AliEMCALLEDEventsCut.h" 


/// \cond CLASSIMP
ClassImp(AliEMCALLEDEventsCut) ;
/// \endcond 

//______________________________________________
/// Default constructor. Initialize parameters.
//______________________________________________
AliEMCALLEDEventsCut::AliEMCALLEDEventsCut() : 
TObject(), fEMCalGeom(0),
fHistogramContainer(), fFillHisto(0),

// Cuts
fLEDHighEnergyCutSM(0),      fLEDHighNCellsCutSM(0), 
fLEDLowEnergyCutSM3(0),      fLEDLowNCellsCutSM3(0),
fLEDMinCellEnergy(0),        fLEDMaxCellEnergy(0),
fRemoveLEDStripEvents(0),    fLEDEventMaxNumberOfStrips(0),
fLEDLowEnergyCutSM3Strip(0), fLEDLowNCellsCutSM3Strip(0),


//Histograms
fhEMCALNSumEnCellsPerSM(0),    
fhEMCALNSumEnCellsPerSMAfter(0), 
fhEMCALNSumEnCellsPerSMAfterStripCut(0),
fhEMCALNSumEnCellsPerStrip(0), 
fhEMCALNSumEnCellsPerStripAfter(0)
{
  // Init default parameters
  //
  fLEDHighEnergyCutSM = 500.; fLEDHighNCellsCutSM = 100;
  fLEDLowEnergyCutSM3 = 20. ; fLEDLowNCellsCutSM3 = 20 ;
  fLEDMinCellEnergy = 0.5;
  fLEDMaxCellEnergy = 15.;
  
  fRemoveLEDStripEvents     = 0  ;
  fLEDEventMaxNumberOfStrips= 0  ; 
  fLEDHighEnergyCutStrip[0] = 80 ; fLEDHighEnergyCutStrip[1] = 55 ; 
  fLEDHighNCellsCutStrip[0] = 24 ; fLEDHighNCellsCutStrip[1] = 15 ;
  fLEDLowEnergyCutSM3Strip  = 100; // open
  fLEDLowNCellsCutSM3Strip  = 100; // open
}

//___________________________________________________
/// Fill the output list of initialized control histograms.
//___________________________________________________
void AliEMCALLEDEventsCut::InitControlHistograms()
{  
  printf("Init EMCal LED event control historams\n");
  fHistogramContainer = new TList();
  fHistogramContainer->SetName("EMCAL_LED") ; 
  fHistogramContainer->SetOwner(kFALSE);
  
  fhEMCALNSumEnCellsPerSM = new TH2F 
  ("hEMCALLEDEventsCut_NSumEnCellsPerSM",
   "Total number of cells and energy in any SM",
   144,0,1152,250,0,5000);
  fhEMCALNSumEnCellsPerSM->SetXTitle("#it{n}_{cells}^{SM}");
  fhEMCALNSumEnCellsPerSM->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
  fHistogramContainer->Add(fhEMCALNSumEnCellsPerSM);
  
  fhEMCALNSumEnCellsPerSMAfter = new TH2F 
  ("hEMCALLEDEventsCut_NSumEnCellsPerSMAfter",
   "Total number of cells and energy in any SM, after LED SM event rejection",
   144,0,1152,250,0,5000);
  fhEMCALNSumEnCellsPerSMAfter->SetXTitle("#it{n}_{cells}^{SM}");
  fhEMCALNSumEnCellsPerSMAfter->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
  fHistogramContainer->Add(fhEMCALNSumEnCellsPerSMAfter);
  
  if ( fRemoveLEDStripEvents )
  {
    fhEMCALNSumEnCellsPerSMAfterStripCut = new TH2F 
    ("hEMCALLEDEventsCut_NSumEnCellsPerSMAfterStripCut",
     "Total number of cells and energy in any SM, after LED SM and strip event rejection ",
     144,0,1152,250,0,5000);
    fhEMCALNSumEnCellsPerSMAfterStripCut->SetXTitle("#it{n}_{cells}^{SM}");
    fhEMCALNSumEnCellsPerSMAfterStripCut->SetYTitle("#Sigma #it{E}_{cells}^{SM} (GeV)");
    fHistogramContainer->Add(fhEMCALNSumEnCellsPerSMAfterStripCut);
  
    fhEMCALNSumEnCellsPerStrip = new TH2F 
    ("hEMCALLEDEventsCut_NSumEnCellsPerStrip",
     "Total number of cells and energy in any strip, after LED SM event rejection",
     48,0,48,100,0,500);
    fhEMCALNSumEnCellsPerStrip->SetXTitle("#it{n}_{cells}^{strip}");
    fhEMCALNSumEnCellsPerStrip->SetYTitle("#Sigma #it{E}_{cells}^{strip} (GeV)");
    fHistogramContainer->Add(fhEMCALNSumEnCellsPerStrip);
    
    fhEMCALNSumEnCellsPerStripAfter = new TH2F 
    ("hEMCALLEDEventsCut_NSumEnCellsPerStripAfter",
     "Total number of cells and energy in any strip, after LED SM event rejection",
     48,0,48,100,0,500);
    fhEMCALNSumEnCellsPerStripAfter->SetXTitle("#it{n}_{cells}^{strip}");
    fhEMCALNSumEnCellsPerStripAfter->SetYTitle("#Sigma #it{E}_{cells}^{strip} (GeV)");
    fHistogramContainer->Add(fhEMCALNSumEnCellsPerStripAfter);
  }
}

//__________________________________________
/// Count total cell energy or multiplicity per SM or per strip
/// Decide if the event is suspicious of having LED events depending on high activity on a given SM or strip.
///   * LHC11a: SM3 LED system started giving problems, check if there is large activity and remove, mostly affected last runs of LHC11a
///   * Run2: SM3 LED system was off, check, if there is low activity in this SM, remove events with large activity on any other SM
///   * A general case not depending on SM3 also considered, but not studied  
//__________________________________________
Bool_t  AliEMCALLEDEventsCut::IsEMCALLEDEvent
 (AliVEvent * event, Int_t runNumber) 
{
  // Init counters per SM
  Int_t   ncellsSM[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Float_t ecellsSM[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
  // Get geometry
  if ( !fEMCalGeom ) 
    fEMCalGeom = AliEMCALGeometry::GetInstanceFromRunNumber(runNumber);
  
  //
  // LHC11a case
  //
  if ( runNumber <= 146860 && runNumber >= 145288 )
  {
    for(Int_t icell = 0; icell < event->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      Int_t absID = event->GetEMCALCells()->GetCellNumber(icell);
      Int_t sm    = fEMCalGeom->GetSuperModuleNumber(absID);
      
      if ( event->GetEMCALCells()->GetAmplitude(icell) > 0.1 && sm == 3) ncellsSM[3]++;
    }
    
    Int_t ncellcut = 21;
    if ( (event->GetFiredTriggerClasses()).Contains("EMC") ) ncellcut = 35;
    
    if ( ncellsSM[3] >= ncellcut )
    {
      AliDebug(1,Form("Reject EMCal LED event with ncells in SM3 %d, cut %d, trig %s",
                      ncellsSM[3],ncellcut,(event->GetFiredTriggerClasses()).Data()));
      return kTRUE;
    }
  }
  //
  // Not LHC11a
  //
  else 
  {
    for(Int_t icell = 0; icell < event->GetEMCALCells()->GetNumberOfCells(); icell++)
    {
      Int_t absID = event->GetEMCALCells()->GetCellNumber(icell);
      Int_t sm    = fEMCalGeom->GetSuperModuleNumber(absID);
      Float_t amp = event->GetEMCALCells()->GetAmplitude(icell);
      
      if ( amp >= fLEDMinCellEnergy && amp <= fLEDMaxCellEnergy ) 
      {
        ncellsSM[sm]++;
        ecellsSM[sm]+=amp;
      }
    }
    
    if ( fFillHisto )
    {
      for(Int_t ism = 0; ism < 20; ism++)
        fhEMCALNSumEnCellsPerSM->Fill(ncellsSM[ism],ecellsSM[ism]);
    }
    //
    // Run2, Low activity on SM3 requiered since LED was off
    //
    if ( runNumber > 246994    ) // ADD LAST RUN WITH SM3 LED ACTIVE, right now last LHC15o
    {
      // if there is some activity in SM3, accept the event
      if ( ncellsSM[3] <= fLEDLowNCellsCutSM3 || ecellsSM[3] <= fLEDLowEnergyCutSM3 ) 
      {      
        for(Int_t ism = 0; ism < 20; ism++)
        {
          if ( ism == 3 ) continue;
          
          if ( ncellsSM[ism] >=  fLEDHighNCellsCutSM )
          {
            AliInfo(Form("Reject EMCal LED event because of SM%d: ",ism));
            for(Int_t jsm = 0; jsm < 20; jsm++){
              if ( ncellsSM[jsm] > 0 ) AliInfo(Form("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]));}
            return kTRUE;
          }
          
          if ( ecellsSM[ism] >= fLEDHighEnergyCutSM )
          {
            AliInfo(Form("Reject EMCal LED event because of SM%d: ",ism));
            for(Int_t jsm = 0; jsm < 20; jsm++) {
              if ( ncellsSM[jsm] > 0 ) AliInfo(Form("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]));}
            return kTRUE;
          }
        }
      } // SM3 activity
      
    } // Low activity on SM3
    //
    // General case without condition on SM3 low activity, NOT STUDIED!
    //
    else 
    {      
      for(Int_t ism = 0; ism < fEMCalGeom->GetNumberOfSuperModules(); ism++)
      {
        if ( ncellsSM[ism] >=  fLEDHighNCellsCutSM )
        {
          AliInfo(Form("Reject EMCal LED event because of SM%d: ",ism));
          for(Int_t jsm = 0; jsm < 20; jsm++){
            if ( ncellsSM[jsm] > 0 ) AliInfo(Form("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]));}
          return kTRUE;
        }
        
        if ( ecellsSM[ism] >= fLEDHighEnergyCutSM )
        {
          AliInfo(Form("Reject EMCal LED event because of SM%d: ",ism));
          for(Int_t jsm = 0; jsm < 20; jsm++) {
            if ( ncellsSM[jsm] > 0 ) AliInfo(Form("\t SM%d: ncells %d; sum E %3.1f \n",jsm,ncellsSM[jsm],ecellsSM[jsm]));}
          return kTRUE;
        }
      } // SM loop
    }  
    
    if ( fFillHisto )
    {
      for(Int_t ism = 0; ism < 20; ism++)
        fhEMCALNSumEnCellsPerSMAfter->Fill(ncellsSM[ism],ecellsSM[ism]);
    } 
  } 
  
  // Check activity inside all the strips 
  // (24 strips, 2x24 cells in full SM, 2x8 cellsin 1/3 SM), 
  // n cells (max 48) and sum of cells energy
  if ( fRemoveLEDStripEvents )
  {    
    Float_t amp1   = 0., amp2   = 0. ;
    Int_t   absId1 = -1, absId2 = -1 ;
    Int_t   eventNStripActiveSM[20];
    Float_t enCellsStrip[20][24];
    Int_t    nCellsStrip[20][24];
    
    for (Int_t ism = 0; ism < 20; ism++)
    {
      eventNStripActiveSM[ism] = 0;
      
      for (Int_t ieta = 0; ieta < 48; ieta=ieta+2)
      {
        enCellsStrip[ism][ieta/2] = 0.; 
        nCellsStrip [ism][ieta/2] = 0 ; 
        
        for (Int_t iphi = 0; iphi < 24; iphi++)
        {
          absId1 = fEMCalGeom->GetAbsCellIdFromCellIndexes(ism, iphi, ieta);
          if ( absId1 < 0 || absId1 > 17664 ) continue;
          
          absId2 = fEMCalGeom->GetAbsCellIdFromCellIndexes(ism, iphi, ieta+1);
          if ( absId2 < 0 || absId2 > 17664 ) continue;   
          
          amp1 = event->GetEMCALCells()->GetCellAmplitude(absId1);
          amp2 = event->GetEMCALCells()->GetCellAmplitude(absId2);
          
          if ( amp1 > fLEDMinCellEnergy && amp1 < fLEDMaxCellEnergy )
          {            
             nCellsStrip[ism][ieta/2]++;
            enCellsStrip[ism][ieta/2]+=amp1;
          }
          
          if ( amp2 > fLEDMinCellEnergy && amp2 < fLEDMaxCellEnergy )
          {            
             nCellsStrip[ism][ieta/2]++;
            enCellsStrip[ism][ieta/2]+=amp2;
          }
        }// iphi
        
        if ( fFillHisto ) fhEMCALNSumEnCellsPerStrip->Fill(nCellsStrip[ism][ieta/2],enCellsStrip[ism][ieta/2]);
        
      } // ieta
    } // ism 
    
    // Count per event over event cut
    // Low activity on SM3 for emin = 0.5
    Bool_t bSM3StripsLowActivity = kTRUE;
    for (Int_t ieta = 0; ieta < 24; ieta++)
    {
      if ( enCellsStrip[3][ieta] > fLEDLowEnergyCutSM3Strip || 
            nCellsStrip[3][ieta] > fLEDLowNCellsCutSM3Strip   ) 
        bSM3StripsLowActivity = kFALSE;
    }
    
    // Count number or active strips, depending on cuts
    //
    Int_t   nStrips = 0;
    
    if ( bSM3StripsLowActivity )
    {
      Int_t   maxNCells = 24;
      Float_t maxECells = 80;
      for (Int_t ism = 0; ism < 20; ism++)
      {
        if ( ism == 3 ) continue ;

        maxNCells = fLEDHighNCellsCutStrip[0];
        maxECells = fLEDHighEnergyCutStrip[0];
        if ( ism == 10 || ism == 11 || 
             ism == 18 || ism == 19   ) 
        {
          maxNCells = fLEDHighNCellsCutStrip[1];
          maxECells = fLEDHighEnergyCutStrip[1];
        }
  
        for (Int_t ieta = 0; ieta < 24; ieta++)
        {
          if( enCellsStrip[ism][ieta] >= maxECells || 
               nCellsStrip[ism][ieta] >= maxNCells   )
          {
            nStrips++;
          }
        } // ieta
      } // ism
    } // bSM03
    
    if ( nStrips > fLEDEventMaxNumberOfStrips ) 
    {
      AliInfo(Form("Reject EMCal LED event because of suspicious nStrips %d > %d \n",nStrips, fLEDEventMaxNumberOfStrips));
      return kTRUE;
    }
    
    if ( fFillHisto )
    {
      for (Int_t ism = 0; ism < 20; ism++)
      {
        fhEMCALNSumEnCellsPerSMAfterStripCut->Fill(ncellsSM[ism],ecellsSM[ism]);
        
        for (Int_t ieta = 0; ieta < 24; ieta++)
        {
          fhEMCALNSumEnCellsPerStripAfter->Fill(nCellsStrip[ism][ieta],enCellsStrip[ism][ieta]);
        }
      }
    }
  } // remove strip LED events
  
  return kFALSE;
}

//______________________________________________
/// Print settings and cuts
//______________________________________________
void AliEMCALLEDEventsCut::PrintParam() const
{
  printf("Remove LED events, %2.1f < Ecell < %1.2f:\n", fLEDMinCellEnergy, fLEDMaxCellEnergy  );
  printf("\t SM - nCell >= %d - Sum E >= %2.0f; \n", fLEDHighNCellsCutSM, fLEDHighEnergyCutSM);
  printf("\t SM3: nCell <= %d - Sum E <= %2.0f  \n", fLEDLowNCellsCutSM3, fLEDLowEnergyCutSM3);
  
  if ( fRemoveLEDStripEvents > 0 )
  {
    printf("Remove LED strip? %d, with n strip > %d\n: "
           "\t Full SM, nCell > %d, Sum E > %2.0f;\n "
           "\t  1/3 SM, nCell > %d, Sum E > %2.0f;\n "
           "\t     SM3, nCell < %d, Sum E < %2.0f\n",
           fRemoveLEDStripEvents    , fLEDEventMaxNumberOfStrips, 
           fLEDHighNCellsCutStrip[0], fLEDHighEnergyCutStrip[0], 
           fLEDHighNCellsCutStrip[1], fLEDHighEnergyCutStrip[1], 
           fLEDLowNCellsCutSM3Strip , fLEDLowEnergyCutSM3Strip);
  }
}
