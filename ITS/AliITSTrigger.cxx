/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

/* $Id$ */


//////////////////////////////////////////////////////////////////////
//
//  ITS SPD Trigger Detector Class
//
//  Several Fast-OR Algoritm are implemented
//
//  Ref: ALICE-INT-2005-025
//  (J. Conrad, J. G. Contreras and C. E. Jorgensen)
//
//////////////////////////////////////////////////////////////////////


#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliITSgeom.h"
#include "AliITSdigitSPD.h"

#include "AliTriggerInput.h"
#include "AliITSTrigger.h"

//______________________________________________________________________
ClassImp(AliITSTrigger)

//______________________________________________________________________
AliITSTrigger::AliITSTrigger()
  : AliTriggerDetector()
{
   SetName("ITS");
   CreateInputs();
}

//______________________________________________________________________
void AliITSTrigger::CreateInputs()
{
   // inputs

   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_GFO_L0",     "Global, Fast OR all detectors", 0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_LAYER_L0",   "Layer, hit in each layer",      0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_SECTOR_L0",  "Sector, hit coincidence in one sector", 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_HSECTOR_L0", "Half Sector, hit coincidence in one half sector", 0x08 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_SW1_L0",     "Sliding Window 1, hit coincidence on phi windows",  0x10 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_SW2_L0",     "Sliding Window 2, hit coincidence on phi windows",  0x20 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_VERTEX_L0",  "Vertex, hit coincidence on theta windows", 0x40 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_SVERTEX_L0", "Sector+Vertex, hit coincidence on one sector and theta windows", 0x80 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_HSVERTEX_L0","Half Sector+Vertex, hit coincidence on one half sector and theta windows", 0x100 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_SWVERTEX_L0","Sliding Window+Vertex, hit coincidence on phi windows and theta windows", 0x200 ) );
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_LUCUT_L0",   "Layer+Upper Cut, hit in each layer and FO occupancy cut", 0x400 ) );
   // Not yet implemented
   fInputs.AddLast( new AliTriggerInput( "ITS_SPD_HMULT_L0",   "High Multiplicity, to be defined!!", 0x600 ) );

}

//______________________________________________________________________
void AliITSTrigger::Trigger()
{

   //  ********** Get run loader for the current event **********
   AliRunLoader* runLoader = gAlice->GetRunLoader();

   AliITSLoader * loader = (AliITSLoader* )runLoader->GetLoader( "ITSLoader" );
   AliITSgeom* geom = loader->GetITSgeom();
   loader->LoadDigits("READ");
   TTree *treeD = loader->TreeD();

   TObjArray *digDet = 0;
   digDet = new TObjArray(3);

   Int_t startSPD = geom->GetStartSPD();
   Int_t lastSPD  = geom->GetLastSPD();

   // Cut on Signal In the Pixel Detector
   TBranch * br = treeD->GetBranch( "ITSDigitsSPD" );
   br->SetAddress( &((*digDet)[0]) );
   ((TClonesArray*)(digDet->At(0)))->Clear();

   Int_t ndig = 0;
   Int_t ndigfo = 0;
   Int_t singhitthreshold = 1; // single hit threshold
   Int_t threshold = 1;

//   const Int_t nlay = 2;      // not used
//   const Int_t nlad = 240;    // not used
//   const Int_t nstave = 40;   // not used
//   const Int_t ndet = 4;      // not used
//   const Int_t nchip = 5;     // not used
//   const Int_t ntotal = 1200;

   Int_t ndigA[5];
   Int_t iFOperlayer[2];
   Int_t iFOperladder[240];
   Int_t iFOperstave[40][2];
//   Int_t iFOperchip[ntotal];   // not used
   Int_t iFOperChipinStave[20][40][2];

   for (Int_t m=startSPD;m<lastSPD;m++) {
      iFOperladder[m] = 0;
   }

   for (Int_t k = 0;k<2;k++){
      iFOperlayer[k] = 0;
      for (Int_t o=0; o<40; o++) {
         iFOperstave[o][k] = 0;
         for (Int_t ich=0; ich<20; ich++) {
            iFOperChipinStave[ich][o][k] = 0;
         }
      }
   }

   // nFO = 0.0;
   Int_t mInStaveCounter = 0;
   Int_t checkStave = 0;

   for (Int_t m=startSPD; m<lastSPD; m++) {
      treeD->GetEvent( m );
      TClonesArray* digits = (TClonesArray*) (digDet->At(0)); // SPD only.
      Int_t lay, stav, det;
      geom->GetModuleId(m,lay,stav,det);
      ndig = digits->GetEntriesFast();
      for( Int_t l=0; l<5; l++ ) {
         ndigA[l] = 0 ;
      }
      for( Int_t dig=0; dig<ndig; dig++) {
         AliITSdigitSPD* dp = (AliITSdigitSPD*) digits->At(dig);
         Int_t column = dp->GetCoord1();
         //   Int_t row = dp->GetCoord2();
         Int_t chip = Int_t(column/32.);
         ndigA[chip]++;
      }

      if (checkStave != stav) {
         mInStaveCounter = 0;
      } else {
         mInStaveCounter += 1;
      }

      // m 0,.., 239
      // stav 1,..,40
      // mInStave 0,..,3
      // chipInStave 0,..,19


      for (Int_t ichip=0; ichip<5; ichip++) {
         //Int_t seq = (m*5+ichip);
         Int_t chipInStave = (mInStaveCounter *5) + ichip;

         if (ndigA[ichip] >= 1) {
            iFOperladder[m]++;
            iFOperlayer[lay-1]++;
            iFOperstave[stav-1][lay-1]++;
            //iFOperHstave[hstav-1][lay-1]++;
            iFOperChipinStave[chipInStave][stav-1][lay-1]++;
        //    nFO++;
         }
      }
      // SIMPLE FO ---> ANY HIT TRIGGERS
      ndigfo += ndig;
      checkStave = stav;
   } // endl loop over SPD's


   //Int_t bit1 = 0;
   Int_t upper_cut = 120;

   if ( ndigfo >= singhitthreshold ) {
      // Set input GLOBAL FO
      SetInput( "ITS_SPD_GFO_L0" );
      // or SetInput( "0x01" );
      // or ((AliTriggerInput*)fInputs.At(0))->Set();
      // bit1 |= (1 << 1);
   }

   //  if ( ndigfo <= upper_cut) { bit1 |= (1 << 10); }

   if (iFOperlayer[0] >= threshold && iFOperlayer[1] >=threshold) {
      // Set input LAYER
      //bit1 |= (1 << 2);
      SetInput( "ITS_SPD_LAYER_L0" );
      if ( ndigfo <= upper_cut) {
         //bit1 |= (1 << 10);
         //  LAYER and UPPER CUT
         SetInput( "ITS_SPD_LUCUT_L0" );
      }
   }

   // Sector coincidence

   // Int_t nsec = 0;
   Int_t finstav = 0;
   Int_t istav;
   Int_t jstav;

   // staves layer 1  1-20
   // staves layer 2: 0-39
   for (istav=1; istav<21; istav++) {
      for (jstav = finstav; jstav<finstav+4; jstav++) {
         if ((iFOperstave[istav-1][0] >= threshold) &&
             (iFOperstave[jstav][1] >= threshold)) {
            // Set input SECTOR
            //bit1 |= (1 << 3);
            SetInput( "ITS_SPD_SECTOR_L0" );
            if ( RequireZ10cm( iFOperChipinStave, istav-1, jstav ) == kTRUE ) {
               // Set input SECTOR & VERTEX
               // bit1 |= (1 << 7);
               SetInput( "ITS_SPD_SVERTEX_L0" );
            }
         }

      }
      if (TMath::Even(istav)) {
         finstav = jstav;
      }

   }

   // half sector coincidence

   finstav = 0;
   for (istav=1; istav<21; istav++){
      for (jstav = finstav; jstav<finstav+2; jstav++) {
         if (iFOperstave[istav-1][0] >= threshold && iFOperstave[jstav][1] >= threshold) {
            // Set input HALF SECTOR
            //bit1 |= (1 << 4);
            SetInput( "ITS_SPD_HSECTOR_L0" );
            if (RequireZ10cm(iFOperChipinStave,istav-1,jstav) == kTRUE) {
               // Set input HALF SECTOR & VERTEX
               SetInput( "ITS_SPD_HSVERTEX_L0" );
               //bit1 |= (1 << 8);
            }
         }
      }
      finstav = jstav;
   }

   finstav = 0;
   for (istav=1; istav<21; istav++) {
      for (jstav = finstav-1; jstav<finstav+3; jstav++) {
         Int_t probe_stav = jstav;
         if (jstav == -1) probe_stav = 39;
         if (jstav == 40) probe_stav = 0;
         if (iFOperstave[istav-1][0] >= threshold && iFOperstave[probe_stav][1] >= threshold) {
            // Set input SLIDING WINDOW1
            SetInput( "ITS_SPD_SW1_L0" );
            //bit1 |= (1 << 5);
         }
      }
      finstav = jstav-1;
   }

  // sliding window coincidence (symmetric): 1 (layer 1), 5 (layer 2)

   finstav = 0;
   for (istav=1;istav<21; istav++) {
      for (jstav = finstav-2; jstav<finstav+3; jstav++) {
         Int_t probe_stav = jstav;
         if (jstav == -2) probe_stav = 38;
         if (jstav == -1) probe_stav = 39;
         if (jstav == 40) probe_stav = 0;
         if (jstav == 41) probe_stav = 1;

         if ((iFOperstave[istav-1][0] >= threshold) && (iFOperstave[probe_stav][1] >= threshold)) {
            // Set input SLIDING WINDOW 2
            SetInput( "ITS_SPD_SW2_L0" );
            // bit1 |= (1 << 6);
            if (RequireZ10cm(iFOperChipinStave,istav-1,probe_stav) == kTRUE) {
               // Set input SLIDING WINDOW & VERTEX
               SetInput( "ITS_SPD_SWVERTEX_L0" );
               //bit1 |= (1 << 9);
            }
         }
      }
      finstav = jstav-1;
   }

}



//______________________________________________________________________
Bool_t AliITSTrigger::RequireZ10cm(Int_t iFOperChipinStave[][40][2], Int_t stave1, Int_t stave2)
{
   // z  sliding window
   // Bool_t zFlag = kFALSE;

   Int_t threshold = 1;
   Int_t start1 = 1;
   Int_t start2 = 3;
   Int_t start3 = 6;
   Int_t start4 = 7;
 //  Int_t i = 1;

   for (Int_t ic=0; ic<=19; ic++) {
      if (ic <= 5) {
         for (Int_t jc=0;jc<=ic-1;jc++) {
            if (iFOperChipinStave[ic][stave1][0] >=threshold) {
               if (iFOperChipinStave[jc][stave2][1] >= threshold){
                  return kTRUE; // zFlag = kTRUE;
               }
            }
         }
      }
      if (ic >= 6 && ic <= 8) {
         for (Int_t jc=(2*start1-1);jc<=(2*start1-1)+5;jc++) {
            if (iFOperChipinStave[ic][stave1][0] >=threshold) {
               if  (iFOperChipinStave[jc][stave2][1] >= threshold) {
                  return kTRUE; // zFlag = kTRUE;
               }
            }
         }
         start1++;
      }
      if (ic >= 9 && ic <= 11) {
         for (Int_t jc=(2*start2); jc<=(2*start2+5); jc++) {
            if (iFOperChipinStave[ic][stave1][0] >= threshold) {
               if (iFOperChipinStave[jc][stave2][1] >= threshold) {
                  return kTRUE; // zFlag = kTRUE;
               }
            }
         }
         start2++;
      }
      if (ic >= 12 && ic <= 13) {
         for (Int_t jc=(2*start3-1);jc<=(2*start3-1)+5;jc++) {
            if (iFOperChipinStave[ic][stave1][0] >= threshold) {
               if  (iFOperChipinStave[jc][stave2][1] >= threshold){
                  return kTRUE; // zFlag = kTRUE;
               }
            }
         }
         start3++;
      }
      if (ic >= 14) {
         for (Int_t jc=(2*start4);jc<=19;jc++) {
            if (iFOperChipinStave[ic][stave1][0] >= threshold) {
               if  (iFOperChipinStave[jc][stave2][1] >= threshold){
                  return kTRUE; // zFlag = kTRUE;
               }
            }
         }
         start4++;
      }
   }

   return kFALSE; // zFlag;
}


