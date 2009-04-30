/* 
ROOT Macro to generate all numbering believed to be relevant for EMCAL
- at least as far as electronics is concerned.
Includes FEE, TRU and LED information. Sorry if this macro does not quite 
adher to ALICE coding conventions (but on the other hand it is not used in AliRoot)

Some further documentation is available at:
http://cern.ch/dsilverm/mapping/emcal_mapping.html

Author: David Silvermyr, ORNL; silvermy@mail.phy.ornl.gov
*/

/* 
First we define a number of constants; the main method EMCALNumbering starts below 
*/

const int kDDLEqIdOffsetEMCAL = 0x1200; /* From AliDAQ; first equipment Id # for EMCAL*/

// global arrays for chip, channel and CSP numbering
// - in the area covered by a single FEC (32 CSPs covering a 4x8 tower area)
const int kNROWS = 8;
const int kNCOLS = 4;
const int kNCSP = 32;

// the way CSPs are populated, as seen from the back where we
// plug in the T-cards
const int kCspMap[kNROWS][kNCOLS] = {
  0, 16,  8, 24, // 7
  1, 17,  9, 25, // |
  2, 18, 10, 26, // |
  3, 19, 11, 27, // row
  4, 20, 12, 28, // |
  5, 21, 13, 29, // |
  6, 22, 14, 30, // |
  7, 23, 15, 31  // 0
  // 0 <-col-> 3
};
// i.e. the highest row comes first, and this map should thus be indexed as [NROWS-1-irow][icol]
// - Csp help array is constructed below.

/*
  The rest of the global Chan/Chip arrays are either fixed
  from the Altro mapping, or a function of the CspMap above
*/

// Altro mapping for chips and channels, high and low gain
const int kChip[kNCSP] = {
  2,   2,   2,   2,   3,   3,   3,   3, 
  0,   0,   0,   0,   4,   4,   4,   4, 
  2,   2,   2,   2,   3,   3,   3,   3, 
  0,   0,   0,   0,   4,   4,   4,   4 
};

const int kChanHigh[kNCSP] = {
  10,  14,   5,   1,   1,   5,  14,  10, 
  10,  14,   5,   1,   1,   5,  14,  10, 
   8,  12,   7,   3,   3,   7,  12,   8, 
   8,  12,   7,   3,   3,   7,  12,   8 
};

const int kChanLow[kNCSP] = {
  11,  15,   4,   0,   0,   4,  15,  11, 
  11,  15,   4,   0,   0,   4,  15,  11, 
   9,  13,   6,   2,   2,   6,  13,   9, 
   9,  13,   6,   2,   2,   6,  13,   9 
};

// Order that CSPs appear in the data 
const int kCspOrder[kNCSP] = { // just from ALTRO mapping of chips/channels to CSP
 11,  27,  10,  26,  24,   8,  25,   9, 
  3,  19,   2,  18,  16,   0,  17,   1, 
  4,  20,   5,  21,  23,   7,  22,   6, 
 12,  28,  13,  29,  31,  15,  30,  14 
};

// LED reference info:
const int kNLED = 24; // per SuperModule; equals number of StripModules per SuperModule
const int kNLEDPerTCard = kNLED / 2;
// CSPs and LED are connected on a special T-card with only 12 connectors
// Half of the StripModules in a SuperModule will be connected to the Top
// (and half to the Bottom) T-card
// First Top
const int kCspMapLEDTop[kNLEDPerTCard] = {
  1, 17, // Strips  0, 1
  2, 18, // Strips  2, 3
  3, 19, // Strips  4, 5
  4, 20, // Strips  6, 7
  5, 21, // Strips  8, 9
  6, 22  // Strips 10,11
};
const int kStripModuleMapLEDTop[kNLEDPerTCard] = {
  0, 1,
  2, 3,
  4, 5,
  6, 7,
  8, 9,
 10,11
};

// Then Bottom
const int kCspMapLEDBottom[kNLEDPerTCard] = {
   9, 25, // Strips 12,13
  10, 26, // Strips 14,15
  11, 27, // Strips 16,17
  12, 28, // Strips 18,19
  13, 29, // Strips 20,21
  14, 30  // Strips 22,23
};
const int kStripModuleMapLEDBottom[kNLEDPerTCard] = {
  12,13,
  14,15,
  16,17,
  18,19,
  20,21,
  22,23
};

// let's make some simpler/normal help index arrays too, that we'll use later on
int ROW[kNCSP];
int COL[kNCSP];
int Csp[kNROWS][kNCOLS];

// Order that Towers, appear in the data 
int towerOrder[kNCSP];

void initTowers()
{
  for(int icol=0; icol<kNCOLS; icol++){
    for(int irow=0; irow<kNROWS; irow++){
      int csp = kCspMap[kNROWS-1-irow][icol];
      COL[csp] = icol;
      ROW[csp] = irow;
      Csp[irow][icol] = csp;

      cout << " icol " << icol
	   << " irow " << irow
	   << " csp " << csp << endl;
    }
  }

  // let's also give the order that Towers appear in the data
  for (int ic=0; ic<kNCSP; ic++) {
    int towerid = ROW[kCspOrder[ic]]*kNCOLS + COL[kCspOrder[ic]];
    towerOrder[ic] = towerid;
  }

}

// help functions for TRU mapping:
int getTRUADC(int iFEC) { // iFEC is a number 0-35, within a SM
  // ADC channel 1-12 on TRU as a function of connected iFEC
  return (iFEC%12 + 1); 
}

int getTRUADCChan(int iCSP) { // iCSP is from 0 to 31
  int bottom = (iCSP%16)/8; // 0 for top, 1 for bottom T-card
  int iADCChan = (iCSP%8)/2 + 1; // within a T-card; 1-4
  return (iADCChan + bottom*4);
}
// ok, done with TRU help methods also; let's do what needs to be done

// HERE STARTS THE MAIN METHOD..
void EMCalNumbering()
{
  /*
    General coord. info: ALICE-INT-2003-038 EDMS doc.

    z goes in the beam direction from RB26 (side C,where the muon arm is,Gex), 
    to RB24 (side A, Bellegarde), i.e. away from muon arm.

    x is horizontal, perpendicular to the beam direction and points to the 
    accelerator (LHC ring) centre. 
    [visual aid.: pos. x = Saleve; Inside/I, negative = Jura; Outside/O]

    y is vertical, perpendicular to z and x. Positive y points upward.    
    [pos. y = Up/U. neg. y = Down/D]

    The usual relation to r, phi(-pi, pi), theta (0,pi) applies:
    x = r * sin(theta) * cos(phi);
    y = r * sin(theta) * sin(phi);
    z = r * cos(theta)
  
    r = sqrt(x*x + y*y +z*z);
    theta = TMath::ACos(z/r);
    phi = TMath::ATan2(y, x);
  */

  /* Numbering rules: ALICE-INT-2003-038 EDMS doc.

  All numbering starts from 0.

  Rotational Numbering: follows phi direction
  [looks counter-clockwise from A, and clockwise from C]

  Linear numbering: follows _reverse_ z-direction from A to C, 
  without interruption at z=0
  [presumably to have a reasonable numbering in the muon system]

  Radial numbering increases outwards, but EMCAL only has one layer.
  so doesn't matter for us.
   */

  // EMCAL specifics

  initTowers(); // prepare setup

  // global info on sectors
  const int kNFullSect = 5;
  const int kNThirdSect = 1;
  const int kNSides = 2; // supermodules for both positive and negative Z
  // let's call side A '0', and side C '1'

  // Number of SuperModules and DDLs total
  const int kNSM = (kNFullSect + kNThirdSect)*kNSides; // 12
  //  const int NDDL = NSides*(NFullSect * 2 + NThirdSect*1); // 22

  // per supermodule info
  const int kNTowersZ = 48;
  const int kNTowersPhi = 24;
  const int kNTowersPerFEC = 32;
  const int kNFEC = 36;
  const int kNGTL = 4; 
  const int kNRCU = 2;
  const int kNTRU = 3;

  const int kNTowersSM = kNTowersZ*kNTowersPhi; // 1152

  const int kNModulesZ = kNTowersZ/2; // 24
  const int kNModulesPhi = kNTowersPhi/2; // 12
  const int kNModulesSM = kNTowersSM/4; // 288

  // per GTL
  const int kNFECPerGTL = kNFEC/kNGTL; // 9
  // per RCU
  const int kNFECPerRCU = kNFEC/kNRCU; // 18
  // per TRU
  const int kNFECPerTRU = kNFEC/kNTRU; // 12

  // and per FrontEndCard:
  // 2 T-cards, each with 16 towers, per FEC
  const int kNTCards = 2;
  // each T-card covers a 2x8 tower area
  const int kNTowersPhiTC = 8; 
  const int kNTowersZTC = 2; 

  // OK, that was all the setup and definitions of constants..

  /* Now, how do we populate the Supermodule and it's GTL space?
     Simplest seems to be to have 3 rows of 12 FECs each,
     where each FEC's 2 T-cards cover a 4(z)*8(phi) tower area. 
     Do the coverage in order..
  */
  const int kNFECPerRow = kNFEC/3; // 12 

  TFile *f = new TFile("map.root","RECREATE");

  // global variables
  int iside = 0; // A=0, C=1
  int isect = 0; // 0-5
  int iSM = 0;   // 0-11, offline SuperModule index
  int iDDLEqId = 0;  // kDDLEqIdOffsetEMCAL = 0x1200 upwards (NDDL) 
 
  // within SM
  int iFEC = 0;  // 0
  int iTRU = 0;  // 0-2, within SM
  int iTRUADC = 0;  // 1-8, within TRU
  int iRCU = 0;  // 0-1, within SM
  int iBranch = 0; // A=0, B=1
  int iGTL = 0;   // address 1-9; TRU is in address slot 0
  // within FEC
  int nTow = kNCSP; // # of towers, per  FEC
  int CSP[kNCSP] = {0};
  int chip[kNCSP] = {0};
  int lowGainChan[kNCSP] = {0};
  int highGainChan[kNCSP] = {0};
  int towerCol[kNCSP] = {0};
  int towerRow[kNCSP] = {0};
  int iTRUADCChan[kNCSP] = {0};

  TTree *t = new TTree("tree","ALICE EMCal tower map");
  t->Branch("iside",&iside,"iside/I");
  t->Branch("isect",&isect,"isect/I");
  t->Branch("iSM",&iSM,"iSM/I");
  t->Branch("iDDLEqId",&iDDLEqId,"iDDLEqId/I");
  t->Branch("iFEC",&iFEC,"iFEC/I");
  t->Branch("iTRU",&iTRU,"iTRU/I");
  t->Branch("iTRUADC",&iTRUADC,"iTRUADC/I");
  t->Branch("iRCU",&iRCU,"iRCU/I");
  t->Branch("iBranch",&iBranch,"iBranch/I");
  t->Branch("iGTL",&iGTL,"iGTL/I");
  t->Branch("nTow",&nTow,"nTow/I");
  t->Branch("CSP",CSP,"CSP[nTow]/I");
  t->Branch("chip",chip,"chip[nTow]/I");
  t->Branch("lowGainChan",lowGainChan,"lowGainChan[nTow]/I");
  t->Branch("highGainChan",highGainChan,"highGainChan[nTow]/I");
  t->Branch("towerCol",towerCol,"towerCol[nTow]/I");
  t->Branch("towerRow",towerRow,"towerRow[nTow]/I");
  t->Branch("towerOrder",towerOrder,"towerOrder[nTow]/I");
  t->Branch("iTRUADCChan",iTRUADCChan,"iTRUADCChan[nTow]/I");

  // LED TTree
  TTree *tLED = new TTree("tLED","ALICE EMCal LED reference map");
  tLED->Branch("iside",&iside,"iside/I");
  tLED->Branch("isect",&isect,"isect/I");
  tLED->Branch("iSM",&iSM,"iSM/I");
  tLED->Branch("iDDLEqId",&iDDLEqId,"iDDLEqId/I");
  tLED->Branch("iRCU",&iRCU,"iRCU/I");
  tLED->Branch("iBranch",&iBranch,"iBranch/I");
  tLED->Branch("iGTL",&iGTL,"iGTL/I");
  int nLED = kNLED; 
  tLED->Branch("nLED",&nLED,"nLED/I");
  int iLEDCSP[kNLED] = {0};
  int iLEDchip[kNLED] = {0};
  int iLEDhighGainChan[kNLED] = {0};
  int iLEDlowGainChan[kNLED] = {0};
  int iLEDStrip[kNLED] = {0};
  tLED->Branch("iLEDCSP",iLEDCSP,"iLEDCSP[nLED]/I");
  tLED->Branch("iLEDchip",iLEDchip,"iLEDchip[nLED]/I");
  tLED->Branch("iLEDlowGainChan",iLEDlowGainChan,"iLEDlowGainChan[nLED]/I");
  tLED->Branch("iLEDhighGainChan",iLEDhighGainChan,"iLEDhighGainChan[nLED]/I");
  tLED->Branch("iLEDStrip",iLEDStrip,"iLEDStrip[nLED]/I");

  // TRU TTree
  TTree *tTRU = new TTree("tTRU","ALICE EMCal TRU fake-altro map");
  tTRU->Branch("iside",&iside,"iside/I");
  tTRU->Branch("isect",&isect,"isect/I");
  tTRU->Branch("iSM",&iSM,"iSM/I");
  tTRU->Branch("iDDLEqId",&iDDLEqId,"iDDLEqId/I");
  tTRU->Branch("iRCU",&iRCU,"iRCU/I");
  tTRU->Branch("iBranch",&iBranch,"iBranch/I");
  tTRU->Branch("iGTL",&iGTL,"iGTL/I");
  // TRU is identified by (GTL==0 && !(Branch==0 && RCU==0))
  int iTRUFirstChan = 0; 
  int iTRULastChan = 127; // maximum allowed number of fake ALTRO channels=128 from TRU
  tTRU->Branch("iTRUFirstChan",&iTRUFirstChan,"iTRUFirstChan/I");
  tTRU->Branch("iTRULastChan",&iTRULastChan,"iTRULastChan/I");

  for (isect = 0; isect<(kNFullSect+kNThirdSect); isect++) { 
    for (iside=0; iside<kNSides; iside++) { // A or C sides
      // half sector only has one third of the FECs
      int MINFEC = 0;
      int MAXFEC = kNFEC;
      int MINTRU = 0;
      int MAXTRU = kNTRU;
      if (isect==kNFullSect) { // meaning last third-size-sector
	if (iside==0) { // A side
	  MAXFEC = kNFEC / 3;
	  MAXTRU = 1;
	}
	else if (iside==1) { // C side
	  MINFEC = 2 * kNFEC / 3;
	  MINTRU = 2;
	}
      }

      iSM = isect*2 + iside;
      for (iFEC=MINFEC; iFEC<MAXFEC; iFEC++) {

	// ok, where does this FEC belong? Use the local iFEC index which starts with 0
	// closest to the crate and revert to global z and phi index at the end..

	iTRU = iFEC / (kNFECPerTRU);
	iTRUADC = getTRUADC(iFEC); 
	iRCU = iFEC / (kNFECPerRCU); 
	iBranch = (iFEC%kNFECPerRCU) / kNFECPerGTL; // local index inside RCU
	iGTL = iFEC % kNFECPerGTL + 1;

	iDDLEqId = kDDLEqIdOffsetEMCAL + iSM*kNRCU + iRCU;

	// tower limits/indices 
	int tcolLow = (iFEC%kNFECPerRow)*kNCOLS;
	int trowLow = (iFEC/kNFECPerRow)*kNROWS;
	/*
	printf("iside %d iSM %d: FEC %02d RCU %d Branch %d iGTL\n",
	       iside, iSM, iFEC, iRCU, iBranch, iGTL);
	*/
	for (int col=0; col<kNCOLS; col++) {
	  for (int row=0; row<kNROWS; row++) {
	    int tcol = col + tcolLow;
	    int trow = row + trowLow;

	    int itow = row*kNCOLS + col;

	    CSP[itow] = Csp[row][col];
	    chip[itow] = kChip[CSP[itow]];
	    lowGainChan[itow] = kChanLow[CSP[itow]];
	    highGainChan[itow] = kChanHigh[CSP[itow]];

	    iTRUADCChan[itow] = getTRUADCChan(CSP[itow]); 

	    // we need to switch to global indices if we are on side C
	    if (iside==1) {
	      tcol =  kNTowersZ-1 - tcol; // flip axis 
	      trow = kNTowersPhi-1 - trow; // just flip axis
	    }

	    towerCol[itow] = tcol;
	    towerRow[itow] = trow;

	  } // row
	}// col
	tree->Fill();
      } // iFEC

      // also handle special LED and TRU mapping
      iGTL = 0; // they are all in GTL/FEC slot 0

      // start with LED..
      iRCU = 0; iBranch = 0; 
      iDDLEqId = kDDLEqIdOffsetEMCAL + iSM*kNRCU + iRCU;

      /* for the special 'half'/third sector, side C, there is
	 no RCU=0.. So we put the LED in RCU=1 instead then.
	 There is only 1 TRU in this sector (RCU=1, branch=1), so no problem..
      */
      if (MINTRU == 2) { // key for this special sector
	iRCU++;
	iDDLEqId++;
      }

      // loop over attached CSPs
      // First Top
      for (int iled = 0; iled<kNLEDPerTCard; iled++) {
	int il = iled;
	iLEDCSP[il] = kCspMapLEDTop[iled];
	iLEDStrip[il] = kStripModuleMapLEDTop[iled];
	/*
	cout << " LED CSP " << iLEDCSP[il]
	     << " Strip " << iLEDStrip[il] << endl;
	*/
	iLEDchip[il] = kChip[iLEDCSP[il]];
	iLEDlowGainChan[il] = kChanLow[iLEDCSP[il]];
	iLEDhighGainChan[il] = kChanHigh[iLEDCSP[il]];
      }
      // Then Bottom
      for (int iled = 0; iled<kNLEDPerTCard; iled++) {
	int il = iled + kNLEDPerTCard; // just a convenient index
	iLEDCSP[il] = kCspMapLEDBottom[iled];
	iLEDStrip[il] = kStripModuleMapLEDBottom[iled];
	/*
	cout << " LED CSP " << iLEDCSP[il]
	     << " Strip " << iLEDStrip[il] << endl;
	*/
	iLEDchip[il] = kChip[iLEDCSP[il]];
	iLEDlowGainChan[il] = kChanLow[iLEDCSP[il]];
	iLEDhighGainChan[il] = kChanHigh[iLEDCSP[il]];
      }
      tLED->Fill();

      // then we also have the TRUs (fake ALTRO)
      for (iTRU=MINTRU; iTRU<MAXTRU; iTRU++) {

	if (iTRU == 0) { iRCU = 0; iBranch = 1; } 
	else if (iTRU == 1) { iRCU = 1; iBranch = 0; } 
	else if (iTRU == 2) { iRCU = 1; iBranch = 1; } 
	/*
	cout << " TRU " << iTRU
	     << " RCU " << iRCU
	     << " Branch " << iBranch << endl;
	*/
	iDDLEqId = kDDLEqIdOffsetEMCAL + iSM*kNRCU + iRCU;
	/* 
	   last sector only has 1 TRU, but it's probably called TRU 0 on A
	   and TRU 2 on C side..
	   This is handled via MINFEC/MAXFEC and MINTRU/MAXTRU at the start
	   of this loop.
	*/
	tTRU->Fill();

      }

    } // isect / iSM
  } // iside

  f->Write();
}
