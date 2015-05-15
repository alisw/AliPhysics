/*
ROOT Macro to generate ascii files for AliCaloAltroMapping; offline decoder
- based on the map file produced by EMCALNumbering.C

Update: 20 Aug 2008 - update to write separate files for A and C sides +
   added TRU and LED channels also.. 

LED channels look like regular FEE data; just from rcu=0, branch=0, FEC=0
TRU data has one block of data from the whole TRU; all other FEC=0 slots
   (rcu, branch) = (0,1), (1,0), (1,1)

Some further documentation is available at:
http://cern.ch/dsilverm/mapping/emcal_mapping.html

Author: David Silvermyr, ORNL; silvermy@mail.phy.ornl.gov
*/

// First we define some constants - the main method WriteRCUs comes later
const int kNTRU = 3; // per SM
const int kNTRUChanBlocks = 128; // max. per TRU 
const int kNLED = 24; // one per StripModule

const int kNGAIN = 2; // low (0) and high (1)
const int kNFEEChannelsPerRCU = 1152;
int NChannelsPerRCU[2];
// RCU 0: FEE + NTRUChanBlocks + NLED*NGAIN 
NChannelsPerRCU[0] = kNFEEChannelsPerRCU + kNTRUChanBlocks + kNLED*kNGAIN; 
// RCU 1: FEE + 2*NTRUChanBlocks
NChannelsPerRCU[1] = kNFEEChannelsPerRCU + 2*kNTRUChanBlocks; 

const int kMaxHWAddress = (1 << 11) // branch
  | (9 << 7) // FEC
  | (4 << 4) // Altro
  | 0xf; // channel

const int kNRCU = 2;
const int kNSides = 2;
char * sideStr[] = {"A", "C"};

const int kFirstFEC = 1;
const int kNFECPerGTL = 9; // NFEC/NGTL

int makeHWAddress(int ibranch, int ifec, int ichip, int ichan) {
  int addr = (ibranch << 11) // branch
    | (ifec << 7) // FEC
    | (ichip << 4) // Altro
    | ichan; // channel
  return addr;
}

// HERE STARTS THE MAIN METHOD..
void WriteRCUs(const char *filename="map.root")
{
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!f) {
    f = new TFile(filename);
  }
  TTree *tree = (TTree*)gDirectory->Get("tree");

  const int kNTOW = 32; // per FEC
  
  //Declaration of leaves types
  Int_t           iside;
  Int_t           isect;
  Int_t           iSM;
  Int_t           iFEC;
  Int_t           iRCU;
  Int_t           iDDLEqId;
  Int_t           iBranch;
  Int_t           iGTL;
  Int_t           nTow;
  Int_t           CSP[kNTOW];
  Int_t           chip[kNTOW];
  Int_t           lowGainChan[kNTOW];
  Int_t           highGainChan[kNTOW];
  Int_t           towerCol[kNTOW];
  Int_t           towerRow[kNTOW];
  Int_t           towerOrder[kNTOW];
  
  // Set branch addresses.
  tree->SetBranchAddress("iside",&iside);
  tree->SetBranchAddress("isect",&isect);
  tree->SetBranchAddress("iSM",&iSM);
  tree->SetBranchAddress("iFEC",&iFEC);
  tree->SetBranchAddress("iRCU",&iRCU);
  tree->SetBranchAddress("iDDLEqId",&iDDLEqId);
  tree->SetBranchAddress("iBranch",&iBranch);
  tree->SetBranchAddress("iGTL",&iGTL);
  tree->SetBranchAddress("nTow",&nTow);
  tree->SetBranchAddress("CSP",CSP);
  tree->SetBranchAddress("chip",chip);
  tree->SetBranchAddress("lowGainChan",lowGainChan);
  tree->SetBranchAddress("highGainChan",highGainChan);
  tree->SetBranchAddress("towerCol",towerCol);
  tree->SetBranchAddress("towerRow",towerRow);
  tree->SetBranchAddress("towerOrder",towerOrder);
  
  Int_t nbytes = 0;

  // also variables and branch setting for LED TTree
  TTree *tLED = (TTree*)gDirectory->Get("tLED");

  //Declaration of leaves types, not already declared
  Int_t           nLED;
  Int_t           iLEDCSP[kNLED];
  Int_t           iLEDchip[kNLED];
  Int_t           iLEDlowGainChan[kNLED];
  Int_t           iLEDhighGainChan[kNLED];
  Int_t           iLEDStrip[kNLED];

  // Set branch addresses.
  tLED->SetBranchAddress("iside",&iside);
  tLED->SetBranchAddress("isect",&isect);
  tLED->SetBranchAddress("iSM",&iSM);
  tLED->SetBranchAddress("iDDLEqId",&iDDLEqId);
  tLED->SetBranchAddress("iRCU",&iRCU);
  tLED->SetBranchAddress("iBranch",&iBranch);
  tLED->SetBranchAddress("iGTL",&iGTL);
  tLED->SetBranchAddress("nLED",&nLED);
  tLED->SetBranchAddress("iLEDCSP",iLEDCSP);
  tLED->SetBranchAddress("iLEDchip",iLEDchip);
  tLED->SetBranchAddress("iLEDlowGainChan",iLEDlowGainChan);
  tLED->SetBranchAddress("iLEDhighGainChan",iLEDhighGainChan);
  tLED->SetBranchAddress("iLEDStrip",iLEDStrip);

  // and TRU TTree
  TTree *tTRU = (TTree*)gDirectory->Get("tTRU");

  //Declaration of leaves types, not already declared
  Int_t           iTRUchip;
  Int_t           iTRUFirstChan;
  Int_t           iTRULastChan;
  
   // Set branch addresses.
  tTRU->SetBranchAddress("iside",&iside);
  tTRU->SetBranchAddress("isect",&isect);
  tTRU->SetBranchAddress("iSM",&iSM);
  tTRU->SetBranchAddress("iDDLEqId",&iDDLEqId);
  tTRU->SetBranchAddress("iRCU",&iRCU);
  tTRU->SetBranchAddress("iBranch",&iBranch);
  tTRU->SetBranchAddress("iGTL",&iGTL);
  tTRU->SetBranchAddress("iTRUFirstChan",&iTRUFirstChan);
  tTRU->SetBranchAddress("iTRULastChan",&iTRULastChan);

   // OK, setup done - let's proceed.. First with the regular data

  int n[kNRCU][kNSides] = {0};
  ofstream out[kNRCU][kNSides];
  ofstream outFinal[kNRCU][kNSides];

  // Open output files, and provide the necessary header
  char fname[100];
  for (iRCU=0; iRCU<kNRCU; iRCU++) {
    for (iside=0; iside<kNSides; iside++) {
      sprintf(fname, "RCU%d%s.data.unsorted",iRCU,sideStr[iside]);
      out[iRCU][iside].open(fname);

      sprintf(fname, "RCU%d%s.data",iRCU,sideStr[iside]);
      outFinal[iRCU][iside].open(fname);

      outFinal[iRCU][iside] << NChannelsPerRCU[iRCU] << endl;
      outFinal[iRCU][iside] << kMaxHWAddress  << endl;

      n[iRCU][iside] = 0;
    }
  }

  int chid[kNGAIN] = {0};
  int ig = 0; // gain flag

  // loop over channels
  for (Long64_t i=0; i<tree->GetEntries();i++) { 
    nbytes += tree->GetEntry(i);

    if (isect==0) { // just the 1st sector; other sectors will look the same internally (or be a subset in case of the 1/3 size ones)

      /*
	FECs should come in order so let's just worry about channel order
	within an FEC
	DS (Aug 2008): in principle we don't really need to do this anymore, 
	since we anyway call a sort command at the end now, but since it was
	already coded I kept this for now..
      */
      // The towerOrder array should have the needed indexing info
      for (int it = 0; it<kNTOW; it++) {
	int itow = towerOrder[it];

	chid[0] = makeHWAddress( iBranch, (iFEC%kNFECPerGTL) + kFirstFEC, 
				 chip[itow], lowGainChan[itow] );
	chid[1] = makeHWAddress( iBranch, (iFEC%kNFECPerGTL) + kFirstFEC, 
				 chip[itow], highGainChan[itow] );

	// cout << " it " << it << " : " << chid[0] << " " << chid[1] << endl;

	// start with the lowest
	if (chid[0] < chid[1]) {
	  ig = 0;
	}
	else {
	  ig = 1;
	}
	
	out[iRCU][iside] << chid[ig] << " " 
			 << towerRow[itow] << " "
			 << towerCol[itow] << " " << ig
			 << endl;
	n[iRCU][iside]++;
	// and then the next
	ig = 1 -ig;
	out[iRCU][iside] << chid[ig] << " " 
			 << towerRow[itow] << " "
			 << towerCol[itow] << " " << ig
			 << endl;
	n[iRCU][iside]++;
      } // loop over towers      
    } // sector==0 selection
  }
  // done, with FEE data

  // report on counts; just for early debugging..
  for (iRCU=0; iRCU<kNRCU; iRCU++) {
    for (iside=0; iside<kNSides; iside++) {
      cout << " post-FEE count: RCU " << iRCU
	   << " iside " << iside
	   << " n = " << n[iRCU][iside] << endl;
    }
  }
  
  // Next, we have the fake ALTRO from the TRU
  int caloFlag = 2; // from enum AliCaloRawStream::kTRUData
  int dummyRow = 0; // need to have the right number of fields in print-out
  int TRUchid = 0;
  for (Long64_t i=0; i<tTRU->GetEntries(); i++) { 
    nbytes += tTRU->GetEntry(i);
    if (isect==0) { // select just the 1st sector; same motivation as for FEE
      for (int ichan=iTRUFirstChan; ichan<=iTRULastChan; ichan++) {
	TRUchid = makeHWAddress( iBranch, iGTL, 
				 ichan/16, ichan%16 );
	out[iRCU][iside] << TRUchid << " " 
			 << dummyRow << " "
			 << ichan << " " // channel # coded as Column.. 
			 << caloFlag << endl;
	n[iRCU][iside]++;
      }
    }
  }

  // report on counts; just for early debugging..
  for (iRCU=0; iRCU<kNRCU; iRCU++) {
    for (iside=0; iside<kNSides; iside++) {
      cout << " post-TRU count: RCU " << iRCU
	   << " iside " << iside
	   << " n = " << n[iRCU][iside] << endl;
    }
  }

  // Then. LED data
  caloFlag = 3; // from enum AliCaloRawStream::kLEDMonData
  int LEDchid = 0; // assigned below

  for (Long64_t i=0; i<tLED->GetEntries(); i++) { 
    nbytes += tLED->GetEntry(i);
    if (isect==0) { // just the 1st sector
      for (int il=0; il<nLED; il++) {

	// first low gain
	ig = 0; 
	LEDchid = makeHWAddress( iBranch, iGTL, 
				 iLEDchip[il], iLEDlowGainChan[il] );
	out[iRCU][iside] << LEDchid << " " 
			 << ig << " " // gain coded as row..
			 << iLEDStrip[il] << " " // strip # coded as Column.. 
			 << caloFlag << endl;
	n[iRCU][iside]++;
	
	// then high gain
	ig = 1; 
	LEDchid = makeHWAddress( iBranch, iGTL, 
				 iLEDchip[il], iLEDhighGainChan[il] );
	out[iRCU][iside] << LEDchid << " " 
			 << ig << " " // gain coded as row..
			 << iLEDStrip[il] << " " // strip # coded as Column.. 
			 << caloFlag << endl;
	n[iRCU][iside]++;
      }
    }
  }

  // Finally, let's close the output files - we should get the counts we expected
  for (iRCU=0; iRCU<kNRCU; iRCU++) {
    for (iside=0; iside<kNSides; iside++) {
      out[iRCU][iside].close();
      outFinal[iRCU][iside].close();
      // report how many entries we encountered
      cout << " final count: RCU " << iRCU
	   << " iside " << iside
	   << " n = " << n[iRCU][iside] 
	   << " expected " << NChannelsPerRCU[iRCU];
      if (n[iRCU][iside] == NChannelsPerRCU[iRCU]) {
	cout << " - OK! " << endl;
      }
      else {
	cout << " - Wrong! " << endl;
      }

    }
  }

  // let's then do a final loop where we sort the lists into the final files
  // I'm not sure if having the lists ordered is really needed but it is at 
  // least more aestethically pleasing..
  char cmd[200];
  sprintf(cmd, "mkdir -p tmp"); // provide a temporary storage space - not to pollute the local dir. too much
  gSystem->Exec(cmd);

  for (iRCU=0; iRCU<kNRCU; iRCU++) {
    for (iside=0; iside<kNSides; iside++) {
      sprintf(cmd, "sort -n RCU%d%s.data.unsorted >> RCU%d%s.data", 
	      iRCU, sideStr[iside], iRCU, sideStr[iside]);
      cout << "executing " << cmd << endl;
      gSystem->Exec(cmd);
      sprintf(cmd, "mv RCU%d%s.data.unsorted tmp/.", 
	      iRCU, sideStr[iside]);
      cout << "executing " << cmd << endl;
      gSystem->Exec(cmd);
    }
  }

}
