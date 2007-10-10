//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Example of usage:                                                    //
// 1. Fix the seed                                                      //
// 2. Randomly store 1% dead pixels all over the SPD                    //
// 3. Print the number of dead pixels temporarily generated             //
// 4. Store the dead pixels in calibration objects for runNrs 1 to 10   //
//                                                                      //
// root [0] .L AliITSStoreDeadSPD.C                                     //
// root [1] SetSeed(11)                                                 //
// root [2] RandomizeDeadPixelsFrac(0.01)                               //
// root [3] PrintNrDead()                                               //
// root [4] StoreCalib(1,10)                                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

AliITSOnlineCalibrationSPDhandler* handler = NULL;
TRandom* rnd = new TRandom();

void StoreCalib(Int_t RunNrStart=0, Int_t RunNrEnd=9999999) {
  if (handler!=NULL) {
    handler->WriteToDB(RunNrStart,RunNrEnd);
  }
}

void ReadCalib(UInt_t runNr) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  else ClearDead();
  handler->ReadFromDB(runNr);
}

void SetSeed(UInt_t seed) {
  rnd->SetSeed(seed);
}

void PrintNrDead() {
  if (handler!=NULL) {
    printf("Nr of dead pixels: %d\n",handler->GetNrDead());
  }
}

void PrintDead() {
  if (handler!=NULL) {
    handler->PrintDead();
  }
}

void ClearDead() {
  if (handler!=NULL) {
    handler->ResetDead();
  }
}

void RandomizeDeadPixelsFrac(Double_t fraction) {
  if (fraction>0 && fraction<=1) {
    UInt_t nrdeadpixels = static_cast<UInt_t>(fraction*256*160*240 + 0.5);
    RandomizeDeadPixels(nrdeadpixels);
  }
}
void RandomizeDeadPixels(UInt_t nrDead) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t nrDeadInserted = 0;
  while (nrDeadInserted<nrDead) {
    UInt_t mod = (UInt_t)(rnd->Rndm()*240);
    UInt_t col = (UInt_t)(rnd->Rndm()*160);
    UInt_t row = (UInt_t)(rnd->Rndm()*256);
    if (handler->SetDeadPixelM(mod,col,row)) nrDeadInserted++;
  }
}

void RandomizeDeadPixelsFrac(Double_t fraction, UInt_t layer) {
  if (fraction>0 && fraction<=1) {
    UInt_t nr_pixels;
    UInt_t nrdeadpixels;
    switch(layer) {
    case 0: 
      nr_pixels = 256*160*80;
      nrdeadpixels = static_cast<UInt_t>(fraction*nr_pixels + 0.5);
      RandomizeDeadPixels(nrdeadpixels,layer);
      break;
    case 1:
      nr_pixels = 256*160*160;
      nrdeadpixels = static_cast<Int_t>(fraction*nr_pixels + 0.5);
      RandomizeDeadPixels(nrdeadpixels,layer);
      break;
    default:
      return;
    }
  }
}
void RandomizeDeadPixels(UInt_t nrDead, UInt_t layer) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t nr_modules;
  UInt_t mod_offset;
  switch(layer) {
  case 0: 
    nr_modules=80;
    mod_offset=0;
    break;
  case 1:
    nr_modules=160;
    mod_offset=80;
    break;
  default:
    return;
  }
  UInt_t nrDeadInserted = 0;
  while (nrDeadInserted<nrDead) {
    UInt_t mod = (UInt_t)(rnd->Rndm()*nr_modules) + mod_offset;
    UInt_t col = (UInt_t)(rnd->Rndm()*160);
    UInt_t row = (UInt_t)(rnd->Rndm()*256);
    if (handler->SetDeadPixelM(mod,col,row)) nrDeadInserted++;
  }
}

void RandomizeDeadChipsFrac(Double_t fraction) {
  if (fraction>0 && fraction<=1) {
    UInt_t nrdeadchips = static_cast<UInt_t>(fraction*240*5 + 0.5);
    RandomizeDeadChips(nrdeadchips);
  }
}
void RandomizeDeadChips(UInt_t nrDeadChips) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t nrDeadChipsInserted = 0;
  AliITSIntMap* chipMap = new AliITSIntMap();
  while (nrDeadChipsInserted<nrDeadChips) {
    UInt_t eq = (UInt_t)(rnd->Rndm()*20);
    UInt_t hs = (UInt_t)(rnd->Rndm()*6);
    UInt_t chip = (UInt_t)(rnd->Rndm()*10);
    if (!chipMap->Find(20*6*chip + 20*hs + eq)) {
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  handler->SetDeadPixel(eq,hs,chip,col,row); 
	}
      }
      chipMap->Insert(20*6*chip + 20*hs + eq, chip);
      nrDeadChipsInserted++;
    }
  }
  delete chipMap;
}

void RandomizeDeadChipsFrac(Double_t fraction, UInt_t layer) {
  if (fraction>0 && fraction<=1) {
    UInt_t nr_chips;
    UInt_t nrdeadchips;
    switch(layer) {
    case 0: 
      nr_chips = 400;
      nrdeadchips = static_cast<UInt_t>(fraction*nr_chips + 0.5);
      RandomizeDeadChips(nrdeadchips,layer);
      break;
    case 1:
      nr_chips = 800;
      nrdeadchips = static_cast<UInt_t>(fraction*nr_chips + 0.5);
      RandomizeDeadChips(nrdeadchips,layer);
      break;
    default:
      return;
    }
  }
}
void RandomizeDeadChips(UInt_t nrDeadChips, UInt_t layer) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t hs_nr;
  UInt_t hs_offset;
  switch(layer) {
  case 0: 
    hs_nr=2;
    hs_offset=0;
    break;
  case 1:
    hs_nr=4;
    hs_offset=2;
    break;
  default:
    return;
  }
  UInt_t nrDeadChipsInserted = 0;
  AliITSIntMap* chipMap = new AliITSIntMap();
  while (nrDeadChipsInserted<nrDeadChips) {
    UInt_t eq = (UInt_t)(rnd->Rndm()*20);
    UInt_t hs = (UInt_t)(rnd->Rndm()*hs_nr+hs_offset);
    UInt_t chip = (UInt_t)(rnd->Rndm()*10);
    if (!chipMap->Find(20*6*chip + 20*hs + eq)) {
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  handler->SetDeadPixel(eq,hs,chip,col,row); 
	}
      }
      chipMap->Insert(20*6*chip + 20*hs + eq, chip);
      nrDeadChipsInserted++;
    }
  }
  delete chipMap;
}

void RandomizeDeadHalfStavesFrac(Double_t fraction) {
  if (fraction>0 && fraction<=1) {
    UInt_t nrdeadhstaves = static_cast<UInt_t>(fraction*120 + 0.5);
    RandomizeDeadHalfStaves(nrdeadhstaves);
  }
}
void RandomizeDeadHalfStaves(UInt_t nrDeadHS) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t nrDeadHSInserted = 0;
  AliITSIntMap* hsMap = new AliITSIntMap();
  while (nrDeadHSInserted<nrDeadHS) {
    UInt_t eq = (UInt_t)(rnd->Rndm()*20);
    UInt_t hs = (UInt_t)(rnd->Rndm()*6);
    if (!hsMap->Find(20*hs + eq)) {
      for (UInt_t chip=0; chip<10; chip++) {
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=0; row<256; row++) {
	    handler->SetDeadPixel(eq,hs,chip,col,row); 
	  }
	}
      }
      hsMap->Insert(20*hs + eq, hs);
      nrDeadHSInserted++;
    }
  }
  delete hsMap;
}

void RandomizeDeadHalfStavesFrac(Double_t fraction, UInt_t layer) {
  if (fraction>0 && fraction<=1) {
    UInt_t nr_hstaves;
    switch(layer) {
    case 0: 
      nr_hstaves=static_cast<UInt_t>(fraction*40 + 0.5);
      RandomizeDeadHalfStaves(nr_hstaves,layer);
      break;
    case 1:
      nr_hstaves=static_cast<UInt_t>(fraction*80 + 0.5);
      RandomizeDeadHalfStaves(nr_hstaves,layer);
      break;
    default:
      return;
    }
  }
}
void RandomizeDeadHalfStaves(UInt_t nrDeadHS, UInt_t layer) {
  if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
  UInt_t hs_nr;
  UInt_t hs_offset;
  switch(layer) {
  case 0: 
    hs_nr=2;
    hs_offset=0;
    break;
  case 1:
    hs_nr=4;
    hs_offset=2;
    break;
  default:
    return;
  }
  UInt_t nrDeadHSInserted = 0;
  AliITSIntMap* hsMap = new AliITSIntMap();
  while (nrDeadHSInserted<nrDeadHS) {
    UInt_t eq = (UInt_t)(rnd->Rndm()*20);
    UInt_t hs = (UInt_t)(rnd->Rndm()*hs_nr+hs_offset);
    if (!hsMap->Find(20*hs + eq)) {
      for (UInt_t chip=0; chip<10; chip++) {
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=0; row<256; row++) {
	    handler->SetDeadPixel(eq,hs,chip,col,row); 
	  }
	}
      }
      hsMap->Insert(20*hs + eq, hs);
      nrDeadHSInserted++;
    }
  }
  delete hsMap;
}

void RandomizeDeadColumnChipsFrac(Double_t chipfraction, Double_t columnfraction) {
  if (chipfraction>0 && chipfraction<=1 && columnfraction>0 && columnfraction<=1) {
    if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
    UInt_t nrColChips = static_cast<UInt_t>(chipfraction*240*5 + 0.5);
    UInt_t nrColChipsUsed = 0;
    AliITSIntMap* colChipMap = new AliITSIntMap();
    while (nrColChipsUsed<nrColChips) {
      UInt_t eq = (UInt_t)(rnd->Rndm()*20);
      UInt_t hs = (UInt_t)(rnd->Rndm()*6);
      UInt_t chip = (UInt_t)(rnd->Rndm()*10);
      if (!colChipMap->Find(20*6*chip + 20*hs + eq)) {
	UInt_t nrCols = static_cast<UInt_t>(columnfraction*32 + 0.5);
	UInt_t nrColsInserted = 0;
	AliITSIntMap* colMap = new AliITSIntMap();
	while (nrColsInserted<nrCols) {
	  UInt_t col = (UInt_t)(rnd->Rndm()*32);
	  if (!colMap->Find(col)) {
	    for (UInt_t row=0; row<256; row++) {
	      handler->SetDeadPixel(eq,hs,chip,col,row);
	    }
	    colMap->Insert(col,col);
	    nrColsInserted++;
	  }
	}
	delete colMap;
	colChipMap->Insert(20*6*chip + 20*hs + eq, chip);
	nrColChipsUsed++;
      }
    }
    delete colChipMap;
  }
}
void RandomizeDeadColumnChipsFrac(Double_t chipfraction, Double_t columnfraction, Int_t layer) {
  if (chipfraction>0 && chipfraction<=1 && columnfraction>0 && columnfraction<=1) {
    if (handler==NULL) handler = new AliITSOnlineCalibrationSPDhandler();
    UInt_t hs_nr;
    UInt_t hs_offset;
    switch(layer) {
    case 0: 
      hs_nr=2;
      hs_offset=0;
      break;
    case 1:
      hs_nr=4;
      hs_offset=2;
      break;
    default:
      return;
    }
    UInt_t nrColChips = static_cast<UInt_t>(chipfraction*240*5 + 0.5);
    UInt_t nrColChipsUsed = 0;
    AliITSIntMap* colChipMap = new AliITSIntMap();
    while (nrColChipsUsed<nrColChips) {
      UInt_t eq = (UInt_t)(rnd->Rndm()*20);
      UInt_t hs = (UInt_t)(rnd->Rndm()*hs_nr+hs_offset);
      UInt_t chip = (UInt_t)(rnd->Rndm()*10);
      if (!colChipMap->Find(20*6*chip + 20*hs + eq)) {
	UInt_t nrCols = static_cast<UInt_t>(columnfraction*32 + 0.5);
	UInt_t nrColsInserted = 0;
	AliITSIntMap* colMap = new AliITSIntMap();
	while (nrColsInserted<nrCols) {
	  UInt_t col = (UInt_t)(rnd->Rndm()*32);
	  if (!colMap->Find(col)) {
	    for (UInt_t row=0; row<256; row++) {
	      handler->SetDeadPixel(eq,hs,chip,col,row);
	    }
	    colMap->Insert(col,col);
	    nrColsInserted++;
	  }
	}
	delete colMap;
	colChipMap->Insert(20*6*chip + 20*hs + eq, chip);
	nrColChipsUsed++;
      }
    }
    delete colChipMap;
  }
}

