#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TArrayI.h>
#include "PixConv.h"
#include "AliRunLoader.h"
#include "AliGeomManager.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliITSUSDigit.h"
#include "AliITSMFTDigitPix.h"
#include "AliITSMFTSensMap.h"
#endif

PixConv converter;

void digits2raw(const char* dataDir="./",int nev=-1,int evStart=0)
{
  const Int_t kMaxROCycleAccept=126;
  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open(Form("%s/galice.root",dataDir));
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  runLoader->LoadHeader();
  runLoader->LoadDigits();
  //
  AliGeomManager::LoadGeometry(Form("%s/geometry.root",dataDir));
  AliITSUGeomTGeo::SetITSsegmentationFileName(Form("%s/%s",dataDir,AliITSUGeomTGeo::GetITSsegmentationFileName()));
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE,kTRUE);
  //
  Int_t nLayers = gm->GetNLayers();
  Int_t nChips = gm->GetNChips();  
  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //DIGITS INIT
  TTree * digTree = 0x0;
  TClonesArray *digArr= new TClonesArray("AliITSMFTDigitPix");

  int nevTot = (Int_t)runLoader->GetNumberOfEvents(), nhitTot=0;
  printf("N Events : %i \n",nevTot);
  evStart = evStart<nevTot ? evStart : nevTot-1;
  if (evStart<0) evStart = 0;
  //
  int lastEv = nev<0 ? nevTot : evStart+nev;
  if (lastEv > nevTot) lastEv = nevTot;
  //
  printf("N Events : %i \n",(Int_t)nevTot);

  //  PixConv converter;
  converter.OpenOutput("pix.bin");
  TArrayI rowcolArr(1000),sortIndArr(1000);
  int *rowcol=rowcolArr.GetArray(),*sortInd=sortIndArr.GetArray();
  //
  for (Int_t iEvent = evStart; iEvent < lastEv; iEvent++) {
    UShort_t cycleID = iEvent;
    printf("\n Event %i \n",iEvent);
    runLoader->GetEvent(iEvent);
    digTree=dl->TreeD();
    //
    digTree->SetBranchAddress("ITSDigitsPix",&digArr);
    //
    UShort_t linkCount = 0; // 1 link per chip for IB, 1 per module for OB
    UShort_t linkID=0;
    //
    for (int ilr=0;ilr<nLayers;ilr++) {
      int nst = gm->GetNStaves(ilr);
      printf("Processing %4d stave of layer %d in event %d\n",nst,ilr,iEvent);
      for (int ist=0;ist<nst;ist++) {
	int nhst = gm->GetNHalfStaves(ilr);
	for (int ihst=0;ihst<nhst;ihst++) {
	  int nmod = gm->GetNModules(ilr); // modules per halfstave
	  for (int imd=0;imd<nmod;imd++) {
	    int nch = gm->GetNChipsPerModule(ilr);
	    //
	    if (ilr>=3) { // new link for OB module
	      linkID = linkCount++;
	      converter.AddToBuffer(converter.MakeLinkHeader(linkID, cycleID),PixConv::kSizeLinkData);
	    }
	    //
	    for (int ich=0;ich<nch;ich++) {
	      //
	      int index = gm->GetChipIndex(ilr,ist,ihst,imd,ich);
	      digTree->GetEntry(index);      
	      int ndig  = digArr->GetEntries();
	      //
#ifdef _DEBUG_PIX_CONV_
	      if (ndig) printf("\n>>> Processing %3d digits of Ev:%d Lr:%d Stave:%d HStave: %d Module %d Chip:%d\n",
			       ndig,iEvent,ilr,ist,ihst,imd,ich);
#endif	    
	      //
	      if (ndig>=rowcolArr.GetSize()) {
		rowcolArr.Set(ndig+100);
		sortIndArr.Set(ndig+100);
		rowcol=rowcolArr.GetArray();
		sortInd=sortIndArr.GetArray();
	      }
	      //printf("processing ChipID%5d: Lr%d Stave:%2d HStave:%2d Mod:%2d Chip:%d | %4d hits\n",index,ilr,ist,ihst,imd,ich,ndig);
	      //
	      converter.ResetMap();
	      //
	      for (int idig=0;idig<ndig;idig++) {
		AliITSMFTDigitPix *pDig = (AliITSMFTDigitPix*)digArr->At(idig);
		rowcol[idig] = (pDig->GetCoord2()<<16)+(pDig->GetCoord1());
		//		printf("#%3d digit, col:%4d/row:%4d ROCycle:%d signal: %.2d\n",idig,pDig->GetCoord1(),pDig->GetCoord2(),pDig->GetROCycle(),pDig->GetSignalPix()); 
	      }
	      // sort
	      if (ndig) {
		TMath::Sort(ndig,rowcol,sortInd,kFALSE);
		for (int idig=0;idig<ndig;idig++) {
		  short row = rowcol[sortInd[idig]]>>16;
		  short col = rowcol[sortInd[idig]]&0xffff;
		  //printf("**AddPixel: col %d row %d\n",col,row);
		  converter.AddPixel(row,col);
		  nhitTot++;
		}
		//
	      } // << digits
	      //
	      if (ilr<3) { // new link for IB chip
		linkID = linkCount++;
		converter.AddToBuffer(converter.MakeLinkHeader(linkID, cycleID),PixConv::kSizeLinkData);
	      }
	      //
	      int chipID = ilr<3 ? 0 : ich; // IB: 1 chip/link, OB: 1module/link
	      //	      if (ndig) return;
#ifdef  _DEBUG_PIX_CONV_
	      converter.Print();
#endif
	      converter.ProcChip(chipID,iEvent,0);
	      if (ilr<3) { // close link for IB chip
		converter.AddToBuffer(converter.MakeLinkTrailer(linkID, cycleID),PixConv::kSizeLinkData);
	      }
	      //
	      converter.FlushBuffer();
	      //
#ifdef _DEBUG_PIX_CONV_
	      if (ndig) printf("\n<<< Processed  %3d digits of Ev:%d Lr:%d Stave:%d HStave: %d Module %d Chip:%d\n",
			       ndig,iEvent,ilr,ist,ihst,imd,ich);
#endif	    
	    } // chips
	    //
	    if (ilr>=3) { // close link for OB module
	      converter.AddToBuffer(converter.MakeLinkTrailer(linkID, cycleID),PixConv::kSizeLinkData);
	    }
	  } // << modules
	} // << hstaves
      } // << staves
    } // << layers
    //
    converter.FlushBuffer();
  }//event loop
  //
  printf("Converted %d hits\n",nhitTot);

}
