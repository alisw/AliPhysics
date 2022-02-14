#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <AliOADBContainer.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TProfile.h>
#include <TSystem.h>
#include "LInfo.h"
#include "TInfo.h"
#include "plotOCDB_Temperature.C"
#include "plotOCDB_LED.C"

class TCalCell : public TObject {
 public:
  TCalCell() : fId(0), fSM(0), fLedM(0), fLedR(0), fMonM(0), fMonR(0), fLocT(0), fSMT(0) {;}
  virtual ~TCalCell() {;}
  Short_t    fId;   //         cell id
  Short_t    fSM;   //         super module index
  Short_t    fRow;  //         row (phi) index
  Short_t    fCol;  //         col (eta) index
  Short_t    fBad;  //         bad cell
  Double32_t fLedM; //[0,0,16] led mean
  Double32_t fLedR; //[0,0,16] led rms
  Double32_t fMonM; //[0,0,16] mon mean
  Double32_t fMonR; //[0,0,16] mon rms
  Double32_t fLocT; //[0,0,16] loc T
  Double32_t fSMT;  //[0,0,16] sm T
  ClassDef(TCalCell, 2); // CalCell class
};

class TCalSM : public TObject {
public:
  TCalSM() : fSM(0), fAvgT(0), fT1(0), fT2(0), fT3(0), fT4(0), fT5(0), fT6(0), fT7(0), fFracLed(0), fFracMon(0) {;}
  virtual ~TCalSM() {;}
  Short_t    fSM;       // super module index
  Double32_t fAvgT;     //[0,0,16] sm T
  Double32_t fT1;       //[0,0,16] sensor T 1
  Double32_t fT2;       //[0,0,16] sensor T 2
  Double32_t fT3;       //[0,0,16] sensor T 3
  Double32_t fT4;       //[0,0,16] sensor T 4
  Double32_t fT5;       //[0,0,16] sensor T 5
  Double32_t fT6;       //[0,0,16] sensor T 6
  Double32_t fT7;       //[0,0,16] sensor T 7
  Double32_t fT8;       //[0,0,16] sensor T 8
  Double32_t fFracLed;  // fraction led info for each SM
  Double32_t fFracMon;  // fraction mon info for each SM
  ClassDef(TCalSM, 1); // CalSM class
};

class TCalInfo : public TObject {
 public:
   TCalInfo() : fRunNo(0), fAvTime(0), fFirstTime(0), fLastTime(0), fMinT(0), fMaxT(0), fFracS(0), fSMs("TCalSM"), fCells("TCalCell") {}
  virtual ~TCalInfo() {;}
  Int_t        fRunNo;     // run number
  UInt_t       fAvTime;    // average start time
  UInt_t       fFirstTime; // first time
  UInt_t       fLastTime;  // last time
  Float_t      fMinT;      // min temperature
  Float_t      fMaxT;      // max temperature
  Float_t      fFracS;     // fraction good sensors
  TClonesArray fSMs;       // array with SM infos
  TClonesArray fCells;     // array with cells
  ClassDef(TCalInfo, 2); // CalInfo class
};
#endif

void createTree(const char *period,
                const char *ofile     = "treeout.root",
                Bool_t doprint        = 0,
                Bool_t appBC          = kFALSE,                   // boolean to switch on bad channel
                TString badpath       = "$ALICE_ROOT/OADB/EMCAL", // location of bad channel map
                Int_t referenceRun    = 286313                    // define reference run to which all are being calibrated
) {
  TDraw td(period);
  td.Compute();
  LDraw ld(period);
  ld.Compute();

  if (appBC) {
    cout << "INFO: will be using bad channel map from: " << badpath.Data() << endl;
  }

  TObjArray *ta = td.GetArray();
  if (!ta) {
    cerr << "No time objects for period " << period << endl;
    return;
  }
  TObjArray *la = ld.GetArray();
  if (!la) {
    cerr << "No led objects for period " << period << endl;
    return;
  }

  const Int_t rns=la->GetEntries();
  cout << "Working on period " << period << " with " << rns << " runs" << endl;

  AliEMCALGeometry *g= AliEMCALGeometry::GetInstanceFromRunNumber(referenceRun);
  const Int_t kSM=g->GetNumberOfSuperModules();
  const Int_t kNcells=g->GetNCells();
  const Int_t gain = 1;

  TCalInfo *info = new TCalInfo;
  TFile *out = TFile::Open(ofile,"recreate");
  out->SetCompressionLevel(9);
  TTree* fTree = new TTree("tcal", "Temp calibration tree");
  fTree->SetDirectory(out);
  fTree->Branch("event", &info, 32000, 99);

  TClonesArray &carr    = info->fCells;
  TClonesArray &cSMarr  = info->fSMs;
  Int_t l = 0;
  Int_t t = 0;

  Int_t runno = -1;
  Int_t idx   = -1;
  TObjArray *arrayBC=0;
  TH2I *badmaps[kSM];
  for (Int_t i=0;i<kSM;++i)
    badmaps[i]=0;

  for (Int_t i=0;i<rns;++i) {
    l++;
    t++;
    LInfo *linfo = dynamic_cast<LInfo*>(la->At(l));
    if (!linfo){
      cout << "skipping due to missing info in LED tree!" << endl;
      t++;
      continue;
    }
    TInfo *tinfo = dynamic_cast<TInfo*>(ta->At(t));
    if (!tinfo){
      cout << "skipping due to missing info in temp tree!" << endl;
      l++;
      continue;
    }
    Int_t runl = linfo->GetRunNo();
    Int_t runt = tinfo->GetRunNo();
    if (runl!=runt) {
      if (runl > runt){
        l--;
      } else {
        t--;
      }
      cout << " Run numbers differ, skipping " << runl << " " << runt << endl;
      continue;
    }
    cout << "Working on run " << runl << endl;
    gSystem->Sleep(0.5);
    info->fRunNo     = runt;
    info->fAvTime    = tinfo->GetAverageTime();
    info->fFirstTime = tinfo->GetFirstTime();
    info->fLastTime  = tinfo->GetLastTime();
    info->fMinT      = tinfo->AbsMinT();
    info->fMaxT      = tinfo->AbsMaxT();
    info->fFracS     = tinfo->Fraction();

    cSMarr.Clear();
    cSMarr.ExpandCreate(kSM);
    for (Int_t sm=0; sm<kSM; sm++) {
      cout << "SM  "<< sm << ":\t"<< tinfo->AvgTempSM(sm) << endl;
      TCalSM *smInfo    = (TCalSM*)cSMarr.At(sm);
      smInfo->fSM       = sm;
      smInfo->fAvgT     = tinfo->AvgTempSM(sm);
      smInfo->fT1       = tinfo->T(sm*8,3);
      smInfo->fT2       = tinfo->T(sm*8+1,3);
      smInfo->fT3       = tinfo->T(sm*8+2,3);
      smInfo->fT4       = tinfo->T(sm*8+3,3);
      smInfo->fT5       = tinfo->T(sm*8+4,3);
      smInfo->fT6       = tinfo->T(sm*8+5,3);
      smInfo->fT7       = tinfo->T(sm*8+6,3);
      smInfo->fT8       = tinfo->T(sm*8+7,3);
      smInfo->fFracLed  = linfo->FracLeds(sm,gain);
      smInfo->fFracMon  = linfo->FracStrips(sm,gain);
    }

    carr.Clear();
    carr.ExpandCreate(kNcells);

    Double_t avg[20];
    for (Int_t sm=0; sm<kSM; ++sm) {
      avg[sm]=tinfo->AvgTempSM(sm);
    }

    if (runt!=runno && appBC) {
      runno=runt;
      AliOADBContainer *contBC = new AliOADBContainer("");
      contBC->InitFromFile(Form("%s/EMCALBadChannels.root",badpath.Data()),"AliEMCALBadChannels");
      if (!contBC) {
        cerr << "Could not load bc map for run " << runno << endl;
      } else {
        Int_t idxCurr     =contBC->GetIndexForRun(runno);
        if ( idx != idxCurr ){
          cout << "INFO: need to switch bad channel map" << endl;
          TObjArray *arrayBC=(TObjArray*)contBC->GetObject(runno);
          if (!arrayBC){
            cout << "WARNING: missing bad channel map for run: " << runno << " Will continue without using bad channel map!"<< endl;
          }
          for (Int_t i=0; i<kSM; ++i) {
            delete badmaps[i];
            if (!arrayBC){
              badmaps[i] = NULL;
            } else {
              badmaps[i] = (TH2I*)arrayBC->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
              badmaps[i]->SetDirectory(0);
            }
          }
          idx = idxCurr;
        }
        delete contBC;
      }
    }

    for (Int_t sm=0; sm<kSM; ++sm) {
      Int_t nrow = g->GetNumberOfCellsInPhiDirection(sm);
      Int_t ncol = g->GetNumberOfCellsInEtaDirection(sm);

      TH1 *hmonm = linfo->GetLedMonHist(sm,gain);
      TH1 *hmonr = linfo->GetLedMonRmsHist(sm,gain);
      TH2 *hledm = linfo->GetLedHist(sm,gain);
      TH2 *hledr = linfo->GetLedRmsHist(sm,gain);

      for (Int_t col=0; col<ncol; ++col) {
        for (Int_t row=0; row<nrow; ++row) {
          Int_t  id = g->GetAbsCellIdFromCellIndexes(sm,row,col);
          TCalCell *cell = (TCalCell*)carr.At(id);
          cell->fId = id;
          Int_t badcell = 0;
          if (appBC && badmaps[sm])
            badcell = badmaps[sm]->GetBinContent(col,row);
          cell->fSM = sm;
          cell->fRow = row;
          cell->fCol = col;
          cell->fBad = badcell;

          Int_t orow=row;
          Int_t ocol=col;
          // shift to online row(phi)/col(eta)
          g->ShiftOfflineToOnlineCellIndexes(sm, orow, ocol);
          Int_t ns = TInfo::SensId(sm,orow,ocol);
          cell->fLedM = hledm->GetBinContent(hledm->FindBin(ocol,orow));
          cell->fLedR = hledm->GetBinError(hledm->FindBin(ocol,orow));
          cell->fMonM = hmonm->GetBinContent(hmonm->FindBin(ocol/2));
          cell->fMonR = hmonm->GetBinError(hmonm->FindBin(ocol/2));
          cell->fLocT = tinfo->T(ns,3);
          cell->fSMT  = avg[sm];
        }
      }
    }
    fTree->Fill();
  }
  out->cd();
  fTree->Write();
  out->Close();
}

void createTree_all()
{
  createTree("lhc16f","tree_lhc16f.root");
  createTree("lhc16g","tree_lhc16g.root");
  createTree("lhc16h","tree_lhc16h.root");
  createTree("lhc16i","tree_lhc16i.root");
  createTree("lhc16j","tree_lhc16j.root");
  createTree("lhc16k","tree_lhc16k.root");
  createTree("lhc16l","tree_lhc16g.root");
  createTree("lhc16o","tree_lhc16g.root");
  createTree("lhc16p","tree_lhc16g.root");
  createTree("lhc17p","tree_lhc17p.root");
  createTree("lhc18d","tree_lhc18d.root");
}
