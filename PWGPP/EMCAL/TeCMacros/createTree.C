#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <TCanvas.h>
#include <TClonesArray.h>
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
  TCalCell() : fLedM(0), fLedR(0), fMonM(0), fMonR(0), fLocT(0), fSMT(0) {;}
  virtual ~TCalCell() {;}
  Short_t  fId;   //         cell id
  Short_t  fSM;   //         super module index
  Short_t  fRow;  //         row (phi) index
  Short_t  fCol;  //         col (eta) index
  Double32_t fLedM; //[0,0,16] led mean
  Double32_t fLedR; //[0,0,16] led rms
  Double32_t fMonM; //[0,0,16] mon mean
  Double32_t fMonR; //[0,0,16] mon rms
  Double32_t fLocT; //[0,0,16] loc T
  Double32_t fSMT;  //[0,0,16] sm T
  ClassDef(TCalCell, 2); // CalCell class
};

class TCalInfo : public TObject {
 public:
  TCalInfo() : fRunNo(0), fAvTime(0), fFirstTime(0), fLastTime(0), fMinT(0), fMaxT(0), fFracS(0), fAvgTemp(160), fFracLed(20), fFracMon(20), fCells("TCalCell") {}
  virtual ~TCalInfo() {;}
  Int_t        fRunNo;     // run number
  UInt_t       fAvTime;    // average start time
  UInt_t       fFirstTime; // first time
  UInt_t       fLastTime;  // last time
  Float_t      fMinT;      // min temperature
  Float_t      fMaxT;      // max temperature
  Float_t      fFracS;     // fraction good sensors
  TArrayF      fAvgTemp;   // average temp for each sensor
  TArrayF      fFracLed;   // fraction led info for each SM
  TArrayF      fFracMon;   // fraction mon info for each SM
  TClonesArray fCells;     // array with cells
  ClassDef(TCalInfo, 1); // CalInfo class
};
#endif

void createTree(const char *period, const char *ofile="treeout.root",Bool_t doprint=0)
{
  TDraw td(period);
  td.Compute();
  LDraw ld(period);
  ld.Compute();

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

  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  const Int_t kSM=g->GetNumberOfSuperModules();
  const Int_t kNcells=g->GetNCells();
  const Int_t gain = 1;

  TCalInfo *info = new TCalInfo;
  TFile *out = TFile::Open(ofile,"recreate");
  out->SetCompressionLevel(9);
  TTree* fTree = new TTree("tcal", "Temp calibration tree");
  fTree->SetDirectory(out);
  fTree->Branch("event", &info, 32000, 99);

  TClonesArray &carr = info->fCells;
  Int_t l = 0;
  Int_t t = 0;
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
    for (Int_t sm=0; sm<kSM; ++sm) {
      info->fFracLed.SetAt(sm,linfo->FracLeds(sm,gain));
      info->fFracMon.SetAt(sm,linfo->FracStrips(sm,gain));
    }
    for (Int_t n=0; n<160; ++n) {
      info->fAvgTemp.SetAt(n,tinfo->T(n,3));
    }

    Int_t ncells = 0;
    carr.Clear();
    carr.ExpandCreate(kNcells);

    Double_t avg[20];
    for (Int_t sm=0; sm<kSM; ++sm) {
      avg[sm]=tinfo->AvgT(sm);
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
          Int_t ns = TInfo::SensId(sm,row,col);
          TCalCell *cell = (TCalCell*)carr.At(id);
          cell->fId = id;
          cell->fSM = sm;
          cell->fRow = row;
          cell->fCol = col;
          cell->fLedM = hledm->GetBinContent(hledm->FindBin(col,row));
          cell->fLedR = hledr->GetBinContent(hledr->FindBin(col,row));
          cell->fMonM = hmonm->GetBinContent(hmonm->FindBin(col/2));
          cell->fMonR = hmonr->GetBinContent(hmonr->FindBin(col/2));
          cell->fLocT = tinfo->T(ns,3);
          cell->fSMT  = avg[sm];
        }
      }
    }
    fTree->Fill();
  }
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
