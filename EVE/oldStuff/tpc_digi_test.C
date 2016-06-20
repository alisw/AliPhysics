// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TTree;

class AliEveTPCData;

AliEveTPCData*  x = 0;
TTree*            tree = 0;

void tpc_digi_test()
{
  gROOT->Macro("alieve_loadlibs.C");
  gSystem->Load("libAlieve");

  TFile* f = new TFile("coctail_1k/TPC.Digits.root");
  tree = (TTree*) gDirectory->Get("Event0/TreeD");

  x = new AliEveTPCData;
  // x->SetSectorBlockSize(8192);
  // x->SetLoadThreshold(5);
  x->CreateAllSectors();
  x->LoadDigits(tree, kFALSE);
  gStyle->SetPalette(1, 0);

  AliEveTPCSector2D* s = new AliEveTPCSector2D();
  s->SetDataSource(x);
  s->SetMainColor(36);
  gEve->AddElement(s);
  gEve->DrawElement(s);
}


void tpc_digi_pad_dump(Int_t s, Int_t r, Int_t p)
{
  if(r >= AliEveTPCSectorData::GetInnSeg().fNRows) {
    r -= AliEveTPCSectorData::GetInnSeg().fNRows;
    s += 36;
  }

  AliSimDigits *digit = 0;
  tree->GetBranch("Segment")->SetAddress(&digit);

  Int_t sbr = (Int_t) tree->GetEntries();
  for (Int_t ent=0; ent<sbr; ent++) {
    tree->GetEntry(ent);
    Int_t sector, row;
    AliEveTPCSectorData::GetParam().AdjustSectorRow(digit->GetID(), sector, row);

    if(sector != s || row != r)
      continue;

    printf("Entry = %d, ID = %d, Sec = %d, Row = %d\n",
	   ent, digit->GetID(), sector, row);

    Int_t pad;
    Short_t time, signal;
    digit->First();
    do {
      pad    = digit->CurrentColumn();
      time   = digit->CurrentRow();
      signal = digit->CurrentDigit();

      if(p == pad)
	printf("%d %d\n", time, signal);

    } while (digit->Next());
  }
}
