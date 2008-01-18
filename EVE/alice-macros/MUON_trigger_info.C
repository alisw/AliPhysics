// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <iomanip.h>

void MUON_trigger_info(Int_t label) {

  AliRunLoader* rl = AliEveEventManager::AssertRunLoader();

  TTree* tt = rl->GetTreeT("MUON", false);

  TClonesArray *tracks = 0;
  tt->SetBranchAddress("AliEveMUONTrack",&tracks);
  tt->GetEntry(0);

  Int_t ntracks = tracks->GetEntriesFast();

  if (ntracks <= label) return;

  AliMUONTrack *mt = tracks->At(label);

  TTree* tr = rl->GetTreeR("MUON", false);

  TClonesArray *lotrig = 0;
  tr->SetBranchAddress("MUONLocalTrigger",&lotrig);
  tr->GetEntry(0);

  Int_t nlotrig = lotrig->GetEntriesFast();

  cout << endl;
  cout << ">>>>>#########################################################################################################################" << endl;
  cout << endl;
  cout << "                   TEveTrack number " << label << endl;
  cout << endl;
  AliMUONTrackParam *mtp = (AliMUONTrackParam*)mt->GetTrackParamAtVertex();
  Float_t pt = TMath::Sqrt(mtp->Px()*mtp->Px()+mtp->Py()*mtp->Py());

  cout << "   Pt = " <<
    setw(8) << setprecision(3) <<
    pt << "  GeV/c" << endl;

  cout << endl;

  if (mt->GetLoTrgNum() >= 0) {
    AliMUONLocalTrigger *lo = (AliMUONLocalTrigger*)lotrig->At(mt->GetLoTrgNum());
    cout << "   Local trigger information" << endl;
    cout << "   -----------------------------------------------------" << endl;
    cout << "   Circuit   " << lo->LoCircuit() << endl;
    cout << "   StripX    " << lo->LoStripX()  << endl;
    cout << "   StripY    " << lo->LoStripY()  << endl;
    cout << "   Dev       " << lo->LoDev()     << endl;
    cout << "   LoLpt     " << lo->LoLpt()     << endl;
    cout << "   LoHpt     " << lo->LoHpt()     << endl;
    cout << "   Pattern X:" << endl;
    printf("   %016b \n",lo->GetX1Pattern());
    printf("   %016b \n",lo->GetX2Pattern());
    printf("   %016b \n",lo->GetX3Pattern());
    printf("   %016b \n",lo->GetX4Pattern());
    cout << "   Pattern Y:" << endl;
    printf("   %016b \n",lo->GetY1Pattern());
    printf("   %016b \n",lo->GetY2Pattern());
    printf("   %016b \n",lo->GetY3Pattern());
    printf("   %016b \n",lo->GetY4Pattern());
    cout << "   Decision:" << endl;
    printf("   %04b  \n",lo->GetLoDecision());
  } else {
    cout << "   The track has no trigger information!" << endl;
  }
  cout << endl;
  cout << "#########################################################################################################################<<<<<" << endl;

}

