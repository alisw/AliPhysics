#include <iomanip.h>

void MUON_track_info(Int_t label) {

  Double_t RADDEG = 180.0/TMath::Pi();

  AliRunLoader* rl = Alieve::Event::AssertRunLoader();

  TTree* tt = rl->GetTreeT("MUON", false);
  
  TClonesArray *tracks = 0;
  tt->SetBranchAddress("MUONTrack",&tracks);
  tt->GetEntry(0);

  Int_t ntracks = tracks->GetEntriesFast();
  
  if (ntracks <= label) return;

  AliMUONTrack *mt = tracks->At(label);  

  cout << endl;
  cout << ">>>>>#########################################################################################################################" << endl;
  cout << endl;
  cout << "                   Track number " << label << endl;
  cout << endl;
  cout << "   Number of clusters       " << mt->GetNTrackHits() << endl;
  cout << "   Match to trigger         " << mt->GetMatchTrigger() << endl;
  if (mt->GetMatchTrigger()) {
    cout << "   Chi2 tracking-trigger    " << mt->GetChi2MatchTrigger() << endl;
    cout << "   Local trigger number     " << mt->GetLoTrgNum() << endl;
  }
 
  TClonesArray *trackParamAtHit = mt->GetTrackParamAtHit();
  Int_t nparam = trackParamAtHit->GetEntries();
  AliMUONTrackParam *mtp;

  cout << endl;
  cout << "   TrackParamAtHit entries  " << nparam << "" << endl;

  cout << "   (slopes [deg], coord [cm], p [GeV/c])" << endl;
  cout << "   -----------------------------------------------------" << endl;
  cout << "   Number   InvBendMom   BendSlope   NonBendSlope   BendCoord   NonBendCoord      Z        Px        Py        Pz         P" << endl;

  for (Int_t i = 0; i < nparam; i++) {

    mtp = (AliMUONTrackParam*)trackParamAtHit->At(i);

    cout << 
      setw(9)<< setprecision(3) << 
      i << "   " << 

      setw(8) << setprecision(3) << 
      mtp->GetInverseBendingMomentum() << "    " << 

      setw(8) << setprecision(3) <<
      mtp->GetBendingSlope()*RADDEG << "    " << 

      setw(8) << setprecision(3) <<
      mtp->GetNonBendingSlope()*RADDEG << "      " << 

      setw(8) << setprecision(3) <<
      mtp->GetBendingCoor() << "    " <<  

      setw(8) << setprecision(3) <<
      mtp->GetNonBendingCoor() << "      " <<  

      setw(10) << setprecision(6) <<
      mtp->GetZ() << "  " <<  

      setw(8) << setprecision(3) <<
      mtp->Px() << "  " <<  

      setw(8) << setprecision(3) <<
      mtp->Py() << "  " <<  

      setw(8) << setprecision(3) <<
      mtp->Pz() << "  " <<  

      setw(8) << setprecision(3) <<
      mtp->P() << "  " <<  

      endl;

  }

  cout << endl;
  
  cout << "   Track parameters at vertex" << endl;
  cout << "   -----------------------------------------------------" << endl;
  cout << "   InvBendMom   BendSlope   NonBendSlope   BendCoord   NonBendCoord         Z          Px        Py        Pz         P" << endl;

  mtp = (AliMUONTrackParam*)mt->GetTrackParamAtVertex();

  cout << "  " <<
    setw(8) << setprecision(3) << 
    mtp->GetInverseBendingMomentum() << "     " << 
    
    setw(8) << setprecision(3) <<
    mtp->GetBendingSlope()*RADDEG << "    " << 
    
    setw(8) << setprecision(3) <<
    mtp->GetNonBendingSlope()*RADDEG << "         " << 
    
    setw(8) << setprecision(3) <<
    mtp->GetBendingCoor() << "    " <<  
    
    setw(8) << setprecision(3) <<
    mtp->GetNonBendingCoor() << "      " <<  
    
    setw(10) << setprecision(3) <<
    mtp->GetZ() << "  " <<  
    
    setw(8) << setprecision(3) <<
    mtp->Px() << "  " <<  
    
    setw(8) << setprecision(3) <<
    mtp->Py() << "  " <<  
    
    setw(8) << setprecision(3) <<
    mtp->Pz() << "  " <<  
    
    setw(8) << setprecision(3) <<
    mtp->P() << "  " <<  
    
    endl;
  
  Float_t pt = TMath::Sqrt(mtp->Px()*mtp->Px()+mtp->Py()*mtp->Py());

  cout << "" << endl;
  cout << "   Pt = " << 
    setw(8) << setprecision(3) <<
    pt << "  GeV/c" << endl;
  
  cout << endl;
  cout << "#########################################################################################################################<<<<<" << endl;
  
}
