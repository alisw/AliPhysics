void emcal_digits()
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();

  rl->LoadgAlice();
  AliEMCAL         * emcal = (AliEMCAL*) rl->GetAliRun()->GetDetector("EMCAL");
  AliEMCALGeometry * geom  = emcal->GetGeometry();

  rl->LoadDigits("EMCAL");
  TTree* dt = rl->GetTreeD("EMCAL", kFALSE);

  gGeoManager = gReve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");

  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
  bbbox->Dump();
  TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();
  sbbox->Dump();

  Reve::RenderElementList* l = new Reve::RenderElementList("EMCAL");
  l->SetTitle("Tooltip");
  gReve->AddRenderElement(l);

  Reve::FrameBox* frame_big = new Reve::FrameBox();
  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());

  Reve::FrameBox* frame_sml = new Reve::FrameBox();
  frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());

  gStyle->SetPalette(1, 0);
  Reve::RGBAPalette* pal = new Reve::RGBAPalette(0, 512);
  pal->SetLimits(0, 1024);

  Reve::QuadSet* smodules[12];

  for (Int_t sm=0; sm<12; ++sm)
  {
    Reve::QuadSet* q = new Reve::QuadSet(Form("SM %d", sm+1));
    q->SetOwnIds(kTRUE);
    q->Reset(Reve::QuadSet::QT_RectangleYZFixedDimX, kFALSE, 32);
    q->SetDefWidth (geom->GetPhiTileSize());
    q->SetDefHeight(geom->GetEtaTileSize());

    q->RefHMTrans().SetFrom(*node->GetDaughter(sm)->GetMatrix());

    q->SetFrame(sm < 10 ? frame_big : frame_sml);
    q->SetPalette(pal);

    gReve->AddRenderElement(q, l);
    smodules[sm] = q;
  }

  TClonesArray *digits = 0;
  dt->SetBranchAddress("EMCAL", &digits);
  dt->GetEntry(0);
  Int_t nEnt = digits->GetEntriesFast();
  AliEMCALDigit * dig;

  Int_t iEvent  = -1 ;
  Float_t amp   = -1 ;
  Float_t time  = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;
  Int_t ieta    =  0 ;
  Double_t x, y, z;

  for(Int_t idig = 0; idig<nEnt; idig++)
  {
    dig = static_cast<AliEMCALDigit *>(digits->At(idig));

    if(dig != 0) {
      id   = dig->GetId() ; //cell (digit) label
      amp  = dig->GetAmp(); //amplitude in cell (digit)
      time = dig->GetTime();//time of creation of digit after collision

      cout<<"Cell ID "<<id<<" Amp "<<amp<<endl;//" time "<<time<<endl;

      //Geometry methods  
      geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta); 
      //Gives SuperModule and Tower numbers
      geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					iIphi, iIeta,iphi,ieta);
      //Gives label of cell in eta-phi position per each supermodule

      cout<< "SModule "<<iSupMod<<"; Tower "<<iTower
	  <<"; Eta "<<iIeta<<"; Phi "<<iIphi
	  <<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<endl;

      geom->RelPosCellInSModule(id, x, y, z);
      cout << x <<" "<< y <<" "<< z <<endl;

      Reve::QuadSet* q = smodules[iSupMod];
      q->AddQuad(y, z);
      q->QuadValue(amp);
      q->QuadId(dig);
    } else {
      cout<<"Digit pointer 0x0"<<endl;
    }
  }

  for (Int_t sm=0; sm<12; ++sm)
  {
    smodules[iSupMod]->RefitPlex();
  }

  gReve->Redraw3D();
}
