MakeCDBEntryProblematic(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), const Char_t *filename = NULL)
{

  /* create object */
  TH1C *obj = new TH1C("hProblematic", "", 157248, 0., 157248.);

  /* update object */
  if (filename) {
    printf("PROBLEMATIC HISTO WILL BE UPDATED ACCORDING TO INPUT LIST\n");
    printf("inputList: %s\n", filename);
    UpdateProblematicHisto(obj, filename);
  }
  else {
    printf("EMPTY PROBLEMATIC HISTO WILL BE GENERATED\n");
  }

  /* create cdb info */
  AliCDBId id("TOF/Calib/Problematic", startRun, endRun);
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("Roberto Preghenella");
  md->SetComment("Problematic");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);

  /* put object in cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->GetDefaultStorage()->Put(obj, id, md);

}

Bool_t
UpdateProblematicHisto(TH1C *histo, const Char_t *filename)
{
  /*
   * this routine updates the problematic histo taking
   * problematic channels from a list specified in the
   * file written in ASCII format according to the
   * following formats based on electronics-oriented indices:
   *
   * # disable single crate
   * # crate[0-71]
   *
   * # disable single TRM
   * # crate[0-71] trm[3-12]
   *
   * # disable single chain
   * # crate[0-71] trm[3-12] chain[0-1]
   *
   * # disable single TDC
   * # crate[0-71] trm[3-12] chain[0-1] tdc[0-14]
   *
   * # disable single channel
   * # crate[0-71] trm[3-12] chain[0-1] tdc[0-14] channel[0-7]
   *
   *
   * EXAMPLE:
   * # this list will set as problematics:
   * # - all channels of crate-39
   * # - all channels of TRM-03 crate-40
   * # - all channels of chain-A TRM-07 crate-04
   * # - all channels of TDC-04 chain-B TRM-09 crate-07
   * # - channel-03 TDC-02 chain-A TRM-04 crate-00
   * 39
   * 40 3
   * 4 7 0
   * 7 9 1 4
   * 0 4 0 2 3
   *
   */

  /* check histo */
  if (!histo) {
    printf("WARNING: NULL histo, will just run in DUMMY mode\n");
  }

  /* open ASCII file */
  ifstream filein;
  filein.open(filename, ifstream::in);
  
  /* loop over lines in file */
  char buf[1024];
  TString str;
  TObjString *ostr;
  TObjArray *oa;
  Int_t crate, trm, chain, tdc, channel;
  Int_t irequest = 0, nflagged;
  printf("processing requests to flag problematic channels:\n");
  while (filein.good()) {
    filein.getline(buf, 1024);
    /* check whether we got to EOF */
    if (filein.eof()) break;
    /* put buffer in a TString */
    str = buf;
    /* check whether commented line */
    if (str.BeginsWith("#")) continue;
    irequest++;
    /* tokenize */
    oa = str.Tokenize(" ");
    switch (oa->GetEntries()) {
    case 1:
      ostr = (TObjString *)oa->At(0);
      crate = atoi(ostr->GetName());
      if (crate < 0 || crate > 71) {
	printf("%d.) invalid crate number: %d\n", irequest, crate);
	break;
      }
      nflagged = FlagAsProblematic(histo, crate);
      printf("%d.) crate flagged as problematic (%d channels): crate-%02d\n", irequest, nflagged, crate);
      break;
    case 2:
      ostr = (TObjString *)oa->At(0);
      crate = atoi(ostr->GetName());
      if (crate < 0 || crate > 71) {
	printf("%d.) invalid crate number: %d\n", irequest, crate);
	break;
      }
      ostr = (TObjString *)oa->At(1);
      trm = atoi(ostr->GetName());
      if (trm < 3 || trm > 12) {
	printf("%d.) invalid TRM number: %d\n", irequest, trm);
	break;
      }
      nflagged = FlagAsProblematic(histo, crate, trm);
      printf("%d.) TRM flagged as problematic (%d channels): crate-%02d TRM-%02d\n", irequest, nflagged, crate, trm);
      break;
    case 3:
      ostr = (TObjString *)oa->At(0);
      crate = atoi(ostr->GetName());
      if (crate < 0 || crate > 71) {
	printf("%d.) invalid crate number: %d\n", irequest, crate);
	break;
      }
      ostr = (TObjString *)oa->At(1);
      trm = atoi(ostr->GetName());
      if (trm < 3 || trm > 12) {
	printf("%d.) invalid TRM number: %d\n", irequest, trm);
	break;
      }
      ostr = (TObjString *)oa->At(2);
      chain = atoi(ostr->GetName());
      if (chain < 0 || chain > 1) {
	printf("%d.) invalid chain number: %d\n", irequest, chain);
	break;
      }
      nflagged = FlagAsProblematic(histo, crate, trm, chain);
      printf("%d.) chain flagged as problematic (%d channels): crate-%02d TRM-%02d chain-%s\n", irequest, nflagged, crate, trm, chain == 0 ? "A" : "B");
      break;
    case 4:
      ostr = (TObjString *)oa->At(0);
      crate = atoi(ostr->GetName());
      if (crate < 0 || crate > 71) {
	printf("%d.) invalid crate number: %d\n", irequest, crate);
	break;
      }
      ostr = (TObjString *)oa->At(1);
      trm = atoi(ostr->GetName());
      if (trm < 3 || trm > 12) {
	printf("%d.) invalid TRM number: %d\n", irequest, trm);
	break;
      }
      ostr = (TObjString *)oa->At(2);
      chain = atoi(ostr->GetName());
      if (chain < 0 || chain > 1) {
	printf("%d.) invalid chain number: %d\n", irequest, chain);
	break;
      }
      ostr = (TObjString *)oa->At(3);
      tdc = atoi(ostr->GetName());
      if (tdc < 0 || tdc > 14) {
	printf("%d.) invalid chain number: %d\n", irequest, chain);
	break;
      }
      nflagged = FlagAsProblematic(histo, crate, trm, chain, tdc);
      printf("%d.) TDC flagged as problematic (%d channels): crate-%02d TRM-%02d chain-%s TDC-%02d\n", irequest, nflagged, crate, trm, chain == 0 ? "A" : "B", tdc);
      break;
    case 5:
      ostr = (TObjString *)oa->At(0);
      crate = atoi(ostr->GetName());
      if (crate < 0 || crate > 71) {
	printf("%d.) invalid crate number: %d\n", irequest, crate);
	break;
      }
      ostr = (TObjString *)oa->At(1);
      trm = atoi(ostr->GetName());
      if (trm < 3 || trm > 12) {
	printf("invalid TRM number: %d\n", irequest, trm);
	break;
      }
      ostr = (TObjString *)oa->At(2);
      chain = atoi(ostr->GetName());
      if (chain < 0 || chain > 1) {
	printf("%d.) invalid chain number: %d\n", irequest, chain);
	break;
      }
      ostr = (TObjString *)oa->At(3);
      tdc = atoi(ostr->GetName());
      if (tdc < 0 || tdc > 14) {
	printf("%d.) invalid chain number: %d\n", irequest, chain);
	break;
      }
      ostr = (TObjString *)oa->At(4);
      channel = atoi(ostr->GetName());
      if (channel < 0 || channel > 7) {
	printf("%d.) invalid channel number: %d\n", irequest, channel);
	break;
      }
      nflagged = FlagAsProblematic(histo, crate, trm, chain, tdc, channel);
      printf("%d.) channel flagged as problematic (%d channels): crate-%02d TRM-%02d chain-%s TDC-%02d, channel-%d\n", irequest, nflagged, crate, trm, chain == 0 ? "A" : "B", tdc, channel);
      break;
    default:
      printf("%d.) invalid format: %s\n", irequest, str.Data());
      break;
    }
  }

  /* close file */
  filein.close();

  return kTRUE;
}

Int_t
FlagAsProblematic(TH1C *histo, Int_t crate = -1, Int_t trm = -1, Int_t chain = -1, Int_t tdc = -1, Int_t channel = -1)
{

  /*
   * flag as problematic according to parameters
   */

  /* loop over everything checking request */
  Int_t det[5], dummy, index, nflagged = 0;
  for (Int_t icrate = 0; icrate < 72; icrate++) {
    if (crate != -1 && icrate != crate) continue;
    for (Int_t itrm = 3; itrm <= 12; itrm++) {
      if (trm != -1 && itrm != trm) continue;
      for (Int_t ichain = 0; ichain < 2; ichain++) {
	if (chain != -1 && ichain != chain) continue;
	for (Int_t itdc = 0; itdc < 15; itdc++) {
	  if (tdc != -1 && itdc != tdc) continue;
	  for (Int_t ichannel = 0; ichannel < 8; ichannel++) {
	    if (channel != -1 && ichannel != channel) continue;
	    AliTOFRawStream::EquipmentId2VolumeId(icrate, itrm, ichain, itdc, ichannel, det);
	    dummy = det[4];
	    det[4] = det[3];
	    det[3] = dummy;
	    if (det[0] < 0 || det[0] > 17 ||
		det[1] < 0 || det[1] > 4 ||
		det[2] < 0 || det[2] > 18 ||
		det[3] < 0 || det[3] > 1 ||
		det[4] < 0 || det[4] > 47) {
	      //	      printf("invalid volume indices: EO = (%d %d %d %d %d), VOL = (%d %d %d %d %d)\n", icrate, itrm, ichain, itdc, ichannel, det[0], det[1], det[2], det[3], det[4]);
	      continue;
	    }
	    index = AliTOFGeometry::GetIndex(det);
	    if (index < 0 || index > 157248) {
	      //	      printf("invalid calib index: EO = (%d %d %d %d %d), VOL = (%d %d %d %d %d), CAL = %d\n", icrate, itrm, ichain, itdc, ichannel, det[0], det[1], det[2], det[3], det[4], index);
	      continue;
	    }
	    nflagged++;
	    if (!histo) continue;
	    histo->SetBinContent(index + 1, 0x1);
	  }
	}
      }
    }
  }

  return nflagged;
}
