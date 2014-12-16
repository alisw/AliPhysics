void ReadAlignmentData()
{
  TFile * f = TFile::Open("$ALICE_ROOT/OCDB/PMD/Align/Data/Run0_999999999_v0_s0.root");

  f->ls();

  AliAlignObjMatrix * aam;

  TGeoHMatrix hh;
  Double_t tr[3];
  if (!AliCDBEntry)
    {
      printf("Something is wrong ************ \n");
    }
  else if(AliCDBEntry)
    {
      AliCDBEntry->PrintId(); 
      AliCDBEntry->PrintMetaData();

      TClonesArray * ncut = 0;
      ncut = (TClonesArray*)AliCDBEntry->GetObject();
      ncut->Print(); 

      Int_t nen = ncut->GetLast();
      cout << nen << endl;

      for (int i=0; i<4; i++)
	{
	  aam = (AliAlignObjMatrix*)ncut->UncheckedAt(i);

	  aam->GetMatrix(hh);
	  //hh.Print();

	  aam->GetTranslation(tr);

	  cout << tr[0] << " " << tr[1] << " " << tr[2] << endl;

	}

    }
}
