void Config(Int_t version)
{
  AliABSO* ABSO = 0;
  switch (version) {
    case 0: ABSO  = new AliABSOv0("ABSO","Muon Absorber"); break;
  }

//=================== ABSO parameters ============================
}  
