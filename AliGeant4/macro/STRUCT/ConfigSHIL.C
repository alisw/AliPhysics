void Config(Int_t version)
{
  AliSHIL* SHIL = 0;
  switch (version) {
    case 0: SHIL  = new AliSHILv0("SHIL","Shielding"); break;
    case 1: SHIL  = new AliSHILvF("SHIL","Shielding"); break;
    case 2: SHIL  = new AliSHILv2("SHIL","Shielding"); break;
  }  

//=================== SHIL parameters ============================
}  
