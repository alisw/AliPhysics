void Config(Int_t version)
{
  AliDIPO* DIPO = 0;
  switch (version) {
    case 1: DIPO  = new AliDIPOv1("DIPO","DIPOv1 module"); break;
    case 2: DIPO  = new AliDIPOv2("DIPO","Dipole version 2"); break;
  }  

//=================== DIPO parameters ============================
}  
