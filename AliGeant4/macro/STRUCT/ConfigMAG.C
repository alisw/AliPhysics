void Config(Int_t version)
{
//=================== All modules in STRUCT parameters =========================

  AliMAG* MAG = 0;
  switch (version) {
    case 0: MAG  = new AliMAG("MAG","Magnet"); break;
  }  

//=================== MAG parameters ============================
// --- Start with Magnet since detector layouts may be depending ---
// --- on the selected Magnet dimensions ---
}  
