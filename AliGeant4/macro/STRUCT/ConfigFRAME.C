void Config(Int_t version)
{
  AliFRAME* FRAME = 0;
  switch (version) {
    case 0: FRAME = new AliFRAMEv0("FRAME","Space Frame");   break;
    case 1: FRAME = new AliFRAMEv1("FRAME","FRAMEv1 module"); break;
  }  

//=================== FRAME parameters ============================
}  
