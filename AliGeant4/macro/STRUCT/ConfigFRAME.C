void Config(Int_t version)
{
  AliFRAME* FRAME = 0;
  switch (version) {
    case 0: FRAME = new AliFRAMEv0("FRAME","FRAMEv0 module");   break;
    case 1: FRAME = new AliFRAMEv1("FRAME","FRAMEv1 module"); break;
    case 2: FRAME = new AliFRAMEv2("FRAME","Space frame"); 
            if (AliRunConfiguration::Holes())
               ((AliFRAMEv2*) FRAME)->SetHoles(1);
	    else 
               ((AliFRAMEv2*) FRAME)->SetHoles(0);
	    break;    
  }  

//=================== FRAME parameters ============================

}  
