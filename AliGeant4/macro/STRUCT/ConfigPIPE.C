void Config(Int_t version)
{
  AliPIPE* PIPE = 0;
  switch (version) {
    case 0: PIPE  = new AliPIPEv0("PIPE","Beam Pipe"); break;
    case 1: PIPE  = new AliPIPEv1("PIPE","PIPEv1 module"); break;
    case 3: PIPE  = new AliPIPEv3("PIPE","PIPEv3 module"); break;
  }

//=================== PIPE parameters ============================
}  
