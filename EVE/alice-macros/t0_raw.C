// $Id$

namespace Alieve {
class Event;
}



void t0_raw()
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();
  Int_t ievt = Alieve::gEvent->GetEventId();
    cout<<ievt<<endl;

  gStyle->SetPalette(1, 0);

  Alieve::T0Module::LoadRaw("raw.root",ievt);


}

