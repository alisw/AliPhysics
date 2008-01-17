// $Id$

namespace Alieve {
class Event;
}

void t0_digits()
{
  AliRunLoader* rl =  Alieve::Event::AssertRunLoader();

  rl->LoadDigits("T0");
  TTree* dt = rl->GetTreeD("T0", false);

  AliT0digit *digits = 0;
  dt->SetBranchAddress("T0", &digits);
  dt->GetEntry(0);

  gStyle->SetPalette(1, 0);

  Alieve::T0Module::MakeModules(digits);
}

