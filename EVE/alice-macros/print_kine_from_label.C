void print_kine_from_label(Int_t label)
{
  AliRunLoader* rl = Alieve::Event::AssertRunLoader();
  rl->LoadKinematics();
  AliStack* stack = rl->Stack();

  printf("Number primaries %d, all particles %d, label %d\n",
	 stack->GetNprimary(), stack->GetNtrack(), label);
  if (label < 0 || label >= stack->GetNtrack()) {
    printf("  Label exceeds available range.\n");
    return;
  }

  TParticle* part = stack->Particle(label);
  if(part != 0) {
    part->Print();
    while(part->GetMother(0) >= 0) {
      part = stack->Particle(part->GetMother(0));
      part->Print();
    }
  }
}
