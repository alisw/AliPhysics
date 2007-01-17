void event_goto(Int_t event=0)
{
  if(Alieve::gEvent == 0) {
    printf("Event not set!\n");
    return;
  }
  Alieve::gEvent->GotoEvent(event);
}
