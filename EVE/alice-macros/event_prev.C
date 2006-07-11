void event_prev()
{
  if(Alieve::gEvent == 0) {
    printf("Event not set!\n");
    return;
  }
  Alieve::gEvent->PrevEvent();
}
