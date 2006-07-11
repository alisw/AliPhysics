void event_next()
{
  if(Alieve::gEvent == 0) {
    printf("Event not set!\n");
    return;
  }
  Alieve::gEvent->NextEvent();
}
