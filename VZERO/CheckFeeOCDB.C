void CheckFeeOCDB(Int_t run)
{
  AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");
  man->SetRun(run);

  AliCDBEntry *ent1 = man->Get("VZERO/Trigger/Data");
  AliVZEROTriggerData *fTriggerData = (AliVZEROTriggerData*)ent1->GetObject();

  Bool_t isRun2 = (AliCDBManager::Instance()->GetRun() >= 215011);

  for (int i=0; i<AliVZEROTriggerData::kNCIUBoards; i++) {
    AliVZEROLogicalSignal clk1BB(fTriggerData->GetClk1Win1(i),(UInt_t)fTriggerData->GetDelayClk1Win1(i),isRun2);
    AliVZEROLogicalSignal clk2BB(fTriggerData->GetClk2Win1(i),(UInt_t)fTriggerData->GetDelayClk2Win1(i),isRun2);
    AliVZEROLogicalSignal bbGate(clk1BB & clk2BB);

    AliVZEROLogicalSignal clk1BG(fTriggerData->GetClk1Win2(i),(UInt_t)fTriggerData->GetDelayClk1Win2(i),isRun2);
    AliVZEROLogicalSignal clk2BG(fTriggerData->GetClk2Win2(i),(UInt_t)fTriggerData->GetDelayClk2Win2(i),isRun2);
    AliVZEROLogicalSignal bgGate(clk1BG & clk2BG);

    printf("Board=%d\n",i);
    printf("  Win1: Clk1=%d DelayClk1=%d (%.2f -> %.2f) Clk2=%d DelayClk2=%d (%.2f -> %.2f) Latch=%d Reset=%d  Start=%.2f Stop=%.2f\n",
	   fTriggerData->GetClk1Win1(i), fTriggerData->GetDelayClk1Win1(i),
	   clk1BB.GetStartTime(), clk1BB.GetStopTime(),
	   fTriggerData->GetClk2Win1(i), fTriggerData->GetDelayClk2Win1(i),
	   clk2BB.GetStartTime(), clk2BB.GetStopTime(),
	   fTriggerData->GetLatchWin1(i), fTriggerData->GetResetWin1(i),
	   bbGate.GetStartTime(), bbGate.GetStopTime());
    printf("  Win2: Clk1=%d DelayClk1=%d (%.2f -> %.2f) Clk2=%d DelayClk2=%d (%.2f -> %.2f) Latch=%d Reset=%d  Start=%.2f Stop=%.2f\n\n",
	   fTriggerData->GetClk1Win2(i), fTriggerData->GetDelayClk1Win2(i),
	   clk1BG.GetStartTime(), clk1BG.GetStopTime(),
	   fTriggerData->GetClk2Win2(i), fTriggerData->GetDelayClk2Win2(i),
	   clk2BG.GetStartTime(), clk2BG.GetStopTime(),
	   fTriggerData->GetLatchWin2(i), fTriggerData->GetResetWin2(i),
	   bgGate.GetStartTime(), bgGate.GetStopTime());
  }

  printf("What will be used in MC:\n");
  AliVZEROTriggerSimulator simulator;
  simulator.Print();
}
