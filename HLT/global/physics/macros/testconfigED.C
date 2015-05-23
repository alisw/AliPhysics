void testconfigED(const char* parent = "GLOBAL-esd-converter", const char* config="")
{
  printf("****testconfigED\n");
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  AliHLTConfiguration zmqsink("ZMQsink" , "ZMQsink" , parent , config);
  printf("****testconfigED, printing task list\n");
  pHLT->PrintTaskList();
}
