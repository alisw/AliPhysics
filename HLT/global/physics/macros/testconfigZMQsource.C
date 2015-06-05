void testconfigZMQsource(const char* parent = "GLOBAL-esd-converter", const char* config="")
{
  printf("****testconfigZMQsource\n");
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  AliHLTConfiguration zmqsink("ZMQsource" , "ZMQsource" , parent , config);
  printf("****testconfigZMQsource, printing task list\n");
  pHLT->PrintTaskList();
}
