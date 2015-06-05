void testconfigZMQsink(const char* parent = "GLOBAL-esd-converter", const char* config="")
{
  printf("****testconfigZMQsink\n");
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  AliHLTConfiguration zmqsink("ZMQsink" , "ZMQsink" , parent , config);
  printf("****testconfigZMQsink, printing task list\n");
  pHLT->PrintTaskList();
}
