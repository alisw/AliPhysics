AliTaskCDBconnect* AddTaskCDBconnect(const char *path/*="raw://"*/, Int_t run=0) {
  return AliTaskCDBconnect::AddTaskCDBconnect(path, run);
}

AliTaskCDBconnect* AddTaskCDBconnect() {
  return AliTaskCDBconnect::AddTaskCDBconnect("cvmfs://");
}
