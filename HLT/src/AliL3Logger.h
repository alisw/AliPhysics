class MLUCLogServer;
class ofstream;

class AliL3Logger{
  public:
  static int kAll;
  static int kDebug;
  static int kInformational;
  static int kWarning;
  static int kError;
  static int kFatal;
  AliL3Logger();
  ~AliL3Logger();
  void Set(int l);
  void UnSet(int l);
  void UseDevNull();
  void UseStdout();
  void UseStderr();
  void UseStream(char *name="AliLevel3.log");
  void NotUseDevNull();
  void NotUseStdout();
  void NotUseStderr();
  void NotUseStream();
  private:
  MLUCLogServer *dn;
  MLUCLogServer *so;
  MLUCLogServer *se;
  MLUCLogServer *sm;
  ofstream *of;
};


