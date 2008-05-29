// $Id$
/**
 * Analyze RCU RAW data ddl files.
 *
 * Usage:
 *   aliroot -b -q -l check-ddl.C'("filename")' | tee check-ddl.C
 *
 * Check an RCU ddl file and print out the channel content by using
 * the AliAltroDecoder.
 *
 * \b Note: The channel data has to be read in reverse order and is:
 * - length of a bunch
 * - end time bin of the bunch
 * - bunch data
 *
 * E.g.
 * <pre>
 * ---------- channel ----------
 * 10       314     3       300
 * 432      3
 * </pre>
 * means two bunches, one at 432 signal 300 and one at 314 signal 10.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_rcu
 */
void check_ddl(const char* filename=NULL)
{
  if (!filename) {
    cout << "usage: aliroot -b -l -q check-ddl.C'(\"filename\")'" << endl;
    return;
  }

  TString param=filename;
  param+="?filetype=raw";
  TFile file(param);
  if (file.IsZombie()) {
    cout << "can not open file " << filename << endl;
    return;
  }

  TArrayC buffer(file.GetSize());
  if (file.ReadBuffer(buffer.GetArray(), buffer.GetSize())) {
    cout << "error reading file " << filename << endl;
    return;
  }

  AliAltroDecoder decoder;
  if (decoder.SetMemory((UChar_t*)buffer.GetArray(), (UInt_t)buffer.GetSize())<0) {
    cout << "error setting up decoder " << endl;
    return;
  }

  if (!decoder.Decode()) {
    cout << "error decoding file " << filename << endl;
    return;    
  }

  cout << "RCU trailer size: " << decoder.GetRCUTrailerSize() << endl;

  AliAltroData channel;
  while (decoder.NextChannel(&channel)) {
    cout << "---------- channel " << channel.GetHadd() << " ----------" << endl;
    decoder.PrintInfo(channel, channel.GetDataSize(), 4);
  }
}
