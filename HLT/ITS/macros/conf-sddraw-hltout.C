// $Id$
/**
 * @file conf-sddraw-hltout.C
 * @brief Configuration macro for sim-sddraw-hltout.C
 *
 * This is the configuration macro for the sim-hlt-rawddl.C example.
 * It defines only one component in the chain, an AliRawReaderPublisher
 * which publishes the raw data blocks according to the selection.
 * The data blocks are published with data type {DDLRAW :ITS }. The
 * AliRawReaderPublisherComponent sets automatically the data specification
 * from the equipment id. 
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_its
 */
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // the configuration
  TString writerInput;
  TString arg;

  // publisher configuration
  // see AliHLTRawReaderPublisherComponent for details
  arg.Form("-detector ITSSDD -skipempty -datatype 'DDL_RAW ' 'ITS ' -verbose");
  AliHLTConfiguration pubconf("publisher", "AliRawReaderPublisher", NULL , arg.Data());
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="publisher";

  // currently, no more components in the chain, the original data is just
  // forwarded to the HLTOUT
  // ----------------------------------------
  // add additional processing here if needed
  // Note: you have to change the name of the chain to run from 'publisher'
  // to the last component in your list
}
