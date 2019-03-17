# construct O2 compatible headers for ZMQ communication

import struct

#this describes the binary format of the O2 header
O2fmtString="<4sIII8s8s8s8s4sI8sQQ"

# this function creates a binary buffer in O2 DataHEader format
# with custom data type, origin and payload size
def make(dataType,origin,payloadSize):
    return struct.pack(
            O2fmtString
            ,"O2O2"                        # 4 magic number
            ,struct.calcsize(O2fmtString)  # 4 header size
            ,0                             # 4 flags
            ,1                             # 4 header version
            ,"DataHDR" + '\0'                    # 8 header type
            ,"NONE    "                    # 8 header serialization method
            ,""                            # 8 data description (1) - unused in HLT
            ,dataType                      # 8 data description (2)
            ,origin                        # 4 origin
            ,0                             # 4 reserved (checksum)
            ,"NONE    "                    # 8 data serialization
            ,0                             # 8 specification
            ,payloadSize                   # 8 payload size
            )

# this function just prints the header values according to the format string
def dump(header):
  magic = struct.unpack("<4s", header[0:4])
  if magic[0]=="O2O2":
    unpacked = struct.unpack(O2fmtString,header[0:80])
    print unpacked
  else:
    print ""

