#!/usr/bin/python
import zmq,sys
import signal
import re
import time
import binascii
import o2header

def Exit_gracefully(signal, frame):
  context.destroy()
  print(" signal caught, exiting")
  sys.exit(0)
signal.signal(signal.SIGINT, Exit_gracefully)

endpoint="SUB>tcp://localhost:60202"
theMessage=[]

if len(sys.argv)==1:
    print "A simple debugging/monitoring tool for ZMQ sockets"
    print "supports multi-part messages (can both send/receive)"
    print "Usage:"
    print sys.argv[0]+" \"PULL@tcp://localhost:60202\""
    print sys.argv[0]+" \"SUB>tcp://localhost:60202\" \"optional subscription\""
    print sys.argv[0]+" \"REQ>tcp://localhost:60202\" INFO \"-reset someoption=value\""
    print sys.argv[0]+" \"REP>tcp://localhost:60202\" INFO \"-reset someoption=value\""
    quit()
if len(sys.argv)>1:
    endpoint=sys.argv[1]
    sepindex=endpoint.find('@')
    if sepindex < 0:
      sepindex=endpoint.find('>')
    if sepindex < 0:
      sepindex=endpoint.find('+')
    if sepindex < 0:
      sepindex=endpoint.find('-')
    if sepindex < 0 or sepindex > 5:
      print("error in endpoint syntax: "+endpoint)
      quit()
    mode=endpoint[0:sepindex]
    endpoint=endpoint[sepindex:]

print "mode: "+mode+" endpoint: "+endpoint

if len(sys.argv)>2:
  theMessage=sys.argv[2:]

if len(theMessage)<1:
  theMessage = [ "","" ]
if len(theMessage)%2==1:
  theMessage.append("")

#  Prepare our context and sockets
context = zmq.Context()

if mode=="PUSH":
    socket = context.socket(zmq.PUSH)
elif mode=="PULL":
    socket = context.socket(zmq.PULL)
elif mode=="SUB":
    socket = context.socket(zmq.SUB)
    socket.setsockopt(zmq.SUBSCRIBE, theMessage[0])
elif mode=="PUB":
    socket = context.socket(zmq.PUB)
elif mode=="REQ":
    socket = context.socket(zmq.REQ)
elif mode=="REP":
    socket = context.socket(zmq.REP)
else:
    print "not a valid mode"
    quit()
#socket.set(zmq.LINGER,10)
if endpoint[0]=='>' or endpoint[0]=='-' or endpoint[0]=='+':
    print "connect to: "+endpoint[1:]
    socket.connect(str(endpoint[1:]))
elif endpoint[0]=='@':
    print "bind to: "+endpoint[1:]
    socket.bind(str(endpoint[1:]))

# fill in the header information
for idx in range(0, len(theMessage)):
    if idx%2:
        continue
    payloadSize=len(theMessage[idx+1])
    origin=theMessage[idx][8:12]
    if origin[:3]=="***":
        origin=origin[:3]+'\0'
    theMessage[idx]=o2header.make(theMessage[idx],origin,payloadSize);

#avoid late subscriber syndrome
if mode=="PUB":
  time.sleep(0.5)

endless=True
if mode=="REQ" or mode=="PUSH" or mode=="PUB":
  endless=False
  socket.send_multipart(theMessage)

while True:
  if mode=="SUB" or mode=="PULL" or mode=="REP" or mode=="REQ":
    #raw_input("press a key...")
    msg = socket.recv_multipart();
    print "###################################################"
    i=0;
    for message in msg:
      dirty = str(message)[0:2000]
      clean = re.sub('[^\s!-~]', '.', dirty)
      if i==0:
        print "topic: "+clean
        o2header.dump(message)
        i=1
      elif i==1:
        print "message size: "+str(sys.getsizeof(message))
        print clean
        print "___________________________________________________"
        i=0
  if mode=="REP":
    socket.send_multipart(theMessage)
  if not endless:
      break
