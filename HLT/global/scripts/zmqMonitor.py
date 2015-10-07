#!/usr/bin/python
import zmq,sys
import signal

def Exit_gracefully(signal, frame):
  context.destroy()
  print(" signal caught, exiting")
  sys.exit(0)
signal.signal(signal.SIGINT, Exit_gracefully)

endpoint="SUB>tcp://localhost:60202"
requestHeader="INFO"
requestBody=""
if len(sys.argv)==1:
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
    if sepindex < 0 or sepindex > 5:
      print("error in endpoint syntax: "+endpoint)
      quit()
    mode=endpoint[0:sepindex]
    endpoint=endpoint[sepindex:]
if len(sys.argv)>2:
    requestHeader=sys.argv[2]
if len(sys.argv)>3:
    requestBody=sys.argv[3]

print "mode: "+mode+" endpoint: "+endpoint

#  Prepare our context and sockets
context = zmq.Context()

if mode=="PUSH":
    socket = context.socket(zmq.PUSH)
elif mode=="PULL":
    socket = context.socket(zmq.PULL)
elif mode=="SUB":
    socket = context.socket(zmq.SUB)
    socket.setsockopt(zmq.SUBSCRIBE, requestHeader)
elif mode=="PUB":
    socket = context.socket(zmq.PUB)
elif mode=="REQ":
    socket = context.socket(zmq.REQ)
elif mode=="REP":
    socket = context.socket(zmq.REP)
else:
    print "not a valid mode"
    quit()
socket.set(zmq.LINGER,10)
if endpoint[0]=='>':
    print "connect to: "+endpoint[1:]
    socket.connect(str(endpoint[1:]))
elif endpoint[0]=='@':
    print "bind to: "+endpoint[1:]
    socket.bind(str(endpoint[1:]))

endless=True
if mode=="REQ" or mode=="PUSH" or mode=="PUB":
    endless=False
    socket.send(requestHeader,zmq.SNDMORE)
    socket.send(requestBody)
    print("sent: "+requestHeader+" "+requestBody)

while True:
    if mode=="SUB" or mode=="PULL" or mode=="REP" or mode=="REQ":
        msg = socket.recv_multipart();
        print "_________________________"
        i=0;
        for message in msg:
            if i==0:
                print "topic: "+str(message)
                i=1
            elif i==1:
                print "  messagesize: "+str(sys.getsizeof(message))
                i=0
    if mode=="REP":
        socket.send(requestHeader,zmq.SNDMORE)
        socket.send(requestBody)
        print "  sent back: "+requestHeader+" "+requestBody
    if not endless:
        break
