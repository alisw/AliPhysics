#!/usr/bin/python
import zmq,time,sys

port="60201"
if (len(sys.argv)>1):
  port=sys.argv[1]

#  Prepare our context and sockets
context = zmq.Context()
socket = context.socket(zmq.PUB)
socket.bind('tcp://*:'+port)
while True:
  print "sending"
  socket.send("1st part",zmq.SNDMORE)
  socket.send("2nd part",zmq.SNDMORE)
  socket.send("3rd part")
  time.sleep(1)

