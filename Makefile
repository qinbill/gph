# Copyright 2014 UNSW CSE. All Rights Reserved.
# Author: jqin@cse.unsw.edu.au (Jianbin Qin)

PROD 	:= DEBUG
#PROD	:= DEBUG
OPT     := -O3
VERSION := \"0.0.1.0_${PROD}\"
TARGETS := mihapp histgenapp histalloapp histdump parthist mihistapp truehist dataproj dimshuff
TESTERS := 
# index_tester search_tester
DEFINES := -DUSEDENSE -DREDUCTION
# -DGREEDYALLO 
# -DGREEDYALLO
# -DREDUCTION
#-DDWORD32BITS
#-DLB_DEBUG 
#-DCEIL_ERROR
#-DEH_DEBUG

SRCS    := array32.cpp bucket_group.cpp linscan.cpp mihasher.cpp sparse_hashtable.cpp hmdata.cpp mihapp.cpp histgram.cpp dense_hashtable.cpp histgenapp.cpp histalloapp.cpp histdump.cpp parthist.cpp histallo.cpp mihistapp.cpp truehist.cpp dataproj.cpp mlregressor.cpp
HEADERS := array32.h bitarray.h bitops.h bucket_group.h linscan.h memusage.h mihasher.h sparse_hashtable.h types.h hmdata.h histgram.h dense_hashtable.h histallo.h mlregressor.h
OBJS    := ${SRCS:.cpp=.o}

CCFLAGS = ${OPT} -Wall -Wno-deprecated -pg -ggdb -D${PROD} ${DEFINES} -DVERSION=${VERSION}  -std=c++11 -O3 -I..
LDFLAGS = ${OPT} -ggdb  ${LIBS} -pg -std=c++11 
TESTFLAGS = -DTESTER_MAIN -DWORD32BITS

LIBS    = -lcrypto
CC	= g++


.PHONY: all clean
all:: ${TARGETS} 

test:: ${TESTERS}

mihapp:mihapp.o array32.o bucket_group.o linscan.o mihasher.o sparse_hashtable.o hmdata.o histgram.o dense_hashtable.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

mihistapp:mihistapp.o array32.o bucket_group.o linscan.o mihasher.o sparse_hashtable.o hmdata.o histgram.o dense_hashtable.o histallo.o mlregressor.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

dataproj:hmdata.o dataproj.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

histgenapp:histgenapp.o hmdata.o histgram.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

dimshuff:dimshuff.cpp hmdata.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

histalloapp:histalloapp.o hmdata.o histgram.o histallo.o mlregressor.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

histdump:histdump.o histgram.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

parthist:parthist.o hmdata.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^

truehist:truehist.o hmdata.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^
mlregressor:miregressor.o
	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^


# mihmaster.so:
# 	${CC} ${CCFLAGS} ${LDFLAGS} ${TESTFLAGS} -o $@ $^ -shared -fPIC

${OBJS}: %.o: %.cpp
	${CC} ${CCFLAGS} -o $@ -c $< 

clean:: 
	-rm -rf *~ *.o ${TARGETS} ${TESTERS} *.dSYM
