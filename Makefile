ifdef PLATFORM
include make.inc.$(PLATFORM)
else
include make.inc
endif

HDRS     = config.h measurement.h orbtimer.h types.h tests.h copyright.h comm.h options.h ugni_utils.h
OBJS     = measurement.o orbtimer.o comm.o net_test.o options.o sysconfidence.o bit_test.o ugni_utils.o

sysconfidence: $(OBJS)
	$(CC) $(CFLAGS) -o sysconfidence $(OBJS) $(LIBS)

sysconfidence.o: sysconfidence.c $(HDRS)
measurement.o:   measurement.c   $(HDRS)
options.o:       options.c       $(HDRS)
comm.o:          comm.c          $(HDRS)
orbtimer.o:      orbtimer.c      $(HDRS)
net_test.o:      net_test.c      $(HDRS)
bit_test.o:      bit_test.c      $(HDRS)
ugni_utils.o:    ugni_utils.c      $(HDRS)

config.h:
	echo No config.h
	echo you need to run scripts/config.sh 
	scripts/config.sh

clean: 
	rm -f *.o *.a sysconfidence read_xdd libxdd.a

distclean:
	make clean
	rm -f config.h
