include ../../makefile.inc

CC=g++

CFLAGS=-w

#INC = -I/usr/local/include/eigen3
INC = -I${EIGENDIR}

LIB = \
-L/usr/lib \
-L/lib/libblas \
-llapack -lblas

SRCS = options.cpp parser.cpp cellinfo.cpp common.cpp ewald.cpp find_bounding_sphere.cpp generate_tvecs.cpp

OBJS = $(SRCS:.cpp=.o)

MAIN = PPP 

LIBNAME = libcommon.a

$(MAIN) : $(OBJS)
	ar -rcs $(LIBNAME) $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INC) $(LIB) -c $< -o $@

clean:
	rm *.o
