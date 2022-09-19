TARGET1   =  PM25-stat
OBJECTS1   = PM25-stat.f90

TARGET2   =  O3-stat
OBJECTS2   = O3-stat.f90

#F_FLAGS   = -O -FR -no-wrap-margin -qopenmp
F_FLAGS   = -no-wrap-margin -qopenmp

FC        = ftn 

LIBS      = ${WGRIB2_LIB}
INC       = -I${WGRIB2_INC}


all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJECTS1)
	$(FC) $(F_FLAGS) -o $@ $(INC) $(OBJECTS1) $(LIBS)


$(TARGET2): $(OBJECTS2)
	$(FC) $(F_FLAGS) -o $@ $(INC) $(OBJECTS2) $(LIBS)

clean:
	rm -f *.o

clobber: clean
	-rm -f $(TARGET1) $(TARGET2)
