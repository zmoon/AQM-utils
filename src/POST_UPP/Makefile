TARGET1   =  PM25-stat

OBJECTS1   = PM25-stat.f90

#F_FLAGS   = -O -FR -no-wrap-margin -qopenmp
F_FLAGS   = -no-wrap-margin -qopenmp

FC        = ftn 

LIBS      = ${WGRIB2_LIB}
INC       = -I${WGRIB2_INC}

$(TARGET1): $(OBJECTS1)
	$(FC) $(F_FLAGS) -o $@ $(INC) $(OBJECTS1) $(LIBS)

clean: 
	rm -f *.o 
# install:
	# -cp $(TARGET1) ../../../bin/

clobber: clean
	-rm -f $(TARGET1)
