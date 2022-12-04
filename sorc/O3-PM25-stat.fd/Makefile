TARGET   =  O3-PM25-stat

OBJECTS   = O3-PM25-stat.f90

#F_FLAGS   = -O -FR -no-wrap-margin -qopenmp
F_FLAGS   = -no-wrap-margin -qopenmp

FC        = ftn 

LIBS      = ${WGRIB2_LIB}
INC       = -I${WGRIB2_INC}

$(TARGET): $(OBJECTS)
	$(FC) $(F_FLAGS) -o $@ $(INC) $(OBJECTS) $(LIBS)

clean: 
	rm -f *.o 
# install:
	# -cp $(TARGET1) ../../../bin/

clobber: clean
	-rm -f $(TARGET)
