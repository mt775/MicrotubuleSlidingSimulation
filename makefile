objects = imain.o rdata.o init.o mkMT.o overlapMap.o overlap.o equations.o\
newpos.o functions.o mprove.o nrutil.o ran1.o ludcmp.o lubksb.o \
iunit.o agr_funcs.o calcav.o write_output.o normalize.o sampleExpDist.o lambert_w.o\
update_av.o gasdev.o ran2.o compact.o dlinmin.o frprmn2.o f1dim.o\
mnbrak.o dbrent.o df1dim.o timecorr.o motorvec.o boundary.o getjMT.o percolation.o \
special_plots.o MTflip.o initialize_variables.o free_all.o actin_overlapMap.o prob_update.o

     MTprog  : $(objects)
	 gcc -mmacosx-version-min=10.6 -framework Accelerate -g -O0 -g3 -Wall -o MT $(objects)
	rm $(objects)

     %.o: %.c
	gcc  -mmacosx-version-min=10.6  -g -O0 -g3 -c $<


     imain.o :   global_var.h
     getjMT.o :  global_var.h
     init.o :   global_var.h nrutil.h
     MTflip.o :   global_var.h nrutil.h
     rdata.o :   global_var.h ran1.h
     mkMT.o :   global_var.h nrutil.h ran1.h sampleExpDist.h lambert_w.h
     overlapMap.o :  global_var.h ran1.h
     overlap.o :  global_var.h
     equations.o :  global_var.h nrutil.h
     newpos.o :  global_var.h
     motorvec.o :  global_var.h
     agr_funcs.o : global_var.h
     calcav.o : global_var.h
     mprove.o : nrutil.h
     ludcmp.o : nrutil.h
     write_output.o : global_var.h nrutil.h
     normalize.o : global_var.h
     update_av.o : global_var.h
     functions.o : global_var.h nrutil.h
     compact.o : global_var.h nrutil.h ran1.h
     timecorr.o : global_var.h
     boundary.o : global_var.h ran1.h
     frprmn2.o : nrutil.h
     dlinmin.o : nrutil.h
     f1dim.o : nrutil.h
     mnbrak.o : nrutil.h
     dbrent.o : nrutil.h
     df1dim.o : nrutil.h
     percolation.o : global_var.h
     special_plots.o : global_var.h
     initialize_variables.o : global_var.h nrutil.h
     free_all.o : global_var.h nrutil.h
     actin_overlapMap.o : global_var.h ran1.h
     prob_update.o : global_var.h

     .PHONY : clean
     clean :
	rm $(objects)
	clean baby
