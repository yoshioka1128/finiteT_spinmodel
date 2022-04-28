# make file
TARGET=magcurve
FC=ifort
OBJS=magcurve_statistical.f90\
~/research/subroutine/fundip.f90\
~/research/subroutine/dipole.f90\
~/research/subroutine/plot_configxyz_open.f90\
~/research/subroutine/grainboundary.f90\
~/research/subroutine/defectR.f90\
~/research/subroutine/magnetization.f90\
~/research/subroutine/reset_position.f90\
~/research/subroutine/Intexch.f90\
~/research/subroutine/forward.f90\
~/research/subroutine/engdens.f90\
~/research/subroutine/plot_engdens.f90\
~/research/subroutine/plot_engdens_unit.f90\
~/research/subroutine/Olm.f90\
~/research/subroutine/momentR_exp.f90\
~/research/subroutine/moment_op.f90\
~/research/subroutine/Wigner3j.f90\
~/research/subroutine/Wigner6j.f90\
~/research/subroutine/Wignersub.f90\
~/research/subroutine/vector_to_angle.f90\
~/research/subroutine/MsT_Kuzmin.f90\
~/research/subroutine/rareearth.f90\
~/research/subroutine/rareearth_SmL3.f90\
~/research/subroutine/rdmatrix.f90\
~/research/subroutine/rdmatrix_SmL3.f90\
~/research/subroutine/stfactor.f90\
~/research/subroutine/numfactor.f90\
~/research/subroutine/CkqL.f90\
~/research/subroutine/CkqS.f90\
~/research/subroutine/CkqJ.f90\
~/research/subroutine/OkqL.f90\
~/research/subroutine/xifactor.f90\
~/research/subroutine/extendedBJ2.f90\
~/research/subroutine/extendedLJ2.f90\
~/research/subroutine/TVJl.f90\
~/research/subroutine/BESSIK.f90\
~/research/subroutine/BESCHB.f90\
~/research/subroutine/CHEBEV.f90\
~/research/subroutine/magcurve_analytical.f90\
~/research/subroutine/Hexch.f90\
~/research/subroutine/tlm.f90

FLAGS=-llapack -lblas -mcmodel=large
.SUFFIXES: .f90 .o
	$(FC) 
${TARGET} : $(OBJS)
	$(FC) -o $@ $(OBJS) ${FLAGS}
