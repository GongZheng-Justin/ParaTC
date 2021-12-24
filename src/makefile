#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP       =  gcc_MPI#intel_MPI

# ============================ CFD part ============================
# Choose pre-processing options
#   -DOverWriteFFT : Overwrite the input DFT or not.
#   -DScalarFlow   : Calculate the scalar field or not.
#   -DSaveNode     : Save the velocity field on grid node or cell center.
#   -DSAVE_SINGLE  : Save binary data in single-precision.

# CFD options

exeName    =  ParaTC

CFD_dir    =  
CFDEdir    =  

CFD_DEFS   = -DOverWriteFFT -DScalarFlow #-DSaveNode
CFDEDEFS   = -DSAVE_SINGLE

CFD_inc    = -I./ThirdParty/fftw/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

# CFD source files
SrcT   := m_Decomp2d.f90
CFDEsrc:= $(addprefix $(CFDEdir), ${SrcT})

SrcT   := m_TypeDef.f90 m_LogInfo.f90 m_Timer.f90 m_Parameters.f90 m_Variables.f90 m_MeshAndMetries.f90 m_BC_and_Halo.f90  \
          m_Tools.f90 m_TScheme.f90  m_FlowCase.f90 m_Poisson.f90 m_IOAndVisu.f90 m_ChannelSystem.f90 main.f90
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
ifeq ($(CMP),intel_MPI)
  FortC = mpiifort
  CFLAGS= -fpp -O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large
else ifeq ($(CMP),gcc_MPI)
  FortC = mpif90
  CFLAGS= -cpp -O3 -funroll-loops -floop-optimize -g -Wall -fbacktrace -ffree-line-length-none
endif

all: $(exeName)
CFDEobj  = $(CFDEsrc:%.f90=%.o)
CFD_obj  = $(CFD_src:%.f90=%.o)
$(exeName):$(CFDEobj) $(CFD_obj)
	         $(FortC) $(CFLAGS) -o  $@ $(CFDEobj) $(CFD_obj) $(CFD_lib)
$(CFDEobj):$(CFDEdir)%.o :$(CFDEdir)%.f90
	         $(FortC) $(CFLAGS) $(CFDEDEFS) -c $<
$(CFD_obj):%.o : %.f90
	         $(FortC) $(CFLAGS) $(CFD_inc) $(CFD_DEFS) -c $<

.PHONY: clean
clean:
	rm -fr  *.o *.mod $(exeName)  $(CFD_dir)*o $(CFDEdir)*.o
