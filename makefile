#EMSOL_HOME =  /home/agr/project/dietprob/EMSOL/v0_2_16_Julian
#EMSOL_HOME =  /home/jajhall/EMSOL/v0_2_16/
#EMSOL_HOME =  /home/agr/bin/EMSOL/v0_2_16_Julian/
EMSOL_HOME =  /home/agr/bin/EMSOL/v0_2_16a_Julian/

FORTRAN   = /sl61/opt/intel_fc_80/bin/ifc
#LLKOPT     = -L$(EMSOL_HOME)/lib -lemsol_fpe -static-libcxa
#LLKOPT     = -L$(EMSOL_HOME)/lib -lEMSOLdeb -static-libcxa
LLKOPT     = -L$(EMSOL_HOME)/lib -lEMSOL -static-libcxa
#LLKOPT     = -L$(EMSOL_HOME)/lib -lemsol -static-libcxa
#LLKOPT     = -L$(EMSOL_HOME)/lib -lEMSOL_ifc_deb
LINC       = -I$(EMSOL_HOME)/install/source -I$(EMSOL_HOME)/inc

OBJ= interface.o common_types.o nonlin.o\
     slpdrive.o rdprobdi.o diet_drv.o register.o rd_prob_da.o\
     solprint.o objfun.o gradient.o confun.o ranprint.o lamprint.o\
     filtersl.o diet_time.o slpmain.o phaseI.o ranging.o prt_mdl.o\
     int2char.o do_range.o\
     getsensitivity.o comfloat.o getindices.o getstatus.o\
     correctmixer.o

INC = msg.inc\
      return.inc\
      pusr.inc\
      SLPCOMM.INC
#MOD = common_types.o

OUT=diet

FFLAGS = -w -check bounds -check arg_temp_created -traceback -fpstkchk -heap-arrays -fpe0 -g
#FFLAGS = -w -g -d5 -CA -CB -CS -quiet 
#FFLAGS = -O

INCLUDE = 

all: $(OUT)

%.o:%.f $(INC)
	$(FORTRAN) $(FFLAGS) $(INCLUDE) $(LINC) -c $< -o $@

$(OUT): $(OBJ)
	$(FORTRAN) $(FFLAGS) -o $(OUT) $(OBJ) $(LLKOPT)

clean:
	rm *.o diet

