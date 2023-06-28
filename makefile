# Variables
f90comp = ifort
switch = -auto-scalar -heap-arrays 1024
#switch = -warn unused -warn uncalled -warn interfaces -O3 -ipo -auto-scalar -heap-arrays 1024
#f90comp = gfortran 
#switch = -fbacktrace -ffree-line-length-0 -Wall -fno-automatic

objects = gear_GlobVARs.o t_dgls.o uawp.o fgauss.o gear.o random_routines.o libdate.o hypevar.o modvar.o worvar.o general_wc.o general_func.o convert.o time.o hypetypes.o readwrite.o hype_indata.o atm_proc.o irrigation.o hype_wbout.o npc_soil_proc.o soil_proc.o regional_groundwater.o sw_proc.o npc_sw_proc.o soilmodel0.o glacier_soilmodel.o model_hype.o compout.o data.o optim.o statedata.o main.o
modfiles = gear_GlobVARs.mod t_dgls.mod uawp.mod fgauss.mod gear_implicit.mod random_routines.mod libdate.mod hypevariables.mod modvar.mod worldvar.mod general_water_concentration.mod general_functions.mod convert.mod timeroutines.mod statetype_module.mod readwrite_routines.mod hype_indata.mod atmospheric_processes.mod irrigation_module.mod hype_waterbalance.mod npc_soil_processes.mod soil_processes.mod regional_groundwater_module.mod surfacewater_processes.mod npc_surfacewater_processes.mod soilmodel_default.mod glacier_soilmodel.mod modelmodule.mod compout.mod datamodule.mod state_datamodule.mod

# Makefile
hype:	$(objects)
	$(f90comp) -o hype $(switch) $(objects)

# All .o files are made from corresponding .f90 files
%.o:	%.f90
	$(f90comp) -c $(switch) $<
%.o:	%.F90
	$(f90comp) -c $(switch) $<
%.mod:	%.f90
	$(f90comp) -c $(switch) $<
%.mod:	%.F90
	$(f90comp) -c $(switch) $<

# Dependencies
modvar.o       : libdate.mod
convert.o      : libdate.mod
atm_proc.o     : hypevariables.mod modvar.mod hype_indata.mod
irrigation.o   : hypevariables.mod modvar.mod statetype_module.mod
npc_soil_proc.o: hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod general_functions.mod
soil_proc.o    : gear_GlobVARs.mod t_dgls.mod uawp.mod fgauss.mod gear_implicit.mod hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod general_functions.mod npc_soil_processes.mod atmospheric_processes.mod hype_indata.mod
npc_sw_proc.o  : hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod general_functions.mod
sw_proc.o      : hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod general_functions.mod soil_processes.mod
soilmodel0.o   : gear_GlobVARs.mod t_dgls.mod uawp.mod fgauss.mod gear_implicit.mod hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod general_functions.mod npc_soil_processes.mod soil_processes.mod atmospheric_processes.mod regional_groundwater_module.mod irrigation_module.mod
glacier_soilmodel.o : hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod npc_soil_processes.mod soil_processes.mod atmospheric_processes.mod regional_groundwater_module.mod
regional_groundwater.o : hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod npc_soil_processes.mod 
model_hype.o   : libdate.mod hypevariables.mod modvar.mod statetype_module.mod general_water_concentration.mod glacier_soilmodel.mod soilmodel_default.mod soil_processes.mod npc_soil_processes.mod surfacewater_processes.mod npc_surfacewater_processes.mod irrigation_module.mod regional_groundwater_module.mod hype_waterbalance.mod hype_indata.mod atmospheric_processes.mod
compout.o      : libdate.mod modvar.mod worldvar.mod timeroutines.mod convert.mod
worldvar.o     : libdate.mod modvar.mod
readwrite.o    : libdate.mod worldvar.mod convert.mod
hype_wbout.o   : libdate.mod readwrite_routines.mod worldvar.mod modvar.mod
hype_indata.o  : libdate.mod readwrite_routines.mod worldvar.mod modvar.mod
timeroutines.o : libdate.mod worldvar.mod modvar.mod
data.o         : libdate.mod modvar.mod worldvar.mod convert.mod timeroutines.mod readwrite_routines.mod compout.mod modelmodule.mod 
optim.o        : random_routines.mod libdate.mod modvar.mod worldvar.mod statetype_module.mod timeroutines.mod modelmodule.mod compout.mod datamodule.mod
statedata.o    : libdate.mod modvar.mod worldvar.mod statetype_module.mod modelmodule.mod readwrite_routines.mod
main.o         : gear_GlobVARs.mod t_dgls.mod uawp.mod fgauss.mod gear_implicit.mod libdate.mod modvar.mod worldvar.mod statetype_module.mod timeroutines.mod readwrite_routines.mod modelmodule.mod compout.mod datamodule.mod state_datamodule.mod


.PHONY : clean
clean:	
	rm -f $(objects)
	rm -f $(modfiles)
