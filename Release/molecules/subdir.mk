################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molecules/MolSphere.cpp \
../molecules/OBMoleculeAnalytic.cpp \
../molecules/PMol.cpp \
../molecules/RDMoleculeAnalytic.cpp 

OBJS += \
./molecules/MolSphere.o \
./molecules/OBMoleculeAnalytic.o \
./molecules/PMol.o \
./molecules/RDMoleculeAnalytic.o 

CPP_DEPS += \
./molecules/MolSphere.d \
./molecules/OBMoleculeAnalytic.d \
./molecules/PMol.d \
./molecules/RDMoleculeAnalytic.d 


# Each subdirectory must supply rules for building sources it contributes
molecules/%.o: ../molecules/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/openbabel-2.0 -I/usr/local/include/rdkit -I/usr/include/eigen3 -I"/home/dkoes/git/shapedb" -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


