################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../GSSTreeCreator.cpp \
../GSSTreeSearcher.cpp \
../GSSTreeStructures.cpp \
../KSamplePartitioner.cpp \
../MGrid.cpp \
../MappableOctTree.cpp \
../MemMapped.cpp \
../ShapeDistance.cpp \
../WorkFile.cpp \
../main.cpp 

OBJS += \
./GSSTreeCreator.o \
./GSSTreeSearcher.o \
./GSSTreeStructures.o \
./KSamplePartitioner.o \
./MGrid.o \
./MappableOctTree.o \
./MemMapped.o \
./ShapeDistance.o \
./WorkFile.o \
./main.o 

CPP_DEPS += \
./GSSTreeCreator.d \
./GSSTreeSearcher.d \
./GSSTreeStructures.d \
./KSamplePartitioner.d \
./MGrid.d \
./MappableOctTree.d \
./MemMapped.d \
./ShapeDistance.d \
./WorkFile.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/openbabel-2.0 -I/usr/local/include/rdkit -I/usr/include/eigen3 -I"/home/dkoes/git/shapedb" -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


