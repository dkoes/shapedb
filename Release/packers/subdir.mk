################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../packers/FullMergePacker.cpp \
../packers/GreedyPacker.cpp \
../packers/MatcherPacker.cpp \
../packers/Packer.cpp \
../packers/SpectralPacker.cpp 

OBJS += \
./packers/FullMergePacker.o \
./packers/GreedyPacker.o \
./packers/MatcherPacker.o \
./packers/Packer.o \
./packers/SpectralPacker.o 

CPP_DEPS += \
./packers/FullMergePacker.d \
./packers/GreedyPacker.d \
./packers/MatcherPacker.d \
./packers/Packer.d \
./packers/SpectralPacker.d 


# Each subdirectory must supply rules for building sources it contributes
packers/%.o: ../packers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/openbabel-2.0 -I/usr/local/include/rdkit -I/usr/include/eigen3 -I"/home/dkoes/git/shapedb" -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


