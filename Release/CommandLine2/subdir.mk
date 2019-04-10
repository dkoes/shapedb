################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../CommandLine2/CommandLine.cpp 

OBJS += \
./CommandLine2/CommandLine.o 

CPP_DEPS += \
./CommandLine2/CommandLine.d 


# Each subdirectory must supply rules for building sources it contributes
CommandLine2/%.o: ../CommandLine2/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/local/include/openbabel-2.0 -I/usr/local/include/rdkit -I/usr/include/eigen3 -I"/home/dkoes/git/shapedb" -O3 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


