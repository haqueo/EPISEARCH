################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../srcf/searchTools/searchTools.cpp 

OBJS += \
./srcf/searchTools/searchTools.o 

CPP_DEPS += \
./srcf/searchTools/searchTools.d 


# Each subdirectory must supply rules for building sources it contributes
srcf/searchTools/%.o: ../srcf/searchTools/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


