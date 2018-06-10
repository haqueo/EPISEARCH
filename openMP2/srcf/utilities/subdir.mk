################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../srcf/utilities/utilities.cpp 

OBJS += \
./srcf/utilities/utilities.o 

CPP_DEPS += \
./srcf/utilities/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
srcf/utilities/%.o: ../srcf/utilities/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


