################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../srcf/entropy/entropy.cpp 

OBJS += \
./srcf/entropy/entropy.o 

CPP_DEPS += \
./srcf/entropy/entropy.d 


# Each subdirectory must supply rules for building sources it contributes
srcf/entropy/%.o: ../srcf/entropy/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


