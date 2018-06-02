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
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Users/Omar/Documents/Year4/M4R/muster/src -I/opt/local/include -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


