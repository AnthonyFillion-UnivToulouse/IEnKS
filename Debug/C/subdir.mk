################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../C/main.cpp 

OBJS += \
./C/main.o 

CPP_DEPS += \
./C/main.d 


# Each subdirectory must supply rules for building sources it contributes
C/%.o: ../C/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++11 -I../C/armadillo-9.200.6/include -I../C/cnpy-master -O2 -g3 -Wall -c -fmessage-length=0 -DARMA_DONT_USE_WRAPPER -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


