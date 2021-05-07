################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../C/Ini/ini.c 

CPP_SRCS += \
../C/Ini/INIReader.cpp 

OBJS += \
./C/Ini/INIReader.o \
./C/Ini/ini.o 

C_DEPS += \
./C/Ini/ini.d 

CPP_DEPS += \
./C/Ini/INIReader.d 


# Each subdirectory must supply rules for building sources it contributes
C/Ini/%.o: ../C/Ini/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I"/home/anthony/workspace/IEnKS/C/armadillo-9.200.6" -O2 -g3 -Wall -c -fmessage-length=0 -DARMA_DONT_USE_WRAPPER -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

C/Ini/%.o: ../C/Ini/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -O2 -g3 -p -pg -ftest-coverage -fprofile-arcs -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


