################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
$(ROOT)/baseRNG.cpp \
$(ROOT)/main.cpp 

OBJS += \
./baseRNG.o \
./main.o 

DEPS += \
${addprefix ./, \
baseRNG.d \
main.d \
}


# Each subdirectory must supply rules for building sources it contributes
%.o: $(ROOT)/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: icc C++ Compiler'
	@echo icc -O2   -c  -o$@ $<
	@icc -O2   -c  -o$@ $< && \
	echo -n $(@:%.o=%.d) $(dir $@) > $(@:%.o=%.d) && \
	icc -MM -MG -P -w -O2   -c  $< >> $(@:%.o=%.d)
	@echo 'Finished building: $<'
	@echo ' '


