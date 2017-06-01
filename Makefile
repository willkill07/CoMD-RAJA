USE_CPLUSPLUS = ON
DOUBLE_PRECISION = ON

# Compiler Flags -- Edit to your desire
CXX := g++
CXXFLAGS := -std=c++11 -fopenmp -O3 -march=native -fno-exceptions -fno-rtti -Wall -Wextra -Wpedantic

CC := gcc
CFLAGS := -std=c99 -fopenmp -O3 -march=native

CPPFLAGS := -MMD

# Check for double precision
ifeq ($(DOUBLE_PRECISION), ON)
CPPFLAGS += -DDOUBLE
else
CPPFLAGS += -DSINGLE
endif

ifeq ($(USE_CPLUSPLUS), ON)
SRCEXT := .cpp
else
SRCEXT := .c
endif

SRCDIR := src
INCDIR := include
OBJDIR := obj

CPPFLAGS += -I$(INCDIR)

BIN := CoMD
SRC := $(wildcard $(SRCDIR)/*$(SRCEXT))
OBJ := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SRC:$(SRCEXT)=.o))
DEP := $(OBJ:.o=.d)
LDLIBS := -lm

ifeq ($(USE_CPLUSPLUS), ON)
LINK := $(LINK.cc)
COMPILE := $(COMPILE.cc)
else
LINK := $(LINK.c)
COMPILE := $(COMPILE.c)
endif

$(BIN) : $(OBJ)
	$(LINK) $^ $(LOADLIBES) $(LDLIBS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%$(SRCEXT) | $(OBJDIR)
	$(COMPILE) $(OUTPUT_OPTION) $<

$(OBJDIR) :
	@mkdir -p $(OBJDIR)

clean :
	-$(RM) $(OBJ) $(BIN)

veryclean : clean
	-$(RM) -r $(OBJDIR)

-include $(DEP)
