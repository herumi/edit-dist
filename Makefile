MCL_USE_OMP=1
include ../mcl/common.mk

SRC=edit.cpp
ALL_SRC=edit.cpp edit_test.cpp

TARGET=edit
all: $(TARGET)

CFLAGS+=-I../mcl/include -std=c++11
ifeq ($(MCL_USE_PROF),2)
  LDFLAGS+=-L /opt/intel/vtune_amplifier/lib64 -ljitprofiling -ldl
endif

MCLLIB=../mcl/lib/libmcl.a

$(MCLLIB):
	$(MAKE) -C ../mcl lib/libmcl.a

%.o: %.cpp edit.hpp
	$(PRE)$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MP -MF $(@:.o=.d)

$(TARGET): edit.o $(MCLLIB)
	$(PRE)$(CXX) $< -o $@ $(MCLLIB) $(LDFLAGS)

edit_test: edit_test.o $(MCLLIB)
	$(PRE)$(CXX) $< -o $@ $(MCLLIB) $(LDFLAGS) -lpthread

test: edit_test
	time ./edit_test -n 4

clean:
	rm -rf $(TARGET) *.o edit_test *.d

DEPEND_FILE=$(ALL_SRC:.cpp=.d)
-include $(DEPEND_FILE)

# don't remove these files automatically
.SECONDARY: $(ALL_SRC:.cpp=.o)

