include ../mcl/common.mk

SRC=edit.cpp
ALL_SRC=edit.cpp edit_test.cpp

TARGET=edit
all: $(TARGET)

CFLAGS+=-I../mcl/include -std=c++14
ifeq ($(OS),mac)
  CFLAGS+=-Xpreprocessor -fopenmp
  LDFLAGS+=-lomp
else
  CFLAGS+=-fopenmp
  LDFLAGS+=-fopenmp
endif
%.o: %.cpp edit.hpp
	$(PRE)$(CXX) $(CFLAGS) -c $< -o $@ -MMD -MP -MF $(@:.o=.d)

$(TARGET): edit.o
	$(PRE)$(CXX) $< -o $@ ../mcl/lib/libmcl.a $(LDFLAGS)

edit_test: edit_test.o
	$(PRE)$(CXX) $< -o $@ ../mcl/lib/libmcl.a $(LDFLAGS) -lpthread

test: edit_test
	time ./edit_test -n 4

clean:
	rm -rf $(TARGET) *.o edit_test *.d

DEPEND_FILE=$(ALL_SRC:.cpp=.d)
-include $(DEPEND_FILE)

# don't remove these files automatically
.SECONDARY: $(ALL_SRC:.cpp=.o)

