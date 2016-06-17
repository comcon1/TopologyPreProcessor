LIBS         = -lboost_program_options -lboost_unit_test_framework -lboost_test_exec_monitor -lopenbabel -lmysqlpp
CXXINC       = -I../include -I/usr/include/openbabel-2.0 -I/usr/include/mysql++ -I.. -I../include -DHAVE_CONFIG_H
CXXFLAGS     = -O0 -pipe -DTPP_UNIT_TEST \
 -DBUILD_DATE=\"`date +"%d-%m-%Y_%X"`\"  -DMYSQLPP_MYSQL_HEADERS_BURIED
CC           = gcc
CXX          = g++ $(CXXFLAGS) $(CXXINC)
LD           = $(CXX)

all: test_matrix 
	
	#test_global topio_test test_align 

test_matrix:  test_matrix.cpp ../src/global.o
	$(LD) -o test_matrix test_matrix.cpp ../src/global.o $(LIBS) -llapack

test_global: $(OBJDIR)/global.o test_global.o test_calc.o $(OBJDIR)/calc.o
	$(LD) -o test_global $(OBJDIR)/global.o $(OBJDIR)/calc.o test_global.o test_calc.o $(LIBS)
test_global.o: test_global.cpp
	$(CXX) -o test_global.o -c test_global.cpp
test_calc.o: test_calc.cpp test_calc.hpp
	$(CXX) -o test_calc.o -c test_calc.cpp

topio_test: $(OBJDIR)/global.o $(OBJDIR)/topio.o $(OBJDIR)/pdbutils.o $(OBJDIR)/strutil.o topio_test.o 
	$(CXX) -o topio_test topio_test.o $(OBJDIR)/global.o $(OBJDIR)/topio.o $(OBJDIR)/strutil.o $(OBJDIR)/pdbutils.o $(LIBS)
topio_test.o: topio_test.cpp
	$(CXX) -o topio_test.o -c topio_test.cpp

test_align: $(OBJDIR)/db_scanner.o $(OBJDIR)/smartfit.o $(OBJDIR)/global.o $(OBJDIR)/topio.o $(OBJDIR)/db_bonds.o $(OBJDIR)/strutil.o test_align.o
	$(CXX) -o test_align test_align.o $(OBJDIR)/global.o $(OBJDIR)/topio.o $(OBJDIR)/smartfit.o $(OBJDIR)/db_scanner.o $(OBJDIR)/db_bonds.o $(OBJDIR)/strutil.o $(LIBS)
test_align.o: test_align.cpp ../include/db_scanner.hpp ../include/global.hpp
	$(CXX) -o test_align.o -c test_align.cpp

$(OBJDIR)/*.o:
	cd .. && make @ && cd tests
clean: 
	rm -f test_global.o topio_test.o test_align.o test_global topio_test test_align test_calc test_calc.o
