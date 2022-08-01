CC=g++ -std=c++17
CFLAGS=-Wall -Wextra -Werror
SOURCES=s21_matrix_oop.cpp
BUILD_PATH=./
LIBO=s21_matrix_oop.o
LIBA=s21_matrix_oop.a
EXE=test.o

all: clean s21_matrix_oop.a test check

s21_matrix_oop.a:
	$(CC) $(CFLAGS) -c $(SOURCES) -o $(BUILD_PATH)$(LIBO)
	ar rcs $(LIBA) $(LIBO)
	ranlib $(LIBA)
	
test:
	$(CC) $(CFLAGS) *.cpp -o $(EXE) -lgtest -lgtest_main
	./test.o

check:
	cppcheck *.cpp
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	python3 ../materials/linters/cpplint.py --extension=cpp *.cpp *.h

clean:
	-rm -rf *.o && rm -rf *.gcno
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.a *.o
	-rm -rf *.out