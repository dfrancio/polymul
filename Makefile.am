AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = gentab.py

include_HEADERS = \
	include/polymul.h \
	include/polymul_tab.h

TESTS = \
	polymul_test \
	benchmark \
	polymul_with_tab_test \
	benchmark_with_tab

check_PROGRAMS = $(TESTS)

AM_CPPFLAGS = -Wall -I$(top_srcdir)/include

polymul_test_SOURCES = tests/polymul_test.cpp
benchmark_SOURCES = tests/benchmark.cpp

polymul_with_tab_test_SOURCES = tests/polymul_test.cpp
polymul_with_tab_test_CPPFLAGS = -DPOLYMUL_USE_TABLES

benchmark_with_tab_SOURCES = tests/benchmark.cpp
benchmark_with_tab_CPPFLAGS = -DPOLYMUL_USE_TABLES
