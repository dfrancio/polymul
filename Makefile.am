AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = gentab.py
AM_CXXFLAGS = -Wall

include_HEADERS = \
	include/polymul.h \
	include/polymul_tab.h

TESTS = \
	polymul_test \
	benchmark \
	polymul_with_tab_test \
	benchmark_with_tab

check_PROGRAMS = $(TESTS)

polymul_test_SOURCES = tests/polymul_test.cpp
polymul_test_CPPFLAGS = -I$(top_srcdir)/include

benchmark_SOURCES = tests/benchmark.cpp
benchmark_CPPFLAGS = -I$(top_srcdir)/include

polymul_with_tab_test_SOURCES = tests/polymul_test.cpp
polymul_with_tab_test_CPPFLAGS = -I$(top_srcdir)/include -DPOLYMUL_USE_TABLES

benchmark_with_tab_SOURCES = tests/benchmark.cpp
benchmark_with_tab_CPPFLAGS = -I$(top_srcdir)/include -DPOLYMUL_USE_TABLES

noinst_PROGRAMS = example

example_SOURCES = examples/example.cpp
example_CPPFLAGS = -I$(top_srcdir)/include
