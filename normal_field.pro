TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11

SOURCES += main.cpp \
    normalfield.cpp \
    grid2D.cpp \
    slae.cpp \
    matrixfem.cpp

HEADERS += \
    normalfield.h \
    grid2D.h \
    slae.h \
    matrixfem.h \
    matrix.h

OTHER_FILES += \
    area.txt \
    grid.txt \
    ku1.txt \
    sigma.txt

