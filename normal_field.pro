TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -std=c++11

SOURCES += main.cpp \
    normalfield.cpp \
    grid2D.cpp \
    slae.cpp \
    matrixfem.cpp \
    anomalousfield.cpp \
    grid3D.cpp \
    hel.cpp \
    mke_3d.cpp

HEADERS += \
    normalfield.h \
    grid2D.h \
    slae.h \
    matrixfem.h \
    matrix.h \
    anomalousfield.h \
    grid3D.h \
    hel.h \
    mke_3d.h \
    search.h

OTHER_FILES += \
    area.txt \
    grid.txt \
    sigma.txt \
    area3D.txt \
    grid3D.txt \
    sigma3D.txt

