TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    grid.cpp \
    normalfield.cpp

HEADERS += \
    grid.h \
    normalfield.h

OTHER_FILES += \
    area.txt \
    grid.txt

