TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    integrator.cpp \
    LA.cpp \
    model.cpp \
    TQuaternion.cpp \
    custom.cpp \
    aes_param.cpp

HEADERS += \
    integrator.h \
    LA.h \
    model.h \
    TQuaternion.h \
    custom.h \
    aes_param.h
