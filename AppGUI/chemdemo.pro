######################################################################
# Automatically generated by qmake (2.01a) Wed Jan 7 16:01:24 2009
######################################################################

CONFIG += qt debug
CONFIG += hdf5 
QT += opengl
QT += xml
TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += .
INCLUDEPATH += /usr/local/hdf5/x86_64/include
LIBPATH += /usr/local/hdf5/lib64
LIBS += /usr/local/hdf5/x86_64/lib/libhdf5.a

# Input
HEADERS += bdcndwindow.h \
           cutplane.h \
           getfile.h \
           glviewer.h \
           mgviewer.h \
           vmergewindow.h \
           grid.h \
           initcndwindow.h \
           qualitydialog.h \
           mainwindow.h \
           pages.h \
           solverwindow.h \
           cutdialog.h \
           importwindow.h \
           physicswindow.h \
           mdlwindow.h \
           fvmadapt.h \
           refdialog.h \
           progressdialog.h \
           stateregion.h \
           varwindow.h \
           pbwindow.h \
           pboperation.h \
           defines.h \
           helpwindow.h
          
SOURCES += bdcndwindow.cpp \
           cutplane.cpp \
           getfile.cpp \
           glviewer.cpp \
           mgviewer.cpp \
           vmergewindow.cpp \
           initcndwindow.cpp \
           qualitydialog.cpp \
           main.cpp \
           mainwindow.cpp \
           make_objects.cpp \
           pages.cpp \
           read_grid.cc \
           solverwindow.cpp \
           cutdialog.cpp \
           importwindow.cpp \
           physicswindow.cpp \
           mdlwindow.cpp \
           fvmadapt.cpp \
           mark_node.cpp \
           refdialog.cpp \
           progressdialog.cpp \
           stateregion.cpp \
           varwindow.cpp \
           pbwindow.cpp \
           pboperation.cpp \
           helpwindow.cpp

RESOURCES     = chemdemo.qrc
#install
target.path = $$[CHEMDEMOPATH]/
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS chemdemo.pro png xml
sources.path = $$[CHEMDEMOPATH]/
INSTALLS += target sources
