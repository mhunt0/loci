/****************************************************************************
** Meta object code from reading C++ file 'qualitydialog.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qualitydialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'qualitydialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_QualityDialog[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_QualityDialog[] = {
    "QualityDialog\0\0buttonToggled(int)\0"
};

const QMetaObject QualityDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_QualityDialog,
      qt_meta_data_QualityDialog, 0 }
};

const QMetaObject *QualityDialog::metaObject() const
{
    return &staticMetaObject;
}

void *QualityDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_QualityDialog))
        return static_cast<void*>(const_cast< QualityDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int QualityDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: buttonToggled((*reinterpret_cast< int(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
