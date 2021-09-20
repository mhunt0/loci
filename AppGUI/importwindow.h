#ifndef IMPORTWINDOW_H
#define IMPORTWINDOW_H

#include <QWidget>
#include <QStringList>
#include <QDomDocument>
#include <QModelIndex>
#include <QGroupBox>
#include "pages.h"
#include "getfile.h"
class QListWidget;
class QListWidgetItem;
class QStackedWidget;
class QCheckBox;
class QButtonGroup;

class VogOption: public QGroupBox{
  Q_OBJECT
  public:
  VogOption( QDomElement& myelem, QWidget *parent=0 );
  QString currentText();
public slots:
  void unitButtonClicked(int);
private:
  QStringList alist;
  QStringList tag;
  QList<QCheckBox*> objs;
  VarGBox* LrefBox;
  QButtonGroup* radioGroup;
};


class XdrOption: public QGroupBox{
  Q_OBJECT
  public:
  XdrOption(QDomElement& myelem, QWidget *parent=0 );
  XdrOption();
  QString currentText();
public slots:
  void update(int);
private:
  QStringList optionList;
  QStringList tag;
  QList<QCheckBox*> objs;
  QList<VarGBox*> options;
  QSignalMapper *signalMapper;
};





class ImportWindow : public GeneralGroup
{
  Q_OBJECT
  
  public:
  ImportWindow(QDomElement& theelem,  QWidget* parent = 0);
  
public slots:
  void changePage(int);
  void helpClicked();
  void convert();
  //  void check();
  void usageButtonClicked();
  void updateFileName(QString);
signals:

private:

  QListWidget *typesWidget;
 
 
  QString importFileName;
  QStringList gridTypes;
  QStringList toXdr;
  QStringList toVog;
  
  int currentRow; //current grid type
  FindFileWindow* getFileWindow;

  QStackedWidget *pagesWidget;
  VogOption* option;
  QList<XdrOption*> otherOptions;
  QPushButton *convertButton;
   
};

#endif




