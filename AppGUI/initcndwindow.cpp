#include <QtGui>
#include <QFile>
#include <QString>
#include <QMessageBox>
#include <QTreeWidgetItem>

#include "initcndwindow.h"
#include "pages.h"
#include "getfile.h"
#include "stateregion.h"

InitCndWindow::InitCndWindow(QDomElement& theelem, QWidget* parent): GeneralGroup(theelem, parent)
{
  typesWidget = new QButtonGroup(this);
  typesWidget->setExclusive(true);
  pagesWidget = new QStackedWidget;
  
  QGroupBox* buttonGroup = new QGroupBox(myelem.attribute("label"));
  QHBoxLayout* buttonLayout = new QHBoxLayout;
  QDomElement elem = myelem.firstChildElement();
  if(elem.isNull()){
    QMessageBox::warning(window(), ".xml",
                         myelem.tagName()+tr(" has no child")
                         );
    return;
  }
  int count = 0;
  int current = 0;
  if(myelem.hasAttribute("current"))current = myelem.attribute("current").toInt();
  
  for (; !elem.isNull(); elem = elem.nextSiblingElement(), count++) {   
      
    QRadioButton *button = new QRadioButton(elem.tagName());
    button->setToolTip(elem.attribute("whatsThis"));
    button->setStatusTip(elem.attribute("whatsThis"));
    button->setWhatsThis(elem.attribute("whatsThis"));
    if(count==current)button->setChecked(true);
    buttonLayout->addWidget(button);
    typesWidget->addButton(button, count);
   
    if(elem.hasAttribute("action") && elem.attribute("action")=="find file"){
      fileSelected="./output/";
      FindFileWindow* getFileWindow = new FindFileWindow(elem, fileSelected);
      connect(this, SIGNAL(directoryChanged(QString)), getFileWindow, SLOT(addDirectory(QString)));
      connect(getFileWindow, SIGNAL(fileNameSelected(QString)), this, SLOT(checkStatus()));
      pagesWidget->addWidget(getFileWindow);
     
    }else if(elem.hasAttribute("element")&&elem.attribute("element")=="panel"){
      VarPanel* bdCndPage = new VarPanel(elem);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), bdCndPage, SLOT(updateShowStatus(const bool&)));
      pagesWidget->addWidget(bdCndPage);
    }else if(elem.hasAttribute("element")&&elem.attribute("element")=="page"){
      VarPage* bdCndPage = new VarPage(elem);
      connect(bdCndPage, SIGNAL(textChanged(const QString&)), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), bdCndPage, SLOT(updateShowStatus(const bool&)));
      pagesWidget->addWidget(bdCndPage);
    }else if(elem.hasAttribute("element")&&elem.attribute("element")=="regionWindow"){
      RegionWindow* bdCndPage = new RegionWindow(elem);
      connect(bdCndPage, SIGNAL(textChanged()), this, SLOT(checkStatus()));
      connect(this, SIGNAL(stateChanged()), bdCndPage, SLOT(changeState()));
      connect(this, SIGNAL(componentsChanged()), bdCndPage, SIGNAL(componentsChanged()));
      connect(this, SIGNAL(showStatus(const bool&)), bdCndPage, SLOT(updateShowStatus(const bool&)));
      connect(bdCndPage, SIGNAL( valueChanged(const QTreeWidgetItem*)),
              this, SIGNAL( valueChanged(const QTreeWidgetItem*)));
      pagesWidget->addWidget(bdCndPage);
    }else{
      QMessageBox::warning(window(), elem.tagName(),
                           tr(" don't know how to handle it yet: ")
                           +elem.attribute("element")
                           );
      
    }
    
  }

  connect(typesWidget,
          SIGNAL(buttonClicked(int)),
          this, SLOT(changePage(int)));
  changePage(current);
  myelem.setAttribute("current", current);
    
  buttonGroup->setLayout(buttonLayout);    
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(buttonGroup);
  mainLayout->addWidget(pagesWidget);
  mainLayout->addStretch(1);
  setLayout(mainLayout);
  setWindowTitle(myelem.attribute("title"));
  checkStatus();
  
}

void  InitCndWindow::setDirectory(QString s){
  QString dir= s+"/output/";
  emit directoryChanged(dir);
}
void InitCndWindow::changePage(int id)
{
  
  myelem.setAttribute("current",id);
  QDomElement elt = myelem.firstChildElement();
  for(int i = 0; i < id; i++)elt = elt.nextSiblingElement();
  
  pagesWidget->setCurrentIndex(id);
  emit showShapes(elt.attribute("element")=="regionWindow");
  checkStatus();
}

void InitCndWindow::checkStatus(){
  int current= myelem.attribute("current").toInt();
  QDomElement elt = myelem.firstChildElement();
  for(int i=0; i < current; i++) elt = elt.nextSiblingElement();
  myelem.setAttribute("status", elt.attribute("status"));
  myelem.setAttribute("currentText", elt.attribute("currentText"));
  
  emit updateStatusTip(myelem.attribute("buttonIndex").toInt());
  
}









