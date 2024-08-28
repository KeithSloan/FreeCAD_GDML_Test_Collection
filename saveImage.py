import FreeCAD
import os
import sys
import glob
import Mesh

def loadAndSaveImage(doc):
    FreeCAD.openDocument(doc)
    Gui.SendMsgToActiveView("ViewFit")
    Gui.activeDocument().activeView().viewIsometric()
    docName = FreeCAD.ActiveDocument.FileName
    name = docName.replace('.FCStd','')
    Gui.activeDocument().activeView().saveImage('{}.png'.format(name),1455,833,'Current')
    FreeCAD.closeDocument(FreeCAD.ActiveDocument.Name)

def loadStlAndSaveImage(doc):
    Mesh.open(doc)
    Gui.SendMsgToActiveView("ViewFit")
    Gui.activeDocument().activeView().viewIsometric()
    docName = FreeCAD.ActiveDocument.FileName
    name = doc.replace('.stl','')
    Gui.activeDocument().activeView().saveImage('{}.png'.format(name),1455,833,'Current')
    FreeCAD.closeDocument(FreeCAD.ActiveDocument.Name)

"""
Usage:
1 - Open FreeCAD
2 - Open python console (View->Panels->Python cosole)
3 - Copy above loadAndSaveImage methos and paste into python console
4 - Change to the desired directory:
   os.chdir('/home/Documents/.....')
5 - get list of FCStd files in current directory:
import glob
docs = glob.glob("./*.FCStd")
# 6 - loop over files:
for doc in docs:
   print(doc)
   loadAndSaveImage(doc)
# or
#  loadStlAndSaveImage(doc)
"""

#os.chdir("needed directory")
#glob.glob("./*.FCStd")
#for doc in docs
#   print(doc)
#  loadAndSaveImage(doc)
