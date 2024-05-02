import logging
import os

import vtk

import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin


#
# FMM
#

class FMM(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "FractionalMyocardialMass"
        self.parent.categories = ["Cardiac"]
        self.parent.dependencies = []  
        self.parent.contributors = ["David Molony (NGHS)"]  #
        self.parent.helpText = """
This module computes the fractional myocardial mass for an input myocardial volume mesh.
"""
        # TODO: replace with organization, grant and thanks
        self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc., Andras Lasso, PerkLab,
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
"""

#
# FMMWidget
#

class FMMWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Uses ScriptedLoadableModuleWidget base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent=None):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)  # needed for parameter node observation
        self.logic = None
        self._parameterNode = None
        self._updatingGUIFromParameterNode = False

    def setup(self):
        """
        Called when the user opens the module the first time and the widget is initialized.
        """
        ScriptedLoadableModuleWidget.setup(self)

        # Load widget from .ui file (created by Qt Designer).
        # Additional widgets can be instantiated manually and added to self.layout.
        uiWidget = slicer.util.loadUI(self.resourcePath('UI/FMM.ui'))
        self.layout.addWidget(uiWidget)
        self.ui = slicer.util.childWidgetVariables(uiWidget)

        # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
        # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
        # "setMRMLScene(vtkMRMLScene*)" slot.
        uiWidget.setMRMLScene(slicer.mrmlScene)

        # Create logic class. Logic implements all computations that should be possible to run
        # in batch mode, without a graphical user interface.
        self.logic = FMMLogic()

        # Connections

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
        # (in the selected parameter node).
        self.ui.inputModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.inputLCAMarkupSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.inputRCAMarkupSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.inputMMARMarkupSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.outputMeshSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
        self.ui.outputTableSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)

        # Buttons
        self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

    def cleanup(self):
        """
        Called when the application closes and the module widget is destroyed.
        """
        self.removeObservers()

    def enter(self):
        """
        Called each time the user opens this module.
        """
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self):
        """
        Called each time the user opens a different module.
        """
        # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
        self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    def onSceneStartClose(self, caller, event):
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event):
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self):
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.

        self.setParameterNode(self.logic.getParameterNode())

        # Select default input nodes if nothing is selected yet to save a few clicks for the user
        if not self._parameterNode.GetNodeReference("InputModel"):
            firstModelNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLModelNode")
            if firstModelNode:
                self._parameterNode.SetNodeReferenceID("InputModel", firstModelNode.GetID())

    def setParameterNode(self, inputParameterNode):
        """
        Set and observe parameter node.
        Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
        """

        if inputParameterNode:
            self.logic.setDefaultParameters(inputParameterNode)

        # Unobserve previously selected parameter node and add an observer to the newly selected.
        # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
        # those are reflected immediately in the GUI.
        if self._parameterNode is not None:
            self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
        self._parameterNode = inputParameterNode
        if self._parameterNode is not None:
            self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

        # Initial GUI update
        self.updateGUIFromParameterNode()

    def updateGUIFromParameterNode(self, caller=None, event=None):
        """
        This method is called whenever parameter node is changed.
        The module GUI is updated to show the current state of the parameter node.
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
        self._updatingGUIFromParameterNode = True

        # Update node selectors and sliders
        self.ui.inputModelSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputModel"))
        self.ui.inputLCAMarkupSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputLCAMarkup"))
        self.ui.inputRCAMarkupSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputRCAMarkup"))
        self.ui.inputMMARMarkupSelector.setCurrentNode(self._parameterNode.GetNodeReference("InputMMARMarkup"))
        self.ui.outputTableSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputTable"))
        self.ui.outputMeshSelector.setCurrentNode(self._parameterNode.GetNodeReference("OutputMesh"))

        # Update buttons states and tooltips
        if self._parameterNode.GetNodeReference("InputModel") and self._parameterNode.GetNodeReference("InputLCAMarkup") and self._parameterNode.GetNodeReference("OutputTable") and self._parameterNode.GetNodeReference("OutputMesh"):
            self.ui.applyButton.toolTip = "Compute output myocardial mass"
            self.ui.applyButton.enabled = True
        else:
            self.ui.applyButton.toolTip = "Select input and output model and mesh nodes"
            self.ui.applyButton.enabled = False

        # All the GUI updates are done
        self._updatingGUIFromParameterNode = False

    def updateParameterNodeFromGUI(self, caller=None, event=None):
        """
        This method is called when the user makes any change in the GUI.
        The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
        """

        if self._parameterNode is None or self._updatingGUIFromParameterNode:
            return

        wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

        self._parameterNode.SetNodeReferenceID("InputModel", self.ui.inputModelSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("InputLCAMarkup", self.ui.inputLCAMarkupSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("InputRCAMarkup", self.ui.inputRCAMarkupSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("InputMMARMarkup", self.ui.inputMMARMarkupSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("OutputMesh", self.ui.outputMeshSelector.currentNodeID)
        self._parameterNode.SetNodeReferenceID("OutputTable", self.ui.outputTableSelector.currentNodeID)

        self._parameterNode.EndModify(wasModified)

    def onApplyButton(self):
        """
        Run processing when user clicks "Apply" button.
        """
        with slicer.util.tryWithErrorDisplay("Failed to compute results.", waitCursor=True):
            inputModelNode = self._parameterNode.GetNodeReference("InputModel")
            outputModelNode = self._parameterNode.GetNodeReference("OutputMesh")
            outputTableNode = self._parameterNode.GetNodeReference("OutputTable")
            slicer.util.showStatusMessage("Segmenting myocardium...")
            slicer.app.processEvents()  # force update
        
            # Compute output
            outputTable, outputMesh = self.logic.segmentMesh(self.ui.inputModelSelector.currentNode(), self.ui.inputLCAMarkupSelector.currentNode(),
                               self.ui.inputRCAMarkupSelector.currentNode(), self.ui.inputMMARMarkupSelector.currentNode(), self.ui.outputTableSelector.currentNode(), self.ui.outputMeshSelector.currentNode())

        outputModelNode.SetAndObserveMesh(outputMesh)
        if not outputModelNode.GetDisplayNode():
            outputModelNode.CreateDefaultDisplayNodes()
            outputModelNode.GetDisplayNode().SetColor(0.0, 1.0, 0.0)
            outputModelNode.GetDisplayNode().SetLineWidth(3)
            outputModelNode.GetDisplayNode().SetScalarVisibility(1)
            inputModelNode.GetDisplayNode().SetVisibility(0)
            if self.ui.inputMMARMarkupSelector.currentNode() is not None:
                outputModelNode.GetDisplayNode().SetActiveScalarName('MMAR')
                MMARNode = slicer.util.getNode('MMAR')
                slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayout3DTableView)
                layoutManager = slicer.app.layoutManager()
                #layoutManager.tableWidget(0).setMRMLTableViewNode()
                slicer.app.applicationLogic().GetSelectionNode().SetReferenceActiveTableID(MMARNode.GetID())
                slicer.app.applicationLogic().PropagateTableSelection()
            else:
                outputModelNode.GetDisplayNode().SetActiveScalarName('Ids')
                slicer.app.layoutManager().setLayout(slicer.vtkMRMLLayoutNode.SlicerLayout3DTableView)
                slicer.app.applicationLogic().GetSelectionNode().SetReferenceActiveTableID(outputTableNode.GetID())
                slicer.app.applicationLogic().PropagateTableSelection()        

# FMMLogic
#

class FMMLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)

    def setDefaultParameters(self, parameterNode):
        """
        Initialize parameter node with default settings.
        """

    def markupToPolyData(self, markups, groupId=0):
        """
        Convert markup data to polydata and assign groupids for every point
        """
        points = vtk.vtkPoints()
        group = vtk.vtkIntArray()
        group.SetName('Ids')
        groupDict = {}
        lengthDict = {}
        childDict = {}

        visited = []
        stack = []

        children = vtk.vtkIdList()
        shNode = slicer.mrmlScene.GetSubjectHierarchyNode()
        itemId = shNode.GetItemByDataNode(markups)
        shNode.GetItemChildren(itemId, children)
        stack.append(itemId)

        while stack:
            currentId = stack.pop()
            if currentId not in visited:
                visited.append(currentId)
                markup = shNode.GetItemDataNode(currentId)
                for j in range(markup.GetNumberOfControlPoints()):
                    point = list(markup.GetNthControlPointPositionWorld(j))
                    points.InsertNextPoint(point[0], point[1], point[2])
                    group.InsertNextTuple([groupId,])

                # Get current nodes children
                shNode.GetItemChildren(currentId,children)
                itemIds = [children.GetId(i) for i in range(children.GetNumberOfIds())]
                stack.extend(itemIds)
                groupDict[markup.GetName()] = groupId
                lengthDict[markup.GetName()] = markup.GetCurveLengthWorld()
                childDict[groupId] = [shNode.GetItemDataNode(itemId).GetName() for itemId in itemIds]
                groupId += 1

        # convert child dictionary names to ids
        for groupId, names in childDict.items():
           childIds = [int(groupDict[name]) for name in names] 
           childDict[groupId] = childIds

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.GetPointData().AddArray(group)
        return polydata, groupDict, lengthDict, childDict
        
    def findMMARCenterline(self, centerlinePolydata, MMARPoint):
        """
        Identify pointid and corresponding groupid for the MMAR fiducial point
        """
        point = vtk.vtkPoints()
        point.InsertNextPoint(MMARPoint[0], MMARPoint[1], MMARPoint[2])
        MMARPolydata = vtk.vtkPolyData()
        MMARPolydata.SetPoints(point)
        locator =  vtk.vtkPointLocator()
        locator.SetDataSet(centerlinePolydata)
        locator.BuildLocator()
        pointId = locator.FindClosestPoint(MMARPoint)
        groupId = centerlinePolydata.GetPointData().GetArray(0).GetTuple(pointId)[0]
        return pointId, groupId
        
    def addMMARArray(self, polydata, childDict, pointId, groupId):
        """
        Add pointdata array indicating whether the vertex is "at risk"
        """

        stack = [groupId]
        childIds = []
        while stack:
            currentId = stack.pop()
            if currentId not in childIds:
                childIds.append(currentId)
                stack.extend(childDict[currentId])

        mmar = vtk.vtkIntArray()
        mmar.SetName("MMAR")

        for i in range(polydata.GetNumberOfPoints()):
            currentGroupId = int(polydata.GetPointData().GetArray(0).GetTuple(i)[0])
            #print(currentGroupId)
            if currentGroupId in childIds:
                if currentGroupId == groupId:
                    if i >= pointId:
                        mmar.InsertNextTuple([1,])
                    else:
                        mmar.InsertNextTuple([0,])
                else:   
                    mmar.InsertNextTuple([1,])
            else:
                mmar.InsertNextTuple([0,])
        
        polydata.GetPointData().AddArray(mmar)

        return polydata
        
    def voronoi(self, inputModel, polydata):
        """
        Assign pointdata to each mesh point using voronoi algorithm
        """
        kernel = vtk.vtkVoronoiKernel()
        
        nullValue = 0.0
        probeFilter = vtk.vtkPointInterpolator()
        probeFilter.SetInputData(inputModel.GetMesh())
        probeFilter.SetSourceData(polydata)
        probeFilter.SetKernel(kernel)
        probeFilter.SetNullPointsStrategyToClosestPoint()
        #probeFilter.SetNullPointsStrategyToNullValue()
        probeFilter.SetNullValue(nullValue)
        probeFilter.Update()
        outputMesh = probeFilter.GetOutput()
        return outputMesh    

    def computeVolume(self, outputMesh, ids, arrayIdx):        
        """Compute the volume associated to each id"""
        from statistics import mode

        vol = {int(idx):0 for idx in ids}

        for i in range(outputMesh.GetNumberOfCells()):
            tetra = outputMesh.GetCell(i)
            pts=tetra.GetPoints()
            cell_vol = tetra.ComputeVolume(pts.GetPoint(0), pts.GetPoint(1), pts.GetPoint(2), pts.GetPoint(3))
            pointIds = [tetra.GetPointIds().GetId(j) for j in range(4)]
            groupIds = [outputMesh.GetPointData().GetArray(arrayIdx).GetTuple(id)[0] for id in pointIds]
            groupId = int(mode(groupIds))
            vol[groupId] += cell_vol
            
        return vol
        
    def segmentMesh(self, inputModel, inputLCAMarkup, inputRCAMarkup, inputMMARMarkup, outputTable, outputModel):
        """
        Run the processing algorithm.
        Can be used without GUI widget.
        :param inputModel: input mesh
        :param inputLCAMarkup: input markup for LCA
        :param inputRCAMarkup: input markup for RCA
        :inputMMARMarkup: input markup for MMAR point
        :param outputTable: output table of myocardial masses for each centerline
        :param outputModel: output mesh of segmented myocardium
        """

        if not inputModel or not inputLCAMarkup or not outputTable or not outputModel:
            raise ValueError("Input model, markup, or output model is invalid")

        firstGroupId = 0
        polydata, groupName, lengthName, childDict = self.markupToPolyData(inputLCAMarkup, firstGroupId)
        if inputRCAMarkup is not None:
            firstGroupId = max(groupName.values()) + 1
            polydataRCA, groupNameRCA, lengthNameRCA, childDictRCA = self.markupToPolyData(inputRCAMarkup, firstGroupId)
            appendFilter = vtk.vtkAppendFilter()
            appendFilter.AddInputData(polydata)
            appendFilter.AddInputData(polydataRCA)
            appendFilter.Update()
            polydata = appendFilter.GetOutput()
            groupName.update(groupNameRCA)
            lengthName.update(lengthNameRCA)
            childDict.update(childDictRCA)

        if inputMMARMarkup is not None:
            MMARPoint = list(inputMMARMarkup.GetNthControlPointPositionWorld(0))
            pointId, groupId = self.findMMARCenterline(polydata, MMARPoint)
            polydata = self.addMMARArray(polydata, childDict, pointId, groupId)

        outputMesh = self.voronoi(inputModel, polydata)

       
        # TODO Allow polydata input containing Radius array and return minDiameterCol
        # TODO Allow surface model input and calculate volume/mass for each section
        
        num_arrays = outputMesh.GetPointData().GetNumberOfArrays() 
        groupIdArray = [i for i in range(num_arrays) if outputMesh.GetPointData().GetArrayName(i) == 'Ids'][0]

        ids = [outputMesh.GetPointData().GetArray(groupIdArray).GetTuple(i)[0] for i in range(outputMesh.GetPointData().GetArray(groupIdArray).GetNumberOfTuples())]
        ids = list(set(ids))
        vol = self.computeVolume(outputMesh, ids, groupIdArray)

        labelCol = vtk.vtkStringArray()
        groupCol = vtk.vtkIntArray()
        volCol = vtk.vtkDoubleArray()
        lengthCol = vtk.vtkDoubleArray()
        labelCol.SetName("Name")
        groupCol.SetName("Id")
        volCol.SetName("Volume (mm^3)")
        lengthCol.SetName("Segment length (mm)")
        for name, groupId in groupName.items():
            labelCol.InsertNextValue(name)
            groupCol.InsertNextTuple([groupId, ])
            lengthCol.InsertNextTuple([lengthName[name],])
            if groupId in vol:
                volCol.InsertNextTuple([vol[groupId], ])
            else:
                volCol.InsertNextTuple([0.0],)
      
        outputTable.AddColumn(labelCol)
        outputTable.AddColumn(groupCol)
        outputTable.AddColumn(volCol)
        outputTable.AddColumn(lengthCol)

        if inputMMARMarkup is not None:
            MMARArray = [i for i in range(num_arrays) if outputMesh.GetPointData().GetArrayName(i) == 'MMAR'][0]
            volMMAR = self.computeVolume(outputMesh, [0, 1], MMARArray)
            labelCol1 = vtk.vtkStringArray()
            labelCol1.InsertNextValue("Total volume (mm^3)")
            labelCol1.InsertNextValue("Myocardial volume at risk (mm^3)")
            volCol1 = vtk.vtkIntArray()
            volCol1.InsertNextTuple([volMMAR[1] + volMMAR[0], ])
            volCol1.InsertNextTuple([volMMAR[1], ])      
            MMARTableNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", "MMAR")
            MMARTableNode.AddColumn(labelCol1)
            MMARTableNode.AddColumn(volCol1)
        
        return outputTable, outputMesh