from numpy import *
#import math
from visit_utils import *
import glob
import os
#OpenDatabase("localhost:/Users/natarajan/Desktop/Research/NGA/AMR/amrex/Tutorials/Amr/FluidsSolverFrameworkRefCrit_bkp/Exec/SingleVortex/movie.visit", 0)
for i in range(1,2,1):
        solnfile="./plt%0.5d/Header"%i
        print solnfile
        #OpenDatabase("/Users/natarajan/Desktop/Research/NGA/AMR/amrex/Tutorials/Amr/FluidsSolverFrameworkRefCrit_bkp/Exec/SingleVortex/plt%0.5d/Header"%i, 0)
        OpenDatabase(solnfile, 0)
        filename='Bubble_AMReX%0.5d'%i
        print filename
	AddPlot("Pseudocolor", "phi11", 1, 1)
	DrawPlots()
	AddOperator("Slice", 1)
	DrawPlots()
	SliceAtts = SliceAttributes()
	SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
	SliceAtts.originPoint = (0, 0, 0)
	SliceAtts.originIntercept = 0
	SliceAtts.originPercent = 0
	SliceAtts.originZone = 0
	SliceAtts.originNode = 0
	SliceAtts.normal = (0, -1, 0)
	SliceAtts.axisType = SliceAtts.YAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
	SliceAtts.upAxis = (0, 0, 1)
	SliceAtts.project2d = 0
	SliceAtts.interactive = 1
	SliceAtts.flip = 0
	SliceAtts.originZoneDomain = 0
	SliceAtts.originNodeDomain = 0
	SliceAtts.meshName = "Mesh"
	SliceAtts.theta = 0
	SliceAtts.phi = 0
	SetOperatorOptions(SliceAtts, 1)
	# Begin spontaneous state
	View2DAtts = View2DAttributes()
	View2DAtts.windowCoords = (-2, 2, -2, 4)
	View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
	View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
	View2DAtts.fullFrameAutoThreshold = 100
	View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
	View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
	View2DAtts.windowValid = 1
	SetView2D(View2DAtts)
	# End spontaneous state
	
	# Begin spontaneous state
	View3DAtts = View3DAttributes()
	View3DAtts.viewNormal = (-0.24266, -0.938392, 0.246042)
	View3DAtts.focus = (0, 0, 1)
	View3DAtts.viewUp = (0.0934594, 0.229829, 0.968733)
	View3DAtts.viewAngle = 30
	View3DAtts.parallelScale = 4.12311
	View3DAtts.nearPlane = -8.24621
	View3DAtts.farPlane = 8.24621
	View3DAtts.imagePan = (0, 0)
	View3DAtts.imageZoom = 1
	View3DAtts.perspective = 1
	View3DAtts.eyeAngle = 2
	View3DAtts.centerOfRotationSet = 0
	View3DAtts.centerOfRotation = (0, 0, 1)
	View3DAtts.axis3DScaleFlag = 0
	View3DAtts.axis3DScales = (1, 1, 1)
	View3DAtts.shear = (0, 0, 1)
	View3DAtts.windowValid = 1
	SetView3D(View3DAtts)
	# End spontaneous state
	
	# Begin spontaneous state
	View3DAtts = View3DAttributes()
	View3DAtts.viewNormal = (-0.25299, -0.938574, 0.234682)
	View3DAtts.focus = (0, 0, 1)
	View3DAtts.viewUp = (0.0616979, 0.226427, 0.972072)
	View3DAtts.viewAngle = 30
	View3DAtts.parallelScale = 4.12311
	View3DAtts.nearPlane = -8.24621
	View3DAtts.farPlane = 8.24621
	View3DAtts.imagePan = (0, 0)
	View3DAtts.imageZoom = 1
	View3DAtts.perspective = 1
	View3DAtts.eyeAngle = 2
	View3DAtts.centerOfRotationSet = 0
	View3DAtts.centerOfRotation = (0, 0, 1)
	View3DAtts.axis3DScaleFlag = 0
	View3DAtts.axis3DScales = (1, 1, 1)
	View3DAtts.shear = (0, 0, 1)
	View3DAtts.windowValid = 1
	SetView3D(View3DAtts)
	# End spontaneous state
	
	AddPlot("Subset", "levels", 1, 0)
	DrawPlots()
	SubsetAtts = SubsetAttributes()
	SubsetAtts.colorType = SubsetAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
	SubsetAtts.colorTableName = "Default"
	SubsetAtts.invertColorTable = 0
	SubsetAtts.legendFlag = 1
	SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
	SubsetAtts.lineWidth = 0
	SubsetAtts.singleColor = (0, 0, 0, 255)
	SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
	SubsetAtts.SetMultiColor(1, (0, 255, 0, 255))
	SubsetAtts.subsetNames = ("1", "2")
	SubsetAtts.opacity = 1
	SubsetAtts.wireframe = 1
	SubsetAtts.drawInternal = 0
	SubsetAtts.smoothingLevel = 0
	SubsetAtts.pointSize = 0.05
	SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
	SubsetAtts.pointSizeVarEnabled = 0
	SubsetAtts.pointSizeVar = "default"
	SubsetAtts.pointSizePixels = 2
	SetPlotOptions(SubsetAtts)
	SubsetAtts = SubsetAttributes()
	SubsetAtts.colorType = SubsetAtts.ColorByMultipleColors  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
	SubsetAtts.colorTableName = "Default"
	SubsetAtts.invertColorTable = 0
	SubsetAtts.legendFlag = 1
	SubsetAtts.lineStyle = SubsetAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
	SubsetAtts.lineWidth = 2
	SubsetAtts.singleColor = (0, 0, 0, 255)
	SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
	SubsetAtts.SetMultiColor(1, (0, 255, 0, 255))
	SubsetAtts.subsetNames = ("1", "2")
	SubsetAtts.opacity = 1
	SubsetAtts.wireframe = 1
	SubsetAtts.drawInternal = 0
	SubsetAtts.smoothingLevel = 0
	SubsetAtts.pointSize = 0.05
	SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
	SubsetAtts.pointSizeVarEnabled = 0
	SubsetAtts.pointSizeVar = "default"
	SubsetAtts.pointSizePixels = 2
	SetPlotOptions(SubsetAtts)
	PseudocolorAtts = PseudocolorAttributes()
        PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
        PseudocolorAtts.skewFactor = 1
        PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, CurrentPlot
        PseudocolorAtts.minFlag = 0
        PseudocolorAtts.min = 0
        PseudocolorAtts.maxFlag = 0
        PseudocolorAtts.max = 1
        PseudocolorAtts.centering = PseudocolorAtts.Natural  # Natural, Nodal, Zonal
        PseudocolorAtts.colorTableName = "hot"
        PseudocolorAtts.invertColorTable = 0
        PseudocolorAtts.opacityType = PseudocolorAtts.Ramp  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
        PseudocolorAtts.opacityVariable = ""
        PseudocolorAtts.opacity = 1
        PseudocolorAtts.opacityVarMin = 0
        PseudocolorAtts.opacityVarMax = 1
        PseudocolorAtts.opacityVarMinFlag = 0
        PseudocolorAtts.opacityVarMaxFlag = 0
        PseudocolorAtts.pointSize = 0.05
        PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
        PseudocolorAtts.pointSizeVarEnabled = 0
        PseudocolorAtts.pointSizeVar = "default"
        PseudocolorAtts.pointSizePixels = 2
        PseudocolorAtts.lineStyle = PseudocolorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
        PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
        PseudocolorAtts.lineWidth = 0
        PseudocolorAtts.tubeResolution = 10
        PseudocolorAtts.tubeRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
        PseudocolorAtts.tubeRadiusAbsolute = 0.125
        PseudocolorAtts.tubeRadiusBBox = 0.005
        PseudocolorAtts.tubeRadiusVarEnabled = 0
        PseudocolorAtts.tubeRadiusVar = ""
        PseudocolorAtts.tubeRadiusVarRatio = 10
        PseudocolorAtts.tailStyle = PseudocolorAtts.None  # None, Spheres, Cones
        PseudocolorAtts.headStyle = PseudocolorAtts.None  # None, Spheres, Cones
        PseudocolorAtts.endPointRadiusSizeType = PseudocolorAtts.FractionOfBBox  # Absolute, FractionOfBBox
        PseudocolorAtts.endPointRadiusAbsolute = 0.125
        PseudocolorAtts.endPointRadiusBBox = 0.05
        PseudocolorAtts.endPointResolution = 10
        PseudocolorAtts.endPointRatio = 5
        PseudocolorAtts.endPointRadiusVarEnabled = 0
        PseudocolorAtts.endPointRadiusVar = ""
        PseudocolorAtts.endPointRadiusVarRatio = 10
        PseudocolorAtts.renderSurfaces = 1
        PseudocolorAtts.renderWireframe = 0
        PseudocolorAtts.renderPoints = 0
        PseudocolorAtts.smoothingLevel = 0
        PseudocolorAtts.legendFlag = 1
        PseudocolorAtts.lightingFlag = 1
        PseudocolorAtts.wireframeColor = (0, 0, 0, 0)
	PseudocolorAtts.pointColor = (0, 0, 0, 0)
        SetPlotOptions(PseudocolorAtts)
	SaveWindowAtts = SaveWindowAttributes()
	SaveWindowAtts.outputToCurrentDirectory = 0
	SaveWindowAtts.outputDirectory = "."
	SaveWindowAtts.fileName = filename
	SaveWindowAtts.family = 0
	SaveWindowAtts.format = SaveWindowAtts.JPEG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
	SaveWindowAtts.width = 1024
	SaveWindowAtts.height = 1024
	SaveWindowAtts.screenCapture = 0
	SaveWindowAtts.saveTiled = 0
	SaveWindowAtts.quality = 100
	SaveWindowAtts.progressive = 0
	SaveWindowAtts.binary = 0
	SaveWindowAtts.stereo = 0
	SaveWindowAtts.compression = SaveWindowAtts.None  # None, PackBits, Jpeg, Deflate
	SaveWindowAtts.forceMerge = 0
	SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
	SaveWindowAtts.advancedMultiWindowSave = 0
	SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win1.layer = 0
	SaveWindowAtts.subWindowAtts.win1.transparency = 0
	SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win2.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win2.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win2.layer = 0
	SaveWindowAtts.subWindowAtts.win2.transparency = 0
	SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win3.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win3.layer = 0
	SaveWindowAtts.subWindowAtts.win3.transparency = 0
	SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win4.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win4.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win4.layer = 0
	SaveWindowAtts.subWindowAtts.win4.transparency = 0
	SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win5.layer = 0
	SaveWindowAtts.subWindowAtts.win5.transparency = 0
	SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win6.layer = 0
	SaveWindowAtts.subWindowAtts.win6.transparency = 0
	SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win7.layer = 0
	SaveWindowAtts.subWindowAtts.win7.transparency = 0
	SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win8.layer = 0
	SaveWindowAtts.subWindowAtts.win8.transparency = 0
	SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win9.layer = 0
	SaveWindowAtts.subWindowAtts.win9.transparency = 0
	SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win10.layer = 0
	SaveWindowAtts.subWindowAtts.win10.transparency = 0
	SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win11.layer = 0
	SaveWindowAtts.subWindowAtts.win11.transparency = 0
	SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win12.layer = 0
	SaveWindowAtts.subWindowAtts.win12.transparency = 0
	SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win13.layer = 0
	SaveWindowAtts.subWindowAtts.win13.transparency = 0
	SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win14.layer = 0
	SaveWindowAtts.subWindowAtts.win14.transparency = 0
	SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win15.layer = 0
	SaveWindowAtts.subWindowAtts.win15.transparency = 0
	SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
	SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
	SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
	SaveWindowAtts.subWindowAtts.win16.layer = 0
	SaveWindowAtts.subWindowAtts.win16.transparency = 0
	SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()
	DeleteActivePlots()
	DeleteActivePlots()
        OpenDatabase(solnfile)
