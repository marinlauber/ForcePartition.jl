# state file generated using paraview version 5.13.3
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1280, 1280]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [33.5, 25.5, 33.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [137.45008547111246, 152.7882546344924, -28.32981985171648]
renderView1.CameraFocalPoint = [19.05338994563121, 8.913358664835044, 47.6880735109266]
renderView1.CameraViewUp = [-0.21583327224792098, -0.31136277146706526, -0.925456224321336]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 52.084066661504075
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1280, 1280)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
dickinsonpvd = PVDReader(registrationName='Dickinson.pvd', FileName='/home/marin/Dickinson.pvd')
dickinsonpvd.PointArrays = ['Q', 'u', 'ϕ1', 'ϕ2', 'λ2', 'p', 'd']

# create a new 'Extract Subset'
extractSubset1 = ExtractSubset(registrationName='ExtractSubset1', Input=dickinsonpvd)
extractSubset1.VOI = [1, 64, 1, 48, 1, 64]

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=extractSubset1)
contour1.ContourBy = ['POINTS', 'd']
contour1.Isosurfaces = [0.0]
contour1.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=extractSubset1)
contour2.ContourBy = ['POINTS', 'ϕ1']
contour2.Isosurfaces = [-5.0, -3.888888888888889, -2.7777777777777777, -1.6666666666666665, -0.5555555555555554, 0.5555555555555554, 1.666666666666667, 2.7777777777777786, 3.8888888888888893, 5.0]
contour2.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour3 = Contour(registrationName='Contour3', Input=extractSubset1)
contour3.ContourBy = ['POINTS', 'ϕ2']
contour3.Isosurfaces = [-5.0, -3.888888888888889, -2.7777777777777777, -1.6666666666666665, -0.5555555555555554, 0.5555555555555554, 1.666666666666667, 2.7777777777777786, 3.8888888888888893, 5.0]
contour3.PointMergeMethod = 'Uniform Binning'

# create a new 'Contour'
contour4 = Contour(registrationName='Contour4', Input=extractSubset1)
contour4.ContourBy = ['POINTS', 'Q']
contour4.Isosurfaces = [0.1]
contour4.PointMergeMethod = 'Uniform Binning'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from extractSubset1
extractSubset1Display = Show(extractSubset1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
extractSubset1Display.Representation = 'Outline'
extractSubset1Display.ColorArrayName = [None, '']
extractSubset1Display.SelectNormalArray = 'None'
extractSubset1Display.SelectTangentArray = 'None'
extractSubset1Display.SelectTCoordArray = 'None'
extractSubset1Display.TextureTransform = 'Transform2'
extractSubset1Display.OSPRayScaleArray = 'Q'
extractSubset1Display.OSPRayScaleFunction = 'Piecewise Function'
extractSubset1Display.Assembly = ''
extractSubset1Display.SelectedBlockSelectors = ['']
extractSubset1Display.SelectOrientationVectors = 'None'
extractSubset1Display.ScaleFactor = 6.300000000000001
extractSubset1Display.SelectScaleArray = 'None'
extractSubset1Display.GlyphType = 'Arrow'
extractSubset1Display.GlyphTableIndexArray = 'None'
extractSubset1Display.GaussianRadius = 0.315
extractSubset1Display.SetScaleArray = ['POINTS', 'Q']
extractSubset1Display.ScaleTransferFunction = 'Piecewise Function'
extractSubset1Display.OpacityArray = ['POINTS', 'Q']
extractSubset1Display.OpacityTransferFunction = 'Piecewise Function'
extractSubset1Display.DataAxesGrid = 'Grid Axes Representation'
extractSubset1Display.PolarAxes = 'Polar Axes Representation'
extractSubset1Display.ScalarOpacityUnitDistance = 1.762960213463563
extractSubset1Display.OpacityArrayName = ['POINTS', 'Q']
extractSubset1Display.ColorArray2Name = ['POINTS', 'Q']
extractSubset1Display.SliceFunction = 'Plane'
extractSubset1Display.Slice = 31
extractSubset1Display.SelectInputVectors = ['POINTS', 'u']
extractSubset1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
extractSubset1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.5202293090518197, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
extractSubset1Display.ScaleTransferFunction.Points = [-2.381096839904785, 0.0, 0.5, 0.0, 1.9338678121566772, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
extractSubset1Display.OpacityTransferFunction.Points = [-2.381096839904785, 0.0, 0.5, 0.0, 1.9338678121566772, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
extractSubset1Display.SliceFunction.Origin = [33.5, 25.5, 33.5]

# show data from contour1
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.AmbientColor = [1.0, 0.0, 0.4980392156862745]
contour1Display.ColorArrayName = ['POINTS', '']
contour1Display.DiffuseColor = [1.0, 0.0, 0.4980392156862745]
contour1Display.Specular = 1.0
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.SelectTCoordArray = 'None'
contour1Display.TextureTransform = 'Transform2'
contour1Display.OSPRayScaleArray = 'd'
contour1Display.OSPRayScaleFunction = 'Piecewise Function'
contour1Display.Assembly = ''
contour1Display.SelectedBlockSelectors = ['']
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 1.828216552734375
contour1Display.SelectScaleArray = 'd'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'd'
contour1Display.GaussianRadius = 0.09141082763671875
contour1Display.SetScaleArray = ['POINTS', 'd']
contour1Display.ScaleTransferFunction = 'Piecewise Function'
contour1Display.OpacityArray = ['POINTS', 'd']
contour1Display.OpacityTransferFunction = 'Piecewise Function'
contour1Display.DataAxesGrid = 'Grid Axes Representation'
contour1Display.PolarAxes = 'Polar Axes Representation'
contour1Display.SelectInputVectors = ['POINTS', 'Normals']
contour1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
contour1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.5202293090518197, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from contour2
contour2Display = Show(contour2, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'a1'
a1TF2D = GetTransferFunction2D('a1')
a1TF2D.ScalarRangeInitialized = 1
a1TF2D.Range = [-5.0, 5.0, 0.0, 1.0]

# get color transfer function/color map for 'a1'
a1LUT = GetColorTransferFunction('a1')
a1LUT.TransferFunction2D = a1TF2D
a1LUT.RGBPoints = [-5.0, 0.498039, 0.231373, 0.031373, -4.37255, 0.62599, 0.30273, 0.026451, -3.7451, 0.746943, 0.387082, 0.037524, -3.117645, 0.85767, 0.490427, 0.071972, -2.490195, 0.936409, 0.617762, 0.236371, -1.8627450000000003, 0.992695, 0.743099, 0.43291, -1.2352950000000003, 0.995156, 0.841523, 0.63714, -0.6078450000000002, 0.985313, 0.913802, 0.813687, 0.01960784999999987, 0.966244, 0.966398, 0.967705, 0.6470600000000006, 0.889965, 0.89504, 0.938178, 1.2745099999999994, 0.806151, 0.804306, 0.894656, 1.9019600000000008, 0.712649, 0.688658, 0.833141, 2.5294100000000004, 0.594233, 0.554325, 0.744637, 3.156865, 0.474894, 0.404229, 0.652364, 3.7843149999999994, 0.366628, 0.217224, 0.563783, 4.411765000000001, 0.266436, 0.089965, 0.434833, 5.0, 0.176471, 0.0, 0.294118]
a1LUT.ColorSpace = 'Lab'
a1LUT.NanColor = [1.0, 0.0, 0.0]
a1LUT.NumberOfTableValues = 11
a1LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour2Display.Representation = 'Surface'
contour2Display.ColorArrayName = ['POINTS', 'ϕ1']
contour2Display.LookupTable = a1LUT
contour2Display.Specular = 1.0
contour2Display.SelectNormalArray = 'Normals'
contour2Display.SelectTangentArray = 'None'
contour2Display.SelectTCoordArray = 'None'
contour2Display.TextureTransform = 'Transform2'
contour2Display.OSPRayScaleArray = 'ϕ1'
contour2Display.OSPRayScaleFunction = 'Piecewise Function'
contour2Display.Assembly = ''
contour2Display.SelectedBlockSelectors = ['']
contour2Display.SelectOrientationVectors = 'None'
contour2Display.ScaleFactor = 1.8667794227600099
contour2Display.SelectScaleArray = 'ϕ1'
contour2Display.GlyphType = 'Arrow'
contour2Display.GlyphTableIndexArray = 'ϕ1'
contour2Display.GaussianRadius = 0.09333897113800049
contour2Display.SetScaleArray = ['POINTS', 'ϕ1']
contour2Display.ScaleTransferFunction = 'Piecewise Function'
contour2Display.OpacityArray = ['POINTS', 'ϕ1']
contour2Display.OpacityTransferFunction = 'Piecewise Function'
contour2Display.DataAxesGrid = 'Grid Axes Representation'
contour2Display.PolarAxes = 'Polar Axes Representation'
contour2Display.SelectInputVectors = ['POINTS', 'Normals']
contour2Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
contour2Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.5202293090518197, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour2Display.ScaleTransferFunction.Points = [-5.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour2Display.OpacityTransferFunction.Points = [-5.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]

# show data from contour3
contour3Display = Show(contour3, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'a2'
a2TF2D = GetTransferFunction2D('a2')
a2TF2D.ScalarRangeInitialized = 1
a2TF2D.Range = [-5.0, 5.0, 0.0, 1.0]

# get color transfer function/color map for 'a2'
a2LUT = GetColorTransferFunction('a2')
a2LUT.TransferFunction2D = a2TF2D
a2LUT.RGBPoints = [-5.0, 0.152941, 0.392157, 0.098039, -4.3725499999999995, 0.246444, 0.505344, 0.117724, -3.7451, 0.351942, 0.614533, 0.161399, -3.1176450000000004, 0.474971, 0.717878, 0.240138, -2.490195, 0.611995, 0.811226, 0.392849, -1.8627450000000008, 0.746328, 0.893118, 0.565321, -1.2352950000000003, 0.859516, 0.94233, 0.747405, -0.6078450000000002, 0.928105, 0.96386, 0.875663, 0.01960784999999987, 0.969089, 0.966859, 0.968012, 0.6470600000000006, 0.983852, 0.910265, 0.948328, 1.2745099999999985, 0.979239, 0.833218, 0.914648, 1.9019600000000008, 0.949712, 0.729873, 0.862976, 2.5294100000000004, 0.905652, 0.58293, 0.763552, 3.156865, 0.85521, 0.410073, 0.652211, 3.7843149999999994, 0.793695, 0.183699, 0.531642, 4.411765000000003, 0.683737, 0.063899, 0.420761, 5.0, 0.556863, 0.003922, 0.321569]
a2LUT.ColorSpace = 'Lab'
a2LUT.NanColor = [1.0, 0.0, 0.0]
a2LUT.NumberOfTableValues = 11
a2LUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
contour3Display.Representation = 'Surface'
contour3Display.ColorArrayName = ['POINTS', 'ϕ2']
contour3Display.LookupTable = a2LUT
contour3Display.Specular = 1.0
contour3Display.SelectNormalArray = 'Normals'
contour3Display.SelectTangentArray = 'None'
contour3Display.SelectTCoordArray = 'None'
contour3Display.TextureTransform = 'Transform2'
contour3Display.OSPRayScaleArray = 'ϕ2'
contour3Display.OSPRayScaleFunction = 'Piecewise Function'
contour3Display.Assembly = ''
contour3Display.SelectedBlockSelectors = ['']
contour3Display.SelectOrientationVectors = 'None'
contour3Display.ScaleFactor = 2.436218643188477
contour3Display.SelectScaleArray = 'ϕ2'
contour3Display.GlyphType = 'Arrow'
contour3Display.GlyphTableIndexArray = 'ϕ2'
contour3Display.GaussianRadius = 0.12181093215942383
contour3Display.SetScaleArray = ['POINTS', 'ϕ2']
contour3Display.ScaleTransferFunction = 'Piecewise Function'
contour3Display.OpacityArray = ['POINTS', 'ϕ2']
contour3Display.OpacityTransferFunction = 'Piecewise Function'
contour3Display.DataAxesGrid = 'Grid Axes Representation'
contour3Display.PolarAxes = 'Polar Axes Representation'
contour3Display.SelectInputVectors = ['POINTS', 'Normals']
contour3Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
contour3Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.5202293090518197, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour3Display.ScaleTransferFunction.Points = [-0.22361451387405396, 0.0, 0.5, 0.0, -0.22358399629592896, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour3Display.OpacityTransferFunction.Points = [-0.22361451387405396, 0.0, 0.5, 0.0, -0.22358399629592896, 1.0, 0.5, 0.0]

# show data from contour4
contour4Display = Show(contour4, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
contour4Display.Representation = 'Surface'
contour4Display.AmbientColor = [0.0, 0.3333333333333333, 0.4980392156862745]
contour4Display.ColorArrayName = ['POINTS', '']
contour4Display.DiffuseColor = [0.0, 0.3333333333333333, 0.4980392156862745]
contour4Display.Opacity = 0.5
contour4Display.Specular = 1.0
contour4Display.SelectNormalArray = 'Normals'
contour4Display.SelectTangentArray = 'None'
contour4Display.SelectTCoordArray = 'None'
contour4Display.TextureTransform = 'Transform2'
contour4Display.OSPRayScaleArray = 'Q'
contour4Display.OSPRayScaleFunction = 'Piecewise Function'
contour4Display.Assembly = ''
contour4Display.SelectedBlockSelectors = ['']
contour4Display.SelectOrientationVectors = 'None'
contour4Display.ScaleFactor = 5.519553637504578
contour4Display.SelectScaleArray = 'Q'
contour4Display.GlyphType = 'Arrow'
contour4Display.GlyphTableIndexArray = 'Q'
contour4Display.GaussianRadius = 0.27597768187522886
contour4Display.SetScaleArray = ['POINTS', 'Q']
contour4Display.ScaleTransferFunction = 'Piecewise Function'
contour4Display.OpacityArray = ['POINTS', 'Q']
contour4Display.OpacityTransferFunction = 'Piecewise Function'
contour4Display.DataAxesGrid = 'Grid Axes Representation'
contour4Display.PolarAxes = 'Polar Axes Representation'
contour4Display.SelectInputVectors = ['POINTS', 'Normals']
contour4Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
contour4Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.5202293090518197, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
contour4Display.ScaleTransferFunction.Points = [0.10000000149011612, 0.0, 0.5, 0.0, 0.10001526027917862, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
contour4Display.OpacityTransferFunction.Points = [0.10000000149011612, 0.0, 0.5, 0.0, 0.10001526027917862, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for a1LUT in view renderView1
a1LUTColorBar = GetScalarBar(a1LUT, renderView1)
a1LUTColorBar.Orientation = 'Horizontal'
a1LUTColorBar.WindowLocation = 'Any Location'
a1LUTColorBar.Position = [0.2129533550396372, 0.023285899094437255]
a1LUTColorBar.Title = 'Potential Wing 1'
a1LUTColorBar.ComponentTitle = ''
a1LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
a1LUTColorBar.TitleFontSize = 14
a1LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
a1LUTColorBar.LabelFontSize = 14
a1LUTColorBar.ScalarBarLength = 0.20812500000000067
a1LUTColorBar.DrawScalarBarOutline = 1
a1LUTColorBar.ScalarBarOutlineColor = [0.0, 0.0, 0.0]
a1LUTColorBar.RangeLabelFormat = '%-#6.1f'

# set color bar visibility
a1LUTColorBar.Visibility = 1

# get color legend/bar for a2LUT in view renderView1
a2LUTColorBar = GetScalarBar(a2LUT, renderView1)
a2LUTColorBar.Orientation = 'Horizontal'
a2LUTColorBar.WindowLocation = 'Any Location'
a2LUTColorBar.Position = [0.7157846121177803, 0.020698576972833106]
a2LUTColorBar.Title = 'Potential Wing 2'
a2LUTColorBar.ComponentTitle = ''
a2LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
a2LUTColorBar.TitleFontSize = 14
a2LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
a2LUTColorBar.LabelFontSize = 14
a2LUTColorBar.ScalarBarLength = 0.20812500000000056
a2LUTColorBar.DrawScalarBarOutline = 1
a2LUTColorBar.ScalarBarOutlineColor = [0.0, 0.0, 0.0]
a2LUTColorBar.RangeLabelFormat = '%-#6.1f'

# set color bar visibility
a2LUTColorBar.Visibility = 1

# show color legend
contour2Display.SetScalarBarVisibility(renderView1, True)

# show color legend
contour3Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'a2'
a2PWF = GetOpacityTransferFunction('a2')
a2PWF.Points = [-5.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
a2PWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'a1'
a1PWF = GetOpacityTransferFunction('a1')
a1PWF.Points = [-5.0, 0.0, 0.5, 0.0, 5.0, 1.0, 0.5, 0.0]
a1PWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 12.0001
animationScene1.EndTime = 12.0001
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(contour2)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------