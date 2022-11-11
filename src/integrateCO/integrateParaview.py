# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

slozka = 'ZZ_cases_flow_isoT/intraTrans_yInf_0.1_R_0.01_T_500_cS_0.0033_k0_100_tort_1_inv_0.1116/intraTrans_yInf_0.1_R_0.01_T_500_cS_0.0033_k0_100_tort_1_inv_0.1116.OpenFOAM'
data = 'ZZ_cases_flow_isoT/intraTrans_yInf_0.1_R_0.01_T_500_cS_0.0033_k0_100_tort_1_inv_0.1116/dataFromInt.csv'

# load state
LoadState('/multipede1/khyrm/reactingHetCatFoam/tutorials/verification/innerSpherePartTransp/integrateCO/integCO.pvsm', LoadStateDataFileOptions='Choose File Names',
    DataDirectory='/multipede1/khyrm/reactingHetCatFoam/tutorials/verification/innerSpherePartTransp/integrateCO',
    intraTrans_yInf_01_R_001_T_500_cS_00033_k0_100_tort_1_inv_0OpenFOAMFileName='/multipede1/khyrm/reactingHetCatFoam/tutorials/verification/innerSpherePartTransp/%s'%slozka)

# find view
spreadSheetView1 = FindViewOrCreate('SpreadSheetView1', viewtype='SpreadSheetView')
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# find source
integrateVariables1 = FindSource('IntegrateVariables1')

# set active source
SetActiveSource(integrateVariables1)

# set active view
SetActiveView(spreadSheetView1)

# set active source
SetActiveSource(integrateVariables1)

# save data
SaveData('/multipede1/khyrm/reactingHetCatFoam/tutorials/verification/innerSpherePartTransp/%s'%data, proxy=integrateVariables1, Precision=8,
    FieldAssociation='Cells')

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).