# python script to compare exact solution to Abaqus solution
# 1/ run with command line
#    abaqus cae -script = comparison.py
#    in folder where patch_XXX.odb file is
# 2/ user only needs to modify the variables below:
#    - element_type (l. 18): type of element used
#    - strain_component (l. 19): strain component for stress-strain plot
#    - stress_component (l. 20): stress component for stress-strain plot
#    - plot_title (l. 21): titple for plots, e.g., loading type
# 3/ generates and saves .png files of stress-strain plot and plot of 
#    strain energies (per deformed unit volume) in the equilibrium (psiE) 
#    and non-equilibrium (psiNE) branches
import csv
from abaqus import *
from abaqusConstants import *
from viewerModules import *

element_type='CPEG4H'
strain_component='LE22'
stress_component='S22'
plot_title='Uniaxial tension'

##################################################################
# do not modify after this line
##################################################################
jobname='patch_'+element_type
SvsLE=stress_component+' vs '+strain_component

# opens odb file and preps Viewport
o2 = session.openOdb(name=jobname+'.odb')
myViewport=session.viewports['Viewport: 1']
myViewport.setValues(displayedObject=o2)
odb = session.odbs[jobname+'.odb']

#instance and element info
instance_name=odb.rootAssembly.instances.keys()[0]
el_label=str(odb.rootAssembly.instances[instance_name].elements[0].label)

# extracts LE, S (Cauchy stress), SENER, PENER
session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_CENTROID, 
    variable=(('LE', INTEGRATION_POINT, ((COMPONENT, strain_component), )), ('S', 
    INTEGRATION_POINT, ((COMPONENT, stress_component), )), ('SENER', INTEGRATION_POINT), ('PENER', INTEGRATION_POINT),), elementLabels=
    ((instance_name, (el_label, )), ))

# renames LE    
xy1 = session.xyDataObjects['LE:'+strain_component+' PI: '+instance_name+' E: '+el_label+' Centroid']
xy1.setValues(
    sourceDescription='Log strain vs time')
tmpName = xy1.name
session.xyDataObjects.changeKey(tmpName, strain_component+' vs t')
    
# renames S 
xy1 = session.xyDataObjects['S:'+stress_component+' PI: '+instance_name+' E: '+el_label+' Centroid']
xy1.setValues(
    sourceDescription='Cauchy stress vs time')
tmpName = xy1.name
session.xyDataObjects.changeKey(tmpName, stress_component+' vs t')

# combines LE and S 
xy1=combine(session.xyDataObjects[strain_component+' vs t'],session.xyDataObjects[stress_component+' vs t'])
xy1.setValues(
    sourceDescription='Cauchy stress vs Log strain') 
tmpName = xy1.name
session.xyDataObjects.changeKey(tmpName, element_type+': '+SvsLE)

# renames SENER    
xy1 = session.xyDataObjects['SENER PI: '+instance_name+' E: '+el_label+' Centroid']
xy1.setValues(
    sourceDescription='Strain energy density per unit deformed volume in equilibrium branch vs time')
tmpName = xy1.name
session.xyDataObjects.changeKey(tmpName, element_type+': psiE vs t')
    
# renames PENER 
xy1 = session.xyDataObjects['PENER PI: '+instance_name+' E: '+el_label+' Centroid']
xy1.setValues(
    sourceDescription='Strain energy density per unit deformed volume in non-equilibrium branch vs time')
tmpName = xy1.name
session.xyDataObjects.changeKey(tmpName, element_type+': psiNE vs t')

# opens exact solution
with open('exact_solution.csv') as f:
	x = csv.reader(f)
	data = list(x)
    
# column indexes in exact solution table
indexes = {'t': 0, 'LE11': 1, 'LE22': 2, 'LE33': 3, 'LE12': 4, 'LE13': 5, 'LE23': 6, 'J': 7,
                   'S11': 8, 'S22': 9, 'S33': 10, 'S12': 11, 'S13': 12, 'S23': 13, 'psiE': 14, 'psiNE': 15}
    
# import S vs LE exact solution    
exact_SvsLE=[(float(row[indexes[strain_component]]),float(row[indexes[stress_component]])) for row in data[1:]]  
session.XYData('Exact: '+SvsLE,exact_SvsLE,'Exact solution '+SvsLE)
xQuantity = visualization.QuantityType(type=STRAIN)
yQuantity = visualization.QuantityType(type=STRESS)
session.xyDataObjects['Exact: '+SvsLE].setValues(
    axis1QuantityType=xQuantity, axis2QuantityType=yQuantity, )

# comparison S vs LE plot 
xyp = session.XYPlot('XYPlot-1')
session.charts['Chart-1'].gridArea.style.setValues(color='#FFFFFF')
session.xyPlots['XYPlot-1'].title.style.setValues(color='#000000')
session.xyPlots['XYPlot-1'].title.setValues(
    text=plot_title)
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects[element_type+': '+SvsLE]
c1 = session.Curve(xyData=xy1)
xy2 = session.xyDataObjects['Exact: '+SvsLE]
c2 = session.Curve(xyData=xy2)
chart.setValues(curvesToPlot=(c1, c2, ), )
session.curves[element_type+': '+SvsLE].lineStyle.setValues(style=DASHED)
session.curves[element_type+': '+SvsLE].lineStyle.setValues(thickness=0.5)
myViewport.setValues(displayedObject=xyp)
session.charts['Chart-1'].axes1[0].axisData.setValues(useSystemTitle=False, 
    title='Log strain '+strain_component)
session.charts['Chart-1'].axes2[0].axisData.setValues(useSystemTitle=False, 
    title='Cauchy stress '+stress_component+' (MPa)')
session.printToFile(fileName=SvsLE, format=PNG,
    canvasObjects=(myViewport,))
    
# import psiE vs t and psiNE vs t exact solutions    
exact_psiEvst=[(float(row[indexes['t']]),float(row[indexes['psiE']])) for row in data[1:]] 
exact_psiNEvst=[(float(row[indexes['t']]),float(row[indexes['psiNE']])) for row in data[1:]]  
session.XYData('Exact: psiE vs t',exact_psiEvst,"Exact solution psiE vs t") 
session.XYData('Exact: psiNE vs t',exact_psiNEvst,"Exact solution psiNE vs t")
xQuantity = visualization.QuantityType(type=TIME)
yQuantity = visualization.QuantityType(type=ENERGY_DENSITY)
session.xyDataObjects['Exact: psiE vs t'].setValues(
    axis1QuantityType=xQuantity, axis2QuantityType=yQuantity, )
session.xyDataObjects['Exact: psiNE vs t'].setValues(
    axis1QuantityType=xQuantity, axis2QuantityType=yQuantity, )

# comparison psiE vs t and psiNE vs t plot 
xyp = session.XYPlot('XYPlot-2')
session.charts['Chart-2'].gridArea.style.setValues(color='#FFFFFF')
session.xyPlots['XYPlot-2'].title.style.setValues(color='#000000')
session.xyPlots['XYPlot-2'].title.setValues(
    text=plot_title)
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
xy1 = session.xyDataObjects[element_type+': psiE vs t']
c1 = session.Curve(xyData=xy1)
xy2 = session.xyDataObjects['Exact: psiE vs t']
c2 = session.Curve(xyData=xy2)
xy3 = session.xyDataObjects[element_type+': psiNE vs t']
c3 = session.Curve(xyData=xy3)
xy4 = session.xyDataObjects['Exact: psiNE vs t']
c4 = session.Curve(xyData=xy4)
chart.setValues(curvesToPlot=(c1, c2, c3, c4), )
session.curves[element_type+': psiE vs t'].lineStyle.setValues(style=DASHED)
session.curves[element_type+': psiE vs t'].lineStyle.setValues(thickness=0.5)
session.curves[element_type+': psiNE vs t'].lineStyle.setValues(style=DASHED)
session.curves[element_type+': psiNE vs t'].lineStyle.setValues(thickness=0.5)
session.curves[element_type+': psiE vs t'].lineStyle.setValues(color='#FF0000')
session.curves['Exact: psiE vs t'].lineStyle.setValues(color='#003366')
session.curves[element_type+': psiNE vs t'].lineStyle.setValues(color='#008000')
session.curves['Exact: psiNE vs t'].lineStyle.setValues(color='#FF7F00')
myViewport.setValues(displayedObject=xyp)
session.charts['Chart-2'].axes2[0].axisData.setValues(useSystemTitle=False, 
    title='Str. en. dens. per unit def. vol. (J/mm^3)')
session.charts['Chart-2'].axes1[0].axisData.setValues(useSystemTitle=False, 
    title='Time (s)')
session.printToFile(fileName='psis vs t', format=PNG,
    canvasObjects=(myViewport,))