# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:09:26 2022

@author: SArden
"""
#OWA Script, locally referenced for packaging with ArcPro project

#import packages
import pandas as pd
import arcpy
import numpy as np
from gekko import GEKKO
m=GEKKO()
import geopandas as gpd
import os

# Set overwrite option
workspace=arcpy.env.workspace
workspace_new=workspace.replace("PRESTo.gdb","Data")
arcpy.env.workspace=workspace_new
arcpy.env.overwriteOutput = True

#eventually we will have an AOI function clip an input feature class
##code to index clipped feature class
##arcpy.management.AddIndex(in_table, fields, {index_name}, {unique}, {ascending})

#define input table of available indicators, set Census Tract (first column) as index
indvaltable=arcpy.GetParameterAsText(0) #csv file, formatted per template
df_indvalues=pd.read_csv(indvaltable,dtype={'GEOID':object}) #convert data to pandas data frame
df_indvalues.set_index('GEOID',inplace=True)

#pull input from the user to select among the available indicators
indicators_raw=arcpy.GetParameterAsText(1) #GetPar... pulls in data as string, need to convert to list for our purposes
indicators=indicators_raw.split(";")
df1=df_indvalues[indicators]

#pull ORness from user input
OR=arcpy.GetParameterAsText(2)

####Weighting subroutine NEEDS TO BE SWAPPED FOR NEW PACKAGE
n=df1.shape[1] #need a 1xn array where n is the number of user-selected indicators
vl=[None]*n

# syntax: x1 = m.Var(value=1,lb=1,ub=5)
vl[0] = m.Var(value=0.1, lb=0, ub=1)

# assigning all values in vl as gekko variables to be solved
for i in range(n):
    vl[i] = m.Var(value=1/n, lb=0, ub=1)

# Use gekko array instead of python list
x = m.Array(m.Var,(n))
# intial guess
ig = [1/n for i in range(n)]
# m.Var definition
i = 0
for xi in x:
    xi.value = ig[i]
    xi.lower = 0.00000001
    xi.upper = 1
    i += 1

# m.Equation() <-- sum of weights == 1
# m.Equation() <-- (1/(n-1))*sum of [(n-rnk)*weight] == ORNess
# m.Maximize(obj) <-- this is the maximum entropy expression to be maximized
m.Equation(sum(x)==1)
m.Equation((1/(n-1))*sum((n-(1+i))*x[i] for i in range(n))==OR)
m.Maximize(-sum(x[i]*m.log(x[i]) for i in range(n)))

#solve the equations, write to 
m.solve()
OW=np.array([elem for singleList in x for elem in singleList])#converts from list objects to values in solution array
####Back to the main program

#calculate results in pandas dataframe
df1_rnk = df1.rank(1,'first',ascending=False) #generate rank matrix corresponding to user-selected indicators
df1_rnkWt = df1_rnk.applymap(lambda x: OW[int(x-1)]) #apply weights
df1_contrib = df1.mul(df1_rnkWt) #these are the individual indicator contributions to the final OWA score
df1_contrib.reset_index(inplace=True) #this addes an object ID column back in
df1_contrib.set_index('GEOID',inplace=True) #set index to GEOID
df1_OWA=df1_contrib.sum(axis=1)
df2_OWA=pd.DataFrame({'GEOID':df1_OWA.index,'OWA':df1_OWA.values}) #adding column headings
df2_OWA.set_index('GEOID',inplace=True)
#read base shapefile, (over)write results shapefile
shapefile_base=os.path.join(workspace_new,"BostonTracts_Base.shp")
base = gpd.read_file(shapefile_base)
base.set_index('GEOID',inplace=True)
base=base.join(df2_OWA,on='GEOID')
base=base.join(df1_contrib,on='GEOID')
shapefile_results=os.path.join(workspace_new,"BostonTracts_Results.shp")
base.to_file(shapefile_results)

#map display code - define project, map, layer, layout, symbology
p=arcpy.mp.ArcGISProject('current') #define p as current project
m=p.listMaps('Tool')[0] #define which map you are working in
lyr=m.listLayers('BostonTracts_Results')[0] #define which layer you want to edit
lyt=p.listLayouts("Tool")[0] #set current layout, point to layout name
sym=lyr.symbology
cp=lyr.connectionProperties #define layer connection properties so we can update parts of it later

#text element cycle
titletext=('ORness = '+ OR)
tableText=lyt.listElements("TEXT_ELEMENT","ORness")[0] #refers to a text element of Name=Title
tableText.text=titletext #replace "Title" text element with whatever string of titles you want

#map symbology
cp_new={'dataset': 'BostonTracts_Results.shp', 'workspace_factory': 'Shape File', 'connection_info': {'database': workspace_new}}#changes the database portion of cp dictionary to the new relative path
lyr.updateConnectionProperties(cp,cp_new) #regardless of the old data source, this should update with the new data source 
fields=base.columns
sym.updateRenderer('SimpleRenderer') #for some reason not resetting the renderer causes updated layer to disappear
sym.updateRenderer('GraduatedColorsRenderer') 
sym.renderer.classificationField=fields[1] #make sure to use the field name and not the alias
sym.renderer.colorRamp=p.listColorRamps("GreenYellowRed")[0]
#class interval definition
classBreakValues=[0.2,0.4,0.6,0.8,1]
classBreakLabels=["0.2","0.4","0.6","0.8","1"]
sym.renderer.breakCount=len(classBreakValues)
count=0
for brk in sym.renderer.classBreaks:
    brk.upperBound=classBreakValues[count]
    brk.label=classBreakLabels[count]
    count+=1
lyr.symbology=sym