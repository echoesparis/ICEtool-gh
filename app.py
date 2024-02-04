from flask import Flask
#import rhinoinside
import ghhops_server as hs
from io import StringIO, BytesIO
import rhino3dm
import numpy as np
from scipy.optimize import root
import pandas as pd
#import csv
import time
#import logging
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import base64



# register hops app as middleware
app = Flask(__name__)


hops = hs.Hops(app)

#### #### #### #### #### #### #### #### #### #### 
#                   RADIANCE                    #
#### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### 
#                   FUNCTIONS                   #
#### #### #### #### #### #### #### #### #### ####

def dict_pandas(dict_tree, key_name: str = None, transp: bool = True, resertindex: bool = True):
    df = pd.DataFrame(dict_tree)
    key_list = [int(x[1:-1]) for x in df.columns]
    df.columns = key_list
    if transp == True:
        df = df.T
    if resertindex == True:
        df.reset_index(inplace=True)
    if key_name is not None:
        df.rename(columns = {'index': key_name}, inplace= True)
    df.sort_values(by =key_name, inplace=True )
    
    return df , key_list

# main temperature calculation function
def params_rad (sh, RAD, alb, em, sigma, Tskyb, hc, TairK, lambd, ep, Tint, Cv, Tsurf0, ETo, kc):
 
    # Evaluate each part of the equation
    #a0= -RAD*(1-alb)
    a0= -(RAD*.8*sh*(1-alb)+ RAD*.2)
    c= (em)*sigma
    a1= c* -(Tskyb)**4 
    a2= -hc*TairK
    b1= hc
    b2= (lambd/ep)
    a3= -b2*Tint
    b3= Cv*ep/3600
    a4= -b3*(Tsurf0)
    a5= ETo*kc

    A= a0+a1+a2+a3+a4+a5
    B= (b1+ b2+ b3)
    C= c
    return A,B,C

def objective1(x,A,B,C):  
# Calculates the equation of the energy balance with the variables given
    return A+ B*x + C*((x)**4)


#### #### #### #### #### #### #### #### #### #### 
#           RADIANCE MAIN COMPONENT             #
#### #### #### ######## #### #### #### #### ####

@hops.component(
    #"/temp_calc_rad",
    #name="temp_calc_rad",
    "/surface_temp",
    name="surface_temp",
    description="Calculate Surface Temperature",
    #icon="examples/pointat.png",
    #icon="https://raw.githubusercontent.com/Art-Ev/ICEtool/main/Icons/icon.png",
    inputs=[
        #hs.HopsPoint("MP", "Mesh_Points", "Points of the mesh", access = hs.HopsParamAccess.LIST),
        
        hs.HopsString("mt", "material_tag", "Material tag for the point", access = hs.HopsParamAccess.LIST),
        hs.HopsNumber("db_t", "db_temp", "Dry bulb temperature per hour of the year in Celsious derees", access = hs.HopsParamAccess.TREE),
        hs.HopsNumber("hc", "hc_data", "hc value per point per hour", access = hs.HopsParamAccess.TREE),
        hs.HopsNumber("pre", "precipitation", "Hourly (Hoy) precipitation data", access = hs.HopsParamAccess.TREE),
        hs.HopsNumber("Tn", "Tint", "Temperature of the ground in Celsious derees", access = hs.HopsParamAccess.ITEM, default=35.0),
        hs.HopsString("mi", "materials_info", "Read file output of the materials info csv file", access = hs.HopsParamAccess.LIST, default="" ),
        hs.HopsNumber("RAD", "Global_radiation", "Global radiation for each point in kWh/m²", access = hs.HopsParamAccess.TREE, default=None),
        hs.HopsNumber("sh", "shadows", "shadows per point per hour", access = hs.HopsParamAccess.TREE, default=10),

    ],
    outputs=[
        
        hs.HopsNumber("stK", "Surf_temp_K", "Surface temperature in kelvin ",  access = hs.HopsParamAccess.LIST),
        hs.HopsNumber("stC", "Surf_temp_C", 'Surface temperature in degrees',  access = hs.HopsParamAccess.LIST),
        hs.HopsInteger("Pth", "Mesh_Path", "Mesh point path for the temperature values",  access = hs.HopsParamAccess.LIST),
        hs.HopsInteger("hp", "HOY_Path", "Hour of the Year path for the temperature values",  access = hs.HopsParamAccess.LIST),
    ]
    
)

#def temp_calc_rad( Points, mt, db_temp, hc, precipitation, Tint, materials_info, GH, sh):
def surface_temp( mt, db_temp, hc, precipitation, Tint, materials_info, GH, sh):
    pre_time = time.time()
    
    # dry bulb temperature

    df_db_temp, temp_hours = dict_pandas(db_temp , key_name = 'Hoy')

    df_sh, sh_hoy = dict_pandas(sh , key_name = 'Hoy')
    #df_sh.to_csv('./df_sh.csv')
    df_db_temp.columns =["Hoy", 'degrees']

    df_db_temp['kelvin'] = df_db_temp['degrees'] + 273.15

    # Estimate sky temperature
    Tskyb = round(df_db_temp['kelvin'].mean() ,2) +1

    # Initial temperature in Kelvin
    Tintk = Tint + 273.15
   
   # wind coeficients
    df_hc, hc_hours = dict_pandas(hc , key_name = 'Hoy')

    # CONSTANT VALUES

    sigma = 5.57e-08 #Stefan - Boltzman constant in W.m-2.K-4
    threshold = 0.01 # Calculation threshold, default = 1%
    

    # Calculation of ETo

    df_precipitation, precip_hours = dict_pandas(precipitation , key_name = 'Hoy')
    df_precipitation.rename(columns = {0: 'milimiters'}, inplace= True)
   
    # radiation 

    Gh_df, gh_hours = dict_pandas(GH , key_name = 'Hoy')
    #Gh_df.to_csv('./gh_df.csv')

    # material info

    str_aux = ''
    for i in materials_info:
        str_aux += i + '\n'
    df_mat = pd.read_csv(StringIO(str_aux))
    #df_mat.to_csv('./mat_df.csv')

    # output matrices
    hours = len(temp_hours)
    min_Tsurf= np.ones((len(mt), hours))
    min_TsurfC= np.ones((len(mt), hours))
    min_TsurfC2= np.ones((len(mt), hours))


    start_time = time.time()

    # Solve the equation at each point during all the hours and fill in the tables
    for m in range (len(mt)): # For each point
        #### material constants for the point ####
        # get line that matchs the material tag string
        point_mat = df_mat.loc[df_mat['ID'] == mt[m]]
        # get values in line:
        alb = point_mat["alb"].values[0]
        em = point_mat["em"].values[0]
        lambd = point_mat["lambd"].values[0]
        ep = point_mat["ep"].values[0]
        Cv = point_mat["Cv"].values[0]
        kc = point_mat["kc"].values[0]
        #GH_col = Gh_df[m]
        hc_col = df_hc[m]
        sh_col = df_sh[m]
        


        Tsurf0 = Tintk # Initialize surface temperature
        fin= False
        count=0
        while fin==False: # Start iterations
            for h in range(len(temp_hours)): # Each hour
                GH_value = Gh_df.iat[h,1] # from kWh/m2 to Wh/m2
                #print(GH_value)
                hc_value = hc_col[h]
                sh_value = sh_col[h]
                temp_value = df_db_temp['kelvin'].iat[h]
                preci_value = df_precipitation['milimiters'].iat[h]

                #### Calculate the different parameters ####
               
                A,B,C = params_rad(sh_value, GH_value, alb, em, sigma, Tskyb, hc_value, temp_value, lambd, ep, Tintk, Cv, Tsurf0, preci_value, kc)
                args= (A,B,C)
                

                # Solve equation A + BT + CT^4 = 0
                sol = root(objective1,[Tsurf0+2],args)

                # Save results
                min_Tsurf[m,h]= sol.x[0] # Fill in the table for the concerned point and hour 
                min_TsurfC[m,h]= (min_Tsurf[m,h] - 273.15) # Result in degrees
                min_TsurfC2[m,h]= round(min_Tsurf[m,h] - 273.15,2) # Round result in degrees
            
                T_previous = Tsurf0 # Save previous value
                Tsurf0 = min_Tsurf[m,h] # Update initial value T0
        
            count=count+1
            # print(abs(T_previous-Tsurf0)/Tsurf0)
            # Iterate at least 3 times, check convergence and stop after 20 times  
            if count>=3:
                if abs(T_previous-Tsurf0)/Tsurf0<=threshold:
                    fin=True
            if count>= 20:
                print("Point number : "+ str(m)+ " at hour "+ str(h)+ "did not respect the threshold convergence value")
                fin =True

    Tsurf_min_K = min_Tsurf.flatten().tolist()
    Tsurf_min_C = min_TsurfC.flatten().tolist()

    pth = []
    for i in range(len(mt)):
        pth_aux = [i]*hours
        pth = pth + pth_aux

    hp = temp_hours * len(mt)

    end_time = time.time()

    timelist = [start_time - pre_time ,  end_time - start_time]
    print("ELAPSED TIME: ", timelist)
    return Tsurf_min_K, Tsurf_min_C, pth, hp


#### #### #### #### #### #### #### #### #### #### 
#                   RUN FLASK APP               #
#### #### #### #### #### #### #### #### #### ####

if __name__ == "__main__":
    app.debug = True
    app.run()
