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
#                   KMEANS                     #
#### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### 
#                   FUNCTIONS                   #
#### #### #### #### #### #### #### #### #### ####

# Plot figures
def plot_figs(coords, x, y , n_clusters, sample_silhouette_values,cluster_labels, silhouette_avg, clusterer):
    # Create a subplot with 1 row and 2 columns
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(18, 7)
    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax1.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax1.set_ylim([0, len(coords) + (n_clusters + 1) * 10])
    y_lower = 10
    img_str_list = []
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )
        # Label the silhouette plots with their cluster numbers at the middle
        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax1.set_title("The silhouette plot for the various clusters.")
    ax1.set_xlabel("The silhouette coefficient values")
    ax1.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax1.set_yticks([])  # Clear the yaxis labels / ticks
    ax1.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])

    # 2nd Plot showing the actual clusters formed
    colors = cm.nipy_spectral(cluster_labels.astype(float) / n_clusters)
    ax2.scatter(
        coords[:, x], coords[:, y], marker=".", s=30, lw=0, alpha=0.7, c=colors, edgecolor="k"
    )

    # Labeling the clusters
    centers = clusterer.cluster_centers_
    # Draw white circles at cluster centers
    ax2.scatter(
        centers[:, 0],
        centers[:, 1],
        marker="o",
        c="white",
        alpha=1,
        s=200,
        edgecolor="k",
    )

    for i, c in enumerate(centers):
        ax2.scatter(c[0], c[1], marker="$%d$" % i, alpha=1, s=50, edgecolor="k")

    ax2.set_title("The visualization of the clustered data.")
    ax2.set_xlabel("Feature space for the {} feature".format(x))
    ax2.set_ylabel("Feature space for the {} feature".format(y))

    plt.suptitle(
        "Silhouette analysis for KMeans clustering on sample data with n_clusters = %d"
        % n_clusters,
        fontsize=14,
        fontweight="bold",
    )

    my_stringIO = BytesIO()
    plt.savefig(my_stringIO, format='png')
    my_stringIO.seek(0)
    string_img = base64.b64encode(my_stringIO.read())
    string_img = string_img.decode('utf-8')
    #plt.savefig(fig_path + '\\Silhouette_n_clusters_' + str(n_clusters) + '.png')
    return string_img

#### #### #### #### #### #### #### #### #### #### 
#               KMENAS MAIN COMPONENT           #
#### #### #### #### #### #### #### #### #### ####
@hops.component(
    "/kmeans",
    name="kmeans",
    description="Get K-means",
    #icon="examples/pointat.png",
    inputs=[
        hs.HopsNumber("d", "data", "data to peform K-means on ", access = hs.HopsParamAccess.TREE),
        hs.HopsInteger( "k", "k_values","index of feature x parameter", access = hs.HopsParamAccess.LIST, default=3),
        hs.HopsInteger( "x", "x_index","index of feature x parameter", access = hs.HopsParamAccess.ITEM, default=1),
        hs.HopsInteger( "y", "y_index", "index of feature y parameter", access = hs.HopsParamAccess.ITEM,  default=2),
        hs.HopsInteger('i', 'max_iter', "Maximum number of iterations of the k-means algorithm for a single run.", access = hs.HopsParamAccess.ITEM, default=50),
        hs.HopsNumber("t", "tol", "Relative tolerance with regards to Frobenius norm of the difference in the cluster centers of two consecutive iterations to declare convergence", access = hs.HopsParamAccess.ITEM,  default= 0.0001),
        hs.HopsInteger("rs", "Random_seed", "Determines random number generation for centroid initialization. If 0 or no input is given, the algorithm remains random"), 
        hs.HopsBoolean('s', 'show', 'Make silhouette plots for all values of k' ),
        


    ],
    outputs=[
        hs.HopsPoint("P", "Points", "Point XY",  access = hs.HopsParamAccess.LIST),
        hs.HopsInteger('KL', 'Cluster_label', 'Index of the cluster',  access = hs.HopsParamAccess.LIST),
        hs.HopsNumber('S', 'Point_scores', 'Point Scores for each K',  access = hs.HopsParamAccess.LIST),
        hs.HopsNumber('GS', "G_scores", "Global silhouette scores for each tested k",  access = hs.HopsParamAccess.LIST),
        hs.HopsInteger('Pth', 'Tree_Path', 'Tree Path for Point Scores and Cluster labels',  access = hs.HopsParamAccess.LIST),
        hs.HopsString('Img_as_string', 'img', 'Image as string of 64 bites do be decoded', access = hs.HopsParamAccess.LIST )
    ]
    
)


def kmenas(data, k, x, y, max_iter, tol, Random_seed=0 , s = False):

    if Random_seed == 0:
        Random_seed = None

# points and coords from data
    P = []
    coords =[]
    KL = []
    S =  []
    GS = []
    Tree_Path =[]
    nK = len(k)
    for i in  data.keys():
        x_coord = float(data[i][x])
        y_coord = float(data[i][y])
        coords.append(data[i])
        p_elem = rhino3dm.Point3d(x_coord, y_coord, 0.0)
        P.append(p_elem)
# corrds as an np array
    coords = np.asarray(coords)
# K-means
    # K values loops
    string_img_list =[]
    for n_clusters in k:
        
        Tree_Path = Tree_Path + ([n_clusters] * len(P))
    # Initialize the clusterer with n_clusters value 

        clusterer = KMeans(n_clusters=n_clusters, max_iter= max_iter, tol= tol, random_state= Random_seed)
        cluster_labels = clusterer.fit_predict(coords)
    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
        silhouette_avg = silhouette_score(coords, cluster_labels)
    # Compute the silhouette scores for each sample
        sample_silhouette_values = silhouette_samples(coords, cluster_labels)
        KL.append(cluster_labels)
        S.append(sample_silhouette_values)
        GS.append(silhouette_avg)
        if s == True:
            string_img = plot_figs(coords,x, y , n_clusters, sample_silhouette_values,cluster_labels, silhouette_avg, clusterer)
            #print('string image ---->  ',string_img)
            string_img_list.append(string_img)
    # Formating values to outputs

    KL =  np.asarray(KL).flatten().tolist()
    S = np.asarray(S).flatten().tolist()
    return P , KL ,S , GS, Tree_Path, string_img_list

#### #### #### #### #### #### #### #### #### #### 
#                   RUN FLASK APP               #
#### #### #### #### #### #### #### #### #### ####

if __name__ == "__main__":
    app.debug = True
    app.run()
