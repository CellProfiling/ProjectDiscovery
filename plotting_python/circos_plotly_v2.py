#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:27:34 2017

@author: devinsullivan
"""

import numpy as np
#from IPython.display import Image
import plotly.plotly as py
import plotly.offline as offline
from plotly.graph_objs import *
#from plotly.tools import FigureFactory as FF
import matplotlib as plt
import csv
import pandas as pd
from math import log
import sys
import argparse

np.set_printoptions(threshold=np.nan)
PI=np.pi

def get_cmap(n, name='plasma'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def check_data(data_matrix):
    L, M=data_matrix.shape
    if L!=M:
        raise ValueError('Data array must have (n,n) shape')
    return L


def moduloAB(x, a, b): #maps a real number onto the unit circle identified with 
                       #the interval [a,b), b-a=2*PI
        if a>=b:
            raise ValueError('Incorrect interval ends')
        y=(x-a)%(b-a)
        return y+b if y<0 else y+a
    
def test_2PI(x):
    return 0<= x <2*PI


def get_ideogram_ends(ideogram_len, gap):
    ideo_ends=[]
    left=0
    for k in range(len(ideogram_len)):
        right=left+ideogram_len[k]
        ideo_ends.append([left, right])   
        left=right+gap
    return ideo_ends 


def make_ideogram_arc(R, phi, a=50):
    # R is the circle radius
    # phi is the list of ends angle coordinates of an arc
    # a is a parameter that controls the number of points to be evaluated on an arc
    if not test_2PI(phi[0]) or not test_2PI(phi[1]):
        phi=[moduloAB(t, 0, 2*PI) for t in phi]
    length=(phi[1]-phi[0])% 2*PI 
    nr=5 if length<=PI/4 else int(a*length/PI)

    if phi[0] < phi[1]:   
        theta=np.linspace(phi[0], phi[1], nr)
    else:
        phi=[moduloAB(t, -PI, PI) for t in phi]
        theta=np.linspace(phi[0], phi[1], nr)
    return R*np.exp(1j*theta)   

def map_data(data_matrix, row_value, ideogram_length,L):
    mapped=np.zeros(data_matrix.shape)
    for j  in range(L):
#        mapped[:, j]=ideogram_length*3/row_value#data_matrix[:,j]/row_value
        mapped[:, j]=ideogram_length*2*data_matrix[:,j]/row_value

    return mapped 

def make_ribbon_ends(mapped_data, ideo_ends,  idx_sort,remove_cols):
       
    L=mapped_data.shape[0]
    ribbon_boundary=np.zeros((L,L+1))

    for k in range(L):
        start=ideo_ends[k][0]
        ribbon_boundary[k][0]=start
        #ideo_linspace = np.linspace(ideo_ends[k][0],ideo_ends[k][1],num_sig[k]+1)
        lincounts = 0
        for j in range(1,L+1):
            J=idx_sort[k][j-1]
            if remove_cols[k] or mapped_data[k][J]==0:
                ribbon_boundary[k][j] = start
            else:
#                ribbon_boundary[k][j]=ideo_linspace[lincounts]#start+mapped_data[k][J]
                ribbon_boundary[k][j]=start+mapped_data[k][J]
                lincounts+=1
            start=ribbon_boundary[k][j]
    
    return [[(ribbon_boundary[k][j],ribbon_boundary[k][j+1] ) for j in range(L)] for k in range(L)]   
#    return [[(ribbon_boundary[k][0],ribbon_boundary[k][1] ) for j in range(L)] for k in range(L)]   

def control_pts(angle, radius):
    #angle is a  3-list containing angular coordinates of the control points b0, b1, b2
    #radius is the distance from b1 to the  origin O(0,0) 

    if len(angle)!=3:
        raise InvalidInputError('angle must have len =3')
    b_cplx=np.array([np.exp(1j*angle[k]) for k in range(3)])
    b_cplx[1]=radius*b_cplx[1]
    return zip(b_cplx.real, b_cplx.imag)

def ctrl_rib_chords(l, r, radius):
    # this function returns a 2-list containing control poligons of the two quadratic Bezier
    #curves that are opposite sides in a ribbon
    #l (r) the list of angular variables of the ribbon arc ends defining 
    #the ribbon starting (ending) arc 
    # radius is a common parameter for both control polygons
    if len(l)!=2 or len(r)!=2:
        raise ValueError('the arc ends must be elements in a list of len 2')
    return [control_pts([l[j], (l[j]+r[j])/2, r[j]], radius) for j in range(2)]

def make_q_bezier(b):# defines the Plotly SVG path for a quadratic Bezier curve defined by the 
                     #list of its control points
    if len(b)!=3:
        raise valueError('control poligon must have 3 points')
    A, B, C=b 
#    print('M '+str(A[0])+','+str(A[1])+' '+'Q '+str(B[0])+', '+str(B[1])+' '+str(C[0])+', '+str(C[1]))
    return 'M '+str(A[0])+',' +str(A[1])+' '+'Q '+\
                str(B[0])+', '+str(B[1])+ ' '+\
                str(C[0])+', '+str(C[1])
        

def make_ribbon_arc(theta0, theta1):

    if test_2PI(theta0) and test_2PI(theta1):
        if theta0 < theta1:
            theta0= moduloAB(theta0, -PI, PI)
            theta1= moduloAB(theta1, -PI, PI)
            if theta0*theta1>0:
                raise ValueError('incorrect angle coordinates for ribbon')
    
        nr=int(40*(theta0-theta1)/PI)
        if nr<=2: nr=3
        theta=np.linspace(theta0, theta1, nr)
        pts=np.exp(1j*theta)# points on arc in polar complex form
    
        string_arc=''
        for k in range(len(theta)):
            string_arc+='L '+str(pts.real[k])+', '+str(pts.imag[k])+' '
        return   string_arc 
    else:
        raise ValueError('the angle coordinates for an arc side of a ribbon must be in [0, 2*pi]')
        





def make_layout(title, plot_size):
    axis=dict(showline=False, # hide axis line, grid, ticklabels and  title
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          title=''
          )

    return Layout(title=title,
                  xaxis=XAxis(axis),
                  yaxis=YAxis(axis),
                  showlegend=False,
                  width=plot_size,
                  height=plot_size,
                  margin=Margin(t=0, b=0, l=0, r=0, pad=3),
                  paper_bgcolor='rgba(0,0,0,0)',
                  plot_bgcolor='rgba(0,100,0,0)',
                  hovermode='closest',
                  font=dict(size=20),
                  shapes=[]# to this list one appends below the dicts defining the ribbon,
                           #respectively the ideogram shapes
                 )  
    
def make_ideo_shape(path, line_color, fill_color):
    #line_color is the color of the shape boundary
    #fill_collor is the color assigned to an ideogram
    return  dict(
                  line=Line(
                  color=line_color, 
                  width=0.45
                 ),

            path=  path,
            type='path',
            fillcolor=fill_color,    
        )   

def make_ribbon(l, r, line_color, fill_color, radius=0.2):
    #l=[l[0], l[1]], r=[r[0], r[1]]  represent the opposite arcs in the ribbon 
    #line_color is the color of the shape boundary
    #fill_color is the fill color for the ribbon shape
    #fill_color = np.linspace(0,1,11)
    

    poligon=ctrl_rib_chords(l,r, radius)
    b,c =poligon 
    
#    ribbon_data = [
#            {
#                type:'path',
#                path:make_q_bezier(b)+make_ribbon_arc(r[0], r[1])+make_q_bezier(c[::-1])+make_ribbon_arc(l[1], l[0]),
#                line: {color:line_color,width:0.5},
#                fillcolor: fill_color}]
#                marker: {gradient:{type:'radial',color='rgba(0,0,0,1)'}}       
            
           
    return  dict(
                line=Line(
                color=line_color, width=0.5
            ),
            path=  make_q_bezier(b)+make_ribbon_arc(r[0], r[1])+
                   make_q_bezier(c[::-1])+make_ribbon_arc(l[1], l[0]),
            type='path', 
            fillcolor=fill_color,            
        )
#    return ribbon_data

def make_self_rel(l, line_color, fill_color, radius):
    #radius is the radius of Bezier control point b_1
    b=control_pts([l[0], (l[0]+l[1])/2, l[1]], radius) 
    return  dict(
                line=Line(
                color=line_color, width=0.5
            ),
            path=  make_q_bezier(b)+make_ribbon_arc(l[1], l[0]),
            type='path',
            fillcolor=fill_color,    
        )

def invPerm(perm):
    # function that returns the inverse of a permutation, perm
    inv = [0] * len(perm)
    for i, s in enumerate(perm):
        inv[s] = i
    return inv

def translate_dictnames(dictnames):
    DICT_HASH = dict()
    DICT_HASH['Nucleus'] = 'NUC'
    DICT_HASH['Nucleoplasm'] = 'NP'
    DICT_HASH['Nuclear bodies'] = 'NB'
    #DICT_HASH['Nuclear bodies (many)'] = 'NB-MANY'
    DICT_HASH['Nuclear speckles'] = 'NS'
    DICT_HASH['Nucleoli'] = 'NUI'
    #DICT_HASH['Nucleoli rim'] = 'NUI-RIM'
    DICT_HASH['Nucleoli fibrillar center'] = 'NUI-FC'
    DICT_HASH['Nuclear membrane'] = 'NM'
    DICT_HASH['Cytosol'] = 'CYT'
    DICT_HASH['Cytoplasm'] = 'CYT'
    DICT_HASH['Aggresome'] = 'AGG'
    DICT_HASH['Mitochondria'] = 'MIT'
    #DICT_HASH['Rods & Rings'] = 'R&R'
    DICT_HASH['Intermediate filaments'] = 'IF'
    DICT_HASH['Microtubule ends'] = 'MT-E'
    DICT_HASH['Microtubule end'] = 'MT-E'
    DICT_HASH['Microtubules'] = 'MT'
    DICT_HASH['Actin filaments'] = 'AF'
    DICT_HASH['Cytokinetic bridge'] = 'CKB'
    DICT_HASH['Microtubule organizing center'] = 'MTOC'
    DICT_HASH['Centrosome'] = 'CEN'
    DICT_HASH['Endoplasmic reticulum'] = 'ER'
    DICT_HASH['Golgi apparatus'] = 'GA'
    DICT_HASH['Vesicles'] = 'VES'
    DICT_HASH['Cell Junctions'] = 'CJ'
    DICT_HASH['Focal adhesion sites'] = 'FA'
    DICT_HASH['Focal Adhesions'] = 'FA'
    DICT_HASH['Plasma membrane'] = 'PM'
    DICT_HASH['Nucleoli (Fibrillar center)'] = 'NUI-FC'
    DICT_HASH['Nucleoli (fib center)'] = 'NUI-FC'
    DICT_HASH['Nucleoli fib center'] = 'NUI-FC'    
    DICT_HASH['Cytoskeleton (Intermediate filaments)'] = 'IF'
    DICT_HASH['Cytoskeleton (Microtubule end)'] = 'MT-E'
    DICT_HASH['Cytoskeleton (Microtubules)'] = 'MT'
    DICT_HASH['Cytoskeleton (Actin filaments)'] = 'AF'
    DICT_HASH['Cytoskeleton (Cytokinetic bridge)'] = 'CKB'
    
    
    dispnames = []
    for name in dictnames:
        try:
            dispnames.append(DICT_HASH[name])
        except:
            dispnames.append('')
            
    return dispnames


def main():
    
    parser = argparse.ArgumentParser(description='Create circos overannotation plot.')
    parser.add_argument('data_path', metavar='data path', type=str,
                    help='path to the input data file')
    parser.add_argument('number_path', metavar='number path', type=str, 
                    help='path to the input counts file')
    parser.add_argument('distance_path', metavar='distance path', type=str,
                    help='path to the input tree distance file')
    parser.add_argument('out_path', metavar='output path', type=str,
                    help='path where we will save the data')
    args = parser.parse_args()
    
#Define data paths 
    #Heatmap data
    data_path = args.data_path#'/Users/devinsullivan/Documents/PD_paper/pd_results_2017/binomial_heatmapv14_v3.txt'
    #number of items in each class
    number_path = args.number_path#'/Users/devinsullivan/Documents/PD_paper/pd_results_2017/tot_hparesult_ml0.txt'
    #distance on heirarchical tree
    distance_path = args.distance_path#'/Users/devinsullivan/Documents/PD_paper/distances.csv'
    out_path = args.out_path
    overann_num = 5
    
    #Read the data
    with open(data_path, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        labels = next(spamreader)
        
    #trim off the blank space
    labels = labels[1:]  
    print(labels)
    
    #Read data in
    data=pd.read_csv(data_path, sep=',')

    #convert data to np array
    npdata = np.array(data)
    npdata = np.delete(npdata,0,1)
    npdata = npdata.astype(float)

    #Get display labels    
    disp_labels = translate_dictnames(labels)
    print(disp_labels)
    
    
    #If the label is not recognized as a display label, remove it. 
    cols_to_remove = [False]*len(disp_labels)
    for i,label in enumerate(disp_labels):
        if label == '':
            cols_to_remove[i] = True
    
    #Remove items from the data matrix
    #matrix_tmp=npdata[np.logical_or(~np.all(npdata == 0, axis=1),~np.all(npdata == 0, axis=0))]
    #matrix=matrix_tmp[:,np.logical_or(~np.all(npdata == 0, axis=1),~np.all(npdata == 0, axis=0))]
    np_keep = ~np.array(cols_to_remove, dtype = bool)
    matrix_tmp=npdata[:,np_keep]
    matrix=matrix_tmp[np_keep,:]
    #Set non-significant enrichment to 0
    matrix[matrix<=2]=0

    
    #Remove items from the label sets
    labels = [d for (d, remove) in zip(labels, cols_to_remove) if not remove]
    disp_labels = [d for (d, remove) in zip(disp_labels, cols_to_remove) if not remove]
    num_labels = len(labels)    
    
        
    #normalize matrix
    matrix_raw = matrix
    matrix = np.divide(matrix,np.amax(matrix))   
    
    #read in how many of each class we have    
    with open(number_path, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        alsolabels = next(spamreader)
        number_data = next(spamreader)
    number_data = [float(d) for d in number_data]
    num_matrix = [d for (d, remove) in zip(number_data, cols_to_remove) if not remove]
    
    #get the distances between each node in the heirarchical tree
    distance_array_wnames = np.genfromtxt(distance_path, delimiter=',',names=True)
    distance_names = distance_array_wnames.dtype.names
    distance_names = [d.replace('_', ' ') for d in distance_names]
    disp_distance_names = translate_dictnames(distance_names)
    distance_array_raw = np.genfromtxt(distance_path, delimiter=',')

    #Go through labels and find each label in distance names
    ind_order = len(labels)*[None]
    for j,label in enumerate(disp_labels):
        for i,dist_name in enumerate(disp_distance_names):
            if dist_name.replace('_',' ') == label:
                ind_order[j] = i
                break
    #fill in distance matrix
    distance_array = np.zeros([len(ind_order),len(ind_order)])
    for i,ind1 in enumerate(ind_order):
        if ind1 == None:
            continue
        for j,ind2 in enumerate(ind_order):
            if ind2 == None:
                continue
            distance_array[i,j] = distance_array_raw[ind1][ind2]
    distance_array = np.divide(distance_array,np.amax(distance_array))


    
    #set up ideograms with gaps 
    L=check_data(matrix)
    row_sum = [log(num,2) for num in num_matrix]

    #set the gap between two consecutive ideograms
    gap=2*PI*0.005
    ideogram_length=2*PI*np.asarray(row_sum)/sum(row_sum)-gap*np.ones(L)
    ideogram_length[ideogram_length<=0] = 0   
    #get endpoints of ideograms
    ideo_ends=get_ideogram_ends(ideogram_length, gap)
        
    
    #IDEA: make the colors here proportional to the precision for each class?
    # Make a list of colors cycling through the default series.
    ideo_cm = get_cmap(num_labels,'Greens')
    ideo_colors_tmp = []
    for i, lab in enumerate(labels):
        ideo_colors_tmp.append(ideo_cm(i))
    ideo_colors_tmp = [tuple(item for item in color) for color in ideo_colors_tmp]
    ideo_colors =[]
    for k in range(L):
        ideo_colors.append('rgba'+str(ideo_colors_tmp[k]))
   
    #
    #z=make_ideogram_arc(2.3, [11*PI/6, PI/17])
    #print z

    #get all significant entries and map them to the ideogram placement    
    bool_mat = matrix>0
    bool_mat = bool_mat.astype(float)
    mapped_data=map_data(bool_mat, row_sum, ideogram_length,L)
    idx_sort=np.argsort(mapped_data, axis=1)

    #find overannotation inds 
#    overann_cols = [False]*len(labels)
    overann_cols = np.sum(bool_mat,1)>overann_num
 
#    for ann in over_annotated:
#        overann_cols = np.logical_or([s==ann for s in labels], overann_cols)

    #define ribbon ends for each connection
    ribbon_ends=make_ribbon_ends(mapped_data, ideo_ends,  idx_sort,overann_cols)
    np.set_printoptions(threshold=np.nan)
    np.set_printoptions(linewidth=120)


    print 'ribbon ends starting from the ideogram[2]\n', ribbon_ends[2]
    print 'num colors: ', len(ideo_colors)
    
    #get continuous colormap
    cont_cmap=plt.cm.plasma
    
    #create a basic layout object for the circos plot    
    layout=make_layout('', 400)
    

    #Create ribbons
    ribbon_info=[]
    for k in range(L):
        
        sigma=idx_sort[k]
        sigma_inv=invPerm(sigma)
        
        for j in range(L):
        
            eta=idx_sort[j]
            eta_inv=invPerm(eta)
            l=ribbon_ends[k][sigma_inv[j]]  
            
            if j==k:
                continue
            else:
                r=ribbon_ends[j][eta_inv[k]]

                if matrix[k][j]<sys.float_info.min:
                    continue
                if matrix[j][k]<sys.float_info.min:
                    lstr = list(r)
                    lstr[1] = lstr[0]#+0.01
                    r =tuple(lstr)
                print('matrix: ',matrix[k][j])  
                print('r: ',r)
                print('l: ',l)
                print('label k: ',disp_labels[k])
                print('label j: ',disp_labels[j])
                print('distance: ',distance_array[k][j])
                zi=0.9*np.exp(1j*(l[0]+l[1])/2)
#                zf=0.9*np.exp(1j*(r[0]+r[1])/2)
                #texti and textf are the strings that will be displayed when hovering the mouse 
                #over the two ribbon ends
                texti=disp_labels[k]+' coannotated p<10^-'+ '{:.2f}'.format(matrix_raw[k][j])+' of '+\
                      disp_labels[j],disp_labels[k],
                if overann_cols[k]:
                    lst = list(l)
                    lst[1]=lst[0]
                    l = tuple(lst)
                    lstr = list(r)
                    lstr[1] = lstr[0]
                    r =tuple(lstr)
                    curr_color = 'rgba(200,200,200,0.3)'
                    
                else:
                    curr_color = 'rgba'+str(cont_cmap(distance_array[k][j],alpha=0.8))


                ribbon_info.append(Scatter(x=zi.real,
                                           y=zi.imag,
                                           mode='markers',
                                           marker=Marker(size=0.3, color=curr_color),
                                           text=texti,
                                           hoverinfo='text'
                                           )
                                  )#,

                r=(r[1], r[0])#IMPORTANT!!!  Reverse these arc ends because otherwise you get
                              # a twisted ribbon
                #append the ribbon shape
                layout['shapes'].append(make_ribbon(l, r, 'rgb(175,175,175)' , curr_color))                
  
    ideograms=[]
    for k in range(len(ideo_ends)):
        print(k)
        z= make_ideogram_arc(1.1, ideo_ends[k])
        zlab_r = 1.25
        zlab= make_ideogram_arc(zlab_r, ideo_ends[k])
        zi=make_ideogram_arc(1.0, ideo_ends[k])
        m=len(z)
        n=len(zi)
        ideograms.append(Scatter(x=z.real,
                                 y=z.imag,
                                 mode='lines',
                                 line=Line(color='rgba'+str(ideo_colors[k]), shape='spline', width=0.25),
                                 text=disp_labels[k]+'<br>'+'{:.2f}'.format(row_sum[k]),
                                 hoverinfo='text'
                                 )
                         )
    
#        ideograms.append(Scatter(x=z.real[0],y=z.imag[0],mode='text',text=labels[k]))
        ideograms.append(Scatter(x=np.mean(zlab.real),y=np.mean(zlab.imag),mode='text',text=disp_labels[k]))
             
#        print(ideograms)
        path='M '
        for s in range(m):
            path+=str(z.real[s])+', '+str(z.imag[s])+' L '
            
        Zi=np.array(zi.tolist()[::-1]) 
    
        for s in range(m):
            path+=str(Zi.real[s])+', '+str(Zi.imag[s])+' L '
        path+=str(z.real[0])+' ,'+str(z.imag[0]) 
       
        layout['shapes'].append(make_ideo_shape(path,'rgb(150,150,150)' , ideo_colors[k]))
#        print(vars(layout))
        lab_x = np.mean(zlab.real)
        lab_y = np.mean(zlab.imag)
        if lab_x>0 and lab_y>0:
            actual_angle = np.rad2deg(np.arcsin(np.mean(zlab.imag)/zlab_r))
        elif lab_x<0:
            actual_angle = 180-np.rad2deg(np.arcsin(np.mean(zlab.imag)/zlab_r))
        else:
            actual_angle = 360+np.rad2deg(np.arcsin(np.mean(zlab.imag)/zlab_r))
        print(actual_angle)
        if 90<k*360/len(ideo_ends)<270:
            layout['annotations'].append(dict(x=np.mean(zlab.real),y=np.mean(zlab.imag),text=disp_labels[k],showarrow=False,textangle=-actual_angle+180))
        else:
            layout['annotations'].append(dict(x=np.mean(zlab.real),y=np.mean(zlab.imag),text=disp_labels[k],showarrow=False,textangle=-actual_angle))
                                    
            

    ideograms = []
    data = Data(ideograms+ribbon_info)
    fig = Figure(data=data, layout=layout)

    offline.init_notebook_mode(connected=True)

    offline.plot({'data': data,'layout': layout}, filename='testing123.html', image='png',image_height=1000,image_width=1000)

    #py.image.save_as(fig,filename=out_path,scale=3)
    #Image('a-simple-plot.png')

    #url = py.plot(fig, filename='chord-diagram-Fb') 
    
    
    
    
    
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
