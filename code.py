# import 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy import stats
import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import os.path
import random
import matplotlib.backends.backend_pdf
import os
from PIL import Image
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline
from heapq import nsmallest


raw_wafer_All = pd.read_csv(path + 'slim_TEOS_109PT_polar.csv')
raw_wafer = raw_wafer_All[['Point', 'x', 'y', 'Tilt43k_8']]

raw_wafer.columns = ['Point', 'x', 'y', 'value']

raw_wafer.head()



def choose_radian_point(raw_wafer, i, arc_divided):    
    #angle 계산하기
    myatan = lambda x,y: np.pi*(1.0-0.5*(1+np.sign(x))*(1-np.sign(y**2))\
    -0.25*(2+np.sign(x))*np.sign(y))\
    -np.sign(x*y)*np.arctan((np.abs(x)-np.abs(y))/(np.abs(x)+np.abs(y)))
    theta = myatan(raw_wafer['x'], raw_wafer['y'])
    raw_wafer['angle'] = (theta * 180.0/np.pi) # converting to degrees
    raw_wafer['angle']= raw_wafer['angle'].fillna(360)
        
    #distance 계산하기
    mydist = lambda x,y: (x**2 + y**2)**(1/2)
    raw_wafer['distance'] = mydist(raw_wafer['x'], raw_wafer['y'])
            
    #360을 더한 angle열 추가
    raw_wafer['angle2'] = raw_wafer['angle'] + 360
                    
    #데이터 추출시작시 cemter 제외 0으로 맞추기
    def split_start(raw_wafer):
        raw_wafer['split'] = 0
        raw_wafer.loc[ (raw_wafer['angle'] == 360), 'split'] = 1
        
        
             #######################빈공간 추출 추가#########3       
    def new_xy(arc_divided, distance1, distance2, distance3, i, zero_pi):
        x1 = distance3
        degree = (180/(arc_divided))*i
        rad = math.pi*degree / 180
        new_y = int(x1 * math.sin(rad))
        new_x = int(x1 * math.cos(rad))
        xy = [[new_x, new_y]]
        xy2 = pd.DataFrame(xy, columns = ['x', 'y'])

        For_angle = raw_wafer[(raw_wafer['distance']>distance1) & (raw_wafer['distance']<distance2)]
        near_angle = nsmallest( 1, For_angle['angle'] , key= lambda x: abs(x- ( zero_pi+(180/arc_divided)*i )) )
        near_angle = pd.DataFrame(near_angle, columns = ['angle'])     
        #제일 근처 각도의 지점에서 값 추출
        xy_value = raw_wafer[(raw_wafer['distance']> distance1) & (raw_wafer['distance']< distance2) &
                                (raw_wafer['angle'] == near_angle.iloc[0,0])]
        xy_value.reset_index(drop=True, inplace=True) #index 제거
        
        #좌표와 새 value 취합
        theta = myatan(xy2['x'], xy2['y'])
        xy2['angle'] = (theta * 180.0/np.pi)
        xy2['distance'] =  mydist(xy2['x'], xy2['y'])
        xy2['value'] = xy_value['value']
        xy2['split'] = 1

        return(xy2)
    
#  value 추출 시작

    split_start(raw_wafer)
    split2 = raw_wafer[raw_wafer['split'] == 1 ]
    
    angles = []
    for k in range(0, 5):
        interval = 150 - 30*k
        interval2 = 120 - 30*k
        newxy = new_xy(arc_divided, interval2, interval, (interval2+interval)/2, i, 0 )
        newxy2 = new_xy(arc_divided, interval2, interval, -(interval2+interval)/2, i, 180 )        

        raw_wafer = pd.concat([raw_wafer, newxy,  newxy2])

    newxy3 = new_xy(arc_divided, 0, 25, 0, i, 180 )
    raw_wafer = pd.concat([raw_wafer, newxy3])
    split2= raw_wafer[raw_wafer['split'] == 1]

    #1, 3분면은 재정렬, 2, 4분면은 y축만 역순 정렬 츠기힘
    coordinate = pd.concat((split2['x'], split2['y'],  split2['distance'],
                            split2['value'], split2['angle']), axis = 1).values
    Temp_coordinate = coordinate.tolist()

    if (i<=arc_divided/2):
        sorted_value = sorted(Temp_coordinate, key = lambda k : [k[0], k[1]])
    else:
        sorted_value = sorted(Temp_coordinate, key = lambda k : [k[0], -k[1]])
        
    pd_sorted_value = pd.DataFrame(sorted_value, columns = ('x', 'y','distance', 'value', 'angle'))     
    result = pd_sorted_value['value'].values
    xy_before = pd_sorted_value
    
    #좌표값 확인시 사용할 부분
    groups = raw_wafer.groupby('split')
    fig, ax = plt.subplots()
    for name, group in groups:
        ax.plot(group.x, group.y, marker = 'o', linestyle = "", label = name)

    
    return (result, xy_before)
  
  def match_coordinate(micro_sub_wafer, i, arc_divided, xy_before):
        #(X1, y1) (x2, y2)정의
    x1 = xy_before.loc[[0],['x']].values
    x2 = xy_before.loc[[len(xy_before)-1],['x']].values
    y1 = xy_before.loc[[0],['y']].values
    y2 = xy_before.loc[[len(xy_before)-1],['y']].values
    
    xy_after = pd.DataFrame(columns = ['x', 'y'])
    for n in range(0,63):
        n_Total = 63
        if x1 < 0:
            x3 = - (abs(x1) - abs(x1- x2) * n/(n_Total - 1))
        else:
            x3 =  (abs(x1) - abs(x1- x2) * n/(n_Total - 1))

        if y1 < 0:
            y3 = - (abs(y1) - abs(y1- y2) * n/(n_Total - 1))
        else:
            y3 =  (abs(y1) - abs(y1- y2) * n/(n_Total - 1))  

        x3y3 = pd.DataFrame(np.concatenate((x3, y3), axis = 1),  columns = ['x', 'y'])
        xy_after = xy_after.append(x3y3)
        
    xy_after2 = pd.DataFrame(xy_after).values
    micro_sub_wafer_temp = pd.DataFrame(micro_sub_wafer).values
    result = pd.DataFrame(np.concatenate((xy_after, micro_sub_wafer_temp), axis = 1),  
                          columns = ['x', 'y', 'value'])

    return(result)
  
  
  def load_data():
    path = '//10.10.10.113/반도체사업본부/기술개발총괄/시뮬레이션팀/개인/백승혜_양지현/WaferMap/'
    csv_data = pd.read_csv(path + 'xscan_data_109PT_removeEdge.csv', header = None)
    
    return csv_data
  
  
  def apply_micro_variation( macro_sub_wafer ):
    dist_data = load_data()
    micro_sub_wafer = np.zeros(63)
    init_idx = 1
    idx_unit = int(63 /len(macro_sub_wafer))
    init_T = []
    
    for i in range(len(macro_sub_wafer)-1):
        if i == 0:
            micro_sub_wafer[0] = macro_sub_wafer[i]
        else:
            init_idx = i * idx_unit + i
            init_T.append(init_idx)
            micro_sub_wafer[init_idx] = macro_sub_wafer[i]
            
    for i in range(len(micro_sub_wafer)):
        if i == 0:
            continue
        else:
            data = dist_data.iloc[i-1]
            estimator = stats.gaussian_kde(data, bw_method='silverman')
            sample_dist = estimator.resample(1000)[0]
            sample_dist1 = np.mean(sample_dist)
            if micro_sub_wafer[i] != 0:
                micro_sub_wafer[i] = micro_sub_wafer[i]
                interpolate_unit = (micro_sub_wafer[i] - micro_sub_wafer[i-5]) / 5                        
            elif micro_sub_wafer[i] == 0:
                temp_num1 = micro_sub_wafer[i-1] * sample_dist1
                micro_sub_wafer[i] = temp_num1
            
    #### spline interpol 추가 ####            
#     k = len(init_T)
    
#     for j in init_T:
#             y = [micro_sub_wafer[j-3],micro_sub_wafer[j-2],micro_sub_wafer[j-1],micro_sub_wafer[j]]
#             x = np.linspace(j-3, j, 4)
#             cs = CubicSpline(x, y)
#             xnew = np.linspace(j-3, j, 8)
#             micro_sub_wafer[j-7:j+1] = cs(xnew)

#     for j in init_T:
#         if max([micro_sub_wafer[j-6],micro_sub_wafer[j-3], micro_sub_wafer[j]]) > micro_sub_wafer[j]:
#             y = [micro_sub_wafer[j-8], micro_sub_wafer[j-7],micro_sub_wafer[j-6], micro_sub_wafer[j]]
#             x = np.linspace(j-8,j, 4)
#             cs = CubicSpline(x, y)
#             xnew = np.linspace(j-8,j, 9)
#             micro_sub_wafer[j-8:j+1] = cs(xnew)
#         else:
#             y = [micro_sub_wafer[j-4],micro_sub_wafer[j-2], micro_sub_wafer[j]]  # 봉우리로부터 3개 깂
#             x = np.linspace(j-4,j, 3)
#             cs = CubicSpline(x, y)
#             xnew = np.linspace(j-4,j, 5)
#             micro_sub_wafer[j-4:j+1] = cs(xnew)
             
    return pd.Series(micro_sub_wafer)
  
  
  import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)



#결과출력
path = '//10.10.10.113/반도체사업본부/기술개발총괄/시뮬레이션팀/개인/백승혜_양지현/WaferMap/'


interpolated_wafer = np.zeros((300,300))
mini_wafer_map = np.zeros((63, 63))
arc_divided = 16

total = pd.DataFrame(columns = ['x', 'y', 'value'])

for i in range(arc_divided):
    macro_sub_wafer, xy_before = choose_radian_point(raw_wafer, i, arc_divided)
    pd_macro_sub_wafer = pd.DataFrame(macro_sub_wafer)
    macro_sub_wafer = pd_macro_sub_wafer.stack()
    # #     temp_macro_wafer = pd.DataFrame(pd_macro_sub_wafer.iloc[::-1][0])
    # #     macro_sub_wafer = ((pd_macro_sub_wafer + temp_macro_wafer)/2).stack()
    micro_sub_wafer = apply_micro_variation(macro_sub_wafer)
    mini_wafer_map = match_coordinate(micro_sub_wafer, i, arc_divided,xy_before)
    total= total.append(mini_wafer_map)
    total.to_csv(path + 'result_ex3.csv')
    print(total)
    
    
    
  def visualize_wafer(mini_wafer_map):
    
    # Interpolation 실행
    data = mini_wafer_map
    x = data['x'].tolist()
    y = data['y'].tolist()
    value = data['value'].tolist()
    
    xx = np.array(x)
    yy = np.array(y)
    zz = np.array(value)
    
    npts = len(xx)            # 현재 PT 개수 (490개)
    ngridx, ngridy = 300,300  # 보간하고자 하는 PT 개수
    
    # Create grid values first.
    xi = np.linspace(-137, 137, ngridx)
    yi = np.linspace(-137, 137, ngridy)
    
    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(xx, yy)
    interpolator = tri.LinearTriInterpolator(triang, zz)  # LinearTriInterpolator
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)
    
    fig = plt.figure(figsize = (10,10))
    ax = plt.axes(projection='3d')
    ax.axes.set_zlim3d(bottom=40000, top=50000) 
    contour = ax.contour3D(Xi, Yi, zi,levels = 100, cmap='jet')
    ax.view_init(35, 50)  # 보는 각도 변경

    fig.colorbar(contour, shrink=0.5, aspect=5)
    
    fig.savefig("3d graph.png")
  
  
  
  def visualize_wafer_2d(mini_wafer_map):
    data = mini_wafer_map
    X = data['x'].tolist()
    Y = data['y'].tolist()
    value = data['value'].tolist()

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(value)
    
    npts = len(X)              # 현재 PT 개수
    ngridx, ngridy = 300,300  # 보간하고자 하는 PT 개수
    lv = 10
    
    fig, ax = plt.subplots(nrows=1)
    
    # Create grid values first.
    xi = np.linspace(-135, 135, ngridx)
    yi = np.linspace(-135, 135, ngridy)

    # Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
    triang = tri.Triangulation(X, Y)
    interpolator = tri.CubicTriInterpolator(triang, Z)  # LinearTriInterpolator
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    s = ax.contour(xi, yi, zi, levels=lv, linewidths=0.5, colors = 'k')
    colors = ['#990033','#CC0033','#FF3300','#FF6600','#FFFF00','#FFFF66','#FFFF99','#FFFFCC','#669933','#339900','#33CC00','#33FF00','#003399','#0033CC','#0033FF','#0066FF','#663399','#663399','#9900CC','#9933FF','#9900FF']
    colors = list(reversed(colors))
    cntr1 = ax.contourf(xi, yi, zi, levels=lv, cmap=matplotlib.colors.ListedColormap(colors))
    fig.colorbar(cntr1, shrink=0.8, aspect=5)   # legend
    fig.savefig('2d graph.png')
