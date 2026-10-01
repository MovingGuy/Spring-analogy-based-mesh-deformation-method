#==========================💦💦版本说明💦💦==============================#
# 运动网格算法-变形网格算法优化（阶段性稳定版本）
# 版本号：v9.12
# 版本特点：不同的模型需要使用不同的参数来适应三种刚度的贡献，可实现较大范围的旋转和位移
# 注意：运行代码时，需要根据模型加载./MovingMesh-backup中对应的代码，其中微粒运动代码需要加载STL几何信息，其余只加载网格信息即可
#==========================================================================#
include("./AddIns.jl")#加载所有相关文件
##💧读取初始网格
include("./methods/readMesh/readMesh.jl")#加载COMSOL网格读取函数
include("./methods/getGeoData/BDYIdentification.jl")#加载边界识别可视化函数
path = "./model/particalCapture"
completePoints,completeMesh,OriginBDY=readInitMesh_direct(path)

# showWholeMesh(completePoints,completeMesh)

#区分几何模型中的多边形
sepBDYLines=getsepBDYLines(OriginBDY)
#取出边界
# ##💧读取几何STL信息(仅微粒运动需要运行这段)
include("./methods/getGeoData/GetGeoInfo.jl")
path2="./modelused"
stl_files = [replace(joinpath(path2, file), "\\" => "/") for file in readdir(path2) if occursin(r"(?i)\.stl$", file)]
mesh = load(stl_files[1])   
GeoPolygon,domainInfo=getDomainBDY(path2);
#==========💦统一管理边界信息============#
partMesh,partPoints=deepcopy(completeMesh),deepcopy(completePoints)
include("./BoundCo.jl")#注：该过程的目的是方便实现类似comsol的效果，用编号进行边界条件的赋值
RegionIndex,BDYs,sepBDYNodes=BoundCo(partMesh,sepBDYLines)

DrawSepBdy(sepBDYLines,showLineIndex=true,color="",linewidth=3.0,linestyle="--")#此时边界与边界节点分组对应，可用line图查看，进行边界设置

#==================💧变形条件设置=======================#
# freeMeshLabel=1#自由变形域
# MoveBDYIndex=[2]#指定运动边界
# MoveBDYPoint=unique(vcat(sepBDYNodes[MoveBDYIndex]...))#运动边界的点
# fixedBDYIndex=filter(i->!in(i,MoveBDYIndex),1:length(sepBDYNodes))
# fixedBDYPoint=unique(vcat(sepBDYNodes[fixedBDYIndex]...))
##
#=============💧变形算法(初始化)=============#
# include("./methods/MovingMesh/movingMesh.jl")#加载弹簧法
# include("./methods/MovingMesh/backup/particalCapturebackup.jl")
# include("./getresult/comper.jl")#加载弹簧法 
# meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)

#=============微粒运动模型====================#
# include("./DataLoading.jl")#加载所有相关文件
# particalMove,traceLine=gettrace()
# N=size(particalMove,1)#计算时间步数量
#=============P30N30机翼旋转========================#
# translation_steps=compute_translation_steps(MoveBDYPoint,Point(1650.0,320.0),-1.1,1.1)#计算P30N30旋转步长
#=============NACA0012机翼旋转========================#
# translation_steps=compute_translation_steps(MoveBDYPoint,Point(0.0,0.0),7.9,7.9)#计算NACA0012旋转步长


#=================绘图==================#

