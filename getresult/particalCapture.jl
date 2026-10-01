#==========初始化算法==========#
include("../methods/MovingMesh/movingMesh.jl")#加载弹簧法
include("../getresult/comper.jl")#加载弹簧法 
include("../methods/generateMesh/GenerateMesh.jl")

#==========数据初始化==========#
meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)
freeMeshLabel=1#自由变形域
MoveBDYIndex=[2]#指定运动边界
MoveBDYPoint=unique(vcat(sepBDYNodes[MoveBDYIndex]...))#运动边界的点
fixedBDYIndex=filter(i->!in(i,MoveBDYIndex),1:length(sepBDYNodes))
fixedBDYPoint=unique(vcat(sepBDYNodes[fixedBDYIndex]...))

showWholeMesh(meshPoint,meshElement,color="#5C565C",cover=true,linewidth=0.3,enlarge=[-9.0,35.0,-4.0,34.0],fontsize=12)
DrawPoint(MoveBDYPoint)

include("../DataLoading.jl")#加载路径文件
particalMove,traceLine=gettrace()
N=size(particalMove,1)#计算时间步数量
#=======单步平移调试==========#
    minimum([TriangleQuality(tri) for tri in meshElement])<0.15
    m=29
    movingPara=[particalMove[m,2],particalMove[m,3]]
    # movingPara=[-11,17]
    showWholeMesh(meshPoint,meshElement,color="#9090C5",enlarge=[-26,10.0,-5.0,15.0])
    Points,Meshs=wholeMovingMesh(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=true)#第一次运动

    # Points,Meshs, fixedBDYPoint2, MoveBDYPoint2 = generate_mesh(GeoPolygon)

    countFlipTri(Points,Meshs)
    minimum([TriangleQuality(tri) for tri in Meshs])
    showWholeMesh(Points,Meshs,color="#9090C5",enlarge=[-26,10.0,-5.0,15.0])
    
    #对比扭转弹簧
    Points,Meshs=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=true)
    countFlipTri(Points,Meshs)
    showWholeMesh(Points,Meshs,color="#609DBF",enlarge=[-6,25.0,3.0,44.0])
#========连续运动最大距离测试==========#
    #🦄J-SPA算法
    p = Progress(N, desc="Moving mesh: ")  
    @showprogress for m in 2:N
        meshPoint_init,meshElement_init,u_init,K1_Init,K2_Init,K3_Init=MovingMesh_init(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,[particalMove[m,2],particalMove[m,3]],Parallel=true)#预变形
        initParams=paramsInit(meshElement,meshElement_init,meshPoint_init,u_init,K1_Init,K2_Init,K3_Init)#初始化参量字典构建
        meshPoint,meshElement=MovingMesh(meshPoint,meshElement,MoveBDYPoint,fixedBDYPoint,[particalMove[m,2],particalMove[m,3]],initParams,Parallel=true)
        showWholeMesh(meshPoint,meshElement,color="#9090C5",linewidth=0.5,fontsize=15)
        sleep(1)
        if minimum([TriangleQuality(tri) for tri in meshElement])<0.15
            @info " - 运动 $m 步"
            break
        end
        # if countFlipTri(meshPoint,meshElement)>0#判断是否出现翻转
        #     @info " - 第 $m 步"
        #     @error "出现翻转网格"
        #     break
        # end
    end# 26%(0.129)

    #💫对比扭转弹簧
    p = Progress(N, desc="Moving mesh: ")  
    @showprogress for m in 2:N
        meshPoint,meshElement=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,[particalMove[m,2],particalMove[m,3]],Parallel=true)
        showWholeMesh(meshPoint,meshElement,color="#609DBF",linewidth=0.5,fontsize=15)
        sleep(1)
        if minimum([TriangleQuality(tri) for tri in meshElement])<0.2
            @info " - 运动 $m 步"
            break
        end
    end#24%(0.11)
#=======================连续运动带网格重构测试===============#
    formMesh=[]
    formMeshPoints=[]
    #🦄J-SPA算法
    m=2
    i=0
    p = Progress(N, desc="Moving mesh: ")  
    Gmsh.initialize()
    while m<N+1#当m大于时间步数量时停止循环
        lastGeoPolygon=deepcopy(GeoPolygon)#储存停止前的最后一个步骤
        # if m == 29
        #     break
        # end
        #====网格变形=====#
        meshPoint,meshElement=wholeMovingMesh(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,[particalMove[m,2],particalMove[m,3]],Parallel=true)#第一次运动
        for p in GeoPolygon[1].Points
            p.x=p.x+particalMove[m,2]
            p.y=p.y+particalMove[m,3]
        end
        #=========质量判断=========#
        if countFlipTri(meshPoint,meshElement)>0#判断是否出现翻转
            showWholeMesh(meshPoint,meshElement,color="#9090C5",linewidth=0.5,fontsize=15)
            @info " - 第 $m 步"
            @error "出现翻转网格"
            break
        end
        if minimum([TriangleQuality(tri) for tri in meshElement])<0.2#判断是否重构
            i=i+1
            @info " - 第 $m 步 - 第 $i 次重构"
            meshPoint,meshElement, fixedBDYPoint, MoveBDYPoint = generate_mesh(GeoPolygon,min_size=0.2,max_size=2.5) #对新几何位置进行重构

            # GeoPolygon=deepcopy(lastGeoPolygon)#GeoPolygon恢复到上一个

            push!(formMesh,[meshElement,particalMove[m,1],2])
            push!(formMeshPoints,[meshPoint,particalMove[m,1],2])
            m=m+1
        else#不需要重构，收集当前运动结构
            push!(formMesh,[meshElement,particalMove[m,1],1])
            push!(formMeshPoints,[meshPoint,particalMove[m,1],1])
            m=m+1
        end
        
        showWholeMesh(meshPoint,meshElement,color="#9090C5",linewidth=0.5,fontsize=15)
        sleep(1)
        next!(p)  # 手动更新进度条
    end# 26%(0.129)
    
#=======================可视化与测试========================#
    #检查是否有翻转三角形
     countFlipTri(Points,Meshs)
    #计算所有单元质量最小值
    minimum([TriangleQuality(tri) for tri in meshElement])
    #绘图
    Color=["#5C565C","#9090C5","#609DBF"]
    showWholeMesh(Points,Meshs,color=Color[3],linewidth=0.5,fontsize=15)
    #局部放大
    showWholeMesh(Points,Meshs,color=Color[2],linewidth=1.0,fontsize=15)
    newLine=[Line(Points[l.PID[1]],Points[l.PID[2]]) for l in sepBDYLines[end]]
    DrawLine(sepBDYLines[end])#原边界位置
    DrawLine(newLine,linewidth=5.0,color=Color[2])#旋转后边界
    enlarge=[-9.0,35.0,-4.0,34.0]
    plt.xlim(enlarge[1], enlarge[2])
    plt.ylim(enlarge[3], enlarge[4])

    #连续运动后局部放大
    showWholeMesh(meshPoint,meshElement,color=Color[2],linewidth=1.0,fontsize=15)
    newLine=[Line(meshPoint[l.PID[1]],meshPoint[l.PID[2]]) for l in sepBDYLines[end]]
    DrawLine(sepBDYLines[end])#原边界位置
    DrawLine(newLine,linewidth=5.0,color=Color[2])#旋转后边界
    enlarge=[-9.0,35.0,-4.0,34.0]
    plt.xlim(enlarge[1], enlarge[2])
    plt.ylim(enlarge[3], enlarge[4])

    

#=========输出comsol使用的轨迹数据============#
    outpath="./result"

    traceLine_X=writedlm(string(outpath, "/traceLine_X.txt"), particalMove[:,1:2])
    traceLine_Y=writedlm(string(outpath, "/traceLine_Y.txt"), particalMove[:,collect([1,3])])

    maximum([p.x for p in MoveBDYPoint])
    maximum([p.y for p in MoveBDYPoint])
#==========输出完整变形构型==================#
    outpath="./result/"
    #输出变形构型
    jldsave(string(outpath, "formMesh.jld2"); formMesh)
    jldsave(string(outpath, "formMeshPoints.jld2"); formMeshPoints)
#========读取网格变形构型=====#
    outpath="./result/particalCapture/"
    formMesh = load(string(outpath, "formMesh.jld2"), "formMesh")
    formMeshPoints = load(string(outpath, "formMeshPoints.jld2"), "formMeshPoints")

    #
    timeStep=length(formMesh)
    p = Progress(timeStep, desc="Moving mesh: ")  
    
    #动态显示网格
    @showprogress for t in 1:timeStep
        ps,ms,label=formMeshPoints[t][1],formMesh[t][1],formMesh[t][end]
        if label==1
            showWholeMesh(ps,ms,color="#5C565C",cover=true,linewidth=0.3)
        else
            showWholeMesh(ps,ms,color="#609DBF",cover=true,linewidth=0.3)
        end
        sleep(1)
    end
    #最小单元质量统计
    minimumQuality=[]#记录最小单元质量
    @showprogress for t in 1:timeStep
        ps,ms,label=formMeshPoints[t][1],formMesh[t][1],formMesh[t][end]
        push!(minimumQuality,[minimum([TriangleQuality(tri) for tri in ms]),formMesh[t][2]])
    end
    formMesh[1][2]
    QualityCurve=Line[]
    for (i,info) in enumerate(minimumQuality) 
        if i==timeStep
            break
        else
            push!(QualityCurve,Line(Point(info[2],info[1]),Point(minimumQuality[i+1][2],minimumQuality[i+1][1])))
        end
    end
    DrawLine(QualityCurve)
#=========输出至matlab==========#
    particalQuality=hcat(hcat(minimumQuality...)'[:,2],hcat(minimumQuality...)'[:,1])#获取按时间节点的质量变化曲线
    outpath="./result/particalCapture"
    writedlm(string(outpath, "/particalQuality.txt"), particalQuality)
#===========按帧数提取数据==========#
    pos2="./result/particalCapture/diedai"
    frame=6#帧数
    gap=floor(Int,timeStep/(frame-1))#取值步长
    index=vcat([1,[1+index*gap for index in 1:frame-2]...,timeStep[end]])#索引

    InterMesh=[m[1] for m in formMesh[index]]
    InterPoint=[m[1] for m in formMeshPoints[index]]
    InterTimeStep=[info[2] for info in formMesh[index]]
    
    q=[TriangleQuality(tri) for tri in InterMesh[1]]
    
    showTriQuality(InterPoint[1],InterMesh[1],q,figerSize=[10.1,4.17],enlarge=[-51.8,49.2,-7.03,38.16],figName="Mesh Quality")
    for i in 1:frame
        q=[TriangleQuality(tri) for tri in InterMesh[i]]
        showTriQuality(InterPoint[i],InterMesh[i],q,figerSize=[10.1,4.17],enlarge=[-51.8,49.2,-7.03,38.16],figName="Mesh Quality")
        plt.savefig(string(pos2, "/时间-",i,".png"),dpi = 300)
        sleep(1)
    end
#=========运动总距离计算===========#
    timeend_laplace=0.46
    timeend_Winslow=0.2473
    timeend_Hyperelastic=0.11
    timeend_Yeoh=0.15
    particalMove
    movestep=findall(s->s<=timeend_laplace,particalMove[:,1])
    sum([sqrt(l[2]^2+l[3]^2) for l in eachrow(particalMove[movestep[2:end],:])])