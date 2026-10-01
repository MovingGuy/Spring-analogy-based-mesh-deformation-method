"""
初始位移量计算
"""
function MovingMesh_Torsion(Points::Vector{Point},Mesh::Vector{MeshUnit},freeMeshLabel,MoveBDYPoint,fixedBDYPoint,movingPara;Parallel=true)
    DeformMesh=filter(t->t.did==freeMeshLabel,Mesh)#取出变形域
    ##计算网格节点运动
    #计算刚度矩阵
    K1=AssembleKForLinear_comper(DeformMesh,Points)
    K2=AssembleKForTorsion_comper(DeformMesh)

    K=K1+K2
    F=zeros(2*length(Points))#2N个F
    movenode =[point.id for point in MoveBDYPoint]#获取节点id
    outnode =[point.id for point in fixedBDYPoint]#获取节点id

    K,F=AddEssBDYForPanMoving(K,F,outnode,[0,0])
    if Parallel
        K,F=AddEssBDYForPanMoving(K,F,movenode,movingPara)#平移运动边界
    else
        K,F=AddEssBDYForMovingVariable(K,F,movenode,movingPara)#边界任意移动
    end
    u_init=Solve(K,F)
    ##更新网格节点位置
    PS=deepcopy(Points)
    MS=deepcopy(Mesh)
    for p in PS
        indices = findfirst(x -> x == p, PS)
        if indices !== nothing
            p.x=p.x+u_init[2*indices-1]
            p.y=p.y+u_init[2*indices]
        end
    end
    #更新网格中所有节点,不检测是否有交叉节点存在
    for ele in MS
        ele.p1=PS[ele.p1.id]
        ele.p2=PS[ele.p2.id]
        ele.p3=PS[ele.p3.id]
    end
    return PS,MS
end
"""
初始化线弹簧刚度
"""
function AssembleKForLinear_comper(Mesh::Vector{MeshUnit},Points::Vector{Point})
    IK=Vector{Int64}()
    JK=Vector{Int64}()
    VK=Vector{Float64}()#稀疏矩阵存储方法
    #取出所有边
    collectLine = Set{Line}()#采用set减少in判断的消耗
    for tri in Mesh 
        push!(collectLine, Line(tri.p1, tri.p2))
        push!(collectLine, Line(tri.p2, tri.p3))
        push!(collectLine, Line(tri.p3, tri.p1))
    end
    
    #按边计算形变
    for line in collectLine
        DOF=line.PID#记录节点在求解域中的序号
        cell=Points[DOF] 
        theK=calculateKForLinear_comper(cell)
        IK,JK,VK=AssembleIJBVK(IK,JK,VK,theK,DOF) 
    end
    KG=sparse(IK,JK,VK)#组装稀疏矩阵，VK中是非零项的值，IK，JK是非零项的J,K坐标
    return KG
end
"""
计算以线为单位的刚度
"""
function calculateKForLinear_comper(seg::Vector{Point})
    l=(1/(getLength(seg[1],seg[2])))#计算线段刚度，将线段长度的倒数作为刚度
    #对二自由度的线段计算单根弹簧刚度矩阵
    K=l * [1 0 -1 0; 0 1 0 -1; -1 0 1 0; 0 -1 0 1]#4*4矩阵
    T=calculateTransformationMatrix_comper(seg)#计算变换矩阵，将轴向坐标映射到全局坐标
    return T' * K * T#刚度矩阵映射  
end
"""
计算线弹簧局部坐标映射矩阵
"""
function calculateTransformationMatrix_comper(seg::Vector{Point})
    # 获取节点坐标
    x1, y1 = seg[1].x, seg[1].y
    x2, y2 = seg[2].x, seg[2].y 

    # 计算边的方向余弦
    L1 = getLength(seg[1], seg[2])#计算长度

    cosα = (x2 - x1) / L1#cos
    sinα = (y2 - y1) / L1#sin
 
    # 构造变换矩阵，该矩阵实际上是T的逆矩阵，[U1x，U1y，U2x，U2y]
    return  [
        cosα sinα     0      0;
       -sinα cosα     0      0;
           0      0    cosα sinα;
           0      0   -sinα cosα;
    ]
end
#=============扭转弹簧部分==============#
#方式二：采用按矩阵组装的方式进行扭转刚度矩阵的计算
"""
扭转弹簧刚度矩阵初始化
"""
function AssembleKForTorsion_comper(Mesh::Vector{MeshUnit})
    IK=Vector{Int64}()
    JK=Vector{Int64}()
    VK=Vector{Float64}()#稀疏矩阵存储方法
    #取出所有边
    for tri in Mesh
        DOF=tri.PID 
        #组装K矩阵，使其可以与线弹簧相加
        TorsionK=calculateKForTorsion_comper(tri)
        IK,JK,VK=AssembleIJBVK(IK,JK,VK,TorsionK,DOF) 
    end
    KG=sparse(IK,JK,VK)#组装稀疏矩阵，VK中是非零项的值，IK，JK是非零项的J,K坐标
    return KG
end
"""
计算单个三角形网格的扭转刚度矩阵
- 输入：tri::MeshUnit参与计算的三角形,localPoints::Vector{Point}节点集
- 输出：K扭转弹簧的刚度矩阵6*6
"""
function calculateKForTorsion_comper(tri::MeshUnit)
    TriPoints=[tri.p1,tri.p2,tri.p3]
    Params=[]
    for point in TriPoints#遍历三角形的三个点i,j,k
        otherPs=filter(p->p!=point,TriPoints)#找到除当前点的另外两个点
        append!(Params,[p-point for p in otherPs])#获得顺序为ij,ik,ji,jk,ki,kj
    end 
    edgeLength=[getLength(TriPoints[1],TriPoints[2]),getLength(TriPoints[2],TriPoints[3]),getLength(TriPoints[3],TriPoints[1])]#计算三个边的边长 
    Rc=[edgeLength[1],edgeLength[3],edgeLength[1],edgeLength[2],edgeLength[3],edgeLength[2]].^2
    #形成旋转增量矩阵
    a=[p.x for p in Params]
    b=[p.y for p in Params]
    a=a./Rc
    b=b./Rc
    Rijk=[b[2]-b[1] a[1]-a[2] b[1] -a[1] -b[2] a[2];
    -b[3] a[3] b[3]-b[4] a[4]-a[3] b[4] -a[4];
    b[5] -a[5] -b[6] a[6] b[6]-b[5] a[5]-a[6]]
    
    C=calculateTorStiff_comper(edgeLength,TriPoints)
    K=Rijk'*C*Rijk
    return K
end
"""
计算三角形三点扭转刚度
输入：edgeLength三角形三边边长
"""
function calculateTorStiff_comper(edgeLength,localPoints::Vector{Point})
    l1,l2,l3=edgeLength[1],edgeLength[2],edgeLength[3]
    S=calculateTriS(localPoints[1],localPoints[2],localPoints[3])^2
    #计算扭转弹簧刚度
    Ci=(l1^2+l3^2)/(4*S)
    Cj=(l1^2+l2^2)/(4*S)
    Ck=(l2^2+l3^2)/(4*S)
    k_torsion = 1  # 扭转系数
    C=k_torsion*[Ci 0 0;0 Cj 0;0 0 Ck]
    return C
end


#=========矩阵组装==========#
#矩阵组装
function AssembleIJBVK(IK,JK,VK,ke,DOF)
    num=size(ke)
    for i=1:num[1]
        for  j=1:num[2]
        ii=(i - 1) ÷ 2 + 1#计算节点在矩阵中的位置
        jj=(j - 1) ÷ 2 + 1
        yux=rem(i,2)#区分x,y的位置
        yuy=rem(j,2) 
        IK=append!(IK,DOF[ii]*2-yux) #在末尾增加DOF
        JK=append!(JK,DOF[jj]*2-yuy)
        VK=append!(VK,ke[i,j])
        end
    end
    return IK,JK,VK
end