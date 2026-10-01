#初始网格
meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)
showWholeMesh(meshPoint,meshElement,color="#5C565C",enlarge=[-20,18,-20.0,17.0],linewidth=0.5,fontsize=15)
#局部放大
showWholeMesh(meshPoint,meshElement,color="#5C565C",enlarge=[-3.5,2.2,-3.2,2.2],linewidth=0.5,fontsize=15)

#最佳性能参数：
#线刚度用长度的倒数不是平方
#扭转弹簧参考角度：90，扭转系数150
#面积能量base系数120
#=======单步旋转测试==========#
translation_steps=compute_translation_steps(MoveBDYPoint,Point(0.0,0.0),2.1,2.1)#计算NACA0012旋转步长(目前最远1.88)
movingPara=translation_steps[1]
#🦄J-SPA算法
meshPoint_init,meshElement_init,u_init,K1_Init,K2_Init,K3_Init=MovingMesh_init(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=false)#预变形
#初始化参量字典构建
initParams=paramsInit(meshElement,meshElement_init,u_init,K1_Init,K2_Init,K3_Init)
Points,Meshs=MovingMesh(meshPoint,meshElement,MoveBDYPoint,fixedBDYPoint,movingPara,initParams,Parallel=false)

countFlipTri(Points,Meshs)

showWholeMesh(Points,Meshs,color=Color[1],enlarge=[-20,18,-20.0,17.0],linewidth=0.5,fontsize=15)
#🩸对比扭转弹簧
Points,Meshs=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=false)
#计算网格质量分布
outpath="./result"
TS_naca0012=GetAllTriQ(Points,Meshs,path=outpath,name="TS_naca0012")#输出网格质量（地址是相对于GetAllTriQ所处文件夹）

JSPA_naca0012=GetAllTriQ(Points,Meshs,path=outpath,name="JSPA_naca0012")#输出网格质量（地址是相对于GetAllTriQ所处文件夹）
#=========单步运动步长连续增加失效网格对比============#
flipNum=[]
for i in 1.6:0.05:2.71#扭转弹簧法
    meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)
    translation_steps=compute_translation_steps(MoveBDYPoint,Point(0.0,0.0),i,i)#计算
    Points,Meshs=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,translation_steps[1],Parallel=false)
    tris=countFlipTri(Points,Meshs)
    push!(flipNum,tris)
end
for i in 1.6:0.05:2.71#J-SPA法
    meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)
    translation_steps=compute_translation_steps(MoveBDYPoint,Point(0.0,0.0),i,i)#计算
    
    meshPoint_init,meshElement_init,u_init,K1_Init,K2_Init,K3_Init=MovingMesh_init(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,translation_steps[1],Parallel=false)#预变形
    #初始化参量字典构建
    initParams=paramsInit(meshElement,meshElement_init,u_init,K1_Init,K2_Init,K3_Init)
    Points,Meshs=MovingMesh(meshPoint,meshElement,MoveBDYPoint,fixedBDYPoint,translation_steps[1],initParams,Parallel=false)

    tris=countFlipTri(Points,Meshs)
    push!(flipNum,tris)
end
#计算网格质量分布
JSPA_naca0012=GetAllTriQ(Points,Meshs,path=outpath,name="JSPA_naca0012")#输出网格质量（地址是相对于GetAllTriQ所处文件夹）
#测试非法网格
tris=countFlipTri(Points,Meshs)
#绘图
Color=["#9090C5","#609DBF"]
showWholeMesh(Points,Meshs,color=Color[1],enlarge=[-20,18,-20.0,17.0],linewidth=0.5,fontsize=15)
#局部放大
showWholeMesh(Points,Meshs,color=Color[1],enlarge=[-3.5,2.2,-3.2,2.2],linewidth=0.5,fontsize=15)
#=========输出至matlab==========#
outpath="./result"
TS_naca0012_flipNum=writedlm(string(outpath, "/TS_naca0012_flipNum.txt"), flipNum)
JSPA_naca0012_flipNum=writedlm(string(outpath, "/JSPA_naca0012_flipNum.txt"), flipNum)