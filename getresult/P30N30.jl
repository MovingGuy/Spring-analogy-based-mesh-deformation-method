#初始网格
meshPoint,meshElement=deepcopy(partPoints),deepcopy(partMesh)
showWholeMesh(meshPoint,meshElement,color="#5C565C",enlarge=[-500.0,2500.0,-1000.0,1500.0],linewidth=0.5,fontsize=15)
#局部放大
showWholeMesh(meshPoint,meshElement,color="#5C565C",enlarge=[1350.0,2250.0,-450.0,450.0],linewidth=0.5,fontsize=15)


#=======单步旋转调试==========#
#🦄J-SPA算法
translation_steps=compute_translation_steps(MoveBDYPoint,Point(1650.0,320.0),-1.38,1.38)#计算P30N30旋转步长
meshPoint_init,meshElement_init,u_init,K1_Init,K2_Init,K3_Init=MovingMesh_init(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,translation_steps[1],Parallel=false)#预变形
initParams=paramsInit(meshElement,meshElement_init,u_init,K1_Init,K2_Init,K3_Init)#初始化参量字典构建
Points,Meshs=MovingMesh(meshPoint,meshElement,MoveBDYPoint,fixedBDYPoint,translation_steps[1],initParams,Parallel=false)

#🪞对比扭转弹簧
translation_steps=compute_translation_steps(MoveBDYPoint,Point(1650.0,320.0),-0.36,0.36)#“计算P30N30旋转步长
Points,Meshs=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,translation_steps[1],Parallel=false)

#=========连续运动最大弧度测试===========#
#🦄J-SPA算法
translation_steps=compute_translation_steps(MoveBDYPoint,Point(1650.0,320.0),-1.8,0.01)
p = Progress(length(translation_steps), desc="Moving mesh: ")  
@showprogress for (i,movingPara) in enumerate(translation_steps)
    meshPoint_init,meshElement_init,u_init,K1_Init,K2_Init,K3_Init=MovingMesh_init(meshPoint,meshElement,freeMeshLabel,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=false)#预变形
    initParams=paramsInit(meshElement,meshElement_init,meshPoint_init,u_init,K1_Init,K2_Init,K3_Init)#初始化参量字典构建
    meshPoint,meshElement=MovingMesh(meshPoint,meshElement,MoveBDYPoint,fixedBDYPoint,movingPara,initParams,Parallel=false)
    # showWholeMesh(meshPoint,meshElement,color="#5C565C",enlarge=[-500.0,2500.0,-1000.0,1500.0],linewidth=0.5,fontsize=15)
    # sleep(1)
    if countFlipTri(meshPoint,meshElement)!=0
        println(i)
        break
    end
end#1.9rad(64%)
#💫对比扭转弹簧
translation_steps=compute_translation_steps(MoveBDYPoint,Point(1650.0,320.0),-2.96,0.01)
p = Progress(length(translation_steps), desc="Moving mesh: ")  
@showprogress for (i,movingPara) in enumerate(translation_steps)
    meshPoint,meshElement=MovingMesh_Torsion(meshPoint,meshElement,1,MoveBDYPoint,fixedBDYPoint,movingPara,Parallel=false)
    if countFlipTri(meshPoint,meshElement)!=0
        println(i)
        break
    end
end#0.94rad(32%)
#=======================可视化与测试========================#
Points=deepcopy(meshPoint)
Meshs=deepcopy(meshElement)
#检查是否有翻转三角形
countFlipTri(Points,Meshs)
#绘图
Color=["#5C565C","#9090C5","#609DBF"]
showWholeMesh(Points,Meshs,color=Color[2],enlarge=[-500.0,2500.0,-1000.0,1500.0],linewidth=0.5,fontsize=15)
#局部放大
plt.figure(figsize=(10.08,6.33))
showWholeMesh(Points,Meshs,color=Color[2],linewidth=0.5,fontsize=22,showticks=false)
newLine=[Line(Points[l.PID[1]],Points[l.PID[2]]) for l in sepBDYLines[end]]
DrawLine(sepBDYLines[end],linewidth=5.0)#原边界位置
DrawLine(newLine,linewidth=5.0,color=Color[2])#旋转后边界
enlarge=[1085.25,2333.90,-300.93,424.158]
plt.xlim(enlarge[1], enlarge[2])
plt.ylim(enlarge[3], enlarge[4])

q=[TriangleQuality(tri) for tri in Meshs]
showTriQuality(Points,Meshs,q,figerSize=[10.08,6.33],enlarge=[1085.25,2333.90,-300.93,424.158],fontsize=22,showticks=false)
#===========输出运动后网格结果=================#
outpath="./result/P30N30/"
#输出变形构型
jldsave(string(outpath, "formMesh.jld2"); meshElement)
jldsave(string(outpath, "formMeshPoints.jld2"); meshPoint)